//! **The self-initialising quadratic sieve** (`big` feature): Contini's method, the one Dario
//! Alpern's page runs, written over the window.
//!
//! One polynomial Q(x) = x² − n grows with x, so most of a wide interval is too big to be
//! smooth. Many polynomials Q_ab(x) = (a·x + b)² − n, each over a short interval [−M, M], keep
//! every value near M·√(n/2). With a = q₁⋯q_s a product of factor-base primes and b² ≡ n (mod a)
//! there are 2^(s−1) choices of b for one a, and stepping through them in Gray-code order moves
//! every prime's two sieve roots by a precomputed amount: that is the "self-initialising".
//! Everything else is the window: `ln p` along two strides per prime into a byte array, a
//! threshold, exact division only at the hits, one large prime kept and merged, a GF(2)
//! dependency, gcd. Threads take different a's. Read next to `qs.rs`; the parts have the
//! same names.

use num_bigint::{BigInt, BigUint};
use num_integer::Integer;
use num_traits::{One, Signed, ToPrimitive, Zero};
use rayon::prelude::*;

use crate::qs::{dependencies_sparse, legendre, sqrt_mod};
use crate::smooth::primes_upto;
use std::sync::atomic::{AtomicU64, Ordering::Relaxed};

/// Nanoseconds spent in each stage, summed over threads, and the candidates seen: the sieve
/// fill, the scan for hits, the exact division of hits, the polynomial switches. Read by
/// `examples/siqs_sweep.rs`; reset with `stats_reset`.
pub static STATS: [AtomicU64; 6] = [AtomicU64::new(0), AtomicU64::new(0), AtomicU64::new(0), AtomicU64::new(0), AtomicU64::new(0), AtomicU64::new(0)];
pub const STAT_NAMES: [&str; 6] = ["sieve fill", "scan", "exact division", "switch", "candidates", "bucket write"];
pub fn stats_reset() { for s in &STATS { s.store(0, Relaxed); } }
/// The clock: `std::time::Instant` natively; in a browser (`wasm32-unknown-unknown`) there is
/// no monotonic clock in `std`, so the stage counters and the progress estimate read zero and
/// the sieve itself runs unchanged.
pub mod clock {
    #[cfg(not(target_arch = "wasm32"))]
    pub use std::time::Instant;
    #[cfg(target_arch = "wasm32")]
    #[derive(Clone, Copy, Debug)]
    pub struct Instant;
    #[cfg(target_arch = "wasm32")]
    impl Instant {
        pub fn now() -> Instant { Instant }
        pub fn elapsed(&self) -> std::time::Duration { std::time::Duration::ZERO }
        pub fn duration_since(&self, _earlier: Instant) -> std::time::Duration { std::time::Duration::ZERO }
    }
}

fn tick(i: usize, t: clock::Instant) { STATS[i].fetch_add(t.elapsed().as_nanos() as u64, Relaxed); }

/// A factor-base prime with the root of n, its log in the byte scale, and the pair for the
/// Granlund–Montgomery exact-divisibility test: p | d iff d·inv (mod 2³²) ≤ limit, no division.
#[derive(Clone, Debug)]
pub struct Prime { pub p: u64, pub root: u64, pub log: u8, pub inv: u32, pub limit: u32 }

impl Prime {
    fn new(p: u64, root: u64) -> Prime {
        let inv = if p % 2 == 1 { inverse_mod_2_32(p as u32) } else { 0 };
        Prime { p, root, log: ((p as f64).log2() + 0.5) as u8, inv, limit: (u32::MAX / p as u32) }
    }
    /// Does p divide d? For odd p by the multiply-and-compare test, for 2 by the bit.
    #[inline]
    pub fn divides(&self, d: u32) -> bool { if self.p == 2 { d & 1 == 0 } else { d.wrapping_mul(self.inv) <= self.limit } }
}

/// The inverse of an odd number modulo 2³², by Newton's iteration.
fn inverse_mod_2_32(p: u32) -> u32 {
    let mut x = p; // correct to 3 bits for odd p
    for _ in 0..5 { x = x.wrapping_mul(2u32.wrapping_sub(p.wrapping_mul(x))); }
    x
}

/// Primes below this are left out of the sieve fill (their strides touch every few cells and
/// cost more than they tell); the threshold is lowered by their expected contribution instead.
/// `SIQS_SMALL` in the environment overrides it, for the sweeps.
/// The sieve block, in cells. Measured on the 60-digit input: the whole 2^17 interval at once
/// is fastest (1.28 s), 32 KiB blocks slower (1.39 s), 8 KiB much slower (2.2 s), so the fill is
/// not cache-latency bound and the default is one block; the blocks stay for a factor base that
/// outgrows one. `SIQS_SEGMENT` in the environment overrides it, for the sweeps.
pub fn segment_len(bucket: bool, primes: usize) -> usize {
    // a bucket entry packs the prime index above the position in the block, so a base past
    // 131,072 primes takes the 16 KB block
    let seg: usize = std::env::var("SIQS_SEGMENT").ok().and_then(|v| v.parse().ok()).unwrap_or(if !bucket { 1 << 20 } else if primes < 1 << 17 { 1 << 15 } else { 1 << 14 });
    seg.next_power_of_two()
}

/// The bucket sieve for primes wider than a block (on unless `SIQS_BUCKET=0`); with it the
/// default block is 32 KB, the L1 size, instead of the whole interval.
/// How many polynomials the large primes are walked through at once (`SIQS_BATCH`, default 16).
pub fn bucket_batch() -> usize { std::env::var("SIQS_BATCH").ok().and_then(|v| v.parse().ok()).filter(|&b| b >= 1).unwrap_or(16) }

/// The bucket sieve for primes wider than a block: `SIQS_BUCKET=1` forces it on, `0` off, and
/// unset leaves it to the rule in the sieve (on once the largest prime is wider than the
/// interval). With it the block is 32 KB, the L1 size, instead of the whole interval.
pub fn bucket_sieve() -> Option<bool> { std::env::var("SIQS_BUCKET").ok().map(|v| v != "0") }

/// Two large primes on (the double-large-prime variation), unless `SIQS_DLP=0`.
pub fn double_large() -> bool { std::env::var("SIQS_DLP").map(|v| v != "0").unwrap_or(true) }

/// Miller–Rabin for u64, deterministic with the first twelve prime bases.
pub fn is_prime_u64(n: u64) -> bool {
    if n < 2 { return false; }
    for p in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] { if n % p == 0 { return n == p; } }
    let (mut d, mut s) = (n - 1, 0);
    while d % 2 == 0 { d /= 2; s += 1; }
    'outer: for a in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        let mut x = pow_mod_u64(a, d, n);
        if x == 1 || x == n - 1 { continue; }
        for _ in 1..s { x = mul_mod_u64(x, x, n); if x == n - 1 { continue 'outer; } }
        return false;
    }
    true
}
fn mul_mod_u64(a: u64, b: u64, m: u64) -> u64 { ((a as u128 * b as u128) % m as u128) as u64 }
fn pow_mod_u64(mut b: u64, mut e: u64, m: u64) -> u64 { let mut r = 1u64; b %= m; while e > 0 { if e & 1 == 1 { r = mul_mod_u64(r, b, m); } b = mul_mod_u64(b, b, m); e >>= 1; } r }

/// A nontrivial factor of a composite u64 by Pollard's rho with Brent's cycle finding.
pub fn rho_u64(n: u64) -> Option<u64> {
    if n % 2 == 0 { return Some(2); }
    for c in 1..20u64 {
        let (mut y, mut r, mut q, mut g) = (2u64, 1u64, 1u64, 1u64);
        let (mut x, mut ys) = (y, y);
        let f = |v: u64| (mul_mod_u64(v, v, n) + c) % n;
        while g == 1 {
            x = y;
            for _ in 0..r { y = f(y); }
            let mut k = 0;
            while k < r && g == 1 {
                ys = y;
                for _ in 0..128.min(r - k) { y = f(y); q = mul_mod_u64(q, x.abs_diff(y), n); }
                g = gcd_u64(q, n);
                k += 128;
            }
            r *= 2;
            if r > 1 << 20 { break; }
        }
        if g == n { loop { ys = f(ys); g = gcd_u64(x.abs_diff(ys), n); if g > 1 { break; } } }
        if g > 1 && g < n { return Some(g); }
    }
    None
}
fn gcd_u64(mut a: u64, mut b: u64) -> u64 { while b != 0 { let t = a % b; a = b; b = t; } a }

pub fn small_bound() -> u64 { std::env::var("SIQS_SMALL").ok().and_then(|v| v.parse().ok()).unwrap_or(30) }

/// The knobs: factor-base size, half the sieve interval, large-prime multiple, threshold slack
/// in largest-prime logs (0.7 everywhere: stricter was faster on every input tried, as also
/// found on Alpern's sieve; 1.3 cost a third more and 2.0 four times), relations past the base.
#[derive(Clone, Copy, Debug)]
pub struct Params { pub base: usize, pub half: usize, pub large: u64, pub slack: f32, pub extra: usize }

impl Params {
    /// By the size of n in digits, from the usual tables, measured here where it says so.
    pub fn for_digits(d: usize) -> Params {
        match d {
            0..=30 => Params { base: 150, half: 1 << 13, large: 30, slack: 0.7, extra: 24 },
            31..=40 => Params { base: 400, half: 1 << 14, large: 40, slack: 0.7, extra: 32 },
            41..=50 => Params { base: 1100, half: 1 << 15, large: 50, slack: 0.7, extra: 48 },
            51..=56 => Params { base: 1800, half: 1 << 15, large: 60, slack: 0.7, extra: 64 },   // measured on the 52-digit case
            57..=62 => Params { base: 3000, half: 1 << 16, large: 60, slack: 0.7, extra: 64 },
            63..=67 => Params { base: 4500, half: 1 << 16, large: 80, slack: 0.7, extra: 80 },
            68..=72 => Params { base: 15000, half: 1 << 16, large: 100, slack: 0.7, extra: 96 },   // 69-digit sweep with the bucket: 7,000/2^17 3.8 s to the count, 15,000/2^16 3.2 s
            73..=77 => Params { base: 30000, half: 1 << 17, large: 120, slack: 0.7, extra: 120 },   // 75-digit sweep: 10,000 26 s to the count, 25,000 13 s, 40,000 12 s
            // with the bucket sieve the exact division no longer grows with the base, and the
            // budgeted sweeps of 2026-11-09 moved to much larger bases: 80 digits 20,000/2^18
            // 217 relations/s → 50,000/2^17 1,200/s (41 s to the count against 92); 90 digits
            // 50,000 at 35/s → 80,000 at 141/s; 100 digits 120,000 at 14.6/s → 200,000 at 27/s
            78..=84 => Params { base: 50000, half: 1 << 17, large: 100, slack: 0.7, extra: 150 },
            85..=92 => Params { base: 80000, half: 1 << 17, large: 200, slack: 0.7, extra: 200 },
            93..=97 => Params { base: 120000, half: 1 << 17, large: 250, slack: 0.7, extra: 250 },
            _ => Params { base: 200000, half: 1 << 17, large: 300, slack: 0.7, extra: 300 },
        }
    }
}

/// **The Knuth–Schroeppel multiplier**: the small square-free k for which k·n has the most
/// small primes in its factor base, scored as Σ g(p)·ln p over the primes to 1000 less ½ ln k,
/// with g = 2/(p−1) for a prime where k·n is a square, 1/(p−1) for a prime dividing k, and the
/// usual cases at 2. The sieve then factors k·n; the gcd at the end is taken with n.
pub fn multiplier(n: &BigUint) -> u64 {
    if let Some(k) = std::env::var("SIQS_MULT").ok().and_then(|v| v.parse().ok()) { return k; }
    let candidates = [1u64, 2, 3, 5, 6, 7, 10, 11, 13, 14, 15, 17, 19, 21, 22, 23, 26, 29, 30, 31, 33, 34, 35, 37, 38, 39, 41, 42, 43, 46, 47, 51, 53, 55, 57, 58, 59, 61, 62, 65, 66, 67, 69, 70, 71, 73];
    let primes = primes_upto(1000);
    let mut best = (f64::MIN, 1u64);
    for &k in &candidates {
        let kn = n * k;
        let mut score = -0.5 * (k as f64).ln();
        let m8 = (&kn % 8u32).to_u64().unwrap();
        score += match m8 { 1 => 2.0 * 2f64.ln(), 5 => 2f64.ln(), 3 | 7 => 0.5 * 2f64.ln(), _ => 0.0 };
        for &p in primes.iter().skip(1) {
            let r = (&kn % p).to_u64().unwrap();
            let lp = (p as f64).ln();
            if k % p == 0 { score += lp / (p - 1) as f64; } else if legendre(r, p) == 1 { score += 2.0 * lp / (p - 1) as f64; }
        }
        if score > best.0 { best = (score, k); }
    }
    best.1
}

/// The first `count` primes modulo which n is a square (2 handled as a stride of its own).
pub fn factor_base(n: &BigUint, count: usize) -> Vec<Prime> {
    let mut out = vec![Prime::new(2, (n % 2u32).to_u64().unwrap())];
    let mut bound = 1000u64;
    let mut checked_through = 2u64;
    loop {
        for p in primes_upto(bound).into_iter().skip(1) {
            if out.len() >= count { return out; }
            if p <= checked_through { continue; }
            let np = (n % p).to_u64().unwrap();
            if np == 0 { out.push(Prime::new(p, 0)); } else if legendre(np, p) == 1 { out.push(Prime::new(p, sqrt_mod(np, p))); }
        }
        checked_through = bound;
        bound *= 2;
    }
}

/// One relation: (a·x + b) mod n, the exponents over the base as `(column, exponent)` pairs
/// (column 0 is the sign, column j + 1 the j-th base prime), the large primes outside the base
/// with their exponents (empty for a full relation).
#[derive(Clone, Debug)]
pub struct Relation { pub y: BigUint, pub exps: Vec<(u32, u32)>, pub large: Vec<(u64, u32)> }

impl Relation {
    /// The product of two relations: exponents summed, large primes summed.
    pub fn times(&self, other: &Relation, n: &BigUint) -> Relation {
        let mut exps: Vec<(u32, u32)> = self.exps.clone();
        for &(c, e) in &other.exps { match exps.iter_mut().find(|(d, _)| *d == c) { Some(t) => t.1 += e, None => exps.push((c, e)) } }
        let mut large = self.large.clone();
        for &(p, e) in &other.large { match large.iter_mut().find(|(q, _)| *q == p) { Some(t) => t.1 += e, None => large.push((p, e)) } }
        Relation { y: &self.y * &other.y % n, exps, large }
    }
    /// Is every exponent even, base and large primes alike?
    pub fn is_square(&self) -> bool { self.exps.iter().all(|(_, e)| e % 2 == 0) && self.large.iter().all(|(_, e)| e % 2 == 0) }
}

/// A polynomial family: a = ∏ q_i, the B_l such that b = Σ ±B_l, and a⁻¹ mod p per prime.
struct Family { a: BigUint, qs: Vec<usize>, bs: Vec<BigUint>, ainv: Vec<u64> }

fn modinv(a: u64, p: u64) -> u64 {
    let (mut t, mut nt, mut r, mut nr) = (0i128, 1i128, p as i128, (a % p) as i128);
    while nr != 0 { let q = r / nr; (t, nt) = (nt, t - q * nt); (r, nr) = (nr, r - q * nr); }
    t.rem_euclid(p as i128) as u64
}

/// Pick a's primes: s of them from the middle of the base so that a ≈ √(2n)/M, chosen by a
/// seeded walk so every thread and every round takes a different set.
fn choose_a(n: &BigUint, base: &[Prime], half: usize, seed: u64) -> Family {
    // The target-A repair is on unless `SIQS_CLOSE_A=0`: measured here, 60 digits
    // 1.27 s → 0.74 s over three alternating pairs, 80 digits 121 s → 91 s
    choose_a_mode(n, base, half, seed, std::env::var("SIQS_CLOSE_A").map(|v| v != "0").unwrap_or(true))
}

fn choose_a_mode(n: &BigUint, base: &[Prime], half: usize, seed: u64, close: bool) -> Family {
    let target = (n.to_f64().unwrap() * 2.0).sqrt() / half as f64;
    let target_log = target.ln();
    // candidates: primes around 2000, the usual choice; s from the target
    let lo = base.iter().position(|q| q.p > 1000).unwrap_or(base.len() / 2).min(base.len().saturating_sub(10));
    let hi = base.len().min(lo + 300).max(lo + 10).min(base.len());
    let pool: Vec<usize> = (lo..hi).collect();
    let mut rng = seed.wrapping_mul(0x9E3779B97F4A7C15).wrapping_add(0xD1B54A32D192ED03);
    let mut next = || { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; rng };
    let mut qs: Vec<usize> = Vec::new();
    let mut log_a = 0.0f64;
    while log_a < target_log - (base[pool[pool.len() / 2]].p as f64).ln() / 2.0 && qs.len() < 20 {
        let i = pool[(next() % pool.len() as u64) as usize];
        if qs.contains(&i) { continue; }
        qs.push(i); log_a += (base[i].p as f64).ln();
    }
    if qs.len() < 2 { qs.push(pool[0]); qs.push(pool[1]); }
    if close { repair_a(&mut qs, base, target_log); }
    qs.sort_unstable();
    let a: BigUint = qs.iter().fold(BigUint::one(), |acc, &i| acc * base[i].p);
    // B_l = (a/q_l) · (t_l · (a/q_l)⁻¹ mod q_l), so B_l ≡ t_l (mod q_l), B_l ≡ 0 mod the others
    let bs: Vec<BigUint> = qs.iter().map(|&i| {
        let q = base[i].p;
        let aq = &a / q;
        let aq_mod = (&aq % q).to_u64().unwrap();
        let gamma = (base[i].root as u128 * modinv(aq_mod, q) as u128 % q as u128) as u64;
        let gamma = if gamma > q / 2 { q - gamma } else { gamma };
        aq * gamma
    }).collect();
    let ainv: Vec<u64> = base.iter().map(|q| { let am = (&a % q.p).to_u64().unwrap(); if am == 0 { 0 } else { modinv(am, q.p) } }).collect();
    Family { a, qs, bs, ainv }
}

/// Keep the number of factors and replace at most one by an unused odd base
/// prime to minimize the logarithmic distance from the usual SIQS target.
/// Selection is approximate; a, b and all relations are subsequently exact.
fn repair_a(qs: &mut [usize], base: &[Prime], target_log: f64) {
    let log_a: f64 = qs.iter().map(|&i| (base[i].p as f64).ln()).sum();
    let mut best = ((log_a - target_log).abs(), None);
    for (slot, &old) in qs.iter().enumerate() {
        let rest = log_a - (base[old].p as f64).ln();
        let wanted = (target_log - rest).exp();
        let at = base.partition_point(|q| (q.p as f64) < wanted);
        // Enough neighbors to step past the already-used factors.
        let radius = qs.len() + 2;
        for (j, q) in base.iter().enumerate().take((at + radius).min(base.len())).skip(at.saturating_sub(radius)) {
            if q.p == 2 || q.root == 0 || (j != old && qs.contains(&j)) { continue; }
            let error = (rest + (q.p as f64).ln() - target_log).abs();
            if error < best.0 { best = (error, Some((slot, j))); }
        }
    }
    if let Some((slot, j)) = best.1 { qs[slot] = j; }
}

/// Move a known reduced residue by a known reduced delta, with one correction.
fn shift_residue(r: u32, delta: u32, p: u32, add: bool) -> u32 {
    debug_assert!(p > 0 && r < p && delta < p);
    if add {
        if r >= p - delta { r - (p - delta) } else { r + delta }
    } else if r >= delta { r - delta } else { p - (delta - r) }
}

/// The same shift without a data-dependent branch: the sum (or difference) and its wrap
/// are both formed and the smaller is the one in range, since the out-of-range one wrapped
/// past 2³². The large-prime bucket walk mispredicted the branches of `shift_residue` on
/// nearly every prime; this is a `min` instead.
#[inline(always)]
fn shift_branchless(r: u32, delta: u32, p: u32, add: bool) -> u32 {
    debug_assert!(p < 1 << 31 && r < p && delta < p);
    if add { let t = r + delta; t.min(t.wrapping_sub(p)) } else { let t = r.wrapping_sub(delta); t.min(t.wrapping_add(p)) }
}

/// Add `lg` at every position of two strides of `p` inside the block `[lo, hi)`, each stride
/// starting at its absolute `next` position, and advance both to their first position at or
/// past `hi`. A root already past the block is left alone (an earlier implementation
/// version clamped it to the block end and lost it). Equal roots, a prime dividing k·n, are
/// one stride.
#[inline]
fn fill_block(block: &mut [u8], lo: usize, hi: usize, p: usize, lg: u8, next: &mut [u32; 2]) {
    let len = hi - lo;
    let (n0, n1) = (next[0] as usize, next[1] as usize);
    if n0 == n1 {
        if n0 < hi { let mut c = n0 - lo; while c < len { block[c] = block[c].wrapping_add(lg); c += p; } next[0] = (lo + c) as u32; next[1] = next[0]; }
        return;
    }
    let (mut c0, mut c1) = (if n0 < hi { n0 - lo } else { len }, if n1 < hi { n1 - lo } else { len });
    // the two roots together: two independent chains of load-add-store in flight
    while c0 < len && c1 < len { block[c0] = block[c0].wrapping_add(lg); block[c1] = block[c1].wrapping_add(lg); c0 += p; c1 += p; }
    while c0 < len { block[c0] = block[c0].wrapping_add(lg); c0 += p; }
    while c1 < len { block[c1] = block[c1].wrapping_add(lg); c1 += p; }
    if n0 < hi { next[0] = (lo + c0) as u32; }
    if n1 < hi { next[1] = (lo + c1) as u32; }
}

/// The bucket store: for each polynomial of a batch and each block, a fixed-capacity run of
/// packed entries, plus one spare run per polynomial that a miss is written into. A hit is then
/// an unconditional store and a conditional length increment, with no branch for the sieve to
/// mispredict (a wide prime hits the interval with probability width/p, the worst case for a
/// predictor). The capacity grows when a chunk of primes could overrun it.
struct Buckets { data: Vec<u32>, len: Vec<u32>, cap: usize, stride: usize }

impl Buckets {
    fn new(batch: usize, blocks: usize) -> Self {
        let stride = blocks + 1;
        Buckets { data: vec![0; batch * stride * 4096], len: vec![0; batch * stride], cap: 4096, stride }
    }
    fn clear(&mut self, count: usize) { self.len[..count * self.stride].iter_mut().for_each(|l| *l = 0); }
    #[inline(always)]
    fn push(&mut self, t: usize, idx: usize, entry: u32, hit: bool) {
        let slot = t * self.stride + idx;
        let l = self.len[slot] as usize;
        self.data[slot * self.cap + l] = entry;
        self.len[slot] = (l + hit as usize) as u32;
    }
    /// Room for `need` more entries in every run of the first `count` polynomials.
    fn reserve(&mut self, count: usize, need: usize) {
        let worst = self.len[..count * self.stride].iter().max().copied().unwrap_or(0) as usize + need;
        if worst <= self.cap { return; }
        let cap = worst.next_power_of_two().max(2 * self.cap);
        let mut data = vec![0u32; self.len.len() * cap];
        for (slot, &l) in self.len.iter().enumerate() { data[slot * cap..slot * cap + l as usize].copy_from_slice(&self.data[slot * self.cap..slot * self.cap + l as usize]); }
        self.data = data;
        self.cap = cap;
    }
    fn run(&self, t: usize, idx: usize) -> &[u32] { let slot = t * self.stride + idx; &self.data[slot * self.cap..slot * self.cap + self.len[slot] as usize] }
}

/// A cofactor between the large-prime bound and its square: two primes each below the bound,
/// or nothing. Below pmax² it is prime and too big; above, it must be composite and split.
fn split_cofactor(v: u128, bound: u128, pmax: u128) -> Option<Vec<(u64, u32)>> {
    if !double_large() || v >= bound * bound || v < pmax * pmax || v > u64::MAX as u128 { return None; }
    let v = v as u64;
    if is_prime_u64(v) { return None; }
    let a = rho_u64(v)?;
    let b = v / a;
    if a as u128 >= bound || b as u128 >= bound { return None; }
    if a == b { Some(vec![(a, 2)]) } else { Some(vec![(a.min(b), 1), (a.max(b), 1)]) }
}

/// Sieve every b of one a. Returns the relations and the number of polynomials used.
fn sieve_family(n: &BigUint, base: &[Prime], fam: &Family, params: &Params, threshold_log: f32) -> (Vec<Relation>, usize) {
    let half = params.half as i64;
    let width = 2 * params.half + 1;
    let s = fam.qs.len();
    let pmax = base.last().unwrap().p;
    let large_bound = BigUint::from(pmax * params.large);
    let nb = BigInt::from(n.clone());
    let ab = BigInt::from(fam.a.clone());
    // start with b = Σ B_l
    let mut b: BigInt = fam.bs.iter().fold(BigInt::zero(), |acc, x| acc + BigInt::from(x.clone()));
    // roots per prime for this b: x ≡ a⁻¹(±t − b) (mod p); and the Gray-code deltas 2·B_l·a⁻¹
    let mut roots: Vec<[u32; 2]> = vec![[0; 2]; base.len()];
    // deltas[v][i]: for Gray bit v, the root shift 2·B_v·a⁻¹ mod p of every prime, contiguous
    // over the primes so a polynomial switch streams through one array (a large-input
    // review: the per-prime layout cost a cache miss a prime and a third of the time at a
    // 60,000-prime base)
    let deltas: Vec<Vec<u32>> = fam.bs.iter().map(|bl| base.iter().enumerate().map(|(i, q)| ((2u128 * (bl % q.p).to_u64().unwrap() as u128 * fam.ainv[i] as u128) % q.p as u128) as u32).collect()).collect();
    let set_roots = |b: &BigInt, roots: &mut Vec<[u32; 2]>| {
        for (i, q) in base.iter().enumerate() {
            if fam.ainv[i] == 0 { continue; }
            let bm = b.mod_floor(&BigInt::from(q.p)).to_u64().unwrap();
            let r1 = ((q.root + q.p - bm) % q.p) as u128 * fam.ainv[i] as u128 % q.p as u128;
            let r2 = ((2 * q.p - q.root - bm) % q.p) as u128 * fam.ainv[i] as u128 % q.p as u128;
            roots[i] = [r1 as u32, r2 as u32];
        }
    };
    set_roots(&b, &mut roots);
    let mut out = Vec::new();
    let mut acc = vec![0u8; width + 8];
    let half_mod: Vec<u32> = base.iter().map(|q| (params.half as u64 % q.p) as u32).collect();
    let mut starts: Vec<[u32; 2]> = vec![[0; 2]; base.len()];
    // **the bucket sieve** (msieve's layout). The interval is sieved in blocks of `seg` bytes
    // that fit the L1 cache; a prime wider than a block hits a block at most once, so walking
    // every such prime through every block would cost more than the hits. Those primes are
    // instead handled once per polynomial: each hit is written as (position, prime index) into
    // the bucket of its block, the bucket is added to the block while it is in cache, and a
    // candidate reads which large primes divide it from its block's bucket instead of testing
    // the whole base. The first version kept one bucket for the whole interval and found a
    // candidate's primes through a position→rank table the size of the interval; at 90 digits
    // that random walk tripled the division time (30 → 95 thread-seconds) and lost overall.
    // the rule: once the largest prime is wider than the interval, or the interval is past 2^17
    // (at 69 digits, 7,000 primes to 77k over 2^17, it won 28% although no prime is wider than
    // the interval; at 66 digits over 2^16 it lost 7%)
    let bucket_on = bucket_sieve().unwrap_or(pmax as usize > width || width > 1 << 17);
    let seg = segment_len(bucket_on, base.len());
    let seg_log = seg.trailing_zeros();
    // The large primes are also walked through a **batch of polynomials at a time** (msieve):
    // one prime's two positions move by a known delta at every Gray-code switch, so a prime-major
    // loop keeps the prime in registers and writes its hits for the next `batch` polynomials in
    // one go, where the per-polynomial layout paid three passes over the base (switch, start,
    // hit) for every polynomial.
    // an entry packs the prime index above the position inside its block into one u32, so the
    // base must stay below 2^(32 − log₂ seg) primes (131,072 at the 32 KB block) for the bucket
    // and it pays only once the base outgrows the interval: at 66 digits (largest prime 100k,
    // interval 131k) it cost 5–10%, at 90 digits (1.3M against 262k) it doubled the sieve rate
    let first_large = if bucket_on && base.len() < 1usize << (32 - seg_log) { base.partition_point(|q| (q.p as usize) <= seg) } else { base.len() };
    let seg_mask = (seg - 1) as u32;
    let batch = bucket_batch();
    let blocks = width.div_ceil(seg);
    let mut buckets = Buckets::new(batch, blocks);
    // the large primes' current positions in the interval (start = root + half mod p), kept
    // `batch` polynomials ahead of the small primes' roots
    let mut lpos: Vec<[u32; 2]> = (first_large..base.len()).map(|i| if fam.ainv[i] == 0 { [0; 2] } else { [shift_residue(roots[i][0], half_mod[i], base[i].p as u32, true), shift_residue(roots[i][1], half_mod[i], base[i].p as u32, true)] }).collect();
    // the Gray-code step from polynomial t−1 to t: which B_v flips, and whether it is added
    let gray_step = |t: usize| { let v = t.trailing_zeros() as usize; (v, (t >> (v + 1)) & 1 == 0) };
    // the primes of a that lie beyond the interval, checked by division on every candidate
    let large_a: Vec<usize> = fam.qs.iter().copied().filter(|&qi| qi >= first_large).collect();
    let polys = 1usize << (s - 1);
    // the small primes are not sieved: take their expected log off the threshold, then start
    // every cell at 128 − threshold so a hit is a byte with its high bit set, found a word at a time
    let small = small_bound();
    let small_log: f32 = base.iter().filter(|q| q.p < small).map(|q| 2.0 * (q.p as f32).log2() / q.p as f32).sum();
    let thr = (threshold_log - small_log).clamp(1.0, 120.0) as u8;
    let init = 128u8 - thr;
    for poly in 0..polys {
        // sieve: position i is x = i − half; p | Q(x) iff x ≡ root (mod p)
        let t0 = clock::Instant::now();
        // in blocks that fit L1: every prime's next position is carried from block to block,
        // so the random byte adds land in cache; the hits of a block are read while it is there
        let mut next: Vec<[u32; 2]> = Vec::with_capacity(first_large);
        for (i, q) in base.iter().enumerate().take(first_large) {
            if fam.ainv[i] == 0 { next.push([u32::MAX; 2]); continue; }
            starts[i] = [shift_residue(roots[i][0], half_mod[i], q.p as u32, true), shift_residue(roots[i][1], half_mod[i], q.p as u32, true)];
            next.push(starts[i]);
        }
        acc[..width].iter_mut().for_each(|v| *v = init);
        let slot = poly % batch;
        let tb = clock::Instant::now();
        if slot == 0 {
            // the large primes into the buckets of the next `batch` polynomials, prime-major; a
            // prime between the block and the interval may land more than once, so each root
            // strides to the end
            let count = batch.min(polys - poly);
            let steps: Vec<(usize, bool)> = (1..count).map(|t| gray_step(poly + t)).collect();
            let dvs: Vec<(&[u32], bool)> = steps.iter().map(|&(v, add)| (&deltas[v][..], add)).collect();
            buckets.clear(count);
            // primes wider than the interval hit at most once per root; the ones between the
            // block and the interval stride on
            let first_wide = base.partition_point(|q| (q.p as usize) < width).max(first_large);
            let last = if poly + count < polys { Some(gray_step(poly + count)) } else { None };
            const CHUNK: usize = 128;
            for i in first_large..base.len() {
                if (i - first_large) % CHUNK == 0 { buckets.reserve(count, 2 * CHUNK * (width / base[i].p as usize + 1)); }
                if fam.ainv[i] == 0 { continue; }
                let p = base[i].p as u32;
                let tag = (i as u32) << seg_log;
                let [mut p0, mut p1] = lpos[i - first_large];
                let two = p0 != p1;
                if i >= first_wide {
                    for t in 0..count {
                        if t > 0 { let (dv, add) = dvs[t - 1]; let d = dv[i]; p0 = shift_branchless(p0, d, p, add); if two { p1 = shift_branchless(p1, d, p, add); } }
                        // a miss goes to the spare run; the last block may be short, so the
                        // test is against the width, not the block count
                        let h0 = (p0 as usize) < width;
                        buckets.push(t, if h0 { p0 as usize >> seg_log } else { blocks }, tag | (p0 & seg_mask), h0);
                        if two { let h1 = (p1 as usize) < width; buckets.push(t, if h1 { p1 as usize >> seg_log } else { blocks }, tag | (p1 & seg_mask), h1); }
                    }
                } else {
                    for t in 0..count {
                        if t > 0 { let (dv, add) = dvs[t - 1]; let d = dv[i]; p0 = shift_branchless(p0, d, p, add); if two { p1 = shift_branchless(p1, d, p, add); } }
                        let mut pos = p0 as usize;
                        while pos < width { buckets.push(t, pos >> seg_log, tag | (pos as u32 & seg_mask), true); pos += p as usize; }
                        if two { let mut pos = p1 as usize; while pos < width { buckets.push(t, pos >> seg_log, tag | (pos as u32 & seg_mask), true); pos += p as usize; } }
                    }
                }
                // leave the prime at the first polynomial of the next batch
                if let Some((v, add)) = last { let d = deltas[v][i]; p0 = shift_branchless(p0, d, p, add); if two { p1 = shift_branchless(p1, d, p, add); } }
                lpos[i - first_large] = [p0, p1];
            }
        }
        let bucket_ns = tb.elapsed().as_nanos() as u64;
        STATS[5].fetch_add(bucket_ns, Relaxed);
        let mut hits: Vec<usize> = Vec::new();
        let mut scan_ns = 0u64;
        let mut lo = 0usize;
        while lo < width {
            let hi = (lo + seg).min(width);
            let block = &mut acc[lo..hi];
            for (i, q) in base.iter().enumerate().take(first_large) {
                if fam.ainv[i] == 0 || q.p < small { continue; }
                let p = q.p as usize;
                let lg = q.log;
                fill_block(block, lo, hi, p, lg, &mut next[i]);
            }
            for &e in buckets.run(slot, lo >> seg_log) { let c = (e & seg_mask) as usize; block[c] = block[c].wrapping_add(base[(e >> seg_log) as usize].log); }
            let ts = clock::Instant::now();
            let full = block.len() / 8 * 8;
            for (w, chunk) in block[..full].chunks_exact(8).enumerate() {
                let word = u64::from_le_bytes(chunk.try_into().unwrap());
                if word & 0x8080_8080_8080_8080 != 0 { for (k, &v) in chunk.iter().enumerate() { if v & 0x80 != 0 { hits.push(lo + w * 8 + k); } } }
            }
            for (k, &v) in block[full..].iter().enumerate() { if v & 0x80 != 0 { hits.push(lo + full + k); } }
            scan_ns += ts.elapsed().as_nanos() as u64;
            lo = hi;
        }
        STATS[0].fetch_add(t0.elapsed().as_nanos() as u64 - scan_ns - bucket_ns, Relaxed);
        STATS[1].fetch_add(scan_ns, Relaxed);

        STATS[4].fetch_add(hits.len() as u64, Relaxed);
        let t2 = clock::Instant::now();
        // c = (b² − n)/a; Q(x)/a = a x² + 2 b x + c. Below about 57 digits this fits i128 and
        // the whole candidate step runs in machine integers; above, the same in BigInt.
        let c = (&b * &b - &nb) / &ab;
        let small = ab.bits() + 2 * (params.half as u64).ilog2() as u64 + 4 < 126 && c.bits() < 125 && b.bits() + (params.half as u64).ilog2() as u64 + 3 < 126;
        let (a128, b128, c128) = if small { (ab.to_i128().unwrap(), b.to_i128().unwrap(), c.to_i128().unwrap()) } else { (0, 0, 0) };
        // the large primes that landed on a candidate are the bucket entries at its position
        let bucket_hits = |i: usize| buckets.run(slot, i >> seg_log).iter().filter(move |&&e| e & seg_mask == i as u32 & seg_mask).map(|&e| (e >> seg_log) as usize);
        for &i in &hits {
            let xi = i as i64 - half;
            let mut exps: Vec<(u32, u32)> = fam.qs.iter().map(|&qi| (qi as u32 + 1, 1)).collect();
            let bump = |exps: &mut Vec<(u32, u32)>, c: u32| { match exps.iter_mut().find(|(d, _)| *d == c) { Some(t) => t.1 += 1, None => exps.push((c, 1)) } };
            let large;
            if small {
                let x = xi as i128;
                let qa = (a128 * x + 2 * b128) * x + c128;
                if qa == 0 { continue; }
                if qa < 0 { exps.push((0, 1)); }
                let mut v = qa.unsigned_abs();
                let iu = i as u32;
                for (j, q) in base.iter().enumerate().take(first_large) {
                    let p = q.p as u128;
                    let divides = if fam.ainv[j] == 0 || q.p < 5 { v % p == 0 } else {
                        let [s0, s1] = starts[j];
                        q.divides(iu.abs_diff(s0)) || q.divides(iu.abs_diff(s1))
                    };
                    if divides { while v % p == 0 { v /= p; bump(&mut exps, j as u32 + 1); } }
                }
                for j in bucket_hits(i) { let p = base[j].p as u128; while v % p == 0 { v /= p; bump(&mut exps, j as u32 + 1); } }
                for &j in &large_a { let p = base[j].p as u128; while v % p == 0 { v /= p; bump(&mut exps, j as u32 + 1); } }
                let lb = (pmax * params.large) as u128;
                if v == 1 { large = vec![]; }
                else if v < lb { large = vec![(v as u64, 1)]; }
                else if let Some(pair) = split_cofactor(v, lb, pmax as u128) { large = pair; }
                else { continue; }
            } else {
                let x = BigInt::from(xi);
                let qa: BigInt = (&ab * &x + BigInt::from(2) * &b) * &x + &c;
                if qa.is_zero() { continue; }
                if qa.is_negative() { exps.push((0, 1)); }
                let mut v = qa.magnitude().clone();
                let iu = i as u32;
                for (j, q) in base.iter().enumerate().take(first_large) {
                    let divides = if fam.ainv[j] == 0 || q.p < 5 { (&v % q.p).is_zero() } else {
                        let [s0, s1] = starts[j];
                        q.divides(iu.abs_diff(s0)) || q.divides(iu.abs_diff(s1))
                    };
                    if divides { while (&v % q.p).is_zero() { v /= q.p; bump(&mut exps, j as u32 + 1); } }
                }
                for j in bucket_hits(i) { while (&v % base[j].p).is_zero() { v /= base[j].p; bump(&mut exps, j as u32 + 1); } }
                for &j in &large_a { while (&v % base[j].p).is_zero() { v /= base[j].p; bump(&mut exps, j as u32 + 1); } }
                let lb = (pmax * params.large) as u128;
                if v.is_one() { large = vec![]; }
                else if v < large_bound { large = vec![(v.to_u64().unwrap(), 1)]; }
                else if v.bits() <= 63 { match split_cofactor(v.to_u128().unwrap(), lb, pmax as u128) { Some(pair) => large = pair, None => continue } }
                else { continue; }
            }
            let y = (&ab * BigInt::from(xi) + &b).mod_floor(&nb).to_biguint().unwrap();
            out.push(Relation { y, exps, large });
        }
        tick(2, t2);
        // next b by the Gray code: flip B_l where l is the lowest set bit of poly+1
        let t3 = clock::Instant::now();
        if poly + 1 < polys {
            let v = (poly + 1).trailing_zeros() as usize;
            let sign_plus = ((poly + 1) >> (v + 1)) & 1 == 0;   // whether B_v is added or subtracted next
            let bl = BigInt::from(fam.bs[v].clone());
            if sign_plus { b -= 2 * &bl; } else { b += 2 * &bl; }
            let dv = &deltas[v];
            for (i, q) in base.iter().enumerate().take(first_large) {
                if fam.ainv[i] == 0 { continue; }
                let d = dv[i];
                let p = q.p as u32;
                for r in roots[i].iter_mut() { *r = shift_residue(*r, d, p, sign_plus); }
            }
        }
        tick(3, t3);
    }
    (out, polys)
}

/// **The cycle finder** for partial relations. A partial with large primes p (and q) is an
/// edge p–q of a graph whose vertex 1 stands for "no second prime"; a cycle multiplies to a
/// relation in which every large prime appears an even number of times. A spanning forest is
/// kept as union-find plus the tree edges; an edge that lands inside one tree closes a cycle,
/// and the path between its ends through the tree edges is the rest of it.
#[derive(Default)]
pub struct Cycles { pub edges: Vec<Relation>, up: std::collections::HashMap<u64, (u64, usize)>, size: std::collections::HashMap<u64, usize> }

impl Cycles {
    fn root(&self, mut v: u64) -> (u64, usize) {
        let mut depth = 0;
        while let Some(&(p, _)) = self.up.get(&v) { v = p; depth += 1; }
        (v, depth)
    }
    /// Make `x` the root of its tree by reversing the parent pointers on its path to the root.
    fn reroot(&mut self, x: u64) {
        let mut path = Vec::new();
        let mut v = x;
        while let Some(&(p, e)) = self.up.get(&v) { path.push((v, p, e)); v = p; }
        for &(child, parent, e) in &path { self.up.remove(&child); self.up.insert(parent, (child, e)); }
        self.up.remove(&x);
    }
    /// Add a partial; when it closes a cycle, the relations of that cycle (this one included).
    pub fn add(&mut self, r: Relation) -> Option<Vec<Relation>> {
        // Only odd exponents need partners. Preserve even powers in the relation for
        // square-root reconstruction, but never put them into the partial forest.
        let mut odd = r.large.iter().filter(|(_, e)| e % 2 == 1).map(|(p, _)| *p);
        let (u, v) = match (odd.next(), odd.next(), odd.next()) {
            (None, None, None) => return Some(vec![r]),
            (Some(p), None, None) => (1u64, p),
            (Some(p), Some(q), None) => (p, q),
            _ => return None,
        };
        let idx = self.edges.len();
        self.edges.push(r);
        let (ru, _) = self.root(u);
        let (rv, _) = self.root(v);
        if ru != rv {
            // join: re-root the smaller tree at its end of the edge and hang it under the other
            let (su, sv) = (*self.size.get(&ru).unwrap_or(&1), *self.size.get(&rv).unwrap_or(&1));
            let (small, big, small_root, big_root) = if su <= sv { (u, v, ru, rv) } else { (v, u, rv, ru) };
            self.reroot(small);
            self.up.insert(small, (big, idx));
            self.size.remove(&small_root);
            self.size.insert(big_root, su + sv);
            return None;
        }
        // the cycle: both paths to the root, cut at the common ancestor
        let mut seen: std::collections::HashMap<u64, usize> = Default::default();
        let mut x = u;
        let mut i = 0;
        loop { seen.insert(x, i); i += 1; match self.up.get(&x) { Some(&(p, _)) => x = p, None => break } }
        let mut path = vec![self.edges[idx].clone()];
        let mut y = v;
        let mut up_from_v: Vec<usize> = Vec::new();
        while !seen.contains_key(&y) { let &(p, e) = self.up.get(&y).unwrap(); up_from_v.push(e); y = p; }
        let meet = y;
        let mut x = u;
        while x != meet { let &(p, e) = self.up.get(&x).unwrap(); path.push(self.edges[e].clone()); x = p; }
        for e in up_from_v { path.push(self.edges[e].clone()); }
        Some(path)
    }
}

/// Collect relations across threads until the base is covered, merging partials.
pub fn relations(n: &BigUint, base: &[Prime], params: &Params, need: usize) -> (Vec<Relation>, usize, usize) {
    // threshold: log2 of the typical value M·√(n/2) less the slack and the large-prime room
    let size = (params.half as f64 * (n.to_f64().unwrap() / 2.0).sqrt()).log2();
    let pmax = base.last().unwrap().p as f64;
    let threshold = (size - params.slack as f64 * pmax.log2() - (pmax * params.large as f64).log2()) as f32;
    let mut full: Vec<Relation> = Vec::new();
    let mut graph = Cycles::default();
    let mut round = 0u64;
    // progress with an estimate, like Alpern's countdown: remaining relations over the rate of
    // the last stretch, printed to stderr every few seconds when `SIQS_PROGRESS` is set
    let progress = std::env::var("SIQS_PROGRESS").map(|v| v != "0").unwrap_or(false);
    // `SIQS_BUDGET_SECS` stops the collection after that many seconds, for sweeps that compare
    // relations a second across parameter rows without running each to the end
    let budget: Option<f64> = std::env::var("SIQS_BUDGET_SECS").ok().and_then(|v| v.parse().ok());
    let started = clock::Instant::now();
    let mut last_print = started;
    let mut mark = (started, 0usize);
    let mut polys = 0usize;
    let threads = rayon::current_num_threads().max(2) as u64;
    while full.len() < need {
        let batch: Vec<(Vec<Relation>, usize)> = (0..threads).into_par_iter().map(|t| {
            let fam = choose_a(n, base, params.half, round * threads + t + 1);
            sieve_family(n, base, &fam, params, threshold)
        }).collect();
        for (rels, np) in batch {
            polys += np;
            for r in rels {
                if r.large.is_empty() { full.push(r); continue; }
                if let Some(cycle) = graph.add(r) {
                    let mut acc = cycle[0].clone();
                    for r in &cycle[1..] { acc = acc.times(r, n); }
                    full.push(acc);
                }
            }
        }
        round += 1;
        if budget.is_some_and(|b| started.elapsed().as_secs_f64() >= b) { break; }
        if progress && last_print.elapsed().as_secs_f64() >= 5.0 {
            let now = clock::Instant::now();
            let elapsed = now.duration_since(started).as_secs_f64();
            let (t0, c0) = mark;
            let recent = (full.len().saturating_sub(c0)) as f64 / now.duration_since(t0).as_secs_f64().max(1e-9);
            let overall = full.len() as f64 / elapsed.max(1e-9);
            let rate = if full.len() > c0 + 20 { recent } else { overall };
            let left = if rate > 0.0 { (need.saturating_sub(full.len())) as f64 / rate } else { f64::INFINITY };
            eprintln!("  {}/{} relations, {} partial edges, {} polynomials, {:.0} s elapsed, about {:.0} s to go", full.len(), need, graph.edges.len(), polys, elapsed, left);
            if now.duration_since(t0).as_secs_f64() >= 20.0 { mark = (now, full.len()); }
            last_print = now;
        }
    }
    (full, polys, graph.edges.len())
}

/// A nontrivial factor of the odd composite n, or None.
pub fn factor(n: &BigUint) -> Option<BigUint> { factor_with(n, Params::for_digits(n.to_string().len())) }

pub fn factor_with(n: &BigUint, params: Params) -> Option<BigUint> {
    if *n < BigUint::from(4u32) { return None; }
    for p in primes_upto(10_000) { if (n % p).is_zero() { return if *n == BigUint::from(p) { None } else { Some(BigUint::from(p)) }; } }
    let r = n.sqrt();
    if &r * &r == *n { return Some(r); }
    if crate::big::is_prime_big(n) { return None; }
    let k = multiplier(n);
    let kn = n * k;
    let base = factor_base(&kn, params.base);
    let need = base.len() + 1 + params.extra;
    let (rels, _, _) = relations(&kn, &base, &params, need);
    let n_only = n;
    let n = &kn;
    let cols = base.len() + 1;
    let odd: Vec<Vec<u32>> = rels.iter().map(|r| { let mut v: Vec<u32> = r.exps.iter().filter(|(_, e)| e % 2 == 1).map(|(c, _)| *c).collect(); v.sort_unstable(); v }).collect();
    let keep = prune_singletons(&odd, cols);
    let kept_rows: Vec<Vec<u32>> = keep.iter().map(|&i| odd[i].clone()).collect();
    for dep in solve(&kept_rows, cols) {
        if dep.is_empty() { continue; }
        let mut x = BigUint::one();
        let mut sum = vec![0u32; cols];
        let mut larges: std::collections::HashMap<u64, u32> = Default::default();
        for &k in &dep {
            let r = &rels[keep[k]];
            x = &x * &r.y % n;
            for &(c, e) in &r.exps { sum[c as usize] += e; }
            for &(p, e) in &r.large { *larges.entry(p).or_default() += e; }
        }
        if sum.iter().any(|e| e % 2 == 1) || larges.values().any(|e| e % 2 == 1) { continue; }
        let mut y = BigUint::one();
        for (j, q) in base.iter().enumerate() { if sum[j + 1] > 0 { y = y * BigUint::from(q.p).modpow(&BigUint::from(sum[j + 1] / 2), n) % n; } }
        for (&p, &e) in &larges { y = y * BigUint::from(p).modpow(&BigUint::from(e / 2), n) % n; }
        let diff = if x >= y { &x - &y } else { &y - &x };
        for g in [diff.gcd(n_only), (&x + &y).gcd(n_only)] { if g > BigUint::one() && g < *n_only { return Some(g); } }
    }
    None
}

/// Rows to keep after removing, repeatedly, every row that is the only one with an odd
/// exponent in some column: such a row can be in no dependency, and dropping it shrinks the
/// matrix the elimination sees. The count-and-XOR form: per column the
/// number of live rows and the XOR of their ids, so when a count reaches one the XOR *is* the
/// remaining row, and the cascade runs in one pass over the nonzeros with a queue.
pub fn prune_singletons(odd: &[Vec<u32>], cols: usize) -> Vec<usize> {
    let mut count = vec![0u32; cols];
    let mut owner = vec![0usize; cols];
    for (i, r) in odd.iter().enumerate() { for &c in r { count[c as usize] += 1; owner[c as usize] ^= i; } }
    let mut alive = vec![true; odd.len()];
    let mut queue: std::collections::VecDeque<usize> = (0..cols).filter(|&c| count[c] == 1).collect();
    while let Some(c) = queue.pop_front() {
        if count[c] != 1 { continue; }
        let i = owner[c];
        if !alive[i] { continue; }
        alive[i] = false;
        for &d in &odd[i] { let d = d as usize; count[d] -= 1; owner[d] ^= i; if count[d] == 1 { queue.push_back(d); } }
    }
    (0..odd.len()).filter(|&i| alive[i]).collect()
}

/// The dependencies of the pruned matrix: dense elimination up to `LANCZOS_FROM` rows, block
/// Lanczos above it (or always with `SIQS_LANCZOS=1`, never with `SIQS_LANCZOS=0`).
pub const LANCZOS_FROM: usize = 6000;
pub fn solve(rows: &[Vec<u32>], cols: usize) -> Vec<Vec<usize>> {
    let force = std::env::var("SIQS_LANCZOS").ok();
    let use_lanczos = match force.as_deref() { Some("1") => true, Some("0") => false, _ => rows.len() >= LANCZOS_FROM };
    if !use_lanczos { return dependencies_sparse(rows, cols); }
    let mut out = Vec::new();
    for seed in 1..=4u64 {
        out.extend(crate::lanczos::dependencies(rows, cols, seed));
        if out.len() >= 32 { break; }
    }
    if out.is_empty() { dependencies_sparse(rows, cols) } else { out }
}

/// The complete factorisation of `n`: trial division to 100,000, perfect squares, primality
/// tests, and the sieve on every composite piece until all are prime. Sorted, with
/// multiplicity. A piece the sieve cannot split (it never happened in the tests) is kept whole.
pub fn factor_completely(n: &BigUint) -> Vec<BigUint> {
    fn split(n: &BigUint, out: &mut Vec<BigUint>) {
        if n.is_one() { return; }
        if crate::big::is_prime_big(n) { out.push(n.clone()); return; }
        let r = n.sqrt();
        if &r * &r == *n { split(&r, out); split(&r, out); return; }
        match factor(n) {
            Some(f) => { split(&f, out); split(&(n / &f), out); }
            None => out.push(n.clone()),
        }
    }
    let mut n = n.clone();
    let mut factors = Vec::new();
    if n.is_zero() { return factors; }
    for p in crate::smooth::primes_upto(100_000) {
        while (&n % p).is_zero() { n /= p; factors.push(BigUint::from(p)); }
    }
    if !n.is_one() { split(&n, &mut factors); }
    factors.sort();
    factors
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn incremental_factor_base_matches_the_original_scan() {
fn original(n: &BigUint, count: usize) -> Vec<Prime> {
    let mut out = vec![Prime::new(2, (n % 2u32).to_u64().unwrap())];
    let mut bound = 1000u64;
    loop {
        for p in primes_upto(bound).into_iter().skip(1) {
            if out.len() >= count { return out; }
            if out.iter().any(|q| q.p == p) { continue; }
            let np = (n % p).to_u64().unwrap();
            if np == 0 { out.push(Prime::new(p, 0)); } else if legendre(np, p) == 1 { out.push(Prime::new(p, sqrt_mod(np, p))); }
        }
        bound *= 2;
    }
}

        for n in [BigUint::from(11u32), BigUint::from(1234567891011121314u64), (BigUint::one()<<180usize)+BigUint::from(37u32)] {
            for count in [1,150,400,1800] {
                let a=factor_base(&n,count);let b=original(&n,count);
                assert_eq!(a.iter().map(|q|(q.p,q.root,q.log,q.inv,q.limit)).collect::<Vec<_>>(),b.iter().map(|q|(q.p,q.root,q.log,q.inv,q.limit)).collect::<Vec<_>>());
            }
        }
    }
    #[test]
    fn repaired_families_keep_exact_crt_contracts() {
        let n = BigUint::parse_bytes(b"13729273257095901138749901421978803731471830970253115531", 10).unwrap();
        let base = factor_base(&n, 1800);
        let half = 1 << 15;
        let target = (n.to_f64().unwrap() * 2.0).sqrt() / half as f64;
        for seed in 1..=64 {
            let old = choose_a_mode(&n, &base, half, seed, false);
            let new = choose_a_mode(&n, &base, half, seed, true);
            assert_eq!(old.qs.len(), new.qs.len());
            assert!(new.qs.windows(2).all(|w| w[0] < w[1]));
            assert!((new.a.to_f64().unwrap() / target).ln().abs() <= (old.a.to_f64().unwrap() / target).ln().abs() + 1e-12);
            assert_eq!(new.a, new.qs.iter().fold(BigUint::one(), |a, &i| a * base[i].p));
            for (slot, bl) in new.bs.iter().enumerate() {
                for (j, &i) in new.qs.iter().enumerate() {
                    if slot == j { assert_eq!((bl * bl) % base[i].p, &n % base[i].p); }
                    else { assert!((bl % base[i].p).is_zero()); }
                }
            }
            let mut b: BigInt = new.bs.iter().fold(BigInt::zero(), |acc, v| acc + BigInt::from(v.clone()));
            let a = BigInt::from(new.a.clone()); let ni = BigInt::from(n.clone());
            for poly in 0..(1usize << (new.qs.len() - 1)) {
                assert!(((&b * &b - &ni) % &a).is_zero());
                if poly + 1 < (1usize << (new.qs.len() - 1)) {
                    let v = (poly + 1).trailing_zeros() as usize;
                    if ((poly + 1) >> (v + 1)) & 1 == 0 { b -= BigInt::from(new.bs[v].clone()) * 2; }
                    else { b += BigInt::from(new.bs[v].clone()) * 2; }
                }
            }
        }
    }
    #[test]
    fn bounded_stride_updates_match_modulo_without_overflow() {
        for p in [1u32, 2, 3, 5, 7, 31, 65537, 1000003, u32::MAX - 4, u32::MAX] {
            for r in [0, p / 2, p - 1] {
                for delta in [0, p / 2, p - 1] {
                    assert_eq!(shift_residue(r, delta, p, true) as u64, (r as u64 + delta as u64) % p as u64);
                    assert_eq!(shift_residue(r, delta, p, false) as u64, (r as u64 + p as u64 - delta as u64) % p as u64);
                }
            }
            let mut state = 20261109u32;
            for _ in 0..10000 {
                state = state.wrapping_mul(1664525).wrapping_add(1013904223);
                let r = state % p;
                state = state.wrapping_mul(1664525).wrapping_add(1013904223);
                let delta = state % p;
                assert_eq!(shift_residue(r, delta, p, true) as u64, (r as u64 + delta as u64) % p as u64);
                assert_eq!(shift_residue(r, delta, p, false) as u64, (r as u64 + p as u64 - delta as u64) % p as u64);
            }
        }
    }

    #[test]
    fn the_divisibility_test_agrees_with_the_remainder() {
        for p in [3u64, 5, 7, 11, 65537, 1_000_003, 4_294_967_291] {
            let q = Prime::new(p, 1);
            assert_eq!(q.inv.wrapping_mul(p as u32), 1, "inverse mod 2^32 of {p}");
            let mut state = 20261109u32;
            for d in (0..5000u32).chain((0..5000).map(|_| { state = state.wrapping_mul(1664525).wrapping_add(1013904223); state })) {
                assert_eq!(q.divides(d), d as u64 % p == 0, "{d} by {p}");
            }
        }
    }
    #[test]
    fn small_prime_tests_and_rho() {
        assert!(is_prime_u64(1_000_000_007) && !is_prime_u64(1_000_000_007 * 3) && is_prime_u64(18446744073709551557));
        let n = 1_000_003u64 * 1_000_033;
        let f = rho_u64(n).unwrap(); assert!(n % f == 0 && f > 1 && f < n);
        assert_eq!(split_cofactor((1_000_003u128) * 1_000_033, 2_000_000, 100_000), Some(vec![(1_000_003, 1), (1_000_033, 1)]));
        assert_eq!(split_cofactor(1_000_003u128 * 1_000_003, 2_000_000, 100_000), Some(vec![(1_000_003, 2)]));
        assert_eq!(split_cofactor(1_000_003u128 * 3_000_017, 2_000_000, 100_000), None, "one prime over the bound");
    }
    #[test]
    fn square_large_cofactors_are_immediately_usable_and_do_not_join_odd_edges() {
        let rel = |large| Relation { y: BigUint::from(101u32), exps: vec![], large };
        let mut graph = Cycles::default();
        assert!(graph.add(rel(vec![(101, 1)])).is_none());
        let square = graph.add(rel(vec![(101, 2)])).expect("an even large exponent needs no partner");
        assert_eq!(square.len(), 1);
        assert_eq!(square[0].large, vec![(101, 2)]);
        let n = BigUint::from(1009u32);
        assert_eq!(&square[0].y * &square[0].y % &n, BigUint::from(101u32).pow(2) % &n);
        assert_eq!(graph.edges.len(), 1, "the square must not enter the forest");
        let pair = graph.add(rel(vec![(101, 1)])).unwrap();
        assert_eq!(pair.len(), 2);
        assert!(pair[0].times(&pair[1], &n).large.iter().all(|(_, e)| e % 2 == 0));
        assert!(graph.add(rel(vec![(103, 2), (107, 1)])).is_none());
        let pair = graph.add(rel(vec![(107, 1)])).unwrap();
        assert!(pair[0].times(&pair[1], &n).large.iter().all(|(_, e)| e % 2 == 0));
    }
    #[test]
    fn cycles_close_and_multiply_to_squares() {
        let n = BigUint::from(1_000_000_007u64);
        let rel = |l: Vec<(u64, u32)>| Relation { y: BigUint::from(3u32), exps: vec![(1, 1)], large: l };
        let mut g = Cycles::default();
        assert!(g.add(rel(vec![(101, 1)])).is_none());
        assert!(g.add(rel(vec![(101, 1), (103, 1)])).is_none());
        let c = g.add(rel(vec![(103, 1)])).expect("1–101–103–1 closes");
        assert_eq!(c.len(), 3);
        let mut acc = c[0].clone(); for r in &c[1..] { acc = acc.times(r, &n); }
        assert!(acc.large.iter().all(|(_, e)| e % 2 == 0), "{:?}", acc.large);
        assert!(g.add(rel(vec![(107, 1), (109, 1)])).is_none());
        assert!(g.add(rel(vec![(107, 1), (109, 1)])).is_some(), "a repeated pair is a two-cycle");
    }
    #[test]
    fn blocked_fill_equals_whole_fill() {
        let mut state = 12345u64;
        let mut rnd = |m: u64| { state ^= state << 13; state ^= state >> 7; state ^= state << 17; state % m };
        for _ in 0..300 {
            let width = 50 + rnd(500) as usize;
            let p = 2 + rnd(60) as usize;
            let r0 = rnd(p as u64) as u32;
            let r1 = if rnd(5) == 0 { r0 } else { rnd(p as u64) as u32 };
            let seg = 1 + rnd(64) as usize;
            let mut whole = vec![0u8; width];
            let mut nx = [r0, r1];
            fill_block(&mut whole, 0, width, p, 3, &mut nx);
            let mut blocked = vec![0u8; width];
            let mut nb = [r0, r1];
            let mut lo = 0;
            while lo < width { let hi = (lo + seg).min(width); fill_block(&mut blocked[lo..hi], lo, hi, p, 3, &mut nb); lo = hi; }
            assert_eq!(whole, blocked, "p={p} roots={r0},{r1} width={width} seg={seg}");
            assert_eq!(nx, nb, "the carried positions agree too");
        }
    }
    #[test]
    fn pruning_keeps_every_dependency() {
        // every 4-by-3 binary matrix: the rows pruned are never part of any dependency
        for bits in 0..4096u32 {
            let odd: Vec<Vec<u32>> = (0..4).map(|r| (0..3).filter(|c| bits >> (3 * r + c) & 1 == 1).collect()).collect();
            let keep = prune_singletons(&odd, 3);
            let kept: std::collections::HashSet<usize> = keep.iter().copied().collect();
            for subset in 1..16u32 {
                let rows: Vec<usize> = (0..4).filter(|r| subset >> r & 1 == 1).collect();
                let mut sum = [0u8; 3];
                for &r in &rows { for &c in &odd[r] { sum[c as usize] ^= 1; } }
                if sum == [0, 0, 0] { assert!(rows.iter().all(|r| kept.contains(r)), "{bits:012b} {rows:?}"); }
            }
        }
        // a 15,000-row chain cascades to nothing in one pass
        let chain: Vec<Vec<u32>> = (0..15000u32).map(|i| vec![i, i + 1]).collect();
        assert!(prune_singletons(&chain, 15001).is_empty());
    }
    #[test]
    fn the_polynomials_are_consistent() {
        let n = BigUint::parse_bytes(b"4534149033143998129927133103", 10).unwrap(); // 25-digit-ish composite
        let base = factor_base(&n, 150);
        assert!(base.iter().skip(1).all(|q| (q.root * q.root) % q.p == (&n % q.p).to_u64().unwrap()));
        let fam = choose_a(&n, &base, 1 << 13, 7);
        // b² ≡ n (mod a) for the first b
        let b: BigUint = fam.bs.iter().fold(BigUint::zero(), |acc, x| acc + x);
        assert_eq!((&b * &b) % &fam.a, &n % &fam.a);
    }
    #[test]
    fn factors_a_few() {
        let (p, q) = (BigUint::from(206158430209u64), BigUint::from(21990232567889u64));
        let f = factor(&(&p * &q)).unwrap(); assert!(f == p || f == q, "25 digits");
        let (p, q) = (BigUint::from(26388279066671u64), BigUint::from(22517998136864981u64));
        let f = factor(&(&p * &q)).unwrap(); assert!(f == p || f == q, "30 digits");
    }
}
