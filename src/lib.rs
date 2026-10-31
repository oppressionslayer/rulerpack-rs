//! rulerpack: XOR rulers over the primes, and the pack that rides on them.
//!
//! Rulers are Mersenne ladders scaled by a power of two. The pack family
//! (L, c, M): K = 2p ^ ruler, X = 2p ^ L*i, N = K*M + X. Decode is mod, div
//! and two XORs. Counting is Meissel with the pi-table pruning; the marks
//! table walks from a stored checkpoint instead of counting.

use std::collections::HashMap;

pub mod bands;
#[cfg(feature = "big")]
pub mod big;
#[cfg(feature = "gmp")]
pub mod gmp;

// ── rulers ──────────────────────────────────────────────────────────────────

/// m ^ (m+1) == 2^(v2(m+1)+1) - 1: 1, 3, 1, 7, 1, 3, 1, 15, ...
pub fn ladder(m: u64) -> u64 { m ^ (m + 1) }
/// The original ruler (4i) ^ (4i+7): 7, 15, 7, 31, ...
pub fn rp(i: u64) -> u64 { (4 * i) ^ (4 * i + 7) }
/// General ruler ((Li) ^ (Li + cL + 1)) ^ 1 == L * (i ^ (i + c)).
pub fn ruler(i: u64, l: u64, c: u64) -> u64 { ((l * i) ^ (l * i + c * l + 1)) ^ 1 }

// ── pack ────────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, Debug)]
pub struct Family { pub l: u64, pub c: u64, pub shift: u32 }

pub const F513: Family = Family { l: 256, c: 2, shift: 32 };
pub const F257: Family = Family { l: 4, c: 64, shift: 32 };
pub const MARKS_FAMILY: Family = Family { l: 256, c: 2, shift: 64 };   // X = 2p ^ 256i; 64 bits covers any u64 index

impl Family {
    /// (K, X) for prime p at index i. Computed in u128 so no size of p or i wraps.
    pub fn encode(&self, i: u64, p: u64) -> (u128, u128) {
        let (i, p, l) = (i as u128, p as u128, self.l as u128);
        let r = ((l * i) ^ (l * i + (self.c as u128) * l + 1)) ^ 1;
        ((2 * p) ^ r, (2 * p) ^ (l * i))
    }
    /// (i, p) from the pair. Panics if the pair did not come from this family.
    pub fn decode(&self, k: u128, x: u128) -> (u64, u64) {
        let t = k ^ x;
        let lb = self.l.trailing_zeros();
        assert!(t & ((1u128 << lb) - 1) == 0 && (t >> lb) >= self.c as u128, "not a pair from this family");
        let i = (t >> lb) - self.c as u128;
        let p = (x ^ ((self.l as u128) * i)) >> 1;
        assert!(i <= u64::MAX as u128 && p <= u64::MAX as u128, "pair does not fit u64");
        (i as u64, p as u64)
    }
    /// One integer for the pair. Panics, never wraps, if X or K overflow the fields.
    /// K and X are each about bits(2p) wide, so a u128 N carries pairs with p < 2^63.
    pub fn pack(&self, i: u64, p: u64) -> u128 {
        let (k, x) = self.encode(i, p);
        assert!(self.shift < 128 && x >> self.shift == 0, "X = {x} needs more than {} bits", self.shift);
        assert!(k >> (128 - self.shift) == 0, "K = {k} does not fit above a {}-bit field in u128", self.shift);
        (k << self.shift) | x
    }
    pub fn unpack(&self, n: u128) -> (u64, u64) {
        let mask = (1u128 << self.shift) - 1;
        self.decode(n >> self.shift, n & mask)
    }
}

// ── primes ──────────────────────────────────────────────────────────────────

pub const SEEDS: [u64; 31] = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61,
    67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131];

/// Seed rule for an odd candidate: a seed, or no seed divides it, and not a square.
pub fn allowed(h: u64) -> bool {
    if SEEDS.contains(&h) { return true; }
    if SEEDS.iter().any(|&q| h % q == 0) { return false; }
    let r = (h as f64).sqrt() as u64;
    !(r * r == h || (r + 1) * (r + 1) == h)
}

fn mulmod(a: u64, b: u64, m: u64) -> u64 { ((a as u128 * b as u128) % m as u128) as u64 }
fn powmod(mut a: u64, mut e: u64, m: u64) -> u64 {
    let mut r = 1u64; a %= m;
    while e > 0 { if e & 1 == 1 { r = mulmod(r, a, m); } a = mulmod(a, a, m); e >>= 1; }
    r
}

/// Deterministic Miller-Rabin, exact for all u64.
pub fn isprime(n: u64) -> bool {
    if n < 2 { return false; }
    for &q in &[2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        if n % q == 0 { return n == q; }
    }
    let s = (n - 1).trailing_zeros();
    let d = (n - 1) >> s;
    'base: for &a in &[2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        let mut x = powmod(a, d, n);
        if x == 1 || x == n - 1 { continue; }
        for _ in 1..s { x = mulmod(x, x, n); if x == n - 1 { continue 'base; } }
        return false;
    }
    true
}

/// Odd-only bit sieve with a rank index: pi(n) in O(1) after O(limit) build.
pub struct PiTable { limit: u64, words: Vec<u64>, prefix: Vec<u32> }

impl PiTable {
    pub fn new(limit: u64) -> Self {
        // bit b of the odd bitset = number 2b+1; word w covers numbers [128w, 128w+128)
        let nbits = (limit / 2 + 1) as usize;
        let nw = nbits.div_ceil(64);
        let mut words = vec![u64::MAX; nw];
        words[0] &= !1; // 1 is not prime
        if nbits % 64 != 0 { words[nw - 1] &= (1u64 << (nbits % 64)) - 1; }
        let mut q = 3u64;
        while q * q <= limit {
            if words[(q / 2) as usize / 64] >> ((q / 2) % 64) & 1 == 1 {
                let mut b = (q * q / 2) as usize;
                while b < nbits { words[b / 64] &= !(1u64 << (b % 64)); b += q as usize; }
            }
            q += 2;
        }
        let mut prefix = Vec::with_capacity(nw + 1);
        let mut acc = 0u32;
        for w in &words { prefix.push(acc); acc += w.count_ones(); }
        prefix.push(acc);
        PiTable { limit, words, prefix }
    }
    pub fn limit(&self) -> u64 { self.limit }
    /// Number of primes <= n.
    pub fn pi(&self, n: u64) -> u64 {
        if n < 2 { return 0; }
        let b = (n / 2) as usize;              // last odd <= n is 2b+1 when n odd, else 2b-1... handle below
        let b = if n % 2 == 0 { b - 1 } else { b };
        let w = b / 64;
        let mask = if b % 64 == 63 { u64::MAX } else { (1u64 << (b % 64 + 1)) - 1 };
        self.prefix[w] as u64 + (self.words[w] & mask).count_ones() as u64 + 1
    }
    pub fn is_prime(&self, n: u64) -> bool {
        if n == 2 { return true; }
        if n < 2 || n % 2 == 0 { return false; }
        let b = (n / 2) as usize;
        self.words[b / 64] >> (b % 64) & 1 == 1
    }
    /// Primes up to lim as a vector.
    pub fn primes_upto(&self, lim: u64) -> Vec<u64> {
        let mut v = vec![2u64];
        let mut n = 3; while n <= lim { if self.is_prime(n) { v.push(n); } n += 2; }
        v
    }
}

// ── counting: Legendre phi and Meissel ──────────────────────────────────────

const WHEEL_PRIMES: [u64; 7] = [2, 3, 5, 7, 11, 13, 17];
const WHEEL_MOD: u64 = 510510;

struct Wheel { cum: Vec<u32>, count: u64 }
impl Wheel {
    fn new() -> Self {
        let mut cum = Vec::with_capacity(WHEEL_MOD as usize + 1);
        let mut acc = 0u32;
        cum.push(0);
        for n in 1..=WHEEL_MOD {
            if WHEEL_PRIMES.iter().all(|&q| n % q != 0) { acc += 1; }
            cum.push(acc);
        }
        Wheel { count: acc as u64, cum }
    }
    fn phi(&self, x: u64) -> u64 { (x / WHEEL_MOD) * self.count + self.cum[(x % WHEEL_MOD) as usize] as u64 }
}

pub struct Counter { pub table: PiTable, primes: Vec<u64>, wheel: Wheel }

impl Counter {
    /// Ready to count any x with x^(2/3) <= limit.
    pub fn new(limit: u64) -> Self {
        let table = PiTable::new(limit);
        let primes = table.primes_upto(((limit as f64).powf(0.75) as u64) + 2); // sqrt(x) for x up to limit^1.5
        Counter { table, primes, wheel: Wheel::new() }
    }
    fn phi(&self, x: u64, a: usize, memo: &mut HashMap<(u64, usize), i64>) -> i64 {
        if a == 0 { return x as i64; }
        if a == 1 { return x.div_ceil(2) as i64; }
        if a == 7 { return self.wheel.phi(x) as i64; }
        let pa = self.primes[a - 1];
        if x < pa { return 1; }
        if x <= self.table.limit && pa * pa > x { return self.table.pi(x) as i64 - a as i64 + 1; }
        if let Some(&v) = memo.get(&(x, a)) { return v; }
        let v = self.phi(x, a - 1, memo) - self.phi(x / pa, a - 1, memo);
        memo.insert((x, a), v);
        v
    }
    fn icbrt(x: u64) -> u64 {
        let mut r = (x as f64).cbrt() as u64;
        while r * r * r > x { r -= 1; }
        while (r + 1) * (r + 1) * (r + 1) <= x { r += 1; }
        r
    }
    /// Meissel: pi(x) = phi(x,a) + a - 1 - sum_{a<i<=b} (pi(x/p_i) - i + 1).
    pub fn primepi(&self, x: u64) -> u64 {
        if x < 2 { return 0; }
        if x <= self.table.limit { return self.table.pi(x); }
        let c = Self::icbrt(x);
        let r = (x as f64).sqrt() as u64;
        let r = if (r + 1) * (r + 1) <= x { r + 1 } else { r };
        assert!(x / c <= self.table.limit, "table too small for x: need {}", x / c);
        let a = self.table.pi(c) as usize;
        let b = self.table.pi(r) as usize;
        let mut memo = HashMap::new();
        let mut total = self.phi(x, a, &mut memo) + a as i64 - 1;
        for i in (a + 1)..=b {
            total -= self.table.pi(x / self.primes[i - 1]) as i64 - i as i64 + 1;
        }
        total as u64
    }
}

// ── li, Riemann R, and the walk ─────────────────────────────────────────────

pub fn li(x: f64) -> f64 {
    if x < 2.0 { return 0.0; }
    let lx = x.ln();
    let (mut total, mut fact, mut pow2) = (0.0f64, 1.0f64, 1.0f64);
    let mut s = 0.0;
    for n in 1..200 {
        fact *= n as f64;
        if n > 1 { pow2 *= 2.0; }
        let k = (n - 1) / 2;
        if n % 2 == 1 { s += 1.0 / (2.0 * k as f64 + 1.0); }
        let term = if n % 2 == 1 { 1.0 } else { -1.0 } * lx.powi(n) / (pow2 * fact) * s;
        total += term;
        if n > 10 && term.abs() < 1e-17 * total.abs() { break; }
    }
    0.5772156649015329 + lx.ln() + x.sqrt() * total
}
const MOBIUS: [i8; 61] = [0,1,-1,-1,0,-1,1,-1,0,0,1,-1,0,-1,1,1,0,-1,0,-1,0,1,1,-1,0,0,1,0,0,-1,-1,-1,0,1,1,1,0,-1,1,1,0,-1,-1,-1,0,0,1,-1,0,0,0,1,0,-1,0,1,0,1,1,-1,0];
pub fn riemann_r(x: f64) -> f64 {
    let mut s = 0.0;
    for (n, &m) in MOBIUS.iter().enumerate().skip(1) {
        if m == 0 { continue; }
        let v = x.powf(1.0 / n as f64); if v < 2.0 { break; }
        s += m as f64 / n as f64 * li(v);
    }
    s
}
pub fn li_inverse(y: f64) -> f64 {
    if y < 2.0 { return 2.0; }
    let mut x = if y > 3.0 { y * y.ln() } else { 3.0 };
    for _ in 0..60 { let step = (li(x) - y) * x.ln(); x -= step; if step.abs() < 0.5 { break; } }
    x
}
pub fn riemann_r_inverse(y: f64) -> f64 {
    let mut x = li_inverse(y);
    for _ in 0..50 { let step = (riemann_r(x) - y) * x.ln(); x -= step; if step.abs() < 0.5 { break; } }
    x
}

pub fn next_prime(mut n: u64) -> u64 { n += if n == 2 { 1 } else { 2 }; while !(allowed(n) && isprime(n)) { n += 2; } n }
pub fn prev_prime(mut n: u64) -> u64 { n -= if n == 3 { 1 } else { 2 }; while !isprime(n) { n -= 2; } n }

impl Counter {
    /// Zero-based index i -> prime: R guesses, Meissel pins, the walk closes.
    pub fn nth_prime(&self, i: u64) -> u64 {
        const SMALL: [u64; 6] = [2, 3, 5, 7, 11, 13];
        if (i as usize) < SMALL.len() { return SMALL[i as usize]; }
        let x = riemann_r_inverse((i + 1) as f64) as u64;
        let x = if x % 2 == 0 { x - 1 } else { x };
        let mut p = if isprime(x) { x } else { prev_prime(x) };
        let mut k = self.primepi(p) as i64 - (i as i64 + 1);
        while k < 0 { p = next_prime(p); k += 1; }
        while k > 0 { p = prev_prime(p); k -= 1; }
        p
    }
    pub fn index_of(&self, p: u64) -> u64 { assert!(isprime(p)); self.primepi(p) - 1 }
}

// ── marks: our own way, no counting ─────────────────────────────────────────

use std::io::{Read, Write};

/// A table of packed checkpoints. Each entry N = pack(i, p) carries its own
/// index and prime; nothing else is stored. Lookups walk from the mark below.
pub struct Marks { pub n: Vec<u128>, pub every: u64, pub family: Family, p: Vec<u64>, last_prime: u64 }

const MAGIC: &[u8; 8] = b"RULRMRK2";   // header: magic, every (u64 LE), flags (u8: bit0 = sparse), count (u64 LE), then u128 LE marks
const MAGIC_V1: &[u8; 8] = b"RULRMRK1";

fn write_marks_file(path: &str, every: u64, sparse: bool, n: &[u128]) -> std::io::Result<()> {
    let mut f = std::fs::File::create(path)?;
    f.write_all(MAGIC)?;
    f.write_all(&every.to_le_bytes())?;
    f.write_all(&[sparse as u8])?;
    f.write_all(&(n.len() as u64).to_le_bytes())?;
    for &x in n { f.write_all(&x.to_le_bytes())?; }
    Ok(())
}

/// (every, sparse, marks). Reads v2, and v1 files (no flags byte, contiguous).
fn read_marks_file(path: &str) -> std::io::Result<(u64, bool, Vec<u128>)> {
    let mut f = std::fs::File::open(path)?;
    let mut buf = Vec::new();
    f.read_to_end(&mut buf)?;
    let bad = |m: &str| std::io::Error::new(std::io::ErrorKind::InvalidData, m.to_string());
    if buf.len() < 24 { return Err(bad("not a marks file")); }
    let (sparse, off) = match &buf[..8] {
        m if m == MAGIC => (buf.get(16).copied().unwrap_or(0) & 1 == 1, 25usize),
        m if m == MAGIC_V1 => (false, 24usize),
        _ => return Err(bad("not a marks file")),
    };
    let every = u64::from_le_bytes(buf[8..16].try_into().unwrap());
    if buf.len() < off { return Err(bad("marks file is truncated")); }
    let count = u64::from_le_bytes(buf[off - 8..off].try_into().unwrap()) as usize;
    if buf.len() != off + 16 * count { return Err(bad("marks file is truncated")); }
    let n = (0..count).map(|k| u128::from_le_bytes(buf[off + 16 * k..off + 16 * k + 16].try_into().unwrap())).collect();
    Ok((every, sparse, n))
}

impl Marks {
    /// Sieve the first n_primes primes once, keep every `every`-th as a packed N.
    pub fn build(n_primes: u64, every: u64) -> Self {
        let nf = n_primes as f64;
        let limit = if n_primes > 10 { (nf * (nf.ln() + nf.ln().ln())) as u64 + 100 } else { 100 };
        let t = PiTable::new(limit);
        let family = MARKS_FAMILY;
        let (mut n, mut p) = (Vec::new(), Vec::new());
        // prime 2 is index 0; then scan the odd bitset word by word
        let mut i = 0u64;
        if n_primes > 0 { n.push(family.pack(0, 2)); p.push(2); i = 1; }
        'outer: for (w, &word) in t.words.iter().enumerate() {
            let mut bits = word;
            while bits != 0 {
                let b = bits.trailing_zeros() as u64;
                bits &= bits - 1;
                if i >= n_primes { break 'outer; }
                if i % every == 0 {
                    let q = 2 * (w as u64 * 64 + b) + 1;
                    n.push(family.pack(i, q));
                    p.push(q);
                }
                i += 1;
            }
        }
        let last_prime = Self::walk_last(&n, every, family);
        Marks { n, every, family, p, last_prime }
    }
    fn walk_last(n: &[u128], every: u64, family: Family) -> u64 {
        let (_, mut p) = family.unpack(*n.last().unwrap());
        for _ in 1..every { p = next_prime(p); }
        p
    }
    pub fn len(&self) -> usize { self.n.len() }
    /// Full check: every mark is prime. Slow on big tables; load() spot-checks.
    pub fn verify(&self) -> bool { self.p.iter().all(|&q| isprime(q)) }
    pub fn is_empty(&self) -> bool { self.n.is_empty() }
    /// Highest index the table can answer without walking past its last mark.
    pub fn last_index(&self) -> u64 { self.family.unpack(*self.n.last().unwrap()).0 + self.every - 1 }

    pub fn save(&self, path: &str) -> std::io::Result<()> {
        write_marks_file(path, self.every, false, &self.n)
    }
    pub fn load(path: &str) -> std::io::Result<Self> {
        let (every, sparse, n) = read_marks_file(path)?;
        let bad = |m: &str| std::io::Error::new(std::io::ErrorKind::InvalidData, m.to_string());
        if sparse { return Err(bad("sparse marks file (LiveMarks::save): load it with LiveMarks::load")); }
        if n.is_empty() { return Err(bad("empty marks file")); }
        let family = MARKS_FAMILY;
        for (k, &x) in n.iter().enumerate() {
            if family.unpack(x).0 != k as u64 * every { return Err(bad("marks are not contiguous at the stated spacing")); }
        }
        let p: Vec<u64> = n.iter().map(|&x| family.unpack(x).1).collect();
        if p.windows(2).any(|w| w[0] >= w[1]) { return Err(bad("marks are not increasing")); }
        let step = (n.len() / 64).max(1);
        if (0..n.len()).step_by(step).any(|k| !isprime(p[k])) { return Err(bad("a mark is not prime")); }
        let last_prime = Self::walk_last(&n, every, family);
        Ok(Marks { n, every, family, p, last_prime })
    }

    /// Prime at zero-based index i: unpack the mark below, walk the rest.
    pub fn nth_prime(&self, i: u64) -> u64 {
        assert!(i <= self.last_index(), "index {i} is past the table (last {})", self.last_index());
        let (mi, mut p) = self.family.unpack(self.n[(i / self.every) as usize]);
        for _ in mi..i { p = next_prime(p); }
        p
    }
    /// Zero-based index of prime p: unpack the mark below, walk up counting.
    pub fn index_of(&self, p: u64) -> u64 {
        assert!(isprime(p), "{p} is not prime");
        assert!(p >= 2 && p <= self.last_prime, "{p} is past the table");
        let k = self.p.partition_point(|&m| m <= p) - 1;
        let (mut i, mut q) = self.family.unpack(self.n[k]);
        while q < p { q = next_prime(q); i += 1; }
        i
    }
}

// ── live marks: the same table, filled on demand ────────────────────────────

/// Marks that compute themselves. Sparse: only the blocks you have touched
/// exist. A miss costs one counting call (Riemann guess, Meissel pin, walk);
/// every later lookup in that block is a walk from the stored mark.
/// Same walk, same pack, same file format as `Marks`; the vector is replaced
/// by a map so gaps are allowed.
pub struct LiveMarks { pub every: u64, pub family: Family, marks: HashMap<u64, u128>, counter: Counter }

impl LiveMarks {
    /// `limit` sizes the counter: it can count any x with x^(2/3) <= limit
    /// (1e8 covers primes to about 1e12).
    pub fn new(every: u64, limit: u64) -> Self {
        LiveMarks { every, family: MARKS_FAMILY, marks: HashMap::new(), counter: Counter::new(limit) }
    }
    /// Seed from a built table so its marks never have to be computed.
    /// Seed from a built table: only its marks that fall on this grid are taken,
    /// so two source marks can never land in one block.
    pub fn with_marks(mut self, m: &Marks) -> Self {
        for &n in &m.n {
            let (i, _) = self.family.unpack(n);
            if i % self.every == 0 { self.marks.insert(i / self.every, n); }
        }
        self
    }
    pub fn len(&self) -> usize { self.marks.len() }
    pub fn is_empty(&self) -> bool { self.marks.is_empty() }
    /// The mark for block k, computing and storing it if absent.
    pub fn mark(&mut self, k: u64) -> (u64, u64) {
        if let Some(&n) = self.marks.get(&k) { return self.family.unpack(n); }
        let i = k * self.every;
        let p = self.counter.nth_prime(i);
        self.marks.insert(k, self.family.pack(i, p));
        (i, p)
    }
    pub fn nth_prime(&mut self, i: u64) -> u64 {
        let (mi, mut p) = self.mark(i / self.every);
        for _ in mi..i { p = next_prime(p); }
        p
    }
    pub fn index_of(&mut self, p: u64) -> u64 {
        assert!(isprime(p), "{p} is not prime");
        self.pi(p) - 1
    }
    /// pi(x) for any x: settle on the block holding x with a few marks, then walk.
    pub fn pi(&mut self, x: u64) -> u64 {
        if x < 2 { return 0; }
        // guess the block from R, then jump by the measured distance (in primes,
        // via the local gap ln x) until the block below x is found
        let lnx = (x as f64).ln().max(1.0);
        let mut k = (riemann_r(x as f64).max(1.0) as u64) / self.every;
        loop {
            let (_, p) = self.mark(k);
            if p > x {
                let back = (((p - x) as f64 / lnx / self.every as f64) as u64).max(1);
                k = k.saturating_sub(back);
                continue;
            }
            let (_, pn) = self.mark(k + 1);
            if pn <= x {
                let fwd = (((x - pn) as f64 / lnx / self.every as f64) as u64).max(1);
                k += fwd;
                continue;
            }
            let (i, mut q) = self.mark(k);
            let mut count = i + 1;
            loop { q = next_prime(q); if q > x { return count; } count += 1; }
        }
    }
    pub fn save(&self, path: &str) -> std::io::Result<()> {
        let mut keys: Vec<_> = self.marks.keys().copied().collect();
        keys.sort_unstable();
        let n: Vec<u128> = keys.iter().map(|k| self.marks[k]).collect();
        write_marks_file(path, self.every, true, &n)
    }
    /// Load a sparse or contiguous marks file. `limit` sizes the counter as in `new`.
    pub fn load(path: &str, limit: u64) -> std::io::Result<Self> {
        let (every, _sparse, n) = read_marks_file(path)?;
        let bad = |m: &str| std::io::Error::new(std::io::ErrorKind::InvalidData, m.to_string());
        let mut lm = LiveMarks::new(every, limit);
        for x in n {
            let (i, p) = lm.family.unpack(x);
            if i % every != 0 { return Err(bad("a mark is off the spacing grid")); }
            if !isprime(p) { return Err(bad("a mark is not prime")); }
            lm.marks.insert(i / every, x);
        }
        Ok(lm)
    }
}

impl Marks {
    /// pi(x) for any x the table covers: walk from the mark below x, counting.
    pub fn pi(&self, x: u64) -> u64 {
        if x < 2 { return 0; }
        assert!(x <= self.last_prime, "{x} is past the table");
        let k = self.p.partition_point(|&m| m <= x) - 1;
        let (i, mut q) = self.family.unpack(self.n[k]);
        let mut count = i + 1;
        loop { q = next_prime(q); if q > x { return count; } count += 1; }
    }
}

// ── the anchored window ─────────────────────────────────────────────────────

/// Offsets d in [0, width) with anchor + d odd and free of every seed factor.
/// One pass per seed: first d is (-anchor) mod q, then every q. No divisions
/// per candidate. Every prime survives. Ascending, so closest to anchor first.
pub fn window(anchor: u64, width: usize) -> Vec<u64> {
    let mut alive = vec![true; width];
    let mut d = (anchor % 2) as usize;               // even candidates out
    while d < width { alive[d] = false; d += 2; }
    for &q in &SEEDS {
        let first = ((q - anchor % q) % q) as usize;
        let mut d = first;
        while d < width { alive[d] = false; d += q as usize; }
        // the seed itself is prime: it sits at offset q - anchor, which is the first stride
        // only when the anchor is not a multiple of q (at anchor 0 it is the second)
        if anchor <= q && q - anchor < width as u64 { alive[(q - anchor) as usize] = true; }
    }
    alive.iter().enumerate().filter(|(_, &a)| a).map(|(d, _)| d as u64).collect()
}

/// Largest prime gap below 2^64 (exhaustively verified). A first window of
/// this plus slack always holds the next prime.
pub const MAX_GAP_BELOW_2_64: u64 = 1550;
pub const WINDOW_WIDTH_U64: usize = MAX_GAP_BELOW_2_64 as usize + 50;

/// Closest prime above anchor: sieve the window, certify survivors in order.
/// The window is the largest known gap plus slack, so the answer is inside.
pub fn next_prime_after(anchor: u64) -> u64 {
    if anchor < 2 { return 2; }
    let mut base = anchor + 1;
    loop {
        for d in window(base, WINDOW_WIDTH_U64) {
            if isprime(base + d) { return base + d; }
        }
        base += WINDOW_WIDTH_U64 as u64;
    }
}

/// Every prime p with lo <= p < hi. No primality test: the seeds are all primes
/// up to sqrt(hi), so every survivor is prime. Sieve of Eratosthenes started
/// at lo, in segments, so memory is one segment regardless of the range.
pub fn primes_between(lo: u64, hi: u64) -> Vec<u64> {
    const SEG: u64 = 1 << 22;
    let mut out = Vec::new();
    if hi <= lo { return out; }
    let root = (hi as f64).sqrt() as u64 + 2;
    let seeds = PiTable::new(root).primes_upto(root);
    let mut base = lo;
    while base < hi {
        let width = (hi - base).min(SEG) as usize;
        let mut alive = vec![true; width];
        for &q in &seeds {
            if q * q >= base + width as u64 { break; }
            let first = if q * q >= base { q * q - base } else { (q - base % q) % q } as usize;
            let mut d = first;
            while d < width { alive[d] = false; d += q as usize; }
        }
        for (d, &a) in alive.iter().enumerate() {
            let n = base + d as u64;
            if a && n >= 2 { out.push(n); }
        }
        base += width as u64;
    }
    out
}

/// The next n primes above `anchor` with exact zero-based indices and packed N.
/// `count` is pi(anchor); pass it from any counter (primecount, Counter::primepi,
/// a mark). Count once, sieve, increment: the sieve returns primes in order so
/// the index just ticks.
pub fn indexed_after(anchor: u64, n: usize, count: u64, family: Family) -> Vec<(u64, u64, u128)> {
    let mut out = Vec::with_capacity(n);
    let mut lo = anchor + 1;
    let width = ((1.5 * n as f64 * (anchor as f64 + 2.0).ln()) as u64).max(4096);
    let mut idx = count;
    while out.len() < n {
        for p in primes_between(lo, lo + width) {
            out.push((idx, p, family.pack(idx, p)));
            idx += 1;
            if out.len() == n { break; }
        }
        lo += width;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pack_round_trips() {
        assert_eq!(F513.encode(78501, 1000039), (2000590, 19669454));
        let n = F513.pack(78501, 1000039);
        assert_eq!(n, 8592468642374094);
        assert_eq!(F513.unpack(n), (78501, 1000039));
        assert_eq!(F257.unpack(F257.pack(78501, 1000039)), (78501, 1000039));
    }

    #[test]
    fn rulers_are_ladders() {
        for i in 0..20000u64 {
            assert_eq!(rp(i), (1u64 << ((i + 1).trailing_zeros() + 3)) - 1);
            assert_eq!(ruler(i, 4, 64), 256 * ladder(i >> 6));
            assert_eq!(ruler(i, 256, 2), 512 * ladder(i >> 1));
        }
    }

    #[test]
    fn isprime_matches_sieve() {
        let t = PiTable::new(200_000);
        for n in 0..200_000u64 { assert_eq!(isprime(n), t.is_prime(n), "{n}"); }
        assert!(isprime((1u64 << 61) - 1));
        assert!(!isprime(3825123056546413051));
        assert_eq!(t.pi(199_999), 17984);
    }

    #[test]
    fn meissel_known_values() {
        let c = Counter::new(1_000_000);
        assert_eq!(c.primepi(1_000_039), 78502);
        assert_eq!(c.primepi(1_000_000_000), 50847534);
        assert_eq!(c.nth_prime(999_999), 15485863);
        assert_eq!(c.index_of(1000039), 78501);
    }

    #[test]
    fn marks_walk_matches_counting() {
        let m = Marks::build(300_000, 64);
        let c = Counter::new(1_000_000);
        for i in [0u64, 1, 2, 63, 64, 1000, 78501, 299_999] {
            let p = m.nth_prime(i);
            assert_eq!(p, c.nth_prime(i), "index {i}");
            assert_eq!(m.index_of(p), i);
        }
        assert_eq!(m.nth_prime(78501), 1000039);
    }

    #[test]
    fn live_marks_fill_on_demand() {
        let mut lm = LiveMarks::new(64, 1_000_000);
        assert_eq!(lm.nth_prime(78501), 1000039);
        assert_eq!(lm.len(), 1);
        assert_eq!(lm.nth_prime(78502), 1000081);     // same block, no new mark
        assert_eq!(lm.len(), 1);
        assert_eq!(lm.index_of(1000039), 78501);
        assert_eq!(lm.pi(1_000_000), 78498);
        assert_eq!(lm.pi(10), 4);
        let m = Marks::build(300_000, 64);
        assert_eq!(m.pi(1_000_000), 78498);
        assert_eq!(m.pi(2), 1);
    }

    #[test]
    fn window_at_anchor_zero_keeps_the_seeds() {
        assert_eq!(window(0, 20), vec![1, 3, 5, 7, 11, 13, 17, 19]);
        assert_eq!(window(3, 10), vec![0, 2, 4, 8]);
    }

    #[test]
    fn window_keeps_every_prime() {
        let t = PiTable::new(2_000_000);
        for anchor in [1_000_000u64, 1 << 20, 1_299_700] {
            let got = window(anchor, 50_000);
            let mut k = 0;
            for n in anchor..anchor + 50_000 {
                if t.is_prime(n) { assert_eq!(anchor + got[k], n); k += 1; }
                else if k < got.len() && anchor + got[k] == n { k += 1; }
            }
        }
        assert_eq!(next_prime_after(1_000_000_000_000), 1_000_000_000_039);
        assert_eq!(next_prime_after(1 << 44), (1 << 44) + 7);
        assert_eq!(next_prime_after(1), 2);
        assert_eq!(next_prime_after(2), 3);
        assert!(window(1_000_000_000_000, 100_000).len() < 20_000);
    }

    #[test]
    fn primes_between_matches_table() {
        let t = PiTable::new(3_000_000);
        for (lo, hi) in [(0u64, 1000u64), (1000, 50_000), (2_000_000, 2_600_000), (131, 132), (1_000_000, 1_000_001), (2, 3)] {
            let want: Vec<u64> = (lo..hi).filter(|&n| n >= 2 && t.is_prime(n)).collect();
            assert_eq!(primes_between(lo, hi), want, "{lo}..{hi}");
        }
        assert_eq!(primes_between(1_000_000_000_000, 1_000_000_000_100), vec![1000000000039, 1000000000061, 1000000000063, 1000000000091]);
        assert!(primes_between(10, 10).is_empty());
    }

    #[test]
    fn review_findings() {
        // 1. no silent overflow: a prime above 2^63 and an index past 1e11 round-trip
        // N holds K and X, each about bits(2p) wide, so u128 carries pairs with p < 2^63
        let big63 = 9223372036854775783u64;                 // largest prime below 2^63
        assert_eq!(MARKS_FAMILY.unpack(MARKS_FAMILY.pack(0, big63)), (0, big63));
        assert_eq!(MARKS_FAMILY.unpack(MARKS_FAMILY.pack(100_000_000_000, big63)), (100_000_000_000, big63));
        // a prime above 2^63 cannot fit: refused loudly, never wrapped
        let big64 = 18446744073709551557u64;
        assert!(std::panic::catch_unwind(|| MARKS_FAMILY.pack(0, big64)).is_err());
        assert_eq!(F513.unpack(F513.pack(78501, 1000039)), (78501, 1000039));
        // 2. seeding from a finer table takes only on-grid marks
        let fine = Marks::build(4000, 16);
        let mut lm = LiveMarks::new(64, 1_000_000).with_marks(&fine);
        assert_eq!(lm.len(), fine.len().div_ceil(4));
        assert_eq!(lm.nth_prime(0), 2);
        assert_eq!(lm.nth_prime(3999), fine.nth_prime(3999));
        // 3. sparse files: LiveMarks::load accepts, Marks::load refuses with an error, no panic
        let path = std::env::temp_dir().join("rulerpack_sparse.marks");
        let s = path.to_str().unwrap();
        let mut live = LiveMarks::new(64, 1_000_000);
        live.nth_prime(78501);
        live.save(s).unwrap();
        assert!(Marks::load(s).is_err());
        let mut back = LiveMarks::load(s, 1_000_000).unwrap();
        assert_eq!(back.len(), 1);
        assert_eq!(back.nth_prime(78501), 1000039);
        // one-mark contiguous file round-trips
        Marks::build(10, 1024).save(s).unwrap();
        assert_eq!(Marks::load(s).unwrap().nth_prime(9), 29);
        std::fs::remove_file(path).ok();
    }

    #[test]
    fn indexed_after_matches_counter() {
        let c = Counter::new(1_000_000);
        let rows = indexed_after(1_000_000, 5, c.primepi(1_000_000), MARKS_FAMILY);
        assert_eq!(rows[0].0, 78498);
        assert_eq!(rows.iter().map(|r| r.1).collect::<Vec<_>>(), vec![1000003, 1000033, 1000037, 1000039, 1000081]);
        for (i, p, n) in &rows { assert_eq!(MARKS_FAMILY.unpack(*n), (*i, *p)); }
        let rows = indexed_after(1_000_000_000_000, 2, 37607912018, MARKS_FAMILY);
        assert_eq!((rows[0].0, rows[0].1), (37607912018, 1000000000039));
    }

    #[test]
    fn marks_save_load() {
        let m = Marks::build(100_000, 256);
        let path = std::env::temp_dir().join("rulerpack_test.marks");
        m.save(path.to_str().unwrap()).unwrap();
        let m2 = Marks::load(path.to_str().unwrap()).unwrap();
        assert_eq!(m2.n, m.n);
        assert_eq!(m2.every, 256);
        assert_eq!(m2.nth_prime(99_999), 1299709);
        std::fs::remove_file(path).ok();
    }
}
