//! **The special number field sieve** (`big` feature), small and readable: 2^k − 1.
//!
//! The polynomial comes from the form: with t = ⌈k/d⌉, m = 2^t and f(x) = x^d − 2^(dt−k),
//! f(m) ≡ 0 (mod 2^k − 1), monic with one tiny coefficient. Then the smooth window runs twice
//! over a plane of coprime pairs (a, b): the rational side a + b·m along a ≡ −b·m (mod p), the
//! algebraic side N(a + b·α) along one stride per ideal (p, α − r), a ≡ −b·r (mod p). A GF(2)
//! dependency makes both products squares; the rational root is read off the exponents and the
//! algebraic root comes by Couveignes' method for odd degree: a square root in each of several
//! prime fields F_q[x]/(f), the sign fixed by the norm, the coefficients put together by the
//! Chinese remainder theorem. Map α → m, gcd(X − Y, n). Same parts as `rulerpack.snfs` in the
//! Python package, in the same order.

use num_bigint::{BigInt, BigUint};
use num_integer::Integer;
use num_traits::{One, Signed, ToPrimitive, Zero};

use crate::qs::{dependencies, legendre};
use crate::smooth::{log_accumulate, primes_upto};

/// `(f, m)` for 2^k − 1: f low to high, monic degree d, f(m) ≡ 0 (mod n).
pub fn polynomial(k: u32, d: u32) -> (Vec<BigInt>, BigUint) {
    let t = k.div_ceil(d);
    let c = BigInt::one() << (d * t - k) as usize;
    let mut f = vec![-c];
    f.extend((1..d).map(|_| BigInt::zero()));
    f.push(BigInt::one());
    (f, BigUint::one() << t as usize)
}

/// f(x) mod q for a small prime q, f's coefficients reduced first.
fn f_mod(f: &[BigInt], q: u64) -> Vec<u64> {
    f.iter().map(|c| { let r = c.mod_floor(&BigInt::from(q)); r.to_u64().unwrap() }).collect()
}
fn eval_mod(fq: &[u64], x: u64, q: u64) -> u64 {
    let mut r = 0u128;
    for c in fq.iter().rev() { r = (r * x as u128 + *c as u128) % q as u128; }
    r as u64
}
/// Roots of f modulo a small prime, by trial.
pub fn roots_mod(fq: &[u64], q: u64) -> Vec<u64> { (0..q).filter(|&r| eval_mod(fq, r, q) == 0).collect() }

/// N(a + b·α) = (−b)^d·f(−a/b) for monic f: for x³ − 4 that is a³ + 4b³.
pub fn norm(f: &[BigInt], a: i64, b: i64) -> BigInt {
    let d = f.len() - 1;
    let (a, b) = (BigInt::from(a), BigInt::from(b));
    let mut s = BigInt::zero();
    for (i, c) in f.iter().enumerate() { s += c * (-&a).pow(i as u32) * b.pow((d - i) as u32); }
    if d % 2 == 1 { -s } else { s }
}

/// A list of `(prime, residue)` pairs: the rational base `(p, m mod p)`, the algebraic base
/// `(p, r)` per ideal, the quadratic characters `(q, s)`.
pub type Base = Vec<(u64, u64)>;

/// The rational base, the algebraic base per ideal, and the quadratic characters beyond it.
pub fn bases(f: &[BigInt], m: &BigUint, rational_bound: u64, algebraic_bound: u64, characters: usize) -> (Base, Base, Base) {
    let rat: Vec<(u64, u64)> = primes_upto(rational_bound).into_iter().map(|p| (p, (m % p).to_u64().unwrap())).collect();
    let mut alg = Vec::new();
    for p in primes_upto(algebraic_bound) { let fq = f_mod(f, p); for r in roots_mod(&fq, p) { alg.push((p, r)); } }
    let mut chars = Vec::new();
    let mut q = algebraic_bound + 1;
    while chars.len() < characters {
        if crate::isprime(q) {
            let fq = f_mod(f, q);
            for s in roots_mod(&fq, q) {
                // a simple root: f'(s) ≠ 0
                let mut fp = 0u128;
                for (i, c) in fq.iter().enumerate().skip(1) { fp = (fp + (i as u128) * (*c as u128) % q as u128 * modpow(s, i as u64 - 1, q) as u128) % q as u128; }
                if fp != 0 { chars.push((q, s)); }
            }
        }
        q += 1;
    }
    chars.truncate(characters);
    (rat, alg, chars)
}

fn modpow(b: u64, mut e: u64, m: u64) -> u64 {
    let (mut r, mut bb, mm) = (1u128, (b % m) as u128, m as u128);
    while e > 0 { if e & 1 == 1 { r = r * bb % mm; } bb = bb * bb % mm; e >>= 1; }
    r as u64
}

/// Exponents of |v| over `primes` and the cofactor, 1 when smooth.
fn divide_out(v: &BigInt, primes: &[u64]) -> (Vec<u32>, BigUint) {
    let mut v = v.magnitude().clone();
    let mut e = vec![0u32; primes.len()];
    for (i, &p) in primes.iter().enumerate() { while (&v % p).is_zero() { v /= p; e[i] += 1; } }
    (e, v)
}

/// One relation, full or merged from two partials: the pairs it is made of, the summed
/// rational and ideal exponents, the rational sign, and the large primes outside the bases with
/// their exponents (always even after a merge: `(p, 2)` on the rational side, `(p, r, 2)` for an
/// ideal on the algebraic side), which the square roots must still account for.
#[derive(Clone, Debug)]
pub struct Relation {
    pub pairs: Vec<(i64, i64)>,
    pub rat: Vec<u32>,
    pub alg: Vec<u32>,
    pub negative: bool,
    pub large_rat: Vec<(u64, u32)>,
    pub large_alg: Vec<(u64, u64, u32)>,
}

/// A pair's sieve result before merging: full, or partial with one large prime on one side.
enum Found { Full(Relation), Partial(LargeKey, Relation) }
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
enum LargeKey { Rat(u64), Alg(u64, u64) }

/// Everything a line of the sieve needs, shared read-only across the threads.
struct Sieve<'a> { f: &'a [BigInt], mb: BigInt, m_f: f64, rat: &'a [(u64, u64)], alg: &'a [(u64, u64)], rat_primes: Vec<u64>, alg_primes: Vec<u64>, half_width: i64, cut_r: f32, cut_a: f32, large: u64, large_log: f32, coupled: bool }

impl Sieve<'_> {
/// Sieve one line of b: every coprime a with |a| ≤ half_width, both sides smooth up to one
/// large prime below `large` on one side.
fn line(&self, b: i64) -> Vec<Found> {
    let Sieve { f, mb, m_f, rat, alg, rat_primes, alg_primes, half_width, cut_r, cut_a, large, large_log, coupled } = self;
    let (f, rat, alg, m_f, half_width, cut_r, cut_a, large, large_log, coupled) = (*f, *rat, *alg, *m_f, *half_width, *cut_r, *cut_a, *large, *large_log, *coupled);
    let d = (f.len() - 1) as i32;
    let width = (2 * half_width + 1) as usize;
    let c = -f[0].to_f64().unwrap();
    let strides_r: Vec<(u64, Vec<usize>)> = rat.iter().map(|&(p, mp)| (p, vec![((-(b as i128) * mp as i128 + half_width as i128).rem_euclid(p as i128)) as usize])).collect();
    let mut strides_a: Vec<(u64, Vec<usize>)> = Vec::new();
    for &(p, r) in alg {
        let s = ((-(b as i128) * r as i128 + half_width as i128).rem_euclid(p as i128)) as usize;
        match strides_a.last_mut() { Some((lp, v)) if *lp == p => v.push(s), _ => strides_a.push((p, vec![s])) }
    }
    let bf = b as f64;
    // the shortfall on each side: how far the logs found fall under the size, less the base
    // slack. Independent screening gives each side the large-prime allowance on its own; the
    // coupled rule (the disjoint-residue construction) gives it to one side or the other, never both,
    // since a relation may carry one large prime in total
    let acc_r = log_accumulate(width, &strides_r);
    let acc_a = log_accumulate(width, &strides_a);
    let short = |i: usize| {
        let a = (i as i64 - half_width) as f64;
        let er = ((a + bf * m_f).abs().max(1.0).ln() as f32) - cut_r - acc_r[i];
        let ea = ((a.powi(d) + c * bf.powi(d)).abs().max(1.0).ln() as f32) - cut_a - acc_a[i];
        (er, ea)
    };
    let large_b = BigUint::from(large);
    let mut out = Vec::new();
    for i in 0..width {
        let (er, ea) = short(i);
        let pass = if coupled { (er <= large_log && ea <= 0.0) || (er <= 0.0 && ea <= large_log) } else { er <= large_log && ea <= large_log };
        if !pass { continue; }
        let a = i as i64 - half_width;
        if a == 0 || a.gcd(&b) != 1 { continue; }
        let rv = BigInt::from(a) + BigInt::from(b) * mb;
        let (re, rr) = divide_out(&rv, rat_primes);
        let rat_large = if rr.is_one() { None } else if rr < large_b { Some(rr.to_u64().unwrap()) } else { continue };
        let nv = norm(f, a, b);
        if nv.is_zero() { continue; }
        let (pe, ar) = divide_out(&nv, alg_primes);
        let alg_large = if ar.is_one() { None } else if ar < large_b && rat_large.is_none() { Some(ar.to_u64().unwrap()) } else { continue };
        let Some(ae) = ideal_exponents_from(a, b, &pe, alg, alg_primes) else { continue };
        let mut rel = Relation { pairs: vec![(a, b)], rat: re, alg: ae, negative: rv.is_negative(), large_rat: vec![], large_alg: vec![] };
        match (rat_large, alg_large) {
            (None, None) => out.push(Found::Full(rel)),
            (Some(p), None) => { rel.large_rat.push((p, 1)); out.push(Found::Partial(LargeKey::Rat(p), rel)); }
            (None, Some(p)) => {
                // the ideal above the large prime: r = −a·b⁻¹ (mod p); p cannot divide b
                if (b as i128).rem_euclid(p as i128) == 0 { continue; }
                let binv = modpow((b as i128).rem_euclid(p as i128) as u64, p - 2, p);
                let r = ((-(a as i128)).rem_euclid(p as i128) as u128 * binv as u128 % p as u128) as u64;
                rel.large_alg.push((p, r, 1));
                out.push(Found::Partial(LargeKey::Alg(p, r), rel));
            }
            (Some(_), Some(_)) => {}
        }
    }
    out
}
}

/// `ideal_exponents` from the norm's prime exponents already in hand.
fn ideal_exponents_from(a: i64, b: i64, pe: &[u32], alg: &[(u64, u64)], alg_primes: &[u64]) -> Option<Vec<u32>> {
    let mut left: std::collections::HashMap<u64, u32> = alg_primes.iter().copied().zip(pe.iter().copied()).collect();
    let mut out = vec![0u32; alg.len()];
    for (i, &(p, r)) in alg.iter().enumerate() {
        let e = left[&p];
        if e > 0 && (a as i128 + b as i128 * r as i128).rem_euclid(p as i128) == 0 { out[i] = e; left.insert(p, 0); }
    }
    if left.values().any(|&e| e > 0) { return None; }
    Some(out)
}

/// Two partials with the same large prime, multiplied: the prime appears squared.
fn merge(x: &Relation, y: &Relation) -> Relation {
    let add = |u: &[u32], v: &[u32]| u.iter().zip(v).map(|(a, b)| a + b).collect::<Vec<u32>>();
    let mut pairs = x.pairs.clone(); pairs.extend(&y.pairs);
    let mut large_rat = x.large_rat.clone();
    for &(p, e) in &y.large_rat { match large_rat.iter_mut().find(|(q, _)| *q == p) { Some(t) => t.1 += e, None => large_rat.push((p, e)) } }
    let mut large_alg = x.large_alg.clone();
    for &(p, r, e) in &y.large_alg { match large_alg.iter_mut().find(|(q, s, _)| *q == p && *s == r) { Some(t) => t.2 += e, None => large_alg.push((p, r, e)) } }
    Relation { pairs, rat: add(&x.rat, &y.rat), alg: add(&x.alg, &y.alg), negative: x.negative ^ y.negative, large_rat, large_alg }
}

/// Relations until `need` are in hand: full ones, plus one per extra partial sharing a large
/// prime. Lines of b run in parallel, `batch` at a time. Returns the relations and the lines sieved.
pub fn relations(f: &[BigInt], m: &BigUint, rat: &[(u64, u64)], alg: &[(u64, u64)], need: usize, bounds: &Bounds, slack: f32) -> (Vec<Relation>, u64) {
    use rayon::prelude::*;
    let rat_primes: Vec<u64> = rat.iter().map(|&(p, _)| p).collect();
    let mut alg_primes: Vec<u64> = alg.iter().map(|&(p, _)| p).collect(); alg_primes.dedup();
    let pmax_r = *rat_primes.last().unwrap();
    let pmax_a = *alg_primes.last().unwrap();
    // one large prime below `large` is allowed on one side; it must be below pmax² to be prime
    let large_factor = bounds.large_factor.max(1);
    let large = (pmax_r.min(pmax_a) * large_factor).min(pmax_r.min(pmax_a).pow(2));
    let lp = if large_factor > 1 { (large as f32).ln() } else { 0.0 };
    let sieve = Sieve { f, mb: BigInt::from(m.clone()), m_f: m.to_f64().unwrap(), rat, alg, rat_primes, alg_primes, half_width: bounds.half_width, cut_r: slack * (pmax_r as f32).ln(), cut_a: slack * (pmax_a as f32).ln(), large, large_log: lp, coupled: bounds.coupled };
    let mut full: Vec<Relation> = Vec::new();
    let mut partials: std::collections::HashMap<LargeKey, Relation> = std::collections::HashMap::new();
    let batch = rayon::current_num_threads().max(4) as i64;
    let mut b0 = 1i64;
    while full.len() < need {
        let found: Vec<Vec<Found>> = (b0..b0 + batch).into_par_iter()
            .map(|b| sieve.line(b))
            .collect();
        for line in found {
            for item in line {
                match item {
                    Found::Full(r) => full.push(r),
                    Found::Partial(key, r) => match partials.get(&key) {
                        Some(first) => full.push(merge(first, &r)),
                        None => { partials.insert(key, r); }
                    },
                }
            }
        }
        b0 += batch;
    }
    full.truncate(need);
    (full, (b0 - 1) as u64)
}

// ---- F_q[x]/(f) with u64 coefficients, q below 2^32 so products fit u128 comfortably

fn mulmod(u: &[u64], v: &[u64], fq: &[u64], q: u64) -> Vec<u64> {
    let d = fq.len() - 1;
    let mut prod = vec![0u128; 2 * d - 1];
    for (i, &ui) in u.iter().enumerate() { if ui != 0 { for (j, &vj) in v.iter().enumerate() { prod[i + j] = (prod[i + j] + ui as u128 * vj as u128) % q as u128; } } }
    for k in (d..2 * d - 1).rev() {
        let c = prod[k];
        if c != 0 { for i in 0..=d { let sub = c * fq[i] as u128 % q as u128; prod[k - d + i] = (prod[k - d + i] + q as u128 - sub) % q as u128; } }
    }
    prod[..d].iter().map(|&x| x as u64).collect()
}
fn powmod(u: &[u64], e: &BigUint, fq: &[u64], q: u64) -> Vec<u64> {
    let d = fq.len() - 1;
    let mut r = vec![0u64; d]; r[0] = 1;
    let mut base = u.to_vec();
    for bit in 0..e.bits() { if e.bit(bit) { r = mulmod(&r, &base, fq, q); } base = mulmod(&base, &base, fq, q); }
    r
}
fn field_norm(u: &[u64], fq: &[u64], q: u64) -> u64 {
    let d = (fq.len() - 1) as u32;
    let e = (BigUint::from(q).pow(d) - BigUint::one()) / BigUint::from(q - 1);
    let r = powmod(u, &e, fq, q);
    assert!(r[1..].iter().all(|&c| c == 0), "the norm is in the base field");
    r[0]
}
fn gcd_poly(mut u: Vec<u64>, mut v: Vec<u64>, q: u64) -> Vec<u64> {
    fn trim(w: &mut Vec<u64>) { while w.last() == Some(&0) { w.pop(); } }
    trim(&mut u); trim(&mut v);
    while !v.is_empty() {
        let inv = modpow(v[v.len() - 1], q - 2, q);
        while u.len() >= v.len() && !u.is_empty() {
            let c = u[u.len() - 1] as u128 * inv as u128 % q as u128;
            let shift = u.len() - v.len();
            for (i, &vc) in v.iter().enumerate() { let sub = c * vc as u128 % q as u128; u[shift + i] = ((u[shift + i] as u128 + q as u128 - sub) % q as u128) as u64; }
            trim(&mut u);
        }
        std::mem::swap(&mut u, &mut v);
    }
    u
}
/// Monic f irreducible over F_q: no root and, for degree 5, no quadratic factor.
pub fn irreducible_mod(fq: &[u64], q: u64) -> bool {
    let d = fq.len() - 1;
    let mut x = vec![0u64; d]; x[1] = 1;
    let mut cur = powmod(&x, &BigUint::from(q), fq, q);
    for _ in 1..=d / 2 {
        let mut mx = cur.clone(); mx[1] = (mx[1] + q - 1) % q;
        if gcd_poly(mx, fq.to_vec(), q).len() > 1 { return false; }
        cur = powmod(&cur, &BigUint::from(q), fq, q);
    }
    true
}

// ---- the square root for any degree: Hensel lifting (the "direct" method)
//
// Compute the product P = f'(α)²·∏(a + b·α) exactly in Z[α], take its square root in the
// field F_q[x]/(f) by Tonelli–Shanks, lift that root to a root modulo q^K by Newton's step
// (with the inverse of 2β lifted alongside) until q^K passes the size of β's coefficients,
// bring the coefficients into the symmetric range, evaluate at m. The factor f'(α)² puts the
// root inside Z[α] whatever the index of Z[α] in the ring of integers; it is cancelled on the
// rational side by f'(m). β(m)² ≡ f'(m)²·∏(a + b·m) (mod n) is the check.

fn poly_mul_exact(u: &[BigInt], v: &[BigInt], f: &[BigInt]) -> Vec<BigInt> {
    let d = f.len() - 1;
    let mut prod = vec![BigInt::zero(); 2 * d - 1];
    for (i, ui) in u.iter().enumerate() { if !ui.is_zero() { for (j, vj) in v.iter().enumerate() { prod[i + j] += ui * vj; } } }
    for k in (d..2 * d - 1).rev() {
        let c = std::mem::take(&mut prod[k]);
        if !c.is_zero() { for i in 0..d { let t = &c * &f[i]; prod[k - d + i] -= t; } }
    }
    prod.truncate(d);
    prod
}
fn poly_mul_mod_big(u: &[BigInt], v: &[BigInt], f: &[BigInt], md: &BigInt) -> Vec<BigInt> {
    poly_mul_exact(u, v, f).into_iter().map(|c| c.mod_floor(md)).collect()
}
fn poly_derivative(f: &[BigInt]) -> Vec<BigInt> {
    let d = f.len() - 1;
    let mut g: Vec<BigInt> = (1..=d).map(|i| &f[i] * BigInt::from(i)).collect();
    g.resize(d, BigInt::zero());
    g
}

/// A square root in F_q[x]/(f), f irreducible: Tonelli–Shanks in the field of q^d elements.
fn field_sqrt(u: &[u64], fq: &[u64], q: u64) -> Option<Vec<u64>> {
    let d = fq.len() - 1;
    let one = { let mut v = vec![0u64; d]; v[0] = 1; v };
    let order = BigUint::from(q).pow(d as u32) - BigUint::one();
    let half = &order / 2u32;
    if powmod(u, &half, fq, q) != one { return None; }
    let (mut t, mut s) = (order.clone(), 0u32);
    while (&t % 2u32).is_zero() { t /= 2u32; s += 1; }
    // a non-residue: x + c for the first c that is one
    let mut z = vec![0u64; d]; z[1] = 1;
    loop { if powmod(&z, &half, fq, q) != one { break; } z[0] += 1; }
    let mut c = powmod(&z, &t, fq, q);
    let mut r = powmod(u, &((&t + BigUint::one()) / 2u32), fq, q);
    let mut w = powmod(u, &t, fq, q);
    let mut m = s;
    while w != one {
        let mut i = 0u32; let mut w2 = w.clone();
        while w2 != one { w2 = mulmod(&w2, &w2, fq, q); i += 1; if i == m { return None; } }
        let b = powmod(&c, &(BigUint::one() << (m - i - 1) as usize), fq, q);
        r = mulmod(&r, &b, fq, q);
        c = mulmod(&b, &b, fq, q);
        w = mulmod(&w, &c, fq, q);
        m = i;
    }
    Some(r)
}

/// β(m) mod n with β² = f'(α)²·∏(a + b·α), any degree, by Hensel lifting; the caller's y must
/// carry f'(m) as well, which `factor` does for even degree.
pub fn hensel_sqrt(pairs: &[(i64, i64)], f: &[BigInt], m: &BigUint, n: &BigUint) -> Option<BigUint> {
    let d = f.len() - 1;
    // the product, exactly
    let mut prod: Vec<BigInt> = vec![BigInt::zero(); d]; prod[0] = BigInt::one();
    for &(a, b) in pairs { let mut e = vec![BigInt::zero(); d]; e[0] = BigInt::from(a); e[1] = BigInt::from(b); prod = poly_mul_exact(&prod, &e, f); }
    let fp = poly_derivative(f);
    let fp2 = poly_mul_exact(&fp, &fp, f);
    prod = poly_mul_exact(&prod, &fp2, f);
    let size_bits = prod.iter().map(|c| c.bits()).max().unwrap_or(0);
    // the root has about half the bits; lift to a modulus comfortably past that
    let want_bits = size_bits / 2 + 64 + 8 * d as u64;
    // a prime with f irreducible, q ≡ 1 mod anything is fine for Tonelli–Shanks
    let mut q = (1u64 << 30) + 3;
    let fq = loop { if crate::isprime(q) { let fq = f_mod(f, q); if irreducible_mod(&fq, q) { break fq; } } q += 2; };
    let prod_q: Vec<u64> = prod.iter().map(|c| c.mod_floor(&BigInt::from(q)).to_u64().unwrap()).collect();
    let root_q = field_sqrt(&prod_q, &fq, q)?;
    // lift: β ← β − (β² − P)/(2β), with t = 1/(2β) lifted by t ← t(2 − 2βt)
    let two = BigInt::from(2);
    let mut modulus = BigInt::from(q);
    let mut beta: Vec<BigInt> = root_q.iter().map(|&c| BigInt::from(c)).collect();
    let two_beta_q: Vec<u64> = beta.iter().map(|c| (c * &two).mod_floor(&BigInt::from(q)).to_u64().unwrap()).collect();
    // inverse of 2β in F_q by Fermat
    let inv_q = powmod(&two_beta_q, &(BigUint::from(q).pow(d as u32) - BigUint::from(2u32)), &fq, q);
    let mut t: Vec<BigInt> = inv_q.iter().map(|&c| BigInt::from(c)).collect();
    while modulus.bits() < want_bits {
        modulus = &modulus * &modulus;
        let beta2 = poly_mul_mod_big(&beta, &beta, f, &modulus);
        let diff: Vec<BigInt> = beta2.iter().zip(&prod).map(|(x, p)| (x - p).mod_floor(&modulus)).collect();
        let step = poly_mul_mod_big(&t, &diff, f, &modulus);
        beta = beta.iter().zip(&step).map(|(x, s)| (x - s).mod_floor(&modulus)).collect();
        let two_beta: Vec<BigInt> = beta.iter().map(|c| (c * &two).mod_floor(&modulus)).collect();
        let bt = poly_mul_mod_big(&two_beta, &t, f, &modulus);
        let mut two_minus = bt.iter().map(|c| (-c).mod_floor(&modulus)).collect::<Vec<_>>();
        two_minus[0] = (&two_minus[0] + &two).mod_floor(&modulus);
        t = poly_mul_mod_big(&t, &two_minus, f, &modulus);
    }
    let half = &modulus / 2;
    let nb = BigInt::from(n.clone());
    let mb = BigInt::from(m.clone());
    let mut value = BigInt::zero();
    for i in (0..d).rev() {
        let c = if beta[i] > half { &beta[i] - &modulus } else { beta[i].clone() };
        value = (value * &mb + c).mod_floor(&nb);
    }
    Some(value.to_biguint().unwrap())
}

/// β(m) mod n with β² = ∏(a + b·α), Couveignes' method, odd degree, monic f. `norm_extra` is
/// the part of |N(β)| from large primes outside the base (1 without the large-prime variation).
pub fn algebraic_sqrt(pairs: &[(i64, i64)], f: &[BigInt], m: &BigUint, n: &BigUint, exps_sum: &[u32], alg: &[(u64, u64)], norm_extra: &BigUint) -> Option<BigUint> {
    let d = f.len() - 1;
    if d % 2 == 0 { return hensel_sqrt(pairs, f, m, n); }
    let mut norm_beta = norm_extra.clone();
    for (&(p, _), &e) in alg.iter().zip(exps_sum) { if e % 2 == 1 { return None; } norm_beta *= BigUint::from(p).pow(e / 2); }
    let bound_bits: u64 = pairs.iter().map(|&(a, b)| (64 - a.abs().max(b.abs()).leading_zeros()) as u64 + 2).sum::<u64>() * d as u64 + 64;
    let (mut primes, mut bits) = (Vec::new(), 0u64);
    let mut q = (1u64 << 20) + 1;
    while bits < bound_bits {
        if crate::isprime(q) && q % 4 == 3 { let fq = f_mod(f, q); if irreducible_mod(&fq, q) { primes.push(q); bits += 20; } }
        q += 2;
    }
    let mut residues = Vec::new();
    for &q in &primes {
        let fq = f_mod(f, q);
        let mut prod = vec![0u64; d]; prod[0] = 1;
        for &(a, b) in pairs {
            let mut e = vec![0u64; d]; e[0] = (a as i128).rem_euclid(q as i128) as u64; e[1] = (b as i128).rem_euclid(q as i128) as u64;
            prod = mulmod(&prod, &e, &fq, q);
        }
        let e = (BigUint::from(q).pow(d as u32) + BigUint::one()) / BigUint::from(4u32);
        let mut root = powmod(&prod, &e, &fq, q);
        if mulmod(&root, &root, &fq, q) != prod { return None; }
        let want = (&norm_beta % q).to_u64().unwrap();
        if field_norm(&root, &fq, q) != want {
            root = root.iter().map(|&c| (q - c) % q).collect();
            if field_norm(&root, &fq, q) != want { return None; }
        }
        residues.push(root);
    }
    // CRT each coefficient into the symmetric range, then evaluate at m mod n
    let big_m: BigUint = primes.iter().fold(BigUint::one(), |acc, &q| acc * q);
    let half = &big_m / 2u32;
    let nb = BigInt::from(n.clone());
    let mb = BigInt::from(m.clone());
    let mut value = BigInt::zero();
    for i in (0..d).rev() {
        let mut x = BigUint::zero();
        for (&q, r) in primes.iter().zip(&residues) {
            let mi = &big_m / q;
            let inv = modpow((&mi % q).to_u64().unwrap(), q - 2, q);
            x += BigUint::from(r[i]) * &mi * inv;
        }
        x %= &big_m;
        let coeff = if x > half { BigInt::from(x) - BigInt::from(big_m.clone()) } else { BigInt::from(x) };
        value = (value * &mb + coeff).mod_floor(&nb);
    }
    Some(value.to_biguint().unwrap())
}

/// The sieve's knobs: the two factor-base bounds, half the width of a line in a, how many
/// relations past the base to collect, the large-prime bound as a multiple of the largest
/// base prime (1 turns the large-prime variation off), and whether the screening gives the
/// large-prime allowance to one side only (`coupled`) or to both independently.
#[derive(Clone, Copy, Debug)]
pub struct Bounds { pub rational: u64, pub algebraic: u64, pub half_width: i64, pub extra: usize, pub large_factor: u64, pub coupled: bool }

impl Bounds {
    /// Measured in the Python package: the bases must grow as bits^1.6, a linear rule
    /// collapses the yield.
    pub fn for_bits(bits: u32) -> Bounds {
        let b = (bits as f64).powf(1.6);
        Bounds { rational: (b as u64).max(300), algebraic: ((1.6 * b) as u64).max(500), half_width: (bits as i64 * 200).max(4000), extra: 16, large_factor: 40, coupled: true }
    }
}

/// A nontrivial factor of n from the polynomial (f, m), or None. `target` is the divisor of n
/// the gcd is taken against: n itself, or n with its small factors removed.
pub fn factor(n: &BigUint, target: &BigUint, f: &[BigInt], m: &BigUint, bounds: Bounds) -> Option<BigUint> {
    let (rat, alg, chars) = bases(f, m, bounds.rational, bounds.algebraic, 24);
    let need = 1 + rat.len() + alg.len() + chars.len() + bounds.extra;
    let (rels, _) = relations(f, m, &rat, &alg, need, &bounds, 1.6);
    let rows: Vec<Vec<u32>> = rels.iter().map(|r| {
        let mut row = vec![r.negative as u32];
        row.extend(&r.rat); row.extend(&r.alg);
        for &(q, s) in &chars {
            let mut bit = 0;
            for &(a, b) in &r.pairs { let v = (a as i128 + b as i128 * s as i128).rem_euclid(q as i128) as u64; if legendre(v, q) != 1 { bit ^= 1; } }
            row.push(bit);
        }
        row
    }).collect();
    for dep in dependencies(&rows) {
        if dep.is_empty() { continue; }
        let pairs: Vec<(i64, i64)> = dep.iter().flat_map(|&i| rels[i].pairs.iter().copied()).collect();
        let mut rat_sum = vec![0u32; rat.len()];
        let mut alg_sum = vec![0u32; alg.len()];
        let mut large_rat: std::collections::HashMap<u64, u32> = Default::default();
        let mut large_alg: std::collections::HashMap<(u64, u64), u32> = Default::default();
        for &i in &dep {
            for (j, e) in rels[i].rat.iter().enumerate() { rat_sum[j] += e; }
            for (j, e) in rels[i].alg.iter().enumerate() { alg_sum[j] += e; }
            for &(p, e) in &rels[i].large_rat { *large_rat.entry(p).or_default() += e; }
            for &(p, r, e) in &rels[i].large_alg { *large_alg.entry((p, r)).or_default() += e; }
        }
        if rat_sum.iter().chain(&alg_sum).any(|e| e % 2 == 1) || large_rat.values().chain(large_alg.values()).any(|e| e % 2 == 1) { continue; }
        let mut y = BigUint::one();
        for (&(p, _), &e) in rat.iter().zip(&rat_sum) { y = y * BigUint::from(p).modpow(&BigUint::from(e / 2), n) % n; }
        for (&p, &e) in &large_rat { y = y * BigUint::from(p).modpow(&BigUint::from(e / 2), n) % n; }
        let mut norm_extra = BigUint::one();
        for (&(p, _), &e) in &large_alg { norm_extra *= BigUint::from(p).pow(e / 2); }
        let Some(x) = algebraic_sqrt(&pairs, f, m, n, &alg_sum, &alg, &norm_extra) else { continue };
        if (f.len() - 1) % 2 == 0 {
            // the even-degree root carries f'(α); match it with f'(m) on the rational side
            let fpm = { let fp = poly_derivative(f); let mut v = BigInt::zero(); for c in fp.iter().rev() { v = (v * BigInt::from(m.clone()) + c).mod_floor(&BigInt::from(n.clone())); } v.to_biguint().unwrap() };
            y = y * fpm % n;
        }
        let diff = if x >= y { &x - &y } else { &y - &x };
        for g in [diff.gcd(target), (&x + &y).gcd(target)] { if g > BigUint::one() && g < *target { return Some(g); } }
    }
    None
}

/// A factor of 2^k − 1 by the special sieve with the polynomial from the form, degree 3.
/// Small factors (to a million) come off by trial division first; the sieve answers for the
/// cofactor, with the polynomial still the one for 2^k − 1.
pub fn factor_mersenne(k: u32) -> Option<BigUint> { factor_mersenne_degree(k, if k < 160 { 3 } else { 4 }) }

/// `factor_mersenne` with the degree chosen: 3 keeps the norms small, 4 keeps m small, which
/// matters once m passes 2^50. Odd degrees take Couveignes' square root, even ones Hensel's.
pub fn factor_mersenne_degree(k: u32, degree: u32) -> Option<BigUint> {
    let n = (BigUint::one() << k as usize) - BigUint::one();
    let mut cofactor = n.clone();
    let mut small = None;
    for p in primes_upto(1_000_000) { while (&cofactor % p).is_zero() { cofactor /= p; small.get_or_insert(p); } }
    if cofactor.is_one() || crate::big::is_prime_big(&cofactor) { return small.map(BigUint::from); }
    let (f, m) = polynomial(k, degree);
    factor(&n, &cofactor, &f, &m, Bounds::for_bits(k))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn the_polynomial_the_norm_and_the_field() {
        let (f, m) = polynomial(67, 3);
        assert_eq!(f[0], BigInt::from(-4));
        let n = (BigUint::one() << 67usize) - BigUint::one();
        let mut v = BigInt::zero(); for c in f.iter().rev() { v = (v * BigInt::from(m.clone()) + c).mod_floor(&BigInt::from(n.clone())); }
        assert!(v.is_zero(), "f(m) ≡ 0 mod n");
        assert_eq!(norm(&f, 3, 1), BigInt::from(31));
        let fq = f_mod(&f, 1048583);
        assert_eq!(mulmod(&[0, 0, 1], &[0, 1, 0], &fq, 1048583), vec![4, 0, 0], "x³ = 4");
        for q in [1048583u64, 1048589, 1048601, 1048609] { assert_eq!(irreducible_mod(&f_mod(&f, q), q), roots_mod(&f_mod(&f, q), q).is_empty(), "{q}"); }
    }
    #[test]
    fn the_field_square_root_and_the_lift() {
        let (f, m) = polynomial(67, 4);
        assert_eq!(f[0], BigInt::from(-2));
        let q = (1u64 << 30) + 3;
        let fq = f_mod(&f, q);
        let mut u = vec![0u64; 4]; u[0] = 5; u[1] = 7;
        let u2 = mulmod(&u, &u, &fq, q);
        if irreducible_mod(&fq, q) {
            let r = field_sqrt(&u2, &fq, q).unwrap();
            assert!(r == u || r == u.iter().map(|&c| (q - c) % q).collect::<Vec<_>>());
        }
        // a planted square: (3 + 2α)² (1 + α)², with f'(α)²: hensel_sqrt(m) ² ≡ f'(m)² ∏ (mod n)
        let n = (BigUint::one() << 67usize) - BigUint::one();
        let pairs = [(3i64, 2i64), (3, 2), (1, 1), (1, 1)];
        let x = hensel_sqrt(&pairs, &f, &m, &n).unwrap();
        let nb = BigInt::from(n.clone()); let mb = BigInt::from(m.clone());
        let fpm = { let fp = poly_derivative(&f); let mut v = BigInt::zero(); for c in fp.iter().rev() { v = (v * &mb + c).mod_floor(&nb); } v };
        let mut want = fpm.clone() * &fpm;
        for &(a, b) in &pairs { want = (want * (BigInt::from(a) + BigInt::from(b) * &mb)).mod_floor(&nb); }
        assert_eq!(BigInt::from(x.clone()) * BigInt::from(x) % &nb, want.mod_floor(&nb));
    }
    #[test]
    fn degree_four_factors_too() {
        let f = factor_mersenne_degree(101, 4).unwrap();
        let n = (BigUint::one() << 101usize) - BigUint::one();
        assert!((&n % &f).is_zero() && f > BigUint::one() && f < n);
    }
    #[test]
    fn cole_s_factor_and_two_more() {
        let cole = factor_mersenne(67).unwrap();
        assert!(cole == BigUint::from(193707721u64) || cole == BigUint::from(761838257287u64), "Cole's two, either way round");
        assert_eq!(factor_mersenne(11), Some(BigUint::from(23u32)), "a small one is trial division");
        assert_eq!(factor_mersenne(31), None, "a Mersenne prime has no factor");
        for k in [71u32, 79] {
            let n = (BigUint::one() << k as usize) - BigUint::one();
            let f = factor_mersenne(k).unwrap();
            assert!((&n % &f).is_zero() && f > BigUint::one() && f < n);
        }
    }
}
