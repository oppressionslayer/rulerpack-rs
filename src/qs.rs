//! **A quadratic sieve on the window** (`big` feature): the smooth sieve, twice per prime,
//! over x² − n. Find x with x² − n smooth over a factor base, combine enough of them so every
//! exponent is even, and X² ≡ Y² (mod n) gives gcd(X − Y, n), a factor more often than not.
//!
//! Each prime has two strides, from the two roots of x² ≡ n (mod p); the log sieve adds ln p
//! along them and the few positions that clear the threshold are divided exactly. Pure Rust
//! with num-bigint. Read it next to `rulerpack.qs` in the Python package: same shape.

use num_bigint::BigUint;
use num_integer::Integer;
use num_traits::{One, ToPrimitive, Zero};

use crate::smooth::{log_sieve, primes_upto};

/// Euler's criterion for an odd prime: 1, −1 or 0.
pub fn legendre(a: u64, p: u64) -> i8 {
    let s = modpow(a % p, (p - 1) / 2, p);
    if s == p - 1 { -1 } else { s as i8 }
}
fn modpow(mut b: u64, mut e: u64, m: u64) -> u64 {
    let mut r = 1u128; let mut bb = (b % m) as u128; let mm = m as u128;
    while e > 0 { if e & 1 == 1 { r = r * bb % mm; } bb = bb * bb % mm; e >>= 1; }
    b = r as u64; b
}

/// A root of x² ≡ n (mod p), odd prime p with `legendre(n, p) == 1` (Tonelli–Shanks).
pub fn sqrt_mod(n: u64, p: u64) -> u64 {
    let n = n % p;
    if p % 4 == 3 { return modpow(n, (p + 1) / 4, p); }
    let (mut q, mut s) = (p - 1, 0u32);
    while q % 2 == 0 { q /= 2; s += 1; }
    let mut z = 2; while legendre(z, p) != -1 { z += 1; }
    let (mut m, mut c, mut t, mut r) = (s, modpow(z, q, p), modpow(n, q, p), modpow(n, q.div_ceil(2), p));
    while t != 1 {
        let (mut i, mut t2) = (0u32, t);
        while t2 != 1 { t2 = ((t2 as u128 * t2 as u128) % p as u128) as u64; i += 1; }
        let b = modpow(c, 1u64 << (m - i - 1), p);
        let b2 = ((b as u128 * b as u128) % p as u128) as u64;
        m = i; c = b2; t = ((t as u128 * b2 as u128) % p as u128) as u64; r = ((r as u128 * b as u128) % p as u128) as u64;
    }
    r
}

/// 2 and the odd primes up to `bound` modulo which n is a square, each with its roots.
pub fn factor_base(n: &BigUint, bound: u64) -> Vec<(u64, Vec<u64>)> {
    let mut base = vec![(2u64, vec![(n % 2u64).to_u64().unwrap()])];
    for p in primes_upto(bound).into_iter().skip(1) {
        let np = (n % p).to_u64().unwrap();
        if legendre(np, p) == 1 {
            let r = sqrt_mod(np, p);
            let mut roots = vec![r, p - r]; roots.sort(); roots.dedup();
            base.push((p, roots));
        }
    }
    base
}

/// Exact: the exponents of `v` over the base and the cofactor left, 1 when smooth.
pub fn divide_out(mut v: BigUint, base: &[(u64, Vec<u64>)]) -> (Vec<u32>, BigUint) {
    let mut exps = vec![0u32; base.len()];
    for (i, (p, _)) in base.iter().enumerate() {
        let pb = BigUint::from(*p);
        while (&v % p).is_zero() { v /= &pb; exps[i] += 1; }
    }
    (exps, v)
}

/// Bound and segment by the size of n, as measured in the Python package.
pub fn parameters(digits: usize) -> (u64, usize) {
    match digits { 0..=20 => (600, 1 << 14), 21..=26 => (1500, 1 << 16), 27..=32 => (4000, 1 << 16), 33..=40 => (8000, 1 << 17), 41..=50 => (20000, 1 << 18), _ => (50000, 1 << 19) }
}

/// Sieve x = m + i upward in segments until `need` relations (x, exponents) are in hand.
pub fn relations(n: &BigUint, base: &[(u64, Vec<u64>)], m: &BigUint, segment: usize, need: usize, slack: f32) -> (Vec<(BigUint, Vec<u32>)>, usize) {
    let mut found = Vec::new();
    let pmax = base.last().unwrap().0 as f32;
    let mut seg = 0usize;
    while found.len() < need {
        let lo = m + BigUint::from(seg * segment);
        let strides: Vec<(u64, Vec<usize>)> = base.iter().map(|(p, roots)| {
            let lp = (&lo % *p).to_u64().unwrap();
            let mut s: Vec<usize> = roots.iter().map(|&r| ((r + p - lp) % p) as usize).collect();
            s.sort(); s.dedup(); (*p, s)
        }).collect();
        // ln(x² − n) with x = lo + i, as D + 2·lo·i + i², three positive terms, no cancellation
        let d = big_f64(&(&lo * &lo - n));
        let lo_f = big_f64(&lo);
        let cut = slack * pmax.ln();
        let hits = log_sieve(segment, &strides, |i| { let i = i as f64; ((d + 2.0 * lo_f * i + i * i).ln() as f32) - cut });
        for i in hits {
            let x = &lo + BigUint::from(i);
            let (exps, rest) = divide_out(&x * &x - n, base);
            if rest.is_one() { found.push((x, exps)); if found.len() >= need { break; } }
        }
        seg += 1;
    }
    (found, seg * segment)
}

fn big_f64(b: &BigUint) -> f64 { b.to_f64().unwrap_or(f64::MAX) }

/// Subsets of `rows` whose exponent sum is even everywhere: GF(2) elimination with bit rows.
pub fn dependencies(rows: &[Vec<u32>]) -> Vec<Vec<usize>> {
    let cols = rows.first().map_or(0, |r| r.len());
    let odd: Vec<Vec<u32>> = rows.iter().map(|r| r.iter().enumerate().filter(|(_, &e)| e & 1 == 1).map(|(j, _)| j as u32).collect()).collect();
    dependencies_sparse(&odd, cols)
}

/// `dependencies` from rows given as the columns with an odd exponent, `cols` columns in all.
pub fn dependencies_sparse(rows: &[Vec<u32>], cols: usize) -> Vec<Vec<usize>> {
    let words = cols.div_ceil(64).max(1);
    let rw = rows.len().div_ceil(64).max(1);
    let mut pivots: Vec<Option<(Vec<u64>, Vec<u64>)>> = vec![None; cols.max(1)];
    let mut out = Vec::new();
    for (i, r) in rows.iter().enumerate() {
        let mut v = vec![0u64; words];
        for &j in r { let j = j as usize; v[j / 64] |= 1 << (j % 64); }
        let mut who = vec![0u64; rw]; who[i / 64] |= 1 << (i % 64);
        // the highest set bit, a word at a time
        let top_bit = |v: &[u64]| v.iter().enumerate().rev().find(|(_, &w)| w != 0).map(|(k, &w)| k * 64 + 63 - w.leading_zeros() as usize);
        while let Some(top) = top_bit(&v) {
            match &pivots[top] {
                Some((pv, pw)) => { for k in 0..words { v[k] ^= pv[k]; } for k in 0..rw { who[k] ^= pw[k]; } }
                None => { pivots[top] = Some((v.clone(), who.clone())); break; }
            }
        }
        if v.iter().all(|&w| w == 0) { out.push((0..rows.len()).filter(|&k| who[k / 64] >> (k % 64) & 1 == 1).collect()); }
    }
    out
}

/// A nontrivial factor of the odd composite `n`, or None (prime, or every dependency trivial).
pub fn factor(n: &BigUint) -> Option<BigUint> { factor_with(n, None, None, 1.5, 12) }

/// `factor` with the knobs: factor-base bound, segment width, log slack, extra relations.
pub fn factor_with(n: &BigUint, bound: Option<u64>, segment: Option<usize>, slack: f32, extra: usize) -> Option<BigUint> {
    if *n < BigUint::from(4u32) { return None; }
    for p in primes_upto(1000) { if (n % p).is_zero() { return if *n == BigUint::from(p) { None } else { Some(BigUint::from(p)) }; } }
    let r = n.sqrt();
    if &r * &r == *n { return Some(r); }
    for k in 3..=n.bits() as u32 { let c = n.nth_root(k); if c.pow(k) == *n { return Some(c); } }
    if crate::big::is_prime_big(n) { return None; }
    let digits = n.to_string().len();
    let (b, s) = parameters(digits);
    let base = factor_base(n, bound.unwrap_or(b));
    let need = base.len() + extra;
    let m = &r + BigUint::one();
    let (rels, _) = relations(n, &base, &m, segment.unwrap_or(s), need, slack);
    if rels.len() <= base.len() { return None; }
    let rows: Vec<Vec<u32>> = rels.iter().map(|(_, e)| e.clone()).collect();
    for dep in dependencies(&rows) {
        let mut x = BigUint::one();
        let mut sum = vec![0u32; base.len()];
        for &i in &dep { x = &x * &rels[i].0 % n; for (j, e) in rels[i].1.iter().enumerate() { sum[j] += e; } }
        let mut y = BigUint::one();
        for (j, (p, _)) in base.iter().enumerate() { y = &y * BigUint::from(*p).modpow(&BigUint::from(sum[j] / 2), n) % n; }
        let diff = if x >= y { &x - &y } else { &y - &x };
        for g in [diff.gcd(n), (&x + &y).gcd(n)] { if g > BigUint::one() && g < *n { return Some(g); } }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn roots_and_base() {
        for p in [3u64, 5, 7, 13, 17, 41, 97, 101, 65537] { for a in 1..p.min(50) { if legendre(a, p) == 1 { let r = sqrt_mod(a, p); assert_eq!(r * r % p, a, "{a} mod {p}"); } } }
        let n = BigUint::from(8051u32);
        assert!(factor_base(&n, 50).len() > 3);
    }
    #[test]
    fn factors_small_and_medium() {
        assert_eq!(factor(&BigUint::from(8051u32)), Some(BigUint::from(83u32)));
        let (p, q) = (BigUint::from(1_000_000_007u64), BigUint::from(998_244_353u64));
        let f = factor(&(&p * &q)).unwrap(); assert!(f == p || f == q);
        // factors of unequal size, so the first x² − n is not a square (the equal-size pair was
        // Fermat-easy and proved nothing)
        let (p, q) = (BigUint::from(206158430209u64), BigUint::from(21990232567889u64));
        let f = factor(&(&p * &q)).unwrap(); assert!(f == p || f == q, "25 digits");
        assert_eq!(factor(&BigUint::from(1_000_000_007u64)), None, "a prime has no factor");
        assert_eq!(factor(&BigUint::from(3u32).pow(7)), Some(BigUint::from(3u32)), "a perfect power");
    }
    #[test]
    fn dependencies_are_even() {
        let rows = vec![vec![1, 0, 1], vec![0, 1, 1], vec![1, 1, 0], vec![1, 1, 2]];
        let deps = dependencies(&rows);
        assert!(!deps.is_empty());
        for d in deps { for j in 0..3 { assert_eq!(d.iter().map(|&i| rows[i][j]).sum::<u32>() % 2, 0); } }
    }
}
