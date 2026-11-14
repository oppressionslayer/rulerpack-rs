//! Big anchors: next prime above any-size integer. Pure Rust bignum (num-bigint)
//! with the window sieve and parallel certification (rayon). Enabled by the
//! `big` feature.

use num_bigint::BigUint;
use num_integer::Integer;
use num_traits::{One, ToPrimitive, Zero};
use rayon::prelude::*;

use crate::PiTable;

/// Primes up to `limit` as u64 (the seed list).
fn seeds_upto(limit: u64) -> Vec<u64> {
    let t = PiTable::new(limit);
    t.primes_upto(limit).into_iter().skip(1).collect()     // odd seeds only
}

/// Survivor offsets d in [0, width): anchor + d odd and free of every seed.
pub fn window_big(anchor: &BigUint, width: usize, seeds: &[u64]) -> Vec<usize> {
    let mut alive = vec![true; width];
    let odd = anchor.is_odd();
    let mut d = if odd { 1 } else { 0 };
    while d < width { alive[d] = false; d += 2; }
    for &q in seeds {
        let r = (anchor % q).to_u64().unwrap();
        let first = ((q - r) % q) as usize;
        let mut d = first;
        while d < width { alive[d] = false; d += q as usize; }
    }
    alive.iter().enumerate().filter(|(_, &a)| a).map(|(d, _)| d).collect()
}

/// Miller-Rabin with the given bases. Base 2 alone is the cheap first pass.
pub fn strong_prp(n: &BigUint, bases: &[u64]) -> bool {
    let one = BigUint::one();
    let two = BigUint::from(2u32);
    if *n < BigUint::from(4u32) { return *n == two || *n == BigUint::from(3u32); }
    if n.is_even() { return false; }
    let nm1 = n - &one;
    let s = nm1.trailing_zeros().unwrap();
    let d = &nm1 >> s;
    'base: for &a in bases {
        let mut x = BigUint::from(a).modpow(&d, n);
        if x == one || x == nm1 { continue; }
        for _ in 1..s {
            x = x.modpow(&two, n);
            if x == nm1 { continue 'base; }
        }
        return false;
    }
    true
}

/// gmpy2.is_prime's shape: base 2 first, then more bases only if it passed.
pub fn is_prime_big(n: &BigUint) -> bool {
    strong_prp(n, &[2]) && strong_prp(n, &[3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97])
}

/// Seed list for an anchor of this many bits (odd primes up to `seed_depth`).
pub fn seeds_for(bits: usize) -> Vec<u64> {
    match seed_depth(bits) {
        None => crate::SEEDS.to_vec(),
        Some(depth) => seeds_upto(depth),
    }
}

pub use crate::bands::{seed_depth, window_width, Certify};

fn passes(c: &BigUint, certify: Certify) -> bool {
    match certify {
        Certify::Full => is_prime_big(c),
        Certify::Prp => strong_prp(c, &[2, 3]),
    }
}

/// Rounds of 2 * threads, lowest candidates first, so every core works the low
/// end of the window together (see gmp.rs for why plain work stealing is slow).
fn rounds(cands: &[BigUint], certify: Certify, lowest: bool) -> Option<BigUint> {
    let per = rayon::current_num_threads() * 2;
    for chunk in cands.chunks(per) {
        let hit = if lowest {
            chunk.par_iter().with_max_len(1).find_first(|c| passes(c, certify))
        } else {
            chunk.par_iter().with_max_len(1).find_any(|c| passes(c, certify))
        };
        if let Some(p) = hit { return Some(p.clone()); }
    }
    None
}

/// Next prime above anchor: window sieve, then certify survivors on every core,
/// keeping the smallest that passes.
pub fn next_prime_after_big(anchor: &BigUint, parallel: bool) -> BigUint {
    next_prime_after_big_with(anchor, parallel, Certify::Full)
}

pub fn next_prime_after_big_with(anchor: &BigUint, parallel: bool, certify: Certify) -> BigUint {
    let bits = anchor.bits() as usize;
    let seeds = seeds_for(bits);
    let width = window_width(bits);
    let mut base = anchor + BigUint::one();
    loop {
        let offs = window_big(&base, width, &seeds);
        let cands: Vec<BigUint> = offs.iter().map(|&d| &base + BigUint::from(d)).collect();
        let hit = if parallel {
            rounds(&cands, certify, true)
        } else {
            cands.iter().find(|c| passes(c, certify)).cloned()
        };
        if let Some(p) = hit { return p; }
        base += BigUint::from(width);
    }
}

/// A prime above anchor: the first one any core certifies, not necessarily the
/// next (rayon's `find_any`). No stop rule to wait on, so every core's pass is
/// an answer. The fast path for "a random prime of this size".
pub fn any_prime_after_big(anchor: &BigUint, certify: Certify) -> BigUint {
    let bits = anchor.bits() as usize;
    let seeds = seeds_for(bits);
    let width = window_width(bits);
    let mut base = anchor + BigUint::one();
    loop {
        let offs = window_big(&base, width, &seeds);
        let cands: Vec<BigUint> = offs.iter().map(|&d| &base + BigUint::from(d)).collect();
        if let Some(p) = rounds(&cands, certify, false) { return p; }
        base += BigUint::from(width);
    }
}

/// **The wheel pass on a big anchor**: for each of `anchor + d`, `d < count`, the first of
/// `[2, seeds...]` that divides it, or `None` when it survives. `2^8191 - 1` with `&SEEDS` is the
/// Python example.
pub fn wheel_pass_big(anchor: &BigUint, count: u64, seeds: &[u64]) -> Vec<Option<u64>> {
    (0..count).map(|d| {
        let n = anchor + BigUint::from(d);
        std::iter::once(2u64).chain(seeds.iter().copied()).find(|&q| (&n % q).is_zero())
    }).collect()
}

#[cfg(test)]
mod tests {
    #[test]
    fn the_big_wheel_pass_matches_the_small_one() {
        use num_bigint::BigUint;
        let big = super::wheel_pass_big(&BigUint::from(8191u64), 100, &crate::SEEDS);
        assert_eq!(big, crate::wheel_pass(8191, 100, &crate::SEEDS));
        let m = (BigUint::from(1u8) << 8191u32) - BigUint::from(1u8);
        let w = super::wheel_pass_big(&m, 100, &crate::SEEDS);
        assert_eq!(w[0], None, "the Mersenne prime survives");
        assert_eq!(w[1], Some(2));
        assert_eq!(w.iter().filter(|q| q.is_none()).count(), 12);
    }
    use super::*;
    #[test]
    fn bands_and_width() {
        assert_eq!(seed_depth(512), None);
        assert_eq!(seed_depth(1024), Some(8191));
        assert_eq!(seed_depth(2048), Some(125 * 2048));
        assert_eq!(seed_depth(8192), Some(1_000_000));
        assert!(window_width(8192) > 39_000 && window_width(8192) < 41_000);
    }
    #[test]
    fn next_and_any_agree_with_known_primes() {
        // 2^256 + 297 is the next prime above 2^256; 2^1024 + 643 above 2^1024
        let a = BigUint::from(1u32) << 256u32;
        assert_eq!(&next_prime_after_big(&a, true) - &a, BigUint::from(297u32));
        assert_eq!(&next_prime_after_big(&a, false) - &a, BigUint::from(297u32));
        assert_eq!(&next_prime_after_big_with(&a, true, Certify::Prp) - &a, BigUint::from(297u32));
        let p = any_prime_after_big(&a, Certify::Full);
        assert!(p > a && is_prime_big(&p));
        let p = any_prime_after_big(&a, Certify::Prp);
        assert!(is_prime_big(&p));
        let b = BigUint::from(1u32) << 1024u32;
        assert_eq!(&next_prime_after_big(&b, true) - &b, BigUint::from(643u32));
    }
}
