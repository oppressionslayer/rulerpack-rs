//! Big anchors on GMP (the `rug` crate) with rayon threads. Same window sieve
//! and the same two searches as `big`, with GMP's modular power instead of
//! num-bigint's. Enabled by the `gmp` feature (needs libgmp-dev).

use rayon::prelude::*;
use rug::{integer::IsPrime, Integer};

pub use crate::bands::{seed_depth, window_width, Certify};
use crate::PiTable;

fn seeds_for(bits: usize) -> Vec<u64> {
    match seed_depth(bits) {
        None => crate::SEEDS.to_vec(),
        Some(depth) => PiTable::new(depth).primes_upto(depth).into_iter().skip(1).collect(),
    }
}

/// Survivor offsets d in [0, width): anchor + d odd and free of every seed.
pub fn window_gmp(anchor: &Integer, width: usize, seeds: &[u64]) -> Vec<usize> {
    let mut alive = vec![true; width];
    let mut d = if anchor.is_odd() { 1 } else { 0 };
    while d < width { alive[d] = false; d += 2; }
    for &q in seeds {
        let r = anchor.mod_u(q as u32) as u64;
        let first = ((q - r) % q) as usize;
        let mut d = first;
        while d < width { alive[d] = false; d += q as usize; }
    }
    alive.iter().enumerate().filter(|(_, &a)| a).map(|(d, _)| d).collect()
}

/// Strong probable-prime test to the given bases (GMP powmod).
pub fn strong_prp(n: &Integer, bases: &[u32]) -> bool {
    if *n < 4 { return *n == 2 || *n == 3; }
    if n.is_even() { return false; }
    let nm1 = Integer::from(n - 1u32);
    let s = nm1.find_one(0).unwrap();
    let d = Integer::from(&nm1 >> s);
    'base: for &a in bases {
        let mut x = Integer::from(a).pow_mod(&d, n).unwrap();
        if x == 1 || x == nm1 { continue; }
        for _ in 1..s {
            x = x.square().modulo(n);
            if x == nm1 { continue 'base; }
        }
        return false;
    }
    true
}

fn passes(c: &Integer, certify: Certify) -> bool {
    match certify {
        // GMP's own test (BPSW plus Miller-Rabin rounds), after the cheap base-2 pass
        Certify::Full => strong_prp(c, &[2]) && c.is_probably_prime(25) != IsPrime::No,
        Certify::Prp => strong_prp(c, &[2, 3]),
    }
}

/// Test candidates in rounds of 2 * threads, lowest first: every round is a
/// parallel map over the next slice, so all cores work the low end of the
/// window together. (Rayon's work stealing alone would spread the threads over
/// the whole window and leave the low candidates to one thread in sequence.)
fn rounds(cands: &[Integer], certify: Certify, lowest: bool) -> Option<Integer> {
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

/// Next prime above anchor: every core, smallest pass wins.
pub fn next_prime_after_gmp(anchor: &Integer, certify: Certify) -> Integer {
    let bits = anchor.significant_bits() as usize;
    let seeds = seeds_for(bits);
    let width = window_width(bits);
    let mut base = Integer::from(anchor + 1u32);
    loop {
        let cands: Vec<Integer> = window_gmp(&base, width, &seeds).into_iter().map(|d| Integer::from(&base + d as u64)).collect();
        if let Some(p) = rounds(&cands, certify, true) { return p; }
        base += width as u64;
    }
}

/// A prime above anchor: every core, first pass from any core wins.
pub fn any_prime_after_gmp(anchor: &Integer, certify: Certify) -> Integer {
    let bits = anchor.significant_bits() as usize;
    let seeds = seeds_for(bits);
    let width = window_width(bits);
    let mut base = Integer::from(anchor + 1u32);
    loop {
        let cands: Vec<Integer> = window_gmp(&base, width, &seeds).into_iter().map(|d| Integer::from(&base + d as u64)).collect();
        if let Some(p) = rounds(&cands, certify, false) { return p; }
        base += width as u64;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn known_primes_above_powers_of_two() {
        let a = Integer::from(Integer::u_pow_u(2, 256));
        assert_eq!(Integer::from(&next_prime_after_gmp(&a, Certify::Full) - &a), 297);
        assert_eq!(Integer::from(&next_prime_after_gmp(&a, Certify::Prp) - &a), 297);
        let b = Integer::from(Integer::u_pow_u(2, 1024));
        assert_eq!(Integer::from(&next_prime_after_gmp(&b, Certify::Full) - &b), 643);
        let p = any_prime_after_gmp(&b, Certify::Full);
        assert!(p > b && p.is_probably_prime(30) != IsPrime::No);
    }
}
