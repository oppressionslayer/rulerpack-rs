//! **The smooth window**: the anchored sieve with the sign flipped.
//!
//! `window` clears a position the moment a seed divides it and keeps the primes. Here each
//! prime walks the same strides, first `d` with `p | anchor + d` and every `p` after, and is
//! divided out instead, as many times as it goes. A position divided down to 1 is B-smooth:
//! every prime factor is at most B, and the exponents came out of the same pass. Exact, no
//! logarithm thresholds. `log_sieve` is the wide-window companion: `log p` added along the
//! strides into a float array, with exact division only at positions that clear a threshold.
//!
//! Lars Rocha's method, 2026-11-09, ported from the Python package.

/// Primes at most `bound`, by a plain sieve.
pub fn primes_upto(bound: u64) -> Vec<u64> {
    if bound < 2 { return vec![]; }
    let n = bound as usize;
    let mut s = vec![true; n + 1];
    s[0] = false; s[1] = false;
    let mut p = 2;
    while p * p <= n {
        if s[p] { let mut m = p * p; while m <= n { s[m] = false; m += p; } }
        p += 1;
    }
    (0..=n).filter(|&i| s[i]).map(|i| i as u64).collect()
}

/// One smooth position: its offset and the `(prime, exponent)` pairs that multiply to it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Smooth { pub offset: usize, pub factors: Vec<(u64, u32)> }

/// Offsets `d` in `[0, width)` with `anchor + d` `bound`-smooth, ascending, each with its
/// factorisation. `anchor + d` must fit a `u64`; `smooth_window_big` takes any size.
pub fn smooth_window(anchor: u64, width: usize, bound: u64) -> Vec<Smooth> {
    smooth_window_seeds(anchor, width, &primes_upto(bound))
}

/// `smooth_window` over a custom seed list (any integers above 1): "smooth" then means
/// "a product of those".
pub fn smooth_window_seeds(anchor: u64, width: usize, seeds: &[u64]) -> Vec<Smooth> {
    assert!(seeds.iter().all(|&p| p > 1), "seeds must be greater than 1");
    let mut rest: Vec<u64> = (0..width as u64).map(|d| anchor + d).collect();
    let mut exps: Vec<Vec<(u64, u32)>> = vec![Vec::new(); width];
    for &p in seeds {
        let first = ((p - anchor % p) % p) as usize;
        let mut d = first;
        while d < width {
            let mut e = 0;
            while rest[d] != 0 && rest[d] % p == 0 { rest[d] /= p; e += 1; }
            if e > 0 { exps[d].push((p, e)); }
            d += p as usize;
        }
    }
    (0..width).filter(|&d| rest[d] == 1 && anchor + d as u64 > 0)
        .map(|d| Smooth { offset: d, factors: std::mem::take(&mut exps[d]) }).collect()
}

/// Trial division: is `n` `bound`-smooth? The oracle the tests hold the sieve against.
pub fn is_smooth(mut n: u64, bound: u64) -> bool {
    if n == 0 { return false; }
    for p in primes_upto(bound) { while n % p == 0 { n /= p; } if n == 1 { return true; } }
    n == 1
}

/// **The log sieve.** `strides[k] = (p, starts)`: the offsets in `[0, p)` where `p` divides
/// the value. Adds `ln p` along every stride and returns the offsets whose total reaches
/// `threshold(d)`. Exact confirmation is the caller's.
pub fn log_sieve(width: usize, strides: &[(u64, Vec<usize>)], threshold: impl Fn(usize) -> f32) -> Vec<usize> {
    let acc = log_accumulate(width, strides);
    (0..width).filter(|&d| acc[d] >= threshold(d)).collect()
}

/// The log sieve's accumulator itself: `ln p` summed along every stride, one entry a position.
pub fn log_accumulate(width: usize, strides: &[(u64, Vec<usize>)]) -> Vec<f32> {
    let mut acc = vec![0f32; width];
    for (p, starts) in strides {
        let lp = (*p as f32).ln();
        for &s in starts { let mut d = s; while d < width { acc[d] += lp; d += *p as usize; } }
    }
    acc
}

#[cfg(test)]
mod tests {
    use super::*;
    fn brute(anchor: u64, width: usize, bound: u64) -> Vec<usize> {
        (0..width).filter(|&d| is_smooth(anchor + d as u64, bound)).collect()
    }
    #[test]
    fn matches_trial_division_and_the_factors_multiply_back() {
        for bound in [2, 3, 7, 13, 50] {
            let got = smooth_window(0, 2000, bound);
            assert_eq!(got.iter().map(|s| s.offset).collect::<Vec<_>>(), brute(0, 2000, bound), "bound {bound}");
            for s in &got {
                let back: u64 = s.factors.iter().map(|&(p, e)| p.pow(e)).product();
                assert_eq!(back, s.offset as u64);
                assert!(s.factors.iter().all(|&(p, _)| p <= bound));
            }
        }
        let five: Vec<usize> = smooth_window(0, 32, 5).iter().map(|s| s.offset).collect();
        assert_eq!(five, vec![1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 15, 16, 18, 20, 24, 25, 27, 30]);
        assert!(smooth_window(0, 1, 5).is_empty(), "zero is not smooth");
    }
    #[test]
    fn a_large_anchor_and_custom_seeds() {
        let anchor = 1_000_000_007u64 * 1_000_003;
        let got: Vec<usize> = smooth_window(anchor, 5000, 1000).iter().map(|s| s.offset).collect();
        assert_eq!(got, brute(anchor, 5000, 1000));
        let planted = 2u64.pow(20) * 3u64.pow(5) * 997u64.pow(2);
        let hit = smooth_window(planted - 10, 20, 1000).into_iter().find(|s| s.offset == 10).unwrap();
        assert_eq!(hit.factors, vec![(2, 20), (3, 5), (997, 2)]);
        let two_three: Vec<usize> = smooth_window_seeds(0, 100, &[2, 3]).iter().map(|s| s.offset).collect();
        assert_eq!(two_three, vec![1, 2, 3, 4, 6, 8, 9, 12, 16, 18, 24, 27, 32, 36, 48, 54, 64, 72, 81, 96]);
    }
    #[test]
    fn the_log_sieve_finds_what_the_exact_one_finds() {
        // seeds 2, 3, 5 over 0..256: the log total equals ln(value) exactly when smooth
        let strides: Vec<(u64, Vec<usize>)> = [2u64, 3, 5].iter().map(|&p| (p, vec![0usize])).collect();
        // one stride per prime only counts p once, so ask for "at least ln p of every prime
        // that divides" by a loose threshold and confirm exactly
        let loose = log_sieve(256, &strides, |_| 0.5);
        let exact: Vec<usize> = smooth_window(0, 256, 5).iter().map(|s| s.offset).filter(|&d| d > 1).collect();
        assert!(exact.iter().all(|d| loose.contains(d)));
    }
}
