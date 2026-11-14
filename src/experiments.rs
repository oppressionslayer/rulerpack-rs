//! **The drawer**: things that were fun in a session and might be more. Nothing else in the crate
//! uses it. Each function reproduces an IPython moment exactly and says what it saw, the twin of
//! `rulerpack.experiments` in the Python package.
//!
//! Experiment 1, the step-10 index (In [475]-[478], [481]):
//!
//! ```
//! use rulerpack::experiments as rx;
//! let primes = rulerpack::PiTable::new(100_000).primes_upto(100_000);
//! let idx = rx::stride_index(primes.len(), 2, 10);        // 2, 12, 22, ...
//! assert_eq!(&rx::stripe(&primes, 2, 10)[..6], &[2, 2, 0, 2, 0, 2]);
//! assert_eq!(&rx::last_before_blocks(&primes, 128)[..4], &[127, 251, 383, 509]);
//! assert_eq!(rx::recover(&primes, 32, 2, 10)[5], primes[5] % 16);
//! let _ = idx;
//! ```

/// `start, start+step, ...` below `len`: every tenth index from the third prime, by default.
pub fn stride_index(len: usize, start: usize, step: usize) -> Vec<u64> {
    (start..len).step_by(step.max(1)).map(|i| i as u64).collect()
}

/// The primes the stride index picks.
pub fn picked(primes: &[u64], start: usize, step: usize) -> Vec<u64> {
    stride_index(primes.len(), start, step).iter().map(|&i| primes[i as usize]).collect()
}

/// `idx ^ (2p)` over the first `len(idx)` primes: the array of In [477]. XOR cancels the top bits
/// the two sides share, so the values sit small and leap only when one side crosses a power of
/// two or a multiple of 256.
pub fn xor_index(primes: &[u64], start: usize, step: usize) -> Vec<u64> {
    stride_index(primes.len(), start, step).iter().zip(primes).map(|(&i, &p)| i ^ (2 * p)).collect()
}

/// `(idx ^ 2p) % 4`: 2, then 0, 2, 0, 2, ... for ever. Bit 1 of idx alternates because the step is
/// 2 mod 4; `2p` is 2 mod 4 for every odd p; the one exception is the first entry, where p = 2.
pub fn stripe(primes: &[u64], start: usize, step: usize) -> Vec<u64> {
    xor_index(primes, start, step).into_iter().map(|z| z % 4).collect()
}

/// Indices of the entries where `idx ^ 2p` has just leapt by more than `more_than`. Two kinds:
/// `2p` crossing a power of two (1024, 2048, 4096, 8192 at 97, 172, 309, 564 on the first 1000
/// primes) and idx crossing a multiple of 256, every 25.6 entries.
pub fn jumps(primes: &[u64], start: usize, step: usize, more_than: u64) -> Vec<usize> {
    let z = xor_index(primes, start, step);
    (1..z.len()).filter(|&i| z[i].abs_diff(z[i - 1]) > more_than).collect()
}

/// The indices where `2p` steps into the next `block`: with 256 that is p crossing a multiple of
/// 128. The spacing is the count of primes in each block, thinning as the prime number theorem says.
pub fn block_crossings(primes: &[u64], block: u64) -> Vec<usize> {
    (0..primes.len().saturating_sub(1)).filter(|&i| (2 * primes[i + 1]) / block != (2 * primes[i]) / block).collect()
}

/// The last prime under each multiple of `width`: 127, 251, 383, 509, 631, ... The largest primes
/// under powers of two show up on the dot: 127, 1021, 2039, 4093, 8191.
pub fn last_before_blocks(primes: &[u64], width: u64) -> Vec<u64> {
    block_crossings(primes, 2 * width).into_iter().map(|i| primes[i]).collect()
}

/// How far each `last_before_blocks` prime sits under its boundary. Always odd.
pub fn gaps_to_blocks(primes: &[u64], width: u64) -> Vec<u64> {
    last_before_blocks(primes, width).iter().enumerate().map(|(k, &p)| width * (k as u64 + 1) - p).collect()
}

/// `z % m` for each modulus, z = idx ^ 2p, over the first `n` entries: one more bit of the primes
/// per row, the stripe under all of them. `(a ^ b) % m == (a % m) ^ (b % m)` for a power of two.
pub fn layers(primes: &[u64], moduli: &[u64], start: usize, step: usize, n: usize) -> Vec<(u64, Vec<u64>)> {
    let z = xor_index(primes, start, step);
    moduli.iter().map(|&m| (m, z.iter().take(n).map(|&v| v % m).collect())).collect()
}

/// `(z >> bit) & 1` for each bit: the layers on their own.
pub fn bit_layers(primes: &[u64], bits: &[u32], start: usize, step: usize, n: usize) -> Vec<(u32, Vec<u64>)> {
    let z = xor_index(primes, start, step);
    bits.iter().map(|&b| (b, z.iter().take(n).map(|&v| (v >> b) & 1).collect())).collect()
}

/// Peel the index back off a layer: `((z % m) ^ (idx % m)) >> 1 == p % (m / 2)`. XOR is its own
/// inverse, so no layer hides anything; the same fact as `unslide512`.
pub fn recover(primes: &[u64], m: u64, start: usize, step: usize) -> Vec<u64> {
    stride_index(primes.len(), start, step).iter().zip(primes).map(|(&i, &p)| (((i ^ (2 * p)) % m) ^ (i % m)) >> 1).collect()
}

/// Rows of `idx % m`, `2p % m`, `z % m` over the first `n` entries: where each repeated value in
/// a layer comes from.
pub fn side_by_side(primes: &[u64], m: u64, start: usize, step: usize, n: usize) -> Vec<[u64; 3]> {
    stride_index(primes.len(), start, step).iter().zip(primes).take(n).map(|(&i, &p)| [i % m, (2 * p) % m, (i ^ (2 * p)) % m]).collect()
}

/// The stripe over odd composites instead of primes: 0, 2, 0, 2, ... with no exception, since
/// there is no even member. The stripe is a fact about bit 1 of an index stepping by 10 and of
/// twice an odd number, which is the control that says it is not about primes.
pub fn stripe_control(n: usize, start: usize, step: usize) -> Vec<u64> {
    let comps: Vec<u64> = (9..).step_by(2).filter(|o| [3u64, 5, 7].iter().any(|q| o % q == 0)).take(n).collect();
    stride_index(comps.len(), start, step).iter().zip(&comps).map(|(&i, &c)| (i ^ (2 * c)) % 4).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn experiment_one_the_step_ten_index() {
        let primes = crate::PiTable::new(120_000).primes_upto(120_000);
        let primes = &primes[..10_000];
        let idx = stride_index(primes.len(), 2, 10);
        assert_eq!(&idx[..4], &[2, 12, 22, 32]);
        assert!(idx.iter().take(100).enumerate().all(|(k, &i)| i % 4 == if k % 2 == 0 { 2 } else { 0 }));
        let x = xor_index(primes, 2, 10);
        assert_eq!(&x[..8], &[6, 10, 28, 46, 60, 46, 28, 110], "In [477]");
        let s = stripe(primes, 2, 10);
        assert_eq!(s[0], 2);
        assert!(s[1..200].iter().enumerate().all(|(k, &v)| v == if k % 2 == 0 { 2 } else { 0 }), "In [478]");
        let j = jumps(primes, 2, 10, 300);
        for i in [97usize, 172, 309, 564] {
            assert!(j.contains(&i), "2p crosses a power of two at {i}");
        }
        assert!(2 * primes[97] >= 1024 && 1024 > 2 * primes[96]);
        let last = last_before_blocks(primes, 128);
        assert_eq!(&last[..4], &[127, 251, 383, 509]);
        for p in [127u64, 1021, 2039, 4093, 8191] {
            assert!(last.contains(&p));
        }
        assert!(gaps_to_blocks(primes, 128).iter().all(|g| g % 2 == 1));
        assert_eq!(picked(primes, 2, 10)[..3], [5, 41, 83]);
    }

    #[test]
    fn experiment_one_continued_the_layers_peel_back_off() {
        let primes = crate::PiTable::new(120_000).primes_upto(120_000);
        let primes = &primes[..10_000];
        let l = layers(primes, &[4, 8, 16, 32], 2, 10, 40);
        assert_eq!(&l[0].1[..12], &[2, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2]);
        assert_eq!(&l[1].1[..12], &[6, 2, 4, 6, 4, 6, 4, 6, 4, 6, 0, 2], "In [481]");
        assert_eq!(&l[3].1[..9], &[6, 10, 28, 14, 28, 14, 28, 14, 28]);
        for k in 0..3 {
            let (lo, hi) = (&l[k], &l[k + 1]);
            assert!(hi.1.iter().zip(&lo.1).all(|(&h, &lw)| h % lo.0 == lw), "every layer keeps the one below");
        }
        let b = bit_layers(primes, &[1, 2, 3, 4], 2, 10, 40);
        assert!(b[0].1.iter().zip(&l[0].1).all(|(&bit, &v)| bit == v / 2));
        for m in [4u64, 8, 16, 32, 256] {
            let r = recover(primes, m, 2, 10);
            assert!(r.iter().zip(primes).all(|(&v, &p)| v == p % (m / 2)), "the % {m} row carries p mod {}", m / 2);
        }
        let cols = side_by_side(primes, 32, 2, 10, 40);
        assert!(cols.iter().all(|c| c[0] ^ c[1] == c[2]));
        assert_eq!(cols[2], [22, 10, 28], "idx 22, 2*5 = 10, and 22 ^ 10 = 28");
        let c = stripe_control(1000, 2, 10);
        assert!(c.iter().enumerate().all(|(k, &v)| v == if k % 2 == 0 { 0 } else { 2 }), "0, 2, 0, 2 with no exception");
    }
}
