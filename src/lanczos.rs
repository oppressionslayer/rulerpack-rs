//! **Block Lanczos over GF(2)** (Montgomery, 1995), for the dependencies of a large sparse
//! relation matrix where dense elimination stops fitting. Written from the paper; nothing
//! here is ported from any GPL siever.
//!
//! The matrix is `rows` relations over `cols` columns, each relation a sorted list of the
//! columns with an odd exponent. A dependency is a set of relations whose columns cancel: a
//! vector z with Bᵀz = 0. The iteration works with the symmetric A = B·Bᵀ applied to blocks of
//! 64 vectors at once (an n×64 block is a `Vec<u64>`, one word a row), solves A·X = A·Y for a
//! random Y so that X − Y lies in the kernel of A, and a final small elimination turns kernel
//! vectors of A into kernel vectors of Bᵀ.

type Block = Vec<u64>;        // n × 64
type Small = [u64; 64];       // 64 × 64, row i is a word

/// The sparse matrix as the black box A·X = B·(Bᵀ·X).
pub struct Sparse<'a> { pub rows: &'a [Vec<u32>], pub cols: usize }

impl Sparse<'_> {
    fn n(&self) -> usize { self.rows.len() }
    /// Bᵀ·X: cols × 64.
    fn bt(&self, x: &[u64]) -> Block {
        let mut out = vec![0u64; self.cols];
        for (r, row) in self.rows.iter().enumerate() { let w = x[r]; if w != 0 { for &c in row { out[c as usize] ^= w; } } }
        out
    }
    /// B·Y: rows × 64.
    fn b(&self, y: &[u64]) -> Block {
        self.rows.iter().map(|row| row.iter().fold(0u64, |acc, &c| acc ^ y[c as usize])).collect()
    }
    fn a(&self, x: &[u64]) -> Block { self.b(&self.bt(x)) }
}

/// Vᵀ·W as a 64×64 matrix: entry (i, j) = Σ_r V[r]_i · W[r]_j.
fn inner(v: &[u64], w: &[u64]) -> Small {
    let mut out = [0u64; 64];
    for (&a, &b) in v.iter().zip(w) {
        let mut bits = a;
        while bits != 0 { let i = bits.trailing_zeros() as usize; out[i] ^= b; bits &= bits - 1; }
    }
    out
}
/// V·M for an n×64 block and a 64×64 matrix: row r becomes XOR of M's rows picked by V[r].
fn mul_small(v: &[u64], m: &Small) -> Block {
    // eight 256-entry tables, one per byte of the input word
    let mut tables = vec![[0u64; 256]; 8];
    for (k, table) in tables.iter_mut().enumerate() {
        for byte in 1..256usize { let low = byte & (byte - 1); table[byte] = table[low] ^ m[k * 8 + (byte ^ low).trailing_zeros() as usize]; }
    }
    v.iter().map(|&w| { let mut acc = 0u64; for (k, table) in tables.iter().enumerate() { acc ^= table[((w >> (8 * k)) & 0xff) as usize]; } acc }).collect()
}
fn small_mul(a: &Small, b: &Small) -> Small {
    let mut out = [0u64; 64];
    for i in 0..64 { let mut bits = a[i]; while bits != 0 { let j = bits.trailing_zeros() as usize; out[i] ^= b[j]; bits &= bits - 1; } }
    out
}
fn small_add(a: &Small, b: &Small) -> Small { let mut o = *a; for i in 0..64 { o[i] ^= b[i]; } o }
fn mask_cols(a: &Small, mask: u64) -> Small { let mut o = *a; for w in o.iter_mut() { *w &= mask; } o }
fn identity() -> Small { let mut o = [0u64; 64]; for (i, w) in o.iter_mut().enumerate() { *w = 1 << i; } o }
fn is_zero(a: &Small) -> bool { a.iter().all(|&w| w == 0) }

/// Montgomery's choice of S_i and W_i⁻¹ from T = Vᵢᵀ·A·Vᵢ, preferring the columns left out
/// of the previous S. Returns (W_inv with zero rows and columns outside S, the mask of S).
fn choose(t: &Small, prev_mask: u64) -> (Small, u64) {
    // M = [T | I] as 64 rows of (left, right); columns ordered: not in prev first
    let mut left = *t;
    let mut right = identity();
    let order: Vec<usize> = (0..64).filter(|&c| prev_mask >> c & 1 == 0).chain((0..64).filter(|&c| prev_mask >> c & 1 == 1)).collect();
    let mut s_mask = 0u64;
    for (j, &cj) in order.iter().enumerate() {
        // a pivot on the left in column cj among the rows order[j..]
        if let Some(k) = (j..64).find(|&k| left[order[k]] >> cj & 1 == 1) {
            let ck = order[k];
            left.swap(cj, ck); right.swap(cj, ck);
            for r in 0..64 { if r != cj && left[r] >> cj & 1 == 1 { left[r] ^= left[cj]; right[r] ^= right[cj]; } }
            s_mask |= 1 << cj;
        } else if let Some(k) = (j..64).find(|&k| right[order[k]] >> cj & 1 == 1) {
            let ck = order[k];
            left.swap(cj, ck); right.swap(cj, ck);
            for r in 0..64 { if r != cj && right[r] >> cj & 1 == 1 { left[r] ^= left[cj]; right[r] ^= right[cj]; } }
            left[cj] = 0; right[cj] = 0;
        }
    }
    // W_inv is the right half on the rows and columns of S
    let mut winv = [0u64; 64];
    for r in 0..64 { if s_mask >> r & 1 == 1 { winv[r] = right[r] & s_mask; } }
    (winv, s_mask)
}

/// Up to 64 dependencies: sets of row indices whose columns cancel. `seed` picks the random
/// start; a second call with another seed gives other dependencies.
pub fn dependencies(rows: &[Vec<u32>], cols: usize, seed: u64) -> Vec<Vec<usize>> {
    let m = Sparse { rows, cols };
    let n = m.n();
    if n == 0 { return vec![]; }
    let mut rng = seed ^ 0x9E3779B97F4A7C15;
    let mut next = || { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; rng };
    let y: Block = (0..n).map(|_| next()).collect();
    let v0 = m.a(&y);
    let mut x: Block = vec![0u64; n];
    let (mut v_prev, mut v_prev2): (Block, Block) = (vec![0u64; n], vec![0u64; n]);
    let mut v = v0.clone();
    let (mut winv_prev, mut winv_prev2) = ([0u64; 64], [0u64; 64]);
    let (mut vav_prev, mut vaav_prev, mut s_prev) = ([0u64; 64], [0u64; 64], u64::MAX);
    let mut iterations = 0usize;
    loop {
        let av = m.a(&v);
        let vav = inner(&v, &av);
        let vaav = inner(&av, &av);
        if is_zero(&vav) { break; }
        let (winv, s) = choose(&vav, s_prev);
        // X += V · Winv · Vᵀ V0
        let vtv0 = inner(&v, &v0);
        let step = mul_small(&v, &small_mul(&winv, &vtv0));
        for (xi, si) in x.iter_mut().zip(&step) { *xi ^= si; }
        // V_{i+1} = A V S Sᵀ + V D + V_{i-1} E + V_{i-2} F
        let vaav_s = mask_cols(&vaav, s);
        let d = small_add(&identity(), &small_mul(&winv, &small_add(&vaav_s, &vav)));
        let e = small_mul(&winv_prev, &mask_cols(&vav, s));
        let inner_f = small_add(&identity(), &small_mul(&vav_prev, &winv_prev));
        let f = small_mul(&winv_prev2, &mask_cols(&small_mul(&inner_f, &small_add(&mask_cols(&vaav_prev, s_prev), &vav_prev)), s));
        let mut v_next: Block = av.iter().map(|&w| w & s).collect();
        for (vn, t) in v_next.iter_mut().zip(mul_small(&v, &d)) { *vn ^= t; }
        for (vn, t) in v_next.iter_mut().zip(mul_small(&v_prev, &e)) { *vn ^= t; }
        for (vn, t) in v_next.iter_mut().zip(mul_small(&v_prev2, &f)) { *vn ^= t; }
        v_prev2 = std::mem::replace(&mut v_prev, std::mem::replace(&mut v, v_next));
        winv_prev2 = std::mem::replace(&mut winv_prev, winv);
        vav_prev = vav; vaav_prev = vaav; s_prev = s;
        iterations += 1;
        if iterations > n / 32 + 100 { break; }   // did not converge: give back what there is
    }
    // kernel of A: columns of Z = [X ⊕ Y | V_last]; make kernel vectors of Bᵀ from them
    let z0: Block = x.iter().zip(&y).map(|(a, b)| a ^ b).collect();
    let z1 = v;
    let (u0, u1) = (m.bt(&z0), m.bt(&z1));
    // the 128 columns of [U0 | U1] as vectors over cols; find their linear dependencies
    let words = cols.div_ceil(64).max(1);
    let mut vectors: Vec<(Vec<u64>, u128)> = (0..128).map(|j| {
        let (src, bit) = if j < 64 { (&u0, j) } else { (&u1, j - 64) };
        let mut vcol = vec![0u64; words];
        for (c, &w) in src.iter().enumerate() { if w >> bit & 1 == 1 { vcol[c / 64] |= 1 << (c % 64); } }
        (vcol, 1u128 << j)
    }).collect();
    let mut combos: Vec<u128> = Vec::new();
    let mut pivots: Vec<(usize, Vec<u64>, u128)> = Vec::new();
    for (mut vcol, mut combo) in vectors.drain(..) {
        loop {
            match vcol.iter().enumerate().rev().find(|(_, &w)| w != 0).map(|(k, &w)| k * 64 + 63 - w.leading_zeros() as usize) {
                None => { combos.push(combo); break; }
                Some(top) => match pivots.iter().find(|(t, _, _)| *t == top) {
                    Some((_, pv, pc)) => { for (a, b) in vcol.iter_mut().zip(pv) { *a ^= b; } combo ^= pc; }
                    None => { pivots.push((top, vcol, combo)); break; }
                },
            }
        }
    }
    let mut out = Vec::new();
    for combo in combos {
        let lo = combo as u64; let hi = (combo >> 64) as u64;
        let zvec: Vec<bool> = (0..n).map(|r| ((z0[r] & lo).count_ones() + (z1[r] & hi).count_ones()) % 2 == 1).collect();
        let dep: Vec<usize> = (0..n).filter(|&r| zvec[r]).collect();
        if dep.is_empty() { continue; }
        // check Bᵀ z = 0 exactly before handing it on
        let mut col_sum = vec![0u8; cols];
        for &r in &dep { for &c in &rows[r] { col_sum[c as usize] ^= 1; } }
        if col_sum.iter().all(|&b| b == 0) { out.push(dep); }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    fn random_rows(n: usize, cols: usize, weight: usize, seed: u64) -> Vec<Vec<u32>> {
        let mut s = seed | 1;
        (0..n).map(|_| { let mut r: Vec<u32> = (0..weight).map(|_| { s ^= s << 13; s ^= s >> 7; s ^= s << 17; (s % cols as u64) as u32 }).collect(); r.sort_unstable(); r.dedup(); r }).collect()
    }
    #[test]
    fn the_small_pieces() {
        let id = identity();
        let t = [0x1234_5678_9abc_def0u64; 64];
        assert_eq!(small_mul(&id, &t), t);
        let v: Vec<u64> = vec![1, 2, 3, u64::MAX];
        assert_eq!(mul_small(&v, &id), v);
        let w = inner(&v, &v);
        assert_eq!(w[0], 1 ^ 3 ^ u64::MAX, "row 0 of VᵀV collects the rows with bit 0");
    }
    #[test]
    fn finds_dependencies_that_check() {
        for (n, cols, seed) in [(300usize, 250usize, 1u64), (2000, 1800, 2), (5000, 4700, 3)] {
            let rows = random_rows(n, cols, 12, seed);
            let deps = dependencies(&rows, cols, seed);
            assert!(!deps.is_empty(), "n={n}: at least one dependency");
            for d in &deps {
                let mut sum = vec![0u8; cols];
                for &r in d { for &c in &rows[r] { sum[c as usize] ^= 1; } }
                assert!(sum.iter().all(|&b| b == 0));
            }
            // the elimination agrees that dependencies exist here
            #[cfg(feature = "big")]
            assert!(!crate::qs::dependencies_sparse(&rows, cols).is_empty());
        }
    }
}
