//! **Square lanes**: an exploratory implementation on Lars Rocha's square family,
//! ported with the same oracle. Over x from ⌈√n⌉ up, Q(x) = x² − n. Fix the phase x ≡ φ (mod 4)
//! and a set of forced odd primes P that must divide Q: by the Chinese remainder theorem the
//! roots fall into lanes x = x₀ + 4P·k, along which C(k) = Q/P is a quadratic in k whose values
//! come from two additions each. Each lane is sieved exactly, the forced primes' exponents kept,
//! and every result is checked against trial division of the full Q.
//!
//! The question the experiment asks is whether sieving only the lanes finds smooth Q faster than
//! sieving every x. `examples/square_lanes.rs` measures it; the answer is in the README.

use crate::smooth::primes_upto;

/// One lane: roots `first + step·k` for k below `count`, with the forced primes.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Lane { pub first: u128, pub step: u128, pub count: usize, pub forced: Vec<u64> }

impl Lane {
    /// `(c, b, a)` with C(k) = Q(x₀ + step·k)/P = c + b·k + a·k², exactly.
    pub fn coefficients(&self, n: u128) -> (i128, i128, i128) {
        let p: u128 = self.forced.iter().map(|&q| q as u128).product();
        let q0 = self.first * self.first - n;
        (((q0 / p) as i128), ((2 * self.first * self.step / p) as i128), ((self.step * self.step / p) as i128))
    }
    /// The cofactors C(k), by two additions each.
    pub fn cofactors(&self, n: u128) -> Vec<i128> {
        let (c, b, a) = self.coefficients(n);
        let (mut value, mut gap) = (c, b + a);
        (0..self.count).map(|_| { let v = value; value += gap; gap += 2 * a; v }).collect()
    }
}

fn roots_of_square(n: u128, p: u64) -> Vec<u64> { (0..p).filter(|&r| ((r as u128 * r as u128) % p as u128) == n % p as u128).collect() }

/// Every lane of roots in `[lo, hi)` with `x ≡ phase (mod 4)` and each forced prime dividing x² − n.
pub fn lanes(n: u128, lo: u128, hi: u128, phase: u128, forced: &[u64]) -> Vec<Lane> {
    assert!(phase < 4 && lo > 0 && hi >= lo);
    let mut residues: Vec<u128> = vec![phase];
    let mut step: u128 = 4;
    for &p in forced {
        let roots = roots_of_square(n, p);
        let inverse = mod_inverse((step % p as u128) as u64, p);
        let mut next = Vec::new();
        for &r in &residues { for &t in &roots {
            let shift = ((t as i128 - (r % p as u128) as i128).rem_euclid(p as i128) as u128 * inverse as u128) % p as u128;
            next.push(r + step * shift);
        } }
        residues = next;
        step *= p as u128;
    }
    residues.sort_unstable();
    let mut out = Vec::new();
    for r in residues {
        let first = lo + (r + step - lo % step) % step;
        if first < hi { out.push(Lane { first, step, count: ((hi - 1 - first) / step + 1) as usize, forced: forced.to_vec() }); }
    }
    out
}

fn mod_inverse(a: u64, p: u64) -> u64 {
    let (mut t, mut nt, mut r, mut nr) = (0i128, 1i128, p as i128, (a % p) as i128);
    while nr != 0 { let qd = r / nr; (t, nt) = (nt, t - qd * nt); (r, nr) = (nr, r - qd * nr); }
    assert!(r == 1, "step and prime are not coprime");
    t.rem_euclid(p as i128) as u64
}

/// A smooth root with its exponent vector `(p, e)`, the forced primes included.
pub type SmoothRoot = (u128, Vec<(u64, u32)>);

/// Exact sieve of one lane's cofactors over `base`: the smooth roots with their full exponent
/// vectors `(p, e)`, the forced primes included. Residues of the quadratic are found by
/// enumeration, which also covers primes dividing the step.
pub fn collect(n: u128, lane: &Lane, base: &[u64]) -> Vec<SmoothRoot> {
    let (c, b, a) = lane.coefficients(n);
    let mut rest = lane.cofactors(n);
    let mut exps: Vec<Vec<(u64, u32)>> = vec![lane.forced.iter().map(|&p| (p, 1)).collect(); rest.len()];
    for &p in base {
        let pm = p as i128;
        let roots: Vec<usize> = (0..p as usize).filter(|&r| { let r = r as i128; ((a % pm * (r * r % pm) + b % pm * r + c % pm) % pm).rem_euclid(pm) == 0 }).collect();
        for r in roots {
            let mut i = r;
            while i < rest.len() {
                let mut e = 0;
                while rest[i] % pm == 0 { rest[i] /= pm; e += 1; }
                if e > 0 { match exps[i].iter_mut().find(|(q, _)| *q == p) { Some(t) => t.1 += e, None => exps[i].push((p, e)) } }
                i += p as usize;
            }
        }
    }
    (0..rest.len()).filter(|&i| rest[i] == 1).map(|i| { let mut v = std::mem::take(&mut exps[i]); v.sort_unstable(); (lane.first + lane.step * i as u128, v) }).collect()
}

/// The oracle: the full x² − n by trial division, or None when it is not smooth.
pub fn trial_factor(value: u128, base: &[u64]) -> Option<Vec<(u64, u32)>> {
    let mut v = value; let mut out = Vec::new();
    for &p in base { let mut e = 0; while v % p as u128 == 0 { v /= p as u128; e += 1; } if e > 0 { out.push((p, e)); } }
    if v == 1 { Some(out) } else { None }
}

/// Run every lane over `[lo, hi)`: `(candidates visited, smooth roots with exponents)`.
pub fn run(n: u128, width: u128, phase: u128, forced: &[u64], bound: u64) -> (usize, Vec<SmoothRoot>) {
    let base = primes_upto(bound);
    let lo = isqrt(n) + 1;
    let mut candidates = 0; let mut got = Vec::new();
    for lane in lanes(n, lo, lo + width, phase, forced) { candidates += lane.count; got.extend(collect(n, &lane, &base)); }
    got.sort_unstable_by_key(|(x, _)| *x);
    (candidates, got)
}

fn isqrt(n: u128) -> u128 { let mut x = (n as f64).sqrt() as u128; while x * x > n { x -= 1; } while (x + 1) * (x + 1) <= n { x += 1; } x }

#[cfg(test)]
mod tests {
    use super::*;
    fn check(n: u128, width: u128, phase: u128, forced: &[u64], bound: u64) -> (usize, usize) {
        let base = primes_upto(bound);
        let lo = isqrt(n) + 1;
        let (candidates, got) = run(n, width, phase, forced, bound);
        let expected: Vec<u128> = (lo..lo + width).filter(|&x| x % 4 == phase && forced.iter().all(|&p| (x * x - n) % p as u128 == 0)).collect();
        assert_eq!(candidates, expected.len(), "every eligible root is visited once");
        let mut want = Vec::new();
        for &x in &expected { if let Some(f) = trial_factor(x * x - n, &base) { want.push((x, f)); } }
        assert_eq!(got, want, "the lanes find exactly what trial division finds, n = {n}");
        for (x, f) in &got { let back: u128 = f.iter().map(|&(p, e)| (p as u128).pow(e)).product(); assert_eq!(back, x * x - n); }
        (candidates, got.len())
    }
    #[test]
    fn the_reference_cases_hold() {
        for &n in &[81u128, 8191, 8051, 10403, 1009 * 1013] { for &width in &[0u128, 1, 257] { for phase in 0..4 { for forced in [vec![], vec![3u64], vec![3, 7]] {
            if n == 81 && width == 0 { continue; }
            check(n, width, phase, &forced, 97);
        } } } }
    }
    #[test]
    fn reference_row_8191() {
        let (c0, s0) = check(8191, 16384, 1, &[], 97);
        let (c1, s1) = check(8191, 16384, 1, &[3, 7], 97);
        assert_eq!((c0, s0), (4096, 87));
        assert_eq!((c1, s1), (780, 45));
        let l = lanes(8191, 91, 16475, 1, &[3, 7]);
        let first = l.iter().find(|l| l.first == 97).expect("the lane at 97");
        assert_eq!(first.step, 84);
        assert_eq!(first.coefficients(8191), (58, 776, 336));
    }
}
