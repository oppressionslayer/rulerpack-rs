//! Seed bands, first-window width and certification modes shared by the
//! `big` (num-bigint) and `gmp` (rug) paths. No dependencies.

/// How deep to seed for an anchor of this many bits. Measured against the cost
/// of one certification: the 31 small seeds through 512 bits, 8191 through
/// 1536, then 125 * bits capped at a million (at 8191 bits a certification
/// costs ~100 ms and the deeper sieve still pays).
pub fn seed_depth(bits: usize) -> Option<u64> {
    if bits <= 512 { None }
    else if bits <= 1536 { Some(8191) }
    else { Some((125 * bits as u64).min(1_000_000)) }
}

/// First-window width above 2^64: 7 ln(anchor). The next prime is inside about
/// 999 times in 1000; the caller slides on a miss.
pub fn window_width(bits: usize) -> usize {
    7 * bits * 693 / 1000 + 1
}

/// How a candidate earns the answer.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Certify {
    /// strong test base 2, then 24 more bases (gmpy2.is_prime's shape)
    Full,
    /// strong tests on bases 2 and 3 only: the sort as the answer (Lars Rocha).
    /// For a random survivor the chance a composite passes both is far below
    /// 2^-100 at these sizes; not for inputs built to fool fixed bases.
    Prp,
}

