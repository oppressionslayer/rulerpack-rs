//! anyprime: a prime above 2^k, the first one any core certifies (not necessarily
//! the next). Pure Rust bignum, window sieve, rayon.
//!   anyprime <k> [--prp] [--next]
use num_bigint::BigUint;
use rulerpack::big::*;
use std::time::Instant;
fn main() {
    let a: Vec<String> = std::env::args().skip(1).collect();
    let k: u64 = a.first().and_then(|s| s.parse().ok()).unwrap_or(2048);
    let certify = if a.iter().any(|s| s == "--prp") { Certify::Prp } else { Certify::Full };
    let next = a.iter().any(|s| s == "--next");
    let anchor = BigUint::from(1u32) << k;
    let t = Instant::now();
    let p = if next { next_prime_after_big_with(&anchor, true, certify) } else { any_prime_after_big(&anchor, certify) };
    let off = &p - &anchor;
    println!("{} prime after 2^{k} = 2^{k} + {off}   {:.2?}   ({:?}, all cores)", if next { "next" } else { "a" }, t.elapsed(), certify);
}
