//! nextbig: next prime above 2^k, any k. Pure Rust bignum, window sieve, all cores.
//!   nextbig <k> [--serial]
use num_bigint::BigUint;
use rulerpack::big::*;
use std::time::Instant;
fn main() {
    let a: Vec<String> = std::env::args().skip(1).collect();
    let k: u64 = a.first().and_then(|s| s.parse().ok()).unwrap_or(2048);
    let parallel = !a.iter().any(|s| s == "--serial");
    let anchor = BigUint::from(1u32) << k;
    let t = Instant::now();
    let p = next_prime_after_big(&anchor, parallel);
    let off = &p - &anchor;
    println!("next prime after 2^{k} = 2^{k} + {off}   {:.2?}   ({})", t.elapsed(), if parallel { "all cores" } else { "one core" });
}
