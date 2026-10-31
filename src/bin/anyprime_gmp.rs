//! anyprime_gmp: a prime above 2^k on GMP with rayon threads.
//!   anyprime_gmp <k> [--prp] [--next]
use rug::Integer;
use rulerpack::gmp::*;
use std::time::Instant;
fn main() {
    let a: Vec<String> = std::env::args().skip(1).collect();
    let k: u32 = a.first().and_then(|s| s.parse().ok()).unwrap_or(2048);
    let certify = if a.iter().any(|s| s == "--prp") { Certify::Prp } else { Certify::Full };
    let next = a.iter().any(|s| s == "--next");
    let anchor = Integer::from(Integer::u_pow_u(2, k));
    let t = Instant::now();
    let p = if next { next_prime_after_gmp(&anchor, certify) } else { any_prime_after_gmp(&anchor, certify) };
    let off = Integer::from(&p - &anchor);
    println!("{} prime after 2^{k} = 2^{k} + {off}   {:.2?}   ({certify:?}, GMP, all cores)", if next { "next" } else { "a" }, t.elapsed());
}
