//! between: every prime in [lo, hi), sieve of Eratosthenes started at lo.
//!   between <lo> <hi>            print them
//!   between <lo> <hi> --count    just how many, and the time
use rulerpack::*;
use std::time::Instant;
fn main() {
    let a: Vec<String> = std::env::args().skip(1).collect();
    if a.len() < 2 { println!("between <lo> <hi> [--count]"); return; }
    let (lo, hi): (u64, u64) = (a[0].parse().unwrap(), a[1].parse().unwrap());
    let t = Instant::now();
    let p = primes_between(lo, hi);
    if a.get(2).map(String::as_str) == Some("--count") {
        println!("{} primes in [{lo}, {hi})  first {:?} last {:?}  {:.2?}", p.len(), p.first(), p.last(), t.elapsed());
    } else {
        for x in p { println!("{x}"); }
    }
}
