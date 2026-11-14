//! `cargo run --release --features big --example factor <n>`: the complete factorisation of n,
//! small primes by trial division, perfect powers, then the self-initialising quadratic sieve
//! on every composite piece until all are prime.
use num_bigint::BigUint;
use num_traits::One;
use std::time::Instant;

fn main() {
    let arg = std::env::args().nth(1).expect("a number");
    // progress lines with an estimate on stderr, unless the caller turned them off
    if std::env::var("SIQS_PROGRESS").is_err() { std::env::set_var("SIQS_PROGRESS", "1"); }
    let n = BigUint::parse_bytes(arg.as_bytes(), 10).expect("a decimal number");
    let t = Instant::now();
    let factors = rulerpack::siqs::factor_completely(&n);
    let shown: Vec<String> = factors.iter().map(|f| f.to_string()).collect();
    println!("{} = {}", arg, shown.join(" × "));
    let back = factors.iter().fold(BigUint::one(), |acc, f| acc * f);
    println!("{} prime factors, product checks: {}, {:.2} s", factors.len(), back.to_string() == arg, t.elapsed().as_secs_f64());
}
