//! `cargo run --release --example square_lanes`: a square-lane experiment, measured.
//! For each anchor: every root in the span against only the lanes the forced primes allow,
//! candidates visited, smooth Q found, and the time, so the yield per second can be compared.
use rulerpack::lanes::run;
use std::time::Instant;
fn main() {
    let cases: Vec<(u128, Vec<u64>)> = vec![(8191, vec![3, 7]), (8051, vec![5, 7]), (10403, vec![7, 13]), (2045243, vec![7, 19]), (200250077, vec![11, 13]), (1100012300027, vec![7, 11])];
    println!("{:>14} {:>8} {:>12} {:>8} {:>9} {:>10}", "anchor", "forced", "candidates", "smooth", "time", "smooth/s");
    for (n, forced) in &cases {
        for (label, f) in [("none", vec![]), ("yes", forced.clone())] {
            let t = Instant::now();
            let (cands, got) = run(*n, 16384, 1, &f, 97);
            let dt = t.elapsed().as_secs_f64();
            println!("{n:>14} {label:>8} {cands:>12} {:>8} {dt:>8.4}s {:>10.0}", got.len(), got.len() as f64 / dt);
        }
    }
    // a wider, more realistic span: a 13-digit anchor, bound 1000, a million roots
    let n: u128 = 1_000_003 * 1_100_009;
    for (label, f) in [("none", vec![]), ("7,11", vec![7u64, 11]), ("7,11,13", vec![7, 11, 13])] {
        let t = Instant::now();
        let (cands, got) = run(n, 1 << 20, 1, &f, 1000);
        let dt = t.elapsed().as_secs_f64();
        println!("{n:>14} {label:>8} {cands:>12} {:>8} {dt:>8.4}s {:>10.0}", got.len(), got.len() as f64 / dt);
    }
}
