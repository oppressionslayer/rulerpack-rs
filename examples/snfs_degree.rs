//! `cargo run --release --features big --example snfs_degree <k> <degree>`: one Mersenne number, one degree.
use std::time::Instant;
fn main() {
    let k: u32 = std::env::args().nth(1).and_then(|a| a.parse().ok()).unwrap_or(173);
    let d: u32 = std::env::args().nth(2).and_then(|a| a.parse().ok()).unwrap_or(5);
    let t = Instant::now();
    let f = rulerpack::snfs::factor_mersenne_degree(k, d);
    println!("2^{k}-1 degree {d}: {:?} in {:.1} s", f.map(|f| f.to_string()), t.elapsed().as_secs_f64());
}
