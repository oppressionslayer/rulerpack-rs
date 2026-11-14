//! `cargo run --release --features big --example snfs_time`: the special sieve on 2^k − 1.
use std::time::Instant;
fn main() {
    for k in [67u32, 71, 79, 101, 113, 137, 149] {
        let t = Instant::now();
        let f = rulerpack::snfs::factor_mersenne(k);
        println!("2^{k}-1 ({} digits): {:?} in {:.1} s", ((k as f64) * std::f64::consts::LOG10_2).ceil(), f.map(|f| f.to_string()), t.elapsed().as_secs_f64());
    }
}
