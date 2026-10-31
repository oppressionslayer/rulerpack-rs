use rulerpack::*;
use std::time::Instant;
fn main() {
    let t = Instant::now();
    let w = window(1_000_000_000_000, 1_000_000);
    println!("window at 1e12, width 1,000,000: {} survivors in {:.2?}", w.len(), t.elapsed());
    for a in [1_000_000_000_000u64, 1 << 44, 1 << 60] {
        let t = Instant::now();
        let p = next_prime_after(a);
        println!("next prime after {a} = {p}  (+{})  {:.1?}", p - a, t.elapsed());
    }
}
