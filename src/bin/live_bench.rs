use rulerpack::*;
use std::time::Instant;
fn main() {
    let t = Instant::now();
    let mut lm = LiveMarks::new(64, 250_000_000);
    println!("counter ready {:.2?}", t.elapsed());
    for i in [999_999_999u64, 999_999_998, 999_999_937, 99_999_999_999] {
        let t = Instant::now();
        let p = lm.nth_prime(i);
        println!("rust live nth_prime({i}) = {p}  {:.3?}  marks={}", t.elapsed(), lm.len());
    }
    let t = Instant::now();
    println!("rust live pi(1e9) = {}  {:.3?}  marks={}", lm.pi(1_000_000_000), t.elapsed(), lm.len());
    let t = Instant::now();
    println!("rust live pi(1e12) = {}  {:.3?}  marks={}", lm.pi(1_000_000_000_000), t.elapsed(), lm.len());
}
