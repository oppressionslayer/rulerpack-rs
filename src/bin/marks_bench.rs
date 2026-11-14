use rulerpack::*;
use std::time::Instant;
fn main() {
    for every in [1024u64, 256, 64, 16] {
        let t = Instant::now();
        let m = Marks::build(50_000_000, every);
        let tb = t.elapsed();
        let idx = [999u64, 999_999, 9_999_999, 49_999_999];
        let mut line = format!("every {every:>5}: {:>8} marks {:>5.1} MB build {:.1?} |", m.n.len(), m.n.len() as f64 * 16.0 / 1e6, tb);
        for &i in &idx {
            let t = Instant::now();
            let mut s = 0u64;
            for d in 0..200 { s = s.wrapping_add(m.nth_prime(i - d)); }
            let per = t.elapsed().as_secs_f64() * 1000.0 / 200.0;
            line += &format!(" nth({i:>10}) {per:6.3} ms");
            let _ = s;
        }
        println!("{line}");
    }
}
