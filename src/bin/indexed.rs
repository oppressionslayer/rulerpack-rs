//! indexed: the next n primes above an anchor with exact indices and packed N.
//!   indexed <anchor> <n> <pi(anchor)>          print index, prime, N
//!   indexed <anchor> <n> <pi(anchor)> --time   just the count and the time
use rulerpack::*;
use std::time::Instant;
fn main() {
    let a: Vec<String> = std::env::args().skip(1).collect();
    if a.len() < 3 { println!("indexed <anchor> <n> <pi(anchor)> [--time]"); return; }
    let (anchor, n, count): (u64, usize, u64) = (a[0].parse().unwrap(), a[1].parse().unwrap(), a[2].parse().unwrap());
    let t = Instant::now();
    let rows = indexed_after(anchor, n, count, MARKS_FAMILY);
    if a.get(3).map(String::as_str) == Some("--time") {
        let (i, p, _) = rows[rows.len() - 1];
        println!("{} indexed primes above {anchor}, last index {i} prime {p}, {:.2?}", rows.len(), t.elapsed());
    } else {
        for (i, p, nn) in rows { println!("{i} {p} {nn}"); }
    }
}
