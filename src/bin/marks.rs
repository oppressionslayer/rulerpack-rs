//! marks: Lars Rocha's ruler checkpoints.
//!
//!   marks build <n_primes> [every] [file]   sieve once, save packed marks (default every 64, file marks.bin)
//!   marks nth <i> [file]                    prime at zero-based index i
//!   marks index <p> [file]                  zero-based index of prime p
//!   marks pack <i> <p>                      one integer for the pair
//!   marks unpack <N>                        the pair back
//!   marks verify [file]                     Miller-Rabin every mark
//!   marks bench [file]                      time lookups across the table
use rulerpack::*;
use std::time::Instant;

fn main() {
    let a: Vec<String> = std::env::args().skip(1).collect();
    let file = |k: usize| a.get(k).cloned().unwrap_or_else(|| "marks.bin".into());
    match a.first().map(String::as_str) {
        Some("build") => {
            let n: u64 = a[1].parse().expect("n_primes");
            let every: u64 = a.get(2).and_then(|s| s.parse().ok()).unwrap_or(64);
            let t = Instant::now();
            let m = Marks::build(n, every);
            let tb = t.elapsed();
            m.save(&file(3)).expect("save");
            println!("{} marks every {every} primes, last index {}, built in {:.2?}, {:.1} MB -> {}",
                m.len(), m.last_index(), tb, m.len() as f64 * 16.0 / 1e6, file(3));
        }
        Some("nth") => {
            let i: u64 = a[1].parse().expect("index");
            let m = Marks::load(&file(2)).expect("load");
            let t = Instant::now();
            let p = m.nth_prime(i);
            println!("{p}   ({:.1?})", t.elapsed());
        }
        Some("index") => {
            let p: u64 = a[1].parse().expect("prime");
            let m = Marks::load(&file(2)).expect("load");
            let t = Instant::now();
            let i = m.index_of(p);
            println!("{i}   ({:.1?})", t.elapsed());
        }
        Some("pack") => {
            let (i, p): (u64, u64) = (a[1].parse().unwrap(), a[2].parse().unwrap());
            let (k, x) = F513.encode(i, p);
            println!("K = {k}\nX = {x}\nN = {}", F513.pack(i, p));
        }
        Some("unpack") => {
            let n: u128 = a[1].parse().unwrap();
            let (i, p) = F513.unpack(n);
            println!("i = {i}\np = {p}");
        }
        Some("verify") => {
            let m = Marks::load(&file(1)).expect("load");
            let t = Instant::now();
            println!("{} marks, all prime: {}  ({:.2?})", m.len(), m.verify(), t.elapsed());
        }
        Some("bench") => {
            let t = Instant::now();
            let m = Marks::load(&file(1)).expect("load");
            println!("loaded {} marks every {} in {:.2?}", m.len(), m.every, t.elapsed());
            let last = m.last_index();
            let mut i = 999u64;
            while i <= last {
                let t = Instant::now();
                let mut s = 0u64;
                for d in 0..200 { s = s.wrapping_add(m.nth_prime(i - d)); }
                let per = t.elapsed().as_secs_f64() * 1e6 / 200.0;
                let p = m.nth_prime(i);
                let t = Instant::now();
                for d in 0..200 { s = s.wrapping_add(m.index_of(m.nth_prime(i - d))); }
                let per2 = t.elapsed().as_secs_f64() * 1e6 / 200.0 - per;
                println!("  nth_prime({i:>12}) = {p:>13}  {per:7.1} us   index_of {per2:7.1} us");
                let _ = s;
                i = i * 10 + 9;
            }
        }
        _ => println!("{}", include_str!("marks.rs").lines().take_while(|l| l.starts_with("//!")).map(|l| l.trim_start_matches("//! ").trim_start_matches("//!")).collect::<Vec<_>>().join("\n")),
    }
}
