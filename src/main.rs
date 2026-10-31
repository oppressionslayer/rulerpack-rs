use rulerpack::*;
use std::time::Instant;

fn main() {
    // pack round trip on the known pair
    let n = F513.pack(78501, 1000039);
    println!("pack(78501, 1000039) = {n}  unpack = {:?}", F513.unpack(n));

    println!("\nMeissel (table to 1e8, enough for x <= 1e12):");
    let t = Instant::now();
    let c = Counter::new(100_000_000);
    println!("  table build {:.2?}", t.elapsed());
    for (x, want) in [(1_000_000_000u64, 50847534u64), (10_000_000_000, 455052511), (100_000_000_000, 4118054813), (1_000_000_000_000, 37607912018)] {
        let t = Instant::now();
        let v = c.primepi(x);
        println!("  pi({x:>14}) = {v:>12}  {}  {:.3?}", if v == want { "OK" } else { "WRONG" }, t.elapsed());
    }
    println!("\nnth_prime by R guess + Meissel + walk:");
    for (i, want) in [(999_999u64, 15485863u64), (9_999_999, 179424673), (99_999_999, 2038074743), (999_999_999, 22801763489)] {
        let t = Instant::now();
        let p = c.nth_prime(i);
        println!("  nth_prime({i:>11}) = {p:>12}  {}  {:.3?}", if p == want { "OK" } else { "WRONG" }, t.elapsed());
    }
    println!("\nmarks: sieve 50,000,000 primes once, keep every 1024th:");
    let t = Instant::now();
    let m = Marks::build(50_000_000, 1024);
    println!("  build {} marks in {:.2?}", m.n.len(), t.elapsed());
    for i in [999u64, 999_999, 9_999_999, 49_999_999] {
        let t = Instant::now();
        let p = m.nth_prime(i);
        let t1 = t.elapsed();
        let t = Instant::now();
        let j = m.index_of(p);
        println!("  nth_prime({i:>10}) = {p:>10} {:>9.1?}   index_of -> {j} {:.1?}", t1, t.elapsed());
    }
    let t = Instant::now();
    let mut ok = true;
    for i in (0..50_000_000u64).step_by(997_331) { ok &= m.index_of(m.nth_prime(i)) == i; }
    println!("  round-trip 51 random indices: {}  {:.1?}", if ok { "OK" } else { "FAILED" }, t.elapsed());
}
