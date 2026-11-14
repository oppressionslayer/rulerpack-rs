//! `cargo run --release --features big --example snfs_coupled <k>`: the same Mersenne number,
//! independent screening against the coupled rule, relations and time for each.
use std::time::Instant;
fn main() {
    let k: u32 = std::env::args().nth(1).and_then(|a| a.parse().ok()).unwrap_or(137);
    let (f, m) = rulerpack::snfs::polynomial(k, 3);
    let mut bounds = rulerpack::snfs::Bounds::for_bits(k);
    let (rat, alg, chars) = rulerpack::snfs::bases(&f, &m, bounds.rational, bounds.algebraic, 24);
    let need = 1 + rat.len() + alg.len() + chars.len() + bounds.extra;
    for coupled in [false, true] {
        bounds.coupled = coupled;
        let t = Instant::now();
        let (rels, lines) = rulerpack::snfs::relations(&f, &m, &rat, &alg, need, &bounds, 1.6);
        let dt = t.elapsed().as_secs_f64();
        let pairs: usize = rels.iter().map(|r| r.pairs.len()).sum();
        println!("2^{k}-1 coupled={coupled}: {} relations ({pairs} pairs) from {lines} lines in {dt:.1} s, {:.0} relations/s", rels.len(), rels.len() as f64 / dt);
    }
}
