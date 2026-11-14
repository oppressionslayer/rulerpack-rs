//! `cargo run --release --features big --example siqs_sweep <n> <base> <half_log2> <large> <slack>`:
//! one number, one parameter set, the time and the counts.
use num_bigint::BigUint;
use std::time::Instant;
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let n = BigUint::parse_bytes(a[1].as_bytes(), 10).unwrap();
    let params = rulerpack::siqs::Params { base: a[2].parse().unwrap(), half: 1 << a[3].parse::<u32>().unwrap(), large: a[4].parse().unwrap(), slack: a[5].parse().unwrap(), extra: 64 };
    rulerpack::siqs::stats_reset();
    let t = Instant::now();
    let k = rulerpack::siqs::multiplier(&n);
    let kn = &n * k;
    let base = rulerpack::siqs::factor_base(&kn, params.base);
    let (rels, polys, partials) = rulerpack::siqs::relations(&kn, &base, &params, base.len() + 1 + params.extra);
    let dt = t.elapsed().as_secs_f64();
    // the stage counters of this collection alone, read before anything else runs
    println!("base {} half 2^{} large {} slack {} multiplier {k}: {} relations ({} partials held) from {} polynomials in {dt:.2} s", a[2], a[3], a[4], a[5], rels.len(), partials, polys);
    let threads = rayon::current_num_threads() as f64;
    let st: Vec<u64> = rulerpack::siqs::STATS.iter().map(|s| s.load(std::sync::atomic::Ordering::Relaxed)).collect();
    println!("  thread-seconds: bucket write {:.2}, sieve fill {:.2}, scan {:.2}, exact division {:.2}, switch {:.2}; {} candidates; per thread ≈ {:.2} s of {dt:.2}", st[5] as f64 / 1e9, st[0] as f64 / 1e9, st[1] as f64 / 1e9, st[2] as f64 / 1e9, st[3] as f64 / 1e9, st[4], (st[0] + st[1] + st[2] + st[3] + st[5]) as f64 / 1e9 / threads);
    rulerpack::siqs::stats_reset();
    if std::env::var("SIQS_BUDGET_SECS").is_ok() {
        println!("  budgeted: {:.1} relations/s, {:.0} partial edges/s, {} needed", rels.len() as f64 / dt, partials as f64 / dt, base.len() + 1 + params.extra);
        return;
    }
    let t2 = Instant::now();
    let f = rulerpack::siqs::factor_with(&n, params);
    let dt2 = t2.elapsed().as_secs_f64();
    match f {
        Some(f) => println!("whole factorisation in {dt2:.2} s: {} = {} × {}", n, f, &n / &f),
        None => println!("whole factorisation in {dt2:.2} s: no factor found (prime, or every dependency trivial)"),
    }
}
