//! Ruler stack: how close can a pile of rulers get to naming the next prime?
//! The Rust twin of spindexinz/sandbox/xor-ruler/ruler_stack.py; the two print the same table.
//! `cargo run --release --example ruler_stack`
use rulerpack::PiTable;
use std::collections::HashMap;

const N: u64 = 2_000_000;
const SCORE: usize = 50_000;

/// run(i) ^ run(i+1): the carry ruler at index i
fn carry(i: u64) -> u64 { ((i + 1) ^ i) ^ ((i + 2) ^ (i + 1)) }

fn score(name: &str, p: &[u64], q: &[u64], idx: &[u64], order: &dyn Fn(u64, u64) -> Vec<u64>) {
    let (mut top1, mut top3) = (0, 0);
    for k in 0..SCORE {
        let c = order(p[k], idx[k]);
        if c[0] == q[k] { top1 += 1; }
        if c.iter().take(3).any(|&x| x == q[k]) { top3 += 1; }
    }
    let m = SCORE as f64;
    println!("{name:34} first guess {:6.1}%   top-3 {:6.1}%", 100.0 * top1 as f64 / m, 100.0 * top3 as f64 / m);
}

fn wheel(p: u64, g: u64, ps: &[u64]) -> bool { ps.iter().all(|&q| (p + g) % q != 0) }

/// gaps ordered by how often they occur, most common first
fn by_frequency(gaps: impl Iterator<Item = u64>) -> Vec<u64> {
    let mut f: HashMap<u64, usize> = HashMap::new();
    for g in gaps { *f.entry(g).or_default() += 1; }
    let mut v: Vec<(u64, usize)> = f.into_iter().collect();
    v.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
    v.into_iter().map(|(g, _)| g).collect()
}

fn main() {
    stack();
    wheels_out_to_23();
}

fn stack() {
    let primes = PiTable::new(N).primes_upto(N - 1);
    let p = &primes[3..primes.len() - 1]; // from 7 on
    let q = &primes[4..];
    let gaps: Vec<u64> = p.iter().zip(q).map(|(a, b)| b - a).collect();
    let idx: Vec<u64> = (3..3 + p.len() as u64).collect();
    let w210 = [3u64, 5, 7];

    score("odds", p, q, &idx, &|p, _| (2..60).step_by(2).map(|g| p + g).collect());
    score("wheel 2·3·5", p, q, &idx, &|p, _| (2..60).step_by(2).filter(|&g| wheel(p, g, &[3, 5])).map(|g| p + g).collect());
    score("wheel 2·3·5·7", p, q, &idx, &|p, _| (2..60).step_by(2).filter(|&g| wheel(p, g, &w210)).map(|g| p + g).collect());

    let g_all = by_frequency(gaps.iter().copied());
    let freq = |p: u64, order: &[u64]| -> Vec<u64> { order.iter().filter(|&&g| wheel(p, g, &w210)).map(|g| p + g).collect() };
    score("gap frequency + wheel 210", p, q, &idx, &|p, _| freq(p, &g_all));

    score("carry twin/sexy + gap freq + wheel", p, q, &idx, &|p, i| {
        let c = carry(i);
        let mut first = Vec::new();
        if p ^ (p + 2) == c { first.push(p + 2); }
        if p ^ (p + 6) == c { first.push(p + 6); }
        let base = freq(p, &g_all);
        first.iter().copied().chain(base.into_iter().filter(|x| !first.contains(x))).collect()
    });

    // how much does the carry know? P(gap=2 | carry matches p^(p+2)) vs P(gap=2)
    for (g, label) in [(2u64, "2"), (6, "6")] {
        let matches: Vec<usize> = (0..p.len()).filter(|&k| p[k] ^ (p[k] + g) == carry(idx[k])).collect();
        let base = gaps.iter().filter(|&&x| x == g).count() as f64 / gaps.len() as f64;
        let cond = matches.iter().filter(|&&k| gaps[k] == g).count() as f64 / matches.len() as f64;
        println!("\nP(gap={label})                        = {base:.3}");
        println!("P(gap={label} | carry matches p^(p+{label})) = {cond:.3}   (matches: {:.1}% of indices)", 100.0 * matches.len() as f64 / p.len() as f64);
    }

    // residue-conditioned gap frequency, learned on the primes past the score set
    println!();
    for m in [6u64, 30, 210] {
        let mut table: HashMap<u64, Vec<u64>> = HashMap::new();
        let mut per: HashMap<u64, Vec<u64>> = HashMap::new();
        for k in SCORE..p.len() { per.entry(p[k] % m).or_default().push(gaps[k]); }
        for (r, gs) in per { table.insert(r, by_frequency(gs.into_iter())); }
        score(&format!("gap freq | p mod {m}, + wheel 210"), p, q, &idx, &|p, _| freq(p, table.get(&(p % m)).unwrap_or(&g_all)));
    }

    println!("\ngap order after p mod 30 (learned on held-out primes):");
    let mut per: HashMap<u64, Vec<u64>> = HashMap::new();
    for k in SCORE..p.len() { per.entry(p[k] % 30).or_default().push(gaps[k]); }
    let mut rs: Vec<u64> = per.keys().copied().collect();
    rs.sort();
    for r in rs {
        let gs = &per[&r];
        let n = gs.len() as f64;
        let line: Vec<String> = by_frequency(gs.iter().copied()).iter().take(4)
            .map(|&g| format!("gap {g} {:5.1}%", 100.0 * gs.iter().filter(|&&x| x == g).count() as f64 / n)).collect();
        println!("  p ≡ {r:2} mod 30: {}", line.join(", "));
    }

    // bigger wheels
    println!();
    for ps in [vec![3u64, 5, 7], vec![3, 5, 7, 11], vec![3, 5, 7, 11, 13]] {
        let name = format!("wheel 2·{}", ps.iter().map(|x| x.to_string()).collect::<Vec<_>>().join("·"));
        score(&name, p, q, &idx, &|p, _| (2..120).step_by(2).filter(|&g| wheel(p, g, &ps)).map(|g| p + g).collect());
    }
}

/// The final ruler: a primality test on each candidate in wheel order. The test makes the answer
/// certain; what the wheel buys is how few candidates need testing. Mean, median and worst count.
fn tests_needed(p: &[u64], q: &[u64], ps: &[u64]) -> (f64, u64, u64) {
    let mut ranks: Vec<u64> = (0..SCORE).map(|k| {
        let (mut g, mut n) = (2u64, 1u64);
        while p[k] + g != q[k] {
            if wheel(p[k], g, ps) { n += 1; }
            g += 2;
        }
        n
    }).collect();
    ranks.sort_unstable();
    (ranks.iter().sum::<u64>() as f64 / SCORE as f64, ranks[SCORE / 2], ranks[SCORE - 1])
}

fn wheels_out_to_23() {
    let primes = PiTable::new(N).primes_upto(N - 1);
    let p = &primes[3..primes.len() - 1];
    let q = &primes[4..];
    let idx: Vec<u64> = (3..3 + p.len() as u64).collect();
    println!();
    for ps in [vec![], vec![3u64], vec![3, 5], vec![3, 5, 7], vec![3, 5, 7, 11], vec![3, 5, 7, 11, 13], vec![3, 5, 7, 11, 13, 17], vec![3, 5, 7, 11, 13, 17, 19, 23]] {
        let name = if ps.is_empty() { "odds".to_string() } else { format!("wheel 2·{}", ps.iter().map(|x| x.to_string()).collect::<Vec<_>>().join("·")) };
        score(&name, p, q, &idx, &|p, _| (2..120).step_by(2).filter(|&g| wheel(p, g, &ps)).map(|g| p + g).collect());
        let (mean, med, mx) = tests_needed(p, q, &ps);
        println!("{:34} candidates tested until the next prime: mean {mean:.2}  median {med}  worst {mx}", "");
    }
}
