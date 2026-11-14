//! `cargo run --release --features big --example qs_time [max_digits]`: time the quadratic sieve.
//! The fixtures are semiprimes whose factors differ in size, so the first x² − n is nowhere
//! near a square: the earlier fixtures with factors of equal size were Fermat-easy, which an earlier review
//! pointed out, and measured nothing.
use num_bigint::BigUint;
use std::time::Instant;
fn main() {
    let cases: Vec<(&str, &str)> = vec![
        ("206158430209", "21990232567889"),                       // 25 digits
        ("26388279066671", "22517998136864981"),                   // 30
        ("27021597764223071", "23058430092136951903"),             // 36
        ("864691128455135281", "23611832414348226080827"),         // 41
        ("885443715538058477627", "24178516392292583494135867"),   // 47
    ];
    let want: Option<usize> = std::env::args().nth(1).and_then(|a| a.parse().ok());
    for (p, q) in cases {
        let (p, q) = (BigUint::parse_bytes(p.as_bytes(), 10).unwrap(), BigUint::parse_bytes(q.as_bytes(), 10).unwrap());
        let n = &p * &q;
        let d = n.to_string().len();
        if want.is_some_and(|w| w < d) { break; }
        let t = Instant::now();
        let f = rulerpack::qs::factor(&n);
        let ok = f.as_ref().is_some_and(|f| *f == p || *f == q);
        println!("{d} digits: {} in {:.1} s", if ok { "factored" } else { "FAILED" }, t.elapsed().as_secs_f64());
    }
}
