//! `cargo run --release --features big --example siqs_time [max_digits]`: the self-initialising
//! quadratic sieve on the honest fixtures and on the 52-digit random semiprime (randomly sampled) that
//! Dario Alpern's page factors in 1.3 s.
use num_bigint::BigUint;
use std::time::Instant;
fn main() {
    let cases: Vec<&str> = vec![
        "4533471823462034378103601",                            // 25 digits
        "594215563283594108826522371",                          // 30 (approx; replaced below by products)
    ];
    let _ = cases;
    let pairs: Vec<(&str, &str)> = vec![
        ("206158430209", "21990232567889"),
        ("26388279066671", "22517998136864981"),
        ("27021597764223071", "23058430092136951903"),
        ("864691128455135281", "23611832414348226080827"),
        ("885443715538058477627", "24178516392292583494135867"),
        ("53535991745037163321506403", "69560032408973211685368869"),   // random semiprime, 52 digits
    ];
    let want: Option<usize> = std::env::args().nth(1).and_then(|a| a.parse().ok());
    for (p, q) in pairs {
        let (p, q) = (BigUint::parse_bytes(p.as_bytes(), 10).unwrap(), BigUint::parse_bytes(q.as_bytes(), 10).unwrap());
        let n = &p * &q;
        let d = n.to_string().len();
        if want.is_some_and(|w| w < d) { break; }
        let t = Instant::now();
        let f = rulerpack::siqs::factor(&n);
        let ok = f.as_ref().is_some_and(|f| *f == p || *f == q);
        println!("{d} digits: {} in {:.2} s", if ok { "factored" } else { "FAILED" }, t.elapsed().as_secs_f64());
    }
}
