//! gmp_bench: one strong test and one full test at 8191 bits, single thread; thread count.
use rug::{integer::IsPrime, Integer};
use rulerpack::gmp::*;
use std::time::Instant;
fn main() {
    let a = Integer::from(Integer::u_pow_u(2, 8191));
    let p = Integer::from(&a + 1911u32);           // the known prime
    let c = Integer::from(&a + 11u32);             // a composite survivor-like number
    println!("rayon threads: {}", rayon::current_num_threads());
    let t = Instant::now(); let r = strong_prp(&c, &[2]); println!("strong_prp base 2, composite: {:.1?} -> {r}", t.elapsed());
    let t = Instant::now(); let r = strong_prp(&p, &[2]); println!("strong_prp base 2, prime:     {:.1?} -> {r}", t.elapsed());
    let t = Instant::now(); let r = p.is_probably_prime(25) != IsPrime::No; println!("is_probably_prime(25), prime: {:.1?} -> {r}", t.elapsed());
    let t = Instant::now(); let _ = Integer::from(2u32).pow_mod(&Integer::from(&p - 1u32), &p); println!("one pow_mod:                  {:.1?}", t.elapsed());
}
