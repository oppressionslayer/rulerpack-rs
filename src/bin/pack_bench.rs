use rulerpack::*;
use std::time::Instant;
fn main() {
    let n = 100_000_000u64;
    let t = Instant::now();
    let mut acc = 0u64;
    for i in 0..n { let (a, b) = MARKS_FAMILY.unpack(MARKS_FAMILY.pack(i, 2 * i + 3)); acc = acc.wrapping_add(a ^ b); }
    let dt = t.elapsed();
    println!("rust pack+unpack {n} pairs: {:.2?}  ({:.2} ns per pair)  [{acc}]", dt, dt.as_nanos() as f64 / n as f64);
    let t = Instant::now();
    let mut acc = 0u64;
    for i in 0..n { acc = acc.wrapping_add(ruler(i, 256, 2)); }
    let dt = t.elapsed();
    println!("rust ruler R513 {n} i: {:.2?}  ({:.2} ns each)  [{acc}]", dt, dt.as_nanos() as f64 / n as f64);
}
