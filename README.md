# rulerpack (Rust)

XOR rulers over the primes, and the pack that rides on them.
Created by Lars Rocha. Rust implementation with Claude assisting.
No dependencies. The Python version is `rulerpack-lib`.

```rust
use rulerpack::*;

let n = F513.pack(78501, 1000039);     // 8592468642374094_u128
F513.unpack(n);                        // (78501, 1000039)
```

## Credit

The ideas are Lars Rocha's: the XOR ruler and the pack, the anchored seed
sieve, "closest to the anchor", treating the survivors as a list to certify in
parallel, counting once and incrementing, and remembering counts as packed
marks. Claude assisted with the Rust and the measurements; that does not imply
ownership. The segmented sieve, Legendre, Meissel, Miller-Rabin and the li / R
estimates are classical. primecount is Kim Walisch's and is the yardstick in
the tables.

## The idea

A ruler is a Mersenne ladder `m ^ (m+1)` scaled by a power of two: the tick
marks of a tape measure, taller at every power of two. For a prime `p` at
index `i`:

```
K = 2p ^ ruler(i)         X = 2p ^ 256 i         N = K << 40 | X
K ^ X = 256 (i + 2)   so   i = (K ^ X)/256 - 2,   p = (X ^ 256 i)/2
```

One integer carries the prime and its index, and comes apart with a shift,
a mask and two XORs.

## What is in the crate

| item | what |
|---|---|
| `ladder`, `rp`, `ruler` | the rulers |
| `Family`, `F513`, `F257`, `MARKS_FAMILY` | pack and unpack on `u128` |
| `isprime` | deterministic Miller-Rabin, exact for all `u64` |
| `PiTable` | odd bit sieve with a rank index: `pi(n)` in O(1) after the build |
| `Counter` | Meissel with pi-table pruning: `primepi`, `nth_prime`, `index_of` |
| `li`, `riemann_r` and their inverses | the guess for `nth_prime` |
| `Marks` | the checkpoint table: build, save, load, `nth_prime`, `index_of`, `pi`, `verify` |
| `LiveMarks` | the same table filled on demand: a missed block costs one counting call, then it is a walk |

## Marks

Sieve once, keep every k-th prime as a packed `N`, walk from the nearest mark
with the seed rule and `isprime`. Each mark carries its own index and prime,
so the file is self-describing and the loader checks it against itself.

```
cargo build --release
target/release/marks build 50000000 64 marks.bin    # 3.2 s, 781,250 marks, 12.5 MB
target/release/marks nth 49999999 marks.bin         # 982451653
target/release/marks index 982451653 marks.bin      # 49999999
target/release/marks bench marks.bin
target/release/marks verify marks.bin               # Miller-Rabin on every mark
```

Measured on one thread, averaged over 200 indices each, against
[primecount](https://github.com/kimwalisch/primecount) 8.1's `nth_prime`:

| index | primecount | marks every 64 | marks every 16 |
|---|---|---|---|
| 999 | 0.001 ms | 0.058 ms | 0.013 ms |
| 999,999 | 0.081 ms | 0.122 ms | 0.027 ms |
| 9,999,999 | 0.271 ms | 0.138 ms | 0.032 ms |
| 49,999,999 | 0.470 ms | 0.153 ms | 0.035 ms |

primecount computes; marks remembers. Its time grows with the prime, ours is
the walk length, which you choose with the spacing. The table at every 16 is
50 MB. Outside the table there is no answer, which is the trade.

## The window: Eratosthenes from any anchor

```rust
window(1_000_000_000_000, 1_000_000);   // survivor offsets, closest first, 2.7 ms
next_prime_after(1 << 60);              // 1152921504606847009 in 24 us
primes_between(10u64.pow(15), 10u64.pow(15) + 10_000_000);   // 289,394 primes, 161 ms
```

`window` strikes with the 31 seeds and leaves `isprime` the rest.
`primes_between` brings every prime up to the square root of the top, so every
survivor is prime with no test. Counts agree with primecount at 10^12, 10^15, 10^18.

```
target/release/between 1000000000000 1000000000100
target/release/between 1000000000000000000 1000000000001000000 --count   # 24,280 primes, 3.7 s
```

### Big anchors (`--features big`)

Same window at any size, pure Rust bignum (num-bigint) and rayon. The first
window is 7 ln(anchor) wide; seed depth follows the cost of a test (31 seeds
through 512 bits, 8191 through 1536, up to a million above). Two searches:

- `next_prime_after_big`: every core, smallest pass wins (`find_first`)
- `any_prime_after_big`: every core, first pass from any core wins (`find_any`).
  A prime above the anchor, not necessarily the next. For "a random prime of
  this size" there is no stop rule to wait on.

`Certify::Full` is base 2 then 24 more bases; `Certify::Prp` is bases 2 and 3
only, the sort as the answer (Lars Rocha): for a random survivor the chance a
composite passes both is far below 2^-100 at these sizes.

```
target/release/anyprime 4096            # a prime above 2^4096, all cores
target/release/anyprime 8191 --prp      # bases 2 and 3 as the answer
target/release/anyprime 8191 --next     # the next prime, all cores
```

### GMP in Rust (`--features gmp`, needs libgmp-dev)

The same two searches on GMP's modular power (the `rug` crate, linked to the
system libgmp) with rayon threads instead of Python processes. Candidates go
out in rounds of 2 x threads, lowest first, so every core works the low end
of the window together; rayon's plain work stealing would spread the threads
over the whole window and leave the low candidates to one thread in sequence
(measured: 5 s instead of 1.2 s at 8191 bits).

```
target/release/anyprime_gmp 8191            # a prime above 2^8191
target/release/anyprime_gmp 12191 --next    # the next prime above 2^12191
```

| anchor | gmpy2.next_prime | Python any_prime_after (16 procs) | Rust GMP any | Rust GMP next | Rust num-bigint next |
|---|---|---|---|---|---|
| 2^2048 | 0.15 s | 0.07 s | 0.03 s | 0.03 s | 0.61 s |
| 2^4096 | 1.29 s | 0.43 s | 0.21 s | 0.21 s | 3.3 s |
| 2^8191 | 7.84 s | 2.91 s | 1.23 s | 1.22 s | 42.9 s |
| 2^12191 | 61.6 s | 14.7 s | 8.6 s | 8.9 s | |

(8 cores / 16 threads, i7-10875H.) Every answer in the `next` column is the
same prime gmpy2 returns. Stated plainly: one GMP modular power is 67 ms at
8191 bits in this build against 103 ms through gmpy2, and threads share the
window without copying it, which is the rest of the gap. The pure-Rust
`big` path (num-bigint) is kept for machines without GMP and is several times
slower, because its modular power is. The u64 path needs neither.

## Live marks

```rust
let mut lm = LiveMarks::new(64, 100_000_000);   // counter sized for primes to about 1e12
lm.nth_prime(999_999_999);    // 22801763489 in 106 ms, one mark computed and stored
lm.nth_prime(999_999_998);    // 22801763477 in 0.9 ms, same mark
lm.pi(1_000_000_000);         // 50847534 in 21 ms, five marks
lm.pi(1_000_000_000_000);     // 37607912018 in 11.6 s, ten marks (each a Meissel at 1e12)
```

Nothing is built up front and the file grows only where you have been.
The counter is the ceiling: past what Meissel reaches here, a mark would be
filled by primecount instead, and the ruler owns it after.

## The fair race: C++ against Rust, same job, one thread

Both sides are given pi(anchor) once (from primecount). Then C++ lists the
range with primesieve; Rust runs `indexed_after`, which sieves the same range
and also attaches the exact index and the packed N to every prime.

| anchor | n | primesieve list | Rust `indexed_after` with indices and packs |
|---|---|---|---|
| 10^12 | 10,000 | 5 ms | 11 ms |
| 10^14 | 1,000,000 | 49 ms | 597 ms |
| 10^15 | 100,000 | 35 ms | 168 ms |

primesieve is 2x to 12x faster at the listing itself; it has fifteen years of
cache tuning in it and this sieve has none. What it does not produce is the
index and the pack, and the only route to indices in those tools, `nth_prime`
per prime, is thousands of times slower than either column.

```
target/release/indexed 1000000000000 10000 37607912018 --time
```

## Baselines, stated exactly

primecount is the yardstick for counting and for `nth_prime` by index, and it
wins at counting. primesieve, from the same author, is the yardstick for
listing primes in a range, and it wins there too: a million primes above 10^14
in 0.04 s against 0.22 s for `between`. The marks table and the pack are what
this crate adds: random access by index with no counting, and one integer that
carries the pair.

## Which to use

| tool | holds | cost per call | reaches |
|---|---|---|---|
| `Counter` (Meissel) | nothing | grows with x, 11 ms at 1e9, 2.1 s at 1e12 | about 1e12 |
| `Marks` | a table built once | flat, the walk length, about 100 µs at every 64 | the end of the table |
| `LiveMarks` | only touched blocks | first touch counts, after that the walk | wherever the counter reaches |
| primecount | nothing | 20 ms at 1e12, one thread | 1e20 and beyond |

## Counting

`Counter` is plain Meissel with one pruning rule and no threads. It is here
so the pack can be checked from nothing, not to compete.

| π(x) | this crate | primecount, 1 thread |
|---|---|---|
| 10^9 | 11 ms | under 1 ms |
| 10^12 | 2.1 s | 20 ms |

If you need counting, use primecount.

## Check

```
cargo test --release
cargo clippy --release --all-targets -- -D warnings
```

MIT license.
