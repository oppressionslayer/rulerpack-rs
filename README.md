# Rulerpack

Prime-number windows, XOR rulers, reversible packing and factorization tools
by Lars Rocha. This directory contains the Rulerpack 0.3.1 source library,
command-line tools, examples and tests.

## Library

```rust
use rulerpack::F513;

let packed = F513.pack(78501, 1000039);
assert_eq!(F513.unpack(packed), (78501, 1000039));
```

The default library includes rulers, packing, prime counting, stored marks and
integer windows. Enable `big` for arbitrary-size integer windows and factoring:

```toml
[dependencies]
rulerpack = { version = "0.3.1", features = ["big"] }
```

## Run an example

From this directory:

```sh
cargo run --release --features big --example factor -- 1009732533765211
```

The separate `square-rulerpack-rs-v2` project uses this library while providing
its own optimized SIQS engine. Use that project's factoring example to reproduce
its performance measurements.

## Features

* Default: rulers, packing, prime counting, marks and fixed-width windows.
* `big`: arbitrary-size integers and factoring using `num-bigint` and Rayon.
* `gmp`: optional GMP-backed tools; requires the system GMP development library.

Use release builds for timings. Probable-prime tests on large integers are not
primality certificates.

## License

MIT. See LICENSE. Source comments identify relevant mathematical algorithms.
