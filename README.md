# MPA-OpenCl

### OpenCL MPA API

===========================================================================

DESCRIPTION:

MPA-OpenCL is a Multiple precision Arithmetic API in OpenCl licensed under the Apache Software License v2.
The API offer a set of helper functions that can be used to carry usual Arbitrary Precision Arithmetic in every openCL enable device like GPUs, Multi core CPUs, Co-processeurs, FPGA and hand held devices that support OpenCL like for example Android Devices supporting OpenCL.
The main motivation behind the development of this API come from the lack of such an API in OpenCL, similar API exists for proprietary GPGPU platform like CUDA.
The API can be used to accelerate applications using multiple precision arithmetics like ECDSA, RSA, Research in physics, Big Data analysis Applications to name a few.

List of supported multiple precision arithmetic operations in the MPA-OpenCl version 1.0-beta

Fifteen operators, spanning thirteen operations: two of the operations have two
operators each, computing the same result by different routes, and the pair
shares a row number below. Operator codes run 1 to 15 with no gaps.

All fifteen are implemented in `mpaKernel_32bits_opt.cl` and verified against
GMP by `mpa_test`. The three original kernels (`mpaKernel_32bits.cl`,
`mpaKernel_16bits.cl`, `mpaKernels_8bits.cl`) implement seven of them.

| # | Operation | Code | Kernels |
|---|---|---|---|
| 1 | Big number comparison | `COMPARE` = 8 | opt |
| 2 | Big number addition | `ADD` = 1 | all |
| 3 | Big number subtraction | `SUBTRACT` = 2 | all |
| 4 | Big number multiplication, operand scanning | `MULTIPLYOPERANDSCANNING` = 5 | all |
|   | Big number multiplication, product scanning | `MULTIPLYPRODUCTSCANNING` = 6 | all |
| 5 | Big number exponentiation | `EXPONENTIATION` = 12 | opt |
| 6 | Big number division | `DIVIDE` = 13 | opt |
| 7 | Big number integer square root | `ISQRT` = 14 | opt |
| 8 | Big number reduction | `REDUCE` = 9 | opt |
| 9 | Modular addition | `ADDMOD` = 3 | all |
| 10 | Modular subtraction | `SUBTRACTMOD` = 4 | all |
| 11 | Modular multiplication | `MODMUL` = 10 | opt |
|    | Modular multiplication, R² variant | `MODMUL_R2` = 15 | opt |
| 12 | Montgomery multiplication | `MONTGOMERYMULTIPLICATION` = 7 | all |
| 13 | Modular exponentiation | `MODEXP` = 11 | opt |

## Data layout and conventions

An operand is `WORDLENGTH_T` words, **most significant word first**, each word
`w` bits wide (`w` = 8, 16 or 32 depending on the kernel). `WORDLENGTH_T` is set
at build time with `-DWORDLENGTH_T=N`, so `BITSLENGTH = N * w`.

The default layout places item `g`'s word `i` at `g*T + i`. The optimized
kernel also supports a word-interleaved layout at `i*N + g`
(`-DMPA_INTERLEAVED=1`) so that adjacent work-items touch adjacent addresses.

Kernel arguments are `(input1, input2, output, params, modulus)`, where
`params` is four ints: `{operator, wordsize, bitslength, m'}` and `m'` is
`-p^-1 mod 2^w`, computed by the host from the modulus in use.

Per-operation notes:

- `ADD`, `SUBTRACT`, `EXPONENTIATION` wrap to the operand width, i.e. the
  result is taken mod `2^(w*T)`. The multiplications instead produce the full
  `2*T`-word product.
- `COMPARE` writes a `T`-word two's-complement value: `0` if `a == b`, `1` if
  `a > b`, and all-ones (`-1`) if `a < b`.
- `DIVIDE` writes `2*T` words: the quotient in words `0..T-1` and the remainder
  in words `T..2T-1`. Read as one big-endian `2T`-word integer this is
  `q * 2^(w*T) + r`. Division by zero yields `q = 0, r = 0`.
- `REDUCE`, `MODMUL`, `MODEXP` and the modular add/subtract use the modulus
  supplied in the `modulus` buffer. `ADDMOD`, `SUBTRACTMOD`, `MODMUL` and
  `MODEXP` require operands already less than the modulus; `REDUCE` does not.
- The modulus must be odd for the Montgomery-based operations. It does not have
  to be prime — RSA-style composite moduli are supported and tested.
- The `modulus` buffer is `2*T` words: the modulus in words `0..T-1`, then
  `R^2 mod p` in words `T..2T-1`, where `R = 2^(w*T)`. Only `MODMUL_R2` reads
  the second half; every other operation ignores it.
- `MODMUL` and `MODMUL_R2` compute the same thing by different routes.
  `MODMUL` enters the Montgomery domain per item with a binary reduction plus
  `w*T` modular doublings. `MODMUL_R2` instead uses the host-supplied `R^2` and
  two CIOS calls, `Mont(Mont(a,b), R^2)`. Both are kept so the cost of the
  domain transfer can be measured rather than assumed; `MODMUL_R2` measured
  17x faster on a CPU OpenCL device.

## Using the kernels from your own host

`mpa_run` is a complete, dependency-free host: standard C99 and an OpenCL ICD,
nothing else. It is the worked reference for what an application actually has to
do to drive the kernels, and it is small enough to read in one sitting.

```sh
./mpa_run <op> <p-hex> <a-hex> <b-hex> [<a-hex> <b-hex> ...]
./mpa_run <op> <p-hex> <a-hex> [<a-hex> ...]          # reduce, isqrt
make run RUNARGS="modmul <p-hex> <a-hex> <b-hex>"
```

Operand pairs beyond the first become additional work-items in the same launch.
The modulus fixes the width — `T = ceil(hex digits / 8)` words — so one binary
drives any size; operators that ignore the modulus still take it as the width
argument. Results print one per line, in hex.

Op names: `add` `sub` `addmod` `submod` `mul` `mulop` `montmul` `compare`
`reduce` `modmul` `modexp` `exp` `div` `isqrt` `modmul_r2`.

```sh
P=ffffffff00000001000000000000000000000000ffffffffffffffffffffffff
./mpa_run modmul $P 1234567890abcdef fedcba0987654321
./mpa_run addmod $P a b  c d  e f                       # three items, one launch
MPA_KERNEL=mpaKernel_32bits.cl ./mpa_run montmul $P a b  # choose the kernel
MPA_BUILD="-DMPA_MULHI=1 -DMPA_FUSED_CIOS=1" ./mpa_run montmul $P a b
```

`mpa_run` assumes the contiguous default layout, so `-DMPA_INTERLEAVED=1` is
rejected at startup rather than answered wrongly.

### What the host has to compute

Everything the kernel reads is a word array, so the only host-side arithmetic is
the two derived values described under *Data layout and conventions*. Neither
needs a bignum library.

`m'` depends only on the modulus mod `2^w`, which makes it a single-word
computation — Hensel lifting from a seed correct to 3 bits:

```c
uint32_t x = p0;                               /* p0 = lowest word, odd */
for (int i = 0; i < 4; i++) x *= 2u - p0 * x;  /* 3 -> 6 -> 12 -> 24 -> 48 bits */
return -x;                                     /* -p^-1 mod 2^32 */
```

`R^2 mod p` is the one genuinely multi-precision quantity, and only `MODMUL_R2`
reads it: `2^(2wT) mod p` by `2wT` modular doublings over word arrays, roughly
fifteen lines. Every other Montgomery path builds the domain in-kernel, so a
host that sticks to `MODMUL` and `MODEXP` computes nothing beyond `m'`.

## Building and testing

```sh
make            # build everything
make test       # GMP-vs-device correctness harness, exits non-zero on failure
make perf       # optimization sweep, each level verified before it is timed
make compare    # head-to-head against GMP and OpenSSL on identical operands
make bench      # the three original timing hosts
make ecdsa      # batched P-256 ECDSA verification, validated against OpenSSL
make run        # minimal host, no GMP and no OpenSSL
```

`ecdsa_bench` generates fresh P-256 key pairs and signatures with OpenSSL,
deliberately corrupts a configurable fraction of them, and checks the device
verdict against OpenSSL's on every signature. It reports false accepts and false
rejects separately: a verifier that accepts everything cannot pass, and a false
accept voids the timing results.

`mpa_compare` reports four OpenCL rows, as a 2x2 over memory path and what is
timed: kernel alone and end-to-end, each with default device buffers
(`clCreateBuffer` plus explicit write/read) and with host-mapped buffers
(`CL_MEM_ALLOC_HOST_PTR` plus map/unmap, no copies). Both kernel rows run
identical code, so comparing them isolates the cost of the memory path from the
cost of the arithmetic. Zero-copy wins on unified memory (Apple silicon,
integrated GPUs) and usually loses on a discrete GPU, where the kernel would
read operands over PCIe instead of from device memory. Every row is verified
against GMP before it is timed.

Timings vary run to run by 15% or more on a loaded machine. The default is
minimum-of-9; raise `--reps` for numbers you intend to publish.

Selecting a device:

```sh
MPA_LIST_DEVICES=1 ./mpa_test      # enumerate platforms and devices
MPA_DEVICE_TYPE=gpu make test      # require a GPU; fail rather than fall back
MPA_DEVICE_TYPE=cpu ./mpa_test     # force the CPU
MPA_DEVICE_INDEX=1 ./mpa_test      # second matching device
```

`mpa_test` covers every operation against five moduli (four primes plus one
RSA-style composite) at three word sizes, with random operands plus directed
edge cases, comparing every output word against GMP.

The composite modulus is not optional. `python3 reach_check.py` exhaustively
evaluates a reference CIOS model over every odd modulus below 64 at five
`(w, T)` settings and shows that the Montgomery post-condition `A == m` is
reachable only for composite moduli — so a suite restricted to primes cannot
detect a wrong conditional subtraction there. It runs in about a second and
exits non-zero if the property fails.

## Kernel dependencies

The kernels use standard OpenCL C as described in the Khronos specifications
and have no external dependencies. `mpaKernel_32bits_opt.cl` is fully
self-contained and needs no `-I`.

## Host dependencies

The host application also uses standard C99 and does not depend on external
libraries. `mpa_run` is the demonstration: it links an OpenCL ICD and libc and
nothing more, and covers the whole host side of the API.

The test and benchmark hosts are a separate matter. Their job is to check the
device against an independent implementation, so they need GMP, and the ones
that generate random vectors or cross-check signatures also need OpenSSL. None
of that is required to *use* the API.

| Host | Purpose | Needs |
|---|---|---|
| `mpa_run` | drive the kernels | C99 + OpenCL |
| `mpa_test` | correctness harness | + GMP |
| `mpa_bench` | optimization sweep | + GMP |
| `mpa_compare` | head-to-head timings | + GMP, OpenSSL |
| `ecdsa_bench` | batched ECDSA | + GMP, OpenSSL |
| `mpa_8bits`, `mpa_16bits`, `mpa_32bits` | original timing hosts | + GMP, OpenSSL |

```sh
# Debian/Ubuntu - mpa_run alone
sudo apt install ocl-icd-opencl-dev pocl-opencl-icd
# ...and the harnesses
sudo apt install libgmp-dev libssl-dev
# macOS (OpenCL ships with the OS)
brew install gmp openssl
```
