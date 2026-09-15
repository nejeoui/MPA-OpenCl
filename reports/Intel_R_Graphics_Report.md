# MPA-OpenCL benchmark report - Intel(R) Graphics

> **Note.** The multi-threaded GMP and OpenSSL baseline columns have been
> removed from this report: they predate the 2026-09-12 timing fix and were
> understated (see `reports/README.md`). The single-threaded GMP column and
> every MPA measurement are unaffected, and every configuration was verified
> word-for-word against GMP before it was timed.


## 1. System under test

1 OpenCL device(s) exercised with the identical kernels and operands.

### Device 0 - Intel(R) Graphics (GPU)

| Property | Value |
|---|---|
| Model | Intel(R) Graphics |
| Type | GPU |
| Vendor | Intel(R) Corporation |
| Device memory | 28.01 GiB |
| Max single allocation | 4.00 GiB |
| Local memory | 64 KiB |
| Global cache | 4096 KiB |
| Compute units | 64 |
| Max clock | 2100 MHz |
| Max work-group size | 1024 |
| OpenCL version | OpenCL 3.0 NEO  |
| Driver | 26.09.37435.12 |

### Host

| Property | Value |
|---|---|
| CPU | Intel(R) Core(TM) Ultra 7 265U |
| Logical cores | 14 |
| OpenMP threads used | 14 |
| RAM | 30.3 GB |
| OS | Ubuntu 26.04 LTS |
| Kernel | 7.0.0-30-generic |
| Arch | x86_64 |
| GMP | 6.3.0 |
| OpenSSL | OpenSSL 3.5.5 27 Jan 2026 |
| CGBN | not measured |

## 2. Method

- Workload auto-sized from the device and host: --min-items from 700 x compute units, --items from ten times that capped by host RAM. Either flag, given explicitly, overrides its half.
- Base workload 150000 items, scaled down per operator by its cost weight and by modulus size. Device rows honour --min-items (15000) so the GPU is not left idle; the CPU libraries keep the smaller count because a full-width MODEXP there costs minutes. Both counts appear in every row as dev/cpu, and throughput is per-second so they remain comparable.
- 5 timed repetitions, **minimum** reported. Two untimed warm-up launches precede them.
- `kernel` times `clEnqueueNDRangeKernel` + `clFinish` only. `e2e` adds the host->device operand writes and the device->host result read.
- Every OpenCL device runs the same kernels on the same operands, so GPU and CPU-OpenCL columns are directly comparable.
- CPU library baselines (GMP, OpenSSL) run those same operands, with temporaries preallocated outside the timed region, so the figure is the arithmetic and not marshalling. The generator is reseeded per modulus and operation so every backend sees identical inputs.
- OpenSSL rows time the nearest BN primitive, which is not always semantically identical (its Montgomery routine expects Montgomery-domain inputs); they measure comparable work, not identical results. Correctness is judged against GMP only.
- Every device cell is checked word-for-word against GMP before it is timed. A cell that mismatches is reported and excluded from the speedup tables.
- Total wall time 624.6 s.

## 3. Correctness

| Device | Kernel | Configs run | Passed | Mismatched | Build/launch failed |
|---|---|---|---|---|---|
| [0] GPU | `mpaKernels_8bits.cl` (w8) | 75 | 40 | 0 | 35 |
| [0] GPU | `mpaKernel_16bits.cl` (w16) | 75 | 0 | 0 | 75 |
| [0] GPU | `mpaKernel_32bits.cl` (w32) | 35 | 0 | 0 | 35 |
| [0] GPU | `mpaKernel_32bits_opt.cl` (w32-opt) | 75 | 0 | 0 | 75 |
| [0] GPU | `mpaKernel_32bits_opt.cl` (w32-o64) | 75 | 0 | 0 | 75 |
| [0] GPU | `mpaKernel_32bits_opt.cl` (w32-il) | 75 | 0 | 0 | 75 |
| [0] GPU | `mpaKernel_32bits_opt.cl` (w32-il64) | 75 | 0 | 0 | 75 |

**FAILURES PRESENT** - 485 configurations, 445 problems.

## 4. Throughput per device

Operations per second, higher is better. Kernel-only timings.

### Device 0 - Intel(R) Graphics (GPU)

#### secp256k1 (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 150000 / 150000 | 68.26 M | build failed | build failed | build failed | build failed | build failed | build failed | 80.33 M | n/a |
| SUBTRACT | 150000 / 150000 | 112.73 M | build failed | build failed | build failed | build failed | build failed | build failed | 126.06 M | n/a |
| ADDMOD | 150000 / 150000 | 60.09 M | build failed | build failed | build failed | build failed | build failed | build failed | 37.83 M | n/a |
| SUBTRACTMOD | 150000 / 150000 | 63.42 M | build failed | build failed | build failed | build failed | build failed | build failed | 44.90 M | n/a |
| MULTIPLYOPERANDSCANNING | 150000 / 150000 | 1.73 M | build failed | build failed | build failed | build failed | build failed | build failed | 89.40 M | n/a |
| MULTIPLYPRODUCTSCANNING | 150000 / 150000 | 12.80 M | build failed | build failed | build failed | build failed | build failed | build failed | 91.16 M | n/a |
| MONTGOMERYMULTIPLICATION | 150000 / 150000 | 721.78 k | build failed | build failed | build failed | build failed | build failed | build failed | 10.56 M | n/a |
| COMPARE | 150000 / 150000 | 107.68 M | build failed | - | build failed | build failed | build failed | build failed | 92.41 M | n/a |
| REDUCE | 18750 / 18750 | 447.20 k | build failed | - | build failed | build failed | build failed | build failed | 82.83 M | n/a |
| MODMUL | 15000 / 9375 | 142.17 k | build failed | - | build failed | build failed | build failed | build failed | 20.57 M | n/a |
| MODEXP | 15000 / 2343 | 3.73 k | build failed | - | build failed | build failed | build failed | build failed | 166.34 k | n/a |
| EXPONENTIATION | 15000 / 2343 | 14.39 k | build failed | - | build failed | build failed | build failed | build failed | 494.42 k | n/a |
| DIVIDE | 18750 / 18750 | 235.27 k | build failed | - | build failed | build failed | build failed | build failed | 50.24 M | n/a |
| ISQRT | 15000 / 4687 | 26.40 k | build failed | - | build failed | build failed | build failed | build failed | 13.19 M | n/a |
| MODMUL_R2 | 150000 / 150000 | 869.29 k | build failed | - | build failed | build failed | build failed | build failed | 20.67 M | n/a |

#### rsa256(composite) (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 150000 / 150000 | 92.94 M | build failed | build failed | build failed | build failed | build failed | build failed | 86.15 M | n/a |
| SUBTRACT | 150000 / 150000 | 114.28 M | build failed | build failed | build failed | build failed | build failed | build failed | 131.53 M | n/a |
| ADDMOD | 150000 / 150000 | 79.19 M | build failed | build failed | build failed | build failed | build failed | build failed | 40.99 M | n/a |
| SUBTRACTMOD | 150000 / 150000 | 74.98 M | build failed | build failed | build failed | build failed | build failed | build failed | 44.96 M | n/a |
| MULTIPLYOPERANDSCANNING | 150000 / 150000 | 1.73 M | build failed | build failed | build failed | build failed | build failed | build failed | 89.32 M | n/a |
| MULTIPLYPRODUCTSCANNING | 150000 / 150000 | 12.72 M | build failed | build failed | build failed | build failed | build failed | build failed | 89.41 M | n/a |
| MONTGOMERYMULTIPLICATION | 150000 / 150000 | 724.90 k | build failed | build failed | build failed | build failed | build failed | build failed | 10.64 M | n/a |
| COMPARE | 150000 / 150000 | 110.47 M | build failed | - | build failed | build failed | build failed | build failed | 92.46 M | n/a |
| REDUCE | 18750 / 18750 | 452.11 k | build failed | - | build failed | build failed | build failed | build failed | 63.41 M | n/a |
| MODMUL | 15000 / 9375 | 142.21 k | build failed | - | build failed | build failed | build failed | build failed | 21.03 M | n/a |
| MODEXP | 15000 / 2343 | 3.75 k | build failed | - | build failed | build failed | build failed | build failed | 183.04 k | n/a |
| EXPONENTIATION | 15000 / 2343 | 14.37 k | build failed | - | build failed | build failed | build failed | build failed | 500.45 k | n/a |
| DIVIDE | 18750 / 18750 | 235.88 k | build failed | - | build failed | build failed | build failed | build failed | 49.50 M | n/a |
| ISQRT | 15000 / 4687 | 24.25 k | build failed | - | build failed | build failed | build failed | build failed | 13.96 M | n/a |
| MODMUL_R2 | 150000 / 150000 | 868.90 k | build failed | - | build failed | build failed | build failed | build failed | 20.80 M | n/a |

#### brainpoolP512r1 (512-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 75000 / 75000 | 35.57 M | build failed | build failed | build failed | build failed | build failed | build failed | 74.49 M | n/a |
| SUBTRACT | 75000 / 75000 | 34.24 M | build failed | build failed | build failed | build failed | build failed | build failed | 91.73 M | n/a |
| ADDMOD | 75000 / 75000 | 26.88 M | build failed | build failed | build failed | build failed | build failed | build failed | 41.27 M | n/a |
| SUBTRACTMOD | 75000 / 75000 | 24.07 M | build failed | build failed | build failed | build failed | build failed | build failed | 42.55 M | n/a |
| MULTIPLYOPERANDSCANNING | 75000 / 75000 | 348.71 k | build failed | build failed | build failed | build failed | build failed | build failed | 38.51 M | n/a |
| MULTIPLYPRODUCTSCANNING | 75000 / 75000 | 1.23 M | build failed | build failed | build failed | build failed | build failed | build failed | 38.91 M | n/a |
| MONTGOMERYMULTIPLICATION | 75000 / 75000 | 107.88 k | build failed | build failed | build failed | build failed | build failed | build failed | 4.44 M | n/a |
| COMPARE | 75000 / 75000 | 33.34 M | build failed | - | build failed | build failed | build failed | build failed | 75.27 M | n/a |
| REDUCE | 15000 / 9375 | 66.88 k | build failed | - | build failed | build failed | build failed | build failed | 57.28 M | n/a |
| MODMUL | 15000 / 4687 | 22.93 k | build failed | - | build failed | build failed | build failed | build failed | 2.98 M | n/a |
| MODEXP | 15000 / 1171 | build failed | build failed | - | build failed | build failed | build failed | build failed | 33.15 k | n/a |
| EXPONENTIATION | 15000 / 1171 | build failed | build failed | - | build failed | build failed | build failed | build failed | 151.37 k | n/a |
| DIVIDE | 15000 / 9375 | build failed | build failed | - | build failed | build failed | build failed | build failed | 47.24 M | n/a |
| ISQRT | 15000 / 2343 | build failed | build failed | - | build failed | build failed | build failed | build failed | 12.41 M | n/a |
| MODMUL_R2 | 75000 / 75000 | build failed | build failed | - | build failed | build failed | build failed | build failed | 9.44 M | n/a |

#### p1024 (1024-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 37500 / 37500 | build failed | build failed | build failed | build failed | build failed | build failed | build failed | 54.98 M | n/a |
| SUBTRACT | 37500 / 37500 | build failed | build failed | build failed | build failed | build failed | build failed | build failed | 49.17 M | n/a |
| ADDMOD | 37500 / 37500 | build failed | build failed | build failed | build failed | build failed | build failed | build failed | 28.82 M | n/a |
| SUBTRACTMOD | 37500 / 37500 | build failed | build failed | build failed | build failed | build failed | build failed | build failed | 34.12 M | n/a |
| MULTIPLYOPERANDSCANNING | 37500 / 37500 | build failed | build failed | build failed | build failed | build failed | build failed | build failed | 10.76 M | n/a |
| MULTIPLYPRODUCTSCANNING | 37500 / 37500 | build failed | build failed | build failed | build failed | build failed | build failed | build failed | 10.85 M | n/a |
| MONTGOMERYMULTIPLICATION | 37500 / 37500 | build failed | build failed | build failed | build failed | build failed | build failed | build failed | 1.44 M | n/a |
| COMPARE | 37500 / 37500 | build failed | build failed | - | build failed | build failed | build failed | build failed | 57.65 M | n/a |
| REDUCE | 15000 / 4687 | build failed | build failed | - | build failed | build failed | build failed | build failed | 73.38 M | n/a |
| MODMUL | 15000 / 2343 | build failed | build failed | - | build failed | build failed | build failed | build failed | 3.36 M | n/a |
| MODEXP | 15000 / 585 | build failed | build failed | - | build failed | build failed | build failed | build failed | 5.10 k | n/a |
| EXPONENTIATION | 15000 / 585 | build failed | build failed | - | build failed | build failed | build failed | build failed | 39.87 k | n/a |
| DIVIDE | 15000 / 4687 | build failed | build failed | - | build failed | build failed | build failed | build failed | 40.67 M | n/a |
| ISQRT | 15000 / 1171 | build failed | build failed | - | build failed | build failed | build failed | build failed | 5.75 M | n/a |
| MODMUL_R2 | 37500 / 37500 | build failed | build failed | - | build failed | build failed | build failed | build failed | 3.25 M | n/a |

#### p2048 (2048-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 18750 / 18750 | build failed | build failed | build failed | build failed | build failed | build failed | build failed | 26.58 M | n/a |
| SUBTRACT | 18750 / 18750 | build failed | build failed | build failed | build failed | build failed | build failed | build failed | 25.77 M | n/a |
| ADDMOD | 18750 / 18750 | build failed | build failed | build failed | build failed | build failed | build failed | build failed | 16.76 M | n/a |
| SUBTRACTMOD | 18750 / 18750 | build failed | build failed | build failed | build failed | build failed | build failed | build failed | 17.81 M | n/a |
| MULTIPLYOPERANDSCANNING | 18750 / 18750 | build failed | build failed | build failed | build failed | build failed | build failed | build failed | 3.04 M | n/a |
| MULTIPLYPRODUCTSCANNING | 18750 / 18750 | build failed | build failed | build failed | build failed | build failed | build failed | build failed | 3.11 M | n/a |
| MONTGOMERYMULTIPLICATION | 18750 / 18750 | build failed | build failed | build failed | build failed | build failed | build failed | build failed | 416.52 k | n/a |
| COMPARE | 18750 / 18750 | build failed | build failed | - | build failed | build failed | build failed | build failed | 27.32 M | n/a |
| REDUCE | 15000 / 2343 | build failed | build failed | - | build failed | build failed | build failed | build failed | 47.07 M | n/a |
| MODMUL | 15000 / 1171 | build failed | build failed | - | build failed | build failed | build failed | build failed | 956.86 k | n/a |
| MODEXP | 15000 / 292 | build failed | build failed | - | build failed | build failed | build failed | build failed | 708.3 | n/a |
| EXPONENTIATION | 15000 / 292 | build failed | build failed | - | build failed | build failed | build failed | build failed | 5.85 k | n/a |
| DIVIDE | 15000 / 2343 | build failed | build failed | - | build failed | build failed | build failed | build failed | 34.28 M | n/a |
| ISQRT | 15000 / 585 | build failed | build failed | - | build failed | build failed | build failed | build failed | 1.70 M | n/a |
| MODMUL_R2 | 18750 / 18750 | build failed | build failed | - | build failed | build failed | build failed | build failed | 941.19 k | n/a |

## 5. Head to head

Best OpenCL GPU result against best OpenCL CPU result and the CPU libraries.
Ratios above 1.00x mean the GPU is faster than that baseline.

### secp256k1 (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w8 | 68.26 M | none | n/a | 80.33 M | n/a | n/a | n/a |
| SUBTRACT | w8 | 112.73 M | none | n/a | 126.06 M | n/a | n/a | n/a |
| ADDMOD | w8 | 60.09 M | none | n/a | 37.83 M | n/a | n/a | n/a |
| SUBTRACTMOD | w8 | 63.42 M | none | n/a | 44.90 M | n/a | n/a | n/a |
| MULTIPLYOPERANDSCANNING | w8 | 1.73 M | none | n/a | 89.40 M | n/a | n/a | n/a |
| MULTIPLYPRODUCTSCANNING | w8 | 12.80 M | none | n/a | 91.16 M | n/a | n/a | n/a |
| MONTGOMERYMULTIPLICATION | w8 | 721.78 k | none | n/a | 10.56 M | n/a | n/a | n/a |
| COMPARE | w8 | 107.68 M | none | n/a | 92.41 M | n/a | n/a | n/a |
| REDUCE | w8 | 447.20 k | none | n/a | 82.83 M | n/a | n/a | n/a |
| MODMUL | w8 | 88.86 k | none | n/a | 20.57 M | n/a | n/a | n/a |
| MODEXP | w8 | 582.8 | none | n/a | 166.34 k | n/a | n/a | n/a |
| EXPONENTIATION | w8 | 2.25 k | none | n/a | 494.42 k | n/a | n/a | n/a |
| DIVIDE | w8 | 235.27 k | none | n/a | 50.24 M | n/a | n/a | n/a |
| ISQRT | w8 | 8.25 k | none | n/a | 13.19 M | n/a | n/a | n/a |
| MODMUL_R2 | w8 | 869.29 k | none | n/a | 20.67 M | n/a | n/a | n/a |

### rsa256(composite) (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w8 | 92.94 M | none | n/a | 86.15 M | n/a | n/a | n/a |
| SUBTRACT | w8 | 114.28 M | none | n/a | 131.53 M | n/a | n/a | n/a |
| ADDMOD | w8 | 79.19 M | none | n/a | 40.99 M | n/a | n/a | n/a |
| SUBTRACTMOD | w8 | 74.98 M | none | n/a | 44.96 M | n/a | n/a | n/a |
| MULTIPLYOPERANDSCANNING | w8 | 1.73 M | none | n/a | 89.32 M | n/a | n/a | n/a |
| MULTIPLYPRODUCTSCANNING | w8 | 12.72 M | none | n/a | 89.41 M | n/a | n/a | n/a |
| MONTGOMERYMULTIPLICATION | w8 | 724.90 k | none | n/a | 10.64 M | n/a | n/a | n/a |
| COMPARE | w8 | 110.47 M | none | n/a | 92.46 M | n/a | n/a | n/a |
| REDUCE | w8 | 452.11 k | none | n/a | 63.41 M | n/a | n/a | n/a |
| MODMUL | w8 | 88.88 k | none | n/a | 21.03 M | n/a | n/a | n/a |
| MODEXP | w8 | 586.4 | none | n/a | 183.04 k | n/a | n/a | n/a |
| EXPONENTIATION | w8 | 2.25 k | none | n/a | 500.45 k | n/a | n/a | n/a |
| DIVIDE | w8 | 235.88 k | none | n/a | 49.50 M | n/a | n/a | n/a |
| ISQRT | w8 | 7.58 k | none | n/a | 13.96 M | n/a | n/a | n/a |
| MODMUL_R2 | w8 | 868.90 k | none | n/a | 20.80 M | n/a | n/a | n/a |

### brainpoolP512r1 (512-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w8 | 35.57 M | none | n/a | 74.49 M | n/a | n/a | n/a |
| SUBTRACT | w8 | 34.24 M | none | n/a | 91.73 M | n/a | n/a | n/a |
| ADDMOD | w8 | 26.88 M | none | n/a | 41.27 M | n/a | n/a | n/a |
| SUBTRACTMOD | w8 | 24.07 M | none | n/a | 42.55 M | n/a | n/a | n/a |
| MULTIPLYOPERANDSCANNING | w8 | 348.71 k | none | n/a | 38.51 M | n/a | n/a | n/a |
| MULTIPLYPRODUCTSCANNING | w8 | 1.23 M | none | n/a | 38.91 M | n/a | n/a | n/a |
| MONTGOMERYMULTIPLICATION | w8 | 107.88 k | none | n/a | 4.44 M | n/a | n/a | n/a |
| COMPARE | w8 | 33.34 M | none | n/a | 75.27 M | n/a | n/a | n/a |
| REDUCE | w8 | 41.80 k | none | n/a | 57.28 M | n/a | n/a | n/a |
| MODMUL | w8 | 7.16 k | none | n/a | 2.98 M | n/a | n/a | n/a |
| MODEXP | none | n/a | none | n/a | 33.15 k | n/a | n/a | n/a |
| EXPONENTIATION | none | n/a | none | n/a | 151.37 k | n/a | n/a | n/a |
| DIVIDE | none | n/a | none | n/a | 47.24 M | n/a | n/a | n/a |
| ISQRT | none | n/a | none | n/a | 12.41 M | n/a | n/a | n/a |
| MODMUL_R2 | none | n/a | none | n/a | 9.44 M | n/a | n/a | n/a |

### p1024 (1024-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | none | n/a | none | n/a | 54.98 M | n/a | n/a | n/a |
| SUBTRACT | none | n/a | none | n/a | 49.17 M | n/a | n/a | n/a |
| ADDMOD | none | n/a | none | n/a | 28.82 M | n/a | n/a | n/a |
| SUBTRACTMOD | none | n/a | none | n/a | 34.12 M | n/a | n/a | n/a |
| MULTIPLYOPERANDSCANNING | none | n/a | none | n/a | 10.76 M | n/a | n/a | n/a |
| MULTIPLYPRODUCTSCANNING | none | n/a | none | n/a | 10.85 M | n/a | n/a | n/a |
| MONTGOMERYMULTIPLICATION | none | n/a | none | n/a | 1.44 M | n/a | n/a | n/a |
| COMPARE | none | n/a | none | n/a | 57.65 M | n/a | n/a | n/a |
| REDUCE | none | n/a | none | n/a | 73.38 M | n/a | n/a | n/a |
| MODMUL | none | n/a | none | n/a | 3.36 M | n/a | n/a | n/a |
| MODEXP | none | n/a | none | n/a | 5.10 k | n/a | n/a | n/a |
| EXPONENTIATION | none | n/a | none | n/a | 39.87 k | n/a | n/a | n/a |
| DIVIDE | none | n/a | none | n/a | 40.67 M | n/a | n/a | n/a |
| ISQRT | none | n/a | none | n/a | 5.75 M | n/a | n/a | n/a |
| MODMUL_R2 | none | n/a | none | n/a | 3.25 M | n/a | n/a | n/a |

### p2048 (2048-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | none | n/a | none | n/a | 26.58 M | n/a | n/a | n/a |
| SUBTRACT | none | n/a | none | n/a | 25.77 M | n/a | n/a | n/a |
| ADDMOD | none | n/a | none | n/a | 16.76 M | n/a | n/a | n/a |
| SUBTRACTMOD | none | n/a | none | n/a | 17.81 M | n/a | n/a | n/a |
| MULTIPLYOPERANDSCANNING | none | n/a | none | n/a | 3.04 M | n/a | n/a | n/a |
| MULTIPLYPRODUCTSCANNING | none | n/a | none | n/a | 3.11 M | n/a | n/a | n/a |
| MONTGOMERYMULTIPLICATION | none | n/a | none | n/a | 416.52 k | n/a | n/a | n/a |
| COMPARE | none | n/a | none | n/a | 27.32 M | n/a | n/a | n/a |
| REDUCE | none | n/a | none | n/a | 47.07 M | n/a | n/a | n/a |
| MODMUL | none | n/a | none | n/a | 956.86 k | n/a | n/a | n/a |
| MODEXP | none | n/a | none | n/a | 708.3 | n/a | n/a | n/a |
| EXPONENTIATION | none | n/a | none | n/a | 5.85 k | n/a | n/a | n/a |
| DIVIDE | none | n/a | none | n/a | 34.28 M | n/a | n/a | n/a |
| ISQRT | none | n/a | none | n/a | 1.70 M | n/a | n/a | n/a |
| MODMUL_R2 | none | n/a | none | n/a | 941.19 k | n/a | n/a | n/a |

## 6. CGBN

Not measured on this run. CGBN is CUDA-only and is not built into this host.
To populate the CGBN columns, produce `cgbn_results.tsv` next to the binary with one
whitespace-separated row per measurement and re-run:

```
# modulus_name  operation_name  items  seconds
secp256k1  MODMUL  20000  0.00123
```

`modulus_name` and `operation_name` must match the spellings used in the tables above.

## 7. Raw data

Also written to `Intel_R_Graphics_Report.csv` for analysis.

```csv
kind,device,device_type,kernel,modulus,bits,operation,items,seconds,ops_per_sec,mismatches
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,secp256k1,256,ADD,150000,0.001867299,80329931.212,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,secp256k1,256,ADD,150000,0.001897254,79061633.357,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,secp256k1,256,ADD,150000,0.001548597,96861868.721,0
opencl-kernel,Intel(R) Graphics,GPU,w8,secp256k1,256,ADD,150000,0.002197517,68258857.773,0
opencl-e2e,Intel(R) Graphics,GPU,w8,secp256k1,256,ADD,150000,0.004895016,30643413.639,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,secp256k1,256,SUBTRACT,150000,0.001189908,126060164.646,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,secp256k1,256,SUBTRACT,150000,0.001333155,112515048.898,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,secp256k1,256,SUBTRACT,150000,0.001526421,98269088.653,0
opencl-kernel,Intel(R) Graphics,GPU,w8,secp256k1,256,SUBTRACT,150000,0.001330605,112730674.333,0
opencl-e2e,Intel(R) Graphics,GPU,w8,secp256k1,256,SUBTRACT,150000,0.003523177,42575209.798,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,secp256k1,256,ADDMOD,150000,0.003965175,37829351.678,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,secp256k1,256,ADDMOD,150000,0.001895131,79150201.631,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,secp256k1,256,ADDMOD,150000,0.003323443,45133916.829,0
opencl-kernel,Intel(R) Graphics,GPU,w8,secp256k1,256,ADDMOD,150000,0.002496101,60093722.379,0
opencl-e2e,Intel(R) Graphics,GPU,w8,secp256k1,256,ADDMOD,150000,0.005232690,28665944.342,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,secp256k1,256,SUBTRACTMOD,150000,0.003340540,44902920.047,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,secp256k1,256,SUBTRACTMOD,150000,0.003120841,48063967.484,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,secp256k1,256,SUBTRACTMOD,150000,0.004252837,35270573.415,0
opencl-kernel,Intel(R) Graphics,GPU,w8,secp256k1,256,SUBTRACTMOD,150000,0.002365180,63420120.239,0
opencl-e2e,Intel(R) Graphics,GPU,w8,secp256k1,256,SUBTRACTMOD,150000,0.004637359,32345996.926,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,secp256k1,256,MULTIPLYOPERANDSCANNING,150000,0.001677853,89399965.468,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,150000,0.001468410,102151307.359,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,150000,0.001865326,80414897.146,0
opencl-kernel,Intel(R) Graphics,GPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,150000,0.086669299,1730716.664,0
opencl-e2e,Intel(R) Graphics,GPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,150000,0.091744623,1634973.202,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,secp256k1,256,MULTIPLYPRODUCTSCANNING,150000,0.001645496,91157924.493,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,150000,0.001111905,134903610.982,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,150000,0.002027371,73987444.967,0
opencl-kernel,Intel(R) Graphics,GPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,150000,0.011720365,12798236.232,0
opencl-e2e,Intel(R) Graphics,GPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,150000,0.016280518,9213466.058,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,secp256k1,256,MONTGOMERYMULTIPLICATION,150000,0.014210056,10555904.925,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,150000,0.005351670,28028634.057,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,150000,0.001807040,83008676.666,0
opencl-kernel,Intel(R) Graphics,GPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,150000,0.207819815,721779.105,0
opencl-e2e,Intel(R) Graphics,GPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,150000,0.212325985,706460.870,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,secp256k1,256,COMPARE,150000,0.001623210,92409484.771,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,secp256k1,256,COMPARE,150000,0.000808432,185544362.476,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,secp256k1,256,COMPARE,150000,0.001494298,100381584.494,0
opencl-kernel,Intel(R) Graphics,GPU,w8,secp256k1,256,COMPARE,150000,0.001393018,107679872.097,0
opencl-e2e,Intel(R) Graphics,GPU,w8,secp256k1,256,COMPARE,150000,0.003373699,44461583.551,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,secp256k1,256,REDUCE,18750,0.000226362,82831925.113,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,secp256k1,256,REDUCE,18750,0.000296230,63295414.194,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,secp256k1,256,REDUCE,18750,0.001730984,10831989.269,0
opencl-kernel,Intel(R) Graphics,GPU,w8,secp256k1,256,REDUCE,18750,0.041927481,447200.727,0
opencl-e2e,Intel(R) Graphics,GPU,w8,secp256k1,256,REDUCE,18750,0.042835182,437724.299,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,secp256k1,256,MODMUL,9375,0.000455734,20571210.406,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,secp256k1,256,MODMUL,9375,0.000560934,16713196.623,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,secp256k1,256,MODMUL,9375,0.001025308,9143593.828,0
opencl-kernel,Intel(R) Graphics,GPU,w8,secp256k1,256,MODMUL,15000,0.105503903,142174.835,0
opencl-e2e,Intel(R) Graphics,GPU,w8,secp256k1,256,MODMUL,15000,0.108323056,138474.675,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,secp256k1,256,MODEXP,2343,0.014085848,166337.163,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,secp256k1,256,MODEXP,2343,0.004856354,482460.710,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,secp256k1,256,MODEXP,2343,0.005008969,467760.931,0
opencl-kernel,Intel(R) Graphics,GPU,w8,secp256k1,256,MODEXP,15000,4.019905329,3731.431,0
opencl-e2e,Intel(R) Graphics,GPU,w8,secp256k1,256,MODEXP,15000,4.016454878,3734.637,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,secp256k1,256,EXPONENTIATION,2343,0.004738916,494416.869,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,secp256k1,256,EXPONENTIATION,2343,0.001476327,1587046.755,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,secp256k1,256,EXPONENTIATION,2343,0.009626171,243398.959,0
opencl-kernel,Intel(R) Graphics,GPU,w8,secp256k1,256,EXPONENTIATION,15000,1.042721131,14385.438,0
opencl-e2e,Intel(R) Graphics,GPU,w8,secp256k1,256,EXPONENTIATION,15000,1.040436017,14417.033,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,secp256k1,256,DIVIDE,18750,0.000373193,50242099.308,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,secp256k1,256,DIVIDE,18750,0.000117523,159543259.373,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,secp256k1,256,DIVIDE,18750,0.000460338,40730941.798,0
opencl-kernel,Intel(R) Graphics,GPU,w8,secp256k1,256,DIVIDE,18750,0.079695606,235270.185,0
opencl-e2e,Intel(R) Graphics,GPU,w8,secp256k1,256,DIVIDE,18750,0.081386760,230381.453,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,secp256k1,256,ISQRT,4687,0.000355214,13194862.928,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,secp256k1,256,ISQRT,4687,0.000093146,50318841.537,0
opencl-kernel,Intel(R) Graphics,GPU,w8,secp256k1,256,ISQRT,15000,0.568279892,26395.444,0
opencl-e2e,Intel(R) Graphics,GPU,w8,secp256k1,256,ISQRT,15000,0.568695017,26376.176,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,secp256k1,256,MODMUL_R2,150000,0.007257648,20667852.763,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,secp256k1,256,MODMUL_R2,150000,0.004058327,36961043.289,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,secp256k1,256,MODMUL_R2,150000,0.005549187,27030986.627,0
opencl-kernel,Intel(R) Graphics,GPU,w8,secp256k1,256,MODMUL_R2,150000,0.172555143,869287.333,0
opencl-e2e,Intel(R) Graphics,GPU,w8,secp256k1,256,MODMUL_R2,150000,0.177123263,846867.867,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,rsa256(composite),256,ADD,150000,0.001741165,86149215.978,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,rsa256(composite),256,ADD,150000,0.001849793,81090153.793,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,rsa256(composite),256,ADD,150000,0.001504601,99694204.080,0
opencl-kernel,Intel(R) Graphics,GPU,w8,rsa256(composite),256,ADD,150000,0.001613969,92938587.326,0
opencl-e2e,Intel(R) Graphics,GPU,w8,rsa256(composite),256,ADD,150000,0.003732105,40191795.360,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACT,150000,0.001140394,131533485.792,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACT,150000,0.001189931,126057728.892,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACT,150000,0.001334954,112363422.685,0
opencl-kernel,Intel(R) Graphics,GPU,w8,rsa256(composite),256,SUBTRACT,150000,0.001312546,114281709.379,0
opencl-e2e,Intel(R) Graphics,GPU,w8,rsa256(composite),256,SUBTRACT,150000,0.003295551,45515909.231,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,rsa256(composite),256,ADDMOD,150000,0.003659673,40987268.467,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,rsa256(composite),256,ADDMOD,150000,0.002720976,55127277.544,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,rsa256(composite),256,ADDMOD,150000,0.004086889,36702733.978,0
opencl-kernel,Intel(R) Graphics,GPU,w8,rsa256(composite),256,ADDMOD,150000,0.001894259,79186637.877,0
opencl-e2e,Intel(R) Graphics,GPU,w8,rsa256(composite),256,ADDMOD,150000,0.003891942,38541170.410,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACTMOD,150000,0.003336436,44958152.968,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACTMOD,150000,0.002327145,64456662.841,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACTMOD,150000,0.003953318,37942811.629,0
opencl-kernel,Intel(R) Graphics,GPU,w8,rsa256(composite),256,SUBTRACTMOD,150000,0.002000643,74975895.302,0
opencl-e2e,Intel(R) Graphics,GPU,w8,rsa256(composite),256,SUBTRACTMOD,150000,0.004176840,35912316.389,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYOPERANDSCANNING,150000,0.001679377,89318837.216,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,150000,0.001055239,142147892.389,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,150000,0.002292915,65418909.575,0
opencl-kernel,Intel(R) Graphics,GPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,150000,0.086737682,1729352.186,0
opencl-e2e,Intel(R) Graphics,GPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,150000,0.091664501,1636402.297,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,150000,0.001677633,89411688.903,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,150000,0.001522440,98526050.440,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,150000,0.002259101,66398094.018,0
opencl-kernel,Intel(R) Graphics,GPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,150000,0.011796345,12715803.101,0
opencl-e2e,Intel(R) Graphics,GPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,150000,0.016981294,8833249.103,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,rsa256(composite),256,MONTGOMERYMULTIPLICATION,150000,0.014099888,10638382.369,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,150000,0.004097934,36603810.626,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,150000,0.002662128,56345901.082,0
opencl-kernel,Intel(R) Graphics,GPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,150000,0.206923946,724904.019,0
opencl-e2e,Intel(R) Graphics,GPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,150000,0.211516554,709164.352,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,rsa256(composite),256,COMPARE,150000,0.001622329,92459667.530,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,rsa256(composite),256,COMPARE,150000,0.001011996,148221928.833,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,rsa256(composite),256,COMPARE,150000,0.001166290,128612953.020,0
opencl-kernel,Intel(R) Graphics,GPU,w8,rsa256(composite),256,COMPARE,150000,0.001357816,110471520.759,0
opencl-e2e,Intel(R) Graphics,GPU,w8,rsa256(composite),256,COMPARE,150000,0.003233496,46389418.999,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,rsa256(composite),256,REDUCE,18750,0.000295707,63407360.270,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,rsa256(composite),256,REDUCE,18750,0.000398754,47021472.273,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,rsa256(composite),256,REDUCE,18750,0.000723230,25925362.171,0
opencl-kernel,Intel(R) Graphics,GPU,w8,rsa256(composite),256,REDUCE,18750,0.041471891,452113.457,0
opencl-e2e,Intel(R) Graphics,GPU,w8,rsa256(composite),256,REDUCE,18750,0.042658752,439534.659,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,rsa256(composite),256,MODMUL,9375,0.000445776,21030741.528,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,rsa256(composite),256,MODMUL,9375,0.000337407,27785432.396,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,rsa256(composite),256,MODMUL,9375,0.000913680,10260703.885,0
opencl-kernel,Intel(R) Graphics,GPU,w8,rsa256(composite),256,MODMUL,15000,0.105474955,142213.855,0
opencl-e2e,Intel(R) Graphics,GPU,w8,rsa256(composite),256,MODMUL,15000,0.106293143,141119.169,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,rsa256(composite),256,MODEXP,2343,0.012800728,183036.465,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,rsa256(composite),256,MODEXP,2343,0.005790837,404604.723,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,rsa256(composite),256,MODEXP,2343,0.005402146,433716.527,0
opencl-kernel,Intel(R) Graphics,GPU,w8,rsa256(composite),256,MODEXP,15000,3.995564952,3754.162,0
opencl-e2e,Intel(R) Graphics,GPU,w8,rsa256(composite),256,MODEXP,15000,3.996963127,3752.849,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,rsa256(composite),256,EXPONENTIATION,2343,0.004681747,500454.211,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,rsa256(composite),256,EXPONENTIATION,2343,0.002156587,1086438.898,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,rsa256(composite),256,EXPONENTIATION,2343,0.009815701,238699.203,0
opencl-kernel,Intel(R) Graphics,GPU,w8,rsa256(composite),256,EXPONENTIATION,15000,1.043506629,14374.609,0
opencl-e2e,Intel(R) Graphics,GPU,w8,rsa256(composite),256,EXPONENTIATION,15000,1.041288446,14405.230,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,rsa256(composite),256,DIVIDE,18750,0.000378767,49502728.676,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,rsa256(composite),256,DIVIDE,18750,0.000106460,176122523.269,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,rsa256(composite),256,DIVIDE,18750,0.000257581,72792634.199,0
opencl-kernel,Intel(R) Graphics,GPU,w8,rsa256(composite),256,DIVIDE,18750,0.079489129,235881.312,0
opencl-e2e,Intel(R) Graphics,GPU,w8,rsa256(composite),256,DIVIDE,18750,0.081068389,231286.205,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,rsa256(composite),256,ISQRT,4687,0.000335845,13955841.917,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,rsa256(composite),256,ISQRT,4687,0.000086464,54207520.843,0
opencl-kernel,Intel(R) Graphics,GPU,w8,rsa256(composite),256,ISQRT,15000,0.618439359,24254.601,0
opencl-e2e,Intel(R) Graphics,GPU,w8,rsa256(composite),256,ISQRT,15000,0.616113324,24346.170,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,rsa256(composite),256,MODMUL_R2,150000,0.007211253,20800823.328,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,rsa256(composite),256,MODMUL_R2,150000,0.003065381,48933558.414,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,rsa256(composite),256,MODMUL_R2,150000,0.005313184,28231659.230,0
opencl-kernel,Intel(R) Graphics,GPU,w8,rsa256(composite),256,MODMUL_R2,150000,0.172632049,868900.073,0
opencl-e2e,Intel(R) Graphics,GPU,w8,rsa256(composite),256,MODMUL_R2,150000,0.177097990,846988.721,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,brainpoolP512r1,512,ADD,75000,0.001006837,74490706.819,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,brainpoolP512r1,512,ADD,75000,0.000841604,89115546.798,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,brainpoolP512r1,512,ADD,75000,0.000793611,94504739.597,0
opencl-kernel,Intel(R) Graphics,GPU,w8,brainpoolP512r1,512,ADD,75000,0.002108537,35569686.793,0
opencl-e2e,Intel(R) Graphics,GPU,w8,brainpoolP512r1,512,ADD,75000,0.004528048,16563428.720,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACT,75000,0.000817659,91725280.146,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACT,75000,0.000858990,87311841.975,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACT,75000,0.000838175,89480119.588,0
opencl-kernel,Intel(R) Graphics,GPU,w8,brainpoolP512r1,512,SUBTRACT,75000,0.002190606,34237101.527,0
opencl-e2e,Intel(R) Graphics,GPU,w8,brainpoolP512r1,512,SUBTRACT,75000,0.004449190,16857000.866,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,brainpoolP512r1,512,ADDMOD,75000,0.001817114,41274240.774,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,brainpoolP512r1,512,ADDMOD,75000,0.001803324,41589863.939,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,brainpoolP512r1,512,ADDMOD,75000,0.004028372,18617942.927,0
opencl-kernel,Intel(R) Graphics,GPU,w8,brainpoolP512r1,512,ADDMOD,75000,0.002789791,26883734.543,0
opencl-e2e,Intel(R) Graphics,GPU,w8,brainpoolP512r1,512,ADDMOD,75000,0.005347280,14025822.534,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACTMOD,75000,0.001762710,42548121.518,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACTMOD,75000,0.002060846,36392820.880,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACTMOD,75000,0.003298600,22736918.689,0
opencl-kernel,Intel(R) Graphics,GPU,w8,brainpoolP512r1,512,SUBTRACTMOD,75000,0.003115693,24071691.345,0
opencl-e2e,Intel(R) Graphics,GPU,w8,brainpoolP512r1,512,SUBTRACTMOD,75000,0.006173785,12148139.289,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,75000,0.001947754,38505889.258,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,75000,0.001782379,42078592.790,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,75000,0.001512231,49595597.987,0
opencl-kernel,Intel(R) Graphics,GPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,75000,0.215077540,348711.446,0
opencl-e2e,Intel(R) Graphics,GPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,75000,0.221293950,338915.727,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,75000,0.001927698,38906509.175,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,75000,0.001423195,52698330.000,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,75000,0.001806045,41527204.927,0
opencl-kernel,Intel(R) Graphics,GPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,75000,0.060773193,1234096.751,0
opencl-e2e,Intel(R) Graphics,GPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,75000,0.065108983,1151914.783,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,75000,0.016889583,4440606.971,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,75000,0.004753224,15778764.060,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,75000,0.004978870,15063658.985,0
opencl-kernel,Intel(R) Graphics,GPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,75000,0.695221048,107879.358,0
opencl-e2e,Intel(R) Graphics,GPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,75000,0.701735393,106877.893,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,brainpoolP512r1,512,COMPARE,75000,0.000996364,75273695.946,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,brainpoolP512r1,512,COMPARE,75000,0.000400006,187497194.440,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,brainpoolP512r1,512,COMPARE,75000,0.000688104,108995153.786,0
opencl-kernel,Intel(R) Graphics,GPU,w8,brainpoolP512r1,512,COMPARE,75000,0.002249706,33337689.270,0
opencl-e2e,Intel(R) Graphics,GPU,w8,brainpoolP512r1,512,COMPARE,75000,0.004427077,16941200.757,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,brainpoolP512r1,512,REDUCE,9375,0.000163681,57276046.455,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,brainpoolP512r1,512,REDUCE,9375,0.000169053,55455979.831,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,brainpoolP512r1,512,REDUCE,9375,0.000556210,16855144.708,0
opencl-kernel,Intel(R) Graphics,GPU,w8,brainpoolP512r1,512,REDUCE,15000,0.224274968,66882.185,0
opencl-e2e,Intel(R) Graphics,GPU,w8,brainpoolP512r1,512,REDUCE,15000,0.225625079,66481.971,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL,4687,0.001574566,2976693.249,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL,4687,0.000584090,8024448.520,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL,4687,0.000888521,5275058.058,0
opencl-kernel,Intel(R) Graphics,GPU,w8,brainpoolP512r1,512,MODMUL,15000,0.654180450,22929.453,0
opencl-e2e,Intel(R) Graphics,GPU,w8,brainpoolP512r1,512,MODMUL,15000,0.656083596,22862.940,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,brainpoolP512r1,512,MODEXP,1171,0.035326811,33147.628,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,brainpoolP512r1,512,MODEXP,1171,0.008535536,137191.150,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,brainpoolP512r1,512,MODEXP,1171,0.008634272,135622.320,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,brainpoolP512r1,512,EXPONENTIATION,1171,0.007736259,151365.149,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,brainpoolP512r1,512,EXPONENTIATION,1171,0.003054501,383368.675,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,brainpoolP512r1,512,EXPONENTIATION,1171,0.015081231,77646.181,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,brainpoolP512r1,512,DIVIDE,9375,0.000198455,47239927.601,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,brainpoolP512r1,512,DIVIDE,9375,0.000077436,121067706.801,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,brainpoolP512r1,512,DIVIDE,9375,0.000192266,48760567.284,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,brainpoolP512r1,512,ISQRT,2343,0.000188829,12408051.799,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,brainpoolP512r1,512,ISQRT,2343,0.000082089,28542187.588,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL_R2,75000,0.007941602,9443938.387,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL_R2,75000,0.005976665,12548804.435,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL_R2,75000,0.006003881,12491919.831,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,p1024,1024,ADD,37500,0.000682007,54984773.083,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,p1024,1024,ADD,37500,0.000615329,60943006.140,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,p1024,1024,ADD,37500,0.000738896,50751391.747,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,p1024,1024,SUBTRACT,37500,0.000762699,49167496.134,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,p1024,1024,SUBTRACT,37500,0.000537474,69770816.895,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,p1024,1024,SUBTRACT,37500,0.000723499,51831445.543,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,p1024,1024,ADDMOD,37500,0.001301014,28823671.804,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,p1024,1024,ADDMOD,37500,0.000870732,43067211.995,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,p1024,1024,ADDMOD,37500,0.002011716,18640802.027,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,p1024,1024,SUBTRACTMOD,37500,0.001099063,34119973.222,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,p1024,1024,SUBTRACTMOD,37500,0.000907699,41313253.607,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,p1024,1024,SUBTRACTMOD,37500,0.003327468,11269830.384,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,p1024,1024,MULTIPLYOPERANDSCANNING,37500,0.003486639,10755343.454,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,p1024,1024,MULTIPLYOPERANDSCANNING,37500,0.002833590,13234095.270,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,p1024,1024,MULTIPLYOPERANDSCANNING,37500,0.004072257,9208652.603,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,p1024,1024,MULTIPLYPRODUCTSCANNING,37500,0.003457805,10845030.313,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,37500,0.002537506,14778290.225,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,37500,0.002592393,14465399.380,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,p1024,1024,MONTGOMERYMULTIPLICATION,37500,0.026066808,1438611.126,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,p1024,1024,MONTGOMERYMULTIPLICATION,37500,0.007024808,5338224.190,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,p1024,1024,MONTGOMERYMULTIPLICATION,37500,0.004324610,8671302.160,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,p1024,1024,COMPARE,37500,0.000650509,57647165.646,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,p1024,1024,COMPARE,37500,0.000322171,116397805.905,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,p1024,1024,COMPARE,37500,0.000533496,70291061.151,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,p1024,1024,REDUCE,4687,0.000063870,73383451.473,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,p1024,1024,REDUCE,4687,0.000128616,36441812.190,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,p1024,1024,REDUCE,4687,0.000376498,12448937.096,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,p1024,1024,MODMUL,2343,0.000697398,3359631.201,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,p1024,1024,MODMUL,2343,0.000646273,3625402.919,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,p1024,1024,MODMUL,2343,0.001251675,1871891.688,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,p1024,1024,MODEXP,585,0.114802133,5095.724,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,p1024,1024,MODEXP,585,0.029700258,19696.799,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,p1024,1024,MODEXP,585,0.030387333,19251.443,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,p1024,1024,EXPONENTIATION,585,0.014671219,39873.987,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,p1024,1024,EXPONENTIATION,585,0.005955679,98225.576,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,p1024,1024,EXPONENTIATION,585,0.045178993,12948.496,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,p1024,1024,DIVIDE,4687,0.000115245,40669877.713,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,p1024,1024,DIVIDE,4687,0.000056170,83443140.899,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,p1024,1024,DIVIDE,4687,0.000168210,27863982.398,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,p1024,1024,ISQRT,1171,0.000203666,5749610.045,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,p1024,1024,ISQRT,1171,0.000090592,12926087.925,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,p1024,1024,MODMUL_R2,37500,0.011523877,3254113.183,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,p1024,1024,MODMUL_R2,37500,0.004879367,7685423.108,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,p1024,1024,MODMUL_R2,37500,0.006786143,5525966.665,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,p2048,2048,ADD,18750,0.000705463,26578289.838,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,p2048,2048,ADD,18750,0.000696984,26901621.718,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,p2048,2048,ADD,18750,0.000581505,32243918.248,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,p2048,2048,SUBTRACT,18750,0.000727703,25766005.867,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,p2048,2048,SUBTRACT,18750,0.000466610,40183452.151,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,p2048,2048,SUBTRACT,18750,0.000724086,25894714.701,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,p2048,2048,ADDMOD,18750,0.001118619,16761739.327,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,p2048,2048,ADDMOD,18750,0.000918897,20404898.163,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,p2048,2048,ADDMOD,18750,0.002708692,6922160.249,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,p2048,2048,SUBTRACTMOD,18750,0.001052977,17806656.858,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,p2048,2048,SUBTRACTMOD,18750,0.000858868,21831061.182,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,p2048,2048,SUBTRACTMOD,18750,0.002530172,7410563.388,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,p2048,2048,MULTIPLYOPERANDSCANNING,18750,0.006174844,3036513.949,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,p2048,2048,MULTIPLYOPERANDSCANNING,18750,0.002870818,6531239.554,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,p2048,2048,MULTIPLYOPERANDSCANNING,18750,0.003478052,5390948.715,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,p2048,2048,MULTIPLYPRODUCTSCANNING,18750,0.006020020,3114607.595,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,18750,0.003324035,5640734.826,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,18750,0.004854690,3862244.534,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,p2048,2048,MONTGOMERYMULTIPLICATION,18750,0.045015957,416518.969,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,p2048,2048,MONTGOMERYMULTIPLICATION,18750,0.009726222,1927778.325,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,p2048,2048,MONTGOMERYMULTIPLICATION,18750,0.004760238,3938878.692,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,p2048,2048,COMPARE,18750,0.000686233,27323081.621,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,p2048,2048,COMPARE,18750,0.000170509,109964862.933,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,p2048,2048,COMPARE,18750,0.000362772,51685356.121,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,p2048,2048,REDUCE,2343,0.000049773,47073721.760,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,p2048,2048,REDUCE,2343,0.000143991,16271847.322,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,p2048,2048,REDUCE,2343,0.000285866,8196147.819,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,p2048,2048,MODMUL,1171,0.001223798,956857.261,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,p2048,2048,MODMUL,1171,0.000680526,1720727.813,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,p2048,2048,MODMUL,1171,0.000715797,1635938.657,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,p2048,2048,MODEXP,292,0.412227085,708.347,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,p2048,2048,MODEXP,292,0.089889729,3248.425,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,p2048,2048,MODEXP,292,0.087966610,3319.441,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,p2048,2048,EXPONENTIATION,292,0.049928353,5848.380,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,p2048,2048,EXPONENTIATION,292,0.011347907,25731.617,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,p2048,2048,EXPONENTIATION,292,0.120286387,2427.540,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,p2048,2048,DIVIDE,2343,0.000068355,34276935.357,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,p2048,2048,DIVIDE,2343,0.000075438,31058616.990,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,p2048,2048,DIVIDE,2343,0.000186149,12586691.183,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,p2048,2048,ISQRT,585,0.000344627,1697487.365,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,p2048,2048,ISQRT,585,0.000133859,4370269.467,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-1t,p2048,2048,MODMUL_R2,18750,0.019921510,941193.714,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,gmp-nt,p2048,2048,MODMUL_R2,18750,0.005244817,3574957.904,0
library,Intel(R) Core(TM) Ultra 7 265U,host-cpu,openssl-nt,p2048,2048,MODMUL_R2,18750,0.007777290,2410865.478,0
```
