# MPA-OpenCL benchmark report - NVIDIA GeForce RTX 3080 Ti

> **Note.** The multi-threaded GMP and OpenSSL baseline columns have been
> removed from this report: they predate the 2026-09-12 timing fix and were
> understated (see `reports/README.md`). The single-threaded GMP column, the
> OpenCL-on-CPU rows and all MPA measurements are unaffected and were verified
> against GMP before timing.


> **Partial report.** The run was interrupted or hit its time budget.
> Rows that never ran are marked `n/a`.

## 1. System under test

2 OpenCL device(s) exercised with the identical kernels and operands.

### Device 0 - NVIDIA GeForce RTX 3080 Ti (GPU)

| Property | Value |
|---|---|
| Model | NVIDIA GeForce RTX 3080 Ti |
| Type | GPU |
| Vendor | NVIDIA Corporation |
| Device memory | 11.63 GiB |
| Max single allocation | 2.91 GiB |
| Local memory | 48 KiB |
| Global cache | 2240 KiB |
| Compute units | 80 |
| Max clock | 1755 MHz |
| Max work-group size | 1024 |
| OpenCL version | OpenCL 3.0 CUDA |
| Driver | 595.58.03 |

### Device 1 - cpu-haswell-AMD Ryzen 9 3900X 12-Core Processor (CPU)

| Property | Value |
|---|---|
| Model | cpu-haswell-AMD Ryzen 9 3900X 12-Core Processor |
| Type | CPU |
| Vendor | AuthenticAMD |
| Device memory | 123.70 GiB |
| Max single allocation | 32.00 GiB |
| Local memory | 512 KiB |
| Global cache | 16384 KiB |
| Compute units | 24 |
| Max clock | 4673 MHz |
| Max work-group size | 4096 |
| OpenCL version | OpenCL 3.0 PoCL HSTR: cpu-x86_64-pc-linux-gnu-haswell |
| Driver | 5.0+debian |

### Host

| Property | Value |
|---|---|
| CPU | AMD Ryzen 9 3900X 12-Core Processor |
| Logical cores | 24 |
| OpenMP threads used | 24 |
| RAM | 125.7 GB |
| OS | Ubuntu 24.04.4 LTS |
| Kernel | 6.8.0-111-generic |
| Arch | x86_64 |
| GMP | 6.3.0 |
| OpenSSL | OpenSSL 3.0.13 30 Jan 2024 |
| CGBN | not measured |

## 2. Method

- Base workload 50000 items, scaled down per operator by its cost weight and by modulus size; the exact count is in every row.
- 5 timed repetitions, **minimum** reported. Two untimed warm-up launches precede them.
- `kernel` times `clEnqueueNDRangeKernel` + `clFinish` only. `e2e` adds the host->device operand writes and the device->host result read.
- Every OpenCL device runs the same kernels on the same operands, so GPU and CPU-OpenCL columns are directly comparable.
- CPU library baselines (GMP, OpenSSL) run those same operands, with temporaries preallocated outside the timed region, so the figure is the arithmetic and not marshalling. The generator is reseeded per modulus and operation so every backend sees identical inputs.
- OpenSSL rows time the nearest BN primitive, which is not always semantically identical (its Montgomery routine expects Montgomery-domain inputs); they measure comparable work, not identical results. Correctness is judged against GMP only.
- Every device cell is checked word-for-word against GMP before it is timed. A cell that mismatches is reported and excluded from the speedup tables.
- Total wall time 452.1 s.

## 3. Correctness

| Device | Kernel | Configs run | Passed | Mismatched | Build/launch failed |
|---|---|---|---|---|---|
| [0] GPU | `mpaKernels_8bits.cl` (w8) | 75 | 75 | 0 | 0 |
| [0] GPU | `mpaKernel_16bits.cl` (w16) | 75 | 75 | 0 | 0 |
| [0] GPU | `mpaKernel_32bits.cl` (w32) | 35 | 35 | 0 | 0 |
| [0] GPU | `mpaKernel_32bits_opt.cl` (w32-opt) | 75 | 75 | 0 | 0 |
| [0] GPU | `mpaKernel_32bits_opt.cl` (w32-o64) | 75 | 75 | 0 | 0 |

**All configurations correct** - 335 configurations, 0 problems.

## 4. Throughput per device

Operations per second, higher is better. Kernel-only timings.

### Device 0 - NVIDIA GeForce RTX 3080 Ti (GPU)

#### secp256k1 (256-bit)

| Operation | items | w8 | w16 | w32 | w32-opt | w32-o64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | 50000 | 1.11 G | 2.06 G | 3.27 G | 3.09 G | 3.33 G | 73.17 M | n/a |
| SUBTRACT | 50000 | 1.11 G | 2.09 G | 3.42 G | 3.07 G | 3.36 G | 94.83 M | n/a |
| ADDMOD | 50000 | 757.19 M | 1.51 G | 3.06 G | 4.08 G | 4.31 G | 22.44 M | n/a |
| SUBTRACTMOD | 50000 | 757.07 M | 1.48 G | 2.95 G | 4.27 G | 4.13 G | 25.50 M | n/a |
| MULTIPLYOPERANDSCANNING | 50000 | 32.33 M | 127.80 M | 503.54 M | 1.77 G | 2.09 G | 58.35 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 | 181.58 M | 599.04 M | 1.77 G | 1.58 G | 1.76 G | 58.70 M | n/a |
| MONTGOMERYMULTIPLICATION | 50000 | 364.46 M | 1.42 G | 3.77 G | 2.72 G | 2.89 G | 7.07 M | n/a |
| COMPARE | 50000 | 1.14 G | 2.15 G | - | 4.24 G | 3.90 G | 164.83 M | n/a |
| REDUCE | 6250 | 78.13 M | 123.31 M | - | 215.78 M | 243.11 M | 59.28 M | n/a |
| MODMUL | 3125 | 21.21 M | 36.18 M | - | 44.27 M | 53.45 M | 12.06 M | n/a |
| MODEXP | 781 | 100.50 k | 597.33 k | - | 859.89 k | 1.50 M | 122.48 k | n/a |
| EXPONENTIATION | 781 | 79.44 k | 301.41 k | - | 4.50 M | 7.92 M | 366.83 k | n/a |
| DIVIDE | 6250 | 20.55 M | 28.61 M | - | 63.58 M | 73.85 M | 30.25 M | n/a |
| ISQRT | 1562 | 864.71 k | 1.05 M | - | 2.90 M | 3.34 M | 12.25 M | n/a |
| MODMUL_R2 | 50000 | 365.67 M | 1.34 G | - | 2.04 G | 2.06 G | 11.97 M | n/a |

#### rsa256(composite) (256-bit)

| Operation | items | w8 | w16 | w32 | w32-opt | w32-o64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | 50000 | 1.23 G | 2.07 G | 3.36 G | 3.09 G | 3.38 G | 71.92 M | n/a |
| SUBTRACT | 50000 | 1.22 G | 2.05 G | 3.38 G | 3.07 G | 3.42 G | 94.74 M | n/a |
| ADDMOD | 50000 | 911.87 M | 1.62 G | 3.12 G | 4.07 G | 4.37 G | 25.81 M | n/a |
| SUBTRACTMOD | 50000 | 842.29 M | 1.48 G | 2.98 G | 4.26 G | 4.09 G | 25.72 M | n/a |
| MULTIPLYOPERANDSCANNING | 50000 | 36.46 M | 125.81 M | 505.43 M | 1.79 G | 2.04 G | 58.43 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 | 204.15 M | 596.25 M | 1.77 G | 1.59 G | 1.79 G | 58.24 M | n/a |
| MONTGOMERYMULTIPLICATION | 50000 | 410.04 M | 1.34 G | 3.67 G | 2.70 G | 2.84 G | 7.00 M | n/a |
| COMPARE | 50000 | 1.27 G | 2.10 G | - | 4.27 G | 3.87 G | 143.36 M | n/a |
| REDUCE | 6250 | 86.63 M | 122.27 M | - | 213.93 M | 239.20 M | 34.52 M | n/a |
| MODMUL | 3125 | 24.42 M | 36.97 M | - | 45.63 M | 55.33 M | 12.24 M | n/a |
| MODEXP | 781 | 115.73 k | 633.88 k | - | 900.38 k | 1.61 M | 129.48 k | n/a |
| EXPONENTIATION | 781 | 78.85 k | 301.22 k | - | 4.55 M | 7.89 M | 370.18 k | n/a |
| DIVIDE | 6250 | 20.15 M | 27.99 M | - | 60.41 M | 70.06 M | 27.61 M | n/a |
| ISQRT | 1562 | 766.62 k | 920.86 k | - | 2.48 M | 2.84 M | 12.19 M | n/a |
| MODMUL_R2 | 50000 | 365.45 M | 1.37 G | - | 2.00 G | 2.09 G | 11.78 M | n/a |

#### brainpoolP512r1 (512-bit)

| Operation | items | w8 | w16 | w32 | w32-opt | w32-o64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | 25000 | 420.73 M | 726.01 M | 1.38 G | 1.25 G | 1.34 G | 66.06 M | n/a |
| SUBTRACT | 25000 | 414.37 M | 743.76 M | 1.38 G | 1.22 G | 1.36 G | 84.61 M | n/a |
| ADDMOD | 25000 | 322.39 M | 577.49 M | 1.20 G | 1.06 G | 1.19 G | 23.54 M | n/a |
| SUBTRACTMOD | 25000 | 275.79 M | 531.60 M | 1.11 G | 1.08 G | 1.21 G | 24.28 M | n/a |
| MULTIPLYOPERANDSCANNING | 25000 | 4.93 M | 18.84 M | 75.96 M | 397.60 M | 433.14 M | 24.44 M | n/a |
| MULTIPLYPRODUCTSCANNING | 25000 | 22.24 M | 80.97 M | 279.77 M | 244.06 M | 277.04 M | 24.17 M | n/a |
| MONTGOMERYMULTIPLICATION | 25000 | 88.92 M | 286.19 M | 1.04 G | 660.84 M | 923.50 M | 3.05 M | n/a |
| COMPARE | 25000 | 430.08 M | 747.32 M | - | 1.36 G | 1.48 G | 162.58 M | n/a |
| REDUCE | 3125 | 25.45 M | 30.96 M | - | 61.38 M | 67.46 M | 33.08 M | n/a |
| MODMUL | 1562 | 2.68 M | 3.77 M | - | 7.19 M | 9.73 M | 5.81 M | n/a |
| MODEXP | 390 | 3.41 k | 34.65 k | - | 38.62 k | 91.50 k | 23.52 k | n/a |
| EXPONENTIATION | 390 | 4.61 k | 18.17 k | - | 51.75 k | 66.48 k | 104.59 k | n/a |
| DIVIDE | 3125 | 3.23 M | 3.24 M | - | 9.19 M | 10.62 M | 24.21 M | n/a |
| ISQRT | 781 | 106.18 k | 132.25 k | - | 387.33 k | 412.45 k | 7.16 M | n/a |
| MODMUL_R2 | 25000 | 60.11 M | 306.29 M | - | 465.46 M | 638.18 M | 6.33 M | n/a |

#### p1024 (1024-bit)

| Operation | items | w8 | w16 | w32 | w32-opt | w32-o64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | 12500 | 144.17 M | 273.55 M | 473.31 M | 484.53 M | 469.57 M | 57.95 M | n/a |
| SUBTRACT | 12500 | 153.44 M | 274.76 M | 485.27 M | 478.58 M | 469.93 M | 72.70 M | n/a |
| ADDMOD | 12500 | 99.97 M | 184.27 M | 364.28 M | 353.25 M | 348.61 M | 17.75 M | n/a |
| SUBTRACTMOD | 12500 | 100.87 M | 189.21 M | 358.01 M | 350.66 M | 355.96 M | 21.26 M | n/a |
| MULTIPLYOPERANDSCANNING | 12500 | 1.07 M | 4.21 M | 17.28 M | 148.13 M | 148.19 M | 6.93 M | n/a |
| MULTIPLYPRODUCTSCANNING | 12500 | 2.97 M | 11.51 M | 44.04 M | 43.56 M | 43.56 M | 6.93 M | n/a |
| MONTGOMERYMULTIPLICATION | 12500 | 16.78 M | 82.63 M | 252.82 M | 187.59 M | 226.11 M | 1.03 M | n/a |
| COMPARE | 12500 | 162.94 M | 295.37 M | - | 494.72 M | 494.50 M | 162.84 M | n/a |
| REDUCE | 1562 | 1.48 M | 2.66 M | - | 9.96 M | 9.88 M | 50.36 M | n/a |
| MODMUL | 781 | 199.31 k | 765.40 k | - | 1.40 M | 1.71 M | 2.01 M | n/a |
| MODEXP | 195 | 200.6 | 1.27 k | - | 2.89 k | 4.74 k | 3.47 k | n/a |
| EXPONENTIATION | 195 | 274.6 | 1.17 k | - | 3.71 k | 4.01 k | 24.99 k | n/a |
| DIVIDE | 1562 | 181.96 k | 183.51 k | - | 1.02 M | 1.04 M | 22.95 M | n/a |
| ISQRT | 390 | 3.56 k | 6.04 k | - | 31.23 k | 32.32 k | 3.67 M | n/a |
| MODMUL_R2 | 12500 | 10.33 M | 70.21 M | - | 123.20 M | 147.63 M | 2.28 M | n/a |

#### p2048 (2048-bit)

| Operation | items | w8 | w16 | w32 | w32-opt | w32-o64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | 6250 | 41.95 M | 80.94 M | 156.15 M | 156.35 M | 155.18 M | 41.21 M | n/a |
| SUBTRACT | 6250 | 41.98 M | 80.69 M | 157.21 M | 156.35 M | 157.30 M | 51.05 M | n/a |
| ADDMOD | 6250 | 31.15 M | 60.90 M | 125.06 M | 143.64 M | 144.07 M | 14.45 M | n/a |
| SUBTRACTMOD | 6250 | 30.45 M | 58.93 M | 119.07 M | 141.97 M | 145.45 M | 17.52 M | n/a |
| MULTIPLYOPERANDSCANNING | 6250 | 145.00 k | 569.75 k | 2.45 M | 41.37 M | 48.56 M | 2.05 M | n/a |
| MULTIPLYPRODUCTSCANNING | 6250 | 375.36 k | 1.45 M | 5.14 M | 5.66 M | 5.69 M | 2.02 M | n/a |
| MONTGOMERYMULTIPLICATION | 6250 | 289.71 k | 8.69 M | 49.89 M | 35.69 M | 43.54 M | 312.65 k | n/a |
| COMPARE | 6250 | 69.41 M | 128.86 M | - | 229.10 M | 232.77 M | 167.24 M | n/a |
| REDUCE | 781 | 24.92 k | 409.02 k | - | 2.38 M | 2.29 M | 32.27 M | n/a |
| MODMUL | 390 | 7.99 k | 55.01 k | - | 138.94 k | 179.50 k | 722.86 k | n/a |
| MODEXP | 97 | 6.3 | 61.9 | - | 245.4 | 210.0 | 500.8 | n/a |
| EXPONENTIATION | 97 | 14.7 | 72.2 | - | 271.0 | 261.9 | 4.18 k | n/a |
| DIVIDE | 781 | 4.98 k | 6.33 k | - | 49.30 k | 68.61 k | 17.90 M | n/a |
| ISQRT | 195 | 129.9 | 187.2 | - | 2.47 k | 7.11 k | 2.36 M | n/a |
| MODMUL_R2 | 6250 | 481.12 k | 8.24 M | - | 19.89 M | 26.32 M | 726.81 k | n/a |

### Device 1 - cpu-haswell-AMD Ryzen 9 3900X 12-Core Processor (CPU)

#### secp256k1 (256-bit)

| Operation | items | w8 | w16 | w32 | w32-opt | w32-o64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | 50000 | - | - | - | - | - | 73.17 M | n/a |
| SUBTRACT | 50000 | - | - | - | - | - | 94.83 M | n/a |
| ADDMOD | 50000 | - | - | - | - | - | 22.44 M | n/a |
| SUBTRACTMOD | 50000 | - | - | - | - | - | 25.50 M | n/a |
| MULTIPLYOPERANDSCANNING | 50000 | - | - | - | - | - | 58.35 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 | - | - | - | - | - | 58.70 M | n/a |
| MONTGOMERYMULTIPLICATION | 50000 | - | - | - | - | - | 7.07 M | n/a |
| COMPARE | 50000 | - | - | - | - | - | 164.83 M | n/a |
| REDUCE | 6250 | - | - | - | - | - | 59.28 M | n/a |
| MODMUL | 3125 | - | - | - | - | - | 12.06 M | n/a |
| MODEXP | 781 | - | - | - | - | - | 122.48 k | n/a |
| EXPONENTIATION | 781 | - | - | - | - | - | 366.83 k | n/a |
| DIVIDE | 6250 | - | - | - | - | - | 30.25 M | n/a |
| ISQRT | 1562 | - | - | - | - | - | 12.25 M | n/a |
| MODMUL_R2 | 50000 | - | - | - | - | - | 11.97 M | n/a |

#### rsa256(composite) (256-bit)

| Operation | items | w8 | w16 | w32 | w32-opt | w32-o64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | 50000 | - | - | - | - | - | 71.92 M | n/a |
| SUBTRACT | 50000 | - | - | - | - | - | 94.74 M | n/a |
| ADDMOD | 50000 | - | - | - | - | - | 25.81 M | n/a |
| SUBTRACTMOD | 50000 | - | - | - | - | - | 25.72 M | n/a |
| MULTIPLYOPERANDSCANNING | 50000 | - | - | - | - | - | 58.43 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 | - | - | - | - | - | 58.24 M | n/a |
| MONTGOMERYMULTIPLICATION | 50000 | - | - | - | - | - | 7.00 M | n/a |
| COMPARE | 50000 | - | - | - | - | - | 143.36 M | n/a |
| REDUCE | 6250 | - | - | - | - | - | 34.52 M | n/a |
| MODMUL | 3125 | - | - | - | - | - | 12.24 M | n/a |
| MODEXP | 781 | - | - | - | - | - | 129.48 k | n/a |
| EXPONENTIATION | 781 | - | - | - | - | - | 370.18 k | n/a |
| DIVIDE | 6250 | - | - | - | - | - | 27.61 M | n/a |
| ISQRT | 1562 | - | - | - | - | - | 12.19 M | n/a |
| MODMUL_R2 | 50000 | - | - | - | - | - | 11.78 M | n/a |

#### brainpoolP512r1 (512-bit)

| Operation | items | w8 | w16 | w32 | w32-opt | w32-o64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | 25000 | - | - | - | - | - | 66.06 M | n/a |
| SUBTRACT | 25000 | - | - | - | - | - | 84.61 M | n/a |
| ADDMOD | 25000 | - | - | - | - | - | 23.54 M | n/a |
| SUBTRACTMOD | 25000 | - | - | - | - | - | 24.28 M | n/a |
| MULTIPLYOPERANDSCANNING | 25000 | - | - | - | - | - | 24.44 M | n/a |
| MULTIPLYPRODUCTSCANNING | 25000 | - | - | - | - | - | 24.17 M | n/a |
| MONTGOMERYMULTIPLICATION | 25000 | - | - | - | - | - | 3.05 M | n/a |
| COMPARE | 25000 | - | - | - | - | - | 162.58 M | n/a |
| REDUCE | 3125 | - | - | - | - | - | 33.08 M | n/a |
| MODMUL | 1562 | - | - | - | - | - | 5.81 M | n/a |
| MODEXP | 390 | - | - | - | - | - | 23.52 k | n/a |
| EXPONENTIATION | 390 | - | - | - | - | - | 104.59 k | n/a |
| DIVIDE | 3125 | - | - | - | - | - | 24.21 M | n/a |
| ISQRT | 781 | - | - | - | - | - | 7.16 M | n/a |
| MODMUL_R2 | 25000 | - | - | - | - | - | 6.33 M | n/a |

#### p1024 (1024-bit)

| Operation | items | w8 | w16 | w32 | w32-opt | w32-o64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | 12500 | - | - | - | - | - | 57.95 M | n/a |
| SUBTRACT | 12500 | - | - | - | - | - | 72.70 M | n/a |
| ADDMOD | 12500 | - | - | - | - | - | 17.75 M | n/a |
| SUBTRACTMOD | 12500 | - | - | - | - | - | 21.26 M | n/a |
| MULTIPLYOPERANDSCANNING | 12500 | - | - | - | - | - | 6.93 M | n/a |
| MULTIPLYPRODUCTSCANNING | 12500 | - | - | - | - | - | 6.93 M | n/a |
| MONTGOMERYMULTIPLICATION | 12500 | - | - | - | - | - | 1.03 M | n/a |
| COMPARE | 12500 | - | - | - | - | - | 162.84 M | n/a |
| REDUCE | 1562 | - | - | - | - | - | 50.36 M | n/a |
| MODMUL | 781 | - | - | - | - | - | 2.01 M | n/a |
| MODEXP | 195 | - | - | - | - | - | 3.47 k | n/a |
| EXPONENTIATION | 195 | - | - | - | - | - | 24.99 k | n/a |
| DIVIDE | 1562 | - | - | - | - | - | 22.95 M | n/a |
| ISQRT | 390 | - | - | - | - | - | 3.67 M | n/a |
| MODMUL_R2 | 12500 | - | - | - | - | - | 2.28 M | n/a |

#### p2048 (2048-bit)

| Operation | items | w8 | w16 | w32 | w32-opt | w32-o64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | 6250 | - | - | - | - | - | 41.21 M | n/a |
| SUBTRACT | 6250 | - | - | - | - | - | 51.05 M | n/a |
| ADDMOD | 6250 | - | - | - | - | - | 14.45 M | n/a |
| SUBTRACTMOD | 6250 | - | - | - | - | - | 17.52 M | n/a |
| MULTIPLYOPERANDSCANNING | 6250 | - | - | - | - | - | 2.05 M | n/a |
| MULTIPLYPRODUCTSCANNING | 6250 | - | - | - | - | - | 2.02 M | n/a |
| MONTGOMERYMULTIPLICATION | 6250 | - | - | - | - | - | 312.65 k | n/a |
| COMPARE | 6250 | - | - | - | - | - | 167.24 M | n/a |
| REDUCE | 781 | - | - | - | - | - | 32.27 M | n/a |
| MODMUL | 390 | - | - | - | - | - | 722.86 k | n/a |
| MODEXP | 97 | - | - | - | - | - | 500.8 | n/a |
| EXPONENTIATION | 97 | - | - | - | - | - | 4.18 k | n/a |
| DIVIDE | 781 | - | - | - | - | - | 17.90 M | n/a |
| ISQRT | 195 | - | - | - | - | - | 2.36 M | n/a |
| MODMUL_R2 | 6250 | - | - | - | - | - | 726.81 k | n/a |

## 5. Head to head

Best OpenCL GPU result against best OpenCL CPU result and the CPU libraries.
Ratios above 1.00x mean the GPU is faster than that baseline.

### secp256k1 (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w32-o64 | 3.33 G | none | n/a | 73.17 M | n/a | n/a | n/a |
| SUBTRACT | w32 | 3.42 G | none | n/a | 94.83 M | n/a | n/a | n/a |
| ADDMOD | w32-o64 | 4.31 G | none | n/a | 22.44 M | n/a | n/a | n/a |
| SUBTRACTMOD | w32-opt | 4.27 G | none | n/a | 25.50 M | n/a | n/a | n/a |
| MULTIPLYOPERANDSCANNING | w32-o64 | 2.09 G | none | n/a | 58.35 M | n/a | n/a | n/a |
| MULTIPLYPRODUCTSCANNING | w32 | 1.77 G | none | n/a | 58.70 M | n/a | n/a | n/a |
| MONTGOMERYMULTIPLICATION | w32 | 3.77 G | none | n/a | 7.07 M | n/a | n/a | n/a |
| COMPARE | w32-opt | 4.24 G | none | n/a | 164.83 M | n/a | n/a | n/a |
| REDUCE | w32-o64 | 243.11 M | none | n/a | 59.28 M | n/a | n/a | n/a |
| MODMUL | w32-o64 | 53.45 M | none | n/a | 12.06 M | n/a | n/a | n/a |
| MODEXP | w32-o64 | 1.50 M | none | n/a | 122.48 k | n/a | n/a | n/a |
| EXPONENTIATION | w32-o64 | 7.92 M | none | n/a | 366.83 k | n/a | n/a | n/a |
| DIVIDE | w32-o64 | 73.85 M | none | n/a | 30.25 M | n/a | n/a | n/a |
| ISQRT | w32-o64 | 3.34 M | none | n/a | 12.25 M | n/a | n/a | n/a |
| MODMUL_R2 | w32-o64 | 2.06 G | none | n/a | 11.97 M | n/a | n/a | n/a |

### rsa256(composite) (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w32-o64 | 3.38 G | none | n/a | 71.92 M | n/a | n/a | n/a |
| SUBTRACT | w32-o64 | 3.42 G | none | n/a | 94.74 M | n/a | n/a | n/a |
| ADDMOD | w32-o64 | 4.37 G | none | n/a | 25.81 M | n/a | n/a | n/a |
| SUBTRACTMOD | w32-opt | 4.26 G | none | n/a | 25.72 M | n/a | n/a | n/a |
| MULTIPLYOPERANDSCANNING | w32-o64 | 2.04 G | none | n/a | 58.43 M | n/a | n/a | n/a |
| MULTIPLYPRODUCTSCANNING | w32-o64 | 1.79 G | none | n/a | 58.24 M | n/a | n/a | n/a |
| MONTGOMERYMULTIPLICATION | w32 | 3.67 G | none | n/a | 7.00 M | n/a | n/a | n/a |
| COMPARE | w32-opt | 4.27 G | none | n/a | 143.36 M | n/a | n/a | n/a |
| REDUCE | w32-o64 | 239.20 M | none | n/a | 34.52 M | n/a | n/a | n/a |
| MODMUL | w32-o64 | 55.33 M | none | n/a | 12.24 M | n/a | n/a | n/a |
| MODEXP | w32-o64 | 1.61 M | none | n/a | 129.48 k | n/a | n/a | n/a |
| EXPONENTIATION | w32-o64 | 7.89 M | none | n/a | 370.18 k | n/a | n/a | n/a |
| DIVIDE | w32-o64 | 70.06 M | none | n/a | 27.61 M | n/a | n/a | n/a |
| ISQRT | w32-o64 | 2.84 M | none | n/a | 12.19 M | n/a | n/a | n/a |
| MODMUL_R2 | w32-o64 | 2.09 G | none | n/a | 11.78 M | n/a | n/a | n/a |

### brainpoolP512r1 (512-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w32 | 1.38 G | none | n/a | 66.06 M | n/a | n/a | n/a |
| SUBTRACT | w32 | 1.38 G | none | n/a | 84.61 M | n/a | n/a | n/a |
| ADDMOD | w32 | 1.20 G | none | n/a | 23.54 M | n/a | n/a | n/a |
| SUBTRACTMOD | w32-o64 | 1.21 G | none | n/a | 24.28 M | n/a | n/a | n/a |
| MULTIPLYOPERANDSCANNING | w32-o64 | 433.14 M | none | n/a | 24.44 M | n/a | n/a | n/a |
| MULTIPLYPRODUCTSCANNING | w32 | 279.77 M | none | n/a | 24.17 M | n/a | n/a | n/a |
| MONTGOMERYMULTIPLICATION | w32 | 1.04 G | none | n/a | 3.05 M | n/a | n/a | n/a |
| COMPARE | w32-o64 | 1.48 G | none | n/a | 162.58 M | n/a | n/a | n/a |
| REDUCE | w32-o64 | 67.46 M | none | n/a | 33.08 M | n/a | n/a | n/a |
| MODMUL | w32-o64 | 9.73 M | none | n/a | 5.81 M | n/a | n/a | n/a |
| MODEXP | w32-o64 | 91.50 k | none | n/a | 23.52 k | n/a | n/a | n/a |
| EXPONENTIATION | w32-o64 | 66.48 k | none | n/a | 104.59 k | n/a | n/a | n/a |
| DIVIDE | w32-o64 | 10.62 M | none | n/a | 24.21 M | n/a | n/a | n/a |
| ISQRT | w32-o64 | 412.45 k | none | n/a | 7.16 M | n/a | n/a | n/a |
| MODMUL_R2 | w32-o64 | 638.18 M | none | n/a | 6.33 M | n/a | n/a | n/a |

### p1024 (1024-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w32-opt | 484.53 M | none | n/a | 57.95 M | n/a | n/a | n/a |
| SUBTRACT | w32 | 485.27 M | none | n/a | 72.70 M | n/a | n/a | n/a |
| ADDMOD | w32 | 364.28 M | none | n/a | 17.75 M | n/a | n/a | n/a |
| SUBTRACTMOD | w32 | 358.01 M | none | n/a | 21.26 M | n/a | n/a | n/a |
| MULTIPLYOPERANDSCANNING | w32-o64 | 148.19 M | none | n/a | 6.93 M | n/a | n/a | n/a |
| MULTIPLYPRODUCTSCANNING | w32 | 44.04 M | none | n/a | 6.93 M | n/a | n/a | n/a |
| MONTGOMERYMULTIPLICATION | w32 | 252.82 M | none | n/a | 1.03 M | n/a | n/a | n/a |
| COMPARE | w32-opt | 494.72 M | none | n/a | 162.84 M | n/a | n/a | n/a |
| REDUCE | w32-opt | 9.96 M | none | n/a | 50.36 M | n/a | n/a | n/a |
| MODMUL | w32-o64 | 1.71 M | none | n/a | 2.01 M | n/a | n/a | n/a |
| MODEXP | w32-o64 | 4.74 k | none | n/a | 3.47 k | n/a | n/a | n/a |
| EXPONENTIATION | w32-o64 | 4.01 k | none | n/a | 24.99 k | n/a | n/a | n/a |
| DIVIDE | w32-o64 | 1.04 M | none | n/a | 22.95 M | n/a | n/a | n/a |
| ISQRT | w32-o64 | 32.32 k | none | n/a | 3.67 M | n/a | n/a | n/a |
| MODMUL_R2 | w32-o64 | 147.63 M | none | n/a | 2.28 M | n/a | n/a | n/a |

### p2048 (2048-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w32-opt | 156.35 M | none | n/a | 41.21 M | n/a | n/a | n/a |
| SUBTRACT | w32-o64 | 157.30 M | none | n/a | 51.05 M | n/a | n/a | n/a |
| ADDMOD | w32-o64 | 144.07 M | none | n/a | 14.45 M | n/a | n/a | n/a |
| SUBTRACTMOD | w32-o64 | 145.45 M | none | n/a | 17.52 M | n/a | n/a | n/a |
| MULTIPLYOPERANDSCANNING | w32-o64 | 48.56 M | none | n/a | 2.05 M | n/a | n/a | n/a |
| MULTIPLYPRODUCTSCANNING | w32-o64 | 5.69 M | none | n/a | 2.02 M | n/a | n/a | n/a |
| MONTGOMERYMULTIPLICATION | w32 | 49.89 M | none | n/a | 312.65 k | n/a | n/a | n/a |
| COMPARE | w32-o64 | 232.77 M | none | n/a | 167.24 M | n/a | n/a | n/a |
| REDUCE | w32-opt | 2.38 M | none | n/a | 32.27 M | n/a | n/a | n/a |
| MODMUL | w32-o64 | 179.50 k | none | n/a | 722.86 k | n/a | n/a | n/a |
| MODEXP | w32-opt | 245.4 | none | n/a | 500.8 | n/a | n/a | n/a |
| EXPONENTIATION | w32-opt | 271.0 | none | n/a | 4.18 k | n/a | n/a | n/a |
| DIVIDE | w32-o64 | 68.61 k | none | n/a | 17.90 M | n/a | n/a | n/a |
| ISQRT | w32-o64 | 7.11 k | none | n/a | 2.36 M | n/a | n/a | n/a |
| MODMUL_R2 | w32-o64 | 26.32 M | none | n/a | 726.81 k | n/a | n/a | n/a |

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

Also written to `NVIDIA_GeForce_RTX_3080_Ti_Report.csv` for analysis.

```csv
kind,device,device_type,kernel,modulus,bits,operation,items,seconds,ops_per_sec,mismatches
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,secp256k1,256,ADD,50000,0.000683382,73165529.387,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,secp256k1,256,ADD,50000,0.007013446,7129163.087,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,secp256k1,256,ADD,50000,0.000070662,707593848.911,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,secp256k1,256,ADD,50000,0.000044984,1111505213.633,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,secp256k1,256,ADD,50000,0.000675337,74037118.753,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,secp256k1,256,ADD,50000,0.000024306,2057114603.470,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,secp256k1,256,ADD,50000,0.000736071,67928227.011,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32,secp256k1,256,ADD,50000,0.000015289,3270362671.134,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32,secp256k1,256,ADD,50000,0.000692199,72233564.527,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,secp256k1,256,ADD,50000,0.000016170,3092129082.793,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,secp256k1,256,ADD,50000,0.000653085,76559732.304,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,secp256k1,256,ADD,50000,0.000015009,3331343015.373,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,secp256k1,256,ADD,50000,0.000686287,72855797.618,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,secp256k1,256,SUBTRACT,50000,0.000527239,94833663.270,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,secp256k1,256,SUBTRACT,50000,0.000059011,847300709.410,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,secp256k1,256,SUBTRACT,50000,0.000062908,794813831.902,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,secp256k1,256,SUBTRACT,50000,0.000045135,1107789741.712,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,secp256k1,256,SUBTRACT,50000,0.000672993,74294975.892,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,secp256k1,256,SUBTRACT,50000,0.000023885,2093349626.655,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,secp256k1,256,SUBTRACT,50000,0.000721213,69327652.178,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32,secp256k1,256,SUBTRACT,50000,0.000014637,3415969916.966,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32,secp256k1,256,SUBTRACT,50000,0.000686528,72830249.049,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,secp256k1,256,SUBTRACT,50000,0.000016310,3065600274.086,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,secp256k1,256,SUBTRACT,50000,0.000657123,76089250.199,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,secp256k1,256,SUBTRACT,50000,0.000014898,3356124912.873,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,secp256k1,256,SUBTRACT,50000,0.000688923,72577069.932,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,secp256k1,256,ADDMOD,50000,0.002228520,22436414.071,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,secp256k1,256,ADDMOD,50000,0.000141335,353769531.274,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,secp256k1,256,ADDMOD,50000,0.000500599,99880312.382,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,secp256k1,256,ADDMOD,50000,0.000066034,757185346.281,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,secp256k1,256,ADDMOD,50000,0.000687610,72715624.988,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,secp256k1,256,ADDMOD,50000,0.000033142,1508657651.886,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,secp256k1,256,ADDMOD,50000,0.000728988,68588226.306,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32,secp256k1,256,ADDMOD,50000,0.000016331,3061666711.339,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32,secp256k1,256,ADDMOD,50000,0.000687630,72713507.547,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,secp256k1,256,ADDMOD,50000,0.000012263,4077319956.711,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,secp256k1,256,ADDMOD,50000,0.000651713,76720888.461,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,secp256k1,256,ADDMOD,50000,0.000011591,4313688704.980,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,secp256k1,256,ADDMOD,50000,0.000676379,73923043.993,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,secp256k1,256,SUBTRACTMOD,50000,0.001961149,25495257.585,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,secp256k1,256,SUBTRACTMOD,50000,0.000196578,254351986.507,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,secp256k1,256,SUBTRACTMOD,50000,0.000728436,68640205.535,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.000066044,757070563.392,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.000686628,72819629.677,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,secp256k1,256,SUBTRACTMOD,50000,0.000033874,1476065660.849,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,secp256k1,256,SUBTRACTMOD,50000,0.000737945,67755719.383,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32,secp256k1,256,SUBTRACTMOD,50000,0.000016971,2946197898.203,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32,secp256k1,256,SUBTRACTMOD,50000,0.000686157,72869617.125,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,secp256k1,256,SUBTRACTMOD,50000,0.000011712,4269181440.102,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,secp256k1,256,SUBTRACTMOD,50000,0.000655881,76233326.092,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,secp256k1,256,SUBTRACTMOD,50000,0.000012093,4134626481.064,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,secp256k1,256,SUBTRACTMOD,50000,0.000663606,75345888.434,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000856928,58347970.874,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000131126,381312550.472,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000242966,205790115.291,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001546641,32328120.320,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.002357512,21208799.660,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000391234,127800771.393,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001318053,37934738.469,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000099297,503540319.408,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000943751,52980072.137,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000028283,1767854559.824,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000838003,59665665.461,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000023915,2090740938.918,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000905820,55198594.308,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000851798,58699371.864,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000127088,393428034.274,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000243858,205037384.209,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000275356,181583116.489,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001062183,47072870.410,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000083467,599039200.196,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000950614,52597584.828,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000028313,1765964645.900,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000869512,57503513.110,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000031609,1581823547.437,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000843773,59257651.139,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000028353,1763484826.935,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000872958,57276525.114,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.007076634,7065505.908,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000876004,57077373.002,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000268463,186245372.927,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000137188,364463347.691,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000791394,63179643.026,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000035226,1419401598.202,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000705664,70855246.595,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000013255,3772147633.936,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000690505,72410756.999,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000018365,2722572673.910,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000660529,75696888.949,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000017302,2889820820.325,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000688612,72609830.110,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,secp256k1,256,COMPARE,50000,0.000303339,164832096.258,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,secp256k1,256,COMPARE,50000,0.000048060,1040366855.283,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,secp256k1,256,COMPARE,50000,0.000072386,690742420.431,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,secp256k1,256,COMPARE,50000,0.000043932,1138119556.726,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,secp256k1,256,COMPARE,50000,0.000713117,70114726.131,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,secp256k1,256,COMPARE,50000,0.000023244,2151097491.786,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,secp256k1,256,COMPARE,50000,0.000687720,72704029.814,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,secp256k1,256,COMPARE,50000,0.000011802,4236587125.412,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,secp256k1,256,COMPARE,50000,0.000670267,74597134.128,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,secp256k1,256,COMPARE,50000,0.000012815,3901678139.535,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,secp256k1,256,COMPARE,50000,0.000685176,72973939.117,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,secp256k1,256,REDUCE,6250,0.000105428,59282139.529,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,secp256k1,256,REDUCE,6250,0.000020448,305655074.411,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,secp256k1,256,REDUCE,6250,0.000090610,68976931.754,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,secp256k1,256,REDUCE,6250,0.000079990,78134638.122,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,secp256k1,256,REDUCE,6250,0.000250970,24903373.572,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,secp256k1,256,REDUCE,6250,0.000050686,123308048.416,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,secp256k1,256,REDUCE,6250,0.000192912,32398181.403,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,secp256k1,256,REDUCE,6250,0.000028964,215784128.617,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,secp256k1,256,REDUCE,6250,0.000167885,37227842.859,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,secp256k1,256,REDUCE,6250,0.000025709,243106219.038,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,secp256k1,256,REDUCE,6250,0.000200596,31157161.080,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,secp256k1,256,MODMUL,3125,0.000259086,12061630.927,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,secp256k1,256,MODMUL,3125,0.000037360,83645598.903,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,secp256k1,256,MODMUL,3125,0.000105679,29570671.179,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,secp256k1,256,MODMUL,3125,0.000147347,21208424.732,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,secp256k1,256,MODMUL,3125,0.000269425,11598780.820,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,secp256k1,256,MODMUL,3125,0.000086382,36176505.089,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,secp256k1,256,MODMUL,3125,0.000200487,15587047.504,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,secp256k1,256,MODMUL,3125,0.000070593,44267790.696,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,secp256k1,256,MODMUL,3125,0.000179026,17455560.917,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,secp256k1,256,MODMUL,3125,0.000058469,53447220.076,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,secp256k1,256,MODMUL,3125,0.000179456,17413731.418,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,secp256k1,256,MODEXP,781,0.006376491,122481.159,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,secp256k1,256,MODEXP,781,0.000927459,842085.757,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,secp256k1,256,MODEXP,781,0.001803603,433022.173,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,secp256k1,256,MODEXP,781,0.007771137,100500.095,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,secp256k1,256,MODEXP,781,0.007792226,100228.099,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,secp256k1,256,MODEXP,781,0.001307493,597326.299,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,secp256k1,256,MODEXP,781,0.001326710,588674.246,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,secp256k1,256,MODEXP,781,0.000908254,859891.666,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,secp256k1,256,MODEXP,781,0.000926479,842976.474,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,secp256k1,256,MODEXP,781,0.000519315,1503904.142,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,secp256k1,256,MODEXP,781,0.000540354,1445348.639,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,secp256k1,256,EXPONENTIATION,781,0.002129023,366834.895,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,secp256k1,256,EXPONENTIATION,781,0.000304602,2564000.821,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,secp256k1,256,EXPONENTIATION,781,0.004510881,173136.914,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,secp256k1,256,EXPONENTIATION,781,0.009831411,79439.258,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,secp256k1,256,EXPONENTIATION,781,0.009848213,79303.727,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,secp256k1,256,EXPONENTIATION,781,0.002591143,301411.370,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,secp256k1,256,EXPONENTIATION,781,0.002612091,298994.189,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,secp256k1,256,EXPONENTIATION,781,0.000173526,4500763.406,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,secp256k1,256,EXPONENTIATION,781,0.000188133,4151318.905,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,secp256k1,256,EXPONENTIATION,781,0.000098655,7916476.584,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,secp256k1,256,EXPONENTIATION,781,0.000113933,6854899.074,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,secp256k1,256,DIVIDE,6250,0.000206587,30253611.696,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,secp256k1,256,DIVIDE,6250,0.000032972,189554320.900,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,secp256k1,256,DIVIDE,6250,0.000089007,70219224.552,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,secp256k1,256,DIVIDE,6250,0.000304200,20545694.432,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,secp256k1,256,DIVIDE,6250,0.000503484,12413504.516,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,secp256k1,256,DIVIDE,6250,0.000218450,28610684.530,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,secp256k1,256,DIVIDE,6250,0.000412724,15143292.588,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,secp256k1,256,DIVIDE,6250,0.000098294,63584833.562,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,secp256k1,256,DIVIDE,6250,0.000265388,23550419.754,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,secp256k1,256,DIVIDE,6250,0.000084629,73851709.727,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,secp256k1,256,DIVIDE,6250,0.000288682,21650115.818,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,secp256k1,256,ISQRT,1562,0.000127489,12252039.346,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,secp256k1,256,ISQRT,1562,0.000018554,84187567.969,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,secp256k1,256,ISQRT,1562,0.001806388,864709.064,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,secp256k1,256,ISQRT,1562,0.001835703,850900.211,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,secp256k1,256,ISQRT,1562,0.001484175,1052436.553,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,secp256k1,256,ISQRT,1562,0.001514141,1031608.033,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,secp256k1,256,ISQRT,1562,0.000538911,2898437.926,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,secp256k1,256,ISQRT,1562,0.000562787,2775472.505,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,secp256k1,256,ISQRT,1562,0.000467367,3342127.246,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,secp256k1,256,ISQRT,1562,0.000490420,3185024.871,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,secp256k1,256,MODMUL_R2,50000,0.004175913,11973429.367,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,secp256k1,256,MODMUL_R2,50000,0.000279555,178855619.149,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,secp256k1,256,MODMUL_R2,50000,0.000727024,68773505.202,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,secp256k1,256,MODMUL_R2,50000,0.000136737,365665439.202,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,secp256k1,256,MODMUL_R2,50000,0.000797997,62656877.473,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,secp256k1,256,MODMUL_R2,50000,0.000037250,1342277950.846,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,secp256k1,256,MODMUL_R2,50000,0.000693151,72134351.564,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,secp256k1,256,MODMUL_R2,50000,0.000024566,2035336601.270,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,secp256k1,256,MODMUL_R2,50000,0.000708490,70572635.092,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,secp256k1,256,MODMUL_R2,50000,0.000024275,2059738776.137,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,secp256k1,256,MODMUL_R2,50000,0.000701637,71261919.597,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,ADD,50000,0.000695234,71918238.906,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,ADD,50000,0.000055815,895815040.630,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,ADD,50000,0.000090951,549747369.768,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,rsa256(composite),256,ADD,50000,0.000040676,1229226710.627,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,rsa256(composite),256,ADD,50000,0.000693311,72117709.356,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,rsa256(composite),256,ADD,50000,0.000024175,2068249027.747,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,rsa256(composite),256,ADD,50000,0.000671139,74500216.755,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32,rsa256(composite),256,ADD,50000,0.000014888,3358381784.061,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32,rsa256(composite),256,ADD,50000,0.000691017,72357105.775,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,rsa256(composite),256,ADD,50000,0.000016171,3091951000.662,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,rsa256(composite),256,ADD,50000,0.000653877,76466989.795,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,rsa256(composite),256,ADD,50000,0.000014808,3376546616.352,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,rsa256(composite),256,ADD,50000,0.000672372,74363589.924,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACT,50000,0.000527770,94738233.642,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACT,50000,0.000108434,461110462.939,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACT,50000,0.000129223,386928413.284,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,rsa256(composite),256,SUBTRACT,50000,0.000041067,1217518594.868,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,rsa256(composite),256,SUBTRACT,50000,0.000695455,71895389.420,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,rsa256(composite),256,SUBTRACT,50000,0.000024346,2053730835.366,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,rsa256(composite),256,SUBTRACT,50000,0.000672792,74317190.190,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32,rsa256(composite),256,SUBTRACT,50000,0.000014787,3381384761.215,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32,rsa256(composite),256,SUBTRACT,50000,0.000687129,72766529.920,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,rsa256(composite),256,SUBTRACT,50000,0.000016311,3065381477.675,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,rsa256(composite),256,SUBTRACT,50000,0.000648577,77091848.879,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,rsa256(composite),256,SUBTRACT,50000,0.000014608,3422775614.032,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,rsa256(composite),256,SUBTRACT,50000,0.000670077,74618284.972,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,ADDMOD,50000,0.001937073,25812140.849,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,ADDMOD,50000,0.000282801,176802776.175,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,ADDMOD,50000,0.003595645,13905710.651,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,rsa256(composite),256,ADDMOD,50000,0.000054832,911874908.918,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,rsa256(composite),256,ADDMOD,50000,0.000728687,68616562.870,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,rsa256(composite),256,ADDMOD,50000,0.000030798,1623474714.425,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,rsa256(composite),256,ADDMOD,50000,0.000685927,72894055.126,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32,rsa256(composite),256,ADDMOD,50000,0.000016010,3123022043.832,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32,rsa256(composite),256,ADDMOD,50000,0.000685416,72948406.757,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,rsa256(composite),256,ADDMOD,50000,0.000012293,4067357945.377,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,rsa256(composite),256,ADDMOD,50000,0.000697538,71680680.182,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,rsa256(composite),256,ADDMOD,50000,0.000011451,4366490408.898,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,rsa256(composite),256,ADDMOD,50000,0.000684624,73032805.940,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACTMOD,50000,0.001943926,25721145.542,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACTMOD,50000,0.000193794,258005890.388,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACTMOD,50000,0.004530368,11036630.679,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,rsa256(composite),256,SUBTRACTMOD,50000,0.000059362,842289180.179,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,rsa256(composite),256,SUBTRACTMOD,50000,0.000721353,69314181.321,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,rsa256(composite),256,SUBTRACTMOD,50000,0.000033834,1477802614.991,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,rsa256(composite),256,SUBTRACTMOD,50000,0.000687370,72741019.269,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32,rsa256(composite),256,SUBTRACTMOD,50000,0.000016781,2979553858.534,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32,rsa256(composite),256,SUBTRACTMOD,50000,0.000689875,72476887.307,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,rsa256(composite),256,SUBTRACTMOD,50000,0.000011742,4258261085.443,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,rsa256(composite),256,SUBTRACTMOD,50000,0.000685546,72934557.353,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,rsa256(composite),256,SUBTRACTMOD,50000,0.000012223,4090678796.884,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,rsa256(composite),256,SUBTRACTMOD,50000,0.000686939,72786655.283,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000855685,58432719.089,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000130595,382863223.278,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000191509,261084301.146,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001371403,36459013.356,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.002206639,22658895.174,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000397426,125809634.250,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001217795,41057811.138,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000098926,505428222.291,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000951926,52525091.727,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000027982,1786874504.289,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000872257,57322544.281,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000024486,2041994226.271,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000933011,53589946.025,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000858470,58243146.797,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000069009,724541704.230,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000239459,208803925.835,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000244919,204148951.251,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001072912,46602143.954,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000083857,596253265.327,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000992954,50354799.406,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000028294,1767170816.567,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000869322,57516080.467,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000031469,1588869062.875,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000926879,53944473.770,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000027982,1786859636.218,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000869171,57526079.751,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.007139883,7002915.870,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000996218,50189826.451,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000266039,187942337.189,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000121938,410044231.269,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000790173,63277285.919,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000037340,1339047257.037,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000688993,72569687.639,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000013625,3669720343.820,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000689043,72564415.491,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000018535,2697606552.188,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000673654,74222075.801,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000017623,2837209205.972,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000683112,73194457.016,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,COMPARE,50000,0.000348765,143363035.221,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,COMPARE,50000,0.000048842,1023708091.050,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,COMPARE,50000,0.000072947,685429101.993,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,rsa256(composite),256,COMPARE,50000,0.000039264,1273428081.453,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,rsa256(composite),256,COMPARE,50000,0.000701706,71254896.982,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,rsa256(composite),256,COMPARE,50000,0.000023834,2097848551.277,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,rsa256(composite),256,COMPARE,50000,0.000647995,77161098.437,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,rsa256(composite),256,COMPARE,50000,0.000011712,4269096570.781,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,rsa256(composite),256,COMPARE,50000,0.000700876,71339307.211,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,rsa256(composite),256,COMPARE,50000,0.000012925,3868503473.123,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,rsa256(composite),256,COMPARE,50000,0.000684764,73017856.887,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,REDUCE,6250,0.000181030,34524661.294,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,REDUCE,6250,0.000027612,226350169.066,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,REDUCE,6250,0.000091873,68028610.818,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,rsa256(composite),256,REDUCE,6250,0.000072145,86631206.351,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,rsa256(composite),256,REDUCE,6250,0.000233729,26740354.572,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,rsa256(composite),256,REDUCE,6250,0.000051116,122271218.588,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,rsa256(composite),256,REDUCE,6250,0.000184726,33833900.014,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,rsa256(composite),256,REDUCE,6250,0.000029215,213930295.351,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,rsa256(composite),256,REDUCE,6250,0.000200646,31149387.770,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,rsa256(composite),256,REDUCE,6250,0.000026129,239198253.477,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,rsa256(composite),256,REDUCE,6250,0.000200196,31219414.888,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MODMUL,3125,0.000255349,12238148.071,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MODMUL,3125,0.000037841,82582819.874,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MODMUL,3125,0.000106720,29282228.122,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,rsa256(composite),256,MODMUL,3125,0.000127950,24423605.882,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,rsa256(composite),256,MODMUL,3125,0.000254037,12301353.245,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,rsa256(composite),256,MODMUL,3125,0.000084539,36965212.123,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,rsa256(composite),256,MODMUL,3125,0.000183945,16988768.642,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,rsa256(composite),256,MODMUL,3125,0.000068488,45628387.851,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,rsa256(composite),256,MODMUL,3125,0.000196579,15896923.842,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,rsa256(composite),256,MODMUL,3125,0.000056476,55333369.613,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,rsa256(composite),256,MODMUL,3125,0.000178154,17541017.689,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MODEXP,781,0.006031804,129480.339,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MODEXP,781,0.000903404,864508.158,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MODEXP,781,0.001821927,428667.004,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,rsa256(composite),256,MODEXP,781,0.006748489,115729.608,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,rsa256(composite),256,MODEXP,781,0.006767855,115398.453,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,rsa256(composite),256,MODEXP,781,0.001232092,633881.218,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,rsa256(composite),256,MODEXP,781,0.001254134,622740.532,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,rsa256(composite),256,MODEXP,781,0.000867408,900383.670,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,rsa256(composite),256,MODEXP,781,0.000884259,883225.503,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,rsa256(composite),256,MODEXP,781,0.000485481,1608713.846,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,rsa256(composite),256,MODEXP,781,0.000506170,1542959.589,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,EXPONENTIATION,781,0.002109807,370176.042,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,EXPONENTIATION,781,0.000315662,2474165.344,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,EXPONENTIATION,781,0.007725712,101091.007,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,rsa256(composite),256,EXPONENTIATION,781,0.009904549,78852.657,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,rsa256(composite),256,EXPONENTIATION,781,0.009921671,78716.578,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,rsa256(composite),256,EXPONENTIATION,781,0.002592825,301215.871,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,rsa256(composite),256,EXPONENTIATION,781,0.002614226,298749.998,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,rsa256(composite),256,EXPONENTIATION,781,0.000171672,4549370.914,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,rsa256(composite),256,EXPONENTIATION,781,0.000188443,4144491.989,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,rsa256(composite),256,EXPONENTIATION,781,0.000098926,7894788.832,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,rsa256(composite),256,EXPONENTIATION,781,0.000113673,6870582.354,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,DIVIDE,6250,0.000226395,27606620.189,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,DIVIDE,6250,0.005245830,1191422.536,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,DIVIDE,6250,0.000049824,125441232.190,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,rsa256(composite),256,DIVIDE,6250,0.000310112,20154007.332,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,rsa256(composite),256,DIVIDE,6250,0.000512191,12202476.903,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,rsa256(composite),256,DIVIDE,6250,0.000223289,27990629.591,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,rsa256(composite),256,DIVIDE,6250,0.000395902,15786736.423,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,rsa256(composite),256,DIVIDE,6250,0.000103464,60407416.259,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,rsa256(composite),256,DIVIDE,6250,0.000306075,20419817.996,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,rsa256(composite),256,DIVIDE,6250,0.000089207,70061793.439,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,rsa256(composite),256,DIVIDE,6250,0.000288712,21647863.524,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,ISQRT,1562,0.000128181,12185897.878,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,ISQRT,1562,0.000012403,125938406.539,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,rsa256(composite),256,ISQRT,1562,0.002037522,766617.520,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,rsa256(composite),256,ISQRT,1562,0.002064453,756616.896,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,rsa256(composite),256,ISQRT,1562,0.001696243,920858.654,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,rsa256(composite),256,ISQRT,1562,0.001719376,908469.081,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,rsa256(composite),256,ISQRT,1562,0.000629441,2481566.528,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,rsa256(composite),256,ISQRT,1562,0.000653066,2391795.106,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,rsa256(composite),256,ISQRT,1562,0.000550463,2837610.746,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,rsa256(composite),256,ISQRT,1562,0.000579057,2697490.145,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MODMUL_R2,50000,0.004244161,11780891.243,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MODMUL_R2,50000,0.007000150,7142704.160,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MODMUL_R2,50000,0.000723217,69135551.075,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,rsa256(composite),256,MODMUL_R2,50000,0.000136817,365451375.963,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,rsa256(composite),256,MODMUL_R2,50000,0.000795542,62850230.579,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,rsa256(composite),256,MODMUL_R2,50000,0.000036478,1370696330.527,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,rsa256(composite),256,MODMUL_R2,50000,0.000672141,74389169.270,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,rsa256(composite),256,MODMUL_R2,50000,0.000024967,2002633189.410,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,rsa256(composite),256,MODMUL_R2,50000,0.000691999,72254441.470,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,rsa256(composite),256,MODMUL_R2,50000,0.000023935,2088991875.486,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,rsa256(composite),256,MODMUL_R2,50000,0.000700384,71389394.397,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,ADD,25000,0.000378470,66055444.794,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,ADD,25000,0.000037841,660658494.026,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,ADD,25000,0.000046427,538478269.635,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,brainpoolP512r1,512,ADD,25000,0.000059421,420727094.052,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,brainpoolP512r1,512,ADD,25000,0.000719680,34737661.897,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,brainpoolP512r1,512,ADD,25000,0.000034435,726006493.708,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,brainpoolP512r1,512,ADD,25000,0.000711906,35116995.643,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32,brainpoolP512r1,512,ADD,25000,0.000018134,1378624669.705,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32,brainpoolP512r1,512,ADD,25000,0.000692991,36075500.727,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,brainpoolP512r1,512,ADD,25000,0.000019958,1252629900.022,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,brainpoolP512r1,512,ADD,25000,0.000676119,36975739.391,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,brainpoolP512r1,512,ADD,25000,0.000018625,1342277950.846,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,brainpoolP512r1,512,ADD,25000,0.000691788,36138236.885,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACT,25000,0.000295484,84606957.719,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACT,25000,0.000061335,407597339.731,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACT,25000,0.000057177,437239364.262,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,brainpoolP512r1,512,SUBTRACT,25000,0.000060333,414367348.955,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,brainpoolP512r1,512,SUBTRACT,25000,0.000722065,34622928.018,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,brainpoolP512r1,512,SUBTRACT,25000,0.000033613,743758493.284,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,brainpoolP512r1,512,SUBTRACT,25000,0.000663615,37672441.947,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32,brainpoolP512r1,512,SUBTRACT,25000,0.000018064,1383955434.685,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32,brainpoolP512r1,512,SUBTRACT,25000,0.000662213,37752205.959,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,brainpoolP512r1,512,SUBTRACT,25000,0.000020489,1220161163.636,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,brainpoolP512r1,512,SUBTRACT,25000,0.000695154,35963250.520,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,brainpoolP512r1,512,SUBTRACT,25000,0.000018404,1358393097.603,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,brainpoolP512r1,512,SUBTRACT,25000,0.000694183,36013549.717,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,ADDMOD,25000,0.001061972,23541110.360,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,ADDMOD,25000,0.000142207,175800185.339,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,ADDMOD,25000,0.000374262,66798127.461,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,brainpoolP512r1,512,ADDMOD,25000,0.000077546,322388840.382,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,brainpoolP512r1,512,ADDMOD,25000,0.000727194,34378726.426,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,brainpoolP512r1,512,ADDMOD,25000,0.000043291,577488570.614,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,brainpoolP512r1,512,ADDMOD,25000,0.000665829,37547175.137,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32,brainpoolP512r1,512,ADDMOD,25000,0.000020779,1203139474.480,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32,brainpoolP512r1,512,ADDMOD,25000,0.000752963,33202157.113,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,brainpoolP512r1,512,ADDMOD,25000,0.000023574,1060485752.099,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,brainpoolP512r1,512,ADDMOD,25000,0.000699102,35760163.590,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,brainpoolP512r1,512,ADDMOD,25000,0.000020959,1192807909.529,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,brainpoolP512r1,512,ADDMOD,25000,0.000739308,33815413.017,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACTMOD,25000,0.001029551,24282429.225,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000140533,177894053.696,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000516810,48373678.031,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000090650,275786546.394,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000744026,33600979.857,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000047028,531597465.146,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000686549,36414000.718,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000022583,1107030222.800,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000693482,36049968.541,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000023084,1083001486.711,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000700273,35700358.220,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000020618,1212527750.299,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000700214,35703373.417,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.001022749,24443922.404,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.002087956,11973433.372,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.000149360,167380907.130,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.005075892,4925242.727,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.005928040,4217245.486,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.001326951,18840183.738,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.002128324,11746331.824,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.000329108,75962910.999,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.001172951,21313764.102,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.000062878,397595275.107,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.000912021,27411652.090,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.000057718,433140305.853,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.000894769,27940169.305,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.001034331,24170211.448,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000154460,161854115.987,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000131867,189585112.048,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.001124129,22239441.834,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.001971909,12678069.329,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000308739,80974539.224,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.001112538,22471144.090,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000089358,279773267.672,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000972455,25708126.895,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000102432,244064404.853,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000938702,26632516.200,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000090239,277041758.219,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000925076,27024800.940,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.008201093,3048374.181,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000785543,31825114.927,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.005971311,4186685.214,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000281148,88921171.889,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000940644,26577537.670,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000087354,286192254.426,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000715633,34934111.048,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000024075,1038424990.087,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000693902,36028147.097,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000037831,660837399.835,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000718138,34812253.099,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000027071,923497943.562,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000709422,35239959.251,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,COMPARE,25000,0.000153768,162582722.719,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,COMPARE,25000,0.010999873,2272753.503,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,COMPARE,25000,0.008600373,2906850.650,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,brainpoolP512r1,512,COMPARE,25000,0.000058129,430078195.320,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,brainpoolP512r1,512,COMPARE,25000,0.000710152,35203726.591,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,brainpoolP512r1,512,COMPARE,25000,0.000033453,747319945.156,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,brainpoolP512r1,512,COMPARE,25000,0.000686989,36390675.782,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,brainpoolP512r1,512,COMPARE,25000,0.000018445,1355375246.462,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,brainpoolP512r1,512,COMPARE,25000,0.000697248,35855252.338,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,brainpoolP512r1,512,COMPARE,25000,0.000016901,1479186973.412,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,brainpoolP512r1,512,COMPARE,25000,0.000698681,35781709.217,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,REDUCE,3125,0.000094458,33083488.541,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,REDUCE,3125,0.007999605,390644.298,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,REDUCE,3125,0.008358869,373854.405,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,brainpoolP512r1,512,REDUCE,3125,0.000122811,25445614.434,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,brainpoolP512r1,512,REDUCE,3125,0.000280527,11139750.376,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,brainpoolP512r1,512,REDUCE,3125,0.000100930,30962056.421,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,brainpoolP512r1,512,REDUCE,3125,0.000241894,12918879.575,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,brainpoolP512r1,512,REDUCE,3125,0.000050915,61376596.747,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,brainpoolP512r1,512,REDUCE,3125,0.000225944,13830862.037,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,brainpoolP512r1,512,REDUCE,3125,0.000046327,67455246.692,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,brainpoolP512r1,512,REDUCE,3125,0.000218040,14332227.555,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL,1562,0.000268644,5814387.304,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL,1562,0.006572198,237667.828,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL,1562,0.008999511,173564.983,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,brainpoolP512r1,512,MODMUL,1562,0.000583414,2677344.101,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,brainpoolP512r1,512,MODMUL,1562,0.000714370,2186541.836,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,brainpoolP512r1,512,MODMUL,1562,0.000413897,3773885.941,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,brainpoolP512r1,512,MODMUL,1562,0.000515688,3028962.011,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,brainpoolP512r1,512,MODMUL,1562,0.000217158,7192923.359,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,brainpoolP512r1,512,MODMUL,1562,0.000342503,4560541.465,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,brainpoolP512r1,512,MODMUL,1562,0.000160492,9732570.805,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,brainpoolP512r1,512,MODMUL,1562,0.000274214,5696276.954,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MODEXP,390,0.016579058,23523.653,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MODEXP,390,0.009002586,43320.886,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MODEXP,390,0.008441203,46201.945,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,brainpoolP512r1,512,MODEXP,390,0.114321433,3411.434,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,brainpoolP512r1,512,MODEXP,390,0.114231854,3414.109,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,brainpoolP512r1,512,MODEXP,390,0.011254870,34651.667,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,brainpoolP512r1,512,MODEXP,390,0.011302168,34506.654,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,brainpoolP512r1,512,MODEXP,390,0.010099202,38616.913,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,brainpoolP512r1,512,MODEXP,390,0.010108809,38580.213,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,brainpoolP512r1,512,MODEXP,390,0.004262259,91500.775,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,brainpoolP512r1,512,MODEXP,390,0.004268070,91376.196,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,EXPONENTIATION,390,0.003728844,104590.055,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,EXPONENTIATION,390,0.008317511,46889.027,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,EXPONENTIATION,390,0.011453824,34049.763,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,brainpoolP512r1,512,EXPONENTIATION,390,0.084640788,4607.708,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,brainpoolP512r1,512,EXPONENTIATION,390,0.084668811,4606.183,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,brainpoolP512r1,512,EXPONENTIATION,390,0.021462974,18170.828,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,brainpoolP512r1,512,EXPONENTIATION,390,0.021358238,18259.933,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,brainpoolP512r1,512,EXPONENTIATION,390,0.007535561,51754.607,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,brainpoolP512r1,512,EXPONENTIATION,390,0.007765753,50220.500,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,brainpoolP512r1,512,EXPONENTIATION,390,0.005866649,66477.474,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,brainpoolP512r1,512,EXPONENTIATION,390,0.005768105,67613.192,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,DIVIDE,3125,0.000129072,24211293.744,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,DIVIDE,3125,0.006021305,518990.486,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,DIVIDE,3125,0.008583861,364055.292,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,brainpoolP512r1,512,DIVIDE,3125,0.000966834,3232199.436,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,brainpoolP512r1,512,DIVIDE,3125,0.001182819,2641993.367,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,brainpoolP512r1,512,DIVIDE,3125,0.000963908,3242010.129,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,brainpoolP512r1,512,DIVIDE,3125,0.001167811,2675946.735,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,brainpoolP512r1,512,DIVIDE,3125,0.000339888,9194204.169,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,brainpoolP512r1,512,DIVIDE,3125,0.000552868,5652343.728,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,brainpoolP512r1,512,DIVIDE,3125,0.000294372,10615820.387,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,brainpoolP512r1,512,DIVIDE,3125,0.000494699,6316972.621,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,ISQRT,781,0.000109105,7158248.275,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,ISQRT,781,0.006448065,121121.609,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,brainpoolP512r1,512,ISQRT,781,0.007355617,106177.359,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,brainpoolP512r1,512,ISQRT,781,0.007381506,105804.967,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,brainpoolP512r1,512,ISQRT,781,0.005905431,132251.143,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,brainpoolP512r1,512,ISQRT,781,0.005930529,131691.457,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,brainpoolP512r1,512,ISQRT,781,0.002016343,387334.904,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,brainpoolP512r1,512,ISQRT,781,0.002013158,387947.682,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,brainpoolP512r1,512,ISQRT,781,0.001893573,412447.796,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,brainpoolP512r1,512,ISQRT,781,0.001923980,405929.384,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL_R2,25000,0.003951872,6326115.825,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL_R2,25000,0.008117557,3079744.298,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL_R2,25000,0.008748360,2857678.475,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,brainpoolP512r1,512,MODMUL_R2,25000,0.000415900,60110600.092,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,brainpoolP512r1,512,MODMUL_R2,25000,0.001074566,23265205.734,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,brainpoolP512r1,512,MODMUL_R2,25000,0.000081623,306286319.345,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,brainpoolP512r1,512,MODMUL_R2,25000,0.000764835,32686791.757,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,brainpoolP512r1,512,MODMUL_R2,25000,0.000053710,465462051.387,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,brainpoolP512r1,512,MODMUL_R2,25000,0.000748745,33389207.593,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,brainpoolP512r1,512,MODMUL_R2,25000,0.000039174,638178568.924,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,brainpoolP512r1,512,MODMUL_R2,25000,0.000716465,34893536.904,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,p1024,1024,ADD,12500,0.000215695,57952199.205,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,p1024,1024,ADD,12500,0.008291753,1507521.990,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,p1024,1024,ADD,12500,0.006387982,1956799.487,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p1024,1024,ADD,12500,0.000086702,144171702.790,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p1024,1024,ADD,12500,0.000742953,16824756.452,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p1024,1024,ADD,12500,0.000045696,273546675.634,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p1024,1024,ADD,12500,0.000724950,17242565.912,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32,p1024,1024,ADD,12500,0.000026410,473305926.122,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32,p1024,1024,ADD,12500,0.000679485,18396284.496,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p1024,1024,ADD,12500,0.000025798,484531788.235,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p1024,1024,ADD,12500,0.000721013,17336715.139,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p1024,1024,ADD,12500,0.000026620,469571871.392,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p1024,1024,ADD,12500,0.000701326,17823381.666,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,p1024,1024,SUBTRACT,12500,0.000171942,72698959.627,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,p1024,1024,SUBTRACT,12500,0.008476389,1474684.544,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,p1024,1024,SUBTRACT,12500,0.005936906,2105473.803,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p1024,1024,SUBTRACT,12500,0.000081463,153443860.055,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p1024,1024,SUBTRACT,12500,0.000738595,16924021.307,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p1024,1024,SUBTRACT,12500,0.000045495,274756222.908,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p1024,1024,SUBTRACT,12500,0.000719941,17362534.155,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32,p1024,1024,SUBTRACT,12500,0.000025759,485267559.701,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32,p1024,1024,SUBTRACT,12500,0.000706396,17695455.295,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p1024,1024,SUBTRACT,12500,0.000026119,478579882.332,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p1024,1024,SUBTRACT,12500,0.000718208,17404433.719,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p1024,1024,SUBTRACT,12500,0.000026600,469925347.058,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p1024,1024,SUBTRACT,12500,0.000696507,17946695.593,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,p1024,1024,ADDMOD,12500,0.000704221,17750110.659,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,p1024,1024,ADDMOD,12500,0.008553504,1461389.389,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,p1024,1024,ADDMOD,12500,0.005841297,2139935.667,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p1024,1024,ADDMOD,12500,0.000125034,99972796.392,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p1024,1024,ADDMOD,12500,0.000795432,15714733.075,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p1024,1024,ADDMOD,12500,0.000067837,184265032.022,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p1024,1024,ADDMOD,12500,0.000745960,16756929.642,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32,p1024,1024,ADDMOD,12500,0.000034314,364284055.178,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32,p1024,1024,ADDMOD,12500,0.000744197,16796626.602,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p1024,1024,ADDMOD,12500,0.000035386,353246379.176,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p1024,1024,ADDMOD,12500,0.000732795,17057974.265,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p1024,1024,ADDMOD,12500,0.000035857,348606156.943,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p1024,1024,ADDMOD,12500,0.000726694,17201187.648,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,p1024,1024,SUBTRACTMOD,12500,0.000588023,21257668.142,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,p1024,1024,SUBTRACTMOD,12500,0.008121614,1539102.921,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,p1024,1024,SUBTRACTMOD,12500,0.008832557,1415218.736,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p1024,1024,SUBTRACTMOD,12500,0.000123922,100869890.144,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p1024,1024,SUBTRACTMOD,12500,0.000775725,16113955.520,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p1024,1024,SUBTRACTMOD,12500,0.000066064,189210942.335,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p1024,1024,SUBTRACTMOD,12500,0.000747612,16719903.306,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32,p1024,1024,SUBTRACTMOD,12500,0.000034915,358011797.891,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32,p1024,1024,SUBTRACTMOD,12500,0.000714451,17495949.610,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p1024,1024,SUBTRACTMOD,12500,0.000035647,350659955.716,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p1024,1024,SUBTRACTMOD,12500,0.000705844,17709300.773,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p1024,1024,SUBTRACTMOD,12500,0.000035116,355963262.654,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p1024,1024,SUBTRACTMOD,12500,0.000709281,17623478.937,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.001803112,6932459.057,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.008713074,1434625.708,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.008060419,1550787.868,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.011645444,1073381.125,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.012454272,1003671.671,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.002969883,4208919.803,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.003783620,3303714.465,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.000723298,17281951.598,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.001565036,7987036.051,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.000084388,148125612.027,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.000943761,13244877.527,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.000084349,148193894.158,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.000925026,13513135.106,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.001802541,6934654.705,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.001606083,7782910.302,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.007001643,1785295.234,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.004205658,2972186.445,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.005036868,2481700.968,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.001085668,11513650.108,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.001889896,6614120.168,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.000283823,44041514.966,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.001122767,11133208.246,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.000286928,43564914.103,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.001112568,11235266.375,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.000286979,43557208.910,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.001173392,10652874.894,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.012128671,1030615.809,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.008819914,1417247.386,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.008614248,1451084.275,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.000744768,16783751.173,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.001406418,8887827.553,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.000151284,82626032.997,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.000780525,16014863.458,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.000049443,252816455.386,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.000719721,17367836.415,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.000066636,187585923.131,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.000738035,16936861.801,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.000055284,226105175.558,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.000734478,17018889.550,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,p1024,1024,COMPARE,12500,0.000076764,162836681.943,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,p1024,1024,COMPARE,12500,0.008151730,1533416.822,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,p1024,1024,COMPARE,12500,0.000525305,23795697.484,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p1024,1024,COMPARE,12500,0.000076714,162942938.656,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p1024,1024,COMPARE,12500,0.000744226,16795975.005,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p1024,1024,COMPARE,12500,0.000042320,295368645.984,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p1024,1024,COMPARE,12500,0.000688622,18152193.613,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p1024,1024,COMPARE,12500,0.000025267,494716148.948,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p1024,1024,COMPARE,12500,0.000717056,17432390.611,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p1024,1024,COMPARE,12500,0.000025278,494497427.442,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p1024,1024,COMPARE,12500,0.000699864,17860612.869,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,p1024,1024,REDUCE,1562,0.000031018,50357968.461,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,p1024,1024,REDUCE,1562,0.006002790,260212.342,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,p1024,1024,REDUCE,1562,0.008600332,181620.899,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p1024,1024,REDUCE,1562,0.001054047,1481907.482,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p1024,1024,REDUCE,1562,0.001222233,1277988.673,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p1024,1024,REDUCE,1562,0.000587242,2659891.688,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p1024,1024,REDUCE,1562,0.000727244,2147834.974,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p1024,1024,REDUCE,1562,0.000156844,9958937.884,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p1024,1024,REDUCE,1562,0.000348524,4481757.898,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p1024,1024,REDUCE,1562,0.000158056,9882563.816,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p1024,1024,REDUCE,1562,0.000332604,4696276.445,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,p1024,1024,MODMUL,781,0.000388680,2009364.896,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,p1024,1024,MODMUL,781,0.005728816,136328.346,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,p1024,1024,MODMUL,781,0.005999784,130171.353,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p1024,1024,MODMUL,781,0.003918471,199312.433,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p1024,1024,MODMUL,781,0.004015292,194506.400,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p1024,1024,MODMUL,781,0.001020385,765397.410,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p1024,1024,MODMUL,781,0.001128016,692366.184,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p1024,1024,MODMUL,781,0.000557145,1401789.219,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p1024,1024,MODMUL,781,0.000704882,1107986.919,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p1024,1024,MODMUL,781,0.000457840,1705836.817,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p1024,1024,MODMUL,781,0.000580600,1345160.186,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,p1024,1024,MODEXP,195,0.056158981,3472.285,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,p1024,1024,MODEXP,195,0.006000746,32495.959,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,p1024,1024,MODEXP,195,0.011998556,16251.956,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p1024,1024,MODEXP,195,0.972214674,200.573,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p1024,1024,MODEXP,195,0.973610364,200.285,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p1024,1024,MODEXP,195,0.153017513,1274.364,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p1024,1024,MODEXP,195,0.153096842,1273.704,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p1024,1024,MODEXP,195,0.067395601,2893.364,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p1024,1024,MODEXP,195,0.067571772,2885.820,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p1024,1024,MODEXP,195,0.041174035,4735.994,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p1024,1024,MODEXP,195,0.041412695,4708.701,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,p1024,1024,EXPONENTIATION,195,0.007804479,24985.653,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,p1024,1024,EXPONENTIATION,195,0.008927526,21842.557,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,p1024,1024,EXPONENTIATION,195,0.016876607,11554.455,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p1024,1024,EXPONENTIATION,195,0.710036739,274.634,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p1024,1024,EXPONENTIATION,195,0.711554132,274.048,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p1024,1024,EXPONENTIATION,195,0.166622542,1170.310,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p1024,1024,EXPONENTIATION,195,0.166826676,1168.878,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p1024,1024,EXPONENTIATION,195,0.052496623,3714.525,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p1024,1024,EXPONENTIATION,195,0.052678243,3701.718,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p1024,1024,EXPONENTIATION,195,0.048588691,4013.280,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p1024,1024,EXPONENTIATION,195,0.048683027,4005.503,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,p1024,1024,DIVIDE,1562,0.000068058,22951003.282,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,p1024,1024,DIVIDE,1562,0.008860761,176282.826,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,p1024,1024,DIVIDE,1562,0.008998639,173581.801,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p1024,1024,DIVIDE,1562,0.008584173,181962.783,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p1024,1024,DIVIDE,1562,0.008648673,180605.740,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p1024,1024,DIVIDE,1562,0.008511832,183509.261,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p1024,1024,DIVIDE,1562,0.008721717,179093.173,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p1024,1024,DIVIDE,1562,0.001537014,1016256.174,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p1024,1024,DIVIDE,1562,0.001758190,888413.600,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p1024,1024,DIVIDE,1562,0.001507088,1036435.776,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p1024,1024,DIVIDE,1562,0.001724937,905540.319,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,p1024,1024,ISQRT,390,0.000106179,3673045.370,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,p1024,1024,ISQRT,390,0.006684890,58340.526,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p1024,1024,ISQRT,390,0.109444872,3563.438,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p1024,1024,ISQRT,390,0.108538072,3593.209,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p1024,1024,ISQRT,390,0.064616651,6035.596,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p1024,1024,ISQRT,390,0.064643142,6033.123,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p1024,1024,ISQRT,390,0.012487223,31231.924,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p1024,1024,ISQRT,390,0.012500458,31198.856,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p1024,1024,ISQRT,390,0.012066083,32322.005,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p1024,1024,ISQRT,390,0.012094156,32246.980,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,p1024,1024,MODMUL_R2,12500,0.005492944,2275646.705,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,p1024,1024,MODMUL_R2,12500,0.006000716,2083084.750,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,p1024,1024,MODMUL_R2,12500,0.008998549,1389112.833,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p1024,1024,MODMUL_R2,12500,0.001209649,10333575.894,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p1024,1024,MODMUL_R2,12500,0.001891448,6608693.691,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p1024,1024,MODMUL_R2,12500,0.000178044,70207378.909,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p1024,1024,MODMUL_R2,12500,0.000865033,14450313.610,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p1024,1024,MODMUL_R2,12500,0.000101461,123200238.657,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p1024,1024,MODMUL_R2,12500,0.000771889,16194038.316,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p1024,1024,MODMUL_R2,12500,0.000084669,147633558.549,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p1024,1024,MODMUL_R2,12500,0.000774132,16147115.197,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,p2048,2048,ADD,6250,0.000151664,41209510.709,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,p2048,2048,ADD,6250,0.005999954,1041674.661,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,p2048,2048,ADD,6250,0.004388982,1424020.434,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p2048,2048,ADD,6250,0.000148989,41949397.956,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p2048,2048,ADD,6250,0.000855485,7305797.490,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p2048,2048,ADD,6250,0.000077215,80942797.525,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p2048,2048,ADD,6250,0.000761098,8211820.339,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32,p2048,2048,ADD,6250,0.000040025,156152464.719,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32,p2048,2048,ADD,6250,0.000711214,8787790.644,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p2048,2048,ADD,6250,0.000039975,156348006.593,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p2048,2048,ADD,6250,0.000700835,8917934.807,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p2048,2048,ADD,6250,0.000040276,155179355.316,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p2048,2048,ADD,6250,0.001080277,5785552.384,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,p2048,2048,SUBTRACT,6250,0.000122419,51054127.725,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,p2048,2048,SUBTRACT,6250,0.008221291,760221.235,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,p2048,2048,SUBTRACT,6250,0.006003361,1041083.493,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p2048,2048,SUBTRACT,6250,0.000148889,41977540.291,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p2048,2048,SUBTRACT,6250,0.000858621,7279114.061,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p2048,2048,SUBTRACT,6250,0.000077456,80690969.757,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p2048,2048,SUBTRACT,6250,0.000742964,8412251.685,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32,p2048,2048,SUBTRACT,6250,0.000039755,157212399.632,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32,p2048,2048,SUBTRACT,6250,0.000711295,8786792.485,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p2048,2048,SUBTRACT,6250,0.000039975,156348006.593,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p2048,2048,SUBTRACT,6250,0.000704001,8877829.095,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p2048,2048,SUBTRACT,6250,0.000039734,157296231.014,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p2048,2048,SUBTRACT,6250,0.001116945,5595620.484,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,p2048,2048,ADDMOD,6250,0.000432411,14453842.173,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,p2048,2048,ADDMOD,6250,0.006000045,1041658.856,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,p2048,2048,ADDMOD,6250,0.000839315,7446549.657,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p2048,2048,ADDMOD,6250,0.000200626,31152496.629,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p2048,2048,ADDMOD,6250,0.000900730,6938816.101,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p2048,2048,ADDMOD,6250,0.000102632,60897197.602,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p2048,2048,ADDMOD,6250,0.000782408,7988159.089,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32,p2048,2048,ADDMOD,6250,0.000049974,125064856.479,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32,p2048,2048,ADDMOD,6250,0.000719360,8688280.365,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p2048,2048,ADDMOD,6250,0.000043512,143639010.713,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p2048,2048,ADDMOD,6250,0.000726564,8602134.483,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p2048,2048,ADDMOD,6250,0.000043381,144072271.361,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p2048,2048,ADDMOD,6250,0.000721544,8661980.092,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,p2048,2048,SUBTRACTMOD,6250,0.000356699,17521770.797,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,p2048,2048,SUBTRACTMOD,6250,0.004423347,1412957.254,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,p2048,2048,SUBTRACTMOD,6250,0.000696036,8979420.175,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p2048,2048,SUBTRACTMOD,6250,0.000205285,30445490.716,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p2048,2048,SUBTRACTMOD,6250,0.000905860,6899521.030,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p2048,2048,SUBTRACTMOD,6250,0.000106059,58929455.567,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p2048,2048,SUBTRACTMOD,6250,0.000753785,8291488.756,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32,p2048,2048,SUBTRACTMOD,6250,0.000052489,119072323.777,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32,p2048,2048,SUBTRACTMOD,6250,0.000727385,8592422.958,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p2048,2048,SUBTRACTMOD,6250,0.000044022,141974505.085,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p2048,2048,SUBTRACTMOD,6250,0.000749707,8336590.282,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p2048,2048,SUBTRACTMOD,6250,0.000042970,145450901.091,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p2048,2048,SUBTRACTMOD,6250,0.000724179,8630463.247,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.003042627,2054145.888,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.008178631,764186.566,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.008997166,694663.181,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.043102765,145002.299,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.044116707,141669.685,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.010969674,569752.575,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.011816733,528910.992,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.002547050,2453819.001,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.003767530,1658911.821,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.000151073,41370702.192,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.001018541,6136227.190,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.000128711,48558445.322,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.000971653,6432337.251,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.003090527,2022308.876,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.007967475,784439.249,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.002202832,2837256.887,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.016650644,375360.854,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.017578235,355553.333,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.004300670,1453261.934,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.005159361,1211390.312,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.001216242,5138780.551,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.002060105,3033825.809,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.001104433,5659011.903,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.001965569,3179741.044,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.001099052,5686719.578,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.001955009,3196916.122,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.019990397,312650.120,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.004239603,1474194.655,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.006004934,1040810.779,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.021573579,289706.220,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.022377918,279293.185,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.000718959,8693125.469,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.001398765,4468227.277,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.000125285,49886257.259,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.000771969,8096181.514,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.000175109,35692074.986,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.000865915,7217797.749,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.000143549,43539087.714,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.000808798,7727516.545,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,p2048,2048,COMPARE,6250,0.000037371,167242211.243,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,p2048,2048,COMPARE,6250,0.000243477,25669794.258,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,p2048,2048,COMPARE,6250,0.001099682,5683460.287,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p2048,2048,COMPARE,6250,0.000090049,69406568.448,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p2048,2048,COMPARE,6250,0.000791124,7900151.507,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p2048,2048,COMPARE,6250,0.000048501,128863451.587,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p2048,2048,COMPARE,6250,0.000726112,8607488.359,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p2048,2048,COMPARE,6250,0.000027281,229097179.336,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p2048,2048,COMPARE,6250,0.000739989,8446070.311,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p2048,2048,COMPARE,6250,0.000026850,232774415.539,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p2048,2048,COMPARE,6250,0.000703581,8883129.016,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,p2048,2048,REDUCE,781,0.000024205,32265962.468,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,p2048,2048,REDUCE,781,0.003884887,201035.444,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,p2048,2048,REDUCE,781,0.008170796,95584.324,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p2048,2048,REDUCE,781,0.031334219,24924.827,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p2048,2048,REDUCE,781,0.030196865,25863.612,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p2048,2048,REDUCE,781,0.001909453,409017.621,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p2048,2048,REDUCE,781,0.002071187,377078.434,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p2048,2048,REDUCE,781,0.000327615,2383895.241,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p2048,2048,REDUCE,781,0.000513053,1522259.784,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p2048,2048,REDUCE,781,0.000340689,2292411.926,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p2048,2048,REDUCE,781,0.000513965,1519558.633,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,p2048,2048,MODMUL,390,0.000539522,722861.885,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,p2048,2048,MODMUL,390,0.005984926,65163.712,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,p2048,2048,MODMUL,390,0.004041541,96497.848,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p2048,2048,MODMUL,390,0.048824488,7987.795,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p2048,2048,MODMUL,390,0.048851248,7983.419,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p2048,2048,MODMUL,390,0.007089584,55010.281,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p2048,2048,MODMUL,390,0.007208878,54099.957,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p2048,2048,MODMUL,390,0.002806867,138944.951,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p2048,2048,MODMUL,390,0.002949334,132233.244,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p2048,2048,MODMUL,390,0.002172717,179498.746,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p2048,2048,MODMUL,390,0.002317148,168310.352,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,p2048,2048,MODEXP,97,0.193704224,500.763,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,p2048,2048,MODEXP,97,0.020001248,4849.697,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,p2048,2048,MODEXP,97,0.015106656,6421.011,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p2048,2048,MODEXP,97,15.475319768,6.268,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p2048,2048,MODEXP,97,15.623051427,6.209,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p2048,2048,MODEXP,97,1.566673162,61.915,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p2048,2048,MODEXP,97,1.567334934,61.888,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p2048,2048,MODEXP,97,0.395226837,245.429,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p2048,2048,MODEXP,97,0.395329510,245.365,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p2048,2048,MODEXP,97,0.461870485,210.016,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p2048,2048,MODEXP,97,0.460972050,210.425,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,p2048,2048,EXPONENTIATION,97,0.023205558,4180.033,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,p2048,2048,EXPONENTIATION,97,0.008998649,10779.396,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,p2048,2048,EXPONENTIATION,97,0.040524565,2393.610,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p2048,2048,EXPONENTIATION,97,6.596189985,14.705,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p2048,2048,EXPONENTIATION,97,6.590429238,14.718,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p2048,2048,EXPONENTIATION,97,1.343384446,72.206,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p2048,2048,EXPONENTIATION,97,1.347497881,71.985,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p2048,2048,EXPONENTIATION,97,0.357910172,271.018,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p2048,2048,EXPONENTIATION,97,0.359078235,270.136,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p2048,2048,EXPONENTIATION,97,0.370353700,261.912,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p2048,2048,EXPONENTIATION,97,0.370657289,261.697,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,p2048,2048,DIVIDE,781,0.000043632,17899707.885,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,p2048,2048,DIVIDE,781,0.008088953,96551.431,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,p2048,2048,DIVIDE,781,0.005174537,150931.379,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p2048,2048,DIVIDE,781,0.156823309,4980.127,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p2048,2048,DIVIDE,781,0.161390990,4839.180,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p2048,2048,DIVIDE,781,0.123446030,6326.651,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p2048,2048,DIVIDE,781,0.123739100,6311.667,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p2048,2048,DIVIDE,781,0.015842740,49297.028,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p2048,2048,DIVIDE,781,0.016089012,48542.446,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p2048,2048,DIVIDE,781,0.011382641,68613.250,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p2048,2048,DIVIDE,781,0.011581984,67432.316,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,p2048,2048,ISQRT,195,0.000082756,2356326.696,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,p2048,2048,ISQRT,195,0.007184576,27141.476,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p2048,2048,ISQRT,195,1.501546080,129.866,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p2048,2048,ISQRT,195,1.502050861,129.823,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p2048,2048,ISQRT,195,1.041814114,187.174,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p2048,2048,ISQRT,195,1.040898294,187.338,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p2048,2048,ISQRT,195,0.078921562,2470.808,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p2048,2048,ISQRT,195,0.078940597,2470.212,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p2048,2048,ISQRT,195,0.027416860,7112.412,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p2048,2048,ISQRT,195,0.027392895,7118.634,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-1t,p2048,2048,MODMUL_R2,6250,0.008599211,726810.847,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,gmp-nt,p2048,2048,MODMUL_R2,6250,0.004782701,1306792.985,0
library,AMD Ryzen 9 3900X 12-Core Processor,host-cpu,openssl-nt,p2048,2048,MODMUL_R2,6250,0.001342508,4655465.704,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p2048,2048,MODMUL_R2,6250,0.012990476,481121.713,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w8,p2048,2048,MODMUL_R2,6250,0.013823940,452114.225,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p2048,2048,MODMUL_R2,6250,0.000758784,8236864.366,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w16,p2048,2048,MODMUL_R2,6250,0.001454028,4298404.644,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p2048,2048,MODMUL_R2,6250,0.000314199,19891858.309,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-opt,p2048,2048,MODMUL_R2,6250,0.000998814,6257420.521,0
opencl-kernel,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p2048,2048,MODMUL_R2,6250,0.000237475,26318569.102,0
opencl-e2e,NVIDIA GeForce RTX 3080 Ti,GPU,w32-o64,p2048,2048,MODMUL_R2,6250,0.000919335,6798393.623,0
```
