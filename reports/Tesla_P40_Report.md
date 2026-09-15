# MPA-OpenCL benchmark report - Tesla P40


> **Note.** Two column groups have been removed from this report: the CGBN
> reference column, which was invalid, and the multi-threaded GMP and OpenSSL
> baselines, which predate the 2026-09-12 timing fix and were understated.
> See `reports/README.md`. The single-threaded GMP column, the OpenCL-on-CPU
> rows and every MPA measurement are unaffected, and every configuration was
> verified word-for-word against GMP before it was timed.


> **Partial report.** The run was interrupted or hit its time budget.
> Rows that never ran are marked `n/a`.

## 1. System under test

2 OpenCL device(s) exercised with the identical kernels and operands.

### Device 0 - Tesla P40 (GPU)

| Property | Value |
|---|---|
| Model | Tesla P40 |
| Type | GPU |
| Vendor | NVIDIA Corporation |
| Device memory | 23.87 GiB |
| Max single allocation | 5.97 GiB |
| Local memory | 48 KiB |
| Global cache | 1440 KiB |
| Compute units | 30 |
| Max clock | 1531 MHz |
| Max work-group size | 1024 |
| OpenCL version | OpenCL 3.0 CUDA |
| Driver | 580.126.09 |

### Device 1 - cpu-skylake-avx512-Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz (CPU)

| Property | Value |
|---|---|
| Model | cpu-skylake-avx512-Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz |
| Type | CPU |
| Vendor | GenuineIntel |
| Device memory | 249.26 GiB |
| Max single allocation | 64.00 GiB |
| Local memory | 1024 KiB |
| Global cache | 25344 KiB |
| Compute units | 48 |
| Max clock | 3200 MHz |
| Max work-group size | 4096 |
| OpenCL version | OpenCL 3.0 PoCL HSTR: cpu-x86_64-pc-linux-gnu-skylake-avx512 |
| Driver | 5.0+debian |

### Host

| Property | Value |
|---|---|
| CPU | Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz |
| Logical cores | 48 |
| OpenMP threads used | 48 |
| RAM | 251.3 GB |
| OS | Ubuntu 24.04.4 LTS |
| Kernel | 5.15.0-171-generic |
| Arch | x86_64 |
| GMP | 6.3.0 |
| OpenSSL | OpenSSL 3.0.13 30 Jan 2024 |

## 2. Method

- Workload auto-sized from the device and host: --min-items from 700 x compute units, --items from ten times that capped by host RAM. Either flag, given explicitly, overrides its half.
- Base workload 50000 items, scaled down per operator by its cost weight and by modulus size. Device rows honour --min-items (21000) so the GPU is not left idle; the CPU libraries keep the smaller count because a full-width MODEXP there costs minutes. Both counts appear in every row as dev/cpu, and throughput is per-second so they remain comparable.
- 5 timed repetitions, **minimum** reported. Two untimed warm-up launches precede them.
- `kernel` times `clEnqueueNDRangeKernel` + `clFinish` only. `e2e` adds the host->device operand writes and the device->host result read.
- Every OpenCL device runs the same kernels on the same operands, so GPU and CPU-OpenCL columns are directly comparable.
- CPU library baselines (GMP, OpenSSL) run those same operands, with temporaries preallocated outside the timed region, so the figure is the arithmetic and not marshalling. The generator is reseeded per modulus and operation so every backend sees identical inputs.
- OpenSSL rows time the nearest BN primitive, which is not always semantically identical (its Montgomery routine expects Montgomery-domain inputs); they measure comparable work, not identical results. Correctness is judged against GMP only.
- Every device cell is checked word-for-word against GMP before it is timed. A cell that mismatches is reported and excluded from the speedup tables.
- Total wall time 5407.7 s.

## 3. Correctness

| Device | Kernel | Configs run | Passed | Mismatched | Build/launch failed |
|---|---|---|---|---|---|
| [0] GPU | `mpaKernels_8bits.cl` (w8) | 75 | 75 | 0 | 0 |
| [0] GPU | `mpaKernel_16bits.cl` (w16) | 74 | 74 | 0 | 0 |

**All configurations correct** - 149 configurations, 0 problems.

## 4. Throughput per device

Operations per second, higher is better. Kernel-only timings.

### Device 0 - Tesla P40 (GPU)

#### secp256k1 (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | 390.28 M | 797.36 M | - | - | - | - | - | 52.62 M |
| SUBTRACT | 50000 / 50000 | 383.62 M | 768.76 M | - | - | - | - | - | 40.61 M |
| ADDMOD | 50000 / 50000 | 286.78 M | 596.71 M | - | - | - | - | - | 16.18 M |
| SUBTRACTMOD | 50000 / 50000 | 284.12 M | 591.66 M | - | - | - | - | - | 23.94 M |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | 8.39 M | 39.34 M | - | - | - | - | - | 48.93 M |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | 15.97 M | 63.07 M | - | - | - | - | - | 49.25 M |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | 138.23 M | 463.44 M | - | - | - | - | - | 6.29 M |
| COMPARE | 50000 / 50000 | 583.30 M | 1.03 G | - | - | - | - | - | 125.73 M |
| REDUCE | 21000 / 6250 | 89.15 M | 161.39 M | - | - | - | - | - | 49.42 M |
| MODMUL | 21000 / 3125 | 31.13 M | 53.39 M | - | - | - | - | - | 8.39 M |
| MODEXP | 21000 / 781 | 715.58 k | 2.93 M | - | - | - | - | - | 95.05 k |
| EXPONENTIATION | 21000 / 781 | 135.30 k | 747.04 k | - | - | - | - | - | 255.10 k |
| DIVIDE | 21000 / 6250 | 46.05 M | 65.01 M | - | - | - | - | - | 17.17 M |
| ISQRT | 21000 / 1562 | 3.26 M | 5.30 M | - | - | - | - | - | 9.00 M |
| MODMUL_R2 | 50000 / 50000 | 128.49 M | 400.26 M | - | - | - | - | - | 11.43 M |

#### rsa256(composite) (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | 459.18 M | 799.42 M | - | - | - | - | - | 49.78 M |
| SUBTRACT | 50000 / 50000 | 458.55 M | 784.52 M | - | - | - | - | - | 70.32 M |
| ADDMOD | 50000 / 50000 | 351.96 M | 627.61 M | - | - | - | - | - | 23.11 M |
| SUBTRACTMOD | 50000 / 50000 | 324.91 M | 589.09 M | - | - | - | - | - | 23.27 M |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | 9.90 M | 39.35 M | - | - | - | - | - | 48.43 M |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | 16.00 M | 63.39 M | - | - | - | - | - | 51.36 M |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | 138.16 M | 459.69 M | - | - | - | - | - | 6.36 M |
| COMPARE | 50000 / 50000 | 587.93 M | 993.18 M | - | - | - | - | - | 129.43 M |
| REDUCE | 21000 / 6250 | 88.29 M | 160.74 M | - | - | - | - | - | 30.04 M |
| MODMUL | 21000 / 3125 | 31.15 M | 53.86 M | - | - | - | - | - | 4.42 M |
| MODEXP | 21000 / 781 | 717.18 k | 2.93 M | - | - | - | - | - | 110.47 k |
| EXPONENTIATION | 21000 / 781 | 135.32 k | 746.79 k | - | - | - | - | - | 261.71 k |
| DIVIDE | 21000 / 6250 | 46.59 M | 75.29 M | - | - | - | - | - | 9.56 M |
| ISQRT | 21000 / 1562 | 3.26 M | 5.31 M | - | - | - | - | - | 10.11 M |
| MODMUL_R2 | 50000 / 50000 | 129.15 M | 396.18 M | - | - | - | - | - | 11.05 M |

#### brainpoolP512r1 (512-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 25000 / 25000 | 204.46 M | 371.28 M | - | - | - | - | - | 48.43 M |
| SUBTRACT | 25000 / 25000 | 205.98 M | 377.45 M | - | - | - | - | - | 65.92 M |
| ADDMOD | 25000 / 25000 | 146.90 M | 284.07 M | - | - | - | - | - | 20.95 M |
| SUBTRACTMOD | 25000 / 25000 | 134.15 M | 267.73 M | - | - | - | - | - | 23.62 M |
| MULTIPLYOPERANDSCANNING | 25000 / 25000 | 2.24 M | 8.67 M | - | - | - | - | - | 22.83 M |
| MULTIPLYPRODUCTSCANNING | 25000 / 25000 | 4.06 M | 16.18 M | - | - | - | - | - | 21.89 M |
| MONTGOMERYMULTIPLICATION | 25000 / 25000 | 25.89 M | 113.01 M | - | - | - | - | - | 2.86 M |
| COMPARE | 25000 / 25000 | 279.30 M | 493.55 M | - | - | - | - | - | 61.85 M |
| REDUCE | 21000 / 3125 | 27.66 M | 48.72 M | - | - | - | - | - | 29.30 M |
| MODMUL | 21000 / 1562 | 8.67 M | 16.81 M | - | - | - | - | - | 4.71 M |
| MODEXP | 21000 / 390 | 40.63 k | 425.59 k | - | - | - | - | - | 19.35 k |
| EXPONENTIATION | 21000 / 390 | 16.62 k | 68.04 k | - | - | - | - | - | 68.65 k |
| DIVIDE | 21000 / 3125 | 11.11 M | 8.69 M | - | - | - | - | - | 15.72 M |
| ISQRT | 21000 / 781 | 437.53 k | 753.33 k | - | - | - | - | - | 5.88 M |
| MODMUL_R2 | 25000 / 25000 | 25.10 M | 114.85 M | - | - | - | - | - | 5.65 M |

#### p1024 (1024-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 21000 / 12500 | 75.79 M | 144.15 M | - | - | - | - | - | 40.19 M |
| SUBTRACT | 21000 / 12500 | 75.22 M | 143.52 M | - | - | - | - | - | 52.20 M |
| ADDMOD | 21000 / 12500 | 52.19 M | 100.35 M | - | - | - | - | - | 5.83 M |
| SUBTRACTMOD | 21000 / 12500 | 52.10 M | 100.84 M | - | - | - | - | - | 19.38 M |
| MULTIPLYOPERANDSCANNING | 21000 / 12500 | 240.82 k | 797.64 k | - | - | - | - | - | 6.30 M |
| MULTIPLYPRODUCTSCANNING | 21000 / 12500 | 1.05 M | 4.16 M | - | - | - | - | - | 6.28 M |
| MONTGOMERYMULTIPLICATION | 21000 / 12500 | 4.92 M | 25.81 M | - | - | - | - | - | 908.01 k |
| COMPARE | 21000 / 12500 | 120.32 M | 229.91 M | - | - | - | - | - | 54.51 M |
| REDUCE | 21000 / 1562 | 4.47 M | 13.67 M | - | - | - | - | - | 40.52 M |
| MODMUL | 21000 / 781 | 1.24 M | 3.85 M | - | - | - | - | - | 1.92 M |
| MODEXP | 21000 / 195 | 4.76 k | 34.74 k | - | - | - | - | - | 2.37 k |
| EXPONENTIATION | 21000 / 195 | 2.02 k | 8.45 k | - | - | - | - | - | 21.93 k |
| DIVIDE | 21000 / 1562 | 159.26 k | 668.89 k | - | - | - | - | - | 10.38 M |
| ISQRT | 21000 / 390 | 14.49 k | 69.18 k | - | - | - | - | - | 3.12 M |
| MODMUL_R2 | 21000 / 12500 | 4.49 M | 33.49 M | - | - | - | - | - | 2.06 M |

#### p2048 (2048-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 21000 / 6250 | 34.51 M | 67.09 M | - | - | - | - | - | 28.38 M |
| SUBTRACT | 21000 / 6250 | 34.61 M | 68.22 M | - | - | - | - | - | 30.49 M |
| ADDMOD | 21000 / 6250 | 25.14 M | 50.33 M | - | - | - | - | - | 11.83 M |
| SUBTRACTMOD | 21000 / 6250 | 24.67 M | 47.93 M | - | - | - | - | - | 14.86 M |
| MULTIPLYOPERANDSCANNING | 21000 / 6250 | 43.81 k | 141.88 k | - | - | - | - | - | 2.01 M |
| MULTIPLYPRODUCTSCANNING | 21000 / 6250 | 265.48 k | 1.05 M | - | - | - | - | - | 1.84 M |
| MONTGOMERYMULTIPLICATION | 21000 / 6250 | 187.09 k | 5.15 M | - | - | - | - | - | 292.06 k |
| COMPARE | 21000 / 6250 | 63.01 M | 120.93 M | - | - | - | - | - | 75.74 M |
| REDUCE | 21000 / 781 | 28.17 k | 2.43 M | - | - | - | - | - | 29.07 M |
| MODMUL | 21000 / 390 | 11.94 k | 651.93 k | - | - | - | - | - | 631.01 k |
| MODEXP | 21000 / 97 | 120.7 | 415.0 | - | - | - | - | - | 366.0 |
| EXPONENTIATION | 21000 / 97 | 228.2 | 1.03 k | - | - | - | - | - | 3.88 k |
| DIVIDE | 21000 / 781 | 11.07 k | 32.62 k | - | - | - | - | - | 13.06 M |
| ISQRT | 21000 / 195 | 795.0 | over budget | - | - | - | - | - | 1.18 M |
| MODMUL_R2 | 21000 / 6250 | 307.53 k | - | - | - | - | - | - | 678.90 k |

### Device 1 - cpu-skylake-avx512-Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz (CPU)

#### secp256k1 (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | - | - | - | - | - | - | - | 52.62 M |
| SUBTRACT | 50000 / 50000 | - | - | - | - | - | - | - | 40.61 M |
| ADDMOD | 50000 / 50000 | - | - | - | - | - | - | - | 16.18 M |
| SUBTRACTMOD | 50000 / 50000 | - | - | - | - | - | - | - | 23.94 M |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | - | - | - | - | - | - | - | 48.93 M |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | - | - | - | - | - | - | - | 49.25 M |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | - | - | - | - | - | - | - | 6.29 M |
| COMPARE | 50000 / 50000 | - | - | - | - | - | - | - | 125.73 M |
| REDUCE | 21000 / 6250 | - | - | - | - | - | - | - | 49.42 M |
| MODMUL | 21000 / 3125 | - | - | - | - | - | - | - | 8.39 M |
| MODEXP | 21000 / 781 | - | - | - | - | - | - | - | 95.05 k |
| EXPONENTIATION | 21000 / 781 | - | - | - | - | - | - | - | 255.10 k |
| DIVIDE | 21000 / 6250 | - | - | - | - | - | - | - | 17.17 M |
| ISQRT | 21000 / 1562 | - | - | - | - | - | - | - | 9.00 M |
| MODMUL_R2 | 50000 / 50000 | - | - | - | - | - | - | - | 11.43 M |

#### rsa256(composite) (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | - | - | - | - | - | - | - | 49.78 M |
| SUBTRACT | 50000 / 50000 | - | - | - | - | - | - | - | 70.32 M |
| ADDMOD | 50000 / 50000 | - | - | - | - | - | - | - | 23.11 M |
| SUBTRACTMOD | 50000 / 50000 | - | - | - | - | - | - | - | 23.27 M |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | - | - | - | - | - | - | - | 48.43 M |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | - | - | - | - | - | - | - | 51.36 M |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | - | - | - | - | - | - | - | 6.36 M |
| COMPARE | 50000 / 50000 | - | - | - | - | - | - | - | 129.43 M |
| REDUCE | 21000 / 6250 | - | - | - | - | - | - | - | 30.04 M |
| MODMUL | 21000 / 3125 | - | - | - | - | - | - | - | 4.42 M |
| MODEXP | 21000 / 781 | - | - | - | - | - | - | - | 110.47 k |
| EXPONENTIATION | 21000 / 781 | - | - | - | - | - | - | - | 261.71 k |
| DIVIDE | 21000 / 6250 | - | - | - | - | - | - | - | 9.56 M |
| ISQRT | 21000 / 1562 | - | - | - | - | - | - | - | 10.11 M |
| MODMUL_R2 | 50000 / 50000 | - | - | - | - | - | - | - | 11.05 M |

#### brainpoolP512r1 (512-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 25000 / 25000 | - | - | - | - | - | - | - | 48.43 M |
| SUBTRACT | 25000 / 25000 | - | - | - | - | - | - | - | 65.92 M |
| ADDMOD | 25000 / 25000 | - | - | - | - | - | - | - | 20.95 M |
| SUBTRACTMOD | 25000 / 25000 | - | - | - | - | - | - | - | 23.62 M |
| MULTIPLYOPERANDSCANNING | 25000 / 25000 | - | - | - | - | - | - | - | 22.83 M |
| MULTIPLYPRODUCTSCANNING | 25000 / 25000 | - | - | - | - | - | - | - | 21.89 M |
| MONTGOMERYMULTIPLICATION | 25000 / 25000 | - | - | - | - | - | - | - | 2.86 M |
| COMPARE | 25000 / 25000 | - | - | - | - | - | - | - | 61.85 M |
| REDUCE | 21000 / 3125 | - | - | - | - | - | - | - | 29.30 M |
| MODMUL | 21000 / 1562 | - | - | - | - | - | - | - | 4.71 M |
| MODEXP | 21000 / 390 | - | - | - | - | - | - | - | 19.35 k |
| EXPONENTIATION | 21000 / 390 | - | - | - | - | - | - | - | 68.65 k |
| DIVIDE | 21000 / 3125 | - | - | - | - | - | - | - | 15.72 M |
| ISQRT | 21000 / 781 | - | - | - | - | - | - | - | 5.88 M |
| MODMUL_R2 | 25000 / 25000 | - | - | - | - | - | - | - | 5.65 M |

#### p1024 (1024-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 21000 / 12500 | - | - | - | - | - | - | - | 40.19 M |
| SUBTRACT | 21000 / 12500 | - | - | - | - | - | - | - | 52.20 M |
| ADDMOD | 21000 / 12500 | - | - | - | - | - | - | - | 5.83 M |
| SUBTRACTMOD | 21000 / 12500 | - | - | - | - | - | - | - | 19.38 M |
| MULTIPLYOPERANDSCANNING | 21000 / 12500 | - | - | - | - | - | - | - | 6.30 M |
| MULTIPLYPRODUCTSCANNING | 21000 / 12500 | - | - | - | - | - | - | - | 6.28 M |
| MONTGOMERYMULTIPLICATION | 21000 / 12500 | - | - | - | - | - | - | - | 908.01 k |
| COMPARE | 21000 / 12500 | - | - | - | - | - | - | - | 54.51 M |
| REDUCE | 21000 / 1562 | - | - | - | - | - | - | - | 40.52 M |
| MODMUL | 21000 / 781 | - | - | - | - | - | - | - | 1.92 M |
| MODEXP | 21000 / 195 | - | - | - | - | - | - | - | 2.37 k |
| EXPONENTIATION | 21000 / 195 | - | - | - | - | - | - | - | 21.93 k |
| DIVIDE | 21000 / 1562 | - | - | - | - | - | - | - | 10.38 M |
| ISQRT | 21000 / 390 | - | - | - | - | - | - | - | 3.12 M |
| MODMUL_R2 | 21000 / 12500 | - | - | - | - | - | - | - | 2.06 M |

#### p2048 (2048-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 21000 / 6250 | - | - | - | - | - | - | - | 28.38 M |
| SUBTRACT | 21000 / 6250 | - | - | - | - | - | - | - | 30.49 M |
| ADDMOD | 21000 / 6250 | - | - | - | - | - | - | - | 11.83 M |
| SUBTRACTMOD | 21000 / 6250 | - | - | - | - | - | - | - | 14.86 M |
| MULTIPLYOPERANDSCANNING | 21000 / 6250 | - | - | - | - | - | - | - | 2.01 M |
| MULTIPLYPRODUCTSCANNING | 21000 / 6250 | - | - | - | - | - | - | - | 1.84 M |
| MONTGOMERYMULTIPLICATION | 21000 / 6250 | - | - | - | - | - | - | - | 292.06 k |
| COMPARE | 21000 / 6250 | - | - | - | - | - | - | - | 75.74 M |
| REDUCE | 21000 / 781 | - | - | - | - | - | - | - | 29.07 M |
| MODMUL | 21000 / 390 | - | - | - | - | - | - | - | 631.01 k |
| MODEXP | 21000 / 97 | - | - | - | - | - | - | - | 366.0 |
| EXPONENTIATION | 21000 / 97 | - | - | - | - | - | - | - | 3.88 k |
| DIVIDE | 21000 / 781 | - | - | - | - | - | - | - | 13.06 M |
| ISQRT | 21000 / 195 | - | - | - | - | - | - | - | 1.18 M |
| MODMUL_R2 | 21000 / 6250 | - | - | - | - | - | - | - | 678.90 k |

## 5. Head to head

Best OpenCL GPU result against best OpenCL CPU result and the CPU libraries.
Ratios above 1.00x mean the GPU is faster than that baseline.

### secp256k1 (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w16 | 797.36 M | none | n/a | 52.62 M | n/a |
| SUBTRACT | w16 | 768.76 M | none | n/a | 40.61 M | n/a |
| ADDMOD | w16 | 596.71 M | none | n/a | 16.18 M | n/a |
| SUBTRACTMOD | w16 | 591.66 M | none | n/a | 23.94 M | n/a |
| MULTIPLYOPERANDSCANNING | w16 | 39.34 M | none | n/a | 48.93 M | n/a |
| MULTIPLYPRODUCTSCANNING | w16 | 63.07 M | none | n/a | 49.25 M | n/a |
| MONTGOMERYMULTIPLICATION | w16 | 463.44 M | none | n/a | 6.29 M | n/a |
| COMPARE | w16 | 1.03 G | none | n/a | 125.73 M | n/a |
| REDUCE | w16 | 48.03 M | none | n/a | 49.42 M | n/a |
| MODMUL | w16 | 7.94 M | none | n/a | 8.39 M | n/a |
| MODEXP | w16 | 108.99 k | none | n/a | 95.05 k | n/a |
| EXPONENTIATION | w16 | 27.78 k | none | n/a | 255.10 k | n/a |
| DIVIDE | w16 | 19.35 M | none | n/a | 17.17 M | n/a |
| ISQRT | w16 | 394.47 k | none | n/a | 9.00 M | n/a |
| MODMUL_R2 | w16 | 400.26 M | none | n/a | 11.43 M | n/a |

### rsa256(composite) (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w16 | 799.42 M | none | n/a | 49.78 M | n/a |
| SUBTRACT | w16 | 784.52 M | none | n/a | 70.32 M | n/a |
| ADDMOD | w16 | 627.61 M | none | n/a | 23.11 M | n/a |
| SUBTRACTMOD | w16 | 589.09 M | none | n/a | 23.27 M | n/a |
| MULTIPLYOPERANDSCANNING | w16 | 39.35 M | none | n/a | 48.43 M | n/a |
| MULTIPLYPRODUCTSCANNING | w16 | 63.39 M | none | n/a | 51.36 M | n/a |
| MONTGOMERYMULTIPLICATION | w16 | 459.69 M | none | n/a | 6.36 M | n/a |
| COMPARE | w16 | 993.18 M | none | n/a | 129.43 M | n/a |
| REDUCE | w16 | 47.84 M | none | n/a | 30.04 M | n/a |
| MODMUL | w16 | 8.01 M | none | n/a | 4.42 M | n/a |
| MODEXP | w16 | 108.92 k | none | n/a | 110.47 k | n/a |
| EXPONENTIATION | w16 | 27.77 k | none | n/a | 261.71 k | n/a |
| DIVIDE | w16 | 22.41 M | none | n/a | 9.56 M | n/a |
| ISQRT | w16 | 394.59 k | none | n/a | 10.11 M | n/a |
| MODMUL_R2 | w16 | 396.18 M | none | n/a | 11.05 M | n/a |

### brainpoolP512r1 (512-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w16 | 371.28 M | none | n/a | 48.43 M | n/a |
| SUBTRACT | w16 | 377.45 M | none | n/a | 65.92 M | n/a |
| ADDMOD | w16 | 284.07 M | none | n/a | 20.95 M | n/a |
| SUBTRACTMOD | w16 | 267.73 M | none | n/a | 23.62 M | n/a |
| MULTIPLYOPERANDSCANNING | w16 | 8.67 M | none | n/a | 22.83 M | n/a |
| MULTIPLYPRODUCTSCANNING | w16 | 16.18 M | none | n/a | 21.89 M | n/a |
| MONTGOMERYMULTIPLICATION | w16 | 113.01 M | none | n/a | 2.86 M | n/a |
| COMPARE | w16 | 493.55 M | none | n/a | 61.85 M | n/a |
| REDUCE | w16 | 7.25 M | none | n/a | 29.30 M | n/a |
| MODMUL | w16 | 1.25 M | none | n/a | 4.71 M | n/a |
| MODEXP | w16 | 7.90 k | none | n/a | 19.35 k | n/a |
| EXPONENTIATION | w16 | 1.26 k | none | n/a | 68.65 k | n/a |
| DIVIDE | w8 | 1.65 M | none | n/a | 15.72 M | n/a |
| ISQRT | w16 | 28.02 k | none | n/a | 5.88 M | n/a |
| MODMUL_R2 | w16 | 114.85 M | none | n/a | 5.65 M | n/a |

### p1024 (1024-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w16 | 85.80 M | none | n/a | 40.19 M | n/a |
| SUBTRACT | w16 | 85.43 M | none | n/a | 52.20 M | n/a |
| ADDMOD | w16 | 59.73 M | none | n/a | 5.83 M | n/a |
| SUBTRACTMOD | w16 | 60.02 M | none | n/a | 19.38 M | n/a |
| MULTIPLYOPERANDSCANNING | w16 | 474.78 k | none | n/a | 6.30 M | n/a |
| MULTIPLYPRODUCTSCANNING | w16 | 2.48 M | none | n/a | 6.28 M | n/a |
| MONTGOMERYMULTIPLICATION | w16 | 15.36 M | none | n/a | 908.01 k | n/a |
| COMPARE | w16 | 136.85 M | none | n/a | 54.51 M | n/a |
| REDUCE | w16 | 1.02 M | none | n/a | 40.52 M | n/a |
| MODMUL | w16 | 143.30 k | none | n/a | 1.92 M | n/a |
| MODEXP | w16 | 322.6 | none | n/a | 2.37 k | n/a |
| EXPONENTIATION | w16 | 78.5 | none | n/a | 21.93 k | n/a |
| DIVIDE | w16 | 49.75 k | none | n/a | 10.38 M | n/a |
| ISQRT | w16 | 1.28 k | none | n/a | 3.12 M | n/a |
| MODMUL_R2 | w16 | 19.93 M | none | n/a | 2.06 M | n/a |

### p2048 (2048-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w16 | 19.97 M | none | n/a | 28.38 M | n/a |
| SUBTRACT | w16 | 20.30 M | none | n/a | 30.49 M | n/a |
| ADDMOD | w16 | 14.98 M | none | n/a | 11.83 M | n/a |
| SUBTRACTMOD | w16 | 14.26 M | none | n/a | 14.86 M | n/a |
| MULTIPLYOPERANDSCANNING | w16 | 42.23 k | none | n/a | 2.01 M | n/a |
| MULTIPLYPRODUCTSCANNING | w16 | 311.54 k | none | n/a | 1.84 M | n/a |
| MONTGOMERYMULTIPLICATION | w16 | 1.53 M | none | n/a | 292.06 k | n/a |
| COMPARE | w16 | 35.99 M | none | n/a | 75.74 M | n/a |
| REDUCE | w16 | 90.31 k | none | n/a | 29.07 M | n/a |
| MODMUL | w16 | 12.11 k | none | n/a | 631.01 k | n/a |
| MODEXP | w16 | 1.9 | none | n/a | 366.0 | n/a |
| EXPONENTIATION | w16 | 4.7 | none | n/a | 3.88 k | n/a |
| DIVIDE | w16 | 1.21 k | none | n/a | 13.06 M | n/a |
| ISQRT | w8 | 7.4 | none | n/a | 1.18 M | n/a |
| MODMUL_R2 | w8 | 91.53 k | none | n/a | 678.90 k | n/a |

## 6. Raw data

Also written to `Tesla_P40_Report.csv` for analysis.

```csv
kind,device,device_type,kernel,modulus,bits,operation,items,seconds,ops_per_sec,mismatches
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,secp256k1,256,ADD,50000,0.000950252,52617624.821,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,secp256k1,256,ADD,50000,0.099992186,500039.073,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,secp256k1,256,ADD,50000,0.052334023,955401.419,0
library,Tesla P40,gpu,cgbn,secp256k1,256,ADD,50000,0.000038752,1290255986.788,0
opencl-kernel,Tesla P40,GPU,w8,secp256k1,256,ADD,50000,0.000128113,390280558.519,0
opencl-e2e,Tesla P40,GPU,w8,secp256k1,256,ADD,50000,0.001243580,40206503.824,0
opencl-kernel,Tesla P40,GPU,w16,secp256k1,256,ADD,50000,0.000062707,797360668.934,0
opencl-e2e,Tesla P40,GPU,w16,secp256k1,256,ADD,50000,0.001188036,42086264.820,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,secp256k1,256,SUBTRACT,50000,0.001231117,40613526.081,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,secp256k1,256,SUBTRACT,50000,0.032914438,1519090.191,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,secp256k1,256,SUBTRACT,50000,0.062065489,805600.679,0
library,Tesla P40,gpu,cgbn,secp256k1,256,SUBTRACT,50000,0.000037792,1323031329.382,0
opencl-kernel,Tesla P40,GPU,w8,secp256k1,256,SUBTRACT,50000,0.000130336,383623766.812,0
opencl-e2e,Tesla P40,GPU,w8,secp256k1,256,SUBTRACT,50000,0.001312469,38096140.809,0
opencl-kernel,Tesla P40,GPU,w16,secp256k1,256,SUBTRACT,50000,0.000065040,768756787.485,0
opencl-e2e,Tesla P40,GPU,w16,secp256k1,256,SUBTRACT,50000,0.001241229,40282655.048,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,secp256k1,256,ADDMOD,50000,0.003090101,16180700.568,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,secp256k1,256,ADDMOD,50000,0.053800250,929363.712,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,secp256k1,256,ADDMOD,50000,0.071094220,703292.053,0
library,Tesla P40,gpu,cgbn,secp256k1,256,ADDMOD,50000,0.000035840,1395089285.714,0
opencl-kernel,Tesla P40,GPU,w8,secp256k1,256,ADDMOD,50000,0.000174349,286780825.401,0
opencl-e2e,Tesla P40,GPU,w8,secp256k1,256,ADDMOD,50000,0.001335362,37443031.202,0
opencl-kernel,Tesla P40,GPU,w16,secp256k1,256,ADDMOD,50000,0.000083793,596708878.318,0
opencl-e2e,Tesla P40,GPU,w16,secp256k1,256,ADDMOD,50000,0.001222024,40915727.592,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,secp256k1,256,SUBTRACTMOD,50000,0.002088882,23936248.819,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,secp256k1,256,SUBTRACTMOD,50000,0.048334619,1034455.241,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,secp256k1,256,SUBTRACTMOD,50000,0.051765002,965903.564,0
library,Tesla P40,gpu,cgbn,secp256k1,256,SUBTRACTMOD,50000,0.000036192,1381520778.073,0
opencl-kernel,Tesla P40,GPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.000175984,284116935.086,0
opencl-e2e,Tesla P40,GPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.001311359,38128384.495,0
opencl-kernel,Tesla P40,GPU,w16,secp256k1,256,SUBTRACTMOD,50000,0.000084508,591660118.085,0
opencl-e2e,Tesla P40,GPU,w16,secp256k1,256,SUBTRACTMOD,50000,0.001181935,42303512.193,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001021963,48925457.721,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.036952193,1353099.670,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.088449618,565293.566,0
opencl-kernel,Tesla P40,GPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.005960913,8387976.608,0
opencl-e2e,Tesla P40,GPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.007664027,6523985.431,0
opencl-kernel,Tesla P40,GPU,w16,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001270819,39344706.370,0
opencl-e2e,Tesla P40,GPU,w16,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.002619784,19085543.952,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001015150,49253805.427,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.032806615,1524082.866,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.043776376,1142168.554,0
library,Tesla P40,gpu,cgbn,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000036704,1362249346.120,0
opencl-kernel,Tesla P40,GPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.003130736,15970685.170,0
opencl-e2e,Tesla P40,GPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.004904941,10193802.624,0
opencl-kernel,Tesla P40,GPU,w16,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000792726,63073481.768,0
opencl-e2e,Tesla P40,GPU,w16,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.002201569,22711075.253,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.007950403,6288989.353,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.022529877,2219275.318,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.050231761,995386.165,0
library,Tesla P40,gpu,cgbn,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000035840,1395089285.714,0
opencl-kernel,Tesla P40,GPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000361719,138228870.647,0
opencl-e2e,Tesla P40,GPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.001757840,28443996.860,0
opencl-kernel,Tesla P40,GPU,w16,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000107888,463444002.806,0
opencl-e2e,Tesla P40,GPU,w16,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.001361523,36723585.631,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,secp256k1,256,COMPARE,50000,0.000397693,125725151.193,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,secp256k1,256,COMPARE,50000,0.070409686,710129.569,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,secp256k1,256,COMPARE,50000,0.047429869,1054188.028,0
library,Tesla P40,gpu,cgbn,secp256k1,256,COMPARE,50000,0.000037856,1320794590.025,0
opencl-kernel,Tesla P40,GPU,w8,secp256k1,256,COMPARE,50000,0.000085719,583300145.317,0
opencl-e2e,Tesla P40,GPU,w8,secp256k1,256,COMPARE,50000,0.001509906,33114648.768,0
opencl-kernel,Tesla P40,GPU,w16,secp256k1,256,COMPARE,50000,0.000048400,1033054954.613,0
opencl-e2e,Tesla P40,GPU,w16,secp256k1,256,COMPARE,50000,0.001152488,43384402.022,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,secp256k1,256,REDUCE,6250,0.000126462,49421974.777,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,secp256k1,256,REDUCE,6250,0.090523096,69043.153,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,secp256k1,256,REDUCE,6250,0.034955174,178800.426,0
library,Tesla P40,gpu,cgbn,secp256k1,256,REDUCE,50000,0.000035840,1395089285.714,0
opencl-kernel,Tesla P40,GPU,w8,secp256k1,256,REDUCE,21000,0.000235562,89148512.021,0
opencl-e2e,Tesla P40,GPU,w8,secp256k1,256,REDUCE,21000,0.000910013,23076593.848,0
opencl-kernel,Tesla P40,GPU,w16,secp256k1,256,REDUCE,21000,0.000130116,161394437.832,0
opencl-e2e,Tesla P40,GPU,w16,secp256k1,256,REDUCE,21000,0.000741218,28331746.472,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,secp256k1,256,MODMUL,3125,0.000372625,8386448.490,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,secp256k1,256,MODMUL,3125,0.082828272,37728.663,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,secp256k1,256,MODMUL,3125,0.088473788,35321.196,0
library,Tesla P40,gpu,cgbn,secp256k1,256,MODMUL,50000,0.000084640,590737240.076,0
opencl-kernel,Tesla P40,GPU,w8,secp256k1,256,MODMUL,21000,0.000674574,31130750.221,0
opencl-e2e,Tesla P40,GPU,w8,secp256k1,256,MODMUL,21000,0.001489175,14101767.588,0
opencl-kernel,Tesla P40,GPU,w16,secp256k1,256,MODMUL,21000,0.000393346,53388117.741,0
opencl-e2e,Tesla P40,GPU,w16,secp256k1,256,MODMUL,21000,0.001077314,19492925.743,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,secp256k1,256,MODEXP,781,0.008216411,95053.669,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,secp256k1,256,MODEXP,781,0.023953097,32605.387,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,secp256k1,256,MODEXP,781,0.097950152,7973.443,0
library,Tesla P40,gpu,cgbn,secp256k1,256,MODEXP,50000,0.020284416,2464946.489,0
opencl-kernel,Tesla P40,GPU,w8,secp256k1,256,MODEXP,21000,0.029347003,715575.625,0
opencl-e2e,Tesla P40,GPU,w8,secp256k1,256,MODEXP,21000,0.030125014,697095.111,0
opencl-kernel,Tesla P40,GPU,w16,secp256k1,256,MODEXP,21000,0.007165670,2930640.151,0
opencl-e2e,Tesla P40,GPU,w16,secp256k1,256,MODEXP,21000,0.007989023,2628606.740,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,secp256k1,256,EXPONENTIATION,781,0.003061485,255104.956,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,secp256k1,256,EXPONENTIATION,781,0.088459156,8828.933,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,secp256k1,256,EXPONENTIATION,781,0.076011186,10274.804,0
opencl-kernel,Tesla P40,GPU,w8,secp256k1,256,EXPONENTIATION,21000,0.155212795,135298.124,0
opencl-e2e,Tesla P40,GPU,w8,secp256k1,256,EXPONENTIATION,21000,0.156007210,134609.163,0
opencl-kernel,Tesla P40,GPU,w16,secp256k1,256,EXPONENTIATION,21000,0.028111122,747035.288,0
opencl-e2e,Tesla P40,GPU,w16,secp256k1,256,EXPONENTIATION,21000,0.028930889,725867.775,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,secp256k1,256,DIVIDE,6250,0.000364032,17168815.214,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,secp256k1,256,DIVIDE,6250,0.100218664,62363.633,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,secp256k1,256,DIVIDE,6250,0.091564133,68258.168,0
library,Tesla P40,gpu,cgbn,secp256k1,256,DIVIDE,50000,0.000036864,1356336805.556,0
opencl-kernel,Tesla P40,GPU,w8,secp256k1,256,DIVIDE,21000,0.000455995,46053127.427,0
opencl-e2e,Tesla P40,GPU,w8,secp256k1,256,DIVIDE,21000,0.001344174,15622978.384,0
opencl-kernel,Tesla P40,GPU,w16,secp256k1,256,DIVIDE,21000,0.000323043,65006784.495,0
opencl-e2e,Tesla P40,GPU,w16,secp256k1,256,DIVIDE,21000,0.001063168,19752290.921,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,secp256k1,256,ISQRT,1562,0.000173577,8998894.600,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,secp256k1,256,ISQRT,1562,0.089685943,17416.330,0
opencl-kernel,Tesla P40,GPU,w8,secp256k1,256,ISQRT,21000,0.006450470,3255576.689,0
opencl-e2e,Tesla P40,GPU,w8,secp256k1,256,ISQRT,21000,0.007211928,2911842.705,0
opencl-kernel,Tesla P40,GPU,w16,secp256k1,256,ISQRT,21000,0.003959726,5303397.153,0
opencl-e2e,Tesla P40,GPU,w16,secp256k1,256,ISQRT,21000,0.004640270,4525598.702,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,secp256k1,256,MODMUL_R2,50000,0.004374592,11429637.064,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,secp256k1,256,MODMUL_R2,50000,0.040578356,1232183.975,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,secp256k1,256,MODMUL_R2,50000,0.100026827,499865.901,0
library,Tesla P40,gpu,cgbn,secp256k1,256,MODMUL_R2,50000,0.000045824,1091131284.916,0
opencl-kernel,Tesla P40,GPU,w8,secp256k1,256,MODMUL_R2,50000,0.000389123,128494041.242,0
opencl-e2e,Tesla P40,GPU,w8,secp256k1,256,MODMUL_R2,50000,0.001695957,29481878.947,0
opencl-kernel,Tesla P40,GPU,w16,secp256k1,256,MODMUL_R2,50000,0.000124920,400255652.641,0
opencl-e2e,Tesla P40,GPU,w16,secp256k1,256,MODMUL_R2,50000,0.001443684,34633621.043,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,rsa256(composite),256,ADD,50000,0.001004452,49778380.743,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,rsa256(composite),256,ADD,50000,0.025347061,1972615.287,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,rsa256(composite),256,ADD,50000,0.027153611,1841375.723,0
library,Tesla P40,gpu,cgbn,rsa256(composite),256,ADD,50000,0.000035840,1395089285.714,0
opencl-kernel,Tesla P40,GPU,w8,rsa256(composite),256,ADD,50000,0.000108891,459175049.767,0
opencl-e2e,Tesla P40,GPU,w8,rsa256(composite),256,ADD,50000,0.001594737,31353133.403,0
opencl-kernel,Tesla P40,GPU,w16,rsa256(composite),256,ADD,50000,0.000062545,799423609.513,0
opencl-e2e,Tesla P40,GPU,w16,rsa256(composite),256,ADD,50000,0.001224439,40835030.886,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACT,50000,0.000711011,70322399.893,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACT,50000,0.063609388,786047.494,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACT,50000,0.043499793,1149430.753,0
library,Tesla P40,gpu,cgbn,rsa256(composite),256,SUBTRACT,50000,0.000038560,1296680497.925,0
opencl-kernel,Tesla P40,GPU,w8,rsa256(composite),256,SUBTRACT,50000,0.000109040,458547551.156,0
opencl-e2e,Tesla P40,GPU,w8,rsa256(composite),256,SUBTRACT,50000,0.001559330,32065055.274,0
opencl-kernel,Tesla P40,GPU,w16,rsa256(composite),256,SUBTRACT,50000,0.000063733,784523363.448,0
opencl-e2e,Tesla P40,GPU,w16,rsa256(composite),256,SUBTRACT,50000,0.001275201,39209502.082,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,rsa256(composite),256,ADDMOD,50000,0.002163506,23110634.955,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,rsa256(composite),256,ADDMOD,50000,0.025325565,1974289.606,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,rsa256(composite),256,ADDMOD,50000,0.073642886,678952.207,0
library,Tesla P40,gpu,cgbn,rsa256(composite),256,ADDMOD,50000,0.000035840,1395089285.714,0
opencl-kernel,Tesla P40,GPU,w8,rsa256(composite),256,ADDMOD,50000,0.000142060,351963990.081,0
opencl-e2e,Tesla P40,GPU,w8,rsa256(composite),256,ADDMOD,50000,0.001625795,30754182.035,0
opencl-kernel,Tesla P40,GPU,w16,rsa256(composite),256,ADDMOD,50000,0.000079667,627612729.457,0
opencl-e2e,Tesla P40,GPU,w16,rsa256(composite),256,ADDMOD,50000,0.001165269,42908549.406,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACTMOD,50000,0.002148890,23267827.333,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACTMOD,50000,0.064464318,775622.882,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACTMOD,50000,0.041747683,1197671.259,0
library,Tesla P40,gpu,cgbn,rsa256(composite),256,SUBTRACTMOD,50000,0.000036864,1356336805.556,0
opencl-kernel,Tesla P40,GPU,w8,rsa256(composite),256,SUBTRACTMOD,50000,0.000153891,324905197.766,0
opencl-e2e,Tesla P40,GPU,w8,rsa256(composite),256,SUBTRACTMOD,50000,0.001727827,28938084.129,0
opencl-kernel,Tesla P40,GPU,w16,rsa256(composite),256,SUBTRACTMOD,50000,0.000084876,589094104.351,0
opencl-e2e,Tesla P40,GPU,w16,rsa256(composite),256,SUBTRACTMOD,50000,0.001208110,41386959.318,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001032416,48430089.714,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.068473373,730210.851,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.033545605,1490508.214,0
opencl-kernel,Tesla P40,GPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.005050814,9899394.411,0
opencl-e2e,Tesla P40,GPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.006831299,7319252.364,0
opencl-kernel,Tesla P40,GPU,w16,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001270761,39346501.359,0
opencl-e2e,Tesla P40,GPU,w16,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.002756699,18137634.980,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000973460,51363177.093,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.040642102,1230251.329,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.047829913,1045370.917,0
library,Tesla P40,gpu,cgbn,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000035840,1395089285.714,0
opencl-kernel,Tesla P40,GPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.003124874,16000645.902,0
opencl-e2e,Tesla P40,GPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.004782213,10455410.720,0
opencl-kernel,Tesla P40,GPU,w16,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000788830,63385014.222,0
opencl-e2e,Tesla P40,GPU,w16,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.002346759,21305980.195,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.007861365,6360218.722,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.029522300,1693634.981,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.034855595,1434489.935,0
library,Tesla P40,gpu,cgbn,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000036800,1358695652.174,0
opencl-kernel,Tesla P40,GPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000361896,138161282.764,0
opencl-e2e,Tesla P40,GPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.001772394,28210432.751,0
opencl-kernel,Tesla P40,GPU,w16,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000108768,459694030.461,0
opencl-e2e,Tesla P40,GPU,w16,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.001309200,38191263.088,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,rsa256(composite),256,COMPARE,50000,0.000386300,129433072.677,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,rsa256(composite),256,COMPARE,50000,0.045214873,1105830.823,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,rsa256(composite),256,COMPARE,50000,0.029838854,1675667.559,0
library,Tesla P40,gpu,cgbn,rsa256(composite),256,COMPARE,50000,0.000036640,1364628820.961,0
opencl-kernel,Tesla P40,GPU,w8,rsa256(composite),256,COMPARE,50000,0.000085044,587931273.254,0
opencl-e2e,Tesla P40,GPU,w8,rsa256(composite),256,COMPARE,50000,0.001565589,31936864.732,0
opencl-kernel,Tesla P40,GPU,w16,rsa256(composite),256,COMPARE,50000,0.000050343,993184619.511,0
opencl-e2e,Tesla P40,GPU,w16,rsa256(composite),256,COMPARE,50000,0.001146797,43599692.049,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,rsa256(composite),256,REDUCE,6250,0.000208077,30036988.060,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,rsa256(composite),256,REDUCE,6250,0.062448208,100082.936,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,rsa256(composite),256,REDUCE,6250,0.082908367,75384.430,0
library,Tesla P40,gpu,cgbn,rsa256(composite),256,REDUCE,50000,0.000037888,1319679054.054,0
opencl-kernel,Tesla P40,GPU,w8,rsa256(composite),256,REDUCE,21000,0.000237857,88288344.436,0
opencl-e2e,Tesla P40,GPU,w8,rsa256(composite),256,REDUCE,21000,0.000931541,22543286.268,0
opencl-kernel,Tesla P40,GPU,w16,rsa256(composite),256,REDUCE,21000,0.000130643,160743231.635,0
opencl-e2e,Tesla P40,GPU,w16,rsa256(composite),256,REDUCE,21000,0.000741076,28337176.232,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,rsa256(composite),256,MODMUL,3125,0.000707248,4418534.985,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,rsa256(composite),256,MODMUL,3125,0.087998973,35511.778,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,rsa256(composite),256,MODMUL,3125,0.087995822,35513.050,0
library,Tesla P40,gpu,cgbn,rsa256(composite),256,MODMUL,50000,0.000084672,590513983.371,0
opencl-kernel,Tesla P40,GPU,w8,rsa256(composite),256,MODMUL,21000,0.000674069,31154073.243,0
opencl-e2e,Tesla P40,GPU,w8,rsa256(composite),256,MODMUL,21000,0.001437387,14609844.484,0
opencl-kernel,Tesla P40,GPU,w16,rsa256(composite),256,MODMUL,21000,0.000389896,53860502.783,0
opencl-e2e,Tesla P40,GPU,w16,rsa256(composite),256,MODMUL,21000,0.000966500,21727888.785,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,rsa256(composite),256,MODEXP,781,0.007069690,110471.602,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,rsa256(composite),256,MODEXP,781,0.071535696,10917.626,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,rsa256(composite),256,MODEXP,781,0.093807269,8325.581,0
library,Tesla P40,gpu,cgbn,rsa256(composite),256,MODEXP,50000,0.020103071,2487182.182,0
opencl-kernel,Tesla P40,GPU,w8,rsa256(composite),256,MODEXP,21000,0.029281274,717181.911,0
opencl-e2e,Tesla P40,GPU,w8,rsa256(composite),256,MODEXP,21000,0.030049303,698851.482,0
opencl-kernel,Tesla P40,GPU,w16,rsa256(composite),256,MODEXP,21000,0.007170147,2928810.234,0
opencl-e2e,Tesla P40,GPU,w16,rsa256(composite),256,MODEXP,21000,0.007842241,2677806.049,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,rsa256(composite),256,EXPONENTIATION,781,0.002984260,261706.422,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,rsa256(composite),256,EXPONENTIATION,781,0.096823986,8066.183,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,rsa256(composite),256,EXPONENTIATION,781,0.095972986,8137.707,0
opencl-kernel,Tesla P40,GPU,w8,rsa256(composite),256,EXPONENTIATION,21000,0.155186883,135320.715,0
opencl-e2e,Tesla P40,GPU,w8,rsa256(composite),256,EXPONENTIATION,21000,0.155993580,134620.925,0
opencl-kernel,Tesla P40,GPU,w16,rsa256(composite),256,EXPONENTIATION,21000,0.028120324,746790.826,0
opencl-e2e,Tesla P40,GPU,w16,rsa256(composite),256,EXPONENTIATION,21000,0.028887889,726948.241,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,rsa256(composite),256,DIVIDE,6250,0.000654022,9556253.803,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,rsa256(composite),256,DIVIDE,6250,0.099994736,62503.290,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,rsa256(composite),256,DIVIDE,6250,0.097546857,64071.772,0
library,Tesla P40,gpu,cgbn,rsa256(composite),256,DIVIDE,50000,0.000035840,1395089285.714,0
opencl-kernel,Tesla P40,GPU,w8,rsa256(composite),256,DIVIDE,21000,0.000450749,46589102.740,0
opencl-e2e,Tesla P40,GPU,w8,rsa256(composite),256,DIVIDE,21000,0.001392831,15077205.041,0
opencl-kernel,Tesla P40,GPU,w16,rsa256(composite),256,DIVIDE,21000,0.000278940,75285018.686,0
opencl-e2e,Tesla P40,GPU,w16,rsa256(composite),256,DIVIDE,21000,0.001141332,18399552.756,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,rsa256(composite),256,ISQRT,1562,0.000154554,10106490.485,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,rsa256(composite),256,ISQRT,1562,0.044039569,35468.104,0
opencl-kernel,Tesla P40,GPU,w8,rsa256(composite),256,ISQRT,21000,0.006442438,3259635.621,0
opencl-e2e,Tesla P40,GPU,w8,rsa256(composite),256,ISQRT,21000,0.007215485,2910407.279,0
opencl-kernel,Tesla P40,GPU,w16,rsa256(composite),256,ISQRT,21000,0.003958519,5305014.526,0
opencl-e2e,Tesla P40,GPU,w16,rsa256(composite),256,ISQRT,21000,0.004683399,4483922.792,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,rsa256(composite),256,MODMUL_R2,50000,0.004525447,11048632.246,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,rsa256(composite),256,MODMUL_R2,50000,0.047959600,1042544.139,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,rsa256(composite),256,MODMUL_R2,50000,0.025091660,1992693.989,0
library,Tesla P40,gpu,cgbn,rsa256(composite),256,MODMUL_R2,50000,0.000044896,1113684960.798,0
opencl-kernel,Tesla P40,GPU,w8,rsa256(composite),256,MODMUL_R2,50000,0.000387155,129147245.700,0
opencl-e2e,Tesla P40,GPU,w8,rsa256(composite),256,MODMUL_R2,50000,0.001683166,29705922.929,0
opencl-kernel,Tesla P40,GPU,w16,rsa256(composite),256,MODMUL_R2,50000,0.000126204,396183994.598,0
opencl-e2e,Tesla P40,GPU,w16,rsa256(composite),256,MODMUL_R2,50000,0.001579079,31664025.368,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,brainpoolP512r1,512,ADD,25000,0.000516232,48427850.815,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,brainpoolP512r1,512,ADD,25000,0.091828341,272247.105,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,brainpoolP512r1,512,ADD,25000,0.074114283,337316.898,0
library,Tesla P40,gpu,cgbn,brainpoolP512r1,512,ADD,50000,0.000050176,996492346.939,0
opencl-kernel,Tesla P40,GPU,w8,brainpoolP512r1,512,ADD,25000,0.000122273,204460719.250,0
opencl-e2e,Tesla P40,GPU,w8,brainpoolP512r1,512,ADD,25000,0.001411282,17714390.282,0
opencl-kernel,Tesla P40,GPU,w16,brainpoolP512r1,512,ADD,25000,0.000067334,371283873.622,0
opencl-e2e,Tesla P40,GPU,w16,brainpoolP512r1,512,ADD,25000,0.001485301,16831605.756,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACT,25000,0.000379227,65923598.894,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACT,25000,0.067256103,371713.478,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACT,25000,0.052507449,476122.922,0
library,Tesla P40,gpu,cgbn,brainpoolP512r1,512,SUBTRACT,50000,0.000051200,976562500.000,0
opencl-kernel,Tesla P40,GPU,w8,brainpoolP512r1,512,SUBTRACT,25000,0.000121369,205983360.830,0
opencl-e2e,Tesla P40,GPU,w8,brainpoolP512r1,512,SUBTRACT,25000,0.001523440,16410229.424,0
opencl-kernel,Tesla P40,GPU,w16,brainpoolP512r1,512,SUBTRACT,25000,0.000066234,377449467.612,0
opencl-e2e,Tesla P40,GPU,w16,brainpoolP512r1,512,SUBTRACT,25000,0.001407021,17768036.684,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,brainpoolP512r1,512,ADDMOD,25000,0.001193461,20947479.494,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,brainpoolP512r1,512,ADDMOD,25000,0.016957426,1474280.356,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,brainpoolP512r1,512,ADDMOD,25000,0.040220094,621579.848,0
library,Tesla P40,gpu,cgbn,brainpoolP512r1,512,ADDMOD,50000,0.000051200,976562500.000,0
opencl-kernel,Tesla P40,GPU,w8,brainpoolP512r1,512,ADDMOD,25000,0.000170187,146897156.022,0
opencl-e2e,Tesla P40,GPU,w8,brainpoolP512r1,512,ADDMOD,25000,0.001449064,17252517.465,0
opencl-kernel,Tesla P40,GPU,w16,brainpoolP512r1,512,ADDMOD,25000,0.000088007,284068453.148,0
opencl-e2e,Tesla P40,GPU,w16,brainpoolP512r1,512,ADDMOD,25000,0.001383210,18073899.693,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACTMOD,25000,0.001058325,23622234.616,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACTMOD,25000,0.036467481,685542.275,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACTMOD,25000,0.041061622,608841.026,0
library,Tesla P40,gpu,cgbn,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000050784,984562066.793,0
opencl-kernel,Tesla P40,GPU,w8,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000186360,134148976.649,0
opencl-e2e,Tesla P40,GPU,w8,brainpoolP512r1,512,SUBTRACTMOD,25000,0.001433787,17436342.437,0
opencl-kernel,Tesla P40,GPU,w16,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000093377,267731990.196,0
opencl-e2e,Tesla P40,GPU,w16,brainpoolP512r1,512,SUBTRACTMOD,25000,0.001383535,18069653.622,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.001095080,22829381.292,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.025174384,993072.958,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.057599198,434033.821,0
opencl-kernel,Tesla P40,GPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.011174738,2237188.928,0
opencl-e2e,Tesla P40,GPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.012802280,1952777.194,0
opencl-kernel,Tesla P40,GPU,w16,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.002882914,8671781.642,0
opencl-e2e,Tesla P40,GPU,w16,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.004510652,5542435.943,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.001142166,21888235.201,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.046236414,540699.372,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.053700760,465542.759,0
library,Tesla P40,gpu,cgbn,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000051200,976562500.000,0
opencl-kernel,Tesla P40,GPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.006157106,4060349.117,0
opencl-e2e,Tesla P40,GPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.007764787,3219663.356,0
opencl-kernel,Tesla P40,GPU,w16,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.001545121,16179960.565,0
opencl-e2e,Tesla P40,GPU,w16,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.003091525,8086623.775,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.008737209,2861325.618,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.027754270,900762.297,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.023771276,1051689.442,0
library,Tesla P40,gpu,cgbn,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000062464,800461065.574,0
opencl-kernel,Tesla P40,GPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000965621,25890072.483,0
opencl-e2e,Tesla P40,GPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.002321288,10769882.975,0
opencl-kernel,Tesla P40,GPU,w16,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000221223,113008087.599,0
opencl-e2e,Tesla P40,GPU,w16,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.001455507,17176145.396,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,brainpoolP512r1,512,COMPARE,25000,0.000404174,61854516.149,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,brainpoolP512r1,512,COMPARE,25000,0.071352946,350370.957,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,brainpoolP512r1,512,COMPARE,25000,0.040584884,615992.889,0
library,Tesla P40,gpu,cgbn,brainpoolP512r1,512,COMPARE,50000,0.000051392,972914072.229,0
opencl-kernel,Tesla P40,GPU,w8,brainpoolP512r1,512,COMPARE,25000,0.000089508,279304594.287,0
opencl-e2e,Tesla P40,GPU,w8,brainpoolP512r1,512,COMPARE,25000,0.001399869,17858815.636,0
opencl-kernel,Tesla P40,GPU,w16,brainpoolP512r1,512,COMPARE,25000,0.000050654,493545058.996,0
opencl-e2e,Tesla P40,GPU,w16,brainpoolP512r1,512,COMPARE,25000,0.001127393,22175053.716,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,brainpoolP512r1,512,REDUCE,3125,0.000106666,29297056.930,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,brainpoolP512r1,512,REDUCE,3125,0.079493401,39311.439,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,brainpoolP512r1,512,REDUCE,3125,0.084721434,36885.589,0
library,Tesla P40,gpu,cgbn,brainpoolP512r1,512,REDUCE,50000,0.000051200,976562500.000,0
opencl-kernel,Tesla P40,GPU,w8,brainpoolP512r1,512,REDUCE,21000,0.000759266,27658296.453,0
opencl-e2e,Tesla P40,GPU,w8,brainpoolP512r1,512,REDUCE,21000,0.001914641,10968114.243,0
opencl-kernel,Tesla P40,GPU,w16,brainpoolP512r1,512,REDUCE,21000,0.000431067,48716316.088,0
opencl-e2e,Tesla P40,GPU,w16,brainpoolP512r1,512,REDUCE,21000,0.001572461,13354863.345,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL,1562,0.000331328,4714361.348,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL,1562,0.092247740,16932.664,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL,1562,0.099982033,15622.807,0
library,Tesla P40,gpu,cgbn,brainpoolP512r1,512,MODMUL,50000,0.000211968,235884661.836,0
opencl-kernel,Tesla P40,GPU,w8,brainpoolP512r1,512,MODMUL,21000,0.002421840,8671092.717,0
opencl-e2e,Tesla P40,GPU,w8,brainpoolP512r1,512,MODMUL,21000,0.003729835,5630275.991,0
opencl-kernel,Tesla P40,GPU,w16,brainpoolP512r1,512,MODMUL,21000,0.001249450,16807394.903,0
opencl-e2e,Tesla P40,GPU,w16,brainpoolP512r1,512,MODMUL,21000,0.002282267,9201377.520,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,brainpoolP512r1,512,MODEXP,390,0.020153182,19351.783,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,brainpoolP512r1,512,MODEXP,390,0.025081169,15549.514,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,brainpoolP512r1,512,MODEXP,390,0.098091459,3975.881,0
library,Tesla P40,gpu,cgbn,brainpoolP512r1,512,MODEXP,50000,0.101150721,494311.850,0
opencl-kernel,Tesla P40,GPU,w8,brainpoolP512r1,512,MODEXP,21000,0.516889159,40627.666,0
opencl-e2e,Tesla P40,GPU,w8,brainpoolP512r1,512,MODEXP,21000,0.518371165,40511.513,0
opencl-kernel,Tesla P40,GPU,w16,brainpoolP512r1,512,MODEXP,21000,0.049343605,425587.064,0
opencl-e2e,Tesla P40,GPU,w16,brainpoolP512r1,512,MODEXP,21000,0.050688557,414294.691,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,brainpoolP512r1,512,EXPONENTIATION,390,0.005681265,68646.684,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,brainpoolP512r1,512,EXPONENTIATION,390,0.017246789,22612.905,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,brainpoolP512r1,512,EXPONENTIATION,390,0.121972468,3197.443,0
opencl-kernel,Tesla P40,GPU,w8,brainpoolP512r1,512,EXPONENTIATION,21000,1.263222017,16624.156,0
opencl-e2e,Tesla P40,GPU,w8,brainpoolP512r1,512,EXPONENTIATION,21000,1.262527809,16633.297,0
opencl-kernel,Tesla P40,GPU,w16,brainpoolP512r1,512,EXPONENTIATION,21000,0.308650542,68038.112,0
opencl-e2e,Tesla P40,GPU,w16,brainpoolP512r1,512,EXPONENTIATION,21000,0.310038782,67733.462,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,brainpoolP512r1,512,DIVIDE,3125,0.000198815,15718128.502,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,brainpoolP512r1,512,DIVIDE,3125,0.100437556,31113.859,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,brainpoolP512r1,512,DIVIDE,3125,0.082853920,37716.984,0
library,Tesla P40,gpu,cgbn,brainpoolP512r1,512,DIVIDE,50000,0.000063104,792342799.189,0
opencl-kernel,Tesla P40,GPU,w8,brainpoolP512r1,512,DIVIDE,21000,0.001890174,11110087.612,0
opencl-e2e,Tesla P40,GPU,w8,brainpoolP512r1,512,DIVIDE,21000,0.003443269,6098855.785,0
opencl-kernel,Tesla P40,GPU,w16,brainpoolP512r1,512,DIVIDE,21000,0.002417028,8688356.292,0
opencl-e2e,Tesla P40,GPU,w16,brainpoolP512r1,512,DIVIDE,21000,0.003968815,5291251.808,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,brainpoolP512r1,512,ISQRT,781,0.000132806,5880762.799,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,brainpoolP512r1,512,ISQRT,781,0.090605600,8619.776,0
opencl-kernel,Tesla P40,GPU,w8,brainpoolP512r1,512,ISQRT,21000,0.047996990,437527.435,0
opencl-e2e,Tesla P40,GPU,w8,brainpoolP512r1,512,ISQRT,21000,0.049420351,424926.161,0
opencl-kernel,Tesla P40,GPU,w16,brainpoolP512r1,512,ISQRT,21000,0.027876265,753329.039,0
opencl-e2e,Tesla P40,GPU,w16,brainpoolP512r1,512,ISQRT,21000,0.029196425,719266.146,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL_R2,25000,0.004422231,5653255.265,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL_R2,25000,0.029386666,850725.975,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL_R2,25000,0.065164088,383646.894,0
library,Tesla P40,gpu,cgbn,brainpoolP512r1,512,MODMUL_R2,50000,0.000117088,427029242.963,0
opencl-kernel,Tesla P40,GPU,w8,brainpoolP512r1,512,MODMUL_R2,25000,0.000996157,25096449.962,0
opencl-e2e,Tesla P40,GPU,w8,brainpoolP512r1,512,MODMUL_R2,25000,0.002264991,11037571.759,0
opencl-kernel,Tesla P40,GPU,w16,brainpoolP512r1,512,MODMUL_R2,25000,0.000217676,114849629.965,0
opencl-e2e,Tesla P40,GPU,w16,brainpoolP512r1,512,MODMUL_R2,25000,0.001564794,15976544.008,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,p1024,1024,ADD,12500,0.000310988,40194455.575,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,p1024,1024,ADD,12500,0.047362756,263920.451,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,p1024,1024,ADD,12500,0.052826808,236622.285,0
library,Tesla P40,gpu,cgbn,p1024,1024,ADD,50000,0.000092000,543478260.870,0
opencl-kernel,Tesla P40,GPU,w8,p1024,1024,ADD,21000,0.000277071,75792839.912,0
opencl-e2e,Tesla P40,GPU,w8,p1024,1024,ADD,21000,0.002152503,9756083.953,0
opencl-kernel,Tesla P40,GPU,w16,p1024,1024,ADD,21000,0.000145681,144150605.672,0
opencl-e2e,Tesla P40,GPU,w16,p1024,1024,ADD,21000,0.002077544,10108088.641,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,p1024,1024,SUBTRACT,12500,0.000239477,52197073.539,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,p1024,1024,SUBTRACT,12500,0.088889262,140624.410,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,p1024,1024,SUBTRACT,12500,0.076357944,163702.679,0
library,Tesla P40,gpu,cgbn,p1024,1024,SUBTRACT,50000,0.000093184,536572802.198,0
opencl-kernel,Tesla P40,GPU,w8,p1024,1024,SUBTRACT,21000,0.000279189,75217839.607,0
opencl-e2e,Tesla P40,GPU,w8,p1024,1024,SUBTRACT,21000,0.002267420,9261627.512,0
opencl-kernel,Tesla P40,GPU,w16,p1024,1024,SUBTRACT,21000,0.000146325,143516166.025,0
opencl-e2e,Tesla P40,GPU,w16,p1024,1024,SUBTRACT,21000,0.002078513,10103377.196,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,p1024,1024,ADDMOD,12500,0.002143586,5831350.382,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,p1024,1024,ADDMOD,12500,0.029518514,423463.052,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,p1024,1024,ADDMOD,12500,0.064888731,192637.455,0
library,Tesla P40,gpu,cgbn,p1024,1024,ADDMOD,50000,0.000091136,548630617.978,0
opencl-kernel,Tesla P40,GPU,w8,p1024,1024,ADDMOD,21000,0.000402392,52187900.564,0
opencl-e2e,Tesla P40,GPU,w8,p1024,1024,ADDMOD,21000,0.002384029,8808617.326,0
opencl-kernel,Tesla P40,GPU,w16,p1024,1024,ADDMOD,21000,0.000209275,100346465.475,0
opencl-e2e,Tesla P40,GPU,w16,p1024,1024,ADDMOD,21000,0.002208886,9507054.818,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,p1024,1024,SUBTRACTMOD,12500,0.000645003,19379753.360,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,p1024,1024,SUBTRACTMOD,12500,0.045700905,273517.560,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,p1024,1024,SUBTRACTMOD,12500,0.024657960,506935.690,0
library,Tesla P40,gpu,cgbn,p1024,1024,SUBTRACTMOD,50000,0.000091136,548630617.978,0
opencl-kernel,Tesla P40,GPU,w8,p1024,1024,SUBTRACTMOD,21000,0.000403089,52097677.877,0
opencl-e2e,Tesla P40,GPU,w8,p1024,1024,SUBTRACTMOD,21000,0.002566610,8181998.832,0
opencl-kernel,Tesla P40,GPU,w16,p1024,1024,SUBTRACTMOD,21000,0.000208258,100836496.187,0
opencl-e2e,Tesla P40,GPU,w16,p1024,1024,SUBTRACTMOD,21000,0.002458086,8543232.402,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.001983803,6301028.887,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.056225671,222318.379,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.042083044,297031.745,0
opencl-kernel,Tesla P40,GPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,21000,0.087202980,240817.458,0
opencl-e2e,Tesla P40,GPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,21000,0.090571774,231860.314,0
opencl-kernel,Tesla P40,GPU,w16,p1024,1024,MULTIPLYOPERANDSCANNING,21000,0.026327818,797635.416,0
opencl-e2e,Tesla P40,GPU,w16,p1024,1024,MULTIPLYOPERANDSCANNING,21000,0.029526736,711219.827,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.001990146,6280946.210,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.018823113,664077.197,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.024840917,503202.036,0
library,Tesla P40,gpu,cgbn,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000128000,390625000.000,0
opencl-kernel,Tesla P40,GPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,21000,0.019980597,1051019.655,0
opencl-e2e,Tesla P40,GPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,21000,0.022573338,930301.043,0
opencl-kernel,Tesla P40,GPU,w16,p1024,1024,MULTIPLYPRODUCTSCANNING,21000,0.005042749,4164395.261,0
opencl-e2e,Tesla P40,GPU,w16,p1024,1024,MULTIPLYPRODUCTSCANNING,21000,0.007761977,2705496.242,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.013766402,908007.777,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.024462307,510990.233,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.062118824,201227.248,0
library,Tesla P40,gpu,cgbn,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000209920,238185975.610,0
opencl-kernel,Tesla P40,GPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,21000,0.004268232,4920069.910,0
opencl-e2e,Tesla P40,GPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,21000,0.006352906,3305573.920,0
opencl-kernel,Tesla P40,GPU,w16,p1024,1024,MONTGOMERYMULTIPLICATION,21000,0.000813690,25808354.083,0
opencl-e2e,Tesla P40,GPU,w16,p1024,1024,MONTGOMERYMULTIPLICATION,21000,0.003061492,6859400.589,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,p1024,1024,COMPARE,12500,0.000229308,54511855.575,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,p1024,1024,COMPARE,12500,0.025087641,498253.306,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,p1024,1024,COMPARE,12500,0.035893704,348250.489,0
library,Tesla P40,gpu,cgbn,p1024,1024,COMPARE,50000,0.000091136,548630617.978,0
opencl-kernel,Tesla P40,GPU,w8,p1024,1024,COMPARE,21000,0.000174529,120323898.761,0
opencl-e2e,Tesla P40,GPU,w8,p1024,1024,COMPARE,21000,0.002061380,10187351.009,0
opencl-kernel,Tesla P40,GPU,w16,p1024,1024,COMPARE,21000,0.000091341,229908064.664,0
opencl-e2e,Tesla P40,GPU,w16,p1024,1024,COMPARE,21000,0.002008679,10454632.867,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,p1024,1024,REDUCE,1562,0.000038552,40516601.741,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,p1024,1024,REDUCE,1562,0.081478349,19170.737,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,p1024,1024,REDUCE,1562,0.085734945,18218.942,0
library,Tesla P40,gpu,cgbn,p1024,1024,REDUCE,50000,0.000092160,542534722.222,0
opencl-kernel,Tesla P40,GPU,w8,p1024,1024,REDUCE,21000,0.004693032,4474719.106,0
opencl-e2e,Tesla P40,GPU,w8,p1024,1024,REDUCE,21000,0.006719952,3125022.329,0
opencl-kernel,Tesla P40,GPU,w16,p1024,1024,REDUCE,21000,0.001536562,13666875.605,0
opencl-e2e,Tesla P40,GPU,w16,p1024,1024,REDUCE,21000,0.003896699,5389176.764,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,p1024,1024,MODMUL,781,0.000405857,1924323.697,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,p1024,1024,MODMUL,781,0.035995527,21697.140,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,p1024,1024,MODMUL,781,0.099943100,7814.446,0
library,Tesla P40,gpu,cgbn,p1024,1024,MODMUL,50000,0.000745312,67085998.884,0
opencl-kernel,Tesla P40,GPU,w8,p1024,1024,MODMUL,21000,0.016961350,1238108.987,0
opencl-e2e,Tesla P40,GPU,w8,p1024,1024,MODMUL,21000,0.018977229,1106589.375,0
opencl-kernel,Tesla P40,GPU,w16,p1024,1024,MODMUL,21000,0.005450154,3853102.163,0
opencl-e2e,Tesla P40,GPU,w16,p1024,1024,MODMUL,21000,0.007426298,2827788.493,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,p1024,1024,MODEXP,195,0.082387942,2366.851,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,p1024,1024,MODEXP,195,0.067724896,2879.296,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,p1024,1024,MODEXP,195,0.106234650,1835.559,0
library,Tesla P40,gpu,cgbn,p1024,1024,MODEXP,50000,0.651071310,76796.503,0
opencl-kernel,Tesla P40,GPU,w8,p1024,1024,MODEXP,21000,4.409977441,4761.929,0
opencl-e2e,Tesla P40,GPU,w8,p1024,1024,MODEXP,21000,4.403779918,4768.631,0
opencl-kernel,Tesla P40,GPU,w16,p1024,1024,MODEXP,21000,0.604520269,34738.289,0
opencl-e2e,Tesla P40,GPU,w16,p1024,1024,MODEXP,21000,0.607743064,34554.076,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,p1024,1024,EXPONENTIATION,195,0.008892343,21928.979,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,p1024,1024,EXPONENTIATION,195,0.053218517,3664.138,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,p1024,1024,EXPONENTIATION,195,0.111961676,1741.667,0
opencl-kernel,Tesla P40,GPU,w8,p1024,1024,EXPONENTIATION,21000,10.376377043,2023.828,0
opencl-e2e,Tesla P40,GPU,w8,p1024,1024,EXPONENTIATION,21000,10.385627565,2022.025,0
opencl-kernel,Tesla P40,GPU,w16,p1024,1024,EXPONENTIATION,21000,2.485098736,8450.368,0
opencl-e2e,Tesla P40,GPU,w16,p1024,1024,EXPONENTIATION,21000,2.487185239,8443.279,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,p1024,1024,DIVIDE,1562,0.000150533,10376467.311,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,p1024,1024,DIVIDE,1562,0.047390319,32960.318,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,p1024,1024,DIVIDE,1562,0.083172354,18780.279,0
library,Tesla P40,gpu,cgbn,p1024,1024,DIVIDE,50000,0.000112320,445156695.157,0
opencl-kernel,Tesla P40,GPU,w8,p1024,1024,DIVIDE,21000,0.131857924,159262.328,0
opencl-e2e,Tesla P40,GPU,w8,p1024,1024,DIVIDE,21000,0.133809834,156939.138,0
opencl-kernel,Tesla P40,GPU,w16,p1024,1024,DIVIDE,21000,0.031395073,668894.765,0
opencl-e2e,Tesla P40,GPU,w16,p1024,1024,DIVIDE,21000,0.033992498,617783.368,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,p1024,1024,ISQRT,390,0.000125049,3118779.561,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,p1024,1024,ISQRT,390,0.096120379,4057.412,0
opencl-kernel,Tesla P40,GPU,w8,p1024,1024,ISQRT,21000,1.449702706,14485.729,0
opencl-e2e,Tesla P40,GPU,w8,p1024,1024,ISQRT,21000,1.451299929,14469.786,0
opencl-kernel,Tesla P40,GPU,w16,p1024,1024,ISQRT,21000,0.303566359,69177.626,0
opencl-e2e,Tesla P40,GPU,w16,p1024,1024,ISQRT,21000,0.305559645,68726.353,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,p1024,1024,MODMUL_R2,12500,0.006062673,2061796.855,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,p1024,1024,MODMUL_R2,12500,0.040179716,311102.247,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,p1024,1024,MODMUL_R2,12500,0.025580916,488645.523,0
library,Tesla P40,gpu,cgbn,p1024,1024,MODMUL_R2,50000,0.000397312,125845682.990,0
opencl-kernel,Tesla P40,GPU,w8,p1024,1024,MODMUL_R2,21000,0.004672076,4494789.812,0
opencl-e2e,Tesla P40,GPU,w8,p1024,1024,MODMUL_R2,21000,0.006738822,3116271.668,0
opencl-kernel,Tesla P40,GPU,w16,p1024,1024,MODMUL_R2,21000,0.000627145,33485082.501,0
opencl-e2e,Tesla P40,GPU,w16,p1024,1024,MODMUL_R2,21000,0.002638313,7959631.816,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,p2048,2048,ADD,6250,0.000220213,28381631.348,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,p2048,2048,ADD,6250,0.055395694,112824.653,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,p2048,2048,ADD,6250,0.069323934,90156.453,0
library,Tesla P40,gpu,cgbn,p2048,2048,ADD,50000,0.000180192,277481797.194,0
opencl-kernel,Tesla P40,GPU,w8,p2048,2048,ADD,21000,0.000608516,34510189.648,0
opencl-e2e,Tesla P40,GPU,w8,p2048,2048,ADD,21000,0.004537322,4628280.640,0
opencl-kernel,Tesla P40,GPU,w16,p2048,2048,ADD,21000,0.000313030,67086202.441,0
opencl-e2e,Tesla P40,GPU,w16,p2048,2048,ADD,21000,0.004083825,5142237.922,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,p2048,2048,SUBTRACT,6250,0.000205005,30487053.345,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,p2048,2048,SUBTRACT,6250,0.067307481,92857.435,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,p2048,2048,SUBTRACT,6250,0.051528149,121292.927,0
library,Tesla P40,gpu,cgbn,p2048,2048,SUBTRACT,50000,0.000183104,273068857.043,0
opencl-kernel,Tesla P40,GPU,w8,p2048,2048,SUBTRACT,21000,0.000606731,34611724.911,0
opencl-e2e,Tesla P40,GPU,w8,p2048,2048,SUBTRACT,21000,0.004398522,4774331.049,0
opencl-kernel,Tesla P40,GPU,w16,p2048,2048,SUBTRACT,21000,0.000307809,68224104.760,0
opencl-e2e,Tesla P40,GPU,w16,p2048,2048,SUBTRACT,21000,0.003875040,5419298.595,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,p2048,2048,ADDMOD,6250,0.000528446,11827132.370,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,p2048,2048,ADDMOD,6250,0.024436550,255764.419,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,p2048,2048,ADDMOD,6250,0.079980085,78144.453,0
library,Tesla P40,gpu,cgbn,p2048,2048,ADDMOD,50000,0.000180224,277432528.409,0
opencl-kernel,Tesla P40,GPU,w8,p2048,2048,ADDMOD,21000,0.000835400,25137663.364,0
opencl-e2e,Tesla P40,GPU,w8,p2048,2048,ADDMOD,21000,0.004562459,4602780.899,0
opencl-kernel,Tesla P40,GPU,w16,p2048,2048,ADDMOD,21000,0.000417232,50331732.261,0
opencl-e2e,Tesla P40,GPU,w16,p2048,2048,ADDMOD,21000,0.003690125,5690864.200,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,p2048,2048,SUBTRACTMOD,6250,0.000420689,14856576.377,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,p2048,2048,SUBTRACTMOD,6250,0.029691619,210497.111,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,p2048,2048,SUBTRACTMOD,6250,0.028906132,216217.098,0
library,Tesla P40,gpu,cgbn,p2048,2048,SUBTRACTMOD,50000,0.000182272,274315308.989,0
opencl-kernel,Tesla P40,GPU,w8,p2048,2048,SUBTRACTMOD,21000,0.000851136,24672909.657,0
opencl-e2e,Tesla P40,GPU,w8,p2048,2048,SUBTRACTMOD,21000,0.004207769,4990768.334,0
opencl-kernel,Tesla P40,GPU,w16,p2048,2048,SUBTRACTMOD,21000,0.000438147,47929120.592,0
opencl-e2e,Tesla P40,GPU,w16,p2048,2048,SUBTRACTMOD,21000,0.004475602,4692106.278,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.003110438,2009363.385,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.049179072,127086.579,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.025245260,247571.227,0
opencl-kernel,Tesla P40,GPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,21000,0.479311020,43812.888,0
opencl-e2e,Tesla P40,GPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,21000,0.486536973,43162.187,0
opencl-kernel,Tesla P40,GPU,w16,p2048,2048,MULTIPLYOPERANDSCANNING,21000,0.148011329,141881.031,0
opencl-e2e,Tesla P40,GPU,w16,p2048,2048,MULTIPLYOPERANDSCANNING,21000,0.153569472,136745.928,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.003397048,1839832.710,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.031273234,199851.412,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.071155012,87836.399,0
library,Tesla P40,gpu,cgbn,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.000447488,111734839.817,0
opencl-kernel,Tesla P40,GPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,21000,0.079103472,265475.073,0
opencl-e2e,Tesla P40,GPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,21000,0.084156183,249536.032,0
opencl-kernel,Tesla P40,GPU,w16,p2048,2048,MULTIPLYPRODUCTSCANNING,21000,0.020061854,1046762.675,0
opencl-e2e,Tesla P40,GPU,w16,p2048,2048,MULTIPLYPRODUCTSCANNING,21000,0.025528608,822606.540,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.021399397,292064.305,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.052149005,119848.883,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.080347490,77787.122,0
library,Tesla P40,gpu,cgbn,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000776192,64417051.451,0
opencl-kernel,Tesla P40,GPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,21000,0.112246063,187088.967,0
opencl-e2e,Tesla P40,GPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,21000,0.115858415,181255.716,0
opencl-kernel,Tesla P40,GPU,w16,p2048,2048,MONTGOMERYMULTIPLICATION,21000,0.004077427,5150306.980,0
opencl-e2e,Tesla P40,GPU,w16,p2048,2048,MONTGOMERYMULTIPLICATION,21000,0.008046199,2609928.035,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,p2048,2048,COMPARE,6250,0.000082523,75736372.921,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,p2048,2048,COMPARE,6250,0.065873304,94879.103,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,p2048,2048,COMPARE,6250,0.066683862,93725.825,0
library,Tesla P40,gpu,cgbn,p2048,2048,COMPARE,50000,0.000181248,275865112.994,0
opencl-kernel,Tesla P40,GPU,w8,p2048,2048,COMPARE,21000,0.000333259,63014074.414,0
opencl-e2e,Tesla P40,GPU,w8,p2048,2048,COMPARE,21000,0.004299915,4883817.550,0
opencl-kernel,Tesla P40,GPU,w16,p2048,2048,COMPARE,21000,0.000173647,120934868.332,0
opencl-e2e,Tesla P40,GPU,w16,p2048,2048,COMPARE,21000,0.004144330,5067164.285,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,p2048,2048,REDUCE,781,0.000026869,29066822.570,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,p2048,2048,REDUCE,781,0.081196378,9618.656,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,p2048,2048,REDUCE,781,0.095932051,8141.179,0
library,Tesla P40,gpu,cgbn,p2048,2048,REDUCE,50000,0.000179200,279017857.143,0
opencl-kernel,Tesla P40,GPU,w8,p2048,2048,REDUCE,21000,0.745434577,28171.486,0
opencl-e2e,Tesla P40,GPU,w8,p2048,2048,REDUCE,21000,0.749813420,28006.967,0
opencl-kernel,Tesla P40,GPU,w16,p2048,2048,REDUCE,21000,0.008647986,2428311.057,0
opencl-e2e,Tesla P40,GPU,w16,p2048,2048,REDUCE,21000,0.012621067,1663884.678,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,p2048,2048,MODMUL,390,0.000618057,631009.667,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,p2048,2048,MODMUL,390,0.087995581,4432.041,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,p2048,2048,MODMUL,390,0.079282054,4919.146,0
library,Tesla P40,gpu,cgbn,p2048,2048,MODMUL,50000,0.002637824,18955017.469,0
opencl-kernel,Tesla P40,GPU,w8,p2048,2048,MODMUL,21000,1.758996411,11938.626,0
opencl-e2e,Tesla P40,GPU,w8,p2048,2048,MODMUL,21000,1.770033887,11864.180,0
opencl-kernel,Tesla P40,GPU,w16,p2048,2048,MODMUL,21000,0.032212003,651930.899,0
opencl-e2e,Tesla P40,GPU,w16,p2048,2048,MODMUL,21000,0.036126473,581291.177,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,p2048,2048,MODEXP,97,0.265019234,366.011,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,p2048,2048,MODEXP,97,0.064131889,1512.508,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,p2048,2048,MODEXP,97,0.119904945,808.974,0
library,Tesla P40,gpu,cgbn,p2048,2048,MODEXP,50000,4.239257812,11794.517,0
opencl-kernel,Tesla P40,GPU,w8,p2048,2048,MODEXP,21000,173.938639484,120.732,0
opencl-e2e,Tesla P40,GPU,w8,p2048,2048,MODEXP,21000,173.962065360,120.716,0
opencl-kernel,Tesla P40,GPU,w16,p2048,2048,MODEXP,21000,50.598355256,415.033,0
opencl-e2e,Tesla P40,GPU,w16,p2048,2048,MODEXP,21000,50.642284648,414.673,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,p2048,2048,EXPONENTIATION,97,0.025023578,3876.344,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,p2048,2048,EXPONENTIATION,97,0.030512018,3179.075,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,p2048,2048,EXPONENTIATION,97,0.144024953,673.494,0
opencl-kernel,Tesla P40,GPU,w8,p2048,2048,EXPONENTIATION,21000,92.038802515,228.165,0
opencl-e2e,Tesla P40,GPU,w8,p2048,2048,EXPONENTIATION,21000,92.009315796,228.238,0
opencl-kernel,Tesla P40,GPU,w16,p2048,2048,EXPONENTIATION,21000,20.461915908,1026.297,0
opencl-e2e,Tesla P40,GPU,w16,p2048,2048,EXPONENTIATION,21000,20.468722883,1025.956,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,p2048,2048,DIVIDE,781,0.000059784,13063712.498,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,p2048,2048,DIVIDE,781,0.104013876,7508.614,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,p2048,2048,DIVIDE,781,0.039969512,19539.893,0
library,Tesla P40,gpu,cgbn,p2048,2048,DIVIDE,50000,0.000183904,271880981.382,0
opencl-kernel,Tesla P40,GPU,w8,p2048,2048,DIVIDE,21000,1.897806430,11065.407,0
opencl-e2e,Tesla P40,GPU,w8,p2048,2048,DIVIDE,21000,1.907009476,11012.006,0
opencl-kernel,Tesla P40,GPU,w16,p2048,2048,DIVIDE,21000,0.643780806,32619.798,0
opencl-e2e,Tesla P40,GPU,w16,p2048,2048,DIVIDE,21000,0.649551456,32330.002,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,p2048,2048,ISQRT,195,0.000165444,1178645.837,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,p2048,2048,ISQRT,195,0.109974664,1773.136,0
opencl-kernel,Tesla P40,GPU,w8,p2048,2048,ISQRT,21000,26.413884440,795.036,0
opencl-e2e,Tesla P40,GPU,w8,p2048,2048,ISQRT,21000,26.431438482,794.508,0
opencl-kernel,Tesla P40,GPU,w16,p2048,2048,ISQRT,21000,0.000000000,inf,0
opencl-e2e,Tesla P40,GPU,w16,p2048,2048,ISQRT,21000,0.000000000,inf,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-1t,p2048,2048,MODMUL_R2,6250,0.009206081,678899.085,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,gmp-nt,p2048,2048,MODMUL_R2,6250,0.047294431,132150.865,0
library,Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz,host-cpu,openssl-nt,p2048,2048,MODMUL_R2,6250,0.068336993,91458.516,0
library,Tesla P40,gpu,cgbn,p2048,2048,MODMUL_R2,50000,0.001475328,33890768.697,0
opencl-kernel,Tesla P40,GPU,w8,p2048,2048,MODMUL_R2,21000,0.068286636,307527.229,0
opencl-e2e,Tesla P40,GPU,w8,p2048,2048,MODMUL_R2,21000,0.071789389,292522.339,0
```
