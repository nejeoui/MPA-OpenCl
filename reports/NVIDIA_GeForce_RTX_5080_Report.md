# MPA-OpenCL benchmark report - NVIDIA GeForce RTX 5080

> **Note.** The multi-threaded GMP and OpenSSL baseline columns have been
> removed from this report: they predate the 2026-09-12 timing fix and were
> understated (see `reports/README.md`). The single-threaded GMP column, the
> OpenCL-on-CPU rows and all MPA measurements are unaffected and were verified
> against GMP before timing.


> **Partial report.** The run was interrupted or hit its time budget.
> Rows that never ran are marked `n/a`.

## 1. System under test

2 OpenCL device(s) exercised with the identical kernels and operands.

### Device 0 - NVIDIA GeForce RTX 5080 (GPU)

| Property | Value |
|---|---|
| Model | NVIDIA GeForce RTX 5080 |
| Type | GPU |
| Vendor | NVIDIA Corporation |
| Device memory | 15.45 GiB |
| Max single allocation | 3.86 GiB |
| Local memory | 48 KiB |
| Global cache | 2688 KiB |
| Compute units | 84 |
| Max clock | 2730 MHz |
| Max work-group size | 1024 |
| OpenCL version | OpenCL 3.0 CUDA |
| Driver | 595.84 |

### Device 1 - cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T (CPU)

| Property | Value |
|---|---|
| Model | cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T |
| Type | CPU |
| Vendor | GenuineIntel |
| Device memory | 29.09 GiB |
| Max single allocation | 8.00 GiB |
| Local memory | 1280 KiB |
| Global cache | 24576 KiB |
| Compute units | 20 |
| Max clock | 4600 MHz |
| Max work-group size | 4096 |
| OpenCL version | OpenCL 3.0 PoCL HSTR: cpu-x86_64-pc-linux-gnu-haswell |
| Driver | 5.0+debian |

### Host

| Property | Value |
|---|---|
| CPU | 13th Gen Intel(R) Core(TM) i5-13500T |
| Logical cores | 20 |
| OpenMP threads used | 20 |
| RAM | 31.1 GB |
| OS | Ubuntu 24.04.4 LTS |
| Kernel | 7.0.0-31-generic |
| Arch | x86_64 |
| GMP | 6.3.0 |
| OpenSSL | OpenSSL 3.0.13 30 Jan 2024 |
| CGBN | cgbn_results.tsv loaded |

## 2. Method

- Workload auto-sized from the device and host: --min-items from 700 x compute units, --items from ten times that capped by host RAM. Either flag, given explicitly, overrides its half.
- Base workload 50000 items, scaled down per operator by its cost weight and by modulus size. Device rows honour --min-items (58800) so the GPU is not left idle; the CPU libraries keep the smaller count because a full-width MODEXP there costs minutes. Both counts appear in every row as dev/cpu, and throughput is per-second so they remain comparable.
- 5 timed repetitions, **minimum** reported. Two untimed warm-up launches precede them.
- `kernel` times `clEnqueueNDRangeKernel` + `clFinish` only. `e2e` adds the host->device operand writes and the device->host result read.
- Every OpenCL device runs the same kernels on the same operands, so GPU and CPU-OpenCL columns are directly comparable.
- CPU library baselines (GMP, OpenSSL) run those same operands, with temporaries preallocated outside the timed region, so the figure is the arithmetic and not marshalling. The generator is reseeded per modulus and operation so every backend sees identical inputs.
- OpenSSL rows time the nearest BN primitive, which is not always semantically identical (its Montgomery routine expects Montgomery-domain inputs); they measure comparable work, not identical results. Correctness is judged against GMP only.
- Every device cell is checked word-for-word against GMP before it is timed. A cell that mismatches is reported and excluded from the speedup tables.
- Total wall time 8692.0 s.

## 3. Correctness

| Device | Kernel | Configs run | Passed | Mismatched | Build/launch failed |
|---|---|---|---|---|---|
| [0] GPU | `mpaKernels_8bits.cl` (w8) | 75 | 75 | 0 | 0 |
| [0] GPU | `mpaKernel_16bits.cl` (w16) | 75 | 75 | 0 | 0 |
| [0] GPU | `mpaKernel_32bits.cl` (w32) | 35 | 35 | 0 | 0 |
| [0] GPU | `mpaKernel_32bits_opt.cl` (w32-opt) | 75 | 75 | 0 | 0 |
| [0] GPU | `mpaKernel_32bits_opt.cl` (w32-o64) | 75 | 75 | 0 | 0 |
| [0] GPU | `mpaKernel_32bits_opt.cl` (w32-il) | 75 | 75 | 0 | 0 |
| [0] GPU | `mpaKernel_32bits_opt.cl` (w32-il64) | 75 | 75 | 0 | 0 |
| [1] CPU | `mpaKernels_8bits.cl` (w8) | 71 | 71 | 0 | 0 |

**All configurations correct** - 556 configurations, 0 problems.

## 4. Throughput per device

Operations per second, higher is better. Kernel-only timings.

### Device 0 - NVIDIA GeForce RTX 5080 (GPU)

#### secp256k1 (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | 1.23 G | 1.92 G | 3.62 G | 3.66 G | 3.76 G | 6.12 G | 5.21 G | 77.99 M | 6.59 G |
| SUBTRACT | 50000 / 50000 | 1.22 G | 2.05 G | 3.51 G | 3.66 G | 3.91 G | 6.42 G | 5.35 G | 107.67 M | 7.23 G |
| ADDMOD | 50000 / 50000 | 817.17 M | 1.51 G | 3.13 G | 4.89 G | 4.87 G | 6.55 G | 6.18 G | 33.24 M | 6.88 G |
| SUBTRACTMOD | 50000 / 50000 | 827.29 M | 1.50 G | 3.23 G | 4.74 G | 4.68 G | 6.44 G | 6.23 G | 36.57 M | 6.15 G |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | 36.35 M | 118.68 M | 585.96 M | 2.48 G | 2.55 G | 5.10 G | 5.09 G | 72.02 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | 247.65 M | 779.14 M | 2.21 G | 2.19 G | 2.27 G | 3.08 G | 2.55 G | 71.92 M | 6.94 G |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | 485.75 M | 1.69 G | 3.80 G | 3.69 G | 4.51 G | 3.92 G | 4.28 G | 8.35 M | 5.39 G |
| COMPARE | 50000 / 50000 | 1.27 G | 2.22 G | - | 4.72 G | 4.90 G | 6.15 G | 6.12 G | 229.66 M | 7.17 G |
| REDUCE | 50000 / 6250 | 294.46 M | 468.96 M | - | 1.41 G | 1.50 G | 1.45 G | 1.19 G | 80.97 M | 4.04 G |
| MODMUL | 50000 / 3125 | 107.63 M | 194.55 M | - | 483.73 M | 475.18 M | 485.09 M | 398.08 M | 15.61 M | 1.34 G |
| MODEXP | 50000 / 781 | 2.34 M | 12.95 M | - | 16.71 M | 26.85 M | 16.98 M | 26.15 M | 125.16 k | 4.43 M |
| EXPONENTIATION | 50000 / 781 | 1.18 M | 3.59 M | - | 81.33 M | 94.48 M | 78.67 M | 87.91 M | 372.36 k | n/a |
| DIVIDE | 50000 / 6250 | 139.57 M | 207.17 M | - | 543.26 M | 553.28 M | 553.96 M | 602.22 M | 35.45 M | 3.09 G |
| ISQRT | 50000 / 1562 | 11.66 M | 16.10 M | - | 77.13 M | 68.31 M | 83.25 M | 53.53 M | 9.09 M | n/a |
| MODMUL_R2 | 50000 / 50000 | 410.28 M | 1.69 G | - | 2.58 G | 3.15 G | 2.75 G | 3.21 G | 14.37 M | 3.31 G |

#### rsa256(composite) (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | 1.22 G | 2.05 G | 3.56 G | 3.61 G | 3.84 G | 6.17 G | 5.37 G | 66.83 M | 6.38 G |
| SUBTRACT | 50000 / 50000 | 1.24 G | 1.95 G | 3.66 G | 3.71 G | 4.03 G | 6.38 G | 5.37 G | 100.96 M | 7.04 G |
| ADDMOD | 50000 / 50000 | 911.08 M | 1.63 G | 3.45 G | 4.51 G | 4.89 G | 6.52 G | 6.14 G | 36.09 M | 6.13 G |
| SUBTRACTMOD | 50000 / 50000 | 826.76 M | 1.50 G | 3.28 G | 4.79 G | 4.82 G | 6.68 G | 6.24 G | 35.42 M | 6.30 G |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | 36.38 M | 119.29 M | 591.21 M | 2.46 G | 2.53 G | 5.07 G | 5.05 G | 72.25 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | 248.03 M | 760.02 M | 2.22 G | 2.20 G | 2.29 G | 3.05 G | 2.53 G | 72.29 M | 7.17 G |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | 487.89 M | 1.68 G | 3.81 G | 3.69 G | 4.59 G | 3.94 G | 4.25 G | 8.23 M | 5.44 G |
| COMPARE | 50000 / 50000 | 1.27 G | 2.23 G | - | 4.92 G | 5.02 G | 6.35 G | 6.19 G | 230.10 M | 6.33 G |
| REDUCE | 50000 / 6250 | 293.84 M | 464.78 M | - | 1.43 G | 1.49 G | 1.43 G | 1.18 G | 48.61 M | 4.09 G |
| MODMUL | 50000 / 3125 | 107.61 M | 194.39 M | - | 485.08 M | 477.14 M | 484.51 M | 397.56 M | 15.83 M | 1.35 G |
| MODEXP | 50000 / 781 | 2.34 M | 12.97 M | - | 16.73 M | 26.95 M | 17.00 M | 26.28 M | 134.03 k | 4.50 M |
| EXPONENTIATION | 50000 / 781 | 1.18 M | 3.61 M | - | 80.69 M | 94.75 M | 77.73 M | 86.30 M | 368.78 k | n/a |
| DIVIDE | 50000 / 6250 | 140.13 M | 199.79 M | - | 494.84 M | 514.09 M | 514.58 M | 556.20 M | 20.84 M | 3.16 G |
| ISQRT | 50000 / 1562 | 11.64 M | 16.10 M | - | 77.18 M | 67.92 M | 82.95 M | 53.57 M | 14.12 M | n/a |
| MODMUL_R2 | 50000 / 50000 | 419.22 M | 1.61 G | - | 2.58 G | 3.25 G | 2.70 G | 3.20 G | 14.41 M | 3.39 G |

#### brainpoolP512r1 (512-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 25000 | 542.09 M | 1.06 G | 1.79 G | 1.98 G | 2.00 G | 3.50 G | 3.44 G | 61.91 M | 5.70 G |
| SUBTRACT | 50000 / 25000 | 546.33 M | 1.06 G | 1.89 G | 1.97 G | 1.97 G | 3.46 G | 3.41 G | 99.43 M | 6.06 G |
| ADDMOD | 50000 / 25000 | 401.45 M | 767.21 M | 1.68 G | 1.94 G | 1.88 G | 3.80 G | 3.73 G | 32.07 M | 5.58 G |
| SUBTRACTMOD | 50000 / 25000 | 341.64 M | 675.99 M | 1.52 G | 1.93 G | 1.95 G | 3.88 G | 3.67 G | 32.60 M | 5.46 G |
| MULTIPLYOPERANDSCANNING | 50000 / 25000 | 6.98 M | 25.02 M | 93.31 M | 677.80 M | 681.12 M | 1.07 G | 1.10 G | 29.81 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 25000 | 35.97 M | 132.13 M | 469.06 M | 360.14 M | 362.25 M | 675.79 M | 717.06 M | 29.84 M | 5.54 G |
| MONTGOMERYMULTIPLICATION | 50000 / 25000 | 145.42 M | 471.39 M | 1.44 G | 1.37 G | 1.80 G | 1.40 G | 1.82 G | 3.39 M | 3.73 G |
| COMPARE | 50000 / 25000 | 516.41 M | 1.02 G | - | 2.44 G | 2.76 G | 5.27 G | 5.19 G | 224.38 M | 5.79 G |
| REDUCE | 50000 / 3125 | 96.24 M | 126.47 M | - | 487.65 M | 504.74 M | 419.64 M | 443.08 M | 46.62 M | 2.44 G |
| MODMUL | 50000 / 1562 | 33.53 M | 48.78 M | - | 130.09 M | 130.64 M | 106.20 M | 110.09 M | 7.09 M | 443.51 M |
| MODEXP | 50000 / 390 | 155.18 k | 1.85 M | - | 2.11 M | 3.73 M | 2.05 M | 3.63 M | 26.11 k | 1.64 M |
| EXPONENTIATION | 50000 / 390 | 150.58 k | 613.25 k | - | 2.23 M | 2.32 M | 2.25 M | 2.30 M | 110.45 k | n/a |
| DIVIDE | 50000 / 3125 | 48.16 M | 51.40 M | - | 162.54 M | 167.37 M | 170.82 M | 174.82 M | 31.81 M | 1.72 G |
| ISQRT | 50000 / 781 | 2.24 M | 2.39 M | - | 12.24 M | 12.51 M | 9.73 M | 9.73 M | 7.32 M | n/a |
| MODMUL_R2 | 50000 / 25000 | 87.98 M | 453.84 M | - | 793.41 M | 1.23 G | 784.44 M | 1.11 G | 6.36 M | 2.12 G |

#### p1024 (1024-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 12500 | 207.00 M | 403.54 M | 1.06 G | 1.04 G | 1.07 G | 2.25 G | 2.24 G | 52.18 M | 3.69 G |
| SUBTRACT | 50000 / 12500 | 207.66 M | 401.08 M | 1.05 G | 1.06 G | 1.07 G | 2.28 G | 2.15 G | 78.43 M | 3.50 G |
| ADDMOD | 50000 / 12500 | 134.26 M | 267.35 M | 862.65 M | 805.30 M | 816.05 M | 2.26 G | 2.25 G | 23.22 M | 3.50 G |
| SUBTRACTMOD | 50000 / 12500 | 136.68 M | 271.02 M | 865.83 M | 806.67 M | 807.57 M | 2.39 G | 2.28 G | 27.29 M | 3.29 G |
| MULTIPLYOPERANDSCANNING | 50000 / 12500 | 1.62 M | 6.04 M | 32.55 M | 249.64 M | 249.10 M | 349.62 M | 347.80 M | 8.20 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 12500 | 4.79 M | 18.72 M | 72.77 M | 72.63 M | 72.62 M | 163.72 M | 174.85 M | 8.24 M | 2.15 G |
| MONTGOMERYMULTIPLICATION | 50000 / 12500 | 15.41 M | 139.37 M | 406.03 M | 341.87 M | 454.17 M | 455.47 M | 660.15 M | 1.14 M | 1.41 G |
| COMPARE | 50000 / 12500 | 214.39 M | 417.00 M | - | 1.35 G | 1.40 G | 3.77 G | 3.50 G | 176.72 M | 3.35 G |
| REDUCE | 50000 / 1562 | 17.80 M | 32.90 M | - | 129.28 M | 136.06 M | 139.87 M | 152.78 M | 42.22 M | 1.63 G |
| MODMUL | 50000 / 781 | 4.13 M | 12.57 M | - | 29.73 M | 29.92 M | 30.30 M | 30.39 M | 2.58 M | 166.92 M |
| MODEXP | 50000 / 195 | 10.02 k | 88.02 k | - | 269.20 k | 442.14 k | 270.48 k | 442.58 k | 4.55 k | 325.58 k |
| EXPONENTIATION | 50000 / 195 | 19.91 k | 82.86 k | - | 313.89 k | 345.24 k | 320.54 k | 333.98 k | 28.91 k | n/a |
| DIVIDE | 50000 / 1562 | 1.98 M | 8.41 M | - | 38.50 M | 41.17 M | 39.40 M | 41.57 M | 14.19 M | 1.19 G |
| ISQRT | 50000 / 390 | 138.89 k | 337.54 k | - | 1.75 M | 1.88 M | 1.81 M | 1.84 M | 4.33 M | n/a |
| MODMUL_R2 | 50000 / 12500 | 8.46 M | 103.93 M | - | 206.01 M | 287.92 M | 243.31 M | 366.73 M | 2.15 M | 757.39 M |

#### p2048 (2048-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 6250 | 101.72 M | 198.17 M | 564.08 M | 568.69 M | 522.11 M | 1.31 G | 1.16 G | 38.83 M | 1.98 G |
| SUBTRACT | 50000 / 6250 | 96.73 M | 198.01 M | 566.16 M | 563.10 M | 530.07 M | 1.31 G | 1.15 G | 54.26 M | 2.02 G |
| ADDMOD | 50000 / 6250 | 73.54 M | 145.78 M | 458.39 M | 421.57 M | 405.93 M | 1.35 G | 1.16 G | 16.88 M | 2.03 G |
| SUBTRACTMOD | 50000 / 6250 | 69.48 M | 133.16 M | 422.46 M | 421.89 M | 403.54 M | 1.36 G | 1.13 G | 19.98 M | 2.05 G |
| MULTIPLYOPERANDSCANNING | 50000 / 6250 | 386.97 k | 1.46 M | 7.53 M | 78.73 M | 79.00 M | 98.50 M | 95.63 M | 2.47 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 6250 | 1.21 M | 4.77 M | 18.87 M | 18.60 M | 18.58 M | 30.70 M | 32.50 M | 2.47 M | 675.53 M |
| MONTGOMERYMULTIPLICATION | 50000 / 6250 | 1.05 M | 25.12 M | 124.86 M | 109.59 M | 131.73 M | 128.74 M | 157.41 M | 341.89 k | 433.67 M |
| COMPARE | 50000 / 6250 | 106.41 M | 221.06 M | - | 783.99 M | 706.63 M | 2.29 G | 1.90 G | 161.25 M | 2.03 G |
| REDUCE | 50000 / 781 | 211.92 k | 8.90 M | - | 42.34 M | 49.13 M | 43.95 M | 51.24 M | 40.06 M | 1.61 G |
| MODMUL | 50000 / 390 | 102.97 k | 2.59 M | - | 8.23 M | 3.83 M | 8.30 M | 3.85 M | 773.40 k | 87.78 M |
| MODEXP | 50000 / 97 | 771.1 | 4.15 k | - | 34.92 k | 17.14 k | 34.85 k | 17.07 k | 654.0 | 60.97 k |
| EXPONENTIATION | 50000 / 97 | 2.07 k | 10.27 k | - | 42.31 k | 45.58 k | 41.99 k | 47.67 k | 5.21 k | n/a |
| DIVIDE | 50000 / 781 | 40.68 k | 497.17 k | - | 3.27 M | 3.42 M | 3.34 M | 3.55 M | 11.25 M | 1.27 G |
| ISQRT | 50000 / 195 | 2.91 k | 17.13 k | - | 280.28 k | 657.86 k | 276.47 k | 665.34 k | 2.73 M | n/a |
| MODMUL_R2 | 50000 / 6250 | 1.48 M | 25.59 M | - | 61.09 M | 76.54 M | 65.93 M | 85.85 M | 663.68 k | 225.18 M |

### Device 1 - cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T (CPU)

#### secp256k1 (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | 258.48 M | - | - | - | - | - | - | 77.99 M | 6.59 G |
| SUBTRACT | 50000 / 50000 | 253.10 M | - | - | - | - | - | - | 107.67 M | 7.23 G |
| ADDMOD | 50000 / 50000 | 182.27 M | - | - | - | - | - | - | 33.24 M | 6.88 G |
| SUBTRACTMOD | 50000 / 50000 | 172.89 M | - | - | - | - | - | - | 36.57 M | 6.15 G |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | 9.79 M | - | - | - | - | - | - | 72.02 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | 10.67 M | - | - | - | - | - | - | 71.92 M | 6.94 G |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | 2.02 M | - | - | - | - | - | - | 8.35 M | 5.39 G |
| COMPARE | 50000 / 50000 | 522.03 M | - | - | - | - | - | - | 229.66 M | 7.17 G |
| REDUCE | 50000 / 6250 | 2.96 M | - | - | - | - | - | - | 80.97 M | 4.04 G |
| MODMUL | 50000 / 3125 | 646.38 k | - | - | - | - | - | - | 15.61 M | 1.34 G |
| MODEXP | 50000 / 781 | 7.95 k | - | - | - | - | - | - | 125.16 k | 4.43 M |
| EXPONENTIATION | 50000 / 781 | 19.46 k | - | - | - | - | - | - | 372.36 k | n/a |
| DIVIDE | 50000 / 6250 | 2.01 M | - | - | - | - | - | - | 35.45 M | 3.09 G |
| ISQRT | 50000 / 1562 | 111.52 k | - | - | - | - | - | - | 9.09 M | n/a |
| MODMUL_R2 | 50000 / 50000 | 1.98 M | - | - | - | - | - | - | 14.37 M | 3.31 G |

#### rsa256(composite) (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | 316.63 M | - | - | - | - | - | - | 66.83 M | 6.38 G |
| SUBTRACT | 50000 / 50000 | 288.92 M | - | - | - | - | - | - | 100.96 M | 7.04 G |
| ADDMOD | 50000 / 50000 | 211.88 M | - | - | - | - | - | - | 36.09 M | 6.13 G |
| SUBTRACTMOD | 50000 / 50000 | 172.67 M | - | - | - | - | - | - | 35.42 M | 6.30 G |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | 10.43 M | - | - | - | - | - | - | 72.25 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | 10.95 M | - | - | - | - | - | - | 72.29 M | 7.17 G |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | 2.14 M | - | - | - | - | - | - | 8.23 M | 5.44 G |
| COMPARE | 50000 / 50000 | 484.31 M | - | - | - | - | - | - | 230.10 M | 6.33 G |
| REDUCE | 50000 / 6250 | 2.96 M | - | - | - | - | - | - | 48.61 M | 4.09 G |
| MODMUL | 50000 / 3125 | 787.75 k | - | - | - | - | - | - | 15.83 M | 1.35 G |
| MODEXP | 50000 / 781 | 7.94 k | - | - | - | - | - | - | 134.03 k | 4.50 M |
| EXPONENTIATION | 50000 / 781 | 19.84 k | - | - | - | - | - | - | 368.78 k | n/a |
| DIVIDE | 50000 / 6250 | 2.07 M | - | - | - | - | - | - | 20.84 M | 3.16 G |
| ISQRT | 50000 / 1562 | 110.36 k | - | - | - | - | - | - | 14.12 M | n/a |
| MODMUL_R2 | 50000 / 50000 | 2.11 M | - | - | - | - | - | - | 14.41 M | 3.39 G |

#### brainpoolP512r1 (512-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 25000 | 161.73 M | - | - | - | - | - | - | 61.91 M | 5.70 G |
| SUBTRACT | 50000 / 25000 | 149.80 M | - | - | - | - | - | - | 99.43 M | 6.06 G |
| ADDMOD | 50000 / 25000 | 119.71 M | - | - | - | - | - | - | 32.07 M | 5.58 G |
| SUBTRACTMOD | 50000 / 25000 | 97.19 M | - | - | - | - | - | - | 32.60 M | 5.46 G |
| MULTIPLYOPERANDSCANNING | 50000 / 25000 | 2.42 M | - | - | - | - | - | - | 29.81 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 25000 | 2.67 M | - | - | - | - | - | - | 29.84 M | 5.54 G |
| MONTGOMERYMULTIPLICATION | 50000 / 25000 | 478.93 k | - | - | - | - | - | - | 3.39 M | 3.73 G |
| COMPARE | 50000 / 25000 | 416.07 M | - | - | - | - | - | - | 224.38 M | 5.79 G |
| REDUCE | 50000 / 3125 | 784.63 k | - | - | - | - | - | - | 46.62 M | 2.44 G |
| MODMUL | 50000 / 1562 | 169.86 k | - | - | - | - | - | - | 7.09 M | 443.51 M |
| MODEXP | 50000 / 390 | 957.9 | - | - | - | - | - | - | 26.11 k | 1.64 M |
| EXPONENTIATION | 50000 / 390 | 2.20 k | - | - | - | - | - | - | 110.45 k | n/a |
| DIVIDE | 50000 / 3125 | 450.24 k | - | - | - | - | - | - | 31.81 M | 1.72 G |
| ISQRT | 50000 / 781 | 25.01 k | - | - | - | - | - | - | 7.32 M | n/a |
| MODMUL_R2 | 50000 / 25000 | 463.81 k | - | - | - | - | - | - | 6.36 M | 2.12 G |

#### p1024 (1024-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 12500 | 80.36 M | - | - | - | - | - | - | 52.18 M | 3.69 G |
| SUBTRACT | 50000 / 12500 | 74.12 M | - | - | - | - | - | - | 78.43 M | 3.50 G |
| ADDMOD | 50000 / 12500 | 46.76 M | - | - | - | - | - | - | 23.22 M | 3.50 G |
| SUBTRACTMOD | 50000 / 12500 | 43.15 M | - | - | - | - | - | - | 27.29 M | 3.29 G |
| MULTIPLYOPERANDSCANNING | 50000 / 12500 | 526.14 k | - | - | - | - | - | - | 8.20 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 12500 | 785.81 k | - | - | - | - | - | - | 8.24 M | 2.15 G |
| MONTGOMERYMULTIPLICATION | 50000 / 12500 | 104.74 k | - | - | - | - | - | - | 1.14 M | 1.41 G |
| COMPARE | 50000 / 12500 | 147.57 M | - | - | - | - | - | - | 176.72 M | 3.35 G |
| REDUCE | 50000 / 1562 | 203.73 k | - | - | - | - | - | - | 42.22 M | 1.63 G |
| MODMUL | 50000 / 781 | 38.72 k | - | - | - | - | - | - | 2.58 M | 166.92 M |
| MODEXP | 50000 / 195 | over budget | - | - | - | - | - | - | 4.55 k | 325.58 k |
| EXPONENTIATION | 50000 / 195 | over budget | - | - | - | - | - | - | 28.91 k | n/a |
| DIVIDE | 50000 / 1562 | 101.42 k | - | - | - | - | - | - | 14.19 M | 1.19 G |
| ISQRT | 50000 / 390 | 5.39 k | - | - | - | - | - | - | 4.33 M | n/a |
| MODMUL_R2 | 50000 / 12500 | 97.39 k | - | - | - | - | - | - | 2.15 M | 757.39 M |

#### p2048 (2048-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 6250 | 25.76 M | - | - | - | - | - | - | 38.83 M | 1.98 G |
| SUBTRACT | 50000 / 6250 | 23.49 M | - | - | - | - | - | - | 54.26 M | 2.02 G |
| ADDMOD | 50000 / 6250 | 22.25 M | - | - | - | - | - | - | 16.88 M | 2.03 G |
| SUBTRACTMOD | 50000 / 6250 | 20.86 M | - | - | - | - | - | - | 19.98 M | 2.05 G |
| MULTIPLYOPERANDSCANNING | 50000 / 6250 | 111.05 k | - | - | - | - | - | - | 2.47 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 6250 | 154.66 k | - | - | - | - | - | - | 2.47 M | 675.53 M |
| MONTGOMERYMULTIPLICATION | 50000 / 6250 | 24.14 k | - | - | - | - | - | - | 341.89 k | 433.67 M |
| COMPARE | 50000 / 6250 | 40.09 M | - | - | - | - | - | - | 161.25 M | 2.03 G |
| REDUCE | 50000 / 781 | 46.90 k | - | - | - | - | - | - | 40.06 M | 1.61 G |
| MODMUL | 50000 / 390 | 9.76 k | - | - | - | - | - | - | 773.40 k | 87.78 M |
| MODEXP | 50000 / 97 | over budget | - | - | - | - | - | - | 654.0 | 60.97 k |
| EXPONENTIATION | 50000 / 97 | - | - | - | - | - | - | - | 5.21 k | n/a |
| DIVIDE | 50000 / 781 | - | - | - | - | - | - | - | 11.25 M | 1.27 G |
| ISQRT | 50000 / 195 | - | - | - | - | - | - | - | 2.73 M | n/a |
| MODMUL_R2 | 50000 / 6250 | - | - | - | - | - | - | - | 663.68 k | 225.18 M |

## 5. Head to head

Best OpenCL GPU result against best OpenCL CPU result and the CPU libraries.
Ratios above 1.00x mean the GPU is faster than that baseline.

### secp256k1 (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w32-il | 6.12 G | w8 | 258.48 M | 77.99 M | 6.59 G | 23.69x | 0.93x |
| SUBTRACT | w32-il | 6.42 G | w8 | 253.10 M | 107.67 M | 7.23 G | 25.35x | 0.89x |
| ADDMOD | w32-il | 6.55 G | w8 | 182.27 M | 33.24 M | 6.88 G | 35.96x | 0.95x |
| SUBTRACTMOD | w32-il | 6.44 G | w8 | 172.89 M | 36.57 M | 6.15 G | 37.22x | 1.05x |
| MULTIPLYOPERANDSCANNING | w32-il | 5.10 G | w8 | 9.79 M | 72.02 M | n/a | 521.24x | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il | 3.08 G | w8 | 10.67 M | 71.92 M | 6.94 G | 288.11x | 0.44x |
| MONTGOMERYMULTIPLICATION | w32-o64 | 4.51 G | w8 | 2.02 M | 8.35 M | 5.39 G | 2232.39x | 0.84x |
| COMPARE | w32-il | 6.15 G | w8 | 522.03 M | 229.66 M | 7.17 G | 11.79x | 0.86x |
| REDUCE | w32-o64 | 187.05 M | w8 | 370.53 k | 80.97 M | 4.04 G | 504.82x | 0.05x |
| MODMUL | w32-il | 30.32 M | w8 | 40.40 k | 15.61 M | 1.34 G | 750.48x | 0.02x |
| MODEXP | w32-o64 | 419.41 k | w8 | 124.2 | 125.16 k | 4.43 M | 3375.97x | 0.09x |
| EXPONENTIATION | w32-o64 | 1.48 M | w8 | 303.9 | 372.36 k | n/a | 4855.67x | n/a |
| DIVIDE | w32-il64 | 75.28 M | w8 | 251.47 k | 35.45 M | 3.09 G | 299.35x | 0.02x |
| ISQRT | w32-il | 2.60 M | w8 | 3.48 k | 9.09 M | n/a | 746.52x | n/a |
| MODMUL_R2 | w32-il64 | 3.21 G | w8 | 1.98 M | 14.37 M | 3.31 G | 1619.70x | 0.97x |

### rsa256(composite) (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w32-il | 6.17 G | w8 | 316.63 M | 66.83 M | 6.38 G | 19.49x | 0.97x |
| SUBTRACT | w32-il | 6.38 G | w8 | 288.92 M | 100.96 M | 7.04 G | 22.07x | 0.91x |
| ADDMOD | w32-il | 6.52 G | w8 | 211.88 M | 36.09 M | 6.13 G | 30.77x | 1.06x |
| SUBTRACTMOD | w32-il | 6.68 G | w8 | 172.67 M | 35.42 M | 6.30 G | 38.67x | 1.06x |
| MULTIPLYOPERANDSCANNING | w32-il | 5.07 G | w8 | 10.43 M | 72.25 M | n/a | 486.50x | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il | 3.05 G | w8 | 10.95 M | 72.29 M | 7.17 G | 278.40x | 0.43x |
| MONTGOMERYMULTIPLICATION | w32-o64 | 4.59 G | w8 | 2.14 M | 8.23 M | 5.44 G | 2149.31x | 0.84x |
| COMPARE | w32-il | 6.35 G | w8 | 484.31 M | 230.10 M | 6.33 G | 13.11x | 1.00x |
| REDUCE | w32-o64 | 186.73 M | w8 | 370.12 k | 48.61 M | 4.09 G | 504.51x | 0.05x |
| MODMUL | w32-opt | 30.32 M | w8 | 49.23 k | 15.83 M | 1.35 G | 615.78x | 0.02x |
| MODEXP | w32-o64 | 421.03 k | w8 | 124.1 | 134.03 k | 4.50 M | 3392.65x | 0.09x |
| EXPONENTIATION | w32-o64 | 1.48 M | w8 | 309.9 | 368.78 k | n/a | 4776.20x | n/a |
| DIVIDE | w32-il64 | 69.52 M | w8 | 259.07 k | 20.84 M | 3.16 G | 268.36x | 0.02x |
| ISQRT | w32-il | 2.59 M | w8 | 3.45 k | 14.12 M | n/a | 751.66x | n/a |
| MODMUL_R2 | w32-o64 | 3.25 G | w8 | 2.11 M | 14.41 M | 3.39 G | 1536.21x | 0.96x |

### brainpoolP512r1 (512-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w32-il | 1.75 G | w8 | 80.86 M | 61.91 M | 5.70 G | 21.62x | 0.31x |
| SUBTRACT | w32-il | 1.73 G | w8 | 74.90 M | 99.43 M | 6.06 G | 23.11x | 0.29x |
| ADDMOD | w32-il | 1.90 G | w8 | 59.86 M | 32.07 M | 5.58 G | 31.71x | 0.34x |
| SUBTRACTMOD | w32-il | 1.94 G | w8 | 48.59 M | 32.60 M | 5.46 G | 39.97x | 0.36x |
| MULTIPLYOPERANDSCANNING | w32-il64 | 548.55 M | w8 | 1.21 M | 29.81 M | n/a | 454.14x | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 358.53 M | w8 | 1.34 M | 29.84 M | 5.54 G | 268.50x | 0.06x |
| MONTGOMERYMULTIPLICATION | w32-il64 | 909.98 M | w8 | 239.46 k | 3.39 M | 3.73 G | 3800.10x | 0.24x |
| COMPARE | w32-il | 2.64 G | w8 | 208.03 M | 224.38 M | 5.79 G | 12.67x | 0.46x |
| REDUCE | w32-o64 | 31.55 M | w8 | 49.04 k | 46.62 M | 2.44 G | 643.29x | 0.01x |
| MODMUL | w32-o64 | 4.08 M | w8 | 5.31 k | 7.09 M | 443.51 M | 769.10x | 0.01x |
| MODEXP | w32-o64 | 29.08 k | w8 | 7.5 | 26.11 k | 1.64 M | 3892.41x | 0.02x |
| EXPONENTIATION | w32-o64 | 18.07 k | w8 | 17.1 | 110.45 k | n/a | 1054.18x | n/a |
| DIVIDE | w32-il64 | 10.93 M | w8 | 28.14 k | 31.81 M | 1.72 G | 388.29x | 0.01x |
| ISQRT | w32-o64 | 195.43 k | w8 | 390.7 | 7.32 M | n/a | 500.19x | n/a |
| MODMUL_R2 | w32-o64 | 615.60 M | w8 | 231.90 k | 6.36 M | 2.12 G | 2654.54x | 0.29x |

### p1024 (1024-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w32-il | 561.44 M | w8 | 20.09 M | 52.18 M | 3.69 G | 27.95x | 0.15x |
| SUBTRACT | w32-il | 570.54 M | w8 | 18.53 M | 78.43 M | 3.50 G | 30.79x | 0.16x |
| ADDMOD | w32-il | 566.07 M | w8 | 11.69 M | 23.22 M | 3.50 G | 48.43x | 0.16x |
| SUBTRACTMOD | w32-il | 597.09 M | w8 | 10.79 M | 27.29 M | 3.29 G | 55.35x | 0.18x |
| MULTIPLYOPERANDSCANNING | w32-il | 87.41 M | w8 | 131.54 k | 8.20 M | n/a | 664.50x | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 43.71 M | w8 | 196.45 k | 8.24 M | 2.15 G | 222.50x | 0.02x |
| MONTGOMERYMULTIPLICATION | w32-il64 | 165.04 M | w8 | 26.18 k | 1.14 M | 1.41 G | 6302.81x | 0.12x |
| COMPARE | w32-il | 943.33 M | w8 | 36.89 M | 176.72 M | 3.35 G | 25.57x | 0.28x |
| REDUCE | w32-il64 | 4.77 M | w8 | 6.36 k | 42.22 M | 1.63 G | 749.93x | 0.00x |
| MODMUL | w32-il64 | 474.67 k | w8 | 604.8 | 2.58 M | 166.92 M | 784.85x | 0.00x |
| MODEXP | w32-il64 | 1.73 k | none | n/a | 4.55 k | 325.58 k | n/a | 0.01x |
| EXPONENTIATION | w32-o64 | 1.35 k | none | n/a | 28.91 k | n/a | n/a | n/a |
| DIVIDE | w32-il64 | 1.30 M | w8 | 3.17 k | 14.19 M | 1.19 G | 409.83x | 0.00x |
| ISQRT | w32-o64 | 14.69 k | w8 | 42.1 | 4.33 M | n/a | 349.39x | n/a |
| MODMUL_R2 | w32-il64 | 91.68 M | w8 | 24.35 k | 2.15 M | 757.39 M | 3765.69x | 0.12x |

### p2048 (2048-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w32-il | 164.02 M | w8 | 3.22 M | 38.83 M | 1.98 G | 50.93x | 0.08x |
| SUBTRACT | w32-il | 163.86 M | w8 | 2.94 M | 54.26 M | 2.02 G | 55.81x | 0.08x |
| ADDMOD | w32-il | 168.49 M | w8 | 2.78 M | 16.88 M | 2.03 G | 60.57x | 0.08x |
| SUBTRACTMOD | w32-il | 169.60 M | w8 | 2.61 M | 19.98 M | 2.05 G | 65.05x | 0.08x |
| MULTIPLYOPERANDSCANNING | w32-il | 12.31 M | w8 | 13.88 k | 2.47 M | n/a | 886.92x | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 4.06 M | w8 | 19.33 k | 2.47 M | 675.53 M | 210.15x | 0.01x |
| MONTGOMERYMULTIPLICATION | w32-il64 | 19.68 M | w8 | 3.02 k | 341.89 k | 433.67 M | 6521.47x | 0.05x |
| COMPARE | w32-il | 286.04 M | w8 | 5.01 M | 161.25 M | 2.03 G | 57.07x | 0.14x |
| REDUCE | w32-il64 | 800.40 k | w8 | 732.6 | 40.06 M | 1.61 G | 1092.58x | 0.00x |
| MODMUL | w32-il | 64.71 k | w8 | 76.1 | 773.40 k | 87.78 M | 850.01x | 0.00x |
| MODEXP | w32-opt | 67.8 | none | n/a | 654.0 | 60.97 k | n/a | 0.00x |
| EXPONENTIATION | w32-il64 | 92.5 | none | n/a | 5.21 k | n/a | n/a | n/a |
| DIVIDE | w32-il64 | 55.40 k | none | n/a | 11.25 M | 1.27 G | n/a | 0.00x |
| ISQRT | w32-il64 | 2.59 k | none | n/a | 2.73 M | n/a | n/a | n/a |
| MODMUL_R2 | w32-il64 | 10.73 M | none | n/a | 663.68 k | 225.18 M | n/a | 0.05x |

## 6. Raw data

Also written to `NVIDIA_GeForce_RTX_5080_Report.csv` for analysis.

```csv
kind,device,device_type,kernel,modulus,bits,operation,items,seconds,ops_per_sec,mismatches
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,secp256k1,256,ADD,50000,0.000641081,77993264.263,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,secp256k1,256,ADD,50000,0.000191078,261673241.908,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,secp256k1,256,ADD,50000,0.000249707,200234676.525,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,secp256k1,256,ADD,50000,0.000007584,6592827004.219,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,secp256k1,256,ADD,50000,0.000040713,1228108982.563,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,secp256k1,256,ADD,50000,0.000727441,68734096.478,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,secp256k1,256,ADD,50000,0.000026070,1917913319.426,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,secp256k1,256,ADD,50000,0.000540795,92456475.918,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32,secp256k1,256,ADD,50000,0.000013817,3618730587.669,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32,secp256k1,256,ADD,50000,0.000489993,102042273.511,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,secp256k1,256,ADD,50000,0.000013643,3664883368.199,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,secp256k1,256,ADD,50000,0.000507757,98472300.301,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,secp256k1,256,ADD,50000,0.000013290,3762226603.393,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,secp256k1,256,ADD,50000,0.000554983,90092849.354,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,secp256k1,256,ADD,50000,0.000008164,6124446704.431,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,secp256k1,256,ADD,50000,0.000486854,102700193.876,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,secp256k1,256,ADD,50000,0.000009604,5206164923.197,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,secp256k1,256,ADD,50000,0.000534227,93593172.680,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,secp256k1,256,ADD,50000,0.000193437,258482093.109,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,secp256k1,256,ADD,50000,0.000573188,87231414.876,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,secp256k1,256,SUBTRACT,50000,0.000464385,107669282.875,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,secp256k1,256,SUBTRACT,50000,0.000324389,154135930.008,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,secp256k1,256,SUBTRACT,50000,0.000331659,150757254.908,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,secp256k1,256,SUBTRACT,50000,0.000006912,7233796296.296,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,secp256k1,256,SUBTRACT,50000,0.000041004,1219393197.076,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,secp256k1,256,SUBTRACT,50000,0.000742739,67318398.762,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,secp256k1,256,SUBTRACT,50000,0.000024449,2045073401.969,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,secp256k1,256,SUBTRACT,50000,0.000537983,92939739.739,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32,secp256k1,256,SUBTRACT,50000,0.000014239,3511482114.580,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32,secp256k1,256,SUBTRACT,50000,0.000487076,102653384.996,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,secp256k1,256,SUBTRACT,50000,0.000013651,3662734654.973,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,secp256k1,256,SUBTRACT,50000,0.000498203,100360696.250,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,secp256k1,256,SUBTRACT,50000,0.000012797,3907166380.829,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,secp256k1,256,SUBTRACT,50000,0.000558029,89601078.392,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,secp256k1,256,SUBTRACT,50000,0.000007794,6415189396.186,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,secp256k1,256,SUBTRACT,50000,0.000489478,102149637.063,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,secp256k1,256,SUBTRACT,50000,0.000009353,5345878460.154,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,secp256k1,256,SUBTRACT,50000,0.000525634,95123222.450,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,secp256k1,256,SUBTRACT,50000,0.000197550,253100478.639,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,secp256k1,256,SUBTRACT,50000,0.000518111,96504417.173,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,secp256k1,256,ADDMOD,50000,0.001504186,33240569.929,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,secp256k1,256,ADDMOD,50000,0.000425234,117582319.263,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,secp256k1,256,ADDMOD,50000,0.000904475,55280687.666,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,secp256k1,256,ADDMOD,50000,0.000007264,6883259911.894,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,secp256k1,256,ADDMOD,50000,0.000061187,817167020.615,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,secp256k1,256,ADDMOD,50000,0.000754565,66263343.734,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,secp256k1,256,ADDMOD,50000,0.000033195,1506250961.382,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,secp256k1,256,ADDMOD,50000,0.000576640,86709212.129,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32,secp256k1,256,ADDMOD,50000,0.000015982,3128519377.141,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32,secp256k1,256,ADDMOD,50000,0.000493582,101300290.578,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,secp256k1,256,ADDMOD,50000,0.000010219,4892847234.847,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,secp256k1,256,ADDMOD,50000,0.000501762,99648837.023,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,secp256k1,256,ADDMOD,50000,0.000010277,4865233474.647,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,secp256k1,256,ADDMOD,50000,0.000485820,102918776.310,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,secp256k1,256,ADDMOD,50000,0.000007629,6553940642.705,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,secp256k1,256,ADDMOD,50000,0.000477871,104630747.843,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,secp256k1,256,ADDMOD,50000,0.000008093,6178178397.698,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,secp256k1,256,ADDMOD,50000,0.000524705,95291639.511,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,secp256k1,256,ADDMOD,50000,0.000274320,182268881.685,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,secp256k1,256,ADDMOD,50000,0.000565873,88359048.669,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,secp256k1,256,SUBTRACTMOD,50000,0.001367062,36574785.949,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,secp256k1,256,SUBTRACTMOD,50000,0.000309056,161782976.716,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,secp256k1,256,SUBTRACTMOD,50000,0.000891428,56089779.506,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,secp256k1,256,SUBTRACTMOD,50000,0.000008128,6151574803.150,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.000060438,827294083.508,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.000801575,62377194.971,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,secp256k1,256,SUBTRACTMOD,50000,0.000033234,1504483419.600,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,secp256k1,256,SUBTRACTMOD,50000,0.000519975,96158468.931,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32,secp256k1,256,SUBTRACTMOD,50000,0.000015485,3228931518.878,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32,secp256k1,256,SUBTRACTMOD,50000,0.000489005,102248443.941,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,secp256k1,256,SUBTRACTMOD,50000,0.000010555,4737090295.681,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,secp256k1,256,SUBTRACTMOD,50000,0.000484842,103126379.438,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,secp256k1,256,SUBTRACTMOD,50000,0.000010678,4682525280.582,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,secp256k1,256,SUBTRACTMOD,50000,0.000486100,102859493.597,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,secp256k1,256,SUBTRACTMOD,50000,0.000007770,6435004437.303,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,secp256k1,256,SUBTRACTMOD,50000,0.000490332,101971725.936,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,secp256k1,256,SUBTRACTMOD,50000,0.000008026,6229753530.870,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,secp256k1,256,SUBTRACTMOD,50000,0.000487947,102470145.369,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.000289197,172892526.666,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.000578611,86413843.140,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000694229,72022344.008,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000221853,225374460.734,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000304308,164307215.267,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001375508,36350206.617,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.002307015,21673027.683,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000421296,118681402.850,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001072675,46612440.834,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000085330,585960384.070,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000775286,64492329.271,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000020173,2478560459.590,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000703405,71082804.408,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000019638,2546084281.367,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000690989,72360052.057,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000009800,5102039846.700,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000668516,74792525.468,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000009817,5093204554.717,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000693623,72085268.392,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.005108105,9788365.743,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.005743229,8705903.947,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000695183,71923507.904,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000263172,189989815.302,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000307967,162355056.310,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000007200,6944444444.444,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000201895,247653483.346,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000947416,52775127.311,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000064173,779143873.277,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000766848,65201969.678,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000022604,2211997950.732,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000725001,68965421.913,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000022811,2191924925.872,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000699895,71439287.203,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000022055,2267059338.586,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000727080,68768223.779,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000016260,3075030589.941,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000686297,72854755.265,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000019572,2554670124.073,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000703251,71098370.260,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.004684683,10673080.759,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.005090699,9821833.897,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.005988899,8348779.963,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000981098,50963308.409,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000378049,132257987.520,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000009280,5387931034.483,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000102933,485752861.426,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000616088,81157237.246,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000029534,1692964063.066,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000539378,92699368.339,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000013156,3800547009.508,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000492325,101558929.971,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000013535,3694126806.276,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000494269,101159490.215,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000011089,4508971539.591,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000486215,102835165.917,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000012745,3923107015.719,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000486731,102726146.110,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000011687,4278258024.503,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000482404,103647564.972,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.024754941,2019798.795,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.022157166,2256606.283,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,secp256k1,256,COMPARE,50000,0.000217717,229655926.670,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,secp256k1,256,COMPARE,50000,0.000144941,344967949.969,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,secp256k1,256,COMPARE,50000,0.000117254,426424686.206,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,secp256k1,256,COMPARE,50000,0.000006976,7167431192.661,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,secp256k1,256,COMPARE,50000,0.000039517,1265278263.654,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,secp256k1,256,COMPARE,50000,0.000607614,82289085.118,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,secp256k1,256,COMPARE,50000,0.000022472,2224991217.046,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,secp256k1,256,COMPARE,50000,0.000517829,96556971.289,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,secp256k1,256,COMPARE,50000,0.000010595,4719207860.370,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,secp256k1,256,COMPARE,50000,0.000494596,101092609.113,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,secp256k1,256,COMPARE,50000,0.000010207,4898600390.779,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,secp256k1,256,COMPARE,50000,0.000479925,104182944.849,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,secp256k1,256,COMPARE,50000,0.000008127,6152331044.485,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,secp256k1,256,COMPARE,50000,0.000496303,100744907.738,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,secp256k1,256,COMPARE,50000,0.000008171,6119203084.914,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,secp256k1,256,COMPARE,50000,0.000476150,105008926.073,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,secp256k1,256,COMPARE,50000,0.000095780,522029649.839,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,secp256k1,256,COMPARE,50000,0.000445336,112274776.400,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,secp256k1,256,REDUCE,6250,0.000077193,80965887.907,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,secp256k1,256,REDUCE,6250,0.000042511,147020767.346,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,secp256k1,256,REDUCE,6250,0.000109214,57227094.911,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,secp256k1,256,REDUCE,50000,0.000012384,4037467700.258,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,secp256k1,256,REDUCE,50000,0.000169804,294457139.517,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,secp256k1,256,REDUCE,50000,0.000758284,65938355.483,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,secp256k1,256,REDUCE,50000,0.000106618,468963975.731,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,secp256k1,256,REDUCE,50000,0.000587195,85150588.588,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,secp256k1,256,REDUCE,50000,0.000035455,1410238357.872,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,secp256k1,256,REDUCE,50000,0.000532901,93826058.061,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,secp256k1,256,REDUCE,50000,0.000033414,1496378662.810,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,secp256k1,256,REDUCE,50000,0.000534929,93470348.327,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,secp256k1,256,REDUCE,50000,0.000034533,1447890539.653,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,secp256k1,256,REDUCE,50000,0.000509840,98069982.793,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,secp256k1,256,REDUCE,50000,0.000041893,1193516729.135,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,secp256k1,256,REDUCE,50000,0.000590111,84729821.779,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,secp256k1,256,REDUCE,50000,0.016867950,2964201.341,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,secp256k1,256,REDUCE,50000,0.017688998,2826615.730,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,secp256k1,256,MODMUL,3125,0.000200137,15614304.390,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,secp256k1,256,MODMUL,3125,0.000059081,52893483.978,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,secp256k1,256,MODMUL,3125,0.000141546,22077628.597,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,secp256k1,256,MODMUL,50000,0.000037376,1337756849.315,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,secp256k1,256,MODMUL,50000,0.000464562,107628260.772,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,secp256k1,256,MODMUL,50000,0.001044230,47882171.514,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,secp256k1,256,MODMUL,50000,0.000257001,194551772.069,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,secp256k1,256,MODMUL,50000,0.000775568,64468879.441,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,secp256k1,256,MODMUL,50000,0.000103364,483727409.476,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,secp256k1,256,MODMUL,50000,0.000614285,81395443.515,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,secp256k1,256,MODMUL,50000,0.000105224,475176763.790,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,secp256k1,256,MODMUL,50000,0.000611532,81761870.280,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,secp256k1,256,MODMUL,50000,0.000103073,485093097.966,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,secp256k1,256,MODMUL,50000,0.000582316,85864032.190,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,secp256k1,256,MODMUL,50000,0.000125603,398079660.772,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,secp256k1,256,MODMUL,50000,0.000713554,70071781.537,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,secp256k1,256,MODMUL,50000,0.077354211,646377.222,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,secp256k1,256,MODMUL,50000,0.067463654,741139.814,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,secp256k1,256,MODEXP,781,0.006239859,125163.085,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,secp256k1,256,MODEXP,781,0.001087410,718220.360,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,secp256k1,256,MODEXP,781,0.001531079,510097.781,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,secp256k1,256,MODEXP,50000,0.011290208,4428616.373,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,secp256k1,256,MODEXP,50000,0.021326678,2344481.405,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,secp256k1,256,MODEXP,50000,0.022052402,2267326.707,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,secp256k1,256,MODEXP,50000,0.003860220,12952629.642,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,secp256k1,256,MODEXP,50000,0.004445019,11248545.845,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,secp256k1,256,MODEXP,50000,0.002992527,16708287.014,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,secp256k1,256,MODEXP,50000,0.003624499,13795010.025,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,secp256k1,256,MODEXP,50000,0.001862153,26850640.095,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,secp256k1,256,MODEXP,50000,0.002486268,20110462.748,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,secp256k1,256,MODEXP,50000,0.002944560,16980465.665,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,secp256k1,256,MODEXP,50000,0.003542520,14114246.354,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,secp256k1,256,MODEXP,50000,0.001912373,26145527.026,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,secp256k1,256,MODEXP,50000,0.002589415,19309380.701,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,secp256k1,256,MODEXP,50000,6.286571423,7953.461,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,secp256k1,256,MODEXP,50000,6.421734047,7786.059,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,secp256k1,256,EXPONENTIATION,781,0.002097429,372360.637,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,secp256k1,256,EXPONENTIATION,781,0.000320769,2434773.943,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,secp256k1,256,EXPONENTIATION,781,0.003544550,220338.266,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,secp256k1,256,EXPONENTIATION,50000,0.042463230,1177489.324,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,secp256k1,256,EXPONENTIATION,50000,0.043202566,1157338.664,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,secp256k1,256,EXPONENTIATION,50000,0.013923177,3591134.409,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,secp256k1,256,EXPONENTIATION,50000,0.014517943,3444014.073,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,secp256k1,256,EXPONENTIATION,50000,0.000614779,81330038.806,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,secp256k1,256,EXPONENTIATION,50000,0.001243319,40214940.743,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,secp256k1,256,EXPONENTIATION,50000,0.000529204,94481523.107,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,secp256k1,256,EXPONENTIATION,50000,0.001094239,45693856.646,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,secp256k1,256,EXPONENTIATION,50000,0.000635531,78674368.398,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,secp256k1,256,EXPONENTIATION,50000,0.001121279,44591934.691,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,secp256k1,256,EXPONENTIATION,50000,0.000568789,87906059.963,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,secp256k1,256,EXPONENTIATION,50000,0.001255344,39829720.006,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,secp256k1,256,EXPONENTIATION,50000,2.569640599,19457.974,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,secp256k1,256,EXPONENTIATION,50000,2.453309329,20380.634,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,secp256k1,256,DIVIDE,6250,0.000176288,35453349.401,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,secp256k1,256,DIVIDE,6250,0.000056842,109953907.227,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,secp256k1,256,DIVIDE,6250,0.000110816,56399798.115,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,secp256k1,256,DIVIDE,50000,0.000016192,3087944664.032,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,secp256k1,256,DIVIDE,50000,0.000358252,139566563.387,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,secp256k1,256,DIVIDE,50000,0.001156998,43215286.400,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,secp256k1,256,DIVIDE,50000,0.000241346,207171447.831,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,secp256k1,256,DIVIDE,50000,0.001043917,47896528.169,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,secp256k1,256,DIVIDE,50000,0.000092037,543259774.272,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,secp256k1,256,DIVIDE,50000,0.000774774,64534948.097,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,secp256k1,256,DIVIDE,50000,0.000090370,553280966.104,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,secp256k1,256,DIVIDE,50000,0.000771306,64825114.773,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,secp256k1,256,DIVIDE,50000,0.000090260,553955237.754,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,secp256k1,256,DIVIDE,50000,0.000735074,68020362.431,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,secp256k1,256,DIVIDE,50000,0.000083026,602220993.705,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,secp256k1,256,DIVIDE,50000,0.000851496,58720181.955,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,secp256k1,256,DIVIDE,50000,0.024853870,2011759.135,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,secp256k1,256,DIVIDE,50000,0.024938178,2004958.021,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,secp256k1,256,ISQRT,1562,0.000171769,9093608.314,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,secp256k1,256,ISQRT,1562,0.000060440,25843811.572,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,secp256k1,256,ISQRT,50000,0.004288696,11658555.421,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,secp256k1,256,ISQRT,50000,0.005002444,9995114.390,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,secp256k1,256,ISQRT,50000,0.003104790,16104148.745,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,secp256k1,256,ISQRT,50000,0.003785045,13209882.581,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,secp256k1,256,ISQRT,50000,0.000648229,77133235.210,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,secp256k1,256,ISQRT,50000,0.001163666,42967655.595,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,secp256k1,256,ISQRT,50000,0.000731989,68307037.478,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,secp256k1,256,ISQRT,50000,0.001267233,39456043.280,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,secp256k1,256,ISQRT,50000,0.000600608,83248974.339,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,secp256k1,256,ISQRT,50000,0.001108741,45096194.608,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,secp256k1,256,ISQRT,50000,0.000933987,53533935.837,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,secp256k1,256,ISQRT,50000,0.001627347,30724854.573,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,secp256k1,256,ISQRT,50000,0.448364097,111516.512,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,secp256k1,256,ISQRT,50000,0.444221222,112556.532,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,secp256k1,256,MODMUL_R2,50000,0.003478303,14374825.875,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,secp256k1,256,MODMUL_R2,50000,0.000660564,75692892.624,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,secp256k1,256,MODMUL_R2,50000,0.001914226,26120217.745,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,secp256k1,256,MODMUL_R2,50000,0.000015104,3310381355.932,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,secp256k1,256,MODMUL_R2,50000,0.000121869,410276604.097,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,secp256k1,256,MODMUL_R2,50000,0.000712112,70213674.278,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,secp256k1,256,MODMUL_R2,50000,0.000029512,1694226242.463,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,secp256k1,256,MODMUL_R2,50000,0.000552155,90554282.591,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,secp256k1,256,MODMUL_R2,50000,0.000019390,2578649085.890,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,secp256k1,256,MODMUL_R2,50000,0.000499595,100081065.605,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,secp256k1,256,MODMUL_R2,50000,0.000015859,3152783710.188,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,secp256k1,256,MODMUL_R2,50000,0.000484971,103098948.500,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,secp256k1,256,MODMUL_R2,50000,0.000018204,2746648768.655,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,secp256k1,256,MODMUL_R2,50000,0.000494801,101050725.617,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,secp256k1,256,MODMUL_R2,50000,0.000015558,3213781049.491,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,secp256k1,256,MODMUL_R2,50000,0.000546158,91548599.429,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,secp256k1,256,MODMUL_R2,50000,0.025199219,1984188.478,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,secp256k1,256,MODMUL_R2,50000,0.023629855,2115967.280,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,rsa256(composite),256,ADD,50000,0.000748192,66827766.119,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,rsa256(composite),256,ADD,50000,0.000316317,158069278.405,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,rsa256(composite),256,ADD,50000,0.000252929,197683932.531,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,rsa256(composite),256,ADD,50000,0.000007840,6377551020.408,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,rsa256(composite),256,ADD,50000,0.000040973,1220315869.874,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,rsa256(composite),256,ADD,50000,0.000590795,84631725.066,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,rsa256(composite),256,ADD,50000,0.000024440,2045826527.832,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,rsa256(composite),256,ADD,50000,0.000546809,91439607.137,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32,rsa256(composite),256,ADD,50000,0.000014028,3564299304.666,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32,rsa256(composite),256,ADD,50000,0.000490687,101897951.636,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,rsa256(composite),256,ADD,50000,0.000013832,3614806444.946,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,rsa256(composite),256,ADD,50000,0.000493567,101303369.817,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,rsa256(composite),256,ADD,50000,0.000013036,3835531724.225,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,rsa256(composite),256,ADD,50000,0.000485574,102970916.729,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,rsa256(composite),256,ADD,50000,0.000008102,6171315176.010,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,rsa256(composite),256,ADD,50000,0.000488383,102378665.792,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,rsa256(composite),256,ADD,50000,0.000009311,5369992883.924,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,rsa256(composite),256,ADD,50000,0.000530612,94230812.696,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,rsa256(composite),256,ADD,50000,0.000157915,316626030.622,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,rsa256(composite),256,ADD,50000,0.000440875,113410830.944,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACT,50000,0.000495268,100955442.466,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACT,50000,0.000211041,236920789.371,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACT,50000,0.000180902,276392743.710,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,rsa256(composite),256,SUBTRACT,50000,0.000007104,7038288288.288,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,rsa256(composite),256,SUBTRACT,50000,0.000040289,1241033531.927,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,rsa256(composite),256,SUBTRACT,50000,0.000593451,84252954.719,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,rsa256(composite),256,SUBTRACT,50000,0.000025706,1945071320.608,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,rsa256(composite),256,SUBTRACT,50000,0.000552970,90420818.015,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32,rsa256(composite),256,SUBTRACT,50000,0.000013676,3656040026.069,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32,rsa256(composite),256,SUBTRACT,50000,0.000493091,101401161.256,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,rsa256(composite),256,SUBTRACT,50000,0.000013466,3713054996.318,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,rsa256(composite),256,SUBTRACT,50000,0.000496036,100799135.878,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,rsa256(composite),256,SUBTRACT,50000,0.000012393,4034535507.412,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,rsa256(composite),256,SUBTRACT,50000,0.000485110,103069407.277,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,rsa256(composite),256,SUBTRACT,50000,0.000007843,6375109861.643,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,rsa256(composite),256,SUBTRACT,50000,0.000492421,101539129.996,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,rsa256(composite),256,SUBTRACT,50000,0.000009309,5371145018.190,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,rsa256(composite),256,SUBTRACT,50000,0.000545669,91630640.311,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,rsa256(composite),256,SUBTRACT,50000,0.000173061,288915465.794,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,rsa256(composite),256,SUBTRACT,50000,0.000459872,108725907.046,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,rsa256(composite),256,ADDMOD,50000,0.001385313,36092926.282,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,rsa256(composite),256,ADDMOD,50000,0.000370151,135080008.066,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,rsa256(composite),256,ADDMOD,50000,0.000743719,67229692.960,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,rsa256(composite),256,ADDMOD,50000,0.000008160,6127450980.392,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,rsa256(composite),256,ADDMOD,50000,0.000054880,911078749.397,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,rsa256(composite),256,ADDMOD,50000,0.000654282,76419647.674,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,rsa256(composite),256,ADDMOD,50000,0.000030588,1634627845.450,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,rsa256(composite),256,ADDMOD,50000,0.000531617,94052673.424,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32,rsa256(composite),256,ADDMOD,50000,0.000014488,3451132071.642,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32,rsa256(composite),256,ADDMOD,50000,0.000490995,101834031.591,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,rsa256(composite),256,ADDMOD,50000,0.000011075,4514672418.407,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,rsa256(composite),256,ADDMOD,50000,0.000494166,101180574.772,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,rsa256(composite),256,ADDMOD,50000,0.000010223,4890931937.879,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,rsa256(composite),256,ADDMOD,50000,0.000485700,102944204.552,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,rsa256(composite),256,ADDMOD,50000,0.000007669,6519753831.499,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,rsa256(composite),256,ADDMOD,50000,0.000488698,102312676.002,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,rsa256(composite),256,ADDMOD,50000,0.000008143,6140242594.683,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,rsa256(composite),256,ADDMOD,50000,0.000542571,92153837.837,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,rsa256(composite),256,ADDMOD,50000,0.000235980,211882361.216,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,rsa256(composite),256,ADDMOD,50000,0.000506681,98681418.717,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACTMOD,50000,0.001411710,35418039.117,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACTMOD,50000,0.000399658,125106965.925,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACTMOD,50000,0.001000844,49957835.568,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,rsa256(composite),256,SUBTRACTMOD,50000,0.000007936,6300403225.806,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,rsa256(composite),256,SUBTRACTMOD,50000,0.000060477,826760595.534,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,rsa256(composite),256,SUBTRACTMOD,50000,0.000645691,77436420.658,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,rsa256(composite),256,SUBTRACTMOD,50000,0.000033410,1496558058.728,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,rsa256(composite),256,SUBTRACTMOD,50000,0.000532437,93907823.926,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32,rsa256(composite),256,SUBTRACTMOD,50000,0.000015259,3276754688.127,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32,rsa256(composite),256,SUBTRACTMOD,50000,0.000488285,102399214.178,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,rsa256(composite),256,SUBTRACTMOD,50000,0.000010444,4787437530.048,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,rsa256(composite),256,SUBTRACTMOD,50000,0.000488861,102278562.155,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,rsa256(composite),256,SUBTRACTMOD,50000,0.000010384,4815100467.360,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,rsa256(composite),256,SUBTRACTMOD,50000,0.000481545,103832455.713,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,rsa256(composite),256,SUBTRACTMOD,50000,0.000007488,6677349672.106,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,rsa256(composite),256,SUBTRACTMOD,50000,0.000489626,102118760.083,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,rsa256(composite),256,SUBTRACTMOD,50000,0.000008011,6241417802.254,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,rsa256(composite),256,SUBTRACTMOD,50000,0.000533289,93757793.691,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,rsa256(composite),256,SUBTRACTMOD,50000,0.000289565,172672802.576,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,rsa256(composite),256,SUBTRACTMOD,50000,0.000573050,87252421.019,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000692054,72248697.361,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000228842,218491360.785,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000335182,149172687.632,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001374322,36381575.771,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.002235292,22368442.257,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000419158,119286760.310,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001184072,42227161.866,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000084573,591205236.680,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000740606,67512280.534,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000020289,2464389506.182,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000663590,75347729.884,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000019764,2529852234.668,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000698942,71536694.045,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000009854,5074081394.586,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000654825,76356278.571,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000009901,5049993385.134,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000714113,70016930.009,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.004793968,10429773.422,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.005298751,9436186.001,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000691632,72292780.179,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000181078,276124099.183,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000390978,127884433.754,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000006976,7167431192.661,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000201591,248026946.268,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000946676,52816380.569,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000065788,760017013.540,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000734341,68088258.746,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000022476,2224595064.806,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000693020,72147989.984,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000022739,2198865640.863,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000666634,75003675.024,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000021817,2291790918.376,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000654315,76415793.481,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000016397,3049338320.761,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000705474,70874333.961,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000019796,2525762801.648,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000712456,70179772.530,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.004564855,10953250.431,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.004418746,11315427.497,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.006077145,8227547.639,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000870360,57447493.012,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000331448,150853225.019,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000009184,5444250871.080,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000102482,487890550.492,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000643530,77696455.249,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000029768,1679656116.047,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000499323,100135583.674,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000013118,3811555987.841,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000482949,103530600.759,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000013541,3692490208.459,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000484619,103173833.763,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000010884,4593899515.438,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000482317,103666260.551,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000012697,3937937852.426,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000530608,94231523.372,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000011759,4252061339.546,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000538343,92877589.231,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.023393054,2137386.593,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.021525281,2322850.048,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,rsa256(composite),256,COMPARE,50000,0.000217298,230098760.017,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,rsa256(composite),256,COMPARE,50000,0.000177305,281999939.496,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,rsa256(composite),256,COMPARE,50000,0.000185947,268893825.439,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,rsa256(composite),256,COMPARE,50000,0.000007904,6325910931.174,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,rsa256(composite),256,COMPARE,50000,0.000039424,1268263057.626,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,rsa256(composite),256,COMPARE,50000,0.000584525,85539540.417,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,rsa256(composite),256,COMPARE,50000,0.000022389,2233239340.799,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,rsa256(composite),256,COMPARE,50000,0.000516603,96786120.342,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,rsa256(composite),256,COMPARE,50000,0.000010164,4919323182.451,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,rsa256(composite),256,COMPARE,50000,0.000495302,100948512.269,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,rsa256(composite),256,COMPARE,50000,0.000009958,5021088934.011,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,rsa256(composite),256,COMPARE,50000,0.000489614,102121263.275,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,rsa256(composite),256,COMPARE,50000,0.000007874,6350011641.715,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,rsa256(composite),256,COMPARE,50000,0.000492793,101462480.487,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,rsa256(composite),256,COMPARE,50000,0.000008080,6188120171.832,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,rsa256(composite),256,COMPARE,50000,0.000525098,95220320.238,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,rsa256(composite),256,COMPARE,50000,0.000103240,484308409.521,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,rsa256(composite),256,COMPARE,50000,0.000379378,131794673.911,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,rsa256(composite),256,REDUCE,6250,0.000128585,48605980.018,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,rsa256(composite),256,REDUCE,6250,0.000057014,109622197.962,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,rsa256(composite),256,REDUCE,6250,0.000103212,60554975.362,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,rsa256(composite),256,REDUCE,50000,0.000012224,4090314136.126,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,rsa256(composite),256,REDUCE,50000,0.000170159,293842816.141,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,rsa256(composite),256,REDUCE,50000,0.000744721,67139237.293,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,rsa256(composite),256,REDUCE,50000,0.000107577,464783358.113,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,rsa256(composite),256,REDUCE,50000,0.000613943,81440785.134,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,rsa256(composite),256,REDUCE,50000,0.000035059,1426167397.797,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,rsa256(composite),256,REDUCE,50000,0.000522556,95683524.522,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,rsa256(composite),256,REDUCE,50000,0.000033471,1493830548.850,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,rsa256(composite),256,REDUCE,50000,0.000511826,97689449.565,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,rsa256(composite),256,REDUCE,50000,0.000035041,1426899878.312,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,rsa256(composite),256,REDUCE,50000,0.000542696,92132611.741,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,rsa256(composite),256,REDUCE,50000,0.000042371,1180052474.888,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,rsa256(composite),256,REDUCE,50000,0.000591940,84468020.530,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,rsa256(composite),256,REDUCE,50000,0.016886609,2960926.022,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,rsa256(composite),256,REDUCE,50000,0.016278431,3071549.094,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,rsa256(composite),256,MODMUL,3125,0.000197472,15825028.079,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,rsa256(composite),256,MODMUL,3125,0.000054630,57203002.228,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,rsa256(composite),256,MODMUL,3125,0.000130282,23986429.084,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,rsa256(composite),256,MODMUL,50000,0.000037088,1348144952.545,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,rsa256(composite),256,MODMUL,50000,0.000464647,107608572.366,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,rsa256(composite),256,MODMUL,50000,0.001060501,47147527.473,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,rsa256(composite),256,MODMUL,50000,0.000257219,194386885.280,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,rsa256(composite),256,MODMUL,50000,0.000780430,64067244.955,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,rsa256(composite),256,MODMUL,50000,0.000103076,485078973.188,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,rsa256(composite),256,MODMUL,50000,0.000614515,81364978.879,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,rsa256(composite),256,MODMUL,50000,0.000104790,477144763.546,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,rsa256(composite),256,MODMUL,50000,0.000614175,81410021.813,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,rsa256(composite),256,MODMUL,50000,0.000103197,484510197.748,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,rsa256(composite),256,MODMUL,50000,0.000596577,83811477.861,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,rsa256(composite),256,MODMUL,50000,0.000125768,397557414.196,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,rsa256(composite),256,MODMUL,50000,0.000710419,70381000.426,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,rsa256(composite),256,MODMUL,50000,0.063472158,787746.968,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,rsa256(composite),256,MODMUL,50000,0.077629215,644087.410,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,rsa256(composite),256,MODEXP,781,0.005827087,134029.233,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,rsa256(composite),256,MODEXP,781,0.001228928,635513.227,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,rsa256(composite),256,MODEXP,781,0.001500595,520460.217,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,rsa256(composite),256,MODEXP,50000,0.011109120,4500806.545,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,rsa256(composite),256,MODEXP,50000,0.021332830,2343805.299,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,rsa256(composite),256,MODEXP,50000,0.022061489,2266392.808,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,rsa256(composite),256,MODEXP,50000,0.003854968,12970276.267,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,rsa256(composite),256,MODEXP,50000,0.004469172,11187754.690,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,rsa256(composite),256,MODEXP,50000,0.002987909,16734110.712,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,rsa256(composite),256,MODEXP,50000,0.003606608,13863441.774,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,rsa256(composite),256,MODEXP,50000,0.001854993,26954279.571,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,rsa256(composite),256,MODEXP,50000,0.002471518,20230481.828,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,rsa256(composite),256,MODEXP,50000,0.002940720,17002638.820,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,rsa256(composite),256,MODEXP,50000,0.003733455,13392420.694,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,rsa256(composite),256,MODEXP,50000,0.001902259,26284538.560,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,rsa256(composite),256,MODEXP,50000,0.002610914,19150381.808,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,rsa256(composite),256,MODEXP,50000,6.293343371,7944.903,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,rsa256(composite),256,MODEXP,50000,6.361847279,7859.352,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,rsa256(composite),256,EXPONENTIATION,781,0.002117819,368775.613,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,rsa256(composite),256,EXPONENTIATION,781,0.000310509,2515225.002,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,rsa256(composite),256,EXPONENTIATION,781,0.003489077,223841.434,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,rsa256(composite),256,EXPONENTIATION,50000,0.042419020,1178716.528,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,rsa256(composite),256,EXPONENTIATION,50000,0.043231658,1156559.852,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,rsa256(composite),256,EXPONENTIATION,50000,0.013834148,3614244.982,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,rsa256(composite),256,EXPONENTIATION,50000,0.014374200,3478454.453,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,rsa256(composite),256,EXPONENTIATION,50000,0.000619680,80686806.271,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,rsa256(composite),256,EXPONENTIATION,50000,0.001124559,44461873.521,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,rsa256(composite),256,EXPONENTIATION,50000,0.000527702,94750445.785,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,rsa256(composite),256,EXPONENTIATION,50000,0.001066780,46870019.992,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,rsa256(composite),256,EXPONENTIATION,50000,0.000643277,77727013.114,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,rsa256(composite),256,EXPONENTIATION,50000,0.001407257,35530112.808,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,rsa256(composite),256,EXPONENTIATION,50000,0.000579359,86302275.603,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,rsa256(composite),256,EXPONENTIATION,50000,0.001277409,39141731.452,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,rsa256(composite),256,EXPONENTIATION,50000,2.520409667,19838.045,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,rsa256(composite),256,EXPONENTIATION,50000,2.461429260,20313.401,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,rsa256(composite),256,DIVIDE,6250,0.000299950,20836806.182,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,rsa256(composite),256,DIVIDE,6250,0.000074642,83733021.923,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,rsa256(composite),256,DIVIDE,6250,0.000137065,45598802.659,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,rsa256(composite),256,DIVIDE,50000,0.000015808,3162955465.587,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,rsa256(composite),256,DIVIDE,50000,0.000356812,140129815.880,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,rsa256(composite),256,DIVIDE,50000,0.001154827,43296528.349,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,rsa256(composite),256,DIVIDE,50000,0.000250268,199785829.634,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,rsa256(composite),256,DIVIDE,50000,0.000924315,54094113.056,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,rsa256(composite),256,DIVIDE,50000,0.000101042,494843718.969,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,rsa256(composite),256,DIVIDE,50000,0.000747722,66869772.565,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,rsa256(composite),256,DIVIDE,50000,0.000097259,514091228.442,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,rsa256(composite),256,DIVIDE,50000,0.000738147,67737185.266,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,rsa256(composite),256,DIVIDE,50000,0.000097166,514583284.789,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,rsa256(composite),256,DIVIDE,50000,0.000892919,55996120.541,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,rsa256(composite),256,DIVIDE,50000,0.000089896,556198280.823,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,rsa256(composite),256,DIVIDE,50000,0.000867211,57656095.365,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,rsa256(composite),256,DIVIDE,50000,0.024124622,2072571.334,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,rsa256(composite),256,DIVIDE,50000,0.024303254,2057337.672,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,rsa256(composite),256,ISQRT,1562,0.000110612,14121433.499,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,rsa256(composite),256,ISQRT,1562,0.000049534,31533896.472,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,rsa256(composite),256,ISQRT,50000,0.004294079,11643940.420,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,rsa256(composite),256,ISQRT,50000,0.005003741,9992523.594,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,rsa256(composite),256,ISQRT,50000,0.003105378,16101099.446,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,rsa256(composite),256,ISQRT,50000,0.003717415,13450206.668,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,rsa256(composite),256,ISQRT,50000,0.000647800,77184315.765,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,rsa256(composite),256,ISQRT,50000,0.001162779,43000432.552,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,rsa256(composite),256,ISQRT,50000,0.000736112,67924446.225,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,rsa256(composite),256,ISQRT,50000,0.001324169,37759530.711,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,rsa256(composite),256,ISQRT,50000,0.000602750,82953131.686,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,rsa256(composite),256,ISQRT,50000,0.001122655,44537279.875,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,rsa256(composite),256,ISQRT,50000,0.000933383,53568578.045,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,rsa256(composite),256,ISQRT,50000,0.001540161,32464138.492,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,rsa256(composite),256,ISQRT,50000,0.453063249,110359.867,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,rsa256(composite),256,ISQRT,50000,0.458620248,109022.661,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,rsa256(composite),256,MODMUL_R2,50000,0.003469103,14412947.664,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,rsa256(composite),256,MODMUL_R2,50000,0.000739652,67599357.412,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,rsa256(composite),256,MODMUL_R2,50000,0.001807847,27657207.741,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,rsa256(composite),256,MODMUL_R2,50000,0.000014752,3389370932.755,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,rsa256(composite),256,MODMUL_R2,50000,0.000119270,419216902.330,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,rsa256(composite),256,MODMUL_R2,50000,0.000714378,69990957.313,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,rsa256(composite),256,MODMUL_R2,50000,0.000031065,1609528253.467,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,rsa256(composite),256,MODMUL_R2,50000,0.000514506,97180596.114,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,rsa256(composite),256,MODMUL_R2,50000,0.000019394,2578117004.935,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,rsa256(composite),256,MODMUL_R2,50000,0.000490985,101836105.079,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,rsa256(composite),256,MODMUL_R2,50000,0.000015403,3246121192.722,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,rsa256(composite),256,MODMUL_R2,50000,0.000487571,102549167.173,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,rsa256(composite),256,MODMUL_R2,50000,0.000018531,2698181687.170,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,rsa256(composite),256,MODMUL_R2,50000,0.000495264,100956257.454,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,rsa256(composite),256,MODMUL_R2,50000,0.000015639,3197135803.053,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,rsa256(composite),256,MODMUL_R2,50000,0.000600423,83274624.481,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,rsa256(composite),256,MODMUL_R2,50000,0.023662197,2113075.130,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,rsa256(composite),256,MODMUL_R2,50000,0.024934292,2005270.493,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,brainpoolP512r1,512,ADD,25000,0.000403813,61909844.540,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,brainpoolP512r1,512,ADD,25000,0.000142122,175905207.709,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,brainpoolP512r1,512,ADD,25000,0.000133487,187284155.850,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,brainpoolP512r1,512,ADD,50000,0.000008768,5702554744.526,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,brainpoolP512r1,512,ADD,50000,0.000092235,542093581.321,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,brainpoolP512r1,512,ADD,50000,0.001043087,47934640.213,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,brainpoolP512r1,512,ADD,50000,0.000047211,1059075204.483,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,brainpoolP512r1,512,ADD,50000,0.000985913,50714413.925,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32,brainpoolP512r1,512,ADD,50000,0.000027936,1789805317.757,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32,brainpoolP512r1,512,ADD,50000,0.000949236,52673940.061,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,brainpoolP512r1,512,ADD,50000,0.000025262,1979257577.940,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,brainpoolP512r1,512,ADD,50000,0.000969977,51547614.121,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,brainpoolP512r1,512,ADD,50000,0.000025002,1999839846.763,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,brainpoolP512r1,512,ADD,50000,0.000961138,52021665.929,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,brainpoolP512r1,512,ADD,50000,0.000014299,3496748194.609,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,brainpoolP512r1,512,ADD,50000,0.000930331,53744312.518,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,brainpoolP512r1,512,ADD,50000,0.000014549,3436662028.371,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,brainpoolP512r1,512,ADD,50000,0.000969067,51596019.802,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,brainpoolP512r1,512,ADD,50000,0.000309163,161726986.097,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,brainpoolP512r1,512,ADD,50000,0.000958505,52164568.622,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACT,25000,0.000251424,99433625.887,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACT,25000,0.000121492,205774869.966,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACT,25000,0.000111583,224048462.443,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,brainpoolP512r1,512,SUBTRACT,50000,0.000008256,6056201550.388,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,brainpoolP512r1,512,SUBTRACT,50000,0.000091519,546334644.491,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,brainpoolP512r1,512,SUBTRACT,50000,0.001043454,47917780.674,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,brainpoolP512r1,512,SUBTRACT,50000,0.000047317,1056702658.738,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,brainpoolP512r1,512,SUBTRACT,50000,0.000961640,51994509.465,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32,brainpoolP512r1,512,SUBTRACT,50000,0.000026525,1885014003.605,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32,brainpoolP512r1,512,SUBTRACT,50000,0.000957120,52240053.552,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,brainpoolP512r1,512,SUBTRACT,50000,0.000025341,1973087103.357,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,brainpoolP512r1,512,SUBTRACT,50000,0.000969029,51598043.045,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,brainpoolP512r1,512,SUBTRACT,50000,0.000025431,1966104375.302,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,brainpoolP512r1,512,SUBTRACT,50000,0.000961108,52023289.863,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,brainpoolP512r1,512,SUBTRACT,50000,0.000014444,3461644527.203,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,brainpoolP512r1,512,SUBTRACT,50000,0.000927345,53917366.163,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,brainpoolP512r1,512,SUBTRACT,50000,0.000014643,3414600995.519,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,brainpoolP512r1,512,SUBTRACT,50000,0.000938787,53260217.629,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,brainpoolP512r1,512,SUBTRACT,50000,0.000333780,149799268.620,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,brainpoolP512r1,512,SUBTRACT,50000,0.000949028,52685484.572,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,brainpoolP512r1,512,ADDMOD,25000,0.000779485,32072458.141,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,brainpoolP512r1,512,ADDMOD,25000,0.000239084,104565758.525,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,brainpoolP512r1,512,ADDMOD,25000,0.000479921,52091906.784,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,brainpoolP512r1,512,ADDMOD,50000,0.000008960,5580357142.857,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,brainpoolP512r1,512,ADDMOD,50000,0.000124547,401454868.660,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,brainpoolP512r1,512,ADDMOD,50000,0.001071344,46670350.566,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,brainpoolP512r1,512,ADDMOD,50000,0.000065171,767212413.574,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,brainpoolP512r1,512,ADDMOD,50000,0.000973046,51385032.126,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32,brainpoolP512r1,512,ADDMOD,50000,0.000029740,1681237380.954,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32,brainpoolP512r1,512,ADDMOD,50000,0.000957105,52240872.239,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,brainpoolP512r1,512,ADDMOD,50000,0.000025763,1940767909.305,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,brainpoolP512r1,512,ADDMOD,50000,0.000935054,53472847.641,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,brainpoolP512r1,512,ADDMOD,50000,0.000026551,1883168321.535,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,brainpoolP512r1,512,ADDMOD,50000,0.000960853,52037096.111,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,brainpoolP512r1,512,ADDMOD,50000,0.000013172,3795931590.738,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,brainpoolP512r1,512,ADDMOD,50000,0.000925285,54037404.718,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,brainpoolP512r1,512,ADDMOD,50000,0.000013387,3734966292.162,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,brainpoolP512r1,512,ADDMOD,50000,0.000952869,52473110.217,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,brainpoolP512r1,512,ADDMOD,50000,0.000417668,119712308.048,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,brainpoolP512r1,512,ADDMOD,50000,0.001092625,45761354.547,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000766845,32601112.280,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000170624,146521005.234,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000476045,52516043.736,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000009152,5463286713.287,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000146353,341639728.250,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,brainpoolP512r1,512,SUBTRACTMOD,50000,0.001078968,46340577.214,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000073966,675986272.771,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000983547,50836411.617,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000032957,1517128241.937,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000955211,52344455.906,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000025898,1930651024.753,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000968217,51641316.028,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000025602,1952972472.200,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000963226,51908897.768,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000012872,3884399609.861,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000930338,53743907.958,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000013627,3669186407.568,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000948354,52722928.472,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000514471,97187207.600,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,brainpoolP512r1,512,SUBTRACTMOD,50000,0.001169522,42752509.209,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.000838609,29811270.893,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.000280752,89046560.816,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.000247639,100953403.943,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.007158891,6984322.013,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.008680675,5759920.743,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.001998453,25019352.458,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.003396647,14720399.264,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000535822,93314570.791,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.001780738,28078246.168,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000073768,677800687.067,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.001364439,36645097.409,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000073408,681124652.569,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.001318508,37921650.765,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000046538,1074390823.207,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.001267517,39447202.678,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000045575,1097092670.090,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.001333859,37485221.459,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.020697536,2415746.493,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.022092309,2263231.064,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000837866,29837706.665,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000217583,114898681.625,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000219995,113638944.720,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000009024,5540780141.844,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.001389970,35971999.425,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.002779062,17991682.081,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000378404,132133910.424,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.001656494,30184232.454,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000106596,469060742.641,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.001377680,36292898.286,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000138833,360144926.608,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.001396726,35798001.848,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000138026,362250590.314,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.001364832,36634545.545,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000073988,675785281.713,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.001296349,38569860.408,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000069729,717061751.737,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.001367771,36555827.026,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.018722318,2670609.484,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.020186193,2476940.550,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.007371159,3391596.898,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.001223307,20436407.257,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000460747,54259713.027,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000013408,3729116945.107,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000343830,145420703.237,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.001260758,39658681.547,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000106070,471386829.556,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.001039556,48097456.895,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000034631,1443793102.285,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000967577,51675473.806,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000036486,1370388548.314,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000970932,51496912.235,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000027764,1800893292.224,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000956171,52291901.794,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000035604,1404336591.227,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000948438,52718258.925,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000027473,1819968633.269,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000971686,51456952.314,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.104400112,478926.689,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.104668107,477700.433,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,brainpoolP512r1,512,COMPARE,25000,0.000111416,224384288.631,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,brainpoolP512r1,512,COMPARE,25000,0.000057804,432496035.544,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,brainpoolP512r1,512,COMPARE,25000,0.000166651,150014101.555,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,brainpoolP512r1,512,COMPARE,50000,0.000008640,5787037037.037,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,brainpoolP512r1,512,COMPARE,50000,0.000096823,516406237.126,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,brainpoolP512r1,512,COMPARE,50000,0.001030870,48502720.991,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,brainpoolP512r1,512,COMPARE,50000,0.000048926,1021951478.737,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,brainpoolP512r1,512,COMPARE,50000,0.001003684,49816476.130,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,brainpoolP512r1,512,COMPARE,50000,0.000020512,2437597775.072,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,brainpoolP512r1,512,COMPARE,50000,0.000957512,52218666.797,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,brainpoolP512r1,512,COMPARE,50000,0.000018098,2762736413.803,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,brainpoolP512r1,512,COMPARE,50000,0.000931141,53697560.127,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,brainpoolP512r1,512,COMPARE,50000,0.000009483,5272591425.989,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,brainpoolP512r1,512,COMPARE,50000,0.000893765,55943116.925,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,brainpoolP512r1,512,COMPARE,50000,0.000009641,5186185768.747,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,brainpoolP512r1,512,COMPARE,50000,0.000944425,52942266.513,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,brainpoolP512r1,512,COMPARE,50000,0.000120173,416066833.339,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,brainpoolP512r1,512,COMPARE,50000,0.000687588,72717964.935,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,brainpoolP512r1,512,REDUCE,3125,0.000067025,46624395.687,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,brainpoolP512r1,512,REDUCE,3125,0.000043735,71453071.170,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,brainpoolP512r1,512,REDUCE,3125,0.000076516,40841131.126,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,brainpoolP512r1,512,REDUCE,50000,0.000020480,2441406250.000,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,brainpoolP512r1,512,REDUCE,50000,0.000519527,96241389.045,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,brainpoolP512r1,512,REDUCE,50000,0.001432078,34914299.380,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,brainpoolP512r1,512,REDUCE,50000,0.000395358,126467656.619,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,brainpoolP512r1,512,REDUCE,50000,0.001310844,38143364.128,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,brainpoolP512r1,512,REDUCE,50000,0.000102533,487647886.706,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,brainpoolP512r1,512,REDUCE,50000,0.001040057,48074288.280,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,brainpoolP512r1,512,REDUCE,50000,0.000099060,504744589.258,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,brainpoolP512r1,512,REDUCE,50000,0.001003271,49836983.332,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,brainpoolP512r1,512,REDUCE,50000,0.000119150,419639118.211,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,brainpoolP512r1,512,REDUCE,50000,0.001001419,49929150.480,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,brainpoolP512r1,512,REDUCE,50000,0.000112846,443081726.092,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,brainpoolP512r1,512,REDUCE,50000,0.001042311,47970327.485,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,brainpoolP512r1,512,REDUCE,50000,0.063724481,784627.810,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,brainpoolP512r1,512,REDUCE,50000,0.051911386,963179.831,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL,1562,0.000220445,7085667.650,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL,1562,0.000056469,27661193.795,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL,1562,0.000104961,14881717.967,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,brainpoolP512r1,512,MODMUL,50000,0.000112736,443514050.525,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,brainpoolP512r1,512,MODMUL,50000,0.001491108,33532111.697,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,brainpoolP512r1,512,MODMUL,50000,0.002476793,20187395.549,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,brainpoolP512r1,512,MODMUL,50000,0.001025007,48780154.631,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,brainpoolP512r1,512,MODMUL,50000,0.001973259,25338792.312,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,brainpoolP512r1,512,MODMUL,50000,0.000384351,130089423.977,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,brainpoolP512r1,512,MODMUL,50000,0.001329501,37608095.084,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,brainpoolP512r1,512,MODMUL,50000,0.000382732,130639716.664,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,brainpoolP512r1,512,MODMUL,50000,0.001315554,38006801.696,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,brainpoolP512r1,512,MODMUL,50000,0.000470794,106203562.756,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,brainpoolP512r1,512,MODMUL,50000,0.001356613,36856494.812,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,brainpoolP512r1,512,MODMUL,50000,0.000454171,110090693.484,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,brainpoolP512r1,512,MODMUL,50000,0.001428757,34995454.063,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,brainpoolP512r1,512,MODMUL,50000,0.294360706,169859.628,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,brainpoolP512r1,512,MODMUL,50000,0.312888170,159801.504,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,brainpoolP512r1,512,MODEXP,390,0.014935877,26111.624,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,brainpoolP512r1,512,MODEXP,390,0.002564852,152055.557,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,brainpoolP512r1,512,MODEXP,390,0.002475629,157535.721,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,brainpoolP512r1,512,MODEXP,50000,0.030459423,1641528.141,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,brainpoolP512r1,512,MODEXP,50000,0.322199942,155183.144,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,brainpoolP512r1,512,MODEXP,50000,0.323566233,154527.868,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,brainpoolP512r1,512,MODEXP,50000,0.027060979,1847678.903,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,brainpoolP512r1,512,MODEXP,50000,0.028251152,1769839.332,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,brainpoolP512r1,512,MODEXP,50000,0.023711413,2108689.178,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,brainpoolP512r1,512,MODEXP,50000,0.024878191,2009792.432,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,brainpoolP512r1,512,MODEXP,50000,0.013410425,3728442.611,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,brainpoolP512r1,512,MODEXP,50000,0.014599409,3424796.168,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,brainpoolP512r1,512,MODEXP,50000,0.024447361,2045210.524,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,brainpoolP512r1,512,MODEXP,50000,0.025493727,1961266.785,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,brainpoolP512r1,512,MODEXP,50000,0.013763017,3632924.380,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,brainpoolP512r1,512,MODEXP,50000,0.015029335,3326827.168,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,brainpoolP512r1,512,MODEXP,50000,52.198876514,957.875,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,brainpoolP512r1,512,MODEXP,50000,52.517536653,952.063,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,brainpoolP512r1,512,EXPONENTIATION,390,0.003531009,110450.016,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,brainpoolP512r1,512,EXPONENTIATION,390,0.000849569,459056.299,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,brainpoolP512r1,512,EXPONENTIATION,390,0.006138757,63530.777,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,brainpoolP512r1,512,EXPONENTIATION,50000,0.332045737,150581.665,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,brainpoolP512r1,512,EXPONENTIATION,50000,0.332099947,150557.085,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,brainpoolP512r1,512,EXPONENTIATION,50000,0.081532992,613248.683,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,brainpoolP512r1,512,EXPONENTIATION,50000,0.082187422,608365.596,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,brainpoolP512r1,512,EXPONENTIATION,50000,0.022386581,2233480.852,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,brainpoolP512r1,512,EXPONENTIATION,50000,0.023540978,2123955.938,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,brainpoolP512r1,512,EXPONENTIATION,50000,0.021582602,2316680.815,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,brainpoolP512r1,512,EXPONENTIATION,50000,0.022918386,2181654.502,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,brainpoolP512r1,512,EXPONENTIATION,50000,0.022190378,2253228.854,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,brainpoolP512r1,512,EXPONENTIATION,50000,0.023669028,2112465.286,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,brainpoolP512r1,512,EXPONENTIATION,50000,0.021760405,2297751.351,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,brainpoolP512r1,512,EXPONENTIATION,50000,0.022993836,2174495.808,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,brainpoolP512r1,512,EXPONENTIATION,50000,22.751972197,2197.612,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,brainpoolP512r1,512,EXPONENTIATION,50000,23.055085257,2168.719,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,brainpoolP512r1,512,DIVIDE,3125,0.000098225,31814710.427,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,brainpoolP512r1,512,DIVIDE,3125,0.000073153,42718684.489,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,brainpoolP512r1,512,DIVIDE,3125,0.000093108,33563174.158,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,brainpoolP512r1,512,DIVIDE,50000,0.000029024,1722712238.148,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,brainpoolP512r1,512,DIVIDE,50000,0.001038141,48163014.498,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,brainpoolP512r1,512,DIVIDE,50000,0.002309324,21651357.733,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,brainpoolP512r1,512,DIVIDE,50000,0.000972776,51399294.353,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,brainpoolP512r1,512,DIVIDE,50000,0.002235268,22368682.424,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,brainpoolP512r1,512,DIVIDE,50000,0.000307621,162537669.281,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,brainpoolP512r1,512,DIVIDE,50000,0.001580930,31626953.713,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,brainpoolP512r1,512,DIVIDE,50000,0.000298734,167372979.938,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,brainpoolP512r1,512,DIVIDE,50000,0.001544864,32365308.519,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,brainpoolP512r1,512,DIVIDE,50000,0.000292709,170818116.846,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,brainpoolP512r1,512,DIVIDE,50000,0.001547303,32314291.387,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,brainpoolP512r1,512,DIVIDE,50000,0.000286003,174823341.058,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,brainpoolP512r1,512,DIVIDE,50000,0.001572102,31804552.141,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,brainpoolP512r1,512,DIVIDE,50000,0.111050843,450244.218,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,brainpoolP512r1,512,DIVIDE,50000,0.106513311,469424.897,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,brainpoolP512r1,512,ISQRT,781,0.000106688,7320410.952,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,brainpoolP512r1,512,ISQRT,781,0.000056227,13890124.028,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,brainpoolP512r1,512,ISQRT,50000,0.022354831,2236653.008,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,brainpoolP512r1,512,ISQRT,50000,0.023447792,2132396.944,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,brainpoolP512r1,512,ISQRT,50000,0.020878526,2394805.074,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,brainpoolP512r1,512,ISQRT,50000,0.022114978,2260911.135,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,brainpoolP512r1,512,ISQRT,50000,0.004084553,12241241.580,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,brainpoolP512r1,512,ISQRT,50000,0.005285030,9460684.238,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,brainpoolP512r1,512,ISQRT,50000,0.003996384,12511310.227,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,brainpoolP512r1,512,ISQRT,50000,0.005182422,9647998.556,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,brainpoolP512r1,512,ISQRT,50000,0.005138899,9729710.587,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,brainpoolP512r1,512,ISQRT,50000,0.006234706,8019624.340,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,brainpoolP512r1,512,ISQRT,50000,0.005139817,9727972.805,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,brainpoolP512r1,512,ISQRT,50000,0.006227969,8028299.435,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,brainpoolP512r1,512,ISQRT,50000,1.998956578,25013.050,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,brainpoolP512r1,512,ISQRT,50000,1.979333186,25261.033,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL_R2,25000,0.003930724,6360151.463,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL_R2,25000,0.000899148,27804099.022,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL_R2,25000,0.001425284,17540363.873,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,brainpoolP512r1,512,MODMUL_R2,50000,0.000023584,2120081411.126,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,brainpoolP512r1,512,MODMUL_R2,50000,0.000568332,87976745.766,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,brainpoolP512r1,512,MODMUL_R2,50000,0.001512130,33065940.052,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,brainpoolP512r1,512,MODMUL_R2,50000,0.000110170,453844059.301,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,brainpoolP512r1,512,MODMUL_R2,50000,0.001015829,49220882.707,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,brainpoolP512r1,512,MODMUL_R2,50000,0.000063019,793411502.465,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,brainpoolP512r1,512,MODMUL_R2,50000,0.001008403,49583351.100,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,brainpoolP512r1,512,MODMUL_R2,50000,0.000040611,1231193593.875,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,brainpoolP512r1,512,MODMUL_R2,50000,0.000972377,51420385.403,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,brainpoolP512r1,512,MODMUL_R2,50000,0.000063740,784436764.932,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,brainpoolP512r1,512,MODMUL_R2,50000,0.000951832,52530278.437,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,brainpoolP512r1,512,MODMUL_R2,50000,0.000044873,1114255817.313,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,brainpoolP512r1,512,MODMUL_R2,50000,0.001003783,49811562.870,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,brainpoolP512r1,512,MODMUL_R2,50000,0.107803719,463805.892,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,brainpoolP512r1,512,MODMUL_R2,50000,0.116525464,429090.761,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,p1024,1024,ADD,12500,0.000239553,52180519.743,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,p1024,1024,ADD,12500,0.000099005,126256245.866,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,p1024,1024,ADD,12500,0.000101933,122629569.041,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,p1024,1024,ADD,50000,0.000013536,3693853427.896,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,p1024,1024,ADD,50000,0.000241546,206999908.642,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,p1024,1024,ADD,50000,0.002147010,23288200.800,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,p1024,1024,ADD,50000,0.000123903,403541483.632,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,p1024,1024,ADD,50000,0.002032101,24605076.206,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32,p1024,1024,ADD,50000,0.000047244,1058335436.611,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32,p1024,1024,ADD,50000,0.001928535,25926415.660,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,p1024,1024,ADD,50000,0.000048005,1041558164.934,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,p1024,1024,ADD,50000,0.001925634,25965474.202,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,p1024,1024,ADD,50000,0.000046632,1072225044.400,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,p1024,1024,ADD,50000,0.001860694,26871694.119,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,p1024,1024,ADD,50000,0.000022264,2245778201.766,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,p1024,1024,ADD,50000,0.001792492,27894127.292,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,p1024,1024,ADD,50000,0.000022334,2238739205.453,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,p1024,1024,ADD,50000,0.002503386,19972948.666,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,p1024,1024,ADD,50000,0.000622218,80357687.908,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,p1024,1024,ADD,50000,0.002244219,22279465.612,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,p1024,1024,SUBTRACT,12500,0.000159380,78428911.692,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,p1024,1024,SUBTRACT,12500,0.000136289,91716866.511,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,p1024,1024,SUBTRACT,12500,0.000157185,79524127.258,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,p1024,1024,SUBTRACT,50000,0.000014304,3495525727.069,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,p1024,1024,SUBTRACT,50000,0.000240775,207662755.824,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,p1024,1024,SUBTRACT,50000,0.002107141,23728834.461,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,p1024,1024,SUBTRACT,50000,0.000124664,401078102.507,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,p1024,1024,SUBTRACT,50000,0.001958274,25532688.496,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32,p1024,1024,SUBTRACT,50000,0.000047501,1052609455.651,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32,p1024,1024,SUBTRACT,50000,0.001907844,26207593.518,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,p1024,1024,SUBTRACT,50000,0.000046966,1064599905.003,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,p1024,1024,SUBTRACT,50000,0.001876012,26652281.540,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,p1024,1024,SUBTRACT,50000,0.000046825,1067805658.760,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,p1024,1024,SUBTRACT,50000,0.001864854,26811750.425,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,p1024,1024,SUBTRACT,50000,0.000021909,2282167169.885,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,p1024,1024,SUBTRACT,50000,0.001813402,27572485.318,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,p1024,1024,SUBTRACT,50000,0.000023296,2146291450.171,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,p1024,1024,SUBTRACT,50000,0.002074625,24100741.085,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,p1024,1024,SUBTRACT,50000,0.000674558,74122610.719,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,p1024,1024,SUBTRACT,50000,0.002526158,19792902.888,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,p1024,1024,ADDMOD,12500,0.000538413,23216378.428,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,p1024,1024,ADDMOD,12500,0.000157366,79432659.625,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,p1024,1024,ADDMOD,12500,0.000300104,41652227.423,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,p1024,1024,ADDMOD,50000,0.000014304,3495525727.069,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,p1024,1024,ADDMOD,50000,0.000372417,134258103.765,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,p1024,1024,ADDMOD,50000,0.002251524,22207180.543,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,p1024,1024,ADDMOD,50000,0.000187019,267352511.234,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,p1024,1024,ADDMOD,50000,0.002021774,24730756.231,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32,p1024,1024,ADDMOD,50000,0.000057961,862649042.348,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32,p1024,1024,ADDMOD,50000,0.001945615,25698815.024,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,p1024,1024,ADDMOD,50000,0.000062089,805295641.210,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,p1024,1024,ADDMOD,50000,0.001929658,25911327.300,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,p1024,1024,ADDMOD,50000,0.000061271,816046754.610,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,p1024,1024,ADDMOD,50000,0.001870831,26726091.267,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,p1024,1024,ADDMOD,50000,0.000022082,2264287633.146,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,p1024,1024,ADDMOD,50000,0.001798675,27798240.387,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,p1024,1024,ADDMOD,50000,0.000022173,2254994750.863,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,p1024,1024,ADDMOD,50000,0.002188637,22845268.563,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,p1024,1024,ADDMOD,50000,0.001069379,46756108.024,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,p1024,1024,ADDMOD,50000,0.002847685,17558121.782,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,p1024,1024,SUBTRACTMOD,12500,0.000458113,27285844.232,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,p1024,1024,SUBTRACTMOD,12500,0.000192138,65057405.774,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,p1024,1024,SUBTRACTMOD,12500,0.000373000,33512064.438,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,p1024,1024,SUBTRACTMOD,50000,0.000015200,3289473684.211,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,p1024,1024,SUBTRACTMOD,50000,0.000365809,136683351.503,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,p1024,1024,SUBTRACTMOD,50000,0.002238043,22340946.961,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,p1024,1024,SUBTRACTMOD,50000,0.000184485,271024742.756,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,p1024,1024,SUBTRACTMOD,50000,0.002032378,24601722.700,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32,p1024,1024,SUBTRACTMOD,50000,0.000057748,865830816.932,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32,p1024,1024,SUBTRACTMOD,50000,0.001894239,26395824.421,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,p1024,1024,SUBTRACTMOD,50000,0.000061983,806672811.873,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,p1024,1024,SUBTRACTMOD,50000,0.001909024,26191394.158,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,p1024,1024,SUBTRACTMOD,50000,0.000061914,807571827.797,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,p1024,1024,SUBTRACTMOD,50000,0.001837396,27212424.524,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,p1024,1024,SUBTRACTMOD,50000,0.000020935,2388344869.477,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,p1024,1024,SUBTRACTMOD,50000,0.001821116,27455692.025,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,p1024,1024,SUBTRACTMOD,50000,0.000021972,2275623660.772,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,p1024,1024,SUBTRACTMOD,50000,0.002234444,22376931.349,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,p1024,1024,SUBTRACTMOD,50000,0.001158741,43150281.315,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,p1024,1024,SUBTRACTMOD,50000,0.002837001,17624244.777,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.001524063,8201760.698,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.000470155,26586976.705,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.000462783,27010499.396,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.030849155,1620789.937,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.033683249,1484417.373,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.008276602,6041126.541,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.011037824,4529878.353,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.001536006,32551956.203,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.004255361,11749884.440,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.000200291,249636777.851,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.002641553,18928259.240,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.000200722,249100747.379,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.002653222,18845011.839,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.000143011,349623449.757,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.002581693,19367136.195,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.000143761,347799470.801,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.002972368,16821604.872,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.095031003,526144.084,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.108392749,461285.468,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.001517624,8236559.246,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.000372807,33529413.295,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.000392027,31885558.970,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000023232,2152203856.749,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.010441447,4788608.322,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.013139089,3805438.871,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.002670661,18721956.845,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.005895420,8481159.952,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000687140,72765375.315,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.003277968,15253352.074,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000688465,72625333.247,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.003293793,15180067.474,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000688519,72619637.456,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.003157167,15836982.963,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000305398,163720784.075,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.002729353,18319359.935,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000285964,174847183.075,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.002946882,16967085.872,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.063628401,785812.612,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.075920488,658583.754,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.010951986,1141345.506,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.001763661,7087529.863,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.001164284,10736212.127,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000035552,1406390639.064,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.003244830,15409127.765,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.005124395,9757249.392,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000358769,139365440.985,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.002209050,22634164.020,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000123144,406028714.023,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.001973801,25331834.392,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000146254,341870994.983,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.002014967,24814302.193,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000110092,454165599.912,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.001898849,26331741.002,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000109776,455472964.826,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.001882882,26555036.333,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000075740,660153145.130,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.002029998,24630566.154,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.477374613,104739.546,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.496452727,100714.524,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,p1024,1024,COMPARE,12500,0.000070733,176720908.130,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,p1024,1024,COMPARE,12500,0.000149074,83850973.546,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,p1024,1024,COMPARE,12500,0.000250574,49885462.727,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,p1024,1024,COMPARE,50000,0.000014944,3345824411.135,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,p1024,1024,COMPARE,50000,0.000233224,214386169.416,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,p1024,1024,COMPARE,50000,0.002038516,24527646.554,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,p1024,1024,COMPARE,50000,0.000119905,416996781.940,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,p1024,1024,COMPARE,50000,0.001958126,25534618.276,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,p1024,1024,COMPARE,50000,0.000036914,1354499655.430,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,p1024,1024,COMPARE,50000,0.001898336,26338856.735,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,p1024,1024,COMPARE,50000,0.000035808,1396335885.465,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,p1024,1024,COMPARE,50000,0.001847557,27062764.511,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,p1024,1024,COMPARE,50000,0.000013251,3773300278.415,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,p1024,1024,COMPARE,50000,0.001783126,28040643.213,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,p1024,1024,COMPARE,50000,0.000014300,3496503558.147,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,p1024,1024,COMPARE,50000,0.001817182,27515130.554,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,p1024,1024,COMPARE,50000,0.000338827,147567932.226,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,p1024,1024,COMPARE,50000,0.002147584,23281976.394,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,p1024,1024,REDUCE,1562,0.000036994,42223064.401,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,p1024,1024,REDUCE,1562,0.000034550,45209839.310,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,p1024,1024,REDUCE,1562,0.000067916,22998998.674,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,p1024,1024,REDUCE,50000,0.000030752,1625910509.886,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,p1024,1024,REDUCE,50000,0.002808546,17802806.150,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,p1024,1024,REDUCE,50000,0.004606917,10853245.235,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,p1024,1024,REDUCE,50000,0.001519689,32901468.618,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,p1024,1024,REDUCE,50000,0.003384046,14775212.858,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,p1024,1024,REDUCE,50000,0.000386747,129283484.919,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,p1024,1024,REDUCE,50000,0.002219014,22532530.206,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,p1024,1024,REDUCE,50000,0.000367475,136063677.605,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,p1024,1024,REDUCE,50000,0.002181427,22920776.167,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,p1024,1024,REDUCE,50000,0.000357464,139874225.279,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,p1024,1024,REDUCE,50000,0.002120986,23573941.567,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,p1024,1024,REDUCE,50000,0.000327269,152779519.003,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,p1024,1024,REDUCE,50000,0.002154402,23208296.335,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,p1024,1024,REDUCE,50000,0.245428765,203725.101,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,p1024,1024,REDUCE,50000,0.263722246,189593.410,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,p1024,1024,MODMUL,781,0.000302165,2584680.548,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,p1024,1024,MODMUL,781,0.000101123,7723267.866,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,p1024,1024,MODMUL,781,0.000139707,5590271.030,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,p1024,1024,MODMUL,50000,0.000299552,166915927.785,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,p1024,1024,MODMUL,50000,0.012105809,4130248.544,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,p1024,1024,MODMUL,50000,0.014255050,3507528.911,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,p1024,1024,MODMUL,50000,0.003977948,12569294.523,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,p1024,1024,MODMUL,50000,0.005858769,8534215.976,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,p1024,1024,MODMUL,50000,0.001681621,29733215.732,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,p1024,1024,MODMUL,50000,0.003555405,14063095.495,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,p1024,1024,MODMUL,50000,0.001671219,29918281.191,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,p1024,1024,MODMUL,50000,0.003487965,14335006.239,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,p1024,1024,MODMUL,50000,0.001650230,30298806.791,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,p1024,1024,MODMUL,50000,0.003414559,14643179.398,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,p1024,1024,MODMUL,50000,0.001645343,30388800.393,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,p1024,1024,MODMUL,50000,0.003474602,14390137.349,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,p1024,1024,MODMUL,50000,1.291344264,38719.342,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,p1024,1024,MODMUL,50000,1.300032063,38460.590,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,p1024,1024,MODEXP,195,0.042834313,4552.425,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,p1024,1024,MODEXP,195,0.006415093,30397.065,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,p1024,1024,MODEXP,195,0.008138394,23960.501,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,p1024,1024,MODEXP,50000,0.153570056,325584.305,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,p1024,1024,MODEXP,50000,4.990554169,10018.927,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,p1024,1024,MODEXP,50000,4.993293023,10013.432,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,p1024,1024,MODEXP,50000,0.568045363,88021.139,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,p1024,1024,MODEXP,50000,0.572662654,87311.438,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,p1024,1024,MODEXP,50000,0.185735888,269199.456,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,p1024,1024,MODEXP,50000,0.188251830,265601.668,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,p1024,1024,MODEXP,50000,0.113085548,442143.146,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,p1024,1024,MODEXP,50000,0.115350417,433461.805,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,p1024,1024,MODEXP,50000,0.184853904,270483.874,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,p1024,1024,MODEXP,50000,0.187253815,267017.257,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,p1024,1024,MODEXP,50000,0.112974767,442576.704,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,p1024,1024,MODEXP,50000,0.115235861,433892.710,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,p1024,1024,MODEXP,50000,0.000000000,inf,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,p1024,1024,MODEXP,50000,0.000000000,inf,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,p1024,1024,EXPONENTIATION,195,0.006745038,28910.141,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,p1024,1024,EXPONENTIATION,195,0.001239415,157332.290,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,p1024,1024,EXPONENTIATION,195,0.015514394,12568.973,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,p1024,1024,EXPONENTIATION,50000,2.511722100,19906.661,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,p1024,1024,EXPONENTIATION,50000,2.506918586,19944.804,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,p1024,1024,EXPONENTIATION,50000,0.603406153,82862.927,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,p1024,1024,EXPONENTIATION,50000,0.605692915,82550.082,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,p1024,1024,EXPONENTIATION,50000,0.159290740,313891.442,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,p1024,1024,EXPONENTIATION,50000,0.161209073,310156.240,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,p1024,1024,EXPONENTIATION,50000,0.144826814,345239.936,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,p1024,1024,EXPONENTIATION,50000,0.146781032,340643.470,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,p1024,1024,EXPONENTIATION,50000,0.155987790,320537.909,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,p1024,1024,EXPONENTIATION,50000,0.159056389,314353.924,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,p1024,1024,EXPONENTIATION,50000,0.149711034,333976.719,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,p1024,1024,EXPONENTIATION,50000,0.152273838,328355.814,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,p1024,1024,EXPONENTIATION,50000,0.000000000,inf,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,p1024,1024,EXPONENTIATION,50000,0.000000000,inf,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,p1024,1024,DIVIDE,1562,0.000110094,14187875.481,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,p1024,1024,DIVIDE,1562,0.000101272,15423808.811,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,p1024,1024,DIVIDE,1562,0.000094128,16594424.373,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,p1024,1024,DIVIDE,50000,0.000042112,1187310030.395,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,p1024,1024,DIVIDE,50000,0.025224160,1982226.564,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,p1024,1024,DIVIDE,50000,0.028240137,1770529.654,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,p1024,1024,DIVIDE,50000,0.005944317,8411395.288,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,p1024,1024,DIVIDE,50000,0.009081019,5505990.021,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,p1024,1024,DIVIDE,50000,0.001298633,38502024.922,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,p1024,1024,DIVIDE,50000,0.003849421,12988966.386,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,p1024,1024,DIVIDE,50000,0.001214611,41165443.146,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,p1024,1024,DIVIDE,50000,0.003980660,12560731.136,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,p1024,1024,DIVIDE,50000,0.001269085,39398464.308,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,p1024,1024,DIVIDE,50000,0.003787704,13200609.124,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,p1024,1024,DIVIDE,50000,0.001202879,41566940.752,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,p1024,1024,DIVIDE,50000,0.003698209,13520057.944,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,p1024,1024,DIVIDE,50000,0.492980606,101423.868,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,p1024,1024,DIVIDE,50000,0.507855734,98453.156,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,p1024,1024,ISQRT,390,0.000090001,4333285.274,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,p1024,1024,ISQRT,390,0.000055598,7014641.214,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,p1024,1024,ISQRT,50000,0.359993052,138891.569,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,p1024,1024,ISQRT,50000,0.360489550,138700.276,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,p1024,1024,ISQRT,50000,0.148128468,337544.840,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,p1024,1024,ISQRT,50000,0.150399063,332448.880,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,p1024,1024,ISQRT,50000,0.028555764,1750959.981,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,p1024,1024,ISQRT,50000,0.030484325,1640187.211,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,p1024,1024,ISQRT,50000,0.026542406,1883777.982,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,p1024,1024,ISQRT,50000,0.028622003,1746907.790,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,p1024,1024,ISQRT,50000,0.027603318,1811376.444,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,p1024,1024,ISQRT,50000,0.029937678,1670136.208,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,p1024,1024,ISQRT,50000,0.027195225,1838558.056,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,p1024,1024,ISQRT,50000,0.029279525,1707677.976,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,p1024,1024,ISQRT,50000,9.273683747,5391.601,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,p1024,1024,ISQRT,50000,9.172362562,5451.158,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,p1024,1024,MODMUL_R2,12500,0.005826330,2145432.888,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,p1024,1024,MODMUL_R2,12500,0.001245742,10034180.448,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,p1024,1024,MODMUL_R2,12500,0.001856220,6734115.561,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,p1024,1024,MODMUL_R2,50000,0.000066016,757392147.358,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,p1024,1024,MODMUL_R2,50000,0.005911296,8458382.058,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,p1024,1024,MODMUL_R2,50000,0.008062973,6201186.585,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,p1024,1024,MODMUL_R2,50000,0.000481096,103929361.700,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,p1024,1024,MODMUL_R2,50000,0.002359197,21193651.922,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,p1024,1024,MODMUL_R2,50000,0.000242703,206013109.591,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,p1024,1024,MODMUL_R2,50000,0.002114641,23644675.385,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,p1024,1024,MODMUL_R2,50000,0.000173659,287920582.653,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,p1024,1024,MODMUL_R2,50000,0.002012910,24839659.980,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,p1024,1024,MODMUL_R2,50000,0.000205495,243314924.604,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,p1024,1024,MODMUL_R2,50000,0.002029709,24634073.188,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,p1024,1024,MODMUL_R2,50000,0.000136339,366732922.639,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,p1024,1024,MODMUL_R2,50000,0.001943999,25720177.850,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,p1024,1024,MODMUL_R2,50000,0.513410748,97387.911,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,p1024,1024,MODMUL_R2,50000,0.525334922,95177.377,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,p2048,2048,ADD,6250,0.000160941,38834107.103,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,p2048,2048,ADD,6250,0.000055207,113210278.385,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,p2048,2048,ADD,6250,0.000128329,48702943.106,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,p2048,2048,ADD,50000,0.000025248,1980354879.594,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,p2048,2048,ADD,50000,0.000491547,101719672.430,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,p2048,2048,ADD,50000,0.004493200,11127926.649,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,p2048,2048,ADD,50000,0.000252303,198174418.856,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,p2048,2048,ADD,50000,0.004240559,11790898.325,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32,p2048,2048,ADD,50000,0.000088640,564079414.778,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32,p2048,2048,ADD,50000,0.004236728,11801560.071,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,p2048,2048,ADD,50000,0.000087922,568685888.534,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,p2048,2048,ADD,50000,0.003856047,12966646.937,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,p2048,2048,ADD,50000,0.000095765,522111413.787,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,p2048,2048,ADD,50000,0.003976149,12574981.480,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,p2048,2048,ADD,50000,0.000038105,1312163838.727,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,p2048,2048,ADD,50000,0.003927242,12731581.095,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,p2048,2048,ADD,50000,0.000043284,1155161194.895,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,p2048,2048,ADD,50000,0.003981381,12558456.476,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,p2048,2048,ADD,50000,0.001940838,25762067.733,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,p2048,2048,ADD,50000,0.006109731,8183666.348,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,p2048,2048,SUBTRACT,6250,0.000115184,54261008.188,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,p2048,2048,SUBTRACT,6250,0.000097977,63790483.071,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,p2048,2048,SUBTRACT,6250,0.000130142,48024465.000,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,p2048,2048,SUBTRACT,50000,0.000024768,2018733850.129,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,p2048,2048,SUBTRACT,50000,0.000516897,96731070.530,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,p2048,2048,SUBTRACT,50000,0.004448349,11240125.264,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,p2048,2048,SUBTRACT,50000,0.000252514,198008822.686,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,p2048,2048,SUBTRACT,50000,0.004410471,11336657.696,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32,p2048,2048,SUBTRACT,50000,0.000088315,566155247.970,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32,p2048,2048,SUBTRACT,50000,0.004030140,12406516.894,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,p2048,2048,SUBTRACT,50000,0.000088794,563101107.417,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,p2048,2048,SUBTRACT,50000,0.003740658,13366632.283,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,p2048,2048,SUBTRACT,50000,0.000094328,530065308.261,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,p2048,2048,SUBTRACT,50000,0.004191567,11928713.054,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,p2048,2048,SUBTRACT,50000,0.000038143,1310856514.217,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,p2048,2048,SUBTRACT,50000,0.003805504,13138864.122,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,p2048,2048,SUBTRACT,50000,0.000043544,1148263784.801,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,p2048,2048,SUBTRACT,50000,0.003985807,12544511.057,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,p2048,2048,SUBTRACT,50000,0.002128714,23488359.631,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,p2048,2048,SUBTRACT,50000,0.006210533,8050838.793,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,p2048,2048,ADDMOD,6250,0.000370242,16880850.907,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,p2048,2048,ADDMOD,6250,0.000118094,52923941.199,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,p2048,2048,ADDMOD,6250,0.000274272,22787597.802,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,p2048,2048,ADDMOD,50000,0.000024672,2026588845.655,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,p2048,2048,ADDMOD,50000,0.000679938,73536116.643,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,p2048,2048,ADDMOD,50000,0.004743892,10539868.957,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,p2048,2048,ADDMOD,50000,0.000342979,145781519.604,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,p2048,2048,ADDMOD,50000,0.004232595,11813083.931,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32,p2048,2048,ADDMOD,50000,0.000109078,458387580.964,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32,p2048,2048,ADDMOD,50000,0.004069630,12286129.211,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,p2048,2048,ADDMOD,50000,0.000118604,421570932.586,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,p2048,2048,ADDMOD,50000,0.003934367,12708524.651,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,p2048,2048,ADDMOD,50000,0.000123175,405926529.550,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,p2048,2048,ADDMOD,50000,0.004357357,11474845.875,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,p2048,2048,ADDMOD,50000,0.000037095,1347890545.747,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,p2048,2048,ADDMOD,50000,0.003926729,12733244.383,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,p2048,2048,ADDMOD,50000,0.000043059,1161197406.684,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,p2048,2048,ADDMOD,50000,0.004004075,12487278.588,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,p2048,2048,ADDMOD,50000,0.002246940,22252485.611,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,p2048,2048,ADDMOD,50000,0.006560484,7621388.912,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,p2048,2048,SUBTRACTMOD,6250,0.000312826,19979157.770,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,p2048,2048,SUBTRACTMOD,6250,0.000103734,60250255.036,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,p2048,2048,SUBTRACTMOD,6250,0.000241549,25874667.464,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,p2048,2048,SUBTRACTMOD,50000,0.000024352,2053219448.095,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,p2048,2048,SUBTRACTMOD,50000,0.000719636,69479570.382,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,p2048,2048,SUBTRACTMOD,50000,0.004642164,10770838.773,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,p2048,2048,SUBTRACTMOD,50000,0.000375480,133162885.485,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,p2048,2048,SUBTRACTMOD,50000,0.004165336,12003833.548,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32,p2048,2048,SUBTRACTMOD,50000,0.000118355,422457855.626,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32,p2048,2048,SUBTRACTMOD,50000,0.004074907,12270218.684,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,p2048,2048,SUBTRACTMOD,50000,0.000118515,421887526.618,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,p2048,2048,SUBTRACTMOD,50000,0.004023252,12427757.444,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,p2048,2048,SUBTRACTMOD,50000,0.000123903,403541471.783,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,p2048,2048,SUBTRACTMOD,50000,0.004355420,11479949.116,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,p2048,2048,SUBTRACTMOD,50000,0.000036852,1356778415.213,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,p2048,2048,SUBTRACTMOD,50000,0.003974499,12580201.935,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,p2048,2048,SUBTRACTMOD,50000,0.000044412,1125821813.357,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,p2048,2048,SUBTRACTMOD,50000,0.004032321,12399806.471,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,p2048,2048,SUBTRACTMOD,50000,0.002397391,20856005.556,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,p2048,2048,SUBTRACTMOD,50000,0.006491857,7701956.468,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.002532887,2467540.005,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.000553763,11286416.781,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.000629337,9931086.168,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.129209249,386969.202,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.134783014,370966.626,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.034220481,1461113.302,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.039896172,1253253.069,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.006638922,7531343.190,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.011971943,4176431.511,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.000635064,78732222.307,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.005916253,8451295.100,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.000632913,78999799.212,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.006098934,8198153.975,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.000507638,98495384.454,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.005885331,8495698.882,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.000522862,95627527.078,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.005962879,8385211.242,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.450234104,111053.338,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.459320824,108856.375,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.002530166,2470193.656,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.000558665,11187384.202,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.000680756,9180969.364,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.000074016,675529615.218,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.041343144,1209390.365,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.046915191,1065752.882,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.010492295,4765401.660,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.015976388,3129618.535,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.002649663,18870324.252,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.007905493,6324716.245,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.002687663,18603522.838,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.008009685,6242442.743,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.002690517,18583788.904,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.008391414,5958471.363,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.001628608,30701064.909,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.007020341,7122161.160,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.001538349,32502377.570,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.007061066,7081083.791,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.323282436,154663.522,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.327129331,152844.748,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.018280592,341892.648,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.002552372,2448702.619,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.001155761,5407692.414,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000115296,433666389.120,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.047678371,1048693.547,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.052002001,961501.462,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.001990731,25116401.991,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.005683266,8797758.189,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000400446,124860779.463,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.004384908,11402747.790,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000456249,109589281.436,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.004381131,11412578.171,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000379562,131730784.576,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.004538064,11017914.255,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000388371,128742876.815,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.004475413,11172153.275,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000317640,157410906.080,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.004336290,11530594.126,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,50000,2.071480246,24137.329,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,50000,2.057316427,24303.505,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,p2048,2048,COMPARE,6250,0.000038760,161248705.064,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,p2048,2048,COMPARE,6250,0.000049607,125990286.676,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,p2048,2048,COMPARE,6250,0.000228512,27350861.342,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,p2048,2048,COMPARE,50000,0.000024672,2026588845.655,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,p2048,2048,COMPARE,50000,0.000469898,106406070.872,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,p2048,2048,COMPARE,50000,0.004408422,11341926.887,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,p2048,2048,COMPARE,50000,0.000226184,221058959.083,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,p2048,2048,COMPARE,50000,0.004127715,12113239.414,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,p2048,2048,COMPARE,50000,0.000063776,783993951.792,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,p2048,2048,COMPARE,50000,0.003960766,12623820.745,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,p2048,2048,COMPARE,50000,0.000070758,706633876.658,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,p2048,2048,COMPARE,50000,0.004730708,10569242.493,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,p2048,2048,COMPARE,50000,0.000021850,2288329621.074,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,p2048,2048,COMPARE,50000,0.003874608,12904531.249,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,p2048,2048,COMPARE,50000,0.000026251,1904689328.060,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,p2048,2048,COMPARE,50000,0.003981730,12557355.716,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,p2048,2048,COMPARE,50000,0.001247062,40094237.541,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,p2048,2048,COMPARE,50000,0.005518370,9060646.534,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,p2048,2048,REDUCE,781,0.000019497,40057441.384,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,p2048,2048,REDUCE,781,0.000031390,24880536.641,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,p2048,2048,REDUCE,781,0.000063162,12365029.862,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,p2048,2048,REDUCE,50000,0.000031136,1605858170.606,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,p2048,2048,REDUCE,50000,0.235940850,211917.521,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,p2048,2048,REDUCE,50000,0.240465417,207930.107,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,p2048,2048,REDUCE,50000,0.005619350,8897826.258,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,p2048,2048,REDUCE,50000,0.009782088,5111383.173,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,p2048,2048,REDUCE,50000,0.001180787,42344639.648,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,p2048,2048,REDUCE,50000,0.005111001,9782819.453,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,p2048,2048,REDUCE,50000,0.001017804,49125371.738,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,p2048,2048,REDUCE,50000,0.005194483,9625597.008,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,p2048,2048,REDUCE,50000,0.001137619,43951446.046,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,p2048,2048,REDUCE,50000,0.005184623,9643902.752,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,p2048,2048,REDUCE,50000,0.000975768,51241688.694,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,p2048,2048,REDUCE,50000,0.004958979,10082720.658,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,p2048,2048,REDUCE,50000,1.066104695,46899.709,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,p2048,2048,REDUCE,50000,1.130038679,44246.273,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,p2048,2048,MODMUL,390,0.000504268,773398.269,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,p2048,2048,MODMUL,390,0.000103474,3769062.699,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,p2048,2048,MODMUL,390,0.000223545,1744615.198,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,p2048,2048,MODMUL,50000,0.000569632,87775967.642,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,p2048,2048,MODMUL,50000,0.485562306,102973.397,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,p2048,2048,MODMUL,50000,0.490851088,101863.887,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,p2048,2048,MODMUL,50000,0.019329646,2586700.242,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,p2048,2048,MODMUL,50000,0.023477635,2129686.401,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,p2048,2048,MODMUL,50000,0.006075082,8230341.581,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,p2048,2048,MODMUL,50000,0.010090450,4955180.394,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,p2048,2048,MODMUL,50000,0.013049275,3831630.493,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,p2048,2048,MODMUL,50000,0.017434531,2867871.811,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,p2048,2048,MODMUL,50000,0.006026698,8296417.043,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,p2048,2048,MODMUL,50000,0.010155918,4923237.860,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,p2048,2048,MODMUL,50000,0.012999981,3846159.468,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,p2048,2048,MODMUL,50000,0.017388896,2875398.185,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,p2048,2048,MODMUL,50000,5.122737479,9760.406,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,p2048,2048,MODMUL,50000,5.135432054,9736.279,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,p2048,2048,MODEXP,97,0.148327541,653.958,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,p2048,2048,MODEXP,97,0.025013297,3877.937,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,p2048,2048,MODEXP,97,0.028378566,3418.073,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,p2048,2048,MODEXP,50000,0.820030451,60973.345,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,p2048,2048,MODEXP,50000,64.846482089,771.052,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,p2048,2048,MODEXP,50000,64.847892070,771.035,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,p2048,2048,MODEXP,50000,12.057345529,4146.850,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,p2048,2048,MODEXP,50000,12.067681075,4143.298,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,p2048,2048,MODEXP,50000,1.431721144,34923.002,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,p2048,2048,MODEXP,50000,1.436246964,34812.954,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,p2048,2048,MODEXP,50000,2.918002941,17135.007,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,p2048,2048,MODEXP,50000,2.919289829,17127.453,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,p2048,2048,MODEXP,50000,1.434516254,34854.955,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,p2048,2048,MODEXP,50000,1.438719743,34753.120,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,p2048,2048,MODEXP,50000,2.928314310,17074.670,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,p2048,2048,MODEXP,50000,2.941887265,16995.893,0
opencl-kernel,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,p2048,2048,MODEXP,50000,0.000000000,inf,0
opencl-e2e,cpu-haswell-13th Gen Intel(R) Core(TM) i5-13500T,CPU,w8,p2048,2048,MODEXP,50000,0.000000000,inf,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,p2048,2048,EXPONENTIATION,97,0.018627066,5207.476,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,p2048,2048,EXPONENTIATION,97,0.003151117,30782.735,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,p2048,2048,EXPONENTIATION,97,0.039928247,2429.358,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,p2048,2048,EXPONENTIATION,50000,24.159079078,2069.615,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,p2048,2048,EXPONENTIATION,50000,24.153785267,2070.069,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,p2048,2048,EXPONENTIATION,50000,4.866796647,10273.698,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,p2048,2048,EXPONENTIATION,50000,4.864582234,10278.375,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,p2048,2048,EXPONENTIATION,50000,1.181776344,42309.190,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,p2048,2048,EXPONENTIATION,50000,1.186961439,42124.368,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,p2048,2048,EXPONENTIATION,50000,1.097086664,45575.251,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,p2048,2048,EXPONENTIATION,50000,1.101118973,45408.354,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,p2048,2048,EXPONENTIATION,50000,1.190769787,41989.644,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,p2048,2048,EXPONENTIATION,50000,1.195587840,41820.432,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,p2048,2048,EXPONENTIATION,50000,1.048837219,47671.840,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,p2048,2048,EXPONENTIATION,50000,1.052103901,47523.823,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,p2048,2048,DIVIDE,781,0.000069401,11253440.518,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,p2048,2048,DIVIDE,781,0.000048061,16250182.166,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,p2048,2048,DIVIDE,781,0.000102205,7641504.949,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,p2048,2048,DIVIDE,50000,0.000039264,1273431132.844,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,p2048,2048,DIVIDE,50000,1.229164252,40678.046,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,p2048,2048,DIVIDE,50000,1.232048170,40582.829,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,p2048,2048,DIVIDE,50000,0.100569485,497168.699,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,p2048,2048,DIVIDE,50000,0.107553953,464882.960,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,p2048,2048,DIVIDE,50000,0.015285019,3271176.830,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,p2048,2048,DIVIDE,50000,0.021029990,2377557.003,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,p2048,2048,DIVIDE,50000,0.014601525,3424299.859,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,p2048,2048,DIVIDE,50000,0.019946151,2506749.297,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,p2048,2048,DIVIDE,50000,0.014989012,3335776.901,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,p2048,2048,DIVIDE,50000,0.020374924,2453996.884,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,p2048,2048,DIVIDE,50000,0.014096664,3546938.482,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,p2048,2048,DIVIDE,50000,0.019590732,2552227.247,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,p2048,2048,ISQRT,195,0.000071475,2728226.569,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,p2048,2048,ISQRT,195,0.000060862,3203969.669,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,p2048,2048,ISQRT,50000,17.188069673,2908.994,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,p2048,2048,ISQRT,50000,17.176441175,2910.964,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,p2048,2048,ISQRT,50000,2.918637365,17131.282,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,p2048,2048,ISQRT,50000,2.924460417,17097.171,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,p2048,2048,ISQRT,50000,0.178391613,280282.235,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,p2048,2048,ISQRT,50000,0.182692434,273684.021,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,p2048,2048,ISQRT,50000,0.076003468,657864.717,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,p2048,2048,ISQRT,50000,0.080313719,622558.644,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,p2048,2048,ISQRT,50000,0.180853480,276466.895,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,p2048,2048,ISQRT,50000,0.184337423,271241.722,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,p2048,2048,ISQRT,50000,0.075149427,665341.068,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,p2048,2048,ISQRT,50000,0.079340809,630192.717,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-1t,p2048,2048,MODMUL_R2,6250,0.009417199,663679.296,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,gmp-nt,p2048,2048,MODMUL_R2,6250,0.001213566,5150111.325,0
library,13th Gen Intel(R) Core(TM) i5-13500T,host-cpu,openssl-nt,p2048,2048,MODMUL_R2,6250,0.002904067,2152154.202,0
library,NVIDIA GeForce RTX 5080,gpu,cgbn,p2048,2048,MODMUL_R2,50000,0.000222048,225176538.406,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w8,p2048,2048,MODMUL_R2,50000,0.033796699,1479434.426,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w8,p2048,2048,MODMUL_R2,50000,0.038193213,1309133.117,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w16,p2048,2048,MODMUL_R2,50000,0.001954212,25585760.412,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w16,p2048,2048,MODMUL_R2,50000,0.005974661,8368675.642,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-opt,p2048,2048,MODMUL_R2,50000,0.000818434,61092281.355,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-opt,p2048,2048,MODMUL_R2,50000,0.004478357,11164808.877,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-o64,p2048,2048,MODMUL_R2,50000,0.000653251,76540258.131,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-o64,p2048,2048,MODMUL_R2,50000,0.004298089,11633076.924,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il,p2048,2048,MODMUL_R2,50000,0.000758412,65927226.941,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il,p2048,2048,MODMUL_R2,50000,0.004414213,11327047.428,0
opencl-kernel,NVIDIA GeForce RTX 5080,GPU,w32-il64,p2048,2048,MODMUL_R2,50000,0.000582407,85850616.573,0
opencl-e2e,NVIDIA GeForce RTX 5080,GPU,w32-il64,p2048,2048,MODMUL_R2,50000,0.004485033,11148189.997,0
```
