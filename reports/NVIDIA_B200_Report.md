# MPA-OpenCL benchmark report - NVIDIA B200

> **Partial report.** The run was interrupted or hit its time budget.
> Rows that never ran are marked `n/a`.

## 1. System under test

2 OpenCL device(s) exercised with the identical kernels and operands.

### Device 0 - NVIDIA B200 (GPU)

| Property | Value |
|---|---|
| Model | NVIDIA B200 |
| Type | GPU |
| Vendor | NVIDIA Corporation |
| Device memory | 178.34 GiB |
| Max single allocation | 44.59 GiB |
| Local memory | 48 KiB |
| Global cache | 4736 KiB |
| Compute units | 148 |
| Max clock | 1965 MHz |
| Max work-group size | 1024 |
| OpenCL version | OpenCL 3.0 CUDA |
| Driver | 595.91.07 |

### Device 1 - cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C (CPU)

| Property | Value |
|---|---|
| Model | cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C |
| Type | CPU |
| Vendor | GenuineIntel |
| Device memory | 1994.00 GiB |
| Max single allocation | 512.00 GiB |
| Local memory | 2048 KiB |
| Global cache | 327680 KiB |
| Compute units | 192 |
| Max clock | 3888 MHz |
| Max work-group size | 4096 |
| OpenCL version | OpenCL 3.0 PoCL HSTR: cpu-x86_64-pc-linux-gnu-skylake-avx512 |
| Driver | 5.0+debian |

### Host

| Property | Value |
|---|---|
| CPU | Intel(R) Xeon(R) Platinum 8559C |
| Logical cores | 192 |
| OpenMP threads used | 192 |
| RAM | 1996.0 GB |
| OS | Ubuntu 24.04.4 LTS |
| Kernel | 7.0.0-1011-aws |
| Arch | x86_64 |
| GMP | 6.3.0 |
| OpenSSL | OpenSSL 3.0.13 30 Jan 2024 |
| CGBN | cgbn_results.tsv loaded |

## 2. Method

- Workload auto-sized from the device and host: --min-items from 700 x compute units, --items from ten times that capped by host RAM. Either flag, given explicitly, overrides its half.
- Base workload 50000 items, scaled down per operator by its cost weight and by modulus size. Device rows honour --min-items (103600) so the GPU is not left idle; the CPU libraries keep the smaller count because a full-width MODEXP there costs minutes. Both counts appear in every row as dev/cpu, and throughput is per-second so they remain comparable.
- 5 timed repetitions, **minimum** reported. Two untimed warm-up launches precede them.
- `kernel` times `clEnqueueNDRangeKernel` + `clFinish` only. `e2e` adds the host->device operand writes and the device->host result read.
- Every OpenCL device runs the same kernels on the same operands, so GPU and CPU-OpenCL columns are directly comparable.
- CPU library baselines (GMP, OpenSSL) run those same operands, with every temporary - including each thread's GMP context, BN_CTX and Montgomery context - allocated outside the timed region, so the figure is the arithmetic and not marshalling. The generator is reseeded per modulus and operation so every backend sees identical inputs.
- Cost weighting drives the wide cells down to a few hundred items, which is tens of microseconds of work - the same order as the cost of entering an OpenMP region. Each baseline pass is therefore repeated until the timed interval reaches 5 ms and the per-pass time is reported; the multi-threaded loop enters one parallel region per interval and partitions the range itself. Without this the multi-threaded GMP figure came out up to 9x slower than the single-threaded one at 2048 bits.
- OpenSSL rows time the nearest BN primitive, which is not always semantically identical (its Montgomery routine expects Montgomery-domain inputs); they measure comparable work, not identical results. Correctness is judged against GMP only.
- Every device cell is checked word-for-word against GMP before it is timed. A cell that mismatches is reported and excluded from the speedup tables.
- Total wall time 6552.5 s.

## 3. Correctness

| Device | Kernel | Configs run | Passed | Mismatched | Launch failed |
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

### Device 0 - NVIDIA B200 (GPU)

#### secp256k1 (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | GMP 192T | OpenSSL 192T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | 2.05 G | 2.87 G | 3.43 G | 3.58 G | 3.65 G | 3.93 G | 3.64 G | 65.93 M | 19.22 M | 23.01 M | 5.11 G |
| SUBTRACT | 50000 / 50000 | 2.16 G | 2.86 G | 3.57 G | 3.54 G | 3.64 G | 4.02 G | 3.95 G | 107.99 M | 28.58 M | 54.31 M | 6.35 G |
| ADDMOD | 50000 / 50000 | 1.63 G | 2.40 G | 3.11 G | 3.56 G | 3.47 G | 4.17 G | 4.14 G | 31.74 M | 35.90 M | 351.46 M | 5.07 G |
| SUBTRACTMOD | 50000 / 50000 | 1.65 G | 2.32 G | 3.37 G | 3.50 G | 3.50 G | 4.06 G | 4.09 G | 36.09 M | 1.87 G | 35.75 M | 4.61 G |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | 29.83 M | 84.20 M | 640.28 M | 2.75 G | 3.02 G | 3.39 G | 3.76 G | 76.07 M | 39.46 M | 43.10 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | 289.55 M | 836.58 M | 2.06 G | 2.14 G | 2.14 G | 2.35 G | 2.47 G | 74.59 M | 61.92 M | 60.95 M | 5.62 G |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | 655.03 M | 1.72 G | 3.37 G | 2.96 G | 3.28 G | 3.03 G | 3.68 G | 8.45 M | 13.99 M | 6.58 M | 4.52 G |
| COMPARE | 50000 / 50000 | 2.12 G | 2.77 G | - | 3.68 G | 3.69 G | 3.92 G | 3.82 G | 169.64 M | 65.63 M | 5.47 G | 5.81 G |
| REDUCE | 50000 / 6250 | 352.22 M | 519.36 M | - | 1.32 G | 1.33 G | 1.33 G | 1.37 G | 81.59 M | 58.87 M | 130.22 M | 3.90 G |
| MODMUL | 50000 / 3125 | 132.09 M | 221.12 M | - | 457.61 M | 577.45 M | 465.83 M | 539.92 M | 16.08 M | 752.45 M | 151.09 M | 1.64 G |
| MODEXP | 50000 / 781 | 2.91 M | 15.53 M | - | 19.61 M | 40.69 M | 19.42 M | 40.38 M | 132.71 k | 7.10 M | 4.43 M | 6.16 M |
| EXPONENTIATION | 50000 / 781 | 2.13 M | 10.64 M | - | 106.02 M | 187.22 M | 105.55 M | 188.64 M | 405.63 k | 15.09 M | 1.50 M | n/a |
| DIVIDE | 50000 / 6250 | 145.82 M | 171.30 M | - | 449.37 M | 457.15 M | 459.07 M | 476.71 M | 43.06 M | 1.62 G | 304.30 M | 3.22 G |
| ISQRT | 50000 / 1562 | 13.46 M | 18.01 M | - | 73.23 M | 79.50 M | 71.47 M | 80.42 M | 18.89 M | 853.55 M | n/a | n/a |
| MODMUL_R2 | 50000 / 50000 | 552.56 M | 1.67 G | - | 2.31 G | 2.83 G | 2.32 G | 3.16 G | 16.07 M | 1.01 G | 27.25 M | 3.50 G |

#### rsa256(composite) (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | GMP 192T | OpenSSL 192T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | 2.13 G | 2.73 G | 3.86 G | 3.75 G | 3.72 G | 3.80 G | 3.64 G | 72.86 M | 3.08 G | 2.96 G | 6.23 G |
| SUBTRACT | 50000 / 50000 | 2.15 G | 2.85 G | 3.61 G | 3.65 G | 3.59 G | 3.77 G | 3.82 G | 107.83 M | 4.53 G | 3.89 G | 5.94 G |
| ADDMOD | 50000 / 50000 | 1.73 G | 2.44 G | 3.47 G | 3.52 G | 3.56 G | 4.16 G | 4.27 G | 36.90 M | 1.90 G | 396.10 M | 4.57 G |
| SUBTRACTMOD | 50000 / 50000 | 1.62 G | 2.36 G | 3.39 G | 3.58 G | 3.43 G | 3.99 G | 4.19 G | 35.76 M | 1.86 G | 249.41 M | 4.40 G |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | 29.06 M | 84.11 M | 618.45 M | 2.82 G | 2.99 G | 3.46 G | 3.59 G | 75.57 M | 3.01 G | 2.09 G | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | 288.68 M | 827.38 M | 2.05 G | 2.09 G | 2.09 G | 2.34 G | 2.51 G | 75.20 M | 2.82 G | 2.10 G | 5.66 G |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | 661.49 M | 1.74 G | 3.20 G | 2.89 G | 3.22 G | 2.97 G | 3.51 G | 8.40 M | 544.70 M | 1.07 G | 4.60 G |
| COMPARE | 50000 / 50000 | 2.13 G | 2.83 G | - | 3.71 G | 3.74 G | 3.91 G | 4.11 G | 169.17 M | 5.47 G | 5.45 G | 5.74 G |
| REDUCE | 50000 / 6250 | 353.59 M | 523.42 M | - | 1.33 G | 1.36 G | 1.31 G | 1.38 G | 48.82 M | 2.30 G | 25.41 M | 3.93 G |
| MODMUL | 50000 / 3125 | 131.30 M | 221.20 M | - | 458.85 M | 573.11 M | 466.19 M | 538.01 M | 16.21 M | 715.28 M | 156.02 M | 1.59 G |
| MODEXP | 50000 / 781 | 2.89 M | 15.56 M | - | 19.73 M | 40.80 M | 19.51 M | 40.62 M | 143.59 k | 7.28 M | 4.41 M | 6.26 M |
| EXPONENTIATION | 50000 / 781 | 2.13 M | 10.58 M | - | 105.77 M | 187.14 M | 105.83 M | 188.95 M | 408.84 k | 15.44 M | 979.95 k | n/a |
| DIVIDE | 50000 / 6250 | 140.32 M | 165.70 M | - | 419.36 M | 429.60 M | 426.67 M | 442.60 M | 42.54 M | 33.11 M | 308.38 M | 3.07 G |
| ISQRT | 50000 / 1562 | 13.45 M | 17.49 M | - | 67.79 M | 71.16 M | 66.77 M | 71.94 M | 18.96 M | 7.92 M | n/a | n/a |
| MODMUL_R2 | 50000 / 50000 | 552.85 M | 1.70 G | - | 2.25 G | 2.87 G | 2.33 G | 3.00 G | 16.28 M | 997.76 M | 158.35 M | 3.33 G |

#### brainpoolP512r1 (512-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | GMP 192T | OpenSSL 192T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 25000 | 857.00 M | 1.50 G | 2.46 G | 2.46 G | 2.50 G | 3.46 G | 2.60 G | 68.50 M | 2.96 G | 2.82 G | 5.70 G |
| SUBTRACT | 50000 / 25000 | 876.65 M | 1.47 G | 2.53 G | 2.29 G | 2.58 G | 3.19 G | 2.96 G | 101.32 M | 3.70 G | 3.28 G | 5.64 G |
| ADDMOD | 50000 / 25000 | 691.28 M | 1.27 G | 2.30 G | 2.39 G | 2.42 G | 3.67 G | 3.48 G | 34.20 M | 34.42 M | 399.36 M | 4.49 G |
| SUBTRACTMOD | 50000 / 25000 | 642.60 M | 1.21 G | 2.20 G | 2.37 G | 2.46 G | 3.58 G | 3.42 G | 32.47 M | 1.69 G | 245.72 M | 4.54 G |
| MULTIPLYOPERANDSCANNING | 50000 / 25000 | 10.25 M | 25.10 M | 59.12 M | 1.35 G | 1.38 G | 1.88 G | 1.69 G | 30.71 M | 1.53 G | 1.57 G | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 25000 | 43.60 M | 156.32 M | 536.23 M | 540.21 M | 543.82 M | 1.02 G | 870.69 M | 30.93 M | 1.39 G | 1.44 G | 4.93 G |
| MONTGOMERYMULTIPLICATION | 50000 / 25000 | 198.04 M | 637.41 M | 1.91 G | 1.30 G | 1.95 G | 1.52 G | 2.27 G | 3.47 M | 232.06 M | 513.26 M | 3.40 G |
| COMPARE | 50000 / 25000 | 1.19 G | 1.84 G | - | 3.11 G | 3.01 G | 3.32 G | 3.31 G | 118.14 M | 5.33 G | 5.40 G | 5.39 G |
| REDUCE | 50000 / 3125 | 124.92 M | 151.90 M | - | 522.48 M | 544.35 M | 543.34 M | 469.73 M | 45.95 M | 2.23 G | 256.07 M | 2.57 G |
| MODMUL | 50000 / 1562 | 43.16 M | 60.83 M | - | 148.06 M | 172.74 M | 150.28 M | 148.36 M | 7.34 M | 372.28 M | 103.58 M | 591.41 M |
| MODEXP | 50000 / 390 | 225.63 k | 2.22 M | - | 2.52 M | 5.86 M | 2.50 M | 5.75 M | 26.73 k | 1.33 M | 1.42 M | 2.05 M |
| EXPONENTIATION | 50000 / 390 | 263.68 k | 1.05 M | - | 4.36 M | 6.15 M | 4.40 M | 4.23 M | 115.50 k | 3.55 M | 365.31 k | n/a |
| DIVIDE | 50000 / 3125 | 42.76 M | 42.33 M | - | 125.34 M | 135.38 M | 129.49 M | 143.86 M | 39.62 M | 775.41 M | 230.02 M | 1.96 G |
| ISQRT | 50000 / 781 | 2.85 M | 2.82 M | - | 14.64 M | 16.03 M | 15.43 M | 12.31 M | 11.12 M | 500.64 M | n/a | n/a |
| MODMUL_R2 | 50000 / 25000 | 126.77 M | 641.77 M | - | 861.76 M | 1.54 G | 942.35 M | 1.58 G | 7.47 M | 496.72 M | 112.27 M | 2.34 G |

#### p1024 (1024-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | GMP 192T | OpenSSL 192T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 12500 | 259.04 M | 494.45 M | 1.08 G | 1.17 G | 1.17 G | 2.20 G | 1.78 G | 48.98 M | 43.39 M | 1.36 G | 4.57 G |
| SUBTRACT | 50000 / 12500 | 254.61 M | 488.24 M | 1.08 G | 1.16 G | 1.19 G | 1.79 G | 1.79 G | 69.03 M | 3.24 G | 209.09 M | 4.49 G |
| ADDMOD | 50000 / 12500 | 201.12 M | 383.61 M | 954.98 M | 1.23 G | 1.23 G | 2.96 G | 2.97 G | 23.45 M | 1.11 G | 116.44 M | 3.98 G |
| SUBTRACTMOD | 50000 / 12500 | 189.99 M | 388.17 M | 973.37 M | 1.25 G | 1.21 G | 2.84 G | 3.07 G | 28.18 M | 1.43 G | 94.05 M | 3.78 G |
| MULTIPLYOPERANDSCANNING | 50000 / 12500 | 2.19 M | 7.28 M | 11.46 M | 511.81 M | 610.58 M | 567.83 M | 767.18 M | 8.46 M | 463.51 M | 239.82 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 12500 | 5.86 M | 22.76 M | 87.94 M | 87.67 M | 88.18 M | 249.76 M | 267.12 M | 8.32 M | 451.97 M | 363.07 M | 2.51 G |
| MONTGOMERYMULTIPLICATION | 50000 / 12500 | 25.96 M | 196.45 M | 730.63 M | 430.20 M | 623.72 M | 505.23 M | 919.98 M | 1.12 M | 87.84 M | 171.23 M | 1.70 G |
| COMPARE | 50000 / 12500 | 575.10 M | 961.19 M | - | 2.22 G | 2.19 G | 2.84 G | 2.69 G | 123.98 M | 4.13 G | 3.17 G | 4.49 G |
| REDUCE | 50000 / 1562 | 22.52 M | 42.81 M | - | 176.93 M | 146.39 M | 155.82 M | 152.11 M | 71.51 M | 1.88 G | 102.64 M | 1.92 G |
| MODMUL | 50000 / 781 | 6.06 M | 16.21 M | - | 34.45 M | 45.70 M | 34.86 M | 45.83 M | 2.62 M | 125.40 M | 34.79 M | 235.10 M |
| MODEXP | 50000 / 195 | 19.65 k | 321.05 k | - | 318.87 k | 615.90 k | 313.65 k | 619.92 k | 4.21 k | 156.30 k | 193.17 k | 367.87 k |
| EXPONENTIATION | 50000 / 195 | 30.69 k | 134.46 k | - | 447.86 k | 519.55 k | 451.55 k | 517.34 k | 30.56 k | 636.93 k | 37.04 k | n/a |
| DIVIDE | 50000 / 1562 | 4.75 M | 5.62 M | - | 31.94 M | 31.41 M | 32.25 M | 31.05 M | 36.71 M | 1.03 G | 176.62 M | 1.47 G |
| ISQRT | 50000 / 390 | 257.67 k | 404.75 k | - | 2.17 M | 2.19 M | 2.13 M | 2.19 M | 6.09 M | 167.34 M | n/a | n/a |
| MODMUL_R2 | 50000 / 12500 | 20.00 M | 213.27 M | - | 268.47 M | 410.67 M | 293.26 M | 499.28 M | 2.55 M | 154.23 M | 53.06 M | 960.36 M |

#### p2048 (2048-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | GMP 192T | OpenSSL 192T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 6250 | 135.70 M | 268.89 M | 669.68 M | 676.53 M | 663.65 M | 1.19 G | 1.17 G | 33.75 M | 2.11 G | 1.90 G | 3.88 G |
| SUBTRACT | 50000 / 6250 | 135.87 M | 267.55 M | 665.49 M | 678.08 M | 656.69 M | 1.20 G | 1.19 G | 38.93 M | 2.27 G | 1.94 G | 3.83 G |
| ADDMOD | 50000 / 6250 | 107.52 M | 209.13 M | 560.63 M | 712.10 M | 700.77 M | 1.88 G | 1.97 G | 16.91 M | 982.24 M | 12.00 M | 3.84 G |
| SUBTRACTMOD | 50000 / 6250 | 104.87 M | 202.42 M | 536.86 M | 713.26 M | 695.41 M | 1.78 G | 1.89 G | 19.24 M | 1.08 G | 120.39 M | 3.55 G |
| MULTIPLYOPERANDSCANNING | 50000 / 6250 | 543.49 k | 2.09 M | 7.03 M | 179.32 M | 218.00 M | 161.14 M | 273.37 M | 2.52 M | 142.22 M | 3.41 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 6250 | 1.48 M | 5.83 M | 23.00 M | 23.03 M | 23.02 M | 59.05 M | 62.96 M | 2.52 M | 145.99 M | 132.53 M | 840.05 M |
| MONTGOMERYMULTIPLICATION | 50000 / 6250 | 1.25 M | 37.66 M | 244.10 M | 135.95 M | 198.11 M | 144.74 M | 219.29 M | 323.77 k | 29.33 M | 67.58 M | 564.08 M |
| COMPARE | 50000 / 6250 | 303.09 M | 576.89 M | - | 1.50 G | 1.43 G | 2.03 G | 1.93 G | 165.07 M | 4.48 G | 4.61 G | 3.80 G |
| REDUCE | 50000 / 781 | 458.87 k | 10.40 M | - | 55.35 M | 51.27 M | 46.82 M | 54.12 M | 48.73 M | 1.39 G | 121.64 M | 1.90 G |
| MODMUL | 50000 / 390 | 277.70 k | 3.25 M | - | 9.19 M | 13.44 M | 9.10 M | 13.94 M | 768.67 k | 1.15 M | 1.30 M | 121.83 M |
| MODEXP | 50000 / 97 | 919.6 | 9.39 k | - | 33.13 k | 51.40 k | 36.37 k | 50.38 k | 590.3 | 16.26 k | 11.22 k | 87.95 k |
| EXPONENTIATION | 50000 / 97 | 3.65 k | 15.61 k | - | 53.02 k | 61.72 k | 52.21 k | 61.64 k | 5.09 k | 28.32 k | 8.46 k | n/a |
| DIVIDE | 50000 / 781 | 112.26 k | 326.21 k | - | 2.98 M | 3.27 M | 2.98 M | 3.27 M | 28.44 M | 686.67 M | 132.08 M | 1.59 G |
| ISQRT | 50000 / 195 | 6.80 k | 11.22 k | - | 324.90 k | 925.21 k | 325.11 k | 923.05 k | 3.67 M | 83.89 M | n/a | n/a |
| MODMUL_R2 | 50000 / 6250 | 1.62 M | 44.28 M | - | 78.50 M | 113.07 M | 79.94 M | 129.16 M | 745.19 k | 57.25 M | 22.31 M | 298.87 M |

### Device 1 - cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C (CPU)

#### secp256k1 (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | GMP 192T | OpenSSL 192T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | 82.80 M | - | - | - | - | - | - | 65.93 M | 19.22 M | 23.01 M | 5.11 G |
| SUBTRACT | 50000 / 50000 | 65.86 M | - | - | - | - | - | - | 107.99 M | 28.58 M | 54.31 M | 6.35 G |
| ADDMOD | 50000 / 50000 | 39.67 M | - | - | - | - | - | - | 31.74 M | 35.90 M | 351.46 M | 5.07 G |
| SUBTRACTMOD | 50000 / 50000 | 48.89 M | - | - | - | - | - | - | 36.09 M | 1.87 G | 35.75 M | 4.61 G |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | 27.34 M | - | - | - | - | - | - | 76.07 M | 39.46 M | 43.10 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | 22.43 M | - | - | - | - | - | - | 74.59 M | 61.92 M | 60.95 M | 5.62 G |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | 16.23 M | - | - | - | - | - | - | 8.45 M | 13.99 M | 6.58 M | 4.52 G |
| COMPARE | 50000 / 50000 | 59.33 M | - | - | - | - | - | - | 169.64 M | 65.63 M | 5.47 G | 5.81 G |
| REDUCE | 50000 / 6250 | 19.43 M | - | - | - | - | - | - | 81.59 M | 58.87 M | 130.22 M | 3.90 G |
| MODMUL | 50000 / 3125 | 6.29 M | - | - | - | - | - | - | 16.08 M | 752.45 M | 151.09 M | 1.64 G |
| MODEXP | 50000 / 781 | 25.08 k | - | - | - | - | - | - | 132.71 k | 7.10 M | 4.43 M | 6.16 M |
| EXPONENTIATION | 50000 / 781 | 68.89 k | - | - | - | - | - | - | 405.63 k | 15.09 M | 1.50 M | n/a |
| DIVIDE | 50000 / 6250 | 13.35 M | - | - | - | - | - | - | 43.06 M | 1.62 G | 304.30 M | 3.22 G |
| ISQRT | 50000 / 1562 | 446.97 k | - | - | - | - | - | - | 18.89 M | 853.55 M | n/a | n/a |
| MODMUL_R2 | 50000 / 50000 | 15.76 M | - | - | - | - | - | - | 16.07 M | 1.01 G | 27.25 M | 3.50 G |

#### rsa256(composite) (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | GMP 192T | OpenSSL 192T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | 209.32 M | - | - | - | - | - | - | 72.86 M | 3.08 G | 2.96 G | 6.23 G |
| SUBTRACT | 50000 / 50000 | 76.24 M | - | - | - | - | - | - | 107.83 M | 4.53 G | 3.89 G | 5.94 G |
| ADDMOD | 50000 / 50000 | 59.86 M | - | - | - | - | - | - | 36.90 M | 1.90 G | 396.10 M | 4.57 G |
| SUBTRACTMOD | 50000 / 50000 | 44.90 M | - | - | - | - | - | - | 35.76 M | 1.86 G | 249.41 M | 4.40 G |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | 33.02 M | - | - | - | - | - | - | 75.57 M | 3.01 G | 2.09 G | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | 30.35 M | - | - | - | - | - | - | 75.20 M | 2.82 G | 2.10 G | 5.66 G |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | 11.18 M | - | - | - | - | - | - | 8.40 M | 544.70 M | 1.07 G | 4.60 G |
| COMPARE | 50000 / 50000 | 83.67 M | - | - | - | - | - | - | 169.17 M | 5.47 G | 5.45 G | 5.74 G |
| REDUCE | 50000 / 6250 | 14.33 M | - | - | - | - | - | - | 48.82 M | 2.30 G | 25.41 M | 3.93 G |
| MODMUL | 50000 / 3125 | 4.51 M | - | - | - | - | - | - | 16.21 M | 715.28 M | 156.02 M | 1.59 G |
| MODEXP | 50000 / 781 | 25.26 k | - | - | - | - | - | - | 143.59 k | 7.28 M | 4.41 M | 6.26 M |
| EXPONENTIATION | 50000 / 781 | 69.58 k | - | - | - | - | - | - | 408.84 k | 15.44 M | 979.95 k | n/a |
| DIVIDE | 50000 / 6250 | 10.41 M | - | - | - | - | - | - | 42.54 M | 33.11 M | 308.38 M | 3.07 G |
| ISQRT | 50000 / 1562 | 459.53 k | - | - | - | - | - | - | 18.96 M | 7.92 M | n/a | n/a |
| MODMUL_R2 | 50000 / 50000 | 11.99 M | - | - | - | - | - | - | 16.28 M | 997.76 M | 158.35 M | 3.33 G |

#### brainpoolP512r1 (512-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | GMP 192T | OpenSSL 192T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 25000 | 62.64 M | - | - | - | - | - | - | 68.50 M | 2.96 G | 2.82 G | 5.70 G |
| SUBTRACT | 50000 / 25000 | 52.11 M | - | - | - | - | - | - | 101.32 M | 3.70 G | 3.28 G | 5.64 G |
| ADDMOD | 50000 / 25000 | 38.06 M | - | - | - | - | - | - | 34.20 M | 34.42 M | 399.36 M | 4.49 G |
| SUBTRACTMOD | 50000 / 25000 | 51.87 M | - | - | - | - | - | - | 32.47 M | 1.69 G | 245.72 M | 4.54 G |
| MULTIPLYOPERANDSCANNING | 50000 / 25000 | 13.45 M | - | - | - | - | - | - | 30.71 M | 1.53 G | 1.57 G | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 25000 | 14.03 M | - | - | - | - | - | - | 30.93 M | 1.39 G | 1.44 G | 4.93 G |
| MONTGOMERYMULTIPLICATION | 50000 / 25000 | 4.79 M | - | - | - | - | - | - | 3.47 M | 232.06 M | 513.26 M | 3.40 G |
| COMPARE | 50000 / 25000 | 60.40 M | - | - | - | - | - | - | 118.14 M | 5.33 G | 5.40 G | 5.39 G |
| REDUCE | 50000 / 3125 | 7.11 M | - | - | - | - | - | - | 45.95 M | 2.23 G | 256.07 M | 2.57 G |
| MODMUL | 50000 / 1562 | 525.48 k | - | - | - | - | - | - | 7.34 M | 372.28 M | 103.58 M | 591.41 M |
| MODEXP | 50000 / 390 | 3.23 k | - | - | - | - | - | - | 26.73 k | 1.33 M | 1.42 M | 2.05 M |
| EXPONENTIATION | 50000 / 390 | 7.32 k | - | - | - | - | - | - | 115.50 k | 3.55 M | 365.31 k | n/a |
| DIVIDE | 50000 / 3125 | 613.33 k | - | - | - | - | - | - | 39.62 M | 775.41 M | 230.02 M | 1.96 G |
| ISQRT | 50000 / 781 | 81.67 k | - | - | - | - | - | - | 11.12 M | 500.64 M | n/a | n/a |
| MODMUL_R2 | 50000 / 25000 | 3.90 M | - | - | - | - | - | - | 7.47 M | 496.72 M | 112.27 M | 2.34 G |

#### p1024 (1024-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | GMP 192T | OpenSSL 192T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 12500 | 42.12 M | - | - | - | - | - | - | 48.98 M | 43.39 M | 1.36 G | 4.57 G |
| SUBTRACT | 50000 / 12500 | 49.23 M | - | - | - | - | - | - | 69.03 M | 3.24 G | 209.09 M | 4.49 G |
| ADDMOD | 50000 / 12500 | 50.35 M | - | - | - | - | - | - | 23.45 M | 1.11 G | 116.44 M | 3.98 G |
| SUBTRACTMOD | 50000 / 12500 | 53.24 M | - | - | - | - | - | - | 28.18 M | 1.43 G | 94.05 M | 3.78 G |
| MULTIPLYOPERANDSCANNING | 50000 / 12500 | 3.99 M | - | - | - | - | - | - | 8.46 M | 463.51 M | 239.82 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 12500 | 4.38 M | - | - | - | - | - | - | 8.32 M | 451.97 M | 363.07 M | 2.51 G |
| MONTGOMERYMULTIPLICATION | 50000 / 12500 | 275.27 k | - | - | - | - | - | - | 1.12 M | 87.84 M | 171.23 M | 1.70 G |
| COMPARE | 50000 / 12500 | 32.62 M | - | - | - | - | - | - | 123.98 M | 4.13 G | 3.17 G | 4.49 G |
| REDUCE | 50000 / 1562 | 506.05 k | - | - | - | - | - | - | 71.51 M | 1.88 G | 102.64 M | 1.92 G |
| MODMUL | 50000 / 781 | 121.27 k | - | - | - | - | - | - | 2.62 M | 125.40 M | 34.79 M | 235.10 M |
| MODEXP | 50000 / 195 | 360.6 | - | - | - | - | - | - | 4.21 k | 156.30 k | 193.17 k | 367.87 k |
| EXPONENTIATION | 50000 / 195 | 932.8 | - | - | - | - | - | - | 30.56 k | 636.93 k | 37.04 k | n/a |
| DIVIDE | 50000 / 1562 | 254.83 k | - | - | - | - | - | - | 36.71 M | 1.03 G | 176.62 M | 1.47 G |
| ISQRT | 50000 / 390 | 17.25 k | - | - | - | - | - | - | 6.09 M | 167.34 M | n/a | n/a |
| MODMUL_R2 | 50000 / 12500 | 256.12 k | - | - | - | - | - | - | 2.55 M | 154.23 M | 53.06 M | 960.36 M |

#### p2048 (2048-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | GMP 192T | OpenSSL 192T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 6250 | 58.91 M | - | - | - | - | - | - | 33.75 M | 2.11 G | 1.90 G | 3.88 G |
| SUBTRACT | 50000 / 6250 | 26.98 M | - | - | - | - | - | - | 38.93 M | 2.27 G | 1.94 G | 3.83 G |
| ADDMOD | 50000 / 6250 | 48.99 M | - | - | - | - | - | - | 16.91 M | 982.24 M | 12.00 M | 3.84 G |
| SUBTRACTMOD | 50000 / 6250 | 54.31 M | - | - | - | - | - | - | 19.24 M | 1.08 G | 120.39 M | 3.55 G |
| MULTIPLYOPERANDSCANNING | 50000 / 6250 | 419.60 k | - | - | - | - | - | - | 2.52 M | 142.22 M | 3.41 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 6250 | 485.60 k | - | - | - | - | - | - | 2.52 M | 145.99 M | 132.53 M | 840.05 M |
| MONTGOMERYMULTIPLICATION | 50000 / 6250 | 69.94 k | - | - | - | - | - | - | 323.77 k | 29.33 M | 67.58 M | 564.08 M |
| COMPARE | 50000 / 6250 | 49.17 M | - | - | - | - | - | - | 165.07 M | 4.48 G | 4.61 G | 3.80 G |
| REDUCE | 50000 / 781 | 128.32 k | - | - | - | - | - | - | 48.73 M | 1.39 G | 121.64 M | 1.90 G |
| MODMUL | 50000 / 390 | 29.27 k | - | - | - | - | - | - | 768.67 k | 1.15 M | 1.30 M | 121.83 M |
| MODEXP | 50000 / 97 | over budget | - | - | - | - | - | - | 590.3 | 16.26 k | 11.22 k | 87.95 k |
| EXPONENTIATION | 50000 / 97 | - | - | - | - | - | - | - | 5.09 k | 28.32 k | 8.46 k | n/a |
| DIVIDE | 50000 / 781 | - | - | - | - | - | - | - | 28.44 M | 686.67 M | 132.08 M | 1.59 G |
| ISQRT | 50000 / 195 | - | - | - | - | - | - | - | 3.67 M | 83.89 M | n/a | n/a |
| MODMUL_R2 | 50000 / 6250 | - | - | - | - | - | - | - | 745.19 k | 57.25 M | 22.31 M | 298.87 M |

## 5. Head to head

Best OpenCL GPU result against best OpenCL CPU result and the CPU libraries.
Ratios above 1.00x mean the GPU is faster than that baseline.

### secp256k1 (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GMP 192T | OpenSSL | CGBN | GPU vs CPU-CL | GPU vs GMP 192T | GPU vs OpenSSL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | w32-il | 3.93 G | w8 | 82.80 M | 65.93 M | 19.22 M | 23.01 M | 5.11 G | 47.43x | 204.29x | 170.68x | 0.77x |
| SUBTRACT | w32-il | 4.02 G | w8 | 65.86 M | 107.99 M | 28.58 M | 54.31 M | 6.35 G | 61.07x | 140.73x | 74.06x | 0.63x |
| ADDMOD | w32-il | 4.17 G | w8 | 39.67 M | 31.74 M | 35.90 M | 351.46 M | 5.07 G | 105.00x | 116.02x | 11.85x | 0.82x |
| SUBTRACTMOD | w32-il64 | 4.09 G | w8 | 48.89 M | 36.09 M | 1.87 G | 35.75 M | 4.61 G | 83.68x | 2.19x | 114.44x | 0.89x |
| MULTIPLYOPERANDSCANNING | w32-il64 | 3.76 G | w8 | 27.34 M | 76.07 M | 39.46 M | 43.10 M | n/a | 137.69x | 95.41x | 87.35x | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 2.47 G | w8 | 22.43 M | 74.59 M | 61.92 M | 60.95 M | 5.62 G | 110.24x | 39.94x | 40.57x | 0.44x |
| MONTGOMERYMULTIPLICATION | w32-il64 | 3.68 G | w8 | 16.23 M | 8.45 M | 13.99 M | 6.58 M | 4.52 G | 226.49x | 262.84x | 558.44x | 0.81x |
| COMPARE | w32-il | 3.92 G | w8 | 59.33 M | 169.64 M | 65.63 M | 5.47 G | 5.81 G | 66.01x | 59.67x | 0.72x | 0.67x |
| REDUCE | w32-il64 | 171.05 M | w8 | 2.43 M | 81.59 M | 58.87 M | 130.22 M | 3.90 G | 70.42x | 2.91x | 1.31x | 0.04x |
| MODMUL | w32-o64 | 36.09 M | w8 | 393.21 k | 16.08 M | 752.45 M | 151.09 M | 1.64 G | 91.78x | 0.05x | 0.24x | 0.02x |
| MODEXP | w32-o64 | 635.64 k | w8 | 391.7 | 132.71 k | 7.10 M | 4.43 M | 6.16 M | 1622.61x | 0.09x | 0.14x | 0.10x |
| EXPONENTIATION | w32-il64 | 2.95 M | w8 | 1.08 k | 405.63 k | 15.09 M | 1.50 M | n/a | 2738.08x | 0.20x | 1.97x | n/a |
| DIVIDE | w32-il64 | 59.59 M | w8 | 1.67 M | 43.06 M | 1.62 G | 304.30 M | 3.22 G | 35.72x | 0.04x | 0.20x | 0.02x |
| ISQRT | w32-il64 | 2.51 M | w8 | 13.96 k | 18.89 M | 853.55 M | n/a | n/a | 179.92x | 0.00x | n/a | n/a |
| MODMUL_R2 | w32-il64 | 3.16 G | w8 | 15.76 M | 16.07 M | 1.01 G | 27.25 M | 3.50 G | 200.46x | 3.13x | 115.95x | 0.90x |

### rsa256(composite) (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GMP 192T | OpenSSL | CGBN | GPU vs CPU-CL | GPU vs GMP 192T | GPU vs OpenSSL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | w32 | 3.86 G | w8 | 209.32 M | 72.86 M | 3.08 G | 2.96 G | 6.23 G | 18.43x | 1.25x | 1.30x | 0.62x |
| SUBTRACT | w32-il64 | 3.82 G | w8 | 76.24 M | 107.83 M | 4.53 G | 3.89 G | 5.94 G | 50.08x | 0.84x | 0.98x | 0.64x |
| ADDMOD | w32-il64 | 4.27 G | w8 | 59.86 M | 36.90 M | 1.90 G | 396.10 M | 4.57 G | 71.27x | 2.25x | 10.77x | 0.93x |
| SUBTRACTMOD | w32-il64 | 4.19 G | w8 | 44.90 M | 35.76 M | 1.86 G | 249.41 M | 4.40 G | 93.29x | 2.25x | 16.79x | 0.95x |
| MULTIPLYOPERANDSCANNING | w32-il64 | 3.59 G | w8 | 33.02 M | 75.57 M | 3.01 G | 2.09 G | n/a | 108.66x | 1.19x | 1.72x | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 2.51 G | w8 | 30.35 M | 75.20 M | 2.82 G | 2.10 G | 5.66 G | 82.58x | 0.89x | 1.19x | 0.44x |
| MONTGOMERYMULTIPLICATION | w32-il64 | 3.51 G | w8 | 11.18 M | 8.40 M | 544.70 M | 1.07 G | 4.60 G | 314.34x | 6.45x | 3.29x | 0.76x |
| COMPARE | w32-il64 | 4.11 G | w8 | 83.67 M | 169.17 M | 5.47 G | 5.45 G | 5.74 G | 49.11x | 0.75x | 0.75x | 0.72x |
| REDUCE | w32-il64 | 172.11 M | w8 | 1.79 M | 48.82 M | 2.30 G | 25.41 M | 3.93 G | 96.10x | 0.07x | 6.77x | 0.04x |
| MODMUL | w32-o64 | 35.82 M | w8 | 281.88 k | 16.21 M | 715.28 M | 156.02 M | 1.59 G | 127.07x | 0.05x | 0.23x | 0.02x |
| MODEXP | w32-o64 | 637.23 k | w8 | 394.6 | 143.59 k | 7.28 M | 4.41 M | 6.26 M | 1614.75x | 0.09x | 0.14x | 0.10x |
| EXPONENTIATION | w32-il64 | 2.95 M | w8 | 1.09 k | 408.84 k | 15.44 M | 979.95 k | n/a | 2715.69x | 0.19x | 3.01x | n/a |
| DIVIDE | w32-il64 | 55.33 M | w8 | 1.30 M | 42.54 M | 33.11 M | 308.38 M | 3.07 G | 42.51x | 1.67x | 0.18x | 0.02x |
| ISQRT | w32-il64 | 2.25 M | w8 | 14.36 k | 18.96 M | 7.92 M | n/a | n/a | 156.56x | 0.28x | n/a | n/a |
| MODMUL_R2 | w32-il64 | 3.00 G | w8 | 11.99 M | 16.28 M | 997.76 M | 158.35 M | 3.33 G | 250.36x | 3.01x | 18.96x | 0.90x |

### brainpoolP512r1 (512-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GMP 192T | OpenSSL | CGBN | GPU vs CPU-CL | GPU vs GMP 192T | GPU vs OpenSSL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | w32-il | 1.73 G | w8 | 31.32 M | 68.50 M | 2.96 G | 2.82 G | 5.70 G | 55.23x | 0.58x | 0.61x | 0.30x |
| SUBTRACT | w32-il | 1.59 G | w8 | 26.05 M | 101.32 M | 3.70 G | 3.28 G | 5.64 G | 61.18x | 0.43x | 0.49x | 0.28x |
| ADDMOD | w32-il | 1.83 G | w8 | 19.03 M | 34.20 M | 34.42 M | 399.36 M | 4.49 G | 96.40x | 53.29x | 4.59x | 0.41x |
| SUBTRACTMOD | w32-il | 1.79 G | w8 | 25.93 M | 32.47 M | 1.69 G | 245.72 M | 4.54 G | 68.94x | 1.06x | 7.28x | 0.39x |
| MULTIPLYOPERANDSCANNING | w32-il | 942.01 M | w8 | 6.72 M | 30.71 M | 1.53 G | 1.57 G | n/a | 140.13x | 0.62x | 0.60x | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il | 509.82 M | w8 | 7.02 M | 30.93 M | 1.39 G | 1.44 G | 4.93 G | 72.65x | 0.37x | 0.35x | 0.10x |
| MONTGOMERYMULTIPLICATION | w32-il64 | 1.14 G | w8 | 2.40 M | 3.47 M | 232.06 M | 513.26 M | 3.40 G | 474.61x | 4.90x | 2.22x | 0.33x |
| COMPARE | w32-il | 1.66 G | w8 | 30.20 M | 118.14 M | 5.33 G | 5.40 G | 5.39 G | 55.00x | 0.31x | 0.31x | 0.31x |
| REDUCE | w32-o64 | 34.02 M | w8 | 444.69 k | 45.95 M | 2.23 G | 256.07 M | 2.57 G | 76.51x | 0.02x | 0.13x | 0.01x |
| MODMUL | w32-o64 | 5.40 M | w8 | 16.42 k | 7.34 M | 372.28 M | 103.58 M | 591.41 M | 328.73x | 0.01x | 0.05x | 0.01x |
| MODEXP | w32-o64 | 45.71 k | w8 | 25.2 | 26.73 k | 1.33 M | 1.42 M | 2.05 M | 1813.14x | 0.03x | 0.03x | 0.02x |
| EXPONENTIATION | w32-o64 | 47.96 k | w8 | 57.1 | 115.50 k | 3.55 M | 365.31 k | n/a | 840.46x | 0.01x | 0.13x | n/a |
| DIVIDE | w32-il64 | 8.99 M | w8 | 38.33 k | 39.62 M | 775.41 M | 230.02 M | 1.96 G | 234.56x | 0.01x | 0.04x | 0.00x |
| ISQRT | w32-o64 | 250.43 k | w8 | 1.28 k | 11.12 M | 500.64 M | n/a | n/a | 196.31x | 0.00x | n/a | n/a |
| MODMUL_R2 | w32-il64 | 792.32 M | w8 | 1.95 M | 7.47 M | 496.72 M | 112.27 M | 2.34 G | 405.88x | 1.60x | 7.06x | 0.34x |

### p1024 (1024-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GMP 192T | OpenSSL | CGBN | GPU vs CPU-CL | GPU vs GMP 192T | GPU vs OpenSSL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | w32-il | 549.43 M | w8 | 10.53 M | 48.98 M | 43.39 M | 1.36 G | 4.57 G | 52.18x | 12.66x | 0.40x | 0.12x |
| SUBTRACT | w32-il64 | 448.06 M | w8 | 12.31 M | 69.03 M | 3.24 G | 209.09 M | 4.49 G | 36.40x | 0.14x | 2.14x | 0.10x |
| ADDMOD | w32-il64 | 743.17 M | w8 | 12.59 M | 23.45 M | 1.11 G | 116.44 M | 3.98 G | 59.04x | 0.67x | 6.38x | 0.19x |
| SUBTRACTMOD | w32-il64 | 767.48 M | w8 | 13.31 M | 28.18 M | 1.43 G | 94.05 M | 3.78 G | 57.66x | 0.54x | 8.16x | 0.20x |
| MULTIPLYOPERANDSCANNING | w32-il64 | 191.79 M | w8 | 998.17 k | 8.46 M | 463.51 M | 239.82 M | n/a | 192.15x | 0.41x | 0.80x | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 66.78 M | w8 | 1.10 M | 8.32 M | 451.97 M | 363.07 M | 2.51 G | 60.98x | 0.15x | 0.18x | 0.03x |
| MONTGOMERYMULTIPLICATION | w32-il64 | 229.99 M | w8 | 68.82 k | 1.12 M | 87.84 M | 171.23 M | 1.70 G | 3342.05x | 2.62x | 1.34x | 0.14x |
| COMPARE | w32-il | 709.22 M | w8 | 8.16 M | 123.98 M | 4.13 G | 3.17 G | 4.49 G | 86.97x | 0.17x | 0.22x | 0.16x |
| REDUCE | w32-opt | 5.53 M | w8 | 15.81 k | 71.51 M | 1.88 G | 102.64 M | 1.92 G | 349.63x | 0.00x | 0.05x | 0.00x |
| MODMUL | w32-il64 | 715.88 k | w8 | 1.89 k | 2.62 M | 125.40 M | 34.79 M | 235.10 M | 377.93x | 0.01x | 0.02x | 0.00x |
| MODEXP | w32-il64 | 2.42 k | w8 | 1.4 | 4.21 k | 156.30 k | 193.17 k | 367.87 k | 1719.35x | 0.02x | 0.01x | 0.01x |
| EXPONENTIATION | w32-o64 | 2.03 k | w8 | 3.6 | 30.56 k | 636.93 k | 37.04 k | n/a | 556.99x | 0.00x | 0.05x | n/a |
| DIVIDE | w32-il | 1.01 M | w8 | 7.96 k | 36.71 M | 1.03 G | 176.62 M | 1.47 G | 126.55x | 0.00x | 0.01x | 0.00x |
| ISQRT | w32-o64 | 17.10 k | w8 | 134.6 | 6.09 M | 167.34 M | n/a | n/a | 127.03x | 0.00x | n/a | n/a |
| MODMUL_R2 | w32-il64 | 124.82 M | w8 | 64.03 k | 2.55 M | 154.23 M | 53.06 M | 960.36 M | 1949.35x | 0.81x | 2.35x | 0.13x |

### p2048 (2048-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GMP 192T | OpenSSL | CGBN | GPU vs CPU-CL | GPU vs GMP 192T | GPU vs OpenSSL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | w32-il | 149.10 M | w8 | 7.36 M | 33.75 M | 2.11 G | 1.90 G | 3.88 G | 20.25x | 0.07x | 0.08x | 0.04x |
| SUBTRACT | w32-il | 149.85 M | w8 | 3.37 M | 38.93 M | 2.27 G | 1.94 G | 3.83 G | 44.43x | 0.07x | 0.08x | 0.04x |
| ADDMOD | w32-il64 | 246.61 M | w8 | 6.12 M | 16.91 M | 982.24 M | 12.00 M | 3.84 G | 40.27x | 0.25x | 20.54x | 0.06x |
| SUBTRACTMOD | w32-il64 | 235.99 M | w8 | 6.79 M | 19.24 M | 1.08 G | 120.39 M | 3.55 G | 34.76x | 0.22x | 1.96x | 0.07x |
| MULTIPLYOPERANDSCANNING | w32-il64 | 34.17 M | w8 | 52.45 k | 2.52 M | 142.22 M | 3.41 M | n/a | 651.48x | 0.24x | 10.03x | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 7.87 M | w8 | 60.70 k | 2.52 M | 145.99 M | 132.53 M | 840.05 M | 129.66x | 0.05x | 0.06x | 0.01x |
| MONTGOMERYMULTIPLICATION | w32 | 30.51 M | w8 | 8.74 k | 323.77 k | 29.33 M | 67.58 M | 564.08 M | 3490.08x | 1.04x | 0.45x | 0.05x |
| COMPARE | w32-il | 253.83 M | w8 | 6.15 M | 165.07 M | 4.48 G | 4.61 G | 3.80 G | 41.30x | 0.06x | 0.06x | 0.07x |
| REDUCE | w32-opt | 864.55 k | w8 | 2.00 k | 48.73 M | 1.39 G | 121.64 M | 1.90 G | 431.32x | 0.00x | 0.01x | 0.00x |
| MODMUL | w32-il64 | 108.76 k | w8 | 228.3 | 768.67 k | 1.15 M | 1.30 M | 121.83 M | 476.41x | 0.09x | 0.08x | 0.00x |
| MODEXP | w32-o64 | 99.7 | none | n/a | 590.3 | 16.26 k | 11.22 k | 87.95 k | n/a | 0.01x | 0.01x | 0.00x |
| EXPONENTIATION | w32-o64 | 119.7 | none | n/a | 5.09 k | 28.32 k | 8.46 k | n/a | n/a | 0.00x | 0.01x | n/a |
| DIVIDE | w32-o64 | 51.11 k | none | n/a | 28.44 M | 686.67 M | 132.08 M | 1.59 G | n/a | 0.00x | 0.00x | 0.00x |
| ISQRT | w32-o64 | 3.61 k | none | n/a | 3.67 M | 83.89 M | n/a | n/a | n/a | 0.00x | n/a | n/a |
| MODMUL_R2 | w32-il64 | 16.15 M | none | n/a | 745.19 k | 57.25 M | 22.31 M | 298.87 M | n/a | 0.28x | 0.72x | 0.05x |

## 6. Raw data

Also written to `NVIDIA_B200_Report.csv` for analysis.

```csv
kind,device,device_type,kernel,modulus,bits,operation,items,seconds,ops_per_sec,mismatches
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,secp256k1,256,ADD,50000,0.000758423,65926243.337,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,secp256k1,256,ADD,50000,0.002601033,19223129.159,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,secp256k1,256,ADD,50000,0.002173063,23009000.307,0
library,NVIDIA B200,gpu,cgbn,secp256k1,256,ADD,50000,0.000009792,5106209150.327,0
opencl-kernel,NVIDIA B200,GPU,w8,secp256k1,256,ADD,50000,0.000024400,2049174501.300,0
opencl-e2e,NVIDIA B200,GPU,w8,secp256k1,256,ADD,50000,0.000428652,116644726.138,0
opencl-kernel,NVIDIA B200,GPU,w16,secp256k1,256,ADD,50000,0.000017443,2866483772.708,0
opencl-e2e,NVIDIA B200,GPU,w16,secp256k1,256,ADD,50000,0.000544663,91799874.664,0
opencl-kernel,NVIDIA B200,GPU,w32,secp256k1,256,ADD,50000,0.000014584,3428403920.943,0
opencl-e2e,NVIDIA B200,GPU,w32,secp256k1,256,ADD,50000,0.000500768,99846644.058,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,secp256k1,256,ADD,50000,0.000013957,3582423301.360,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,secp256k1,256,ADD,50000,0.000629042,79485957.824,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,secp256k1,256,ADD,50000,0.000013699,3649917394.814,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,secp256k1,256,ADD,50000,0.000613162,81544525.758,0
opencl-kernel,NVIDIA B200,GPU,w32-il,secp256k1,256,ADD,50000,0.000012732,3927114482.431,0
opencl-e2e,NVIDIA B200,GPU,w32-il,secp256k1,256,ADD,50000,0.000517545,96609950.068,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,secp256k1,256,ADD,50000,0.000013733,3640851852.229,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,secp256k1,256,ADD,50000,0.000507157,98588787.649,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,secp256k1,256,ADD,50000,0.000603841,82803265.577,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,secp256k1,256,ADD,50000,0.001097838,45544058.364,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,secp256k1,256,SUBTRACT,50000,0.000462987,107994455.304,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,secp256k1,256,SUBTRACT,50000,0.001749407,28581117.679,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,secp256k1,256,SUBTRACT,50000,0.000920653,54309294.114,0
library,NVIDIA B200,gpu,cgbn,secp256k1,256,SUBTRACT,50000,0.000007872,6351626016.260,0
opencl-kernel,NVIDIA B200,GPU,w8,secp256k1,256,SUBTRACT,50000,0.000023170,2157960546.453,0
opencl-e2e,NVIDIA B200,GPU,w8,secp256k1,256,SUBTRACT,50000,0.000396282,126172793.002,0
opencl-kernel,NVIDIA B200,GPU,w16,secp256k1,256,SUBTRACT,50000,0.000017471,2861880590.371,0
opencl-e2e,NVIDIA B200,GPU,w16,secp256k1,256,SUBTRACT,50000,0.000503726,99260321.153,0
opencl-kernel,NVIDIA B200,GPU,w32,secp256k1,256,SUBTRACT,50000,0.000014012,3568374815.972,0
opencl-e2e,NVIDIA B200,GPU,w32,secp256k1,256,SUBTRACT,50000,0.000545339,91686075.928,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,secp256k1,256,SUBTRACT,50000,0.000014114,3542591675.877,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,secp256k1,256,SUBTRACT,50000,0.000559671,89338205.738,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,secp256k1,256,SUBTRACT,50000,0.000013753,3635551047.081,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,secp256k1,256,SUBTRACT,50000,0.000505597,98892995.490,0
opencl-kernel,NVIDIA B200,GPU,w32-il,secp256k1,256,SUBTRACT,50000,0.000012431,4022182854.788,0
opencl-e2e,NVIDIA B200,GPU,w32-il,secp256k1,256,SUBTRACT,50000,0.000521681,95844004.850,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,secp256k1,256,SUBTRACT,50000,0.000012665,3947906808.468,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,secp256k1,256,SUBTRACT,50000,0.000616564,81094588.844,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,secp256k1,256,SUBTRACT,50000,0.000759195,65859236.441,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,secp256k1,256,SUBTRACT,50000,0.001398876,35742983.241,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,secp256k1,256,ADDMOD,50000,0.001575398,31738011.858,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,secp256k1,256,ADDMOD,50000,0.001392565,35904957.525,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,secp256k1,256,ADDMOD,50000,0.000142262,351464950.017,0
library,NVIDIA B200,gpu,cgbn,secp256k1,256,ADDMOD,50000,0.000009856,5073051948.052,0
opencl-kernel,NVIDIA B200,GPU,w8,secp256k1,256,ADDMOD,50000,0.000030585,1634789225.153,0
opencl-e2e,NVIDIA B200,GPU,w8,secp256k1,256,ADDMOD,50000,0.000400778,124757339.927,0
opencl-kernel,NVIDIA B200,GPU,w16,secp256k1,256,ADDMOD,50000,0.000020803,2403505000.672,0
opencl-e2e,NVIDIA B200,GPU,w16,secp256k1,256,ADDMOD,50000,0.000583740,85654579.394,0
opencl-kernel,NVIDIA B200,GPU,w32,secp256k1,256,ADDMOD,50000,0.000016099,3105790985.545,0
opencl-e2e,NVIDIA B200,GPU,w32,secp256k1,256,ADDMOD,50000,0.000497927,100416334.577,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,secp256k1,256,ADDMOD,50000,0.000014055,3557468501.048,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,secp256k1,256,ADDMOD,50000,0.000492767,101467816.880,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,secp256k1,256,ADDMOD,50000,0.000014397,3472953849.388,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,secp256k1,256,ADDMOD,50000,0.000688815,72588428.353,0
opencl-kernel,NVIDIA B200,GPU,w32-il,secp256k1,256,ADDMOD,50000,0.000012003,4165624650.599,0
opencl-e2e,NVIDIA B200,GPU,w32-il,secp256k1,256,ADDMOD,50000,0.000569151,87850146.268,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,secp256k1,256,ADDMOD,50000,0.000012065,4144201254.366,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,secp256k1,256,ADDMOD,50000,0.000698088,71624210.935,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,secp256k1,256,ADDMOD,50000,0.001260337,39671930.459,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,secp256k1,256,ADDMOD,50000,0.002378526,21021422.270,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,secp256k1,256,SUBTRACTMOD,50000,0.001385536,36087117.031,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,secp256k1,256,SUBTRACTMOD,50000,0.000026766,1868067457.974,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,secp256k1,256,SUBTRACTMOD,50000,0.001398660,35748492.963,0
library,NVIDIA B200,gpu,cgbn,secp256k1,256,SUBTRACTMOD,50000,0.000010848,4609144542.773,0
opencl-kernel,NVIDIA B200,GPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.000030373,1646199452.668,0
opencl-e2e,NVIDIA B200,GPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.000548831,91102713.863,0
opencl-kernel,NVIDIA B200,GPU,w16,secp256k1,256,SUBTRACTMOD,50000,0.000021538,2321478458.462,0
opencl-e2e,NVIDIA B200,GPU,w16,secp256k1,256,SUBTRACTMOD,50000,0.000605449,82583328.001,0
opencl-kernel,NVIDIA B200,GPU,w32,secp256k1,256,SUBTRACTMOD,50000,0.000014833,3370849033.473,0
opencl-e2e,NVIDIA B200,GPU,w32,secp256k1,256,SUBTRACTMOD,50000,0.000524587,95313075.234,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,secp256k1,256,SUBTRACTMOD,50000,0.000014281,3501151268.820,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,secp256k1,256,SUBTRACTMOD,50000,0.000610723,81870170.953,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,secp256k1,256,SUBTRACTMOD,50000,0.000014302,3496021502.120,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,secp256k1,256,SUBTRACTMOD,50000,0.000550686,90795838.369,0
opencl-kernel,NVIDIA B200,GPU,w32-il,secp256k1,256,SUBTRACTMOD,50000,0.000012301,4064701931.576,0
opencl-e2e,NVIDIA B200,GPU,w32-il,secp256k1,256,SUBTRACTMOD,50000,0.000566655,88237117.766,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,secp256k1,256,SUBTRACTMOD,50000,0.000012222,4090990509.211,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,secp256k1,256,SUBTRACTMOD,50000,0.000672129,74390483.488,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.001022742,48888184.256,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.001248861,40036478.646,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000657290,76069884.509,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001267165,39458170.692,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001160176,43096909.241,0
opencl-kernel,NVIDIA B200,GPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001676169,29829927.806,0
opencl-e2e,NVIDIA B200,GPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.002307451,21668932.694,0
opencl-kernel,NVIDIA B200,GPU,w16,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000593799,84203570.219,0
opencl-e2e,NVIDIA B200,GPU,w16,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001222211,40909468.666,0
opencl-kernel,NVIDIA B200,GPU,w32,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000078091,640279205.004,0
opencl-e2e,NVIDIA B200,GPU,w32,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000797679,62681859.717,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000018211,2745581016.672,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000653678,76490263.330,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000016558,3019691276.225,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000594080,84163754.637,0
opencl-kernel,NVIDIA B200,GPU,w32-il,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000014736,3393058433.730,0
opencl-e2e,NVIDIA B200,GPU,w32-il,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000725740,68895198.072,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000013282,3764510168.199,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000620810,80539944.197,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001828751,27341064.108,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.002597143,19251924.562,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000670352,74587680.392,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000807535,61916857.368,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000820349,60949669.005,0
library,NVIDIA B200,gpu,cgbn,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000008896,5620503597.122,0
opencl-kernel,NVIDIA B200,GPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000172680,289552911.799,0
opencl-e2e,NVIDIA B200,GPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000785020,63692649.056,0
opencl-kernel,NVIDIA B200,GPU,w16,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000059767,836581442.359,0
opencl-e2e,NVIDIA B200,GPU,w16,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000667805,74872159.321,0
opencl-kernel,NVIDIA B200,GPU,w32,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000024301,2057528501.897,0
opencl-e2e,NVIDIA B200,GPU,w32,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000613710,81471701.144,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000023326,2143528837.295,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000613843,81454055.978,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000023355,2140868364.753,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000672213,74381181.885,0
opencl-kernel,NVIDIA B200,GPU,w32-il,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000021252,2352710593.030,0
opencl-e2e,NVIDIA B200,GPU,w32-il,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000686535,72829495.713,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000020220,2472806006.172,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000637345,78450459.135,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.002229005,22431533.534,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.003605834,13866417.683,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.005916322,8451196.538,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.003575098,13985630.595,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.007595871,6582523.590,0
library,NVIDIA B200,gpu,cgbn,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000011072,4515895953.757,0
opencl-kernel,NVIDIA B200,GPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000076332,655034162.084,0
opencl-e2e,NVIDIA B200,GPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000447778,111662467.814,0
opencl-kernel,NVIDIA B200,GPU,w16,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000028990,1724734078.114,0
opencl-e2e,NVIDIA B200,GPU,w16,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000579302,86310754.764,0
opencl-kernel,NVIDIA B200,GPU,w32,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000014849,3367228499.749,0
opencl-e2e,NVIDIA B200,GPU,w32,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000550870,90765521.519,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000016889,2960515110.115,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000628846,79510722.940,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000015227,3283639244.948,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000501701,99660948.227,0
opencl-kernel,NVIDIA B200,GPU,w32-il,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000016506,3029211338.294,0
opencl-e2e,NVIDIA B200,GPU,w32-il,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000558723,89489783.418,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000013602,3675939144.129,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000691056,72353034.555,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.003080687,16230145.963,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.004541264,11010150.400,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,secp256k1,256,COMPARE,50000,0.000294738,169642307.399,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,secp256k1,256,COMPARE,50000,0.000761795,65634489.792,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,secp256k1,256,COMPARE,50000,0.000009133,5474643471.167,0
library,NVIDIA B200,gpu,cgbn,secp256k1,256,COMPARE,50000,0.000008608,5808550185.874,0
opencl-kernel,NVIDIA B200,GPU,w8,secp256k1,256,COMPARE,50000,0.000023584,2120081593.405,0
opencl-e2e,NVIDIA B200,GPU,w8,secp256k1,256,COMPARE,50000,0.000387270,129108889.252,0
opencl-kernel,NVIDIA B200,GPU,w16,secp256k1,256,COMPARE,50000,0.000018026,2773774102.634,0
opencl-e2e,NVIDIA B200,GPU,w16,secp256k1,256,COMPARE,50000,0.000593089,84304374.596,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,secp256k1,256,COMPARE,50000,0.000013594,3678111257.076,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,secp256k1,256,COMPARE,50000,0.000644835,77539222.074,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,secp256k1,256,COMPARE,50000,0.000013566,3685654838.156,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,secp256k1,256,COMPARE,50000,0.000705381,70883674.400,0
opencl-kernel,NVIDIA B200,GPU,w32-il,secp256k1,256,COMPARE,50000,0.000012767,3916335937.557,0
opencl-e2e,NVIDIA B200,GPU,w32-il,secp256k1,256,COMPARE,50000,0.000513782,97317534.374,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,secp256k1,256,COMPARE,50000,0.000013084,3821485270.932,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,secp256k1,256,COMPARE,50000,0.000520025,96149217.917,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,secp256k1,256,COMPARE,50000,0.000842747,59329783.778,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,secp256k1,256,COMPARE,50000,0.001843834,27117409.402,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,secp256k1,256,REDUCE,6250,0.000076605,81587763.072,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,secp256k1,256,REDUCE,6250,0.000106167,58869534.537,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,secp256k1,256,REDUCE,6250,0.000047997,130217319.687,0
library,NVIDIA B200,gpu,cgbn,secp256k1,256,REDUCE,50000,0.000012832,3896508728.180,0
opencl-kernel,NVIDIA B200,GPU,w8,secp256k1,256,REDUCE,50000,0.000141956,352222033.276,0
opencl-e2e,NVIDIA B200,GPU,w8,secp256k1,256,REDUCE,50000,0.000574977,86959997.003,0
opencl-kernel,NVIDIA B200,GPU,w16,secp256k1,256,REDUCE,50000,0.000096272,519361923.165,0
opencl-e2e,NVIDIA B200,GPU,w16,secp256k1,256,REDUCE,50000,0.000746374,66990541.955,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,secp256k1,256,REDUCE,50000,0.000037959,1317211996.381,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,secp256k1,256,REDUCE,50000,0.000659586,75805134.244,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,secp256k1,256,REDUCE,50000,0.000037473,1334292445.579,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,secp256k1,256,REDUCE,50000,0.000550502,90826194.685,0
opencl-kernel,NVIDIA B200,GPU,w32-il,secp256k1,256,REDUCE,50000,0.000037458,1334827386.702,0
opencl-e2e,NVIDIA B200,GPU,w32-il,secp256k1,256,REDUCE,50000,0.000558609,89508060.288,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,secp256k1,256,REDUCE,50000,0.000036538,1368438469.504,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,secp256k1,256,REDUCE,50000,0.000520263,96105242.105,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,secp256k1,256,REDUCE,50000,0.002572986,19432675.113,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,secp256k1,256,REDUCE,50000,0.003819028,13092336.474,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,secp256k1,256,MODMUL,3125,0.000194320,16081756.664,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,secp256k1,256,MODMUL,3125,0.000004153,752445650.235,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,secp256k1,256,MODMUL,3125,0.000020683,151088400.438,0
library,NVIDIA B200,gpu,cgbn,secp256k1,256,MODMUL,50000,0.000030400,1644736842.105,0
opencl-kernel,NVIDIA B200,GPU,w8,secp256k1,256,MODMUL,50000,0.000378534,132088502.615,0
opencl-e2e,NVIDIA B200,GPU,w8,secp256k1,256,MODMUL,50000,0.000886761,56384978.507,0
opencl-kernel,NVIDIA B200,GPU,w16,secp256k1,256,MODMUL,50000,0.000226118,221123464.853,0
opencl-e2e,NVIDIA B200,GPU,w16,secp256k1,256,MODMUL,50000,0.000824291,60658187.957,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,secp256k1,256,MODMUL,50000,0.000109264,457607074.585,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,secp256k1,256,MODMUL,50000,0.000736819,67859270.903,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,secp256k1,256,MODMUL,50000,0.000086588,577447420.424,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,secp256k1,256,MODMUL,50000,0.000584439,85552141.111,0
opencl-kernel,NVIDIA B200,GPU,w32-il,secp256k1,256,MODMUL,50000,0.000107336,465827046.808,0
opencl-e2e,NVIDIA B200,GPU,w32-il,secp256k1,256,MODMUL,50000,0.000645573,77450572.742,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,secp256k1,256,MODMUL,50000,0.000092607,539916038.457,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,secp256k1,256,MODMUL,50000,0.000604971,82648594.857,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,secp256k1,256,MODMUL,50000,0.007947397,6291368.142,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,secp256k1,256,MODMUL,50000,0.008096659,6175386.682,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,secp256k1,256,MODEXP,781,0.005884845,132713.777,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,secp256k1,256,MODEXP,781,0.000110057,7096326.053,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,secp256k1,256,MODEXP,781,0.000176147,4433797.383,0
library,NVIDIA B200,gpu,cgbn,secp256k1,256,MODEXP,50000,0.008116383,6160379.568,0
opencl-kernel,NVIDIA B200,GPU,w8,secp256k1,256,MODEXP,50000,0.017207982,2905628.327,0
opencl-e2e,NVIDIA B200,GPU,w8,secp256k1,256,MODEXP,50000,0.017736248,2819085.529,0
opencl-kernel,NVIDIA B200,GPU,w16,secp256k1,256,MODEXP,50000,0.003220264,15526677.388,0
opencl-e2e,NVIDIA B200,GPU,w16,secp256k1,256,MODEXP,50000,0.003953284,12647712.630,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,secp256k1,256,MODEXP,50000,0.002549374,19612657.772,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,secp256k1,256,MODEXP,50000,0.003249582,15386594.244,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,secp256k1,256,MODEXP,50000,0.001228675,40694243.355,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,secp256k1,256,MODEXP,50000,0.001938267,25796239.246,0
opencl-kernel,NVIDIA B200,GPU,w32-il,secp256k1,256,MODEXP,50000,0.002574028,19424808.261,0
opencl-e2e,NVIDIA B200,GPU,w32-il,secp256k1,256,MODEXP,50000,0.003347053,14938514.228,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,secp256k1,256,MODEXP,50000,0.001238119,40383841.188,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,secp256k1,256,MODEXP,50000,0.001870909,26724976.899,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,secp256k1,256,MODEXP,50000,1.993655698,25079.556,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,secp256k1,256,MODEXP,50000,1.975283113,25312.827,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,secp256k1,256,EXPONENTIATION,781,0.001925397,405630.528,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,secp256k1,256,EXPONENTIATION,781,0.000051751,15091502.778,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,secp256k1,256,EXPONENTIATION,781,0.000522033,1496073.945,0
opencl-kernel,NVIDIA B200,GPU,w8,secp256k1,256,EXPONENTIATION,50000,0.023444199,2132723.749,0
opencl-e2e,NVIDIA B200,GPU,w8,secp256k1,256,EXPONENTIATION,50000,0.023964585,2086412.092,0
opencl-kernel,NVIDIA B200,GPU,w16,secp256k1,256,EXPONENTIATION,50000,0.004697283,10644451.200,0
opencl-e2e,NVIDIA B200,GPU,w16,secp256k1,256,EXPONENTIATION,50000,0.005465483,9148322.198,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,secp256k1,256,EXPONENTIATION,50000,0.000471590,106024281.251,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,secp256k1,256,EXPONENTIATION,50000,0.001168216,42800304.178,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,secp256k1,256,EXPONENTIATION,50000,0.000267065,187220329.972,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,secp256k1,256,EXPONENTIATION,50000,0.000908813,55016820.237,0
opencl-kernel,NVIDIA B200,GPU,w32-il,secp256k1,256,EXPONENTIATION,50000,0.000473693,105553594.011,0
opencl-e2e,NVIDIA B200,GPU,w32-il,secp256k1,256,EXPONENTIATION,50000,0.001178558,42424726.366,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,secp256k1,256,EXPONENTIATION,50000,0.000265061,188635794.694,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,secp256k1,256,EXPONENTIATION,50000,0.000775578,64468051.162,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,secp256k1,256,EXPONENTIATION,50000,0.725758990,68893.394,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,secp256k1,256,EXPONENTIATION,50000,0.701901275,71235.089,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,secp256k1,256,DIVIDE,6250,0.000145151,43058516.942,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,secp256k1,256,DIVIDE,6250,0.000003859,1619783200.891,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,secp256k1,256,DIVIDE,6250,0.000020539,304299226.534,0
library,NVIDIA B200,gpu,cgbn,secp256k1,256,DIVIDE,50000,0.000015520,3221649484.536,0
opencl-kernel,NVIDIA B200,GPU,w8,secp256k1,256,DIVIDE,50000,0.000342900,145815085.237,0
opencl-e2e,NVIDIA B200,GPU,w8,secp256k1,256,DIVIDE,50000,0.000866835,57681103.876,0
opencl-kernel,NVIDIA B200,GPU,w16,secp256k1,256,DIVIDE,50000,0.000291891,171296860.557,0
opencl-e2e,NVIDIA B200,GPU,w16,secp256k1,256,DIVIDE,50000,0.000891257,56100543.598,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,secp256k1,256,DIVIDE,50000,0.000111268,449365421.721,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,secp256k1,256,DIVIDE,50000,0.000694248,72020372.044,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,secp256k1,256,DIVIDE,50000,0.000109373,457151175.670,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,secp256k1,256,DIVIDE,50000,0.000709817,70440698.009,0
opencl-kernel,NVIDIA B200,GPU,w32-il,secp256k1,256,DIVIDE,50000,0.000108915,459073455.091,0
opencl-e2e,NVIDIA B200,GPU,w32-il,secp256k1,256,DIVIDE,50000,0.000842797,59326268.034,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,secp256k1,256,DIVIDE,50000,0.000104885,476712743.256,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,secp256k1,256,DIVIDE,50000,0.000730421,68453679.314,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,secp256k1,256,DIVIDE,50000,0.003746361,13346284.546,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,secp256k1,256,DIVIDE,50000,0.005469974,9140811.249,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,secp256k1,256,ISQRT,1562,0.000082669,18894706.679,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,secp256k1,256,ISQRT,1562,0.000001830,853551592.014,0
opencl-kernel,NVIDIA B200,GPU,w8,secp256k1,256,ISQRT,50000,0.003715724,13456327.809,0
opencl-e2e,NVIDIA B200,GPU,w8,secp256k1,256,ISQRT,50000,0.004255414,11749738.073,0
opencl-kernel,NVIDIA B200,GPU,w16,secp256k1,256,ISQRT,50000,0.002775656,18013759.245,0
opencl-e2e,NVIDIA B200,GPU,w16,secp256k1,256,ISQRT,50000,0.003269195,15294285.104,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,secp256k1,256,ISQRT,50000,0.000682747,73233570.593,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,secp256k1,256,ISQRT,50000,0.001399096,35737362.220,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,secp256k1,256,ISQRT,50000,0.000628953,79497198.068,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,secp256k1,256,ISQRT,50000,0.001249602,40012737.991,0
opencl-kernel,NVIDIA B200,GPU,w32-il,secp256k1,256,ISQRT,50000,0.000699565,71472981.942,0
opencl-e2e,NVIDIA B200,GPU,w32-il,secp256k1,256,ISQRT,50000,0.001386223,36069235.454,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,secp256k1,256,ISQRT,50000,0.000621756,80417389.012,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,secp256k1,256,ISQRT,50000,0.001271049,39337585.696,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,secp256k1,256,ISQRT,50000,0.111865059,446967.091,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,secp256k1,256,ISQRT,50000,0.113345486,441129.169,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,secp256k1,256,MODMUL_R2,50000,0.003111946,16067114.498,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,secp256k1,256,MODMUL_R2,50000,0.000049573,1008614335.583,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,secp256k1,256,MODMUL_R2,50000,0.001834716,27252174.066,0
library,NVIDIA B200,gpu,cgbn,secp256k1,256,MODMUL_R2,50000,0.000014272,3503363228.700,0
opencl-kernel,NVIDIA B200,GPU,w8,secp256k1,256,MODMUL_R2,50000,0.000090488,552559456.365,0
opencl-e2e,NVIDIA B200,GPU,w8,secp256k1,256,MODMUL_R2,50000,0.000474306,105417165.432,0
opencl-kernel,NVIDIA B200,GPU,w16,secp256k1,256,MODMUL_R2,50000,0.000029882,1673244364.104,0
opencl-e2e,NVIDIA B200,GPU,w16,secp256k1,256,MODMUL_R2,50000,0.000528655,94579638.115,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,secp256k1,256,MODMUL_R2,50000,0.000021616,2313101732.012,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,secp256k1,256,MODMUL_R2,50000,0.000504765,99055981.147,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,secp256k1,256,MODMUL_R2,50000,0.000017684,2827422119.233,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,secp256k1,256,MODMUL_R2,50000,0.000490122,102015414.682,0
opencl-kernel,NVIDIA B200,GPU,w32-il,secp256k1,256,MODMUL_R2,50000,0.000021538,2321478458.462,0
opencl-e2e,NVIDIA B200,GPU,w32-il,secp256k1,256,MODMUL_R2,50000,0.000527922,94710970.481,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,secp256k1,256,MODMUL_R2,50000,0.000015823,3159946215.025,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,secp256k1,256,MODMUL_R2,50000,0.000632220,79086400.497,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,secp256k1,256,MODMUL_R2,50000,0.003171959,15763129.543,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,secp256k1,256,MODMUL_R2,50000,0.004481144,11157865.004,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,rsa256(composite),256,ADD,50000,0.000686255,72859179.905,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,rsa256(composite),256,ADD,50000,0.000016235,3079830455.775,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,rsa256(composite),256,ADD,50000,0.000016875,2962904631.686,0
library,NVIDIA B200,gpu,cgbn,rsa256(composite),256,ADD,50000,0.000008032,6225099601.594,0
opencl-kernel,NVIDIA B200,GPU,w8,rsa256(composite),256,ADD,50000,0.000023520,2125853092.781,0
opencl-e2e,NVIDIA B200,GPU,w8,rsa256(composite),256,ADD,50000,0.000393144,127179834.425,0
opencl-kernel,NVIDIA B200,GPU,w16,rsa256(composite),256,ADD,50000,0.000018319,2729406831.513,0
opencl-e2e,NVIDIA B200,GPU,w16,rsa256(composite),256,ADD,50000,0.000635720,78650967.167,0
opencl-kernel,NVIDIA B200,GPU,w32,rsa256(composite),256,ADD,50000,0.000012963,3857142994.674,0
opencl-e2e,NVIDIA B200,GPU,w32,rsa256(composite),256,ADD,50000,0.000504770,99055021.644,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,rsa256(composite),256,ADD,50000,0.000013348,3745894133.859,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,rsa256(composite),256,ADD,50000,0.000622302,80346848.406,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,rsa256(composite),256,ADD,50000,0.000013447,3718296666.061,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,rsa256(composite),256,ADD,50000,0.000519202,96301636.871,0
opencl-kernel,NVIDIA B200,GPU,w32-il,rsa256(composite),256,ADD,50000,0.000013143,3804289975.022,0
opencl-e2e,NVIDIA B200,GPU,w32-il,rsa256(composite),256,ADD,50000,0.000525178,95205808.951,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,rsa256(composite),256,ADD,50000,0.000013750,3636351341.100,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,rsa256(composite),256,ADD,50000,0.000566393,88277923.959,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,rsa256(composite),256,ADD,50000,0.000238867,209321495.702,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,rsa256(composite),256,ADD,50000,0.001132810,44138027.559,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACT,50000,0.000463702,107827873.569,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACT,50000,0.000011026,4534813586.595,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACT,50000,0.000012838,3894623380.251,0
library,NVIDIA B200,gpu,cgbn,rsa256(composite),256,SUBTRACT,50000,0.000008416,5941064638.783,0
opencl-kernel,NVIDIA B200,GPU,w8,rsa256(composite),256,SUBTRACT,50000,0.000023287,2147118637.832,0
opencl-e2e,NVIDIA B200,GPU,w8,rsa256(composite),256,SUBTRACT,50000,0.000392343,127439499.781,0
opencl-kernel,NVIDIA B200,GPU,w16,rsa256(composite),256,SUBTRACT,50000,0.000017520,2853874718.265,0
opencl-e2e,NVIDIA B200,GPU,w16,rsa256(composite),256,SUBTRACT,50000,0.000563610,88713827.586,0
opencl-kernel,NVIDIA B200,GPU,w32,rsa256(composite),256,SUBTRACT,50000,0.000013844,3611674581.858,0
opencl-e2e,NVIDIA B200,GPU,w32,rsa256(composite),256,SUBTRACT,50000,0.000691427,72314210.347,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,rsa256(composite),256,SUBTRACT,50000,0.000013682,3654420475.121,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,rsa256(composite),256,SUBTRACT,50000,0.000677224,73830813.486,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,rsa256(composite),256,SUBTRACT,50000,0.000013945,3585503682.370,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,rsa256(composite),256,SUBTRACT,50000,0.000531613,94053378.841,0
opencl-kernel,NVIDIA B200,GPU,w32-il,rsa256(composite),256,SUBTRACT,50000,0.000013268,3768440753.869,0
opencl-e2e,NVIDIA B200,GPU,w32-il,rsa256(composite),256,SUBTRACT,50000,0.000542461,92172515.127,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,rsa256(composite),256,SUBTRACT,50000,0.000013096,3817952331.680,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,rsa256(composite),256,SUBTRACT,50000,0.000371334,134649658.859,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,rsa256(composite),256,SUBTRACT,50000,0.000655854,76236478.945,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,rsa256(composite),256,SUBTRACT,50000,0.000912874,54772071.158,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,rsa256(composite),256,ADDMOD,50000,0.001354838,36904772.154,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,rsa256(composite),256,ADDMOD,50000,0.000026350,1897508977.782,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,rsa256(composite),256,ADDMOD,50000,0.000126229,396104463.180,0
library,NVIDIA B200,gpu,cgbn,rsa256(composite),256,ADDMOD,50000,0.000010944,4568713450.292,0
opencl-kernel,NVIDIA B200,GPU,w8,rsa256(composite),256,ADDMOD,50000,0.000028846,1733344349.335,0
opencl-e2e,NVIDIA B200,GPU,w8,rsa256(composite),256,ADDMOD,50000,0.000443546,112727883.211,0
opencl-kernel,NVIDIA B200,GPU,w16,rsa256(composite),256,ADDMOD,50000,0.000020526,2435934876.387,0
opencl-e2e,NVIDIA B200,GPU,w16,rsa256(composite),256,ADDMOD,50000,0.000500507,99898688.607,0
opencl-kernel,NVIDIA B200,GPU,w32,rsa256(composite),256,ADDMOD,50000,0.000014407,3470540419.377,0
opencl-e2e,NVIDIA B200,GPU,w32,rsa256(composite),256,ADDMOD,50000,0.000678915,73646927.528,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,rsa256(composite),256,ADDMOD,50000,0.000014215,3517408887.360,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,rsa256(composite),256,ADDMOD,50000,0.000544489,91829217.650,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,rsa256(composite),256,ADDMOD,50000,0.000014032,3563312367.567,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,rsa256(composite),256,ADDMOD,50000,0.000538670,92821209.127,0
opencl-kernel,NVIDIA B200,GPU,w32-il,rsa256(composite),256,ADDMOD,50000,0.000012026,4157680679.948,0
opencl-e2e,NVIDIA B200,GPU,w32-il,rsa256(composite),256,ADDMOD,50000,0.000703700,71053005.218,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,rsa256(composite),256,ADDMOD,50000,0.000011720,4266213020.244,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,rsa256(composite),256,ADDMOD,50000,0.000733510,68165395.681,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,rsa256(composite),256,ADDMOD,50000,0.000835330,59856584.149,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,rsa256(composite),256,ADDMOD,50000,0.001317313,37956053.148,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACTMOD,50000,0.001398196,35760357.123,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACTMOD,50000,0.000026835,1863215903.034,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACTMOD,50000,0.000200474,249408497.413,0
library,NVIDIA B200,gpu,cgbn,rsa256(composite),256,SUBTRACTMOD,50000,0.000011360,4401408450.704,0
opencl-kernel,NVIDIA B200,GPU,w8,rsa256(composite),256,SUBTRACTMOD,50000,0.000030879,1619227023.766,0
opencl-e2e,NVIDIA B200,GPU,w8,rsa256(composite),256,SUBTRACTMOD,50000,0.000450014,111107649.412,0
opencl-kernel,NVIDIA B200,GPU,w16,rsa256(composite),256,SUBTRACTMOD,50000,0.000021190,2359599879.135,0
opencl-e2e,NVIDIA B200,GPU,w16,rsa256(composite),256,SUBTRACTMOD,50000,0.000520616,96040062.450,0
opencl-kernel,NVIDIA B200,GPU,w32,rsa256(composite),256,SUBTRACTMOD,50000,0.000014762,3387064623.635,0
opencl-e2e,NVIDIA B200,GPU,w32,rsa256(composite),256,SUBTRACTMOD,50000,0.000487873,102485691.767,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,rsa256(composite),256,SUBTRACTMOD,50000,0.000013949,3584486272.022,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,rsa256(composite),256,SUBTRACTMOD,50000,0.000510921,97862486.915,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,rsa256(composite),256,SUBTRACTMOD,50000,0.000014557,3434764799.591,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,rsa256(composite),256,SUBTRACTMOD,50000,0.000529896,94358137.798,0
opencl-kernel,NVIDIA B200,GPU,w32-il,rsa256(composite),256,SUBTRACTMOD,50000,0.000012538,3987899067.781,0
opencl-e2e,NVIDIA B200,GPU,w32-il,rsa256(composite),256,SUBTRACTMOD,50000,0.000668497,74794657.862,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,rsa256(composite),256,SUBTRACTMOD,50000,0.000011937,4188659127.348,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,rsa256(composite),256,SUBTRACTMOD,50000,0.000371292,134664899.632,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,rsa256(composite),256,SUBTRACTMOD,50000,0.001113625,44898416.908,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,rsa256(composite),256,SUBTRACTMOD,50000,0.001638120,30522793.672,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000661603,75573992.512,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000016634,3005980368.944,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000023976,2085375939.094,0
opencl-kernel,NVIDIA B200,GPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001720712,29057739.935,0
opencl-e2e,NVIDIA B200,GPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.002283989,21891524.129,0
opencl-kernel,NVIDIA B200,GPU,w16,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000594454,84110797.108,0
opencl-e2e,NVIDIA B200,GPU,w16,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001422433,35151041.747,0
opencl-kernel,NVIDIA B200,GPU,w32,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000080847,618451641.033,0
opencl-e2e,NVIDIA B200,GPU,w32,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000698148,71618048.219,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000017707,2823741499.783,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000782298,63914259.536,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000016702,2993655281.629,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000643648,77682214.490,0
opencl-kernel,NVIDIA B200,GPU,w32-il,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000014455,3458996920.302,0
opencl-e2e,NVIDIA B200,GPU,w32-il,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000759092,65868165.068,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000013937,3587570203.311,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000655264,76305121.472,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001514376,33016898.831,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.003100184,16128074.847,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000664859,75203946.005,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000017706,2823862176.067,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000023827,2098417403.073,0
library,NVIDIA B200,gpu,cgbn,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000008832,5661231884.058,0
opencl-kernel,NVIDIA B200,GPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000173205,288675195.655,0
opencl-e2e,NVIDIA B200,GPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000685609,72927857.516,0
opencl-kernel,NVIDIA B200,GPU,w16,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000060432,827376108.586,0
opencl-e2e,NVIDIA B200,GPU,w16,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000687550,72721978.053,0
opencl-kernel,NVIDIA B200,GPU,w32,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000024341,2054153197.473,0
opencl-e2e,NVIDIA B200,GPU,w32,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000636916,78503284.570,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000023950,2087681997.190,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000940604,53157331.875,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000023949,2087773330.741,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000646305,77362864.426,0
opencl-kernel,NVIDIA B200,GPU,w32-il,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000021324,2344785635.281,0
opencl-e2e,NVIDIA B200,GPU,w32-il,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000661698,75563166.203,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000019950,2506268517.643,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000645857,77416523.281,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001647467,30349623.227,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.002879336,17365114.670,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.005952723,8399517.365,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000091794,544698225.500,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000046865,1066895025.176,0
library,NVIDIA B200,gpu,cgbn,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000010880,4595588235.294,0
opencl-kernel,NVIDIA B200,GPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000075587,661488784.022,0
opencl-e2e,NVIDIA B200,GPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000474569,105358774.404,0
opencl-kernel,NVIDIA B200,GPU,w16,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000028723,1740763062.307,0
opencl-e2e,NVIDIA B200,GPU,w16,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000518706,96393709.831,0
opencl-kernel,NVIDIA B200,GPU,w32,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000015620,3201019039.314,0
opencl-e2e,NVIDIA B200,GPU,w32,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000628650,79535503.493,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000017302,2889840264.293,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000526566,94954847.930,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000015537,3218144099.025,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000558504,89524870.407,0
opencl-kernel,NVIDIA B200,GPU,w32-il,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000016835,2969993704.534,0
opencl-e2e,NVIDIA B200,GPU,w32-il,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000536841,93137428.000,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000014233,3512978321.610,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000649152,77023566.091,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.004473995,11175694.327,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.005585386,8951932.853,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,rsa256(composite),256,COMPARE,50000,0.000295556,169172551.881,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,rsa256(composite),256,COMPARE,50000,0.000009148,5465544911.978,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,rsa256(composite),256,COMPARE,50000,0.000009181,5445918631.397,0
library,NVIDIA B200,gpu,cgbn,rsa256(composite),256,COMPARE,50000,0.000008704,5744485294.118,0
opencl-kernel,NVIDIA B200,GPU,w8,rsa256(composite),256,COMPARE,50000,0.000023449,2132280488.120,0
opencl-e2e,NVIDIA B200,GPU,w8,rsa256(composite),256,COMPARE,50000,0.000405053,123440628.227,0
opencl-kernel,NVIDIA B200,GPU,w16,rsa256(composite),256,COMPARE,50000,0.000017678,2828371712.116,0
opencl-e2e,NVIDIA B200,GPU,w16,rsa256(composite),256,COMPARE,50000,0.000531306,94107722.633,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,rsa256(composite),256,COMPARE,50000,0.000013468,3712511384.833,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,rsa256(composite),256,COMPARE,50000,0.000555182,90060557.977,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,rsa256(composite),256,COMPARE,50000,0.000013352,3744751027.526,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,rsa256(composite),256,COMPARE,50000,0.000528837,94547096.061,0
opencl-kernel,NVIDIA B200,GPU,w32-il,rsa256(composite),256,COMPARE,50000,0.000012791,3908993297.777,0
opencl-e2e,NVIDIA B200,GPU,w32-il,rsa256(composite),256,COMPARE,50000,0.000593326,84270696.727,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,rsa256(composite),256,COMPARE,50000,0.000012168,4109112153.306,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,rsa256(composite),256,COMPARE,50000,0.000658728,75903855.609,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,rsa256(composite),256,COMPARE,50000,0.000597551,83674866.333,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,rsa256(composite),256,COMPARE,50000,0.001200047,41665035.745,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,rsa256(composite),256,REDUCE,6250,0.000128027,48817794.834,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,rsa256(composite),256,REDUCE,6250,0.000002717,2300723325.066,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,rsa256(composite),256,REDUCE,6250,0.000245963,25410312.899,0
library,NVIDIA B200,gpu,cgbn,rsa256(composite),256,REDUCE,50000,0.000012736,3925879396.985,0
opencl-kernel,NVIDIA B200,GPU,w8,rsa256(composite),256,REDUCE,50000,0.000141405,353594198.870,0
opencl-e2e,NVIDIA B200,GPU,w8,rsa256(composite),256,REDUCE,50000,0.000697973,71636001.877,0
opencl-kernel,NVIDIA B200,GPU,w16,rsa256(composite),256,REDUCE,50000,0.000095526,523417766.848,0
opencl-e2e,NVIDIA B200,GPU,w16,rsa256(composite),256,REDUCE,50000,0.000697703,71663732.452,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,rsa256(composite),256,REDUCE,50000,0.000037571,1330815443.154,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,rsa256(composite),256,REDUCE,50000,0.000539319,92709508.269,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,rsa256(composite),256,REDUCE,50000,0.000036772,1359730551.589,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,rsa256(composite),256,REDUCE,50000,0.000541516,92333375.598,0
opencl-kernel,NVIDIA B200,GPU,w32-il,rsa256(composite),256,REDUCE,50000,0.000038053,1313955969.444,0
opencl-e2e,NVIDIA B200,GPU,w32-il,rsa256(composite),256,REDUCE,50000,0.000701071,71319453.073,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,rsa256(composite),256,REDUCE,50000,0.000036314,1376878931.829,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,rsa256(composite),256,REDUCE,50000,0.000543560,91986162.263,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,rsa256(composite),256,REDUCE,50000,0.003489599,14328293.959,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,rsa256(composite),256,REDUCE,50000,0.004042088,12369844.331,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,rsa256(composite),256,MODMUL,3125,0.000192812,16207503.442,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,rsa256(composite),256,MODMUL,3125,0.000004369,715278821.019,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,rsa256(composite),256,MODMUL,3125,0.000020029,156020928.355,0
library,NVIDIA B200,gpu,cgbn,rsa256(composite),256,MODMUL,50000,0.000031424,1591140529.532,0
opencl-kernel,NVIDIA B200,GPU,w8,rsa256(composite),256,MODMUL,50000,0.000380810,131299072.707,0
opencl-e2e,NVIDIA B200,GPU,w8,rsa256(composite),256,MODMUL,50000,0.000882929,56629690.243,0
opencl-kernel,NVIDIA B200,GPU,w16,rsa256(composite),256,MODMUL,50000,0.000226041,221198855.216,0
opencl-e2e,NVIDIA B200,GPU,w16,rsa256(composite),256,MODMUL,50000,0.000976709,51192322.929,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,rsa256(composite),256,MODMUL,50000,0.000108968,458850301.059,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,rsa256(composite),256,MODMUL,50000,0.000842093,59375862.940,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,rsa256(composite),256,MODMUL,50000,0.000087244,573105507.726,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,rsa256(composite),256,MODMUL,50000,0.000601633,83107145.295,0
opencl-kernel,NVIDIA B200,GPU,w32-il,rsa256(composite),256,MODMUL,50000,0.000107253,466187554.923,0
opencl-e2e,NVIDIA B200,GPU,w32-il,rsa256(composite),256,MODMUL,50000,0.000907943,55069536.185,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,rsa256(composite),256,MODMUL,50000,0.000092935,538010822.430,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,rsa256(composite),256,MODMUL,50000,0.000590970,84606674.385,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,rsa256(composite),256,MODMUL,50000,0.011086192,4510114.885,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,rsa256(composite),256,MODMUL,50000,0.009924215,5038181.904,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,rsa256(composite),256,MODEXP,781,0.005438965,143593.498,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,rsa256(composite),256,MODEXP,781,0.000107242,7282592.652,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,rsa256(composite),256,MODEXP,781,0.000177054,4411081.534,0
library,NVIDIA B200,gpu,cgbn,rsa256(composite),256,MODEXP,50000,0.007989472,6258235.838,0
opencl-kernel,NVIDIA B200,GPU,w8,rsa256(composite),256,MODEXP,50000,0.017277273,2893975.224,0
opencl-e2e,NVIDIA B200,GPU,w8,rsa256(composite),256,MODEXP,50000,0.017779773,2812184.387,0
opencl-kernel,NVIDIA B200,GPU,w16,rsa256(composite),256,MODEXP,50000,0.003214069,15556604.689,0
opencl-e2e,NVIDIA B200,GPU,w16,rsa256(composite),256,MODEXP,50000,0.003932155,12715673.817,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,rsa256(composite),256,MODEXP,50000,0.002533743,19733650.642,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,rsa256(composite),256,MODEXP,50000,0.003017899,16567817.554,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,rsa256(composite),256,MODEXP,50000,0.001225611,40795979.756,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,rsa256(composite),256,MODEXP,50000,0.001837793,27206545.533,0
opencl-kernel,NVIDIA B200,GPU,w32-il,rsa256(composite),256,MODEXP,50000,0.002562968,19508631.720,0
opencl-e2e,NVIDIA B200,GPU,w32-il,rsa256(composite),256,MODEXP,50000,0.003301527,15144507.366,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,rsa256(composite),256,MODEXP,50000,0.001230940,40619364.396,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,rsa256(composite),256,MODEXP,50000,0.002002839,24964562.815,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,rsa256(composite),256,MODEXP,50000,1.979051890,25264.623,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,rsa256(composite),256,MODEXP,50000,1.940473147,25766.912,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,rsa256(composite),256,EXPONENTIATION,781,0.001910281,408840.371,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,rsa256(composite),256,EXPONENTIATION,781,0.000050576,15442108.063,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,rsa256(composite),256,EXPONENTIATION,781,0.000796976,979954.830,0
opencl-kernel,NVIDIA B200,GPU,w8,rsa256(composite),256,EXPONENTIATION,50000,0.023432676,2133772.517,0
opencl-e2e,NVIDIA B200,GPU,w8,rsa256(composite),256,EXPONENTIATION,50000,0.023958455,2086945.927,0
opencl-kernel,NVIDIA B200,GPU,w16,rsa256(composite),256,EXPONENTIATION,50000,0.004727258,10576956.033,0
opencl-e2e,NVIDIA B200,GPU,w16,rsa256(composite),256,EXPONENTIATION,50000,0.005308423,9418993.247,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,rsa256(composite),256,EXPONENTIATION,50000,0.000472733,105767937.615,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,rsa256(composite),256,EXPONENTIATION,50000,0.001260895,39654374.603,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,rsa256(composite),256,EXPONENTIATION,50000,0.000267176,187142587.441,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,rsa256(composite),256,EXPONENTIATION,50000,0.000771027,64848576.982,0
opencl-kernel,NVIDIA B200,GPU,w32-il,rsa256(composite),256,EXPONENTIATION,50000,0.000472434,105834893.254,0
opencl-e2e,NVIDIA B200,GPU,w32-il,rsa256(composite),256,EXPONENTIATION,50000,0.001120043,44641143.753,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,rsa256(composite),256,EXPONENTIATION,50000,0.000264627,188945246.198,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,rsa256(composite),256,EXPONENTIATION,50000,0.000991006,50453779.651,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,rsa256(composite),256,EXPONENTIATION,50000,0.718644670,69575.413,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,rsa256(composite),256,EXPONENTIATION,50000,0.719928336,69451.357,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,rsa256(composite),256,DIVIDE,6250,0.000146927,42538076.657,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,rsa256(composite),256,DIVIDE,6250,0.000188739,33114436.856,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,rsa256(composite),256,DIVIDE,6250,0.000020267,308376473.173,0
library,NVIDIA B200,gpu,cgbn,rsa256(composite),256,DIVIDE,50000,0.000016288,3069744597.250,0
opencl-kernel,NVIDIA B200,GPU,w8,rsa256(composite),256,DIVIDE,50000,0.000356339,140315809.220,0
opencl-e2e,NVIDIA B200,GPU,w8,rsa256(composite),256,DIVIDE,50000,0.000961352,52010084.535,0
opencl-kernel,NVIDIA B200,GPU,w16,rsa256(composite),256,DIVIDE,50000,0.000301747,165701729.831,0
opencl-e2e,NVIDIA B200,GPU,w16,rsa256(composite),256,DIVIDE,50000,0.001148605,43531065.160,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,rsa256(composite),256,DIVIDE,50000,0.000119228,419364464.767,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,rsa256(composite),256,DIVIDE,50000,0.000864994,57803868.245,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,rsa256(composite),256,DIVIDE,50000,0.000116387,429601122.673,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,rsa256(composite),256,DIVIDE,50000,0.000740531,67519116.482,0
opencl-kernel,NVIDIA B200,GPU,w32-il,rsa256(composite),256,DIVIDE,50000,0.000117187,426668345.140,0
opencl-e2e,NVIDIA B200,GPU,w32-il,rsa256(composite),256,DIVIDE,50000,0.000972665,51405158.388,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,rsa256(composite),256,DIVIDE,50000,0.000112968,442603548.481,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,rsa256(composite),256,DIVIDE,50000,0.000948490,52715369.034,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,rsa256(composite),256,DIVIDE,50000,0.004802507,10411228.970,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,rsa256(composite),256,DIVIDE,50000,0.005740128,8710607.163,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,rsa256(composite),256,ISQRT,1562,0.000082403,18955691.676,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,rsa256(composite),256,ISQRT,1562,0.000197323,7915948.173,0
opencl-kernel,NVIDIA B200,GPU,w8,rsa256(composite),256,ISQRT,50000,0.003718837,13445063.852,0
opencl-e2e,NVIDIA B200,GPU,w8,rsa256(composite),256,ISQRT,50000,0.004239963,11792555.686,0
opencl-kernel,NVIDIA B200,GPU,w16,rsa256(composite),256,ISQRT,50000,0.002858488,17491764.843,0
opencl-e2e,NVIDIA B200,GPU,w16,rsa256(composite),256,ISQRT,50000,0.003595490,13906310.823,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,rsa256(composite),256,ISQRT,50000,0.000737546,67792380.678,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,rsa256(composite),256,ISQRT,50000,0.001476168,33871482.305,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,rsa256(composite),256,ISQRT,50000,0.000702638,71160403.830,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,rsa256(composite),256,ISQRT,50000,0.001367387,36566092.511,0
opencl-kernel,NVIDIA B200,GPU,w32-il,rsa256(composite),256,ISQRT,50000,0.000748796,66773867.852,0
opencl-e2e,NVIDIA B200,GPU,w32-il,rsa256(composite),256,ISQRT,50000,0.001489918,33558894.284,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,rsa256(composite),256,ISQRT,50000,0.000695001,71942332.062,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,rsa256(composite),256,ISQRT,50000,0.001439999,34722246.960,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,rsa256(composite),256,ISQRT,50000,0.108807956,459525.221,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,rsa256(composite),256,ISQRT,50000,0.115098149,434411.852,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,rsa256(composite),256,MODMUL_R2,50000,0.003070495,16284016.706,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,rsa256(composite),256,MODMUL_R2,50000,0.000050112,997764548.075,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,rsa256(composite),256,MODMUL_R2,50000,0.000315749,158353381.294,0
library,NVIDIA B200,gpu,cgbn,rsa256(composite),256,MODMUL_R2,50000,0.000015008,3331556503.198,0
opencl-kernel,NVIDIA B200,GPU,w8,rsa256(composite),256,MODMUL_R2,50000,0.000090441,552846090.853,0
opencl-e2e,NVIDIA B200,GPU,w8,rsa256(composite),256,MODMUL_R2,50000,0.000565958,88345781.976,0
opencl-kernel,NVIDIA B200,GPU,w16,rsa256(composite),256,MODMUL_R2,50000,0.000029482,1695946398.999,0
opencl-e2e,NVIDIA B200,GPU,w16,rsa256(composite),256,MODMUL_R2,50000,0.000516963,96718728.446,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,rsa256(composite),256,MODMUL_R2,50000,0.000022222,2250016918.999,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,rsa256(composite),256,MODMUL_R2,50000,0.000518976,96343566.601,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,rsa256(composite),256,MODMUL_R2,50000,0.000017432,2868283221.584,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,rsa256(composite),256,MODMUL_R2,50000,0.000513542,97363024.166,0
opencl-kernel,NVIDIA B200,GPU,w32-il,rsa256(composite),256,MODMUL_R2,50000,0.000021447,2331319876.892,0
opencl-e2e,NVIDIA B200,GPU,w32-il,rsa256(composite),256,MODMUL_R2,50000,0.000600199,83305706.067,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,rsa256(composite),256,MODMUL_R2,50000,0.000016654,3002276921.786,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,rsa256(composite),256,MODMUL_R2,50000,0.000504799,99049333.539,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,rsa256(composite),256,MODMUL_R2,50000,0.004169448,11991994.971,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,rsa256(composite),256,MODMUL_R2,50000,0.004270375,11708573.582,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,brainpoolP512r1,512,ADD,25000,0.000364964,68499908.044,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,brainpoolP512r1,512,ADD,25000,0.000008439,2962264324.732,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,brainpoolP512r1,512,ADD,25000,0.000008874,2817326501.438,0
library,NVIDIA B200,gpu,cgbn,brainpoolP512r1,512,ADD,50000,0.000008768,5702554744.526,0
opencl-kernel,NVIDIA B200,GPU,w8,brainpoolP512r1,512,ADD,50000,0.000058343,857000077.021,0
opencl-e2e,NVIDIA B200,GPU,w8,brainpoolP512r1,512,ADD,50000,0.000796217,62796950.633,0
opencl-kernel,NVIDIA B200,GPU,w16,brainpoolP512r1,512,ADD,50000,0.000033399,1497046768.167,0
opencl-e2e,NVIDIA B200,GPU,w16,brainpoolP512r1,512,ADD,50000,0.000894100,55922160.135,0
opencl-kernel,NVIDIA B200,GPU,w32,brainpoolP512r1,512,ADD,50000,0.000020284,2464986194.824,0
opencl-e2e,NVIDIA B200,GPU,w32,brainpoolP512r1,512,ADD,50000,0.000922448,54203597.278,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,brainpoolP512r1,512,ADD,50000,0.000020359,2455923020.094,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,brainpoolP512r1,512,ADD,50000,0.000928583,53845484.423,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,brainpoolP512r1,512,ADD,50000,0.000019988,2501495256.733,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,brainpoolP512r1,512,ADD,50000,0.000748203,66826792.286,0
opencl-kernel,NVIDIA B200,GPU,w32-il,brainpoolP512r1,512,ADD,50000,0.000014454,3459247655.023,0
opencl-e2e,NVIDIA B200,GPU,w32-il,brainpoolP512r1,512,ADD,50000,0.000912256,54809178.968,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,brainpoolP512r1,512,ADD,50000,0.000019201,2604020526.750,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,brainpoolP512r1,512,ADD,50000,0.001011640,49424698.982,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,brainpoolP512r1,512,ADD,50000,0.000798243,62637569.206,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,brainpoolP512r1,512,ADD,50000,0.002383684,20975934.369,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACT,25000,0.000246746,101318765.086,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACT,25000,0.000006761,3697744570.430,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACT,25000,0.000007618,3281701069.522,0
library,NVIDIA B200,gpu,cgbn,brainpoolP512r1,512,SUBTRACT,50000,0.000008864,5640794223.827,0
opencl-kernel,NVIDIA B200,GPU,w8,brainpoolP512r1,512,SUBTRACT,50000,0.000057035,876654541.595,0
opencl-e2e,NVIDIA B200,GPU,w8,brainpoolP512r1,512,SUBTRACT,50000,0.000820205,60960372.277,0
opencl-kernel,NVIDIA B200,GPU,w16,brainpoolP512r1,512,SUBTRACT,50000,0.000034029,1469334397.055,0
opencl-e2e,NVIDIA B200,GPU,w16,brainpoolP512r1,512,SUBTRACT,50000,0.000889704,56198467.229,0
opencl-kernel,NVIDIA B200,GPU,w32,brainpoolP512r1,512,SUBTRACT,50000,0.000019784,2527298738.989,0
opencl-e2e,NVIDIA B200,GPU,w32,brainpoolP512r1,512,SUBTRACT,50000,0.000776472,64393819.346,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,brainpoolP512r1,512,SUBTRACT,50000,0.000021875,2285712086.427,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,brainpoolP512r1,512,SUBTRACT,50000,0.000777821,64282146.983,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,brainpoolP512r1,512,SUBTRACT,50000,0.000019400,2577330894.602,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,brainpoolP512r1,512,SUBTRACT,50000,0.000741825,67401343.706,0
opencl-kernel,NVIDIA B200,GPU,w32-il,brainpoolP512r1,512,SUBTRACT,50000,0.000015683,3188164209.151,0
opencl-e2e,NVIDIA B200,GPU,w32-il,brainpoolP512r1,512,SUBTRACT,50000,0.000979458,51048642.131,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,brainpoolP512r1,512,SUBTRACT,50000,0.000016868,2964192895.545,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,brainpoolP512r1,512,SUBTRACT,50000,0.000739335,67628344.790,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,brainpoolP512r1,512,SUBTRACT,50000,0.000959537,52108465.245,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,brainpoolP512r1,512,SUBTRACT,50000,0.002205318,22672467.639,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,brainpoolP512r1,512,ADDMOD,25000,0.000730906,34204143.991,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,brainpoolP512r1,512,ADDMOD,25000,0.000726280,34422004.194,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,brainpoolP512r1,512,ADDMOD,25000,0.000062600,399362288.978,0
library,NVIDIA B200,gpu,cgbn,brainpoolP512r1,512,ADDMOD,50000,0.000011136,4489942528.736,0
opencl-kernel,NVIDIA B200,GPU,w8,brainpoolP512r1,512,ADDMOD,50000,0.000072330,691276061.225,0
opencl-e2e,NVIDIA B200,GPU,w8,brainpoolP512r1,512,ADDMOD,50000,0.000831940,60100489.567,0
opencl-kernel,NVIDIA B200,GPU,w16,brainpoolP512r1,512,ADDMOD,50000,0.000039403,1268939635.831,0
opencl-e2e,NVIDIA B200,GPU,w16,brainpoolP512r1,512,ADDMOD,50000,0.000903267,55354621.495,0
opencl-kernel,NVIDIA B200,GPU,w32,brainpoolP512r1,512,ADDMOD,50000,0.000021742,2299688533.596,0
opencl-e2e,NVIDIA B200,GPU,w32,brainpoolP512r1,512,ADDMOD,50000,0.000753535,66353922.250,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,brainpoolP512r1,512,ADDMOD,50000,0.000020913,2390848073.390,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,brainpoolP512r1,512,ADDMOD,50000,0.000909841,54954654.656,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,brainpoolP512r1,512,ADDMOD,50000,0.000020700,2415467713.471,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,brainpoolP512r1,512,ADDMOD,50000,0.000744504,67158813.367,0
opencl-kernel,NVIDIA B200,GPU,w32-il,brainpoolP512r1,512,ADDMOD,50000,0.000013629,3668654585.213,0
opencl-e2e,NVIDIA B200,GPU,w32-il,brainpoolP512r1,512,ADDMOD,50000,0.000783558,63811485.263,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,brainpoolP512r1,512,ADDMOD,50000,0.000014384,3476101957.801,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,brainpoolP512r1,512,ADDMOD,50000,0.000897599,55704166.937,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,brainpoolP512r1,512,ADDMOD,50000,0.001313803,38057455.552,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,brainpoolP512r1,512,ADDMOD,50000,0.001919086,26054069.831,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000769951,32469598.366,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000014785,1690857448.805,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000101742,245720556.332,0
library,NVIDIA B200,gpu,cgbn,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000011008,4542151162.791,0
opencl-kernel,NVIDIA B200,GPU,w8,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000077809,642599397.343,0
opencl-e2e,NVIDIA B200,GPU,w8,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000935114,53469417.203,0
opencl-kernel,NVIDIA B200,GPU,w16,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000041272,1211474341.935,0
opencl-e2e,NVIDIA B200,GPU,w16,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000988650,50574014.085,0
opencl-kernel,NVIDIA B200,GPU,w32,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000022716,2201091224.780,0
opencl-e2e,NVIDIA B200,GPU,w32,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000856076,58406036.089,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000021062,2373946250.573,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000763360,65499899.896,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000020350,2457004831.670,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000846362,59076378.099,0
opencl-kernel,NVIDIA B200,GPU,w32-il,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000013983,3575772227.819,0
opencl-e2e,NVIDIA B200,GPU,w32-il,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000746308,66996466.958,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000014630,3417655204.902,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000879509,56849899.847,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000964021,51866092.154,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,brainpoolP512r1,512,SUBTRACTMOD,50000,0.002461271,20314707.227,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.000814126,30707769.910,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.000016340,1529987010.507,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.000015885,1573849961.670,0
opencl-kernel,NVIDIA B200,GPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.004877874,10250367.175,0
opencl-e2e,NVIDIA B200,GPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.006033201,8287474.588,0
opencl-kernel,NVIDIA B200,GPU,w16,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.001992140,25098638.083,0
opencl-e2e,NVIDIA B200,GPU,w16,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.003173527,15755341.188,0
opencl-kernel,NVIDIA B200,GPU,w32,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000845774,59117442.077,0
opencl-e2e,NVIDIA B200,GPU,w32,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.002027234,24664147.532,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000037155,1345713984.566,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.001150120,43473722.602,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000036126,1384044630.059,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.001073137,46592376.800,0
opencl-kernel,NVIDIA B200,GPU,w32-il,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000026539,1884022010.107,0
opencl-e2e,NVIDIA B200,GPU,w32-il,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.001272147,39303632.214,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000029516,1693993198.733,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.001133728,44102290.739,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.003718851,13445013.346,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.005958063,8391989.218,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000808308,30928803.699,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000018039,1385886852.045,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000017338,1441898047.800,0
library,NVIDIA B200,gpu,cgbn,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000010144,4929022082.019,0
opencl-kernel,NVIDIA B200,GPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.001146886,43596310.882,0
opencl-e2e,NVIDIA B200,GPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.002239245,22328954.725,0
opencl-kernel,NVIDIA B200,GPU,w16,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000319862,156317392.453,0
opencl-e2e,NVIDIA B200,GPU,w16,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.001492795,33494215.680,0
opencl-kernel,NVIDIA B200,GPU,w32,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000093243,536232794.974,0
opencl-e2e,NVIDIA B200,GPU,w32,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.001301863,38406499.055,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000092556,540213482.926,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.001196941,41773152.993,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000091942,543820943.660,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.001138901,43901972.599,0
opencl-kernel,NVIDIA B200,GPU,w32-il,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000049037,1019637318.773,0
opencl-e2e,NVIDIA B200,GPU,w32-il,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.001298668,38500987.775,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000057426,870685017.617,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.001039752,48088388.542,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.003562738,14034150.318,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.005486630,9113062.277,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.007202395,3471067.608,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000107733,232055028.852,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000048708,513262136.860,0
library,NVIDIA B200,gpu,cgbn,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000014720,3396739130.435,0
opencl-kernel,NVIDIA B200,GPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000252477,198037833.399,0
opencl-e2e,NVIDIA B200,GPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.001234140,40514041.698,0
opencl-kernel,NVIDIA B200,GPU,w16,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000078442,637413298.665,0
opencl-e2e,NVIDIA B200,GPU,w16,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000888112,56299205.435,0
opencl-kernel,NVIDIA B200,GPU,w32,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000026193,1908908289.925,0
opencl-e2e,NVIDIA B200,GPU,w32,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000889437,56215341.028,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000038355,1303610770.123,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000894917,55871106.730,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000025650,1949316166.512,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000786803,63548311.673,0
opencl-kernel,NVIDIA B200,GPU,w32-il,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000032826,1523183884.982,0
opencl-e2e,NVIDIA B200,GPU,w32-il,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000759991,65790252.583,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000021988,2273961379.953,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000749104,66746410.135,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.010435680,4791254.599,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.011132682,4491280.727,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,brainpoolP512r1,512,COMPARE,25000,0.000211606,118144156.444,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,brainpoolP512r1,512,COMPARE,25000,0.000004687,5334130618.267,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,brainpoolP512r1,512,COMPARE,25000,0.000004626,5403888433.703,0
library,NVIDIA B200,gpu,cgbn,brainpoolP512r1,512,COMPARE,50000,0.000009280,5387931034.483,0
opencl-kernel,NVIDIA B200,GPU,w8,brainpoolP512r1,512,COMPARE,50000,0.000041847,1194828757.341,0
opencl-e2e,NVIDIA B200,GPU,w8,brainpoolP512r1,512,COMPARE,50000,0.000948608,52708809.114,0
opencl-kernel,NVIDIA B200,GPU,w16,brainpoolP512r1,512,COMPARE,50000,0.000027211,1837489912.339,0
opencl-e2e,NVIDIA B200,GPU,w16,brainpoolP512r1,512,COMPARE,50000,0.000830135,60231170.197,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,brainpoolP512r1,512,COMPARE,50000,0.000016102,3105207169.143,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,brainpoolP512r1,512,COMPARE,50000,0.000731244,68376619.953,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,brainpoolP512r1,512,COMPARE,50000,0.000016613,3009682418.976,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,brainpoolP512r1,512,COMPARE,50000,0.000725641,68904604.114,0
opencl-kernel,NVIDIA B200,GPU,w32-il,brainpoolP512r1,512,COMPARE,50000,0.000015050,3322246688.171,0
opencl-e2e,NVIDIA B200,GPU,w32-il,brainpoolP512r1,512,COMPARE,50000,0.000906362,55165598.188,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,brainpoolP512r1,512,COMPARE,50000,0.000015084,3314785286.718,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,brainpoolP512r1,512,COMPARE,50000,0.000730459,68450111.856,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,brainpoolP512r1,512,COMPARE,50000,0.000827827,60399099.702,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,brainpoolP512r1,512,COMPARE,50000,0.001526468,32755352.384,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,brainpoolP512r1,512,REDUCE,3125,0.000068011,45948242.513,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,brainpoolP512r1,512,REDUCE,3125,0.000001399,2233672404.927,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,brainpoolP512r1,512,REDUCE,3125,0.000012204,256072585.147,0
library,NVIDIA B200,gpu,cgbn,brainpoolP512r1,512,REDUCE,50000,0.000019456,2569901315.789,0
opencl-kernel,NVIDIA B200,GPU,w8,brainpoolP512r1,512,REDUCE,50000,0.000400262,124918193.114,0
opencl-e2e,NVIDIA B200,GPU,w8,brainpoolP512r1,512,REDUCE,50000,0.001443171,34645927.588,0
opencl-kernel,NVIDIA B200,GPU,w16,brainpoolP512r1,512,REDUCE,50000,0.000329169,151897613.368,0
opencl-e2e,NVIDIA B200,GPU,w16,brainpoolP512r1,512,REDUCE,50000,0.001185981,42159192.320,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,brainpoolP512r1,512,REDUCE,50000,0.000095698,522476679.480,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,brainpoolP512r1,512,REDUCE,50000,0.000816666,61224535.617,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,brainpoolP512r1,512,REDUCE,50000,0.000091853,544348215.171,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,brainpoolP512r1,512,REDUCE,50000,0.000908852,55014459.447,0
opencl-kernel,NVIDIA B200,GPU,w32-il,brainpoolP512r1,512,REDUCE,50000,0.000092023,543342116.609,0
opencl-e2e,NVIDIA B200,GPU,w32-il,brainpoolP512r1,512,REDUCE,50000,0.000826626,60486848.848,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,brainpoolP512r1,512,REDUCE,50000,0.000106444,469730561.373,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,brainpoolP512r1,512,REDUCE,50000,0.001016683,49179540.423,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,brainpoolP512r1,512,REDUCE,50000,0.007027433,7114973.540,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,brainpoolP512r1,512,REDUCE,50000,0.009880248,5060601.722,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL,1562,0.000212663,7344956.526,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL,1562,0.000004196,372281791.012,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL,1562,0.000015080,103581772.378,0
library,NVIDIA B200,gpu,cgbn,brainpoolP512r1,512,MODMUL,50000,0.000084544,591408024.224,0
opencl-kernel,NVIDIA B200,GPU,w8,brainpoolP512r1,512,MODMUL,50000,0.001158470,43160377.607,0
opencl-e2e,NVIDIA B200,GPU,w8,brainpoolP512r1,512,MODMUL,50000,0.002161749,23129419.138,0
opencl-kernel,NVIDIA B200,GPU,w16,brainpoolP512r1,512,MODMUL,50000,0.000822004,60826959.734,0
opencl-e2e,NVIDIA B200,GPU,w16,brainpoolP512r1,512,MODMUL,50000,0.001808299,27650294.795,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,brainpoolP512r1,512,MODMUL,50000,0.000337694,148063007.270,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,brainpoolP512r1,512,MODMUL,50000,0.001366501,36589801.906,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,brainpoolP512r1,512,MODMUL,50000,0.000289456,172737832.890,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,brainpoolP512r1,512,MODMUL,50000,0.001114159,44876897.876,0
opencl-kernel,NVIDIA B200,GPU,w32-il,brainpoolP512r1,512,MODMUL,50000,0.000332709,150281453.448,0
opencl-e2e,NVIDIA B200,GPU,w32-il,brainpoolP512r1,512,MODMUL,50000,0.001219764,40991535.980,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,brainpoolP512r1,512,MODMUL,50000,0.000337017,148360467.766,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,brainpoolP512r1,512,MODMUL,50000,0.001365596,36614050.832,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,brainpoolP512r1,512,MODMUL,50000,0.095151770,525476.299,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,brainpoolP512r1,512,MODMUL,50000,0.100690250,496572.409,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,brainpoolP512r1,512,MODEXP,390,0.014591250,26728.347,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,brainpoolP512r1,512,MODEXP,390,0.000293492,1328826.724,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,brainpoolP512r1,512,MODEXP,390,0.000273967,1423528.875,0
library,NVIDIA B200,gpu,cgbn,brainpoolP512r1,512,MODEXP,50000,0.024437215,2046059.668,0
opencl-kernel,NVIDIA B200,GPU,w8,brainpoolP512r1,512,MODEXP,50000,0.221602253,225629.475,0
opencl-e2e,NVIDIA B200,GPU,w8,brainpoolP512r1,512,MODEXP,50000,0.222837134,224379.120,0
opencl-kernel,NVIDIA B200,GPU,w16,brainpoolP512r1,512,MODEXP,50000,0.022489820,2223228.107,0
opencl-e2e,NVIDIA B200,GPU,w16,brainpoolP512r1,512,MODEXP,50000,0.023492128,2128372.544,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,brainpoolP512r1,512,MODEXP,50000,0.019811588,2523775.466,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,brainpoolP512r1,512,MODEXP,50000,0.020680119,2417781.066,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,brainpoolP512r1,512,MODEXP,50000,0.008532473,5859965.802,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,brainpoolP512r1,512,MODEXP,50000,0.009399720,5319307.383,0
opencl-kernel,NVIDIA B200,GPU,w32-il,brainpoolP512r1,512,MODEXP,50000,0.019978409,2502701.800,0
opencl-e2e,NVIDIA B200,GPU,w32-il,brainpoolP512r1,512,MODEXP,50000,0.020868410,2395965.959,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,brainpoolP512r1,512,MODEXP,50000,0.008695586,5750043.756,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,brainpoolP512r1,512,MODEXP,50000,0.009724349,5141732.420,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,brainpoolP512r1,512,MODEXP,50000,15.470551177,3231.947,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,brainpoolP512r1,512,MODEXP,50000,15.511474692,3223.420,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,brainpoolP512r1,512,EXPONENTIATION,390,0.003376679,115498.114,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,brainpoolP512r1,512,EXPONENTIATION,390,0.000109766,3553028.827,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,brainpoolP512r1,512,EXPONENTIATION,390,0.001067594,365307.575,0
opencl-kernel,NVIDIA B200,GPU,w8,brainpoolP512r1,512,EXPONENTIATION,50000,0.189625680,263677.367,0
opencl-e2e,NVIDIA B200,GPU,w8,brainpoolP512r1,512,EXPONENTIATION,50000,0.191178753,261535.339,0
opencl-kernel,NVIDIA B200,GPU,w16,brainpoolP512r1,512,EXPONENTIATION,50000,0.047443527,1053884.547,0
opencl-e2e,NVIDIA B200,GPU,w16,brainpoolP512r1,512,EXPONENTIATION,50000,0.048424801,1032528.766,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,brainpoolP512r1,512,EXPONENTIATION,50000,0.011459971,4363012.800,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,brainpoolP512r1,512,EXPONENTIATION,50000,0.012463125,4011834.897,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,brainpoolP512r1,512,EXPONENTIATION,50000,0.008132254,6148356.896,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,brainpoolP512r1,512,EXPONENTIATION,50000,0.009263139,5397738.291,0
opencl-kernel,NVIDIA B200,GPU,w32-il,brainpoolP512r1,512,EXPONENTIATION,50000,0.011357277,4402463.750,0
opencl-e2e,NVIDIA B200,GPU,w32-il,brainpoolP512r1,512,EXPONENTIATION,50000,0.012201611,4097819.558,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,brainpoolP512r1,512,EXPONENTIATION,50000,0.011810177,4233636.801,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,brainpoolP512r1,512,EXPONENTIATION,50000,0.012805640,3904529.567,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,brainpoolP512r1,512,EXPONENTIATION,50000,6.834823998,7315.477,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,brainpoolP512r1,512,EXPONENTIATION,50000,7.053126300,7089.055,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,brainpoolP512r1,512,DIVIDE,3125,0.000078883,39615466.291,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,brainpoolP512r1,512,DIVIDE,3125,0.000004030,775406686.181,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,brainpoolP512r1,512,DIVIDE,3125,0.000013586,230023502.682,0
library,NVIDIA B200,gpu,cgbn,brainpoolP512r1,512,DIVIDE,50000,0.000025504,1960476787.955,0
opencl-kernel,NVIDIA B200,GPU,w8,brainpoolP512r1,512,DIVIDE,50000,0.001169316,42760040.731,0
opencl-e2e,NVIDIA B200,GPU,w8,brainpoolP512r1,512,DIVIDE,50000,0.002806719,17814394.397,0
opencl-kernel,NVIDIA B200,GPU,w16,brainpoolP512r1,512,DIVIDE,50000,0.001181273,42327217.367,0
opencl-e2e,NVIDIA B200,GPU,w16,brainpoolP512r1,512,DIVIDE,50000,0.002367324,21120894.640,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,brainpoolP512r1,512,DIVIDE,50000,0.000398919,125338732.210,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,brainpoolP512r1,512,DIVIDE,50000,0.001570417,31838677.980,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,brainpoolP512r1,512,DIVIDE,50000,0.000369335,135378428.804,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,brainpoolP512r1,512,DIVIDE,50000,0.001588639,31473479.314,0
opencl-kernel,NVIDIA B200,GPU,w32-il,brainpoolP512r1,512,DIVIDE,50000,0.000386140,129486728.067,0
opencl-e2e,NVIDIA B200,GPU,w32-il,brainpoolP512r1,512,DIVIDE,50000,0.001364501,36643433.353,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,brainpoolP512r1,512,DIVIDE,50000,0.000347554,143862542.598,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,brainpoolP512r1,512,DIVIDE,50000,0.001622169,30822928.185,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,brainpoolP512r1,512,DIVIDE,50000,0.081522818,613325.217,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,brainpoolP512r1,512,DIVIDE,50000,0.013160407,3799274.608,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,brainpoolP512r1,512,ISQRT,781,0.000070206,11124414.670,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,brainpoolP512r1,512,ISQRT,781,0.000001560,500641169.325,0
opencl-kernel,NVIDIA B200,GPU,w8,brainpoolP512r1,512,ISQRT,50000,0.017514417,2854791.007,0
opencl-e2e,NVIDIA B200,GPU,w8,brainpoolP512r1,512,ISQRT,50000,0.018854634,2651867.967,0
opencl-kernel,NVIDIA B200,GPU,w16,brainpoolP512r1,512,ISQRT,50000,0.017701686,2824589.720,0
opencl-e2e,NVIDIA B200,GPU,w16,brainpoolP512r1,512,ISQRT,50000,0.018706042,2672933.167,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,brainpoolP512r1,512,ISQRT,50000,0.003414194,14644744.923,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,brainpoolP512r1,512,ISQRT,50000,0.004247874,11770593.889,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,brainpoolP512r1,512,ISQRT,50000,0.003118604,16032814.909,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,brainpoolP512r1,512,ISQRT,50000,0.003976805,12572907.085,0
opencl-kernel,NVIDIA B200,GPU,w32-il,brainpoolP512r1,512,ISQRT,50000,0.003241116,15426785.342,0
opencl-e2e,NVIDIA B200,GPU,w32-il,brainpoolP512r1,512,ISQRT,50000,0.004147600,12055164.336,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,brainpoolP512r1,512,ISQRT,50000,0.004061833,12309713.300,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,brainpoolP512r1,512,ISQRT,50000,0.004972359,10055589.209,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,brainpoolP512r1,512,ISQRT,50000,0.612200705,81672.562,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,brainpoolP512r1,512,ISQRT,50000,0.619646360,80691.186,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL_R2,25000,0.003348895,7465147.916,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL_R2,25000,0.000050330,496721531.510,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL_R2,25000,0.000222678,112269732.137,0
library,NVIDIA B200,gpu,cgbn,brainpoolP512r1,512,MODMUL_R2,50000,0.000021376,2339071856.287,0
opencl-kernel,NVIDIA B200,GPU,w8,brainpoolP512r1,512,MODMUL_R2,50000,0.000394408,126772273.527,0
opencl-e2e,NVIDIA B200,GPU,w8,brainpoolP512r1,512,MODMUL_R2,50000,0.001635819,30565729.315,0
opencl-kernel,NVIDIA B200,GPU,w16,brainpoolP512r1,512,MODMUL_R2,50000,0.000077910,641765952.525,0
opencl-e2e,NVIDIA B200,GPU,w16,brainpoolP512r1,512,MODMUL_R2,50000,0.000894768,55880411.314,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,brainpoolP512r1,512,MODMUL_R2,50000,0.000058021,861757978.796,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,brainpoolP512r1,512,MODMUL_R2,50000,0.000775848,64445618.602,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,brainpoolP512r1,512,MODMUL_R2,50000,0.000032502,1538367167.879,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,brainpoolP512r1,512,MODMUL_R2,50000,0.001073953,46556972.358,0
opencl-kernel,NVIDIA B200,GPU,w32-il,brainpoolP512r1,512,MODMUL_R2,50000,0.000053059,942347900.380,0
opencl-e2e,NVIDIA B200,GPU,w32-il,brainpoolP512r1,512,MODMUL_R2,50000,0.001007046,49650160.754,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,brainpoolP512r1,512,MODMUL_R2,50000,0.000031553,1584636580.848,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,brainpoolP512r1,512,MODMUL_R2,50000,0.000809020,61803168.232,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,brainpoolP512r1,512,MODMUL_R2,50000,0.012806668,3904216.164,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,brainpoolP512r1,512,MODMUL_R2,50000,0.080395276,621927.090,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,p1024,1024,ADD,12500,0.000255199,48981382.440,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,p1024,1024,ADD,12500,0.000288077,43391251.260,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,p1024,1024,ADD,12500,0.000009170,1363079248.349,0
library,NVIDIA B200,gpu,cgbn,p1024,1024,ADD,50000,0.000010944,4568713450.292,0
opencl-kernel,NVIDIA B200,GPU,w8,p1024,1024,ADD,50000,0.000193024,259035229.421,0
opencl-e2e,NVIDIA B200,GPU,w8,p1024,1024,ADD,50000,0.001968956,25394167.161,0
opencl-kernel,NVIDIA B200,GPU,w16,p1024,1024,ADD,50000,0.000101122,494451884.803,0
opencl-e2e,NVIDIA B200,GPU,w16,p1024,1024,ADD,50000,0.001936749,25816458.764,0
opencl-kernel,NVIDIA B200,GPU,w32,p1024,1024,ADD,50000,0.000046457,1076262970.007,0
opencl-e2e,NVIDIA B200,GPU,w32,p1024,1024,ADD,50000,0.001508580,33143749.867,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,p1024,1024,ADD,50000,0.000042864,1166479891.580,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,p1024,1024,ADD,50000,0.001614274,30973676.256,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,p1024,1024,ADD,50000,0.000042786,1168606367.389,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,p1024,1024,ADD,50000,0.001382602,36163699.419,0
opencl-kernel,NVIDIA B200,GPU,w32-il,p1024,1024,ADD,50000,0.000022751,2197701118.559,0
opencl-e2e,NVIDIA B200,GPU,w32-il,p1024,1024,ADD,50000,0.001342727,37237648.546,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,p1024,1024,ADD,50000,0.000028060,1781892718.870,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,p1024,1024,ADD,50000,0.001298581,38503566.076,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,p1024,1024,ADD,50000,0.001187167,42117073.176,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,p1024,1024,ADD,50000,0.003564034,14029046.835,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,p1024,1024,SUBTRACT,12500,0.000181084,69028845.879,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,p1024,1024,SUBTRACT,12500,0.000003860,3238405211.648,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,p1024,1024,SUBTRACT,12500,0.000059782,209092793.063,0
library,NVIDIA B200,gpu,cgbn,p1024,1024,SUBTRACT,50000,0.000011136,4489942528.736,0
opencl-kernel,NVIDIA B200,GPU,w8,p1024,1024,SUBTRACT,50000,0.000196378,254611031.766,0
opencl-e2e,NVIDIA B200,GPU,w8,p1024,1024,SUBTRACT,50000,0.001986502,25169871.827,0
opencl-kernel,NVIDIA B200,GPU,w16,p1024,1024,SUBTRACT,50000,0.000102409,488238122.878,0
opencl-e2e,NVIDIA B200,GPU,w16,p1024,1024,SUBTRACT,50000,0.001931110,25891843.947,0
opencl-kernel,NVIDIA B200,GPU,w32,p1024,1024,SUBTRACT,50000,0.000046505,1075152962.263,0
opencl-e2e,NVIDIA B200,GPU,w32,p1024,1024,SUBTRACT,50000,0.001424938,35089246.612,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,p1024,1024,SUBTRACT,50000,0.000043245,1156202162.209,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,p1024,1024,SUBTRACT,50000,0.001430678,34948463.407,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,p1024,1024,SUBTRACT,50000,0.000042139,1186550109.263,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,p1024,1024,SUBTRACT,50000,0.001602501,31201228.066,0
opencl-kernel,NVIDIA B200,GPU,w32-il,p1024,1024,SUBTRACT,50000,0.000027903,1791921603.772,0
opencl-e2e,NVIDIA B200,GPU,w32-il,p1024,1024,SUBTRACT,50000,0.001513971,33025731.324,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,p1024,1024,SUBTRACT,50000,0.000027898,1792243136.011,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,p1024,1024,SUBTRACT,50000,0.001294005,38639728.133,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,p1024,1024,SUBTRACT,50000,0.001015574,49233242.814,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,p1024,1024,SUBTRACT,50000,0.003634550,13756861.574,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,p1024,1024,ADDMOD,12500,0.000532999,23452223.186,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,p1024,1024,ADDMOD,12500,0.000011261,1110009535.558,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,p1024,1024,ADDMOD,12500,0.000107350,116441017.481,0
library,NVIDIA B200,gpu,cgbn,p1024,1024,ADDMOD,50000,0.000012576,3975826972.010,0
opencl-kernel,NVIDIA B200,GPU,w8,p1024,1024,ADDMOD,50000,0.000248608,201119878.475,0
opencl-e2e,NVIDIA B200,GPU,w8,p1024,1024,ADDMOD,50000,0.002058565,24288763.853,0
opencl-kernel,NVIDIA B200,GPU,w16,p1024,1024,ADDMOD,50000,0.000130341,383609376.038,0
opencl-e2e,NVIDIA B200,GPU,w16,p1024,1024,ADDMOD,50000,0.001920558,26034101.453,0
opencl-kernel,NVIDIA B200,GPU,w32,p1024,1024,ADDMOD,50000,0.000052357,954982578.050,0
opencl-e2e,NVIDIA B200,GPU,w32,p1024,1024,ADDMOD,50000,0.001485583,33656822.364,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,p1024,1024,ADDMOD,50000,0.000040537,1233441686.339,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,p1024,1024,ADDMOD,50000,0.001576871,31708363.637,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,p1024,1024,ADDMOD,50000,0.000040778,1226152591.070,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,p1024,1024,ADDMOD,50000,0.001488022,33601652.823,0
opencl-kernel,NVIDIA B200,GPU,w32-il,p1024,1024,ADDMOD,50000,0.000016868,2964192895.545,0
opencl-e2e,NVIDIA B200,GPU,w32-il,p1024,1024,ADDMOD,50000,0.001303178,38367746.763,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,p1024,1024,ADDMOD,50000,0.000016820,2972666004.070,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,p1024,1024,ADDMOD,50000,0.001520961,32873951.478,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,p1024,1024,ADDMOD,50000,0.000992995,50352721.401,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,p1024,1024,ADDMOD,50000,0.003299369,15154412.808,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,p1024,1024,SUBTRACTMOD,12500,0.000443621,28177226.976,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,p1024,1024,SUBTRACTMOD,12500,0.000008726,1432464311.856,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,p1024,1024,SUBTRACTMOD,12500,0.000132903,94053243.821,0
library,NVIDIA B200,gpu,cgbn,p1024,1024,SUBTRACTMOD,50000,0.000013216,3783292978.208,0
opencl-kernel,NVIDIA B200,GPU,w8,p1024,1024,SUBTRACTMOD,50000,0.000263172,189989838.939,0
opencl-e2e,NVIDIA B200,GPU,w8,p1024,1024,SUBTRACTMOD,50000,0.002030839,24620366.512,0
opencl-kernel,NVIDIA B200,GPU,w16,p1024,1024,SUBTRACTMOD,50000,0.000128811,388165629.683,0
opencl-e2e,NVIDIA B200,GPU,w16,p1024,1024,SUBTRACTMOD,50000,0.001810623,27614802.882,0
opencl-kernel,NVIDIA B200,GPU,w32,p1024,1024,SUBTRACTMOD,50000,0.000051368,973368105.011,0
opencl-e2e,NVIDIA B200,GPU,w32,p1024,1024,SUBTRACTMOD,50000,0.001569823,31850724.315,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,p1024,1024,SUBTRACTMOD,50000,0.000040046,1248566041.071,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,p1024,1024,SUBTRACTMOD,50000,0.001303646,38353969.885,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,p1024,1024,SUBTRACTMOD,50000,0.000041257,1211915319.347,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,p1024,1024,SUBTRACTMOD,50000,0.001625340,30762794.885,0
opencl-kernel,NVIDIA B200,GPU,w32-il,p1024,1024,SUBTRACTMOD,50000,0.000017577,2844612941.597,0
opencl-e2e,NVIDIA B200,GPU,w32-il,p1024,1024,SUBTRACTMOD,50000,0.001559767,32056071.130,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,p1024,1024,SUBTRACTMOD,50000,0.000016287,3069938883.806,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,p1024,1024,SUBTRACTMOD,50000,0.001349485,37051170.998,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,p1024,1024,SUBTRACTMOD,50000,0.000939153,53239460.471,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,p1024,1024,SUBTRACTMOD,50000,0.003083131,16217280.193,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.001476785,8464333.143,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.000026968,463512356.395,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.000052123,239819048.087,0
opencl-kernel,NVIDIA B200,GPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.022864408,2186804.926,0
opencl-e2e,NVIDIA B200,GPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.025032931,1997368.985,0
opencl-kernel,NVIDIA B200,GPU,w16,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.006871182,7276768.387,0
opencl-e2e,NVIDIA B200,GPU,w16,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.009201855,5433686.983,0
opencl-kernel,NVIDIA B200,GPU,w32,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.004363819,11457853.708,0
opencl-e2e,NVIDIA B200,GPU,w32,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.006300188,7936271.084,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.000097692,511812568.639,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.001961138,25495401.360,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.000081889,610582751.439,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.001817555,27509483.132,0
opencl-kernel,NVIDIA B200,GPU,w32-il,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.000088054,567833450.471,0
opencl-e2e,NVIDIA B200,GPU,w32-il,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.001807964,27655417.036,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.000065174,767177639.326,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.001838251,27199767.353,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.012522941,3992672.332,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.017764997,2814523.418,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.001501646,8324198.861,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.000027657,451970511.304,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.000034428,363073211.671,0
library,NVIDIA B200,gpu,cgbn,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000019904,2512057877.814,0
opencl-kernel,NVIDIA B200,GPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.008534379,5858657.037,0
opencl-e2e,NVIDIA B200,GPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.010768098,4643345.582,0
opencl-kernel,NVIDIA B200,GPU,w16,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.002196994,22758368.644,0
opencl-e2e,NVIDIA B200,GPU,w16,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.004461512,11206963.097,0
opencl-kernel,NVIDIA B200,GPU,w32,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000568579,87938535.021,0
opencl-e2e,NVIDIA B200,GPU,w32,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.002474179,20208723.572,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000570331,87668388.804,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.002271781,22009163.958,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000567008,88182170.682,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.002297627,21761582.322,0
opencl-kernel,NVIDIA B200,GPU,w32-il,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000200194,249757642.875,0
opencl-e2e,NVIDIA B200,GPU,w32-il,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.002115063,23639957.546,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000187184,267116734.073,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.001944568,25712653.290,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.011414673,4380326.991,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.013477415,3709910.256,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.011181190,1117948.981,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.000142307,87838252.367,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.000073002,171228296.140,0
library,NVIDIA B200,gpu,cgbn,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000029408,1700217627.856,0
opencl-kernel,NVIDIA B200,GPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.001926351,25955809.858,0
opencl-e2e,NVIDIA B200,GPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.003411990,14654204.695,0
opencl-kernel,NVIDIA B200,GPU,w16,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000254513,196453621.315,0
opencl-e2e,NVIDIA B200,GPU,w16,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.001964383,25453284.044,0
opencl-kernel,NVIDIA B200,GPU,w32,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000068434,730631698.600,0
opencl-e2e,NVIDIA B200,GPU,w32,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.001493187,33485425.862,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000116225,430200106.775,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.001553323,32189055.263,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000080164,623721839.144,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.001633614,30606986.385,0
opencl-kernel,NVIDIA B200,GPU,w32-il,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000098964,505234397.137,0
opencl-e2e,NVIDIA B200,GPU,w32-il,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.001492534,33500075.510,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000054349,919978857.675,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.001345707,37155190.897,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.181637295,275273.864,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.116282914,429985.785,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,p1024,1024,COMPARE,12500,0.000100823,123979170.470,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,p1024,1024,COMPARE,12500,0.000003027,4129287573.068,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,p1024,1024,COMPARE,12500,0.000003942,3170610336.376,0
library,NVIDIA B200,gpu,cgbn,p1024,1024,COMPARE,50000,0.000011136,4489942528.736,0
opencl-kernel,NVIDIA B200,GPU,w8,p1024,1024,COMPARE,50000,0.000086942,575096112.883,0
opencl-e2e,NVIDIA B200,GPU,w8,p1024,1024,COMPARE,50000,0.001501109,33308707.420,0
opencl-kernel,NVIDIA B200,GPU,w16,p1024,1024,COMPARE,50000,0.000052019,961188995.589,0
opencl-e2e,NVIDIA B200,GPU,w16,p1024,1024,COMPARE,50000,0.001765532,28320074.642,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,p1024,1024,COMPARE,50000,0.000022513,2220941280.865,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,p1024,1024,COMPARE,50000,0.001512117,33066224.692,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,p1024,1024,COMPARE,50000,0.000022879,2185411463.957,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,p1024,1024,COMPARE,50000,0.001304823,38319370.933,0
opencl-kernel,NVIDIA B200,GPU,w32-il,p1024,1024,COMPARE,50000,0.000017625,2836890622.668,0
opencl-e2e,NVIDIA B200,GPU,w32-il,p1024,1024,COMPARE,50000,0.001450649,34467331.590,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,p1024,1024,COMPARE,50000,0.000018613,2686285327.579,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,p1024,1024,COMPARE,50000,0.001442007,34673897.751,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,p1024,1024,COMPARE,50000,0.001532771,32620658.019,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,p1024,1024,COMPARE,50000,0.003715958,13455480.460,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,p1024,1024,REDUCE,1562,0.000021844,71505685.761,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,p1024,1024,REDUCE,1562,0.000000830,1881896517.445,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,p1024,1024,REDUCE,1562,0.000015218,102644607.410,0
library,NVIDIA B200,gpu,cgbn,p1024,1024,REDUCE,50000,0.000026080,1917177914.110,0
opencl-kernel,NVIDIA B200,GPU,w8,p1024,1024,REDUCE,50000,0.002220078,22521730.993,0
opencl-e2e,NVIDIA B200,GPU,w8,p1024,1024,REDUCE,50000,0.003929054,12725709.498,0
opencl-kernel,NVIDIA B200,GPU,w16,p1024,1024,REDUCE,50000,0.001168072,42805580.814,0
opencl-e2e,NVIDIA B200,GPU,w16,p1024,1024,REDUCE,50000,0.003020568,16553177.791,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,p1024,1024,REDUCE,50000,0.000282598,176929797.512,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,p1024,1024,REDUCE,50000,0.001834696,27252472.354,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,p1024,1024,REDUCE,50000,0.000341556,146388864.136,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,p1024,1024,REDUCE,50000,0.001966524,25425572.543,0
opencl-kernel,NVIDIA B200,GPU,w32-il,p1024,1024,REDUCE,50000,0.000320885,155819072.037,0
opencl-e2e,NVIDIA B200,GPU,w32-il,p1024,1024,REDUCE,50000,0.001734736,28822829.089,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,p1024,1024,REDUCE,50000,0.000328708,152110699.750,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,p1024,1024,REDUCE,50000,0.001614796,30963663.565,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,p1024,1024,REDUCE,50000,0.098803653,506054.164,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,p1024,1024,REDUCE,50000,0.101485989,492678.847,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,p1024,1024,MODMUL,781,0.000298578,2615728.787,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,p1024,1024,MODMUL,781,0.000006228,125404873.509,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,p1024,1024,MODMUL,781,0.000022447,34793345.013,0
library,NVIDIA B200,gpu,cgbn,p1024,1024,MODMUL,50000,0.000212672,235103821.848,0
opencl-kernel,NVIDIA B200,GPU,w8,p1024,1024,MODMUL,50000,0.008250133,6060508.345,0
opencl-e2e,NVIDIA B200,GPU,w8,p1024,1024,MODMUL,50000,0.009693144,5158285.069,0
opencl-kernel,NVIDIA B200,GPU,w16,p1024,1024,MODMUL,50000,0.003085219,16206304.606,0
opencl-e2e,NVIDIA B200,GPU,w16,p1024,1024,MODMUL,50000,0.004893787,10217036.475,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,p1024,1024,MODMUL,50000,0.001451465,34447952.685,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,p1024,1024,MODMUL,50000,0.002837677,17620046.128,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,p1024,1024,MODMUL,50000,0.001094033,45702460.099,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,p1024,1024,MODMUL,50000,0.002452245,20389479.993,0
opencl-kernel,NVIDIA B200,GPU,w32-il,p1024,1024,MODMUL,50000,0.001434221,34862129.544,0
opencl-e2e,NVIDIA B200,GPU,w32-il,p1024,1024,MODMUL,50000,0.002863923,17458570.116,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,p1024,1024,MODMUL,50000,0.001090965,45830984.457,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,p1024,1024,MODMUL,50000,0.002462041,20308354.069,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,p1024,1024,MODMUL,50000,0.412303886,121269.776,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,p1024,1024,MODMUL,50000,0.407912650,122575.262,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,p1024,1024,MODEXP,195,0.046328750,4209.049,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,p1024,1024,MODEXP,195,0.001247576,156303.106,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,p1024,1024,MODEXP,195,0.001009479,193168.942,0
library,NVIDIA B200,gpu,cgbn,p1024,1024,MODEXP,50000,0.135919064,367865.982,0
opencl-kernel,NVIDIA B200,GPU,w8,p1024,1024,MODEXP,50000,2.545114651,19645.480,0
opencl-e2e,NVIDIA B200,GPU,w8,p1024,1024,MODEXP,50000,2.592670893,19285.132,0
opencl-kernel,NVIDIA B200,GPU,w16,p1024,1024,MODEXP,50000,0.155739292,321049.360,0
opencl-e2e,NVIDIA B200,GPU,w16,p1024,1024,MODEXP,50000,0.157711364,317034.858,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,p1024,1024,MODEXP,50000,0.156801257,318874.995,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,p1024,1024,MODEXP,50000,0.158281011,315893.863,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,p1024,1024,MODEXP,50000,0.081181383,615904.757,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,p1024,1024,MODEXP,50000,0.082645828,604991.216,0
opencl-kernel,NVIDIA B200,GPU,w32-il,p1024,1024,MODEXP,50000,0.159411356,313653.941,0
opencl-e2e,NVIDIA B200,GPU,w32-il,p1024,1024,MODEXP,50000,0.161054814,310453.309,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,p1024,1024,MODEXP,50000,0.080655214,619922.724,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,p1024,1024,MODEXP,50000,0.082004966,609719.173,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,p1024,1024,MODEXP,50000,138.674670182,360.556,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,p1024,1024,MODEXP,50000,138.782773579,360.275,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,p1024,1024,EXPONENTIATION,195,0.006380928,30559.818,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,p1024,1024,EXPONENTIATION,195,0.000306155,636932.134,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,p1024,1024,EXPONENTIATION,195,0.005264667,37039.380,0
opencl-kernel,NVIDIA B200,GPU,w8,p1024,1024,EXPONENTIATION,50000,1.629350840,30687.068,0
opencl-e2e,NVIDIA B200,GPU,w8,p1024,1024,EXPONENTIATION,50000,1.632659738,30624.875,0
opencl-kernel,NVIDIA B200,GPU,w16,p1024,1024,EXPONENTIATION,50000,0.371844025,134464.982,0
opencl-e2e,NVIDIA B200,GPU,w16,p1024,1024,EXPONENTIATION,50000,0.373308733,133937.397,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,p1024,1024,EXPONENTIATION,50000,0.111641164,447863.478,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,p1024,1024,EXPONENTIATION,50000,0.113073769,442189.205,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,p1024,1024,EXPONENTIATION,50000,0.096236938,519551.027,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,p1024,1024,EXPONENTIATION,50000,0.097621081,512184.454,0
opencl-kernel,NVIDIA B200,GPU,w32-il,p1024,1024,EXPONENTIATION,50000,0.110729438,451551.104,0
opencl-e2e,NVIDIA B200,GPU,w32-il,p1024,1024,EXPONENTIATION,50000,0.112234410,445496.172,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,p1024,1024,EXPONENTIATION,50000,0.096649077,517335.515,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,p1024,1024,EXPONENTIATION,50000,0.098115442,509603.779,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,p1024,1024,EXPONENTIATION,50000,53.602803859,932.787,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,p1024,1024,EXPONENTIATION,50000,53.895404631,927.723,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,p1024,1024,DIVIDE,1562,0.000042552,36707695.636,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,p1024,1024,DIVIDE,1562,0.000001515,1031178047.609,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,p1024,1024,DIVIDE,1562,0.000008844,176616904.189,0
library,NVIDIA B200,gpu,cgbn,p1024,1024,DIVIDE,50000,0.000034048,1468515037.594,0
opencl-kernel,NVIDIA B200,GPU,w8,p1024,1024,DIVIDE,50000,0.010516125,4754603.045,0
opencl-e2e,NVIDIA B200,GPU,w8,p1024,1024,DIVIDE,50000,0.012356908,4046319.658,0
opencl-kernel,NVIDIA B200,GPU,w16,p1024,1024,DIVIDE,50000,0.008904322,5615250.703,0
opencl-e2e,NVIDIA B200,GPU,w16,p1024,1024,DIVIDE,50000,0.010758792,4647361.875,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,p1024,1024,DIVIDE,50000,0.001565230,31944185.504,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,p1024,1024,DIVIDE,50000,0.003365007,14858810.117,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,p1024,1024,DIVIDE,50000,0.001591809,31410804.126,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,p1024,1024,DIVIDE,50000,0.003532253,14155271.079,0
opencl-kernel,NVIDIA B200,GPU,w32-il,p1024,1024,DIVIDE,50000,0.001550388,32249990.434,0
opencl-e2e,NVIDIA B200,GPU,w32-il,p1024,1024,DIVIDE,50000,0.003399045,14710013.746,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,p1024,1024,DIVIDE,50000,0.001610468,31046875.327,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,p1024,1024,DIVIDE,50000,0.003310258,15104562.689,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,p1024,1024,DIVIDE,50000,0.196209027,254830.273,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,p1024,1024,DIVIDE,50000,0.200490707,249388.117,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,p1024,1024,ISQRT,390,0.000064031,6090796.988,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,p1024,1024,ISQRT,390,0.000002331,167341796.747,0
opencl-kernel,NVIDIA B200,GPU,w8,p1024,1024,ISQRT,50000,0.194049357,257666.404,0
opencl-e2e,NVIDIA B200,GPU,w8,p1024,1024,ISQRT,50000,0.198012073,252509.856,0
opencl-kernel,NVIDIA B200,GPU,w16,p1024,1024,ISQRT,50000,0.123533640,404748.051,0
opencl-e2e,NVIDIA B200,GPU,w16,p1024,1024,ISQRT,50000,0.125434814,398613.418,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,p1024,1024,ISQRT,50000,0.023092022,2165249.979,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,p1024,1024,ISQRT,50000,0.024693767,2024802.457,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,p1024,1024,ISQRT,50000,0.022812495,2191781.302,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,p1024,1024,ISQRT,50000,0.024390753,2049957.210,0
opencl-kernel,NVIDIA B200,GPU,w32-il,p1024,1024,ISQRT,50000,0.023509284,2126819.344,0
opencl-e2e,NVIDIA B200,GPU,w32-il,p1024,1024,ISQRT,50000,0.025092903,1992595.281,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,p1024,1024,ISQRT,50000,0.022820831,2190980.695,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,p1024,1024,ISQRT,50000,0.024263717,2060690.042,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,p1024,1024,ISQRT,50000,2.897829244,17254.295,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,p1024,1024,ISQRT,50000,2.845134633,17573.861,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,p1024,1024,MODMUL_R2,12500,0.004896758,2552709.360,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,p1024,1024,MODMUL_R2,12500,0.000081049,154227602.750,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,p1024,1024,MODMUL_R2,12500,0.000235571,53062545.045,0
library,NVIDIA B200,gpu,cgbn,p1024,1024,MODMUL_R2,50000,0.000052064,960356484.327,0
opencl-kernel,NVIDIA B200,GPU,w8,p1024,1024,MODMUL_R2,50000,0.002499380,20004960.970,0
opencl-e2e,NVIDIA B200,GPU,w8,p1024,1024,MODMUL_R2,50000,0.004267128,11717482.860,0
opencl-kernel,NVIDIA B200,GPU,w16,p1024,1024,MODMUL_R2,50000,0.000234445,213269553.715,0
opencl-e2e,NVIDIA B200,GPU,w16,p1024,1024,MODMUL_R2,50000,0.002045572,24443040.430,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,p1024,1024,MODMUL_R2,50000,0.000186238,268473713.504,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,p1024,1024,MODMUL_R2,50000,0.001758056,28440504.832,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,p1024,1024,MODMUL_R2,50000,0.000121753,410667619.257,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,p1024,1024,MODMUL_R2,50000,0.001506774,33183477.021,0
opencl-kernel,NVIDIA B200,GPU,w32-il,p1024,1024,MODMUL_R2,50000,0.000170499,293256703.374,0
opencl-e2e,NVIDIA B200,GPU,w32-il,p1024,1024,MODMUL_R2,50000,0.001451804,34439911.727,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,p1024,1024,MODMUL_R2,50000,0.000100145,499276049.563,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,p1024,1024,MODMUL_R2,50000,0.001357347,36836564.106,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,p1024,1024,MODMUL_R2,50000,0.195217622,256124.419,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,p1024,1024,MODMUL_R2,50000,0.199546985,250567.554,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,p2048,2048,ADD,6250,0.000185199,33747567.343,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,p2048,2048,ADD,6250,0.000002961,2111064067.430,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,p2048,2048,ADD,6250,0.000003289,1900238082.191,0
library,NVIDIA B200,gpu,cgbn,p2048,2048,ADD,50000,0.000012896,3877171215.881,0
opencl-kernel,NVIDIA B200,GPU,w8,p2048,2048,ADD,50000,0.000368461,135699567.275,0
opencl-e2e,NVIDIA B200,GPU,w8,p2048,2048,ADD,50000,0.003710947,13473649.720,0
opencl-kernel,NVIDIA B200,GPU,w16,p2048,2048,ADD,50000,0.000185950,268889543.216,0
opencl-e2e,NVIDIA B200,GPU,w16,p2048,2048,ADD,50000,0.003521698,14197696.901,0
opencl-kernel,NVIDIA B200,GPU,w32,p2048,2048,ADD,50000,0.000074663,669676042.099,0
opencl-e2e,NVIDIA B200,GPU,w32,p2048,2048,ADD,50000,0.002603076,19208044.274,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,p2048,2048,ADD,50000,0.000073907,676526219.489,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,p2048,2048,ADD,50000,0.002583843,19351021.260,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,p2048,2048,ADD,50000,0.000075341,663649543.399,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,p2048,2048,ADD,50000,0.002584615,19345240.796,0
opencl-kernel,NVIDIA B200,GPU,w32-il,p2048,2048,ADD,50000,0.000041918,1192804596.846,0
opencl-e2e,NVIDIA B200,GPU,w32-il,p2048,2048,ADD,50000,0.002569615,19458167.661,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,p2048,2048,ADD,50000,0.000042589,1174011184.216,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,p2048,2048,ADD,50000,0.002550609,19603161.870,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,p2048,2048,ADD,50000,0.000848718,58912384.012,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,p2048,2048,ADD,50000,0.004786828,10445330.248,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,p2048,2048,SUBTRACT,6250,0.000160556,38927330.725,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,p2048,2048,SUBTRACT,6250,0.000002751,2271685038.416,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,p2048,2048,SUBTRACT,6250,0.000003215,1943853717.230,0
library,NVIDIA B200,gpu,cgbn,p2048,2048,SUBTRACT,50000,0.000013056,3829656862.745,0
opencl-kernel,NVIDIA B200,GPU,w8,p2048,2048,SUBTRACT,50000,0.000367987,135874376.333,0
opencl-e2e,NVIDIA B200,GPU,w8,p2048,2048,SUBTRACT,50000,0.003685734,13565818.904,0
opencl-kernel,NVIDIA B200,GPU,w16,p2048,2048,SUBTRACT,50000,0.000186884,267545532.672,0
opencl-e2e,NVIDIA B200,GPU,w16,p2048,2048,SUBTRACT,50000,0.003574617,13987512.609,0
opencl-kernel,NVIDIA B200,GPU,w32,p2048,2048,SUBTRACT,50000,0.000075133,665486079.072,0
opencl-e2e,NVIDIA B200,GPU,w32,p2048,2048,SUBTRACT,50000,0.002648201,18880742.371,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,p2048,2048,SUBTRACT,50000,0.000073738,678076001.295,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,p2048,2048,SUBTRACT,50000,0.002624698,19049810.804,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,p2048,2048,SUBTRACT,50000,0.000076140,656685676.193,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,p2048,2048,SUBTRACT,50000,0.002613220,19133482.723,0
opencl-kernel,NVIDIA B200,GPU,w32-il,p2048,2048,SUBTRACT,50000,0.000041709,1198780638.499,0
opencl-e2e,NVIDIA B200,GPU,w32-il,p2048,2048,SUBTRACT,50000,0.002796746,17877919.710,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,p2048,2048,SUBTRACT,50000,0.000042152,1186183083.989,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,p2048,2048,SUBTRACT,50000,0.002675080,18691029.987,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,p2048,2048,SUBTRACT,50000,0.001853171,26980781.660,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,p2048,2048,SUBTRACT,50000,0.005644836,8857653.265,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,p2048,2048,ADDMOD,6250,0.000369522,16913730.467,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,p2048,2048,ADDMOD,6250,0.000006363,982239924.422,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,p2048,2048,ADDMOD,6250,0.000520662,12003940.223,0
library,NVIDIA B200,gpu,cgbn,p2048,2048,ADDMOD,50000,0.000013024,3839066339.066,0
opencl-kernel,NVIDIA B200,GPU,w8,p2048,2048,ADDMOD,50000,0.000465028,107520410.158,0
opencl-e2e,NVIDIA B200,GPU,w8,p2048,2048,ADDMOD,50000,0.003746694,13345098.536,0
opencl-kernel,NVIDIA B200,GPU,w16,p2048,2048,ADDMOD,50000,0.000239085,209130594.213,0
opencl-e2e,NVIDIA B200,GPU,w16,p2048,2048,ADDMOD,50000,0.002721566,18371775.605,0
opencl-kernel,NVIDIA B200,GPU,w32,p2048,2048,ADDMOD,50000,0.000089185,560632625.021,0
opencl-e2e,NVIDIA B200,GPU,w32,p2048,2048,ADDMOD,50000,0.002711164,18442263.515,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,p2048,2048,ADDMOD,50000,0.000070215,712098858.312,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,p2048,2048,ADDMOD,50000,0.002571295,19445454.410,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,p2048,2048,ADDMOD,50000,0.000071350,700771800.532,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,p2048,2048,ADDMOD,50000,0.002587323,19324993.733,0
opencl-kernel,NVIDIA B200,GPU,w32-il,p2048,2048,ADDMOD,50000,0.000026536,1884236909.389,0
opencl-e2e,NVIDIA B200,GPU,w32-il,p2048,2048,ADDMOD,50000,0.002615212,19118908.927,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,p2048,2048,ADDMOD,50000,0.000025344,1972847212.729,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,p2048,2048,ADDMOD,50000,0.002548335,19620655.138,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,p2048,2048,ADDMOD,50000,0.001020520,48994630.467,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,p2048,2048,ADDMOD,50000,0.006032259,8288768.648,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,p2048,2048,SUBTRACTMOD,6250,0.000324878,19237980.677,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,p2048,2048,SUBTRACTMOD,6250,0.000005770,1083110731.830,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,p2048,2048,SUBTRACTMOD,6250,0.000051913,120393922.270,0
library,NVIDIA B200,gpu,cgbn,p2048,2048,SUBTRACTMOD,50000,0.000014080,3551136363.636,0
opencl-kernel,NVIDIA B200,GPU,w8,p2048,2048,SUBTRACTMOD,50000,0.000476769,104872603.747,0
opencl-e2e,NVIDIA B200,GPU,w8,p2048,2048,SUBTRACTMOD,50000,0.003783435,13215503.975,0
opencl-kernel,NVIDIA B200,GPU,w16,p2048,2048,SUBTRACTMOD,50000,0.000247015,202416828.012,0
opencl-e2e,NVIDIA B200,GPU,w16,p2048,2048,SUBTRACTMOD,50000,0.002685636,18617563.618,0
opencl-kernel,NVIDIA B200,GPU,w32,p2048,2048,SUBTRACTMOD,50000,0.000093135,536855477.405,0
opencl-e2e,NVIDIA B200,GPU,w32,p2048,2048,SUBTRACTMOD,50000,0.002644005,18910705.549,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,p2048,2048,SUBTRACTMOD,50000,0.000070101,713256592.650,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,p2048,2048,SUBTRACTMOD,50000,0.002789743,17922797.831,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,p2048,2048,SUBTRACTMOD,50000,0.000071900,695410626.668,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,p2048,2048,SUBTRACTMOD,50000,0.002756113,18141491.601,0
opencl-kernel,NVIDIA B200,GPU,w32-il,p2048,2048,SUBTRACTMOD,50000,0.000028153,1776005464.928,0
opencl-e2e,NVIDIA B200,GPU,w32-il,p2048,2048,SUBTRACTMOD,50000,0.002495650,20034860.905,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,p2048,2048,SUBTRACTMOD,50000,0.000026484,1887939205.697,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,p2048,2048,SUBTRACTMOD,50000,0.002693516,18563097.664,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,p2048,2048,SUBTRACTMOD,50000,0.000920572,54314058.748,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,p2048,2048,SUBTRACTMOD,50000,0.004863675,10280291.951,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.002484331,2515767.331,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.000043945,142224806.844,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.001834506,3406910.870,0
opencl-kernel,NVIDIA B200,GPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.091998481,543487.235,0
opencl-e2e,NVIDIA B200,GPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.095720597,522353.617,0
opencl-kernel,NVIDIA B200,GPU,w16,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.023954527,2087288.137,0
opencl-e2e,NVIDIA B200,GPU,w16,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.027168754,1840349.394,0
opencl-kernel,NVIDIA B200,GPU,w32,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.007110679,7031677.344,0
opencl-e2e,NVIDIA B200,GPU,w32,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.010609194,4712893.392,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.000278836,179316874.577,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.003617047,13823431.168,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.000229354,218003594.473,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.003565539,14023125.162,0
opencl-kernel,NVIDIA B200,GPU,w32-il,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.000310286,161141622.842,0
opencl-e2e,NVIDIA B200,GPU,w32-il,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.003581488,13960678.097,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.000182905,273365952.664,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.003499285,14288633.377,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.119159704,419604.936,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.198023659,252495.082,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.002484130,2515970.881,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.000042810,145993966.307,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.000047159,132531758.369,0
library,NVIDIA B200,gpu,cgbn,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.000059520,840053763.441,0
opencl-kernel,NVIDIA B200,GPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.033731553,1482291.669,0
opencl-e2e,NVIDIA B200,GPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.037784258,1323302.418,0
opencl-kernel,NVIDIA B200,GPU,w16,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.008569414,5834704.664,0
opencl-e2e,NVIDIA B200,GPU,w16,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.012002370,4165843.882,0
opencl-kernel,NVIDIA B200,GPU,w32,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.002174191,22997060.188,0
opencl-e2e,NVIDIA B200,GPU,w32,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.005487255,9112024.240,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.002171335,23027308.859,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.005403641,9253020.474,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.002171585,23024657.242,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.005436128,9197723.076,0
opencl-kernel,NVIDIA B200,GPU,w32-il,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.000846771,59047837.719,0
opencl-e2e,NVIDIA B200,GPU,w32-il,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.004199965,11904860.993,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.000794144,62960873.748,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.004136961,12086166.429,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.102965514,485599.480,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.119091578,419844.970,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.019303917,323768.486,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.000213072,29332816.398,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.000092484,67579300.405,0
library,NVIDIA B200,gpu,cgbn,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000088640,564079422.383,0
opencl-kernel,NVIDIA B200,GPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.039934931,1252036.720,0
opencl-e2e,NVIDIA B200,GPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.043705986,1144008.054,0
opencl-kernel,NVIDIA B200,GPU,w16,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.001327660,37660242.888,0
opencl-e2e,NVIDIA B200,GPU,w16,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.003845176,13003305.699,0
opencl-kernel,NVIDIA B200,GPU,w32,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000204830,244104770.679,0
opencl-e2e,NVIDIA B200,GPU,w32,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.002718714,18391047.675,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000367793,135946026.832,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.003033714,16481448.045,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000252381,198113195.962,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.002943342,16987492.625,0
opencl-kernel,NVIDIA B200,GPU,w32-il,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000345444,144741229.558,0
opencl-e2e,NVIDIA B200,GPU,w32-il,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.002936089,17029456.172,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000228009,219289633.827,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.002777167,18003958.621,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.714873948,69942.401,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.703952951,71027.474,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,p2048,2048,COMPARE,6250,0.000037863,165068235.802,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,p2048,2048,COMPARE,6250,0.000001396,4476872744.684,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,p2048,2048,COMPARE,6250,0.000001355,4611232388.072,0
library,NVIDIA B200,gpu,cgbn,p2048,2048,COMPARE,50000,0.000013152,3801703163.017,0
opencl-kernel,NVIDIA B200,GPU,w8,p2048,2048,COMPARE,50000,0.000164969,303087229.522,0
opencl-e2e,NVIDIA B200,GPU,w8,p2048,2048,COMPARE,50000,0.003518191,14211849.484,0
opencl-kernel,NVIDIA B200,GPU,w16,p2048,2048,COMPARE,50000,0.000086671,576894403.477,0
opencl-e2e,NVIDIA B200,GPU,w16,p2048,2048,COMPARE,50000,0.002572273,19438061.045,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,p2048,2048,COMPARE,50000,0.000033297,1501637057.678,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,p2048,2048,COMPARE,50000,0.002571341,19445106.662,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,p2048,2048,COMPARE,50000,0.000034949,1430659072.846,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,p2048,2048,COMPARE,50000,0.002621280,19074650.338,0
opencl-kernel,NVIDIA B200,GPU,w32-il,p2048,2048,COMPARE,50000,0.000024623,2030621387.168,0
opencl-e2e,NVIDIA B200,GPU,w32-il,p2048,2048,COMPARE,50000,0.002740440,18245245.563,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,p2048,2048,COMPARE,50000,0.000025871,1932658943.172,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,p2048,2048,COMPARE,50000,0.002769704,18052470.740,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,p2048,2048,COMPARE,50000,0.001016947,49166766.310,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,p2048,2048,COMPARE,50000,0.005743182,8705975.142,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,p2048,2048,REDUCE,781,0.000016028,48727913.744,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,p2048,2048,REDUCE,781,0.000000562,1390619071.912,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,p2048,2048,REDUCE,781,0.000006421,121640704.338,0
library,NVIDIA B200,gpu,cgbn,p2048,2048,REDUCE,50000,0.000026368,1896237864.078,0
opencl-kernel,NVIDIA B200,GPU,w8,p2048,2048,REDUCE,50000,0.108963885,458867.633,0
opencl-e2e,NVIDIA B200,GPU,w8,p2048,2048,REDUCE,50000,0.113104931,442067.376,0
opencl-kernel,NVIDIA B200,GPU,w16,p2048,2048,REDUCE,50000,0.004808530,10398188.266,0
opencl-e2e,NVIDIA B200,GPU,w16,p2048,2048,REDUCE,50000,0.007257370,6889548.099,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,p2048,2048,REDUCE,50000,0.000903363,55348729.238,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,p2048,2048,REDUCE,50000,0.003451665,14485762.953,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,p2048,2048,REDUCE,50000,0.000975197,51271688.370,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,p2048,2048,REDUCE,50000,0.003653723,13684671.592,0
opencl-kernel,NVIDIA B200,GPU,w32-il,p2048,2048,REDUCE,50000,0.001067820,46824370.924,0
opencl-e2e,NVIDIA B200,GPU,w32-il,p2048,2048,REDUCE,50000,0.003570165,14004954.877,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,p2048,2048,REDUCE,50000,0.000923839,54121981.618,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,p2048,2048,REDUCE,50000,0.003462186,14441743.014,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,p2048,2048,REDUCE,50000,0.389642516,128322.752,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,p2048,2048,REDUCE,50000,0.403473734,123923.804,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,p2048,2048,MODMUL,390,0.000507371,768668.734,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,p2048,2048,MODMUL,390,0.000340560,1145173.332,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,p2048,2048,MODMUL,390,0.000298865,1304937.877,0
library,NVIDIA B200,gpu,cgbn,p2048,2048,MODMUL,50000,0.000410400,121832358.674,0
opencl-kernel,NVIDIA B200,GPU,w8,p2048,2048,MODMUL,50000,0.180053435,277695.341,0
opencl-e2e,NVIDIA B200,GPU,w8,p2048,2048,MODMUL,50000,0.183460964,272537.541,0
opencl-kernel,NVIDIA B200,GPU,w16,p2048,2048,MODMUL,50000,0.015398938,3246977.167,0
opencl-e2e,NVIDIA B200,GPU,w16,p2048,2048,MODMUL,50000,0.018147278,2755234.147,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,p2048,2048,MODMUL,50000,0.005439628,9191805.170,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,p2048,2048,MODMUL,50000,0.008016694,6236984.981,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,p2048,2048,MODMUL,50000,0.003718906,13444814.271,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,p2048,2048,MODMUL,50000,0.006317389,7914662.111,0
opencl-kernel,NVIDIA B200,GPU,w32-il,p2048,2048,MODMUL,50000,0.005494802,9099509.104,0
opencl-e2e,NVIDIA B200,GPU,w32-il,p2048,2048,MODMUL,50000,0.008019490,6234810.403,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,p2048,2048,MODMUL,50000,0.003585841,13943730.195,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,p2048,2048,MODMUL,50000,0.006371682,7847221.592,0
opencl-kernel,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,p2048,2048,MODMUL,50000,1.708317896,29268.557,0
opencl-e2e,cpu-skylake-avx512-Intel(R) Xeon(R) Platinum 8559C,CPU,w8,p2048,2048,MODMUL,50000,1.719275875,29082.011,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,p2048,2048,MODEXP,97,0.164325335,590.292,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,p2048,2048,MODEXP,97,0.005965829,16259.266,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,p2048,2048,MODEXP,97,0.008644351,11221.201,0
library,NVIDIA B200,gpu,cgbn,p2048,2048,MODEXP,50000,0.568473935,87954.780,0
opencl-kernel,NVIDIA B200,GPU,w8,p2048,2048,MODEXP,50000,54.369481465,919.634,0
opencl-e2e,NVIDIA B200,GPU,w8,p2048,2048,MODEXP,50000,54.347270832,920.009,0
opencl-kernel,NVIDIA B200,GPU,w16,p2048,2048,MODEXP,50000,5.322783562,9393.581,0
opencl-e2e,NVIDIA B200,GPU,w16,p2048,2048,MODEXP,50000,5.296765810,9439.723,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,p2048,2048,MODEXP,50000,1.509343314,33126.989,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,p2048,2048,MODEXP,50000,1.513836685,33028.662,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,p2048,2048,MODEXP,50000,0.972694497,51403.601,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,p2048,2048,MODEXP,50000,0.974268040,51320.579,0
opencl-kernel,NVIDIA B200,GPU,w32-il,p2048,2048,MODEXP,50000,1.374828277,36368.178,0
opencl-e2e,NVIDIA B200,GPU,w32-il,p2048,2048,MODEXP,50000,1.377631162,36294.185,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,p2048,2048,MODEXP,50000,0.992403211,50382.747,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,p2048,2048,MODEXP,50000,0.992558422,50374.869,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,p2048,2048,EXPONENTIATION,97,0.019058353,5089.632,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,p2048,2048,EXPONENTIATION,97,0.003425054,28320.721,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,p2048,2048,EXPONENTIATION,97,0.011465978,8459.810,0
opencl-kernel,NVIDIA B200,GPU,w8,p2048,2048,EXPONENTIATION,50000,13.683773244,3653.963,0
opencl-e2e,NVIDIA B200,GPU,w8,p2048,2048,EXPONENTIATION,50000,13.688342862,3652.743,0
opencl-kernel,NVIDIA B200,GPU,w16,p2048,2048,EXPONENTIATION,50000,3.202554483,15612.537,0
opencl-e2e,NVIDIA B200,GPU,w16,p2048,2048,EXPONENTIATION,50000,3.209493248,15578.783,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,p2048,2048,EXPONENTIATION,50000,0.943071801,53018.232,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,p2048,2048,EXPONENTIATION,50000,0.944433360,52941.798,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,p2048,2048,EXPONENTIATION,50000,0.810076200,61722.589,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,p2048,2048,EXPONENTIATION,50000,0.815874246,61283.954,0
opencl-kernel,NVIDIA B200,GPU,w32-il,p2048,2048,EXPONENTIATION,50000,0.957685337,52209.215,0
opencl-e2e,NVIDIA B200,GPU,w32-il,p2048,2048,EXPONENTIATION,50000,0.959174474,52128.160,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,p2048,2048,EXPONENTIATION,50000,0.811183117,61638.364,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,p2048,2048,EXPONENTIATION,50000,0.813502191,61462.649,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,p2048,2048,DIVIDE,781,0.000027463,28438188.612,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,p2048,2048,DIVIDE,781,0.000001137,686671254.320,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,p2048,2048,DIVIDE,781,0.000005913,132076423.581,0
library,NVIDIA B200,gpu,cgbn,p2048,2048,DIVIDE,50000,0.000031360,1594387755.102,0
opencl-kernel,NVIDIA B200,GPU,w8,p2048,2048,DIVIDE,50000,0.445404737,112257.450,0
opencl-e2e,NVIDIA B200,GPU,w8,p2048,2048,DIVIDE,50000,0.450883748,110893.329,0
opencl-kernel,NVIDIA B200,GPU,w16,p2048,2048,DIVIDE,50000,0.153276575,326207.706,0
opencl-e2e,NVIDIA B200,GPU,w16,p2048,2048,DIVIDE,50000,0.156060170,320389.245,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,p2048,2048,DIVIDE,50000,0.016772343,2981098.095,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,p2048,2048,DIVIDE,50000,0.020060612,2492446.391,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,p2048,2048,DIVIDE,50000,0.015282165,3271787.729,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,p2048,2048,DIVIDE,50000,0.018602030,2687878.688,0
opencl-kernel,NVIDIA B200,GPU,w32-il,p2048,2048,DIVIDE,50000,0.016759002,2983471.212,0
opencl-e2e,NVIDIA B200,GPU,w32-il,p2048,2048,DIVIDE,50000,0.020119048,2485207.059,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,p2048,2048,DIVIDE,50000,0.015300883,3267785.269,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,p2048,2048,DIVIDE,50000,0.018786111,2661540.748,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,p2048,2048,ISQRT,195,0.000053126,3670488.343,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,p2048,2048,ISQRT,195,0.000002325,83888997.480,0
opencl-kernel,NVIDIA B200,GPU,w8,p2048,2048,ISQRT,50000,7.352760116,6800.167,0
opencl-e2e,NVIDIA B200,GPU,w8,p2048,2048,ISQRT,50000,7.302886179,6846.608,0
opencl-kernel,NVIDIA B200,GPU,w16,p2048,2048,ISQRT,50000,4.455297622,11222.595,0
opencl-e2e,NVIDIA B200,GPU,w16,p2048,2048,ISQRT,50000,4.434906692,11274.194,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,p2048,2048,ISQRT,50000,0.153893562,324899.881,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,p2048,2048,ISQRT,50000,0.156388153,319717.312,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,p2048,2048,ISQRT,50000,0.054041785,925210.001,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,p2048,2048,ISQRT,50000,0.056615809,883145.554,0
opencl-kernel,NVIDIA B200,GPU,w32-il,p2048,2048,ISQRT,50000,0.153792389,325113.618,0
opencl-e2e,NVIDIA B200,GPU,w32-il,p2048,2048,ISQRT,50000,0.156318825,319859.108,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,p2048,2048,ISQRT,50000,0.054167957,923054.934,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,p2048,2048,ISQRT,50000,0.056776096,880652.308,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-1t,p2048,2048,MODMUL_R2,6250,0.008387175,745185.354,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,gmp-nt,p2048,2048,MODMUL_R2,6250,0.000109175,57247483.971,0
library,Intel(R) Xeon(R) Platinum 8559C,host-cpu,openssl-nt,p2048,2048,MODMUL_R2,6250,0.000280094,22313938.759,0
library,NVIDIA B200,gpu,cgbn,p2048,2048,MODMUL_R2,50000,0.000167296,298871461.362,0
opencl-kernel,NVIDIA B200,GPU,w8,p2048,2048,MODMUL_R2,50000,0.030953125,1615345.780,0
opencl-e2e,NVIDIA B200,GPU,w8,p2048,2048,MODMUL_R2,50000,0.033788118,1479810.149,0
opencl-kernel,NVIDIA B200,GPU,w16,p2048,2048,MODMUL_R2,50000,0.001129149,44281135.866,0
opencl-e2e,NVIDIA B200,GPU,w16,p2048,2048,MODMUL_R2,50000,0.003701824,13506854.587,0
opencl-kernel,NVIDIA B200,GPU,w32-opt,p2048,2048,MODMUL_R2,50000,0.000636941,78500214.045,0
opencl-e2e,NVIDIA B200,GPU,w32-opt,p2048,2048,MODMUL_R2,50000,0.003206414,15593744.025,0
opencl-kernel,NVIDIA B200,GPU,w32-o64,p2048,2048,MODMUL_R2,50000,0.000442220,113065891.091,0
opencl-e2e,NVIDIA B200,GPU,w32-o64,p2048,2048,MODMUL_R2,50000,0.002958535,16900255.804,0
opencl-kernel,NVIDIA B200,GPU,w32-il,p2048,2048,MODMUL_R2,50000,0.000625457,79941557.456,0
opencl-e2e,NVIDIA B200,GPU,w32-il,p2048,2048,MODMUL_R2,50000,0.003196493,15642142.693,0
opencl-kernel,NVIDIA B200,GPU,w32-il64,p2048,2048,MODMUL_R2,50000,0.000387111,129161926.563,0
opencl-e2e,NVIDIA B200,GPU,w32-il64,p2048,2048,MODMUL_R2,50000,0.002889518,17303923.782,0
```
