# MPA-OpenCL benchmark report - NVIDIA H100 80GB HBM3

> **Partial report.** The run was interrupted or hit its time budget.
> Rows that never ran are marked `n/a`.

## 1. System under test

2 OpenCL device(s) exercised with the identical kernels and operands.

### Device 0 - NVIDIA H100 80GB HBM3 (GPU)

| Property | Value |
|---|---|
| Model | NVIDIA H100 80GB HBM3 |
| Type | GPU |
| Vendor | NVIDIA Corporation |
| Device memory | 79.19 GiB |
| Max single allocation | 19.80 GiB |
| Local memory | 48 KiB |
| Global cache | 4224 KiB |
| Compute units | 132 |
| Max clock | 1980 MHz |
| Max work-group size | 1024 |
| OpenCL version | OpenCL 3.0 CUDA |
| Driver | 570.172.08 |

### Device 1 - cpu-skylake-avx512-AMD EPYC 9554 64-Core Processor (CPU)

| Property | Value |
|---|---|
| Model | cpu-skylake-avx512-AMD EPYC 9554 64-Core Processor |
| Type | CPU |
| Vendor | AuthenticAMD |
| Device memory | 1509.53 GiB |
| Max single allocation | 512.00 GiB |
| Local memory | 1024 KiB |
| Global cache | 32768 KiB |
| Compute units | 224 |
| Max clock | 3764 MHz |
| Max work-group size | 4096 |
| OpenCL version | OpenCL 3.0 PoCL HSTR: cpu-x86_64-pc-linux-gnu-skylake-avx512 |
| Driver | 5.0+debian |

### Host

| Property | Value |
|---|---|
| CPU | AMD EPYC 9554 64-Core Processor |
| Logical cores | 224 |
| OpenMP threads used | 224 |
| RAM | 1511.5 GB |
| OS | Ubuntu 24.04.4 LTS |
| Kernel | 5.15.0-171-generic |
| Arch | x86_64 |
| GMP | 6.3.0 |
| OpenSSL | OpenSSL 3.0.13 30 Jan 2024 |
| CGBN | not measured (`cgbn_results_NVIDIA_H100_80GB_HBM3.tsv` absent) |

## 2. Method

- Workload auto-sized from the device and host: --min-items from 700 x compute units, --items from ten times that capped by host RAM. Either flag, given explicitly, overrides its half.
- Base workload 50000 items, scaled down per operator by its cost weight and by modulus size. Device rows honour --min-items (92400) so the GPU is not left idle; the CPU libraries keep the smaller count because a full-width MODEXP there costs minutes. Both counts appear in every row as dev/cpu, and throughput is per-second so they remain comparable.
- 5 timed repetitions, **minimum** reported. Two untimed warm-up launches precede them.
- `kernel` times `clEnqueueNDRangeKernel` + `clFinish` only. `e2e` adds the host->device operand writes and the device->host result read.
- Every OpenCL device runs the same kernels on the same operands, so GPU and CPU-OpenCL columns are directly comparable.
- CPU library baselines (GMP, OpenSSL) run those same operands, with every temporary - including each thread's GMP context, BN_CTX and Montgomery context - allocated outside the timed region, so the figure is the arithmetic and not marshalling. The generator is reseeded per modulus and operation so every backend sees identical inputs.
- Cost weighting drives the wide cells down to a few hundred items, which is tens of microseconds of work - the same order as the cost of entering an OpenMP region. Each baseline pass is therefore repeated until the timed interval reaches 5 ms and the per-pass time is reported; the multi-threaded loop enters one parallel region per interval and partitions the range itself. Without this the multi-threaded GMP figure came out up to 9x slower than the single-threaded one at 2048 bits.
- OpenSSL rows time the nearest BN primitive, which is not always semantically identical (its Montgomery routine expects Montgomery-domain inputs); they measure comparable work, not identical results. Correctness is judged against GMP only.
- Every device cell is checked word-for-word against GMP before it is timed. A cell that mismatches is reported and excluded from the speedup tables.
- Total wall time 2308.9 s.

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

**All configurations correct** - 485 configurations, 0 problems.

## 4. Throughput per device

Operations per second, higher is better. Kernel-only timings.

### Device 0 - NVIDIA H100 80GB HBM3 (GPU)

#### secp256k1 (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | GMP 224T | OpenSSL 224T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | 1.79 G | 2.62 G | 3.87 G | 3.92 G | 3.85 G | 5.32 G | 5.45 G | 71.27 M | 9.06 M | 17.57 M | n/a |
| SUBTRACT | 50000 / 50000 | 1.90 G | 2.68 G | 3.73 G | 3.87 G | 3.92 G | 5.58 G | 5.62 G | 92.79 M | 43.98 M | 43.93 M | n/a |
| ADDMOD | 50000 / 50000 | 1.29 G | 2.05 G | 3.71 G | 4.80 G | 4.71 G | 5.89 G | 5.96 G | 25.78 M | 17.59 M | 8.94 M | n/a |
| SUBTRACTMOD | 50000 / 50000 | 1.31 G | 2.07 G | 3.63 G | 4.86 G | 4.75 G | 5.90 G | 5.92 G | 33.20 M | 25.80 M | 25.80 M | n/a |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | 54.80 M | 195.16 M | 673.65 M | 2.56 G | 2.53 G | 4.76 G | 5.30 G | 57.25 M | 34.10 M | 34.54 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | 295.50 M | 888.65 M | 2.41 G | 2.31 G | 2.45 G | 2.87 G | 3.09 G | 57.05 M | 34.94 M | 34.85 M | n/a |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | 635.43 M | 1.92 G | 4.28 G | 3.71 G | 4.04 G | 3.92 G | 4.95 G | 7.12 M | 8.68 M | 4.53 M | n/a |
| COMPARE | 50000 / 50000 | 1.80 G | 2.74 G | - | 5.22 G | 5.20 G | 5.52 G | 5.52 G | 154.78 M | 27.76 M | 52.55 M | n/a |
| REDUCE | 50000 / 6250 | 362.01 M | 544.67 M | - | 1.42 G | 1.46 G | 1.49 G | 1.52 G | 67.25 M | 33.26 M | 16.63 M | n/a |
| MODMUL | 50000 / 3125 | 135.34 M | 229.09 M | - | 483.53 M | 613.40 M | 489.41 M | 602.37 M | 12.43 M | 3.67 M | 7.25 M | n/a |
| MODEXP | 50000 / 781 | 2.97 M | 15.73 M | - | 19.71 M | 41.27 M | 19.96 M | 42.06 M | 117.10 k | 4.23 M | 3.99 M | n/a |
| EXPONENTIATION | 50000 / 781 | 1.87 M | 7.07 M | - | 108.36 M | 131.85 M | 107.76 M | 122.97 M | 375.20 k | 13.53 M | 1.60 M | n/a |
| DIVIDE | 50000 / 6250 | 152.17 M | 181.12 M | - | 502.21 M | 511.83 M | 507.36 M | 530.38 M | 40.53 M | 788.51 M | 417.10 M | n/a |
| ISQRT | 50000 / 1562 | 13.71 M | 18.41 M | - | 75.70 M | 83.59 M | 73.60 M | 84.48 M | 20.47 M | 639.33 M | n/a | n/a |
| MODMUL_R2 | 50000 / 50000 | 530.83 M | 1.82 G | - | 2.81 G | 3.75 G | 2.82 G | 3.99 G | 12.30 M | 375.30 M | 203.21 M | n/a |

#### rsa256(composite) (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | GMP 224T | OpenSSL 224T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | 1.83 G | 2.70 G | 3.94 G | 3.84 G | 3.98 G | 5.40 G | 5.44 G | 69.96 M | 18.00 M | 33.18 M | n/a |
| SUBTRACT | 50000 / 50000 | 1.68 G | 2.75 G | 3.93 G | 3.81 G | 3.83 G | 5.60 G | 5.51 G | 94.88 M | 44.13 M | 43.49 M | n/a |
| ADDMOD | 50000 / 50000 | 1.40 G | 2.23 G | 3.68 G | 4.89 G | 4.81 G | 5.79 G | 5.88 G | 30.38 M | 886.45 M | 25.70 M | n/a |
| SUBTRACTMOD | 50000 / 50000 | 1.27 G | 2.08 G | 3.56 G | 4.97 G | 4.73 G | 5.79 G | 6.01 G | 33.13 M | 26.66 M | 25.88 M | n/a |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | 54.36 M | 192.37 M | 671.74 M | 2.57 G | 2.62 G | 4.67 G | 5.39 G | 55.94 M | 21.24 M | 34.65 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | 294.62 M | 905.59 M | 2.40 G | 2.31 G | 2.50 G | 2.96 G | 3.08 G | 56.99 M | 34.77 M | 30.68 M | n/a |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | 636.87 M | 1.95 G | 4.38 G | 3.70 G | 4.09 G | 3.90 G | 4.86 G | 7.10 M | 8.64 M | 7.10 M | n/a |
| COMPARE | 50000 / 50000 | 1.81 G | 2.80 G | - | 5.35 G | 5.18 G | 5.21 G | 5.55 G | 154.71 M | 52.82 M | 27.41 M | n/a |
| REDUCE | 50000 / 6250 | 357.21 M | 544.55 M | - | 1.41 G | 1.48 G | 1.47 G | 1.51 G | 42.86 M | 25.28 M | 24.76 M | n/a |
| MODMUL | 50000 / 3125 | 134.82 M | 229.14 M | - | 485.50 M | 611.37 M | 489.27 M | 603.52 M | 12.37 M | 9.31 M | 9.31 M | n/a |
| MODEXP | 50000 / 781 | 2.97 M | 15.74 M | - | 19.75 M | 41.34 M | 20.01 M | 41.90 M | 124.84 k | 71.29 k | 70.50 k | n/a |
| EXPONENTIATION | 50000 / 781 | 1.87 M | 7.06 M | - | 108.53 M | 131.16 M | 107.76 M | 123.30 M | 379.01 k | 18.65 M | 1.78 M | n/a |
| DIVIDE | 50000 / 6250 | 151.21 M | 180.41 M | - | 482.82 M | 502.97 M | 490.71 M | 511.95 M | 40.02 M | 764.45 M | 315.92 M | n/a |
| ISQRT | 50000 / 1562 | 13.71 M | 17.95 M | - | 72.87 M | 77.69 M | 71.18 M | 77.01 M | 20.63 M | 481.54 M | n/a | n/a |
| MODMUL_R2 | 50000 / 50000 | 532.81 M | 1.79 G | - | 2.73 G | 3.81 G | 2.79 G | 3.92 G | 12.36 M | 656.47 M | 211.34 M | n/a |

#### brainpoolP512r1 (512-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | GMP 224T | OpenSSL 224T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 25000 | 827.94 M | 1.50 G | 2.53 G | 2.48 G | 2.51 G | 3.41 G | 3.50 G | 67.14 M | 1.47 G | 1.34 G | n/a |
| SUBTRACT | 50000 / 25000 | 849.34 M | 1.52 G | 2.45 G | 2.46 G | 2.45 G | 3.40 G | 3.53 G | 91.18 M | 1.34 G | 39.52 M | n/a |
| ADDMOD | 50000 / 25000 | 656.80 M | 1.19 G | 2.31 G | 2.05 G | 2.08 G | 4.18 G | 4.28 G | 27.38 M | 17.26 M | 17.26 M | n/a |
| SUBTRACTMOD | 50000 / 25000 | 569.65 M | 1.10 G | 2.05 G | 2.04 G | 2.10 G | 4.21 G | 4.25 G | 29.28 M | 22.10 M | 21.39 M | n/a |
| MULTIPLYOPERANDSCANNING | 50000 / 25000 | 10.45 M | 38.66 M | 148.61 M | 831.38 M | 814.30 M | 1.58 G | 1.60 G | 23.23 M | 17.48 M | 16.50 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 25000 | 44.03 M | 159.17 M | 563.16 M | 559.88 M | 560.64 M | 838.49 M | 905.56 M | 23.28 M | 17.56 M | 688.61 M | n/a |
| MONTGOMERYMULTIPLICATION | 50000 / 25000 | 185.75 M | 615.44 M | 1.87 G | 1.33 G | 1.92 G | 1.60 G | 2.56 G | 2.99 M | 6.76 M | 4.37 M | n/a |
| COMPARE | 50000 / 25000 | 829.45 M | 1.51 G | - | 2.65 G | 2.65 G | 4.52 G | 4.65 G | 152.77 M | 48.16 M | 48.75 M | n/a |
| REDUCE | 50000 / 3125 | 119.78 M | 154.59 M | - | 556.07 M | 542.89 M | 504.90 M | 466.41 M | 41.86 M | 23.64 M | 11.95 M | n/a |
| MODMUL | 50000 / 1562 | 43.30 M | 61.17 M | - | 151.03 M | 187.96 M | 125.50 M | 151.69 M | 6.26 M | 4.15 M | 2.10 M | n/a |
| MODEXP | 50000 / 390 | 222.74 k | 2.24 M | - | 2.57 M | 5.76 M | 2.40 M | 5.84 M | 20.37 k | 64.70 k | 41.21 k | n/a |
| EXPONENTIATION | 50000 / 390 | 254.72 k | 959.36 k | - | 3.24 M | 3.50 M | 3.24 M | 3.36 M | 110.51 k | 136.28 k | 66.18 k | n/a |
| DIVIDE | 50000 / 3125 | 43.58 M | 44.35 M | - | 131.90 M | 141.31 M | 149.35 M | 154.66 M | 37.51 M | 11.91 M | 15.41 M | n/a |
| ISQRT | 50000 / 781 | 2.91 M | 2.93 M | - | 14.89 M | 16.34 M | 12.05 M | 12.81 M | 11.27 M | 5.72 M | n/a | n/a |
| MODMUL_R2 | 50000 / 25000 | 123.66 M | 610.47 M | - | 909.03 M | 1.45 G | 938.95 M | 1.68 G | 6.20 M | 4.59 M | 8.48 M | n/a |

#### p1024 (1024-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | GMP 224T | OpenSSL 224T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 12500 | 262.13 M | 514.32 M | 1.14 G | 1.16 G | 1.17 G | 2.28 G | 2.31 G | 60.06 M | 22.03 M | 10.90 M | n/a |
| SUBTRACT | 50000 / 12500 | 263.80 M | 513.84 M | 1.14 G | 1.15 G | 1.17 G | 2.13 G | 2.33 G | 53.82 M | 10.36 M | 14.64 M | n/a |
| ADDMOD | 50000 / 12500 | 187.93 M | 370.14 M | 909.37 M | 850.64 M | 862.08 M | 2.96 G | 3.06 G | 17.75 M | 12.96 M | 6.60 M | n/a |
| SUBTRACTMOD | 50000 / 12500 | 187.64 M | 372.10 M | 901.15 M | 829.86 M | 832.64 M | 3.02 G | 2.98 G | 22.77 M | 6.74 M | 12.93 M | n/a |
| MULTIPLYOPERANDSCANNING | 50000 / 12500 | 2.20 M | 8.55 M | 38.87 M | 334.91 M | 333.88 M | 532.29 M | 561.45 M | 6.36 M | 6.51 M | 3.46 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 12500 | 5.90 M | 22.97 M | 89.60 M | 89.63 M | 89.61 M | 256.51 M | 277.51 M | 6.36 M | 6.46 M | 6.47 M | n/a |
| MONTGOMERYMULTIPLICATION | 50000 / 12500 | 27.24 M | 183.28 M | 671.21 M | 385.40 M | 538.56 M | 530.20 M | 932.49 M | 943.28 k | 2.13 M | 1.83 M | n/a |
| COMPARE | 50000 / 12500 | 345.21 M | 721.02 M | - | 1.22 G | 1.22 G | 3.58 G | 3.52 G | 143.66 M | 38.52 M | 39.41 M | n/a |
| REDUCE | 50000 / 1562 | 22.74 M | 41.85 M | - | 155.06 M | 156.49 M | 158.69 M | 181.30 M | 53.36 M | 19.48 M | 19.23 M | n/a |
| MODMUL | 50000 / 781 | 5.92 M | 16.29 M | - | 34.66 M | 45.36 M | 35.40 M | 46.43 M | 2.20 M | 913.60 k | 1.28 M | n/a |
| MODEXP | 50000 / 195 | 18.77 k | 319.75 k | - | 320.83 k | 612.47 k | 317.91 k | 607.79 k | 3.36 k | 31.71 k | 23.26 k | n/a |
| EXPONENTIATION | 50000 / 195 | 30.90 k | 127.41 k | - | 442.87 k | 480.17 k | 453.05 k | 481.69 k | 24.65 k | 32.99 k | 15.94 k | n/a |
| DIVIDE | 50000 / 1562 | 3.89 M | 7.16 M | - | 32.97 M | 32.63 M | 33.79 M | 32.62 M | 32.43 M | 13.95 M | 7.15 M | n/a |
| ISQRT | 50000 / 390 | 256.52 k | 408.83 k | - | 2.19 M | 2.23 M | 2.23 M | 2.29 M | 5.48 M | 2.75 M | n/a | n/a |
| MODMUL_R2 | 50000 / 12500 | 18.14 M | 198.40 M | - | 249.41 M | 362.98 M | 303.20 M | 501.61 M | 2.17 M | 2.18 M | 1.67 M | n/a |

#### p2048 (2048-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | GMP 224T | OpenSSL 224T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 6250 | 131.41 M | 266.26 M | 635.00 M | 610.93 M | 578.42 M | 733.43 M | 671.39 M | 43.30 M | 18.61 M | 18.59 M | n/a |
| SUBTRACT | 50000 / 6250 | 131.27 M | 266.93 M | 643.44 M | 623.67 M | 585.96 M | 745.24 M | 672.65 M | 43.55 M | 20.17 M | 20.74 M | n/a |
| ADDMOD | 50000 / 6250 | 104.33 M | 206.67 M | 519.03 M | 448.24 M | 444.08 M | 1.28 G | 1.19 G | 15.23 M | 5.76 M | 10.55 M | n/a |
| SUBTRACTMOD | 50000 / 6250 | 98.85 M | 193.63 M | 460.22 M | 452.83 M | 434.96 M | 1.33 G | 1.29 G | 17.55 M | 11.97 M | 6.15 M | n/a |
| MULTIPLYOPERANDSCANNING | 50000 / 6250 | 545.23 k | 2.16 M | 9.65 M | 115.53 M | 115.13 M | 146.82 M | 157.76 M | 1.95 M | 1.13 M | 2.20 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 6250 | 1.49 M | 5.88 M | 23.17 M | 22.95 M | 22.96 M | 53.27 M | 56.97 M | 1.94 M | 1.13 M | 2.18 M | n/a |
| MONTGOMERYMULTIPLICATION | 50000 / 6250 | 1.17 M | 37.62 M | 220.31 M | 125.00 M | 171.19 M | 146.35 M | 191.91 M | 286.85 k | 1.06 M | 1.08 M | n/a |
| COMPARE | 50000 / 6250 | 172.06 M | 358.86 M | - | 670.85 M | 673.83 M | 2.30 G | 2.10 G | 143.61 M | 32.13 M | 16.82 M | n/a |
| REDUCE | 50000 / 781 | 357.83 k | 10.71 M | - | 49.93 M | 52.11 M | 56.07 M | 61.10 M | 37.45 M | 12.95 M | 12.90 M | n/a |
| MODMUL | 50000 / 390 | 160.22 k | 3.36 M | - | 9.01 M | 13.21 M | 9.12 M | 13.83 M | 682.53 k | 25.95 M | 9.08 M | n/a |
| MODEXP | 50000 / 97 | 1.21 k | 8.35 k | - | 40.70 k | 30.89 k | 39.93 k | 30.46 k | 466.0 | 11.75 k | 15.84 k | n/a |
| EXPONENTIATION | 50000 / 97 | 3.65 k | 15.81 k | - | 53.02 k | 61.92 k | 54.61 k | 60.71 k | 3.94 k | 100.28 k | 2.07 k | n/a |
| DIVIDE | 50000 / 781 | 73.52 k | 342.93 k | - | 4.11 M | 4.15 M | 4.39 M | 4.46 M | 24.11 M | 10.26 M | 10.23 M | n/a |
| ISQRT | 50000 / 195 | 5.12 k | 12.07 k | - | 818.67 k | 958.64 k | 874.11 k | 960.19 k | 3.23 M | 1.79 M | n/a | n/a |
| MODMUL_R2 | 50000 / 6250 | 2.13 M | 40.88 M | - | 71.00 M | 102.10 M | 77.56 M | 112.47 M | 661.60 k | 1.08 M | 1.06 M | n/a |

### Device 1 - cpu-skylake-avx512-AMD EPYC 9554 64-Core Processor (CPU)

#### secp256k1 (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | GMP 224T | OpenSSL 224T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | - | - | - | - | - | - | - | 71.27 M | 9.06 M | 17.57 M | n/a |
| SUBTRACT | 50000 / 50000 | - | - | - | - | - | - | - | 92.79 M | 43.98 M | 43.93 M | n/a |
| ADDMOD | 50000 / 50000 | - | - | - | - | - | - | - | 25.78 M | 17.59 M | 8.94 M | n/a |
| SUBTRACTMOD | 50000 / 50000 | - | - | - | - | - | - | - | 33.20 M | 25.80 M | 25.80 M | n/a |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | - | - | - | - | - | - | - | 57.25 M | 34.10 M | 34.54 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | - | - | - | - | - | - | - | 57.05 M | 34.94 M | 34.85 M | n/a |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | - | - | - | - | - | - | - | 7.12 M | 8.68 M | 4.53 M | n/a |
| COMPARE | 50000 / 50000 | - | - | - | - | - | - | - | 154.78 M | 27.76 M | 52.55 M | n/a |
| REDUCE | 50000 / 6250 | - | - | - | - | - | - | - | 67.25 M | 33.26 M | 16.63 M | n/a |
| MODMUL | 50000 / 3125 | - | - | - | - | - | - | - | 12.43 M | 3.67 M | 7.25 M | n/a |
| MODEXP | 50000 / 781 | - | - | - | - | - | - | - | 117.10 k | 4.23 M | 3.99 M | n/a |
| EXPONENTIATION | 50000 / 781 | - | - | - | - | - | - | - | 375.20 k | 13.53 M | 1.60 M | n/a |
| DIVIDE | 50000 / 6250 | - | - | - | - | - | - | - | 40.53 M | 788.51 M | 417.10 M | n/a |
| ISQRT | 50000 / 1562 | - | - | - | - | - | - | - | 20.47 M | 639.33 M | n/a | n/a |
| MODMUL_R2 | 50000 / 50000 | - | - | - | - | - | - | - | 12.30 M | 375.30 M | 203.21 M | n/a |

#### rsa256(composite) (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | GMP 224T | OpenSSL 224T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | - | - | - | - | - | - | - | 69.96 M | 18.00 M | 33.18 M | n/a |
| SUBTRACT | 50000 / 50000 | - | - | - | - | - | - | - | 94.88 M | 44.13 M | 43.49 M | n/a |
| ADDMOD | 50000 / 50000 | - | - | - | - | - | - | - | 30.38 M | 886.45 M | 25.70 M | n/a |
| SUBTRACTMOD | 50000 / 50000 | - | - | - | - | - | - | - | 33.13 M | 26.66 M | 25.88 M | n/a |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | - | - | - | - | - | - | - | 55.94 M | 21.24 M | 34.65 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | - | - | - | - | - | - | - | 56.99 M | 34.77 M | 30.68 M | n/a |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | - | - | - | - | - | - | - | 7.10 M | 8.64 M | 7.10 M | n/a |
| COMPARE | 50000 / 50000 | - | - | - | - | - | - | - | 154.71 M | 52.82 M | 27.41 M | n/a |
| REDUCE | 50000 / 6250 | - | - | - | - | - | - | - | 42.86 M | 25.28 M | 24.76 M | n/a |
| MODMUL | 50000 / 3125 | - | - | - | - | - | - | - | 12.37 M | 9.31 M | 9.31 M | n/a |
| MODEXP | 50000 / 781 | - | - | - | - | - | - | - | 124.84 k | 71.29 k | 70.50 k | n/a |
| EXPONENTIATION | 50000 / 781 | - | - | - | - | - | - | - | 379.01 k | 18.65 M | 1.78 M | n/a |
| DIVIDE | 50000 / 6250 | - | - | - | - | - | - | - | 40.02 M | 764.45 M | 315.92 M | n/a |
| ISQRT | 50000 / 1562 | - | - | - | - | - | - | - | 20.63 M | 481.54 M | n/a | n/a |
| MODMUL_R2 | 50000 / 50000 | - | - | - | - | - | - | - | 12.36 M | 656.47 M | 211.34 M | n/a |

#### brainpoolP512r1 (512-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | GMP 224T | OpenSSL 224T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 25000 | - | - | - | - | - | - | - | 67.14 M | 1.47 G | 1.34 G | n/a |
| SUBTRACT | 50000 / 25000 | - | - | - | - | - | - | - | 91.18 M | 1.34 G | 39.52 M | n/a |
| ADDMOD | 50000 / 25000 | - | - | - | - | - | - | - | 27.38 M | 17.26 M | 17.26 M | n/a |
| SUBTRACTMOD | 50000 / 25000 | - | - | - | - | - | - | - | 29.28 M | 22.10 M | 21.39 M | n/a |
| MULTIPLYOPERANDSCANNING | 50000 / 25000 | - | - | - | - | - | - | - | 23.23 M | 17.48 M | 16.50 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 25000 | - | - | - | - | - | - | - | 23.28 M | 17.56 M | 688.61 M | n/a |
| MONTGOMERYMULTIPLICATION | 50000 / 25000 | - | - | - | - | - | - | - | 2.99 M | 6.76 M | 4.37 M | n/a |
| COMPARE | 50000 / 25000 | - | - | - | - | - | - | - | 152.77 M | 48.16 M | 48.75 M | n/a |
| REDUCE | 50000 / 3125 | - | - | - | - | - | - | - | 41.86 M | 23.64 M | 11.95 M | n/a |
| MODMUL | 50000 / 1562 | - | - | - | - | - | - | - | 6.26 M | 4.15 M | 2.10 M | n/a |
| MODEXP | 50000 / 390 | - | - | - | - | - | - | - | 20.37 k | 64.70 k | 41.21 k | n/a |
| EXPONENTIATION | 50000 / 390 | - | - | - | - | - | - | - | 110.51 k | 136.28 k | 66.18 k | n/a |
| DIVIDE | 50000 / 3125 | - | - | - | - | - | - | - | 37.51 M | 11.91 M | 15.41 M | n/a |
| ISQRT | 50000 / 781 | - | - | - | - | - | - | - | 11.27 M | 5.72 M | n/a | n/a |
| MODMUL_R2 | 50000 / 25000 | - | - | - | - | - | - | - | 6.20 M | 4.59 M | 8.48 M | n/a |

#### p1024 (1024-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | GMP 224T | OpenSSL 224T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 12500 | - | - | - | - | - | - | - | 60.06 M | 22.03 M | 10.90 M | n/a |
| SUBTRACT | 50000 / 12500 | - | - | - | - | - | - | - | 53.82 M | 10.36 M | 14.64 M | n/a |
| ADDMOD | 50000 / 12500 | - | - | - | - | - | - | - | 17.75 M | 12.96 M | 6.60 M | n/a |
| SUBTRACTMOD | 50000 / 12500 | - | - | - | - | - | - | - | 22.77 M | 6.74 M | 12.93 M | n/a |
| MULTIPLYOPERANDSCANNING | 50000 / 12500 | - | - | - | - | - | - | - | 6.36 M | 6.51 M | 3.46 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 12500 | - | - | - | - | - | - | - | 6.36 M | 6.46 M | 6.47 M | n/a |
| MONTGOMERYMULTIPLICATION | 50000 / 12500 | - | - | - | - | - | - | - | 943.28 k | 2.13 M | 1.83 M | n/a |
| COMPARE | 50000 / 12500 | - | - | - | - | - | - | - | 143.66 M | 38.52 M | 39.41 M | n/a |
| REDUCE | 50000 / 1562 | - | - | - | - | - | - | - | 53.36 M | 19.48 M | 19.23 M | n/a |
| MODMUL | 50000 / 781 | - | - | - | - | - | - | - | 2.20 M | 913.60 k | 1.28 M | n/a |
| MODEXP | 50000 / 195 | - | - | - | - | - | - | - | 3.36 k | 31.71 k | 23.26 k | n/a |
| EXPONENTIATION | 50000 / 195 | - | - | - | - | - | - | - | 24.65 k | 32.99 k | 15.94 k | n/a |
| DIVIDE | 50000 / 1562 | - | - | - | - | - | - | - | 32.43 M | 13.95 M | 7.15 M | n/a |
| ISQRT | 50000 / 390 | - | - | - | - | - | - | - | 5.48 M | 2.75 M | n/a | n/a |
| MODMUL_R2 | 50000 / 12500 | - | - | - | - | - | - | - | 2.17 M | 2.18 M | 1.67 M | n/a |

#### p2048 (2048-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | GMP 224T | OpenSSL 224T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 6250 | - | - | - | - | - | - | - | 43.30 M | 18.61 M | 18.59 M | n/a |
| SUBTRACT | 50000 / 6250 | - | - | - | - | - | - | - | 43.55 M | 20.17 M | 20.74 M | n/a |
| ADDMOD | 50000 / 6250 | - | - | - | - | - | - | - | 15.23 M | 5.76 M | 10.55 M | n/a |
| SUBTRACTMOD | 50000 / 6250 | - | - | - | - | - | - | - | 17.55 M | 11.97 M | 6.15 M | n/a |
| MULTIPLYOPERANDSCANNING | 50000 / 6250 | - | - | - | - | - | - | - | 1.95 M | 1.13 M | 2.20 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 6250 | - | - | - | - | - | - | - | 1.94 M | 1.13 M | 2.18 M | n/a |
| MONTGOMERYMULTIPLICATION | 50000 / 6250 | - | - | - | - | - | - | - | 286.85 k | 1.06 M | 1.08 M | n/a |
| COMPARE | 50000 / 6250 | - | - | - | - | - | - | - | 143.61 M | 32.13 M | 16.82 M | n/a |
| REDUCE | 50000 / 781 | - | - | - | - | - | - | - | 37.45 M | 12.95 M | 12.90 M | n/a |
| MODMUL | 50000 / 390 | - | - | - | - | - | - | - | 682.53 k | 25.95 M | 9.08 M | n/a |
| MODEXP | 50000 / 97 | - | - | - | - | - | - | - | 466.0 | 11.75 k | 15.84 k | n/a |
| EXPONENTIATION | 50000 / 97 | - | - | - | - | - | - | - | 3.94 k | 100.28 k | 2.07 k | n/a |
| DIVIDE | 50000 / 781 | - | - | - | - | - | - | - | 24.11 M | 10.26 M | 10.23 M | n/a |
| ISQRT | 50000 / 195 | - | - | - | - | - | - | - | 3.23 M | 1.79 M | n/a | n/a |
| MODMUL_R2 | 50000 / 6250 | - | - | - | - | - | - | - | 661.60 k | 1.08 M | 1.06 M | n/a |

## 5. Head to head

Best OpenCL GPU result against best OpenCL CPU result and the CPU libraries.
Ratios above 1.00x mean the GPU is faster than that baseline.

### secp256k1 (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GMP 224T | OpenSSL | CGBN | GPU vs CPU-CL | GPU vs GMP 224T | GPU vs OpenSSL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | w32-il64 | 5.45 G | none | n/a | 71.27 M | 9.06 M | 17.57 M | n/a | n/a | 601.23x | 310.14x | n/a |
| SUBTRACT | w32-il64 | 5.62 G | none | n/a | 92.79 M | 43.98 M | 43.93 M | n/a | n/a | 127.70x | 127.85x | n/a |
| ADDMOD | w32-il64 | 5.96 G | none | n/a | 25.78 M | 17.59 M | 8.94 M | n/a | n/a | 338.75x | 666.68x | n/a |
| SUBTRACTMOD | w32-il64 | 5.92 G | none | n/a | 33.20 M | 25.80 M | 25.80 M | n/a | n/a | 229.35x | 229.35x | n/a |
| MULTIPLYOPERANDSCANNING | w32-il64 | 5.30 G | none | n/a | 57.25 M | 34.10 M | 34.54 M | n/a | n/a | 155.45x | 153.46x | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 3.09 G | none | n/a | 57.05 M | 34.94 M | 34.85 M | n/a | n/a | 88.46x | 88.69x | n/a |
| MONTGOMERYMULTIPLICATION | w32-il64 | 4.95 G | none | n/a | 7.12 M | 8.68 M | 4.53 M | n/a | n/a | 570.23x | 1092.19x | n/a |
| COMPARE | w32-il | 5.52 G | none | n/a | 154.78 M | 27.76 M | 52.55 M | n/a | n/a | 198.91x | 105.08x | n/a |
| REDUCE | w32-il64 | 190.43 M | none | n/a | 67.25 M | 33.26 M | 16.63 M | n/a | n/a | 5.73x | 11.45x | n/a |
| MODMUL | w32-o64 | 38.34 M | none | n/a | 12.43 M | 3.67 M | 7.25 M | n/a | n/a | 10.45x | 5.29x | n/a |
| MODEXP | w32-il64 | 656.92 k | none | n/a | 117.10 k | 4.23 M | 3.99 M | n/a | n/a | 0.16x | 0.16x | n/a |
| EXPONENTIATION | w32-o64 | 2.06 M | none | n/a | 375.20 k | 13.53 M | 1.60 M | n/a | n/a | 0.15x | 1.28x | n/a |
| DIVIDE | w32-il64 | 66.30 M | none | n/a | 40.53 M | 788.51 M | 417.10 M | n/a | n/a | 0.08x | 0.16x | n/a |
| ISQRT | w32-il64 | 2.64 M | none | n/a | 20.47 M | 639.33 M | n/a | n/a | n/a | 0.00x | n/a | n/a |
| MODMUL_R2 | w32-il64 | 3.99 G | none | n/a | 12.30 M | 375.30 M | 203.21 M | n/a | n/a | 10.64x | 19.65x | n/a |

### rsa256(composite) (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GMP 224T | OpenSSL | CGBN | GPU vs CPU-CL | GPU vs GMP 224T | GPU vs OpenSSL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | w32-il64 | 5.44 G | none | n/a | 69.96 M | 18.00 M | 33.18 M | n/a | n/a | 302.38x | 164.08x | n/a |
| SUBTRACT | w32-il | 5.60 G | none | n/a | 94.88 M | 44.13 M | 43.49 M | n/a | n/a | 126.87x | 128.71x | n/a |
| ADDMOD | w32-il64 | 5.88 G | none | n/a | 30.38 M | 886.45 M | 25.70 M | n/a | n/a | 6.63x | 228.85x | n/a |
| SUBTRACTMOD | w32-il64 | 6.01 G | none | n/a | 33.13 M | 26.66 M | 25.88 M | n/a | n/a | 225.64x | 232.42x | n/a |
| MULTIPLYOPERANDSCANNING | w32-il64 | 5.39 G | none | n/a | 55.94 M | 21.24 M | 34.65 M | n/a | n/a | 253.80x | 155.58x | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 3.08 G | none | n/a | 56.99 M | 34.77 M | 30.68 M | n/a | n/a | 88.48x | 100.28x | n/a |
| MONTGOMERYMULTIPLICATION | w32-il64 | 4.86 G | none | n/a | 7.10 M | 8.64 M | 7.10 M | n/a | n/a | 562.91x | 684.33x | n/a |
| COMPARE | w32-il64 | 5.55 G | none | n/a | 154.71 M | 52.82 M | 27.41 M | n/a | n/a | 105.16x | 202.61x | n/a |
| REDUCE | w32-il64 | 188.66 M | none | n/a | 42.86 M | 25.28 M | 24.76 M | n/a | n/a | 7.46x | 7.62x | n/a |
| MODMUL | w32-o64 | 38.21 M | none | n/a | 12.37 M | 9.31 M | 9.31 M | n/a | n/a | 4.10x | 4.10x | n/a |
| MODEXP | w32-il64 | 654.45 k | none | n/a | 124.84 k | 71.29 k | 70.50 k | n/a | n/a | 9.18x | 9.28x | n/a |
| EXPONENTIATION | w32-o64 | 2.05 M | none | n/a | 379.01 k | 18.65 M | 1.78 M | n/a | n/a | 0.11x | 1.15x | n/a |
| DIVIDE | w32-il64 | 63.99 M | none | n/a | 40.02 M | 764.45 M | 315.92 M | n/a | n/a | 0.08x | 0.20x | n/a |
| ISQRT | w32-o64 | 2.43 M | none | n/a | 20.63 M | 481.54 M | n/a | n/a | n/a | 0.01x | n/a | n/a |
| MODMUL_R2 | w32-il64 | 3.92 G | none | n/a | 12.36 M | 656.47 M | 211.34 M | n/a | n/a | 5.97x | 18.56x | n/a |

### brainpoolP512r1 (512-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GMP 224T | OpenSSL | CGBN | GPU vs CPU-CL | GPU vs GMP 224T | GPU vs OpenSSL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | w32-il64 | 1.75 G | none | n/a | 67.14 M | 1.47 G | 1.34 G | n/a | n/a | 1.19x | 1.31x | n/a |
| SUBTRACT | w32-il64 | 1.76 G | none | n/a | 91.18 M | 1.34 G | 39.52 M | n/a | n/a | 1.32x | 44.64x | n/a |
| ADDMOD | w32-il64 | 2.14 G | none | n/a | 27.38 M | 17.26 M | 17.26 M | n/a | n/a | 123.86x | 123.82x | n/a |
| SUBTRACTMOD | w32-il64 | 2.12 G | none | n/a | 29.28 M | 22.10 M | 21.39 M | n/a | n/a | 96.04x | 99.23x | n/a |
| MULTIPLYOPERANDSCANNING | w32-il64 | 799.06 M | none | n/a | 23.23 M | 17.48 M | 16.50 M | n/a | n/a | 45.71x | 48.44x | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 452.78 M | none | n/a | 23.28 M | 17.56 M | 688.61 M | n/a | n/a | 25.78x | 0.66x | n/a |
| MONTGOMERYMULTIPLICATION | w32-il64 | 1.28 G | none | n/a | 2.99 M | 6.76 M | 4.37 M | n/a | n/a | 188.93x | 292.33x | n/a |
| COMPARE | w32-il64 | 2.33 G | none | n/a | 152.77 M | 48.16 M | 48.75 M | n/a | n/a | 48.31x | 47.72x | n/a |
| REDUCE | w32-opt | 34.75 M | none | n/a | 41.86 M | 23.64 M | 11.95 M | n/a | n/a | 1.47x | 2.91x | n/a |
| MODMUL | w32-o64 | 5.87 M | none | n/a | 6.26 M | 4.15 M | 2.10 M | n/a | n/a | 1.42x | 2.80x | n/a |
| MODEXP | w32-il64 | 45.52 k | none | n/a | 20.37 k | 64.70 k | 41.21 k | n/a | n/a | 0.70x | 1.10x | n/a |
| EXPONENTIATION | w32-o64 | 27.28 k | none | n/a | 110.51 k | 136.28 k | 66.18 k | n/a | n/a | 0.20x | 0.41x | n/a |
| DIVIDE | w32-il64 | 9.67 M | none | n/a | 37.51 M | 11.91 M | 15.41 M | n/a | n/a | 0.81x | 0.63x | n/a |
| ISQRT | w32-o64 | 255.28 k | none | n/a | 11.27 M | 5.72 M | n/a | n/a | n/a | 0.04x | n/a | n/a |
| MODMUL_R2 | w32-il64 | 841.86 M | none | n/a | 6.20 M | 4.59 M | 8.48 M | n/a | n/a | 183.61x | 99.25x | n/a |

### p1024 (1024-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GMP 224T | OpenSSL | CGBN | GPU vs CPU-CL | GPU vs GMP 224T | GPU vs OpenSSL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | w32-il64 | 578.38 M | none | n/a | 60.06 M | 22.03 M | 10.90 M | n/a | n/a | 26.26x | 53.07x | n/a |
| SUBTRACT | w32-il64 | 581.84 M | none | n/a | 53.82 M | 10.36 M | 14.64 M | n/a | n/a | 56.18x | 39.73x | n/a |
| ADDMOD | w32-il64 | 765.21 M | none | n/a | 17.75 M | 12.96 M | 6.60 M | n/a | n/a | 59.03x | 115.92x | n/a |
| SUBTRACTMOD | w32-il | 754.54 M | none | n/a | 22.77 M | 6.74 M | 12.93 M | n/a | n/a | 111.97x | 58.35x | n/a |
| MULTIPLYOPERANDSCANNING | w32-il64 | 140.36 M | none | n/a | 6.36 M | 6.51 M | 3.46 M | n/a | n/a | 21.56x | 40.61x | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 69.38 M | none | n/a | 6.36 M | 6.46 M | 6.47 M | n/a | n/a | 10.75x | 10.73x | n/a |
| MONTGOMERYMULTIPLICATION | w32-il64 | 233.12 M | none | n/a | 943.28 k | 2.13 M | 1.83 M | n/a | n/a | 109.25x | 127.62x | n/a |
| COMPARE | w32-il | 895.98 M | none | n/a | 143.66 M | 38.52 M | 39.41 M | n/a | n/a | 23.26x | 22.73x | n/a |
| REDUCE | w32-il64 | 5.66 M | none | n/a | 53.36 M | 19.48 M | 19.23 M | n/a | n/a | 0.29x | 0.29x | n/a |
| MODMUL | w32-il64 | 725.27 k | none | n/a | 2.20 M | 913.60 k | 1.28 M | n/a | n/a | 0.79x | 0.57x | n/a |
| MODEXP | w32-o64 | 2.39 k | none | n/a | 3.36 k | 31.71 k | 23.26 k | n/a | n/a | 0.08x | 0.10x | n/a |
| EXPONENTIATION | w32-il64 | 1.88 k | none | n/a | 24.65 k | 32.99 k | 15.94 k | n/a | n/a | 0.06x | 0.12x | n/a |
| DIVIDE | w32-il | 1.06 M | none | n/a | 32.43 M | 13.95 M | 7.15 M | n/a | n/a | 0.08x | 0.15x | n/a |
| ISQRT | w32-il64 | 17.89 k | none | n/a | 5.48 M | 2.75 M | n/a | n/a | n/a | 0.01x | n/a | n/a |
| MODMUL_R2 | w32-il64 | 125.40 M | none | n/a | 2.17 M | 2.18 M | 1.67 M | n/a | n/a | 57.51x | 75.00x | n/a |

### p2048 (2048-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GMP 224T | OpenSSL | CGBN | GPU vs CPU-CL | GPU vs GMP 224T | GPU vs OpenSSL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | w32-il | 91.68 M | none | n/a | 43.30 M | 18.61 M | 18.59 M | n/a | n/a | 4.93x | 4.93x | n/a |
| SUBTRACT | w32-il | 93.16 M | none | n/a | 43.55 M | 20.17 M | 20.74 M | n/a | n/a | 4.62x | 4.49x | n/a |
| ADDMOD | w32-il | 159.48 M | none | n/a | 15.23 M | 5.76 M | 10.55 M | n/a | n/a | 27.69x | 15.12x | n/a |
| SUBTRACTMOD | w32-il | 166.33 M | none | n/a | 17.55 M | 11.97 M | 6.15 M | n/a | n/a | 13.89x | 27.06x | n/a |
| MULTIPLYOPERANDSCANNING | w32-il64 | 19.72 M | none | n/a | 1.95 M | 1.13 M | 2.20 M | n/a | n/a | 17.39x | 8.97x | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 7.12 M | none | n/a | 1.94 M | 1.13 M | 2.18 M | n/a | n/a | 6.32x | 3.26x | n/a |
| MONTGOMERYMULTIPLICATION | w32 | 27.54 M | none | n/a | 286.85 k | 1.06 M | 1.08 M | n/a | n/a | 25.97x | 25.47x | n/a |
| COMPARE | w32-il | 288.12 M | none | n/a | 143.61 M | 32.13 M | 16.82 M | n/a | n/a | 8.97x | 17.13x | n/a |
| REDUCE | w32-il64 | 954.43 k | none | n/a | 37.45 M | 12.95 M | 12.90 M | n/a | n/a | 0.07x | 0.07x | n/a |
| MODMUL | w32-il64 | 107.91 k | none | n/a | 682.53 k | 25.95 M | 9.08 M | n/a | n/a | 0.00x | 0.01x | n/a |
| MODEXP | w32-opt | 79.0 | none | n/a | 466.0 | 11.75 k | 15.84 k | n/a | n/a | 0.01x | 0.00x | n/a |
| EXPONENTIATION | w32-o64 | 120.1 | none | n/a | 3.94 k | 100.28 k | 2.07 k | n/a | n/a | 0.00x | 0.06x | n/a |
| DIVIDE | w32-il64 | 69.65 k | none | n/a | 24.11 M | 10.26 M | 10.23 M | n/a | n/a | 0.01x | 0.01x | n/a |
| ISQRT | w32-il64 | 3.74 k | none | n/a | 3.23 M | 1.79 M | n/a | n/a | n/a | 0.00x | n/a | n/a |
| MODMUL_R2 | w32-il64 | 14.06 M | none | n/a | 661.60 k | 1.08 M | 1.06 M | n/a | n/a | 12.96x | 13.20x | n/a |

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

Also written to `NVIDIA_H100_80GB_HBM3_Report.csv` for analysis.

```csv
kind,device,device_type,kernel,modulus,bits,operation,items,seconds,ops_per_sec,mismatches
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,secp256k1,256,ADD,50000,0.000701555,71270244.488,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,secp256k1,256,ADD,50000,0.005516557,9063624.248,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,secp256k1,256,ADD,50000,0.002845631,17570795.487,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,ADD,50000,0.000027992,1786235400.586,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,ADD,50000,0.000352694,141766061.970,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,ADD,50000,0.000019109,2616585008.285,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,ADD,50000,0.000325864,153438158.985,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,secp256k1,256,ADD,50000,0.000012930,3866831691.155,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,secp256k1,256,ADD,50000,0.000342667,145914211.167,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,ADD,50000,0.000012770,3915336289.382,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,ADD,50000,0.000333373,149982375.487,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,ADD,50000,0.000013001,3845780171.920,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,ADD,50000,0.000339122,147439351.880,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,ADD,50000,0.000009406,5315553584.158,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,ADD,50000,0.000364812,137056860.873,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,ADD,50000,0.000009175,5449359642.712,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,ADD,50000,0.000336368,148646880.711,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,secp256k1,256,SUBTRACT,50000,0.000538841,92791834.291,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,secp256k1,256,SUBTRACT,50000,0.001136992,43975678.306,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,secp256k1,256,SUBTRACT,50000,0.001138299,43925206.105,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,SUBTRACT,50000,0.000026299,1901235611.587,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,SUBTRACT,50000,0.000348896,143309267.959,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,SUBTRACT,50000,0.000018669,2678194712.162,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,SUBTRACT,50000,0.000346694,144219339.172,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,secp256k1,256,SUBTRACT,50000,0.000013392,3733455577.191,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,secp256k1,256,SUBTRACT,50000,0.000346014,144502708.260,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,SUBTRACT,50000,0.000012910,3872968633.675,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,SUBTRACT,50000,0.000332654,150306540.047,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,SUBTRACT,50000,0.000012770,3915336289.382,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,SUBTRACT,50000,0.000337830,148003515.446,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,SUBTRACT,50000,0.000008954,5584261618.473,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,SUBTRACT,50000,0.000367224,136156597.955,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,SUBTRACT,50000,0.000008903,5615804518.828,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,SUBTRACT,50000,0.000345383,144766891.374,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,secp256k1,256,ADDMOD,50000,0.001939426,25780821.468,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,secp256k1,256,ADDMOD,50000,0.002843128,17586266.704,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,secp256k1,256,ADDMOD,50000,0.005595461,8935813.482,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,ADDMOD,50000,0.000038778,1289377280.369,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,ADDMOD,50000,0.000362778,137825305.368,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,ADDMOD,50000,0.000024356,2052886631.998,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,ADDMOD,50000,0.000344422,145170870.153,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,secp256k1,256,ADDMOD,50000,0.000013480,3709209009.258,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,secp256k1,256,ADDMOD,50000,0.000344761,145028124.392,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,ADDMOD,50000,0.000010407,4804643923.394,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,ADDMOD,50000,0.000332823,150229991.661,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,ADDMOD,50000,0.000010626,4705266538.124,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,ADDMOD,50000,0.000342900,145815134.742,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,ADDMOD,50000,0.000008492,5888033691.599,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,ADDMOD,50000,0.000351161,142384927.438,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,ADDMOD,50000,0.000008393,5957289303.151,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,ADDMOD,50000,0.000339424,147308278.137,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,secp256k1,256,SUBTRACTMOD,50000,0.001505843,33203991.539,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,secp256k1,256,SUBTRACTMOD,50000,0.001938226,25796789.281,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,secp256k1,256,SUBTRACTMOD,50000,0.001938213,25796962.818,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.000038067,1313477790.282,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.000363547,137533664.655,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,SUBTRACTMOD,50000,0.000024127,2072380575.928,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,SUBTRACTMOD,50000,0.000345761,144608577.324,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,secp256k1,256,SUBTRACTMOD,50000,0.000013771,3630940835.926,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,secp256k1,256,SUBTRACTMOD,50000,0.000344820,145003055.249,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,SUBTRACTMOD,50000,0.000010297,4855923589.001,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,SUBTRACTMOD,50000,0.000333413,149963941.899,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,SUBTRACTMOD,50000,0.000010535,4746029985.856,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,SUBTRACTMOD,50000,0.000355948,140470049.922,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,SUBTRACTMOD,50000,0.000008471,5902274758.135,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,SUBTRACTMOD,50000,0.000342056,146174829.013,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,SUBTRACTMOD,50000,0.000008451,5916584879.877,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,SUBTRACTMOD,50000,0.000351662,142182055.859,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000873420,57246232.596,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001466269,34100162.475,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001447533,34541531.832,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000912460,54796948.597,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001368828,36527595.648,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000256198,195161915.010,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000654018,76450547.529,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000074223,673648504.316,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000464270,107695977.982,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000019519,2561651455.292,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000412764,121134586.938,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000019770,2529069681.553,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000425743,117441759.819,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000010507,4758650168.410,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000430891,116038548.759,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000009432,5300858135.861,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000414617,120593117.576,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000876494,57045441.935,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001430907,34942864.897,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001434528,34854663.435,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000169205,295500331.348,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000568178,88000660.901,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000056265,888653146.622,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000444362,112520887.808,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000020761,2408356863.449,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000417631,119722879.035,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000021635,2311110253.982,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000414386,120660332.986,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000020372,2454379226.479,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000423601,118035632.907,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000017416,2870967443.850,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000423631,118027329.127,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000016175,3091149884.846,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000425503,117508079.146,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.007022494,7119977.147,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.005762137,8677335.492,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.011036428,4530451.465,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000078687,635425389.987,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000405913,123179038.468,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000025989,1923854769.584,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000346363,144357391.155,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000011686,4278537711.189,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000346012,144503486.144,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000013491,3706136352.340,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000333864,149761471.084,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000012368,4042702650.602,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000339855,147121779.259,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000012750,3921628283.419,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000368237,135781935.902,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000010105,4948119004.608,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000340834,146698867.661,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,secp256k1,256,COMPARE,50000,0.000323033,154782818.147,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,secp256k1,256,COMPARE,50000,0.001801032,27761866.501,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,secp256k1,256,COMPARE,50000,0.000951451,52551299.954,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,COMPARE,50000,0.000027841,1795915273.968,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,COMPARE,50000,0.000355558,140623847.367,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,COMPARE,50000,0.000018267,2737182175.997,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,COMPARE,50000,0.000345182,144851258.917,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,COMPARE,50000,0.000009574,5222479688.716,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,COMPARE,50000,0.000334254,149587049.390,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,COMPARE,50000,0.000009624,5195189781.304,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,COMPARE,50000,0.000335800,148898362.002,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,COMPARE,50000,0.000009054,5522227031.475,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,COMPARE,50000,0.000358291,139551381.813,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,COMPARE,50000,0.000009062,5517686659.815,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,COMPARE,50000,0.000346055,144485596.947,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,secp256k1,256,REDUCE,6250,0.000092933,67252425.511,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,secp256k1,256,REDUCE,6250,0.000187895,33263267.532,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,secp256k1,256,REDUCE,6250,0.000375761,16632903.386,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,REDUCE,50000,0.000138119,362007034.200,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,REDUCE,50000,0.000463180,107949336.465,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,REDUCE,50000,0.000091799,544670594.919,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,REDUCE,50000,0.000411103,121624154.883,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,REDUCE,50000,0.000035252,1418342259.326,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,REDUCE,50000,0.000358252,139566618.660,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,REDUCE,50000,0.000034181,1462783804.697,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,REDUCE,50000,0.000375008,133330414.094,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,REDUCE,50000,0.000033522,1491556681.669,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,REDUCE,50000,0.000388457,128714537.931,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,REDUCE,50000,0.000032820,1523470238.365,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,REDUCE,50000,0.000364421,137203972.461,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,secp256k1,256,MODMUL,3125,0.000251347,12433022.352,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,secp256k1,256,MODMUL,3125,0.000851649,3669350.489,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,secp256k1,256,MODMUL,3125,0.000430903,7252216.376,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,MODMUL,50000,0.000369430,135343788.319,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,MODMUL,50000,0.000700418,71385953.398,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,MODMUL,50000,0.000218252,229093268.927,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,MODMUL,50000,0.000560576,89193956.612,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,MODMUL,50000,0.000103407,483528092.802,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,MODMUL,50000,0.000422798,118259756.463,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,MODMUL,50000,0.000081513,613398510.123,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,MODMUL,50000,0.000411522,121500292.395,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,MODMUL,50000,0.000102164,489408113.183,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,MODMUL,50000,0.000462148,108190370.555,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,MODMUL,50000,0.000083005,602372946.166,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,MODMUL,50000,0.000423711,118005018.507,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,secp256k1,256,MODEXP,781,0.006669369,117102.540,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,secp256k1,256,MODEXP,781,0.000184698,4228523.707,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,secp256k1,256,MODEXP,781,0.000195645,3991928.313,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,MODEXP,50000,0.016818859,2972853.327,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,MODEXP,50000,0.017239286,2900352.199,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,MODEXP,50000,0.003179274,15726858.442,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,MODEXP,50000,0.003506878,14257694.287,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,MODEXP,50000,0.002537204,19706732.660,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,MODEXP,50000,0.002862936,17464589.496,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,MODEXP,50000,0.001211451,41272808.287,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,MODEXP,50000,0.001540648,32453880.350,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,MODEXP,50000,0.002504455,19964423.990,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,MODEXP,50000,0.002885360,17328858.551,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,MODEXP,50000,0.001188885,42056197.456,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,MODEXP,50000,0.001517475,32949479.555,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,secp256k1,256,EXPONENTIATION,781,0.002081572,375197.137,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,secp256k1,256,EXPONENTIATION,781,0.000057710,13533104.679,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,secp256k1,256,EXPONENTIATION,781,0.000486909,1603995.484,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,EXPONENTIATION,50000,0.026735194,1870194.020,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,EXPONENTIATION,50000,0.027152816,1841429.644,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,EXPONENTIATION,50000,0.007067621,7074516.551,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,EXPONENTIATION,50000,0.007399714,6757018.253,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,EXPONENTIATION,50000,0.000461439,108356761.513,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,EXPONENTIATION,50000,0.000775302,64491012.572,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,EXPONENTIATION,50000,0.000379233,131845174.092,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,EXPONENTIATION,50000,0.000736011,67933749.554,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,EXPONENTIATION,50000,0.000464002,107758232.741,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,EXPONENTIATION,50000,0.000826878,60468377.909,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,EXPONENTIATION,50000,0.000406606,122969127.099,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,EXPONENTIATION,50000,0.000747891,66854649.196,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,secp256k1,256,DIVIDE,6250,0.000154212,40528530.944,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,secp256k1,256,DIVIDE,6250,0.000007926,788507221.643,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,secp256k1,256,DIVIDE,6250,0.000014985,417097262.190,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,DIVIDE,50000,0.000328576,152171706.830,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,DIVIDE,50000,0.000739748,67590616.116,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,DIVIDE,50000,0.000276066,181115871.859,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,DIVIDE,50000,0.000652047,76681603.702,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,DIVIDE,50000,0.000099560,502208482.535,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,DIVIDE,50000,0.000479054,104372431.276,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,DIVIDE,50000,0.000097688,511832086.336,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,DIVIDE,50000,0.000491964,101633514.942,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,DIVIDE,50000,0.000098549,507362697.513,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,DIVIDE,50000,0.000502318,99538510.828,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,DIVIDE,50000,0.000094272,530379072.157,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,DIVIDE,50000,0.000465604,107387498.450,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,secp256k1,256,ISQRT,1562,0.000076306,20470332.486,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,secp256k1,256,ISQRT,1562,0.000002443,639333441.838,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,ISQRT,50000,0.003647253,13708947.839,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,ISQRT,50000,0.004011132,12465309.595,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,ISQRT,50000,0.002715284,18414280.021,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,ISQRT,50000,0.003041446,16439547.690,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,ISQRT,50000,0.000660527,75697075.727,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,ISQRT,50000,0.001024267,48815410.830,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,ISQRT,50000,0.000598164,83589077.559,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,ISQRT,50000,0.000965409,51791521.513,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,ISQRT,50000,0.000679366,73597999.627,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,ISQRT,50000,0.001034664,48324867.277,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,ISQRT,50000,0.000591835,84483005.234,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,ISQRT,50000,0.000924958,54056513.413,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,secp256k1,256,MODMUL_R2,50000,0.004064097,12302855.751,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,secp256k1,256,MODMUL_R2,50000,0.000133227,375300355.817,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,secp256k1,256,MODMUL_R2,50000,0.000246047,203213172.288,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,MODMUL_R2,50000,0.000094192,530830065.851,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,MODMUL_R2,50000,0.000429578,116393263.611,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,MODMUL_R2,50000,0.000027481,1819407997.831,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,MODMUL_R2,50000,0.000357730,139770095.025,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,MODMUL_R2,50000,0.000017809,2807608576.509,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,MODMUL_R2,50000,0.000340946,146650781.233,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,MODMUL_R2,50000,0.000013320,3753817032.583,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,MODMUL_R2,50000,0.000342648,145922143.098,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,MODMUL_R2,50000,0.000017727,2820589009.142,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,MODMUL_R2,50000,0.000373524,133860320.943,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,MODMUL_R2,50000,0.000012519,3993980895.700,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,MODMUL_R2,50000,0.000341406,146453157.002,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,ADD,50000,0.000714668,69962594.503,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,ADD,50000,0.002777292,18003145.842,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,ADD,50000,0.001507057,33177237.909,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,ADD,50000,0.000027381,1826091537.415,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,ADD,50000,0.000361776,138207075.226,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,ADD,50000,0.000018509,2701373211.231,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,ADD,50000,0.000345722,144624938.580,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,rsa256(composite),256,ADD,50000,0.000012698,3937735895.555,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,rsa256(composite),256,ADD,50000,0.000340864,146686041.530,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,ADD,50000,0.000013009,3843577548.683,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,ADD,50000,0.000339173,147417490.060,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,ADD,50000,0.000012558,3981540433.106,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,ADD,50000,0.000349807,142936116.421,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,ADD,50000,0.000009265,5396772336.148,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,ADD,50000,0.000356609,140209584.598,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,ADD,50000,0.000009185,5443834029.609,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,ADD,50000,0.000346795,144177510.420,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACT,50000,0.000526974,94881278.339,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACT,50000,0.001133081,44127488.583,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACT,50000,0.001149559,43494925.344,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,SUBTRACT,50000,0.000029674,1684988111.230,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,SUBTRACT,50000,0.000358522,139461479.634,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,SUBTRACT,50000,0.000018207,2746142772.379,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,SUBTRACT,50000,0.000343438,145586584.374,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,rsa256(composite),256,SUBTRACT,50000,0.000012709,3934273134.985,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,rsa256(composite),256,SUBTRACT,50000,0.000348385,143519207.862,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,SUBTRACT,50000,0.000013119,3811379468.976,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,SUBTRACT,50000,0.000333333,149999975.413,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,SUBTRACT,50000,0.000013052,3830961267.304,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,SUBTRACT,50000,0.000342628,145930869.218,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,SUBTRACT,50000,0.000008931,5598236830.031,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,SUBTRACT,50000,0.000358341,139531796.469,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,SUBTRACT,50000,0.000009075,5509758949.097,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,SUBTRACT,50000,0.000348447,143493890.555,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,ADDMOD,50000,0.001645657,30383003.460,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,ADDMOD,50000,0.000056405,886452202.629,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,ADDMOD,50000,0.001945884,25695258.747,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,ADDMOD,50000,0.000035824,1395702469.714,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,ADDMOD,50000,0.000373553,133849641.486,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,ADDMOD,50000,0.000022434,2228789903.686,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,ADDMOD,50000,0.000347637,143828337.513,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,rsa256(composite),256,ADDMOD,50000,0.000013581,3681737155.397,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,rsa256(composite),256,ADDMOD,50000,0.000337221,148270838.028,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,ADDMOD,50000,0.000010226,4889534717.668,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,ADDMOD,50000,0.000304470,164219878.748,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,ADDMOD,50000,0.000010386,4814122238.164,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,ADDMOD,50000,0.000350779,142539921.305,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,ADDMOD,50000,0.000008643,5785246896.552,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,ADDMOD,50000,0.000358053,139644305.720,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,ADDMOD,50000,0.000008503,5880294764.513,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,ADDMOD,50000,0.000341877,146251283.616,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACTMOD,50000,0.001509371,33126383.909,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACTMOD,50000,0.001875725,26656364.608,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACTMOD,50000,0.001932140,25878036.074,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,SUBTRACTMOD,50000,0.000039380,1269678630.215,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,SUBTRACTMOD,50000,0.000363769,137449861.493,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,SUBTRACTMOD,50000,0.000024049,2079122112.927,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,SUBTRACTMOD,50000,0.000347747,143782884.384,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,rsa256(composite),256,SUBTRACTMOD,50000,0.000014061,3555907484.435,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,rsa256(composite),256,SUBTRACTMOD,50000,0.000352753,141742107.792,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,SUBTRACTMOD,50000,0.000010056,4971947694.017,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,SUBTRACTMOD,50000,0.000333354,149990755.889,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,SUBTRACTMOD,50000,0.000010576,4727641000.352,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,SUBTRACTMOD,50000,0.000347804,143759013.742,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,SUBTRACTMOD,50000,0.000008633,5791487723.840,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,SUBTRACTMOD,50000,0.000356490,140256470.330,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,SUBTRACTMOD,50000,0.000008313,6014686444.096,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,SUBTRACTMOD,50000,0.000339853,147122585.596,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000893826,55939320.106,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.002353785,21242385.306,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001442848,34653679.192,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000919761,54361940.327,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001353966,36928547.197,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000259912,192373068.461,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000651676,76725219.370,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000074433,671743590.000,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000468338,106760522.914,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000019468,2568268809.797,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000406556,122984338.520,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000019100,2617860893.310,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000421166,118717916.785,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000010706,4670067084.203,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000420395,118935682.195,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000009274,5391352801.767,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000413025,121058106.530,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000877336,56990729.827,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001438089,34768369.625,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001629923,30676286.900,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000169707,294624640.274,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000591423,84541807.394,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000055213,905591579.516,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000449279,111289341.432,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000020821,2401462300.948,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000427946,116837048.644,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000021623,2312304729.089,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000406304,123060452.109,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000019960,2504996789.847,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000418141,119576750.650,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000016894,2959597089.305,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000434946,114956728.192,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000016253,3076271556.269,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000417029,119895599.198,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.007042082,7100172.905,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.005789800,8635877.460,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.007038737,7103547.415,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000078509,636872656.528,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000419773,119111950.445,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000025649,1949422338.417,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000354795,140926535.734,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000011427,4375476055.420,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000344321,145213277.362,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000013530,3695422026.432,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000334134,149640417.645,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000012217,4092627778.625,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000344982,144934942.309,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000012828,3897712443.735,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000357229,139966137.257,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000010286,4861199855.125,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000341376,146465942.436,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,COMPARE,50000,0.000323190,154707736.286,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,COMPARE,50000,0.000946697,52815200.269,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,COMPARE,50000,0.001823958,27412907.112,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,COMPARE,50000,0.000027662,1807524449.532,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,COMPARE,50000,0.000362998,137741853.317,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,COMPARE,50000,0.000017857,2799994325.649,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,COMPARE,50000,0.000343651,145496626.485,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,COMPARE,50000,0.000009354,5345190282.756,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,COMPARE,50000,0.000333626,149868494.925,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,COMPARE,50000,0.000009654,5179152151.264,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,COMPARE,50000,0.000337189,148284761.942,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,COMPARE,50000,0.000009594,5211327043.293,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,COMPARE,50000,0.000355106,140803088.447,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,COMPARE,50000,0.000009002,5554220070.350,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,COMPARE,50000,0.000351651,142186574.572,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,REDUCE,6250,0.000145829,42858392.782,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,REDUCE,6250,0.000247238,25279322.663,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,REDUCE,6250,0.000252471,24755348.040,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,REDUCE,50000,0.000139972,357213802.281,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,REDUCE,50000,0.000484923,103109174.506,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,REDUCE,50000,0.000091819,544549053.657,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,REDUCE,50000,0.000417532,119751185.978,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,REDUCE,50000,0.000035573,1405568415.541,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,REDUCE,50000,0.000358222,139578229.920,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,REDUCE,50000,0.000033682,1484463064.757,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,REDUCE,50000,0.000367405,136089641.012,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,REDUCE,50000,0.000033921,1474029191.148,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,REDUCE,50000,0.000380825,131293814.746,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,REDUCE,50000,0.000033129,1509251411.222,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,REDUCE,50000,0.000373164,133989276.284,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MODMUL,3125,0.000252714,12365770.882,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MODMUL,3125,0.000335582,9312185.034,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MODMUL,3125,0.000335596,9311792.836,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,MODMUL,50000,0.000370871,134817666.616,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,MODMUL,50000,0.000708390,70582585.968,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,MODMUL,50000,0.000218211,229136290.770,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,MODMUL,50000,0.000548137,91218012.838,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,MODMUL,50000,0.000102986,485504532.465,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,MODMUL,50000,0.000426965,117105662.995,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,MODMUL,50000,0.000081783,611372801.603,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,MODMUL,50000,0.000416219,120128998.998,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,MODMUL,50000,0.000102194,489265389.593,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,MODMUL,50000,0.000455298,109818298.450,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,MODMUL,50000,0.000082847,603524115.293,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,MODMUL,50000,0.000416590,120022112.629,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MODEXP,781,0.006256234,124835.486,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MODEXP,781,0.010955406,71289.003,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MODEXP,781,0.011077341,70504.286,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,MODEXP,50000,0.016841531,2968851.264,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,MODEXP,50000,0.017274316,2894470.545,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,MODEXP,50000,0.003177362,15736326.827,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,MODEXP,50000,0.003565267,14024196.171,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,MODEXP,50000,0.002531124,19754067.557,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,MODEXP,50000,0.002930589,17061416.546,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,MODEXP,50000,0.001209527,41338464.950,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,MODEXP,50000,0.001541849,32428592.345,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,MODEXP,50000,0.002499148,20006816.325,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,MODEXP,50000,0.002910359,17180012.122,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,MODEXP,50000,0.001193373,41898065.046,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,MODEXP,50000,0.001523182,32826022.559,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,EXPONENTIATION,781,0.002060631,379010.065,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,EXPONENTIATION,781,0.000041877,18650029.606,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,EXPONENTIATION,781,0.000438052,1782893.684,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,EXPONENTIATION,50000,0.026757328,1868646.985,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,EXPONENTIATION,50000,0.027131684,1842863.855,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,EXPONENTIATION,50000,0.007086249,7055919.079,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,EXPONENTIATION,50000,0.007412903,6744995.844,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,EXPONENTIATION,50000,0.000460686,108533757.060,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,EXPONENTIATION,50000,0.000785457,63657210.883,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,EXPONENTIATION,50000,0.000381226,131155893.662,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,EXPONENTIATION,50000,0.000712927,70133365.382,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,EXPONENTIATION,50000,0.000463992,107760395.658,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,EXPONENTIATION,50000,0.000841221,59437424.246,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,EXPONENTIATION,50000,0.000405502,123304083.564,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,EXPONENTIATION,50000,0.000750413,66629961.154,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,DIVIDE,6250,0.000156172,40020075.378,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,DIVIDE,6250,0.000008176,764453948.967,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,DIVIDE,6250,0.000019783,315922050.424,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,DIVIDE,50000,0.000330659,151213352.786,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,DIVIDE,50000,0.000685066,72985670.161,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,DIVIDE,50000,0.000277149,180408659.009,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,DIVIDE,50000,0.000664756,75215601.533,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,DIVIDE,50000,0.000103557,482823634.369,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,DIVIDE,50000,0.000481127,103922701.633,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,DIVIDE,50000,0.000099409,502970687.652,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,DIVIDE,50000,0.000500185,99962930.717,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,DIVIDE,50000,0.000101894,490705352.442,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,DIVIDE,50000,0.000497742,100453724.417,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,DIVIDE,50000,0.000097666,511949223.786,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,DIVIDE,50000,0.000500707,99858808.474,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,ISQRT,1562,0.000075722,20627969.112,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,ISQRT,1562,0.000003244,481538468.943,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,ISQRT,50000,0.003646852,13710453.248,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,ISQRT,50000,0.003995828,12513050.077,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,ISQRT,50000,0.002786212,17945513.813,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,ISQRT,50000,0.003114147,16055760.440,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,ISQRT,50000,0.000686187,72866402.819,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,ISQRT,50000,0.001042176,47976539.640,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,ISQRT,50000,0.000643572,77691404.425,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,ISQRT,50000,0.000985852,50717577.063,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,ISQRT,50000,0.000702461,71178305.682,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,ISQRT,50000,0.001080342,46281655.989,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,ISQRT,50000,0.000649272,77009382.773,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,ISQRT,50000,0.000976916,51181452.034,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MODMUL_R2,50000,0.004044668,12361954.616,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MODMUL_R2,50000,0.000076164,656473890.023,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MODMUL_R2,50000,0.000236588,211338211.421,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,MODMUL_R2,50000,0.000093842,532810892.995,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,MODMUL_R2,50000,0.000382477,130726672.218,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,MODMUL_R2,50000,0.000027902,1791958985.314,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,MODMUL_R2,50000,0.000357540,139844366.069,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,MODMUL_R2,50000,0.000018286,2734393969.644,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,MODMUL_R2,50000,0.000340704,146755008.365,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,MODMUL_R2,50000,0.000013130,3808135281.600,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,MODMUL_R2,50000,0.000342427,146016599.307,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,MODMUL_R2,50000,0.000017926,2789229592.685,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,MODMUL_R2,50000,0.000366474,136435486.839,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,MODMUL_R2,50000,0.000012748,3922201285.798,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,MODMUL_R2,50000,0.000354005,141240933.414,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,ADD,25000,0.000372340,67142963.033,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,ADD,25000,0.000016973,1472954606.882,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,ADD,25000,0.000018662,1339593777.714,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,ADD,50000,0.000060391,827942310.777,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,ADD,50000,0.000548620,91137800.895,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,ADD,50000,0.000033431,1495628794.295,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,ADD,50000,0.000608910,82113958.141,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,brainpoolP512r1,512,ADD,50000,0.000019789,2526689156.627,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,brainpoolP512r1,512,ADD,50000,0.000634279,78829650.599,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,ADD,50000,0.000020141,2482525256.636,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,ADD,50000,0.000588529,84957607.829,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,ADD,50000,0.000019940,2507570817.375,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,ADD,50000,0.000596011,83891061.032,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,ADD,50000,0.000014663,3410003252.033,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,ADD,50000,0.000635469,78682003.136,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,ADD,50000,0.000014283,3500723213.354,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,ADD,50000,0.000595281,83993959.742,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACT,25000,0.000274187,91178732.301,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACT,25000,0.000018672,1338934084.107,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACT,25000,0.000632528,39523928.548,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,SUBTRACT,50000,0.000058869,849344901.123,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,SUBTRACT,50000,0.000545233,91703831.648,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,SUBTRACT,50000,0.000032838,1522606103.233,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,SUBTRACT,50000,0.000601158,83172872.617,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,brainpoolP512r1,512,SUBTRACT,50000,0.000020390,2452137170.001,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,brainpoolP512r1,512,SUBTRACT,50000,0.000634627,78786385.022,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,SUBTRACT,50000,0.000020351,2456850228.812,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,SUBTRACT,50000,0.000591213,84571905.294,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,SUBTRACT,50000,0.000020420,2448558387.303,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,SUBTRACT,50000,0.000602720,82957218.395,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,SUBTRACT,50000,0.000014711,3398777614.586,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,SUBTRACT,50000,0.000584112,85599951.530,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,SUBTRACT,50000,0.000014171,3528331440.589,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,SUBTRACT,50000,0.000585154,85447635.539,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,ADDMOD,25000,0.000912964,27383325.717,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,ADDMOD,25000,0.001448604,17257991.278,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,ADDMOD,25000,0.001448174,17263124.399,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,ADDMOD,50000,0.000076126,656803171.030,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,ADDMOD,50000,0.000563633,88710254.529,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,ADDMOD,50000,0.000042073,1188398512.485,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,ADDMOD,50000,0.000611095,81820371.313,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,brainpoolP512r1,512,ADDMOD,50000,0.000021623,2312304729.089,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,brainpoolP512r1,512,ADDMOD,50000,0.000635000,78740164.206,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,ADDMOD,50000,0.000024417,2047718788.619,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,ADDMOD,50000,0.000584092,85602954.251,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,ADDMOD,50000,0.000024077,2076709391.923,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,ADDMOD,50000,0.000594558,84096057.343,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,ADDMOD,50000,0.000011958,4181237632.399,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,ADDMOD,50000,0.000594709,84074722.659,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,ADDMOD,50000,0.000011696,4275130689.600,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,ADDMOD,50000,0.000597874,83629702.693,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000853844,29279340.809,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACTMOD,25000,0.001131067,22103014.769,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACTMOD,25000,0.001168626,21392640.820,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000087773,569648485.877,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000563713,88697650.351,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000045599,1096505273.477,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000617793,80933279.867,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000024406,2048656460.353,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000640368,78080093.079,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000024546,2036996934.284,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000584774,85503158.485,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000023795,2101256015.656,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000596451,83829233.301,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000011878,4209431644.974,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000596901,83765928.247,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000011778,4245381243.081,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000590172,84721112.465,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.001076302,23227690.950,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.001430069,17481678.475,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.001515418,16497100.391,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.004782869,10453977.458,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.005439362,9192255.281,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.001293406,38657625.088,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.002019191,24762390.491,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000336459,148606557.976,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.001080101,46291951.886,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000060141,831378394.450,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000765536,65313716.502,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000061402,814304431.973,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000781631,63968795.739,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000031559,1584344307.384,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000746278,66999152.887,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000031287,1598115473.001,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000753848,66326378.549,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.001073910,23279409.766,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.001423446,17563006.882,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000036305,688613862.808,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.001135476,44034398.724,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.001731647,28874242.178,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000314126,159171903.110,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.001023125,48869888.546,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000088785,563159182.646,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000822565,60785500.361,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000089305,559882064.866,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000793830,62985812.701,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000089183,560642138.680,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000793589,63004883.407,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000059631,838493958.893,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000774058,64594677.659,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000055214,905561029.585,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000767820,65119464.364,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.008361347,2989948.861,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.003697136,6761990.140,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.005720517,4370234.623,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000269186,185745343.833,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000829784,60256630.609,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000081243,615437687.140,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000649972,76926403.994,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000026690,1873371875.218,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000642441,77828132.719,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000037577,1330601050.858,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000604805,82671328.172,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000025989,1923854769.584,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000605265,82608488.101,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000031298,1597544819.377,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000607368,82322467.630,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000019569,2555068113.459,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000605525,82572912.649,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,COMPARE,25000,0.000163647,152767603.722,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,COMPARE,25000,0.000519071,48162989.012,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,COMPARE,25000,0.000512797,48752263.039,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,COMPARE,50000,0.000060281,829451707.196,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,COMPARE,50000,0.000530841,94190192.742,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,COMPARE,50000,0.000033133,1509081718.012,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,COMPARE,50000,0.000602970,82922878.934,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,COMPARE,50000,0.000018889,2647031417.020,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,COMPARE,50000,0.000594178,84149837.146,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,COMPARE,50000,0.000018839,2654097844.572,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,COMPARE,50000,0.000583492,85690945.540,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,COMPARE,50000,0.000011066,4518354755.092,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,COMPARE,50000,0.000581270,86018533.200,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,COMPARE,50000,0.000010746,4653067360.028,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,COMPARE,50000,0.000608921,82112451.057,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,REDUCE,3125,0.000074654,41859782.839,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,REDUCE,3125,0.000132172,23643504.004,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,REDUCE,3125,0.000261479,11951231.071,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,REDUCE,50000,0.000417430,119780575.263,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,REDUCE,50000,0.000986401,50689324.546,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,REDUCE,50000,0.000323428,154593988.678,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,REDUCE,50000,0.000905419,55223065.545,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,REDUCE,50000,0.000089917,556066321.415,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,REDUCE,50000,0.000662720,75446663.388,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,REDUCE,50000,0.000092100,542886089.876,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,REDUCE,50000,0.000678286,73715222.215,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,REDUCE,50000,0.000099029,504900605.650,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,REDUCE,50000,0.000667648,74889718.531,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,REDUCE,50000,0.000107201,466414358.939,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,REDUCE,50000,0.000689812,72483516.768,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL,1562,0.000249552,6259223.204,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL,1562,0.000376455,4149235.793,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL,1562,0.000744428,2098254.286,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,MODMUL,50000,0.001154846,43295831.794,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,MODMUL,50000,0.001730366,28895626.292,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,MODMUL,50000,0.000817357,61172807.675,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,MODMUL,50000,0.001402399,35653202.321,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,MODMUL,50000,0.000331052,151033835.213,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,MODMUL,50000,0.000909584,54970205.742,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,MODMUL,50000,0.000266012,187961583.598,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,MODMUL,50000,0.000854403,58520408.804,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,MODMUL,50000,0.000398392,125504573.018,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,MODMUL,50000,0.000978049,51122188.746,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,MODMUL,50000,0.000329629,151685873.152,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,MODMUL,50000,0.000910938,54888490.480,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MODEXP,390,0.019147336,20368.368,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MODEXP,390,0.006027419,64704.311,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MODEXP,390,0.009463994,41208.818,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,MODEXP,50000,0.224474479,222742.471,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,MODEXP,50000,0.225202620,222022.284,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,MODEXP,50000,0.022341885,2237949.069,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,MODEXP,50000,0.022999369,2173972.723,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,MODEXP,50000,0.019462051,2569102.318,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,MODEXP,50000,0.020090250,2488769.391,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,MODEXP,50000,0.008686766,5755881.975,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,MODEXP,50000,0.009338232,5354332.737,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,MODEXP,50000,0.020807886,2402935.061,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,MODEXP,50000,0.021521073,2323304.251,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,MODEXP,50000,0.008568190,5835538.168,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,MODEXP,50000,0.009265672,5396262.436,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,EXPONENTIATION,390,0.003529013,110512.482,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,EXPONENTIATION,390,0.002861685,136283.338,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,EXPONENTIATION,390,0.005892972,66180.534,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,EXPONENTIATION,50000,0.196294200,254719.702,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,EXPONENTIATION,50000,0.197149683,253614.408,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,EXPONENTIATION,50000,0.052118313,959355.695,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,EXPONENTIATION,50000,0.052880293,945531.830,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,EXPONENTIATION,50000,0.015435757,3239232.104,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,EXPONENTIATION,50000,0.016066531,3112059.528,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,EXPONENTIATION,50000,0.014298439,3496885.213,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,EXPONENTIATION,50000,0.015040139,3324437.419,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,EXPONENTIATION,50000,0.015445022,3237289.030,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,EXPONENTIATION,50000,0.016195675,3087243.903,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,EXPONENTIATION,50000,0.014870953,3362259.340,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,EXPONENTIATION,50000,0.015476430,3230719.264,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,DIVIDE,3125,0.000083318,37507059.929,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,DIVIDE,3125,0.000262275,11914977.510,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,DIVIDE,3125,0.000202828,15407146.502,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,DIVIDE,50000,0.001147294,43580792.566,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,DIVIDE,50000,0.001827952,27353023.007,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,DIVIDE,50000,0.001127493,44346187.624,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,DIVIDE,50000,0.001822535,27434316.164,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,DIVIDE,50000,0.000379082,131897648.365,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,DIVIDE,50000,0.001049807,47627788.670,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,DIVIDE,50000,0.000353836,141308593.207,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,DIVIDE,50000,0.001036447,48241754.877,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,DIVIDE,50000,0.000334775,149354010.160,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,DIVIDE,50000,0.001023598,48847300.651,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,DIVIDE,50000,0.000323288,154660791.408,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,DIVIDE,50000,0.001027934,48641243.196,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,ISQRT,781,0.000069316,11267212.881,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,ISQRT,781,0.000136519,5720799.364,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,ISQRT,50000,0.017172311,2911664.100,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,ISQRT,50000,0.017788885,2810743.874,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,ISQRT,50000,0.017068041,2929451.515,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,ISQRT,50000,0.017705886,2823919.733,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,ISQRT,50000,0.003357012,14894198.355,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,ISQRT,50000,0.004004503,12485944.981,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,ISQRT,50000,0.003059354,16343321.250,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,ISQRT,50000,0.003654933,13680142.817,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,ISQRT,50000,0.004149843,12048648.566,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,ISQRT,50000,0.004815068,10384069.135,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,ISQRT,50000,0.003902618,12811912.933,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,ISQRT,50000,0.004557917,10969923.576,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL_R2,25000,0.004031999,6200398.307,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL_R2,25000,0.005452380,4585153.989,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL_R2,25000,0.002947404,8482040.325,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,MODMUL_R2,50000,0.000404321,123664223.193,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,MODMUL_R2,50000,0.000970928,51497124.484,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,MODMUL_R2,50000,0.000081904,610469062.130,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,MODMUL_R2,50000,0.000658635,75914575.144,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,MODMUL_R2,50000,0.000055004,909026264.815,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,MODMUL_R2,50000,0.000625955,79877954.038,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,MODMUL_R2,50000,0.000034573,1446233802.058,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,MODMUL_R2,50000,0.000610203,81940004.884,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,MODMUL_R2,50000,0.000053251,938946643.814,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,MODMUL_R2,50000,0.000632826,79010630.294,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,MODMUL_R2,50000,0.000029696,1683719851.973,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,MODMUL_R2,50000,0.000607317,82329284.682,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,p1024,1024,ADD,12500,0.000208114,60063370.569,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,p1024,1024,ADD,12500,0.000567463,22027880.791,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,p1024,1024,ADD,12500,0.001146913,10898818.406,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,ADD,50000,0.000190748,262126081.225,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,ADD,50000,0.001283329,38961171.119,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,ADD,50000,0.000097215,514322992.029,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,ADD,50000,0.001116728,44773643.750,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,p1024,1024,ADD,50000,0.000043947,1137727625.668,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,p1024,1024,ADD,50000,0.001095967,45621796.316,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,ADD,50000,0.000043044,1161605677.442,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,ADD,50000,0.001055956,47350462.858,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,ADD,50000,0.000042865,1166451379.655,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,ADD,50000,0.001095816,45628077.619,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,ADD,50000,0.000021882,2284945999.319,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,ADD,50000,0.001071338,46670628.858,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,ADD,50000,0.000021612,2313500439.541,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,ADD,50000,0.001075985,46469053.171,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,p1024,1024,SUBTRACT,12500,0.000232240,53823683.749,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,p1024,1024,SUBTRACT,12500,0.001206943,10356743.031,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,p1024,1024,SUBTRACT,12500,0.000853596,14643936.189,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,SUBTRACT,50000,0.000189537,263800481.539,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,SUBTRACT,50000,0.001296740,38558229.888,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,SUBTRACT,50000,0.000097306,513840577.324,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,SUBTRACT,50000,0.001110217,45036256.780,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,p1024,1024,SUBTRACT,50000,0.000044037,1135417714.237,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,p1024,1024,SUBTRACT,50000,0.001095435,45643982.473,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,SUBTRACT,50000,0.000043295,1154859129.238,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,SUBTRACT,50000,0.001050178,47610978.166,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,SUBTRACT,50000,0.000042733,1170061267.544,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,SUBTRACT,50000,0.001105040,45247218.102,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,SUBTRACT,50000,0.000023495,2128075598.541,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,SUBTRACT,50000,0.001074545,46531318.860,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,SUBTRACT,50000,0.000021484,2327340523.669,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,SUBTRACT,50000,0.001073603,46572167.957,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,p1024,1024,ADDMOD,12500,0.000704374,17746242.159,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,p1024,1024,ADDMOD,12500,0.000964274,12963120.319,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,p1024,1024,ADDMOD,12500,0.001893662,6600967.096,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,ADDMOD,50000,0.000266053,187932633.230,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,ADDMOD,50000,0.001370171,36491793.286,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,ADDMOD,50000,0.000135083,370143481.978,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,ADDMOD,50000,0.001152251,43393326.328,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,p1024,1024,ADDMOD,50000,0.000054983,909365005.590,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,p1024,1024,ADDMOD,50000,0.001077177,46417626.685,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,ADDMOD,50000,0.000058779,850636803.245,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,ADDMOD,50000,0.001072550,46617864.927,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,ADDMOD,50000,0.000057999,862083165.264,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,ADDMOD,50000,0.001123007,44523304.639,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,ADDMOD,50000,0.000016915,2956012069.155,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,ADDMOD,50000,0.001067042,46858496.257,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,ADDMOD,50000,0.000016335,3060837582.668,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,ADDMOD,50000,0.001083508,46146400.231,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,p1024,1024,SUBTRACTMOD,12500,0.000548900,22772826.885,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,p1024,1024,SUBTRACTMOD,12500,0.001854940,6738762.414,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,p1024,1024,SUBTRACTMOD,12500,0.000966612,12931762.457,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,SUBTRACTMOD,50000,0.000266472,187637061.114,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,SUBTRACTMOD,50000,0.001369702,36504298.752,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,SUBTRACTMOD,50000,0.000134373,372098329.660,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,SUBTRACTMOD,50000,0.001154656,43302955.795,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,p1024,1024,SUBTRACTMOD,50000,0.000055484,901153001.209,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,p1024,1024,SUBTRACTMOD,50000,0.001078989,46339660.060,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,SUBTRACTMOD,50000,0.000060251,829861984.110,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,SUBTRACTMOD,50000,0.001068896,46777249.464,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,SUBTRACTMOD,50000,0.000060050,832642005.025,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,SUBTRACTMOD,50000,0.001116376,44787762.743,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,SUBTRACTMOD,50000,0.000016566,3018163436.024,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,SUBTRACTMOD,50000,0.001080673,46267456.764,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,SUBTRACTMOD,50000,0.000016775,2980629091.717,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,SUBTRACTMOD,50000,0.001075836,46475489.495,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.001963915,6364838.611,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.001920215,6509687.159,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.003616146,3456719.315,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.022727778,2199951.073,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.024392914,2049775.601,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.005850701,8545984.984,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.007215204,6929811.265,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.001286404,38868032.081,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.002581250,19370460.630,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.000149295,334907994.810,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.001434138,34864147.272,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.000149755,333879097.999,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.001449741,34488908.989,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.000093933,532293190.561,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.001427155,35034737.235,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.000089055,561451247.621,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.001422089,35159553.517,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.001965571,6359476.564,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.001936366,6455390.645,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.001932400,6468640.139,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.008471582,5902085.289,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.009923708,5038439.416,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.002177199,22965283.371,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.003525417,14182718.178,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000558024,89601837.196,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.001847092,27069576.564,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000557873,89626069.595,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.001860760,26870738.788,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000557993,89606921.921,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.001886280,26507195.793,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000194924,256510292.502,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.001519855,32897872.705,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000180174,277510034.116,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.001517212,32955183.188,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.013251686,943276.166,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.005858123,2133789.160,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.006842921,1826705.410,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.001835555,27239720.089,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.003022987,16539930.029,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000272812,183276179.292,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.001293786,38646271.502,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000074493,671206101.068,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.001110520,45023944.073,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000129737,385395187.504,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.001141736,43792950.005,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000092840,538561996.670,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.001157109,43211152.285,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000094304,530200983.626,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.001148827,43522639.809,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000053620,932488470.490,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.001108246,45116340.076,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,p1024,1024,COMPARE,12500,0.000087009,143662756.086,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,p1024,1024,COMPARE,12500,0.000324479,38523258.143,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,p1024,1024,COMPARE,12500,0.000317143,39414439.480,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,COMPARE,50000,0.000144839,345210205.761,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,COMPARE,50000,0.001223039,40881774.680,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,COMPARE,50000,0.000069346,721019221.058,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,COMPARE,50000,0.001064299,46979293.689,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,COMPARE,50000,0.000041012,1219163666.091,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,COMPARE,50000,0.001009636,49522816.598,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,COMPARE,50000,0.000041053,1217946715.064,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,COMPARE,50000,0.001049107,47659583.655,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,COMPARE,50000,0.000013951,3583917970.628,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,COMPARE,50000,0.001051469,47552529.570,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,COMPARE,50000,0.000014212,3518158007.864,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,COMPARE,50000,0.001070729,46697177.490,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,p1024,1024,REDUCE,1562,0.000029272,53360821.403,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,p1024,1024,REDUCE,1562,0.000080198,19476738.301,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,p1024,1024,REDUCE,1562,0.000081222,19231204.095,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,REDUCE,50000,0.002198532,22742445.829,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,REDUCE,50000,0.003366958,14850198.548,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,REDUCE,50000,0.001194857,41846009.565,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,REDUCE,50000,0.002225444,22467429.543,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,REDUCE,50000,0.000322457,155059240.518,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,REDUCE,50000,0.001348309,37083481.519,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,REDUCE,50000,0.000319513,156488370.439,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,REDUCE,50000,0.001376713,36318398.305,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,REDUCE,50000,0.000315085,158687311.421,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,REDUCE,50000,0.001361929,36712636.560,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,REDUCE,50000,0.000275787,181299358.377,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,REDUCE,50000,0.001335358,37443135.658,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,p1024,1024,MODMUL,781,0.000354892,2200667.271,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,p1024,1024,MODMUL,781,0.000854860,913600.360,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,p1024,1024,MODMUL,781,0.000612275,1275570.107,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,MODMUL,50000,0.008449402,5917578.751,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,MODMUL,50000,0.009682933,5163724.739,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,MODMUL,50000,0.003070099,16286118.282,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,MODMUL,50000,0.004155669,12031756.055,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,MODMUL,50000,0.001442570,34660357.352,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,MODMUL,50000,0.002470374,20239849.443,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,MODMUL,50000,0.001102185,45364440.261,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,MODMUL,50000,0.002167074,23072585.327,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,MODMUL,50000,0.001412304,35403139.701,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,MODMUL,50000,0.002525026,19801776.466,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,MODMUL,50000,0.001076847,46431837.932,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,MODMUL,50000,0.002132352,23448281.356,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,p1024,1024,MODEXP,195,0.057963520,3364.185,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,p1024,1024,MODEXP,195,0.006150275,31705.897,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,p1024,1024,MODEXP,195,0.008383980,23258.643,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,MODEXP,50000,2.663202697,18774.388,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,MODEXP,50000,2.671068354,18719.102,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,MODEXP,50000,0.156373492,319747.289,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,MODEXP,50000,0.157608720,317241.331,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,MODEXP,50000,0.155846801,320827.889,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,MODEXP,50000,0.157068357,318332.737,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,MODEXP,50000,0.081636101,612474.131,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,MODEXP,50000,0.082968315,602639.694,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,MODEXP,50000,0.157279342,317905.704,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,MODEXP,50000,0.158658098,315143.069,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,MODEXP,50000,0.082264671,607794.322,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,MODEXP,50000,0.083658438,597668.342,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,p1024,1024,EXPONENTIATION,195,0.007910063,24652.141,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,p1024,1024,EXPONENTIATION,195,0.005910223,32993.677,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,p1024,1024,EXPONENTIATION,195,0.012236293,15936.198,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,EXPONENTIATION,50000,1.618205179,30898.430,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,EXPONENTIATION,50000,1.622654101,30813.714,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,EXPONENTIATION,50000,0.392424453,127413.059,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,EXPONENTIATION,50000,0.393568903,127042.557,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,EXPONENTIATION,50000,0.112899948,442870.001,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,EXPONENTIATION,50000,0.114107033,438185.084,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,EXPONENTIATION,50000,0.104130616,480166.177,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,EXPONENTIATION,50000,0.105821902,472491.980,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,EXPONENTIATION,50000,0.110363955,453046.469,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,EXPONENTIATION,50000,0.111754747,447408.287,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,EXPONENTIATION,50000,0.103800688,481692.376,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,EXPONENTIATION,50000,0.105218649,475200.932,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,p1024,1024,DIVIDE,1562,0.000048168,32428091.918,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,p1024,1024,DIVIDE,1562,0.000111996,13946981.835,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,p1024,1024,DIVIDE,1562,0.000218354,7153530.281,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,DIVIDE,50000,0.012850299,3890959.902,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,DIVIDE,50000,0.014689730,3403738.419,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,DIVIDE,50000,0.006980389,7162924.222,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,DIVIDE,50000,0.008328727,6003318.744,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,DIVIDE,50000,0.001516543,32969714.120,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,DIVIDE,50000,0.002785219,17951910.482,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,DIVIDE,50000,0.001532236,32632046.418,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,DIVIDE,50000,0.002828525,17677055.460,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,DIVIDE,50000,0.001479696,33790713.565,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,DIVIDE,50000,0.002778148,17997599.475,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,DIVIDE,50000,0.001532957,32616701.823,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,DIVIDE,50000,0.002883777,17338372.423,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,p1024,1024,ISQRT,390,0.000071172,5479679.658,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,p1024,1024,ISQRT,390,0.000141707,2752154.594,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,ISQRT,50000,0.194914633,256522.556,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,ISQRT,50000,0.198208513,252259.599,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,ISQRT,50000,0.122298958,408834.228,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,ISQRT,50000,0.123418700,405124.993,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,ISQRT,50000,0.022800406,2192943.364,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,ISQRT,50000,0.023951966,2087511.286,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,ISQRT,50000,0.022459364,2226242.959,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,ISQRT,50000,0.023661338,2113151.861,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,ISQRT,50000,0.022423316,2229821.872,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,ISQRT,50000,0.023638634,2115181.443,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,ISQRT,50000,0.021800436,2293532.103,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,ISQRT,50000,0.023044717,2169694.726,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,p1024,1024,MODMUL_R2,12500,0.005758764,2170604.581,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,p1024,1024,MODMUL_R2,12500,0.005732434,2180574.653,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,p1024,1024,MODMUL_R2,12500,0.007476179,1671977.105,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,MODMUL_R2,50000,0.002755735,18143979.279,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,MODMUL_R2,50000,0.003966864,12604414.343,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,MODMUL_R2,50000,0.000252010,198404588.424,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,MODMUL_R2,50000,0.001287056,38848344.460,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,MODMUL_R2,50000,0.000200475,249408111.197,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,MODMUL_R2,50000,0.001250440,39985916.806,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,MODMUL_R2,50000,0.000137748,362981158.303,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,MODMUL_R2,50000,0.001182478,42284087.594,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,MODMUL_R2,50000,0.000164909,303196991.021,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,MODMUL_R2,50000,0.001220243,40975443.320,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,MODMUL_R2,50000,0.000099679,501607878.165,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,MODMUL_R2,50000,0.001155706,43263593.595,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,p2048,2048,ADD,6250,0.000144340,43300576.759,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,p2048,2048,ADD,6250,0.000335860,18608955.717,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,p2048,2048,ADD,6250,0.000336236,18588144.067,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,ADD,50000,0.000380475,131414653.442,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,ADD,50000,0.002337141,21393660.375,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,ADD,50000,0.000187784,266262751.944,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,ADD,50000,0.002152022,23233963.499,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,p2048,2048,ADD,50000,0.000078740,635004508.788,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,p2048,2048,ADD,50000,0.002187235,22859909.237,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,ADD,50000,0.000081843,610927549.557,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,ADD,50000,0.002070550,24148174.643,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,ADD,50000,0.000086442,578424961.214,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,ADD,50000,0.002102157,23785091.978,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,ADD,50000,0.000068173,733430207.650,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,ADD,50000,0.002139233,23372862.501,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,ADD,50000,0.000074472,671390765.845,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,ADD,50000,0.002159024,23158616.084,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,p2048,2048,SUBTRACT,6250,0.000143510,43550958.582,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,p2048,2048,SUBTRACT,6250,0.000309798,20174456.426,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,p2048,2048,SUBTRACT,6250,0.000301413,20735638.471,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,SUBTRACT,50000,0.000380907,131265565.433,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,SUBTRACT,50000,0.002334129,21421266.178,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,SUBTRACT,50000,0.000187313,266932625.319,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,SUBTRACT,50000,0.002167726,23065646.436,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,p2048,2048,SUBTRACT,50000,0.000077708,643436937.606,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,p2048,2048,SUBTRACT,50000,0.002113944,23652472.300,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,SUBTRACT,50000,0.000080170,623673836.574,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,SUBTRACT,50000,0.002064379,24220359.558,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,SUBTRACT,50000,0.000085330,585962882.277,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,SUBTRACT,50000,0.002124120,23539164.733,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,SUBTRACT,50000,0.000067092,745240022.210,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,SUBTRACT,50000,0.002155086,23200929.989,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,SUBTRACT,50000,0.000074333,672652557.195,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,SUBTRACT,50000,0.002150279,23252801.489,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,p2048,2048,ADDMOD,6250,0.000410295,15232938.752,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,p2048,2048,ADDMOD,6250,0.001085150,5759573.103,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,p2048,2048,ADDMOD,6250,0.000592586,10546998.510,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,ADDMOD,50000,0.000479234,104333081.734,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,ADDMOD,50000,0.002464466,20288372.459,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,ADDMOD,50000,0.000241935,206666863.760,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,ADDMOD,50000,0.002229389,22427671.623,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,p2048,2048,ADDMOD,50000,0.000096334,519026771.593,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,p2048,2048,ADDMOD,50000,0.001783846,28029326.247,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,ADDMOD,50000,0.000111548,448236605.607,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,ADDMOD,50000,0.002095196,23864111.304,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,ADDMOD,50000,0.000112591,444084000.860,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,ADDMOD,50000,0.002147745,23280227.429,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,ADDMOD,50000,0.000039190,1275833916.350,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,ADDMOD,50000,0.002098501,23826534.496,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,ADDMOD,50000,0.000041854,1194639323.543,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,ADDMOD,50000,0.002101706,23790193.248,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,p2048,2048,SUBTRACTMOD,6250,0.000356112,17550665.867,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,p2048,2048,SUBTRACTMOD,6250,0.000522033,11972433.533,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,p2048,2048,SUBTRACTMOD,6250,0.001016708,6147288.464,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,SUBTRACTMOD,50000,0.000505814,98850501.368,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,SUBTRACTMOD,50000,0.002505237,19958189.695,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,SUBTRACTMOD,50000,0.000258220,193633066.197,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,SUBTRACTMOD,50000,0.002229810,22423437.584,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,p2048,2048,SUBTRACTMOD,50000,0.000108644,460217144.425,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,p2048,2048,SUBTRACTMOD,50000,0.002099784,23811972.008,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,SUBTRACTMOD,50000,0.000110416,452833981.680,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,SUBTRACTMOD,50000,0.002103768,23766875.972,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,SUBTRACTMOD,50000,0.000114953,434959825.002,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,SUBTRACTMOD,50000,0.002155898,23192190.344,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,SUBTRACTMOD,50000,0.000037577,1330601050.858,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,SUBTRACTMOD,50000,0.002108727,23710991.649,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,SUBTRACTMOD,50000,0.000038899,1285364183.107,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,SUBTRACTMOD,50000,0.002107875,23720566.921,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.003207318,1948668.910,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.005510829,1134130.556,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.002843743,2197807.461,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.091704814,545227.649,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.094951781,526583.069,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.023124136,2162242.945,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.025926551,1928524.943,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.005181639,9649457.362,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.007826068,6388904.766,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.000432774,115533629.446,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.003027324,16516238.782,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.000434287,115131266.620,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.003087135,16196246.515,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.000340555,146819221.808,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.002989037,16727795.818,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.000316940,157758442.841,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.002967864,16847131.596,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.003213727,1944782.542,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.005547574,1126618.633,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.002863658,2182523.451,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.033495732,1492727.507,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.036242064,1379612.372,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.008509709,5875641.874,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.011135669,4490075.818,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.002157683,23173010.289,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.004749287,10527896.260,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.002178613,22950380.717,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.004793314,10431195.881,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.002177641,22960627.895,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.004842158,10325973.625,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.000938559,53273158.756,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.003307709,15116200.636,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.000877665,56969321.700,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.003523454,14190620.625,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.021788379,286850.158,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.005894037,1060393.757,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.005779594,1081390.804,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.042640463,1172595.161,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.044664493,1119457.467,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.001329189,37616918.090,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.003379537,14794926.958,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000226952,220310771.150,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.002160475,23143062.481,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000400005,124998466.131,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.002365176,21140079.099,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000292072,171190622.748,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.002290461,21829664.580,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000341656,146346166.553,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.002150571,23249639.566,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000260543,191906844.536,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.002349941,21277129.394,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,p2048,2048,COMPARE,6250,0.000043520,143611505.113,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,p2048,2048,COMPARE,6250,0.000194497,32134104.578,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,p2048,2048,COMPARE,6250,0.000371640,16817354.386,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,COMPARE,50000,0.000290599,172058569.103,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,COMPARE,50000,0.002145384,23305856.429,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,COMPARE,50000,0.000139331,358856537.839,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,COMPARE,50000,0.002072571,24124627.686,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,COMPARE,50000,0.000074532,670853841.156,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,COMPARE,50000,0.001965880,25433897.277,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,COMPARE,50000,0.000074202,673834515.651,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,COMPARE,50000,0.002055896,24320291.080,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,COMPARE,50000,0.000021692,2304958406.320,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,COMPARE,50000,0.002112534,23668259.273,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,COMPARE,50000,0.000023816,2099448271.547,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,COMPARE,50000,0.002123058,23550936.298,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,p2048,2048,REDUCE,781,0.000020857,37446002.945,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,p2048,2048,REDUCE,781,0.000060308,12950115.140,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,p2048,2048,REDUCE,781,0.000060532,12902330.239,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,REDUCE,50000,0.139730405,357831.926,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,REDUCE,50000,0.142788716,350167.727,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,REDUCE,50000,0.004669247,10708364.306,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,REDUCE,50000,0.006739914,7418491.824,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,REDUCE,50000,0.001001464,49926896.888,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,REDUCE,50000,0.002982236,16765940.905,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,REDUCE,50000,0.000959590,52105582.558,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,REDUCE,50000,0.003066273,16306438.917,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,REDUCE,50000,0.000891777,56067843.992,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,REDUCE,50000,0.003023509,16537076.971,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,REDUCE,50000,0.000818286,61103323.818,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,REDUCE,50000,0.002991442,16714349.154,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,p2048,2048,MODMUL,390,0.000571402,682531.711,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,p2048,2048,MODMUL,390,0.000015026,25954278.804,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,p2048,2048,MODMUL,390,0.000042940,9082435.478,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,MODMUL,50000,0.312067669,160221.660,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,MODMUL,50000,0.313985968,159242.785,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,MODMUL,50000,0.014859717,3364801.554,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,MODMUL,50000,0.017018478,2937982.993,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,MODMUL,50000,0.005551530,9006526.051,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,MODMUL,50000,0.007591993,6585886.278,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,MODMUL,50000,0.003784480,13211856.222,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,MODMUL,50000,0.005933316,8426990.307,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,MODMUL,50000,0.005483449,9118349.511,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,MODMUL,50000,0.007700367,6493197.097,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,MODMUL,50000,0.003614252,13834119.482,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,MODMUL,50000,0.005806405,8611180.096,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,p2048,2048,MODEXP,97,0.208140874,466.031,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,p2048,2048,MODEXP,97,0.008253083,11753.184,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,p2048,2048,MODEXP,97,0.006124627,15837.700,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,MODEXP,50000,41.174704803,1214.338,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,MODEXP,50000,41.176832391,1214.275,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,MODEXP,50000,5.985886863,8352.981,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,MODEXP,50000,5.990351891,8346.755,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,MODEXP,50000,1.228589462,40697.077,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,MODEXP,50000,1.235523360,40468.680,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,MODEXP,50000,1.618831165,30886.482,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,MODEXP,50000,1.740665555,28724.645,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,MODEXP,50000,1.252043365,39934.719,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,MODEXP,50000,1.254966434,39841.703,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,MODEXP,50000,1.641353380,30462.666,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,MODEXP,50000,1.637912832,30526.655,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,p2048,2048,EXPONENTIATION,97,0.024642453,3936.297,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,p2048,2048,EXPONENTIATION,97,0.000967270,100282.262,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,p2048,2048,EXPONENTIATION,97,0.046805019,2072.427,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,EXPONENTIATION,50000,13.687958580,3652.846,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,EXPONENTIATION,50000,13.691817531,3651.816,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,EXPONENTIATION,50000,3.162340090,15811.076,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,EXPONENTIATION,50000,3.165997395,15792.811,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,EXPONENTIATION,50000,0.942989884,53022.838,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,EXPONENTIATION,50000,0.945207806,52898.421,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,EXPONENTIATION,50000,0.807524275,61917.643,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,EXPONENTIATION,50000,0.809451288,61770.240,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,EXPONENTIATION,50000,0.915632101,54607.085,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,EXPONENTIATION,50000,0.918348525,54445.560,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,EXPONENTIATION,50000,0.823549511,60712.804,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,EXPONENTIATION,50000,0.824798649,60620.856,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,p2048,2048,DIVIDE,781,0.000032396,24108085.603,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,p2048,2048,DIVIDE,781,0.000076151,10255946.422,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,p2048,2048,DIVIDE,781,0.000076333,10231544.147,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,DIVIDE,50000,0.680083023,73520.435,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,DIVIDE,50000,0.681787860,73336.595,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,DIVIDE,50000,0.145802287,342930.148,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,DIVIDE,50000,0.148261722,337241.462,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,DIVIDE,50000,0.012155840,4113249.211,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,DIVIDE,50000,0.014812484,3375530.985,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,DIVIDE,50000,0.012049239,4149639.613,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,DIVIDE,50000,0.014248462,3509150.575,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,DIVIDE,50000,0.011389052,4390180.836,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,DIVIDE,50000,0.014321975,3491138.517,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,DIVIDE,50000,0.011213236,4459016.289,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,DIVIDE,50000,0.013869820,3604949.409,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,p2048,2048,ISQRT,195,0.000060338,3231798.617,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,p2048,2048,ISQRT,195,0.000109162,1786344.070,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,ISQRT,50000,9.765763927,5119.927,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,ISQRT,50000,9.776821485,5114.137,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,ISQRT,50000,4.141952042,12071.603,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,ISQRT,50000,4.180727566,11959.641,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,ISQRT,50000,0.061074888,818667.072,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,ISQRT,50000,0.063040389,793142.317,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,ISQRT,50000,0.052157173,958640.915,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,ISQRT,50000,0.054553924,916524.353,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,ISQRT,50000,0.057201326,874105.612,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,ISQRT,50000,0.059488861,840493.485,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,ISQRT,50000,0.052072914,960192.082,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,ISQRT,50000,0.054388884,919305.489,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-1t,p2048,2048,MODMUL_R2,6250,0.009446857,661595.668,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,gmp-nt,p2048,2048,MODMUL_R2,6250,0.005760437,1084987.152,0
library,AMD EPYC 9554 64-Core Processor,host-cpu,openssl-nt,p2048,2048,MODMUL_R2,6250,0.005869050,1064908.367,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,MODMUL_R2,50000,0.023438791,2133215.847,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,MODMUL_R2,50000,0.025829133,1935798.646,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,MODMUL_R2,50000,0.001223199,40876420.892,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,MODMUL_R2,50000,0.003283393,15228151.340,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,MODMUL_R2,50000,0.000704225,71000020.102,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,MODMUL_R2,50000,0.002723416,18359293.712,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,MODMUL_R2,50000,0.000489701,102103206.469,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,MODMUL_R2,50000,0.002528301,19776130.181,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,MODMUL_R2,50000,0.000644665,77559637.332,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,MODMUL_R2,50000,0.002777036,18004806.192,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,MODMUL_R2,50000,0.000444563,112469971.886,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,MODMUL_R2,50000,0.002467278,20265244.563,0
```
