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

### Device 1 - cpu-skylake-avx512-AMD Eng Sample: 100-000000897-03 (CPU)

| Property | Value |
|---|---|
| Model | cpu-skylake-avx512-AMD Eng Sample: 100-000000897-03 |
| Type | CPU |
| Vendor | AuthenticAMD |
| Device memory | 753.22 GiB |
| Max single allocation | 256.00 GiB |
| Local memory | 1024 KiB |
| Global cache | 32768 KiB |
| Compute units | 64 |
| Max clock | 2550 MHz |
| Max work-group size | 4096 |
| OpenCL version | OpenCL 3.0 PoCL HSTR: cpu-x86_64-pc-linux-gnu-skylake-avx512 |
| Driver | 5.0+debian |

### Host

| Property | Value |
|---|---|
| CPU | AMD Eng Sample: 100-000000897-03 |
| Logical cores | 64 |
| OpenMP threads used | 64 |
| RAM | 755.2 GB |
| OS | Ubuntu 24.04.4 LTS |
| Kernel | 6.8.0-85-generic |
| Arch | x86_64 |
| GMP | 6.3.0 |
| OpenSSL | OpenSSL 3.0.13 30 Jan 2024 |
| CGBN | cgbn_results.tsv loaded |

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
- Total wall time 2336.6 s.

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

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | GMP 64T | OpenSSL 64T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | 1.75 G | 2.65 G | 3.77 G | 3.63 G | 3.73 G | 5.16 G | 4.96 G | 66.80 M | 1.48 G | 1.23 G | 6.57 G |
| SUBTRACT | 50000 / 50000 | 1.73 G | 2.68 G | 3.82 G | 3.86 G | 3.82 G | 5.26 G | 5.32 G | 88.15 M | 1.39 G | 953.28 M | 6.94 G |
| ADDMOD | 50000 / 50000 | 1.27 G | 1.98 G | 3.58 G | 4.68 G | 4.69 G | 5.23 G | 5.52 G | 24.29 M | 404.12 M | 92.23 M | 5.85 G |
| SUBTRACTMOD | 50000 / 50000 | 1.27 G | 2.08 G | 3.52 G | 4.65 G | 4.69 G | 5.52 G | 5.47 G | 30.59 M | 459.58 M | 70.94 M | 5.43 G |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | 56.45 M | 198.91 M | 692.54 M | 2.65 G | 2.66 G | 4.44 G | 4.92 G | 53.06 M | 718.00 M | 39.39 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | 294.11 M | 883.96 M | 2.39 G | 2.31 G | 2.43 G | 2.87 G | 2.99 G | 53.02 M | 54.20 M | 243.93 M | 6.98 G |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | 629.81 M | 1.95 G | 3.98 G | 3.61 G | 3.95 G | 3.81 G | 4.64 G | 6.61 M | 102.22 M | 320.20 M | 5.17 G |
| COMPARE | 50000 / 50000 | 1.77 G | 2.67 G | - | 4.96 G | 4.96 G | 5.18 G | 5.06 G | 143.62 M | 1.60 G | 2.21 G | 6.94 G |
| REDUCE | 50000 / 6250 | 362.67 M | 541.49 M | - | 1.38 G | 1.42 G | 1.46 G | 1.47 G | 62.28 M | 817.38 M | 85.04 M | 4.14 G |
| MODMUL | 50000 / 3125 | 135.02 M | 227.60 M | - | 481.77 M | 607.22 M | 487.70 M | 598.19 M | 11.53 M | 163.09 M | 51.52 M | 1.54 G |
| MODEXP | 50000 / 781 | 2.97 M | 15.69 M | - | 19.70 M | 41.20 M | 19.91 M | 41.58 M | 108.99 k | 1.79 M | 885.95 k | 5.39 M |
| EXPONENTIATION | 50000 / 781 | 1.94 M | 6.70 M | - | 108.04 M | 125.74 M | 106.90 M | 123.61 M | 351.28 k | 5.27 M | 337.25 k | n/a |
| DIVIDE | 50000 / 6250 | 151.17 M | 181.26 M | - | 498.70 M | 505.00 M | 505.98 M | 521.46 M | 37.46 M | 485.23 M | 97.46 M | 3.24 G |
| ISQRT | 50000 / 1562 | 13.72 M | 18.37 M | - | 75.51 M | 83.60 M | 73.65 M | 83.45 M | 19.03 M | 263.03 M | n/a | n/a |
| MODMUL_R2 | 50000 / 50000 | 532.93 M | 1.75 G | - | 2.75 G | 3.60 G | 2.70 G | 3.89 G | 11.46 M | 186.87 M | 51.67 M | 3.43 G |

#### rsa256(composite) (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | GMP 64T | OpenSSL 64T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | 1.76 G | 2.71 G | 3.69 G | 3.77 G | 3.86 G | 5.09 G | 5.21 G | 65.37 M | 1.21 G | 1.09 G | 6.85 G |
| SUBTRACT | 50000 / 50000 | 1.72 G | 2.71 G | 3.79 G | 3.79 G | 3.83 G | 5.28 G | 5.33 G | 88.05 M | 1.15 G | 65.56 M | 6.85 G |
| ADDMOD | 50000 / 50000 | 1.36 G | 2.21 G | 3.63 G | 4.66 G | 4.48 G | 5.44 G | 5.71 G | 28.61 M | 128.73 M | 100.25 M | 5.70 G |
| SUBTRACTMOD | 50000 / 50000 | 1.32 G | 2.06 G | 3.46 G | 4.73 G | 4.61 G | 5.59 G | 5.47 G | 30.59 M | 440.31 M | 74.79 M | 5.48 G |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | 54.63 M | 198.59 M | 677.77 M | 2.63 G | 2.52 G | 4.56 G | 4.95 G | 52.95 M | 728.38 M | 420.04 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | 294.89 M | 880.05 M | 2.35 G | 2.28 G | 2.37 G | 2.85 G | 2.98 G | 52.90 M | 723.97 M | 420.20 M | 6.88 G |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | 635.59 M | 1.95 G | 4.17 G | 3.51 G | 3.98 G | 3.71 G | 4.68 G | 6.62 M | 107.56 M | 192.88 M | 5.11 G |
| COMPARE | 50000 / 50000 | 1.76 G | 2.78 G | - | 5.09 G | 4.96 G | 5.28 G | 5.25 G | 123.49 M | 1.14 G | 2.07 G | 6.76 G |
| REDUCE | 50000 / 6250 | 358.92 M | 536.89 M | - | 1.35 G | 1.45 G | 1.42 G | 1.46 G | 39.22 M | 594.75 M | 79.62 M | 4.07 G |
| MODMUL | 50000 / 3125 | 134.16 M | 228.34 M | - | 481.90 M | 608.25 M | 483.68 M | 597.26 M | 11.52 M | 173.78 M | 51.11 M | 1.55 G |
| MODEXP | 50000 / 781 | 2.97 M | 15.73 M | - | 19.71 M | 41.30 M | 19.99 M | 41.86 M | 116.20 k | 1.82 M | 885.17 k | 5.48 M |
| EXPONENTIATION | 50000 / 781 | 1.95 M | 6.71 M | - | 108.07 M | 127.58 M | 106.95 M | 124.08 M | 351.99 k | 5.29 M | 335.88 k | n/a |
| DIVIDE | 50000 / 6250 | 150.45 M | 177.94 M | - | 477.02 M | 499.76 M | 482.00 M | 504.40 M | 37.01 M | 510.38 M | 91.59 M | 3.26 G |
| ISQRT | 50000 / 1562 | 13.69 M | 17.88 M | - | 72.46 M | 77.52 M | 70.79 M | 76.73 M | 19.11 M | 274.56 M | n/a | n/a |
| MODMUL_R2 | 50000 / 50000 | 530.27 M | 1.75 G | - | 2.67 G | 3.57 G | 2.70 G | 3.81 G | 11.47 M | 179.01 M | 44.26 M | 3.41 G |

#### brainpoolP512r1 (512-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | GMP 64T | OpenSSL 64T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 25000 | 846.77 M | 1.48 G | 2.50 G | 2.57 G | 2.52 G | 3.43 G | 3.44 G | 58.49 M | 50.74 M | 921.47 M | 6.62 G |
| SUBTRACT | 50000 / 25000 | 824.54 M | 1.51 G | 2.50 G | 2.52 G | 2.59 G | 3.36 G | 3.41 G | 85.03 M | 1.01 G | 758.68 M | 6.51 G |
| ADDMOD | 50000 / 25000 | 661.88 M | 1.22 G | 2.28 G | 2.02 G | 2.09 G | 4.09 G | 4.10 G | 25.69 M | 366.48 M | 95.18 M | 5.24 G |
| SUBTRACTMOD | 50000 / 25000 | 564.12 M | 1.05 G | 2.04 G | 2.02 G | 2.15 G | 4.00 G | 4.05 G | 27.49 M | 410.35 M | 23.51 M | 5.26 G |
| MULTIPLYOPERANDSCANNING | 50000 / 25000 | 10.43 M | 38.74 M | 144.72 M | 819.11 M | 830.55 M | 1.51 G | 1.55 G | 21.74 M | 266.88 M | 274.64 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 25000 | 43.99 M | 158.85 M | 558.32 M | 559.38 M | 551.17 M | 828.63 M | 900.53 M | 21.81 M | 271.74 M | 292.52 M | 5.37 G |
| MONTGOMERYMULTIPLICATION | 50000 / 25000 | 184.77 M | 595.05 M | 1.87 G | 1.30 G | 1.81 G | 1.55 G | 2.47 G | 2.79 M | 45.65 M | 151.78 M | 3.86 G |
| COMPARE | 50000 / 25000 | 839.93 M | 1.48 G | - | 2.57 G | 2.64 G | 4.39 G | 4.51 G | 142.82 M | 1.67 G | 1.61 G | 6.62 G |
| REDUCE | 50000 / 3125 | 119.32 M | 154.32 M | - | 555.35 M | 542.19 M | 506.44 M | 463.22 M | 38.32 M | 443.68 M | 86.76 M | 2.93 G |
| MODMUL | 50000 / 1562 | 43.27 M | 61.19 M | - | 150.77 M | 186.94 M | 125.03 M | 151.33 M | 5.84 M | 90.56 M | 32.85 M | 584.77 M |
| MODEXP | 50000 / 390 | 222.72 k | 2.24 M | - | 2.57 M | 5.77 M | 2.40 M | 5.83 M | 19.07 k | 242.66 k | 271.12 k | 2.28 M |
| EXPONENTIATION | 50000 / 390 | 258.30 k | 991.85 k | - | 3.00 M | 3.20 M | 3.36 M | 3.50 M | 103.00 k | 205.85 k | 117.14 k | n/a |
| DIVIDE | 50000 / 3125 | 43.31 M | 44.36 M | - | 130.15 M | 141.32 M | 148.59 M | 153.55 M | 35.60 M | 998.28 M | 154.88 M | 2.05 G |
| ISQRT | 50000 / 781 | 2.91 M | 2.93 M | - | 14.87 M | 16.28 M | 12.04 M | 12.81 M | 10.44 M | 328.61 M | n/a | n/a |
| MODMUL_R2 | 50000 / 25000 | 124.61 M | 603.91 M | - | 898.26 M | 1.46 G | 913.87 M | 1.65 G | 5.80 M | 135.92 M | 45.69 M | 2.38 G |

#### p1024 (1024-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | GMP 64T | OpenSSL 64T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 12500 | 261.61 M | 510.07 M | 1.13 G | 1.15 G | 1.15 G | 2.21 G | 2.26 G | 59.27 M | 788.69 M | 741.93 M | 4.75 G |
| SUBTRACT | 50000 / 12500 | 261.10 M | 510.85 M | 1.11 G | 1.16 G | 1.16 G | 2.26 G | 2.28 G | 70.74 M | 874.49 M | 509.78 M | 4.71 G |
| ADDMOD | 50000 / 12500 | 187.72 M | 369.27 M | 871.60 M | 836.55 M | 828.91 M | 2.93 G | 2.88 G | 18.27 M | 299.84 M | 35.25 M | 4.33 G |
| SUBTRACTMOD | 50000 / 12500 | 189.15 M | 371.00 M | 891.36 M | 833.75 M | 835.86 M | 2.96 G | 2.92 G | 22.64 M | 340.33 M | 32.47 M | 4.17 G |
| MULTIPLYOPERANDSCANNING | 50000 / 12500 | 2.19 M | 8.54 M | 38.97 M | 333.61 M | 334.24 M | 546.94 M | 552.76 M | 5.96 M | 76.02 M | 70.70 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 12500 | 5.90 M | 22.96 M | 89.60 M | 89.27 M | 89.37 M | 251.37 M | 274.74 M | 5.96 M | 74.87 M | 73.78 M | 2.26 G |
| MONTGOMERYMULTIPLICATION | 50000 / 12500 | 27.34 M | 183.25 M | 668.70 M | 380.70 M | 530.39 M | 528.93 M | 921.47 M | 880.01 k | 14.04 M | 45.00 M | 1.49 G |
| COMPARE | 50000 / 12500 | 354.43 M | 686.16 M | - | 1.24 G | 1.23 G | 3.48 G | 3.39 G | 140.54 M | 2.00 G | 1.65 G | 4.88 G |
| REDUCE | 50000 / 1562 | 22.70 M | 41.88 M | - | 154.00 M | 156.53 M | 157.83 M | 180.72 M | 49.59 M | 690.24 M | 35.48 M | 2.01 G |
| MODMUL | 50000 / 781 | 5.91 M | 16.29 M | - | 34.63 M | 45.26 M | 35.37 M | 46.36 M | 2.06 M | 28.42 M | 10.83 M | 238.55 M |
| MODEXP | 50000 / 195 | 18.74 k | 319.63 k | - | 320.66 k | 612.56 k | 317.82 k | 607.65 k | 3.14 k | 28.76 k | 36.95 k | 442.40 k |
| EXPONENTIATION | 50000 / 195 | 31.03 k | 128.90 k | - | 442.83 k | 494.63 k | 453.21 k | 485.03 k | 22.89 k | 242.53 k | 19.04 k | n/a |
| DIVIDE | 50000 / 1562 | 3.86 M | 7.08 M | - | 32.77 M | 32.60 M | 33.62 M | 32.48 M | 31.74 M | 561.75 M | 97.51 M | 1.52 G |
| ISQRT | 50000 / 390 | 251.78 k | 408.86 k | - | 2.19 M | 2.23 M | 2.23 M | 2.29 M | 5.12 M | 96.58 M | n/a | n/a |
| MODMUL_R2 | 50000 / 12500 | 18.18 M | 198.30 M | - | 246.64 M | 364.29 M | 299.94 M | 496.08 M | 2.03 M | 45.24 M | 13.10 M | 824.97 M |

#### p2048 (2048-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | GMP 64T | OpenSSL 64T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 6250 | 130.31 M | 265.80 M | 635.26 M | 622.04 M | 583.10 M | 718.65 M | 671.39 M | 44.04 M | 13.99 M | 16.64 M | 2.15 G |
| SUBTRACT | 50000 / 6250 | 130.96 M | 266.27 M | 636.97 M | 615.30 M | 584.47 M | 730.00 M | 669.78 M | 44.98 M | 37.53 M | 757.76 M | 2.11 G |
| ADDMOD | 50000 / 6250 | 103.87 M | 204.98 M | 517.63 M | 443.86 M | 444.85 M | 1.29 G | 1.20 G | 14.78 M | 316.27 M | 16.62 M | 2.14 G |
| SUBTRACTMOD | 50000 / 6250 | 98.60 M | 194.46 M | 455.65 M | 454.19 M | 443.54 M | 1.29 G | 1.25 G | 17.29 M | 243.14 M | 41.04 M | 2.10 G |
| MULTIPLYOPERANDSCANNING | 50000 / 6250 | 542.82 k | 2.16 M | 9.63 M | 116.01 M | 118.97 M | 146.22 M | 164.38 M | 1.82 M | 24.14 M | 22.96 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 6250 | 1.49 M | 5.87 M | 23.12 M | 23.18 M | 23.15 M | 53.54 M | 57.31 M | 1.82 M | 25.29 M | 23.96 M | 776.59 M |
| MONTGOMERYMULTIPLICATION | 50000 / 6250 | 1.16 M | 37.55 M | 220.11 M | 124.38 M | 171.11 M | 145.23 M | 187.78 M | 267.17 k | 4.28 M | 1.81 M | 509.46 M |
| COMPARE | 50000 / 6250 | 185.99 M | 377.62 M | - | 668.52 M | 652.61 M | 2.22 G | 2.03 G | 122.07 M | 1.98 G | 1.99 G | 2.13 G |
| REDUCE | 50000 / 781 | 363.03 k | 10.69 M | - | 49.61 M | 52.28 M | 55.96 M | 60.96 M | 34.13 M | 391.58 M | 45.29 M | 1.73 G |
| MODMUL | 50000 / 390 | 162.00 k | 3.36 M | - | 9.01 M | 13.19 M | 9.10 M | 13.83 M | 635.48 k | 8.60 M | 3.96 M | 115.66 M |
| MODEXP | 50000 / 97 | 1.21 k | 8.36 k | - | 40.73 k | 31.08 k | 39.51 k | 29.38 k | 433.5 | 7.22 k | 13.47 k | 79.25 k |
| EXPONENTIATION | 50000 / 97 | 3.63 k | 15.78 k | - | 52.97 k | 61.26 k | 54.21 k | 61.46 k | 3.64 k | 49.85 k | 4.84 k | n/a |
| DIVIDE | 50000 / 781 | 73.78 k | 339.46 k | - | 4.27 M | 4.45 M | 4.42 M | 4.55 M | 23.77 M | 310.84 M | 60.58 M | 1.48 G |
| ISQRT | 50000 / 195 | 5.15 k | 11.83 k | - | 817.92 k | 957.82 k | 867.11 k | 960.73 k | 3.00 M | 48.26 M | n/a | n/a |
| MODMUL_R2 | 50000 / 6250 | 2.10 M | 40.79 M | - | 70.93 M | 101.94 M | 76.52 M | 111.86 M | 618.11 k | 10.21 M | 4.58 M | 269.82 M |

### Device 1 - cpu-skylake-avx512-AMD Eng Sample: 100-000000897-03 (CPU)

#### secp256k1 (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | GMP 64T | OpenSSL 64T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | - | - | - | - | - | - | - | 66.80 M | 1.48 G | 1.23 G | 6.57 G |
| SUBTRACT | 50000 / 50000 | - | - | - | - | - | - | - | 88.15 M | 1.39 G | 953.28 M | 6.94 G |
| ADDMOD | 50000 / 50000 | - | - | - | - | - | - | - | 24.29 M | 404.12 M | 92.23 M | 5.85 G |
| SUBTRACTMOD | 50000 / 50000 | - | - | - | - | - | - | - | 30.59 M | 459.58 M | 70.94 M | 5.43 G |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | - | - | - | - | - | - | - | 53.06 M | 718.00 M | 39.39 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | - | - | - | - | - | - | - | 53.02 M | 54.20 M | 243.93 M | 6.98 G |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | - | - | - | - | - | - | - | 6.61 M | 102.22 M | 320.20 M | 5.17 G |
| COMPARE | 50000 / 50000 | - | - | - | - | - | - | - | 143.62 M | 1.60 G | 2.21 G | 6.94 G |
| REDUCE | 50000 / 6250 | - | - | - | - | - | - | - | 62.28 M | 817.38 M | 85.04 M | 4.14 G |
| MODMUL | 50000 / 3125 | - | - | - | - | - | - | - | 11.53 M | 163.09 M | 51.52 M | 1.54 G |
| MODEXP | 50000 / 781 | - | - | - | - | - | - | - | 108.99 k | 1.79 M | 885.95 k | 5.39 M |
| EXPONENTIATION | 50000 / 781 | - | - | - | - | - | - | - | 351.28 k | 5.27 M | 337.25 k | n/a |
| DIVIDE | 50000 / 6250 | - | - | - | - | - | - | - | 37.46 M | 485.23 M | 97.46 M | 3.24 G |
| ISQRT | 50000 / 1562 | - | - | - | - | - | - | - | 19.03 M | 263.03 M | n/a | n/a |
| MODMUL_R2 | 50000 / 50000 | - | - | - | - | - | - | - | 11.46 M | 186.87 M | 51.67 M | 3.43 G |

#### rsa256(composite) (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | GMP 64T | OpenSSL 64T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | - | - | - | - | - | - | - | 65.37 M | 1.21 G | 1.09 G | 6.85 G |
| SUBTRACT | 50000 / 50000 | - | - | - | - | - | - | - | 88.05 M | 1.15 G | 65.56 M | 6.85 G |
| ADDMOD | 50000 / 50000 | - | - | - | - | - | - | - | 28.61 M | 128.73 M | 100.25 M | 5.70 G |
| SUBTRACTMOD | 50000 / 50000 | - | - | - | - | - | - | - | 30.59 M | 440.31 M | 74.79 M | 5.48 G |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | - | - | - | - | - | - | - | 52.95 M | 728.38 M | 420.04 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | - | - | - | - | - | - | - | 52.90 M | 723.97 M | 420.20 M | 6.88 G |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | - | - | - | - | - | - | - | 6.62 M | 107.56 M | 192.88 M | 5.11 G |
| COMPARE | 50000 / 50000 | - | - | - | - | - | - | - | 123.49 M | 1.14 G | 2.07 G | 6.76 G |
| REDUCE | 50000 / 6250 | - | - | - | - | - | - | - | 39.22 M | 594.75 M | 79.62 M | 4.07 G |
| MODMUL | 50000 / 3125 | - | - | - | - | - | - | - | 11.52 M | 173.78 M | 51.11 M | 1.55 G |
| MODEXP | 50000 / 781 | - | - | - | - | - | - | - | 116.20 k | 1.82 M | 885.17 k | 5.48 M |
| EXPONENTIATION | 50000 / 781 | - | - | - | - | - | - | - | 351.99 k | 5.29 M | 335.88 k | n/a |
| DIVIDE | 50000 / 6250 | - | - | - | - | - | - | - | 37.01 M | 510.38 M | 91.59 M | 3.26 G |
| ISQRT | 50000 / 1562 | - | - | - | - | - | - | - | 19.11 M | 274.56 M | n/a | n/a |
| MODMUL_R2 | 50000 / 50000 | - | - | - | - | - | - | - | 11.47 M | 179.01 M | 44.26 M | 3.41 G |

#### brainpoolP512r1 (512-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | GMP 64T | OpenSSL 64T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 25000 | - | - | - | - | - | - | - | 58.49 M | 50.74 M | 921.47 M | 6.62 G |
| SUBTRACT | 50000 / 25000 | - | - | - | - | - | - | - | 85.03 M | 1.01 G | 758.68 M | 6.51 G |
| ADDMOD | 50000 / 25000 | - | - | - | - | - | - | - | 25.69 M | 366.48 M | 95.18 M | 5.24 G |
| SUBTRACTMOD | 50000 / 25000 | - | - | - | - | - | - | - | 27.49 M | 410.35 M | 23.51 M | 5.26 G |
| MULTIPLYOPERANDSCANNING | 50000 / 25000 | - | - | - | - | - | - | - | 21.74 M | 266.88 M | 274.64 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 25000 | - | - | - | - | - | - | - | 21.81 M | 271.74 M | 292.52 M | 5.37 G |
| MONTGOMERYMULTIPLICATION | 50000 / 25000 | - | - | - | - | - | - | - | 2.79 M | 45.65 M | 151.78 M | 3.86 G |
| COMPARE | 50000 / 25000 | - | - | - | - | - | - | - | 142.82 M | 1.67 G | 1.61 G | 6.62 G |
| REDUCE | 50000 / 3125 | - | - | - | - | - | - | - | 38.32 M | 443.68 M | 86.76 M | 2.93 G |
| MODMUL | 50000 / 1562 | - | - | - | - | - | - | - | 5.84 M | 90.56 M | 32.85 M | 584.77 M |
| MODEXP | 50000 / 390 | - | - | - | - | - | - | - | 19.07 k | 242.66 k | 271.12 k | 2.28 M |
| EXPONENTIATION | 50000 / 390 | - | - | - | - | - | - | - | 103.00 k | 205.85 k | 117.14 k | n/a |
| DIVIDE | 50000 / 3125 | - | - | - | - | - | - | - | 35.60 M | 998.28 M | 154.88 M | 2.05 G |
| ISQRT | 50000 / 781 | - | - | - | - | - | - | - | 10.44 M | 328.61 M | n/a | n/a |
| MODMUL_R2 | 50000 / 25000 | - | - | - | - | - | - | - | 5.80 M | 135.92 M | 45.69 M | 2.38 G |

#### p1024 (1024-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | GMP 64T | OpenSSL 64T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 12500 | - | - | - | - | - | - | - | 59.27 M | 788.69 M | 741.93 M | 4.75 G |
| SUBTRACT | 50000 / 12500 | - | - | - | - | - | - | - | 70.74 M | 874.49 M | 509.78 M | 4.71 G |
| ADDMOD | 50000 / 12500 | - | - | - | - | - | - | - | 18.27 M | 299.84 M | 35.25 M | 4.33 G |
| SUBTRACTMOD | 50000 / 12500 | - | - | - | - | - | - | - | 22.64 M | 340.33 M | 32.47 M | 4.17 G |
| MULTIPLYOPERANDSCANNING | 50000 / 12500 | - | - | - | - | - | - | - | 5.96 M | 76.02 M | 70.70 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 12500 | - | - | - | - | - | - | - | 5.96 M | 74.87 M | 73.78 M | 2.26 G |
| MONTGOMERYMULTIPLICATION | 50000 / 12500 | - | - | - | - | - | - | - | 880.01 k | 14.04 M | 45.00 M | 1.49 G |
| COMPARE | 50000 / 12500 | - | - | - | - | - | - | - | 140.54 M | 2.00 G | 1.65 G | 4.88 G |
| REDUCE | 50000 / 1562 | - | - | - | - | - | - | - | 49.59 M | 690.24 M | 35.48 M | 2.01 G |
| MODMUL | 50000 / 781 | - | - | - | - | - | - | - | 2.06 M | 28.42 M | 10.83 M | 238.55 M |
| MODEXP | 50000 / 195 | - | - | - | - | - | - | - | 3.14 k | 28.76 k | 36.95 k | 442.40 k |
| EXPONENTIATION | 50000 / 195 | - | - | - | - | - | - | - | 22.89 k | 242.53 k | 19.04 k | n/a |
| DIVIDE | 50000 / 1562 | - | - | - | - | - | - | - | 31.74 M | 561.75 M | 97.51 M | 1.52 G |
| ISQRT | 50000 / 390 | - | - | - | - | - | - | - | 5.12 M | 96.58 M | n/a | n/a |
| MODMUL_R2 | 50000 / 12500 | - | - | - | - | - | - | - | 2.03 M | 45.24 M | 13.10 M | 824.97 M |

#### p2048 (2048-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | GMP 64T | OpenSSL 64T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 6250 | - | - | - | - | - | - | - | 44.04 M | 13.99 M | 16.64 M | 2.15 G |
| SUBTRACT | 50000 / 6250 | - | - | - | - | - | - | - | 44.98 M | 37.53 M | 757.76 M | 2.11 G |
| ADDMOD | 50000 / 6250 | - | - | - | - | - | - | - | 14.78 M | 316.27 M | 16.62 M | 2.14 G |
| SUBTRACTMOD | 50000 / 6250 | - | - | - | - | - | - | - | 17.29 M | 243.14 M | 41.04 M | 2.10 G |
| MULTIPLYOPERANDSCANNING | 50000 / 6250 | - | - | - | - | - | - | - | 1.82 M | 24.14 M | 22.96 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 6250 | - | - | - | - | - | - | - | 1.82 M | 25.29 M | 23.96 M | 776.59 M |
| MONTGOMERYMULTIPLICATION | 50000 / 6250 | - | - | - | - | - | - | - | 267.17 k | 4.28 M | 1.81 M | 509.46 M |
| COMPARE | 50000 / 6250 | - | - | - | - | - | - | - | 122.07 M | 1.98 G | 1.99 G | 2.13 G |
| REDUCE | 50000 / 781 | - | - | - | - | - | - | - | 34.13 M | 391.58 M | 45.29 M | 1.73 G |
| MODMUL | 50000 / 390 | - | - | - | - | - | - | - | 635.48 k | 8.60 M | 3.96 M | 115.66 M |
| MODEXP | 50000 / 97 | - | - | - | - | - | - | - | 433.5 | 7.22 k | 13.47 k | 79.25 k |
| EXPONENTIATION | 50000 / 97 | - | - | - | - | - | - | - | 3.64 k | 49.85 k | 4.84 k | n/a |
| DIVIDE | 50000 / 781 | - | - | - | - | - | - | - | 23.77 M | 310.84 M | 60.58 M | 1.48 G |
| ISQRT | 50000 / 195 | - | - | - | - | - | - | - | 3.00 M | 48.26 M | n/a | n/a |
| MODMUL_R2 | 50000 / 6250 | - | - | - | - | - | - | - | 618.11 k | 10.21 M | 4.58 M | 269.82 M |

## 5. Head to head

Best OpenCL GPU result against best OpenCL CPU result and the CPU libraries.
Ratios above 1.00x mean the GPU is faster than that baseline.

### secp256k1 (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GMP 64T | OpenSSL | CGBN | GPU vs CPU-CL | GPU vs GMP 64T | GPU vs OpenSSL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | w32-il | 5.16 G | none | n/a | 66.80 M | 1.48 G | 1.23 G | 6.57 G | n/a | 3.49x | 4.19x | 0.79x |
| SUBTRACT | w32-il64 | 5.32 G | none | n/a | 88.15 M | 1.39 G | 953.28 M | 6.94 G | n/a | 3.82x | 5.58x | 0.77x |
| ADDMOD | w32-il64 | 5.52 G | none | n/a | 24.29 M | 404.12 M | 92.23 M | 5.85 G | n/a | 13.65x | 59.81x | 0.94x |
| SUBTRACTMOD | w32-il | 5.52 G | none | n/a | 30.59 M | 459.58 M | 70.94 M | 5.43 G | n/a | 12.02x | 77.85x | 1.02x |
| MULTIPLYOPERANDSCANNING | w32-il64 | 4.92 G | none | n/a | 53.06 M | 718.00 M | 39.39 M | n/a | n/a | 6.86x | 124.99x | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 2.99 G | none | n/a | 53.02 M | 54.20 M | 243.93 M | 6.98 G | n/a | 55.12x | 12.25x | 0.43x |
| MONTGOMERYMULTIPLICATION | w32-il64 | 4.64 G | none | n/a | 6.61 M | 102.22 M | 320.20 M | 5.17 G | n/a | 45.39x | 14.49x | 0.90x |
| COMPARE | w32-il | 5.18 G | none | n/a | 143.62 M | 1.60 G | 2.21 G | 6.94 G | n/a | 3.24x | 2.35x | 0.75x |
| REDUCE | w32-il64 | 184.04 M | none | n/a | 62.28 M | 817.38 M | 85.04 M | 4.14 G | n/a | 0.23x | 2.16x | 0.04x |
| MODMUL | w32-o64 | 37.95 M | none | n/a | 11.53 M | 163.09 M | 51.52 M | 1.54 G | n/a | 0.23x | 0.74x | 0.02x |
| MODEXP | w32-il64 | 649.41 k | none | n/a | 108.99 k | 1.79 M | 885.95 k | 5.39 M | n/a | 0.36x | 0.73x | 0.12x |
| EXPONENTIATION | w32-o64 | 1.96 M | none | n/a | 351.28 k | 5.27 M | 337.25 k | n/a | n/a | 0.37x | 5.82x | n/a |
| DIVIDE | w32-il64 | 65.18 M | none | n/a | 37.46 M | 485.23 M | 97.46 M | 3.24 G | n/a | 0.13x | 0.67x | 0.02x |
| ISQRT | w32-o64 | 2.61 M | none | n/a | 19.03 M | 263.03 M | n/a | n/a | n/a | 0.01x | n/a | n/a |
| MODMUL_R2 | w32-il64 | 3.89 G | none | n/a | 11.46 M | 186.87 M | 51.67 M | 3.43 G | n/a | 20.81x | 75.26x | 1.13x |

### rsa256(composite) (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GMP 64T | OpenSSL | CGBN | GPU vs CPU-CL | GPU vs GMP 64T | GPU vs OpenSSL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | w32-il64 | 5.21 G | none | n/a | 65.37 M | 1.21 G | 1.09 G | 6.85 G | n/a | 4.32x | 4.76x | 0.76x |
| SUBTRACT | w32-il64 | 5.33 G | none | n/a | 88.05 M | 1.15 G | 65.56 M | 6.85 G | n/a | 4.62x | 81.27x | 0.78x |
| ADDMOD | w32-il64 | 5.71 G | none | n/a | 28.61 M | 128.73 M | 100.25 M | 5.70 G | n/a | 44.38x | 56.98x | 1.00x |
| SUBTRACTMOD | w32-il | 5.59 G | none | n/a | 30.59 M | 440.31 M | 74.79 M | 5.48 G | n/a | 12.70x | 74.75x | 1.02x |
| MULTIPLYOPERANDSCANNING | w32-il64 | 4.95 G | none | n/a | 52.95 M | 728.38 M | 420.04 M | n/a | n/a | 6.79x | 11.78x | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 2.98 G | none | n/a | 52.90 M | 723.97 M | 420.20 M | 6.88 G | n/a | 4.11x | 7.08x | 0.43x |
| MONTGOMERYMULTIPLICATION | w32-il64 | 4.68 G | none | n/a | 6.62 M | 107.56 M | 192.88 M | 5.11 G | n/a | 43.50x | 24.26x | 0.92x |
| COMPARE | w32-il | 5.28 G | none | n/a | 123.49 M | 1.14 G | 2.07 G | 6.76 G | n/a | 4.63x | 2.55x | 0.78x |
| REDUCE | w32-il64 | 182.21 M | none | n/a | 39.22 M | 594.75 M | 79.62 M | 4.07 G | n/a | 0.31x | 2.29x | 0.04x |
| MODMUL | w32-o64 | 38.02 M | none | n/a | 11.52 M | 173.78 M | 51.11 M | 1.55 G | n/a | 0.22x | 0.74x | 0.02x |
| MODEXP | w32-il64 | 653.93 k | none | n/a | 116.20 k | 1.82 M | 885.17 k | 5.48 M | n/a | 0.36x | 0.74x | 0.12x |
| EXPONENTIATION | w32-o64 | 1.99 M | none | n/a | 351.99 k | 5.29 M | 335.88 k | n/a | n/a | 0.38x | 5.93x | n/a |
| DIVIDE | w32-il64 | 63.05 M | none | n/a | 37.01 M | 510.38 M | 91.59 M | 3.26 G | n/a | 0.12x | 0.69x | 0.02x |
| ISQRT | w32-o64 | 2.42 M | none | n/a | 19.11 M | 274.56 M | n/a | n/a | n/a | 0.01x | n/a | n/a |
| MODMUL_R2 | w32-il64 | 3.81 G | none | n/a | 11.47 M | 179.01 M | 44.26 M | 3.41 G | n/a | 21.27x | 86.05x | 1.12x |

### brainpoolP512r1 (512-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GMP 64T | OpenSSL | CGBN | GPU vs CPU-CL | GPU vs GMP 64T | GPU vs OpenSSL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | w32-il64 | 1.72 G | none | n/a | 58.49 M | 50.74 M | 921.47 M | 6.62 G | n/a | 33.93x | 1.87x | 0.26x |
| SUBTRACT | w32-il64 | 1.70 G | none | n/a | 85.03 M | 1.01 G | 758.68 M | 6.51 G | n/a | 1.69x | 2.25x | 0.26x |
| ADDMOD | w32-il64 | 2.05 G | none | n/a | 25.69 M | 366.48 M | 95.18 M | 5.24 G | n/a | 5.59x | 21.53x | 0.39x |
| SUBTRACTMOD | w32-il64 | 2.02 G | none | n/a | 27.49 M | 410.35 M | 23.51 M | 5.26 G | n/a | 4.93x | 86.09x | 0.38x |
| MULTIPLYOPERANDSCANNING | w32-il64 | 776.90 M | none | n/a | 21.74 M | 266.88 M | 274.64 M | n/a | n/a | 2.91x | 2.83x | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 450.26 M | none | n/a | 21.81 M | 271.74 M | 292.52 M | 5.37 G | n/a | 1.66x | 1.54x | 0.08x |
| MONTGOMERYMULTIPLICATION | w32-il64 | 1.23 G | none | n/a | 2.79 M | 45.65 M | 151.78 M | 3.86 G | n/a | 27.01x | 8.12x | 0.32x |
| COMPARE | w32-il64 | 2.25 G | none | n/a | 142.82 M | 1.67 G | 1.61 G | 6.62 G | n/a | 1.35x | 1.40x | 0.34x |
| REDUCE | w32-opt | 34.71 M | none | n/a | 38.32 M | 443.68 M | 86.76 M | 2.93 G | n/a | 0.08x | 0.40x | 0.01x |
| MODMUL | w32-o64 | 5.84 M | none | n/a | 5.84 M | 90.56 M | 32.85 M | 584.77 M | n/a | 0.06x | 0.18x | 0.01x |
| MODEXP | w32-il64 | 45.44 k | none | n/a | 19.07 k | 242.66 k | 271.12 k | 2.28 M | n/a | 0.19x | 0.17x | 0.02x |
| EXPONENTIATION | w32-il64 | 27.30 k | none | n/a | 103.00 k | 205.85 k | 117.14 k | n/a | n/a | 0.13x | 0.23x | n/a |
| DIVIDE | w32-il64 | 9.60 M | none | n/a | 35.60 M | 998.28 M | 154.88 M | 2.05 G | n/a | 0.01x | 0.06x | 0.00x |
| ISQRT | w32-o64 | 254.23 k | none | n/a | 10.44 M | 328.61 M | n/a | n/a | n/a | 0.00x | n/a | n/a |
| MODMUL_R2 | w32-il64 | 827.37 M | none | n/a | 5.80 M | 135.92 M | 45.69 M | 2.38 G | n/a | 6.09x | 18.11x | 0.35x |

### p1024 (1024-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GMP 64T | OpenSSL | CGBN | GPU vs CPU-CL | GPU vs GMP 64T | GPU vs OpenSSL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | w32-il64 | 564.26 M | none | n/a | 59.27 M | 788.69 M | 741.93 M | 4.75 G | n/a | 0.72x | 0.76x | 0.12x |
| SUBTRACT | w32-il64 | 570.72 M | none | n/a | 70.74 M | 874.49 M | 509.78 M | 4.71 G | n/a | 0.65x | 1.12x | 0.12x |
| ADDMOD | w32-il | 733.36 M | none | n/a | 18.27 M | 299.84 M | 35.25 M | 4.33 G | n/a | 2.45x | 20.81x | 0.17x |
| SUBTRACTMOD | w32-il | 739.87 M | none | n/a | 22.64 M | 340.33 M | 32.47 M | 4.17 G | n/a | 2.17x | 22.79x | 0.18x |
| MULTIPLYOPERANDSCANNING | w32-il64 | 138.19 M | none | n/a | 5.96 M | 76.02 M | 70.70 M | n/a | n/a | 1.82x | 1.95x | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 68.68 M | none | n/a | 5.96 M | 74.87 M | 73.78 M | 2.26 G | n/a | 0.92x | 0.93x | 0.03x |
| MONTGOMERYMULTIPLICATION | w32-il64 | 230.37 M | none | n/a | 880.01 k | 14.04 M | 45.00 M | 1.49 G | n/a | 16.41x | 5.12x | 0.15x |
| COMPARE | w32-il | 871.02 M | none | n/a | 140.54 M | 2.00 G | 1.65 G | 4.88 G | n/a | 0.43x | 0.53x | 0.18x |
| REDUCE | w32-il64 | 5.65 M | none | n/a | 49.59 M | 690.24 M | 35.48 M | 2.01 G | n/a | 0.01x | 0.16x | 0.00x |
| MODMUL | w32-il64 | 724.15 k | none | n/a | 2.06 M | 28.42 M | 10.83 M | 238.55 M | n/a | 0.03x | 0.07x | 0.00x |
| MODEXP | w32-o64 | 2.39 k | none | n/a | 3.14 k | 28.76 k | 36.95 k | 442.40 k | n/a | 0.08x | 0.06x | 0.01x |
| EXPONENTIATION | w32-o64 | 1.93 k | none | n/a | 22.89 k | 242.53 k | 19.04 k | n/a | n/a | 0.01x | 0.10x | n/a |
| DIVIDE | w32-il | 1.05 M | none | n/a | 31.74 M | 561.75 M | 97.51 M | 1.52 G | n/a | 0.00x | 0.01x | 0.00x |
| ISQRT | w32-il64 | 17.88 k | none | n/a | 5.12 M | 96.58 M | n/a | n/a | n/a | 0.00x | n/a | n/a |
| MODMUL_R2 | w32-il64 | 124.02 M | none | n/a | 2.03 M | 45.24 M | 13.10 M | 824.97 M | n/a | 2.74x | 9.47x | 0.15x |

### p2048 (2048-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GMP 64T | OpenSSL | CGBN | GPU vs CPU-CL | GPU vs GMP 64T | GPU vs OpenSSL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | w32-il | 89.83 M | none | n/a | 44.04 M | 13.99 M | 16.64 M | 2.15 G | n/a | 6.42x | 5.40x | 0.04x |
| SUBTRACT | w32-il | 91.25 M | none | n/a | 44.98 M | 37.53 M | 757.76 M | 2.11 G | n/a | 2.43x | 0.12x | 0.04x |
| ADDMOD | w32-il | 161.13 M | none | n/a | 14.78 M | 316.27 M | 16.62 M | 2.14 G | n/a | 0.51x | 9.69x | 0.08x |
| SUBTRACTMOD | w32-il | 161.51 M | none | n/a | 17.29 M | 243.14 M | 41.04 M | 2.10 G | n/a | 0.66x | 3.94x | 0.08x |
| MULTIPLYOPERANDSCANNING | w32-il64 | 20.55 M | none | n/a | 1.82 M | 24.14 M | 22.96 M | n/a | n/a | 0.85x | 0.89x | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 7.16 M | none | n/a | 1.82 M | 25.29 M | 23.96 M | 776.59 M | n/a | 0.28x | 0.30x | 0.01x |
| MONTGOMERYMULTIPLICATION | w32 | 27.51 M | none | n/a | 267.17 k | 4.28 M | 1.81 M | 509.46 M | n/a | 6.42x | 15.17x | 0.05x |
| COMPARE | w32-il | 277.73 M | none | n/a | 122.07 M | 1.98 G | 1.99 G | 2.13 G | n/a | 0.14x | 0.14x | 0.13x |
| REDUCE | w32-il64 | 952.27 k | none | n/a | 34.13 M | 391.58 M | 45.29 M | 1.73 G | n/a | 0.00x | 0.02x | 0.00x |
| MODMUL | w32-il64 | 107.85 k | none | n/a | 635.48 k | 8.60 M | 3.96 M | 115.66 M | n/a | 0.01x | 0.03x | 0.00x |
| MODEXP | w32-opt | 79.0 | none | n/a | 433.5 | 7.22 k | 13.47 k | 79.25 k | n/a | 0.01x | 0.01x | 0.00x |
| EXPONENTIATION | w32-il64 | 119.2 | none | n/a | 3.64 k | 49.85 k | 4.84 k | n/a | n/a | 0.00x | 0.02x | n/a |
| DIVIDE | w32-il64 | 71.04 k | none | n/a | 23.77 M | 310.84 M | 60.58 M | 1.48 G | n/a | 0.00x | 0.00x | 0.00x |
| ISQRT | w32-il64 | 3.75 k | none | n/a | 3.00 M | 48.26 M | n/a | n/a | n/a | 0.00x | n/a | n/a |
| MODMUL_R2 | w32-il64 | 13.98 M | none | n/a | 618.11 k | 10.21 M | 4.58 M | 269.82 M | n/a | 1.37x | 3.06x | 0.05x |

## 6. Raw data

Also written to `NVIDIA_H100_80GB_HBM3_Report.csv` for analysis.

```csv
kind,device,device_type,kernel,modulus,bits,operation,items,seconds,ops_per_sec,mismatches
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,secp256k1,256,ADD,50000,0.000748534,66797234.029,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,secp256k1,256,ADD,50000,0.000033830,1477960260.288,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,secp256k1,256,ADD,50000,0.000040616,1231041189.605,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,secp256k1,256,ADD,50000,0.000007616,6565126050.420,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,ADD,50000,0.000028613,1747456016.665,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,ADD,50000,0.000411896,121389855.066,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,ADD,50000,0.000018878,2648565814.432,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,ADD,50000,0.000419217,119269992.552,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,secp256k1,256,ADD,50000,0.000013270,3767911794.224,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,secp256k1,256,ADD,50000,0.000383163,130492810.718,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,ADD,50000,0.000013771,3630818057.011,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,ADD,50000,0.000420369,118943060.208,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,ADD,50000,0.000013390,3734169691.700,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,ADD,50000,0.000415461,120348267.168,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,ADD,50000,0.000009685,5162592610.044,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,ADD,50000,0.000403744,123840797.939,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,ADD,50000,0.000010075,4962755703.457,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,ADD,50000,0.000412928,121086502.262,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,secp256k1,256,SUBTRACT,50000,0.000567216,88149774.401,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,secp256k1,256,SUBTRACT,50000,0.000035918,1392067229.418,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,secp256k1,256,SUBTRACT,50000,0.000052451,953277990.852,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,secp256k1,256,SUBTRACT,50000,0.000007200,6944444444.444,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,SUBTRACT,50000,0.000028823,1734723531.027,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,SUBTRACT,50000,0.000399658,125106969.341,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,SUBTRACT,50000,0.000018628,2684119699.526,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,SUBTRACT,50000,0.000410404,121831223.637,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,secp256k1,256,SUBTRACT,50000,0.000013089,3820057720.222,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,secp256k1,256,SUBTRACT,50000,0.000378195,132206943.026,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,SUBTRACT,50000,0.000012939,3864257189.642,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,SUBTRACT,50000,0.000365807,136684020.291,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,SUBTRACT,50000,0.000013080,3822641688.917,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,SUBTRACT,50000,0.000406618,122965536.062,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,SUBTRACT,50000,0.000009505,5260345992.553,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,SUBTRACT,50000,0.000411796,121419299.116,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,SUBTRACT,50000,0.000009404,5316869641.000,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,SUBTRACT,50000,0.000414590,120601109.033,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,secp256k1,256,ADDMOD,50000,0.002058368,24291084.035,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,secp256k1,256,ADDMOD,50000,0.000123725,404121914.583,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,secp256k1,256,ADDMOD,50000,0.000542096,92234668.602,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,secp256k1,256,ADDMOD,50000,0.000008544,5852059925.094,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,ADDMOD,50000,0.000039288,1272650777.227,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,ADDMOD,50000,0.000410154,121905431.773,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,ADDMOD,50000,0.000025248,1980342722.243,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,ADDMOD,50000,0.000402311,124282003.250,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,secp256k1,256,ADDMOD,50000,0.000013971,3578841176.569,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,secp256k1,256,ADDMOD,50000,0.000384555,130020424.896,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,ADDMOD,50000,0.000010676,4683409260.027,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,ADDMOD,50000,0.000406677,122947724.806,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,ADDMOD,50000,0.000010656,4692209775.603,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,ADDMOD,50000,0.000412006,121357476.243,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,ADDMOD,50000,0.000009565,5227437618.364,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,ADDMOD,50000,0.000396974,125952784.971,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,ADDMOD,50000,0.000009064,5516411025.200,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,ADDMOD,50000,0.000415842,120237960.883,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,secp256k1,256,SUBTRACTMOD,50000,0.001634431,30591680.339,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,secp256k1,256,SUBTRACTMOD,50000,0.000108796,459575649.474,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,secp256k1,256,SUBTRACTMOD,50000,0.000704824,70939696.273,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,secp256k1,256,SUBTRACTMOD,50000,0.000009216,5425347222.222,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.000039249,1273919102.110,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.000410824,121706594.366,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,SUBTRACTMOD,50000,0.000024096,2075023816.334,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,SUBTRACTMOD,50000,0.000416914,119928809.940,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,secp256k1,256,SUBTRACTMOD,50000,0.000014211,3518446216.106,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,secp256k1,256,SUBTRACTMOD,50000,0.000375081,133304508.799,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,SUBTRACTMOD,50000,0.000010746,4652966541.720,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,SUBTRACTMOD,50000,0.000417554,119744975.998,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,SUBTRACTMOD,50000,0.000010656,4692209775.603,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,SUBTRACTMOD,50000,0.000407820,122603078.602,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,SUBTRACTMOD,50000,0.000009054,5522511052.821,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,SUBTRACTMOD,50000,0.000395211,126514646.531,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,SUBTRACTMOD,50000,0.000009144,5468091686.400,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,SUBTRACTMOD,50000,0.000410604,121771880.717,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000942410,53055483.356,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000069638,718001832.660,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001269295,39391944.819,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000885764,56448450.620,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001392943,35895221.060,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000251366,198913092.425,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000715490,69882182.072,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000072198,692540068.626,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000516993,96713109.489,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000018848,2652786401.818,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000486408,102794384.201,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000018818,2657020461.997,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000484044,103296401.987,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000011257,4441721783.735,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000475020,105258744.777,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000010155,4923614380.044,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000482081,103717016.675,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000942968,53024081.400,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000922424,54204979.122,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000204979,243927177.727,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000007168,6975446428.571,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000170004,294110573.818,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000620548,80573955.373,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000056564,883956387.585,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000514890,97108151.404,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000020881,2394526998.428,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000452867,110407688.855,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000021672,2307137567.684,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000490193,102000635.900,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000020550,2433078389.341,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000481831,103770843.626,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000017396,2874233618.417,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000477214,104774772.677,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000016735,2987761767.488,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000486457,102784003.002,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.007565074,6609320.761,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000489131,102222086.993,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000156154,320196883.145,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000009664,5173841059.603,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000079389,629810468.278,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000454590,109989231.370,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000025628,1950998580.916,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000417344,119805300.149,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000012578,3975202043.612,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000383383,130417920.379,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000013841,3612434013.491,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000420669,118858268.281,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000012659,3949758410.888,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000406157,123105106.937,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000013119,3811244184.148,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000404955,123470472.250,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000010776,4639897258.173,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000421651,118581432.132,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,secp256k1,256,COMPARE,50000,0.000348144,143618725.528,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,secp256k1,256,COMPARE,50000,0.000031240,1600511507.747,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,secp256k1,256,COMPARE,50000,0.000022654,2207131373.365,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,secp256k1,256,COMPARE,50000,0.000007200,6944444444.444,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,COMPARE,50000,0.000028312,1766037260.175,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,COMPARE,50000,0.000398356,125515869.652,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,COMPARE,50000,0.000018758,2665529259.604,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,COMPARE,50000,0.000411175,121602736.153,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,COMPARE,50000,0.000010085,4957829038.439,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,COMPARE,50000,0.000413359,120960188.896,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,COMPARE,50000,0.000010084,4958401403.833,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,COMPARE,50000,0.000400990,124691384.181,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,COMPARE,50000,0.000009645,5184028118.286,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,COMPARE,50000,0.000403304,123975922.105,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,COMPARE,50000,0.000009885,5058258504.299,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,COMPARE,50000,0.000418696,119418359.927,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,secp256k1,256,REDUCE,6250,0.000100349,62282632.831,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,secp256k1,256,REDUCE,6250,0.000007646,817375501.114,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,secp256k1,256,REDUCE,6250,0.000073496,85038686.383,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,secp256k1,256,REDUCE,50000,0.000012064,4144562334.218,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,REDUCE,50000,0.000137866,362671588.117,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,REDUCE,50000,0.000524414,95344517.016,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,REDUCE,50000,0.000092338,541488454.798,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,REDUCE,50000,0.000477213,104775028.274,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,REDUCE,50000,0.000036304,1377263056.361,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,REDUCE,50000,0.000444895,112386096.570,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,REDUCE,50000,0.000035253,1418314156.831,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,REDUCE,50000,0.000423273,118127051.658,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,REDUCE,50000,0.000034301,1457680216.125,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,REDUCE,50000,0.000420178,118997237.581,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,REDUCE,50000,0.000033961,1472280903.052,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,REDUCE,50000,0.000444966,112368160.560,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,secp256k1,256,MODMUL,3125,0.000271017,11530631.109,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,secp256k1,256,MODMUL,3125,0.000019161,163090692.634,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,secp256k1,256,MODMUL,3125,0.000060652,51523453.919,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,secp256k1,256,MODMUL,50000,0.000032384,1543972332.016,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,MODMUL,50000,0.000370314,135020509.316,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,MODMUL,50000,0.000776201,64416302.803,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,MODMUL,50000,0.000219679,227604611.267,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,MODMUL,50000,0.000603933,82790624.213,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,MODMUL,50000,0.000103785,481765383.071,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,MODMUL,50000,0.000514389,97202697.707,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,MODMUL,50000,0.000082343,607217001.640,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,MODMUL,50000,0.000480198,104123751.257,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,MODMUL,50000,0.000102523,487696476.295,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,MODMUL,50000,0.000490594,101917276.617,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,MODMUL,50000,0.000083585,598194857.853,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,MODMUL,50000,0.000489231,102201168.077,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,secp256k1,256,MODEXP,781,0.007165837,108989.360,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,secp256k1,256,MODEXP,781,0.000435692,1792550.597,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,secp256k1,256,MODEXP,781,0.000881538,885951.804,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,secp256k1,256,MODEXP,50000,0.009278368,5388878.734,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,MODEXP,50000,0.016824573,2971843.616,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,MODEXP,50000,0.017286333,2892458.449,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,MODEXP,50000,0.003187254,15687484.864,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,MODEXP,50000,0.003655544,13677854.725,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,MODEXP,50000,0.002537454,19704790.608,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,MODEXP,50000,0.003005084,16638470.653,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,MODEXP,50000,0.001213585,41200244.652,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,MODEXP,50000,0.001614144,30976173.739,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,MODEXP,50000,0.002510885,19913297.427,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,MODEXP,50000,0.003002420,16653232.683,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,MODEXP,50000,0.001202629,41575583.302,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,MODEXP,50000,0.001593363,31380168.751,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,secp256k1,256,EXPONENTIATION,781,0.002223315,351277.267,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,secp256k1,256,EXPONENTIATION,781,0.000148082,5274120.205,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,secp256k1,256,EXPONENTIATION,781,0.002315824,337245.034,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,EXPONENTIATION,50000,0.025718424,1944131.569,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,EXPONENTIATION,50000,0.026090641,1916395.995,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,EXPONENTIATION,50000,0.007465570,6697412.253,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,EXPONENTIATION,50000,0.007938036,6298787.259,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,EXPONENTIATION,50000,0.000462782,108042261.725,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,EXPONENTIATION,50000,0.000862249,57987893.818,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,EXPONENTIATION,50000,0.000397655,125737150.151,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,EXPONENTIATION,50000,0.000791273,63189310.082,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,EXPONENTIATION,50000,0.000467710,106903858.873,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,EXPONENTIATION,50000,0.000871914,57345106.896,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,EXPONENTIATION,50000,0.000404485,123613967.362,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,EXPONENTIATION,50000,0.000795550,62849605.179,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,secp256k1,256,DIVIDE,6250,0.000166854,37457852.365,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,secp256k1,256,DIVIDE,6250,0.000012881,485226040.158,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,secp256k1,256,DIVIDE,6250,0.000064130,97458377.264,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,secp256k1,256,DIVIDE,50000,0.000015424,3241701244.813,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,DIVIDE,50000,0.000330755,151169284.702,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,DIVIDE,50000,0.000788349,63423689.955,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,DIVIDE,50000,0.000275843,181262631.337,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,DIVIDE,50000,0.000720828,69364668.034,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,DIVIDE,50000,0.000100260,498703858.917,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,DIVIDE,50000,0.000552506,90496762.446,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,DIVIDE,50000,0.000099009,505003903.198,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,DIVIDE,50000,0.000556782,89801771.711,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,DIVIDE,50000,0.000098818,505980785.071,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,DIVIDE,50000,0.000544454,91835127.770,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,DIVIDE,50000,0.000095884,521464250.713,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,DIVIDE,50000,0.000554199,90220284.027,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,secp256k1,256,ISQRT,1562,0.000082068,19033054.456,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,secp256k1,256,ISQRT,1562,0.000005939,263026556.713,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,ISQRT,50000,0.003645282,13716359.919,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,ISQRT,50000,0.004082546,12247259.120,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,ISQRT,50000,0.002722069,18368379.767,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,ISQRT,50000,0.003189627,15675814.877,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,ISQRT,50000,0.000662180,75508155.224,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,ISQRT,50000,0.001056440,47328756.709,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,ISQRT,50000,0.000598064,83603103.091,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,ISQRT,50000,0.000994347,50284255.090,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,ISQRT,50000,0.000678855,73653431.744,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,ISQRT,50000,0.001097271,45567595.093,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,ISQRT,50000,0.000599136,83453523.437,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,ISQRT,50000,0.001002991,49850890.684,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,secp256k1,256,MODMUL_R2,50000,0.004363566,11458518.263,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,secp256k1,256,MODMUL_R2,50000,0.000267560,186874045.325,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,secp256k1,256,MODMUL_R2,50000,0.000967767,51665323.751,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,secp256k1,256,MODMUL_R2,50000,0.000014560,3434065934.066,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,MODMUL_R2,50000,0.000093821,532929895.423,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,secp256k1,256,MODMUL_R2,50000,0.000474589,105354329.201,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,MODMUL_R2,50000,0.000028643,1745637821.492,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,secp256k1,256,MODMUL_R2,50000,0.000344935,144954801.972,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,MODMUL_R2,50000,0.000018177,2750715573.204,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,secp256k1,256,MODMUL_R2,50000,0.000406417,123026330.299,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,MODMUL_R2,50000,0.000013881,3602072608.944,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,secp256k1,256,MODMUL_R2,50000,0.000413689,120863721.978,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,MODMUL_R2,50000,0.000018538,2697166098.970,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,secp256k1,256,MODMUL_R2,50000,0.000396082,126236503.936,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,MODMUL_R2,50000,0.000012859,3888326147.495,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,secp256k1,256,MODMUL_R2,50000,0.000405135,123415621.462,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,rsa256(composite),256,ADD,50000,0.000764834,65373663.117,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,rsa256(composite),256,ADD,50000,0.000041460,1205995851.567,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,rsa256(composite),256,ADD,50000,0.000045703,1094013901.597,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,rsa256(composite),256,ADD,50000,0.000007296,6853070175.439,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,ADD,50000,0.000028423,1759136642.747,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,ADD,50000,0.000408461,122410682.159,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,ADD,50000,0.000018448,2710339944.215,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,ADD,50000,0.000391705,127647093.844,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,rsa256(composite),256,ADD,50000,0.000013540,3692753117.584,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,rsa256(composite),256,ADD,50000,0.000378706,132028529.938,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,ADD,50000,0.000013270,3767911794.224,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,ADD,50000,0.000408150,122504044.404,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,ADD,50000,0.000012939,3864257189.642,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,ADD,50000,0.000412968,121074760.104,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,ADD,50000,0.000009815,5094256074.013,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,ADD,50000,0.000390244,128124908.149,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,ADD,50000,0.000009604,5206147174.477,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,ADD,50000,0.000421180,118714110.364,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACT,50000,0.000567845,88052165.001,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACT,50000,0.000043327,1154021316.911,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACT,50000,0.000762623,65563227.307,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,rsa256(composite),256,SUBTRACT,50000,0.000007296,6853070175.439,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,SUBTRACT,50000,0.000029023,1722769326.049,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,SUBTRACT,50000,0.000402231,124306678.807,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,SUBTRACT,50000,0.000018467,2707537852.865,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,SUBTRACT,50000,0.000398096,125597867.832,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,rsa256(composite),256,SUBTRACT,50000,0.000013190,3790725049.867,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,rsa256(composite),256,SUBTRACT,50000,0.000396282,126172718.871,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,SUBTRACT,50000,0.000013199,3788117212.912,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,SUBTRACT,50000,0.000406759,122922952.509,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,SUBTRACT,50000,0.000013040,3834381401.993,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,SUBTRACT,50000,0.000413668,120869912.478,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,SUBTRACT,50000,0.000009474,5277669324.158,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,SUBTRACT,50000,0.000388040,128852704.067,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,SUBTRACT,50000,0.000009384,5328214688.368,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,SUBTRACT,50000,0.000414860,120522594.889,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,rsa256(composite),256,ADDMOD,50000,0.001747577,28611047.923,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,rsa256(composite),256,ADDMOD,50000,0.000388418,128727397.235,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,rsa256(composite),256,ADDMOD,50000,0.000498729,100254789.064,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,rsa256(composite),256,ADDMOD,50000,0.000008768,5702554744.526,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,ADDMOD,50000,0.000036855,1356668192.127,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,ADDMOD,50000,0.000415882,120226382.712,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,ADDMOD,50000,0.000022604,2211985134.523,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,ADDMOD,50000,0.000391535,127702505.783,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,rsa256(composite),256,ADDMOD,50000,0.000013790,3625791260.890,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,rsa256(composite),256,ADDMOD,50000,0.000389482,128375595.359,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,ADDMOD,50000,0.000010726,4661551723.539,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,ADDMOD,50000,0.000406398,123032109.935,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,ADDMOD,50000,0.000011156,4481954435.029,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,ADDMOD,50000,0.000406959,122862541.543,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,ADDMOD,50000,0.000009184,5444248061.858,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,ADDMOD,50000,0.000391966,127562095.919,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,ADDMOD,50000,0.000008753,5712304218.758,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,ADDMOD,50000,0.000418076,119595463.413,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACTMOD,50000,0.001634381,30592621.673,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACTMOD,50000,0.000113557,440308655.403,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACTMOD,50000,0.000668510,74793194.737,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,rsa256(composite),256,SUBTRACTMOD,50000,0.000009120,5482456140.351,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,SUBTRACTMOD,50000,0.000037967,1316933314.527,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,SUBTRACTMOD,50000,0.000416643,120006753.288,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,SUBTRACTMOD,50000,0.000024226,2063895865.449,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,SUBTRACTMOD,50000,0.000400239,124925387.301,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,rsa256(composite),256,SUBTRACTMOD,50000,0.000014441,3462343041.403,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,rsa256(composite),256,SUBTRACTMOD,50000,0.000382252,130803750.122,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,SUBTRACTMOD,50000,0.000010576,4727641000.352,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,SUBTRACTMOD,50000,0.000401550,124517503.206,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,SUBTRACTMOD,50000,0.000010846,4610015773.995,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,SUBTRACTMOD,50000,0.000409042,122236837.198,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,SUBTRACTMOD,50000,0.000008943,5590949356.938,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,SUBTRACTMOD,50000,0.000385396,129736703.352,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,SUBTRACTMOD,50000,0.000009134,5474085261.280,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,SUBTRACTMOD,50000,0.000411546,121493074.870,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000944320,52948166.445,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000068645,728382156.112,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000119035,420042825.721,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000915168,54634786.201,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001397390,35780993.791,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000251777,198588429.905,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000709450,70477142.850,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000073771,677773177.968,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000511295,97790917.010,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000018978,2634625994.357,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000476593,104911336.324,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000019810,2523986751.760,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000492185,101587792.531,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000010966,4559510070.278,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000450583,110967303.694,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000010105,4948119004.608,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000476532,104924714.930,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000945163,52900901.151,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000069063,723973794.478,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000118991,420201453.253,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000007264,6883259911.894,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000169554,294891256.459,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000627628,79665016.645,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000056815,880047720.874,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000489241,102199076.656,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000021312,2346079257.115,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000468110,106812508.822,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000021932,2279755035.139,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000474880,105289761.001,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000021141,2365070096.916,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000494650,101081600.090,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000017566,2846385027.702,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000471665,106007480.952,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000016805,2975301893.955,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000475060,105249871.616,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.007549982,6622532.299,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000464855,107560423.273,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000259228,192880360.992,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000009792,5106209150.327,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000078667,635590888.857,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000459237,108876238.993,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000025658,1948714744.102,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000396092,126233313.152,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000011989,4170438016.818,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000382943,130567787.116,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000014231,3513438120.480,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000407479,122705708.502,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000012549,3984347560.206,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000421540,118612739.464,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000013470,3711901766.516,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000387089,129169229.400,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000010687,4678511683.841,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000405487,123308544.041,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,rsa256(composite),256,COMPARE,50000,0.000404894,123489230.622,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,rsa256(composite),256,COMPARE,50000,0.000043840,1140521647.505,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,rsa256(composite),256,COMPARE,50000,0.000024130,2072092629.759,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,rsa256(composite),256,COMPARE,50000,0.000007392,6764069264.069,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,COMPARE,50000,0.000028483,1755441007.741,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,COMPARE,50000,0.000406127,123114211.186,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,COMPARE,50000,0.000017977,2781354290.895,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,COMPARE,50000,0.000388811,128597148.505,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,COMPARE,50000,0.000009815,5094256074.013,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,COMPARE,50000,0.000404164,123712168.012,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,COMPARE,50000,0.000010085,4957829038.439,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,COMPARE,50000,0.000420849,118807437.983,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,COMPARE,50000,0.000009464,5283252510.640,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,COMPARE,50000,0.000393820,126961557.580,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,COMPARE,50000,0.000009515,5254811089.632,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,COMPARE,50000,0.000397294,125851438.966,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,rsa256(composite),256,REDUCE,6250,0.000159372,39216378.321,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,rsa256(composite),256,REDUCE,6250,0.000010509,594752195.684,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,rsa256(composite),256,REDUCE,6250,0.000078495,79623385.199,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,rsa256(composite),256,REDUCE,50000,0.000012288,4069010416.667,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,REDUCE,50000,0.000139308,358917114.669,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,REDUCE,50000,0.000515851,96927221.117,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,REDUCE,50000,0.000093129,536889703.140,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,REDUCE,50000,0.000477264,104763834.322,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,REDUCE,50000,0.000036945,1353359412.143,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,REDUCE,50000,0.000430684,116094379.530,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,REDUCE,50000,0.000034531,1447969555.660,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,REDUCE,50000,0.000442242,113060265.819,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,REDUCE,50000,0.000035202,1420377964.297,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,REDUCE,50000,0.000421180,118714044.739,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,REDUCE,50000,0.000034301,1457680216.125,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,REDUCE,50000,0.000422332,118390279.519,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,rsa256(composite),256,MODMUL,3125,0.000271350,11516475.566,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,rsa256(composite),256,MODMUL,3125,0.000017982,173781443.318,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,rsa256(composite),256,MODMUL,3125,0.000061146,51107199.756,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,rsa256(composite),256,MODMUL,50000,0.000032224,1551638530.288,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,MODMUL,50000,0.000372687,134160877.328,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,MODMUL,50000,0.000774167,64585547.043,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,MODMUL,50000,0.000218967,228344935.648,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,MODMUL,50000,0.000600648,83243414.824,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,MODMUL,50000,0.000103755,481904845.340,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,MODMUL,50000,0.000496433,100718502.329,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,MODMUL,50000,0.000082203,608250645.926,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,MODMUL,50000,0.000490483,101940353.782,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,MODMUL,50000,0.000103374,483680560.736,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,MODMUL,50000,0.000488030,102452712.668,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,MODMUL,50000,0.000083715,597264839.398,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,MODMUL,50000,0.000471675,106005178.522,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,rsa256(composite),256,MODEXP,781,0.006721442,116195.307,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,rsa256(composite),256,MODEXP,781,0.000430023,1816181.651,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,rsa256(composite),256,MODEXP,781,0.000882320,885166.507,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,rsa256(composite),256,MODEXP,50000,0.009120096,5482398.431,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,MODEXP,50000,0.016840727,2968992.954,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,MODEXP,50000,0.017339813,2883537.431,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,MODEXP,50000,0.003179123,15727608.260,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,MODEXP,50000,0.003648143,13705602.989,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,MODEXP,50000,0.002536413,19712877.780,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,MODEXP,50000,0.002935249,17034329.238,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,MODEXP,50000,0.001210601,41301797.329,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,MODEXP,50000,0.001620223,30859948.626,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,MODEXP,50000,0.002501741,19986081.280,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,MODEXP,50000,0.002991203,16715682.699,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,MODEXP,50000,0.001194316,41864968.375,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,MODEXP,50000,0.001584329,31559102.795,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,rsa256(composite),256,EXPONENTIATION,781,0.002218819,351989.129,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,rsa256(composite),256,EXPONENTIATION,781,0.000147650,5289518.726,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,rsa256(composite),256,EXPONENTIATION,781,0.002325207,335884.007,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,EXPONENTIATION,50000,0.025700177,1945511.905,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,EXPONENTIATION,50000,0.026207725,1907834.423,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,EXPONENTIATION,50000,0.007453922,6707878.126,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,EXPONENTIATION,50000,0.007921671,6311799.596,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,EXPONENTIATION,50000,0.000462682,108065586.019,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,EXPONENTIATION,50000,0.000858904,58213717.072,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,EXPONENTIATION,50000,0.000391906,127581572.506,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,EXPONENTIATION,50000,0.000806697,61981144.543,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,EXPONENTIATION,50000,0.000467520,106947302.243,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,EXPONENTIATION,50000,0.000865434,57774476.535,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,EXPONENTIATION,50000,0.000402952,124084305.628,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,EXPONENTIATION,50000,0.000812986,61501679.044,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,rsa256(composite),256,DIVIDE,6250,0.000168867,37011330.350,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,rsa256(composite),256,DIVIDE,6250,0.000012246,510384669.694,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,rsa256(composite),256,DIVIDE,6250,0.000068235,91594689.029,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,rsa256(composite),256,DIVIDE,50000,0.000015328,3262004175.365,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,DIVIDE,50000,0.000332347,150445219.519,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,DIVIDE,50000,0.000785044,63690703.370,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,DIVIDE,50000,0.000280990,177942254.958,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,DIVIDE,50000,0.000728098,68672076.522,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,DIVIDE,50000,0.000104817,477021419.591,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,DIVIDE,50000,0.000557844,89630782.972,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,DIVIDE,50000,0.000100049,499755332.820,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,DIVIDE,50000,0.000575471,86885355.691,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,DIVIDE,50000,0.000103735,481996783.222,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,DIVIDE,50000,0.000562291,88921908.289,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,DIVIDE,50000,0.000099128,504398967.472,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,DIVIDE,50000,0.000567248,88144853.697,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,rsa256(composite),256,ISQRT,1562,0.000081734,19110890.379,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,rsa256(composite),256,ISQRT,1562,0.000005689,274557297.142,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,ISQRT,50000,0.003651080,13694577.340,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,ISQRT,50000,0.004106611,12175489.677,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,ISQRT,50000,0.002795759,17884230.294,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,ISQRT,50000,0.003269868,15291136.179,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,ISQRT,50000,0.000690042,72459353.216,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,ISQRT,50000,0.001080676,46267337.144,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,ISQRT,50000,0.000644984,77521308.071,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,ISQRT,50000,0.001054827,47401142.713,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,ISQRT,50000,0.000706266,70794841.679,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,ISQRT,50000,0.001124492,44464526.004,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,ISQRT,50000,0.000651624,76731360.230,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,ISQRT,50000,0.001059144,47207932.312,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,rsa256(composite),256,MODMUL_R2,50000,0.004359660,11468784.338,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,rsa256(composite),256,MODMUL_R2,50000,0.000279308,179013807.551,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,rsa256(composite),256,MODMUL_R2,50000,0.001129779,44256441.767,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,rsa256(composite),256,MODMUL_R2,50000,0.000014656,3411572052.402,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,MODMUL_R2,50000,0.000094291,530272990.318,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,rsa256(composite),256,MODMUL_R2,50000,0.000409202,122189055.533,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,MODMUL_R2,50000,0.000028582,1749349251.786,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,rsa256(composite),256,MODMUL_R2,50000,0.000407710,122636195.628,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,MODMUL_R2,50000,0.000018728,2669804127.505,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,rsa256(composite),256,MODMUL_R2,50000,0.000415832,120240855.774,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,MODMUL_R2,50000,0.000013991,3573778745.216,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,rsa256(composite),256,MODMUL_R2,50000,0.000425226,117584516.287,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,MODMUL_R2,50000,0.000018537,2697335487.031,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,rsa256(composite),256,MODMUL_R2,50000,0.000418636,119435495.332,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,MODMUL_R2,50000,0.000013129,3808337881.502,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,rsa256(composite),256,MODMUL_R2,50000,0.000418956,119344295.938,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,brainpoolP512r1,512,ADD,25000,0.000427455,58485687.271,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,brainpoolP512r1,512,ADD,25000,0.000492751,50735550.547,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,brainpoolP512r1,512,ADD,25000,0.000027131,921470191.277,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,brainpoolP512r1,512,ADD,50000,0.000007552,6620762711.864,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,ADD,50000,0.000059048,846769494.773,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,ADD,50000,0.000702941,71129715.688,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,ADD,50000,0.000033691,1484073231.883,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,ADD,50000,0.000659145,75855849.405,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,brainpoolP512r1,512,ADD,50000,0.000019981,2502369721.970,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,brainpoolP512r1,512,ADD,50000,0.000648190,77137872.128,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,ADD,50000,0.000019449,2570820691.224,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,ADD,50000,0.000739916,67575238.734,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,ADD,50000,0.000019819,2522830347.031,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,ADD,50000,0.000780567,64056014.406,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,ADD,50000,0.000014582,3428841845.761,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,ADD,50000,0.000670623,74557534.482,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,ADD,50000,0.000014521,3443301181.715,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,ADD,50000,0.000774068,64593803.339,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACT,25000,0.000294025,85026885.052,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACT,25000,0.000024799,1008085232.788,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACT,25000,0.000032952,758684717.907,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,brainpoolP512r1,512,SUBTRACT,50000,0.000007680,6510416666.667,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,SUBTRACT,50000,0.000060640,824537678.683,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,SUBTRACT,50000,0.000694068,72039059.658,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,SUBTRACT,50000,0.000033020,1514231877.027,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,SUBTRACT,50000,0.000667488,74907717.106,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,brainpoolP512r1,512,SUBTRACT,50000,0.000020000,2500010067.638,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,brainpoolP512r1,512,SUBTRACT,50000,0.000651914,76697241.557,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,SUBTRACT,50000,0.000019820,2522711801.330,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,SUBTRACT,50000,0.000731254,68375694.683,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,SUBTRACT,50000,0.000019269,2594864181.539,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,SUBTRACT,50000,0.000785505,63653323.959,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,SUBTRACT,50000,0.000014883,3359537636.494,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,SUBTRACT,50000,0.000678094,73736078.385,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,SUBTRACT,50000,0.000014671,3408055049.832,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,SUBTRACT,50000,0.000769841,64948474.765,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,brainpoolP512r1,512,ADDMOD,25000,0.000973005,25693591.957,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,brainpoolP512r1,512,ADDMOD,25000,0.000068217,366477294.106,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,brainpoolP512r1,512,ADDMOD,25000,0.000262665,95178191.685,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,brainpoolP512r1,512,ADDMOD,50000,0.000009536,5243288590.604,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,ADDMOD,50000,0.000075542,661884311.296,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,ADDMOD,50000,0.000708790,70542776.118,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,ADDMOD,50000,0.000041112,1216187822.738,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,ADDMOD,50000,0.000677453,73805844.918,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,brainpoolP512r1,512,ADDMOD,50000,0.000021892,2283925348.308,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,brainpoolP512r1,512,ADDMOD,50000,0.000657262,76073185.548,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,ADDMOD,50000,0.000024727,2022093622.471,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,ADDMOD,50000,0.000734067,68113668.321,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,ADDMOD,50000,0.000023936,2088910594.918,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,ADDMOD,50000,0.000779826,64116870.285,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,ADDMOD,50000,0.000012228,4088965227.822,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,ADDMOD,50000,0.000673137,74279068.950,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,ADDMOD,50000,0.000012198,4099033494.942,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,ADDMOD,50000,0.000770763,64870781.494,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000909508,27487394.876,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000060923,410350983.751,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACTMOD,25000,0.001063192,23514086.200,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000009504,5260942760.943,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000088633,564123729.366,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000728850,68601219.210,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000047582,1050818224.434,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000674337,74146912.064,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000024566,2035336601.270,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000642871,77776099.033,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000024707,2023713340.118,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000739025,67656714.244,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000023245,2150989761.311,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000775840,64446266.500,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000012499,4000304841.384,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000665986,75076655.615,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000012349,4048876577.613,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000770652,64880130.154,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.001150213,21735098.344,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.000093673,266884854.992,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.000091029,274636753.903,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.004791656,10434805.743,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.005640926,8863792.817,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.001290780,38736274.186,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.002063475,24230968.621,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000345496,144719477.806,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.001109359,45071069.023,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000061042,819109385.024,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000894327,55907965.234,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000060201,830552035.303,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000952074,52516922.267,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000033050,1512855777.779,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000806787,61974240.106,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000032179,1553805603.149,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000917492,54496389.837,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.001146399,21807407.070,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000092000,271739124.405,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000085465,292516711.447,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000009312,5369415807.560,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.001136670,43988145.089,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.001878409,26618271.838,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000314771,158845654.697,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.001070099,46724650.960,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000089554,558321628.986,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000843271,59292909.701,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000089384,559384749.714,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000921929,54234109.992,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000090716,551171044.756,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000969740,51560219.265,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000060341,828625974.487,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000849300,58872007.799,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000055523,900529480.981,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000934758,53489767.448,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.008960721,2789954.041,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000547699,45645515.355,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000164716,151776355.078,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000012960,3858024691.358,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000270605,184771101.322,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000908919,55010407.818,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000084026,595053783.296,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000712555,70170031.630,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000026740,1869864818.411,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000665676,75111606.659,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000038408,1301812942.453,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000750793,66596239.447,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000027662,1807524449.532,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000784463,63737867.764,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000032309,1547557505.441,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000681860,73328845.858,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000020280,2465481444.743,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000621760,80416877.075,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,brainpoolP512r1,512,COMPARE,25000,0.000175041,142823846.676,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,brainpoolP512r1,512,COMPARE,25000,0.000014927,1674782879.416,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,brainpoolP512r1,512,COMPARE,25000,0.000015529,1609870585.660,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,brainpoolP512r1,512,COMPARE,50000,0.000007552,6620762711.864,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,COMPARE,50000,0.000059529,839927113.719,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,COMPARE,50000,0.000681790,73336358.377,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,COMPARE,50000,0.000033671,1484955778.061,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,COMPARE,50000,0.000651083,76795129.703,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,COMPARE,50000,0.000019419,2574796949.787,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,COMPARE,50000,0.000739546,67609044.160,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,COMPARE,50000,0.000018949,2638672541.623,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,COMPARE,50000,0.000783492,63816870.724,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,COMPARE,50000,0.000011387,4390953540.393,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,COMPARE,50000,0.000666177,75055112.996,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,COMPARE,50000,0.000011087,4509814876.727,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,COMPARE,50000,0.000697533,71681206.563,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,brainpoolP512r1,512,REDUCE,3125,0.000081549,38320559.581,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,brainpoolP512r1,512,REDUCE,3125,0.000007043,443683968.166,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,brainpoolP512r1,512,REDUCE,3125,0.000036020,86757434.788,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,brainpoolP512r1,512,REDUCE,50000,0.000017088,2926029962.547,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,REDUCE,50000,0.000419056,119315783.251,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,REDUCE,50000,0.001019455,49045812.404,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,REDUCE,50000,0.000323994,154323915.826,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,REDUCE,50000,0.000951613,52542376.710,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,REDUCE,50000,0.000090034,555345881.099,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,REDUCE,50000,0.000799956,62503435.110,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,REDUCE,50000,0.000092218,542193900.619,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,REDUCE,50000,0.000852645,58641058.502,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,REDUCE,50000,0.000098728,506441381.587,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,REDUCE,50000,0.000765725,65297610.317,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,REDUCE,50000,0.000107941,463216081.000,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,REDUCE,50000,0.000814488,61388246.755,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL,1562,0.000267452,5840300.281,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL,1562,0.000017248,90558774.599,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL,1562,0.000047552,32847920.871,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,brainpoolP512r1,512,MODMUL,50000,0.000085504,584767964.072,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,MODMUL,50000,0.001155428,43274002.957,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,MODMUL,50000,0.001764990,28328773.626,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,MODMUL,50000,0.000817161,61187448.656,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,MODMUL,50000,0.001453452,34400865.290,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,MODMUL,50000,0.000331625,150772766.194,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,MODMUL,50000,0.001046885,47760736.084,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,MODMUL,50000,0.000267470,186936999.619,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,MODMUL,50000,0.001023521,48850978.653,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,MODMUL,50000,0.000399919,125025319.523,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,MODMUL,50000,0.001061207,47116164.840,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,MODMUL,50000,0.000330414,151325341.603,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,MODMUL,50000,0.001052694,47497186.597,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,brainpoolP512r1,512,MODEXP,390,0.020450515,19070.424,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,brainpoolP512r1,512,MODEXP,390,0.001607213,242656.051,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,brainpoolP512r1,512,MODEXP,390,0.001438502,271115.371,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,brainpoolP512r1,512,MODEXP,50000,0.021893024,2283832.512,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,MODEXP,50000,0.224496395,222720.726,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,MODEXP,50000,0.225280538,221945.493,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,MODEXP,50000,0.022337801,2238358.216,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,MODEXP,50000,0.023070787,2167242.922,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,MODEXP,50000,0.019464444,2568786.432,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,MODEXP,50000,0.020205172,2474613.924,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,MODEXP,50000,0.008664627,5770588.948,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,MODEXP,50000,0.009491404,5267924.505,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,MODEXP,50000,0.020815840,2402016.927,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,MODEXP,50000,0.021553344,2319825.628,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,MODEXP,50000,0.008581995,5826151.123,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,MODEXP,50000,0.009293810,5379924.918,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,brainpoolP512r1,512,EXPONENTIATION,390,0.003786318,103002.455,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,brainpoolP512r1,512,EXPONENTIATION,390,0.001894549,205853.792,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,brainpoolP512r1,512,EXPONENTIATION,390,0.003329425,117137.360,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,EXPONENTIATION,50000,0.193571014,258303.136,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,EXPONENTIATION,50000,0.194375737,257233.751,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,EXPONENTIATION,50000,0.050410795,991851.050,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,EXPONENTIATION,50000,0.051080065,978855.453,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,EXPONENTIATION,50000,0.016673311,2998804.516,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,EXPONENTIATION,50000,0.017367769,2878895.954,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,EXPONENTIATION,50000,0.015638607,3197215.693,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,EXPONENTIATION,50000,0.016478632,3034232.476,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,EXPONENTIATION,50000,0.014888426,3358313.404,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,EXPONENTIATION,50000,0.015624517,3200098.947,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,EXPONENTIATION,50000,0.014286456,3499818.249,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,EXPONENTIATION,50000,0.014944030,3345817.701,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,brainpoolP512r1,512,DIVIDE,3125,0.000087792,35595341.771,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,brainpoolP512r1,512,DIVIDE,3125,0.000003130,998277731.912,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,brainpoolP512r1,512,DIVIDE,3125,0.000020177,154882288.935,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,brainpoolP512r1,512,DIVIDE,50000,0.000024352,2053219448.095,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,DIVIDE,50000,0.001154526,43307819.970,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,DIVIDE,50000,0.001891048,26440366.368,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,DIVIDE,50000,0.001127225,44356712.229,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,DIVIDE,50000,0.001899920,26316895.183,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,DIVIDE,50000,0.000384175,130148946.558,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,DIVIDE,50000,0.001115818,44810173.437,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,DIVIDE,50000,0.000353800,141322821.161,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,DIVIDE,50000,0.001239153,40350147.018,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,DIVIDE,50000,0.000336503,148587021.664,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,DIVIDE,50000,0.001084953,46084952.070,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,DIVIDE,50000,0.000325617,153554676.464,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,DIVIDE,50000,0.001063541,47012760.159,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,brainpoolP512r1,512,ISQRT,781,0.000074825,10437643.717,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,brainpoolP512r1,512,ISQRT,781,0.000002377,328605412.247,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,ISQRT,50000,0.017183209,2909817.357,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,ISQRT,50000,0.017903637,2792728.629,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,ISQRT,50000,0.017063346,2930257.642,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,ISQRT,50000,0.017818916,2806006.847,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,ISQRT,50000,0.003361356,14874949.369,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,ISQRT,50000,0.004050477,12344224.603,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,ISQRT,50000,0.003072024,16275914.033,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,ISQRT,50000,0.003913102,12777586.770,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,ISQRT,50000,0.004154333,12035625.993,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,ISQRT,50000,0.004849122,10311145.187,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,ISQRT,50000,0.003902977,12810734.397,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,ISQRT,50000,0.004586969,10900444.193,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL_R2,25000,0.004307322,5804070.410,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL_R2,25000,0.000183925,135924945.313,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL_R2,25000,0.000547144,45691851.129,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,brainpoolP512r1,512,MODMUL_R2,50000,0.000020992,2381859756.098,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,MODMUL_R2,50000,0.000401250,124610564.967,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,brainpoolP512r1,512,MODMUL_R2,50000,0.001038232,48158798.270,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,MODMUL_R2,50000,0.000082793,603914477.747,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,brainpoolP512r1,512,MODMUL_R2,50000,0.000721969,69255056.561,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,MODMUL_R2,50000,0.000055663,898261875.342,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,brainpoolP512r1,512,MODMUL_R2,50000,0.000667088,74952633.693,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,MODMUL_R2,50000,0.000034191,1462375397.858,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,brainpoolP512r1,512,MODMUL_R2,50000,0.000793687,62997139.160,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,MODMUL_R2,50000,0.000054712,913873383.634,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,brainpoolP512r1,512,MODMUL_R2,50000,0.000680668,73457245.779,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,MODMUL_R2,50000,0.000030216,1654749029.489,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,brainpoolP512r1,512,MODMUL_R2,50000,0.000642171,77860865.580,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,p1024,1024,ADD,12500,0.000210898,59270305.483,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,p1024,1024,ADD,12500,0.000015849,788688346.054,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,p1024,1024,ADD,12500,0.000016848,741932205.238,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,p1024,1024,ADD,50000,0.000010528,4749240121.581,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,ADD,50000,0.000191126,261607500.244,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,ADD,50000,0.001267455,39449130.680,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,ADD,50000,0.000098026,510069319.602,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,ADD,50000,0.001253183,39898402.714,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,p1024,1024,ADD,50000,0.000044426,1125468349.335,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,p1024,1024,ADD,50000,0.001157561,43194264.804,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,ADD,50000,0.000043505,1149296581.251,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,ADD,50000,0.001132864,44135922.989,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,ADD,50000,0.000043445,1150879524.103,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,ADD,50000,0.001190791,41988899.751,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,ADD,50000,0.000022594,2212988095.631,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,ADD,50000,0.001105283,45237286.340,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,ADD,50000,0.000022153,2257040388.456,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,ADD,50000,0.001147997,43554125.829,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,p1024,1024,SUBTRACT,12500,0.000176704,70739657.297,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,p1024,1024,SUBTRACT,12500,0.000014294,874486174.258,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,p1024,1024,SUBTRACT,12500,0.000024520,509782462.271,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,p1024,1024,SUBTRACT,50000,0.000010624,4706325301.205,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,SUBTRACT,50000,0.000191497,261100807.929,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,SUBTRACT,50000,0.001263871,39560996.668,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,SUBTRACT,50000,0.000097876,510850730.064,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,SUBTRACT,50000,0.001261586,39632657.337,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,p1024,1024,SUBTRACT,50000,0.000044847,1114904082.735,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,p1024,1024,SUBTRACT,50000,0.001213705,41196166.377,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,SUBTRACT,50000,0.000042924,1164844296.423,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,SUBTRACT,50000,0.001149508,43496872.313,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,SUBTRACT,50000,0.000042955,1164010866.714,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,SUBTRACT,50000,0.001131092,44205070.437,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,SUBTRACT,50000,0.000022163,2256020809.127,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,SUBTRACT,50000,0.001101898,45376249.584,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,SUBTRACT,50000,0.000021902,2282881340.293,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,SUBTRACT,50000,0.001130841,44214872.699,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,p1024,1024,ADDMOD,12500,0.000684060,18273255.001,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,p1024,1024,ADDMOD,12500,0.000041689,299839383.720,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,p1024,1024,ADDMOD,12500,0.000354634,35247640.532,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,p1024,1024,ADDMOD,50000,0.000011552,4328254847.645,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,ADDMOD,50000,0.000266358,187717266.186,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,ADDMOD,50000,0.001356619,36856326.988,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,ADDMOD,50000,0.000135402,369270866.843,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,ADDMOD,50000,0.001294365,38628975.676,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,p1024,1024,ADDMOD,50000,0.000057366,871596748.179,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,p1024,1024,ADDMOD,50000,0.001145543,43647420.626,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,ADDMOD,50000,0.000059769,836553741.634,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,ADDMOD,50000,0.001149669,43490785.311,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,ADDMOD,50000,0.000060320,828913833.992,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,ADDMOD,50000,0.001222448,40901528.900,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,ADDMOD,50000,0.000017045,2933440310.353,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,ADDMOD,50000,0.001098933,45498672.382,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,ADDMOD,50000,0.000017335,2884348042.389,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,ADDMOD,50000,0.001137892,43940900.584,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,p1024,1024,SUBTRACTMOD,12500,0.000552066,22642237.921,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,p1024,1024,SUBTRACTMOD,12500,0.000036729,340329289.031,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,p1024,1024,SUBTRACTMOD,12500,0.000385006,32467038.412,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,p1024,1024,SUBTRACTMOD,50000,0.000012000,4166666666.667,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,SUBTRACTMOD,50000,0.000264345,189146780.605,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,SUBTRACTMOD,50000,0.001375367,36353934.922,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,SUBTRACTMOD,50000,0.000134772,370997156.041,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,SUBTRACTMOD,50000,0.001291871,38703559.619,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,p1024,1024,SUBTRACTMOD,50000,0.000056094,891364243.051,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,p1024,1024,SUBTRACTMOD,50000,0.001151973,43403798.220,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,SUBTRACTMOD,50000,0.000059970,833750819.392,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,SUBTRACTMOD,50000,0.001163480,42974520.828,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,SUBTRACTMOD,50000,0.000059819,835856939.125,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,SUBTRACTMOD,50000,0.001146775,43600532.996,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,SUBTRACTMOD,50000,0.000016895,2959474729.545,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,SUBTRACTMOD,50000,0.001098983,45496609.560,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,SUBTRACTMOD,50000,0.000017105,2923138430.545,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,SUBTRACTMOD,50000,0.001140606,43836350.828,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.002096981,5960950.406,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.000164431,76019757.429,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.000176815,70695559.001,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.022783855,2194536.440,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.024361936,2052382.051,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.005854091,8541035.796,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.007353171,6799787.459,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.001283179,38965723.837,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.002659326,18801756.864,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.000149874,333613531.685,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.001490699,33541314.202,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.000149594,334238179.863,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.001486563,33634634.838,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.000091417,546944257.869,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.001426904,35040899.821,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.000090456,552755729.785,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.001525211,32782351.283,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.002098093,5957792.397,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.000166960,74868256.942,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.000169429,73777420.605,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000022080,2264492753.623,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.008477308,5898098.883,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.009986335,5006841.760,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.002177254,22964706.245,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.003593151,13915362.898,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000558064,89595444.720,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.001983426,25208908.585,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000560088,89271709.736,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.001903606,26265937.926,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000559446,89374103.828,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.001906101,26231556.821,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000198907,251373775.224,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.001526183,32761471.301,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000181992,274737241.476,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.001638440,30516834.024,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.014204374,880010.633,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.000890291,14040352.103,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.000277796,44997038.198,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000033632,1486679352.997,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.001829155,27335025.950,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.003028019,16512445.388,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000272858,183245526.813,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.001375777,36343100.572,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000074772,668700126.735,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.001244531,40175776.638,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000131336,380702811.638,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.001225092,40813261.950,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000094271,530385621.809,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.001176820,42487384.722,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000094531,528927739.670,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.001167927,42810897.156,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000054261,921473015.546,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.001200715,41641852.290,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,p1024,1024,COMPARE,12500,0.000088942,140541019.835,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,p1024,1024,COMPARE,12500,0.000006235,2004761510.644,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,p1024,1024,COMPARE,12500,0.000007557,1654187517.717,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,p1024,1024,COMPARE,50000,0.000010240,4882812500.000,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,COMPARE,50000,0.000141071,354431065.397,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,COMPARE,50000,0.001235608,40465907.905,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,COMPARE,50000,0.000072869,686162778.541,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,COMPARE,50000,0.001233044,40550050.719,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,COMPARE,50000,0.000040280,1241305677.391,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,COMPARE,50000,0.001129128,44281957.651,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,COMPARE,50000,0.000040551,1233009685.013,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,COMPARE,50000,0.001129329,44274078.901,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,COMPARE,50000,0.000014351,3484082041.631,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,COMPARE,50000,0.001172784,42633596.006,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,COMPARE,50000,0.000014752,3389336565.657,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,COMPARE,50000,0.001158032,43176704.709,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,p1024,1024,REDUCE,1562,0.000031497,49592810.961,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,p1024,1024,REDUCE,1562,0.000002263,690240664.191,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,p1024,1024,REDUCE,1562,0.000044022,35482223.272,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,p1024,1024,REDUCE,50000,0.000024864,2010939510.940,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,REDUCE,50000,0.002202374,22702775.098,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,REDUCE,50000,0.003349249,14928719.760,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,REDUCE,50000,0.001194015,41875523.897,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,REDUCE,50000,0.002321520,21537614.318,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,REDUCE,50000,0.000324686,153995019.630,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,REDUCE,50000,0.001425181,35083261.881,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,REDUCE,50000,0.000319437,156525326.173,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,REDUCE,50000,0.001404670,35595547.431,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,REDUCE,50000,0.000316803,157826733.076,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,REDUCE,50000,0.001468686,34044038.467,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,REDUCE,50000,0.000276674,180718068.242,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,REDUCE,50000,0.001425031,35086953.371,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,p1024,1024,MODMUL,781,0.000379676,2057018.175,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,p1024,1024,MODMUL,781,0.000027480,28420573.570,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,p1024,1024,MODMUL,781,0.000072117,10829640.043,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,p1024,1024,MODMUL,50000,0.000209600,238549618.321,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,MODMUL,50000,0.008461335,5909233.016,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,MODMUL,50000,0.009643683,5184741.153,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,MODMUL,50000,0.003068767,16293188.644,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,MODMUL,50000,0.004244204,11780772.328,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,MODMUL,50000,0.001443659,34634213.122,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,MODMUL,50000,0.002560168,19529968.048,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,MODMUL,50000,0.001104812,45256562.894,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,MODMUL,50000,0.002197706,22750994.356,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,MODMUL,50000,0.001413684,35368580.438,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,MODMUL,50000,0.002509844,19921556.689,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,MODMUL,50000,0.001078504,46360508.275,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,MODMUL,50000,0.002231337,22408088.552,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,p1024,1024,MODEXP,195,0.062050672,3142.593,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,p1024,1024,MODEXP,195,0.006781091,28756.435,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,p1024,1024,MODEXP,195,0.005277162,36951.680,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,p1024,1024,MODEXP,50000,0.113018684,442404.727,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,MODEXP,50000,2.668123824,18739.760,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,MODEXP,50000,2.667858135,18741.626,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,MODEXP,50000,0.156429298,319633.219,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,MODEXP,50000,0.157483613,317493.351,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,MODEXP,50000,0.155928718,320659.341,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,MODEXP,50000,0.157563482,317332.413,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,MODEXP,50000,0.081624613,612560.330,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,MODEXP,50000,0.083316113,600124.011,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,MODEXP,50000,0.157319575,317824.403,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,MODEXP,50000,0.158680539,315098.501,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,MODEXP,50000,0.082283773,607653.225,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,MODEXP,50000,0.083577937,598244.006,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,p1024,1024,EXPONENTIATION,195,0.008518220,22892.106,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,p1024,1024,EXPONENTIATION,195,0.000804013,242533.430,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,p1024,1024,EXPONENTIATION,195,0.010239975,19043.015,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,EXPONENTIATION,50000,1.611529812,31026.420,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,EXPONENTIATION,50000,1.613971929,30979.473,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,EXPONENTIATION,50000,0.387889197,128902.791,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,EXPONENTIATION,50000,0.389160417,128481.721,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,EXPONENTIATION,50000,0.112911325,442825.377,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,EXPONENTIATION,50000,0.114532109,436558.799,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,EXPONENTIATION,50000,0.101085436,494631.096,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,EXPONENTIATION,50000,0.102582194,487414.024,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,EXPONENTIATION,50000,0.110324573,453208.190,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,EXPONENTIATION,50000,0.111903123,446815.054,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,EXPONENTIATION,50000,0.103085523,485034.159,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,EXPONENTIATION,50000,0.103432461,483407.234,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,p1024,1024,DIVIDE,1562,0.000049217,31737178.042,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,p1024,1024,DIVIDE,1562,0.000002781,561748411.211,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,p1024,1024,DIVIDE,1562,0.000016019,97506931.199,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,p1024,1024,DIVIDE,50000,0.000032832,1522904483.431,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,DIVIDE,50000,0.012967985,3855649.156,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,DIVIDE,50000,0.014316212,3492543.940,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,DIVIDE,50000,0.007058040,7084119.665,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,DIVIDE,50000,0.008473266,5900912.410,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,DIVIDE,50000,0.001525701,32771820.445,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,DIVIDE,50000,0.003247566,15396145.650,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,DIVIDE,50000,0.001533563,32603811.828,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,DIVIDE,50000,0.003189910,15674423.562,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,DIVIDE,50000,0.001487084,33622849.245,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,DIVIDE,50000,0.003095589,16152015.640,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,DIVIDE,50000,0.001539612,32475715.586,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,DIVIDE,50000,0.003259545,15339564.456,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,p1024,1024,ISQRT,390,0.000076154,5121190.389,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,p1024,1024,ISQRT,390,0.000004038,96575785.219,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,ISQRT,50000,0.198586034,251780.042,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,ISQRT,50000,0.199198550,251005.843,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,ISQRT,50000,0.122290257,408863.316,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,ISQRT,50000,0.123460206,404988.795,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,ISQRT,50000,0.022803578,2192638.360,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,ISQRT,50000,0.024450240,2044969.686,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,ISQRT,50000,0.022468740,2225313.936,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,ISQRT,50000,0.024039980,2079868.617,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,ISQRT,50000,0.022423225,2229830.902,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,ISQRT,50000,0.023950578,2087632.295,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,ISQRT,50000,0.021812202,2292294.933,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,ISQRT,50000,0.023405364,2136262.438,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,p1024,1024,MODMUL_R2,12500,0.006167274,2026827.422,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,p1024,1024,MODMUL_R2,12500,0.000276313,45238563.309,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,p1024,1024,MODMUL_R2,12500,0.000954096,13101404.164,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,p1024,1024,MODMUL_R2,50000,0.000060608,824973600.845,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,MODMUL_R2,50000,0.002749753,18183451.705,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p1024,1024,MODMUL_R2,50000,0.003914755,12772191.893,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,MODMUL_R2,50000,0.000252147,198297046.898,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p1024,1024,MODMUL_R2,50000,0.001335626,37435629.368,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,MODMUL_R2,50000,0.000202723,246641871.897,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p1024,1024,MODMUL_R2,50000,0.001591960,31407824.944,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,MODMUL_R2,50000,0.000137255,364285291.074,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p1024,1024,MODMUL_R2,50000,0.001169699,42746039.378,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,MODMUL_R2,50000,0.000166699,299941707.928,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p1024,1024,MODMUL_R2,50000,0.001529828,32683413.912,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,MODMUL_R2,50000,0.000100791,496076094.379,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p1024,1024,MODMUL_R2,50000,0.001519152,32913099.682,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,p2048,2048,ADD,6250,0.000141910,44042000.704,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,p2048,2048,ADD,6250,0.000446799,13988382.329,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,p2048,2048,ADD,6250,0.000375567,16641490.775,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,p2048,2048,ADD,50000,0.000023296,2146291208.791,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,ADD,50000,0.000383714,130305390.103,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,ADD,50000,0.002468622,20254216.144,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,ADD,50000,0.000188111,265800706.497,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,ADD,50000,0.002316623,21583140.680,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,p2048,2048,ADD,50000,0.000078708,635259977.281,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,p2048,2048,ADD,50000,0.002257816,22145296.401,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,ADD,50000,0.000080381,622037124.965,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,ADD,50000,0.002778635,17994447.596,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,ADD,50000,0.000085748,583103750.628,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,ADD,50000,0.002927458,17079663.482,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,ADD,50000,0.000069575,718647375.871,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,ADD,50000,0.003109921,16077578.805,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,ADD,50000,0.000074472,671394963.952,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,ADD,50000,0.002778134,17997692.992,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,p2048,2048,SUBTRACT,6250,0.000138961,44976755.632,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,p2048,2048,SUBTRACT,6250,0.000166532,37530335.418,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,p2048,2048,SUBTRACT,6250,0.000008248,757760165.179,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,p2048,2048,SUBTRACT,50000,0.000023744,2105795148.248,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,SUBTRACT,50000,0.000381802,130957939.344,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,SUBTRACT,50000,0.002441341,20480548.843,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,SUBTRACT,50000,0.000187781,266267373.908,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,SUBTRACT,50000,0.002271976,22007273.862,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,p2048,2048,SUBTRACT,50000,0.000078497,636967108.221,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,p2048,2048,SUBTRACT,50000,0.002263324,22091400.850,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,SUBTRACT,50000,0.000081261,615301907.952,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,SUBTRACT,50000,0.002823392,17709195.624,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,SUBTRACT,50000,0.000085548,584466980.562,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,SUBTRACT,50000,0.002928680,17072537.555,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,SUBTRACT,50000,0.000068493,730002089.912,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,SUBTRACT,50000,0.003102260,16117282.534,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,SUBTRACT,50000,0.000074651,669782564.000,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,SUBTRACT,50000,0.002786578,17943155.217,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,p2048,2048,ADDMOD,6250,0.000422818,14781762.991,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,p2048,2048,ADDMOD,6250,0.000019762,316270033.484,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,p2048,2048,ADDMOD,6250,0.000375941,16624935.745,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,p2048,2048,ADDMOD,50000,0.000023328,2143347050.754,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,ADDMOD,50000,0.000481380,103868063.967,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,ADDMOD,50000,0.002661420,18786963.122,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,ADDMOD,50000,0.000243924,204981801.868,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,ADDMOD,50000,0.002362892,21160509.803,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,p2048,2048,ADDMOD,50000,0.000096594,517630583.222,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,p2048,2048,ADDMOD,50000,0.002285367,21878323.078,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,ADDMOD,50000,0.000112648,443860040.511,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,ADDMOD,50000,0.002963151,16873928.658,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,ADDMOD,50000,0.000112398,444847528.100,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,ADDMOD,50000,0.002949511,16951962.926,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,ADDMOD,50000,0.000038788,1289052215.566,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,ADDMOD,50000,0.003109961,16077371.775,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,ADDMOD,50000,0.000041562,1203018154.929,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,ADDMOD,50000,0.002768109,18062871.957,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,p2048,2048,SUBTRACTMOD,6250,0.000361477,17290173.196,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,p2048,2048,SUBTRACTMOD,6250,0.000025706,243136445.723,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,p2048,2048,SUBTRACTMOD,6250,0.000152285,41041445.718,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,p2048,2048,SUBTRACTMOD,50000,0.000023840,2097315436.242,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,SUBTRACTMOD,50000,0.000507098,98600240.040,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,SUBTRACTMOD,50000,0.002681350,18647322.463,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,SUBTRACTMOD,50000,0.000257124,194458734.863,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,SUBTRACTMOD,50000,0.002333168,21430090.451,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,p2048,2048,SUBTRACTMOD,50000,0.000109734,455647235.754,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,p2048,2048,SUBTRACTMOD,50000,0.002273440,21993101.723,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,SUBTRACTMOD,50000,0.000110085,454193981.540,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,SUBTRACTMOD,50000,0.002989640,16724421.704,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,SUBTRACTMOD,50000,0.000112729,443541928.302,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,SUBTRACTMOD,50000,0.002958535,16900256.469,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,SUBTRACTMOD,50000,0.000038698,1292061446.639,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,SUBTRACTMOD,50000,0.003118053,16035648.680,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,SUBTRACTMOD,50000,0.000040050,1248442645.614,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,SUBTRACTMOD,50000,0.002764205,18088383.573,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.003430290,1822003.124,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.000258942,24136677.554,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.000272227,22958781.943,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.092110865,542824.128,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.094645356,528287.938,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.023158429,2159041.081,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.025991945,1923672.904,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.005191533,9631066.348,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.007935125,6301098.034,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.000431004,116008209.392,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.002814849,17762942.754,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.000420288,118965990.591,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.003792271,13184711.684,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.000341951,146219816.051,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.003887133,12862950.748,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.000304175,164379017.500,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.003513204,14232023.231,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.003431532,1821343.963,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.000247134,25289874.787,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.000260825,23962423.599,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.000064384,776590457.256,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.033502928,1492406.873,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.036325060,1376460.224,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.008511464,5874430.147,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.010890310,4591237.560,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.002162243,23124135.374,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.005225383,9568676.206,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.002157236,23177807.331,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.004620339,10821716.708,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.002159719,23151158.629,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.005495267,9098738.988,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.000933966,53535145.074,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.004476345,11169827.390,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.000872515,57305611.029,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.004080303,12253991.882,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.023393546,267167.705,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.001459133,4283365.202,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.003445704,1813852.843,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000098144,509455493.968,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.043045861,1161551.857,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.045342166,1102726.324,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.001331511,37551324.580,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.003167705,15784298.535,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000227159,220110250.502,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.003064622,16315224.914,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000402002,124377450.530,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.002588681,19314856.120,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000292207,171111644.366,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.003139745,15924860.249,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000344274,145233213.381,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.003465162,14429339.293,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000266268,187780789.988,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.002982610,16763840.289,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,p2048,2048,COMPARE,6250,0.000051202,122065463.765,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,p2048,2048,COMPARE,6250,0.000003149,1984988851.542,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,p2048,2048,COMPARE,6250,0.000003138,1991504210.356,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,p2048,2048,COMPARE,50000,0.000023520,2125850340.136,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,COMPARE,50000,0.000268832,185989706.441,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,COMPARE,50000,0.002324807,21507161.744,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,COMPARE,50000,0.000132408,377620707.312,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,COMPARE,50000,0.002401570,20819715.000,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,COMPARE,50000,0.000074792,668521101.146,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,COMPARE,50000,0.002926836,17083293.907,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,COMPARE,50000,0.000076615,652613558.055,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,COMPARE,50000,0.002915180,17151600.602,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,COMPARE,50000,0.000022504,2221825944.089,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,COMPARE,50000,0.003099466,16131811.220,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,COMPARE,50000,0.000024687,2025354756.201,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,COMPARE,50000,0.002758646,18124832.712,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,p2048,2048,REDUCE,781,0.000022880,34134396.637,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,p2048,2048,REDUCE,781,0.000001994,391584338.744,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,p2048,2048,REDUCE,781,0.000017243,45292774.484,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,p2048,2048,REDUCE,50000,0.000028864,1732261640.798,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,REDUCE,50000,0.137730537,363027.700,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,REDUCE,50000,0.140036646,357049.397,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,REDUCE,50000,0.004678093,10688115.566,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,REDUCE,50000,0.006468490,7729779.285,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,REDUCE,50000,0.001007877,49609225.284,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,REDUCE,50000,0.003839801,13021508.713,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,REDUCE,50000,0.000956350,52282113.952,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,REDUCE,50000,0.003813823,13110203.888,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,REDUCE,50000,0.000893566,55955572.058,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,REDUCE,50000,0.003985100,12546736.719,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,REDUCE,50000,0.000820147,60964672.811,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,REDUCE,50000,0.003547325,14095127.893,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,p2048,2048,MODMUL,390,0.000613711,635478.563,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,p2048,2048,MODMUL,390,0.000045326,8604256.108,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,p2048,2048,MODMUL,390,0.000098559,3957025.705,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,p2048,2048,MODMUL,50000,0.000432320,115655070.318,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,MODMUL,50000,0.308650013,161995.781,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,MODMUL,50000,0.311279936,160627.121,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,MODMUL,50000,0.014862818,3364099.606,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,MODMUL,50000,0.016676089,2998304.933,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,MODMUL,50000,0.005552352,9005192.851,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,MODMUL,50000,0.007918460,6314359.242,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,MODMUL,50000,0.003789537,13194224.382,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,MODMUL,50000,0.006658277,7509450.258,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,MODMUL,50000,0.005491913,9104295.611,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,MODMUL,50000,0.008601454,5812970.790,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,MODMUL,50000,0.003616288,13826332.143,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,MODMUL,50000,0.006364238,7856400.112,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,p2048,2048,MODEXP,97,0.223744308,433.531,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,p2048,2048,MODEXP,97,0.013428724,7223.322,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,p2048,2048,MODEXP,97,0.007201191,13469.994,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,p2048,2048,MODEXP,50000,0.630893409,79252.690,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,MODEXP,50000,41.177744162,1214.248,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,MODEXP,50000,41.170890865,1214.450,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,MODEXP,50000,5.982879414,8357.180,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,MODEXP,50000,5.869589012,8518.484,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,MODEXP,50000,1.227683828,40727.098,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,MODEXP,50000,1.239477790,40339.569,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,MODEXP,50000,1.608675414,31081.472,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,MODEXP,50000,1.607744418,31099.470,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,MODEXP,50000,1.265466806,39511.111,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,MODEXP,50000,1.266663235,39473.791,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,MODEXP,50000,1.701812338,29380.443,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,MODEXP,50000,1.654067649,30228.510,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,p2048,2048,EXPONENTIATION,97,0.026668645,3637.230,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,p2048,2048,EXPONENTIATION,97,0.001946020,49845.327,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,p2048,2048,EXPONENTIATION,97,0.020040482,4840.203,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,EXPONENTIATION,50000,13.769309211,3631.264,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,EXPONENTIATION,50000,13.764768149,3632.462,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,EXPONENTIATION,50000,3.168514098,15780.267,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,EXPONENTIATION,50000,3.165133402,15797.122,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,EXPONENTIATION,50000,0.943866312,52973.604,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,EXPONENTIATION,50000,0.945367729,52889.472,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,EXPONENTIATION,50000,0.816244479,61256.157,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,EXPONENTIATION,50000,0.819034400,61047.497,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,EXPONENTIATION,50000,0.922355657,54209.024,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,EXPONENTIATION,50000,0.926634566,53958.704,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,EXPONENTIATION,50000,0.813492040,61463.416,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,EXPONENTIATION,50000,0.815551340,61308.219,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,p2048,2048,DIVIDE,781,0.000032852,23773517.467,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,p2048,2048,DIVIDE,781,0.000002513,310838435.394,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,p2048,2048,DIVIDE,781,0.000012891,60584214.912,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,p2048,2048,DIVIDE,50000,0.000033888,1475448536.355,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,DIVIDE,50000,0.677662082,73783.086,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,DIVIDE,50000,0.680329410,73493.809,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,DIVIDE,50000,0.147293724,339457.776,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,DIVIDE,50000,0.149842175,333684.425,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,DIVIDE,50000,0.011718254,4266847.237,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,DIVIDE,50000,0.014746352,3390669.073,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,DIVIDE,50000,0.011225608,4454101.751,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,DIVIDE,50000,0.013966777,3579923.954,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,DIVIDE,50000,0.011305939,4422454.452,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,DIVIDE,50000,0.015096328,3312063.704,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,DIVIDE,50000,0.010994033,4547921.508,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,DIVIDE,50000,0.013633270,3667498.749,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,p2048,2048,ISQRT,195,0.000065035,2998396.150,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,p2048,2048,ISQRT,195,0.000004041,48257553.939,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,ISQRT,50000,9.709492530,5149.600,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,ISQRT,50000,9.719505388,5144.295,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,ISQRT,50000,4.225232380,11833.669,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,ISQRT,50000,4.213225696,11867.392,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,ISQRT,50000,0.061130766,817918.757,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,ISQRT,50000,0.063863873,782915.248,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,ISQRT,50000,0.052201978,957818.113,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,ISQRT,50000,0.054501618,917403.956,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,ISQRT,50000,0.057662489,867114.843,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,ISQRT,50000,0.060918428,820769.702,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,ISQRT,50000,0.052043607,960732.795,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,ISQRT,50000,0.054484197,917697.292,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-1t,p2048,2048,MODMUL_R2,6250,0.010111392,618114.696,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,gmp-nt,p2048,2048,MODMUL_R2,6250,0.000612115,10210499.297,0
library,AMD Eng Sample: 100-000000897-03,host-cpu,openssl-nt,p2048,2048,MODMUL_R2,6250,0.001366114,4575020.789,0
library,NVIDIA H100 80GB HBM3,gpu,cgbn,p2048,2048,MODMUL_R2,50000,0.000185312,269815230.530,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,MODMUL_R2,50000,0.023851393,2096313.618,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w8,p2048,2048,MODMUL_R2,50000,0.026195196,1908746.934,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,MODMUL_R2,50000,0.001225673,40793918.345,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w16,p2048,2048,MODMUL_R2,50000,0.003301036,15146760.695,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,MODMUL_R2,50000,0.000704944,70927629.751,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-opt,p2048,2048,MODMUL_R2,50000,0.003448176,14500420.080,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,MODMUL_R2,50000,0.000490483,101940353.782,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-o64,p2048,2048,MODMUL_R2,50000,0.002669021,18733460.544,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,MODMUL_R2,50000,0.000653387,76524320.772,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il,p2048,2048,MODMUL_R2,50000,0.003905801,12801471.879,0
opencl-kernel,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,MODMUL_R2,50000,0.000446999,111857078.537,0
opencl-e2e,NVIDIA H100 80GB HBM3,GPU,w32-il64,p2048,2048,MODMUL_R2,50000,0.002501982,19984156.310,0
```
