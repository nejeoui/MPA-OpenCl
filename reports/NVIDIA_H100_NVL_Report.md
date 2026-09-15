# MPA-OpenCL benchmark report - NVIDIA H100 NVL

## 1. System under test

1 OpenCL device(s) exercised with the identical kernels and operands.

### Device 0 - NVIDIA H100 NVL (GPU)

| Property | Value |
|---|---|
| Model | NVIDIA H100 NVL |
| Type | GPU |
| Vendor | NVIDIA Corporation |
| Device memory | 93.09 GiB |
| Max single allocation | 23.27 GiB |
| Local memory | 48 KiB |
| Global cache | 4224 KiB |
| Compute units | 132 |
| Max clock | 1785 MHz |
| Max work-group size | 1024 |
| OpenCL version | OpenCL 3.0 CUDA |
| Driver | 595.71.05 |

### Host

| Property | Value |
|---|---|
| CPU | AMD EPYC 9V84 96-Core Processor |
| Logical cores | 40 |
| OpenMP threads used | 40 |
| RAM | 314.7 GB |
| OS | Ubuntu 24.04.4 LTS |
| Kernel | 6.17.0-1022-azure |
| Arch | x86_64 |
| GMP | 6.3.0 |
| OpenSSL | OpenSSL 3.0.13 30 Jan 2024 |
| CGBN | 60 rows from `cgbn_results_NVIDIA_H100_NVL.tsv` |

## 2. Method

- Base workload 500000 items, scaled down per operator by its cost weight and by modulus size. Device rows honour --min-items (500000) so the GPU is not left idle; the CPU libraries keep the smaller count because a full-width MODEXP there costs minutes. Both counts appear in every row as dev/cpu, and throughput is per-second so they remain comparable.
- 5 timed repetitions, **minimum** reported. Two untimed warm-up launches precede them.
- `kernel` times `clEnqueueNDRangeKernel` + `clFinish` only. `e2e` adds the host->device operand writes and the device->host result read.
- Every OpenCL device runs the same kernels on the same operands, so GPU and CPU-OpenCL columns are directly comparable.
- CPU library baselines (GMP, OpenSSL) run those same operands, with every temporary - including each thread's GMP context, BN_CTX and Montgomery context - allocated outside the timed region, so the figure is the arithmetic and not marshalling. The generator is reseeded per modulus and operation so every backend sees identical inputs.
- Cost weighting drives the wide cells down to a few hundred items, which is tens of microseconds of work - the same order as the cost of entering an OpenMP region. Each baseline pass is therefore repeated until the timed interval reaches 5 ms and the per-pass time is reported; the multi-threaded loop enters one parallel region per interval and partitions the range itself. Without this the multi-threaded GMP figure came out up to 9x slower than the single-threaded one at 2048 bits.
- OpenSSL rows time the nearest BN primitive, which is not always semantically identical (its Montgomery routine expects Montgomery-domain inputs); they measure comparable work, not identical results. Correctness is judged against GMP only.
- Every device cell is checked word-for-word against GMP before it is timed. A cell that mismatches is reported and excluded from the speedup tables.
- Total wall time 7784.8 s.

## 3. Correctness

| Device | Kernel | Configs run | Passed | Mismatched | Launch failed |
|---|---|---|---|---|---|
| [0] GPU | `mpaKernel_32bits_opt.cl` (w32-opt) | 75 | 75 | 0 | 0 |
| [0] GPU | `mpaKernel_32bits_opt.cl` (w32-o64) | 75 | 75 | 0 | 0 |
| [0] GPU | `mpaKernel_32bits_opt.cl` (w32-il) | 75 | 75 | 0 | 0 |
| [0] GPU | `mpaKernel_32bits_opt.cl` (w32-il64) | 75 | 75 | 0 | 0 |

**All configurations correct** - 300 configurations, 0 problems.

## 4. Throughput per device

Operations per second, higher is better. Kernel-only timings.

### Device 0 - NVIDIA H100 NVL (GPU)

#### secp256k1 (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | GMP 40T | OpenSSL 40T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 500000 / 500000 | - | - | - | 9.58 G | 9.78 G | 13.70 G | 13.92 G | 68.07 M | 2.79 G | 2.64 G | 12.29 G |
| SUBTRACT | 500000 / 500000 | - | - | - | 9.82 G | 9.69 G | 13.59 G | 13.87 G | 88.49 M | 3.53 G | 1.85 G | 12.18 G |
| ADDMOD | 500000 / 500000 | - | - | - | 11.18 G | 11.17 G | 16.19 G | 16.23 G | 25.41 M | 983.83 M | 157.72 M | 9.84 G |
| SUBTRACTMOD | 500000 / 500000 | - | - | - | 10.93 G | 11.06 G | 16.33 G | 16.36 G | 31.68 M | 1.19 G | 94.62 M | 9.08 G |
| MULTIPLYOPERANDSCANNING | 500000 / 500000 | - | - | - | 4.82 G | 4.71 G | 13.63 G | 16.32 G | 55.34 M | 2.12 G | 488.47 M | n/a |
| MULTIPLYPRODUCTSCANNING | 500000 / 500000 | - | - | - | 3.97 G | 4.33 G | 5.24 G | 5.43 G | 55.42 M | 2.12 G | 497.35 M | 12.33 G |
| MONTGOMERYMULTIPLICATION | 500000 / 500000 | - | - | - | 9.45 G | 12.95 G | 9.30 G | 13.95 G | 6.98 M | 253.02 M | 1.01 G | 8.08 G |
| COMPARE | 500000 / 500000 | - | - | - | 17.16 G | 15.71 G | 35.54 G | 36.71 G | 129.17 M | 2.28 G | 2.16 G | 12.12 G |
| REDUCE | 500000 / 62500 | - | - | - | 7.00 G | 12.66 G | 6.87 G | 12.83 G | 65.74 M | 997.82 M | 111.16 M | 5.80 G |
| MODMUL | 500000 / 31250 | - | - | - | 670.92 M | 878.57 M | 667.17 M | 857.48 M | 12.14 M | 215.75 M | 92.47 M | 1.79 G |
| MODEXP | 500000 / 7812 | - | - | - | 25.20 M | 57.04 M | 25.21 M | 56.04 M | 113.45 k | 1.42 M | 936.10 k | 5.29 M |
| EXPONENTIATION | 500000 / 7812 | - | - | - | 128.02 M | 246.11 M | 129.14 M | 245.28 M | 358.81 k | 13.20 M | 1.08 M | n/a |
| DIVIDE | 500000 / 62500 | - | - | - | 1.50 G | 1.61 G | 1.53 G | 1.70 G | 35.62 M | 559.99 M | 147.01 M | 4.44 G |
| ISQRT | 500000 / 15625 | - | - | - | 97.52 M | 135.71 M | 97.03 M | 129.78 M | 15.40 M | 330.71 M | n/a | n/a |
| MODMUL_R2 | 500000 / 500000 | - | - | - | 5.77 G | 9.53 G | 5.66 G | 9.58 G | 12.12 M | 477.23 M | 89.38 M | 4.89 G |

#### rsa256(composite) (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | GMP 40T | OpenSSL 40T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 500000 / 500000 | - | - | - | 9.70 G | 9.77 G | 13.80 G | 13.78 G | 63.51 M | 1.66 G | 2.61 G | 12.19 G |
| SUBTRACT | 500000 / 500000 | - | - | - | 9.68 G | 9.79 G | 13.68 G | 13.75 G | 82.29 M | 3.51 G | 1.88 G | 12.25 G |
| ADDMOD | 500000 / 500000 | - | - | - | 11.40 G | 11.15 G | 16.43 G | 16.20 G | 29.72 M | 1.18 G | 173.64 M | 9.21 G |
| SUBTRACTMOD | 500000 / 500000 | - | - | - | 11.32 G | 11.11 G | 16.16 G | 16.32 G | 31.54 M | 1.24 G | 113.20 M | 9.03 G |
| MULTIPLYOPERANDSCANNING | 500000 / 500000 | - | - | - | 4.79 G | 4.64 G | 13.99 G | 16.17 G | 55.32 M | 2.12 G | 1.34 G | n/a |
| MULTIPLYPRODUCTSCANNING | 500000 / 500000 | - | - | - | 4.02 G | 4.34 G | 5.26 G | 5.41 G | 55.46 M | 2.12 G | 494.74 M | 12.40 G |
| MONTGOMERYMULTIPLICATION | 500000 / 500000 | - | - | - | 9.44 G | 13.20 G | 9.39 G | 13.69 G | 6.97 M | 252.64 M | 1.02 G | 8.13 G |
| COMPARE | 500000 / 500000 | - | - | - | 18.20 G | 18.85 G | 35.82 G | 35.09 G | 123.50 M | 2.84 G | 4.14 G | 12.17 G |
| REDUCE | 500000 / 62500 | - | - | - | 6.86 G | 12.69 G | 6.14 G | 12.54 G | 40.90 M | 697.55 M | 122.30 M | 5.80 G |
| MODMUL | 500000 / 31250 | - | - | - | 670.60 M | 878.40 M | 668.29 M | 858.66 M | 12.16 M | 215.30 M | 90.19 M | 1.79 G |
| MODEXP | 500000 / 7812 | - | - | - | 25.21 M | 57.09 M | 25.21 M | 56.06 M | 120.94 k | 4.78 M | 2.44 M | 5.38 M |
| EXPONENTIATION | 500000 / 7812 | - | - | - | 127.99 M | 246.29 M | 129.05 M | 245.26 M | 362.96 k | 13.16 M | 1.11 M | n/a |
| DIVIDE | 500000 / 62500 | - | - | - | 1.49 G | 1.61 G | 1.53 G | 1.68 G | 35.46 M | 38.86 M | 134.34 M | 4.43 G |
| ISQRT | 500000 / 15625 | - | - | - | 97.47 M | 135.71 M | 96.98 M | 129.76 M | 15.36 M | 330.53 M | n/a | n/a |
| MODMUL_R2 | 500000 / 500000 | - | - | - | 5.68 G | 9.33 G | 5.64 G | 9.63 G | 12.12 M | 476.53 M | 56.56 M | 4.91 G |

#### brainpoolP512r1 (512-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | GMP 40T | OpenSSL 40T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 500000 / 250000 | - | - | - | 4.54 G | 4.55 G | 6.09 G | 5.95 G | 58.55 M | 1.52 G | 925.90 M | 9.19 G |
| SUBTRACT | 500000 / 250000 | - | - | - | 4.53 G | 4.59 G | 6.14 G | 6.01 G | 74.89 M | 1.16 G | 1.30 G | 9.30 G |
| ADDMOD | 500000 / 250000 | - | - | - | 3.97 G | 3.79 G | 9.55 G | 9.34 G | 26.79 M | 1.04 G | 184.00 M | 7.37 G |
| SUBTRACTMOD | 500000 / 250000 | - | - | - | 3.96 G | 3.96 G | 9.58 G | 9.28 G | 28.72 M | 1.10 G | 116.64 M | 7.30 G |
| MULTIPLYOPERANDSCANNING | 500000 / 250000 | - | - | - | 1.13 G | 1.13 G | 2.48 G | 2.35 G | 22.88 M | 885.09 M | 816.89 M | n/a |
| MULTIPLYPRODUCTSCANNING | 500000 / 250000 | - | - | - | 683.57 M | 683.61 M | 1.45 G | 1.54 G | 22.90 M | 895.07 M | 823.03 M | 7.58 G |
| MONTGOMERYMULTIPLICATION | 500000 / 250000 | - | - | - | 2.36 G | 4.50 G | 2.60 G | 5.89 G | 2.94 M | 117.08 M | 70.02 M | 5.25 G |
| COMPARE | 500000 / 250000 | - | - | - | 4.97 G | 4.97 G | 24.89 G | 25.11 G | 95.53 M | 68.00 M | 74.97 M | 9.30 G |
| REDUCE | 500000 / 31250 | - | - | - | 1.66 G | 3.94 G | 1.71 G | 4.12 G | 37.37 M | 41.24 M | 39.75 M | 3.34 G |
| MODMUL | 500000 / 15625 | - | - | - | 185.54 M | 237.94 M | 185.43 M | 238.53 M | 6.13 M | 229.93 M | 57.05 M | 566.00 M |
| MODEXP | 500000 / 3906 | - | - | - | 2.95 M | 7.16 M | 2.96 M | 7.13 M | 19.66 k | 488.66 k | 843.26 k | 2.20 M |
| EXPONENTIATION | 500000 / 3906 | - | - | - | 4.25 M | 4.67 M | 4.42 M | 4.58 M | 114.78 k | 3.04 M | 306.26 k | n/a |
| DIVIDE | 500000 / 31250 | - | - | - | 443.91 M | 416.12 M | 505.87 M | 470.42 M | 33.28 M | 988.30 M | 139.20 M | 2.34 G |
| ISQRT | 500000 / 7812 | - | - | - | 23.94 M | 26.10 M | 24.51 M | 28.30 M | 8.62 M | 193.84 M | n/a | n/a |
| MODMUL_R2 | 500000 / 250000 | - | - | - | 1.34 G | 2.90 G | 1.39 G | 3.33 G | 6.14 M | 54.01 M | 33.53 M | 2.93 G |

#### p1024 (1024-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | GMP 40T | OpenSSL 40T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 500000 / 125000 | - | - | - | 1.67 G | 1.68 G | 1.84 G | 1.85 G | 52.29 M | 71.80 M | 72.34 M | 5.23 G |
| SUBTRACT | 500000 / 125000 | - | - | - | 1.67 G | 1.67 G | 1.83 G | 1.85 G | 56.57 M | 878.03 M | 70.72 M | 5.21 G |
| ADDMOD | 500000 / 125000 | - | - | - | 1.34 G | 1.34 G | 4.53 G | 5.09 G | 19.14 M | 34.48 M | 34.77 M | 4.89 G |
| SUBTRACTMOD | 500000 / 125000 | - | - | - | 1.33 G | 1.33 G | 4.65 G | 4.75 G | 23.67 M | 144.55 M | 34.11 M | 4.75 G |
| MULTIPLYOPERANDSCANNING | 500000 / 125000 | - | - | - | 436.84 M | 438.42 M | 686.68 M | 743.91 M | 6.29 M | 247.24 M | 130.03 M | n/a |
| MULTIPLYPRODUCTSCANNING | 500000 / 125000 | - | - | - | 101.97 M | 101.98 M | 271.65 M | 295.90 M | 6.29 M | 247.32 M | 131.52 M | 2.52 G |
| MONTGOMERYMULTIPLICATION | 500000 / 125000 | - | - | - | 499.68 M | 737.62 M | 691.11 M | 1.20 G | 927.73 k | 36.18 M | 135.83 M | 1.62 G |
| COMPARE | 500000 / 125000 | - | - | - | 1.81 G | 1.82 G | 12.59 G | 12.81 G | 122.35 M | 2.40 G | 1.90 G | 5.19 G |
| REDUCE | 500000 / 15625 | - | - | - | 380.46 M | 590.75 M | 439.35 M | 757.86 M | 52.67 M | 830.46 M | 44.17 M | 2.18 G |
| MODMUL | 500000 / 7812 | - | - | - | 40.20 M | 51.77 M | 40.25 M | 52.22 M | 2.14 M | 50.27 M | 20.90 M | 230.40 M |
| MODEXP | 500000 / 1953 | - | - | - | 366.27 k | 676.13 k | 364.14 k | 678.49 k | 3.12 k | 123.21 k | 95.01 k | 413.64 k |
| EXPONENTIATION | 500000 / 1953 | - | - | - | 504.47 k | 577.59 k | 498.46 k | 577.67 k | 23.61 k | 702.78 k | 40.60 k | n/a |
| DIVIDE | 500000 / 15625 | - | - | - | 126.47 M | 135.57 M | 137.92 M | 152.34 M | 29.53 M | 1.17 G | 64.12 M | 1.64 G |
| ISQRT | 500000 / 3906 | - | - | - | 3.43 M | 4.28 M | 3.72 M | 4.27 M | 4.49 M | 87.07 M | n/a | n/a |
| MODMUL_R2 | 500000 / 125000 | - | - | - | 309.20 M | 465.65 M | 354.16 M | 605.25 M | 2.13 M | 84.14 M | 26.07 M | 839.33 M |

#### p2048 (2048-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | GMP 40T | OpenSSL 40T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 500000 / 62500 | - | - | - | 855.67 M | 818.55 M | 935.78 M | 910.11 M | 36.24 M | 715.89 M | 520.23 M | 2.73 G |
| SUBTRACT | 500000 / 62500 | - | - | - | 855.04 M | 817.71 M | 937.92 M | 896.82 M | 34.77 M | 50.87 M | 1.14 G | 2.74 G |
| ADDMOD | 500000 / 62500 | - | - | - | 669.55 M | 655.93 M | 2.09 G | 1.87 G | 15.45 M | 564.44 M | 85.16 M | 2.75 G |
| SUBTRACTMOD | 500000 / 62500 | - | - | - | 678.55 M | 656.53 M | 1.86 G | 2.04 G | 17.83 M | 526.25 M | 70.30 M | 2.73 G |
| MULTIPLYOPERANDSCANNING | 500000 / 62500 | - | - | - | 149.41 M | 148.45 M | 177.42 M | 209.59 M | 1.92 M | 75.19 M | 58.45 M | n/a |
| MULTIPLYPRODUCTSCANNING | 500000 / 62500 | - | - | - | 26.26 M | 26.32 M | 67.23 M | 70.59 M | 1.92 M | 75.32 M | 58.58 M | 769.21 M |
| MONTGOMERYMULTIPLICATION | 500000 / 62500 | - | - | - | 154.36 M | 213.85 M | 174.76 M | 246.59 M | 282.21 k | 11.20 M | 38.37 M | 490.27 M |
| COMPARE | 500000 / 62500 | - | - | - | 962.30 M | 916.15 M | 6.96 G | 5.80 G | 129.29 M | 2.34 G | 2.10 G | 2.74 G |
| REDUCE | 500000 / 7812 | - | - | - | 110.54 M | 148.25 M | 113.82 M | 160.63 M | 35.93 M | 573.63 M | 77.53 M | 1.83 G |
| MODMUL | 500000 / 3906 | - | - | - | 10.54 M | 15.78 M | 11.36 M | 15.86 M | 653.27 k | 13.60 M | 11.55 M | 107.17 M |
| MODEXP | 500000 / 976 | - | - | - | 46.53 k | 69.05 k | 45.00 k | 69.44 k | 413.9 | 15.95 k | 18.92 k | 71.87 k |
| EXPONENTIATION | 500000 / 976 | - | - | - | 60.77 k | 70.66 k | 59.99 k | 70.71 k | 3.90 k | 108.48 k | 9.12 k | n/a |
| DIVIDE | 500000 / 7812 | - | - | - | 7.30 M | 7.41 M | 7.08 M | 7.70 M | 23.08 M | 24.10 M | 25.40 M | 1.60 G |
| ISQRT | 500000 / 1953 | - | - | - | 362.89 k | 1.07 M | 363.72 k | 1.09 M | 2.73 M | 3.81 M | n/a | n/a |
| MODMUL_R2 | 500000 / 62500 | - | - | - | 89.45 M | 120.36 M | 95.48 M | 133.51 M | 650.77 k | 25.79 M | 11.70 M | 252.97 M |

## 5. Head to head

Best OpenCL GPU result against best OpenCL CPU result and the CPU libraries.
Ratios above 1.00x mean the GPU is faster than that baseline.

### secp256k1 (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GMP 40T | OpenSSL | CGBN | GPU vs CPU-CL | GPU vs GMP 40T | GPU vs OpenSSL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | w32-il64 | 13.92 G | none | n/a | 68.07 M | 2.79 G | 2.64 G | 12.29 G | n/a | 4.98x | 5.27x | 1.13x |
| SUBTRACT | w32-il64 | 13.87 G | none | n/a | 88.49 M | 3.53 G | 1.85 G | 12.18 G | n/a | 3.92x | 7.49x | 1.14x |
| ADDMOD | w32-il64 | 16.23 G | none | n/a | 25.41 M | 983.83 M | 157.72 M | 9.84 G | n/a | 16.50x | 102.92x | 1.65x |
| SUBTRACTMOD | w32-il64 | 16.36 G | none | n/a | 31.68 M | 1.19 G | 94.62 M | 9.08 G | n/a | 13.80x | 172.86x | 1.80x |
| MULTIPLYOPERANDSCANNING | w32-il64 | 16.32 G | none | n/a | 55.34 M | 2.12 G | 488.47 M | n/a | n/a | 7.71x | 33.42x | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 5.43 G | none | n/a | 55.42 M | 2.12 G | 497.35 M | 12.33 G | n/a | 2.56x | 10.91x | 0.44x |
| MONTGOMERYMULTIPLICATION | w32-il64 | 13.95 G | none | n/a | 6.98 M | 253.02 M | 1.01 G | 8.08 G | n/a | 55.12x | 13.78x | 1.73x |
| COMPARE | w32-il64 | 36.71 G | none | n/a | 129.17 M | 2.28 G | 2.16 G | 12.12 G | n/a | 16.11x | 16.96x | 3.03x |
| REDUCE | w32-il64 | 12.83 G | none | n/a | 65.74 M | 997.82 M | 111.16 M | 5.80 G | n/a | 12.86x | 115.42x | 2.21x |
| MODMUL | w32-o64 | 878.57 M | none | n/a | 12.14 M | 215.75 M | 92.47 M | 1.79 G | n/a | 4.07x | 9.50x | 0.49x |
| MODEXP | w32-o64 | 57.04 M | none | n/a | 113.45 k | 1.42 M | 936.10 k | 5.29 M | n/a | 40.21x | 60.94x | 10.78x |
| EXPONENTIATION | w32-o64 | 246.11 M | none | n/a | 358.81 k | 13.20 M | 1.08 M | n/a | n/a | 18.64x | 227.27x | n/a |
| DIVIDE | w32-il64 | 1.70 G | none | n/a | 35.62 M | 559.99 M | 147.01 M | 4.44 G | n/a | 3.04x | 11.57x | 0.38x |
| ISQRT | w32-o64 | 135.71 M | none | n/a | 15.40 M | 330.71 M | n/a | n/a | n/a | 0.41x | n/a | n/a |
| MODMUL_R2 | w32-il64 | 9.58 G | none | n/a | 12.12 M | 477.23 M | 89.38 M | 4.89 G | n/a | 20.07x | 107.17x | 1.96x |

### rsa256(composite) (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GMP 40T | OpenSSL | CGBN | GPU vs CPU-CL | GPU vs GMP 40T | GPU vs OpenSSL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | w32-il | 13.80 G | none | n/a | 63.51 M | 1.66 G | 2.61 G | 12.19 G | n/a | 8.29x | 5.29x | 1.13x |
| SUBTRACT | w32-il64 | 13.75 G | none | n/a | 82.29 M | 3.51 G | 1.88 G | 12.25 G | n/a | 3.91x | 7.31x | 1.12x |
| ADDMOD | w32-il | 16.43 G | none | n/a | 29.72 M | 1.18 G | 173.64 M | 9.21 G | n/a | 13.92x | 94.60x | 1.78x |
| SUBTRACTMOD | w32-il64 | 16.32 G | none | n/a | 31.54 M | 1.24 G | 113.20 M | 9.03 G | n/a | 13.14x | 144.21x | 1.81x |
| MULTIPLYOPERANDSCANNING | w32-il64 | 16.17 G | none | n/a | 55.32 M | 2.12 G | 1.34 G | n/a | n/a | 7.61x | 12.05x | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 5.41 G | none | n/a | 55.46 M | 2.12 G | 494.74 M | 12.40 G | n/a | 2.55x | 10.94x | 0.44x |
| MONTGOMERYMULTIPLICATION | w32-il64 | 13.69 G | none | n/a | 6.97 M | 252.64 M | 1.02 G | 8.13 G | n/a | 54.18x | 13.48x | 1.68x |
| COMPARE | w32-il | 35.82 G | none | n/a | 123.50 M | 2.84 G | 4.14 G | 12.17 G | n/a | 12.62x | 8.65x | 2.94x |
| REDUCE | w32-o64 | 12.69 G | none | n/a | 40.90 M | 697.55 M | 122.30 M | 5.80 G | n/a | 18.19x | 103.76x | 2.19x |
| MODMUL | w32-o64 | 878.40 M | none | n/a | 12.16 M | 215.30 M | 90.19 M | 1.79 G | n/a | 4.08x | 9.74x | 0.49x |
| MODEXP | w32-o64 | 57.09 M | none | n/a | 120.94 k | 4.78 M | 2.44 M | 5.38 M | n/a | 11.94x | 23.43x | 10.61x |
| EXPONENTIATION | w32-o64 | 246.29 M | none | n/a | 362.96 k | 13.16 M | 1.11 M | n/a | n/a | 18.71x | 222.01x | n/a |
| DIVIDE | w32-il64 | 1.68 G | none | n/a | 35.46 M | 38.86 M | 134.34 M | 4.43 G | n/a | 43.11x | 12.47x | 0.38x |
| ISQRT | w32-o64 | 135.71 M | none | n/a | 15.36 M | 330.53 M | n/a | n/a | n/a | 0.41x | n/a | n/a |
| MODMUL_R2 | w32-il64 | 9.63 G | none | n/a | 12.12 M | 476.53 M | 56.56 M | 4.91 G | n/a | 20.21x | 170.27x | 1.96x |

### brainpoolP512r1 (512-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GMP 40T | OpenSSL | CGBN | GPU vs CPU-CL | GPU vs GMP 40T | GPU vs OpenSSL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | w32-il | 6.09 G | none | n/a | 58.55 M | 1.52 G | 925.90 M | 9.19 G | n/a | 4.02x | 6.58x | 0.66x |
| SUBTRACT | w32-il | 6.14 G | none | n/a | 74.89 M | 1.16 G | 1.30 G | 9.30 G | n/a | 5.32x | 4.72x | 0.66x |
| ADDMOD | w32-il | 9.55 G | none | n/a | 26.79 M | 1.04 G | 184.00 M | 7.37 G | n/a | 9.22x | 51.90x | 1.30x |
| SUBTRACTMOD | w32-il | 9.58 G | none | n/a | 28.72 M | 1.10 G | 116.64 M | 7.30 G | n/a | 8.68x | 82.13x | 1.31x |
| MULTIPLYOPERANDSCANNING | w32-il | 2.48 G | none | n/a | 22.88 M | 885.09 M | 816.89 M | n/a | n/a | 2.80x | 3.03x | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 1.54 G | none | n/a | 22.90 M | 895.07 M | 823.03 M | 7.58 G | n/a | 1.72x | 1.88x | 0.20x |
| MONTGOMERYMULTIPLICATION | w32-il64 | 5.89 G | none | n/a | 2.94 M | 117.08 M | 70.02 M | 5.25 G | n/a | 50.30x | 84.10x | 1.12x |
| COMPARE | w32-il64 | 25.11 G | none | n/a | 95.53 M | 68.00 M | 74.97 M | 9.30 G | n/a | 369.34x | 334.97x | 2.70x |
| REDUCE | w32-il64 | 4.12 G | none | n/a | 37.37 M | 41.24 M | 39.75 M | 3.34 G | n/a | 99.99x | 103.73x | 1.24x |
| MODMUL | w32-il64 | 238.53 M | none | n/a | 6.13 M | 229.93 M | 57.05 M | 566.00 M | n/a | 1.04x | 4.18x | 0.42x |
| MODEXP | w32-o64 | 7.16 M | none | n/a | 19.66 k | 488.66 k | 843.26 k | 2.20 M | n/a | 14.64x | 8.49x | 3.25x |
| EXPONENTIATION | w32-o64 | 4.67 M | none | n/a | 114.78 k | 3.04 M | 306.26 k | n/a | n/a | 1.54x | 15.25x | n/a |
| DIVIDE | w32-il | 505.87 M | none | n/a | 33.28 M | 988.30 M | 139.20 M | 2.34 G | n/a | 0.51x | 3.63x | 0.22x |
| ISQRT | w32-il64 | 28.30 M | none | n/a | 8.62 M | 193.84 M | n/a | n/a | n/a | 0.15x | n/a | n/a |
| MODMUL_R2 | w32-il64 | 3.33 G | none | n/a | 6.14 M | 54.01 M | 33.53 M | 2.93 G | n/a | 61.66x | 99.31x | 1.14x |

### p1024 (1024-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GMP 40T | OpenSSL | CGBN | GPU vs CPU-CL | GPU vs GMP 40T | GPU vs OpenSSL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | w32-il64 | 1.85 G | none | n/a | 52.29 M | 71.80 M | 72.34 M | 5.23 G | n/a | 25.75x | 25.56x | 0.35x |
| SUBTRACT | w32-il64 | 1.85 G | none | n/a | 56.57 M | 878.03 M | 70.72 M | 5.21 G | n/a | 2.10x | 26.10x | 0.35x |
| ADDMOD | w32-il64 | 5.09 G | none | n/a | 19.14 M | 34.48 M | 34.77 M | 4.89 G | n/a | 147.70x | 146.47x | 1.04x |
| SUBTRACTMOD | w32-il64 | 4.75 G | none | n/a | 23.67 M | 144.55 M | 34.11 M | 4.75 G | n/a | 32.87x | 139.30x | 1.00x |
| MULTIPLYOPERANDSCANNING | w32-il64 | 743.91 M | none | n/a | 6.29 M | 247.24 M | 130.03 M | n/a | n/a | 3.01x | 5.72x | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 295.90 M | none | n/a | 6.29 M | 247.32 M | 131.52 M | 2.52 G | n/a | 1.20x | 2.25x | 0.12x |
| MONTGOMERYMULTIPLICATION | w32-il64 | 1.20 G | none | n/a | 927.73 k | 36.18 M | 135.83 M | 1.62 G | n/a | 33.11x | 8.82x | 0.74x |
| COMPARE | w32-il64 | 12.81 G | none | n/a | 122.35 M | 2.40 G | 1.90 G | 5.19 G | n/a | 5.33x | 6.74x | 2.47x |
| REDUCE | w32-il64 | 757.86 M | none | n/a | 52.67 M | 830.46 M | 44.17 M | 2.18 G | n/a | 0.91x | 17.16x | 0.35x |
| MODMUL | w32-il64 | 52.22 M | none | n/a | 2.14 M | 50.27 M | 20.90 M | 230.40 M | n/a | 1.04x | 2.50x | 0.23x |
| MODEXP | w32-il64 | 678.49 k | none | n/a | 3.12 k | 123.21 k | 95.01 k | 413.64 k | n/a | 5.51x | 7.14x | 1.64x |
| EXPONENTIATION | w32-il64 | 577.67 k | none | n/a | 23.61 k | 702.78 k | 40.60 k | n/a | n/a | 0.82x | 14.23x | n/a |
| DIVIDE | w32-il64 | 152.34 M | none | n/a | 29.53 M | 1.17 G | 64.12 M | 1.64 G | n/a | 0.13x | 2.38x | 0.09x |
| ISQRT | w32-o64 | 4.28 M | none | n/a | 4.49 M | 87.07 M | n/a | n/a | n/a | 0.05x | n/a | n/a |
| MODMUL_R2 | w32-il64 | 605.25 M | none | n/a | 2.13 M | 84.14 M | 26.07 M | 839.33 M | n/a | 7.19x | 23.22x | 0.72x |

### p2048 (2048-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GMP 40T | OpenSSL | CGBN | GPU vs CPU-CL | GPU vs GMP 40T | GPU vs OpenSSL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | w32-il | 935.78 M | none | n/a | 36.24 M | 715.89 M | 520.23 M | 2.73 G | n/a | 1.31x | 1.80x | 0.34x |
| SUBTRACT | w32-il | 937.92 M | none | n/a | 34.77 M | 50.87 M | 1.14 G | 2.74 G | n/a | 18.44x | 0.82x | 0.34x |
| ADDMOD | w32-il | 2.09 G | none | n/a | 15.45 M | 564.44 M | 85.16 M | 2.75 G | n/a | 3.70x | 24.50x | 0.76x |
| SUBTRACTMOD | w32-il64 | 2.04 G | none | n/a | 17.83 M | 526.25 M | 70.30 M | 2.73 G | n/a | 3.88x | 29.03x | 0.75x |
| MULTIPLYOPERANDSCANNING | w32-il64 | 209.59 M | none | n/a | 1.92 M | 75.19 M | 58.45 M | n/a | n/a | 2.79x | 3.59x | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 70.59 M | none | n/a | 1.92 M | 75.32 M | 58.58 M | 769.21 M | n/a | 0.94x | 1.21x | 0.09x |
| MONTGOMERYMULTIPLICATION | w32-il64 | 246.59 M | none | n/a | 282.21 k | 11.20 M | 38.37 M | 490.27 M | n/a | 22.02x | 6.43x | 0.50x |
| COMPARE | w32-il | 6.96 G | none | n/a | 129.29 M | 2.34 G | 2.10 G | 2.74 G | n/a | 2.97x | 3.31x | 2.54x |
| REDUCE | w32-il64 | 160.63 M | none | n/a | 35.93 M | 573.63 M | 77.53 M | 1.83 G | n/a | 0.28x | 2.07x | 0.09x |
| MODMUL | w32-il64 | 15.86 M | none | n/a | 653.27 k | 13.60 M | 11.55 M | 107.17 M | n/a | 1.17x | 1.37x | 0.15x |
| MODEXP | w32-il64 | 69.44 k | none | n/a | 413.9 | 15.95 k | 18.92 k | 71.87 k | n/a | 4.35x | 3.67x | 0.97x |
| EXPONENTIATION | w32-il64 | 70.71 k | none | n/a | 3.90 k | 108.48 k | 9.12 k | n/a | n/a | 0.65x | 7.75x | n/a |
| DIVIDE | w32-il64 | 7.70 M | none | n/a | 23.08 M | 24.10 M | 25.40 M | 1.60 G | n/a | 0.32x | 0.30x | 0.00x |
| ISQRT | w32-il64 | 1.09 M | none | n/a | 2.73 M | 3.81 M | n/a | n/a | n/a | 0.29x | n/a | n/a |
| MODMUL_R2 | w32-il64 | 133.51 M | none | n/a | 650.77 k | 25.79 M | 11.70 M | 252.97 M | n/a | 5.18x | 11.42x | 0.53x |

## 6. Raw data

Also written to `NVIDIA_H100_NVL_Report.csv` for analysis.

```csv
kind,device,device_type,kernel,modulus,bits,operation,items,seconds,ops_per_sec,mismatches
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,secp256k1,256,ADD,500000,0.007345750,68066569.136,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,secp256k1,256,ADD,500000,0.000178899,2794872554.793,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,secp256k1,256,ADD,500000,0.000189369,2640347394.904,0
library,NVIDIA H100 NVL,gpu,cgbn,secp256k1,256,ADD,500000,0.000040672,12293469708.891,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,secp256k1,256,ADD,500000,0.000052209,9576890633.083,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,secp256k1,256,ADD,500000,0.002155365,231979274.992,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,secp256k1,256,ADD,500000,0.000051100,9784739326.408,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,secp256k1,256,ADD,500000,0.002154700,232050866.369,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,secp256k1,256,ADD,500000,0.000036509,13695249819.840,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,secp256k1,256,ADD,500000,0.002304058,217008426.181,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,secp256k1,256,ADD,500000,0.000035930,13915940907.832,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,secp256k1,256,ADD,500000,0.002281835,219121889.813,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,secp256k1,256,SUBTRACT,500000,0.005650650,88485395.307,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,secp256k1,256,SUBTRACT,500000,0.000141530,3532819730.640,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,secp256k1,256,SUBTRACT,500000,0.000270008,1851796631.946,0
library,NVIDIA H100 NVL,gpu,cgbn,secp256k1,256,SUBTRACT,500000,0.000041056,12178487918.940,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,secp256k1,256,SUBTRACT,500000,0.000050939,9815654650.873,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,secp256k1,256,SUBTRACT,500000,0.002164875,230960218.964,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,secp256k1,256,SUBTRACT,500000,0.000051600,9689925212.159,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,secp256k1,256,SUBTRACT,500000,0.002168999,230521084.185,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,secp256k1,256,SUBTRACT,500000,0.000036780,13594335901.348,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,secp256k1,256,SUBTRACT,500000,0.002277968,219493860.475,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,secp256k1,256,SUBTRACT,500000,0.000036060,13865781109.465,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,secp256k1,256,SUBTRACT,500000,0.002291395,218207687.645,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,secp256k1,256,ADDMOD,500000,0.019680530,25405819.854,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,secp256k1,256,ADDMOD,500000,0.000508216,983833577.614,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,secp256k1,256,ADDMOD,500000,0.003170079,157724774.562,0
library,NVIDIA H100 NVL,gpu,cgbn,secp256k1,256,ADDMOD,500000,0.000050816,9839420654.912,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,secp256k1,256,ADDMOD,500000,0.000044719,11180938102.168,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,secp256k1,256,ADDMOD,500000,0.002141945,233432697.213,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,secp256k1,256,ADDMOD,500000,0.000044779,11165953580.997,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,secp256k1,256,ADDMOD,500000,0.002155169,232000364.289,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,secp256k1,256,ADDMOD,500000,0.000030890,16186442464.202,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,secp256k1,256,ADDMOD,500000,0.002293408,218016154.435,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,secp256k1,256,ADDMOD,500000,0.000030800,16233765339.986,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,secp256k1,256,ADDMOD,500000,0.002361125,211763458.268,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,secp256k1,256,SUBTRACTMOD,500000,0.015781276,31683116.146,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,secp256k1,256,SUBTRACTMOD,500000,0.000421927,1185039063.406,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,secp256k1,256,SUBTRACTMOD,500000,0.005284326,94619445.645,0
library,NVIDIA H100 NVL,gpu,cgbn,secp256k1,256,SUBTRACTMOD,500000,0.000055072,9079023823.359,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,secp256k1,256,SUBTRACTMOD,500000,0.000045760,10926570839.444,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,secp256k1,256,SUBTRACTMOD,500000,0.002164584,230991266.420,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,secp256k1,256,SUBTRACTMOD,500000,0.000045200,11061933898.325,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,secp256k1,256,SUBTRACTMOD,500000,0.002159369,231549123.467,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,secp256k1,256,SUBTRACTMOD,500000,0.000030610,16334524853.768,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,secp256k1,256,SUBTRACTMOD,500000,0.002280087,219289880.147,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,secp256k1,256,SUBTRACTMOD,500000,0.000030570,16355923197.027,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,secp256k1,256,SUBTRACTMOD,500000,0.002302856,217121694.907,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,secp256k1,256,MULTIPLYOPERANDSCANNING,500000,0.009035420,55337770.907,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,500000,0.000236268,2116241052.285,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,500000,0.001023613,488465867.832,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,secp256k1,256,MULTIPLYOPERANDSCANNING,500000,0.000103629,4824903314.155,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,secp256k1,256,MULTIPLYOPERANDSCANNING,500000,0.002871469,174126900.841,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,secp256k1,256,MULTIPLYOPERANDSCANNING,500000,0.000106169,4709474005.315,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,secp256k1,256,MULTIPLYOPERANDSCANNING,500000,0.002878376,173709064.222,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,secp256k1,256,MULTIPLYOPERANDSCANNING,500000,0.000036689,13628067674.422,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,secp256k1,256,MULTIPLYOPERANDSCANNING,500000,0.003086283,162007178.840,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,secp256k1,256,MULTIPLYOPERANDSCANNING,500000,0.000030630,16323877685.906,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,secp256k1,256,MULTIPLYOPERANDSCANNING,500000,0.003086360,162003136.548,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,secp256k1,256,MULTIPLYPRODUCTSCANNING,500000,0.009021631,55422351.159,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,500000,0.000235478,2123340208.163,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,500000,0.001005334,497347138.961,0
library,NVIDIA H100 NVL,gpu,cgbn,secp256k1,256,MULTIPLYPRODUCTSCANNING,500000,0.000040544,12332280978.690,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,secp256k1,256,MULTIPLYPRODUCTSCANNING,500000,0.000125969,3969230417.197,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,secp256k1,256,MULTIPLYPRODUCTSCANNING,500000,0.002896010,172651337.376,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,secp256k1,256,MULTIPLYPRODUCTSCANNING,500000,0.000115389,4333168509.580,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,secp256k1,256,MULTIPLYPRODUCTSCANNING,500000,0.002899005,172472968.010,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,secp256k1,256,MULTIPLYPRODUCTSCANNING,500000,0.000095379,5242243277.336,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,secp256k1,256,MULTIPLYPRODUCTSCANNING,500000,0.003149242,158768365.671,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,secp256k1,256,MULTIPLYPRODUCTSCANNING,500000,0.000092129,5427173344.857,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,secp256k1,256,MULTIPLYPRODUCTSCANNING,500000,0.003140520,159209302.239,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,secp256k1,256,MONTGOMERYMULTIPLICATION,500000,0.071624617,6980840.123,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,500000,0.001976106,253022860.981,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,500000,0.000493856,1012440930.222,0
library,NVIDIA H100 NVL,gpu,cgbn,secp256k1,256,MONTGOMERYMULTIPLICATION,500000,0.000061888,8079110651.499,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,secp256k1,256,MONTGOMERYMULTIPLICATION,500000,0.000052889,9453762705.448,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,secp256k1,256,MONTGOMERYMULTIPLICATION,500000,0.002162375,231227236.246,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,secp256k1,256,MONTGOMERYMULTIPLICATION,500000,0.000038609,12950339955.797,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,secp256k1,256,MONTGOMERYMULTIPLICATION,500000,0.002143020,233315602.987,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,secp256k1,256,MONTGOMERYMULTIPLICATION,500000,0.000053750,9302325265.481,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,secp256k1,256,MONTGOMERYMULTIPLICATION,500000,0.002306537,216775191.160,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,secp256k1,256,MONTGOMERYMULTIPLICATION,500000,0.000035849,13947393395.802,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,secp256k1,256,MONTGOMERYMULTIPLICATION,500000,0.002305945,216830846.192,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,secp256k1,256,COMPARE,500000,0.003870815,129171762.278,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,secp256k1,256,COMPARE,500000,0.000219479,2278122587.797,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,secp256k1,256,COMPARE,500000,0.000230949,2164979933.331,0
library,NVIDIA H100 NVL,gpu,cgbn,secp256k1,256,COMPARE,500000,0.000041248,12121799844.841,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,secp256k1,256,COMPARE,500000,0.000029130,17164421204.916,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,secp256k1,256,COMPARE,500000,0.002125665,235220505.075,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,secp256k1,256,COMPARE,500000,0.000031820,15713400380.125,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,secp256k1,256,COMPARE,500000,0.002139618,233686576.144,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,secp256k1,256,COMPARE,500000,0.000014069,35539213795.443,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,secp256k1,256,COMPARE,500000,0.002267627,220494818.063,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,secp256k1,256,COMPARE,500000,0.000013620,36710691020.984,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,secp256k1,256,COMPARE,500000,0.002271406,220127973.281,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,secp256k1,256,REDUCE,62500,0.000950678,65742554.299,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,secp256k1,256,REDUCE,62500,0.000062636,997823117.333,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,secp256k1,256,REDUCE,62500,0.000562240,111162427.479,0
library,NVIDIA H100 NVL,gpu,cgbn,secp256k1,256,REDUCE,500000,0.000086272,5795623145.401,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,secp256k1,256,REDUCE,500000,0.000071399,7002901939.876,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,secp256k1,256,REDUCE,500000,0.002179875,229370946.111,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,secp256k1,256,REDUCE,500000,0.000039489,12661751185.114,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,secp256k1,256,REDUCE,500000,0.002138258,233835210.635,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,secp256k1,256,REDUCE,500000,0.000072730,6874747070.607,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,secp256k1,256,REDUCE,500000,0.002334236,214202843.708,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,secp256k1,256,REDUCE,500000,0.000038970,12830391967.725,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,secp256k1,256,REDUCE,500000,0.002287415,218587359.956,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,secp256k1,256,MODMUL,31250,0.002573301,12143934.983,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,secp256k1,256,MODMUL,31250,0.000144844,215748618.545,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,secp256k1,256,MODMUL,31250,0.000337953,92468484.606,0
library,NVIDIA H100 NVL,gpu,cgbn,secp256k1,256,MODMUL,500000,0.000279744,1787348432.853,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,secp256k1,256,MODMUL,500000,0.000745245,670920280.942,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,secp256k1,256,MODMUL,500000,0.002874990,173913648.766,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,secp256k1,256,MODMUL,500000,0.000569107,878569390.941,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,secp256k1,256,MODMUL,500000,0.002670215,187250836.875,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,secp256k1,256,MODMUL,500000,0.000749436,667168400.179,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,secp256k1,256,MODMUL,500000,0.003004002,166444627.559,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,secp256k1,256,MODMUL,500000,0.000583106,857477096.702,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,secp256k1,256,MODMUL,500000,0.002830562,176643365.783,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,secp256k1,256,MODEXP,7812,0.068860875,113446.133,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,secp256k1,256,MODEXP,7812,0.005505990,1418818.413,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,secp256k1,256,MODEXP,7812,0.008345235,936103.055,0
library,NVIDIA H100 NVL,gpu,cgbn,secp256k1,256,MODEXP,500000,0.094490685,5291526.884,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,secp256k1,256,MODEXP,500000,0.019840700,25200723.729,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,secp256k1,256,MODEXP,500000,0.022092384,22632233.810,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,secp256k1,256,MODEXP,500000,0.008765219,57043640.623,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,secp256k1,256,MODEXP,500000,0.010947597,45672123.380,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,secp256k1,256,MODEXP,500000,0.019831257,25212723.531,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,secp256k1,256,MODEXP,500000,0.022296443,22425101.619,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,secp256k1,256,MODEXP,500000,0.008921802,56042489.795,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,secp256k1,256,MODEXP,500000,0.011293851,44271878.511,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,secp256k1,256,EXPONENTIATION,7812,0.021771842,358812.084,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,secp256k1,256,EXPONENTIATION,7812,0.000591627,13204264.831,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,secp256k1,256,EXPONENTIATION,7812,0.007214191,1082865.699,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,secp256k1,256,EXPONENTIATION,500000,0.003905692,128018286.456,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,secp256k1,256,EXPONENTIATION,500000,0.006010077,83193609.380,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,secp256k1,256,EXPONENTIATION,500000,0.002031629,246107922.326,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,secp256k1,256,EXPONENTIATION,500000,0.004126266,121174931.692,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,secp256k1,256,EXPONENTIATION,500000,0.003871801,129138867.840,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,secp256k1,256,EXPONENTIATION,500000,0.006191756,80752535.430,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,secp256k1,256,EXPONENTIATION,500000,0.002038491,245279478.877,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,secp256k1,256,EXPONENTIATION,500000,0.004326641,115563088.287,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,secp256k1,256,DIVIDE,62500,0.001754740,35617802.462,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,secp256k1,256,DIVIDE,62500,0.000111610,559988059.060,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,secp256k1,256,DIVIDE,62500,0.000425133,147012984.751,0
library,NVIDIA H100 NVL,gpu,cgbn,secp256k1,256,DIVIDE,500000,0.000112576,4441444002.274,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,secp256k1,256,DIVIDE,500000,0.000333858,1497642649.282,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,secp256k1,256,DIVIDE,500000,0.003114548,160536938.887,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,secp256k1,256,DIVIDE,500000,0.000310068,1612549407.839,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,secp256k1,256,DIVIDE,500000,0.003067843,162980962.180,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,secp256k1,256,DIVIDE,500000,0.000326029,1533606062.343,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,secp256k1,256,DIVIDE,500000,0.003338602,149763283.154,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,secp256k1,256,DIVIDE,500000,0.000293898,1701270573.994,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,secp256k1,256,DIVIDE,500000,0.003294396,151772890.335,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,secp256k1,256,ISQRT,15625,0.001014744,15397964.723,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,secp256k1,256,ISQRT,15625,0.000047247,330707108.463,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,secp256k1,256,ISQRT,500000,0.005127344,97516375.150,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,secp256k1,256,ISQRT,500000,0.007291588,68572168.358,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,secp256k1,256,ISQRT,500000,0.003684378,135708116.949,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,secp256k1,256,ISQRT,500000,0.005830747,85752305.584,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,secp256k1,256,ISQRT,500000,0.005152878,97033153.081,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,secp256k1,256,ISQRT,500000,0.007483323,66815236.613,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,secp256k1,256,ISQRT,500000,0.003852663,129780361.036,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,secp256k1,256,ISQRT,500000,0.006141624,81411691.682,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,secp256k1,256,MODMUL_R2,500000,0.041252015,12120620.052,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,secp256k1,256,MODMUL_R2,500000,0.001047704,477234019.011,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,secp256k1,256,MODMUL_R2,500000,0.005593999,89381496.168,0
library,NVIDIA H100 NVL,gpu,cgbn,secp256k1,256,MODMUL_R2,500000,0.000102240,4890453834.116,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,secp256k1,256,MODMUL_R2,500000,0.000086650,5770339150.739,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,secp256k1,256,MODMUL_R2,500000,0.002178584,229506873.129,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,secp256k1,256,MODMUL_R2,500000,0.000052470,9529252022.347,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,secp256k1,256,MODMUL_R2,500000,0.002157507,231748958.552,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,secp256k1,256,MODMUL_R2,500000,0.000088280,5663796455.188,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,secp256k1,256,MODMUL_R2,500000,0.002364494,211461731.413,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,secp256k1,256,MODMUL_R2,500000,0.000052199,9578727471.420,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,secp256k1,256,MODMUL_R2,500000,0.002308850,216558031.994,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,ADD,500000,0.007873007,63508136.253,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,ADD,500000,0.000300368,1664624709.924,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,ADD,500000,0.000191689,2608392465.545,0
library,NVIDIA H100 NVL,gpu,cgbn,rsa256(composite),256,ADD,500000,0.000041024,12187987519.501,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,rsa256(composite),256,ADD,500000,0.000051570,9695557834.048,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,rsa256(composite),256,ADD,500000,0.002163845,231070155.424,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,rsa256(composite),256,ADD,500000,0.000051200,9765628355.326,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,rsa256(composite),256,ADD,500000,0.002154928,232026314.523,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,rsa256(composite),256,ADD,500000,0.000036220,13804524820.210,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,rsa256(composite),256,ADD,500000,0.002294375,217924273.301,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,rsa256(composite),256,ADD,500000,0.000036280,13781690254.554,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,rsa256(composite),256,ADD,500000,0.002282820,219027343.209,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACT,500000,0.006076277,82287229.488,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACT,500000,0.000142339,3512739848.350,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACT,500000,0.000265758,1881410773.579,0
library,NVIDIA H100 NVL,gpu,cgbn,rsa256(composite),256,SUBTRACT,500000,0.000040832,12245297805.643,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,rsa256(composite),256,SUBTRACT,500000,0.000051679,9675114820.497,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,rsa256(composite),256,SUBTRACT,500000,0.002160225,231457370.765,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,rsa256(composite),256,SUBTRACT,500000,0.000051090,9786645618.193,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,rsa256(composite),256,SUBTRACT,500000,0.002155757,231937082.848,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,rsa256(composite),256,SUBTRACT,500000,0.000036560,13676149174.648,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,rsa256(composite),256,SUBTRACT,500000,0.002330257,214568609.947,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,rsa256(composite),256,SUBTRACT,500000,0.000036360,13751376095.796,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,rsa256(composite),256,SUBTRACT,500000,0.002288906,218444972.542,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,ADDMOD,500000,0.016825099,29717507.245,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,ADDMOD,500000,0.000423687,1180116491.986,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,ADDMOD,500000,0.002879544,173638604.696,0
library,NVIDIA H100 NVL,gpu,cgbn,rsa256(composite),256,ADDMOD,500000,0.000054272,9212853773.585,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,rsa256(composite),256,ADDMOD,500000,0.000043850,11402512543.523,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,rsa256(composite),256,ADDMOD,500000,0.002129215,234828329.807,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,rsa256(composite),256,ADDMOD,500000,0.000044860,11145785903.263,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,rsa256(composite),256,ADDMOD,500000,0.002148427,232728408.945,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,rsa256(composite),256,ADDMOD,500000,0.000030440,16425731021.348,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,rsa256(composite),256,ADDMOD,500000,0.002274077,219869420.425,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,rsa256(composite),256,ADDMOD,500000,0.000030860,16202196257.434,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,rsa256(composite),256,ADDMOD,500000,0.002297926,217587511.932,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACTMOD,500000,0.015854054,31537674.815,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACTMOD,500000,0.000402487,1242276149.534,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACTMOD,500000,0.004417106,113196287.719,0
library,NVIDIA H100 NVL,gpu,cgbn,rsa256(composite),256,SUBTRACTMOD,500000,0.000055360,9031791907.514,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,rsa256(composite),256,SUBTRACTMOD,500000,0.000044179,11317600903.307,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,rsa256(composite),256,SUBTRACTMOD,500000,0.002151944,232348049.967,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,rsa256(composite),256,SUBTRACTMOD,500000,0.000045000,11111112739.201,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,rsa256(composite),256,SUBTRACTMOD,500000,0.002147047,232877998.490,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,rsa256(composite),256,SUBTRACTMOD,500000,0.000030940,16160314915.953,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,rsa256(composite),256,SUBTRACTMOD,500000,0.002343207,213382768.428,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,rsa256(composite),256,SUBTRACTMOD,500000,0.000030630,16323877685.906,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,rsa256(composite),256,SUBTRACTMOD,500000,0.002275106,219769982.533,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYOPERANDSCANNING,500000,0.009038761,55317316.369,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,500000,0.000235519,2122970766.050,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,500000,0.000372758,1341352976.679,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,500000,0.000104410,4788813709.826,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,500000,0.002865920,174464045.357,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,500000,0.000107779,4639122951.736,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,500000,0.002881014,173550005.774,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,rsa256(composite),256,MULTIPLYOPERANDSCANNING,500000,0.000035730,13993836461.471,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,rsa256(composite),256,MULTIPLYOPERANDSCANNING,500000,0.003098802,161352677.544,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,500000,0.000030930,16165545844.444,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,500000,0.003092721,161669934.633,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,500000,0.009014991,55463172.425,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,500000,0.000235799,2120450554.000,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,500000,0.001010635,494738430.664,0
library,NVIDIA H100 NVL,gpu,cgbn,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,500000,0.000040320,12400793650.794,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,500000,0.000124389,4019647565.474,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,500000,0.002907730,171955442.983,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,500000,0.000115089,4344463771.752,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,500000,0.002892573,172856487.113,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,500000,0.000095139,5255466945.981,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,500000,0.003135793,159449298.988,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,500000,0.000092389,5411899460.129,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,500000,0.003146261,158918791.047,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MONTGOMERYMULTIPLICATION,500000,0.071686520,6974812.000,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,500000,0.001979129,252636389.381,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,500000,0.000492298,1015645024.277,0
library,NVIDIA H100 NVL,gpu,cgbn,rsa256(composite),256,MONTGOMERYMULTIPLICATION,500000,0.000061472,8133784487.246,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,500000,0.000052950,9442871393.488,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,500000,0.002156105,231899657.685,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,500000,0.000037870,13203055926.665,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,500000,0.002141528,233478151.375,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,rsa256(composite),256,MONTGOMERYMULTIPLICATION,500000,0.000053229,9393378429.183,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,rsa256(composite),256,MONTGOMERYMULTIPLICATION,500000,0.002298927,217492764.627,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,500000,0.000036530,13687371975.888,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,500000,0.002949221,169536294.415,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,COMPARE,500000,0.004048508,123502287.455,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,COMPARE,500000,0.000176159,2838345444.140,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,COMPARE,500000,0.000120700,4142502491.790,0
library,NVIDIA H100 NVL,gpu,cgbn,rsa256(composite),256,COMPARE,500000,0.000041088,12169003115.265,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,rsa256(composite),256,COMPARE,500000,0.000027479,18195698648.128,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,rsa256(composite),256,COMPARE,500000,0.002138405,233819132.583,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,rsa256(composite),256,COMPARE,500000,0.000026530,18846585850.253,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,rsa256(composite),256,COMPARE,500000,0.002115259,236377672.316,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,rsa256(composite),256,COMPARE,500000,0.000013960,35816614999.729,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,rsa256(composite),256,COMPARE,500000,0.002275658,219716669.540,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,rsa256(composite),256,COMPARE,500000,0.000014250,35087738772.042,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,rsa256(composite),256,COMPARE,500000,0.002944892,169785511.078,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,REDUCE,62500,0.001528168,40898635.821,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,REDUCE,62500,0.000089599,697549753.947,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,REDUCE,62500,0.000511020,122304328.367,0
library,NVIDIA H100 NVL,gpu,cgbn,rsa256(composite),256,REDUCE,500000,0.000086240,5797773654.917,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,rsa256(composite),256,REDUCE,500000,0.000072939,6855045872.497,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,rsa256(composite),256,REDUCE,500000,0.002171494,230256219.009,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,rsa256(composite),256,REDUCE,500000,0.000039399,12690671165.807,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,rsa256(composite),256,REDUCE,500000,0.002141209,233512932.841,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,rsa256(composite),256,REDUCE,500000,0.000081489,6135794664.318,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,rsa256(composite),256,REDUCE,500000,0.002318327,215672765.231,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,rsa256(composite),256,REDUCE,500000,0.000039869,12541075755.098,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,rsa256(composite),256,REDUCE,500000,0.002979522,167812153.045,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MODMUL,31250,0.002568901,12164734.934,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MODMUL,31250,0.000145149,215295991.684,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MODMUL,31250,0.000346478,90193312.431,0
library,NVIDIA H100 NVL,gpu,cgbn,rsa256(composite),256,MODMUL,500000,0.000279840,1786735277.301,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,rsa256(composite),256,MODMUL,500000,0.000745604,670597270.328,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,rsa256(composite),256,MODMUL,500000,0.002847360,175601259.858,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,rsa256(composite),256,MODMUL,500000,0.000569217,878399589.982,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,rsa256(composite),256,MODMUL,500000,0.002669577,187295588.791,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,rsa256(composite),256,MODMUL,500000,0.000748176,668291942.641,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,rsa256(composite),256,MODMUL,500000,0.003023953,165346485.081,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,rsa256(composite),256,MODMUL,500000,0.000582306,858655063.240,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,rsa256(composite),256,MODMUL,500000,0.003510299,142438010.299,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MODEXP,7812,0.064595168,120937.839,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MODEXP,7812,0.001633491,4782395.443,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MODEXP,7812,0.003206393,2436382.577,0
library,NVIDIA H100 NVL,gpu,cgbn,rsa256(composite),256,MODEXP,500000,0.092926592,5380591.166,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,rsa256(composite),256,MODEXP,500000,0.019835449,25207395.106,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,rsa256(composite),256,MODEXP,500000,0.022081893,22642986.274,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,rsa256(composite),256,MODEXP,500000,0.008758490,57087466.161,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,rsa256(composite),256,MODEXP,500000,0.010957278,45631771.114,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,rsa256(composite),256,MODEXP,500000,0.019832769,25210801.375,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,rsa256(composite),256,MODEXP,500000,0.022324409,22397009.435,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,rsa256(composite),256,MODEXP,500000,0.008919025,56059939.195,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,rsa256(composite),256,MODEXP,500000,0.011299801,44248566.852,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,EXPONENTIATION,7812,0.021522853,362963.033,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,EXPONENTIATION,7812,0.000593427,13164214.261,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,EXPONENTIATION,7812,0.007041712,1109389.308,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,rsa256(composite),256,EXPONENTIATION,500000,0.003906613,127988104.972,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,rsa256(composite),256,EXPONENTIATION,500000,0.005991637,83449648.120,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,rsa256(composite),256,EXPONENTIATION,500000,0.002030088,246294743.824,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,rsa256(composite),256,EXPONENTIATION,500000,0.004091407,122207348.297,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,rsa256(composite),256,EXPONENTIATION,500000,0.003874399,129052273.630,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,rsa256(composite),256,EXPONENTIATION,500000,0.006200206,80642481.874,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,rsa256(composite),256,EXPONENTIATION,500000,0.002038647,245260703.177,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,rsa256(composite),256,EXPONENTIATION,500000,0.004353043,114862178.470,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,DIVIDE,62500,0.001762675,35457462.254,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,DIVIDE,62500,0.001608527,38855436.736,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,DIVIDE,62500,0.000465247,134337103.308,0
library,NVIDIA H100 NVL,gpu,cgbn,rsa256(composite),256,DIVIDE,500000,0.000112992,4425092041.914,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,rsa256(composite),256,DIVIDE,500000,0.000334557,1494513509.009,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,rsa256(composite),256,DIVIDE,500000,0.003121308,160189254.360,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,rsa256(composite),256,DIVIDE,500000,0.000310838,1608554992.026,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,rsa256(composite),256,DIVIDE,500000,0.003078302,162427207.230,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,rsa256(composite),256,DIVIDE,500000,0.000326268,1532482651.861,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,rsa256(composite),256,DIVIDE,500000,0.003330042,150148256.224,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,rsa256(composite),256,DIVIDE,500000,0.000298498,1675052965.812,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,rsa256(composite),256,DIVIDE,500000,0.003279310,152471098.599,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,ISQRT,15625,0.001017442,15357144.672,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,ISQRT,15625,0.000047272,330532110.419,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,rsa256(composite),256,ISQRT,500000,0.005129614,97473221.718,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,rsa256(composite),256,ISQRT,500000,0.007275478,68724007.023,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,rsa256(composite),256,ISQRT,500000,0.003684379,135708080.501,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,rsa256(composite),256,ISQRT,500000,0.005836267,85671200.087,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,rsa256(composite),256,ISQRT,500000,0.005155912,96976054.646,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,rsa256(composite),256,ISQRT,500000,0.007451739,67098431.766,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,rsa256(composite),256,ISQRT,500000,0.003853276,129759715.383,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,rsa256(composite),256,ISQRT,500000,0.006143353,81388778.975,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MODMUL_R2,500000,0.041255186,12119688.430,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MODMUL_R2,500000,0.001049244,476533597.872,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MODMUL_R2,500000,0.008840402,56558513.929,0
library,NVIDIA H100 NVL,gpu,cgbn,rsa256(composite),256,MODMUL_R2,500000,0.000101856,4908890983.349,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,rsa256(composite),256,MODMUL_R2,500000,0.000088049,5678657335.704,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,rsa256(composite),256,MODMUL_R2,500000,0.002190414,228267350.244,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,rsa256(composite),256,MODMUL_R2,500000,0.000053590,9330100820.164,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,rsa256(composite),256,MODMUL_R2,500000,0.002159407,231545047.769,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,rsa256(composite),256,MODMUL_R2,500000,0.000088679,5638315759.902,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,rsa256(composite),256,MODMUL_R2,500000,0.002334458,214182478.628,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,rsa256(composite),256,MODMUL_R2,500000,0.000051919,9630392575.460,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,rsa256(composite),256,MODMUL_R2,500000,0.002306096,216816649.283,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,ADD,250000,0.004269627,58553123.656,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,ADD,250000,0.000164979,1515344718.322,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,ADD,250000,0.000270008,925898515.576,0
library,NVIDIA H100 NVL,gpu,cgbn,brainpoolP512r1,512,ADD,500000,0.000054400,9191176470.588,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,brainpoolP512r1,512,ADD,500000,0.000110139,4539717061.962,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,brainpoolP512r1,512,ADD,500000,0.004450968,112335114.493,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,brainpoolP512r1,512,ADD,500000,0.000109860,4551245663.700,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,brainpoolP512r1,512,ADD,500000,0.004462205,112052224.091,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,brainpoolP512r1,512,ADD,500000,0.000082099,6090208509.341,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,brainpoolP512r1,512,ADD,500000,0.004525115,110494430.479,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,brainpoolP512r1,512,ADD,500000,0.000084089,5946081019.605,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,brainpoolP512r1,512,ADD,500000,0.004499202,111130818.011,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACT,250000,0.003338132,74892185.313,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACT,250000,0.000216289,1155860823.071,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACT,250000,0.000191979,1302225806.158,0
library,NVIDIA H100 NVL,gpu,cgbn,brainpoolP512r1,512,SUBTRACT,500000,0.000053760,9300595238.095,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,brainpoolP512r1,512,SUBTRACT,500000,0.000110480,4525705884.832,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,brainpoolP512r1,512,SUBTRACT,500000,0.004442080,112559882.307,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,brainpoolP512r1,512,SUBTRACT,500000,0.000108939,4589723958.328,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,brainpoolP512r1,512,SUBTRACT,500000,0.004451484,112322094.204,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,brainpoolP512r1,512,SUBTRACT,500000,0.000081379,6144093751.431,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,brainpoolP512r1,512,SUBTRACT,500000,0.004515605,110727134.283,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,brainpoolP512r1,512,SUBTRACT,500000,0.000083249,6006077855.863,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,brainpoolP512r1,512,SUBTRACT,500000,0.004486943,111434444.018,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,ADDMOD,250000,0.009331950,26789684.938,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,ADDMOD,250000,0.000241439,1035458340.506,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,ADDMOD,250000,0.001358663,184004426.765,0
library,NVIDIA H100 NVL,gpu,cgbn,brainpoolP512r1,512,ADDMOD,500000,0.000067840,7370283018.868,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,brainpoolP512r1,512,ADDMOD,500000,0.000125869,3972383902.612,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,brainpoolP512r1,512,ADDMOD,500000,0.004442210,112556587.364,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,brainpoolP512r1,512,ADDMOD,500000,0.000131900,3790750142.651,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,brainpoolP512r1,512,ADDMOD,500000,0.004477854,111660631.665,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,brainpoolP512r1,512,ADDMOD,500000,0.000052360,9549273729.403,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,brainpoolP512r1,512,ADDMOD,500000,0.004468716,111888963.621,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,brainpoolP512r1,512,ADDMOD,500000,0.000053519,9342480759.672,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,brainpoolP512r1,512,ADDMOD,500000,0.004489113,111380577.989,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACTMOD,250000,0.008704090,28722129.398,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACTMOD,250000,0.000226578,1103372903.319,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACTMOD,250000,0.002143285,116643376.632,0
library,NVIDIA H100 NVL,gpu,cgbn,brainpoolP512r1,512,SUBTRACTMOD,500000,0.000068512,7297991592.714,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,brainpoolP512r1,512,SUBTRACTMOD,500000,0.000126139,3963882064.612,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,brainpoolP512r1,512,SUBTRACTMOD,500000,0.004472550,111793048.590,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,brainpoolP512r1,512,SUBTRACTMOD,500000,0.000126149,3963565645.268,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,brainpoolP512r1,512,SUBTRACTMOD,500000,0.004481504,111569686.548,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,brainpoolP512r1,512,SUBTRACTMOD,500000,0.000052190,9580383365.473,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,brainpoolP512r1,512,SUBTRACTMOD,500000,0.004485915,111459980.664,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,brainpoolP512r1,512,SUBTRACTMOD,500000,0.000053900,9276437337.135,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,brainpoolP512r1,512,SUBTRACTMOD,500000,0.004493942,111260892.537,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,250000,0.010924714,22883894.186,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,250000,0.000282457,885090530.510,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,250000,0.000306037,816894631.012,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,500000,0.000443207,1128141068.351,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,500000,0.006304878,79303675.853,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,500000,0.000443717,1126844364.686,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,500000,0.006261565,79852241.106,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,500000,0.000201749,2478327336.159,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,500000,0.006123256,81655903.377,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,500000,0.000212319,2354946808.228,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,500000,0.006155023,81234464.897,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,250000,0.010918004,22897958.329,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,250000,0.000279307,895072395.340,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,250000,0.000303757,823026225.514,0
library,NVIDIA H100 NVL,gpu,cgbn,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,500000,0.000065920,7584951456.311,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,500000,0.000731456,683568083.223,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,500000,0.006578956,75999900.589,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,500000,0.000731415,683606435.149,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,500000,0.006568133,76125133.089,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,500000,0.000345368,1447731130.518,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,500000,0.006250195,79997504.431,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,500000,0.000323838,1543981999.801,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,500000,0.006261772,79849601.545,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,250000,0.085139102,2936371.116,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,250000,0.002135268,117081321.598,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,250000,0.003570481,70018576.763,0
library,NVIDIA H100 NVL,gpu,cgbn,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,500000,0.000095296,5246809939.557,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,500000,0.000211489,2364188733.605,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,500000,0.004541400,110098207.477,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,500000,0.000111049,4502517598.478,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,500000,0.004459514,112119841.073,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,500000,0.000191989,2604315797.069,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,500000,0.004632624,107930191.527,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,500000,0.000084909,5888655300.584,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,500000,0.004524443,110510842.020,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,COMPARE,250000,0.002616986,95529743.284,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,COMPARE,250000,0.003676739,67995035.865,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,COMPARE,250000,0.003334622,74971015.862,0
library,NVIDIA H100 NVL,gpu,cgbn,brainpoolP512r1,512,COMPARE,500000,0.000053760,9300595238.095,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,brainpoolP512r1,512,COMPARE,500000,0.000100509,4974680243.767,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,brainpoolP512r1,512,COMPARE,500000,0.004286211,116653147.003,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,brainpoolP512r1,512,COMPARE,500000,0.000100620,4969189398.133,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,brainpoolP512r1,512,COMPARE,500000,0.004299655,116288400.570,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,brainpoolP512r1,512,COMPARE,500000,0.000020090,24887973112.672,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,brainpoolP512r1,512,COMPARE,500000,0.004456475,112196299.409,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,brainpoolP512r1,512,COMPARE,500000,0.000019910,25113022888.400,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,brainpoolP512r1,512,COMPARE,500000,0.004449213,112379425.763,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,REDUCE,31250,0.000836145,37373882.557,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,REDUCE,31250,0.000757848,41235192.814,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,REDUCE,31250,0.000786224,39746962.561,0
library,NVIDIA H100 NVL,gpu,cgbn,brainpoolP512r1,512,REDUCE,500000,0.000149792,3337961973.937,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,brainpoolP512r1,512,REDUCE,500000,0.000301747,1657017298.312,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,brainpoolP512r1,512,REDUCE,500000,0.004641688,107719432.607,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,brainpoolP512r1,512,REDUCE,500000,0.000126969,3937968768.323,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,brainpoolP512r1,512,REDUCE,500000,0.004481615,111566924.605,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,brainpoolP512r1,512,REDUCE,500000,0.000292989,1706548605.715,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,brainpoolP512r1,512,REDUCE,500000,0.004743354,105410644.306,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,brainpoolP512r1,512,REDUCE,500000,0.000121269,4123064373.093,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,brainpoolP512r1,512,REDUCE,500000,0.004535122,110250616.859,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL,15625,0.002548966,6129936.600,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL,15625,0.000067955,229933248.505,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL,15625,0.000273878,57050849.205,0
library,NVIDIA H100 NVL,gpu,cgbn,brainpoolP512r1,512,MODMUL,500000,0.000883392,566000144.896,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,brainpoolP512r1,512,MODMUL,500000,0.002694842,185539634.732,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,brainpoolP512r1,512,MODMUL,500000,0.007072952,70691841.318,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,brainpoolP512r1,512,MODMUL,500000,0.002101388,237937972.062,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,brainpoolP512r1,512,MODMUL,500000,0.006477544,77189749.386,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,brainpoolP512r1,512,MODMUL,500000,0.002696425,185430708.237,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,brainpoolP512r1,512,MODMUL,500000,0.007145001,69978996.789,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,brainpoolP512r1,512,MODMUL,500000,0.002096187,238528333.730,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,brainpoolP512r1,512,MODMUL,500000,0.006587321,75903390.807,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MODEXP,3906,0.198694477,19658.322,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MODEXP,3906,0.007993336,488657.051,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MODEXP,3906,0.004632015,843261.526,0
library,NVIDIA H100 NVL,gpu,cgbn,brainpoolP512r1,512,MODEXP,500000,0.226987049,2202768.846,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,brainpoolP512r1,512,MODEXP,500000,0.169435127,2950981.941,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,brainpoolP512r1,512,MODEXP,500000,0.174063937,2872507.704,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,brainpoolP512r1,512,MODEXP,500000,0.069878203,7155307.073,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,brainpoolP512r1,512,MODEXP,500000,0.074452398,6715700.413,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,brainpoolP512r1,512,MODEXP,500000,0.168718213,2963521.194,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,brainpoolP512r1,512,MODEXP,500000,0.173514138,2881609.567,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,brainpoolP512r1,512,MODEXP,500000,0.070137242,7128880.257,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,brainpoolP512r1,512,MODEXP,500000,0.076380901,6546139.069,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,EXPONENTIATION,3906,0.034030813,114778.333,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,EXPONENTIATION,3906,0.001286923,3035146.596,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,EXPONENTIATION,3906,0.012754030,306256.140,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,brainpoolP512r1,512,EXPONENTIATION,500000,0.117570492,4252767.779,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,brainpoolP512r1,512,EXPONENTIATION,500000,0.122712888,4074551.648,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,brainpoolP512r1,512,EXPONENTIATION,500000,0.107058902,4670326.245,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,brainpoolP512r1,512,EXPONENTIATION,500000,0.111448568,4486374.379,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,brainpoolP512r1,512,EXPONENTIATION,500000,0.113069521,4422058.177,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,brainpoolP512r1,512,EXPONENTIATION,500000,0.118048375,4235551.739,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,brainpoolP512r1,512,EXPONENTIATION,500000,0.109116242,4582269.245,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,brainpoolP512r1,512,EXPONENTIATION,500000,0.115430211,4331621.641,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,DIVIDE,31250,0.000939127,33275575.915,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,DIVIDE,31250,0.000031620,988298478.998,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,DIVIDE,31250,0.000224491,139203645.886,0
library,NVIDIA H100 NVL,gpu,cgbn,brainpoolP512r1,512,DIVIDE,500000,0.000213632,2340473337.328,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,brainpoolP512r1,512,DIVIDE,500000,0.001126352,443910962.462,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,brainpoolP512r1,512,DIVIDE,500000,0.007964376,62779557.283,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,brainpoolP512r1,512,DIVIDE,500000,0.001201573,416121200.219,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,brainpoolP512r1,512,DIVIDE,500000,0.006824832,73261876.606,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,brainpoolP512r1,512,DIVIDE,500000,0.000988395,505870622.750,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,brainpoolP512r1,512,DIVIDE,500000,0.006746474,74112788.254,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,brainpoolP512r1,512,DIVIDE,500000,0.001062887,470416908.763,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,brainpoolP512r1,512,DIVIDE,500000,0.008296782,60264329.162,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,ISQRT,7812,0.000906671,8616135.276,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,ISQRT,7812,0.000040302,193837510.616,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,brainpoolP512r1,512,ISQRT,500000,0.020882118,23943931.393,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,brainpoolP512r1,512,ISQRT,500000,0.026621940,18781501.258,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,brainpoolP512r1,512,ISQRT,500000,0.019158261,26098402.128,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,brainpoolP512r1,512,ISQRT,500000,0.023469056,21304648.979,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,brainpoolP512r1,512,ISQRT,500000,0.020397520,24512783.857,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,brainpoolP512r1,512,ISQRT,500000,0.024959985,20032063.327,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,brainpoolP512r1,512,ISQRT,500000,0.017670710,28295410.805,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,brainpoolP512r1,512,ISQRT,500000,0.023685149,21110274.650,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL_R2,250000,0.040727076,6138422.507,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL_R2,250000,0.004628795,54009736.590,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL_R2,250000,0.007454909,33534949.766,0
library,NVIDIA H100 NVL,gpu,cgbn,brainpoolP512r1,512,MODMUL_R2,500000,0.000170464,2933170640.135,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,brainpoolP512r1,512,MODMUL_R2,500000,0.000371888,1344490644.215,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,brainpoolP512r1,512,MODMUL_R2,500000,0.006158999,81182023.343,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,brainpoolP512r1,512,MODMUL_R2,500000,0.000172409,2900080687.004,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,brainpoolP512r1,512,MODMUL_R2,500000,0.004500844,111090275.710,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,brainpoolP512r1,512,MODMUL_R2,500000,0.000359428,1391099195.683,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,brainpoolP512r1,512,MODMUL_R2,500000,0.004784844,104496614.710,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,brainpoolP512r1,512,MODMUL_R2,500000,0.000150129,3330468588.479,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,brainpoolP512r1,512,MODMUL_R2,500000,0.006060462,82501960.560,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,p1024,1024,ADD,125000,0.002390452,52291365.372,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,p1024,1024,ADD,125000,0.001740950,71799878.790,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,p1024,1024,ADD,125000,0.001727965,72339407.835,0
library,NVIDIA H100 NVL,gpu,cgbn,p1024,1024,ADD,500000,0.000095552,5232752846.618,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,p1024,1024,ADD,500000,0.000298668,1674099724.483,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,p1024,1024,ADD,500000,0.012370396,40419077.941,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,p1024,1024,ADD,500000,0.000298058,1677525670.244,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,p1024,1024,ADD,500000,0.009493236,52669079.207,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,p1024,1024,ADD,500000,0.000272238,1836628333.511,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,p1024,1024,ADD,500000,0.009634678,51895870.210,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,p1024,1024,ADD,500000,0.000270469,1848640709.746,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,p1024,1024,ADD,500000,0.012472603,40087862.960,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,p1024,1024,SUBTRACT,125000,0.002209477,56574462.228,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,p1024,1024,SUBTRACT,125000,0.000142365,878027874.275,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,p1024,1024,SUBTRACT,125000,0.001767610,70716937.160,0
library,NVIDIA H100 NVL,gpu,cgbn,p1024,1024,SUBTRACT,500000,0.000096000,5208333333.333,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,p1024,1024,SUBTRACT,500000,0.000298968,1672419833.142,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,p1024,1024,SUBTRACT,500000,0.012395506,40337199.731,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,p1024,1024,SUBTRACT,500000,0.000299379,1670123841.306,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,p1024,1024,SUBTRACT,500000,0.009444367,52941610.457,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,p1024,1024,SUBTRACT,500000,0.000273098,1830844571.835,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,p1024,1024,SUBTRACT,500000,0.009586198,52158321.664,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,p1024,1024,SUBTRACT,500000,0.000270859,1845978978.459,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,p1024,1024,SUBTRACT,500000,0.012541132,39868809.324,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,p1024,1024,ADDMOD,125000,0.006529144,19144929.396,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,p1024,1024,ADDMOD,125000,0.003625010,34482663.823,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,p1024,1024,ADDMOD,125000,0.003594720,34773223.148,0
library,NVIDIA H100 NVL,gpu,cgbn,p1024,1024,ADDMOD,500000,0.000102176,4893517068.588,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,p1024,1024,ADDMOD,500000,0.000373068,1340238332.076,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,p1024,1024,ADDMOD,500000,0.012490125,40031625.094,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,p1024,1024,ADDMOD,500000,0.000374107,1336516090.205,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,p1024,1024,ADDMOD,500000,0.009538406,52419659.834,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,p1024,1024,ADDMOD,500000,0.000110419,4528206137.743,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,p1024,1024,ADDMOD,500000,0.009471199,52791626.759,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,p1024,1024,ADDMOD,500000,0.000098169,5093256271.008,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,p1024,1024,ADDMOD,500000,0.012339323,40520861.631,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,p1024,1024,SUBTRACTMOD,125000,0.005279979,23674336.457,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,p1024,1024,SUBTRACTMOD,125000,0.000864745,144551283.686,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,p1024,1024,SUBTRACTMOD,125000,0.003664669,34109492.575,0
library,NVIDIA H100 NVL,gpu,cgbn,p1024,1024,SUBTRACTMOD,500000,0.000105312,4747797022.182,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,p1024,1024,SUBTRACTMOD,500000,0.000376067,1329550497.378,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,p1024,1024,SUBTRACTMOD,500000,0.012488485,40036881.902,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,p1024,1024,SUBTRACTMOD,500000,0.000376738,1327182194.920,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,p1024,1024,SUBTRACTMOD,500000,0.009518686,52528258.718,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,p1024,1024,SUBTRACTMOD,500000,0.000107499,4651205446.771,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,p1024,1024,SUBTRACTMOD,500000,0.009394429,53223032.593,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,p1024,1024,SUBTRACTMOD,500000,0.000105229,4751540441.227,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,p1024,1024,SUBTRACTMOD,500000,0.012367033,40430069.226,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,p1024,1024,MULTIPLYOPERANDSCANNING,125000,0.019886494,6285673.089,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,p1024,1024,MULTIPLYOPERANDSCANNING,125000,0.000505587,247237356.001,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,p1024,1024,MULTIPLYOPERANDSCANNING,125000,0.000961304,130031706.664,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,p1024,1024,MULTIPLYOPERANDSCANNING,500000,0.001144592,436836877.816,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,p1024,1024,MULTIPLYOPERANDSCANNING,500000,0.016048596,31155373.327,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,p1024,1024,MULTIPLYOPERANDSCANNING,500000,0.001140453,438422283.904,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,p1024,1024,MULTIPLYOPERANDSCANNING,500000,0.013213446,37840242.557,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,p1024,1024,MULTIPLYOPERANDSCANNING,500000,0.000728136,686684917.542,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,p1024,1024,MULTIPLYOPERANDSCANNING,500000,0.012928470,38674336.622,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,p1024,1024,MULTIPLYOPERANDSCANNING,500000,0.000672125,743909280.135,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,p1024,1024,MULTIPLYOPERANDSCANNING,500000,0.015685840,31875882.872,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,p1024,1024,MULTIPLYPRODUCTSCANNING,125000,0.019877024,6288667.762,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,125000,0.000505427,247315628.739,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,125000,0.000950394,131524404.133,0
library,NVIDIA H100 NVL,gpu,cgbn,p1024,1024,MULTIPLYPRODUCTSCANNING,500000,0.000198560,2518130539.887,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,p1024,1024,MULTIPLYPRODUCTSCANNING,500000,0.004903615,101965590.969,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,p1024,1024,MULTIPLYPRODUCTSCANNING,500000,0.019829870,25214487.000,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,p1024,1024,MULTIPLYPRODUCTSCANNING,500000,0.004902972,101978962.447,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,p1024,1024,MULTIPLYPRODUCTSCANNING,500000,0.017023184,29371708.628,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,p1024,1024,MULTIPLYPRODUCTSCANNING,500000,0.001840630,271646117.868,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,p1024,1024,MULTIPLYPRODUCTSCANNING,500000,0.014040707,35610742.321,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,p1024,1024,MULTIPLYPRODUCTSCANNING,500000,0.001689766,295898957.074,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,p1024,1024,MULTIPLYPRODUCTSCANNING,500000,0.016757581,29837242.066,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,p1024,1024,MONTGOMERYMULTIPLICATION,125000,0.134737275,927731.395,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,p1024,1024,MONTGOMERYMULTIPLICATION,125000,0.003454660,36183010.698,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,p1024,1024,MONTGOMERYMULTIPLICATION,125000,0.000920284,135827642.107,0
library,NVIDIA H100 NVL,gpu,cgbn,p1024,1024,MONTGOMERYMULTIPLICATION,500000,0.000308064,1623039368.443,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,p1024,1024,MONTGOMERYMULTIPLICATION,500000,0.001000633,499683701.415,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,p1024,1024,MONTGOMERYMULTIPLICATION,500000,0.013114977,38124351.896,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,p1024,1024,MONTGOMERYMULTIPLICATION,500000,0.000677856,737619815.144,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,p1024,1024,MONTGOMERYMULTIPLICATION,500000,0.009871864,50648996.011,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,p1024,1024,MONTGOMERYMULTIPLICATION,500000,0.000723476,691107931.043,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,p1024,1024,MONTGOMERYMULTIPLICATION,500000,0.010038927,49806119.604,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,p1024,1024,MONTGOMERYMULTIPLICATION,500000,0.000417357,1198015072.232,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,p1024,1024,MONTGOMERYMULTIPLICATION,500000,0.012619927,39619880.522,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,p1024,1024,COMPARE,125000,0.001021619,122354813.257,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,p1024,1024,COMPARE,125000,0.000052040,2402021010.434,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,p1024,1024,COMPARE,125000,0.000065715,1902153411.723,0
library,NVIDIA H100 NVL,gpu,cgbn,p1024,1024,COMPARE,500000,0.000096256,5194481382.979,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,p1024,1024,COMPARE,500000,0.000275518,1814763539.581,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,p1024,1024,COMPARE,500000,0.012071754,41419001.854,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,p1024,1024,COMPARE,500000,0.000274408,1822104135.211,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,p1024,1024,COMPARE,500000,0.009181939,54454728.964,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,p1024,1024,COMPARE,500000,0.000039710,12591298001.208,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,p1024,1024,COMPARE,500000,0.009351371,53468095.904,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,p1024,1024,COMPARE,500000,0.000039019,12814276005.794,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,p1024,1024,COMPARE,500000,0.012248319,40821928.327,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,p1024,1024,REDUCE,15625,0.000296646,52672154.613,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,p1024,1024,REDUCE,15625,0.000018815,830458634.763,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,p1024,1024,REDUCE,15625,0.000353737,44171247.539,0
library,NVIDIA H100 NVL,gpu,cgbn,p1024,1024,REDUCE,500000,0.000229088,2182567397.681,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,p1024,1024,REDUCE,500000,0.001314211,380456411.743,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,p1024,1024,REDUCE,500000,0.013425240,37243282.051,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,p1024,1024,REDUCE,500000,0.000846375,590754680.202,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,p1024,1024,REDUCE,500000,0.010040034,49800628.085,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,p1024,1024,REDUCE,500000,0.001138054,439346482.322,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,p1024,1024,REDUCE,500000,0.010492515,47653017.442,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,p1024,1024,REDUCE,500000,0.000659755,757857056.131,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,p1024,1024,REDUCE,500000,0.012859825,38880777.791,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,p1024,1024,MODMUL,7812,0.003652374,2138883.111,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,p1024,1024,MODMUL,7812,0.000155399,50270593.599,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,p1024,1024,MODMUL,7812,0.000373728,20902903.135,0
library,NVIDIA H100 NVL,gpu,cgbn,p1024,1024,MODMUL,500000,0.002170112,230402854.784,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,p1024,1024,MODMUL,500000,0.012437830,40199938.288,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,p1024,1024,MODMUL,500000,0.024743431,20207383.515,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,p1024,1024,MODMUL,500000,0.009657736,51771968.118,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,p1024,1024,MODMUL,500000,0.019006123,26307311.664,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,p1024,1024,MODMUL,500000,0.012421685,40252187.885,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,p1024,1024,MODMUL,500000,0.021910575,22820030.921,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,p1024,1024,MODMUL,500000,0.009574051,52224497.155,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,p1024,1024,MODMUL,500000,0.022023539,22702981.568,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,p1024,1024,MODEXP,1953,0.625928962,3120.162,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,p1024,1024,MODEXP,1953,0.015851468,123206.254,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,p1024,1024,MODEXP,1953,0.020556050,95008.525,0
library,NVIDIA H100 NVL,gpu,cgbn,p1024,1024,MODEXP,500000,1.208781123,413639.815,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,p1024,1024,MODEXP,500000,1.365102335,366272.907,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,p1024,1024,MODEXP,500000,1.377453672,362988.615,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,p1024,1024,MODEXP,500000,0.739499807,676132.698,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,p1024,1024,MODEXP,500000,0.749025977,667533.591,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,p1024,1024,MODEXP,500000,1.373088800,364142.508,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,p1024,1024,MODEXP,500000,1.382922621,361553.128,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,p1024,1024,MODEXP,500000,0.736928246,678492.109,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,p1024,1024,MODEXP,500000,0.749645872,666981.596,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,p1024,1024,EXPONENTIATION,1953,0.082709748,23612.694,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,p1024,1024,EXPONENTIATION,1953,0.002778952,702782.917,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,p1024,1024,EXPONENTIATION,1953,0.048103927,40599.596,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,p1024,1024,EXPONENTIATION,500000,0.991129641,504474.873,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,p1024,1024,EXPONENTIATION,500000,1.003393738,498308.870,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,p1024,1024,EXPONENTIATION,500000,0.865665216,577590.494,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,p1024,1024,EXPONENTIATION,500000,0.874845534,571529.465,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,p1024,1024,EXPONENTIATION,500000,1.003094347,498457.599,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,p1024,1024,EXPONENTIATION,500000,1.012827094,493667.678,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,p1024,1024,EXPONENTIATION,500000,0.865547301,577669.180,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,p1024,1024,EXPONENTIATION,500000,0.878166463,569368.134,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,p1024,1024,DIVIDE,15625,0.000529082,29532260.472,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,p1024,1024,DIVIDE,15625,0.000013301,1174698060.329,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,p1024,1024,DIVIDE,15625,0.000243701,64115378.790,0
library,NVIDIA H100 NVL,gpu,cgbn,p1024,1024,DIVIDE,500000,0.000304512,1641971416.562,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,p1024,1024,DIVIDE,500000,0.003953508,126469960.829,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,p1024,1024,DIVIDE,500000,0.018288780,27339166.412,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,p1024,1024,DIVIDE,500000,0.003688099,135571198.827,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,p1024,1024,DIVIDE,500000,0.015285655,32710407.278,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,p1024,1024,DIVIDE,500000,0.003625353,137917604.794,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,p1024,1024,DIVIDE,500000,0.015512125,32232850.141,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,p1024,1024,DIVIDE,500000,0.003282069,152342927.783,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,p1024,1024,DIVIDE,500000,0.018056454,27690929.849,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,p1024,1024,ISQRT,3906,0.000869691,4491252.388,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,p1024,1024,ISQRT,3906,0.000044860,87071670.329,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,p1024,1024,ISQRT,500000,0.145902710,3426941.144,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,p1024,1024,ISQRT,500000,0.158175673,3161042.343,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,p1024,1024,ISQRT,500000,0.116699413,4284511.697,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,p1024,1024,ISQRT,500000,0.125992407,3968493.118,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,p1024,1024,ISQRT,500000,0.134488646,3717785.961,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,p1024,1024,ISQRT,500000,0.143995930,3472320.363,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,p1024,1024,ISQRT,500000,0.117192394,4266488.490,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,p1024,1024,ISQRT,500000,0.129627089,3857218.455,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,p1024,1024,MODMUL_R2,125000,0.058612481,2132651.577,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,p1024,1024,MODMUL_R2,125000,0.001485601,84141033.906,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,p1024,1024,MODMUL_R2,125000,0.004795700,26065016.465,0
library,NVIDIA H100 NVL,gpu,cgbn,p1024,1024,MODMUL_R2,500000,0.000595712,839331757.628,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,p1024,1024,MODMUL_R2,500000,0.001617091,309197192.425,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,p1024,1024,MODMUL_R2,500000,0.013706944,36477861.122,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,p1024,1024,MODMUL_R2,500000,0.001073757,465654700.695,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,p1024,1024,MODMUL_R2,500000,0.010213311,48955720.604,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,p1024,1024,MODMUL_R2,500000,0.001411790,354160314.216,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,p1024,1024,MODMUL_R2,500000,0.010704837,46707857.213,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,p1024,1024,MODMUL_R2,500000,0.000826105,605249930.808,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,p1024,1024,MODMUL_R2,500000,0.013075741,38238750.653,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,p2048,2048,ADD,62500,0.001724799,36236096.188,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,p2048,2048,ADD,62500,0.000087304,715885392.127,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,p2048,2048,ADD,62500,0.000120139,520230683.306,0
library,NVIDIA H100 NVL,gpu,cgbn,p2048,2048,ADD,500000,0.000183264,2728304522.438,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,p2048,2048,ADD,500000,0.000584337,855670633.429,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,p2048,2048,ADD,500000,0.025335317,19735296.816,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,p2048,2048,ADD,500000,0.000610837,818548955.318,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,p2048,2048,ADD,500000,0.019409949,25759985.229,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,p2048,2048,ADD,500000,0.000534316,935775831.286,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,p2048,2048,ADD,500000,0.019339960,25853207.542,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,p2048,2048,ADD,500000,0.000549387,910105339.085,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,p2048,2048,ADD,500000,0.025278368,19779757.908,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,p2048,2048,SUBTRACT,62500,0.001797665,34767316.641,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,p2048,2048,SUBTRACT,62500,0.001228596,50871089.002,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,p2048,2048,SUBTRACT,62500,0.000054943,1137542808.309,0
library,NVIDIA H100 NVL,gpu,cgbn,p2048,2048,SUBTRACT,500000,0.000182624,2737865778.868,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,p2048,2048,SUBTRACT,500000,0.000584767,855041457.583,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,p2048,2048,SUBTRACT,500000,0.025292347,19768825.716,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,p2048,2048,SUBTRACT,500000,0.000611466,817706944.926,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,p2048,2048,SUBTRACT,500000,0.019587288,25526759.972,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,p2048,2048,SUBTRACT,500000,0.000533096,937917326.380,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,p2048,2048,SUBTRACT,500000,0.019267486,25950453.524,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,p2048,2048,SUBTRACT,500000,0.000557526,896819138.078,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,p2048,2048,SUBTRACT,500000,0.025275288,19782168.278,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,p2048,2048,ADDMOD,62500,0.004044509,15453048.166,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,p2048,2048,ADDMOD,62500,0.000110730,564438512.268,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,p2048,2048,ADDMOD,62500,0.000733900,85161466.100,0
library,NVIDIA H100 NVL,gpu,cgbn,p2048,2048,ADDMOD,500000,0.000181792,2750396057.032,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,p2048,2048,ADDMOD,500000,0.000746766,669553792.548,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,p2048,2048,ADDMOD,500000,0.025643476,19498136.685,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,p2048,2048,ADDMOD,500000,0.000762277,655929544.536,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,p2048,2048,ADDMOD,500000,0.019653458,25440815.589,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,p2048,2048,ADDMOD,500000,0.000239618,2086654538.229,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,p2048,2048,ADDMOD,500000,0.019008446,26304096.584,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,p2048,2048,ADDMOD,500000,0.000267648,1868125553.782,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,p2048,2048,ADDMOD,500000,0.024992750,20005801.708,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,p2048,2048,SUBTRACTMOD,62500,0.003505938,17826898.312,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,p2048,2048,SUBTRACTMOD,62500,0.000118764,526251550.402,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,p2048,2048,SUBTRACTMOD,62500,0.000889084,70297032.117,0
library,NVIDIA H100 NVL,gpu,cgbn,p2048,2048,SUBTRACTMOD,500000,0.000183200,2729257641.921,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,p2048,2048,SUBTRACTMOD,500000,0.000736866,678549396.639,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,p2048,2048,SUBTRACTMOD,500000,0.019498831,25642562.866,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,p2048,2048,SUBTRACTMOD,500000,0.000761576,656533296.810,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,p2048,2048,SUBTRACTMOD,500000,0.019522448,25611542.143,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,p2048,2048,SUBTRACTMOD,500000,0.000269418,1855852242.790,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,p2048,2048,SUBTRACTMOD,500000,0.018987696,26332842.099,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,p2048,2048,SUBTRACTMOD,500000,0.000245018,2040666491.027,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,p2048,2048,SUBTRACTMOD,500000,0.024968540,20025199.707,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,p2048,2048,MULTIPLYOPERANDSCANNING,62500,0.032594835,1917481.713,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,p2048,2048,MULTIPLYOPERANDSCANNING,62500,0.000831255,75187518.504,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,p2048,2048,MULTIPLYOPERANDSCANNING,62500,0.001069373,58445463.694,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,p2048,2048,MULTIPLYOPERANDSCANNING,500000,0.003346481,149410681.554,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,p2048,2048,MULTIPLYOPERANDSCANNING,500000,0.028064383,17816176.456,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,p2048,2048,MULTIPLYOPERANDSCANNING,500000,0.003368246,148445216.204,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,p2048,2048,MULTIPLYOPERANDSCANNING,500000,0.028101548,17792614.111,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,p2048,2048,MULTIPLYOPERANDSCANNING,500000,0.002818130,177422616.186,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,p2048,2048,MULTIPLYOPERANDSCANNING,500000,0.027500813,18181280.669,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,p2048,2048,MULTIPLYOPERANDSCANNING,500000,0.002385555,209594831.367,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,p2048,2048,MULTIPLYOPERANDSCANNING,500000,0.033053751,15126876.226,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,p2048,2048,MULTIPLYPRODUCTSCANNING,62500,0.032582235,1918223.228,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,62500,0.000829825,75317084.640,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,62500,0.001066963,58577475.691,0
library,NVIDIA H100 NVL,gpu,cgbn,p2048,2048,MULTIPLYPRODUCTSCANNING,500000,0.000650016,769211834.786,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,p2048,2048,MULTIPLYPRODUCTSCANNING,500000,0.019036967,26264688.045,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,p2048,2048,MULTIPLYPRODUCTSCANNING,500000,0.043738886,11431475.411,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,p2048,2048,MULTIPLYPRODUCTSCANNING,500000,0.018994820,26322965.927,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,p2048,2048,MULTIPLYPRODUCTSCANNING,500000,0.043658108,11452626.398,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,p2048,2048,MULTIPLYPRODUCTSCANNING,500000,0.007437044,67231012.771,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,p2048,2048,MULTIPLYPRODUCTSCANNING,500000,0.032225839,15515499.827,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,p2048,2048,MULTIPLYPRODUCTSCANNING,500000,0.007082916,70592394.618,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,p2048,2048,MULTIPLYPRODUCTSCANNING,500000,0.037998241,13158503.849,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,p2048,2048,MONTGOMERYMULTIPLICATION,62500,0.221467565,282208.368,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,p2048,2048,MONTGOMERYMULTIPLICATION,62500,0.005580315,11200084.665,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,p2048,2048,MONTGOMERYMULTIPLICATION,62500,0.001628940,38368508.361,0
library,NVIDIA H100 NVL,gpu,cgbn,p2048,2048,MONTGOMERYMULTIPLICATION,500000,0.001019840,490272983.997,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,p2048,2048,MONTGOMERYMULTIPLICATION,500000,0.003239212,154358527.562,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,p2048,2048,MONTGOMERYMULTIPLICATION,500000,0.028108652,17788117.349,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,p2048,2048,MONTGOMERYMULTIPLICATION,500000,0.002338034,213854889.763,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,p2048,2048,MONTGOMERYMULTIPLICATION,500000,0.021090368,23707504.790,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,p2048,2048,MONTGOMERYMULTIPLICATION,500000,0.002860989,174764739.580,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,p2048,2048,MONTGOMERYMULTIPLICATION,500000,0.021674608,23068467.954,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,p2048,2048,MONTGOMERYMULTIPLICATION,500000,0.002027617,246594893.093,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,p2048,2048,MONTGOMERYMULTIPLICATION,500000,0.026739692,18698794.310,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,p2048,2048,COMPARE,62500,0.000483422,129286626.962,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,p2048,2048,COMPARE,62500,0.000026660,2344358157.337,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,p2048,2048,COMPARE,62500,0.000029732,2102094179.431,0
library,NVIDIA H100 NVL,gpu,cgbn,p2048,2048,COMPARE,500000,0.000182560,2738825591.586,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,p2048,2048,COMPARE,500000,0.000519587,962302718.543,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,p2048,2048,COMPARE,500000,0.024633061,20297923.977,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,p2048,2048,COMPARE,500000,0.000545760,916153562.585,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,p2048,2048,COMPARE,500000,0.018701665,26735587.465,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,p2048,2048,COMPARE,500000,0.000071889,6955170461.271,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,p2048,2048,COMPARE,500000,0.018904235,26449099.842,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,p2048,2048,COMPARE,500000,0.000086190,5801136728.904,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,p2048,2048,COMPARE,500000,0.024752158,20200258.925,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,p2048,2048,REDUCE,7812,0.000217424,35929854.755,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,p2048,2048,REDUCE,7812,0.000013618,573634378.142,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,p2048,2048,REDUCE,7812,0.000100763,77528512.954,0
library,NVIDIA H100 NVL,gpu,cgbn,p2048,2048,REDUCE,500000,0.000272992,1831555503.458,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,p2048,2048,REDUCE,500000,0.004523285,110539131.916,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,p2048,2048,REDUCE,500000,0.029261824,17087109.826,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,p2048,2048,REDUCE,500000,0.003372602,148253485.429,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,p2048,2048,REDUCE,500000,0.022360211,22361148.554,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,p2048,2048,REDUCE,500000,0.004392899,113820052.413,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,p2048,2048,REDUCE,500000,0.023475142,21299125.669,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,p2048,2048,REDUCE,500000,0.003112711,160631680.234,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,p2048,2048,REDUCE,500000,0.027858778,17947664.492,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,p2048,2048,MODMUL,3906,0.005979152,653269.896,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,p2048,2048,MODMUL,3906,0.000287168,13601793.173,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,p2048,2048,MODMUL,3906,0.000338208,11549104.280,0
library,NVIDIA H100 NVL,gpu,cgbn,p2048,2048,MODMUL,500000,0.004665568,107168087.573,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,p2048,2048,MODMUL,500000,0.047459717,10535250.355,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,p2048,2048,MODMUL,500000,0.072528976,6893796.486,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,p2048,2048,MODMUL,500000,0.031691362,15777169.800,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,p2048,2048,MODMUL,500000,0.050770760,9848188.200,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,p2048,2048,MODMUL,500000,0.043998795,11363947.586,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,p2048,2048,MODMUL,500000,0.063122183,7921145.570,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,p2048,2048,MODMUL,500000,0.031518157,15863871.712,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,p2048,2048,MODMUL,500000,0.056544713,8842559.696,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,p2048,2048,MODEXP,976,2.358257644,413.865,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,p2048,2048,MODEXP,976,0.061185500,15951.492,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,p2048,2048,MODEXP,976,0.051573156,18924.574,0
library,NVIDIA H100 NVL,gpu,cgbn,p2048,2048,MODEXP,500000,6.957280636,71867.160,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,p2048,2048,MODEXP,500000,10.745799155,46529.811,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,p2048,2048,MODEXP,500000,10.761588281,46461.543,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,p2048,2048,MODEXP,500000,7.240888903,69052.295,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,p2048,2048,MODEXP,500000,7.261649981,68854.875,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,p2048,2048,MODEXP,500000,11.110361513,45003.036,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,p2048,2048,MODEXP,500000,11.129459804,44925.810,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,p2048,2048,MODEXP,500000,7.200777227,69436.949,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,p2048,2048,MODEXP,500000,7.221244849,69240.139,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,p2048,2048,EXPONENTIATION,976,0.250527088,3895.786,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,p2048,2048,EXPONENTIATION,976,0.008996663,108484.668,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,p2048,2048,EXPONENTIATION,976,0.107010780,9120.576,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,p2048,2048,EXPONENTIATION,500000,8.227823265,60769.414,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,p2048,2048,EXPONENTIATION,500000,8.245878758,60636.351,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,p2048,2048,EXPONENTIATION,500000,7.076301531,70658.380,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,p2048,2048,EXPONENTIATION,500000,7.095226257,70469.916,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,p2048,2048,EXPONENTIATION,500000,8.334121002,59994.329,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,p2048,2048,EXPONENTIATION,500000,8.349462000,59884.098,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,p2048,2048,EXPONENTIATION,500000,7.071025205,70711.104,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,p2048,2048,EXPONENTIATION,500000,7.096105519,70461.184,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,p2048,2048,DIVIDE,7812,0.000338526,23076508.058,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,p2048,2048,DIVIDE,7812,0.000324098,24103827.010,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,p2048,2048,DIVIDE,7812,0.000307505,25404434.982,0
library,NVIDIA H100 NVL,gpu,cgbn,p2048,2048,DIVIDE,500000,0.000312416,1600430195.637,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,p2048,2048,DIVIDE,500000,0.068455697,7303994.000,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,p2048,2048,DIVIDE,500000,0.091980179,5435953.760,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,p2048,2048,DIVIDE,500000,0.067513685,7405905.931,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,p2048,2048,DIVIDE,500000,0.091026917,5492880.750,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,p2048,2048,DIVIDE,500000,0.070670610,7075076.895,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,p2048,2048,DIVIDE,500000,0.094732770,5278004.641,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,p2048,2048,DIVIDE,500000,0.064936588,7699819.400,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,p2048,2048,DIVIDE,500000,0.089456147,5589330.827,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,p2048,2048,ISQRT,1953,0.000715470,2729675.280,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,p2048,2048,ISQRT,1953,0.000512982,3807147.955,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,p2048,2048,ISQRT,500000,1.377843422,362885.936,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,p2048,2048,ISQRT,500000,1.396302269,358088.654,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,p2048,2048,ISQRT,500000,0.466383826,1072078.344,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,p2048,2048,ISQRT,500000,0.484952724,1031028.336,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,p2048,2048,ISQRT,500000,1.374702456,363715.070,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,p2048,2048,ISQRT,500000,1.393803928,358730.514,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,p2048,2048,ISQRT,500000,0.460766456,1085148.438,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,p2048,2048,ISQRT,500000,0.479719977,1042274.710,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-1t,p2048,2048,MODMUL_R2,62500,0.096039876,650771.353,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,gmp-nt,p2048,2048,MODMUL_R2,62500,0.002423254,25791765.774,0
library,AMD EPYC 9V84 96-Core Processor,host-cpu,openssl-nt,p2048,2048,MODMUL_R2,62500,0.005343906,11695565.006,0
library,NVIDIA H100 NVL,gpu,cgbn,p2048,2048,MODMUL_R2,500000,0.001976480,252974985.833,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-opt,p2048,2048,MODMUL_R2,500000,0.005589639,89451214.826,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-opt,p2048,2048,MODMUL_R2,500000,0.024584967,20337631.527,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-o64,p2048,2048,MODMUL_R2,500000,0.004154171,120360956.598,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-o64,p2048,2048,MODMUL_R2,500000,0.023194964,21556403.323,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il,p2048,2048,MODMUL_R2,500000,0.005236824,95477717.304,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il,p2048,2048,MODMUL_R2,500000,0.024339829,20542461.454,0
opencl-kernel,NVIDIA H100 NVL,GPU,w32-il64,p2048,2048,MODMUL_R2,500000,0.003745080,133508495.409,0
opencl-e2e,NVIDIA H100 NVL,GPU,w32-il64,p2048,2048,MODMUL_R2,500000,0.022548139,22174779.057,0
```
