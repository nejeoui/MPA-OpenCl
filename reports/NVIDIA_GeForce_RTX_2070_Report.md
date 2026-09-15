# MPA-OpenCL benchmark report - NVIDIA GeForce RTX 2070

> **Note.** The CGBN reference column has been removed from this report:
> it was invalid (see `reports/README.md`). Every other column is
> unaffected, and every configuration was verified word-for-word against
> GMP before it was timed.


> **Partial report.** The run was interrupted or hit its time budget.
> Rows that never ran are marked `n/a`.

## 1. System under test

2 OpenCL device(s) exercised with the identical kernels and operands.

### Device 0 - NVIDIA GeForce RTX 2070 (GPU)

| Property | Value |
|---|---|
| Model | NVIDIA GeForce RTX 2070 |
| Type | GPU |
| Vendor | NVIDIA Corporation |
| Device memory | 7.60 GiB |
| Max single allocation | 1.90 GiB |
| Local memory | 48 KiB |
| Global cache | 1152 KiB |
| Compute units | 36 |
| Max clock | 1620 MHz |
| Max work-group size | 1024 |
| OpenCL version | OpenCL 3.0 CUDA |
| Driver | 580.159.03 |

### Device 1 - cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz (CPU)

| Property | Value |
|---|---|
| Model | cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz |
| Type | CPU |
| Vendor | GenuineIntel |
| Device memory | 29.27 GiB |
| Max single allocation | 8.00 GiB |
| Local memory | 256 KiB |
| Global cache | 35840 KiB |
| Compute units | 56 |
| Max clock | 3300 MHz |
| Max work-group size | 4096 |
| OpenCL version | OpenCL 3.0 PoCL HSTR: cpu-x86_64-pc-linux-gnu-haswell |
| Driver | 5.0+debian |

### Host

| Property | Value |
|---|---|
| CPU | Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz |
| Logical cores | 56 |
| OpenMP threads used | 56 |
| RAM | 31.3 GB |
| OS | Ubuntu 24.04.4 LTS |
| Kernel | 5.15.0-185-generic |
| Arch | x86_64 |
| GMP | 6.3.0 |
| OpenSSL | OpenSSL 3.0.13 30 Jan 2024 |

## 2. Method

- Workload auto-sized from the device and host: --min-items from 700 x compute units, --items from ten times that capped by host RAM. Either flag, given explicitly, overrides its half.
- Base workload 50000 items, scaled down per operator by its cost weight and by modulus size. Device rows honour --min-items (25200) so the GPU is not left idle; the CPU libraries keep the smaller count because a full-width MODEXP there costs minutes. Both counts appear in every row as dev/cpu, and throughput is per-second so they remain comparable.
- 5 timed repetitions, **minimum** reported. Two untimed warm-up launches precede them.
- `kernel` times `clEnqueueNDRangeKernel` + `clFinish` only. `e2e` adds the host->device operand writes and the device->host result read.
- Every OpenCL device runs the same kernels on the same operands, so GPU and CPU-OpenCL columns are directly comparable.
- CPU library baselines (GMP, OpenSSL) run those same operands, with every temporary - including each thread's GMP context, BN_CTX and Montgomery context - allocated outside the timed region, so the figure is the arithmetic and not marshalling. The generator is reseeded per modulus and operation so every backend sees identical inputs.
- Cost weighting drives the wide cells down to a few hundred items, which is tens of microseconds of work - the same order as the cost of entering an OpenMP region. Each baseline pass is therefore repeated until the timed interval reaches 5 ms and the per-pass time is reported; the multi-threaded loop enters one parallel region per interval and partitions the range itself. Without this the multi-threaded GMP figure came out up to 9x slower than the single-threaded one at 2048 bits.
- OpenSSL rows time the nearest BN primitive, which is not always semantically identical (its Montgomery routine expects Montgomery-domain inputs); they measure comparable work, not identical results. Correctness is judged against GMP only.
- Every device cell is checked word-for-word against GMP before it is timed. A cell that mismatches is reported and excluded from the speedup tables.
- Total wall time 5506.0 s.

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
| [1] CPU | `mpaKernels_8bits.cl` (w8) | 56 | 56 | 0 | 0 |

**All configurations correct** - 541 configurations, 0 problems.

## 4. Throughput per device

Operations per second, higher is better. Kernel-only timings.

### Device 0 - NVIDIA GeForce RTX 2070 (GPU)

#### secp256k1 (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | GMP 56T | OpenSSL 56T |
|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | 467.98 M | 1.09 G | 1.48 G | 1.67 G | 1.54 G | 2.41 G | 2.56 G | 49.31 M | 728.38 M | 774.49 M |
| SUBTRACT | 50000 / 50000 | 463.40 M | 1.10 G | 1.50 G | 1.71 G | 1.55 G | 2.36 G | 2.52 G | 61.02 M | 38.37 M | 72.55 M |
| ADDMOD | 50000 / 50000 | 312.74 M | 761.93 M | 1.17 G | 1.94 G | 2.09 G | 2.54 G | 2.66 G | 15.25 M | 146.31 M | 27.66 M |
| SUBTRACTMOD | 50000 / 50000 | 307.40 M | 755.64 M | 1.20 G | 2.04 G | 1.99 G | 2.77 G | 2.74 G | 21.22 M | 179.04 M | 23.48 M |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | 11.24 M | 56.48 M | 203.91 M | 929.20 M | 941.18 M | 1.93 G | 2.07 G | 39.97 M | 343.16 M | 111.35 M |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | 38.42 M | 173.30 M | 496.69 M | 583.21 M | 548.64 M | 1.05 G | 1.12 G | 39.96 M | 343.64 M | 126.52 M |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | 148.32 M | 804.05 M | 1.82 G | 1.37 G | 1.74 G | 1.68 G | 1.72 G | 5.14 M | 44.54 M | 196.20 M |
| COMPARE | 50000 / 50000 | 492.97 M | 1.18 G | - | 1.90 G | 2.02 G | 4.67 G | 4.47 G | 114.70 M | 1.05 G | 855.98 M |
| REDUCE | 25200 / 6250 | 82.37 M | 175.41 M | - | 369.21 M | 516.87 M | 497.12 M | 538.92 M | 41.47 M | 414.70 M | 27.15 M |
| MODMUL | 25200 / 3125 | 30.73 M | 72.06 M | - | 126.74 M | 198.52 M | 157.92 M | 197.78 M | 8.89 M | 74.88 M | 14.02 M |
| MODEXP | 25200 / 781 | 898.32 k | 4.90 M | - | 5.29 M | 12.97 M | 7.06 M | 12.75 M | 86.19 k | 866.78 k | 441.12 k |
| EXPONENTIATION | 25200 / 781 | 544.80 k | 1.83 M | - | 28.94 M | 46.18 M | 37.46 M | 46.61 M | 244.10 k | 2.29 M | 118.08 k |
| DIVIDE | 25200 / 6250 | 56.06 M | 86.37 M | - | 163.68 M | 218.62 M | 217.96 M | 227.15 M | 19.29 M | 167.11 M | 30.22 M |
| ISQRT | 25200 / 1562 | 4.69 M | 5.90 M | - | 19.48 M | 29.75 M | 25.71 M | 29.57 M | 8.44 M | 94.87 M | n/a |
| MODMUL_R2 | 50000 / 50000 | 164.99 M | 784.87 M | - | 679.75 M | 1.15 G | 831.09 M | 1.22 G | 8.99 M | 77.12 M | 14.81 M |

#### rsa256(composite) (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | GMP 56T | OpenSSL 56T |
|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | 604.44 M | 1.11 G | 1.50 G | 1.41 G | 1.53 G | 2.29 G | 2.32 G | 49.29 M | 451.30 M | 493.81 M |
| SUBTRACT | 50000 / 50000 | 609.90 M | 1.10 G | 1.55 G | 1.41 G | 1.56 G | 2.34 G | 2.34 G | 44.27 M | 528.45 M | 506.38 M |
| ADDMOD | 50000 / 50000 | 439.28 M | 834.29 M | 1.34 G | 1.71 G | 2.16 G | 2.57 G | 2.61 G | 14.48 M | 157.17 M | 28.06 M |
| SUBTRACTMOD | 50000 / 50000 | 406.98 M | 771.11 M | 1.21 G | 1.73 G | 1.62 G | 2.52 G | 2.61 G | 20.90 M | 176.58 M | 24.46 M |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | 14.85 M | 56.55 M | 204.67 M | 745.52 M | 737.12 M | 1.65 G | 1.99 G | 40.10 M | 344.26 M | 188.76 M |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | 50.36 M | 172.93 M | 496.98 M | 464.78 M | 443.66 M | 851.67 M | 933.52 M | 39.89 M | 344.48 M | 181.15 M |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | 196.67 M | 813.27 M | 1.73 G | 1.41 G | 1.42 G | 1.53 G | 1.57 G | 5.18 M | 42.33 M | 187.42 M |
| COMPARE | 50000 / 50000 | 644.34 M | 1.17 G | - | 1.75 G | 1.73 G | 4.18 G | 4.04 G | 116.00 M | 1.05 G | 860.03 M |
| REDUCE | 25200 / 6250 | 107.68 M | 173.98 M | - | 366.82 M | 399.44 M | 392.95 M | 417.28 M | 26.10 M | 231.98 M | 27.77 M |
| MODMUL | 25200 / 3125 | 40.63 M | 71.78 M | - | 124.85 M | 153.89 M | 126.30 M | 154.63 M | 8.74 M | 75.45 M | 15.32 M |
| MODEXP | 25200 / 781 | 901.63 k | 4.91 M | - | 5.25 M | 9.80 M | 5.31 M | 9.64 M | 89.66 k | 890.27 k | 445.85 k |
| EXPONENTIATION | 25200 / 781 | 541.34 k | 1.82 M | - | 29.14 M | 35.68 M | 28.84 M | 35.94 M | 246.37 k | 2.16 M | 124.77 k |
| DIVIDE | 25200 / 6250 | 55.46 M | 85.60 M | - | 159.23 M | 169.74 M | 167.23 M | 168.78 M | 18.95 M | 221.26 M | 32.46 M |
| ISQRT | 25200 / 1562 | 4.68 M | 5.90 M | - | 19.27 M | 22.76 M | 19.54 M | 22.52 M | 8.66 M | 96.61 M | n/a |
| MODMUL_R2 | 50000 / 50000 | 165.47 M | 750.69 M | - | 725.28 M | 1.03 G | 837.68 M | 1.16 G | 8.85 M | 77.25 M | 15.24 M |

#### brainpoolP512r1 (512-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | GMP 56T | OpenSSL 56T |
|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 25200 / 25000 | 222.43 M | 401.84 M | 385.66 M | 512.94 M | 505.98 M | 779.29 M | 783.36 M | 46.05 M | 423.08 M | 460.34 M |
| SUBTRACT | 25200 / 25000 | 223.26 M | 404.70 M | 409.91 M | 513.92 M | 515.11 M | 773.08 M | 776.51 M | 54.78 M | 501.31 M | 481.30 M |
| ADDMOD | 25200 / 25000 | 165.57 M | 319.21 M | 327.51 M | 426.66 M | 426.98 M | 1.00 G | 1.02 G | 17.92 M | 147.19 M | 22.45 M |
| SUBTRACTMOD | 25200 / 25000 | 144.07 M | 298.03 M | 304.54 M | 424.86 M | 418.17 M | 934.20 M | 975.04 M | 13.89 M | 19.67 M | 8.47 M |
| MULTIPLYOPERANDSCANNING | 25200 / 25000 | 2.09 M | 8.06 M | 24.66 M | 178.15 M | 177.53 M | 347.27 M | 367.73 M | 12.70 M | 10.53 M | 22.30 M |
| MULTIPLYPRODUCTSCANNING | 25200 / 25000 | 7.23 M | 27.19 M | 64.70 M | 74.18 M | 73.82 M | 165.59 M | 184.67 M | 17.63 M | 19.80 M | 18.60 M |
| MONTGOMERYMULTIPLICATION | 25200 / 25000 | 51.34 M | 174.45 M | 454.22 M | 331.16 M | 434.11 M | 400.48 M | 557.37 M | 1.69 M | 22.06 M | 91.61 M |
| COMPARE | 25200 / 25000 | 265.37 M | 471.17 M | - | 643.33 M | 654.02 M | 1.91 G | 1.95 G | 54.82 M | 18.61 M | 766.40 M |
| REDUCE | 25200 / 3125 | 36.45 M | 46.73 M | - | 133.42 M | 139.93 M | 135.26 M | 132.95 M | 17.57 M | 192.68 M | 24.58 M |
| MODMUL | 25200 / 1562 | 13.15 M | 18.53 M | - | 36.86 M | 45.17 M | 30.70 M | 38.01 M | 3.29 M | 41.26 M | 9.98 M |
| MODEXP | 25200 / 390 | 62.36 k | 664.68 k | - | 839.20 k | 1.39 M | 784.55 k | 1.42 M | 16.47 k | 155.93 k | 154.15 k |
| EXPONENTIATION | 25200 / 390 | 69.90 k | 277.82 k | - | 898.95 k | 958.36 k | 882.91 k | 942.98 k | 67.56 k | 402.41 k | 33.57 k |
| DIVIDE | 25200 / 3125 | 18.52 M | 20.54 M | - | 64.00 M | 66.87 M | 67.22 M | 69.75 M | 8.88 M | 150.23 M | 28.43 M |
| ISQRT | 25200 / 781 | 786.61 k | 913.46 k | - | 4.46 M | 4.96 M | 3.78 M | 4.10 M | 3.82 M | 54.40 M | n/a |
| MODMUL_R2 | 25200 / 25000 | 35.67 M | 165.02 M | - | 264.20 M | 266.76 M | 282.64 M | 299.24 M | 4.48 M | 43.15 M | 10.43 M |

#### p1024 (1024-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | GMP 56T | OpenSSL 56T |
|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 25200 / 12500 | 54.82 M | 104.68 M | 74.68 M | 199.14 M | 204.67 M | 496.37 M | 509.31 M | 35.38 M | 370.52 M | 177.48 M |
| SUBTRACT | 25200 / 12500 | 54.94 M | 105.33 M | 74.12 M | 200.80 M | 202.15 M | 461.70 M | 501.77 M | 32.42 M | 303.18 M | 281.54 M |
| ADDMOD | 25200 / 12500 | 40.89 M | 78.78 M | 69.21 M | 196.19 M | 195.67 M | 573.44 M | 578.63 M | 8.59 M | 97.65 M | 21.31 M |
| SUBTRACTMOD | 25200 / 12500 | 40.72 M | 78.79 M | 69.41 M | 196.23 M | 193.96 M | 568.37 M | 570.05 M | 16.12 M | 116.47 M | 18.01 M |
| MULTIPLYOPERANDSCANNING | 25200 / 12500 | 400.45 k | 1.60 M | 2.93 M | 83.12 M | 86.09 M | 132.06 M | 136.96 M | 5.23 M | 53.16 M | 36.20 M |
| MULTIPLYPRODUCTSCANNING | 25200 / 12500 | 930.78 k | 3.67 M | 13.67 M | 13.91 M | 14.18 M | 48.82 M | 51.49 M | 5.23 M | 53.20 M | 8.79 M |
| MONTGOMERYMULTIPLICATION | 25200 / 12500 | 9.37 M | 46.46 M | 174.85 M | 106.46 M | 130.46 M | 152.73 M | 213.31 M | 577.27 k | 8.27 M | 30.30 M |
| COMPARE | 25200 / 12500 | 101.19 M | 193.41 M | - | 344.88 M | 343.82 M | 1.83 G | 1.73 G | 66.39 M | 950.50 M | 28.13 M |
| REDUCE | 25200 / 1562 | 6.91 M | 12.74 M | - | 49.90 M | 50.52 M | 56.69 M | 58.48 M | 22.48 M | 352.89 M | 19.81 M |
| MODMUL | 25200 / 781 | 1.66 M | 4.81 M | - | 10.90 M | 14.36 M | 8.59 M | 15.23 M | 1.28 M | 17.09 M | 5.01 M |
| MODEXP | 25200 / 195 | 7.13 k | 50.36 k | - | 97.48 k | 163.51 k | 97.03 k | 163.30 k | 2.59 k | 20.67 k | 29.04 k |
| EXPONENTIATION | 25200 / 195 | 7.48 k | 35.20 k | - | 115.52 k | 117.24 k | 116.22 k | 118.19 k | 15.63 k | 68.07 k | 10.11 k |
| DIVIDE | 25200 / 1562 | 265.50 k | 1.81 M | - | 12.74 M | 13.65 M | 12.83 M | 13.46 M | 11.10 M | 121.06 M | 34.67 M |
| ISQRT | 25200 / 390 | 19.39 k | 113.49 k | - | 673.86 k | 717.19 k | 671.85 k | 722.27 k | 2.09 M | 28.31 M | n/a |
| MODMUL_R2 | 25200 / 12500 | 6.58 M | 40.65 M | - | 70.20 M | 91.89 M | 85.15 M | 116.21 M | 1.27 M | 17.77 M | 5.80 M |

#### p2048 (2048-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | GMP 56T | OpenSSL 56T |
|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 25200 / 6250 | 28.04 M | 56.60 M | 110.82 M | 109.02 M | 109.93 M | 271.63 M | 273.97 M | 19.37 M | 285.11 M | 278.12 M |
| SUBTRACT | 25200 / 6250 | 28.28 M | 56.16 M | 110.81 M | 108.91 M | 108.93 M | 272.08 M | 273.43 M | 28.30 M | 224.05 M | 263.43 M |
| ADDMOD | 25200 / 6250 | 16.87 M | 43.66 M | 89.44 M | 100.66 M | 104.90 M | 283.58 M | 293.32 M | 9.85 M | 88.10 M | 22.69 M |
| SUBTRACTMOD | 25200 / 6250 | 16.18 M | 31.14 M | 83.55 M | 103.63 M | 102.51 M | 271.15 M | 282.75 M | 12.21 M | 101.38 M | 16.71 M |
| MULTIPLYOPERANDSCANNING | 25200 / 6250 | 91.91 k | 342.66 k | 1.58 M | 21.39 M | 26.40 M | 29.93 M | 29.36 M | 1.58 M | 16.41 M | 13.00 M |
| MULTIPLYPRODUCTSCANNING | 25200 / 6250 | 234.43 k | 928.00 k | 3.68 M | 2.77 M | 2.76 M | 7.54 M | 7.81 M | 1.16 M | 15.82 M | 12.32 M |
| MONTGOMERYMULTIPLICATION | 25200 / 6250 | 275.30 k | 11.11 M | 61.27 M | 36.18 M | 44.63 M | 33.60 M | 41.92 M | 180.40 k | 1.64 M | 1.69 M |
| COMPARE | 25200 / 6250 | 51.61 M | 101.82 M | - | 197.89 M | 191.99 M | 555.31 M | 569.92 M | 57.51 M | 36.89 M | 40.35 M |
| REDUCE | 25200 / 781 | 32.17 k | 2.37 M | - | 14.16 M | 15.72 M | 13.02 M | 12.52 M | 14.77 M | 41.51 M | 19.38 M |
| MODMUL | 25200 / 390 | 19.35 k | 1.05 M | - | 2.70 M | 3.56 M | 2.13 M | 3.47 M | 433.65 k | 5.83 M | 1.83 M |
| MODEXP | 25200 / 97 | 151.0 | 1.03 k | - | 10.31 k | 5.00 k | 10.56 k | 5.06 k | 370.9 | 3.17 k | 7.48 k |
| EXPONENTIATION | 25200 / 97 | 627.6 | 3.18 k | - | 12.29 k | 12.52 k | 12.44 k | 12.76 k | 2.27 k | 28.10 k | 2.44 k |
| DIVIDE | 25200 / 781 | 13.46 k | 65.11 k | - | 884.16 k | 998.97 k | 876.73 k | 1.01 M | 5.48 M | 94.30 M | 19.77 M |
| ISQRT | 25200 / 195 | 1.03 k | 2.84 k | - | 274.55 k | 321.56 k | 274.69 k | 320.63 k | 1.11 M | 14.32 M | n/a |
| MODMUL_R2 | 25200 / 6250 | 378.91 k | 10.61 M | - | 20.68 M | 27.69 M | 23.44 M | 30.74 M | 423.96 k | 6.03 M | 1.95 M |

### Device 1 - cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz (CPU)

#### secp256k1 (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | GMP 56T | OpenSSL 56T |
|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | 127.93 M | - | - | - | - | - | - | 49.31 M | 728.38 M | 774.49 M |
| SUBTRACT | 50000 / 50000 | 191.23 M | - | - | - | - | - | - | 61.02 M | 38.37 M | 72.55 M |
| ADDMOD | 50000 / 50000 | 186.07 M | - | - | - | - | - | - | 15.25 M | 146.31 M | 27.66 M |
| SUBTRACTMOD | 50000 / 50000 | 123.41 M | - | - | - | - | - | - | 21.22 M | 179.04 M | 23.48 M |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | 6.84 M | - | - | - | - | - | - | 39.97 M | 343.16 M | 111.35 M |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | 18.16 M | - | - | - | - | - | - | 39.96 M | 343.64 M | 126.52 M |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | 655.28 k | - | - | - | - | - | - | 5.14 M | 44.54 M | 196.20 M |
| COMPARE | 50000 / 50000 | 163.59 M | - | - | - | - | - | - | 114.70 M | 1.05 G | 855.98 M |
| REDUCE | 25200 / 6250 | 2.80 M | - | - | - | - | - | - | 41.47 M | 414.70 M | 27.15 M |
| MODMUL | 25200 / 3125 | 711.84 k | - | - | - | - | - | - | 8.89 M | 74.88 M | 14.02 M |
| MODEXP | 25200 / 781 | 8.17 k | - | - | - | - | - | - | 86.19 k | 866.78 k | 441.12 k |
| EXPONENTIATION | 25200 / 781 | 18.70 k | - | - | - | - | - | - | 244.10 k | 2.29 M | 118.08 k |
| DIVIDE | 25200 / 6250 | 1.45 M | - | - | - | - | - | - | 19.29 M | 167.11 M | 30.22 M |
| ISQRT | 25200 / 1562 | 98.34 k | - | - | - | - | - | - | 8.44 M | 94.87 M | n/a |
| MODMUL_R2 | 50000 / 50000 | 3.03 M | - | - | - | - | - | - | 8.99 M | 77.12 M | 14.81 M |

#### rsa256(composite) (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | GMP 56T | OpenSSL 56T |
|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | 139.86 M | - | - | - | - | - | - | 49.29 M | 451.30 M | 493.81 M |
| SUBTRACT | 50000 / 50000 | 126.20 M | - | - | - | - | - | - | 44.27 M | 528.45 M | 506.38 M |
| ADDMOD | 50000 / 50000 | 106.33 M | - | - | - | - | - | - | 14.48 M | 157.17 M | 28.06 M |
| SUBTRACTMOD | 50000 / 50000 | 172.64 M | - | - | - | - | - | - | 20.90 M | 176.58 M | 24.46 M |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | 13.30 M | - | - | - | - | - | - | 40.10 M | 344.26 M | 188.76 M |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | 10.15 M | - | - | - | - | - | - | 39.89 M | 344.48 M | 181.15 M |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | 3.22 M | - | - | - | - | - | - | 5.18 M | 42.33 M | 187.42 M |
| COMPARE | 50000 / 50000 | 167.76 M | - | - | - | - | - | - | 116.00 M | 1.05 G | 860.03 M |
| REDUCE | 25200 / 6250 | 2.49 M | - | - | - | - | - | - | 26.10 M | 231.98 M | 27.77 M |
| MODMUL | 25200 / 3125 | 928.31 k | - | - | - | - | - | - | 8.74 M | 75.45 M | 15.32 M |
| MODEXP | 25200 / 781 | 8.35 k | - | - | - | - | - | - | 89.66 k | 890.27 k | 445.85 k |
| EXPONENTIATION | 25200 / 781 | 18.45 k | - | - | - | - | - | - | 246.37 k | 2.16 M | 124.77 k |
| DIVIDE | 25200 / 6250 | 1.68 M | - | - | - | - | - | - | 18.95 M | 221.26 M | 32.46 M |
| ISQRT | 25200 / 1562 | 106.96 k | - | - | - | - | - | - | 8.66 M | 96.61 M | n/a |
| MODMUL_R2 | 50000 / 50000 | 3.46 M | - | - | - | - | - | - | 8.85 M | 77.25 M | 15.24 M |

#### brainpoolP512r1 (512-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | GMP 56T | OpenSSL 56T |
|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 25200 / 25000 | 70.37 M | - | - | - | - | - | - | 46.05 M | 423.08 M | 460.34 M |
| SUBTRACT | 25200 / 25000 | 62.88 M | - | - | - | - | - | - | 54.78 M | 501.31 M | 481.30 M |
| ADDMOD | 25200 / 25000 | 58.24 M | - | - | - | - | - | - | 17.92 M | 147.19 M | 22.45 M |
| SUBTRACTMOD | 25200 / 25000 | 54.52 M | - | - | - | - | - | - | 13.89 M | 19.67 M | 8.47 M |
| MULTIPLYOPERANDSCANNING | 25200 / 25000 | 1.45 M | - | - | - | - | - | - | 12.70 M | 10.53 M | 22.30 M |
| MULTIPLYPRODUCTSCANNING | 25200 / 25000 | 2.46 M | - | - | - | - | - | - | 17.63 M | 19.80 M | 18.60 M |
| MONTGOMERYMULTIPLICATION | 25200 / 25000 | 895.60 k | - | - | - | - | - | - | 1.69 M | 22.06 M | 91.61 M |
| COMPARE | 25200 / 25000 | 111.33 M | - | - | - | - | - | - | 54.82 M | 18.61 M | 766.40 M |
| REDUCE | 25200 / 3125 | 791.45 k | - | - | - | - | - | - | 17.57 M | 192.68 M | 24.58 M |
| MODMUL | 25200 / 1562 | 144.16 k | - | - | - | - | - | - | 3.29 M | 41.26 M | 9.98 M |
| MODEXP | 25200 / 390 | over budget | - | - | - | - | - | - | 16.47 k | 155.93 k | 154.15 k |
| EXPONENTIATION | 25200 / 390 | 2.41 k | - | - | - | - | - | - | 67.56 k | 402.41 k | 33.57 k |
| DIVIDE | 25200 / 3125 | 287.28 k | - | - | - | - | - | - | 8.88 M | 150.23 M | 28.43 M |
| ISQRT | 25200 / 781 | 24.58 k | - | - | - | - | - | - | 3.82 M | 54.40 M | n/a |
| MODMUL_R2 | 25200 / 25000 | 257.04 k | - | - | - | - | - | - | 4.48 M | 43.15 M | 10.43 M |

#### p1024 (1024-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | GMP 56T | OpenSSL 56T |
|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 25200 / 12500 | 50.62 M | - | - | - | - | - | - | 35.38 M | 370.52 M | 177.48 M |
| SUBTRACT | 25200 / 12500 | 44.05 M | - | - | - | - | - | - | 32.42 M | 303.18 M | 281.54 M |
| ADDMOD | 25200 / 12500 | 34.38 M | - | - | - | - | - | - | 8.59 M | 97.65 M | 21.31 M |
| SUBTRACTMOD | 25200 / 12500 | 36.98 M | - | - | - | - | - | - | 16.12 M | 116.47 M | 18.01 M |
| MULTIPLYOPERANDSCANNING | 25200 / 12500 | 249.14 k | - | - | - | - | - | - | 5.23 M | 53.16 M | 36.20 M |
| MULTIPLYPRODUCTSCANNING | 25200 / 12500 | 690.38 k | - | - | - | - | - | - | 5.23 M | 53.20 M | 8.79 M |
| MONTGOMERYMULTIPLICATION | 25200 / 12500 | 127.21 k | - | - | - | - | - | - | 577.27 k | 8.27 M | 30.30 M |
| COMPARE | 25200 / 12500 | 68.99 M | - | - | - | - | - | - | 66.39 M | 950.50 M | 28.13 M |
| REDUCE | 25200 / 1562 | 143.08 k | - | - | - | - | - | - | 22.48 M | 352.89 M | 19.81 M |
| MODMUL | 25200 / 781 | 39.59 k | - | - | - | - | - | - | 1.28 M | 17.09 M | 5.01 M |
| MODEXP | 25200 / 195 | over budget | - | - | - | - | - | - | 2.59 k | 20.67 k | 29.04 k |
| EXPONENTIATION | 25200 / 195 | - | - | - | - | - | - | - | 15.63 k | 68.07 k | 10.11 k |
| DIVIDE | 25200 / 1562 | - | - | - | - | - | - | - | 11.10 M | 121.06 M | 34.67 M |
| ISQRT | 25200 / 390 | - | - | - | - | - | - | - | 2.09 M | 28.31 M | n/a |
| MODMUL_R2 | 25200 / 12500 | - | - | - | - | - | - | - | 1.27 M | 17.77 M | 5.80 M |

#### p2048 (2048-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | GMP 56T | OpenSSL 56T |
|---|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 25200 / 6250 | - | - | - | - | - | - | - | 19.37 M | 285.11 M | 278.12 M |
| SUBTRACT | 25200 / 6250 | - | - | - | - | - | - | - | 28.30 M | 224.05 M | 263.43 M |
| ADDMOD | 25200 / 6250 | - | - | - | - | - | - | - | 9.85 M | 88.10 M | 22.69 M |
| SUBTRACTMOD | 25200 / 6250 | - | - | - | - | - | - | - | 12.21 M | 101.38 M | 16.71 M |
| MULTIPLYOPERANDSCANNING | 25200 / 6250 | - | - | - | - | - | - | - | 1.58 M | 16.41 M | 13.00 M |
| MULTIPLYPRODUCTSCANNING | 25200 / 6250 | - | - | - | - | - | - | - | 1.16 M | 15.82 M | 12.32 M |
| MONTGOMERYMULTIPLICATION | 25200 / 6250 | - | - | - | - | - | - | - | 180.40 k | 1.64 M | 1.69 M |
| COMPARE | 25200 / 6250 | - | - | - | - | - | - | - | 57.51 M | 36.89 M | 40.35 M |
| REDUCE | 25200 / 781 | - | - | - | - | - | - | - | 14.77 M | 41.51 M | 19.38 M |
| MODMUL | 25200 / 390 | - | - | - | - | - | - | - | 433.65 k | 5.83 M | 1.83 M |
| MODEXP | 25200 / 97 | - | - | - | - | - | - | - | 370.9 | 3.17 k | 7.48 k |
| EXPONENTIATION | 25200 / 97 | - | - | - | - | - | - | - | 2.27 k | 28.10 k | 2.44 k |
| DIVIDE | 25200 / 781 | - | - | - | - | - | - | - | 5.48 M | 94.30 M | 19.77 M |
| ISQRT | 25200 / 195 | - | - | - | - | - | - | - | 1.11 M | 14.32 M | n/a |
| MODMUL_R2 | 25200 / 6250 | - | - | - | - | - | - | - | 423.96 k | 6.03 M | 1.95 M |

## 5. Head to head

Best OpenCL GPU result against best OpenCL CPU result and the CPU libraries.
Ratios above 1.00x mean the GPU is faster than that baseline.

### secp256k1 (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GMP 56T | OpenSSL | GPU vs CPU-CL | GPU vs GMP 56T | GPU vs OpenSSL |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | w32-il64 | 2.56 G | w8 | 127.93 M | 49.31 M | 728.38 M | 774.49 M | 19.99x | 3.51x | 3.30x |
| SUBTRACT | w32-il64 | 2.52 G | w8 | 191.23 M | 61.02 M | 38.37 M | 72.55 M | 13.19x | 65.73x | 34.77x |
| ADDMOD | w32-il64 | 2.66 G | w8 | 186.07 M | 15.25 M | 146.31 M | 27.66 M | 14.29x | 18.18x | 96.15x |
| SUBTRACTMOD | w32-il | 2.77 G | w8 | 123.41 M | 21.22 M | 179.04 M | 23.48 M | 22.42x | 15.46x | 117.85x |
| MULTIPLYOPERANDSCANNING | w32-il64 | 2.07 G | w8 | 6.84 M | 39.97 M | 343.16 M | 111.35 M | 302.27x | 6.03x | 18.58x |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 1.12 G | w8 | 18.16 M | 39.96 M | 343.64 M | 126.52 M | 61.72x | 3.26x | 8.86x |
| MONTGOMERYMULTIPLICATION | w32 | 1.82 G | w8 | 655.28 k | 5.14 M | 44.54 M | 196.20 M | 2784.30x | 40.97x | 9.30x |
| COMPARE | w32-il | 4.67 G | w8 | 163.59 M | 114.70 M | 1.05 G | 855.98 M | 28.55x | 4.46x | 5.46x |
| REDUCE | w32-il64 | 133.66 M | w8 | 694.38 k | 41.47 M | 414.70 M | 27.15 M | 192.49x | 0.32x | 4.92x |
| MODMUL | w32-o64 | 24.62 M | w8 | 88.27 k | 8.89 M | 74.88 M | 14.02 M | 278.88x | 0.33x | 1.76x |
| MODEXP | w32-o64 | 402.07 k | w8 | 253.3 | 86.19 k | 866.78 k | 441.12 k | 1587.36x | 0.46x | 0.91x |
| EXPONENTIATION | w32-il64 | 1.44 M | w8 | 579.4 | 244.10 k | 2.29 M | 118.08 k | 2493.27x | 0.63x | 12.23x |
| DIVIDE | w32-il64 | 56.34 M | w8 | 360.70 k | 19.29 M | 167.11 M | 30.22 M | 156.19x | 0.34x | 1.86x |
| ISQRT | w32-o64 | 1.84 M | w8 | 6.10 k | 8.44 M | 94.87 M | n/a | 302.56x | 0.02x | n/a |
| MODMUL_R2 | w32-il64 | 1.22 G | w8 | 3.03 M | 8.99 M | 77.12 M | 14.81 M | 403.99x | 15.85x | 82.55x |

### rsa256(composite) (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GMP 56T | OpenSSL | GPU vs CPU-CL | GPU vs GMP 56T | GPU vs OpenSSL |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | w32-il64 | 2.32 G | w8 | 139.86 M | 49.29 M | 451.30 M | 493.81 M | 16.60x | 5.14x | 4.70x |
| SUBTRACT | w32-il | 2.34 G | w8 | 126.20 M | 44.27 M | 528.45 M | 506.38 M | 18.57x | 4.43x | 4.63x |
| ADDMOD | w32-il64 | 2.61 G | w8 | 106.33 M | 14.48 M | 157.17 M | 28.06 M | 24.54x | 16.60x | 93.01x |
| SUBTRACTMOD | w32-il64 | 2.61 G | w8 | 172.64 M | 20.90 M | 176.58 M | 24.46 M | 15.14x | 14.80x | 106.88x |
| MULTIPLYOPERANDSCANNING | w32-il64 | 1.99 G | w8 | 13.30 M | 40.10 M | 344.26 M | 188.76 M | 149.40x | 5.77x | 10.52x |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 933.52 M | w8 | 10.15 M | 39.89 M | 344.48 M | 181.15 M | 91.98x | 2.71x | 5.15x |
| MONTGOMERYMULTIPLICATION | w32 | 1.73 G | w8 | 3.22 M | 5.18 M | 42.33 M | 187.42 M | 538.23x | 40.94x | 9.25x |
| COMPARE | w32-il | 4.18 G | w8 | 167.76 M | 116.00 M | 1.05 G | 860.03 M | 24.93x | 4.00x | 4.86x |
| REDUCE | w32-il64 | 103.49 M | w8 | 616.84 k | 26.10 M | 231.98 M | 27.77 M | 167.78x | 0.45x | 3.73x |
| MODMUL | w32-il64 | 19.18 M | w8 | 115.12 k | 8.74 M | 75.45 M | 15.32 M | 166.57x | 0.25x | 1.25x |
| MODEXP | w32-o64 | 303.63 k | w8 | 258.9 | 89.66 k | 890.27 k | 445.85 k | 1172.71x | 0.34x | 0.68x |
| EXPONENTIATION | w32-il64 | 1.11 M | w8 | 571.9 | 246.37 k | 2.16 M | 124.77 k | 1947.56x | 0.52x | 8.93x |
| DIVIDE | w32-o64 | 42.10 M | w8 | 417.75 k | 18.95 M | 221.26 M | 32.46 M | 100.77x | 0.19x | 1.30x |
| ISQRT | w32-o64 | 1.41 M | w8 | 6.63 k | 8.66 M | 96.61 M | n/a | 212.79x | 0.01x | n/a |
| MODMUL_R2 | w32-il64 | 1.16 G | w8 | 3.46 M | 8.85 M | 77.25 M | 15.24 M | 334.20x | 14.96x | 75.79x |

### brainpoolP512r1 (512-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GMP 56T | OpenSSL | GPU vs CPU-CL | GPU vs GMP 56T | GPU vs OpenSSL |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | w32-il64 | 777.15 M | w8 | 69.81 M | 46.05 M | 423.08 M | 460.34 M | 11.13x | 1.84x | 1.69x |
| SUBTRACT | w32-il64 | 770.35 M | w8 | 62.38 M | 54.78 M | 501.31 M | 481.30 M | 12.35x | 1.54x | 1.60x |
| ADDMOD | w32-il64 | 1.01 G | w8 | 57.78 M | 17.92 M | 147.19 M | 22.45 M | 17.47x | 6.86x | 44.97x |
| SUBTRACTMOD | w32-il64 | 967.31 M | w8 | 54.09 M | 13.89 M | 19.67 M | 8.47 M | 17.88x | 49.17x | 114.16x |
| MULTIPLYOPERANDSCANNING | w32-il64 | 364.81 M | w8 | 1.44 M | 12.70 M | 10.53 M | 22.30 M | 253.66x | 34.65x | 16.36x |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 183.20 M | w8 | 2.44 M | 17.63 M | 19.80 M | 18.60 M | 75.22x | 9.25x | 9.85x |
| MONTGOMERYMULTIPLICATION | w32-il64 | 552.95 M | w8 | 888.50 k | 1.69 M | 22.06 M | 91.61 M | 622.34x | 25.07x | 6.04x |
| COMPARE | w32-il64 | 1.93 G | w8 | 110.45 M | 54.82 M | 18.61 M | 766.40 M | 17.52x | 103.99x | 2.52x |
| REDUCE | w32-o64 | 17.35 M | w8 | 98.15 k | 17.57 M | 192.68 M | 24.58 M | 176.80x | 0.09x | 0.71x |
| MODMUL | w32-o64 | 2.80 M | w8 | 8.94 k | 3.29 M | 41.26 M | 9.98 M | 313.32x | 0.07x | 0.28x |
| MODEXP | w32-il64 | 21.94 k | none | n/a | 16.47 k | 155.93 k | 154.15 k | n/a | 0.14x | 0.14x |
| EXPONENTIATION | w32-o64 | 14.83 k | w8 | 37.3 | 67.56 k | 402.41 k | 33.57 k | 397.49x | 0.04x | 0.44x |
| DIVIDE | w32-il64 | 8.65 M | w8 | 35.62 k | 8.88 M | 150.23 M | 28.43 M | 242.80x | 0.06x | 0.30x |
| ISQRT | w32-o64 | 153.83 k | w8 | 761.9 | 3.82 M | 54.40 M | n/a | 201.90x | 0.00x | n/a |
| MODMUL_R2 | w32-il64 | 296.86 M | w8 | 255.00 k | 4.48 M | 43.15 M | 10.43 M | 1164.16x | 6.88x | 28.47x |

### p1024 (1024-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GMP 56T | OpenSSL | GPU vs CPU-CL | GPU vs GMP 56T | GPU vs OpenSSL |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | w32-il64 | 252.63 M | w8 | 25.11 M | 35.38 M | 370.52 M | 177.48 M | 10.06x | 0.68x | 1.42x |
| SUBTRACT | w32-il64 | 248.89 M | w8 | 21.85 M | 32.42 M | 303.18 M | 281.54 M | 11.39x | 0.82x | 0.88x |
| ADDMOD | w32-il64 | 287.02 M | w8 | 17.06 M | 8.59 M | 97.65 M | 21.31 M | 16.83x | 2.94x | 13.47x |
| SUBTRACTMOD | w32-il64 | 282.76 M | w8 | 18.34 M | 16.12 M | 116.47 M | 18.01 M | 15.41x | 2.43x | 15.70x |
| MULTIPLYOPERANDSCANNING | w32-il64 | 67.94 M | w8 | 123.58 k | 5.23 M | 53.16 M | 36.20 M | 549.73x | 1.28x | 1.88x |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 25.54 M | w8 | 342.45 k | 5.23 M | 53.20 M | 8.79 M | 74.58x | 0.48x | 2.91x |
| MONTGOMERYMULTIPLICATION | w32-il64 | 105.81 M | w8 | 63.10 k | 577.27 k | 8.27 M | 30.30 M | 1676.78x | 12.80x | 3.49x |
| COMPARE | w32-il | 905.73 M | w8 | 34.22 M | 66.39 M | 950.50 M | 28.13 M | 26.47x | 0.95x | 32.20x |
| REDUCE | w32-il64 | 3.62 M | w8 | 8.87 k | 22.48 M | 352.89 M | 19.81 M | 408.72x | 0.01x | 0.18x |
| MODMUL | w32-il64 | 472.15 k | w8 | 1.23 k | 1.28 M | 17.09 M | 5.01 M | 384.77x | 0.03x | 0.09x |
| MODEXP | w32-o64 | 1.27 k | none | n/a | 2.59 k | 20.67 k | 29.04 k | n/a | 0.06x | 0.04x |
| EXPONENTIATION | w32-il64 | 914.6 | none | n/a | 15.63 k | 68.07 k | 10.11 k | n/a | 0.01x | 0.09x |
| DIVIDE | w32-o64 | 846.31 k | none | n/a | 11.10 M | 121.06 M | 34.67 M | n/a | 0.01x | 0.02x |
| ISQRT | w32-il64 | 11.18 k | none | n/a | 2.09 M | 28.31 M | n/a | n/a | 0.00x | n/a |
| MODMUL_R2 | w32-il64 | 57.65 M | none | n/a | 1.27 M | 17.77 M | 5.80 M | n/a | 3.24x | 9.94x |

### p2048 (2048-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GMP 56T | OpenSSL | GPU vs CPU-CL | GPU vs GMP 56T | GPU vs OpenSSL |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | w32-il64 | 67.95 M | none | n/a | 19.37 M | 285.11 M | 278.12 M | n/a | 0.24x | 0.24x |
| SUBTRACT | w32-il64 | 67.82 M | none | n/a | 28.30 M | 224.05 M | 263.43 M | n/a | 0.30x | 0.26x |
| ADDMOD | w32-il64 | 72.75 M | none | n/a | 9.85 M | 88.10 M | 22.69 M | n/a | 0.83x | 3.21x |
| SUBTRACTMOD | w32-il64 | 70.13 M | none | n/a | 12.21 M | 101.38 M | 16.71 M | n/a | 0.69x | 4.20x |
| MULTIPLYOPERANDSCANNING | w32-il | 7.42 M | none | n/a | 1.58 M | 16.41 M | 13.00 M | n/a | 0.45x | 0.57x |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 1.94 M | none | n/a | 1.16 M | 15.82 M | 12.32 M | n/a | 0.12x | 0.16x |
| MONTGOMERYMULTIPLICATION | w32 | 15.20 M | none | n/a | 180.40 k | 1.64 M | 1.69 M | n/a | 9.28x | 8.98x |
| COMPARE | w32-il64 | 141.35 M | none | n/a | 57.51 M | 36.89 M | 40.35 M | n/a | 3.83x | 3.50x |
| REDUCE | w32-o64 | 487.21 k | none | n/a | 14.77 M | 41.51 M | 19.38 M | n/a | 0.01x | 0.03x |
| MODMUL | w32-o64 | 55.07 k | none | n/a | 433.65 k | 5.83 M | 1.83 M | n/a | 0.01x | 0.03x |
| MODEXP | w32-il | 40.6 | none | n/a | 370.9 | 3.17 k | 7.48 k | n/a | 0.01x | 0.01x |
| EXPONENTIATION | w32-il64 | 49.1 | none | n/a | 2.27 k | 28.10 k | 2.44 k | n/a | 0.00x | 0.02x |
| DIVIDE | w32-il64 | 31.26 k | none | n/a | 5.48 M | 94.30 M | 19.77 M | n/a | 0.00x | 0.00x |
| ISQRT | w32-o64 | 2.49 k | none | n/a | 1.11 M | 14.32 M | n/a | n/a | 0.00x | n/a |
| MODMUL_R2 | w32-il64 | 7.62 M | none | n/a | 423.96 k | 6.03 M | 1.95 M | n/a | 1.26x | 3.91x |

## 6. Raw data

Also written to `NVIDIA_GeForce_RTX_2070_Report.csv` for analysis.

```csv
kind,device,device_type,kernel,modulus,bits,operation,items,seconds,ops_per_sec,mismatches
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,secp256k1,256,ADD,50000,0.001013913,49313919.702,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,secp256k1,256,ADD,50000,0.000068646,728379840.005,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,secp256k1,256,ADD,50000,0.000064558,774491411.961,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,secp256k1,256,ADD,50000,0.000036480,1370614035.088,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,secp256k1,256,ADD,50000,0.000106843,467976316.187,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,secp256k1,256,ADD,50000,0.001102500,45351473.385,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,secp256k1,256,ADD,50000,0.000046062,1085493649.624,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,secp256k1,256,ADD,50000,0.001047012,47754945.062,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32,secp256k1,256,ADD,50000,0.000033718,1482887041.126,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32,secp256k1,256,ADD,50000,0.001051343,47558218.939,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,secp256k1,256,ADD,50000,0.000029942,1669895595.155,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,secp256k1,256,ADD,50000,0.001105003,45248745.902,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,secp256k1,256,ADD,50000,0.000032487,1539076921.644,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,secp256k1,256,ADD,50000,0.001029464,48568963.631,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,secp256k1,256,ADD,50000,0.000020753,2409289041.718,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,secp256k1,256,ADD,50000,0.001234443,40504097.332,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,secp256k1,256,ADD,50000,0.000019551,2557414875.098,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,secp256k1,256,ADD,50000,0.001091567,45805709.323,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,secp256k1,256,ADD,50000,0.000390838,127930246.239,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,secp256k1,256,ADD,50000,0.001701321,29388928.001,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,secp256k1,256,SUBTRACT,50000,0.000819375,61022144.994,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,secp256k1,256,SUBTRACT,50000,0.001302937,38374850.924,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,secp256k1,256,SUBTRACT,50000,0.000689209,72546899.277,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,secp256k1,256,SUBTRACT,50000,0.000034816,1436121323.529,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,secp256k1,256,SUBTRACT,50000,0.000107897,463404875.283,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,secp256k1,256,SUBTRACT,50000,0.001124965,44445827.288,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,secp256k1,256,SUBTRACT,50000,0.000045417,1100909518.646,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,secp256k1,256,SUBTRACT,50000,0.001049417,47645501.634,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32,secp256k1,256,SUBTRACT,50000,0.000033243,1504075759.770,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32,secp256k1,256,SUBTRACT,50000,0.001023704,48842242.914,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,secp256k1,256,SUBTRACT,50000,0.000029255,1709109802.088,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,secp256k1,256,SUBTRACT,50000,0.001045785,47810974.517,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,secp256k1,256,SUBTRACT,50000,0.000032271,1549378862.451,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,secp256k1,256,SUBTRACT,50000,0.001030389,48525362.007,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,secp256k1,256,SUBTRACT,50000,0.000021186,2360047199.105,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,secp256k1,256,SUBTRACT,50000,0.001292708,38678494.776,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,secp256k1,256,SUBTRACT,50000,0.000019822,2522448817.832,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,secp256k1,256,SUBTRACT,50000,0.001065101,46943904.460,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,secp256k1,256,SUBTRACT,50000,0.000261468,191227991.183,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,secp256k1,256,SUBTRACT,50000,0.001249585,40013284.102,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,secp256k1,256,ADDMOD,50000,0.003279674,15245417.762,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,secp256k1,256,ADDMOD,50000,0.000341735,146312190.896,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,secp256k1,256,ADDMOD,50000,0.001807582,27661262.660,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,secp256k1,256,ADDMOD,50000,0.000034816,1436121323.529,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,secp256k1,256,ADDMOD,50000,0.000159875,312744341.278,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,secp256k1,256,ADDMOD,50000,0.001190503,41999053.777,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,secp256k1,256,ADDMOD,50000,0.000065623,761928020.050,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,secp256k1,256,ADDMOD,50000,0.001108229,45117029.644,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32,secp256k1,256,ADDMOD,50000,0.000042587,1174067346.328,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32,secp256k1,256,ADDMOD,50000,0.001091278,45817839.960,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,secp256k1,256,ADDMOD,50000,0.000025752,1941596730.679,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,secp256k1,256,ADDMOD,50000,0.001053733,47450350.575,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,secp256k1,256,ADDMOD,50000,0.000023935,2088991875.486,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,secp256k1,256,ADDMOD,50000,0.001016516,49187617.041,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,secp256k1,256,ADDMOD,50000,0.000019721,2535366145.373,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,secp256k1,256,ADDMOD,50000,0.001039738,48089037.346,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,secp256k1,256,ADDMOD,50000,0.000018800,2659574802.891,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,secp256k1,256,ADDMOD,50000,0.001112108,44959662.093,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,secp256k1,256,ADDMOD,50000,0.000268710,186074192.223,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,secp256k1,256,ADDMOD,50000,0.001181314,42325749.093,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,secp256k1,256,SUBTRACTMOD,50000,0.002356304,21219668.376,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,secp256k1,256,SUBTRACTMOD,50000,0.000279269,179038852.967,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,secp256k1,256,SUBTRACTMOD,50000,0.002129442,23480335.041,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,secp256k1,256,SUBTRACTMOD,50000,0.000034816,1436121323.529,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.000162653,307402904.926,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.001181457,42320626.573,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,secp256k1,256,SUBTRACTMOD,50000,0.000066169,755641034.524,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,secp256k1,256,SUBTRACTMOD,50000,0.001083512,46146235.375,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32,secp256k1,256,SUBTRACTMOD,50000,0.000041533,1203862000.767,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32,secp256k1,256,SUBTRACTMOD,50000,0.001051825,47536424.534,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,secp256k1,256,SUBTRACTMOD,50000,0.000024566,2035334189.965,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,secp256k1,256,SUBTRACTMOD,50000,0.001033029,48401350.616,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,secp256k1,256,SUBTRACTMOD,50000,0.000025099,1992111408.809,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,secp256k1,256,SUBTRACTMOD,50000,0.001015605,49231738.838,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,secp256k1,256,SUBTRACTMOD,50000,0.000018069,2767170793.401,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,secp256k1,256,SUBTRACTMOD,50000,0.000878931,56887285.443,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,secp256k1,256,SUBTRACTMOD,50000,0.000018269,2736876939.008,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,secp256k1,256,SUBTRACTMOD,50000,0.001165335,42906117.887,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.000405149,123411383.731,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.001690525,29576610.666,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001250973,39968904.012,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000145705,343159135.149,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000449022,111353118.133,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.004447108,11243261.867,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.005669394,8819284.767,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000885237,56482048.908,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.002122248,23559923.262,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000245209,203907688.981,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001315194,38017204.850,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000053810,929195578.986,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001220970,40951047.855,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000053125,941176522.562,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001259351,39702989.757,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000025898,1930652109.574,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001118544,44700968.626,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000024168,2068851764.559,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001352174,36977490.274,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.007305153,6844483.635,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.007153081,6989994.932,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001251332,39957437.293,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000145499,343644978.510,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000395185,126523018.133,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000034816,1436121323.529,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001301411,38419838.563,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.002559726,19533340.601,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000288513,173302419.605,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001516453,32971678.195,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000100667,496687065.800,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001281323,39022166.972,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000085732,583212820.460,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001302911,38375605.831,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000091135,548636770.440,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001334466,37468170.479,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000047429,1054207361.590,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001162608,43006757.980,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000044596,1121176708.323,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001375307,36355518.934,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.002752627,18164466.267,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.004122478,12128627.534,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.009729830,5138835.917,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.001122704,44535335.776,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000254847,196196154.526,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000033728,1482447817.837,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000337117,148316458.832,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.001349644,37046806.201,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000062185,804052172.244,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.001059548,47189933.909,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000027405,1824483866.151,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.001035223,48298771.571,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000036462,1371291512.934,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.001054073,47435044.168,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000028777,1737497275.304,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.001024206,48818304.401,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000029707,1683104576.061,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000897422,55715148.730,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000029078,1719513846.696,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.001393136,35890250.089,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.076303673,655276.451,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.021880536,2285135.978,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,secp256k1,256,COMPARE,50000,0.000435933,114696436.350,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,secp256k1,256,COMPARE,50000,0.000047748,1047164143.559,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,secp256k1,256,COMPARE,50000,0.000058413,855978702.211,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,secp256k1,256,COMPARE,50000,0.000036000,1388888888.889,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,secp256k1,256,COMPARE,50000,0.000101426,492970219.379,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,secp256k1,256,COMPARE,50000,0.001100436,45436536.125,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,secp256k1,256,COMPARE,50000,0.000042253,1183347948.549,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,secp256k1,256,COMPARE,50000,0.001045851,47807956.991,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,secp256k1,256,COMPARE,50000,0.000026277,1902804407.717,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,secp256k1,256,COMPARE,50000,0.001076027,46467235.735,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,secp256k1,256,COMPARE,50000,0.000024808,2015477434.106,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,secp256k1,256,COMPARE,50000,0.001028122,48632360.655,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,secp256k1,256,COMPARE,50000,0.000010705,4670714609.552,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,secp256k1,256,COMPARE,50000,0.001023313,48860905.634,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,secp256k1,256,COMPARE,50000,0.000011178,4473073916.677,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,secp256k1,256,COMPARE,50000,0.001092698,45758296.615,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,secp256k1,256,COMPARE,50000,0.000305639,163591689.403,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,secp256k1,256,COMPARE,50000,0.001570806,31830792.710,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,secp256k1,256,REDUCE,6250,0.000150709,41470648.277,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,secp256k1,256,REDUCE,6250,0.000015071,414700065.387,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,secp256k1,256,REDUCE,6250,0.000230208,27149329.173,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,secp256k1,256,REDUCE,50000,0.000034816,1436121323.529,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,secp256k1,256,REDUCE,25200,0.000305931,82371516.075,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,secp256k1,256,REDUCE,25200,0.000981299,25680246.210,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,secp256k1,256,REDUCE,25200,0.000143662,175411720.838,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,secp256k1,256,REDUCE,25200,0.000824038,30581114.070,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,secp256k1,256,REDUCE,25200,0.000068253,369214417.665,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,secp256k1,256,REDUCE,25200,0.000803535,31361422.138,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,secp256k1,256,REDUCE,25200,0.000048755,516870078.225,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,secp256k1,256,REDUCE,25200,0.000708911,35547479.822,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,secp256k1,256,REDUCE,25200,0.000050692,497119820.408,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,secp256k1,256,REDUCE,25200,0.000606566,41545355.454,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,secp256k1,256,REDUCE,25200,0.000046760,538922278.485,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,secp256k1,256,REDUCE,25200,0.000787944,31981967.894,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,secp256k1,256,REDUCE,25200,0.009000789,2799754.547,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,secp256k1,256,REDUCE,25200,0.006729449,3744734.524,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,secp256k1,256,MODMUL,3125,0.000351708,8885211.568,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,secp256k1,256,MODMUL,3125,0.000041732,74883304.875,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,secp256k1,256,MODMUL,3125,0.000222824,14024496.577,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,secp256k1,256,MODMUL,50000,0.000082112,608924395.947,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,secp256k1,256,MODMUL,25200,0.000820039,30730245.202,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,secp256k1,256,MODMUL,25200,0.001500435,16795129.425,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,secp256k1,256,MODMUL,25200,0.000349732,72055173.789,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,secp256k1,256,MODMUL,25200,0.001025940,24562839.919,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,secp256k1,256,MODMUL,25200,0.000198828,126742713.411,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,secp256k1,256,MODMUL,25200,0.000885763,28450048.091,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,secp256k1,256,MODMUL,25200,0.000126939,198520541.443,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,secp256k1,256,MODMUL,25200,0.000783847,32149130.954,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,secp256k1,256,MODMUL,25200,0.000159577,157917493.499,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,secp256k1,256,MODMUL,25200,0.001014787,24832798.008,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,secp256k1,256,MODMUL,25200,0.000127415,197778914.330,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,secp256k1,256,MODMUL,25200,0.000856986,29405380.930,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,secp256k1,256,MODMUL,25200,0.035401152,711841.242,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,secp256k1,256,MODMUL,25200,0.034722638,725751.310,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,secp256k1,256,MODEXP,781,0.009061894,86185.073,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,secp256k1,256,MODEXP,781,0.000901035,866780.976,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,secp256k1,256,MODEXP,781,0.001770490,441120.821,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,secp256k1,256,MODEXP,50000,0.129795253,385222.101,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,secp256k1,256,MODEXP,25200,0.028052231,898324.272,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,secp256k1,256,MODEXP,25200,0.028676706,878762.016,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,secp256k1,256,MODEXP,25200,0.005143585,4899306.609,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,secp256k1,256,MODEXP,25200,0.005848243,4308986.485,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,secp256k1,256,MODEXP,25200,0.004764668,5288930.968,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,secp256k1,256,MODEXP,25200,0.005307180,4748284.381,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,secp256k1,256,MODEXP,25200,0.001942452,12973293.454,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,secp256k1,256,MODEXP,25200,0.002614347,9639118.288,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,secp256k1,256,MODEXP,25200,0.003568794,7061208.920,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,secp256k1,256,MODEXP,25200,0.004245451,5935765.113,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,secp256k1,256,MODEXP,25200,0.001976401,12750448.973,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,secp256k1,256,MODEXP,25200,0.002847574,8849638.380,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,secp256k1,256,MODEXP,25200,3.083374583,8172.864,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,secp256k1,256,MODEXP,25200,3.186029143,7909.532,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,secp256k1,256,EXPONENTIATION,781,0.003199521,244099.076,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,secp256k1,256,EXPONENTIATION,781,0.000340993,2290369.649,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,secp256k1,256,EXPONENTIATION,781,0.006614023,118082.444,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,secp256k1,256,EXPONENTIATION,25200,0.046255616,544798.712,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,secp256k1,256,EXPONENTIATION,25200,0.046655485,540129.419,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,secp256k1,256,EXPONENTIATION,25200,0.013750429,1832670.093,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,secp256k1,256,EXPONENTIATION,25200,0.014503939,1737459.042,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,secp256k1,256,EXPONENTIATION,25200,0.000870816,28938375.938,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,secp256k1,256,EXPONENTIATION,25200,0.001464364,17208836.368,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,secp256k1,256,EXPONENTIATION,25200,0.000545632,46184973.910,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,secp256k1,256,EXPONENTIATION,25200,0.001209548,20834228.896,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,secp256k1,256,EXPONENTIATION,25200,0.000672738,37458861.608,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,secp256k1,256,EXPONENTIATION,25200,0.001334185,18887935.552,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,secp256k1,256,EXPONENTIATION,25200,0.000540601,46614785.939,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,secp256k1,256,EXPONENTIATION,25200,0.001290541,19526694.228,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,secp256k1,256,EXPONENTIATION,25200,1.347861791,18696.279,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,secp256k1,256,EXPONENTIATION,25200,1.396261475,18048.195,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,secp256k1,256,DIVIDE,6250,0.000324009,19289607.507,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,secp256k1,256,DIVIDE,6250,0.000037401,167107328.137,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,secp256k1,256,DIVIDE,6250,0.000206796,30222989.159,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,secp256k1,256,DIVIDE,50000,0.000049216,1015929778.934,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,secp256k1,256,DIVIDE,25200,0.000449490,56063537.948,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,secp256k1,256,DIVIDE,25200,0.001249889,20161790.381,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,secp256k1,256,DIVIDE,25200,0.000291764,86371174.165,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,secp256k1,256,DIVIDE,25200,0.001092045,23075971.828,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,secp256k1,256,DIVIDE,25200,0.000153962,163676752.601,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,secp256k1,256,DIVIDE,25200,0.000861992,29234609.667,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,secp256k1,256,DIVIDE,25200,0.000115268,218620970.695,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,secp256k1,256,DIVIDE,25200,0.000908022,27752631.858,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,secp256k1,256,DIVIDE,25200,0.000115618,217959153.964,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,secp256k1,256,DIVIDE,25200,0.000899505,28015407.986,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,secp256k1,256,DIVIDE,25200,0.000110939,227151859.111,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,secp256k1,256,DIVIDE,25200,0.000929605,27108288.426,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,secp256k1,256,DIVIDE,25200,0.017327433,1454341.218,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,secp256k1,256,DIVIDE,25200,0.018477089,1363851.202,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,secp256k1,256,ISQRT,1562,0.000185088,8439242.282,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,secp256k1,256,ISQRT,1562,0.000016464,94872652.096,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,secp256k1,256,ISQRT,25200,0.005369841,4692876.377,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,secp256k1,256,ISQRT,25200,0.006049673,4165514.390,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,secp256k1,256,ISQRT,25200,0.004273434,5896896.976,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,secp256k1,256,ISQRT,25200,0.004959447,5081211.657,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,secp256k1,256,ISQRT,25200,0.001293943,19475355.974,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,secp256k1,256,ISQRT,25200,0.001945911,12950232.535,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,secp256k1,256,ISQRT,25200,0.000846929,29754560.153,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,secp256k1,256,ISQRT,25200,0.001511708,16669886.102,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,secp256k1,256,ISQRT,25200,0.000980118,25711189.915,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,secp256k1,256,ISQRT,25200,0.001607920,15672421.376,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,secp256k1,256,ISQRT,25200,0.000852234,29569343.902,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,secp256k1,256,ISQRT,25200,0.001829613,13773404.538,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,secp256k1,256,ISQRT,25200,0.256244108,98343.725,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,secp256k1,256,ISQRT,25200,0.260376381,96782.972,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,secp256k1,256,MODMUL_R2,50000,0.005560400,8992158.848,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,secp256k1,256,MODMUL_R2,50000,0.000648341,77119911.126,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,secp256k1,256,MODMUL_R2,50000,0.003377199,14805168.432,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,secp256k1,256,MODMUL_R2,50000,0.000034816,1436121323.529,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,secp256k1,256,MODMUL_R2,50000,0.000303051,164988716.568,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,secp256k1,256,MODMUL_R2,50000,0.001306800,38261401.982,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,secp256k1,256,MODMUL_R2,50000,0.000063705,784867797.378,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,secp256k1,256,MODMUL_R2,50000,0.001083016,46167370.148,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,secp256k1,256,MODMUL_R2,50000,0.000073556,679754178.126,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,secp256k1,256,MODMUL_R2,50000,0.001064863,46954397.050,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,secp256k1,256,MODMUL_R2,50000,0.000043341,1153642373.782,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,secp256k1,256,MODMUL_R2,50000,0.001033315,48387955.341,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,secp256k1,256,MODMUL_R2,50000,0.000060162,831089222.993,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,secp256k1,256,MODMUL_R2,50000,0.001126964,44366988.612,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,secp256k1,256,MODMUL_R2,50000,0.000040912,1222135850.469,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,secp256k1,256,MODMUL_R2,50000,0.001140208,43851647.634,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,secp256k1,256,MODMUL_R2,50000,0.016528110,3025149.278,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,secp256k1,256,MODMUL_R2,50000,0.018705994,2672940.024,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,rsa256(composite),256,ADD,50000,0.001014394,49290536.324,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,rsa256(composite),256,ADD,50000,0.000110791,451300228.292,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,rsa256(composite),256,ADD,50000,0.000101254,493805231.846,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,rsa256(composite),256,ADD,50000,0.000033376,1498082454.458,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,rsa256(composite),256,ADD,50000,0.000082721,604441418.133,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,rsa256(composite),256,ADD,50000,0.001087993,45956179.043,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,rsa256(composite),256,ADD,50000,0.000045247,1105045696.625,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,rsa256(composite),256,ADD,50000,0.001036185,48253930.813,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32,rsa256(composite),256,ADD,50000,0.000033293,1501816896.181,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32,rsa256(composite),256,ADD,50000,0.001050842,47580892.520,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,rsa256(composite),256,ADD,50000,0.000035360,1414027555.146,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,rsa256(composite),256,ADD,50000,0.001054399,47420378.519,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,rsa256(composite),256,ADD,50000,0.000032613,1533131161.708,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,rsa256(composite),256,ADD,50000,0.001022887,48881255.697,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,rsa256(composite),256,ADD,50000,0.000021803,2293263544.975,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,rsa256(composite),256,ADD,50000,0.001209909,41325421.685,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,rsa256(composite),256,ADD,50000,0.000021535,2321801610.678,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,rsa256(composite),256,ADD,50000,0.001157065,43212784.803,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,rsa256(composite),256,ADD,50000,0.000357496,139861705.061,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,rsa256(composite),256,ADD,50000,0.001717335,29114878.720,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACT,50000,0.001129323,44274312.804,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACT,50000,0.000094616,528454611.740,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACT,50000,0.000098740,506382940.247,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,rsa256(composite),256,SUBTRACT,50000,0.000032960,1516990291.262,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,rsa256(composite),256,SUBTRACT,50000,0.000081981,609897352.713,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,rsa256(composite),256,SUBTRACT,50000,0.001070465,46708672.485,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,rsa256(composite),256,SUBTRACT,50000,0.000045336,1102875663.239,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,rsa256(composite),256,SUBTRACT,50000,0.001038814,48131812.196,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32,rsa256(composite),256,SUBTRACT,50000,0.000032160,1554725222.758,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32,rsa256(composite),256,SUBTRACT,50000,0.001060625,47142014.683,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,rsa256(composite),256,SUBTRACT,50000,0.000035457,1410159065.512,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,rsa256(composite),256,SUBTRACT,50000,0.001057894,47263714.591,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,rsa256(composite),256,SUBTRACT,50000,0.000031985,1563233094.934,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,rsa256(composite),256,SUBTRACT,50000,0.001022087,48919514.553,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,rsa256(composite),256,SUBTRACT,50000,0.000021338,2343237727.779,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,rsa256(composite),256,SUBTRACT,50000,0.001071811,46650016.405,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,rsa256(composite),256,SUBTRACT,50000,0.000021405,2335900284.851,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,rsa256(composite),256,SUBTRACT,50000,0.001111604,44980046.391,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,rsa256(composite),256,SUBTRACT,50000,0.000396194,126200793.044,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,rsa256(composite),256,SUBTRACT,50000,0.001732327,28862910.738,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,rsa256(composite),256,ADDMOD,50000,0.003452138,14483777.939,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,rsa256(composite),256,ADDMOD,50000,0.000318128,157169439.039,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,rsa256(composite),256,ADDMOD,50000,0.001782159,28055858.078,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,rsa256(composite),256,ADDMOD,50000,0.000034144,1464386129.335,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,rsa256(composite),256,ADDMOD,50000,0.000113822,439282404.074,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,rsa256(composite),256,ADDMOD,50000,0.001109905,45048899.889,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,rsa256(composite),256,ADDMOD,50000,0.000059931,834292559.394,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,rsa256(composite),256,ADDMOD,50000,0.001051279,47561114.169,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32,rsa256(composite),256,ADDMOD,50000,0.000037368,1338043978.523,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32,rsa256(composite),256,ADDMOD,50000,0.001058550,47234424.513,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,rsa256(composite),256,ADDMOD,50000,0.000029168,1714207091.763,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,rsa256(composite),256,ADDMOD,50000,0.001069119,46767478.802,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,rsa256(composite),256,ADDMOD,50000,0.000023131,2161598888.743,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,rsa256(composite),256,ADDMOD,50000,0.001008251,49590826.909,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,rsa256(composite),256,ADDMOD,50000,0.000019435,2572680142.741,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,rsa256(composite),256,ADDMOD,50000,0.001058384,47241833.264,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,rsa256(composite),256,ADDMOD,50000,0.000019161,2609466936.223,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,rsa256(composite),256,ADDMOD,50000,0.001099036,45494415.608,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,rsa256(composite),256,ADDMOD,50000,0.000470223,106332522.427,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,rsa256(composite),256,ADDMOD,50000,0.001891202,26438212.279,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACTMOD,50000,0.002392555,20898156.806,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACTMOD,50000,0.000283160,176578304.559,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACTMOD,50000,0.002044414,24456891.907,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,rsa256(composite),256,SUBTRACTMOD,50000,0.000033856,1476843100.189,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,rsa256(composite),256,SUBTRACTMOD,50000,0.000122857,406977157.321,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,rsa256(composite),256,SUBTRACTMOD,50000,0.001149576,43494296.581,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,rsa256(composite),256,SUBTRACTMOD,50000,0.000064842,771105201.633,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,rsa256(composite),256,SUBTRACTMOD,50000,0.001061833,47088384.046,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32,rsa256(composite),256,SUBTRACTMOD,50000,0.000041467,1205777471.582,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32,rsa256(composite),256,SUBTRACTMOD,50000,0.001049947,47621451.659,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,rsa256(composite),256,SUBTRACTMOD,50000,0.000028854,1732861804.127,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,rsa256(composite),256,SUBTRACTMOD,50000,0.001111553,44982110.923,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,rsa256(composite),256,SUBTRACTMOD,50000,0.000030838,1621375630.105,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,rsa256(composite),256,SUBTRACTMOD,50000,0.001017518,49139178.814,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,rsa256(composite),256,SUBTRACTMOD,50000,0.000019814,2523467721.106,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,rsa256(composite),256,SUBTRACTMOD,50000,0.001145631,43644070.858,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,rsa256(composite),256,SUBTRACTMOD,50000,0.000019129,2613834131.943,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,rsa256(composite),256,SUBTRACTMOD,50000,0.001101943,45374396.781,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,rsa256(composite),256,SUBTRACTMOD,50000,0.000289625,172637034.543,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,rsa256(composite),256,SUBTRACTMOD,50000,0.001391264,35938542.248,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001246978,40096922.290,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000145239,344261326.073,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000264888,188759020.127,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.003366636,14851620.449,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.004604257,10859515.437,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000884189,56548995.203,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.002115928,23630293.491,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000244301,204665542.435,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001497562,33387598.993,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000067067,745523319.001,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001332470,37524296.641,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000067832,737115154.053,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001279309,39083598.620,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000030311,1649566211.096,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001305568,38297507.208,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000025171,1986412876.863,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001330576,37577710.902,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.003760474,13296196.148,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.005511388,9072124.831,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001253602,39885082.689,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000145144,344484292.941,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000276008,181154497.377,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000032928,1518464528.669,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000992759,50364691.467,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.002219724,22525322.935,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000289136,172928990.306,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001502306,33282166.973,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000100608,496978451.683,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001326502,37693121.254,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000107577,464783310.960,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001343451,37217583.566,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000112698,443663687.555,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001326359,37697184.732,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000058708,851672482.540,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001426591,35048587.513,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000053561,933515030.576,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001384439,36115711.720,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.004926367,10149467.156,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.006475900,7720934.531,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.009660627,5175647.505,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.001181270,42327325.823,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000266775,187423863.211,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000033184,1506750241.080,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000254236,196667650.592,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.001260715,39660033.781,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000061480,813272884.196,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000953297,52449552.211,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000028853,1732921233.591,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.001044913,47850874.105,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000035492,1408767981.527,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000942426,53054563.093,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000035254,1418279030.279,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.001039709,48090379.430,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000032667,1530596044.653,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.001090280,45859778.586,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000031843,1570204026.078,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.001104009,45289486.337,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.015529596,3219658.770,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.021769518,2296789.483,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,rsa256(composite),256,COMPARE,50000,0.000431023,116003089.100,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,rsa256(composite),256,COMPARE,50000,0.000047829,1045398034.215,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,rsa256(composite),256,COMPARE,50000,0.000058137,860032521.825,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,rsa256(composite),256,COMPARE,50000,0.000033760,1481042654.028,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,rsa256(composite),256,COMPARE,50000,0.000077599,644338281.824,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,rsa256(composite),256,COMPARE,50000,0.001087365,45982720.815,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,rsa256(composite),256,COMPARE,50000,0.000042869,1166344472.347,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,rsa256(composite),256,COMPARE,50000,0.000975591,51250984.557,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,rsa256(composite),256,COMPARE,50000,0.000028626,1746663642.076,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,rsa256(composite),256,COMPARE,50000,0.001055534,47369388.488,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,rsa256(composite),256,COMPARE,50000,0.000028948,1727236268.618,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,rsa256(composite),256,COMPARE,50000,0.001031752,48461257.832,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,rsa256(composite),256,COMPARE,50000,0.000011957,4181644723.980,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,rsa256(composite),256,COMPARE,50000,0.001090392,45855068.438,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,rsa256(composite),256,COMPARE,50000,0.000012361,4044977569.746,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,rsa256(composite),256,COMPARE,50000,0.001098097,45533317.640,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,rsa256(composite),256,COMPARE,50000,0.000298046,167759348.690,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,rsa256(composite),256,COMPARE,50000,0.001667794,29979721.450,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,rsa256(composite),256,REDUCE,6250,0.000239503,26095688.561,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,rsa256(composite),256,REDUCE,6250,0.000026942,231979826.297,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,rsa256(composite),256,REDUCE,6250,0.000225046,27772151.457,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,rsa256(composite),256,REDUCE,50000,0.000032768,1525878906.250,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,rsa256(composite),256,REDUCE,25200,0.000234032,107677573.042,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,rsa256(composite),256,REDUCE,25200,0.000882780,28546183.493,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,rsa256(composite),256,REDUCE,25200,0.000144844,173980280.283,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,rsa256(composite),256,REDUCE,25200,0.000830588,30339952.333,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,rsa256(composite),256,REDUCE,25200,0.000068699,366817670.004,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,rsa256(composite),256,REDUCE,25200,0.000752400,33492823.158,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,rsa256(composite),256,REDUCE,25200,0.000063089,399435631.909,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,rsa256(composite),256,REDUCE,25200,0.000733079,34375557.945,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,rsa256(composite),256,REDUCE,25200,0.000064131,392945757.284,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,rsa256(composite),256,REDUCE,25200,0.000698589,36072712.682,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,rsa256(composite),256,REDUCE,25200,0.000060391,417280712.549,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,rsa256(composite),256,REDUCE,25200,0.000728464,34593335.356,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,rsa256(composite),256,REDUCE,25200,0.010132238,2487110.941,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,rsa256(composite),256,REDUCE,25200,0.011158101,2258448.817,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,rsa256(composite),256,MODMUL,3125,0.000357473,8741925.653,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,rsa256(composite),256,MODMUL,3125,0.000041419,75448138.134,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,rsa256(composite),256,MODMUL,3125,0.000204008,15318033.609,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,rsa256(composite),256,MODMUL,50000,0.000061856,808329022.245,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,rsa256(composite),256,MODMUL,25200,0.000620281,40626747.961,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,rsa256(composite),256,MODMUL,25200,0.001285872,19597595.804,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,rsa256(composite),256,MODMUL,25200,0.000351065,71781582.184,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,rsa256(composite),256,MODMUL,25200,0.001018093,24752158.557,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,rsa256(composite),256,MODMUL,25200,0.000201842,124850135.911,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,rsa256(composite),256,MODMUL,25200,0.000890306,28304875.164,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,rsa256(composite),256,MODMUL,25200,0.000163752,153891246.825,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,rsa256(composite),256,MODMUL,25200,0.000831866,30293340.638,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,rsa256(composite),256,MODMUL,25200,0.000199528,126298061.165,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,rsa256(composite),256,MODMUL,25200,0.000864887,29136754.318,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,rsa256(composite),256,MODMUL,25200,0.000162971,154628735.462,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,rsa256(composite),256,MODMUL,25200,0.000813922,30961197.725,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,rsa256(composite),256,MODMUL,25200,0.027146245,928305.185,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,rsa256(composite),256,MODMUL,25200,0.071998699,350006.324,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,rsa256(composite),256,MODEXP,781,0.008710443,89662.489,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,rsa256(composite),256,MODEXP,781,0.000877266,890265.872,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,rsa256(composite),256,MODEXP,781,0.001751724,445846.494,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,rsa256(composite),256,MODEXP,50000,0.125685096,397819.643,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,rsa256(composite),256,MODEXP,25200,0.027949351,901630.954,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,rsa256(composite),256,MODEXP,25200,0.028662844,879187.006,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,rsa256(composite),256,MODEXP,25200,0.005136798,4905779.808,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,rsa256(composite),256,MODEXP,25200,0.005970055,4221066.637,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,rsa256(composite),256,MODEXP,25200,0.004797329,5252923.043,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,rsa256(composite),256,MODEXP,25200,0.005539256,4549347.419,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,rsa256(composite),256,MODEXP,25200,0.002572244,9796893.281,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,rsa256(composite),256,MODEXP,25200,0.003237243,7784401.693,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,rsa256(composite),256,MODEXP,25200,0.004749584,5305727.815,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,rsa256(composite),256,MODEXP,25200,0.005477191,4600898.520,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,rsa256(composite),256,MODEXP,25200,0.002614539,9638410.442,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,rsa256(composite),256,MODEXP,25200,0.003310980,7611039.601,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,rsa256(composite),256,MODEXP,25200,3.016494413,8354.068,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,rsa256(composite),256,MODEXP,25200,3.078868015,8184.826,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,rsa256(composite),256,EXPONENTIATION,781,0.003170085,246365.596,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,rsa256(composite),256,EXPONENTIATION,781,0.000362412,2155005.905,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,rsa256(composite),256,EXPONENTIATION,781,0.006259489,124770.558,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,rsa256(composite),256,EXPONENTIATION,25200,0.046551547,541335.393,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,rsa256(composite),256,EXPONENTIATION,25200,0.047200628,533891.202,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,rsa256(composite),256,EXPONENTIATION,25200,0.013860849,1818070.452,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,rsa256(composite),256,EXPONENTIATION,25200,0.014552015,1731718.941,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,rsa256(composite),256,EXPONENTIATION,25200,0.000864869,29137361.239,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,rsa256(composite),256,EXPONENTIATION,25200,0.001536513,16400772.299,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,rsa256(composite),256,EXPONENTIATION,25200,0.000706365,35675606.205,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,rsa256(composite),256,EXPONENTIATION,25200,0.001368972,18407973.186,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,rsa256(composite),256,EXPONENTIATION,25200,0.000873814,28839089.731,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,rsa256(composite),256,EXPONENTIATION,25200,0.001531651,16452834.114,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,rsa256(composite),256,EXPONENTIATION,25200,0.000701183,35939263.285,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,rsa256(composite),256,EXPONENTIATION,25200,0.001372669,18358395.489,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,rsa256(composite),256,EXPONENTIATION,25200,1.365594824,18453.497,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,rsa256(composite),256,EXPONENTIATION,25200,1.488851627,16925.797,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,rsa256(composite),256,DIVIDE,6250,0.000329734,18954701.282,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,rsa256(composite),256,DIVIDE,6250,0.000028247,221261643.809,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,rsa256(composite),256,DIVIDE,6250,0.000192563,32456893.524,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,rsa256(composite),256,DIVIDE,50000,0.000037088,1348144952.545,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,rsa256(composite),256,DIVIDE,25200,0.000454421,55455182.453,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,rsa256(composite),256,DIVIDE,25200,0.001245601,20231197.853,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,rsa256(composite),256,DIVIDE,25200,0.000294390,85600728.458,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,rsa256(composite),256,DIVIDE,25200,0.001081024,23311230.713,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,rsa256(composite),256,DIVIDE,25200,0.000158262,159229629.346,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,rsa256(composite),256,DIVIDE,25200,0.000991200,25423728.567,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,rsa256(composite),256,DIVIDE,25200,0.000148464,169738124.889,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,rsa256(composite),256,DIVIDE,25200,0.000960930,26224595.233,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,rsa256(composite),256,DIVIDE,25200,0.000150688,167232966.408,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,rsa256(composite),256,DIVIDE,25200,0.000960792,26228361.400,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,rsa256(composite),256,DIVIDE,25200,0.000149307,168779777.522,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,rsa256(composite),256,DIVIDE,25200,0.000948484,26568713.118,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,rsa256(composite),256,DIVIDE,25200,0.014961128,1684364.978,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,rsa256(composite),256,DIVIDE,25200,0.013450653,1873514.988,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,rsa256(composite),256,ISQRT,1562,0.000180376,8659685.935,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,rsa256(composite),256,ISQRT,1562,0.000016168,96609588.848,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,rsa256(composite),256,ISQRT,25200,0.005381081,4683073.910,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,rsa256(composite),256,ISQRT,25200,0.006048711,4166176.901,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,rsa256(composite),256,ISQRT,25200,0.004271265,5899891.480,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,rsa256(composite),256,ISQRT,25200,0.004937720,5103570.079,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,rsa256(composite),256,ISQRT,25200,0.001307469,19273879.862,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,rsa256(composite),256,ISQRT,25200,0.001894226,13303586.748,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,rsa256(composite),256,ISQRT,25200,0.001107273,22758614.834,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,rsa256(composite),256,ISQRT,25200,0.001768247,14251402.728,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,rsa256(composite),256,ISQRT,25200,0.001289333,19544989.342,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,rsa256(composite),256,ISQRT,25200,0.001975268,12757762.432,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,rsa256(composite),256,ISQRT,25200,0.001119034,22519422.689,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,rsa256(composite),256,ISQRT,25200,0.001787574,14097318.478,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,rsa256(composite),256,ISQRT,25200,0.235612283,106955.375,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,rsa256(composite),256,ISQRT,25200,0.204174580,123423.788,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,rsa256(composite),256,MODMUL_R2,50000,0.005648958,8851189.898,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,rsa256(composite),256,MODMUL_R2,50000,0.000647249,77250023.827,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,rsa256(composite),256,MODMUL_R2,50000,0.003280002,15243893.085,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,rsa256(composite),256,MODMUL_R2,50000,0.000033536,1490935114.504,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,rsa256(composite),256,MODMUL_R2,50000,0.000302174,165467571.615,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,rsa256(composite),256,MODMUL_R2,50000,0.001316401,37982347.280,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,rsa256(composite),256,MODMUL_R2,50000,0.000066605,750694514.994,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,rsa256(composite),256,MODMUL_R2,50000,0.001095088,45658430.096,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,rsa256(composite),256,MODMUL_R2,50000,0.000068939,725278870.752,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,rsa256(composite),256,MODMUL_R2,50000,0.001081127,46248035.444,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,rsa256(composite),256,MODMUL_R2,50000,0.000048451,1031970868.166,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,rsa256(composite),256,MODMUL_R2,50000,0.001055847,47355345.570,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,rsa256(composite),256,MODMUL_R2,50000,0.000059689,837675456.507,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,rsa256(composite),256,MODMUL_R2,50000,0.001169829,42741289.938,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,rsa256(composite),256,MODMUL_R2,50000,0.000043276,1155374832.812,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,rsa256(composite),256,MODMUL_R2,50000,0.001117487,44743249.751,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,rsa256(composite),256,MODMUL_R2,50000,0.014462782,3457149.535,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,rsa256(composite),256,MODMUL_R2,50000,0.015471346,3231780.866,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,brainpoolP512r1,512,ADD,25000,0.000542861,46052304.193,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,brainpoolP512r1,512,ADD,25000,0.000059090,423083419.584,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,brainpoolP512r1,512,ADD,25000,0.000054308,460341571.331,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,brainpoolP512r1,512,ADD,50000,0.000046752,1069472963.723,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,brainpoolP512r1,512,ADD,25200,0.000113296,222426206.310,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,brainpoolP512r1,512,ADD,25200,0.001137292,22157897.639,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,brainpoolP512r1,512,ADD,25200,0.000062711,401843293.359,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,brainpoolP512r1,512,ADD,25200,0.001060740,23756999.750,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32,brainpoolP512r1,512,ADD,25200,0.000065343,385657290.178,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32,brainpoolP512r1,512,ADD,25200,0.001105097,22803427.990,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,brainpoolP512r1,512,ADD,25200,0.000049129,512935207.797,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,brainpoolP512r1,512,ADD,25200,0.001070725,23535454.862,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,brainpoolP512r1,512,ADD,25200,0.000049804,505983280.598,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,brainpoolP512r1,512,ADD,25200,0.001085338,23218573.211,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,brainpoolP512r1,512,ADD,25200,0.000032337,779292971.389,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,brainpoolP512r1,512,ADD,25200,0.001086060,23203138.229,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,brainpoolP512r1,512,ADD,25200,0.000032169,783363165.056,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,brainpoolP512r1,512,ADD,25200,0.001109706,22708717.340,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,brainpoolP512r1,512,ADD,25200,0.000358120,70367471.755,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,brainpoolP512r1,512,ADD,25200,0.001685573,14950405.661,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACT,25000,0.000456333,54784506.354,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACT,25000,0.000049870,501305422.800,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACT,25000,0.000051943,481296830.398,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,brainpoolP512r1,512,SUBTRACT,50000,0.000046656,1071673525.377,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,brainpoolP512r1,512,SUBTRACT,25200,0.000112874,223257798.123,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,brainpoolP512r1,512,SUBTRACT,25200,0.001133309,22235771.379,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,brainpoolP512r1,512,SUBTRACT,25200,0.000062269,404695670.099,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,brainpoolP512r1,512,SUBTRACT,25200,0.001071272,23523437.329,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32,brainpoolP512r1,512,SUBTRACT,25200,0.000061477,409909520.361,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32,brainpoolP512r1,512,SUBTRACT,25200,0.001084073,23245667.442,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,brainpoolP512r1,512,SUBTRACT,25200,0.000049035,513918559.661,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,brainpoolP512r1,512,SUBTRACT,25200,0.001061849,23732187.797,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,brainpoolP512r1,512,SUBTRACT,25200,0.000048922,515105715.334,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,brainpoolP512r1,512,SUBTRACT,25200,0.001050795,23981842.504,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,brainpoolP512r1,512,SUBTRACT,25200,0.000032597,773077547.333,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,brainpoolP512r1,512,SUBTRACT,25200,0.001084852,23228974.942,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,brainpoolP512r1,512,SUBTRACT,25200,0.000032453,776507975.128,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,brainpoolP512r1,512,SUBTRACT,25200,0.001106844,22767435.615,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,brainpoolP512r1,512,SUBTRACT,25200,0.000400750,62882096.759,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,brainpoolP512r1,512,SUBTRACT,25200,0.001681568,14986012.932,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,brainpoolP512r1,512,ADDMOD,25000,0.001395167,17919001.724,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,brainpoolP512r1,512,ADDMOD,25000,0.000169848,147190139.271,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,brainpoolP512r1,512,ADDMOD,25000,0.001113687,22447952.020,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,brainpoolP512r1,512,ADDMOD,50000,0.000047104,1061480978.261,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,brainpoolP512r1,512,ADDMOD,25200,0.000152198,165573774.941,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,brainpoolP512r1,512,ADDMOD,25200,0.001184967,21266414.541,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,brainpoolP512r1,512,ADDMOD,25200,0.000078945,319209522.797,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,brainpoolP512r1,512,ADDMOD,25200,0.001076961,23399176.136,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32,brainpoolP512r1,512,ADDMOD,25200,0.000076944,327510880.980,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32,brainpoolP512r1,512,ADDMOD,25200,0.001097794,22955126.142,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,brainpoolP512r1,512,ADDMOD,25200,0.000059064,426655744.129,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,brainpoolP512r1,512,ADDMOD,25200,0.000941525,26765088.278,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,brainpoolP512r1,512,ADDMOD,25200,0.000059019,426981015.068,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,brainpoolP512r1,512,ADDMOD,25200,0.001071115,23526885.639,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,brainpoolP512r1,512,ADDMOD,25200,0.000025180,1000794526.529,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,brainpoolP512r1,512,ADDMOD,25200,0.001055748,23869332.495,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,brainpoolP512r1,512,ADDMOD,25200,0.000024763,1017646340.985,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,brainpoolP512r1,512,ADDMOD,25200,0.001110024,22702211.976,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,brainpoolP512r1,512,ADDMOD,25200,0.000432708,58237886.612,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,brainpoolP512r1,512,ADDMOD,25200,0.001813669,13894486.661,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACTMOD,25000,0.001799886,13889764.768,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACTMOD,25000,0.001270707,19674086.973,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACTMOD,25000,0.002950414,8473387.134,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000047104,1061480978.261,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,brainpoolP512r1,512,SUBTRACTMOD,25200,0.000174919,144066684.105,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,brainpoolP512r1,512,SUBTRACTMOD,25200,0.001177257,21405691.636,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,brainpoolP512r1,512,SUBTRACTMOD,25200,0.000084554,298034424.841,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,brainpoolP512r1,512,SUBTRACTMOD,25200,0.001084148,23244059.324,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32,brainpoolP512r1,512,SUBTRACTMOD,25200,0.000082747,304542684.765,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32,brainpoolP512r1,512,SUBTRACTMOD,25200,0.001089898,23121429.791,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,brainpoolP512r1,512,SUBTRACTMOD,25200,0.000059314,424857438.111,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,brainpoolP512r1,512,SUBTRACTMOD,25200,0.001086895,23185312.185,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,brainpoolP512r1,512,SUBTRACTMOD,25200,0.000060262,418173883.480,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,brainpoolP512r1,512,SUBTRACTMOD,25200,0.001061809,23733082.223,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,brainpoolP512r1,512,SUBTRACTMOD,25200,0.000026975,934198273.810,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,brainpoolP512r1,512,SUBTRACTMOD,25200,0.001101579,22876253.228,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,brainpoolP512r1,512,SUBTRACTMOD,25200,0.000025845,975044009.781,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,brainpoolP512r1,512,SUBTRACTMOD,25200,0.001052821,23935692.844,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,brainpoolP512r1,512,SUBTRACTMOD,25200,0.000462189,54523150.021,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,brainpoolP512r1,512,SUBTRACTMOD,25200,0.001832642,13750639.927,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.001968735,12698509.358,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.002374482,10528612.107,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.001121007,22301377.376,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25200,0.012039789,2093059.937,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25200,0.013325419,1891122.521,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25200,0.003124957,8064110.957,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25200,0.004370672,5765703.741,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25200,0.001021889,24660212.382,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25200,0.002275891,11072586.561,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25200,0.000141457,178146027.739,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25200,0.001279575,19694039.136,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25200,0.000141947,177531042.391,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25200,0.001384759,18198112.429,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25200,0.000072565,347274775.499,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25200,0.001360598,18521267.975,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25200,0.000068528,367732921.066,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25200,0.001403129,17959859.839,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25200,0.017382529,1449731.508,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25200,0.019152236,1315773.258,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.001417737,17633736.160,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.001262791,19797409.177,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.001344342,18596465.144,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000046592,1073145604.396,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25200,0.003483101,7234932.325,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25200,0.004710863,5349338.333,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25200,0.000926921,27186782.949,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25200,0.002164322,11643369.251,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25200,0.000389513,64696170.779,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25200,0.001456494,17301821.969,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25200,0.000339725,74177644.426,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25200,0.001588260,15866419.873,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25200,0.000341376,73818897.922,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25200,0.001570516,16045681.761,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25200,0.000152183,165590113.933,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25200,0.001418712,17762590.344,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25200,0.000136463,184665439.540,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25200,0.001449481,17385533.142,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25200,0.010264369,2455094.904,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25200,0.011929620,2112389.163,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.014795682,1689682.164,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.001133485,22055871.960,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000272909,91605634.973,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000044768,1116869192.280,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25200,0.000490871,51337317.132,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25200,0.001499167,16809334.748,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25200,0.000144454,174449985.489,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25200,0.001225097,20569799.972,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25200,0.000055480,454217803.948,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25200,0.000962624,26178445.413,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25200,0.000076097,331156288.711,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25200,0.001128846,22323682.895,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25200,0.000058050,434108486.272,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25200,0.001069835,23555034.634,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25200,0.000062924,400483155.296,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25200,0.001154092,21835347.710,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25200,0.000045212,557374324.657,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25200,0.001161947,21687736.143,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25200,0.028137423,895604.405,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25200,0.090723344,277767.539,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,brainpoolP512r1,512,COMPARE,25000,0.000456055,54817950.107,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,brainpoolP512r1,512,COMPARE,25000,0.001343610,18606589.737,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,brainpoolP512r1,512,COMPARE,25000,0.000032620,766400980.239,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,brainpoolP512r1,512,COMPARE,50000,0.000046752,1069472963.723,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,brainpoolP512r1,512,COMPARE,25200,0.000094960,265374956.134,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,brainpoolP512r1,512,COMPARE,25200,0.001112906,22643421.809,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,brainpoolP512r1,512,COMPARE,25200,0.000053484,471169011.019,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,brainpoolP512r1,512,COMPARE,25200,0.001052833,23935419.577,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,brainpoolP512r1,512,COMPARE,25200,0.000039171,643333226.991,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,brainpoolP512r1,512,COMPARE,25200,0.001069161,23569883.451,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,brainpoolP512r1,512,COMPARE,25200,0.000038531,654018881.026,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,brainpoolP512r1,512,COMPARE,25200,0.001040848,24211027.905,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,brainpoolP512r1,512,COMPARE,25200,0.000013193,1910103962.148,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,brainpoolP512r1,512,COMPARE,25200,0.001119186,22516364.657,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,brainpoolP512r1,512,COMPARE,25200,0.000012921,1950314231.564,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,brainpoolP512r1,512,COMPARE,25200,0.001116998,22560470.139,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,brainpoolP512r1,512,COMPARE,25200,0.000226346,111333980.295,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,brainpoolP512r1,512,COMPARE,25200,0.001328168,18973503.347,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,brainpoolP512r1,512,REDUCE,3125,0.000177845,17571530.157,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,brainpoolP512r1,512,REDUCE,3125,0.000016218,192684172.467,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,brainpoolP512r1,512,REDUCE,3125,0.000127132,24580775.029,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,brainpoolP512r1,512,REDUCE,50000,0.000044928,1112891737.892,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,brainpoolP512r1,512,REDUCE,25200,0.000691373,36449211.090,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,brainpoolP512r1,512,REDUCE,25200,0.001705174,14778550.588,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,brainpoolP512r1,512,REDUCE,25200,0.000539245,46732005.226,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,brainpoolP512r1,512,REDUCE,25200,0.001552318,16233786.950,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,brainpoolP512r1,512,REDUCE,25200,0.000188881,133417331.552,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,brainpoolP512r1,512,REDUCE,25200,0.001135073,22201215.428,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,brainpoolP512r1,512,REDUCE,25200,0.000180094,139926920.755,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,brainpoolP512r1,512,REDUCE,25200,0.001184935,21266989.633,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,brainpoolP512r1,512,REDUCE,25200,0.000186311,135257697.668,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,brainpoolP512r1,512,REDUCE,25200,0.001243337,20268036.754,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,brainpoolP512r1,512,REDUCE,25200,0.000189546,132949257.061,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,brainpoolP512r1,512,REDUCE,25200,0.001310791,19225032.745,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,brainpoolP512r1,512,REDUCE,25200,0.031840405,791447.219,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,brainpoolP512r1,512,REDUCE,25200,0.034253829,735684.177,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL,1562,0.000474175,3294145.236,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL,1562,0.000037855,41262497.134,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL,1562,0.000156436,9984939.652,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,brainpoolP512r1,512,MODMUL,50000,0.000188416,265370244.565,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,brainpoolP512r1,512,MODMUL,25200,0.001916810,13146842.752,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,brainpoolP512r1,512,MODMUL,25200,0.002930628,8598839.595,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,brainpoolP512r1,512,MODMUL,25200,0.001360187,18526864.500,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,brainpoolP512r1,512,MODMUL,25200,0.002359860,10678599.602,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,brainpoolP512r1,512,MODMUL,25200,0.000683630,36862045.165,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,brainpoolP512r1,512,MODMUL,25200,0.001720288,14648709.884,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,brainpoolP512r1,512,MODMUL,25200,0.000557924,45167441.381,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,brainpoolP512r1,512,MODMUL,25200,0.001586832,15880697.958,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,brainpoolP512r1,512,MODMUL,25200,0.000820871,30699098.660,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,brainpoolP512r1,512,MODMUL,25200,0.001876908,13426337.322,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,brainpoolP512r1,512,MODMUL,25200,0.000663062,38005495.451,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,brainpoolP512r1,512,MODMUL,25200,0.001752863,14376480.096,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,brainpoolP512r1,512,MODMUL,25200,0.174810158,144156.382,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,brainpoolP512r1,512,MODMUL,25200,0.109938380,229219.314,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,brainpoolP512r1,512,MODEXP,390,0.023682176,16468.081,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,brainpoolP512r1,512,MODEXP,390,0.002501188,155925.905,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,brainpoolP512r1,512,MODEXP,390,0.002530026,154148.614,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,brainpoolP512r1,512,MODEXP,50000,1.160501957,43084.805,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,brainpoolP512r1,512,MODEXP,25200,0.404135209,62355.369,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,brainpoolP512r1,512,MODEXP,25200,0.405182505,62194.196,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,brainpoolP512r1,512,MODEXP,25200,0.037912705,664684.834,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,brainpoolP512r1,512,MODEXP,25200,0.039488458,638161.156,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,brainpoolP512r1,512,MODEXP,25200,0.030028638,839198.901,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,brainpoolP512r1,512,MODEXP,25200,0.031872879,790640.846,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,brainpoolP512r1,512,MODEXP,25200,0.018097940,1392423.669,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,brainpoolP512r1,512,MODEXP,25200,0.019078167,1320881.614,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,brainpoolP512r1,512,MODEXP,25200,0.032120299,784550.605,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,brainpoolP512r1,512,MODEXP,25200,0.033140996,760387.527,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,brainpoolP512r1,512,MODEXP,25200,0.017776076,1417635.702,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,brainpoolP512r1,512,MODEXP,25200,0.018983059,1327499.431,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,brainpoolP512r1,512,EXPONENTIATION,390,0.005772633,67560.159,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,brainpoolP512r1,512,EXPONENTIATION,390,0.000969171,402405.767,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,brainpoolP512r1,512,EXPONENTIATION,390,0.011619023,33565.645,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,brainpoolP512r1,512,EXPONENTIATION,25200,0.360518771,69899.273,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,brainpoolP512r1,512,EXPONENTIATION,25200,0.361881645,69636.027,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,brainpoolP512r1,512,EXPONENTIATION,25200,0.090707829,277815.049,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,brainpoolP512r1,512,EXPONENTIATION,25200,0.091740873,274686.726,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,brainpoolP512r1,512,EXPONENTIATION,25200,0.028032857,898945.120,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,brainpoolP512r1,512,EXPONENTIATION,25200,0.029078543,866618.386,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,brainpoolP512r1,512,EXPONENTIATION,25200,0.026294839,958362.971,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,brainpoolP512r1,512,EXPONENTIATION,25200,0.027732791,908671.615,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,brainpoolP512r1,512,EXPONENTIATION,25200,0.028541885,882912.954,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,brainpoolP512r1,512,EXPONENTIATION,25200,0.029599778,851357.736,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,brainpoolP512r1,512,EXPONENTIATION,25200,0.026723855,942977.725,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,brainpoolP512r1,512,EXPONENTIATION,25200,0.027628714,912094.569,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,brainpoolP512r1,512,EXPONENTIATION,25200,10.451844374,2411.058,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,brainpoolP512r1,512,EXPONENTIATION,25200,11.279578649,2234.126,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,brainpoolP512r1,512,DIVIDE,3125,0.000351762,8883858.838,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,brainpoolP512r1,512,DIVIDE,3125,0.000020801,150234779.490,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,brainpoolP512r1,512,DIVIDE,3125,0.000109917,28430657.854,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,brainpoolP512r1,512,DIVIDE,50000,0.000092192,542346407.497,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,brainpoolP512r1,512,DIVIDE,25200,0.001360965,18516273.509,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,brainpoolP512r1,512,DIVIDE,25200,0.002445672,10303916.476,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,brainpoolP512r1,512,DIVIDE,25200,0.001226592,20544728.909,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,brainpoolP512r1,512,DIVIDE,25200,0.002461867,10236133.798,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,brainpoolP512r1,512,DIVIDE,25200,0.000393722,64004551.008,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,brainpoolP512r1,512,DIVIDE,25200,0.001465515,17195320.408,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,brainpoolP512r1,512,DIVIDE,25200,0.000376832,66873301.262,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,brainpoolP512r1,512,DIVIDE,25200,0.001603386,15716739.297,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,brainpoolP512r1,512,DIVIDE,25200,0.000374890,67219719.437,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,brainpoolP512r1,512,DIVIDE,25200,0.001638050,15384145.902,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,brainpoolP512r1,512,DIVIDE,25200,0.000361283,69751415.717,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,brainpoolP512r1,512,DIVIDE,25200,0.001681161,14989641.117,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,brainpoolP512r1,512,DIVIDE,25200,0.087719590,287279.044,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,brainpoolP512r1,512,DIVIDE,25200,0.086317023,291947.047,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,brainpoolP512r1,512,ISQRT,781,0.000204392,3821097.060,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,brainpoolP512r1,512,ISQRT,781,0.000014358,54396650.583,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,brainpoolP512r1,512,ISQRT,25200,0.032036139,786611.645,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,brainpoolP512r1,512,ISQRT,25200,0.033056326,762335.173,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,brainpoolP512r1,512,ISQRT,25200,0.027587409,913460.194,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,brainpoolP512r1,512,ISQRT,25200,0.028598179,881174.987,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,brainpoolP512r1,512,ISQRT,25200,0.005649709,4460406.724,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,brainpoolP512r1,512,ISQRT,25200,0.006482432,3887429.912,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,brainpoolP512r1,512,ISQRT,25200,0.005077193,4963372.481,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,brainpoolP512r1,512,ISQRT,25200,0.006108792,4125201.838,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,brainpoolP512r1,512,ISQRT,25200,0.006664351,3781313.452,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,brainpoolP512r1,512,ISQRT,25200,0.007951135,3169358.844,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,brainpoolP512r1,512,ISQRT,25200,0.006144775,4101045.203,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,brainpoolP512r1,512,ISQRT,25200,0.007240747,3480303.904,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,brainpoolP512r1,512,ISQRT,25200,1.025075132,24583.564,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,brainpoolP512r1,512,ISQRT,25200,1.155611173,21806.643,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL_R2,25000,0.005585699,4475715.579,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL_R2,25000,0.000579379,43149650.013,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL_R2,25000,0.002397823,10426124.021,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,brainpoolP512r1,512,MODMUL_R2,50000,0.000063552,786757301.108,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,brainpoolP512r1,512,MODMUL_R2,25200,0.000706417,35672979.659,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,brainpoolP512r1,512,MODMUL_R2,25200,0.001566585,16085945.124,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,brainpoolP512r1,512,MODMUL_R2,25200,0.000152706,165022991.862,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,brainpoolP512r1,512,MODMUL_R2,25200,0.001143309,22041285.491,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,brainpoolP512r1,512,MODMUL_R2,25200,0.000095381,264203579.103,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,brainpoolP512r1,512,MODMUL_R2,25200,0.000905250,27837614.060,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,brainpoolP512r1,512,MODMUL_R2,25200,0.000094466,266762566.246,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,brainpoolP512r1,512,MODMUL_R2,25200,0.001116076,22579107.138,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,brainpoolP512r1,512,MODMUL_R2,25200,0.000089159,282641116.271,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,brainpoolP512r1,512,MODMUL_R2,25200,0.001237406,20365183.015,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,brainpoolP512r1,512,MODMUL_R2,25200,0.000084214,299237657.163,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,brainpoolP512r1,512,MODMUL_R2,25200,0.001175950,21429482.723,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,brainpoolP512r1,512,MODMUL_R2,25200,0.098038623,257041.554,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,brainpoolP512r1,512,MODMUL_R2,25200,0.099534243,253179.200,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,p1024,1024,ADD,12500,0.000353262,35384516.425,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,p1024,1024,ADD,12500,0.000033736,370519362.936,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,p1024,1024,ADD,12500,0.000070432,177476143.530,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,p1024,1024,ADD,50000,0.000072384,690760389.036,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,p1024,1024,ADD,25200,0.000459661,54823009.823,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,p1024,1024,ADD,25200,0.001893932,13305651.947,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,p1024,1024,ADD,25200,0.000240744,104675504.080,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,p1024,1024,ADD,25200,0.002032937,12395858.877,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32,p1024,1024,ADD,25200,0.000337434,74681277.165,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32,p1024,1024,ADD,25200,0.002023337,12454672.643,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,p1024,1024,ADD,25200,0.000126547,199135490.419,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,p1024,1024,ADD,25200,0.001551919,16237960.815,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,p1024,1024,ADD,25200,0.000123125,204670040.435,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,p1024,1024,ADD,25200,0.001937896,13003793.716,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,p1024,1024,ADD,25200,0.000050769,496365766.576,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,p1024,1024,ADD,25200,0.002211118,11396949.414,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,p1024,1024,ADD,25200,0.000049479,509306827.345,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,p1024,1024,ADD,25200,0.001979123,12732912.597,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,p1024,1024,ADD,25200,0.000497852,50617453.955,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,p1024,1024,ADD,25200,0.002672231,9430322.353,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,p1024,1024,SUBTRACT,12500,0.000385532,32422743.985,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,p1024,1024,SUBTRACT,12500,0.000041229,303183435.917,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,p1024,1024,SUBTRACT,12500,0.000044399,281536800.330,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,p1024,1024,SUBTRACT,50000,0.000073152,683508311.461,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,p1024,1024,SUBTRACT,25200,0.000458650,54943856.980,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,p1024,1024,SUBTRACT,25200,0.001902173,13248006.483,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,p1024,1024,SUBTRACT,25200,0.000239254,105327390.640,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,p1024,1024,SUBTRACT,25200,0.002086938,12075107.111,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32,p1024,1024,SUBTRACT,25200,0.000339971,74123966.940,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32,p1024,1024,SUBTRACT,25200,0.002032936,12395864.911,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,p1024,1024,SUBTRACT,25200,0.000125495,200804783.439,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,p1024,1024,SUBTRACT,25200,0.001579715,15952244.579,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,p1024,1024,SUBTRACT,25200,0.000124657,202154708.082,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,p1024,1024,SUBTRACT,25200,0.001925345,13088563.441,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,p1024,1024,SUBTRACT,25200,0.000054581,461699096.493,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,p1024,1024,SUBTRACT,25200,0.002059600,12235385.536,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,p1024,1024,SUBTRACT,25200,0.000050222,501772067.856,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,p1024,1024,SUBTRACT,25200,0.001969822,12793034.007,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,p1024,1024,SUBTRACT,25200,0.000572128,44046086.229,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,p1024,1024,SUBTRACT,25200,0.003075991,8192481.784,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,p1024,1024,ADDMOD,12500,0.001454349,8594912.823,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,p1024,1024,ADDMOD,12500,0.000128014,97645567.369,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,p1024,1024,ADDMOD,12500,0.000586528,21311843.501,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,p1024,1024,ADDMOD,50000,0.000071680,697544642.857,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,p1024,1024,ADDMOD,25200,0.000616280,40890505.400,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,p1024,1024,ADDMOD,25200,0.002034166,12388369.406,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,p1024,1024,ADDMOD,25200,0.000319874,78781019.673,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,p1024,1024,ADDMOD,25200,0.002150193,11719878.156,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32,p1024,1024,ADDMOD,25200,0.000364088,69214038.507,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32,p1024,1024,ADDMOD,25200,0.002103390,11980659.806,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,p1024,1024,ADDMOD,25200,0.000128450,196185303.589,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,p1024,1024,ADDMOD,25200,0.001579572,15953688.913,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,p1024,1024,ADDMOD,25200,0.000128789,195668900.006,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,p1024,1024,ADDMOD,25200,0.001937838,13004182.950,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,p1024,1024,ADDMOD,25200,0.000043945,573444344.651,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,p1024,1024,ADDMOD,25200,0.001992036,12650373.748,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,p1024,1024,ADDMOD,25200,0.000043551,578631935.473,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,p1024,1024,ADDMOD,25200,0.001975891,12753739.921,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,p1024,1024,ADDMOD,25200,0.000732901,34383908.072,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,p1024,1024,ADDMOD,25200,0.003294985,7647986.226,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,p1024,1024,SUBTRACTMOD,12500,0.000775283,16123139.439,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,p1024,1024,SUBTRACTMOD,12500,0.000107321,116472741.808,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,p1024,1024,SUBTRACTMOD,12500,0.000694116,18008523.656,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,p1024,1024,SUBTRACTMOD,50000,0.000071680,697544642.857,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,p1024,1024,SUBTRACTMOD,25200,0.000618932,40715295.883,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,p1024,1024,SUBTRACTMOD,25200,0.002014517,12509201.968,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,p1024,1024,SUBTRACTMOD,25200,0.000319818,78794813.155,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,p1024,1024,SUBTRACTMOD,25200,0.002157346,11681019.041,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32,p1024,1024,SUBTRACTMOD,25200,0.000363045,69412883.901,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32,p1024,1024,SUBTRACTMOD,25200,0.002204978,11428685.485,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,p1024,1024,SUBTRACTMOD,25200,0.000128421,196229586.882,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,p1024,1024,SUBTRACTMOD,25200,0.001997688,12614582.461,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,p1024,1024,SUBTRACTMOD,25200,0.000129927,193955055.213,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,p1024,1024,SUBTRACTMOD,25200,0.001965870,12818751.977,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,p1024,1024,SUBTRACTMOD,25200,0.000044337,568374312.888,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,p1024,1024,SUBTRACTMOD,25200,0.001965885,12818654.244,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,p1024,1024,SUBTRACTMOD,25200,0.000044207,570045450.662,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,p1024,1024,SUBTRACTMOD,25200,0.001964624,12826882.007,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,p1024,1024,SUBTRACTMOD,25200,0.000681411,36982086.777,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,p1024,1024,SUBTRACTMOD,25200,0.002873955,8768404.505,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.002391303,5227276.751,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.000235139,53160045.893,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.000345276,36202973.673,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,25200,0.062929518,400448.006,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,25200,0.065007286,387648.855,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,p1024,1024,MULTIPLYOPERANDSCANNING,25200,0.015795083,1595433.213,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,p1024,1024,MULTIPLYOPERANDSCANNING,25200,0.018109296,1391550.504,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32,p1024,1024,MULTIPLYOPERANDSCANNING,25200,0.008598995,2930575.031,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32,p1024,1024,MULTIPLYOPERANDSCANNING,25200,0.010886361,2314823.106,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,p1024,1024,MULTIPLYOPERANDSCANNING,25200,0.000303181,83118663.796,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,p1024,1024,MULTIPLYOPERANDSCANNING,25200,0.002397054,10512904.628,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,p1024,1024,MULTIPLYOPERANDSCANNING,25200,0.000292711,86091745.598,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,p1024,1024,MULTIPLYOPERANDSCANNING,25200,0.002560515,9841770.167,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,p1024,1024,MULTIPLYOPERANDSCANNING,25200,0.000190822,132060244.522,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,p1024,1024,MULTIPLYOPERANDSCANNING,25200,0.002563224,9831368.622,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,p1024,1024,MULTIPLYOPERANDSCANNING,25200,0.000183993,136961739.639,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,p1024,1024,MULTIPLYOPERANDSCANNING,25200,0.002858157,8816870.428,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,25200,0.101146607,249143.305,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,25200,0.114560303,219971.485,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.002391103,5227711.818,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.000234975,53197260.797,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.001421981,8790556.541,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000080224,623254886.318,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,25200,0.027073930,930784.707,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,25200,0.029267938,861010.434,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,p1024,1024,MULTIPLYPRODUCTSCANNING,25200,0.006870213,3668008.550,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,p1024,1024,MULTIPLYPRODUCTSCANNING,25200,0.009095575,2770578.004,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32,p1024,1024,MULTIPLYPRODUCTSCANNING,25200,0.001843344,13670807.090,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32,p1024,1024,MULTIPLYPRODUCTSCANNING,25200,0.004150060,6072201.358,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,p1024,1024,MULTIPLYPRODUCTSCANNING,25200,0.001811467,13911376.842,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,p1024,1024,MULTIPLYPRODUCTSCANNING,25200,0.003667930,6870360.130,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,p1024,1024,MULTIPLYPRODUCTSCANNING,25200,0.001777190,14179688.259,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,p1024,1024,MULTIPLYPRODUCTSCANNING,25200,0.004060602,6205976.326,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,p1024,1024,MULTIPLYPRODUCTSCANNING,25200,0.000516202,48818097.419,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,p1024,1024,MULTIPLYPRODUCTSCANNING,25200,0.002806316,8979744.313,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,p1024,1024,MULTIPLYPRODUCTSCANNING,25200,0.000489421,51489411.505,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,p1024,1024,MULTIPLYPRODUCTSCANNING,25200,0.003267690,7711869.885,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,25200,0.036501865,690375.684,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,25200,0.087908035,286663.216,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.021653580,577271.749,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.001512227,8265954.792,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.000412577,30297376.453,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000112544,444270685.243,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,25200,0.002690747,9365428.999,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,25200,0.004173561,6038009.269,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,p1024,1024,MONTGOMERYMULTIPLICATION,25200,0.000542419,46458549.387,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,p1024,1024,MONTGOMERYMULTIPLICATION,25200,0.002324959,10838900.809,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32,p1024,1024,MONTGOMERYMULTIPLICATION,25200,0.000144127,174845796.137,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32,p1024,1024,MONTGOMERYMULTIPLICATION,25200,0.001989180,12668536.852,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,p1024,1024,MONTGOMERYMULTIPLICATION,25200,0.000236718,106455789.990,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,p1024,1024,MONTGOMERYMULTIPLICATION,25200,0.001659188,15188152.432,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,p1024,1024,MONTGOMERYMULTIPLICATION,25200,0.000193166,130457745.366,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,p1024,1024,MONTGOMERYMULTIPLICATION,25200,0.002025435,12441771.711,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,p1024,1024,MONTGOMERYMULTIPLICATION,25200,0.000164999,152728211.041,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,p1024,1024,MONTGOMERYMULTIPLICATION,25200,0.002019879,12475994.845,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,p1024,1024,MONTGOMERYMULTIPLICATION,25200,0.000118140,213306220.068,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,p1024,1024,MONTGOMERYMULTIPLICATION,25200,0.002094557,12031183.769,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,25200,0.198094235,127212.183,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,25200,0.286566427,87937.726,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,p1024,1024,COMPARE,12500,0.000188293,66385722.509,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,p1024,1024,COMPARE,12500,0.000013151,950497891.184,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,p1024,1024,COMPARE,12500,0.000444438,28125400.038,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,p1024,1024,COMPARE,50000,0.000072384,690760389.036,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,p1024,1024,COMPARE,25200,0.000249040,101188566.962,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,p1024,1024,COMPARE,25200,0.002123836,11865322.983,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,p1024,1024,COMPARE,25200,0.000130291,193413191.640,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,p1024,1024,COMPARE,25200,0.001933520,13033224.443,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,p1024,1024,COMPARE,25200,0.000073068,344884129.946,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,p1024,1024,COMPARE,25200,0.001496859,16835253.221,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,p1024,1024,COMPARE,25200,0.000073295,343816088.914,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,p1024,1024,COMPARE,25200,0.001903325,13239987.922,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,p1024,1024,COMPARE,25200,0.000013801,1825957525.914,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,p1024,1024,COMPARE,25200,0.001908922,13201167.885,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,p1024,1024,COMPARE,25200,0.000014528,1734582467.323,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,p1024,1024,COMPARE,25200,0.001898547,13273308.523,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,p1024,1024,COMPARE,25200,0.000365257,68992517.361,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,p1024,1024,COMPARE,25200,0.002894907,8704942.864,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,p1024,1024,REDUCE,1562,0.000069490,22478037.384,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,p1024,1024,REDUCE,1562,0.000004426,352889416.963,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,p1024,1024,REDUCE,1562,0.000078835,19813620.645,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,p1024,1024,REDUCE,50000,0.000073280,682314410.480,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,p1024,1024,REDUCE,25200,0.003647346,6909133.397,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,p1024,1024,REDUCE,25200,0.005757771,4376693.692,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,p1024,1024,REDUCE,25200,0.001977524,12743208.034,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,p1024,1024,REDUCE,25200,0.003789111,6650636.552,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,p1024,1024,REDUCE,25200,0.000505013,49899707.141,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,p1024,1024,REDUCE,25200,0.001951353,12914116.425,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,p1024,1024,REDUCE,25200,0.000498798,50521452.921,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,p1024,1024,REDUCE,25200,0.002330832,10811589.984,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,p1024,1024,REDUCE,25200,0.000444510,56691639.104,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,p1024,1024,REDUCE,25200,0.002363366,10662758.042,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,p1024,1024,REDUCE,25200,0.000430922,58479261.447,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,p1024,1024,REDUCE,25200,0.002283391,11036217.675,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,p1024,1024,REDUCE,25200,0.176128535,143077.327,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,p1024,1024,REDUCE,25200,0.226360959,111326.618,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,p1024,1024,MODMUL,781,0.000610796,1278660.116,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,p1024,1024,MODMUL,781,0.000045706,17087471.580,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,p1024,1024,MODMUL,781,0.000155998,5006482.466,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,p1024,1024,MODMUL,50000,0.000559264,89403215.655,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,p1024,1024,MODMUL,25200,0.015198370,1658072.543,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,p1024,1024,MODMUL,25200,0.017292426,1457285.404,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,p1024,1024,MODMUL,25200,0.005239217,4809879.043,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,p1024,1024,MODMUL,25200,0.007062731,3568024.880,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,p1024,1024,MODMUL,25200,0.002311303,10902940.946,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,p1024,1024,MODMUL,25200,0.004175763,6034825.259,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,p1024,1024,MODMUL,25200,0.001754515,14362943.692,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,p1024,1024,MODMUL,25200,0.003576445,7046103.072,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,p1024,1024,MODMUL,25200,0.002934496,8587505.299,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,p1024,1024,MODMUL,25200,0.004898891,5144021.373,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,p1024,1024,MODMUL,25200,0.001654141,15234493.307,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,p1024,1024,MODMUL,25200,0.003538237,7122191.115,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,p1024,1024,MODMUL,25200,0.636462545,39593.846,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,CPU,w8,p1024,1024,MODMUL,25200,0.708873219,35549.375,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,p1024,1024,MODEXP,195,0.075400048,2586.205,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,p1024,1024,MODEXP,195,0.009435278,20667.118,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,p1024,1024,MODEXP,195,0.006715797,29036.018,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,p1024,1024,MODEXP,50000,2.689949274,18587.711,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,p1024,1024,MODEXP,25200,3.536048670,7126.599,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,p1024,1024,MODEXP,25200,3.546370271,7105.857,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,p1024,1024,MODEXP,25200,0.500413621,50358.341,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,p1024,1024,MODEXP,25200,0.502493710,50149.881,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,p1024,1024,MODEXP,25200,0.258515252,97479.742,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,p1024,1024,MODEXP,25200,0.259967876,96935.054,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,p1024,1024,MODEXP,25200,0.154121407,163507.461,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,p1024,1024,MODEXP,25200,0.155456913,162102.794,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,p1024,1024,MODEXP,25200,0.259708620,97031.820,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,p1024,1024,MODEXP,25200,0.262178660,96117.663,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,p1024,1024,MODEXP,25200,0.154313031,163304.420,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,p1024,1024,MODEXP,25200,0.155737785,161810.443,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,p1024,1024,EXPONENTIATION,195,0.012475711,15630.372,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,p1024,1024,EXPONENTIATION,195,0.002864542,68073.709,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,p1024,1024,EXPONENTIATION,195,0.019282515,10112.789,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,p1024,1024,EXPONENTIATION,25200,3.367435705,7483.439,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,p1024,1024,EXPONENTIATION,25200,3.372107468,7473.071,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,p1024,1024,EXPONENTIATION,25200,0.715951234,35197.928,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,p1024,1024,EXPONENTIATION,25200,0.719651489,35016.950,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,p1024,1024,EXPONENTIATION,25200,0.218137193,115523.628,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,p1024,1024,EXPONENTIATION,25200,0.219808465,114645.266,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,p1024,1024,EXPONENTIATION,25200,0.214936263,117244.059,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,p1024,1024,EXPONENTIATION,25200,0.216713662,116282.470,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,p1024,1024,EXPONENTIATION,25200,0.216830071,116220.042,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,p1024,1024,EXPONENTIATION,25200,0.218971749,115083.339,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,p1024,1024,EXPONENTIATION,25200,0.213207596,118194.663,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,p1024,1024,EXPONENTIATION,25200,0.215287540,117052.757,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,p1024,1024,DIVIDE,1562,0.000140752,11097554.704,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,p1024,1024,DIVIDE,1562,0.000012903,121060537.882,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,p1024,1024,DIVIDE,1562,0.000045056,34667654.408,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,p1024,1024,DIVIDE,50000,0.000135424,369210775.047,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,p1024,1024,DIVIDE,25200,0.094915204,265500.140,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,p1024,1024,DIVIDE,25200,0.097523763,258398.561,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,p1024,1024,DIVIDE,25200,0.013891925,1814003.460,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,p1024,1024,DIVIDE,25200,0.016163608,1559057.854,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,p1024,1024,DIVIDE,25200,0.001978693,12735679.523,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,p1024,1024,DIVIDE,25200,0.004342116,5803622.002,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,p1024,1024,DIVIDE,25200,0.001845650,13653726.355,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,p1024,1024,DIVIDE,25200,0.004106570,6136508.074,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,p1024,1024,DIVIDE,25200,0.001963552,12833884.714,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,p1024,1024,DIVIDE,25200,0.004318020,5836008.189,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,p1024,1024,DIVIDE,25200,0.001872158,13460402.394,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,p1024,1024,DIVIDE,25200,0.004216125,5977052.385,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,p1024,1024,ISQRT,390,0.000186453,2091676.877,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,p1024,1024,ISQRT,390,0.000013775,28312438.987,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,p1024,1024,ISQRT,25200,1.299885815,19386.318,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,p1024,1024,ISQRT,25200,1.303116788,19338.251,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,p1024,1024,ISQRT,25200,0.222037316,113494.436,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,p1024,1024,ISQRT,25200,0.224052344,112473.717,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,p1024,1024,ISQRT,25200,0.037396365,673862.286,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,p1024,1024,ISQRT,25200,0.039392552,639714.838,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,p1024,1024,ISQRT,25200,0.035136987,717192.968,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,p1024,1024,ISQRT,25200,0.036985584,681346.548,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,p1024,1024,ISQRT,25200,0.037508116,671854.592,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,p1024,1024,ISQRT,25200,0.039331332,640710.566,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,p1024,1024,ISQRT,25200,0.034890032,722269.329,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,p1024,1024,ISQRT,25200,0.037052441,680117.134,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,p1024,1024,MODMUL_R2,12500,0.009835645,1270887.674,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,p1024,1024,MODMUL_R2,12500,0.000703602,17765725.242,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,p1024,1024,MODMUL_R2,12500,0.002155618,5798801.047,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,p1024,1024,MODMUL_R2,50000,0.000206848,241723391.089,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,p1024,1024,MODMUL_R2,25200,0.003829862,6579871.562,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,p1024,1024,MODMUL_R2,25200,0.005631275,4475007.873,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,p1024,1024,MODMUL_R2,25200,0.000619996,40645422.277,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,p1024,1024,MODMUL_R2,25200,0.002412242,10446713.058,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,p1024,1024,MODMUL_R2,25200,0.000358960,70202806.317,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,p1024,1024,MODMUL_R2,25200,0.002306682,10924782.872,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,p1024,1024,MODMUL_R2,25200,0.000274250,91886966.861,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,p1024,1024,MODMUL_R2,25200,0.002117186,11902591.465,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,p1024,1024,MODMUL_R2,25200,0.000295945,85150955.189,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,p1024,1024,MODMUL_R2,25200,0.002433540,10355284.857,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,p1024,1024,MODMUL_R2,25200,0.000216841,116214185.015,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,p1024,1024,MODMUL_R2,25200,0.002042538,12337591.811,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,p2048,2048,ADD,6250,0.000322725,19366313.594,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,p2048,2048,ADD,6250,0.000021921,285114705.894,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,p2048,2048,ADD,6250,0.000022472,278122541.372,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,p2048,2048,ADD,50000,0.000136896,365240766.713,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,p2048,2048,ADD,25200,0.000898605,28043468.007,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,p2048,2048,ADD,25200,0.004209873,5985928.785,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,p2048,2048,ADD,25200,0.000445192,56604793.123,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,p2048,2048,ADD,25200,0.003891608,6475472.336,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32,p2048,2048,ADD,25200,0.000227396,110819892.226,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32,p2048,2048,ADD,25200,0.003533574,7131589.725,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,p2048,2048,ADD,25200,0.000231152,109019176.870,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,p2048,2048,ADD,25200,0.003604410,6991435.467,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,p2048,2048,ADD,25200,0.000229236,109930377.212,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,p2048,2048,ADD,25200,0.003522792,7153416.916,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,p2048,2048,ADD,25200,0.000092773,271630755.287,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,p2048,2048,ADD,25200,0.003641799,6919657.006,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,p2048,2048,ADD,25200,0.000091980,273972579.775,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,p2048,2048,ADD,25200,0.003367047,7484303.010,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,p2048,2048,SUBTRACT,6250,0.000220823,28303171.735,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,p2048,2048,SUBTRACT,6250,0.000027896,224046447.479,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,p2048,2048,SUBTRACT,6250,0.000023725,263432642.155,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,p2048,2048,SUBTRACT,50000,0.000137088,364729225.023,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,p2048,2048,SUBTRACT,25200,0.000891066,28280733.071,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,p2048,2048,SUBTRACT,25200,0.004212464,5982246.967,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,p2048,2048,SUBTRACT,25200,0.000448751,56155863.546,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,p2048,2048,SUBTRACT,25200,0.003182738,7917711.101,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32,p2048,2048,SUBTRACT,25200,0.000227426,110805270.912,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32,p2048,2048,SUBTRACT,25200,0.003655632,6893472.816,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,p2048,2048,SUBTRACT,25200,0.000231392,108906106.692,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,p2048,2048,SUBTRACT,25200,0.003606911,6986587.708,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,p2048,2048,SUBTRACT,25200,0.000231332,108934345.390,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,p2048,2048,SUBTRACT,25200,0.003537710,7123252.049,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,p2048,2048,SUBTRACT,25200,0.000092620,272079461.712,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,p2048,2048,SUBTRACT,25200,0.003976692,6336925.296,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,p2048,2048,SUBTRACT,25200,0.000092162,273431584.070,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,p2048,2048,SUBTRACT,25200,0.003408340,7393628.549,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,p2048,2048,ADDMOD,6250,0.000634440,9851204.763,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,p2048,2048,ADDMOD,6250,0.000070945,88096827.307,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,p2048,2048,ADDMOD,6250,0.000275393,22694853.557,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,p2048,2048,ADDMOD,50000,0.000137088,364729225.023,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,p2048,2048,ADDMOD,25200,0.001493548,16872574.738,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,p2048,2048,ADDMOD,25200,0.004803657,5246003.193,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,p2048,2048,ADDMOD,25200,0.000577222,43657379.437,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,p2048,2048,ADDMOD,25200,0.003947859,6383206.712,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32,p2048,2048,ADDMOD,25200,0.000281765,89436238.615,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32,p2048,2048,ADDMOD,25200,0.003772287,6680297.674,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,p2048,2048,ADDMOD,25200,0.000250359,100655455.696,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,p2048,2048,ADDMOD,25200,0.003597339,7005178.043,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,p2048,2048,ADDMOD,25200,0.000240236,104896840.076,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,p2048,2048,ADDMOD,25200,0.003532924,7132901.771,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,p2048,2048,ADDMOD,25200,0.000088864,283579475.635,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,p2048,2048,ADDMOD,25200,0.003702726,6805796.627,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,p2048,2048,ADDMOD,25200,0.000085912,293323385.056,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,p2048,2048,ADDMOD,25200,0.003404184,7402655.060,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,p2048,2048,SUBTRACTMOD,6250,0.000511736,12213328.799,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,p2048,2048,SUBTRACTMOD,6250,0.000061647,101383448.889,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,p2048,2048,SUBTRACTMOD,6250,0.000373996,16711421.492,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,p2048,2048,SUBTRACTMOD,50000,0.000137120,364644107.351,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,p2048,2048,SUBTRACTMOD,25200,0.001557166,16183245.824,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,p2048,2048,SUBTRACTMOD,25200,0.004856156,5189289.655,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,p2048,2048,SUBTRACTMOD,25200,0.000809279,31138828.184,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,p2048,2048,SUBTRACTMOD,25200,0.003826906,6584954.003,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32,p2048,2048,SUBTRACTMOD,25200,0.000301613,83550773.464,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32,p2048,2048,SUBTRACTMOD,25200,0.003736726,6743871.492,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,p2048,2048,SUBTRACTMOD,25200,0.000243184,103625246.133,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,p2048,2048,SUBTRACTMOD,25200,0.003089197,8157459.742,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,p2048,2048,SUBTRACTMOD,25200,0.000245820,102514033.271,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,p2048,2048,SUBTRACTMOD,25200,0.003512835,7173693.074,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,p2048,2048,SUBTRACTMOD,25200,0.000092936,271154397.541,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,p2048,2048,SUBTRACTMOD,25200,0.003770542,6683389.317,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,p2048,2048,SUBTRACTMOD,25200,0.000089124,282752058.067,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,p2048,2048,SUBTRACTMOD,25200,0.003408240,7393845.548,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.003945074,1584254.385,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.000380786,16413418.353,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.000480891,12996721.521,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,25200,0.274186929,91908.101,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,25200,0.280221990,89928.703,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,p2048,2048,MULTIPLYOPERANDSCANNING,25200,0.073543150,342655.978,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,p2048,2048,MULTIPLYOPERANDSCANNING,25200,0.077358970,325754.079,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32,p2048,2048,MULTIPLYOPERANDSCANNING,25200,0.015929718,1581948.908,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32,p2048,2048,MULTIPLYOPERANDSCANNING,25200,0.020476401,1230685.021,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,p2048,2048,MULTIPLYOPERANDSCANNING,25200,0.001178082,21390701.313,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,p2048,2048,MULTIPLYOPERANDSCANNING,25200,0.005222564,4825216.110,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,p2048,2048,MULTIPLYOPERANDSCANNING,25200,0.000954689,26396029.947,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,p2048,2048,MULTIPLYOPERANDSCANNING,25200,0.005008798,5031147.188,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,p2048,2048,MULTIPLYOPERANDSCANNING,25200,0.000841834,29934642.686,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,p2048,2048,MULTIPLYOPERANDSCANNING,25200,0.005473845,4603710.921,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,p2048,2048,MULTIPLYOPERANDSCANNING,25200,0.000858428,29355986.078,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,p2048,2048,MULTIPLYOPERANDSCANNING,25200,0.005063159,4977129.878,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.005380868,1161522.637,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.000395032,15821503.263,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.000507128,12324304.994,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.000274944,181855214.153,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,25200,0.107495848,234427.659,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,25200,0.111757168,225488.892,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,p2048,2048,MULTIPLYPRODUCTSCANNING,25200,0.027155032,928004.798,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,p2048,2048,MULTIPLYPRODUCTSCANNING,25200,0.031643008,796384.465,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32,p2048,2048,MULTIPLYPRODUCTSCANNING,25200,0.006845121,3681454.274,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32,p2048,2048,MULTIPLYPRODUCTSCANNING,25200,0.011203245,2249348.292,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,p2048,2048,MULTIPLYPRODUCTSCANNING,25200,0.009099131,2769495.237,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,p2048,2048,MULTIPLYPRODUCTSCANNING,25200,0.013326559,1890960.748,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,p2048,2048,MULTIPLYPRODUCTSCANNING,25200,0.009146831,2755052.543,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,p2048,2048,MULTIPLYPRODUCTSCANNING,25200,0.013295868,1895325.675,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,p2048,2048,MULTIPLYPRODUCTSCANNING,25200,0.003343387,7537266.822,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,p2048,2048,MULTIPLYPRODUCTSCANNING,25200,0.008150695,3091760.895,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,p2048,2048,MULTIPLYPRODUCTSCANNING,25200,0.003226980,7809159.064,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,p2048,2048,MULTIPLYPRODUCTSCANNING,25200,0.007432034,3390727.219,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.034645507,180398.572,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.003816645,1637563.877,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.003695167,1691398.513,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000522080,95770763.101,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,25200,0.091536440,275300.197,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,25200,0.094982798,265311.199,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,p2048,2048,MONTGOMERYMULTIPLICATION,25200,0.002267849,11111850.821,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,p2048,2048,MONTGOMERYMULTIPLICATION,25200,0.005720300,4405363.364,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32,p2048,2048,MONTGOMERYMULTIPLICATION,25200,0.000411262,61274806.576,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32,p2048,2048,MONTGOMERYMULTIPLICATION,25200,0.003822225,6593018.464,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,p2048,2048,MONTGOMERYMULTIPLICATION,25200,0.000696500,36180905.692,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,p2048,2048,MONTGOMERYMULTIPLICATION,25200,0.004007195,6288688.258,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,p2048,2048,MONTGOMERYMULTIPLICATION,25200,0.000564691,44626176.129,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,p2048,2048,MONTGOMERYMULTIPLICATION,25200,0.003872509,6507409.024,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,p2048,2048,MONTGOMERYMULTIPLICATION,25200,0.000750016,33599283.909,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,p2048,2048,MONTGOMERYMULTIPLICATION,25200,0.004291665,5871846.923,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,p2048,2048,MONTGOMERYMULTIPLICATION,25200,0.000601107,41922653.216,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,p2048,2048,MONTGOMERYMULTIPLICATION,25200,0.003170440,7948423.556,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,p2048,2048,COMPARE,6250,0.000108673,57511796.124,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,p2048,2048,COMPARE,6250,0.000169417,36891298.537,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,p2048,2048,COMPARE,6250,0.000154912,40345486.503,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,p2048,2048,COMPARE,50000,0.000137184,364473991.136,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,p2048,2048,COMPARE,25200,0.000488257,51612162.578,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,p2048,2048,COMPARE,25200,0.003808513,6616755.686,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,p2048,2048,COMPARE,25200,0.000247496,101819825.082,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,p2048,2048,COMPARE,25200,0.004090590,6160480.525,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,p2048,2048,COMPARE,25200,0.000127344,197889160.473,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,p2048,2048,COMPARE,25200,0.003429753,7347467.863,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,p2048,2048,COMPARE,25200,0.000131259,191986873.434,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,p2048,2048,COMPARE,25200,0.003441155,7323122.582,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,p2048,2048,COMPARE,25200,0.000045380,555310683.615,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,p2048,2048,COMPARE,25200,0.003373583,7469802.893,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,p2048,2048,COMPARE,25200,0.000044217,569916754.783,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,p2048,2048,COMPARE,25200,0.002593258,9717505.976,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,p2048,2048,REDUCE,781,0.000052895,14765110.842,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,p2048,2048,REDUCE,781,0.000018815,41508632.254,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,p2048,2048,REDUCE,781,0.000040296,19381602.538,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,p2048,2048,REDUCE,50000,0.000137216,364388992.537,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,p2048,2048,REDUCE,25200,0.783238189,32174.121,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,p2048,2048,REDUCE,25200,0.787070225,32017.473,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,p2048,2048,REDUCE,25200,0.010642840,2367789.055,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,p2048,2048,REDUCE,25200,0.015522453,1623454.746,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,p2048,2048,REDUCE,25200,0.001780236,14155426.629,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,p2048,2048,REDUCE,25200,0.005081546,4959120.710,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,p2048,2048,REDUCE,25200,0.001603006,15720465.108,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,p2048,2048,REDUCE,25200,0.004935570,5105793.249,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,p2048,2048,REDUCE,25200,0.001935118,13022461.717,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,p2048,2048,REDUCE,25200,0.006077920,4146155.264,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,p2048,2048,REDUCE,25200,0.002012181,12523724.205,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,p2048,2048,REDUCE,25200,0.005768942,4368218.639,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,p2048,2048,MODMUL,390,0.000899335,433653.996,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,p2048,2048,MODMUL,390,0.000066901,5829538.075,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,p2048,2048,MODMUL,390,0.000212675,1833781.105,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,p2048,2048,MODMUL,50000,0.002570848,19448835.559,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,p2048,2048,MODMUL,25200,1.302109756,19353.207,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,p2048,2048,MODMUL,25200,1.305643464,19300.828,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,p2048,2048,MODMUL,25200,0.023941614,1052560.617,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,p2048,2048,MODMUL,25200,0.028518523,883636.225,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,p2048,2048,MODMUL,25200,0.009324723,2702493.152,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,p2048,2048,MODMUL,25200,0.012688533,1986045.196,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,p2048,2048,MODMUL,25200,0.007082201,3558215.863,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,p2048,2048,MODMUL,25200,0.010398469,2423433.682,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,p2048,2048,MODMUL,25200,0.011852026,2126218.760,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,p2048,2048,MODMUL,25200,0.013642410,1847180.961,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,p2048,2048,MODMUL,25200,0.007264280,3469029.281,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,p2048,2048,MODMUL,25200,0.011395977,2211306.670,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,p2048,2048,MODEXP,97,0.261496785,370.941,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,p2048,2048,MODEXP,97,0.030565020,3173.562,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,p2048,2048,MODEXP,97,0.012974926,7475.958,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,p2048,2048,MODEXP,50000,3.618140221,13819.254,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,p2048,2048,MODEXP,25200,166.915239916,150.975,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,p2048,2048,MODEXP,25200,166.921411666,150.969,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,p2048,2048,MODEXP,25200,24.569531317,1025.661,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,p2048,2048,MODEXP,25200,24.745639354,1018.361,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,p2048,2048,MODEXP,25200,2.443615774,10312.587,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,p2048,2048,MODEXP,25200,2.447266133,10297.205,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,p2048,2048,MODEXP,25200,5.042539623,4997.482,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,p2048,2048,MODEXP,25200,5.057097644,4983.095,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,p2048,2048,MODEXP,25200,2.386695814,10558.530,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,p2048,2048,MODEXP,25200,2.429002734,10374.628,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,p2048,2048,MODEXP,25200,4.982085668,5058.123,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,p2048,2048,MODEXP,25200,5.004292110,5035.677,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,p2048,2048,EXPONENTIATION,97,0.042759016,2268.527,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,p2048,2048,EXPONENTIATION,97,0.003452520,28095.420,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,p2048,2048,EXPONENTIATION,97,0.039794993,2437.493,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,p2048,2048,EXPONENTIATION,25200,40.154482551,627.576,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,p2048,2048,EXPONENTIATION,25200,40.229617843,626.404,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,p2048,2048,EXPONENTIATION,25200,7.922459739,3180.830,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,p2048,2048,EXPONENTIATION,25200,7.933568211,3176.377,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,p2048,2048,EXPONENTIATION,25200,2.051151400,12285.783,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,p2048,2048,EXPONENTIATION,25200,2.061092078,12226.528,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,p2048,2048,EXPONENTIATION,25200,2.012448529,12522.059,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,p2048,2048,EXPONENTIATION,25200,2.016521815,12496.765,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,p2048,2048,EXPONENTIATION,25200,2.026036769,12438.076,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,p2048,2048,EXPONENTIATION,25200,2.044285358,12327.046,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,p2048,2048,EXPONENTIATION,25200,1.975529807,12756.072,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,p2048,2048,EXPONENTIATION,25200,1.975061778,12759.095,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,p2048,2048,DIVIDE,781,0.000142538,5479240.619,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,p2048,2048,DIVIDE,781,0.000008282,94302869.207,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,p2048,2048,DIVIDE,781,0.000039496,19774045.512,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,p2048,2048,DIVIDE,50000,0.000198656,251691365.979,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,p2048,2048,DIVIDE,25200,1.871957548,13461.844,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,p2048,2048,DIVIDE,25200,1.876919798,13426.253,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,p2048,2048,DIVIDE,25200,0.387007945,65114.942,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,p2048,2048,DIVIDE,25200,0.390832994,64477.668,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,p2048,2048,DIVIDE,25200,0.028501575,884161.665,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,p2048,2048,DIVIDE,25200,0.032685429,770985.750,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,p2048,2048,DIVIDE,25200,0.025225862,998974.782,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,p2048,2048,DIVIDE,25200,0.029573180,852123.444,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,p2048,2048,DIVIDE,25200,0.028743281,876726.634,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,p2048,2048,DIVIDE,25200,0.032820253,767818.579,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,p2048,2048,DIVIDE,25200,0.024982536,1008704.640,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,p2048,2048,DIVIDE,25200,0.030436197,827961.522,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,p2048,2048,ISQRT,195,0.000175320,1112254.888,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,p2048,2048,ISQRT,195,0.000013616,14321442.599,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,p2048,2048,ISQRT,25200,24.358499084,1034.547,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,p2048,2048,ISQRT,25200,24.382106123,1033.545,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,p2048,2048,ISQRT,25200,8.883979022,2836.567,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,p2048,2048,ISQRT,25200,8.894815092,2833.111,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,p2048,2048,ISQRT,25200,0.091787798,274546.296,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,p2048,2048,ISQRT,25200,0.095053688,265113.333,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,p2048,2048,ISQRT,25200,0.078369075,321555.409,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,p2048,2048,ISQRT,25200,0.082159931,306718.855,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,p2048,2048,ISQRT,25200,0.091738568,274693.627,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,p2048,2048,ISQRT,25200,0.095588684,263629.532,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,p2048,2048,ISQRT,25200,0.078594099,320634.759,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,p2048,2048,ISQRT,25200,0.081870012,307805.012,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-1t,p2048,2048,MODMUL_R2,6250,0.014742030,423957.894,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,gmp-nt,p2048,2048,MODMUL_R2,6250,0.001037000,6027001.081,0
library,Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz,host-cpu,openssl-nt,p2048,2048,MODMUL_R2,6250,0.003201297,1952333.701,0
library,NVIDIA GeForce RTX 2070,gpu,cgbn,p2048,2048,MODMUL_R2,50000,0.000751616,66523331.063,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w8,p2048,2048,MODMUL_R2,25200,0.066506062,378912.828,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w8,p2048,2048,MODMUL_R2,25200,0.069759041,361243.498,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w16,p2048,2048,MODMUL_R2,25200,0.002375393,10608770.844,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w16,p2048,2048,MODMUL_R2,25200,0.005585192,4511930.831,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-opt,p2048,2048,MODMUL_R2,25200,0.001218774,20676515.795,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-opt,p2048,2048,MODMUL_R2,25200,0.004541119,5549293.017,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-o64,p2048,2048,MODMUL_R2,25200,0.000910178,27686891.969,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-o64,p2048,2048,MODMUL_R2,25200,0.005110605,4930923.043,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il,p2048,2048,MODMUL_R2,25200,0.001074958,23442776.138,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il,p2048,2048,MODMUL_R2,25200,0.004459342,5651057.952,0
opencl-kernel,NVIDIA GeForce RTX 2070,GPU,w32-il64,p2048,2048,MODMUL_R2,25200,0.000819781,30739916.774,0
opencl-e2e,NVIDIA GeForce RTX 2070,GPU,w32-il64,p2048,2048,MODMUL_R2,25200,0.004478516,5626863.900,0
```
