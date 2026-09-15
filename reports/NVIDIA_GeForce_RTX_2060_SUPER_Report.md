# MPA-OpenCL benchmark report - NVIDIA GeForce RTX 2060 SUPER


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

### Device 0 - NVIDIA GeForce RTX 2060 SUPER (GPU)

| Property | Value |
|---|---|
| Model | NVIDIA GeForce RTX 2060 SUPER |
| Type | GPU |
| Vendor | NVIDIA Corporation |
| Device memory | 7.60 GiB |
| Max single allocation | 1.90 GiB |
| Local memory | 48 KiB |
| Global cache | 1088 KiB |
| Compute units | 34 |
| Max clock | 1650 MHz |
| Max work-group size | 1024 |
| OpenCL version | OpenCL 3.0 CUDA |
| Driver | 595.58.03 |

### Device 1 - cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz (CPU)

| Property | Value |
|---|---|
| Model | cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz |
| Type | CPU |
| Vendor | GenuineIntel |
| Device memory | 501.74 GiB |
| Max single allocation | 128.00 GiB |
| Local memory | 256 KiB |
| Global cache | 51200 KiB |
| Compute units | 80 |
| Max clock | 3600 MHz |
| Max work-group size | 4096 |
| OpenCL version | OpenCL 3.0 PoCL HSTR: cpu-x86_64-pc-linux-gnu-haswell |
| Driver | 5.0+debian |

### Host

| Property | Value |
|---|---|
| CPU | Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz |
| Logical cores | 80 |
| OpenMP threads used | 80 |
| RAM | 503.7 GB |
| OS | Ubuntu 24.04.4 LTS |
| Kernel | 6.8.0-110-generic |
| Arch | x86_64 |
| GMP | 6.3.0 |
| OpenSSL | OpenSSL 3.0.13 30 Jan 2024 |

## 2. Method

- Workload auto-sized from the device and host: --min-items from 700 x compute units, --items from ten times that capped by host RAM. Either flag, given explicitly, overrides its half.
- Base workload 50000 items, scaled down per operator by its cost weight and by modulus size. Device rows honour --min-items (23800) so the GPU is not left idle; the CPU libraries keep the smaller count because a full-width MODEXP there costs minutes. Both counts appear in every row as dev/cpu, and throughput is per-second so they remain comparable.
- 5 timed repetitions, **minimum** reported. Two untimed warm-up launches precede them.
- `kernel` times `clEnqueueNDRangeKernel` + `clFinish` only. `e2e` adds the host->device operand writes and the device->host result read.
- Every OpenCL device runs the same kernels on the same operands, so GPU and CPU-OpenCL columns are directly comparable.
- CPU library baselines (GMP, OpenSSL) run those same operands, with temporaries preallocated outside the timed region, so the figure is the arithmetic and not marshalling. The generator is reseeded per modulus and operation so every backend sees identical inputs.
- OpenSSL rows time the nearest BN primitive, which is not always semantically identical (its Montgomery routine expects Montgomery-domain inputs); they measure comparable work, not identical results. Correctness is judged against GMP only.
- Every device cell is checked word-for-word against GMP before it is timed. A cell that mismatches is reported and excluded from the speedup tables.
- Total wall time 5648.0 s.

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
| [1] CPU | `mpaKernels_8bits.cl` (w8) | 56 | 56 | 0 | 0 |

**All configurations correct** - 541 configurations, 0 problems.

## 4. Throughput per device

Operations per second, higher is better. Kernel-only timings.

### Device 0 - NVIDIA GeForce RTX 2060 SUPER (GPU)

#### secp256k1 (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | 453.56 M | 1.07 G | 1.49 G | 1.61 G | 1.43 G | 2.38 G | 2.38 G | 34.66 M |
| SUBTRACT | 50000 / 50000 | 478.25 M | 1.07 G | 1.53 G | 1.65 G | 1.52 G | 2.39 G | 2.46 G | 55.58 M |
| ADDMOD | 50000 / 50000 | 315.22 M | 730.84 M | 1.21 G | 1.87 G | 2.00 G | 2.66 G | 2.67 G | 10.17 M |
| SUBTRACTMOD | 50000 / 50000 | 316.85 M | 748.27 M | 1.22 G | 1.83 G | 1.89 G | 2.72 G | 2.66 G | 14.19 M |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | 11.61 M | 54.72 M | 198.10 M | 911.28 M | 909.35 M | 1.86 G | 2.09 G | 20.08 M |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | 39.19 M | 166.56 M | 500.72 M | 566.38 M | 536.60 M | 1.02 G | 1.06 G | 28.06 M |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | 154.93 M | 777.02 M | 1.71 G | 1.57 G | 1.61 G | 1.68 G | 1.58 G | 3.61 M |
| COMPARE | 50000 / 50000 | 501.61 M | 1.14 G | - | 2.24 G | 2.01 G | 4.47 G | 4.59 G | 81.82 M |
| REDUCE | 23800 / 6250 | 80.36 M | 162.46 M | - | 448.06 M | 468.95 M | 471.53 M | 498.25 M | 28.43 M |
| MODMUL | 23800 / 3125 | 29.44 M | 67.10 M | - | 154.30 M | 184.04 M | 154.07 M | 185.66 M | 6.22 M |
| MODEXP | 23800 / 781 | 858.84 k | 4.55 M | - | 6.54 M | 12.04 M | 6.54 M | 11.92 M | 60.54 k |
| EXPONENTIATION | 23800 / 781 | 583.03 k | 1.78 M | - | 34.99 M | 45.08 M | 35.03 M | 45.09 M | 171.15 k |
| DIVIDE | 23800 / 6250 | 54.02 M | 82.55 M | - | 196.31 M | 208.54 M | 202.78 M | 207.26 M | 11.87 M |
| ISQRT | 23800 / 1562 | 4.40 M | 5.53 M | - | 23.60 M | 27.68 M | 23.90 M | 27.53 M | 6.87 M |
| MODMUL_R2 | 50000 / 50000 | 165.09 M | 731.04 M | - | 786.54 M | 1.17 G | 952.73 M | 1.17 G | 6.24 M |

#### rsa256(composite) (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | 602.64 M | 1.06 G | 1.54 G | 1.36 G | 1.51 G | 2.37 G | 2.42 G | 23.43 M |
| SUBTRACT | 50000 / 50000 | 604.33 M | 1.07 G | 1.50 G | 1.38 G | 1.51 G | 2.43 G | 2.40 G | 42.05 M |
| ADDMOD | 50000 / 50000 | 439.72 M | 789.32 M | 1.32 G | 1.62 G | 2.05 G | 2.66 G | 2.39 G | 13.01 M |
| SUBTRACTMOD | 50000 / 50000 | 400.24 M | 748.61 M | 1.20 G | 1.65 G | 1.91 G | 2.73 G | 2.33 G | 21.50 M |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | 14.86 M | 55.51 M | 199.46 M | 732.96 M | 920.76 M | 1.88 G | 1.92 G | 41.43 M |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | 49.71 M | 167.26 M | 499.64 M | 465.69 M | 442.06 M | 901.87 M | 946.29 M | 40.86 M |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | 196.53 M | 780.73 M | 1.77 G | 1.46 G | 1.49 G | 1.57 G | 1.68 G | 3.61 M |
| COMPARE | 50000 / 50000 | 653.33 M | 1.13 G | - | 1.98 G | 1.74 G | 3.82 G | 4.13 G | 81.71 M |
| REDUCE | 23800 / 6250 | 101.73 M | 161.12 M | - | 356.22 M | 379.91 M | 384.89 M | 399.79 M | 14.93 M |
| MODMUL | 23800 / 3125 | 37.60 M | 67.25 M | - | 123.12 M | 147.62 M | 124.63 M | 147.67 M | 3.55 M |
| MODEXP | 23800 / 781 | 853.19 k | 4.57 M | - | 5.22 M | 9.63 M | 5.24 M | 9.39 M | 63.16 k |
| EXPONENTIATION | 23800 / 781 | 583.36 k | 1.77 M | - | 28.05 M | 36.31 M | 27.86 M | 36.44 M | 114.79 k |
| DIVIDE | 23800 / 6250 | 53.79 M | 80.38 M | - | 154.47 M | 162.48 M | 161.78 M | 165.49 M | 11.44 M |
| ISQRT | 23800 / 1562 | 4.41 M | 5.48 M | - | 18.72 M | 22.12 M | 18.96 M | 22.00 M | 6.95 M |
| MODMUL_R2 | 50000 / 50000 | 165.91 M | 732.09 M | - | 714.09 M | 1.05 G | 856.31 M | 1.12 G | 6.25 M |

#### brainpoolP512r1 (512-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 25000 / 25000 | 221.04 M | 393.39 M | 423.48 M | 513.59 M | 521.43 M | 748.26 M | 751.99 M | 32.43 M |
| SUBTRACT | 25000 / 25000 | 219.95 M | 396.67 M | 414.49 M | 519.52 M | 525.31 M | 770.80 M | 760.15 M | 36.88 M |
| ADDMOD | 25000 / 25000 | 162.63 M | 311.69 M | 345.08 M | 430.59 M | 436.19 M | 935.59 M | 953.84 M | 11.68 M |
| SUBTRACTMOD | 25000 / 25000 | 142.27 M | 286.95 M | 323.22 M | 428.29 M | 427.92 M | 920.23 M | 939.14 M | 12.92 M |
| MULTIPLYOPERANDSCANNING | 25000 / 25000 | 2.07 M | 7.90 M | 31.37 M | 178.35 M | 180.06 M | 340.08 M | 367.67 M | 12.33 M |
| MULTIPLYPRODUCTSCANNING | 25000 / 25000 | 7.21 M | 26.69 M | 65.54 M | 76.37 M | 76.96 M | 163.85 M | 174.28 M | 11.25 M |
| MONTGOMERYMULTIPLICATION | 25000 / 25000 | 51.37 M | 169.87 M | 467.44 M | 340.54 M | 404.17 M | 408.34 M | 563.84 M | 1.52 M |
| COMPARE | 25000 / 25000 | 262.78 M | 453.75 M | - | 660.55 M | 652.06 M | 1.85 G | 1.91 G | 80.58 M |
| REDUCE | 23800 / 3125 | 34.78 M | 43.81 M | - | 131.00 M | 137.39 M | 124.91 M | 133.03 M | 10.42 M |
| MODMUL | 23800 / 1562 | 12.37 M | 17.27 M | - | 35.86 M | 43.95 M | 29.85 M | 37.10 M | 3.17 M |
| MODEXP | 23800 / 390 | 63.41 k | 642.72 k | - | 820.34 k | 1.37 M | 775.66 k | 1.30 M | 11.64 k |
| EXPONENTIATION | 23800 / 390 | 70.26 k | 276.82 k | - | 917.12 k | 967.36 k | 885.72 k | 930.03 k | 47.37 k |
| DIVIDE | 23800 / 3125 | 17.85 M | 19.02 M | - | 60.48 M | 62.72 M | 62.61 M | 65.35 M | 7.61 M |
| ISQRT | 23800 / 781 | 741.15 k | 843.59 k | - | 4.25 M | 4.61 M | 3.57 M | 3.83 M | 3.99 M |
| MODMUL_R2 | 25000 / 25000 | 35.73 M | 162.58 M | - | 267.42 M | 260.18 M | 280.29 M | 307.03 M | 4.88 M |

#### p1024 (1024-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 23800 / 12500 | 54.89 M | 102.41 M | 82.81 M | 202.07 M | 202.21 M | 433.21 M | 461.22 M | 25.99 M |
| SUBTRACT | 23800 / 12500 | 54.11 M | 102.23 M | 81.75 M | 203.46 M | 198.60 M | 467.04 M | 452.63 M | 29.34 M |
| ADDMOD | 23800 / 12500 | 40.72 M | 78.02 M | 79.45 M | 190.41 M | 192.80 M | 589.55 M | 570.27 M | 8.17 M |
| SUBTRACTMOD | 23800 / 12500 | 40.53 M | 77.85 M | 79.01 M | 194.75 M | 192.66 M | 560.55 M | 571.83 M | 11.49 M |
| MULTIPLYOPERANDSCANNING | 23800 / 12500 | 404.63 k | 1.58 M | 3.35 M | 84.44 M | 84.87 M | 133.25 M | 137.29 M | 2.89 M |
| MULTIPLYPRODUCTSCANNING | 23800 / 12500 | 931.70 k | 3.64 M | 11.15 M | 14.11 M | 14.23 M | 47.41 M | 50.06 M | 3.49 M |
| MONTGOMERYMULTIPLICATION | 23800 / 12500 | 8.84 M | 43.98 M | 133.60 M | 102.93 M | 130.65 M | 145.38 M | 200.66 M | 553.45 k |
| COMPARE | 23800 / 12500 | 101.34 M | 191.53 M | - | 342.73 M | 332.12 M | 1.62 G | 1.51 G | 58.68 M |
| REDUCE | 23800 / 1562 | 6.58 M | 11.95 M | - | 47.66 M | 47.26 M | 52.70 M | 54.30 M | 22.91 M |
| MODMUL | 23800 / 781 | 1.56 M | 4.50 M | - | 10.40 M | 13.57 M | 10.64 M | 14.19 M | 1.03 M |
| MODEXP | 23800 / 195 | 6.71 k | 47.69 k | - | 96.53 k | 161.56 k | 95.63 k | 162.22 k | 2.79 k |
| EXPONENTIATION | 23800 / 195 | 6.83 k | 33.51 k | - | 114.99 k | 115.16 k | 114.55 k | 112.66 k | 12.61 k |
| DIVIDE | 23800 / 1562 | 269.43 k | 1.81 M | - | 11.91 M | 12.76 M | 11.98 M | 12.61 M | 11.51 M |
| ISQRT | 23800 / 390 | 18.43 k | 103.25 k | - | 627.73 k | 671.76 k | 627.34 k | 674.59 k | 1.22 M |
| MODMUL_R2 | 23800 / 12500 | 6.18 M | 38.69 M | - | 66.50 M | 88.05 M | 81.12 M | 111.62 M | 1.41 M |

#### p2048 (2048-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 23800 / 6250 | 28.10 M | 55.89 M | 88.53 M | 108.12 M | 108.83 M | 256.14 M | 254.56 M | 19.39 M |
| SUBTRACT | 23800 / 6250 | 28.13 M | 55.76 M | 87.43 M | 108.39 M | 108.31 M | 255.31 M | 254.49 M | 22.01 M |
| ADDMOD | 23800 / 6250 | 17.31 M | 43.08 M | 70.63 M | 100.02 M | 102.13 M | 288.99 M | 289.70 M | 8.53 M |
| SUBTRACTMOD | 23800 / 6250 | 16.68 M | 32.04 M | 65.15 M | 101.43 M | 102.00 M | 270.06 M | 280.19 M | 8.78 M |
| MULTIPLYOPERANDSCANNING | 23800 / 6250 | 99.65 k | 391.45 k | 1.34 M | 27.56 M | 28.16 M | 32.45 M | 37.57 M | 1.11 M |
| MULTIPLYPRODUCTSCANNING | 23800 / 6250 | 232.00 k | 916.88 k | 3.64 M | 2.89 M | 2.89 M | 7.86 M | 8.16 M | 1.12 M |
| MONTGOMERYMULTIPLICATION | 23800 / 6250 | 274.80 k | 10.40 M | 58.84 M | 27.56 M | 43.00 M | 32.02 M | 40.97 M | 268.05 k |
| COMPARE | 23800 / 6250 | 50.98 M | 99.82 M | - | 154.87 M | 187.69 M | 589.91 M | 597.28 M | 85.54 M |
| REDUCE | 23800 / 781 | 30.72 k | 2.91 M | - | 10.40 M | 14.75 M | 12.64 M | 12.07 M | 15.96 M |
| MODMUL | 23800 / 390 | 18.58 k | 1.00 M | - | 2.01 M | 3.36 M | 2.07 M | 2.67 M | 397.91 k |
| MODEXP | 23800 / 97 | 136.5 | 917.7 | - | 10.17 k | 5.36 k | 10.19 k | 5.30 k | 397.8 |
| EXPONENTIATION | 23800 / 97 | 562.1 | 2.74 k | - | 11.98 k | 10.86 k | 12.07 k | 11.00 k | 2.17 k |
| DIVIDE | 23800 / 781 | 13.00 k | 63.80 k | - | 918.37 k | 1.02 M | 950.67 k | 1.01 M | 9.59 M |
| ISQRT | 23800 / 195 | 1.01 k | 2.68 k | - | 256.08 k | 299.27 k | 257.92 k | 298.83 k | 1.11 M |
| MODMUL_R2 | 23800 / 6250 | 362.64 k | 10.28 M | - | 19.76 M | 26.47 M | 22.14 M | 29.23 M | 404.80 k |

### Device 1 - cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz (CPU)

#### secp256k1 (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | 109.69 M | - | - | - | - | - | - | 34.66 M |
| SUBTRACT | 50000 / 50000 | 103.41 M | - | - | - | - | - | - | 55.58 M |
| ADDMOD | 50000 / 50000 | 87.25 M | - | - | - | - | - | - | 10.17 M |
| SUBTRACTMOD | 50000 / 50000 | 104.24 M | - | - | - | - | - | - | 14.19 M |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | 8.59 M | - | - | - | - | - | - | 20.08 M |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | 10.07 M | - | - | - | - | - | - | 28.06 M |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | 535.00 k | - | - | - | - | - | - | 3.61 M |
| COMPARE | 50000 / 50000 | 105.26 M | - | - | - | - | - | - | 81.82 M |
| REDUCE | 23800 / 6250 | 2.03 M | - | - | - | - | - | - | 28.43 M |
| MODMUL | 23800 / 3125 | 239.83 k | - | - | - | - | - | - | 6.22 M |
| MODEXP | 23800 / 781 | 5.64 k | - | - | - | - | - | - | 60.54 k |
| EXPONENTIATION | 23800 / 781 | 16.01 k | - | - | - | - | - | - | 171.15 k |
| DIVIDE | 23800 / 6250 | 324.46 k | - | - | - | - | - | - | 11.87 M |
| ISQRT | 23800 / 1562 | 64.40 k | - | - | - | - | - | - | 6.87 M |
| MODMUL_R2 | 50000 / 50000 | 566.64 k | - | - | - | - | - | - | 6.24 M |

#### rsa256(composite) (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | 100.10 M | - | - | - | - | - | - | 23.43 M |
| SUBTRACT | 50000 / 50000 | 114.38 M | - | - | - | - | - | - | 42.05 M |
| ADDMOD | 50000 / 50000 | 98.14 M | - | - | - | - | - | - | 13.01 M |
| SUBTRACTMOD | 50000 / 50000 | 92.57 M | - | - | - | - | - | - | 21.50 M |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | 8.56 M | - | - | - | - | - | - | 41.43 M |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | 9.94 M | - | - | - | - | - | - | 40.86 M |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | 509.52 k | - | - | - | - | - | - | 3.61 M |
| COMPARE | 50000 / 50000 | 108.77 M | - | - | - | - | - | - | 81.71 M |
| REDUCE | 23800 / 6250 | 2.01 M | - | - | - | - | - | - | 14.93 M |
| MODMUL | 23800 / 3125 | 240.61 k | - | - | - | - | - | - | 3.55 M |
| MODEXP | 23800 / 781 | 5.85 k | - | - | - | - | - | - | 63.16 k |
| EXPONENTIATION | 23800 / 781 | 14.73 k | - | - | - | - | - | - | 114.79 k |
| DIVIDE | 23800 / 6250 | 1.22 M | - | - | - | - | - | - | 11.44 M |
| ISQRT | 23800 / 1562 | 71.35 k | - | - | - | - | - | - | 6.95 M |
| MODMUL_R2 | 50000 / 50000 | 1.85 M | - | - | - | - | - | - | 6.25 M |

#### brainpoolP512r1 (512-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 25000 / 25000 | 39.95 M | - | - | - | - | - | - | 32.43 M |
| SUBTRACT | 25000 / 25000 | 43.38 M | - | - | - | - | - | - | 36.88 M |
| ADDMOD | 25000 / 25000 | 37.51 M | - | - | - | - | - | - | 11.68 M |
| SUBTRACTMOD | 25000 / 25000 | 34.77 M | - | - | - | - | - | - | 12.92 M |
| MULTIPLYOPERANDSCANNING | 25000 / 25000 | 1.10 M | - | - | - | - | - | - | 12.33 M |
| MULTIPLYPRODUCTSCANNING | 25000 / 25000 | 1.50 M | - | - | - | - | - | - | 11.25 M |
| MONTGOMERYMULTIPLICATION | 25000 / 25000 | 151.75 k | - | - | - | - | - | - | 1.52 M |
| COMPARE | 25000 / 25000 | 72.04 M | - | - | - | - | - | - | 80.58 M |
| REDUCE | 23800 / 3125 | 249.55 k | - | - | - | - | - | - | 10.42 M |
| MODMUL | 23800 / 1562 | 86.84 k | - | - | - | - | - | - | 3.17 M |
| MODEXP | 23800 / 390 | over budget | - | - | - | - | - | - | 11.64 k |
| EXPONENTIATION | 23800 / 390 | over budget | - | - | - | - | - | - | 47.37 k |
| DIVIDE | 23800 / 3125 | 203.95 k | - | - | - | - | - | - | 7.61 M |
| ISQRT | 23800 / 781 | 17.67 k | - | - | - | - | - | - | 3.99 M |
| MODMUL_R2 | 25000 / 25000 | 184.83 k | - | - | - | - | - | - | 4.88 M |

#### p1024 (1024-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 23800 / 12500 | 30.45 M | - | - | - | - | - | - | 25.99 M |
| SUBTRACT | 23800 / 12500 | 32.72 M | - | - | - | - | - | - | 29.34 M |
| ADDMOD | 23800 / 12500 | 26.81 M | - | - | - | - | - | - | 8.17 M |
| SUBTRACTMOD | 23800 / 12500 | 25.82 M | - | - | - | - | - | - | 11.49 M |
| MULTIPLYOPERANDSCANNING | 23800 / 12500 | 173.78 k | - | - | - | - | - | - | 2.89 M |
| MULTIPLYPRODUCTSCANNING | 23800 / 12500 | 242.69 k | - | - | - | - | - | - | 3.49 M |
| MONTGOMERYMULTIPLICATION | 23800 / 12500 | 61.88 k | - | - | - | - | - | - | 553.45 k |
| COMPARE | 23800 / 12500 | 41.51 M | - | - | - | - | - | - | 58.68 M |
| REDUCE | 23800 / 1562 | 84.47 k | - | - | - | - | - | - | 22.91 M |
| MODMUL | 23800 / 781 | 23.73 k | - | - | - | - | - | - | 1.03 M |
| MODEXP | 23800 / 195 | over budget | - | - | - | - | - | - | 2.79 k |
| EXPONENTIATION | 23800 / 195 | - | - | - | - | - | - | - | 12.61 k |
| DIVIDE | 23800 / 1562 | - | - | - | - | - | - | - | 11.51 M |
| ISQRT | 23800 / 390 | - | - | - | - | - | - | - | 1.22 M |
| MODMUL_R2 | 23800 / 12500 | - | - | - | - | - | - | - | 1.41 M |

#### p2048 (2048-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 23800 / 6250 | - | - | - | - | - | - | - | 19.39 M |
| SUBTRACT | 23800 / 6250 | - | - | - | - | - | - | - | 22.01 M |
| ADDMOD | 23800 / 6250 | - | - | - | - | - | - | - | 8.53 M |
| SUBTRACTMOD | 23800 / 6250 | - | - | - | - | - | - | - | 8.78 M |
| MULTIPLYOPERANDSCANNING | 23800 / 6250 | - | - | - | - | - | - | - | 1.11 M |
| MULTIPLYPRODUCTSCANNING | 23800 / 6250 | - | - | - | - | - | - | - | 1.12 M |
| MONTGOMERYMULTIPLICATION | 23800 / 6250 | - | - | - | - | - | - | - | 268.05 k |
| COMPARE | 23800 / 6250 | - | - | - | - | - | - | - | 85.54 M |
| REDUCE | 23800 / 781 | - | - | - | - | - | - | - | 15.96 M |
| MODMUL | 23800 / 390 | - | - | - | - | - | - | - | 397.91 k |
| MODEXP | 23800 / 97 | - | - | - | - | - | - | - | 397.8 |
| EXPONENTIATION | 23800 / 97 | - | - | - | - | - | - | - | 2.17 k |
| DIVIDE | 23800 / 781 | - | - | - | - | - | - | - | 9.59 M |
| ISQRT | 23800 / 195 | - | - | - | - | - | - | - | 1.11 M |
| MODMUL_R2 | 23800 / 6250 | - | - | - | - | - | - | - | 404.80 k |

## 5. Head to head

Best OpenCL GPU result against best OpenCL CPU result and the CPU libraries.
Ratios above 1.00x mean the GPU is faster than that baseline.

### secp256k1 (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w32-il64 | 2.38 G | w8 | 109.69 M | 34.66 M | 21.72x |
| SUBTRACT | w32-il64 | 2.46 G | w8 | 103.41 M | 55.58 M | 23.75x |
| ADDMOD | w32-il64 | 2.67 G | w8 | 87.25 M | 10.17 M | 30.61x |
| SUBTRACTMOD | w32-il | 2.72 G | w8 | 104.24 M | 14.19 M | 26.13x |
| MULTIPLYOPERANDSCANNING | w32-il64 | 2.09 G | w8 | 8.59 M | 20.08 M | 242.85x |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 1.06 G | w8 | 10.07 M | 28.06 M | 105.48x |
| MONTGOMERYMULTIPLICATION | w32 | 1.71 G | w8 | 535.00 k | 3.61 M | 3190.20x |
| COMPARE | w32-il64 | 4.59 G | w8 | 105.26 M | 81.82 M | 43.63x |
| REDUCE | w32-il64 | 130.84 M | w8 | 533.15 k | 28.43 M | 245.41x |
| MODMUL | w32-il64 | 24.38 M | w8 | 31.49 k | 6.22 M | 774.13x |
| MODEXP | w32-o64 | 395.10 k | w8 | 185.1 | 60.54 k | 2134.98x |
| EXPONENTIATION | w32-il64 | 1.48 M | w8 | 525.3 | 171.15 k | 2816.33x |
| DIVIDE | w32-o64 | 54.76 M | w8 | 85.20 k | 11.87 M | 642.74x |
| ISQRT | w32-o64 | 1.82 M | w8 | 4.23 k | 6.87 M | 429.76x |
| MODMUL_R2 | w32-il64 | 1.17 G | w8 | 566.64 k | 6.24 M | 2067.53x |

### rsa256(composite) (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w32-il64 | 2.42 G | w8 | 100.10 M | 23.43 M | 24.13x |
| SUBTRACT | w32-il | 2.43 G | w8 | 114.38 M | 42.05 M | 21.28x |
| ADDMOD | w32-il | 2.66 G | w8 | 98.14 M | 13.01 M | 27.05x |
| SUBTRACTMOD | w32-il | 2.73 G | w8 | 92.57 M | 21.50 M | 29.48x |
| MULTIPLYOPERANDSCANNING | w32-il64 | 1.92 G | w8 | 8.56 M | 41.43 M | 223.95x |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 946.29 M | w8 | 9.94 M | 40.86 M | 95.22x |
| MONTGOMERYMULTIPLICATION | w32 | 1.77 G | w8 | 509.52 k | 3.61 M | 3464.51x |
| COMPARE | w32-il64 | 4.13 G | w8 | 108.77 M | 81.71 M | 37.93x |
| REDUCE | w32-il64 | 104.99 M | w8 | 528.06 k | 14.93 M | 198.82x |
| MODMUL | w32-il64 | 19.39 M | w8 | 31.59 k | 3.55 M | 613.71x |
| MODEXP | w32-o64 | 315.90 k | w8 | 191.9 | 63.16 k | 1646.45x |
| EXPONENTIATION | w32-il64 | 1.20 M | w8 | 483.3 | 114.79 k | 2474.44x |
| DIVIDE | w32-il64 | 43.46 M | w8 | 321.21 k | 11.44 M | 135.30x |
| ISQRT | w32-o64 | 1.45 M | w8 | 4.68 k | 6.95 M | 310.10x |
| MODMUL_R2 | w32-il64 | 1.12 G | w8 | 1.85 M | 6.25 M | 606.63x |

### brainpoolP512r1 (512-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w32-il64 | 751.99 M | w8 | 39.95 M | 32.43 M | 18.83x |
| SUBTRACT | w32-il | 770.80 M | w8 | 43.38 M | 36.88 M | 17.77x |
| ADDMOD | w32-il64 | 953.84 M | w8 | 37.51 M | 11.68 M | 25.43x |
| SUBTRACTMOD | w32-il64 | 939.14 M | w8 | 34.77 M | 12.92 M | 27.01x |
| MULTIPLYOPERANDSCANNING | w32-il64 | 367.67 M | w8 | 1.10 M | 12.33 M | 333.81x |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 174.28 M | w8 | 1.50 M | 11.25 M | 116.33x |
| MONTGOMERYMULTIPLICATION | w32-il64 | 563.84 M | w8 | 151.75 k | 1.52 M | 3715.51x |
| COMPARE | w32-il64 | 1.91 G | w8 | 72.04 M | 80.58 M | 26.56x |
| REDUCE | w32-o64 | 18.04 M | w8 | 32.77 k | 10.42 M | 550.56x |
| MODMUL | w32-o64 | 2.88 M | w8 | 5.70 k | 3.17 M | 506.16x |
| MODEXP | w32-o64 | 22.50 k | none | n/a | 11.64 k | n/a |
| EXPONENTIATION | w32-o64 | 15.85 k | none | n/a | 47.37 k | n/a |
| DIVIDE | w32-il64 | 8.58 M | w8 | 26.78 k | 7.61 M | 320.40x |
| ISQRT | w32-o64 | 151.43 k | w8 | 579.9 | 3.99 M | 261.15x |
| MODMUL_R2 | w32-il64 | 307.03 M | w8 | 184.83 k | 4.88 M | 1661.17x |

### p1024 (1024-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w32-il64 | 242.24 M | w8 | 15.99 M | 25.99 M | 15.15x |
| SUBTRACT | w32-il | 245.29 M | w8 | 17.18 M | 29.34 M | 14.27x |
| ADDMOD | w32-il | 309.64 M | w8 | 14.08 M | 8.17 M | 21.99x |
| SUBTRACTMOD | w32-il64 | 300.33 M | w8 | 13.56 M | 11.49 M | 22.15x |
| MULTIPLYOPERANDSCANNING | w32-il64 | 72.11 M | w8 | 91.27 k | 2.89 M | 790.07x |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 26.29 M | w8 | 127.46 k | 3.49 M | 206.26x |
| MONTGOMERYMULTIPLICATION | w32-il64 | 105.39 M | w8 | 32.50 k | 553.45 k | 3242.58x |
| COMPARE | w32-il | 851.91 M | w8 | 21.80 M | 58.68 M | 39.07x |
| REDUCE | w32-il64 | 3.56 M | w8 | 5.54 k | 22.91 M | 642.81x |
| MODMUL | w32-il64 | 465.50 k | w8 | 778.8 | 1.03 M | 597.73x |
| MODEXP | w32-il64 | 1.33 k | none | n/a | 2.79 k | n/a |
| EXPONENTIATION | w32-o64 | 943.5 | none | n/a | 12.61 k | n/a |
| DIVIDE | w32-o64 | 837.48 k | none | n/a | 11.51 M | n/a |
| ISQRT | w32-il64 | 11.05 k | none | n/a | 1.22 M | n/a |
| MODMUL_R2 | w32-il64 | 58.62 M | none | n/a | 1.41 M | n/a |

### p2048 (2048-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w32-il | 67.26 M | none | n/a | 19.39 M | n/a |
| SUBTRACT | w32-il | 67.04 M | none | n/a | 22.01 M | n/a |
| ADDMOD | w32-il64 | 76.08 M | none | n/a | 8.53 M | n/a |
| SUBTRACTMOD | w32-il64 | 73.58 M | none | n/a | 8.78 M | n/a |
| MULTIPLYOPERANDSCANNING | w32-il64 | 9.87 M | none | n/a | 1.11 M | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 2.14 M | none | n/a | 1.12 M | n/a |
| MONTGOMERYMULTIPLICATION | w32 | 15.45 M | none | n/a | 268.05 k | n/a |
| COMPARE | w32-il64 | 156.85 M | none | n/a | 85.54 M | n/a |
| REDUCE | w32-o64 | 484.02 k | none | n/a | 15.96 M | n/a |
| MODMUL | w32-o64 | 55.10 k | none | n/a | 397.91 k | n/a |
| MODEXP | w32-il | 41.5 | none | n/a | 397.8 | n/a |
| EXPONENTIATION | w32-il | 49.2 | none | n/a | 2.17 k | n/a |
| DIVIDE | w32-o64 | 33.39 k | none | n/a | 9.59 M | n/a |
| ISQRT | w32-o64 | 2.45 k | none | n/a | 1.11 M | n/a |
| MODMUL_R2 | w32-il64 | 7.68 M | none | n/a | 404.80 k | n/a |

## 6. Raw data

Also written to `NVIDIA_GeForce_RTX_2060_SUPER_Report.csv` for analysis.

```csv
kind,device,device_type,kernel,modulus,bits,operation,items,seconds,ops_per_sec,mismatches
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,secp256k1,256,ADD,50000,0.001442691,34657459.811,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,secp256k1,256,ADD,50000,0.007409031,6748520.792,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,secp256k1,256,ADD,50000,0.005556150,8999036.553,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,secp256k1,256,ADD,50000,0.000034816,1436121323.529,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,secp256k1,256,ADD,50000,0.000110238,453564686.103,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,secp256k1,256,ADD,50000,0.001300490,38447047.012,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,secp256k1,256,ADD,50000,0.000046904,1066013228.096,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,secp256k1,256,ADD,50000,0.001117965,44724111.009,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,secp256k1,256,ADD,50000,0.000033449,1494795946.096,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,secp256k1,256,ADD,50000,0.001109317,45072781.248,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,secp256k1,256,ADD,50000,0.000031125,1606411968.702,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,secp256k1,256,ADD,50000,0.000996037,50198942.390,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,secp256k1,256,ADD,50000,0.000034906,1432419722.519,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,secp256k1,256,ADD,50000,0.001109137,45080104.608,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,secp256k1,256,ADD,50000,0.000021038,2376639199.628,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,secp256k1,256,ADD,50000,0.001081742,46221727.541,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,secp256k1,256,ADD,50000,0.000020988,2382334148.344,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,secp256k1,256,ADD,50000,0.001080670,46267596.320,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,secp256k1,256,ADD,50000,0.000455832,109689508.546,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,secp256k1,256,ADD,50000,0.002188719,22844409.061,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,secp256k1,256,SUBTRACT,50000,0.000899606,55579863.450,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,secp256k1,256,SUBTRACT,50000,0.005877577,8506907.178,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,secp256k1,256,SUBTRACT,50000,0.006352522,7870889.195,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,secp256k1,256,SUBTRACT,50000,0.000032896,1519941634.241,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,secp256k1,256,SUBTRACT,50000,0.000104547,478253742.100,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,secp256k1,256,SUBTRACT,50000,0.001215182,41146091.638,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,secp256k1,256,SUBTRACT,50000,0.000046750,1069517230.938,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,secp256k1,256,SUBTRACT,50000,0.001082519,46188562.786,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,secp256k1,256,SUBTRACT,50000,0.000032724,1527914370.687,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,secp256k1,256,SUBTRACT,50000,0.001125792,44413194.304,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,secp256k1,256,SUBTRACT,50000,0.000030358,1647020115.657,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,secp256k1,256,SUBTRACT,50000,0.001103030,45329680.644,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,secp256k1,256,SUBTRACT,50000,0.000032975,1516306079.392,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,secp256k1,256,SUBTRACT,50000,0.001103947,45292031.916,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,secp256k1,256,SUBTRACT,50000,0.000020960,2385456820.403,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,secp256k1,256,SUBTRACT,50000,0.001095276,45650599.829,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,secp256k1,256,SUBTRACT,50000,0.000020361,2455670266.438,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,secp256k1,256,SUBTRACT,50000,0.001077822,46389851.573,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,secp256k1,256,SUBTRACT,50000,0.000483498,103413048.972,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,secp256k1,256,SUBTRACT,50000,0.002252571,22196860.394,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,secp256k1,256,ADDMOD,50000,0.004916154,10170552.335,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,secp256k1,256,ADDMOD,50000,0.009052312,5523450.804,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,secp256k1,256,ADDMOD,50000,0.006624805,7547392.178,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,secp256k1,256,ADDMOD,50000,0.000033056,1512584704.743,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,secp256k1,256,ADDMOD,50000,0.000158619,315220479.579,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,secp256k1,256,ADDMOD,50000,0.001284390,38928979.007,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,secp256k1,256,ADDMOD,50000,0.000068414,730844296.819,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,secp256k1,256,ADDMOD,50000,0.001121248,44593162.897,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,secp256k1,256,ADDMOD,50000,0.000041319,1210095370.329,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,secp256k1,256,ADDMOD,50000,0.001112110,44959563.260,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,secp256k1,256,ADDMOD,50000,0.000026743,1869620629.973,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,secp256k1,256,ADDMOD,50000,0.001084942,46085416.896,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,secp256k1,256,ADDMOD,50000,0.000024971,2002315755.711,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,secp256k1,256,ADDMOD,50000,0.001114833,44849759.950,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,secp256k1,256,ADDMOD,50000,0.000018829,2655476255.719,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,secp256k1,256,ADDMOD,50000,0.001106958,45168835.937,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,secp256k1,256,ADDMOD,50000,0.000018720,2670933121.067,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,secp256k1,256,ADDMOD,50000,0.001067359,46844594.911,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,secp256k1,256,ADDMOD,50000,0.000573078,87248132.246,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,secp256k1,256,ADDMOD,50000,0.002300553,21733904.623,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,secp256k1,256,SUBTRACTMOD,50000,0.003522857,14193025.346,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,secp256k1,256,SUBTRACTMOD,50000,0.008282059,6037145.830,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,secp256k1,256,SUBTRACTMOD,50000,0.018398907,2717552.728,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,secp256k1,256,SUBTRACTMOD,50000,0.000032768,1525878906.250,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.000157802,316852965.530,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.001225975,40783854.522,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,secp256k1,256,SUBTRACTMOD,50000,0.000066821,748267785.389,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,secp256k1,256,SUBTRACTMOD,50000,0.001179285,42398576.266,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,secp256k1,256,SUBTRACTMOD,50000,0.000041080,1217132164.273,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,secp256k1,256,SUBTRACTMOD,50000,0.001177717,42455021.019,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,secp256k1,256,SUBTRACTMOD,50000,0.000027262,1834045305.321,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,secp256k1,256,SUBTRACTMOD,50000,0.001081920,46214147.942,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,secp256k1,256,SUBTRACTMOD,50000,0.000026392,1894526473.287,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,secp256k1,256,SUBTRACTMOD,50000,0.001116887,44767296.839,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,secp256k1,256,SUBTRACTMOD,50000,0.000018359,2723435864.658,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,secp256k1,256,SUBTRACTMOD,50000,0.001087093,45994233.661,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,secp256k1,256,SUBTRACTMOD,50000,0.000018810,2658171570.035,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,secp256k1,256,SUBTRACTMOD,50000,0.001072369,46625739.515,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.000479649,104242831.970,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.002204926,22676500.525,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.002489969,20080570.201,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.006109027,8184609.174,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.006854905,7294047.306,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.004305600,11612784.099,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.005734088,8719782.619,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000913699,54722634.983,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.002331998,21440844.166,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000252402,198096748.329,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001666892,29995947.171,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000054868,911279003.293,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001421891,35164447.204,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000054984,909349602.805,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001504212,33239992.063,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000026904,1858456494.046,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001468287,34053291.436,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000023978,2085259504.389,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001409200,35481123.535,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.005823034,8586589.257,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.008033440,6223983.761,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001782028,28057920.455,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.007415317,6742800.073,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.007960531,6280988.095,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000032768,1525878906.250,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001275860,39189259.707,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.002711896,18437286.203,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000300198,166556815.735,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001574459,31756939.599,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000099857,500716665.190,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001452557,34422061.576,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000088280,566382259.638,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001509853,33115802.843,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000093179,536599929.036,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001550449,32248721.574,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000049111,1018102521.216,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001414937,35337263.185,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000047053,1062626748.214,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001482978,33715941.898,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.004963046,10074458.429,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.006953789,7190324.184,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.013851925,3609606.538,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.005977766,8364328.576,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.006497894,7694800.869,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000031136,1605858170.606,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000322721,154932604.561,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.001450713,34465815.877,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000064349,777015243.002,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.001161633,43042869.381,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000029295,1706763243.312,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.001102701,45343195.133,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000031803,1572188450.275,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.001136869,43980433.586,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000031051,1610266528.696,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.001108266,45115524.943,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000029826,1676385730.121,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.001104563,45266770.291,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000031602,1582173173.212,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.001125170,44437732.651,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.093457514,535002.461,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.100442355,497797.966,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,secp256k1,256,COMPARE,50000,0.000611115,81817628.090,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,secp256k1,256,COMPARE,50000,0.005108006,9788555.776,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,secp256k1,256,COMPARE,50000,0.007394464,6761815.124,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,secp256k1,256,COMPARE,50000,0.000033376,1498082454.458,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,secp256k1,256,COMPARE,50000,0.000099679,501610221.481,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,secp256k1,256,COMPARE,50000,0.001188304,42076781.708,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,secp256k1,256,COMPARE,50000,0.000043987,1136691816.815,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,secp256k1,256,COMPARE,50000,0.001116759,44772430.231,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,secp256k1,256,COMPARE,50000,0.000022303,2241866215.680,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,secp256k1,256,COMPARE,50000,0.001117034,44761380.900,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,secp256k1,256,COMPARE,50000,0.000024816,2014827411.244,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,secp256k1,256,COMPARE,50000,0.001098107,45532900.087,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,secp256k1,256,COMPARE,50000,0.000011196,4465922821.611,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,secp256k1,256,COMPARE,50000,0.001109221,45076679.174,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,secp256k1,256,COMPARE,50000,0.000010887,4592761982.976,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,secp256k1,256,COMPARE,50000,0.001153836,43333713.664,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,secp256k1,256,COMPARE,50000,0.000475014,105260034.605,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,secp256k1,256,COMPARE,50000,0.001977955,25278633.853,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,secp256k1,256,REDUCE,6250,0.000219841,28429634.933,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,secp256k1,256,REDUCE,6250,0.004743883,1317486.140,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,secp256k1,256,REDUCE,6250,0.005265082,1187065.980,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,secp256k1,256,REDUCE,50000,0.000031936,1565631262.525,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,secp256k1,256,REDUCE,23800,0.000296184,80355365.879,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,secp256k1,256,REDUCE,23800,0.001022966,23265677.816,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,secp256k1,256,REDUCE,23800,0.000146493,162464750.605,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,secp256k1,256,REDUCE,23800,0.000818697,29070591.058,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,secp256k1,256,REDUCE,23800,0.000053118,448055253.504,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,secp256k1,256,REDUCE,23800,0.000776938,30633068.552,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,secp256k1,256,REDUCE,23800,0.000050752,468947424.257,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,secp256k1,256,REDUCE,23800,0.000767121,31025088.821,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,secp256k1,256,REDUCE,23800,0.000050474,471525936.383,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,secp256k1,256,REDUCE,23800,0.000737451,32273329.963,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,secp256k1,256,REDUCE,23800,0.000047767,498251209.530,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,secp256k1,256,REDUCE,23800,0.000776019,30669354.248,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,secp256k1,256,REDUCE,23800,0.011722685,2030251.630,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,secp256k1,256,REDUCE,23800,0.012665891,1879062.423,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,secp256k1,256,MODMUL,3125,0.000502561,6218147.922,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,secp256k1,256,MODMUL,3125,0.011134325,280663.622,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,secp256k1,256,MODMUL,3125,0.007351137,425104.321,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,secp256k1,256,MODMUL,50000,0.000084512,591631957.592,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,secp256k1,256,MODMUL,23800,0.000808398,29440934.978,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,secp256k1,256,MODMUL,23800,0.001507936,15783162.735,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,secp256k1,256,MODMUL,23800,0.000354684,67101903.587,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,secp256k1,256,MODMUL,23800,0.000959218,24811877.646,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,secp256k1,256,MODMUL,23800,0.000154241,154303990.648,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,secp256k1,256,MODMUL,23800,0.000887091,26829257.330,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,secp256k1,256,MODMUL,23800,0.000129317,184043955.919,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,secp256k1,256,MODMUL,23800,0.000855420,27822597.073,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,secp256k1,256,MODMUL,23800,0.000154474,154071415.978,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,secp256k1,256,MODMUL,23800,0.000894144,26617645.113,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,secp256k1,256,MODMUL,23800,0.000128191,185660511.255,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,secp256k1,256,MODMUL,23800,0.000852642,27913235.314,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,secp256k1,256,MODMUL,23800,0.099235886,239832.595,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,secp256k1,256,MODMUL,23800,0.102235260,232796.395,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,secp256k1,256,MODEXP,781,0.012901203,60536.991,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,secp256k1,256,MODEXP,781,0.010969883,71194.926,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,secp256k1,256,MODEXP,781,0.015468755,50488.873,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,secp256k1,256,MODEXP,50000,0.136671394,365841.004,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,secp256k1,256,MODEXP,23800,0.027711872,858837.699,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,secp256k1,256,MODEXP,23800,0.028267214,841964.827,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,secp256k1,256,MODEXP,23800,0.005228211,4552226.224,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,secp256k1,256,MODEXP,23800,0.005942592,4004986.422,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,secp256k1,256,MODEXP,23800,0.003640182,6538134.900,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,secp256k1,256,MODEXP,23800,0.004340128,5483709.183,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,secp256k1,256,MODEXP,23800,0.001976730,12040084.566,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,secp256k1,256,MODEXP,23800,0.002687801,8854821.748,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,secp256k1,256,MODEXP,23800,0.003640654,6537286.926,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,secp256k1,256,MODEXP,23800,0.004354642,5465431.761,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,secp256k1,256,MODEXP,23800,0.001996340,11921818.332,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,secp256k1,256,MODEXP,23800,0.002719908,8750295.897,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,secp256k1,256,MODEXP,23800,4.220276189,5639.441,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,secp256k1,256,MODEXP,23800,3.943783917,6034.813,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,secp256k1,256,EXPONENTIATION,781,0.004563249,171149.979,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,secp256k1,256,EXPONENTIATION,781,0.007766600,100558.806,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,secp256k1,256,EXPONENTIATION,781,0.015109615,51688.941,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,secp256k1,256,EXPONENTIATION,23800,0.040821099,583031.826,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,secp256k1,256,EXPONENTIATION,23800,0.041369868,575297.944,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,secp256k1,256,EXPONENTIATION,23800,0.013395867,1776667.406,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,secp256k1,256,EXPONENTIATION,23800,0.014156730,1681179.154,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,secp256k1,256,EXPONENTIATION,23800,0.000680107,34994485.400,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,secp256k1,256,EXPONENTIATION,23800,0.001392040,17097210.636,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,secp256k1,256,EXPONENTIATION,23800,0.000527953,45079747.870,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,secp256k1,256,EXPONENTIATION,23800,0.001276986,18637632.775,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,secp256k1,256,EXPONENTIATION,23800,0.000679508,35025349.514,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,secp256k1,256,EXPONENTIATION,23800,0.001400360,16995625.181,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,secp256k1,256,EXPONENTIATION,23800,0.000527866,45087184.384,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,secp256k1,256,EXPONENTIATION,23800,0.001212647,19626491.068,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,secp256k1,256,EXPONENTIATION,23800,1.486647103,16009.179,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,secp256k1,256,EXPONENTIATION,23800,1.688509164,14095.274,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,secp256k1,256,DIVIDE,6250,0.000526595,11868705.244,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,secp256k1,256,DIVIDE,6250,0.005942492,1051747.347,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,secp256k1,256,DIVIDE,6250,0.007309144,855093.303,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,secp256k1,256,DIVIDE,50000,0.000049152,1017252604.167,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,secp256k1,256,DIVIDE,23800,0.000440608,54016295.556,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,secp256k1,256,DIVIDE,23800,0.001119823,21253363.926,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,secp256k1,256,DIVIDE,23800,0.000288306,82551234.345,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,secp256k1,256,DIVIDE,23800,0.001156727,20575298.272,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,secp256k1,256,DIVIDE,23800,0.000121236,196311573.648,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,secp256k1,256,DIVIDE,23800,0.000993856,23947125.659,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,secp256k1,256,DIVIDE,23800,0.000114126,208541197.395,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,secp256k1,256,DIVIDE,23800,0.000968249,24580452.750,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,secp256k1,256,DIVIDE,23800,0.000117368,202780884.531,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,secp256k1,256,DIVIDE,23800,0.001000065,23798444.617,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,secp256k1,256,DIVIDE,23800,0.000114829,207264210.899,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,secp256k1,256,DIVIDE,23800,0.001010094,23562167.301,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,secp256k1,256,DIVIDE,23800,0.073353423,324456.571,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,secp256k1,256,DIVIDE,23800,0.020612263,1154652.471,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,secp256k1,256,ISQRT,1562,0.000227319,6871398.666,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,secp256k1,256,ISQRT,1562,0.007119553,219395.791,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,secp256k1,256,ISQRT,23800,0.005404453,4403776.092,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,secp256k1,256,ISQRT,23800,0.005935295,4009910.194,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,secp256k1,256,ISQRT,23800,0.004305261,5528119.890,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,secp256k1,256,ISQRT,23800,0.005026099,4735282.953,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,secp256k1,256,ISQRT,23800,0.001008587,23597359.476,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,secp256k1,256,ISQRT,23800,0.001734602,13720728.301,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,secp256k1,256,ISQRT,23800,0.000859883,27678194.343,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,secp256k1,256,ISQRT,23800,0.001579812,15065085.634,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,secp256k1,256,ISQRT,23800,0.000996020,23895098.744,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,secp256k1,256,ISQRT,23800,0.001712455,13898175.844,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,secp256k1,256,ISQRT,23800,0.000864378,27534244.435,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,secp256k1,256,ISQRT,23800,0.001585384,15012136.813,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,secp256k1,256,ISQRT,23800,0.369539473,64404.487,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,secp256k1,256,ISQRT,23800,0.410021898,58045.680,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,secp256k1,256,MODMUL_R2,50000,0.008007843,6243878.621,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,secp256k1,256,MODMUL_R2,50000,0.005448160,9177409.972,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,secp256k1,256,MODMUL_R2,50000,0.008315986,6012516.118,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,secp256k1,256,MODMUL_R2,50000,0.000033664,1485266159.696,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,secp256k1,256,MODMUL_R2,50000,0.000302867,165088948.681,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,secp256k1,256,MODMUL_R2,50000,0.001188674,42063677.312,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,secp256k1,256,MODMUL_R2,50000,0.000068396,731038354.019,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,secp256k1,256,MODMUL_R2,50000,0.001136967,43976650.893,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,secp256k1,256,MODMUL_R2,50000,0.000063569,786543375.771,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,secp256k1,256,MODMUL_R2,50000,0.001145932,43632592.850,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,secp256k1,256,MODMUL_R2,50000,0.000042721,1170380109.654,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,secp256k1,256,MODMUL_R2,50000,0.001126761,44374979.553,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,secp256k1,256,MODMUL_R2,50000,0.000052481,952726503.523,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,secp256k1,256,MODMUL_R2,50000,0.001123802,44491812.574,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,secp256k1,256,MODMUL_R2,50000,0.000042679,1171542163.837,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,secp256k1,256,MODMUL_R2,50000,0.001130459,44229825.650,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,secp256k1,256,MODMUL_R2,50000,0.088239715,566638.273,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,secp256k1,256,MODMUL_R2,50000,0.023088326,2165596.611,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,rsa256(composite),256,ADD,50000,0.002133955,23430674.549,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,rsa256(composite),256,ADD,50000,0.006990466,7152598.710,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,rsa256(composite),256,ADD,50000,0.007494356,6671687.318,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,rsa256(composite),256,ADD,50000,0.000035200,1420454545.455,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,rsa256(composite),256,ADD,50000,0.000082968,602640031.879,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,rsa256(composite),256,ADD,50000,0.000983590,50834175.597,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,rsa256(composite),256,ADD,50000,0.000047016,1063468716.201,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,rsa256(composite),256,ADD,50000,0.001121664,44576630.832,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,rsa256(composite),256,ADD,50000,0.000032475,1539657614.821,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,rsa256(composite),256,ADD,50000,0.001170428,42719415.311,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,rsa256(composite),256,ADD,50000,0.000036895,1355204180.182,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,rsa256(composite),256,ADD,50000,0.001155590,43267934.580,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,rsa256(composite),256,ADD,50000,0.000033017,1514381371.733,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,rsa256(composite),256,ADD,50000,0.001087647,45970800.410,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,rsa256(composite),256,ADD,50000,0.000021070,2373067438.725,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,rsa256(composite),256,ADD,50000,0.001169708,42745707.542,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,rsa256(composite),256,ADD,50000,0.000020702,2415236799.604,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,rsa256(composite),256,ADD,50000,0.001079489,46318211.187,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,rsa256(composite),256,ADD,50000,0.000499521,100095908.598,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,rsa256(composite),256,ADD,50000,0.002267546,22050269.123,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACT,50000,0.001189139,42047238.312,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACT,50000,0.009298729,5377079.081,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACT,50000,0.007578249,6597830.107,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,rsa256(composite),256,SUBTRACT,50000,0.000034816,1436121323.529,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,rsa256(composite),256,SUBTRACT,50000,0.000082737,604325752.493,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,rsa256(composite),256,SUBTRACT,50000,0.000973489,51361641.519,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,rsa256(composite),256,SUBTRACT,50000,0.000046609,1072754889.502,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,rsa256(composite),256,SUBTRACT,50000,0.001114608,44858810.087,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,rsa256(composite),256,SUBTRACT,50000,0.000033232,1504577627.689,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,rsa256(composite),256,SUBTRACT,50000,0.001175492,42535378.375,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,rsa256(composite),256,SUBTRACT,50000,0.000036260,1378925648.533,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,rsa256(composite),256,SUBTRACT,50000,0.001124393,44468429.928,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,rsa256(composite),256,SUBTRACT,50000,0.000033150,1508297383.023,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,rsa256(composite),256,SUBTRACT,50000,0.001098415,45520121.382,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,rsa256(composite),256,SUBTRACT,50000,0.000020541,2434181551.087,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,rsa256(composite),256,SUBTRACT,50000,0.001088705,45926126.872,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,rsa256(composite),256,SUBTRACT,50000,0.000020839,2399369453.197,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,rsa256(composite),256,SUBTRACT,50000,0.001077029,46424008.626,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,rsa256(composite),256,SUBTRACT,50000,0.000437144,114378858.046,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,rsa256(composite),256,SUBTRACT,50000,0.002191887,22811392.463,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,rsa256(composite),256,ADDMOD,50000,0.003843169,13010097.536,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,rsa256(composite),256,ADDMOD,50000,0.006171728,8101458.757,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,rsa256(composite),256,ADDMOD,50000,0.007916103,6316239.096,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,rsa256(composite),256,ADDMOD,50000,0.000034784,1437442502.300,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,rsa256(composite),256,ADDMOD,50000,0.000113710,439715723.003,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,rsa256(composite),256,ADDMOD,50000,0.001032216,48439473.498,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,rsa256(composite),256,ADDMOD,50000,0.000063345,789324519.786,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,rsa256(composite),256,ADDMOD,50000,0.001205169,41487968.842,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,rsa256(composite),256,ADDMOD,50000,0.000037768,1323874712.105,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,rsa256(composite),256,ADDMOD,50000,0.001101290,45401307.229,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,rsa256(composite),256,ADDMOD,50000,0.000030906,1617811999.397,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,rsa256(composite),256,ADDMOD,50000,0.001097935,45540026.075,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,rsa256(composite),256,ADDMOD,50000,0.000024413,2048070314.914,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,rsa256(composite),256,ADDMOD,50000,0.001112209,44955572.628,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,rsa256(composite),256,ADDMOD,50000,0.000018832,2655082277.886,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,rsa256(composite),256,ADDMOD,50000,0.001094445,45685250.966,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,rsa256(composite),256,ADDMOD,50000,0.000020949,2386729403.396,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,rsa256(composite),256,ADDMOD,50000,0.001078940,46341780.039,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,rsa256(composite),256,ADDMOD,50000,0.000509464,98142326.331,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,rsa256(composite),256,ADDMOD,50000,0.002010486,24869607.545,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACTMOD,50000,0.002325209,21503440.356,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACTMOD,50000,0.007677573,6512474.898,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACTMOD,50000,0.008118799,6158546.390,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,rsa256(composite),256,SUBTRACTMOD,50000,0.000034816,1436121323.529,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,rsa256(composite),256,SUBTRACTMOD,50000,0.000124926,400237749.184,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,rsa256(composite),256,SUBTRACTMOD,50000,0.001019265,49054954.485,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,rsa256(composite),256,SUBTRACTMOD,50000,0.000066791,748606882.704,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,rsa256(composite),256,SUBTRACTMOD,50000,0.001143022,43743690.940,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,rsa256(composite),256,SUBTRACTMOD,50000,0.000041697,1199122022.693,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,rsa256(composite),256,SUBTRACTMOD,50000,0.001106663,45180866.814,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,rsa256(composite),256,SUBTRACTMOD,50000,0.000030321,1649018373.929,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,rsa256(composite),256,SUBTRACTMOD,50000,0.001105346,45234713.563,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,rsa256(composite),256,SUBTRACTMOD,50000,0.000026180,1909858992.192,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,rsa256(composite),256,SUBTRACTMOD,50000,0.001082367,46195040.884,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,rsa256(composite),256,SUBTRACTMOD,50000,0.000018321,2729111996.747,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,rsa256(composite),256,SUBTRACTMOD,50000,0.001089305,45900820.300,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,rsa256(composite),256,SUBTRACTMOD,50000,0.000021445,2331535022.691,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,rsa256(composite),256,SUBTRACTMOD,50000,0.001247756,40071946.215,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,rsa256(composite),256,SUBTRACTMOD,50000,0.000540125,92571132.588,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,rsa256(composite),256,SUBTRACTMOD,50000,0.002218943,22533250.286,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001206849,41430195.675,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.007455347,6706595.855,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.007879913,6345247.839,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.003365261,14857688.537,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.004567607,10946650.325,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000900725,55510844.485,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.002325875,21497288.248,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000250680,199457550.095,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001709559,29247308.791,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000068216,732964594.895,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001476766,33857765.840,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000054303,920757899.070,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001492568,33499304.707,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000026625,1877926131.137,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001396404,35806253.729,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000026081,1917122239.680,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001427156,35034702.941,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.005840677,8560651.002,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.007852645,6367281.241,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001223585,40863524.645,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.007170102,6973401.403,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.006829871,7320782.631,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000035232,1419164396.004,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001005831,49710154.560,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.002168811,23054102.438,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000298943,167255757.848,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001741670,28708077.840,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000100073,499635571.066,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001523810,32812490.324,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000107367,465692189.323,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001553171,32192208.632,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000113108,442056437.091,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001552572,32204625.444,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000055440,901872064.641,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001502417,33279708.010,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000052838,946286495.871,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001399817,35718956.090,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.005031343,9937705.260,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.007140342,7002465.565,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.013865293,3606126.462,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.007459893,6702509.101,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.007252958,6893738.948,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000034400,1453488372.093,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000254418,196527153.258,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.001134861,44058248.998,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000064043,780727126.248,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.001164231,42946804.113,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000028325,1765238831.440,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.001154327,43315288.648,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000034330,1456454327.686,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.001110043,45043303.735,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000033649,1485921623.005,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.001109050,45083644.138,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000031829,1570900374.532,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.001095235,45652307.845,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000029771,1679479805.421,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.001082224,46201162.966,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.098131482,509520.482,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.024571888,2034845.662,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,rsa256(composite),256,COMPARE,50000,0.000611889,81714143.599,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,rsa256(composite),256,COMPARE,50000,0.007925555,6308706.369,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,rsa256(composite),256,COMPARE,50000,0.006305710,7929321.269,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,rsa256(composite),256,COMPARE,50000,0.000034976,1429551692.589,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,rsa256(composite),256,COMPARE,50000,0.000076531,653330305.630,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,rsa256(composite),256,COMPARE,50000,0.000978206,51113987.512,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,rsa256(composite),256,COMPARE,50000,0.000044169,1132006182.198,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,rsa256(composite),256,COMPARE,50000,0.001233253,40543191.297,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,rsa256(composite),256,COMPARE,50000,0.000025300,1976296817.655,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,rsa256(composite),256,COMPARE,50000,0.001084998,46083023.636,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,rsa256(composite),256,COMPARE,50000,0.000028739,1739782918.807,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,rsa256(composite),256,COMPARE,50000,0.001131473,44190197.871,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,rsa256(composite),256,COMPARE,50000,0.000013089,3819921818.635,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,rsa256(composite),256,COMPARE,50000,0.001084661,46097367.139,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,rsa256(composite),256,COMPARE,50000,0.000012120,4125492081.300,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,rsa256(composite),256,COMPARE,50000,0.001084451,46106274.494,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,rsa256(composite),256,COMPARE,50000,0.000459685,108770138.356,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,rsa256(composite),256,COMPARE,50000,0.001883999,26539292.516,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,rsa256(composite),256,REDUCE,6250,0.000418720,14926448.353,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,rsa256(composite),256,REDUCE,6250,0.008791341,710926.810,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,rsa256(composite),256,REDUCE,6250,0.006281276,995020.794,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,rsa256(composite),256,REDUCE,50000,0.000032992,1515518913.676,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,rsa256(composite),256,REDUCE,23800,0.000233961,101726239.028,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,rsa256(composite),256,REDUCE,23800,0.000757819,31405933.762,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,rsa256(composite),256,REDUCE,23800,0.000147720,161115768.639,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,rsa256(composite),256,REDUCE,23800,0.000880050,27043917.736,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,rsa256(composite),256,REDUCE,23800,0.000066813,356217666.730,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,rsa256(composite),256,REDUCE,23800,0.000795472,29919345.453,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,rsa256(composite),256,REDUCE,23800,0.000062647,379907612.425,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,rsa256(composite),256,REDUCE,23800,0.000788966,30166058.048,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,rsa256(composite),256,REDUCE,23800,0.000061836,384888478.390,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,rsa256(composite),256,REDUCE,23800,0.000805891,29532525.475,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,rsa256(composite),256,REDUCE,23800,0.000059531,399791233.103,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,rsa256(composite),256,REDUCE,23800,0.000803562,29618112.513,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,rsa256(composite),256,REDUCE,23800,0.011835712,2010863.392,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,rsa256(composite),256,REDUCE,23800,0.012757887,1865512.663,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,rsa256(composite),256,MODMUL,3125,0.000879109,3554734.031,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,rsa256(composite),256,MODMUL,3125,0.005151401,606631.057,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,rsa256(composite),256,MODMUL,3125,0.005181797,603072.604,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,rsa256(composite),256,MODMUL,50000,0.000065536,762939453.125,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,rsa256(composite),256,MODMUL,23800,0.000632914,37603857.954,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,rsa256(composite),256,MODMUL,23800,0.001158042,20551933.760,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,rsa256(composite),256,MODMUL,23800,0.000353919,67247048.915,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,rsa256(composite),256,MODMUL,23800,0.001078396,22069812.665,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,rsa256(composite),256,MODMUL,23800,0.000193301,123123961.597,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,rsa256(composite),256,MODMUL,23800,0.000924688,25738415.886,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,rsa256(composite),256,MODMUL,23800,0.000161224,147620660.558,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,rsa256(composite),256,MODMUL,23800,0.000895715,26570942.246,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,rsa256(composite),256,MODMUL,23800,0.000190971,124626282.464,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,rsa256(composite),256,MODMUL,23800,0.000914946,26012471.720,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,rsa256(composite),256,MODMUL,23800,0.000161174,147666722.974,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,rsa256(composite),256,MODMUL,23800,0.000894601,26604039.382,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,rsa256(composite),256,MODMUL,23800,0.098913976,240613.116,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,rsa256(composite),256,MODMUL,23800,0.102014852,233299.363,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,rsa256(composite),256,MODEXP,781,0.012366129,63156.384,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,rsa256(composite),256,MODEXP,781,0.007937707,98391.135,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,rsa256(composite),256,MODEXP,781,0.008984498,86927.508,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,rsa256(composite),256,MODEXP,50000,0.132571131,377156.019,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,rsa256(composite),256,MODEXP,23800,0.027895458,853185.489,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,rsa256(composite),256,MODEXP,23800,0.028382049,838558.204,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,rsa256(composite),256,MODEXP,23800,0.005202471,4574749.356,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,rsa256(composite),256,MODEXP,23800,0.005946806,4002148.582,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,rsa256(composite),256,MODEXP,23800,0.004561481,5217603.593,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,rsa256(composite),256,MODEXP,23800,0.005272205,4514240.257,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,rsa256(composite),256,MODEXP,23800,0.002472330,9626547.067,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,rsa256(composite),256,MODEXP,23800,0.003196498,7445648.533,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,rsa256(composite),256,MODEXP,23800,0.004544088,5237574.987,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,rsa256(composite),256,MODEXP,23800,0.005286459,4502068.500,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,rsa256(composite),256,MODEXP,23800,0.002533563,9393885.409,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,rsa256(composite),256,MODEXP,23800,0.003255817,7309992.918,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,rsa256(composite),256,MODEXP,23800,4.070579233,5846.834,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,rsa256(composite),256,MODEXP,23800,3.998289400,5952.546,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,rsa256(composite),256,EXPONENTIATION,781,0.006803718,114790.171,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,rsa256(composite),256,EXPONENTIATION,781,0.010256291,76148.386,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,rsa256(composite),256,EXPONENTIATION,781,0.010067899,77573.283,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,rsa256(composite),256,EXPONENTIATION,23800,0.040797982,583362.182,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,rsa256(composite),256,EXPONENTIATION,23800,0.041389059,575031.196,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,rsa256(composite),256,EXPONENTIATION,23800,0.013464998,1767545.752,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,rsa256(composite),256,EXPONENTIATION,23800,0.014196216,1676503.130,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,rsa256(composite),256,EXPONENTIATION,23800,0.000848399,28052838.030,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,rsa256(composite),256,EXPONENTIATION,23800,0.001563759,15219736.010,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,rsa256(composite),256,EXPONENTIATION,23800,0.000655507,36307762.728,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,rsa256(composite),256,EXPONENTIATION,23800,0.001366565,17415929.671,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,rsa256(composite),256,EXPONENTIATION,23800,0.000854403,27855714.591,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,rsa256(composite),256,EXPONENTIATION,23800,0.001575059,15110542.595,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,rsa256(composite),256,EXPONENTIATION,23800,0.000653042,36444797.285,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,rsa256(composite),256,EXPONENTIATION,23800,0.001400855,16989625.356,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,rsa256(composite),256,EXPONENTIATION,23800,1.615916391,14728.485,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,rsa256(composite),256,EXPONENTIATION,23800,1.748993221,13607.829,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,rsa256(composite),256,DIVIDE,6250,0.000546515,11436090.158,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,rsa256(composite),256,DIVIDE,6250,0.005276479,1184502.052,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,rsa256(composite),256,DIVIDE,6250,0.005858578,1066811.730,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,rsa256(composite),256,DIVIDE,50000,0.000038912,1284950657.895,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,rsa256(composite),256,DIVIDE,23800,0.000442479,53787886.988,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,rsa256(composite),256,DIVIDE,23800,0.001099216,21651798.999,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,rsa256(composite),256,DIVIDE,23800,0.000296089,80381272.859,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,rsa256(composite),256,DIVIDE,23800,0.001180426,20162216.862,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,rsa256(composite),256,DIVIDE,23800,0.000154074,154470945.932,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,rsa256(composite),256,DIVIDE,23800,0.001016117,23422495.731,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,rsa256(composite),256,DIVIDE,23800,0.000146476,162484377.316,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,rsa256(composite),256,DIVIDE,23800,0.001025280,23213171.401,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,rsa256(composite),256,DIVIDE,23800,0.000147117,161776181.579,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,rsa256(composite),256,DIVIDE,23800,0.000977415,24349940.339,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,rsa256(composite),256,DIVIDE,23800,0.000143815,165490044.464,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,rsa256(composite),256,DIVIDE,23800,0.001036776,22955773.211,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,rsa256(composite),256,DIVIDE,23800,0.019457562,1223174.832,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,rsa256(composite),256,DIVIDE,23800,0.075582184,314889.022,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,rsa256(composite),256,ISQRT,1562,0.000224675,6952262.777,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,rsa256(composite),256,ISQRT,1562,0.009649753,161869.426,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,rsa256(composite),256,ISQRT,23800,0.005401407,4406259.031,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,rsa256(composite),256,ISQRT,23800,0.005919574,4020559.704,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,rsa256(composite),256,ISQRT,23800,0.004339793,5484132.834,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,rsa256(composite),256,ISQRT,23800,0.005058927,4704554.942,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,rsa256(composite),256,ISQRT,23800,0.001271461,18718621.682,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,rsa256(composite),256,ISQRT,23800,0.001994286,11934094.526,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,rsa256(composite),256,ISQRT,23800,0.001075708,22124966.537,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,rsa256(composite),256,ISQRT,23800,0.001767896,13462329.482,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,rsa256(composite),256,ISQRT,23800,0.001255264,18960154.390,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,rsa256(composite),256,ISQRT,23800,0.001977796,12033595.773,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,rsa256(composite),256,ISQRT,23800,0.001081717,22002053.760,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,rsa256(composite),256,ISQRT,23800,0.001800836,13216083.616,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,rsa256(composite),256,ISQRT,23800,0.333578126,71347.604,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,rsa256(composite),256,ISQRT,23800,0.373941857,63646.258,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,rsa256(composite),256,MODMUL_R2,50000,0.007998742,6250982.840,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,rsa256(composite),256,MODMUL_R2,50000,0.006370686,7848448.353,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,rsa256(composite),256,MODMUL_R2,50000,0.007277050,6870916.047,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,rsa256(composite),256,MODMUL_R2,50000,0.000034784,1437442502.300,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,rsa256(composite),256,MODMUL_R2,50000,0.000301360,165914436.791,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,rsa256(composite),256,MODMUL_R2,50000,0.001190644,41994089.065,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,rsa256(composite),256,MODMUL_R2,50000,0.000068298,732085051.374,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,rsa256(composite),256,MODMUL_R2,50000,0.001154425,43311602.031,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,rsa256(composite),256,MODMUL_R2,50000,0.000070020,714085514.012,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,rsa256(composite),256,MODMUL_R2,50000,0.001132403,44153890.856,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,rsa256(composite),256,MODMUL_R2,50000,0.000047568,1051121685.332,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,rsa256(composite),256,MODMUL_R2,50000,0.001281220,39025304.000,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,rsa256(composite),256,MODMUL_R2,50000,0.000058390,856314906.173,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,rsa256(composite),256,MODMUL_R2,50000,0.001177610,42458882.248,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,rsa256(composite),256,MODMUL_R2,50000,0.000044545,1122468167.135,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,rsa256(composite),256,MODMUL_R2,50000,0.001142337,43769921.285,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,rsa256(composite),256,MODMUL_R2,50000,0.027022101,1850337.249,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,rsa256(composite),256,MODMUL_R2,50000,0.090998211,549461.349,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,ADD,25000,0.000770810,32433431.261,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,ADD,25000,0.011421669,2188821.914,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,ADD,25000,0.007964328,3138996.630,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,brainpoolP512r1,512,ADD,50000,0.000047104,1061480978.261,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,brainpoolP512r1,512,ADD,25000,0.000113102,221039138.687,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,brainpoolP512r1,512,ADD,25000,0.000994491,25138478.049,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,brainpoolP512r1,512,ADD,25000,0.000063550,393389836.818,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,brainpoolP512r1,512,ADD,25000,0.001129914,22125576.228,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,brainpoolP512r1,512,ADD,25000,0.000059034,423483267.206,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,brainpoolP512r1,512,ADD,25000,0.001137446,21979071.548,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,brainpoolP512r1,512,ADD,25000,0.000048677,513589882.621,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,brainpoolP512r1,512,ADD,25000,0.001131626,22092107.641,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,brainpoolP512r1,512,ADD,25000,0.000047945,521431330.310,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,brainpoolP512r1,512,ADD,25000,0.001194406,20930909.131,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,brainpoolP512r1,512,ADD,25000,0.000033411,748262570.907,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,brainpoolP512r1,512,ADD,25000,0.001147005,21795898.787,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,brainpoolP512r1,512,ADD,25000,0.000033245,751993769.697,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,brainpoolP512r1,512,ADD,25000,0.001105114,22622102.853,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,brainpoolP512r1,512,ADD,25000,0.000625841,39946227.931,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,brainpoolP512r1,512,ADD,25000,0.002078042,12030556.266,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACT,25000,0.000677787,36884756.681,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACT,25000,0.005339206,4682344.390,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACT,25000,0.004952554,5047900.885,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,brainpoolP512r1,512,SUBTRACT,50000,0.000046720,1070205479.452,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,brainpoolP512r1,512,SUBTRACT,25000,0.000113662,219950637.273,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,brainpoolP512r1,512,SUBTRACT,25000,0.001018342,24549706.910,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,brainpoolP512r1,512,SUBTRACT,25000,0.000063025,396665518.006,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,brainpoolP512r1,512,SUBTRACT,25000,0.001121965,22282339.552,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,brainpoolP512r1,512,SUBTRACT,25000,0.000060315,414488914.967,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,brainpoolP512r1,512,SUBTRACT,25000,0.001117981,22361738.820,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,brainpoolP512r1,512,SUBTRACT,25000,0.000048121,519518978.130,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,brainpoolP512r1,512,SUBTRACT,25000,0.001189910,21009994.553,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,brainpoolP512r1,512,SUBTRACT,25000,0.000047591,525308863.905,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,brainpoolP512r1,512,SUBTRACT,25000,0.001178648,21210745.727,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,brainpoolP512r1,512,SUBTRACT,25000,0.000032434,770801441.472,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,brainpoolP512r1,512,SUBTRACT,25000,0.001093008,22872657.790,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,brainpoolP512r1,512,SUBTRACT,25000,0.000032888,760149676.469,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,brainpoolP512r1,512,SUBTRACT,25000,0.001087150,22995915.081,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,brainpoolP512r1,512,SUBTRACT,25000,0.000576341,43377093.488,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,brainpoolP512r1,512,SUBTRACT,25000,0.002407864,10382645.326,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,ADDMOD,25000,0.002140077,11681823.586,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,ADDMOD,25000,0.006870482,3638755.078,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,ADDMOD,25000,0.007414329,3371849.355,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,brainpoolP512r1,512,ADDMOD,50000,0.000047104,1061480978.261,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,brainpoolP512r1,512,ADDMOD,25000,0.000153722,162630987.895,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,brainpoolP512r1,512,ADDMOD,25000,0.001041788,23997201.531,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,brainpoolP512r1,512,ADDMOD,25000,0.000080207,311693892.930,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,brainpoolP512r1,512,ADDMOD,25000,0.001146122,21812688.855,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,brainpoolP512r1,512,ADDMOD,25000,0.000072447,345079292.192,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,brainpoolP512r1,512,ADDMOD,25000,0.001177186,21237091.488,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,brainpoolP512r1,512,ADDMOD,25000,0.000058060,430592156.045,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,brainpoolP512r1,512,ADDMOD,25000,0.001176308,21252938.797,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,brainpoolP512r1,512,ADDMOD,25000,0.000057315,436186079.312,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,brainpoolP512r1,512,ADDMOD,25000,0.001155023,21644590.676,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,brainpoolP512r1,512,ADDMOD,25000,0.000026721,935592269.487,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,brainpoolP512r1,512,ADDMOD,25000,0.001092783,22877375.142,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,brainpoolP512r1,512,ADDMOD,25000,0.000026210,953843674.158,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,brainpoolP512r1,512,ADDMOD,25000,0.001142493,21881972.460,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,brainpoolP512r1,512,ADDMOD,25000,0.000666433,37513147.572,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,brainpoolP512r1,512,ADDMOD,25000,0.002280981,10960197.404,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACTMOD,25000,0.001935060,12919495.128,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACTMOD,25000,0.009087161,2751134.357,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACTMOD,25000,0.005282237,4732843.544,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000047072,1062202583.277,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000175722,142270222.599,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,brainpoolP512r1,512,SUBTRACTMOD,25000,0.001073261,23293499.768,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000087122,286954034.624,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,brainpoolP512r1,512,SUBTRACTMOD,25000,0.001229713,20329951.212,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000077346,323223446.258,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,brainpoolP512r1,512,SUBTRACTMOD,25000,0.001172135,21328599.073,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000058371,428294079.824,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,brainpoolP512r1,512,SUBTRACTMOD,25000,0.001129300,22137600.829,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000058422,427921976.726,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,brainpoolP512r1,512,SUBTRACTMOD,25000,0.001210970,20644609.385,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000027167,920229190.449,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,brainpoolP512r1,512,SUBTRACTMOD,25000,0.001113536,22451006.231,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000026620,939143742.784,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,brainpoolP512r1,512,SUBTRACTMOD,25000,0.001081786,23109928.639,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000719080,34766635.950,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,brainpoolP512r1,512,SUBTRACTMOD,25000,0.002521245,9915736.826,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.002027026,12333340.118,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.006423140,3892177.380,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.006510060,3840210.217,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.012057582,2073384.214,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.013258466,1885587.665,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.003163732,7902059.666,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.004571801,5468304.442,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.000797006,31367380.597,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.002224166,11240166.143,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.000140176,178347616.311,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.001548159,16148212.832,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.000138841,180061950.838,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.001572035,15902955.965,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.000073512,340080138.852,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.001525363,16389541.874,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.000067995,367674473.010,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.001428281,17503559.000,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.022697146,1101460.067,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.024887083,1004537.169,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.002221694,11252673.666,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.006501823,3845075.333,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.006951284,3596457.796,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000045760,1092657342.657,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.003465873,7213190.441,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.004669436,5353965.904,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000936678,26690077.653,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.002355522,10613358.690,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000381426,65543520.863,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.001802257,13871497.141,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000327370,76366269.571,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.001774248,14090475.042,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000324845,76959816.743,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.001789249,13972346.112,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000152580,163848500.441,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.001607128,15555699.811,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000143448,174278748.137,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.001582958,15793214.816,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.016687465,1498130.492,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.018931882,1320523.761,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.016473356,1517602.105,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.010759880,2323445.921,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.008732452,2862884.312,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000043712,1143850658.858,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000486702,51366113.367,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.001397882,17884197.528,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000147170,169871667.642,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.001216364,20553056.557,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000053483,467441783.845,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.001134835,22029630.700,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000073412,340543930.581,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.001236409,20219847.181,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000061856,404166788.623,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.001146027,21814488.059,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000061224,408335167.861,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.001212058,20626073.553,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000044339,563839348.016,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.001215381,20569680.051,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.164741491,151752.906,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.117629154,212532.346,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,COMPARE,25000,0.000310260,80577492.612,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,COMPARE,25000,0.005852065,4271996.259,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,COMPARE,25000,0.011039212,2264654.326,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,brainpoolP512r1,512,COMPARE,50000,0.000046944,1065098841.172,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,brainpoolP512r1,512,COMPARE,25000,0.000095135,262784280.057,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,brainpoolP512r1,512,COMPARE,25000,0.000985848,25358884.356,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,brainpoolP512r1,512,COMPARE,25000,0.000055096,453752524.552,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,brainpoolP512r1,512,COMPARE,25000,0.001140865,21913196.835,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,brainpoolP512r1,512,COMPARE,25000,0.000037847,660552822.481,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,brainpoolP512r1,512,COMPARE,25000,0.001178970,21204948.379,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,brainpoolP512r1,512,COMPARE,25000,0.000038340,652064653.728,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,brainpoolP512r1,512,COMPARE,25000,0.001127577,22171436.406,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,brainpoolP512r1,512,COMPARE,25000,0.000013495,1852492709.016,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,brainpoolP512r1,512,COMPARE,25000,0.001112839,22465069.799,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,brainpoolP512r1,512,COMPARE,25000,0.000013067,1913228010.406,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,brainpoolP512r1,512,COMPARE,25000,0.001097546,22778089.422,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,brainpoolP512r1,512,COMPARE,25000,0.000347023,72041258.838,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,brainpoolP512r1,512,COMPARE,25000,0.001972417,12674805.329,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,REDUCE,3125,0.000299791,10423933.047,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,REDUCE,3125,0.005037110,620395.388,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,REDUCE,3125,0.008551062,365451.679,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,brainpoolP512r1,512,REDUCE,50000,0.000044928,1112891737.892,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,brainpoolP512r1,512,REDUCE,23800,0.000684267,34781749.130,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,brainpoolP512r1,512,REDUCE,23800,0.001548993,15364821.526,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,brainpoolP512r1,512,REDUCE,23800,0.000543197,43814700.329,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,brainpoolP512r1,512,REDUCE,23800,0.001633581,14569217.264,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,brainpoolP512r1,512,REDUCE,23800,0.000181673,131004538.941,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,brainpoolP512r1,512,REDUCE,23800,0.001236160,19253166.660,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,brainpoolP512r1,512,REDUCE,23800,0.000173228,137391293.702,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,brainpoolP512r1,512,REDUCE,23800,0.001305676,18228104.682,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,brainpoolP512r1,512,REDUCE,23800,0.000190530,124914729.745,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,brainpoolP512r1,512,REDUCE,23800,0.001309839,18170170.987,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,brainpoolP512r1,512,REDUCE,23800,0.000178901,133034461.922,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,brainpoolP512r1,512,REDUCE,23800,0.001249003,19055195.008,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,brainpoolP512r1,512,REDUCE,23800,0.095371972,249549.208,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,brainpoolP512r1,512,REDUCE,23800,0.095214550,249961.796,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL,1562,0.000492202,3173496.051,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL,1562,0.007315833,213509.518,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL,1562,0.004939302,316238.994,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,brainpoolP512r1,512,MODMUL,50000,0.000197344,253364682.990,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,brainpoolP512r1,512,MODMUL,23800,0.001924130,12369225.748,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,brainpoolP512r1,512,MODMUL,23800,0.002784760,8546516.850,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,brainpoolP512r1,512,MODMUL,23800,0.001377980,17271661.534,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,brainpoolP512r1,512,MODMUL,23800,0.002469926,9635917.484,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,brainpoolP512r1,512,MODMUL,23800,0.000663749,35856905.805,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,brainpoolP512r1,512,MODMUL,23800,0.001769154,13452758.658,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,brainpoolP512r1,512,MODMUL,23800,0.000541490,43952793.909,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,brainpoolP512r1,512,MODMUL,23800,0.001631945,14583821.501,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,brainpoolP512r1,512,MODMUL,23800,0.000797191,29854838.870,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,brainpoolP512r1,512,MODMUL,23800,0.001937580,12283364.866,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,brainpoolP512r1,512,MODMUL,23800,0.000641532,37098680.987,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,brainpoolP512r1,512,MODMUL,23800,0.001755500,13557393.013,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,brainpoolP512r1,512,MODMUL,23800,0.274080300,86835.865,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,brainpoolP512r1,512,MODMUL,23800,0.273116562,87142.280,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,MODEXP,390,0.033512007,11637.620,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,MODEXP,390,0.007275999,53600.887,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,MODEXP,390,0.007337364,53152.601,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,brainpoolP512r1,512,MODEXP,50000,1.179690123,42384.012,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,brainpoolP512r1,512,MODEXP,23800,0.375358629,63406.029,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,brainpoolP512r1,512,MODEXP,23800,0.376145936,63273.314,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,brainpoolP512r1,512,MODEXP,23800,0.037030405,642715.081,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,brainpoolP512r1,512,MODEXP,23800,0.038231343,622525.868,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,brainpoolP512r1,512,MODEXP,23800,0.029012456,820337.317,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,brainpoolP512r1,512,MODEXP,23800,0.029837184,797662.404,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,brainpoolP512r1,512,MODEXP,23800,0.017332402,1373150.699,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,brainpoolP512r1,512,MODEXP,23800,0.018301029,1300473.340,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,brainpoolP512r1,512,MODEXP,23800,0.030683674,775656.789,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,brainpoolP512r1,512,MODEXP,23800,0.031807513,748250.893,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,brainpoolP512r1,512,MODEXP,23800,0.018263831,1303122.015,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,brainpoolP512r1,512,MODEXP,23800,0.018578793,1281030.479,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,brainpoolP512r1,512,MODEXP,23800,0.000000000,inf,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,brainpoolP512r1,512,MODEXP,23800,0.000000000,inf,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,EXPONENTIATION,390,0.008233075,47369.907,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,EXPONENTIATION,390,0.007635932,51074.316,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,EXPONENTIATION,390,0.016550521,23564.213,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,brainpoolP512r1,512,EXPONENTIATION,23800,0.338764616,70255.271,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,brainpoolP512r1,512,EXPONENTIATION,23800,0.339775194,70046.314,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,brainpoolP512r1,512,EXPONENTIATION,23800,0.085977990,276815.031,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,brainpoolP512r1,512,EXPONENTIATION,23800,0.087301590,272618.173,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,brainpoolP512r1,512,EXPONENTIATION,23800,0.025950856,917118.108,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,brainpoolP512r1,512,EXPONENTIATION,23800,0.026941194,883405.547,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,brainpoolP512r1,512,EXPONENTIATION,23800,0.024602960,967363.282,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,brainpoolP512r1,512,EXPONENTIATION,23800,0.025975966,916231.569,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,brainpoolP512r1,512,EXPONENTIATION,23800,0.026870671,885724.076,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,brainpoolP512r1,512,EXPONENTIATION,23800,0.028199053,843999.963,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,brainpoolP512r1,512,EXPONENTIATION,23800,0.025590587,930029.470,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,brainpoolP512r1,512,EXPONENTIATION,23800,0.026751843,889658.322,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,brainpoolP512r1,512,EXPONENTIATION,23800,0.000000000,inf,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,brainpoolP512r1,512,EXPONENTIATION,23800,0.000000000,inf,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,DIVIDE,3125,0.000410526,7612192.873,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,DIVIDE,3125,0.007109519,439551.559,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,DIVIDE,3125,0.004576239,682875.203,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,brainpoolP512r1,512,DIVIDE,50000,0.000096416,518586126.784,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,brainpoolP512r1,512,DIVIDE,23800,0.001333034,17854012.698,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,brainpoolP512r1,512,DIVIDE,23800,0.002714148,8768866.632,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,brainpoolP512r1,512,DIVIDE,23800,0.001251323,19019867.111,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,brainpoolP512r1,512,DIVIDE,23800,0.002620385,9082634.705,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,brainpoolP512r1,512,DIVIDE,23800,0.000393528,60478538.889,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,brainpoolP512r1,512,DIVIDE,23800,0.001799651,13224786.685,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,brainpoolP512r1,512,DIVIDE,23800,0.000379463,62720180.763,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,brainpoolP512r1,512,DIVIDE,23800,0.001905107,12492737.026,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,brainpoolP512r1,512,DIVIDE,23800,0.000380101,62614988.879,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,brainpoolP512r1,512,DIVIDE,23800,0.001758828,13531739.504,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,brainpoolP512r1,512,DIVIDE,23800,0.000364218,65345496.556,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,brainpoolP512r1,512,DIVIDE,23800,0.001749629,13602882.608,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,brainpoolP512r1,512,DIVIDE,23800,0.116695839,203949.003,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,brainpoolP512r1,512,DIVIDE,23800,0.132686953,179369.557,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,ISQRT,781,0.000195826,3988235.713,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,ISQRT,781,0.008762983,89124.900,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,brainpoolP512r1,512,ISQRT,23800,0.032112206,741151.203,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,brainpoolP512r1,512,ISQRT,23800,0.033239659,716012.161,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,brainpoolP512r1,512,ISQRT,23800,0.028212689,843592.049,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,brainpoolP512r1,512,ISQRT,23800,0.029304101,812173.009,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,brainpoolP512r1,512,ISQRT,23800,0.005604672,4246457.408,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,brainpoolP512r1,512,ISQRT,23800,0.006720869,3541208.556,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,brainpoolP512r1,512,ISQRT,23800,0.005157389,4614738.004,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,brainpoolP512r1,512,ISQRT,23800,0.006353983,3745681.927,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,brainpoolP512r1,512,ISQRT,23800,0.006674021,3566066.085,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,brainpoolP512r1,512,ISQRT,23800,0.007774236,3061394.014,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,brainpoolP512r1,512,ISQRT,23800,0.006212067,3831253.034,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,brainpoolP512r1,512,ISQRT,23800,0.007268495,3274405.533,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,brainpoolP512r1,512,ISQRT,23800,1.346871842,17670.575,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,brainpoolP512r1,512,ISQRT,23800,1.582828301,15036.375,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL_R2,25000,0.005127993,4875201.713,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL_R2,25000,0.009594617,2605627.733,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL_R2,25000,0.008605477,2905126.449,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,brainpoolP512r1,512,MODMUL_R2,50000,0.000067232,743693479.296,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,brainpoolP512r1,512,MODMUL_R2,25000,0.000699684,35730414.226,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,brainpoolP512r1,512,MODMUL_R2,25000,0.001855226,13475448.882,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,brainpoolP512r1,512,MODMUL_R2,25000,0.000153770,162580753.321,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,brainpoolP512r1,512,MODMUL_R2,25000,0.001251753,19971992.014,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,brainpoolP512r1,512,MODMUL_R2,25000,0.000093487,267417930.773,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,brainpoolP512r1,512,MODMUL_R2,25000,0.001174307,21289152.492,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,brainpoolP512r1,512,MODMUL_R2,25000,0.000096086,260183729.034,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,brainpoolP512r1,512,MODMUL_R2,25000,0.001546921,16161138.259,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,brainpoolP512r1,512,MODMUL_R2,25000,0.000089192,280293262.469,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,brainpoolP512r1,512,MODMUL_R2,25000,0.001149886,21741289.022,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,brainpoolP512r1,512,MODMUL_R2,25000,0.000081425,307030757.353,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,brainpoolP512r1,512,MODMUL_R2,25000,0.001176555,21248472.255,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,brainpoolP512r1,512,MODMUL_R2,25000,0.135260779,184828.153,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,brainpoolP512r1,512,MODMUL_R2,25000,0.109272649,228785.522,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,p1024,1024,ADD,12500,0.000480898,25993027.739,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,p1024,1024,ADD,12500,0.006297575,1984890.912,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,p1024,1024,ADD,12500,0.007657624,1632360.148,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,p1024,1024,ADD,50000,0.000072224,692290651.307,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p1024,1024,ADD,23800,0.000433601,54889175.679,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p1024,1024,ADD,23800,0.002351010,10123309.771,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p1024,1024,ADD,23800,0.000232400,102409474.353,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p1024,1024,ADD,23800,0.002211917,10759896.478,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,p1024,1024,ADD,23800,0.000287390,82814337.855,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,p1024,1024,ADD,23800,0.002197042,10832746.408,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p1024,1024,ADD,23800,0.000117779,202073755.466,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p1024,1024,ADD,23800,0.002045316,11636342.836,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p1024,1024,ADD,23800,0.000117701,202208066.238,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p1024,1024,ADD,23800,0.002631826,9043151.508,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p1024,1024,ADD,23800,0.000054939,433206285.948,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p1024,1024,ADD,23800,0.002048084,11620619.500,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p1024,1024,ADD,23800,0.000051603,461215987.063,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p1024,1024,ADD,23800,0.002076430,11461981.630,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,p1024,1024,ADD,23800,0.000781681,30447205.888,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,p1024,1024,ADD,23800,0.003660357,6502097.857,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,p1024,1024,SUBTRACT,12500,0.000426092,29336406.765,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,p1024,1024,SUBTRACT,12500,0.003284755,3805458.666,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,p1024,1024,SUBTRACT,12500,0.005854235,2135206.382,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,p1024,1024,SUBTRACT,50000,0.000071904,695371606.587,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p1024,1024,SUBTRACT,23800,0.000439821,54112889.154,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p1024,1024,SUBTRACT,23800,0.002354247,10109389.413,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p1024,1024,SUBTRACT,23800,0.000232800,102233716.495,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p1024,1024,SUBTRACT,23800,0.002146282,11088942.372,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,p1024,1024,SUBTRACT,23800,0.000291139,81747925.233,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,p1024,1024,SUBTRACT,23800,0.002203365,10801660.891,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p1024,1024,SUBTRACT,23800,0.000116977,203458957.280,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p1024,1024,SUBTRACT,23800,0.002034212,11699862.176,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p1024,1024,SUBTRACT,23800,0.000119840,198598470.676,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p1024,1024,SUBTRACT,23800,0.002662504,8938953.398,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p1024,1024,SUBTRACT,23800,0.000050959,467040506.811,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p1024,1024,SUBTRACT,23800,0.001715505,13873469.354,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p1024,1024,SUBTRACT,23800,0.000052581,452634331.610,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p1024,1024,SUBTRACT,23800,0.002046778,11628030.072,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,p1024,1024,SUBTRACT,23800,0.000727406,32718993.570,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,p1024,1024,SUBTRACT,23800,0.003310507,7189230.781,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,p1024,1024,ADDMOD,12500,0.001529267,8173850.301,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,p1024,1024,ADDMOD,12500,0.007881914,1585909.250,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,p1024,1024,ADDMOD,12500,0.007331476,1704977.368,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,p1024,1024,ADDMOD,50000,0.000071744,696922390.723,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p1024,1024,ADDMOD,23800,0.000584547,40715292.955,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p1024,1024,ADDMOD,23800,0.002484649,9578817.950,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p1024,1024,ADDMOD,23800,0.000305031,78024866.571,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p1024,1024,ADDMOD,23800,0.002197868,10828677.144,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,p1024,1024,ADDMOD,23800,0.000299547,79453342.944,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,p1024,1024,ADDMOD,23800,0.002204278,10797186.109,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p1024,1024,ADDMOD,23800,0.000124991,190413801.049,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p1024,1024,ADDMOD,23800,0.002081793,11432453.788,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p1024,1024,ADDMOD,23800,0.000123447,192795589.673,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p1024,1024,ADDMOD,23800,0.002853965,8339275.590,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p1024,1024,ADDMOD,23800,0.000040370,589546114.176,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p1024,1024,ADDMOD,23800,0.001978893,12026927.180,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p1024,1024,ADDMOD,23800,0.000041735,570266229.539,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p1024,1024,ADDMOD,23800,0.002000729,11895663.509,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,p1024,1024,ADDMOD,23800,0.000887890,26805125.858,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,p1024,1024,ADDMOD,23800,0.003635216,6547066.184,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,p1024,1024,SUBTRACTMOD,12500,0.001087964,11489350.257,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,p1024,1024,SUBTRACTMOD,12500,0.005903325,2117450.743,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,p1024,1024,SUBTRACTMOD,12500,0.007801804,1602193.593,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,p1024,1024,SUBTRACTMOD,50000,0.000071936,695062277.580,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p1024,1024,SUBTRACTMOD,23800,0.000587226,40529548.064,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p1024,1024,SUBTRACTMOD,23800,0.002495552,9536968.374,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p1024,1024,SUBTRACTMOD,23800,0.000305698,77854550.197,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p1024,1024,SUBTRACTMOD,23800,0.002217080,10734840.525,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,p1024,1024,SUBTRACTMOD,23800,0.000301236,79007868.033,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,p1024,1024,SUBTRACTMOD,23800,0.002179484,10920014.183,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p1024,1024,SUBTRACTMOD,23800,0.000122205,194754150.843,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p1024,1024,SUBTRACTMOD,23800,0.002062693,11538315.505,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p1024,1024,SUBTRACTMOD,23800,0.000123534,192659688.648,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p1024,1024,SUBTRACTMOD,23800,0.002343570,10155447.115,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p1024,1024,SUBTRACTMOD,23800,0.000042458,560553102.968,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p1024,1024,SUBTRACTMOD,23800,0.001965699,12107653.608,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p1024,1024,SUBTRACTMOD,23800,0.000041621,571829389.376,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p1024,1024,SUBTRACTMOD,23800,0.001974341,12054655.236,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,p1024,1024,SUBTRACTMOD,23800,0.000921791,25819303.068,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,p1024,1024,SUBTRACTMOD,23800,0.003712230,6411240.862,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.004326971,2888856.728,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.007751940,1612499.619,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.007368587,1696390.376,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,23800,0.058819345,404628.782,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,23800,0.061391252,387677.384,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p1024,1024,MULTIPLYOPERANDSCANNING,23800,0.015091444,1577052.544,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p1024,1024,MULTIPLYOPERANDSCANNING,23800,0.017574247,1354254.291,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,p1024,1024,MULTIPLYOPERANDSCANNING,23800,0.007104510,3349984.739,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,p1024,1024,MULTIPLYOPERANDSCANNING,23800,0.009627694,2472035.411,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p1024,1024,MULTIPLYOPERANDSCANNING,23800,0.000281850,84441976.745,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p1024,1024,MULTIPLYOPERANDSCANNING,23800,0.002721034,8746675.008,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p1024,1024,MULTIPLYOPERANDSCANNING,23800,0.000280417,84873587.985,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p1024,1024,MULTIPLYOPERANDSCANNING,23800,0.003599535,6611965.630,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p1024,1024,MULTIPLYOPERANDSCANNING,23800,0.000178615,133247415.283,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p1024,1024,MULTIPLYOPERANDSCANNING,23800,0.002680793,8877970.266,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p1024,1024,MULTIPLYOPERANDSCANNING,23800,0.000173350,137294598.009,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p1024,1024,MULTIPLYOPERANDSCANNING,23800,0.002655011,8964180.857,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,23800,0.136958236,173775.603,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,23800,0.135757884,175312.102,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.003577474,3494085.453,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.006954012,1797523.509,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.007688322,1625842.490,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000083968,595464939.024,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,23800,0.025544824,931695.592,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,23800,0.028005167,849843.162,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p1024,1024,MULTIPLYPRODUCTSCANNING,23800,0.006546023,3635795.185,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p1024,1024,MULTIPLYPRODUCTSCANNING,23800,0.009011703,2641010.228,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,p1024,1024,MULTIPLYPRODUCTSCANNING,23800,0.002135092,11147061.024,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,p1024,1024,MULTIPLYPRODUCTSCANNING,23800,0.004621601,5149730.879,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p1024,1024,MULTIPLYPRODUCTSCANNING,23800,0.001687201,14106201.395,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p1024,1024,MULTIPLYPRODUCTSCANNING,23800,0.004141456,5746771.209,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p1024,1024,MULTIPLYPRODUCTSCANNING,23800,0.001673080,14225260.964,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p1024,1024,MULTIPLYPRODUCTSCANNING,23800,0.004185898,5685756.962,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p1024,1024,MULTIPLYPRODUCTSCANNING,23800,0.000502000,47410349.387,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p1024,1024,MULTIPLYPRODUCTSCANNING,23800,0.002948021,8073212.366,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p1024,1024,MULTIPLYPRODUCTSCANNING,23800,0.000475447,50058188.059,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p1024,1024,MULTIPLYPRODUCTSCANNING,23800,0.002940828,8092958.907,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,23800,0.098067535,242689.898,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,23800,0.101478577,234532.260,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.022585763,553445.988,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.006470011,1931990.509,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.009996731,1250408.779,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000117568,425285792.052,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,23800,0.002692756,8838528.989,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,23800,0.004638393,5131087.488,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p1024,1024,MONTGOMERYMULTIPLICATION,23800,0.000541122,43982674.446,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p1024,1024,MONTGOMERYMULTIPLICATION,23800,0.002479736,9597796.816,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,p1024,1024,MONTGOMERYMULTIPLICATION,23800,0.000178142,133601295.552,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,p1024,1024,MONTGOMERYMULTIPLICATION,23800,0.002096462,11352459.181,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p1024,1024,MONTGOMERYMULTIPLICATION,23800,0.000231233,102926503.311,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p1024,1024,MONTGOMERYMULTIPLICATION,23800,0.002150271,11068371.860,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p1024,1024,MONTGOMERYMULTIPLICATION,23800,0.000182160,130654243.307,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p1024,1024,MONTGOMERYMULTIPLICATION,23800,0.002112821,11264561.872,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p1024,1024,MONTGOMERYMULTIPLICATION,23800,0.000163713,145376355.906,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p1024,1024,MONTGOMERYMULTIPLICATION,23800,0.002085013,11414797.758,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p1024,1024,MONTGOMERYMULTIPLICATION,23800,0.000118611,200656069.813,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p1024,1024,MONTGOMERYMULTIPLICATION,23800,0.002060968,11547971.839,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,23800,0.384605262,61881.629,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,23800,0.368990515,64500.303,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,p1024,1024,COMPARE,12500,0.000213009,58683010.627,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,p1024,1024,COMPARE,12500,0.005431009,2301598.105,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,p1024,1024,COMPARE,12500,0.005387802,2320055.702,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,p1024,1024,COMPARE,50000,0.000072192,692597517.730,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p1024,1024,COMPARE,23800,0.000234855,101339178.824,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p1024,1024,COMPARE,23800,0.002133711,11154276.518,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p1024,1024,COMPARE,23800,0.000124262,191530520.112,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p1024,1024,COMPARE,23800,0.002052882,11593456.147,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p1024,1024,COMPARE,23800,0.000069443,342726455.276,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p1024,1024,COMPARE,23800,0.001997832,11912912.353,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p1024,1024,COMPARE,23800,0.000071660,332123224.028,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p1024,1024,COMPARE,23800,0.001614480,14741590.324,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p1024,1024,COMPARE,23800,0.000014673,1622028271.101,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p1024,1024,COMPARE,23800,0.001947745,12219260.253,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p1024,1024,COMPARE,23800,0.000015720,1514014776.420,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p1024,1024,COMPARE,23800,0.001979600,12022632.615,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,p1024,1024,COMPARE,23800,0.000573308,41513480.946,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,p1024,1024,COMPARE,23800,0.003068317,7756694.794,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,p1024,1024,REDUCE,1562,0.000068179,22910325.300,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,p1024,1024,REDUCE,1562,0.004627584,337541.144,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,p1024,1024,REDUCE,1562,0.007240391,215734.209,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,p1024,1024,REDUCE,50000,0.000073856,676993067.591,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p1024,1024,REDUCE,23800,0.003616049,6581769.301,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p1024,1024,REDUCE,23800,0.005515854,4314835.117,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p1024,1024,REDUCE,23800,0.001991713,11949513.083,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p1024,1024,REDUCE,23800,0.003919684,6071918.255,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p1024,1024,REDUCE,23800,0.000499390,47658136.667,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p1024,1024,REDUCE,23800,0.002414217,9858267.854,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p1024,1024,REDUCE,23800,0.000503550,47264395.979,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p1024,1024,REDUCE,23800,0.002439538,9755945.525,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p1024,1024,REDUCE,23800,0.000451583,52703495.605,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p1024,1024,REDUCE,23800,0.002375909,10017220.100,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p1024,1024,REDUCE,23800,0.000438330,54296962.439,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p1024,1024,REDUCE,23800,0.002367703,10051937.353,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,p1024,1024,REDUCE,23800,0.281765078,84467.529,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,p1024,1024,REDUCE,23800,0.310438358,76665.784,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,p1024,1024,MODMUL,781,0.000761372,1025779.726,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,p1024,1024,MODMUL,781,0.007160935,109063.969,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,p1024,1024,MODMUL,781,0.007378279,105851.239,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,p1024,1024,MODMUL,50000,0.000590464,84679167.570,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p1024,1024,MODMUL,23800,0.015260006,1559632.384,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p1024,1024,MODMUL,23800,0.017213023,1382674.053,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p1024,1024,MODMUL,23800,0.005292159,4497219.343,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p1024,1024,MODMUL,23800,0.007193802,3308403.574,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p1024,1024,MODMUL,23800,0.002289159,10396832.426,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p1024,1024,MODMUL,23800,0.004231951,5623883.995,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p1024,1024,MODMUL,23800,0.001754098,13568226.268,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p1024,1024,MODMUL,23800,0.003690270,6449392.984,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p1024,1024,MODMUL,23800,0.002237577,10636505.511,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p1024,1024,MODMUL,23800,0.004163756,5715992.916,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p1024,1024,MODMUL,23800,0.001677779,14185420.002,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p1024,1024,MODMUL,23800,0.003575755,6655936.471,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,p1024,1024,MODMUL,23800,1.002867103,23731.958,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,p1024,1024,MODMUL,23800,1.028748281,23134.911,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,p1024,1024,MODEXP,195,0.069792181,2794.009,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,p1024,1024,MODEXP,195,0.009088467,21455.764,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,p1024,1024,MODEXP,195,0.008378872,23272.822,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,p1024,1024,MODEXP,50000,2.783334494,17964.064,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p1024,1024,MODEXP,23800,3.549183641,6705.767,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p1024,1024,MODEXP,23800,3.557758815,6689.605,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p1024,1024,MODEXP,23800,0.499028730,47692.645,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p1024,1024,MODEXP,23800,0.501072490,47498.118,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p1024,1024,MODEXP,23800,0.246561968,96527.458,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p1024,1024,MODEXP,23800,0.249254238,95484.836,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p1024,1024,MODEXP,23800,0.147317810,161555.483,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p1024,1024,MODEXP,23800,0.149658519,159028.702,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p1024,1024,MODEXP,23800,0.248867422,95633.249,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p1024,1024,MODEXP,23800,0.251249981,94726.375,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p1024,1024,MODEXP,23800,0.146711327,162223.330,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p1024,1024,MODEXP,23800,0.149226660,159488.928,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,p1024,1024,MODEXP,23800,0.000000000,inf,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,CPU,w8,p1024,1024,MODEXP,23800,0.000000000,inf,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,p1024,1024,EXPONENTIATION,195,0.015469211,12605.685,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,p1024,1024,EXPONENTIATION,195,0.007647313,25499.152,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,p1024,1024,EXPONENTIATION,195,0.098680849,1976.067,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p1024,1024,EXPONENTIATION,23800,3.483985253,6831.257,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p1024,1024,EXPONENTIATION,23800,3.493376620,6812.893,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p1024,1024,EXPONENTIATION,23800,0.710220582,33510.716,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p1024,1024,EXPONENTIATION,23800,0.725062718,32824.747,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p1024,1024,EXPONENTIATION,23800,0.206978664,114987.698,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p1024,1024,EXPONENTIATION,23800,0.209453136,113629.237,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p1024,1024,EXPONENTIATION,23800,0.206674739,115156.792,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p1024,1024,EXPONENTIATION,23800,0.208903888,113927.990,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p1024,1024,EXPONENTIATION,23800,0.207765154,114552.414,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p1024,1024,EXPONENTIATION,23800,0.209451463,113630.144,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p1024,1024,EXPONENTIATION,23800,0.211258887,112657.983,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p1024,1024,EXPONENTIATION,23800,0.213197817,111633.413,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,p1024,1024,DIVIDE,1562,0.000135712,11509679.412,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,p1024,1024,DIVIDE,1562,0.005364892,291152.201,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,p1024,1024,DIVIDE,1562,0.005844803,267245.959,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,p1024,1024,DIVIDE,50000,0.000142240,351518560.180,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p1024,1024,DIVIDE,23800,0.088333622,269433.082,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p1024,1024,DIVIDE,23800,0.090943356,261701.360,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p1024,1024,DIVIDE,23800,0.013170727,1807037.674,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p1024,1024,DIVIDE,23800,0.015641786,1521565.406,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p1024,1024,DIVIDE,23800,0.001998748,11907453.086,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p1024,1024,DIVIDE,23800,0.004488988,5301863.316,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p1024,1024,DIVIDE,23800,0.001865113,12760623.070,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p1024,1024,DIVIDE,23800,0.004320965,5508028.689,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p1024,1024,DIVIDE,23800,0.001986120,11983160.969,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p1024,1024,DIVIDE,23800,0.004388905,5422764.743,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p1024,1024,DIVIDE,23800,0.001887794,12607310.108,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p1024,1024,DIVIDE,23800,0.004349969,5471303.533,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,p1024,1024,ISQRT,390,0.000318845,1223165.618,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,p1024,1024,ISQRT,390,0.004967843,78504.894,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p1024,1024,ISQRT,23800,1.291496363,18428.236,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p1024,1024,ISQRT,23800,1.292848145,18408.968,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p1024,1024,ISQRT,23800,0.230506817,103250.743,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p1024,1024,ISQRT,23800,0.231834534,102659.425,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p1024,1024,ISQRT,23800,0.037914281,627731.807,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p1024,1024,ISQRT,23800,0.039792269,598106.130,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p1024,1024,ISQRT,23800,0.035429465,671757.248,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p1024,1024,ISQRT,23800,0.037374141,636803.931,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p1024,1024,ISQRT,23800,0.037937678,627344.669,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p1024,1024,ISQRT,23800,0.039944653,595824.428,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p1024,1024,ISQRT,23800,0.035280639,674590.959,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p1024,1024,ISQRT,23800,0.037273508,638523.210,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,p1024,1024,MODMUL_R2,12500,0.008866171,1409853.194,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,p1024,1024,MODMUL_R2,12500,0.006458579,1935410.221,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,p1024,1024,MODMUL_R2,12500,0.006034624,2071380.147,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,p1024,1024,MODMUL_R2,50000,0.000221184,226056134.259,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p1024,1024,MODMUL_R2,23800,0.003852978,6177040.389,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p1024,1024,MODMUL_R2,23800,0.005781116,4116852.204,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p1024,1024,MODMUL_R2,23800,0.000615188,38687365.603,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p1024,1024,MODMUL_R2,23800,0.002515155,9462637.903,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p1024,1024,MODMUL_R2,23800,0.000357917,66495855.334,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p1024,1024,MODMUL_R2,23800,0.002333110,10200975.512,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p1024,1024,MODMUL_R2,23800,0.000270304,88048923.505,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p1024,1024,MODMUL_R2,23800,0.002177597,10929480.906,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p1024,1024,MODMUL_R2,23800,0.000293388,81121236.898,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p1024,1024,MODMUL_R2,23800,0.002197212,10831908.436,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p1024,1024,MODMUL_R2,23800,0.000213232,111615574.220,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p1024,1024,MODMUL_R2,23800,0.002132877,11158638.073,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,p2048,2048,ADD,6250,0.000322407,19385428.474,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,p2048,2048,ADD,6250,0.007144323,874820.516,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,p2048,2048,ADD,6250,0.005975146,1045999.572,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,p2048,2048,ADD,50000,0.000137216,364388992.537,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p2048,2048,ADD,23800,0.000846895,28102659.889,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p2048,2048,ADD,23800,0.004352347,5468314.000,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p2048,2048,ADD,23800,0.000425834,55890357.182,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p2048,2048,ADD,23800,0.003997568,5953619.740,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,p2048,2048,ADD,23800,0.000268830,88531790.348,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,p2048,2048,ADD,23800,0.003744729,6355600.330,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p2048,2048,ADD,23800,0.000220130,108117833.211,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p2048,2048,ADD,23800,0.003718311,6400755.669,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p2048,2048,ADD,23800,0.000218680,108834763.214,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p2048,2048,ADD,23800,0.003827393,6218331.814,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p2048,2048,ADD,23800,0.000092917,256142242.693,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p2048,2048,ADD,23800,0.003590595,6628427.921,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p2048,2048,ADD,23800,0.000093494,254561582.365,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p2048,2048,ADD,23800,0.003612779,6587726.364,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,p2048,2048,SUBTRACT,6250,0.000283981,22008498.553,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,p2048,2048,SUBTRACT,6250,0.005354982,1167137.472,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,p2048,2048,SUBTRACT,6250,0.010579221,590780.742,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,p2048,2048,SUBTRACT,50000,0.000135968,367733584.373,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p2048,2048,SUBTRACT,23800,0.000845929,28134744.172,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p2048,2048,SUBTRACT,23800,0.004475484,5317860.431,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p2048,2048,SUBTRACT,23800,0.000426825,55760539.846,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p2048,2048,SUBTRACT,23800,0.003900891,6101169.954,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,p2048,2048,SUBTRACT,23800,0.000272219,87429733.610,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,p2048,2048,SUBTRACT,23800,0.003900355,6102009.089,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p2048,2048,SUBTRACT,23800,0.000219571,108393215.253,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p2048,2048,SUBTRACT,23800,0.003762433,6325694.308,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p2048,2048,SUBTRACT,23800,0.000219731,108314424.420,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p2048,2048,SUBTRACT,23800,0.003707795,6418909.459,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p2048,2048,SUBTRACT,23800,0.000093221,255306736.179,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p2048,2048,SUBTRACT,23800,0.003625784,6564097.460,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p2048,2048,SUBTRACT,23800,0.000093520,254490600.760,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p2048,2048,SUBTRACT,23800,0.003585159,6638478.520,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,p2048,2048,ADDMOD,6250,0.000732917,8527564.666,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,p2048,2048,ADDMOD,6250,0.006799917,919128.838,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,p2048,2048,ADDMOD,6250,0.007409493,843512.508,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,p2048,2048,ADDMOD,50000,0.000135808,368166823.751,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p2048,2048,ADDMOD,23800,0.001374889,17310486.258,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p2048,2048,ADDMOD,23800,0.004934165,4823511.469,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p2048,2048,ADDMOD,23800,0.000552442,43081415.204,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p2048,2048,ADDMOD,23800,0.004043288,5886298.386,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,p2048,2048,ADDMOD,23800,0.000336972,70629000.006,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,p2048,2048,ADDMOD,23800,0.003823682,6224367.417,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p2048,2048,ADDMOD,23800,0.000237949,100021547.977,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p2048,2048,ADDMOD,23800,0.003734573,6372883.424,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p2048,2048,ADDMOD,23800,0.000233031,102132183.702,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p2048,2048,ADDMOD,23800,0.003723816,6391293.207,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p2048,2048,ADDMOD,23800,0.000082355,288992801.049,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p2048,2048,ADDMOD,23800,0.003592705,6624535.201,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p2048,2048,ADDMOD,23800,0.000082155,289697158.142,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p2048,2048,ADDMOD,23800,0.003631897,6553048.786,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,p2048,2048,SUBTRACTMOD,6250,0.000712079,8777116.072,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,p2048,2048,SUBTRACTMOD,6250,0.006922950,902794.324,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,p2048,2048,SUBTRACTMOD,6250,0.008611652,725760.859,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,p2048,2048,SUBTRACTMOD,50000,0.000135616,368688060.406,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p2048,2048,SUBTRACTMOD,23800,0.001426787,16680834.680,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p2048,2048,SUBTRACTMOD,23800,0.004968369,4790304.507,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p2048,2048,SUBTRACTMOD,23800,0.000742715,32044599.530,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p2048,2048,SUBTRACTMOD,23800,0.004225381,5632628.533,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,p2048,2048,SUBTRACTMOD,23800,0.000365298,65152326.628,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,p2048,2048,SUBTRACTMOD,23800,0.003845050,6189776.466,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p2048,2048,SUBTRACTMOD,23800,0.000234656,101425049.259,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p2048,2048,SUBTRACTMOD,23800,0.003750793,6345325.342,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p2048,2048,SUBTRACTMOD,23800,0.000233324,102004176.798,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p2048,2048,SUBTRACTMOD,23800,0.003754890,6338401.293,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p2048,2048,SUBTRACTMOD,23800,0.000088129,270059500.792,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p2048,2048,SUBTRACTMOD,23800,0.003609356,6593974.096,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p2048,2048,SUBTRACTMOD,23800,0.000084941,280193579.422,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p2048,2048,SUBTRACTMOD,23800,0.003579284,6649374.448,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.005641929,1107777.170,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.005793655,1078766.373,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.007038496,887973.798,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,23800,0.238827128,99653.671,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,23800,0.244030693,97528.715,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p2048,2048,MULTIPLYOPERANDSCANNING,23800,0.060799957,391447.647,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p2048,2048,MULTIPLYOPERANDSCANNING,23800,0.065475832,363492.901,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,p2048,2048,MULTIPLYOPERANDSCANNING,23800,0.017819636,1335605.281,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,p2048,2048,MULTIPLYOPERANDSCANNING,23800,0.020403249,1166480.877,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p2048,2048,MULTIPLYOPERANDSCANNING,23800,0.000863523,27561520.309,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p2048,2048,MULTIPLYOPERANDSCANNING,23800,0.005474617,4347336.258,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p2048,2048,MULTIPLYOPERANDSCANNING,23800,0.000845175,28159840.673,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p2048,2048,MULTIPLYOPERANDSCANNING,23800,0.005453613,4364079.276,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p2048,2048,MULTIPLYOPERANDSCANNING,23800,0.000733485,32447835.550,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p2048,2048,MULTIPLYOPERANDSCANNING,23800,0.005397185,4409705.993,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p2048,2048,MULTIPLYOPERANDSCANNING,23800,0.000633415,37574112.090,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p2048,2048,MULTIPLYOPERANDSCANNING,23800,0.005238651,4543154.077,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.005600937,1115884.682,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.007313909,854536.198,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.006862883,910695.987,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.000294912,169542100.694,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,23800,0.102588222,231995.443,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,23800,0.107415550,221569.410,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p2048,2048,MULTIPLYPRODUCTSCANNING,23800,0.025957723,916875.485,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p2048,2048,MULTIPLYPRODUCTSCANNING,23800,0.030628877,777044.481,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,p2048,2048,MULTIPLYPRODUCTSCANNING,23800,0.006544335,3636732.989,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,p2048,2048,MULTIPLYPRODUCTSCANNING,23800,0.011092527,2145588.588,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p2048,2048,MULTIPLYPRODUCTSCANNING,23800,0.008233148,2890753.445,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p2048,2048,MULTIPLYPRODUCTSCANNING,23800,0.012977834,1833896.119,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p2048,2048,MULTIPLYPRODUCTSCANNING,23800,0.008237731,2889145.019,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p2048,2048,MULTIPLYPRODUCTSCANNING,23800,0.012817447,1856844.024,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p2048,2048,MULTIPLYPRODUCTSCANNING,23800,0.003026528,7863795.665,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p2048,2048,MULTIPLYPRODUCTSCANNING,23800,0.007654309,3109359.598,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p2048,2048,MULTIPLYPRODUCTSCANNING,23800,0.002916803,8159619.569,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p2048,2048,MULTIPLYPRODUCTSCANNING,23800,0.007538988,3156922.476,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.023316489,268050.648,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.007828070,798408.818,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.007964952,784687.725,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000536576,93183444.656,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,23800,0.086609321,274797.212,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,23800,0.090548170,262843.524,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p2048,2048,MONTGOMERYMULTIPLICATION,23800,0.002289185,10396713.991,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p2048,2048,MONTGOMERYMULTIPLICATION,23800,0.005780973,4116954.010,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,p2048,2048,MONTGOMERYMULTIPLICATION,23800,0.000404470,58842416.209,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32,p2048,2048,MONTGOMERYMULTIPLICATION,23800,0.003864651,6158383.350,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p2048,2048,MONTGOMERYMULTIPLICATION,23800,0.000863452,27563779.634,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p2048,2048,MONTGOMERYMULTIPLICATION,23800,0.004444010,5355523.394,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p2048,2048,MONTGOMERYMULTIPLICATION,23800,0.000553534,42996463.229,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p2048,2048,MONTGOMERYMULTIPLICATION,23800,0.003246694,7330534.057,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p2048,2048,MONTGOMERYMULTIPLICATION,23800,0.000743366,32016516.767,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p2048,2048,MONTGOMERYMULTIPLICATION,23800,0.004318882,5510685.685,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p2048,2048,MONTGOMERYMULTIPLICATION,23800,0.000580891,40971539.490,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p2048,2048,MONTGOMERYMULTIPLICATION,23800,0.004066759,5852326.443,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,p2048,2048,COMPARE,6250,0.000073068,85536943.000,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,p2048,2048,COMPARE,6250,0.007314211,854500.944,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,p2048,2048,COMPARE,6250,0.005128564,1218664.754,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,p2048,2048,COMPARE,50000,0.000136672,365839381.878,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p2048,2048,COMPARE,23800,0.000466883,50976348.883,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p2048,2048,COMPARE,23800,0.003990025,5964875.215,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p2048,2048,COMPARE,23800,0.000238433,99818391.508,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p2048,2048,COMPARE,23800,0.003738563,6366082.297,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p2048,2048,COMPARE,23800,0.000153679,154868329.851,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p2048,2048,COMPARE,23800,0.003768266,6315902.028,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p2048,2048,COMPARE,23800,0.000126805,187689528.271,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p2048,2048,COMPARE,23800,0.003645156,6529212.936,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p2048,2048,COMPARE,23800,0.000040345,589913559.815,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p2048,2048,COMPARE,23800,0.003566482,6673242.513,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p2048,2048,COMPARE,23800,0.000039847,597283084.484,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p2048,2048,COMPARE,23800,0.003574484,6658303.636,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,p2048,2048,REDUCE,781,0.000048934,15960267.679,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,p2048,2048,REDUCE,781,0.007144363,109316.952,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,p2048,2048,REDUCE,781,0.004858198,160759.191,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,p2048,2048,REDUCE,50000,0.000137504,363625785.432,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p2048,2048,REDUCE,23800,0.774783386,30718.263,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p2048,2048,REDUCE,23800,0.779054722,30549.844,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p2048,2048,REDUCE,23800,0.008165640,2914652.024,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p2048,2048,REDUCE,23800,0.011707272,2032924.508,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p2048,2048,REDUCE,23800,0.002288397,10400293.597,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p2048,2048,REDUCE,23800,0.005789465,4110915.396,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p2048,2048,REDUCE,23800,0.001613586,14749754.216,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p2048,2048,REDUCE,23800,0.005158069,4614129.753,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p2048,2048,REDUCE,23800,0.001882422,12643287.600,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p2048,2048,REDUCE,23800,0.005389628,4415889.185,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p2048,2048,REDUCE,23800,0.001972067,12068557.290,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p2048,2048,REDUCE,23800,0.005453884,4363862.416,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,p2048,2048,MODMUL,390,0.000980124,397908.697,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,p2048,2048,MODMUL,390,0.007270919,53638.339,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,p2048,2048,MODMUL,390,0.007445336,52381.784,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,p2048,2048,MODMUL,50000,0.002639680,18941689.902,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p2048,2048,MODMUL,23800,1.280775795,18582.487,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p2048,2048,MODMUL,23800,1.284676414,18526.066,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p2048,2048,MODMUL,23800,0.023721653,1003302.762,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p2048,2048,MODMUL,23800,0.027294918,871957.191,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p2048,2048,MODMUL,23800,0.011835335,2010927.477,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p2048,2048,MODMUL,23800,0.015356987,1549783.192,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p2048,2048,MODMUL,23800,0.007078481,3362303.341,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p2048,2048,MODMUL,23800,0.010579614,2249609.491,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p2048,2048,MODMUL,23800,0.011479571,2073248.208,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p2048,2048,MODMUL,23800,0.012806304,1858459.738,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p2048,2048,MODMUL,23800,0.008917650,2668864.493,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p2048,2048,MODMUL,23800,0.012496228,1904574.745,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,p2048,2048,MODEXP,97,0.243860572,397.768,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,p2048,2048,MODEXP,97,0.064705384,1499.102,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,p2048,2048,MODEXP,97,0.100806141,962.243,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,p2048,2048,MODEXP,50000,3.654630423,13681.274,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p2048,2048,MODEXP,23800,174.304460130,136.543,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p2048,2048,MODEXP,23800,174.352885406,136.505,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p2048,2048,MODEXP,23800,25.934900045,917.682,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p2048,2048,MODEXP,23800,26.202306868,908.317,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p2048,2048,MODEXP,23800,2.340470542,10168.895,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p2048,2048,MODEXP,23800,2.354775711,10107.120,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p2048,2048,MODEXP,23800,4.438617268,5362.030,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p2048,2048,MODEXP,23800,4.455675978,5341.502,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p2048,2048,MODEXP,23800,2.335415437,10190.906,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p2048,2048,MODEXP,23800,2.337203730,10183.109,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p2048,2048,MODEXP,23800,4.492484751,5297.736,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p2048,2048,MODEXP,23800,4.487108624,5304.084,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,p2048,2048,EXPONENTIATION,97,0.044613049,2174.252,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,p2048,2048,EXPONENTIATION,97,0.009603532,10100.451,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,p2048,2048,EXPONENTIATION,97,0.127117927,763.071,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p2048,2048,EXPONENTIATION,23800,42.343252960,562.073,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p2048,2048,EXPONENTIATION,23800,42.367453668,561.752,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p2048,2048,EXPONENTIATION,23800,8.675393295,2743.391,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p2048,2048,EXPONENTIATION,23800,8.706604424,2733.557,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p2048,2048,EXPONENTIATION,23800,1.986578079,11980.400,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p2048,2048,EXPONENTIATION,23800,2.009737312,11842.344,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p2048,2048,EXPONENTIATION,23800,2.190963780,10862.799,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p2048,2048,EXPONENTIATION,23800,2.214176372,10748.918,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p2048,2048,EXPONENTIATION,23800,1.972008512,12068.913,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p2048,2048,EXPONENTIATION,23800,1.989379230,11963.531,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p2048,2048,EXPONENTIATION,23800,2.163959480,10998.358,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p2048,2048,EXPONENTIATION,23800,2.178982468,10922.529,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,p2048,2048,DIVIDE,781,0.000081458,9587802.646,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,p2048,2048,DIVIDE,781,0.008364628,93369.364,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,p2048,2048,DIVIDE,781,0.005101572,153090.073,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,p2048,2048,DIVIDE,50000,0.000204128,244944348.644,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p2048,2048,DIVIDE,23800,1.830125789,13004.571,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p2048,2048,DIVIDE,23800,1.834940291,12970.449,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p2048,2048,DIVIDE,23800,0.373037573,63800.544,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p2048,2048,DIVIDE,23800,0.377864736,62985.502,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p2048,2048,DIVIDE,23800,0.025915403,918372.756,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p2048,2048,DIVIDE,23800,0.030543090,779226.979,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p2048,2048,DIVIDE,23800,0.023388142,1017609.683,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p2048,2048,DIVIDE,23800,0.028129502,846086.782,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p2048,2048,DIVIDE,23800,0.025035071,950666.360,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p2048,2048,DIVIDE,23800,0.029669736,802164.197,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p2048,2048,DIVIDE,23800,0.023645692,1006525.838,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p2048,2048,DIVIDE,23800,0.028441701,836799.456,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,p2048,2048,ISQRT,195,0.000175593,1110522.912,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,p2048,2048,ISQRT,195,0.005246395,37168.381,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p2048,2048,ISQRT,23800,23.560343521,1010.172,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p2048,2048,ISQRT,23800,23.566774937,1009.896,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p2048,2048,ISQRT,23800,8.883395259,2679.156,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p2048,2048,ISQRT,23800,8.907601469,2671.875,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p2048,2048,ISQRT,23800,0.092941134,256076.067,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p2048,2048,ISQRT,23800,0.096877842,245670.213,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p2048,2048,ISQRT,23800,0.079525999,299273.198,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p2048,2048,ISQRT,23800,0.083193505,286080.026,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p2048,2048,ISQRT,23800,0.092277610,257917.386,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p2048,2048,ISQRT,23800,0.095993189,247934.258,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p2048,2048,ISQRT,23800,0.079643156,298832.960,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p2048,2048,ISQRT,23800,0.083271232,285812.993,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-1t,p2048,2048,MODMUL_R2,6250,0.015439611,404802.947,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,gmp-nt,p2048,2048,MODMUL_R2,6250,0.007111495,878858.818,0
library,Intel(R) Xeon(R) CPU E5-2673 v4 @ 2.30GHz,host-cpu,openssl-nt,p2048,2048,MODMUL_R2,6250,0.009512466,657032.578,0
library,NVIDIA GeForce RTX 2060 SUPER,gpu,cgbn,p2048,2048,MODMUL_R2,50000,0.000818144,61113935.933,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p2048,2048,MODMUL_R2,23800,0.065629751,362640.413,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w8,p2048,2048,MODMUL_R2,23800,0.068964669,345104.245,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p2048,2048,MODMUL_R2,23800,0.002315641,10277929.992,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w16,p2048,2048,MODMUL_R2,23800,0.005830785,4081782.974,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p2048,2048,MODMUL_R2,23800,0.001204680,19756280.780,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-opt,p2048,2048,MODMUL_R2,23800,0.004744113,5016743.955,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p2048,2048,MODMUL_R2,23800,0.000899296,26465152.301,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-o64,p2048,2048,MODMUL_R2,23800,0.004391058,5420106.188,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p2048,2048,MODMUL_R2,23800,0.001075003,22139476.590,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il,p2048,2048,MODMUL_R2,23800,0.004585522,5190248.530,0
opencl-kernel,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p2048,2048,MODMUL_R2,23800,0.000814115,29234193.249,0
opencl-e2e,NVIDIA GeForce RTX 2060 SUPER,GPU,w32-il64,p2048,2048,MODMUL_R2,23800,0.004313874,5517082.714,0
```
