# MPA-OpenCL benchmark report - NVIDIA GB10


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

### Device 0 - NVIDIA GB10 (GPU)

| Property | Value |
|---|---|
| Model | NVIDIA GB10 |
| Type | GPU |
| Vendor | NVIDIA Corporation |
| Device memory | 121.63 GiB |
| Max single allocation | 30.41 GiB |
| Local memory | 48 KiB |
| Global cache | 1536 KiB |
| Compute units | 48 |
| Max clock | 2418 MHz |
| Max work-group size | 1024 |
| OpenCL version | OpenCL 3.0 CUDA |
| Driver | 580.173.02 |

### Device 1 - cpu--0xd87 (CPU)

| Property | Value |
|---|---|
| Model | cpu--0xd87 |
| Type | CPU |
| Vendor | ARM |
| Device memory | 119.63 GiB |
| Max single allocation | 32.00 GiB |
| Local memory | 512 KiB |
| Global cache | 8192 KiB |
| Compute units | 20 |
| Max clock | 2808 MHz |
| Max work-group size | 4096 |
| OpenCL version | OpenCL 3.0 PoCL HSTR: cpu-aarch64-unknown-linux-gnu-(null) |
| Driver | 5.0+debian |

### Host

| Property | Value |
|---|---|
| CPU | unknown |
| Logical cores | 20 |
| OpenMP threads used | 20 |
| RAM | 121.6 GB |
| OS | Ubuntu 24.04.4 LTS |
| Kernel | 6.17.0-1029-nvidia |
| Arch | aarch64 |
| GMP | 6.3.0 |
| OpenSSL | OpenSSL 3.0.13 30 Jan 2024 |

## 2. Method

- Workload auto-sized from the device and host: --min-items from 700 x compute units, --items from ten times that capped by host RAM. Either flag, given explicitly, overrides its half.
- Base workload 50000 items, scaled down per operator by its cost weight and by modulus size. Device rows honour --min-items (33600) so the GPU is not left idle; the CPU libraries keep the smaller count because a full-width MODEXP there costs minutes. Both counts appear in every row as dev/cpu, and throughput is per-second so they remain comparable.
- 5 timed repetitions, **minimum** reported. Two untimed warm-up launches precede them.
- `kernel` times `clEnqueueNDRangeKernel` + `clFinish` only. `e2e` adds the host->device operand writes and the device->host result read.
- Every OpenCL device runs the same kernels on the same operands, so GPU and CPU-OpenCL columns are directly comparable.
- CPU library baselines (GMP, OpenSSL) run those same operands, with temporaries preallocated outside the timed region, so the figure is the arithmetic and not marshalling. The generator is reseeded per modulus and operation so every backend sees identical inputs.
- OpenSSL rows time the nearest BN primitive, which is not always semantically identical (its Montgomery routine expects Montgomery-domain inputs); they measure comparable work, not identical results. Correctness is judged against GMP only.
- Every device cell is checked word-for-word against GMP before it is timed. A cell that mismatches is reported and excluded from the speedup tables.
- Total wall time 3173.5 s.

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

**All configurations correct** - 485 configurations, 0 problems.

## 4. Throughput per device

Operations per second, higher is better. Kernel-only timings.

### Device 0 - NVIDIA GB10 (GPU)

#### secp256k1 (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | 664.62 M | 1.20 G | 1.94 G | 2.09 G | 2.07 G | 4.56 G | 4.46 G | 93.99 M |
| SUBTRACT | 50000 / 50000 | 662.37 M | 1.22 G | 1.98 G | 2.09 G | 2.03 G | 4.60 G | 4.66 G | 131.28 M |
| ADDMOD | 50000 / 50000 | 427.03 M | 820.64 M | 1.56 G | 2.44 G | 3.07 G | 4.90 G | 4.84 G | 39.26 M |
| SUBTRACTMOD | 50000 / 50000 | 427.50 M | 815.50 M | 1.57 G | 2.48 G | 2.98 G | 4.94 G | 4.84 G | 50.37 M |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | 18.04 M | 60.15 M | 209.41 M | 1.02 G | 1.03 G | 3.49 G | 3.76 G | 70.06 M |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | 157.53 M | 486.31 M | 951.87 M | 968.09 M | 1.05 G | 2.08 G | 1.97 G | 70.64 M |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | 286.96 M | 1.07 G | 3.04 G | 2.74 G | 2.97 G | 2.79 G | 3.23 G | 9.08 M |
| COMPARE | 50000 / 50000 | 686.36 M | 1.26 G | - | 2.89 G | 3.80 G | 5.48 G | 5.44 G | 142.23 M |
| REDUCE | 33600 / 6250 | 203.08 M | 320.71 M | - | 972.70 M | 959.78 M | 975.38 M | 798.76 M | 86.08 M |
| MODMUL | 33600 / 3125 | 68.33 M | 126.70 M | - | 318.57 M | 320.81 M | 323.08 M | 276.06 M | 17.29 M |
| MODEXP | 33600 / 781 | 1.35 M | 8.88 M | - | 10.98 M | 19.08 M | 11.13 M | 18.58 M | 121.35 k |
| EXPONENTIATION | 33600 / 781 | 639.84 k | 2.13 M | - | 48.22 M | 50.18 M | 43.32 M | 43.90 M | 486.83 k |
| DIVIDE | 33600 / 6250 | 90.32 M | 142.39 M | - | 363.83 M | 374.47 M | 379.07 M | 408.71 M | 52.62 M |
| ISQRT | 33600 / 1562 | 7.97 M | 11.03 M | - | 52.95 M | 47.61 M | 57.23 M | 37.55 M | 26.98 M |
| MODMUL_R2 | 50000 / 50000 | 249.46 M | 1.12 G | - | 1.77 G | 1.93 G | 1.81 G | 2.08 G | 17.24 M |

#### rsa256(composite) (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | 664.05 M | 1.22 G | 1.97 G | 2.09 G | 2.07 G | 4.65 G | 4.36 G | 92.85 M |
| SUBTRACT | 50000 / 50000 | 655.00 M | 1.21 G | 1.94 G | 2.07 G | 2.07 G | 4.58 G | 4.30 G | 131.81 M |
| ADDMOD | 50000 / 50000 | 471.34 M | 907.90 M | 1.76 G | 2.51 G | 3.05 G | 4.79 G | 4.73 G | 46.85 M |
| SUBTRACTMOD | 50000 / 50000 | 427.73 M | 817.63 M | 1.54 G | 2.48 G | 2.96 G | 4.82 G | 4.71 G | 50.46 M |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | 17.99 M | 61.97 M | 209.28 M | 1.04 G | 1.03 G | 3.43 G | 3.68 G | 70.66 M |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | 157.58 M | 485.63 M | 957.71 M | 948.69 M | 1.06 G | 2.12 G | 1.93 G | 70.20 M |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | 289.65 M | 1.10 G | 3.01 G | 2.81 G | 3.02 G | 2.79 G | 3.31 G | 9.09 M |
| COMPARE | 50000 / 50000 | 694.61 M | 1.29 G | - | 3.10 G | 3.60 G | 5.31 G | 5.31 G | 143.20 M |
| REDUCE | 33600 / 6250 | 199.79 M | 320.02 M | - | 968.63 M | 952.84 M | 969.98 M | 790.94 M | 58.28 M |
| MODMUL | 33600 / 3125 | 68.43 M | 126.71 M | - | 319.39 M | 321.25 M | 322.88 M | 276.31 M | 17.27 M |
| MODEXP | 33600 / 781 | 1.35 M | 8.89 M | - | 11.00 M | 19.16 M | 11.15 M | 18.66 M | 127.44 k |
| EXPONENTIATION | 33600 / 781 | 646.56 k | 2.19 M | - | 47.60 M | 49.76 M | 43.98 M | 44.02 M | 488.08 k |
| DIVIDE | 33600 / 6250 | 91.59 M | 137.73 M | - | 342.19 M | 349.94 M | 347.22 M | 382.02 M | 51.84 M |
| ISQRT | 33600 / 1562 | 7.98 M | 11.03 M | - | 52.73 M | 47.36 M | 57.27 M | 37.54 M | 27.07 M |
| MODMUL_R2 | 50000 / 50000 | 250.62 M | 1.17 G | - | 1.72 G | 2.03 G | 1.77 G | 2.08 G | 17.24 M |

#### brainpoolP512r1 (512-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 33600 / 25000 | 288.94 M | 569.72 M | 877.19 M | 861.03 M | 840.02 M | 2.34 G | 2.24 G | 86.96 M |
| SUBTRACT | 33600 / 25000 | 286.22 M | 555.56 M | 880.87 M | 887.95 M | 859.62 M | 2.24 G | 2.24 G | 126.90 M |
| ADDMOD | 33600 / 25000 | 213.81 M | 412.65 M | 744.68 M | 858.54 M | 886.82 M | 2.08 G | 2.09 G | 42.33 M |
| SUBTRACTMOD | 33600 / 25000 | 185.17 M | 361.82 M | 629.69 M | 871.39 M | 883.84 M | 2.12 G | 2.08 G | 44.14 M |
| MULTIPLYOPERANDSCANNING | 33600 / 25000 | 3.74 M | 12.82 M | 40.10 M | 248.08 M | 246.77 M | 594.06 M | 588.88 M | 29.22 M |
| MULTIPLYPRODUCTSCANNING | 33600 / 25000 | 26.38 M | 94.80 M | 312.13 M | 241.96 M | 241.88 M | 474.90 M | 521.86 M | 29.18 M |
| MONTGOMERYMULTIPLICATION | 33600 / 25000 | 82.97 M | 248.93 M | 891.72 M | 839.66 M | 951.52 M | 833.33 M | 1.18 G | 3.94 M |
| COMPARE | 33600 / 25000 | 284.79 M | 536.82 M | - | 1.09 G | 1.10 G | 3.25 G | 3.43 G | 140.49 M |
| REDUCE | 33600 / 3125 | 65.03 M | 84.25 M | - | 313.57 M | 330.04 M | 292.77 M | 307.69 M | 56.53 M |
| MODMUL | 33600 / 1562 | 21.83 M | 31.03 M | - | 85.01 M | 89.93 M | 72.37 M | 74.77 M | 8.17 M |
| MODEXP | 33600 / 390 | 120.26 k | 1.25 M | - | 1.37 M | 2.60 M | 1.29 M | 2.50 M | 28.01 k |
| EXPONENTIATION | 33600 / 390 | 84.01 k | 328.19 k | - | 1.05 M | 1.09 M | 1.01 M | 1.05 M | 157.30 k |
| DIVIDE | 33600 / 3125 | 31.10 M | 35.47 M | - | 106.24 M | 111.55 M | 116.38 M | 121.29 M | 47.85 M |
| ISQRT | 33600 / 781 | 1.51 M | 1.63 M | - | 8.38 M | 8.76 M | 6.60 M | 6.77 M | 14.67 M |
| MODMUL_R2 | 33600 / 25000 | 57.21 M | 280.04 M | - | 471.81 M | 615.31 M | 487.47 M | 708.01 M | 8.11 M |

#### p1024 (1024-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 33600 / 12500 | 97.20 M | 193.80 M | 351.88 M | 350.18 M | 350.64 M | 1.47 G | 1.43 G | 73.06 M |
| SUBTRACT | 33600 / 12500 | 103.44 M | 199.43 M | 353.00 M | 350.06 M | 351.82 M | 1.48 G | 1.41 G | 101.62 M |
| ADDMOD | 33600 / 12500 | 69.55 M | 137.04 M | 283.94 M | 306.98 M | 311.80 M | 1.26 G | 1.30 G | 29.03 M |
| SUBTRACTMOD | 33600 / 12500 | 71.69 M | 135.83 M | 276.83 M | 306.66 M | 314.04 M | 1.27 G | 1.30 G | 38.69 M |
| MULTIPLYOPERANDSCANNING | 33600 / 12500 | 922.63 k | 3.32 M | 13.14 M | 100.99 M | 97.73 M | 181.88 M | 179.39 M | 8.52 M |
| MULTIPLYPRODUCTSCANNING | 33600 / 12500 | 3.51 M | 13.69 M | 52.45 M | 52.09 M | 52.57 M | 115.20 M | 126.44 M | 8.51 M |
| MONTGOMERYMULTIPLICATION | 33600 / 12500 | 10.52 M | 74.68 M | 200.61 M | 180.37 M | 229.41 M | 269.47 M | 415.10 M | 1.26 M |
| COMPARE | 33600 / 12500 | 104.80 M | 225.28 M | - | 419.41 M | 411.21 M | 2.35 G | 2.36 G | 159.31 M |
| REDUCE | 33600 / 1562 | 12.22 M | 22.81 M | - | 86.68 M | 89.48 M | 105.92 M | 108.10 M | 81.29 M |
| MODMUL | 33600 / 781 | 2.78 M | 8.16 M | - | 20.06 M | 20.74 M | 20.62 M | 21.00 M | 2.98 M |
| MODEXP | 33600 / 195 | 5.99 k | 176.84 k | - | 170.35 k | 301.09 k | 170.33 k | 298.33 k | 4.89 k |
| EXPONENTIATION | 33600 / 195 | 9.93 k | 40.68 k | - | 143.15 k | 154.68 k | 147.06 k | 159.76 k | 40.67 k |
| DIVIDE | 33600 / 1562 | 653.45 k | 5.88 M | - | 26.57 M | 28.19 M | 27.21 M | 28.48 M | 42.82 M |
| ISQRT | 33600 / 390 | 78.57 k | 228.78 k | - | 1.18 M | 1.32 M | 1.23 M | 1.28 M | 6.89 M |
| MODMUL_R2 | 33600 / 12500 | 5.20 M | 91.33 M | - | 117.12 M | 156.95 M | 141.24 M | 206.14 M | 2.92 M |

#### p2048 (2048-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 33600 / 6250 | 43.64 M | 97.68 M | 188.73 M | 186.02 M | 176.61 M | 315.32 M | 256.56 M | 56.93 M |
| SUBTRACT | 33600 / 6250 | 44.72 M | 96.97 M | 189.87 M | 189.58 M | 176.12 M | 302.16 M | 257.76 M | 74.38 M |
| ADDMOD | 33600 / 6250 | 35.65 M | 72.74 M | 163.02 M | 132.68 M | 123.38 M | 273.55 M | 216.16 M | 24.95 M |
| SUBTRACTMOD | 33600 / 6250 | 33.06 M | 67.25 M | 140.41 M | 131.91 M | 117.69 M | 266.53 M | 218.86 M | 31.31 M |
| MULTIPLYOPERANDSCANNING | 33600 / 6250 | 222.02 k | 818.20 k | 3.16 M | 34.67 M | 33.64 M | 46.38 M | 44.62 M | 2.78 M |
| MULTIPLYPRODUCTSCANNING | 33600 / 6250 | 878.73 k | 3.45 M | 13.61 M | 13.43 M | 13.44 M | 21.88 M | 22.95 M | 2.78 M |
| MONTGOMERYMULTIPLICATION | 33600 / 6250 | 579.05 k | 8.72 M | 68.19 M | 52.68 M | 55.99 M | 61.22 M | 70.57 M | 421.62 k |
| COMPARE | 33600 / 6250 | 45.38 M | 104.37 M | - | 220.50 M | 212.51 M | 855.75 M | 602.92 M | 179.35 M |
| REDUCE | 33600 / 781 | 42.23 k | 5.00 M | - | 28.30 M | 32.76 M | 29.39 M | 35.10 M | 60.34 M |
| MODMUL | 33600 / 390 | 20.17 k | 1.21 M | - | 5.14 M | 2.64 M | 5.25 M | 2.67 M | 1.01 M |
| MODEXP | 33600 / 97 | 448.4 | 1.66 k | - | 20.66 k | 8.20 k | 20.68 k | 8.20 k | 720.3 |
| EXPONENTIATION | 33600 / 97 | 1.17 k | 5.19 k | - | 20.64 k | 21.44 k | 20.99 k | 20.85 k | 6.37 k |
| DIVIDE | 33600 / 781 | 13.15 k | 154.32 k | - | 1.65 M | 1.77 M | 1.69 M | 1.84 M | 36.43 M |
| ISQRT | 33600 / 195 | 930.0 | 8.30 k | - | 190.94 k | 444.81 k | 189.76 k | 437.89 k | 2.37 M |
| MODMUL_R2 | 33600 / 6250 | 801.27 k | 9.47 M | - | 29.83 M | 33.73 M | 32.93 M | 38.55 M | 966.14 k |

### Device 1 - cpu--0xd87 (CPU)

#### secp256k1 (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | - | - | - | - | - | - | - | 93.99 M |
| SUBTRACT | 50000 / 50000 | - | - | - | - | - | - | - | 131.28 M |
| ADDMOD | 50000 / 50000 | - | - | - | - | - | - | - | 39.26 M |
| SUBTRACTMOD | 50000 / 50000 | - | - | - | - | - | - | - | 50.37 M |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | - | - | - | - | - | - | - | 70.06 M |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | - | - | - | - | - | - | - | 70.64 M |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | - | - | - | - | - | - | - | 9.08 M |
| COMPARE | 50000 / 50000 | - | - | - | - | - | - | - | 142.23 M |
| REDUCE | 33600 / 6250 | - | - | - | - | - | - | - | 86.08 M |
| MODMUL | 33600 / 3125 | - | - | - | - | - | - | - | 17.29 M |
| MODEXP | 33600 / 781 | - | - | - | - | - | - | - | 121.35 k |
| EXPONENTIATION | 33600 / 781 | - | - | - | - | - | - | - | 486.83 k |
| DIVIDE | 33600 / 6250 | - | - | - | - | - | - | - | 52.62 M |
| ISQRT | 33600 / 1562 | - | - | - | - | - | - | - | 26.98 M |
| MODMUL_R2 | 50000 / 50000 | - | - | - | - | - | - | - | 17.24 M |

#### rsa256(composite) (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | - | - | - | - | - | - | - | 92.85 M |
| SUBTRACT | 50000 / 50000 | - | - | - | - | - | - | - | 131.81 M |
| ADDMOD | 50000 / 50000 | - | - | - | - | - | - | - | 46.85 M |
| SUBTRACTMOD | 50000 / 50000 | - | - | - | - | - | - | - | 50.46 M |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | - | - | - | - | - | - | - | 70.66 M |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | - | - | - | - | - | - | - | 70.20 M |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | - | - | - | - | - | - | - | 9.09 M |
| COMPARE | 50000 / 50000 | - | - | - | - | - | - | - | 143.20 M |
| REDUCE | 33600 / 6250 | - | - | - | - | - | - | - | 58.28 M |
| MODMUL | 33600 / 3125 | - | - | - | - | - | - | - | 17.27 M |
| MODEXP | 33600 / 781 | - | - | - | - | - | - | - | 127.44 k |
| EXPONENTIATION | 33600 / 781 | - | - | - | - | - | - | - | 488.08 k |
| DIVIDE | 33600 / 6250 | - | - | - | - | - | - | - | 51.84 M |
| ISQRT | 33600 / 1562 | - | - | - | - | - | - | - | 27.07 M |
| MODMUL_R2 | 50000 / 50000 | - | - | - | - | - | - | - | 17.24 M |

#### brainpoolP512r1 (512-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 33600 / 25000 | - | - | - | - | - | - | - | 86.96 M |
| SUBTRACT | 33600 / 25000 | - | - | - | - | - | - | - | 126.90 M |
| ADDMOD | 33600 / 25000 | - | - | - | - | - | - | - | 42.33 M |
| SUBTRACTMOD | 33600 / 25000 | - | - | - | - | - | - | - | 44.14 M |
| MULTIPLYOPERANDSCANNING | 33600 / 25000 | - | - | - | - | - | - | - | 29.22 M |
| MULTIPLYPRODUCTSCANNING | 33600 / 25000 | - | - | - | - | - | - | - | 29.18 M |
| MONTGOMERYMULTIPLICATION | 33600 / 25000 | - | - | - | - | - | - | - | 3.94 M |
| COMPARE | 33600 / 25000 | - | - | - | - | - | - | - | 140.49 M |
| REDUCE | 33600 / 3125 | - | - | - | - | - | - | - | 56.53 M |
| MODMUL | 33600 / 1562 | - | - | - | - | - | - | - | 8.17 M |
| MODEXP | 33600 / 390 | - | - | - | - | - | - | - | 28.01 k |
| EXPONENTIATION | 33600 / 390 | - | - | - | - | - | - | - | 157.30 k |
| DIVIDE | 33600 / 3125 | - | - | - | - | - | - | - | 47.85 M |
| ISQRT | 33600 / 781 | - | - | - | - | - | - | - | 14.67 M |
| MODMUL_R2 | 33600 / 25000 | - | - | - | - | - | - | - | 8.11 M |

#### p1024 (1024-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 33600 / 12500 | - | - | - | - | - | - | - | 73.06 M |
| SUBTRACT | 33600 / 12500 | - | - | - | - | - | - | - | 101.62 M |
| ADDMOD | 33600 / 12500 | - | - | - | - | - | - | - | 29.03 M |
| SUBTRACTMOD | 33600 / 12500 | - | - | - | - | - | - | - | 38.69 M |
| MULTIPLYOPERANDSCANNING | 33600 / 12500 | - | - | - | - | - | - | - | 8.52 M |
| MULTIPLYPRODUCTSCANNING | 33600 / 12500 | - | - | - | - | - | - | - | 8.51 M |
| MONTGOMERYMULTIPLICATION | 33600 / 12500 | - | - | - | - | - | - | - | 1.26 M |
| COMPARE | 33600 / 12500 | - | - | - | - | - | - | - | 159.31 M |
| REDUCE | 33600 / 1562 | - | - | - | - | - | - | - | 81.29 M |
| MODMUL | 33600 / 781 | - | - | - | - | - | - | - | 2.98 M |
| MODEXP | 33600 / 195 | - | - | - | - | - | - | - | 4.89 k |
| EXPONENTIATION | 33600 / 195 | - | - | - | - | - | - | - | 40.67 k |
| DIVIDE | 33600 / 1562 | - | - | - | - | - | - | - | 42.82 M |
| ISQRT | 33600 / 390 | - | - | - | - | - | - | - | 6.89 M |
| MODMUL_R2 | 33600 / 12500 | - | - | - | - | - | - | - | 2.92 M |

#### p2048 (2048-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 33600 / 6250 | - | - | - | - | - | - | - | 56.93 M |
| SUBTRACT | 33600 / 6250 | - | - | - | - | - | - | - | 74.38 M |
| ADDMOD | 33600 / 6250 | - | - | - | - | - | - | - | 24.95 M |
| SUBTRACTMOD | 33600 / 6250 | - | - | - | - | - | - | - | 31.31 M |
| MULTIPLYOPERANDSCANNING | 33600 / 6250 | - | - | - | - | - | - | - | 2.78 M |
| MULTIPLYPRODUCTSCANNING | 33600 / 6250 | - | - | - | - | - | - | - | 2.78 M |
| MONTGOMERYMULTIPLICATION | 33600 / 6250 | - | - | - | - | - | - | - | 421.62 k |
| COMPARE | 33600 / 6250 | - | - | - | - | - | - | - | 179.35 M |
| REDUCE | 33600 / 781 | - | - | - | - | - | - | - | 60.34 M |
| MODMUL | 33600 / 390 | - | - | - | - | - | - | - | 1.01 M |
| MODEXP | 33600 / 97 | - | - | - | - | - | - | - | 720.3 |
| EXPONENTIATION | 33600 / 97 | - | - | - | - | - | - | - | 6.37 k |
| DIVIDE | 33600 / 781 | - | - | - | - | - | - | - | 36.43 M |
| ISQRT | 33600 / 195 | - | - | - | - | - | - | - | 2.37 M |
| MODMUL_R2 | 33600 / 6250 | - | - | - | - | - | - | - | 966.14 k |

## 5. Head to head

Best OpenCL GPU result against best OpenCL CPU result and the CPU libraries.
Ratios above 1.00x mean the GPU is faster than that baseline.

### secp256k1 (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w32-il | 4.56 G | none | n/a | 93.99 M | n/a |
| SUBTRACT | w32-il64 | 4.66 G | none | n/a | 131.28 M | n/a |
| ADDMOD | w32-il | 4.90 G | none | n/a | 39.26 M | n/a |
| SUBTRACTMOD | w32-il | 4.94 G | none | n/a | 50.37 M | n/a |
| MULTIPLYOPERANDSCANNING | w32-il64 | 3.76 G | none | n/a | 70.06 M | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il | 2.08 G | none | n/a | 70.64 M | n/a |
| MONTGOMERYMULTIPLICATION | w32-il64 | 3.23 G | none | n/a | 9.08 M | n/a |
| COMPARE | w32-il | 5.48 G | none | n/a | 142.23 M | n/a |
| REDUCE | w32-il | 181.43 M | none | n/a | 86.08 M | n/a |
| MODMUL | w32-il | 30.05 M | none | n/a | 17.29 M | n/a |
| MODEXP | w32-o64 | 443.44 k | none | n/a | 121.35 k | n/a |
| EXPONENTIATION | w32-o64 | 1.17 M | none | n/a | 486.83 k | n/a |
| DIVIDE | w32-il64 | 76.03 M | none | n/a | 52.62 M | n/a |
| ISQRT | w32-il | 2.66 M | none | n/a | 26.98 M | n/a |
| MODMUL_R2 | w32-il64 | 2.08 G | none | n/a | 17.24 M | n/a |

### rsa256(composite) (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w32-il | 4.65 G | none | n/a | 92.85 M | n/a |
| SUBTRACT | w32-il | 4.58 G | none | n/a | 131.81 M | n/a |
| ADDMOD | w32-il | 4.79 G | none | n/a | 46.85 M | n/a |
| SUBTRACTMOD | w32-il | 4.82 G | none | n/a | 50.46 M | n/a |
| MULTIPLYOPERANDSCANNING | w32-il64 | 3.68 G | none | n/a | 70.66 M | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il | 2.12 G | none | n/a | 70.20 M | n/a |
| MONTGOMERYMULTIPLICATION | w32-il64 | 3.31 G | none | n/a | 9.09 M | n/a |
| COMPARE | w32-il | 5.31 G | none | n/a | 143.20 M | n/a |
| REDUCE | w32-il | 180.43 M | none | n/a | 58.28 M | n/a |
| MODMUL | w32-il | 30.03 M | none | n/a | 17.27 M | n/a |
| MODEXP | w32-o64 | 445.26 k | none | n/a | 127.44 k | n/a |
| EXPONENTIATION | w32-o64 | 1.16 M | none | n/a | 488.08 k | n/a |
| DIVIDE | w32-il64 | 71.06 M | none | n/a | 51.84 M | n/a |
| ISQRT | w32-il | 2.66 M | none | n/a | 27.07 M | n/a |
| MODMUL_R2 | w32-il64 | 2.08 G | none | n/a | 17.24 M | n/a |

### brainpoolP512r1 (512-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w32-il | 1.74 G | none | n/a | 86.96 M | n/a |
| SUBTRACT | w32-il | 1.67 G | none | n/a | 126.90 M | n/a |
| ADDMOD | w32-il64 | 1.56 G | none | n/a | 42.33 M | n/a |
| SUBTRACTMOD | w32-il | 1.58 G | none | n/a | 44.14 M | n/a |
| MULTIPLYOPERANDSCANNING | w32-il | 442.01 M | none | n/a | 29.22 M | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 388.29 M | none | n/a | 29.18 M | n/a |
| MONTGOMERYMULTIPLICATION | w32-il64 | 877.29 M | none | n/a | 3.94 M | n/a |
| COMPARE | w32-il64 | 2.55 G | none | n/a | 140.49 M | n/a |
| REDUCE | w32-o64 | 30.70 M | none | n/a | 56.53 M | n/a |
| MODMUL | w32-o64 | 4.18 M | none | n/a | 8.17 M | n/a |
| MODEXP | w32-o64 | 30.19 k | none | n/a | 28.01 k | n/a |
| EXPONENTIATION | w32-o64 | 12.65 k | none | n/a | 157.30 k | n/a |
| DIVIDE | w32-il64 | 11.28 M | none | n/a | 47.85 M | n/a |
| ISQRT | w32-o64 | 203.72 k | none | n/a | 14.67 M | n/a |
| MODMUL_R2 | w32-il64 | 526.79 M | none | n/a | 8.11 M | n/a |

### p1024 (1024-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w32-il | 547.86 M | none | n/a | 73.06 M | n/a |
| SUBTRACT | w32-il | 550.95 M | none | n/a | 101.62 M | n/a |
| ADDMOD | w32-il64 | 482.55 M | none | n/a | 29.03 M | n/a |
| SUBTRACTMOD | w32-il64 | 482.55 M | none | n/a | 38.69 M | n/a |
| MULTIPLYOPERANDSCANNING | w32-il | 67.66 M | none | n/a | 8.52 M | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 47.04 M | none | n/a | 8.51 M | n/a |
| MONTGOMERYMULTIPLICATION | w32-il64 | 154.43 M | none | n/a | 1.26 M | n/a |
| COMPARE | w32-il64 | 876.82 M | none | n/a | 159.31 M | n/a |
| REDUCE | w32-il64 | 5.03 M | none | n/a | 81.29 M | n/a |
| MODMUL | w32-il64 | 488.09 k | none | n/a | 2.98 M | n/a |
| MODEXP | w32-o64 | 1.75 k | none | n/a | 4.89 k | n/a |
| EXPONENTIATION | w32-il64 | 927.2 | none | n/a | 40.67 k | n/a |
| DIVIDE | w32-il64 | 1.32 M | none | n/a | 42.82 M | n/a |
| ISQRT | w32-o64 | 15.35 k | none | n/a | 6.89 M | n/a |
| MODMUL_R2 | w32-il64 | 76.69 M | none | n/a | 2.92 M | n/a |

### p2048 (2048-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w32-il | 58.65 M | none | n/a | 56.93 M | n/a |
| SUBTRACT | w32-il | 56.21 M | none | n/a | 74.38 M | n/a |
| ADDMOD | w32-il | 50.88 M | none | n/a | 24.95 M | n/a |
| SUBTRACTMOD | w32-il | 49.58 M | none | n/a | 31.31 M | n/a |
| MULTIPLYOPERANDSCANNING | w32-il | 8.63 M | none | n/a | 2.78 M | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 4.27 M | none | n/a | 2.78 M | n/a |
| MONTGOMERYMULTIPLICATION | w32-il64 | 13.13 M | none | n/a | 421.62 k | n/a |
| COMPARE | w32-il | 159.18 M | none | n/a | 179.35 M | n/a |
| REDUCE | w32-il64 | 815.76 k | none | n/a | 60.34 M | n/a |
| MODMUL | w32-il | 60.98 k | none | n/a | 1.01 M | n/a |
| MODEXP | w32-il | 59.7 | none | n/a | 720.3 | n/a |
| EXPONENTIATION | w32-o64 | 61.9 | none | n/a | 6.37 k | n/a |
| DIVIDE | w32-il64 | 42.87 k | none | n/a | 36.43 M | n/a |
| ISQRT | w32-o64 | 2.58 k | none | n/a | 2.37 M | n/a |
| MODMUL_R2 | w32-il64 | 7.17 M | none | n/a | 966.14 k | n/a |

## 6. Raw data

Also written to `NVIDIA_GB10_Report.csv` for analysis.

```csv
kind,device,device_type,kernel,modulus,bits,operation,items,seconds,ops_per_sec,mismatches
library,unknown,host-cpu,gmp-1t,secp256k1,256,ADD,50000,0.000531981,93988319.012,0
library,unknown,host-cpu,gmp-nt,secp256k1,256,ADD,50000,0.000049104,1018247042.008,0
library,unknown,host-cpu,openssl-nt,secp256k1,256,ADD,50000,0.000664652,75227336.439,0
library,NVIDIA GB10,gpu,cgbn,secp256k1,256,ADD,700000,0.000095840,7303839732.888,0
opencl-kernel,NVIDIA GB10,GPU,w8,secp256k1,256,ADD,50000,0.000075231,664619760.867,0
opencl-e2e,NVIDIA GB10,GPU,w8,secp256k1,256,ADD,50000,0.000163503,305804775.699,0
opencl-kernel,NVIDIA GB10,GPU,w16,secp256k1,256,ADD,50000,0.000041632,1200999338.606,0
opencl-e2e,NVIDIA GB10,GPU,w16,secp256k1,256,ADD,50000,0.000127775,391312860.864,0
opencl-kernel,NVIDIA GB10,GPU,w32,secp256k1,256,ADD,50000,0.000025712,1944617503.841,0
opencl-e2e,NVIDIA GB10,GPU,w32,secp256k1,256,ADD,50000,0.000116111,430622430.461,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,secp256k1,256,ADD,50000,0.000023903,2091787199.090,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,secp256k1,256,ADD,50000,0.000109535,456475092.266,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,secp256k1,256,ADD,50000,0.000024128,2072281835.371,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,secp256k1,256,ADD,50000,0.000109391,457076010.457,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,secp256k1,256,ADD,50000,0.000010959,4562464595.218,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,secp256k1,256,ADD,50000,0.000095999,520838762.462,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,secp256k1,256,ADD,50000,0.000011201,4463892134.084,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,secp256k1,256,ADD,50000,0.000096833,516352825.671,0
library,unknown,host-cpu,gmp-1t,secp256k1,256,SUBTRACT,50000,0.000380878,131275630.866,0
library,unknown,host-cpu,gmp-nt,secp256k1,256,SUBTRACT,50000,0.000758972,65878582.754,0
library,unknown,host-cpu,openssl-nt,secp256k1,256,SUBTRACT,50000,0.000688684,72602238.651,0
library,NVIDIA GB10,gpu,cgbn,secp256k1,256,SUBTRACT,700000,0.000096128,7281957390.146,0
opencl-kernel,NVIDIA GB10,GPU,w8,secp256k1,256,SUBTRACT,50000,0.000075487,662365723.277,0
opencl-e2e,NVIDIA GB10,GPU,w8,secp256k1,256,SUBTRACT,50000,0.000162975,306795540.537,0
opencl-kernel,NVIDIA GB10,GPU,w16,secp256k1,256,SUBTRACT,50000,0.000040928,1221657868.056,0
opencl-e2e,NVIDIA GB10,GPU,w16,secp256k1,256,SUBTRACT,50000,0.000127552,391996975.896,0
opencl-kernel,NVIDIA GB10,GPU,w32,secp256k1,256,SUBTRACT,50000,0.000025199,1984205830.468,0
opencl-e2e,NVIDIA GB10,GPU,w32,secp256k1,256,SUBTRACT,50000,0.000116687,428496778.805,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,secp256k1,256,SUBTRACT,50000,0.000023952,2087508231.808,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,secp256k1,256,SUBTRACT,50000,0.000108960,458884005.474,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,secp256k1,256,SUBTRACT,50000,0.000024591,2033263784.136,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,secp256k1,256,SUBTRACT,50000,0.000108896,459153696.429,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,secp256k1,256,SUBTRACT,50000,0.000010880,4595587662.973,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,secp256k1,256,SUBTRACT,50000,0.000096016,520746564.206,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,secp256k1,256,SUBTRACT,50000,0.000010720,4664177778.743,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,secp256k1,256,SUBTRACT,50000,0.000095393,524147477.674,0
library,unknown,host-cpu,gmp-1t,secp256k1,256,ADDMOD,50000,0.001273417,39264435.724,0
library,unknown,host-cpu,gmp-nt,secp256k1,256,ADDMOD,50000,0.000841547,59414388.528,0
library,unknown,host-cpu,openssl-nt,secp256k1,256,ADDMOD,50000,0.001184569,42209444.870,0
library,NVIDIA GB10,gpu,cgbn,secp256k1,256,ADDMOD,700000,0.000090880,7702464788.732,0
opencl-kernel,NVIDIA GB10,GPU,w8,secp256k1,256,ADDMOD,50000,0.000117087,427032856.775,0
opencl-e2e,NVIDIA GB10,GPU,w8,secp256k1,256,ADDMOD,50000,0.000201678,247919962.861,0
opencl-kernel,NVIDIA GB10,GPU,w16,secp256k1,256,ADDMOD,50000,0.000060928,820640810.902,0
opencl-e2e,NVIDIA GB10,GPU,w16,secp256k1,256,ADDMOD,50000,0.000152079,328776469.224,0
opencl-kernel,NVIDIA GB10,GPU,w32,secp256k1,256,ADDMOD,50000,0.000032000,1562500181.581,0
opencl-e2e,NVIDIA GB10,GPU,w32,secp256k1,256,ADDMOD,50000,0.000129263,386808277.940,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,secp256k1,256,ADDMOD,50000,0.000020480,2441406047.997,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,secp256k1,256,ADDMOD,50000,0.000106543,469294037.189,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,secp256k1,256,ADDMOD,50000,0.000016288,3069744150.422,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,secp256k1,256,ADDMOD,50000,0.000104448,478707101.582,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,secp256k1,256,ADDMOD,50000,0.000010208,4898115062.104,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,secp256k1,256,ADDMOD,50000,0.000095648,522750124.268,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,secp256k1,256,ADDMOD,50000,0.000010337,4836988088.756,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,secp256k1,256,ADDMOD,50000,0.000096353,518925260.269,0
library,unknown,host-cpu,gmp-1t,secp256k1,256,SUBTRACTMOD,50000,0.000992570,50374280.925,0
library,unknown,host-cpu,gmp-nt,secp256k1,256,SUBTRACTMOD,50000,0.000881147,56744220.667,0
library,unknown,host-cpu,openssl-nt,secp256k1,256,SUBTRACTMOD,50000,0.001784101,28025319.003,0
library,NVIDIA GB10,gpu,cgbn,secp256k1,256,SUBTRACTMOD,700000,0.000096768,7233796296.296,0
opencl-kernel,NVIDIA GB10,GPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.000116959,427500249.995,0
opencl-e2e,NVIDIA GB10,GPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.000201631,247977738.185,0
opencl-kernel,NVIDIA GB10,GPU,w16,secp256k1,256,SUBTRACTMOD,50000,0.000061312,815501151.774,0
opencl-e2e,NVIDIA GB10,GPU,w16,secp256k1,256,SUBTRACTMOD,50000,0.000152319,328258456.856,0
opencl-kernel,NVIDIA GB10,GPU,w32,secp256k1,256,SUBTRACTMOD,50000,0.000031888,1567987724.698,0
opencl-e2e,NVIDIA GB10,GPU,w32,secp256k1,256,SUBTRACTMOD,50000,0.000119151,419635607.549,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,secp256k1,256,SUBTRACTMOD,50000,0.000020192,2476227627.877,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,secp256k1,256,SUBTRACTMOD,50000,0.000107248,466209187.973,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,secp256k1,256,SUBTRACTMOD,50000,0.000016768,2981870710.100,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,secp256k1,256,SUBTRACTMOD,50000,0.000104032,480621320.642,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,secp256k1,256,SUBTRACTMOD,50000,0.000010128,4936807856.423,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,secp256k1,256,SUBTRACTMOD,50000,0.000097215,514323915.885,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,secp256k1,256,SUBTRACTMOD,50000,0.000010336,4837457973.630,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,secp256k1,256,SUBTRACTMOD,50000,0.000096242,519523691.247,0
library,unknown,host-cpu,gmp-1t,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000713660,70061373.431,0
library,unknown,host-cpu,gmp-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000619565,80701783.847,0
library,unknown,host-cpu,openssl-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000787772,63470140.670,0
opencl-kernel,NVIDIA GB10,GPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.002772320,18035436.095,0
opencl-e2e,NVIDIA GB10,GPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.002906495,17202850.816,0
opencl-kernel,NVIDIA GB10,GPU,w16,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000831292,60147336.132,0
opencl-e2e,NVIDIA GB10,GPU,w16,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000917852,54475013.540,0
opencl-kernel,NVIDIA GB10,GPU,w32,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000238767,209409176.811,0
opencl-e2e,NVIDIA GB10,GPU,w32,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000351103,142408349.802,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000048912,1022244044.119,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000159567,313347990.714,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000048768,1025262296.319,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000159951,312595733.109,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000014336,3487724695.837,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000130047,384476372.628,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000013296,3760526036.369,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000131266,380905939.661,0
library,unknown,host-cpu,gmp-1t,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000707820,70639428.137,0
library,unknown,host-cpu,gmp-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000566828,88210179.163,0
library,unknown,host-cpu,openssl-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000494109,101192249.828,0
library,NVIDIA GB10,gpu,cgbn,secp256k1,256,MULTIPLYPRODUCTSCANNING,700000,0.000094528,7405213270.142,0
opencl-kernel,NVIDIA GB10,GPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000317407,157526456.224,0
opencl-e2e,NVIDIA GB10,GPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000427901,116849461.407,0
opencl-kernel,NVIDIA GB10,GPU,w16,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000102816,486305614.099,0
opencl-e2e,NVIDIA GB10,GPU,w16,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000214575,233018769.398,0
opencl-kernel,NVIDIA GB10,GPU,w32,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000052528,951873200.659,0
opencl-e2e,NVIDIA GB10,GPU,w32,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000166991,299417337.511,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000051648,968091843.426,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000168079,297479159.790,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000047791,1046221936.786,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000159632,313220392.954,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000024048,2079174953.860,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000134944,370524071.987,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000025408,1967884654.965,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000142305,351358003.584,0
library,unknown,host-cpu,gmp-1t,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.005506592,9080026.263,0
library,unknown,host-cpu,gmp-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.001552999,32195771.089,0
library,unknown,host-cpu,openssl-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000583292,85720359.666,0
library,NVIDIA GB10,gpu,cgbn,secp256k1,256,MONTGOMERYMULTIPLICATION,700000,0.000096992,7217090069.284,0
opencl-kernel,NVIDIA GB10,GPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000174239,286962158.941,0
opencl-e2e,NVIDIA GB10,GPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000263374,189844104.912,0
opencl-kernel,NVIDIA GB10,GPU,w16,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000046512,1074991501.621,0
opencl-e2e,NVIDIA GB10,GPU,w16,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000132447,377509475.088,0
opencl-kernel,NVIDIA GB10,GPU,w32,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000016464,3036930305.454,0
opencl-e2e,NVIDIA GB10,GPU,w32,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000104431,478785013.413,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000018240,2741226442.346,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000104608,477974917.366,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000016848,2967713208.238,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000103696,482178666.631,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000017920,2790178623.788,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000103792,481732691.392,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000015488,3228306560.672,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000104673,477678154.636,0
library,unknown,host-cpu,gmp-1t,secp256k1,256,COMPARE,50000,0.000351550,142227279.466,0
library,unknown,host-cpu,gmp-nt,secp256k1,256,COMPARE,50000,0.000695692,71870884.203,0
library,unknown,host-cpu,openssl-nt,secp256k1,256,COMPARE,50000,0.000753756,66334463.553,0
library,NVIDIA GB10,gpu,cgbn,secp256k1,256,COMPARE,700000,0.000091456,7653953813.856,0
opencl-kernel,NVIDIA GB10,GPU,w8,secp256k1,256,COMPARE,50000,0.000072848,686360701.835,0
opencl-e2e,NVIDIA GB10,GPU,w8,secp256k1,256,COMPARE,50000,0.000157391,317680177.096,0
opencl-kernel,NVIDIA GB10,GPU,w16,secp256k1,256,COMPARE,50000,0.000039840,1255020058.851,0
opencl-e2e,NVIDIA GB10,GPU,w16,secp256k1,256,COMPARE,50000,0.000128095,390335270.196,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,secp256k1,256,COMPARE,50000,0.000017295,2891009813.858,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,secp256k1,256,COMPARE,50000,0.000104400,478927227.847,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,secp256k1,256,COMPARE,50000,0.000013168,3797084131.913,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,secp256k1,256,COMPARE,50000,0.000105471,474063956.185,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,secp256k1,256,COMPARE,50000,0.000009120,5482452884.692,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,secp256k1,256,COMPARE,50000,0.000094303,530205810.736,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,secp256k1,256,COMPARE,50000,0.000009184,5444248061.858,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,secp256k1,256,COMPARE,50000,0.000094737,527776903.625,0
library,unknown,host-cpu,gmp-1t,secp256k1,256,REDUCE,6250,0.000072607,86079855.377,0
library,unknown,host-cpu,gmp-nt,secp256k1,256,REDUCE,6250,0.000395998,15782908.197,0
library,unknown,host-cpu,openssl-nt,secp256k1,256,REDUCE,6250,0.000608109,10277762.772,0
library,NVIDIA GB10,gpu,cgbn,secp256k1,256,REDUCE,700000,0.000092704,7550914739.386,0
opencl-kernel,NVIDIA GB10,GPU,w8,secp256k1,256,REDUCE,33600,0.000165455,203076365.592,0
opencl-e2e,NVIDIA GB10,GPU,w8,secp256k1,256,REDUCE,33600,0.000226639,148253394.397,0
opencl-kernel,NVIDIA GB10,GPU,w16,secp256k1,256,REDUCE,33600,0.000104768,320708619.494,0
opencl-e2e,NVIDIA GB10,GPU,w16,secp256k1,256,REDUCE,33600,0.000167135,201035086.917,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,secp256k1,256,REDUCE,33600,0.000034543,972700644.723,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,secp256k1,256,REDUCE,33600,0.000104272,322234169.332,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,secp256k1,256,REDUCE,33600,0.000035008,959780797.741,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,secp256k1,256,REDUCE,33600,0.000104143,322633278.200,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,secp256k1,256,REDUCE,33600,0.000034448,975382994.493,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,secp256k1,256,REDUCE,33600,0.000103792,323724368.615,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,secp256k1,256,REDUCE,33600,0.000042065,798763759.569,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,secp256k1,256,REDUCE,33600,0.000104673,320999675.289,0
library,unknown,host-cpu,gmp-1t,secp256k1,256,MODMUL,3125,0.000180751,17288977.849,0
library,unknown,host-cpu,gmp-nt,secp256k1,256,MODMUL,3125,0.000795147,3930090.897,0
library,unknown,host-cpu,openssl-nt,secp256k1,256,MODMUL,3125,0.000682924,4575911.808,0
library,NVIDIA GB10,gpu,cgbn,secp256k1,256,MODMUL,700000,0.000221216,3164328077.535,0
opencl-kernel,NVIDIA GB10,GPU,w8,secp256k1,256,MODMUL,33600,0.000491709,68333099.832,0
opencl-e2e,NVIDIA GB10,GPU,w8,secp256k1,256,MODMUL,33600,0.000558093,60205019.484,0
opencl-kernel,NVIDIA GB10,GPU,w16,secp256k1,256,MODMUL,33600,0.000265183,126704954.687,0
opencl-e2e,NVIDIA GB10,GPU,w16,secp256k1,256,MODMUL,33600,0.000327167,102699844.629,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,secp256k1,256,MODMUL,33600,0.000105472,318567989.747,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,secp256k1,256,MODMUL,33600,0.000169391,198357660.980,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,secp256k1,256,MODMUL,33600,0.000104736,320806604.595,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,secp256k1,256,MODMUL,33600,0.000167440,200668882.102,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,secp256k1,256,MODMUL,33600,0.000103999,323080023.335,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,secp256k1,256,MODMUL,33600,0.000167615,200459387.713,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,secp256k1,256,MODMUL,33600,0.000121714,276056998.969,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,secp256k1,256,MODMUL,33600,0.000190354,176513223.451,0
library,unknown,host-cpu,gmp-1t,secp256k1,256,MODEXP,781,0.006436170,121345.458,0
library,unknown,host-cpu,gmp-nt,secp256k1,256,MODEXP,781,0.001269321,615289.591,0
library,unknown,host-cpu,openssl-nt,secp256k1,256,MODEXP,781,0.000781675,999136.469,0
library,NVIDIA GB10,gpu,cgbn,secp256k1,256,MODEXP,700000,0.063936286,10948399.474,0
opencl-kernel,NVIDIA GB10,GPU,w8,secp256k1,256,MODEXP,33600,0.024914081,1348634.935,0
opencl-e2e,NVIDIA GB10,GPU,w8,secp256k1,256,MODEXP,33600,0.024978592,1345151.880,0
opencl-kernel,NVIDIA GB10,GPU,w16,secp256k1,256,MODEXP,33600,0.003783470,8880736.455,0
opencl-e2e,NVIDIA GB10,GPU,w16,secp256k1,256,MODEXP,33600,0.003851518,8723833.064,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,secp256k1,256,MODEXP,33600,0.003061249,10975912.090,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,secp256k1,256,MODEXP,33600,0.003128081,10741409.805,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,secp256k1,256,MODEXP,33600,0.001761223,19077652.287,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,secp256k1,256,MODEXP,33600,0.001826311,18397742.813,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,secp256k1,256,MODEXP,33600,0.003018065,11132960.996,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,secp256k1,256,MODEXP,33600,0.003082097,10901668.576,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,secp256k1,256,MODEXP,33600,0.001808552,18578398.644,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,secp256k1,256,MODEXP,33600,0.001873528,17934079.422,0
library,unknown,host-cpu,gmp-1t,secp256k1,256,EXPONENTIATION,781,0.001604247,486832.761,0
library,unknown,host-cpu,gmp-nt,secp256k1,256,EXPONENTIATION,781,0.000595852,1310728.156,0
library,unknown,host-cpu,openssl-nt,secp256k1,256,EXPONENTIATION,781,0.002548129,306499.395,0
opencl-kernel,NVIDIA GB10,GPU,w8,secp256k1,256,EXPONENTIATION,33600,0.052513362,639837.152,0
opencl-e2e,NVIDIA GB10,GPU,w8,secp256k1,256,EXPONENTIATION,33600,0.052344723,641898.516,0
opencl-kernel,NVIDIA GB10,GPU,w16,secp256k1,256,EXPONENTIATION,33600,0.015773351,2130175.129,0
opencl-e2e,NVIDIA GB10,GPU,w16,secp256k1,256,EXPONENTIATION,33600,0.016043222,2094342.396,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,secp256k1,256,EXPONENTIATION,33600,0.000696844,48217391.601,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,secp256k1,256,EXPONENTIATION,33600,0.000779692,43093939.864,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,secp256k1,256,EXPONENTIATION,33600,0.000669581,50180635.301,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,secp256k1,256,EXPONENTIATION,33600,0.000733245,45823701.236,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,secp256k1,256,EXPONENTIATION,33600,0.000775676,43317054.926,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,secp256k1,256,EXPONENTIATION,33600,0.000840076,39996381.376,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,secp256k1,256,EXPONENTIATION,33600,0.000765354,43901253.870,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,secp256k1,256,EXPONENTIATION,33600,0.000828778,40541616.331,0
library,unknown,host-cpu,gmp-1t,secp256k1,256,DIVIDE,6250,0.000118767,52624046.003,0
library,unknown,host-cpu,gmp-nt,secp256k1,256,DIVIDE,6250,0.000535037,11681435.012,0
library,unknown,host-cpu,openssl-nt,secp256k1,256,DIVIDE,6250,0.000720492,8674627.934,0
library,NVIDIA GB10,gpu,cgbn,secp256k1,256,DIVIDE,700000,0.000090560,7729681978.799,0
opencl-kernel,NVIDIA GB10,GPU,w8,secp256k1,256,DIVIDE,33600,0.000371998,90323065.924,0
opencl-e2e,NVIDIA GB10,GPU,w8,secp256k1,256,DIVIDE,33600,0.000454254,73967425.984,0
opencl-kernel,NVIDIA GB10,GPU,w16,secp256k1,256,DIVIDE,33600,0.000235967,142392790.166,0
opencl-e2e,NVIDIA GB10,GPU,w16,secp256k1,256,DIVIDE,33600,0.000318862,105374737.910,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,secp256k1,256,DIVIDE,33600,0.000092352,363825370.448,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,secp256k1,256,DIVIDE,33600,0.000173712,193423609.714,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,secp256k1,256,DIVIDE,33600,0.000089728,374465062.595,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,secp256k1,256,DIVIDE,33600,0.000171743,195641163.046,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,secp256k1,256,DIVIDE,33600,0.000088639,379065629.535,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,secp256k1,256,DIVIDE,33600,0.000169359,198395122.771,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,secp256k1,256,DIVIDE,33600,0.000082209,408714407.835,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,secp256k1,256,DIVIDE,33600,0.000162226,207118461.261,0
library,unknown,host-cpu,gmp-1t,secp256k1,256,ISQRT,1562,0.000057887,26983601.581,0
library,unknown,host-cpu,gmp-nt,secp256k1,256,ISQRT,1562,0.000484621,3223137.211,0
opencl-kernel,NVIDIA GB10,GPU,w8,secp256k1,256,ISQRT,33600,0.004213223,7974892.379,0
opencl-e2e,NVIDIA GB10,GPU,w8,secp256k1,256,ISQRT,33600,0.004271656,7865801.931,0
opencl-kernel,NVIDIA GB10,GPU,w16,secp256k1,256,ISQRT,33600,0.003045298,11033402.974,0
opencl-e2e,NVIDIA GB10,GPU,w16,secp256k1,256,ISQRT,33600,0.003113906,10790306.447,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,secp256k1,256,ISQRT,33600,0.000634541,52951661.714,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,secp256k1,256,ISQRT,33600,0.000698989,48069426.310,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,secp256k1,256,ISQRT,33600,0.000705804,47605283.753,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,secp256k1,256,ISQRT,33600,0.000767500,43778501.956,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,secp256k1,256,ISQRT,33600,0.000587133,57227237.461,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,secp256k1,256,ISQRT,33600,0.000649197,51756246.724,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,secp256k1,256,ISQRT,33600,0.000894827,37549157.674,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,secp256k1,256,ISQRT,33600,0.000964076,34852023.722,0
library,unknown,host-cpu,gmp-1t,secp256k1,256,MODMUL_R2,50000,0.002900896,17236053.894,0
library,unknown,host-cpu,gmp-nt,secp256k1,256,MODMUL_R2,50000,0.000725820,68887602.441,0
library,unknown,host-cpu,openssl-nt,secp256k1,256,MODMUL_R2,50000,0.001703638,29348957.913,0
library,NVIDIA GB10,gpu,cgbn,secp256k1,256,MODMUL_R2,700000,0.000106624,6565126050.420,0
opencl-kernel,NVIDIA GB10,GPU,w8,secp256k1,256,MODMUL_R2,50000,0.000200430,249463666.226,0
opencl-e2e,NVIDIA GB10,GPU,w8,secp256k1,256,MODMUL_R2,50000,0.000287391,173979005.435,0
opencl-kernel,NVIDIA GB10,GPU,w16,secp256k1,256,MODMUL_R2,50000,0.000044768,1116869023.291,0
opencl-e2e,NVIDIA GB10,GPU,w16,secp256k1,256,MODMUL_R2,50000,0.000130160,384142597.680,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,secp256k1,256,MODMUL_R2,50000,0.000028256,1769536161.076,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,secp256k1,256,MODMUL_R2,50000,0.000118752,421045555.834,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,secp256k1,256,MODMUL_R2,50000,0.000025855,1933861999.547,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,secp256k1,256,MODMUL_R2,50000,0.000116239,430148242.882,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,secp256k1,256,MODMUL_R2,50000,0.000027695,1805379945.418,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,secp256k1,256,MODMUL_R2,50000,0.000125568,398190633.572,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,secp256k1,256,MODMUL_R2,50000,0.000024001,2083246734.928,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,secp256k1,256,MODMUL_R2,50000,0.000116593,428842209.178,0
library,unknown,host-cpu,gmp-1t,rsa256(composite),256,ADD,50000,0.000538525,92846198.285,0
library,unknown,host-cpu,gmp-nt,rsa256(composite),256,ADD,50000,0.000728748,68610823.264,0
library,unknown,host-cpu,openssl-nt,rsa256(composite),256,ADD,50000,0.000665628,75117032.108,0
library,NVIDIA GB10,gpu,cgbn,rsa256(composite),256,ADD,700000,0.000085824,8156226696.495,0
opencl-kernel,NVIDIA GB10,GPU,w8,rsa256(composite),256,ADD,50000,0.000075296,664045862.959,0
opencl-e2e,NVIDIA GB10,GPU,w8,rsa256(composite),256,ADD,50000,0.000163519,305774867.250,0
opencl-kernel,NVIDIA GB10,GPU,w16,rsa256(composite),256,ADD,50000,0.000041136,1215480097.706,0
opencl-e2e,NVIDIA GB10,GPU,w16,rsa256(composite),256,ADD,50000,0.000127583,391901742.413,0
opencl-kernel,NVIDIA GB10,GPU,w32,rsa256(composite),256,ADD,50000,0.000025360,1971608753.642,0
opencl-e2e,NVIDIA GB10,GPU,w32,rsa256(composite),256,ADD,50000,0.000116367,429675044.744,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,rsa256(composite),256,ADD,50000,0.000023936,2088902975.190,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,rsa256(composite),256,ADD,50000,0.000116896,427730629.829,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,rsa256(composite),256,ADD,50000,0.000024191,2066884206.562,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,rsa256(composite),256,ADD,50000,0.000117247,426450169.588,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,rsa256(composite),256,ADD,50000,0.000010752,4650296447.558,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,rsa256(composite),256,ADD,50000,0.000095695,522493284.678,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,rsa256(composite),256,ADD,50000,0.000011472,4358437035.327,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,rsa256(composite),256,ADD,50000,0.000096097,520307570.866,0
library,unknown,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACT,50000,0.000379326,131812742.722,0
library,unknown,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACT,50000,0.000526061,95046010.857,0
library,unknown,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACT,50000,0.001165481,42900742.518,0
library,NVIDIA GB10,gpu,cgbn,rsa256(composite),256,SUBTRACT,700000,0.000096096,7284382284.382,0
opencl-kernel,NVIDIA GB10,GPU,w8,rsa256(composite),256,SUBTRACT,50000,0.000076336,654998949.020,0
opencl-e2e,NVIDIA GB10,GPU,w8,rsa256(composite),256,SUBTRACT,50000,0.000163727,305386426.252,0
opencl-kernel,NVIDIA GB10,GPU,w16,rsa256(composite),256,SUBTRACT,50000,0.000041264,1211710173.647,0
opencl-e2e,NVIDIA GB10,GPU,w16,rsa256(composite),256,SUBTRACT,50000,0.000128175,390091668.471,0
opencl-kernel,NVIDIA GB10,GPU,w32,rsa256(composite),256,SUBTRACT,50000,0.000025728,1943407627.650,0
opencl-e2e,NVIDIA GB10,GPU,w32,rsa256(composite),256,SUBTRACT,50000,0.000117199,426624798.115,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,rsa256(composite),256,SUBTRACT,50000,0.000024208,2065433029.988,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,rsa256(composite),256,SUBTRACT,50000,0.000109215,457812611.288,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,rsa256(composite),256,SUBTRACT,50000,0.000024144,2070907943.825,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,rsa256(composite),256,SUBTRACT,50000,0.000109408,457004942.474,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,rsa256(composite),256,SUBTRACT,50000,0.000010928,4575399234.321,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,rsa256(composite),256,SUBTRACT,50000,0.000096224,519620879.058,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,rsa256(composite),256,SUBTRACT,50000,0.000011616,4304410095.647,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,rsa256(composite),256,SUBTRACT,50000,0.000097025,515331118.033,0
library,unknown,host-cpu,gmp-1t,rsa256(composite),256,ADDMOD,50000,0.001067210,46851135.054,0
library,unknown,host-cpu,gmp-nt,rsa256(composite),256,ADDMOD,50000,0.000711068,70316762.786,0
library,unknown,host-cpu,openssl-nt,rsa256(composite),256,ADDMOD,50000,0.000835019,59878877.785,0
library,NVIDIA GB10,gpu,cgbn,rsa256(composite),256,ADDMOD,700000,0.000090208,7759843916.282,0
opencl-kernel,NVIDIA GB10,GPU,w8,rsa256(composite),256,ADDMOD,50000,0.000106080,471342389.066,0
opencl-e2e,NVIDIA GB10,GPU,w8,rsa256(composite),256,ADDMOD,50000,0.000188735,264921704.231,0
opencl-kernel,NVIDIA GB10,GPU,w16,rsa256(composite),256,ADDMOD,50000,0.000055072,907902387.807,0
opencl-e2e,NVIDIA GB10,GPU,w16,rsa256(composite),256,ADDMOD,50000,0.000140271,356452855.929,0
opencl-kernel,NVIDIA GB10,GPU,w32,rsa256(composite),256,ADDMOD,50000,0.000028416,1759571757.532,0
opencl-e2e,NVIDIA GB10,GPU,w32,rsa256(composite),256,ADDMOD,50000,0.000117151,426799622.734,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,rsa256(composite),256,ADDMOD,50000,0.000019904,2512056895.159,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,rsa256(composite),256,ADDMOD,50000,0.000106592,469078384.578,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,rsa256(composite),256,ADDMOD,50000,0.000016383,3051944596.303,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,rsa256(composite),256,ADDMOD,50000,0.000104160,480030701.634,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,rsa256(composite),256,ADDMOD,50000,0.000010432,4792941985.989,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,rsa256(composite),256,ADDMOD,50000,0.000095535,523368415.755,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,rsa256(composite),256,ADDMOD,50000,0.000010576,4727686534.934,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,rsa256(composite),256,ADDMOD,50000,0.000095729,522307747.832,0
library,unknown,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACTMOD,50000,0.000990971,50455562.970,0
library,unknown,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACTMOD,50000,0.000721548,69295459.207,0
library,unknown,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACTMOD,50000,0.000980171,51011506.594,0
library,NVIDIA GB10,gpu,cgbn,rsa256(composite),256,SUBTRACTMOD,700000,0.000092896,7535308301.757,0
opencl-kernel,NVIDIA GB10,GPU,w8,rsa256(composite),256,SUBTRACTMOD,50000,0.000116895,427734303.868,0
opencl-e2e,NVIDIA GB10,GPU,w8,rsa256(composite),256,SUBTRACTMOD,50000,0.000202911,246413447.679,0
opencl-kernel,NVIDIA GB10,GPU,w16,rsa256(composite),256,SUBTRACTMOD,50000,0.000061152,817634637.513,0
opencl-e2e,NVIDIA GB10,GPU,w16,rsa256(composite),256,SUBTRACTMOD,50000,0.000152255,328396437.999,0
opencl-kernel,NVIDIA GB10,GPU,w32,rsa256(composite),256,SUBTRACTMOD,50000,0.000032416,1542447867.603,0
opencl-e2e,NVIDIA GB10,GPU,w32,rsa256(composite),256,SUBTRACTMOD,50000,0.000118864,420648802.684,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,rsa256(composite),256,SUBTRACTMOD,50000,0.000020128,2484102876.188,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,rsa256(composite),256,SUBTRACTMOD,50000,0.000105519,473848342.111,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,rsa256(composite),256,SUBTRACTMOD,50000,0.000016912,2956480072.691,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,rsa256(composite),256,SUBTRACTMOD,50000,0.000104351,479152033.448,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,rsa256(composite),256,SUBTRACTMOD,50000,0.000010384,4815100467.360,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,rsa256(composite),256,SUBTRACTMOD,50000,0.000095695,522493284.678,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,rsa256(composite),256,SUBTRACTMOD,50000,0.000010624,4706329948.019,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,rsa256(composite),256,SUBTRACTMOD,50000,0.000095537,523357414.680,0
library,unknown,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000707659,70655499.602,0
library,unknown,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000698332,71599182.988,0
library,unknown,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000865323,57781891.623,0
opencl-kernel,NVIDIA GB10,GPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.002779968,17985818.570,0
opencl-e2e,NVIDIA GB10,GPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.002872816,17404525.670,0
opencl-kernel,NVIDIA GB10,GPU,w16,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000806845,61969771.377,0
opencl-e2e,NVIDIA GB10,GPU,w16,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000914971,54646541.013,0
opencl-kernel,NVIDIA GB10,GPU,w32,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000238911,209282966.296,0
opencl-e2e,NVIDIA GB10,GPU,w32,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000351710,142162578.081,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000048112,1039241915.281,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000160511,311505128.523,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000048544,1029993494.085,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000167359,298758952.885,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000014560,3434064496.940,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000130240,383906661.792,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000013569,3684868263.670,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000131074,381463892.913,0
library,unknown,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000712284,70196720.387,0
library,unknown,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000675420,74028012.163,0
library,unknown,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000708828,70538974.564,0
library,NVIDIA GB10,gpu,cgbn,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,700000,0.000093920,7453151618.399,0
opencl-kernel,NVIDIA GB10,GPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000317294,157582554.983,0
opencl-e2e,NVIDIA GB10,GPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000427549,116945658.274,0
opencl-kernel,NVIDIA GB10,GPU,w16,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000102959,485630175.220,0
opencl-e2e,NVIDIA GB10,GPU,w16,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000215375,232153225.863,0
opencl-kernel,NVIDIA GB10,GPU,w32,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000052208,957707481.910,0
opencl-e2e,NVIDIA GB10,GPU,w32,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000166831,299704491.982,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000052704,948694675.314,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000167983,297649163.284,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000047296,1057171974.761,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000159039,314388290.685,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000023632,2115776098.445,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000136079,367433649.629,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000025937,1927748078.437,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000142962,349743282.486,0
library,unknown,host-cpu,gmp-1t,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.005503216,9085596.492,0
library,unknown,host-cpu,gmp-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.001291016,38729186.925,0
library,unknown,host-cpu,openssl-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000833691,59974259.341,0
library,NVIDIA GB10,gpu,cgbn,rsa256(composite),256,MONTGOMERYMULTIPLICATION,700000,0.000089792,7795794725.588,0
opencl-kernel,NVIDIA GB10,GPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000172623,289648521.723,0
opencl-e2e,NVIDIA GB10,GPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000261999,190840425.020,0
opencl-kernel,NVIDIA GB10,GPU,w16,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000045263,1104655256.138,0
opencl-e2e,NVIDIA GB10,GPU,w16,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000129887,384949980.971,0
opencl-kernel,NVIDIA GB10,GPU,w32,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000016623,3007880283.352,0
opencl-e2e,NVIDIA GB10,GPU,w32,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000104048,480547447.401,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000017776,2812779939.061,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000105151,475506632.482,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000016560,3019322437.093,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000104032,480621320.642,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000017936,2787688510.595,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000104000,480769202.552,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000015120,3306879860.794,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000104641,477824164.208,0
library,unknown,host-cpu,gmp-1t,rsa256(composite),256,COMPARE,50000,0.000349166,143198362.992,0
library,unknown,host-cpu,gmp-nt,rsa256(composite),256,COMPARE,50000,0.000679116,73625124.684,0
library,unknown,host-cpu,openssl-nt,rsa256(composite),256,COMPARE,50000,0.000793820,62986571.289,0
library,NVIDIA GB10,gpu,cgbn,rsa256(composite),256,COMPARE,700000,0.000090432,7740622788.393,0
opencl-kernel,NVIDIA GB10,GPU,w8,rsa256(composite),256,COMPARE,50000,0.000071983,694608462.401,0
opencl-e2e,NVIDIA GB10,GPU,w8,rsa256(composite),256,COMPARE,50000,0.000163727,305386399.110,0
opencl-kernel,NVIDIA GB10,GPU,w16,rsa256(composite),256,COMPARE,50000,0.000038640,1293995655.064,0
opencl-e2e,NVIDIA GB10,GPU,w16,rsa256(composite),256,COMPARE,50000,0.000127744,391407808.326,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,rsa256(composite),256,COMPARE,50000,0.000016144,3097126516.058,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,rsa256(composite),256,COMPARE,50000,0.000104735,477395291.007,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,rsa256(composite),256,COMPARE,50000,0.000013872,3604381346.450,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,rsa256(composite),256,COMPARE,50000,0.000103855,481440419.725,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,rsa256(composite),256,COMPARE,50000,0.000009424,5305605573.708,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,rsa256(composite),256,COMPARE,50000,0.000095680,522575213.807,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,rsa256(composite),256,COMPARE,50000,0.000009424,5305605573.708,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,rsa256(composite),256,COMPARE,50000,0.000094817,527331566.360,0
library,unknown,host-cpu,gmp-1t,rsa256(composite),256,REDUCE,6250,0.000107247,58276686.190,0
library,unknown,host-cpu,gmp-nt,rsa256(composite),256,REDUCE,6250,0.000751308,8318825.403,0
library,unknown,host-cpu,openssl-nt,rsa256(composite),256,REDUCE,6250,0.000503277,12418608.474,0
library,NVIDIA GB10,gpu,cgbn,rsa256(composite),256,REDUCE,700000,0.000084288,8304859529.233,0
opencl-kernel,NVIDIA GB10,GPU,w8,rsa256(composite),256,REDUCE,33600,0.000168175,199791883.456,0
opencl-e2e,NVIDIA GB10,GPU,w8,rsa256(composite),256,REDUCE,33600,0.000229614,146332540.061,0
opencl-kernel,NVIDIA GB10,GPU,w16,rsa256(composite),256,REDUCE,33600,0.000104992,320024396.182,0
opencl-e2e,NVIDIA GB10,GPU,w16,rsa256(composite),256,REDUCE,33600,0.000168175,199791883.456,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,rsa256(composite),256,REDUCE,33600,0.000034688,968634771.282,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,rsa256(composite),256,REDUCE,33600,0.000097279,345398286.233,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,rsa256(composite),256,REDUCE,33600,0.000035263,952840051.224,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,rsa256(composite),256,REDUCE,33600,0.000103728,323924103.886,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,rsa256(composite),256,REDUCE,33600,0.000034640,969976772.529,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,rsa256(composite),256,REDUCE,33600,0.000104144,322630167.600,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,rsa256(composite),256,REDUCE,33600,0.000042481,790941857.535,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,rsa256(composite),256,REDUCE,33600,0.000105105,319680346.340,0
library,unknown,host-cpu,gmp-1t,rsa256(composite),256,MODMUL,3125,0.000180911,17273687.417,0
library,unknown,host-cpu,gmp-nt,rsa256(composite),256,MODMUL,3125,0.000720188,4339144.736,0
library,unknown,host-cpu,openssl-nt,rsa256(composite),256,MODMUL,3125,0.000692604,4511957.720,0
library,NVIDIA GB10,gpu,cgbn,rsa256(composite),256,MODMUL,700000,0.000218752,3199970743.125,0
opencl-kernel,NVIDIA GB10,GPU,w8,rsa256(composite),256,MODMUL,33600,0.000491037,68426616.925,0
opencl-e2e,NVIDIA GB10,GPU,w8,rsa256(composite),256,MODMUL,33600,0.000557085,60313955.203,0
opencl-kernel,NVIDIA GB10,GPU,w16,rsa256(composite),256,MODMUL,33600,0.000265166,126713083.187,0
opencl-e2e,NVIDIA GB10,GPU,w16,rsa256(composite),256,MODMUL,33600,0.000327822,102494647.861,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,rsa256(composite),256,MODMUL,33600,0.000105199,319394681.353,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,rsa256(composite),256,MODMUL,33600,0.000168064,199923840.929,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,rsa256(composite),256,MODMUL,33600,0.000104591,321251341.170,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,rsa256(composite),256,MODMUL,33600,0.000167168,200995396.844,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,rsa256(composite),256,MODMUL,33600,0.000104064,322878211.944,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,rsa256(composite),256,MODMUL,33600,0.000167455,200650920.739,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,rsa256(composite),256,MODMUL,33600,0.000121602,276311271.666,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,rsa256(composite),256,MODMUL,33600,0.000190050,176795587.606,0
library,unknown,host-cpu,gmp-1t,rsa256(composite),256,MODEXP,781,0.006128157,127444.515,0
library,unknown,host-cpu,gmp-nt,rsa256(composite),256,MODEXP,781,0.000973098,802591.304,0
library,unknown,host-cpu,openssl-nt,rsa256(composite),256,MODEXP,781,0.000897211,870475.290,0
library,NVIDIA GB10,gpu,cgbn,rsa256(composite),256,MODEXP,700000,0.062912092,11126636.832,0
opencl-kernel,NVIDIA GB10,GPU,w8,rsa256(composite),256,MODEXP,33600,0.024913168,1348684.359,0
opencl-e2e,NVIDIA GB10,GPU,w8,rsa256(composite),256,MODEXP,33600,0.024974929,1345349.170,0
opencl-kernel,NVIDIA GB10,GPU,w16,rsa256(composite),256,MODEXP,33600,0.003777871,8893898.187,0
opencl-e2e,NVIDIA GB10,GPU,w16,rsa256(composite),256,MODEXP,33600,0.003841438,8746724.543,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,rsa256(composite),256,MODEXP,33600,0.003055345,10997121.420,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,rsa256(composite),256,MODEXP,33600,0.003118256,10775253.842,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,rsa256(composite),256,MODEXP,33600,0.001754039,19155788.481,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,rsa256(composite),256,MODEXP,33600,0.001824439,18416620.184,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,rsa256(composite),256,MODEXP,33600,0.003014257,11147025.661,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,rsa256(composite),256,MODEXP,33600,0.003077553,10917764.850,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,rsa256(composite),256,MODEXP,33600,0.001800631,18660125.273,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,rsa256(composite),256,MODEXP,33600,0.001863912,18026602.160,0
library,unknown,host-cpu,gmp-1t,rsa256(composite),256,EXPONENTIATION,781,0.001600151,488078.935,0
library,unknown,host-cpu,gmp-nt,rsa256(composite),256,EXPONENTIATION,781,0.000758524,1029631.239,0
library,unknown,host-cpu,openssl-nt,rsa256(composite),256,EXPONENTIATION,781,0.003111822,250978.365,0
opencl-kernel,NVIDIA GB10,GPU,w8,rsa256(composite),256,EXPONENTIATION,33600,0.051967381,646559.425,0
opencl-e2e,NVIDIA GB10,GPU,w8,rsa256(composite),256,EXPONENTIATION,33600,0.052784384,636551.901,0
opencl-kernel,NVIDIA GB10,GPU,w16,rsa256(composite),256,EXPONENTIATION,33600,0.015311753,2194392.765,0
opencl-e2e,NVIDIA GB10,GPU,w16,rsa256(composite),256,EXPONENTIATION,33600,0.015849303,2119967.042,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,rsa256(composite),256,EXPONENTIATION,33600,0.000705868,47600967.498,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,rsa256(composite),256,EXPONENTIATION,33600,0.000768060,43746581.720,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,rsa256(composite),256,EXPONENTIATION,33600,0.000675260,49758612.421,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,rsa256(composite),256,EXPONENTIATION,33600,0.000743148,45213066.308,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,rsa256(composite),256,EXPONENTIATION,33600,0.000763980,43980208.769,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,rsa256(composite),256,EXPONENTIATION,33600,0.000841323,39937099.294,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,rsa256(composite),256,EXPONENTIATION,33600,0.000763306,44019043.401,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,rsa256(composite),256,EXPONENTIATION,33600,0.000816523,41150095.689,0
library,unknown,host-cpu,gmp-1t,rsa256(composite),256,DIVIDE,6250,0.000120560,51841407.490,0
library,unknown,host-cpu,gmp-nt,rsa256(composite),256,DIVIDE,6250,0.000547565,11414170.092,0
library,unknown,host-cpu,openssl-nt,rsa256(composite),256,DIVIDE,6250,0.000662284,9437039.081,0
library,NVIDIA GB10,gpu,cgbn,rsa256(composite),256,DIVIDE,700000,0.000089856,7790242165.242,0
opencl-kernel,NVIDIA GB10,GPU,w8,rsa256(composite),256,DIVIDE,33600,0.000366862,91587571.256,0
opencl-e2e,NVIDIA GB10,GPU,w8,rsa256(composite),256,DIVIDE,33600,0.000451278,74455213.913,0
opencl-kernel,NVIDIA GB10,GPU,w16,rsa256(composite),256,DIVIDE,33600,0.000243951,137732578.013,0
opencl-e2e,NVIDIA GB10,GPU,w16,rsa256(composite),256,DIVIDE,33600,0.000330319,101719851.892,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,rsa256(composite),256,DIVIDE,33600,0.000098192,342186779.785,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,rsa256(composite),256,DIVIDE,33600,0.000183007,183599538.518,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,rsa256(composite),256,DIVIDE,33600,0.000096016,349941638.110,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,rsa256(composite),256,DIVIDE,33600,0.000176687,190166795.038,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,rsa256(composite),256,DIVIDE,33600,0.000096768,347222239.181,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,rsa256(composite),256,DIVIDE,33600,0.000182063,184551502.775,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,rsa256(composite),256,DIVIDE,33600,0.000087953,382022225.093,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,rsa256(composite),256,DIVIDE,33600,0.000171346,196094459.835,0
library,unknown,host-cpu,gmp-1t,rsa256(composite),256,ISQRT,1562,0.000057696,27072933.519,0
library,unknown,host-cpu,gmp-nt,rsa256(composite),256,ISQRT,1562,0.000533485,2927917.377,0
opencl-kernel,NVIDIA GB10,GPU,w8,rsa256(composite),256,ISQRT,33600,0.004211912,7977374.665,0
opencl-e2e,NVIDIA GB10,GPU,w8,rsa256(composite),256,ISQRT,33600,0.004281079,7848488.684,0
opencl-kernel,NVIDIA GB10,GPU,w16,rsa256(composite),256,ISQRT,33600,0.003047186,11026566.777,0
opencl-e2e,NVIDIA GB10,GPU,w16,rsa256(composite),256,ISQRT,33600,0.003138417,10706034.285,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,rsa256(composite),256,ISQRT,33600,0.000637181,52732269.054,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,rsa256(composite),256,ISQRT,33600,0.000706124,47583710.305,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,rsa256(composite),256,ISQRT,33600,0.000709485,47358295.055,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,rsa256(composite),256,ISQRT,33600,0.000775964,43300977.918,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,rsa256(composite),256,ISQRT,33600,0.000586685,57270936.432,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,rsa256(composite),256,ISQRT,33600,0.000648092,51844490.186,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,rsa256(composite),256,ISQRT,33600,0.000895052,37539718.419,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,rsa256(composite),256,ISQRT,33600,0.000963276,34880968.522,0
library,unknown,host-cpu,gmp-1t,rsa256(composite),256,MODMUL_R2,50000,0.002900623,17237676.162,0
library,unknown,host-cpu,gmp-nt,rsa256(composite),256,MODMUL_R2,50000,0.001022938,48878817.745,0
library,unknown,host-cpu,openssl-nt,rsa256(composite),256,MODMUL_R2,50000,0.001851365,27007100.177,0
library,NVIDIA GB10,gpu,cgbn,rsa256(composite),256,MODMUL_R2,700000,0.000101280,6911532385.466,0
opencl-kernel,NVIDIA GB10,GPU,w8,rsa256(composite),256,MODMUL_R2,50000,0.000199503,250622793.398,0
opencl-e2e,NVIDIA GB10,GPU,w8,rsa256(composite),256,MODMUL_R2,50000,0.000287102,174154135.227,0
opencl-kernel,NVIDIA GB10,GPU,w16,rsa256(composite),256,MODMUL_R2,50000,0.000042896,1165609718.192,0
opencl-e2e,NVIDIA GB10,GPU,w16,rsa256(composite),256,MODMUL_R2,50000,0.000129711,385472321.261,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,rsa256(composite),256,MODMUL_R2,50000,0.000029023,1722771917.401,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,rsa256(composite),256,MODMUL_R2,50000,0.000118175,423101341.276,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,rsa256(composite),256,MODMUL_R2,50000,0.000024688,2025274770.932,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,rsa256(composite),256,MODMUL_R2,50000,0.000115920,431331993.315,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,rsa256(composite),256,MODMUL_R2,50000,0.000028239,1770601218.921,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,rsa256(composite),256,MODMUL_R2,50000,0.000116704,428434319.596,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,rsa256(composite),256,MODMUL_R2,50000,0.000024096,2075033841.427,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,rsa256(composite),256,MODMUL_R2,50000,0.000112498,444452328.218,0
library,unknown,host-cpu,gmp-1t,brainpoolP512r1,512,ADD,25000,0.000287502,86955916.831,0
library,unknown,host-cpu,gmp-nt,brainpoolP512r1,512,ADD,25000,0.000698476,35792210.643,0
library,unknown,host-cpu,openssl-nt,brainpoolP512r1,512,ADD,25000,0.000877707,28483309.241,0
library,NVIDIA GB10,gpu,cgbn,brainpoolP512r1,512,ADD,700000,0.000146560,4776200873.362,0
opencl-kernel,NVIDIA GB10,GPU,w8,brainpoolP512r1,512,ADD,33600,0.000116288,288937791.884,0
opencl-e2e,NVIDIA GB10,GPU,w8,brainpoolP512r1,512,ADD,33600,0.000235439,142712123.979,0
opencl-kernel,NVIDIA GB10,GPU,w16,brainpoolP512r1,512,ADD,33600,0.000058976,569723257.582,0
opencl-e2e,NVIDIA GB10,GPU,w16,brainpoolP512r1,512,ADD,33600,0.000179471,187216879.441,0
opencl-kernel,NVIDIA GB10,GPU,w32,brainpoolP512r1,512,ADD,33600,0.000038304,877192928.102,0
opencl-e2e,NVIDIA GB10,GPU,w32,brainpoolP512r1,512,ADD,33600,0.000156383,214857112.233,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,brainpoolP512r1,512,ADD,33600,0.000039023,861030719.384,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,brainpoolP512r1,512,ADD,33600,0.000158784,211608230.133,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,brainpoolP512r1,512,ADD,33600,0.000039999,840021107.476,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,brainpoolP512r1,512,ADD,33600,0.000156352,214899705.999,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,brainpoolP512r1,512,ADD,33600,0.000014384,2335928699.594,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,brainpoolP512r1,512,ADD,33600,0.000131904,254730709.294,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,brainpoolP512r1,512,ADD,33600,0.000014992,2241193988.534,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,brainpoolP512r1,512,ADD,33600,0.000131506,255501636.195,0
library,unknown,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACT,25000,0.000197007,126899043.759,0
library,unknown,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACT,25000,0.000681756,36670010.644,0
library,unknown,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACT,25000,0.001021370,24476928.032,0
library,NVIDIA GB10,gpu,cgbn,brainpoolP512r1,512,SUBTRACT,700000,0.000148640,4709364908.504,0
opencl-kernel,NVIDIA GB10,GPU,w8,brainpoolP512r1,512,SUBTRACT,33600,0.000117391,286222932.119,0
opencl-e2e,NVIDIA GB10,GPU,w8,brainpoolP512r1,512,SUBTRACT,33600,0.000241199,139304060.715,0
opencl-kernel,NVIDIA GB10,GPU,w16,brainpoolP512r1,512,SUBTRACT,33600,0.000060479,555564772.696,0
opencl-e2e,NVIDIA GB10,GPU,w16,brainpoolP512r1,512,SUBTRACT,33600,0.000180144,186517458.989,0
opencl-kernel,NVIDIA GB10,GPU,w32,brainpoolP512r1,512,SUBTRACT,33600,0.000038144,880872389.334,0
opencl-e2e,NVIDIA GB10,GPU,w32,brainpoolP512r1,512,SUBTRACT,33600,0.000156223,215077161.610,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,brainpoolP512r1,512,SUBTRACT,33600,0.000037840,887949259.937,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,brainpoolP512r1,512,SUBTRACT,33600,0.000157407,213459389.325,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,brainpoolP512r1,512,SUBTRACT,33600,0.000039087,859620905.709,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,brainpoolP512r1,512,SUBTRACT,33600,0.000155856,215583610.938,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,brainpoolP512r1,512,SUBTRACT,33600,0.000014992,2241196163.940,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,brainpoolP512r1,512,SUBTRACT,33600,0.000132432,253715107.010,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,brainpoolP512r1,512,SUBTRACT,33600,0.000014992,2241196163.940,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,brainpoolP512r1,512,SUBTRACT,33600,0.000133154,252339374.573,0
library,unknown,host-cpu,gmp-1t,brainpoolP512r1,512,ADDMOD,25000,0.000590652,42326106.767,0
library,unknown,host-cpu,gmp-nt,brainpoolP512r1,512,ADDMOD,25000,0.000799867,31255196.573,0
library,unknown,host-cpu,openssl-nt,brainpoolP512r1,512,ADDMOD,25000,0.000772107,32378931.950,0
library,NVIDIA GB10,gpu,cgbn,brainpoolP512r1,512,ADDMOD,700000,0.000147872,4733823847.652,0
opencl-kernel,NVIDIA GB10,GPU,w8,brainpoolP512r1,512,ADDMOD,33600,0.000157151,213807112.395,0
opencl-e2e,NVIDIA GB10,GPU,w8,brainpoolP512r1,512,ADDMOD,33600,0.000278350,120711330.545,0
opencl-kernel,NVIDIA GB10,GPU,w16,brainpoolP512r1,512,ADDMOD,33600,0.000081424,412654721.469,0
opencl-e2e,NVIDIA GB10,GPU,w16,brainpoolP512r1,512,ADDMOD,33600,0.000203967,164732536.273,0
opencl-kernel,NVIDIA GB10,GPU,w32,brainpoolP512r1,512,ADDMOD,33600,0.000045120,744680800.707,0
opencl-e2e,NVIDIA GB10,GPU,w32,brainpoolP512r1,512,ADDMOD,33600,0.000168127,199848931.777,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,brainpoolP512r1,512,ADDMOD,33600,0.000039136,858544703.505,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,brainpoolP512r1,512,ADDMOD,33600,0.000157455,213394307.265,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,brainpoolP512r1,512,ADDMOD,33600,0.000037888,886824165.796,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,brainpoolP512r1,512,ADDMOD,33600,0.000155935,215474388.437,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,brainpoolP512r1,512,ADDMOD,33600,0.000016144,2081269018.791,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,brainpoolP512r1,512,ADDMOD,33600,0.000134303,250180558.721,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,brainpoolP512r1,512,ADDMOD,33600,0.000016049,2093588161.829,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,brainpoolP512r1,512,ADDMOD,33600,0.000133809,251104171.141,0
library,unknown,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000566412,44137483.254,0
library,unknown,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000726763,34399109.525,0
library,unknown,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000995802,25105392.574,0
library,NVIDIA GB10,gpu,cgbn,brainpoolP512r1,512,SUBTRACTMOD,700000,0.000142688,4905808477.237,0
opencl-kernel,NVIDIA GB10,GPU,w8,brainpoolP512r1,512,SUBTRACTMOD,33600,0.000181455,185169886.498,0
opencl-e2e,NVIDIA GB10,GPU,w8,brainpoolP512r1,512,SUBTRACTMOD,33600,0.000302447,111093842.298,0
opencl-kernel,NVIDIA GB10,GPU,w16,brainpoolP512r1,512,SUBTRACTMOD,33600,0.000092864,361819405.821,0
opencl-e2e,NVIDIA GB10,GPU,w16,brainpoolP512r1,512,SUBTRACTMOD,33600,0.000216735,155028027.439,0
opencl-kernel,NVIDIA GB10,GPU,w32,brainpoolP512r1,512,SUBTRACTMOD,33600,0.000053360,629685205.269,0
opencl-e2e,NVIDIA GB10,GPU,w32,brainpoolP512r1,512,SUBTRACTMOD,33600,0.000170671,196869998.226,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,brainpoolP512r1,512,SUBTRACTMOD,33600,0.000038559,871391663.882,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,brainpoolP512r1,512,SUBTRACTMOD,33600,0.000156879,214177801.670,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,brainpoolP512r1,512,SUBTRACTMOD,33600,0.000038016,883838596.166,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,brainpoolP512r1,512,SUBTRACTMOD,33600,0.000156015,215363910.590,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,brainpoolP512r1,512,SUBTRACTMOD,33600,0.000015840,2121213117.979,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,brainpoolP512r1,512,SUBTRACTMOD,33600,0.000133071,252496773.312,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,brainpoolP512r1,512,SUBTRACTMOD,33600,0.000016192,2075100245.195,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,brainpoolP512r1,512,SUBTRACTMOD,33600,0.000133938,250862346.788,0
library,unknown,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.000855643,29217793.050,0
library,unknown,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.000735339,33997924.298,0
library,unknown,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.000704316,35495431.011,0
opencl-kernel,NVIDIA GB10,GPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,33600,0.008982077,3740782.891,0
opencl-e2e,NVIDIA GB10,GPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,33600,0.009174348,3662385.599,0
opencl-kernel,NVIDIA GB10,GPU,w16,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,33600,0.002620516,12821902.201,0
opencl-e2e,NVIDIA GB10,GPU,w16,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,33600,0.002762947,12160928.174,0
opencl-kernel,NVIDIA GB10,GPU,w32,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,33600,0.000837900,40100250.318,0
opencl-e2e,NVIDIA GB10,GPU,w32,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,33600,0.000990491,33922569.333,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,33600,0.000135439,248082141.383,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,33600,0.000289487,116067389.117,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,33600,0.000136159,246770299.762,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,33600,0.000290447,115683761.399,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,33600,0.000056560,594059498.285,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,33600,0.000209311,160526684.804,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,33600,0.000057057,588884827.065,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,33600,0.000215314,156051166.854,0
library,unknown,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000856619,29184503.481,0
library,unknown,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000870155,28730513.183,0
library,unknown,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000904827,27629591.061,0
library,NVIDIA GB10,gpu,cgbn,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,700000,0.000153184,4569667850.428,0
opencl-kernel,NVIDIA GB10,GPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,33600,0.001273800,26377767.451,0
opencl-e2e,NVIDIA GB10,GPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,33600,0.001431416,23473259.928,0
opencl-kernel,NVIDIA GB10,GPU,w16,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,33600,0.000354446,94795820.855,0
opencl-e2e,NVIDIA GB10,GPU,w16,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,33600,0.000509454,65952960.089,0
opencl-kernel,NVIDIA GB10,GPU,w32,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,33600,0.000107648,312128386.131,0
opencl-e2e,NVIDIA GB10,GPU,w32,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,33600,0.000265039,126773798.301,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,33600,0.000138864,241963370.419,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,33600,0.000291598,115227127.627,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,33600,0.000138911,241881473.028,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,33600,0.000290062,115837302.941,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,33600,0.000070752,474898179.248,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,33600,0.000228559,147007992.467,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,33600,0.000064385,521860741.948,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,33600,0.000218163,154013285.733,0
library,unknown,host-cpu,gmp-1t,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.006344236,3940584.801,0
library,unknown,host-cpu,gmp-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.001404360,17801703.331,0
library,unknown,host-cpu,openssl-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.001012187,24698993.367,0
library,NVIDIA GB10,gpu,cgbn,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,700000,0.002374688,294775566.306,0
opencl-kernel,NVIDIA GB10,GPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,33600,0.000404958,82971567.005,0
opencl-e2e,NVIDIA GB10,GPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,33600,0.000525677,63917577.421,0
opencl-kernel,NVIDIA GB10,GPU,w16,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,33600,0.000134975,248934971.118,0
opencl-e2e,NVIDIA GB10,GPU,w16,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,33600,0.000252895,132861464.773,0
opencl-kernel,NVIDIA GB10,GPU,w32,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,33600,0.000037680,891719705.073,0
opencl-e2e,NVIDIA GB10,GPU,w32,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,33600,0.000156223,215077161.610,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,33600,0.000040016,839664006.427,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,33600,0.000158351,212186848.013,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,33600,0.000035312,951517961.318,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,33600,0.000156111,215231454.522,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,33600,0.000040320,833333424.161,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,33600,0.000157391,213481079.009,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,33600,0.000028497,1179071438.187,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,33600,0.000146705,229031044.945,0
library,unknown,host-cpu,gmp-1t,brainpoolP512r1,512,COMPARE,25000,0.000177951,140488118.842,0
library,unknown,host-cpu,gmp-nt,brainpoolP512r1,512,COMPARE,25000,0.000955642,26160424.166,0
library,unknown,host-cpu,openssl-nt,brainpoolP512r1,512,COMPARE,25000,0.000558220,44785210.487,0
library,NVIDIA GB10,gpu,cgbn,brainpoolP512r1,512,COMPARE,700000,0.000247200,2831715210.356,0
opencl-kernel,NVIDIA GB10,GPU,w8,brainpoolP512r1,512,COMPARE,33600,0.000117983,284786792.151,0
opencl-e2e,NVIDIA GB10,GPU,w8,brainpoolP512r1,512,COMPARE,33600,0.000241119,139350283.621,0
opencl-kernel,NVIDIA GB10,GPU,w16,brainpoolP512r1,512,COMPARE,33600,0.000062591,536818488.320,0
opencl-e2e,NVIDIA GB10,GPU,w16,brainpoolP512r1,512,COMPARE,33600,0.000179919,186750697.458,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,brainpoolP512r1,512,COMPARE,33600,0.000030912,1086956759.604,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,brainpoolP512r1,512,COMPARE,33600,0.000156831,214243362.588,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,brainpoolP512r1,512,COMPARE,33600,0.000030640,1096605622.875,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,brainpoolP512r1,512,COMPARE,33600,0.000148495,226270238.091,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,brainpoolP512r1,512,COMPARE,33600,0.000010352,3245749719.321,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,brainpoolP512r1,512,COMPARE,33600,0.000132271,254023933.118,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,brainpoolP512r1,512,COMPARE,33600,0.000009809,3425427394.420,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,brainpoolP512r1,512,COMPARE,33600,0.000131650,255222154.980,0
library,unknown,host-cpu,gmp-1t,brainpoolP512r1,512,REDUCE,3125,0.000055280,56530400.767,0
library,unknown,host-cpu,gmp-nt,brainpoolP512r1,512,REDUCE,3125,0.000671132,4656312.008,0
library,unknown,host-cpu,openssl-nt,brainpoolP512r1,512,REDUCE,3125,0.000568029,5501479.629,0
library,NVIDIA GB10,gpu,cgbn,brainpoolP512r1,512,REDUCE,700000,0.000243744,2871865563.870,0
opencl-kernel,NVIDIA GB10,GPU,w8,brainpoolP512r1,512,REDUCE,33600,0.000516701,65027935.661,0
opencl-e2e,NVIDIA GB10,GPU,w8,brainpoolP512r1,512,REDUCE,33600,0.000636524,52786698.005,0
opencl-kernel,NVIDIA GB10,GPU,w16,brainpoolP512r1,512,REDUCE,33600,0.000398798,84253180.450,0
opencl-e2e,NVIDIA GB10,GPU,w16,brainpoolP512r1,512,REDUCE,33600,0.000511630,65672458.393,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,brainpoolP512r1,512,REDUCE,33600,0.000107152,313573214.973,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,brainpoolP512r1,512,REDUCE,33600,0.000230623,145692319.294,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,brainpoolP512r1,512,REDUCE,33600,0.000101807,330036234.110,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,brainpoolP512r1,512,REDUCE,33600,0.000218575,153722974.571,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,brainpoolP512r1,512,REDUCE,33600,0.000114767,292767091.540,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,brainpoolP512r1,512,REDUCE,33600,0.000232031,144808235.982,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,brainpoolP512r1,512,REDUCE,33600,0.000109202,307686684.645,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,brainpoolP512r1,512,REDUCE,33600,0.000230003,146085049.287,0
library,unknown,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL,1562,0.000191231,8168131.895,0
library,unknown,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL,1562,0.000823947,1895752.993,0
library,unknown,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL,1562,0.000772028,2023242.687,0
library,NVIDIA GB10,gpu,cgbn,brainpoolP512r1,512,MODMUL,700000,0.002380640,294038577.861,0
opencl-kernel,NVIDIA GB10,GPU,w8,brainpoolP512r1,512,MODMUL,33600,0.001539431,21826246.354,0
opencl-e2e,NVIDIA GB10,GPU,w8,brainpoolP512r1,512,MODMUL,33600,0.001666982,20156186.409,0
opencl-kernel,NVIDIA GB10,GPU,w16,brainpoolP512r1,512,MODMUL,33600,0.001082651,31034931.934,0
opencl-e2e,NVIDIA GB10,GPU,w16,brainpoolP512r1,512,MODMUL,33600,0.001201995,27953527.410,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,brainpoolP512r1,512,MODMUL,33600,0.000395230,85013788.783,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,brainpoolP512r1,512,MODMUL,33600,0.000514029,65365961.671,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,brainpoolP512r1,512,MODMUL,33600,0.000373614,89932392.150,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,brainpoolP512r1,512,MODMUL,33600,0.000491229,68399870.200,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,brainpoolP512r1,512,MODMUL,33600,0.000464253,72374330.905,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,brainpoolP512r1,512,MODMUL,33600,0.000587949,57147813.314,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,brainpoolP512r1,512,MODMUL,33600,0.000449365,74772177.704,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,brainpoolP512r1,512,MODMUL,33600,0.000573832,58553722.706,0
library,unknown,host-cpu,gmp-1t,brainpoolP512r1,512,MODEXP,390,0.013924400,28008.388,0
library,unknown,host-cpu,gmp-nt,brainpoolP512r1,512,MODEXP,390,0.002195315,177651.043,0
library,unknown,host-cpu,openssl-nt,brainpoolP512r1,512,MODEXP,390,0.002044900,190718.373,0
library,NVIDIA GB10,gpu,cgbn,brainpoolP512r1,512,MODEXP,700000,0.211561188,3308735.438,0
opencl-kernel,NVIDIA GB10,GPU,w8,brainpoolP512r1,512,MODEXP,33600,0.279386457,120263.524,0
opencl-e2e,NVIDIA GB10,GPU,w8,brainpoolP512r1,512,MODEXP,33600,0.279473545,120226.049,0
opencl-kernel,NVIDIA GB10,GPU,w16,brainpoolP512r1,512,MODEXP,33600,0.026788501,1254269.509,0
opencl-e2e,NVIDIA GB10,GPU,w16,brainpoolP512r1,512,MODEXP,33600,0.026959364,1246320.203,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,brainpoolP512r1,512,MODEXP,33600,0.024534743,1369486.528,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,brainpoolP512r1,512,MODEXP,33600,0.024618216,1364843.009,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,brainpoolP512r1,512,MODEXP,33600,0.012919776,2600664.285,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,brainpoolP512r1,512,MODEXP,33600,0.013050671,2574580.265,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,brainpoolP512r1,512,MODEXP,33600,0.026013742,1291625.019,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,brainpoolP512r1,512,MODEXP,33600,0.026123630,1286191.850,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,brainpoolP512r1,512,MODEXP,33600,0.013446475,2498796.153,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,brainpoolP512r1,512,MODEXP,33600,0.013581068,2474032.235,0
library,unknown,host-cpu,gmp-1t,brainpoolP512r1,512,EXPONENTIATION,390,0.002479282,157303.607,0
library,unknown,host-cpu,gmp-nt,brainpoolP512r1,512,EXPONENTIATION,390,0.000924715,421751.566,0
library,unknown,host-cpu,openssl-nt,brainpoolP512r1,512,EXPONENTIATION,390,0.003822970,102014.926,0
opencl-kernel,NVIDIA GB10,GPU,w8,brainpoolP512r1,512,EXPONENTIATION,33600,0.399943716,84011.821,0
opencl-e2e,NVIDIA GB10,GPU,w8,brainpoolP512r1,512,EXPONENTIATION,33600,0.397664817,84493.268,0
opencl-kernel,NVIDIA GB10,GPU,w16,brainpoolP512r1,512,EXPONENTIATION,33600,0.102378359,328194.360,0
opencl-e2e,NVIDIA GB10,GPU,w16,brainpoolP512r1,512,EXPONENTIATION,33600,0.102455446,327947.428,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,brainpoolP512r1,512,EXPONENTIATION,33600,0.031951683,1051587.799,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,brainpoolP512r1,512,EXPONENTIATION,33600,0.032009603,1049684.996,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,brainpoolP512r1,512,EXPONENTIATION,33600,0.030831992,1089777.138,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,brainpoolP512r1,512,EXPONENTIATION,33600,0.031321797,1072735.386,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,brainpoolP512r1,512,EXPONENTIATION,33600,0.033159643,1013279.908,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,brainpoolP512r1,512,EXPONENTIATION,33600,0.032905740,1021098.447,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,brainpoolP512r1,512,EXPONENTIATION,33600,0.032025366,1049168.337,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,brainpoolP512r1,512,EXPONENTIATION,33600,0.032463164,1035019.261,0
library,unknown,host-cpu,gmp-1t,brainpoolP512r1,512,DIVIDE,3125,0.000065312,47847250.539,0
library,unknown,host-cpu,gmp-nt,brainpoolP512r1,512,DIVIDE,3125,0.000631245,4950534.272,0
library,unknown,host-cpu,openssl-nt,brainpoolP512r1,512,DIVIDE,3125,0.000702924,4445715.369,0
library,NVIDIA GB10,gpu,cgbn,brainpoolP512r1,512,DIVIDE,700000,0.002377216,294462093.474,0
opencl-kernel,NVIDIA GB10,GPU,w8,brainpoolP512r1,512,DIVIDE,33600,0.001080362,31100686.758,0
opencl-e2e,NVIDIA GB10,GPU,w8,brainpoolP512r1,512,DIVIDE,33600,0.001226105,27403851.718,0
opencl-kernel,NVIDIA GB10,GPU,w16,brainpoolP512r1,512,DIVIDE,33600,0.000947180,35473721.982,0
opencl-e2e,NVIDIA GB10,GPU,w16,brainpoolP512r1,512,DIVIDE,33600,0.001101579,30501670.567,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,brainpoolP512r1,512,DIVIDE,33600,0.000316271,106238007.813,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,brainpoolP512r1,512,DIVIDE,33600,0.000475405,70676580.669,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,brainpoolP512r1,512,DIVIDE,33600,0.000301198,111554525.346,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,brainpoolP512r1,512,DIVIDE,33600,0.000454973,73850534.979,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,brainpoolP512r1,512,DIVIDE,33600,0.000288719,116376133.966,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,brainpoolP512r1,512,DIVIDE,33600,0.000441390,76123155.982,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,brainpoolP512r1,512,DIVIDE,33600,0.000277028,121287380.175,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,brainpoolP512r1,512,DIVIDE,33600,0.000428245,78459758.357,0
library,unknown,host-cpu,gmp-1t,brainpoolP512r1,512,ISQRT,781,0.000053247,14667492.915,0
library,unknown,host-cpu,gmp-nt,brainpoolP512r1,512,ISQRT,781,0.000796155,980964.762,0
opencl-kernel,NVIDIA GB10,GPU,w8,brainpoolP512r1,512,ISQRT,33600,0.022180689,1514831.212,0
opencl-e2e,NVIDIA GB10,GPU,w8,brainpoolP512r1,512,ISQRT,33600,0.022355376,1502994.179,0
opencl-kernel,NVIDIA GB10,GPU,w16,brainpoolP512r1,512,ISQRT,33600,0.020581553,1632529.869,0
opencl-e2e,NVIDIA GB10,GPU,w16,brainpoolP512r1,512,ISQRT,33600,0.020689968,1623975.445,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,brainpoolP512r1,512,ISQRT,33600,0.004010716,8377556.522,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,brainpoolP512r1,512,ISQRT,33600,0.004126444,8142604.154,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,brainpoolP512r1,512,ISQRT,33600,0.003833757,8764248.750,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,brainpoolP512r1,512,ISQRT,33600,0.003953645,8498486.827,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,brainpoolP512r1,512,ISQRT,33600,0.005093991,6596006.943,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,brainpoolP512r1,512,ISQRT,33600,0.005210742,6448217.925,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,brainpoolP512r1,512,ISQRT,33600,0.004965631,6766511.647,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,brainpoolP512r1,512,ISQRT,33600,0.005084913,6607782.663,0
library,unknown,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL_R2,25000,0.003081998,8111621.092,0
library,unknown,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL_R2,25000,0.000965259,25899784.564,0
library,unknown,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL_R2,25000,0.001776246,14074627.013,0
library,NVIDIA GB10,gpu,cgbn,brainpoolP512r1,512,MODMUL_R2,700000,0.000170240,4111842105.263,0
opencl-kernel,NVIDIA GB10,GPU,w8,brainpoolP512r1,512,MODMUL_R2,33600,0.000587325,57208531.156,0
opencl-e2e,NVIDIA GB10,GPU,w8,brainpoolP512r1,512,MODMUL_R2,33600,0.000709948,47327409.962,0
opencl-kernel,NVIDIA GB10,GPU,w16,brainpoolP512r1,512,MODMUL_R2,33600,0.000119984,280037328.010,0
opencl-e2e,NVIDIA GB10,GPU,w16,brainpoolP512r1,512,MODMUL_R2,33600,0.000242559,138522994.852,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,brainpoolP512r1,512,MODMUL_R2,33600,0.000071215,471810669.278,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,brainpoolP512r1,512,MODMUL_R2,33600,0.000192479,174564500.092,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,brainpoolP512r1,512,MODMUL_R2,33600,0.000054607,615305642.495,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,brainpoolP512r1,512,MODMUL_R2,33600,0.000173599,193549509.183,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,brainpoolP512r1,512,MODMUL_R2,33600,0.000068928,487465190.574,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,brainpoolP512r1,512,MODMUL_R2,33600,0.000186383,180273955.296,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,brainpoolP512r1,512,MODMUL_R2,33600,0.000047457,708009400.877,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,brainpoolP512r1,512,MODMUL_R2,33600,0.000169762,197924149.408,0
library,unknown,host-cpu,gmp-1t,p1024,1024,ADD,12500,0.000171087,73062243.909,0
library,unknown,host-cpu,gmp-nt,p1024,1024,ADD,12500,0.000822555,15196552.256,0
library,unknown,host-cpu,openssl-nt,p1024,1024,ADD,12500,0.001075385,11623744.131,0
library,NVIDIA GB10,gpu,cgbn,p1024,1024,ADD,700000,0.000460992,1518464528.669,0
opencl-kernel,NVIDIA GB10,GPU,w8,p1024,1024,ADD,33600,0.000345662,97204787.005,0
opencl-e2e,NVIDIA GB10,GPU,w8,p1024,1024,ADD,33600,0.000562253,59759574.354,0
opencl-kernel,NVIDIA GB10,GPU,w16,p1024,1024,ADD,33600,0.000173375,193799555.685,0
opencl-e2e,NVIDIA GB10,GPU,w16,p1024,1024,ADD,33600,0.000389823,86192963.147,0
opencl-kernel,NVIDIA GB10,GPU,w32,p1024,1024,ADD,33600,0.000095487,351880396.879,0
opencl-e2e,NVIDIA GB10,GPU,w32,p1024,1024,ADD,33600,0.000314190,106941656.512,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,p1024,1024,ADD,33600,0.000095951,350178764.271,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,p1024,1024,ADD,33600,0.000316847,106044873.614,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,p1024,1024,ADD,33600,0.000095824,350642852.605,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,p1024,1024,ADD,33600,0.000314591,106805343.617,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,p1024,1024,ADD,33600,0.000022816,1472650378.039,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,p1024,1024,ADD,33600,0.000252062,133300535.207,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,p1024,1024,ADD,33600,0.000023520,1428571510.621,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,p1024,1024,ADD,33600,0.000263283,127619329.725,0
library,unknown,host-cpu,gmp-1t,p1024,1024,SUBTRACT,12500,0.000123008,101619411.608,0
library,unknown,host-cpu,gmp-nt,p1024,1024,SUBTRACT,12500,0.000595612,20986816.910,0
library,unknown,host-cpu,openssl-nt,p1024,1024,SUBTRACT,12500,0.000799324,15638214.297,0
library,NVIDIA GB10,gpu,cgbn,p1024,1024,SUBTRACT,700000,0.004714080,148491328.106,0
opencl-kernel,NVIDIA GB10,GPU,w8,p1024,1024,SUBTRACT,33600,0.000324814,103443817.617,0
opencl-e2e,NVIDIA GB10,GPU,w8,p1024,1024,SUBTRACT,33600,0.000540237,62194925.482,0
opencl-kernel,NVIDIA GB10,GPU,w16,p1024,1024,SUBTRACT,33600,0.000168479,199431378.334,0
opencl-e2e,NVIDIA GB10,GPU,w16,p1024,1024,SUBTRACT,33600,0.000387550,86698491.198,0
opencl-kernel,NVIDIA GB10,GPU,w32,p1024,1024,SUBTRACT,33600,0.000095183,353004210.463,0
opencl-e2e,NVIDIA GB10,GPU,w32,p1024,1024,SUBTRACT,33600,0.000317791,105729867.605,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,p1024,1024,SUBTRACT,33600,0.000095984,350058303.520,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,p1024,1024,SUBTRACT,33600,0.000316766,106071988.918,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,p1024,1024,SUBTRACT,33600,0.000095504,351817773.490,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,p1024,1024,SUBTRACT,33600,0.000315903,106361762.257,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,p1024,1024,SUBTRACT,33600,0.000022688,1480959561.037,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,p1024,1024,SUBTRACT,33600,0.000253695,132442497.197,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,p1024,1024,SUBTRACT,33600,0.000023840,1409396253.082,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,p1024,1024,SUBTRACT,33600,0.000252707,132960304.335,0
library,unknown,host-cpu,gmp-1t,p1024,1024,ADDMOD,12500,0.000430653,29025688.714,0
library,unknown,host-cpu,gmp-nt,p1024,1024,ADDMOD,12500,0.000848604,14730074.292,0
library,unknown,host-cpu,openssl-nt,p1024,1024,ADDMOD,12500,0.001002698,12466365.724,0
library,NVIDIA GB10,gpu,cgbn,p1024,1024,ADDMOD,700000,0.000246464,2840171384.056,0
opencl-kernel,NVIDIA GB10,GPU,w8,p1024,1024,ADDMOD,33600,0.000483085,69552978.097,0
opencl-e2e,NVIDIA GB10,GPU,w8,p1024,1024,ADDMOD,33600,0.000696172,48263935.338,0
opencl-kernel,NVIDIA GB10,GPU,w16,p1024,1024,ADDMOD,33600,0.000245183,137040491.927,0
opencl-e2e,NVIDIA GB10,GPU,w16,p1024,1024,ADDMOD,33600,0.000461341,72831157.569,0
opencl-kernel,NVIDIA GB10,GPU,w32,p1024,1024,ADDMOD,33600,0.000118335,283939672.583,0
opencl-e2e,NVIDIA GB10,GPU,w32,p1024,1024,ADDMOD,33600,0.000339518,98963826.470,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,p1024,1024,ADDMOD,33600,0.000109455,306975482.681,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,p1024,1024,ADDMOD,33600,0.000338815,99169165.286,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,p1024,1024,ADDMOD,33600,0.000107760,311804001.038,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,p1024,1024,ADDMOD,33600,0.000340159,98777334.313,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,p1024,1024,ADDMOD,33600,0.000026608,1262777999.270,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,p1024,1024,ADDMOD,33600,0.000256255,131119396.597,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,p1024,1024,ADDMOD,33600,0.000025904,1297096481.863,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,p1024,1024,ADDMOD,33600,0.000255923,131289486.764,0
library,unknown,host-cpu,gmp-1t,p1024,1024,SUBTRACTMOD,12500,0.000323054,38693223.026,0
library,unknown,host-cpu,gmp-nt,p1024,1024,SUBTRACTMOD,12500,0.000843580,14817800.374,0
library,unknown,host-cpu,openssl-nt,p1024,1024,SUBTRACTMOD,12500,0.000603437,20714672.724,0
library,NVIDIA GB10,gpu,cgbn,p1024,1024,SUBTRACTMOD,700000,0.000246464,2840171384.056,0
opencl-kernel,NVIDIA GB10,GPU,w8,p1024,1024,SUBTRACTMOD,33600,0.000468670,71692233.609,0
opencl-e2e,NVIDIA GB10,GPU,w8,p1024,1024,SUBTRACTMOD,33600,0.000688124,48828409.050,0
opencl-kernel,NVIDIA GB10,GPU,w16,p1024,1024,SUBTRACTMOD,33600,0.000247359,135834960.890,0
opencl-e2e,NVIDIA GB10,GPU,w16,p1024,1024,SUBTRACTMOD,33600,0.000471662,71237455.777,0
opencl-kernel,NVIDIA GB10,GPU,w32,p1024,1024,SUBTRACTMOD,33600,0.000121376,276825738.524,0
opencl-e2e,NVIDIA GB10,GPU,w32,p1024,1024,SUBTRACTMOD,33600,0.000340286,98740470.662,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,p1024,1024,SUBTRACTMOD,33600,0.000109568,306658904.041,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,p1024,1024,SUBTRACTMOD,33600,0.000332094,101176172.600,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,p1024,1024,SUBTRACTMOD,33600,0.000106992,314042182.266,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,p1024,1024,SUBTRACTMOD,33600,0.000328751,102205010.931,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,p1024,1024,SUBTRACTMOD,33600,0.000026464,1269648799.005,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,p1024,1024,SUBTRACTMOD,33600,0.000255903,131299752.201,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,p1024,1024,SUBTRACTMOD,33600,0.000025904,1297097210.525,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,p1024,1024,SUBTRACTMOD,33600,0.000255364,131576889.755,0
library,unknown,host-cpu,gmp-1t,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.001467223,8519495.657,0
library,unknown,host-cpu,gmp-nt,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.000993546,12581199.102,0
library,unknown,host-cpu,openssl-nt,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.000660572,18922994.089,0
opencl-kernel,NVIDIA GB10,GPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,33600,0.036417631,922630.031,0
opencl-e2e,NVIDIA GB10,GPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,33600,0.036536990,919615.984,0
opencl-kernel,NVIDIA GB10,GPU,w16,p1024,1024,MULTIPLYOPERANDSCANNING,33600,0.010133761,3315649.542,0
opencl-e2e,NVIDIA GB10,GPU,w16,p1024,1024,MULTIPLYOPERANDSCANNING,33600,0.010570288,3178721.338,0
opencl-kernel,NVIDIA GB10,GPU,w32,p1024,1024,MULTIPLYOPERANDSCANNING,33600,0.002557268,13139021.787,0
opencl-e2e,NVIDIA GB10,GPU,w32,p1024,1024,MULTIPLYOPERANDSCANNING,33600,0.002855762,11765686.316,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,p1024,1024,MULTIPLYOPERANDSCANNING,33600,0.000332702,100991279.021,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,p1024,1024,MULTIPLYOPERANDSCANNING,33600,0.000631021,53247040.463,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,p1024,1024,MULTIPLYOPERANDSCANNING,33600,0.000343806,97729535.248,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,p1024,1024,MULTIPLYOPERANDSCANNING,33600,0.000630029,53330878.677,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,p1024,1024,MULTIPLYOPERANDSCANNING,33600,0.000184735,181882148.855,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,p1024,1024,MULTIPLYOPERANDSCANNING,33600,0.000476206,70557700.118,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,p1024,1024,MULTIPLYOPERANDSCANNING,33600,0.000187298,179393274.063,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,p1024,1024,MULTIPLYOPERANDSCANNING,33600,0.000477862,70313185.779,0
library,unknown,host-cpu,gmp-1t,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.001469544,8506039.945,0
library,unknown,host-cpu,gmp-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.000861019,14517681.990,0
library,unknown,host-cpu,openssl-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.001181977,10575501.886,0
library,NVIDIA GB10,gpu,cgbn,p1024,1024,MULTIPLYPRODUCTSCANNING,700000,0.000449600,1556939501.779,0
opencl-kernel,NVIDIA GB10,GPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,33600,0.009573209,3509794.888,0
opencl-e2e,NVIDIA GB10,GPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,33600,0.009867623,3405075.367,0
opencl-kernel,NVIDIA GB10,GPU,w16,p1024,1024,MULTIPLYPRODUCTSCANNING,33600,0.002453877,13692617.798,0
opencl-e2e,NVIDIA GB10,GPU,w16,p1024,1024,MULTIPLYPRODUCTSCANNING,33600,0.002747764,12228124.413,0
opencl-kernel,NVIDIA GB10,GPU,w32,p1024,1024,MULTIPLYPRODUCTSCANNING,33600,0.000640653,52446488.108,0
opencl-e2e,NVIDIA GB10,GPU,w32,p1024,1024,MULTIPLYPRODUCTSCANNING,33600,0.000939419,35766787.834,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,p1024,1024,MULTIPLYPRODUCTSCANNING,33600,0.000644989,52093912.361,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,p1024,1024,MULTIPLYPRODUCTSCANNING,33600,0.000936123,35892719.070,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,p1024,1024,MULTIPLYPRODUCTSCANNING,33600,0.000639197,52565954.289,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,p1024,1024,MULTIPLYPRODUCTSCANNING,33600,0.000935148,35930141.109,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,p1024,1024,MULTIPLYPRODUCTSCANNING,33600,0.000291662,115201843.366,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,p1024,1024,MULTIPLYPRODUCTSCANNING,33600,0.000591805,56775459.430,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,p1024,1024,MULTIPLYPRODUCTSCANNING,33600,0.000265747,126436051.744,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,p1024,1024,MULTIPLYPRODUCTSCANNING,33600,0.000567784,59177432.730,0
library,unknown,host-cpu,gmp-1t,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.009957479,1255337.823,0
library,unknown,host-cpu,gmp-nt,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.001762198,7093414.048,0
library,unknown,host-cpu,openssl-nt,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.000926602,13490150.155,0
library,NVIDIA GB10,gpu,cgbn,p1024,1024,MONTGOMERYMULTIPLICATION,700000,0.008761887,79891466.302,0
opencl-kernel,NVIDIA GB10,GPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,33600,0.003194557,10517890.247,0
opencl-e2e,NVIDIA GB10,GPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,33600,0.003384573,9927397.068,0
opencl-kernel,NVIDIA GB10,GPU,w16,p1024,1024,MONTGOMERYMULTIPLICATION,33600,0.000449917,74680440.611,0
opencl-e2e,NVIDIA GB10,GPU,w16,p1024,1024,MONTGOMERYMULTIPLICATION,33600,0.000669229,50207028.952,0
opencl-kernel,NVIDIA GB10,GPU,w32,p1024,1024,MONTGOMERYMULTIPLICATION,33600,0.000167487,200612584.859,0
opencl-e2e,NVIDIA GB10,GPU,w32,p1024,1024,MONTGOMERYMULTIPLICATION,33600,0.000387422,86727135.130,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,p1024,1024,MONTGOMERYMULTIPLICATION,33600,0.000186287,180366855.589,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,p1024,1024,MONTGOMERYMULTIPLICATION,33600,0.000414430,81075213.854,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,p1024,1024,MONTGOMERYMULTIPLICATION,33600,0.000146464,229407919.155,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,p1024,1024,MONTGOMERYMULTIPLICATION,33600,0.000379198,88608062.753,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,p1024,1024,MONTGOMERYMULTIPLICATION,33600,0.000124688,269472607.243,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,p1024,1024,MONTGOMERYMULTIPLICATION,33600,0.000353278,95109235.272,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,p1024,1024,MONTGOMERYMULTIPLICATION,33600,0.000080945,415096599.268,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,p1024,1024,MONTGOMERYMULTIPLICATION,33600,0.000311732,107784896.789,0
library,unknown,host-cpu,gmp-1t,p1024,1024,COMPARE,12500,0.000078463,159310766.893,0
library,unknown,host-cpu,gmp-nt,p1024,1024,COMPARE,12500,0.000647244,19312654.888,0
library,unknown,host-cpu,openssl-nt,p1024,1024,COMPARE,12500,0.000953530,13109183.700,0
library,NVIDIA GB10,gpu,cgbn,p1024,1024,COMPARE,700000,0.000458752,1525878906.250,0
opencl-kernel,NVIDIA GB10,GPU,w8,p1024,1024,COMPARE,33600,0.000320622,104796302.521,0
opencl-e2e,NVIDIA GB10,GPU,w8,p1024,1024,COMPARE,33600,0.000540109,62209666.528,0
opencl-kernel,NVIDIA GB10,GPU,w16,p1024,1024,COMPARE,33600,0.000149151,225275049.671,0
opencl-e2e,NVIDIA GB10,GPU,w16,p1024,1024,COMPARE,33600,0.000373070,90063527.605,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,p1024,1024,COMPARE,33600,0.000080112,419412800.738,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,p1024,1024,COMPARE,33600,0.000303726,110626023.665,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,p1024,1024,COMPARE,33600,0.000081711,411205284.417,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,p1024,1024,COMPARE,33600,0.000307566,109244842.847,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,p1024,1024,COMPARE,33600,0.000014288,2351623461.426,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,p1024,1024,COMPARE,33600,0.000239742,140150660.189,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,p1024,1024,COMPARE,33600,0.000014256,2356902020.925,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,p1024,1024,COMPARE,33600,0.000239571,140250696.088,0
library,unknown,host-cpu,gmp-1t,p1024,1024,REDUCE,1562,0.000019216,81286456.598,0
library,unknown,host-cpu,gmp-nt,p1024,1024,REDUCE,1562,0.000783500,1993618.380,0
library,unknown,host-cpu,openssl-nt,p1024,1024,REDUCE,1562,0.000735948,2122432.541,0
library,NVIDIA GB10,gpu,cgbn,p1024,1024,REDUCE,700000,0.000454048,1541687222.496,0
opencl-kernel,NVIDIA GB10,GPU,w8,p1024,1024,REDUCE,33600,0.002750000,12218181.820,0
opencl-e2e,NVIDIA GB10,GPU,w8,p1024,1024,REDUCE,33600,0.002971791,11306313.291,0
opencl-kernel,NVIDIA GB10,GPU,w16,p1024,1024,REDUCE,33600,0.001473066,22809568.550,0
opencl-e2e,NVIDIA GB10,GPU,w16,p1024,1024,REDUCE,33600,0.001691272,19866703.853,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,p1024,1024,REDUCE,33600,0.000387646,86677017.409,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,p1024,1024,REDUCE,33600,0.000613117,54801938.836,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,p1024,1024,REDUCE,33600,0.000375503,89479976.389,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,p1024,1024,REDUCE,33600,0.000605981,55447281.978,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,p1024,1024,REDUCE,33600,0.000317230,105916845.609,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,p1024,1024,REDUCE,33600,0.000541758,62020312.114,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,p1024,1024,REDUCE,33600,0.000310820,108101150.573,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,p1024,1024,REDUCE,33600,0.000535511,62743808.683,0
library,unknown,host-cpu,gmp-1t,p1024,1024,MODMUL,781,0.000261903,2982020.080,0
library,unknown,host-cpu,gmp-nt,p1024,1024,MODMUL,781,0.000720956,1083283.864,0
library,unknown,host-cpu,openssl-nt,p1024,1024,MODMUL,781,0.000866299,901536.305,0
library,NVIDIA GB10,gpu,cgbn,p1024,1024,MODMUL,700000,0.004740128,147675337.037,0
opencl-kernel,NVIDIA GB10,GPU,w8,p1024,1024,MODMUL,33600,0.012086043,2780066.229,0
opencl-e2e,NVIDIA GB10,GPU,w8,p1024,1024,MODMUL,33600,0.012319929,2727288.446,0
opencl-kernel,NVIDIA GB10,GPU,w16,p1024,1024,MODMUL,33600,0.004119117,8157088.032,0
opencl-e2e,NVIDIA GB10,GPU,w16,p1024,1024,MODMUL,33600,0.004348620,7726589.123,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,p1024,1024,MODMUL,33600,0.001674984,20059893.113,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,p1024,1024,MODMUL,33600,0.001904215,17645066.393,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,p1024,1024,MODMUL,33600,0.001619672,20744940.825,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,p1024,1024,MODMUL,33600,0.001847190,18189791.013,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,p1024,1024,MODMUL,33600,0.001629208,20623517.731,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,p1024,1024,MODMUL,33600,0.001857319,18090591.843,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,p1024,1024,MODMUL,33600,0.001600116,20998477.465,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,p1024,1024,MODMUL,33600,0.001824295,18418073.804,0
library,unknown,host-cpu,gmp-1t,p1024,1024,MODEXP,195,0.039880619,4889.593,0
library,unknown,host-cpu,gmp-nt,p1024,1024,MODEXP,195,0.005969421,32666.485,0
library,unknown,host-cpu,openssl-nt,p1024,1024,MODEXP,195,0.005064467,38503.558,0
library,NVIDIA GB10,gpu,cgbn,p1024,1024,MODEXP,700000,1.212423444,577356.041,0
opencl-kernel,NVIDIA GB10,GPU,w8,p1024,1024,MODEXP,33600,5.609375797,5989.971,0
opencl-e2e,NVIDIA GB10,GPU,w8,p1024,1024,MODEXP,33600,5.607541333,5991.931,0
opencl-kernel,NVIDIA GB10,GPU,w16,p1024,1024,MODEXP,33600,0.190004481,176837.935,0
opencl-e2e,NVIDIA GB10,GPU,w16,p1024,1024,MODEXP,33600,0.191068267,175853.377,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,p1024,1024,MODEXP,33600,0.197245415,170346.165,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,p1024,1024,MODEXP,33600,0.197688261,169964.569,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,p1024,1024,MODEXP,33600,0.111594520,301090.054,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,p1024,1024,MODEXP,33600,0.112197684,299471.422,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,p1024,1024,MODEXP,33600,0.197258377,170334.971,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,p1024,1024,MODEXP,33600,0.197419848,170195.653,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,p1024,1024,MODEXP,33600,0.112627541,298328.452,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,p1024,1024,MODEXP,33600,0.113114954,297042.953,0
library,unknown,host-cpu,gmp-1t,p1024,1024,EXPONENTIATION,195,0.004794949,40667.794,0
library,unknown,host-cpu,gmp-nt,p1024,1024,EXPONENTIATION,195,0.001381736,141126.814,0
library,unknown,host-cpu,openssl-nt,p1024,1024,EXPONENTIATION,195,0.013743985,14188.025,0
opencl-kernel,NVIDIA GB10,GPU,w8,p1024,1024,EXPONENTIATION,33600,3.384674711,9927.099,0
opencl-e2e,NVIDIA GB10,GPU,w8,p1024,1024,EXPONENTIATION,33600,3.356873452,10009.314,0
opencl-kernel,NVIDIA GB10,GPU,w16,p1024,1024,EXPONENTIATION,33600,0.825917655,40682.022,0
opencl-e2e,NVIDIA GB10,GPU,w16,p1024,1024,EXPONENTIATION,33600,0.831163487,40425.260,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,p1024,1024,EXPONENTIATION,33600,0.234713710,143153.120,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,p1024,1024,EXPONENTIATION,33600,0.234895982,143042.038,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,p1024,1024,EXPONENTIATION,33600,0.217217949,154683.350,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,p1024,1024,EXPONENTIATION,33600,0.216358321,155297.933,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,p1024,1024,EXPONENTIATION,33600,0.228470669,147064.830,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,p1024,1024,EXPONENTIATION,33600,0.228311118,147167.603,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,p1024,1024,EXPONENTIATION,33600,0.210321037,159755.774,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,p1024,1024,EXPONENTIATION,33600,0.210769475,159415.874,0
library,unknown,host-cpu,gmp-1t,p1024,1024,DIVIDE,1562,0.000036480,42817991.190,0
library,unknown,host-cpu,gmp-nt,p1024,1024,DIVIDE,1562,0.000712252,2193044.068,0
library,unknown,host-cpu,openssl-nt,p1024,1024,DIVIDE,1562,0.000625004,2499183.996,0
library,NVIDIA GB10,gpu,cgbn,p1024,1024,DIVIDE,700000,0.004723648,148190551.032,0
opencl-kernel,NVIDIA GB10,GPU,w8,p1024,1024,DIVIDE,33600,0.051419434,653449.433,0
opencl-e2e,NVIDIA GB10,GPU,w8,p1024,1024,DIVIDE,33600,0.051847703,648051.853,0
opencl-kernel,NVIDIA GB10,GPU,w16,p1024,1024,DIVIDE,33600,0.005712181,5882166.539,0
opencl-e2e,NVIDIA GB10,GPU,w16,p1024,1024,DIVIDE,33600,0.006008564,5592018.326,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,p1024,1024,DIVIDE,33600,0.001264537,26570990.000,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,p1024,1024,DIVIDE,33600,0.001549625,21682665.189,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,p1024,1024,DIVIDE,33600,0.001191770,28193359.385,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,p1024,1024,DIVIDE,33600,0.001498633,22420432.396,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,p1024,1024,DIVIDE,33600,0.001234810,27210664.011,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,p1024,1024,DIVIDE,33600,0.001525129,22030923.251,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,p1024,1024,DIVIDE,33600,0.001179726,28481189.493,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,p1024,1024,DIVIDE,33600,0.001472322,22821094.938,0
library,unknown,host-cpu,gmp-1t,p1024,1024,ISQRT,390,0.000056576,6893382.558,0
library,unknown,host-cpu,gmp-nt,p1024,1024,ISQRT,390,0.000652941,597297.457,0
opencl-kernel,NVIDIA GB10,GPU,w8,p1024,1024,ISQRT,33600,0.427648256,78569.244,0
opencl-e2e,NVIDIA GB10,GPU,w8,p1024,1024,ISQRT,33600,0.429547877,78221.781,0
opencl-kernel,NVIDIA GB10,GPU,w16,p1024,1024,ISQRT,33600,0.146866038,228779.917,0
opencl-e2e,NVIDIA GB10,GPU,w16,p1024,1024,ISQRT,33600,0.146926709,228685.446,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,p1024,1024,ISQRT,33600,0.028392469,1183412.404,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,p1024,1024,ISQRT,33600,0.028618323,1174072.989,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,p1024,1024,ISQRT,33600,0.025406387,1322502.093,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,p1024,1024,ISQRT,33600,0.025667744,1309035.963,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,p1024,1024,ISQRT,33600,0.027308584,1230382.359,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,p1024,1024,ISQRT,33600,0.027655430,1214951.277,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,p1024,1024,ISQRT,33600,0.026217378,1281592.690,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,p1024,1024,ISQRT,33600,0.026478325,1268962.444,0
library,unknown,host-cpu,gmp-1t,p1024,1024,MODMUL_R2,12500,0.004276776,2922762.378,0
library,unknown,host-cpu,gmp-nt,p1024,1024,MODMUL_R2,12500,0.001141241,10952988.944,0
library,unknown,host-cpu,openssl-nt,p1024,1024,MODMUL_R2,12500,0.002114356,5911965.627,0
library,NVIDIA GB10,gpu,cgbn,p1024,1024,MODMUL_R2,700000,0.000556064,1258847902.400,0
opencl-kernel,NVIDIA GB10,GPU,w8,p1024,1024,MODMUL_R2,33600,0.006456123,5204361.816,0
opencl-e2e,NVIDIA GB10,GPU,w8,p1024,1024,MODMUL_R2,33600,0.006713433,5004890.939,0
opencl-kernel,NVIDIA GB10,GPU,w16,p1024,1024,MODMUL_R2,33600,0.000367886,91332638.989,0
opencl-e2e,NVIDIA GB10,GPU,w16,p1024,1024,MODMUL_R2,33600,0.000587214,57219343.959,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,p1024,1024,MODMUL_R2,33600,0.000286879,117122548.538,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,p1024,1024,MODMUL_R2,33600,0.000520781,64518482.482,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,p1024,1024,MODMUL_R2,33600,0.000214079,156951414.966,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,p1024,1024,MODMUL_R2,33600,0.000441646,76079031.679,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,p1024,1024,MODMUL_R2,33600,0.000237887,141243532.995,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,p1024,1024,MODMUL_R2,33600,0.000467950,71802543.316,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,p1024,1024,MODMUL_R2,33600,0.000162994,206142564.593,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,p1024,1024,MODMUL_R2,33600,0.000391013,85930646.030,0
library,unknown,host-cpu,gmp-1t,p2048,2048,ADD,6250,0.000109775,56934635.704,0
library,unknown,host-cpu,gmp-nt,p2048,2048,ADD,6250,0.000850219,7351047.174,0
library,unknown,host-cpu,openssl-nt,p2048,2048,ADD,6250,0.001046090,5974629.369,0
library,NVIDIA GB10,gpu,cgbn,p2048,2048,ADD,700000,0.000459456,1523540883.131,0
opencl-kernel,NVIDIA GB10,GPU,w8,p2048,2048,ADD,33600,0.000769868,43643845.534,0
opencl-e2e,NVIDIA GB10,GPU,w8,p2048,2048,ADD,33600,0.001201961,27954317.978,0
opencl-kernel,NVIDIA GB10,GPU,w16,p2048,2048,ADD,33600,0.000343967,97683790.485,0
opencl-e2e,NVIDIA GB10,GPU,w16,p2048,2048,ADD,33600,0.000776972,43244801.840,0
opencl-kernel,NVIDIA GB10,GPU,w32,p2048,2048,ADD,33600,0.000178031,188731178.542,0
opencl-e2e,NVIDIA GB10,GPU,w32,p2048,2048,ADD,33600,0.000613693,54750502.075,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,p2048,2048,ADD,33600,0.000180624,186021787.617,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,p2048,2048,ADD,33600,0.000616109,54535804.660,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,p2048,2048,ADD,33600,0.000190255,176605083.279,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,p2048,2048,ADD,33600,0.000625068,53754151.404,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,p2048,2048,ADD,33600,0.000106560,315315291.244,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,p2048,2048,ADD,33600,0.000550958,60984683.994,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,p2048,2048,ADD,33600,0.000130962,256562975.667,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,p2048,2048,ADD,33600,0.000564711,59499459.908,0
library,unknown,host-cpu,gmp-1t,p2048,2048,SUBTRACT,6250,0.000084032,74376428.917,0
library,unknown,host-cpu,gmp-nt,p2048,2048,SUBTRACT,6250,0.000799036,7821925.394,0
library,unknown,host-cpu,openssl-nt,p2048,2048,SUBTRACT,6250,0.000761196,8210763.033,0
library,NVIDIA GB10,gpu,cgbn,p2048,2048,SUBTRACT,700000,0.000459776,1522480512.249,0
opencl-kernel,NVIDIA GB10,GPU,w8,p2048,2048,SUBTRACT,33600,0.000751292,44722957.353,0
opencl-e2e,NVIDIA GB10,GPU,w8,p2048,2048,SUBTRACT,33600,0.001186585,28316555.251,0
opencl-kernel,NVIDIA GB10,GPU,w16,p2048,2048,SUBTRACT,33600,0.000346511,96966621.121,0
opencl-e2e,NVIDIA GB10,GPU,w16,p2048,2048,SUBTRACT,33600,0.000779740,43091286.670,0
opencl-kernel,NVIDIA GB10,GPU,w32,p2048,2048,SUBTRACT,33600,0.000176959,189874504.500,0
opencl-e2e,NVIDIA GB10,GPU,w32,p2048,2048,SUBTRACT,33600,0.000613485,54769065.499,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,p2048,2048,SUBTRACT,33600,0.000177231,189583079.964,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,p2048,2048,SUBTRACT,33600,0.000614941,54639387.837,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,p2048,2048,SUBTRACT,33600,0.000190783,176116319.698,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,p2048,2048,SUBTRACT,33600,0.000622125,54008439.405,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,p2048,2048,SUBTRACT,33600,0.000111200,302158265.989,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,p2048,2048,SUBTRACT,33600,0.000550478,61037861.825,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,p2048,2048,SUBTRACT,33600,0.000130354,257759658.428,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,p2048,2048,SUBTRACT,33600,0.000560999,59893154.469,0
library,unknown,host-cpu,gmp-1t,p2048,2048,ADDMOD,6250,0.000250527,24947409.913,0
library,unknown,host-cpu,gmp-nt,p2048,2048,ADDMOD,6250,0.000826827,7559017.893,0
library,unknown,host-cpu,openssl-nt,p2048,2048,ADDMOD,6250,0.000992043,6300130.123,0
library,NVIDIA GB10,gpu,cgbn,p2048,2048,ADDMOD,700000,0.000464160,1508100654.947,0
opencl-kernel,NVIDIA GB10,GPU,w8,p2048,2048,ADDMOD,33600,0.000942427,35652628.964,0
opencl-e2e,NVIDIA GB10,GPU,w8,p2048,2048,ADDMOD,33600,0.001363272,24646585.433,0
opencl-kernel,NVIDIA GB10,GPU,w16,p2048,2048,ADDMOD,33600,0.000461934,72737665.579,0
opencl-e2e,NVIDIA GB10,GPU,w16,p2048,2048,ADDMOD,33600,0.000895628,37515575.974,0
opencl-kernel,NVIDIA GB10,GPU,w32,p2048,2048,ADDMOD,33600,0.000206111,163018956.606,0
opencl-e2e,NVIDIA GB10,GPU,w32,p2048,2048,ADDMOD,33600,0.000637357,52717707.097,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,p2048,2048,ADDMOD,33600,0.000253247,132676796.280,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,p2048,2048,ADDMOD,33600,0.000699901,48006789.020,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,p2048,2048,ADDMOD,33600,0.000272319,123384701.534,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,p2048,2048,ADDMOD,33600,0.000727084,46211992.033,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,p2048,2048,ADDMOD,33600,0.000122831,273546581.783,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,p2048,2048,ADDMOD,33600,0.000556718,60353715.332,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,p2048,2048,ADDMOD,33600,0.000155442,216157793.224,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,p2048,2048,ADDMOD,33600,0.000587431,57198207.021,0
library,unknown,host-cpu,gmp-1t,p2048,2048,SUBTRACTMOD,6250,0.000199599,31312781.756,0
library,unknown,host-cpu,gmp-nt,p2048,2048,SUBTRACTMOD,6250,0.000892091,7006011.687,0
library,unknown,host-cpu,openssl-nt,p2048,2048,SUBTRACTMOD,6250,0.000882747,7080171.321,0
library,NVIDIA GB10,gpu,cgbn,p2048,2048,SUBTRACTMOD,700000,0.009610528,72836788.988,0
opencl-kernel,NVIDIA GB10,GPU,w8,p2048,2048,SUBTRACTMOD,33600,0.001016474,33055444.320,0
opencl-e2e,NVIDIA GB10,GPU,w8,p2048,2048,SUBTRACTMOD,33600,0.001436696,23386993.573,0
opencl-kernel,NVIDIA GB10,GPU,w16,p2048,2048,SUBTRACTMOD,33600,0.000499614,67251917.795,0
opencl-e2e,NVIDIA GB10,GPU,w16,p2048,2048,SUBTRACTMOD,33600,0.000927435,36228953.998,0
opencl-kernel,NVIDIA GB10,GPU,w32,p2048,2048,SUBTRACTMOD,33600,0.000239295,140412463.190,0
opencl-e2e,NVIDIA GB10,GPU,w32,p2048,2048,SUBTRACTMOD,33600,0.000675981,49705538.637,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,p2048,2048,SUBTRACTMOD,33600,0.000254719,131910068.576,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,p2048,2048,SUBTRACTMOD,33600,0.000712124,47182794.636,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,p2048,2048,SUBTRACTMOD,33600,0.000285503,117687028.158,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,p2048,2048,SUBTRACTMOD,33600,0.000735148,45705082.409,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,p2048,2048,SUBTRACTMOD,33600,0.000126063,266533381.376,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,p2048,2048,SUBTRACTMOD,33600,0.000559677,60034626.462,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,p2048,2048,SUBTRACTMOD,33600,0.000153522,218861150.950,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,p2048,2048,SUBTRACTMOD,33600,0.000586935,57246543.595,0
library,unknown,host-cpu,gmp-1t,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.002244419,2784685.021,0
library,unknown,host-cpu,gmp-nt,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.001036602,6029315.014,0
library,unknown,host-cpu,openssl-nt,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.001050858,5947520.969,0
opencl-kernel,NVIDIA GB10,GPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,33600,0.151335566,222023.156,0
opencl-e2e,NVIDIA GB10,GPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,33600,0.153142499,219403.498,0
opencl-kernel,NVIDIA GB10,GPU,w16,p2048,2048,MULTIPLYOPERANDSCANNING,33600,0.041065714,818200.799,0
opencl-e2e,NVIDIA GB10,GPU,w16,p2048,2048,MULTIPLYOPERANDSCANNING,33600,0.040877170,821974.711,0
opencl-kernel,NVIDIA GB10,GPU,w32,p2048,2048,MULTIPLYOPERANDSCANNING,33600,0.010618956,3164152.859,0
opencl-e2e,NVIDIA GB10,GPU,w32,p2048,2048,MULTIPLYOPERANDSCANNING,33600,0.011081642,3032041.644,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,p2048,2048,MULTIPLYOPERANDSCANNING,33600,0.000969228,34666765.713,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,p2048,2048,MULTIPLYOPERANDSCANNING,33600,0.001543688,21766056.477,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,p2048,2048,MULTIPLYOPERANDSCANNING,33600,0.000998811,33639998.094,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,p2048,2048,MULTIPLYOPERANDSCANNING,33600,0.001574216,21343957.840,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,p2048,2048,MULTIPLYOPERANDSCANNING,33600,0.000724413,46382381.601,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,p2048,2048,MULTIPLYOPERANDSCANNING,33600,0.001296873,25908473.677,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,p2048,2048,MULTIPLYOPERANDSCANNING,33600,0.000753017,44620506.759,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,p2048,2048,MULTIPLYOPERANDSCANNING,33600,0.001321744,25420958.992,0
library,unknown,host-cpu,gmp-1t,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.002244851,2784149.150,0
library,unknown,host-cpu,gmp-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.001167130,5355016.210,0
library,unknown,host-cpu,openssl-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.000867563,7204087.820,0
library,NVIDIA GB10,gpu,cgbn,p2048,2048,MULTIPLYPRODUCTSCANNING,700000,0.009617024,72787590.007,0
opencl-kernel,NVIDIA GB10,GPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,33600,0.038236885,878732.669,0
opencl-e2e,NVIDIA GB10,GPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,33600,0.038807746,865806.532,0
opencl-kernel,NVIDIA GB10,GPU,w16,p2048,2048,MULTIPLYPRODUCTSCANNING,33600,0.009735027,3451454.220,0
opencl-e2e,NVIDIA GB10,GPU,w16,p2048,2048,MULTIPLYPRODUCTSCANNING,33600,0.010307344,3259811.646,0
opencl-kernel,NVIDIA GB10,GPU,w32,p2048,2048,MULTIPLYPRODUCTSCANNING,33600,0.002468948,13609035.075,0
opencl-e2e,NVIDIA GB10,GPU,w32,p2048,2048,MULTIPLYPRODUCTSCANNING,33600,0.003046545,11028886.782,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,p2048,2048,MULTIPLYPRODUCTSCANNING,33600,0.002502483,13426664.660,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,p2048,2048,MULTIPLYPRODUCTSCANNING,33600,0.003076256,10922367.950,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,p2048,2048,MULTIPLYPRODUCTSCANNING,33600,0.002499364,13443419.948,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,p2048,2048,MULTIPLYPRODUCTSCANNING,33600,0.003073984,10930440.756,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,p2048,2048,MULTIPLYPRODUCTSCANNING,33600,0.001535417,21883305.929,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,p2048,2048,MULTIPLYPRODUCTSCANNING,33600,0.002116181,15877658.858,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,p2048,2048,MULTIPLYPRODUCTSCANNING,33600,0.001464050,22950036.068,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,p2048,2048,MULTIPLYPRODUCTSCANNING,33600,0.002030553,16547216.524,0
library,unknown,host-cpu,gmp-1t,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.014823931,421615.562,0
library,unknown,host-cpu,gmp-nt,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.002358050,2650495.114,0
library,unknown,host-cpu,openssl-nt,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.001227817,5090335.109,0
library,NVIDIA GB10,gpu,cgbn,p2048,2048,MONTGOMERYMULTIPLICATION,700000,0.000977280,716273739.358,0
opencl-kernel,NVIDIA GB10,GPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,33600,0.058025652,579054.243,0
opencl-e2e,NVIDIA GB10,GPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,33600,0.058352338,575812.404,0
opencl-kernel,NVIDIA GB10,GPU,w16,p2048,2048,MONTGOMERYMULTIPLICATION,33600,0.003853982,8718255.573,0
opencl-e2e,NVIDIA GB10,GPU,w16,p2048,2048,MONTGOMERYMULTIPLICATION,33600,0.004326636,7765848.570,0
opencl-kernel,NVIDIA GB10,GPU,w32,p2048,2048,MONTGOMERYMULTIPLICATION,33600,0.000492734,68190950.194,0
opencl-e2e,NVIDIA GB10,GPU,w32,p2048,2048,MONTGOMERYMULTIPLICATION,33600,0.000925307,36312272.496,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,p2048,2048,MONTGOMERYMULTIPLICATION,33600,0.000637789,52682000.265,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,p2048,2048,MONTGOMERYMULTIPLICATION,33600,0.001072907,31316786.851,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,p2048,2048,MONTGOMERYMULTIPLICATION,33600,0.000600141,55986842.397,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,p2048,2048,MONTGOMERYMULTIPLICATION,33600,0.001032619,32538622.827,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,p2048,2048,MONTGOMERYMULTIPLICATION,33600,0.000548862,61217572.908,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,p2048,2048,MONTGOMERYMULTIPLICATION,33600,0.000977291,34380752.580,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,p2048,2048,MONTGOMERYMULTIPLICATION,33600,0.000476134,70568370.134,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,p2048,2048,MONTGOMERYMULTIPLICATION,33600,0.000916459,36662851.216,0
library,unknown,host-cpu,gmp-1t,p2048,2048,COMPARE,6250,0.000034848,179350270.865,0
library,unknown,host-cpu,gmp-nt,p2048,2048,COMPARE,6250,0.001037706,6022900.456,0
library,unknown,host-cpu,openssl-nt,p2048,2048,COMPARE,6250,0.000710636,8794938.631,0
library,NVIDIA GB10,gpu,cgbn,p2048,2048,COMPARE,700000,0.000457248,1530897893.484,0
opencl-kernel,NVIDIA GB10,GPU,w8,p2048,2048,COMPARE,33600,0.000740363,45383143.093,0
opencl-e2e,NVIDIA GB10,GPU,w8,p2048,2048,COMPARE,33600,0.001174506,28607772.028,0
opencl-kernel,NVIDIA GB10,GPU,w16,p2048,2048,COMPARE,33600,0.000321919,104374082.771,0
opencl-e2e,NVIDIA GB10,GPU,w16,p2048,2048,COMPARE,33600,0.000756188,44433394.591,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,p2048,2048,COMPARE,33600,0.000152383,220497037.565,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,p2048,2048,COMPARE,33600,0.000588574,57087129.655,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,p2048,2048,COMPARE,33600,0.000158111,212508937.239,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,p2048,2048,COMPARE,33600,0.000598173,56171041.581,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,p2048,2048,COMPARE,33600,0.000039264,855745573.656,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,p2048,2048,COMPARE,33600,0.000489150,68690585.803,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,p2048,2048,COMPARE,33600,0.000055729,602917646.733,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,p2048,2048,COMPARE,33600,0.000497334,67560231.633,0
library,unknown,host-cpu,gmp-1t,p2048,2048,REDUCE,781,0.000012944,60336897.114,0
library,unknown,host-cpu,gmp-nt,p2048,2048,REDUCE,781,0.000618924,1261867.353,0
library,unknown,host-cpu,openssl-nt,p2048,2048,REDUCE,781,0.000577229,1353015.869,0
library,NVIDIA GB10,gpu,cgbn,p2048,2048,REDUCE,700000,0.000461600,1516464471.404,0
opencl-kernel,NVIDIA GB10,GPU,w8,p2048,2048,REDUCE,33600,0.795573864,42233.665,0
opencl-e2e,NVIDIA GB10,GPU,w8,p2048,2048,REDUCE,33600,0.795670007,42228.562,0
opencl-kernel,NVIDIA GB10,GPU,w16,p2048,2048,REDUCE,33600,0.006720929,4999308.876,0
opencl-e2e,NVIDIA GB10,GPU,w16,p2048,2048,REDUCE,33600,0.007156334,4695141.393,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,p2048,2048,REDUCE,33600,0.001187275,28300098.877,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,p2048,2048,REDUCE,33600,0.001635880,20539403.796,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,p2048,2048,REDUCE,33600,0.001025659,32759425.793,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,p2048,2048,REDUCE,33600,0.001480873,22689318.998,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,p2048,2048,REDUCE,33600,0.001143162,29392159.815,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,p2048,2048,REDUCE,33600,0.001595896,21054003.528,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,p2048,2048,REDUCE,33600,0.000957388,35095489.186,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,p2048,2048,REDUCE,33600,0.001404449,23923973.126,0
library,unknown,host-cpu,gmp-1t,p2048,2048,MODMUL,390,0.000386253,1009700.897,0
library,unknown,host-cpu,gmp-nt,p2048,2048,MODMUL,390,0.000821643,474658.702,0
library,unknown,host-cpu,openssl-nt,p2048,2048,MODMUL,390,0.000918219,424735.277,0
library,NVIDIA GB10,gpu,cgbn,p2048,2048,MODMUL,700000,0.009618752,72774513.783,0
opencl-kernel,NVIDIA GB10,GPU,w8,p2048,2048,MODMUL,33600,1.665851931,20169.860,0
opencl-e2e,NVIDIA GB10,GPU,w8,p2048,2048,MODMUL,33600,1.666248278,20165.062,0
opencl-kernel,NVIDIA GB10,GPU,w16,p2048,2048,MODMUL,33600,0.027845854,1206642.827,0
opencl-e2e,NVIDIA GB10,GPU,w16,p2048,2048,MODMUL,33600,0.028263901,1188795.560,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,p2048,2048,MODMUL,33600,0.006543280,5135039.305,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,p2048,2048,MODMUL,33600,0.006984030,4810975.899,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,p2048,2048,MODMUL,33600,0.012728657,2639712.895,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,p2048,2048,MODMUL,33600,0.013155087,2554145.021,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,p2048,2048,MODMUL,33600,0.006395504,5253690.722,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,p2048,2048,MODMUL,33600,0.006850669,4904630.482,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,p2048,2048,MODMUL,33600,0.012560793,2674990.346,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,p2048,2048,MODMUL,33600,0.013012559,2582120.856,0
library,unknown,host-cpu,gmp-1t,p2048,2048,MODEXP,97,0.134665897,720.301,0
library,unknown,host-cpu,gmp-nt,p2048,2048,MODEXP,97,0.019415312,4996.057,0
library,unknown,host-cpu,openssl-nt,p2048,2048,MODEXP,97,0.016331170,5939.562,0
library,NVIDIA GB10,gpu,cgbn,p2048,2048,MODEXP,700000,6.987298489,100181.780,0
opencl-kernel,NVIDIA GB10,GPU,w8,p2048,2048,MODEXP,33600,74.936678228,448.379,0
opencl-e2e,NVIDIA GB10,GPU,w8,p2048,2048,MODEXP,33600,75.019208644,447.885,0
opencl-kernel,NVIDIA GB10,GPU,w16,p2048,2048,MODEXP,33600,20.216655012,1661.996,0
opencl-e2e,NVIDIA GB10,GPU,w16,p2048,2048,MODEXP,33600,20.039069863,1676.725,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,p2048,2048,MODEXP,33600,1.626154148,20662.248,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,p2048,2048,MODEXP,33600,1.626999747,20651.509,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,p2048,2048,MODEXP,33600,4.096316203,8202.492,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,p2048,2048,MODEXP,33600,4.118985821,8157.348,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,p2048,2048,MODEXP,33600,1.624559274,20682.533,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,p2048,2048,MODEXP,33600,1.625393028,20671.923,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,p2048,2048,MODEXP,33600,4.098514503,8198.092,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,p2048,2048,MODEXP,33600,4.119858283,8155.620,0
library,unknown,host-cpu,gmp-1t,p2048,2048,EXPONENTIATION,97,0.015230984,6368.597,0
library,unknown,host-cpu,gmp-nt,p2048,2048,EXPONENTIATION,97,0.004347431,22312.027,0
library,unknown,host-cpu,openssl-nt,p2048,2048,EXPONENTIATION,97,0.028399053,3415.607,0
opencl-kernel,NVIDIA GB10,GPU,w8,p2048,2048,EXPONENTIATION,33600,28.625600252,1173.775,0
opencl-e2e,NVIDIA GB10,GPU,w8,p2048,2048,EXPONENTIATION,33600,28.621325246,1173.950,0
opencl-kernel,NVIDIA GB10,GPU,w16,p2048,2048,EXPONENTIATION,33600,6.478926717,5186.044,0
opencl-e2e,NVIDIA GB10,GPU,w16,p2048,2048,EXPONENTIATION,33600,6.454730663,5205.484,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,p2048,2048,EXPONENTIATION,33600,1.627526245,20644.828,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,p2048,2048,EXPONENTIATION,33600,1.622246063,20712.024,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,p2048,2048,EXPONENTIATION,33600,1.567091282,21440.997,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,p2048,2048,EXPONENTIATION,33600,1.570137426,21399.401,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,p2048,2048,EXPONENTIATION,33600,1.600969211,20987.287,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,p2048,2048,EXPONENTIATION,33600,1.602815806,20963.107,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,p2048,2048,EXPONENTIATION,33600,1.611731560,20847.144,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,p2048,2048,EXPONENTIATION,33600,1.610514846,20862.894,0
library,unknown,host-cpu,gmp-1t,p2048,2048,DIVIDE,781,0.000021439,36428934.498,0
library,unknown,host-cpu,gmp-nt,p2048,2048,DIVIDE,781,0.000782684,997848.428,0
library,unknown,host-cpu,openssl-nt,p2048,2048,DIVIDE,781,0.000715932,1090885.714,0
library,NVIDIA GB10,gpu,cgbn,p2048,2048,DIVIDE,700000,0.000884864,791082019.384,0
opencl-kernel,NVIDIA GB10,GPU,w8,p2048,2048,DIVIDE,33600,2.554980799,13150.784,0
opencl-e2e,NVIDIA GB10,GPU,w8,p2048,2048,DIVIDE,33600,2.556112233,13144.963,0
opencl-kernel,NVIDIA GB10,GPU,w16,p2048,2048,DIVIDE,33600,0.217723768,154323.987,0
opencl-e2e,NVIDIA GB10,GPU,w16,p2048,2048,DIVIDE,33600,0.217574633,154429.768,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,p2048,2048,DIVIDE,33600,0.020308220,1654502.462,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,p2048,2048,DIVIDE,33600,0.020979144,1601590.608,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,p2048,2048,DIVIDE,33600,0.018940210,1774003.562,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,p2048,2048,DIVIDE,33600,0.019369295,1734704.335,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,p2048,2048,DIVIDE,33600,0.019868608,1691109.916,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,p2048,2048,DIVIDE,33600,0.020401767,1646916.171,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,p2048,2048,DIVIDE,33600,0.018218735,1844255.378,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,p2048,2048,DIVIDE,33600,0.018728546,1794052.779,0
library,unknown,host-cpu,gmp-1t,p2048,2048,ISQRT,195,0.000082111,2374834.047,0
library,unknown,host-cpu,gmp-nt,p2048,2048,ISQRT,195,0.000757883,257295.647,0
opencl-kernel,NVIDIA GB10,GPU,w8,p2048,2048,ISQRT,33600,36.127179921,930.048,0
opencl-e2e,NVIDIA GB10,GPU,w8,p2048,2048,ISQRT,33600,36.111542230,930.450,0
opencl-kernel,NVIDIA GB10,GPU,w16,p2048,2048,ISQRT,33600,4.050270313,8295.743,0
opencl-e2e,NVIDIA GB10,GPU,w16,p2048,2048,ISQRT,33600,4.055599958,8284.841,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,p2048,2048,ISQRT,33600,0.175973786,190937.530,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,p2048,2048,ISQRT,33600,0.176002586,190906.286,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,p2048,2048,ISQRT,33600,0.075537384,444812.863,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,p2048,2048,ISQRT,33600,0.076273220,440521.588,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,p2048,2048,ISQRT,33600,0.177067408,189758.242,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,p2048,2048,ISQRT,33600,0.178502039,188233.144,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,p2048,2048,ISQRT,33600,0.076732184,437886.663,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,p2048,2048,ISQRT,33600,0.076546230,438950.422,0
library,unknown,host-cpu,gmp-1t,p2048,2048,MODMUL_R2,6250,0.006469066,966136.378,0
library,unknown,host-cpu,gmp-nt,p2048,2048,MODMUL_R2,6250,0.001520775,4109746.684,0
library,unknown,host-cpu,openssl-nt,p2048,2048,MODMUL_R2,6250,0.003045758,2052034.338,0
library,NVIDIA GB10,gpu,cgbn,p2048,2048,MODMUL_R2,700000,0.001894976,369397818.231,0
opencl-kernel,NVIDIA GB10,GPU,w8,p2048,2048,MODMUL_R2,33600,0.041933231,801273.816,0
opencl-e2e,NVIDIA GB10,GPU,w8,p2048,2048,MODMUL_R2,33600,0.042204109,796131.012,0
opencl-kernel,NVIDIA GB10,GPU,w16,p2048,2048,MODMUL_R2,33600,0.003549070,9467268.902,0
opencl-e2e,NVIDIA GB10,GPU,w16,p2048,2048,MODMUL_R2,33600,0.003976413,8449826.504,0
opencl-kernel,NVIDIA GB10,GPU,w32-opt,p2048,2048,MODMUL_R2,33600,0.001126394,29829704.456,0
opencl-e2e,NVIDIA GB10,GPU,w32-opt,p2048,2048,MODMUL_R2,33600,0.001552824,21637996.269,0
opencl-kernel,NVIDIA GB10,GPU,w32-o64,p2048,2048,MODMUL_R2,33600,0.000996139,33730232.347,0
opencl-e2e,NVIDIA GB10,GPU,w32-o64,p2048,2048,MODMUL_R2,33600,0.001431577,23470620.245,0
opencl-kernel,NVIDIA GB10,GPU,w32-il,p2048,2048,MODMUL_R2,33600,0.001020349,32929909.449,0
opencl-e2e,NVIDIA GB10,GPU,w32-il,p2048,2048,MODMUL_R2,33600,0.001451603,23146824.484,0
opencl-kernel,NVIDIA GB10,GPU,w32-il64,p2048,2048,MODMUL_R2,33600,0.000871703,38545238.477,0
opencl-e2e,NVIDIA GB10,GPU,w32-il64,p2048,2048,MODMUL_R2,33600,0.001312395,25602048.070,0
```
