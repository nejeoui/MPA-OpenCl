# MPA-OpenCL benchmark report - NVIDIA GeForce GTX 1060 3GB


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

### Device 0 - NVIDIA GeForce GTX 1060 3GB (GPU)

| Property | Value |
|---|---|
| Model | NVIDIA GeForce GTX 1060 3GB |
| Type | GPU |
| Vendor | NVIDIA Corporation |
| Device memory | 2.94 GiB |
| Max single allocation | 0.73 GiB |
| Local memory | 48 KiB |
| Global cache | 432 KiB |
| Compute units | 9 |
| Max clock | 1708 MHz |
| Max work-group size | 1024 |
| OpenCL version | OpenCL 3.0 CUDA |
| Driver | 570.86.16 |

### Device 1 - cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz (CPU)

| Property | Value |
|---|---|
| Model | cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz |
| Type | CPU |
| Vendor | GenuineIntel |
| Device memory | 60.71 GiB |
| Max single allocation | 16.00 GiB |
| Local memory | 256 KiB |
| Global cache | 40960 KiB |
| Compute units | 32 |
| Max clock | 2600 MHz |
| Max work-group size | 4096 |
| OpenCL version | OpenCL 3.0 PoCL HSTR: cpu-x86_64-pc-linux-gnu-haswell |
| Driver | 5.0+debian |

### Host

| Property | Value |
|---|---|
| CPU | Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz |
| Logical cores | 32 |
| OpenMP threads used | 32 |
| RAM | 62.7 GB |
| OS | Ubuntu 24.04.4 LTS |
| Kernel | 6.8.0-110-generic |
| Arch | x86_64 |
| GMP | 6.3.0 |
| OpenSSL | OpenSSL 3.0.13 30 Jan 2024 |

## 2. Method

- Workload auto-sized from the device and host: --min-items from 700 x compute units, --items from ten times that capped by host RAM. Either flag, given explicitly, overrides its half.
- Base workload 50000 items, scaled down per operator by its cost weight and by modulus size. Device rows honour --min-items (6300) so the GPU is not left idle; the CPU libraries keep the smaller count because a full-width MODEXP there costs minutes. Both counts appear in every row as dev/cpu, and throughput is per-second so they remain comparable.
- 5 timed repetitions, **minimum** reported. Two untimed warm-up launches precede them.
- `kernel` times `clEnqueueNDRangeKernel` + `clFinish` only. `e2e` adds the host->device operand writes and the device->host result read.
- Every OpenCL device runs the same kernels on the same operands, so GPU and CPU-OpenCL columns are directly comparable.
- CPU library baselines (GMP, OpenSSL) run those same operands, with temporaries preallocated outside the timed region, so the figure is the arithmetic and not marshalling. The generator is reseeded per modulus and operation so every backend sees identical inputs.
- OpenSSL rows time the nearest BN primitive, which is not always semantically identical (its Montgomery routine expects Montgomery-domain inputs); they measure comparable work, not identical results. Correctness is judged against GMP only.
- Every device cell is checked word-for-word against GMP before it is timed. A cell that mismatches is reported and excluded from the speedup tables.
- Total wall time 5457.1 s.

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
| [1] CPU | `mpaKernels_8bits.cl` (w8) | 72 | 72 | 0 | 0 |

**All configurations correct** - 557 configurations, 0 problems.

## 4. Throughput per device

Operations per second, higher is better. Kernel-only timings.

### Device 0 - NVIDIA GeForce GTX 1060 3GB (GPU)

#### secp256k1 (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | 190.07 M | 442.84 M | 719.88 M | 709.04 M | 715.33 M | 1.08 G | 1.17 G | 34.55 M |
| SUBTRACT | 50000 / 50000 | 189.83 M | 440.72 M | 716.82 M | 709.68 M | 694.51 M | 1.09 G | 1.19 G | 47.10 M |
| ADDMOD | 50000 / 50000 | 137.20 M | 316.95 M | 548.93 M | 701.85 M | 674.24 M | 1.26 G | 1.26 G | 11.66 M |
| SUBTRACTMOD | 50000 / 50000 | 137.10 M | 315.02 M | 548.78 M | 687.66 M | 668.55 M | 1.27 G | 1.28 G | 16.47 M |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | 5.49 M | 22.77 M | 79.50 M | 419.50 M | 410.33 M | 955.22 M | 897.96 M | 31.64 M |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | 9.90 M | 46.19 M | 172.94 M | 175.26 M | 173.07 M | 272.01 M | 316.41 M | 31.59 M |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | 51.30 M | 206.82 M | 629.67 M | 575.60 M | 602.90 M | 754.39 M | 771.88 M | 4.07 M |
| COMPARE | 50000 / 50000 | 286.03 M | 559.35 M | - | 945.26 M | 947.83 M | 2.46 G | 2.51 G | 92.14 M |
| REDUCE | 6300 / 6250 | 32.06 M | 59.98 M | - | 127.84 M | 125.71 M | 133.04 M | 138.48 M | 33.17 M |
| MODMUL | 6300 / 3125 | 10.94 M | 23.60 M | - | 52.06 M | 58.40 M | 53.26 M | 60.46 M | 7.08 M |
| MODEXP | 6300 / 781 | 260.20 k | 1.07 M | - | 1.60 M | 1.99 M | 1.60 M | 1.99 M | 68.28 k |
| EXPONENTIATION | 6300 / 781 | 52.77 k | 286.71 k | - | 9.64 M | 10.14 M | 9.63 M | 10.23 M | 193.19 k |
| DIVIDE | 6300 / 6250 | 17.73 M | 28.03 M | - | 56.14 M | 53.50 M | 58.48 M | 60.09 M | 14.35 M |
| ISQRT | 6300 / 1562 | 1.18 M | 2.34 M | - | 8.03 M | 9.57 M | 8.86 M | 9.51 M | 7.14 M |
| MODMUL_R2 | 50000 / 50000 | 53.05 M | 200.77 M | - | 183.48 M | 204.87 M | 199.82 M | 241.28 M | 7.05 M |

#### rsa256(composite) (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | 228.75 M | 441.71 M | 719.84 M | 717.17 M | 712.47 M | 1.08 G | 1.17 G | 38.96 M |
| SUBTRACT | 50000 / 50000 | 227.20 M | 441.94 M | 716.59 M | 721.85 M | 710.99 M | 1.09 G | 1.19 G | 47.63 M |
| ADDMOD | 50000 / 50000 | 176.51 M | 341.17 M | 586.40 M | 679.94 M | 681.29 M | 1.26 G | 1.27 G | 15.01 M |
| SUBTRACTMOD | 50000 / 50000 | 164.63 M | 317.38 M | 557.64 M | 692.67 M | 677.70 M | 1.27 G | 1.28 G | 16.45 M |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | 6.09 M | 22.71 M | 79.46 M | 426.76 M | 415.81 M | 936.42 M | 892.79 M | 31.55 M |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | 10.92 M | 46.26 M | 172.80 M | 175.34 M | 173.27 M | 272.48 M | 317.95 M | 31.63 M |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | 54.31 M | 209.86 M | 625.48 M | 587.49 M | 595.54 M | 759.49 M | 767.81 M | 4.07 M |
| COMPARE | 50000 / 50000 | 300.11 M | 567.54 M | - | 954.51 M | 954.13 M | 2.40 G | 2.46 G | 92.24 M |
| REDUCE | 6300 / 6250 | 33.41 M | 59.28 M | - | 127.12 M | 124.51 M | 132.24 M | 136.68 M | 20.32 M |
| MODMUL | 6300 / 3125 | 11.51 M | 23.70 M | - | 52.51 M | 59.02 M | 52.54 M | 61.59 M | 7.10 M |
| MODEXP | 6300 / 781 | 259.40 k | 1.07 M | - | 1.65 M | 1.99 M | 1.63 M | 2.01 M | 71.43 k |
| EXPONENTIATION | 6300 / 781 | 52.25 k | 286.75 k | - | 9.63 M | 10.15 M | 9.62 M | 10.26 M | 193.97 k |
| DIVIDE | 6300 / 6250 | 17.65 M | 27.82 M | - | 56.35 M | 54.83 M | 57.94 M | 59.74 M | 13.60 M |
| ISQRT | 6300 / 1562 | 1.18 M | 2.31 M | - | 7.68 M | 9.03 M | 8.39 M | 8.98 M | 6.95 M |
| MODMUL_R2 | 50000 / 50000 | 52.67 M | 202.32 M | - | 182.52 M | 202.39 M | 199.06 M | 246.49 M | 7.03 M |

#### brainpoolP512r1 (512-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 25000 / 25000 | 102.57 M | 189.54 M | 328.01 M | 345.08 M | 342.12 M | 324.03 M | 479.85 M | 36.36 M |
| SUBTRACT | 25000 / 25000 | 102.75 M | 192.40 M | 329.48 M | 341.35 M | 340.93 M | 323.00 M | 478.61 M | 43.01 M |
| ADDMOD | 25000 / 25000 | 75.74 M | 154.08 M | 270.50 M | 260.41 M | 261.67 M | 502.86 M | 575.25 M | 13.47 M |
| SUBTRACTMOD | 25000 / 25000 | 68.63 M | 141.92 M | 247.59 M | 266.23 M | 262.13 M | 489.51 M | 553.73 M | 15.01 M |
| MULTIPLYOPERANDSCANNING | 25000 / 25000 | 1.12 M | 4.37 M | 5.33 M | 116.99 M | 117.12 M | 170.00 M | 216.17 M | 13.90 M |
| MULTIPLYPRODUCTSCANNING | 25000 / 25000 | 3.10 M | 12.19 M | 42.20 M | 45.62 M | 44.89 M | 34.43 M | 62.57 M | 13.94 M |
| MONTGOMERYMULTIPLICATION | 25000 / 25000 | 10.53 M | 48.70 M | 201.82 M | 155.12 M | 60.64 M | 154.55 M | 77.78 M | 1.82 M |
| COMPARE | 25000 / 25000 | 140.89 M | 263.28 M | - | 461.44 M | 463.59 M | 1.16 G | 1.32 G | 85.62 M |
| REDUCE | 6300 / 3125 | 10.31 M | 17.96 M | - | 43.27 M | 41.35 M | 35.31 M | 40.09 M | 19.19 M |
| MODMUL | 6300 / 1562 | 3.17 M | 6.12 M | - | 13.34 M | 11.34 M | 10.86 M | 12.21 M | 3.64 M |
| MODEXP | 6300 / 390 | 14.71 k | 154.83 k | - | 168.71 k | 110.23 k | 165.35 k | 111.67 k | 13.13 k |
| EXPONENTIATION | 6300 / 390 | 6.46 k | 26.04 k | - | 260.12 k | 268.64 k | 236.23 k | 267.71 k | 53.43 k |
| DIVIDE | 6300 / 3125 | 4.01 M | 3.16 M | - | 16.24 M | 18.24 M | 20.46 M | 18.75 M | 10.94 M |
| ISQRT | 6300 / 781 | 164.97 k | 273.38 k | - | 1.38 M | 1.58 M | 1.25 M | 1.71 M | 3.64 M |
| MODMUL_R2 | 25000 / 25000 | 10.76 M | 54.48 M | - | 89.93 M | 31.84 M | 81.18 M | 33.52 M | 3.61 M |

#### p1024 (1024-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 12500 / 12500 | 34.13 M | 66.53 M | 96.03 M | 123.36 M | 124.62 M | 163.01 M | 162.71 M | 29.24 M |
| SUBTRACT | 12500 / 12500 | 34.29 M | 66.58 M | 96.33 M | 123.67 M | 124.76 M | 161.74 M | 162.08 M | 36.32 M |
| ADDMOD | 12500 / 12500 | 24.12 M | 47.41 M | 79.01 M | 98.55 M | 97.85 M | 250.96 M | 266.79 M | 9.24 M |
| SUBTRACTMOD | 12500 / 12500 | 24.05 M | 47.45 M | 78.19 M | 98.64 M | 100.29 M | 247.03 M | 259.91 M | 12.99 M |
| MULTIPLYOPERANDSCANNING | 12500 / 12500 | 225.78 k | 897.14 k | 981.68 k | 40.13 M | 41.19 M | 55.33 M | 60.58 M | 4.15 M |
| MULTIPLYPRODUCTSCANNING | 12500 / 12500 | 545.30 k | 2.15 M | 8.07 M | 7.37 M | 7.48 M | 8.74 M | 9.06 M | 4.15 M |
| MONTGOMERYMULTIPLICATION | 12500 / 12500 | 1.76 M | 9.94 M | 49.71 M | 32.21 M | 36.08 M | 40.81 M | 47.23 M | 621.64 k |
| COMPARE | 12500 / 12500 | 52.46 M | 101.72 M | - | 187.28 M | 186.89 M | 547.69 M | 540.59 M | 75.76 M |
| REDUCE | 6300 / 1562 | 1.64 M | 4.98 M | - | 14.29 M | 13.44 M | 14.66 M | 15.10 M | 26.01 M |
| MODMUL | 6300 / 781 | 453.70 k | 1.40 M | - | 3.04 M | 3.19 M | 3.03 M | 3.23 M | 1.40 M |
| MODEXP | 6300 / 195 | 1.75 k | 12.88 k | - | 20.76 k | 22.56 k | 20.34 k | 23.52 k | 2.05 k |
| EXPONENTIATION | 6300 / 195 | 789.7 | 3.25 k | - | 27.65 k | 28.45 k | 26.65 k | 29.80 k | 14.01 k |
| DIVIDE | 6300 / 1562 | 81.75 k | 257.19 k | - | 2.32 M | 5.41 M | 2.45 M | 5.47 M | 8.01 M |
| ISQRT | 6300 / 390 | 6.12 k | 25.79 k | - | 225.84 k | 283.08 k | 284.83 k | 296.61 k | 1.40 M |
| MODMUL_R2 | 12500 / 12500 | 1.66 M | 13.31 M | - | 18.75 M | 21.76 M | 20.83 M | 24.81 M | 1.37 M |

#### p2048 (2048-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 6300 / 6250 | 15.61 M | 30.83 M | 58.92 M | 58.77 M | 58.15 M | 82.32 M | 82.41 M | 21.41 M |
| SUBTRACT | 6300 / 6250 | 15.64 M | 31.02 M | 58.90 M | 58.84 M | 58.03 M | 82.42 M | 82.18 M | 24.23 M |
| ADDMOD | 6300 / 6250 | 11.58 M | 23.35 M | 44.64 M | 46.81 M | 47.38 M | 89.62 M | 114.87 M | 7.86 M |
| SUBTRACTMOD | 6300 / 6250 | 11.31 M | 22.76 M | 44.14 M | 48.50 M | 47.61 M | 117.80 M | 111.63 M | 9.89 M |
| MULTIPLYOPERANDSCANNING | 6300 / 6250 | 51.04 k | 188.28 k | 660.90 k | 8.89 M | 11.52 M | 10.00 M | 12.63 M | 1.26 M |
| MULTIPLYPRODUCTSCANNING | 6300 / 6250 | 133.31 k | 532.01 k | 2.06 M | 1.87 M | 1.88 M | 1.83 M | 1.96 M | 1.26 M |
| MONTGOMERYMULTIPLICATION | 6300 / 6250 | 77.64 k | 1.96 M | 8.95 M | 8.28 M | 9.34 M | 8.45 M | 10.02 M | 195.51 k |
| COMPARE | 6300 / 6250 | 25.37 M | 50.45 M | - | 91.65 M | 92.46 M | 185.15 M | 189.53 M | 69.65 M |
| REDUCE | 6300 / 781 | 16.15 k | 841.90 k | - | 4.29 M | 4.51 M | 4.88 M | 4.62 M | 17.50 M |
| MODMUL | 6300 / 390 | 6.47 k | 240.36 k | - | 575.45 k | 882.52 k | 528.89 k | 907.05 k | 472.99 k |
| MODEXP | 6300 / 97 | 58.8 | 383.0 | - | 1.19 k | 2.05 k | 1.21 k | 2.09 k | 292.6 |
| EXPONENTIATION | 6300 / 97 | 91.8 | 395.0 | - | 3.29 k | 3.79 k | 3.32 k | 3.73 k | 2.45 k |
| DIVIDE | 6300 / 781 | 5.67 k | 16.35 k | - | 128.46 k | 1.31 M | 128.34 k | 1.32 M | 5.72 M |
| ISQRT | 6300 / 195 | 388.8 | 680.7 | - | 19.80 k | 71.46 k | 19.80 k | 71.62 k | 774.36 k |
| MODMUL_R2 | 6300 / 6250 | 145.94 k | 1.04 M | - | 3.68 M | 4.68 M | 3.32 M | 4.54 M | 457.81 k |

### Device 1 - cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz (CPU)

#### secp256k1 (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | 97.88 M | - | - | - | - | - | - | 34.55 M |
| SUBTRACT | 50000 / 50000 | 91.83 M | - | - | - | - | - | - | 47.10 M |
| ADDMOD | 50000 / 50000 | 79.91 M | - | - | - | - | - | - | 11.66 M |
| SUBTRACTMOD | 50000 / 50000 | 85.56 M | - | - | - | - | - | - | 16.47 M |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | 4.54 M | - | - | - | - | - | - | 31.64 M |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | 5.30 M | - | - | - | - | - | - | 31.59 M |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | 1.65 M | - | - | - | - | - | - | 4.07 M |
| COMPARE | 50000 / 50000 | 172.85 M | - | - | - | - | - | - | 92.14 M |
| REDUCE | 6300 / 6250 | 1.21 M | - | - | - | - | - | - | 33.17 M |
| MODMUL | 6300 / 3125 | 344.60 k | - | - | - | - | - | - | 7.08 M |
| MODEXP | 6300 / 781 | 9.93 k | - | - | - | - | - | - | 68.28 k |
| EXPONENTIATION | 6300 / 781 | 23.50 k | - | - | - | - | - | - | 193.19 k |
| DIVIDE | 6300 / 6250 | 738.52 k | - | - | - | - | - | - | 14.35 M |
| ISQRT | 6300 / 1562 | 108.74 k | - | - | - | - | - | - | 7.14 M |
| MODMUL_R2 | 50000 / 50000 | 1.40 M | - | - | - | - | - | - | 7.05 M |

#### rsa256(composite) (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | 126.69 M | - | - | - | - | - | - | 38.96 M |
| SUBTRACT | 50000 / 50000 | 82.99 M | - | - | - | - | - | - | 47.63 M |
| ADDMOD | 50000 / 50000 | 75.63 M | - | - | - | - | - | - | 15.01 M |
| SUBTRACTMOD | 50000 / 50000 | 75.22 M | - | - | - | - | - | - | 16.45 M |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | 5.74 M | - | - | - | - | - | - | 31.55 M |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | 5.32 M | - | - | - | - | - | - | 31.63 M |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | 1.37 M | - | - | - | - | - | - | 4.07 M |
| COMPARE | 50000 / 50000 | 183.43 M | - | - | - | - | - | - | 92.24 M |
| REDUCE | 6300 / 6250 | 1.20 M | - | - | - | - | - | - | 20.32 M |
| MODMUL | 6300 / 3125 | 343.86 k | - | - | - | - | - | - | 7.10 M |
| MODEXP | 6300 / 781 | 9.93 k | - | - | - | - | - | - | 71.43 k |
| EXPONENTIATION | 6300 / 781 | 24.66 k | - | - | - | - | - | - | 193.97 k |
| DIVIDE | 6300 / 6250 | 749.66 k | - | - | - | - | - | - | 13.60 M |
| ISQRT | 6300 / 1562 | 88.88 k | - | - | - | - | - | - | 6.95 M |
| MODMUL_R2 | 50000 / 50000 | 1.39 M | - | - | - | - | - | - | 7.03 M |

#### brainpoolP512r1 (512-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 25000 / 25000 | 53.49 M | - | - | - | - | - | - | 36.36 M |
| SUBTRACT | 25000 / 25000 | 43.41 M | - | - | - | - | - | - | 43.01 M |
| ADDMOD | 25000 / 25000 | 44.98 M | - | - | - | - | - | - | 13.47 M |
| SUBTRACTMOD | 25000 / 25000 | 45.40 M | - | - | - | - | - | - | 15.01 M |
| MULTIPLYOPERANDSCANNING | 25000 / 25000 | 1.11 M | - | - | - | - | - | - | 13.90 M |
| MULTIPLYPRODUCTSCANNING | 25000 / 25000 | 1.50 M | - | - | - | - | - | - | 13.94 M |
| MONTGOMERYMULTIPLICATION | 25000 / 25000 | 401.97 k | - | - | - | - | - | - | 1.82 M |
| COMPARE | 25000 / 25000 | 108.86 M | - | - | - | - | - | - | 85.62 M |
| REDUCE | 6300 / 3125 | 340.09 k | - | - | - | - | - | - | 19.19 M |
| MODMUL | 6300 / 1562 | 123.96 k | - | - | - | - | - | - | 3.64 M |
| MODEXP | 6300 / 390 | 1.19 k | - | - | - | - | - | - | 13.13 k |
| EXPONENTIATION | 6300 / 390 | 2.71 k | - | - | - | - | - | - | 53.43 k |
| DIVIDE | 6300 / 3125 | 214.52 k | - | - | - | - | - | - | 10.94 M |
| ISQRT | 6300 / 781 | 27.20 k | - | - | - | - | - | - | 3.64 M |
| MODMUL_R2 | 25000 / 25000 | 501.46 k | - | - | - | - | - | - | 3.61 M |

#### p1024 (1024-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 12500 / 12500 | 29.85 M | - | - | - | - | - | - | 29.24 M |
| SUBTRACT | 12500 / 12500 | 27.46 M | - | - | - | - | - | - | 36.32 M |
| ADDMOD | 12500 / 12500 | 20.43 M | - | - | - | - | - | - | 9.24 M |
| SUBTRACTMOD | 12500 / 12500 | 22.10 M | - | - | - | - | - | - | 12.99 M |
| MULTIPLYOPERANDSCANNING | 12500 / 12500 | 321.42 k | - | - | - | - | - | - | 4.15 M |
| MULTIPLYPRODUCTSCANNING | 12500 / 12500 | 475.77 k | - | - | - | - | - | - | 4.15 M |
| MONTGOMERYMULTIPLICATION | 12500 / 12500 | 147.81 k | - | - | - | - | - | - | 621.64 k |
| COMPARE | 12500 / 12500 | 79.43 M | - | - | - | - | - | - | 75.76 M |
| REDUCE | 6300 / 1562 | 109.76 k | - | - | - | - | - | - | 26.01 M |
| MODMUL | 6300 / 781 | 42.10 k | - | - | - | - | - | - | 1.40 M |
| MODEXP | 6300 / 195 | 155.1 | - | - | - | - | - | - | 2.05 k |
| EXPONENTIATION | 6300 / 195 | 376.2 | - | - | - | - | - | - | 14.01 k |
| DIVIDE | 6300 / 1562 | 80.23 k | - | - | - | - | - | - | 8.01 M |
| ISQRT | 6300 / 390 | 6.40 k | - | - | - | - | - | - | 1.40 M |
| MODMUL_R2 | 12500 / 12500 | 155.32 k | - | - | - | - | - | - | 1.37 M |

#### p2048 (2048-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 6300 / 6250 | 14.30 M | - | - | - | - | - | - | 21.41 M |
| SUBTRACT | 6300 / 6250 | 11.93 M | - | - | - | - | - | - | 24.23 M |
| ADDMOD | 6300 / 6250 | 10.19 M | - | - | - | - | - | - | 7.86 M |
| SUBTRACTMOD | 6300 / 6250 | 9.99 M | - | - | - | - | - | - | 9.89 M |
| MULTIPLYOPERANDSCANNING | 6300 / 6250 | 94.09 k | - | - | - | - | - | - | 1.26 M |
| MULTIPLYPRODUCTSCANNING | 6300 / 6250 | 123.47 k | - | - | - | - | - | - | 1.26 M |
| MONTGOMERYMULTIPLICATION | 6300 / 6250 | 33.47 k | - | - | - | - | - | - | 195.51 k |
| COMPARE | 6300 / 6250 | 31.88 M | - | - | - | - | - | - | 69.65 M |
| REDUCE | 6300 / 781 | 48.57 k | - | - | - | - | - | - | 17.50 M |
| MODMUL | 6300 / 390 | 11.57 k | - | - | - | - | - | - | 472.99 k |
| MODEXP | 6300 / 97 | over budget | - | - | - | - | - | - | 292.6 |
| EXPONENTIATION | 6300 / 97 | over budget | - | - | - | - | - | - | 2.45 k |
| DIVIDE | 6300 / 781 | - | - | - | - | - | - | - | 5.72 M |
| ISQRT | 6300 / 195 | - | - | - | - | - | - | - | 774.36 k |
| MODMUL_R2 | 6300 / 6250 | - | - | - | - | - | - | - | 457.81 k |

## 5. Head to head

Best OpenCL GPU result against best OpenCL CPU result and the CPU libraries.
Ratios above 1.00x mean the GPU is faster than that baseline.

### secp256k1 (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w32-il64 | 1.17 G | w8 | 97.88 M | 34.55 M | 11.91x |
| SUBTRACT | w32-il64 | 1.19 G | w8 | 91.83 M | 47.10 M | 12.95x |
| ADDMOD | w32-il | 1.26 G | w8 | 79.91 M | 11.66 M | 15.80x |
| SUBTRACTMOD | w32-il64 | 1.28 G | w8 | 85.56 M | 16.47 M | 14.99x |
| MULTIPLYOPERANDSCANNING | w32-il | 955.22 M | w8 | 4.54 M | 31.64 M | 210.44x |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 316.41 M | w8 | 5.30 M | 31.59 M | 59.72x |
| MONTGOMERYMULTIPLICATION | w32-il64 | 771.88 M | w8 | 1.65 M | 4.07 M | 467.68x |
| COMPARE | w32-il64 | 2.51 G | w8 | 172.85 M | 92.14 M | 14.52x |
| REDUCE | w32-il64 | 137.38 M | w8 | 1.20 M | 33.17 M | 114.75x |
| MODMUL | w32-il64 | 29.99 M | w8 | 170.93 k | 7.08 M | 175.44x |
| MODEXP | w32-o64 | 246.63 k | w8 | 1.23 k | 68.28 k | 200.41x |
| EXPONENTIATION | w32-il64 | 1.27 M | w8 | 2.91 k | 193.19 k | 435.18x |
| DIVIDE | w32-il64 | 59.61 M | w8 | 732.66 k | 14.35 M | 81.36x |
| ISQRT | w32-o64 | 2.37 M | w8 | 26.96 k | 7.14 M | 88.05x |
| MODMUL_R2 | w32-il64 | 241.28 M | w8 | 1.40 M | 7.05 M | 172.79x |

### rsa256(composite) (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w32-il64 | 1.17 G | w8 | 126.69 M | 38.96 M | 9.27x |
| SUBTRACT | w32-il64 | 1.19 G | w8 | 82.99 M | 47.63 M | 14.30x |
| ADDMOD | w32-il64 | 1.27 G | w8 | 75.63 M | 15.01 M | 16.84x |
| SUBTRACTMOD | w32-il64 | 1.28 G | w8 | 75.22 M | 16.45 M | 17.02x |
| MULTIPLYOPERANDSCANNING | w32-il | 936.42 M | w8 | 5.74 M | 31.55 M | 163.13x |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 317.95 M | w8 | 5.32 M | 31.63 M | 59.82x |
| MONTGOMERYMULTIPLICATION | w32-il64 | 767.81 M | w8 | 1.37 M | 4.07 M | 561.76x |
| COMPARE | w32-il64 | 2.46 G | w8 | 183.43 M | 92.24 M | 13.42x |
| REDUCE | w32-il64 | 135.60 M | w8 | 1.19 M | 20.32 M | 114.01x |
| MODMUL | w32-il64 | 30.55 M | w8 | 170.57 k | 7.10 M | 179.12x |
| MODEXP | w32-il64 | 248.80 k | w8 | 1.23 k | 71.43 k | 202.12x |
| EXPONENTIATION | w32-il64 | 1.27 M | w8 | 3.06 k | 193.97 k | 415.84x |
| DIVIDE | w32-il64 | 59.27 M | w8 | 743.71 k | 13.60 M | 79.69x |
| ISQRT | w32-o64 | 2.24 M | w8 | 22.04 k | 6.95 M | 101.55x |
| MODMUL_R2 | w32-il64 | 246.49 M | w8 | 1.39 M | 7.03 M | 177.45x |

### brainpoolP512r1 (512-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w32-il64 | 479.85 M | w8 | 53.49 M | 36.36 M | 8.97x |
| SUBTRACT | w32-il64 | 478.61 M | w8 | 43.41 M | 43.01 M | 11.03x |
| ADDMOD | w32-il64 | 575.25 M | w8 | 44.98 M | 13.47 M | 12.79x |
| SUBTRACTMOD | w32-il64 | 553.73 M | w8 | 45.40 M | 15.01 M | 12.20x |
| MULTIPLYOPERANDSCANNING | w32-il64 | 216.17 M | w8 | 1.11 M | 13.90 M | 195.32x |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 62.57 M | w8 | 1.50 M | 13.94 M | 41.83x |
| MONTGOMERYMULTIPLICATION | w32 | 201.82 M | w8 | 401.97 k | 1.82 M | 502.08x |
| COMPARE | w32-il64 | 1.32 G | w8 | 108.86 M | 85.62 M | 12.14x |
| REDUCE | w32-opt | 21.46 M | w8 | 168.70 k | 19.19 M | 127.23x |
| MODMUL | w32-opt | 3.31 M | w8 | 30.74 k | 3.64 M | 107.65x |
| MODEXP | w32-opt | 10.44 k | w8 | 73.6 | 13.13 k | 141.91x |
| EXPONENTIATION | w32-o64 | 16.63 k | w8 | 168.0 | 53.43 k | 99.01x |
| DIVIDE | w32-il | 10.15 M | w8 | 106.41 k | 10.94 M | 95.37x |
| ISQRT | w32-il64 | 211.45 k | w8 | 3.37 k | 3.64 M | 62.71x |
| MODMUL_R2 | w32-opt | 89.93 M | w8 | 501.46 k | 3.61 M | 179.34x |

### p1024 (1024-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w32-il | 163.01 M | w8 | 29.85 M | 29.24 M | 5.46x |
| SUBTRACT | w32-il64 | 162.08 M | w8 | 27.46 M | 36.32 M | 5.90x |
| ADDMOD | w32-il64 | 266.79 M | w8 | 20.43 M | 9.24 M | 13.06x |
| SUBTRACTMOD | w32-il64 | 259.91 M | w8 | 22.10 M | 12.99 M | 11.76x |
| MULTIPLYOPERANDSCANNING | w32-il64 | 60.58 M | w8 | 321.42 k | 4.15 M | 188.49x |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 9.06 M | w8 | 475.77 k | 4.15 M | 19.05x |
| MONTGOMERYMULTIPLICATION | w32 | 49.71 M | w8 | 147.81 k | 621.64 k | 336.31x |
| COMPARE | w32-il | 547.69 M | w8 | 79.43 M | 75.76 M | 6.90x |
| REDUCE | w32-il64 | 3.74 M | w8 | 27.21 k | 26.01 M | 137.61x |
| MODMUL | w32-il64 | 400.16 k | w8 | 5.22 k | 1.40 M | 76.68x |
| MODEXP | w32-il64 | 727.8 | w8 | 4.8 | 2.05 k | 151.59x |
| EXPONENTIATION | w32-il64 | 922.3 | w8 | 11.6 | 14.01 k | 79.22x |
| DIVIDE | w32-il64 | 1.36 M | w8 | 19.89 k | 8.01 M | 68.12x |
| ISQRT | w32-il64 | 18.36 k | w8 | 395.9 | 1.40 M | 46.38x |
| MODMUL_R2 | w32-il64 | 24.81 M | w8 | 155.32 k | 1.37 M | 159.75x |

### p2048 (2048-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w32-il64 | 81.76 M | w8 | 14.18 M | 21.41 M | 5.76x |
| SUBTRACT | w32-il | 81.77 M | w8 | 11.84 M | 24.23 M | 6.91x |
| ADDMOD | w32-il64 | 113.96 M | w8 | 10.10 M | 7.86 M | 11.28x |
| SUBTRACTMOD | w32-il | 116.87 M | w8 | 9.91 M | 9.89 M | 11.79x |
| MULTIPLYOPERANDSCANNING | w32-il64 | 12.53 M | w8 | 93.35 k | 1.26 M | 134.19x |
| MULTIPLYPRODUCTSCANNING | w32 | 2.05 M | w8 | 122.49 k | 1.26 M | 16.71x |
| MONTGOMERYMULTIPLICATION | w32-il64 | 9.94 M | w8 | 33.21 k | 195.51 k | 299.32x |
| COMPARE | w32-il64 | 188.02 M | w8 | 31.63 M | 69.65 M | 5.94x |
| REDUCE | w32-il | 605.48 k | w8 | 6.02 k | 17.50 M | 100.56x |
| MODMUL | w32-il64 | 56.15 k | w8 | 716.5 | 472.99 k | 78.37x |
| MODEXP | w32-il64 | 32.2 | none | n/a | 292.6 | n/a |
| EXPONENTIATION | w32-o64 | 58.4 | none | n/a | 2.45 k | n/a |
| DIVIDE | w32-il64 | 163.55 k | none | n/a | 5.72 M | n/a |
| ISQRT | w32-il64 | 2.22 k | none | n/a | 774.36 k | n/a |
| MODMUL_R2 | w32-o64 | 4.64 M | none | n/a | 457.81 k | n/a |

## 6. Raw data

Also written to `NVIDIA_GeForce_GTX_1060_3GB_Report.csv` for analysis.

```csv
kind,device,device_type,kernel,modulus,bits,operation,items,seconds,ops_per_sec,mismatches
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,secp256k1,256,ADD,50000,0.001447355,34545765.936,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,secp256k1,256,ADD,50000,0.000126613,394903208.533,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,secp256k1,256,ADD,50000,0.000149122,335297038.434,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,secp256k1,256,ADD,50000,0.000052800,946969696.970,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,secp256k1,256,ADD,50000,0.000263060,190071059.060,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,secp256k1,256,ADD,50000,0.001217220,41077209.445,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,secp256k1,256,ADD,50000,0.000112908,442838570.038,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,secp256k1,256,ADD,50000,0.001152970,43366266.664,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32,secp256k1,256,ADD,50000,0.000069456,719878398.455,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32,secp256k1,256,ADD,50000,0.001063790,47001750.251,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,secp256k1,256,ADD,50000,0.000070518,709040006.339,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,secp256k1,256,ADD,50000,0.001051493,47551434.503,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,secp256k1,256,ADD,50000,0.000069898,715331919.203,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,secp256k1,256,ADD,50000,0.001086367,46024949.635,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,secp256k1,256,ADD,50000,0.000046307,1079745207.353,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,secp256k1,256,ADD,50000,0.001060881,47130651.946,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,secp256k1,256,ADD,50000,0.000042874,1166198001.564,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,secp256k1,256,ADD,50000,0.001015980,49213577.046,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,secp256k1,256,ADD,50000,0.000510816,97882693.388,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,secp256k1,256,ADD,50000,0.001563167,31986341.587,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,secp256k1,256,SUBTRACT,50000,0.001061603,47098566.877,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,secp256k1,256,SUBTRACT,50000,0.000142755,350250461.241,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,secp256k1,256,SUBTRACT,50000,0.000196904,253930921.750,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,secp256k1,256,SUBTRACT,50000,0.000052768,947543966.040,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,secp256k1,256,SUBTRACT,50000,0.000263393,189830459.380,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,secp256k1,256,SUBTRACT,50000,0.001323869,37768074.434,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,secp256k1,256,SUBTRACT,50000,0.000113450,440722821.304,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,secp256k1,256,SUBTRACT,50000,0.001090720,45841266.719,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32,secp256k1,256,SUBTRACT,50000,0.000069752,716821875.668,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32,secp256k1,256,SUBTRACT,50000,0.001129707,44259251.089,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,secp256k1,256,SUBTRACT,50000,0.000070455,709677345.671,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,secp256k1,256,SUBTRACT,50000,0.001047131,47749532.351,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,secp256k1,256,SUBTRACT,50000,0.000071993,694511024.294,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,secp256k1,256,SUBTRACT,50000,0.001089791,45880363.781,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,secp256k1,256,SUBTRACT,50000,0.000045775,1092311112.920,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,secp256k1,256,SUBTRACT,50000,0.001058228,47248783.029,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,secp256k1,256,SUBTRACT,50000,0.000042031,1189609820.519,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,secp256k1,256,SUBTRACT,50000,0.001015171,49252766.157,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,secp256k1,256,SUBTRACT,50000,0.000544496,91828059.277,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,secp256k1,256,SUBTRACT,50000,0.001666168,30008983.172,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,secp256k1,256,ADDMOD,50000,0.004286705,11663971.475,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,secp256k1,256,ADDMOD,50000,0.000636471,78558120.941,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,secp256k1,256,ADDMOD,50000,0.001762642,28366515.448,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,secp256k1,256,ADDMOD,50000,0.000052224,957414215.686,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,secp256k1,256,ADDMOD,50000,0.000364443,137195557.577,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,secp256k1,256,ADDMOD,50000,0.001396984,35791394.133,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,secp256k1,256,ADDMOD,50000,0.000157751,316954914.278,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,secp256k1,256,ADDMOD,50000,0.001137901,43940549.937,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32,secp256k1,256,ADDMOD,50000,0.000091087,548925311.848,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32,secp256k1,256,ADDMOD,50000,0.001073312,46584776.218,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,secp256k1,256,ADDMOD,50000,0.000071241,701847088.661,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,secp256k1,256,ADDMOD,50000,0.001055317,47379128.756,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,secp256k1,256,ADDMOD,50000,0.000074157,674240715.344,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,secp256k1,256,ADDMOD,50000,0.001089519,45891815.616,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,secp256k1,256,ADDMOD,50000,0.000039592,1262869100.489,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,secp256k1,256,ADDMOD,50000,0.001052491,47506327.913,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,secp256k1,256,ADDMOD,50000,0.000039632,1261622672.369,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,secp256k1,256,ADDMOD,50000,0.001014756,49272926.781,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,secp256k1,256,ADDMOD,50000,0.000625668,79914575.339,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,secp256k1,256,ADDMOD,50000,0.001742937,28687215.570,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,secp256k1,256,SUBTRACTMOD,50000,0.003036486,16466402.364,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,secp256k1,256,SUBTRACTMOD,50000,0.000493286,101361040.060,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,secp256k1,256,SUBTRACTMOD,50000,0.001729330,28912930.796,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,secp256k1,256,SUBTRACTMOD,50000,0.000052224,957414215.686,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.000364710,137095359.598,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.001346912,37121943.583,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,secp256k1,256,SUBTRACTMOD,50000,0.000158722,315017023.224,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,secp256k1,256,SUBTRACTMOD,50000,0.001145974,43631014.888,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32,secp256k1,256,SUBTRACTMOD,50000,0.000091111,548779425.534,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32,secp256k1,256,SUBTRACTMOD,50000,0.001077935,46384981.822,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,secp256k1,256,SUBTRACTMOD,50000,0.000072710,687661276.770,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,secp256k1,256,SUBTRACTMOD,50000,0.001052285,47515661.989,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,secp256k1,256,SUBTRACTMOD,50000,0.000074789,668548157.003,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,secp256k1,256,SUBTRACTMOD,50000,0.001086842,46004835.697,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,secp256k1,256,SUBTRACTMOD,50000,0.000039466,1266922106.853,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,secp256k1,256,SUBTRACTMOD,50000,0.001057593,47277159.377,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,secp256k1,256,SUBTRACTMOD,50000,0.000038993,1282294143.499,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,secp256k1,256,SUBTRACTMOD,50000,0.001021801,48933228.091,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.000584388,85559571.749,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.001558850,32074935.506,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001580179,31641990.117,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000255361,195801085.371,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000507975,98430041.398,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.009113029,5486650.152,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.010122972,4939261.111,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.002196014,22768525.907,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.003398674,14711621.578,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000628954,79497095.067,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001821734,27446377.843,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000119191,419495946.242,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001321113,37846883.614,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000121852,410332557.820,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001287980,38820478.424,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000052344,955218333.215,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001219839,40989020.547,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000055682,897957636.984,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001205856,41464321.451,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.011015134,4539209.447,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.012383658,4037579.060,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001582762,31590341.976,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000256117,195222946.575,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000506025,98809384.916,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000052224,957414215.686,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.005051045,9898942.200,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.006232260,8022772.099,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001082530,46188085.943,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.002287500,21857927.256,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000289120,172938703.775,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001423791,35117512.481,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000285298,175255572.966,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001419382,35226594.403,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000288906,173066926.276,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001460060,34245156.813,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000183815,272012419.314,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001353713,36935457.621,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000158021,316413187.642,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001356989,36846278.523,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.009436658,5298486.163,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.010909909,4582989.487,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.012285752,4069754.875,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.001867609,26772197.693,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000446048,112095651.230,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000064352,776976628.543,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000974676,51299117.290,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.001932582,25872124.925,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000241760,206816536.974,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.001254329,39861935.543,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000079406,629671966.409,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.001048895,47669232.012,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000086866,575597083.798,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.001102680,45344056.811,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000082932,602900583.955,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.001114594,44859391.069,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000066279,754392423.348,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.001088275,45944284.680,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000064777,771876401.070,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.001091955,45789423.085,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.030294795,1650448.557,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.038973592,1282919.976,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,secp256k1,256,COMPARE,50000,0.000542657,92139157.059,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,secp256k1,256,COMPARE,50000,0.000082081,609153007.920,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,secp256k1,256,COMPARE,50000,0.000129852,385052437.100,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,secp256k1,256,COMPARE,50000,0.000051904,963316892.725,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,secp256k1,256,COMPARE,50000,0.000174807,286029106.330,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,secp256k1,256,COMPARE,50000,0.001037667,48185034.815,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,secp256k1,256,COMPARE,50000,0.000089390,559345410.598,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,secp256k1,256,COMPARE,50000,0.001075532,46488609.024,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,secp256k1,256,COMPARE,50000,0.000052895,945261835.341,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,secp256k1,256,COMPARE,50000,0.001039132,48117060.300,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,secp256k1,256,COMPARE,50000,0.000052752,947831842.096,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,secp256k1,256,COMPARE,50000,0.001072785,46607666.330,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,secp256k1,256,COMPARE,50000,0.000020292,2464066972.645,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,secp256k1,256,COMPARE,50000,0.001036528,48237940.490,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,secp256k1,256,COMPARE,50000,0.000019928,2508977063.277,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,secp256k1,256,COMPARE,50000,0.000986444,50687123.130,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,secp256k1,256,COMPARE,50000,0.000289265,172851843.552,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,secp256k1,256,COMPARE,50000,0.001354469,36914835.631,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,secp256k1,256,REDUCE,6250,0.000188448,33165729.649,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,secp256k1,256,REDUCE,6250,0.000038404,162743389.272,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,secp256k1,256,REDUCE,6250,0.000215501,29002240.354,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,secp256k1,256,REDUCE,50000,0.000055296,904224537.037,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,secp256k1,256,REDUCE,6300,0.000196477,32064756.839,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,secp256k1,256,REDUCE,6300,0.000487557,12921574.547,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,secp256k1,256,REDUCE,6300,0.000105036,59979194.297,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,secp256k1,256,REDUCE,6300,0.000364142,17300964.954,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,secp256k1,256,REDUCE,6300,0.000049280,127840902.052,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,secp256k1,256,REDUCE,6300,0.000294641,21381977.606,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,secp256k1,256,REDUCE,6300,0.000050116,125707527.897,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,secp256k1,256,REDUCE,6300,0.000256279,24582537.453,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,secp256k1,256,REDUCE,6300,0.000047356,133035192.952,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,secp256k1,256,REDUCE,6300,0.000252109,24989189.107,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,secp256k1,256,REDUCE,6300,0.000045493,138482097.347,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,secp256k1,256,REDUCE,6300,0.000248283,25374255.384,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,secp256k1,256,REDUCE,6300,0.005220538,1206772.174,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,secp256k1,256,REDUCE,6300,0.005527735,1139707.363,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,secp256k1,256,MODMUL,3125,0.000441125,7084165.235,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,secp256k1,256,MODMUL,3125,0.000081548,38320769.283,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,secp256k1,256,MODMUL,3125,0.000260223,12008944.498,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,secp256k1,256,MODMUL,50000,0.000217088,230321344.340,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,secp256k1,256,MODMUL,6300,0.000575846,10940423.883,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,secp256k1,256,MODMUL,6300,0.000872543,7220273.855,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,secp256k1,256,MODMUL,6300,0.000266938,23601026.757,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,secp256k1,256,MODMUL,6300,0.000564506,11160194.761,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,secp256k1,256,MODMUL,6300,0.000121022,52056804.297,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,secp256k1,256,MODMUL,6300,0.000394881,15954182.762,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,secp256k1,256,MODMUL,6300,0.000107868,58404910.045,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,secp256k1,256,MODMUL,6300,0.000389656,16168105.096,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,secp256k1,256,MODMUL,6300,0.000118291,53258487.184,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,secp256k1,256,MODMUL,6300,0.000392571,16048048.708,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,secp256k1,256,MODMUL,6300,0.000104208,60456274.722,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,secp256k1,256,MODMUL,6300,0.000358313,17582377.244,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,secp256k1,256,MODMUL,6300,0.018282281,344595.946,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,secp256k1,256,MODMUL,6300,0.018552309,339580.374,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,secp256k1,256,MODEXP,781,0.011438826,68276.237,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,secp256k1,256,MODEXP,781,0.001588132,491772.652,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,secp256k1,256,MODEXP,781,0.002839707,275028.390,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,secp256k1,256,MODEXP,50000,0.055261184,904794.222,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,secp256k1,256,MODEXP,6300,0.024212465,260196.559,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,secp256k1,256,MODEXP,6300,0.024510939,257028.100,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,secp256k1,256,MODEXP,6300,0.005890895,1069447.060,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,secp256k1,256,MODEXP,6300,0.006197657,1016513.104,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,secp256k1,256,MODEXP,6300,0.003935836,1600676.347,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,secp256k1,256,MODEXP,6300,0.004155658,1516005.340,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,secp256k1,256,MODEXP,6300,0.003166649,1989484.475,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,secp256k1,256,MODEXP,6300,0.003491253,1804509.903,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,secp256k1,256,MODEXP,6300,0.003935063,1600990.781,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,secp256k1,256,MODEXP,6300,0.004272627,1474502.788,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,secp256k1,256,MODEXP,6300,0.003167959,1988662.145,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,secp256k1,256,MODEXP,6300,0.003476202,1812322.505,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,secp256k1,256,MODEXP,6300,0.634623285,9927.149,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,secp256k1,256,MODEXP,6300,0.633215841,9949.214,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,secp256k1,256,EXPONENTIATION,781,0.004042570,193193.955,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,secp256k1,256,EXPONENTIATION,781,0.000536941,1454535.995,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,secp256k1,256,EXPONENTIATION,781,0.009761838,80005.423,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,secp256k1,256,EXPONENTIATION,6300,0.119391607,52767.528,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,secp256k1,256,EXPONENTIATION,6300,0.119720075,52622.753,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,secp256k1,256,EXPONENTIATION,6300,0.021973701,286706.365,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,secp256k1,256,EXPONENTIATION,6300,0.022274500,282834.635,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,secp256k1,256,EXPONENTIATION,6300,0.000653651,9638176.555,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,secp256k1,256,EXPONENTIATION,6300,0.000922062,6832516.031,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,secp256k1,256,EXPONENTIATION,6300,0.000621205,10141576.475,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,secp256k1,256,EXPONENTIATION,6300,0.000845261,7453319.874,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,secp256k1,256,EXPONENTIATION,6300,0.000654478,9625997.511,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,secp256k1,256,EXPONENTIATION,6300,0.000904433,6965686.801,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,secp256k1,256,EXPONENTIATION,6300,0.000616061,10226267.202,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,secp256k1,256,EXPONENTIATION,6300,0.000923367,6822854.313,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,secp256k1,256,EXPONENTIATION,6300,0.268100051,23498.690,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,secp256k1,256,EXPONENTIATION,6300,0.281094281,22412.409,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,secp256k1,256,DIVIDE,6250,0.000435663,14345938.135,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,secp256k1,256,DIVIDE,6250,0.000057004,109640674.422,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,secp256k1,256,DIVIDE,6250,0.000172788,36171435.347,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,secp256k1,256,DIVIDE,50000,0.000070432,709904588.823,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,secp256k1,256,DIVIDE,6300,0.000355426,17725197.547,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,secp256k1,256,DIVIDE,6300,0.000688087,9155816.362,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,secp256k1,256,DIVIDE,6300,0.000224730,28033640.381,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,secp256k1,256,DIVIDE,6300,0.000553818,11375569.640,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,secp256k1,256,DIVIDE,6300,0.000112219,56140334.715,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,secp256k1,256,DIVIDE,6300,0.000400692,15722790.747,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,secp256k1,256,DIVIDE,6300,0.000117755,53501111.147,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,secp256k1,256,DIVIDE,6300,0.000411095,15324921.255,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,secp256k1,256,DIVIDE,6300,0.000107726,58481658.954,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,secp256k1,256,DIVIDE,6300,0.000387402,16262166.722,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,secp256k1,256,DIVIDE,6300,0.000104850,60085746.311,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,secp256k1,256,DIVIDE,6300,0.000428498,14702526.193,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,secp256k1,256,DIVIDE,6300,0.008530604,738517.484,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,secp256k1,256,DIVIDE,6300,0.008869469,710301.810,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,secp256k1,256,ISQRT,1562,0.000218775,7139751.431,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,secp256k1,256,ISQRT,1562,0.000032054,48729871.843,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,secp256k1,256,ISQRT,6300,0.005351823,1177169.015,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,secp256k1,256,ISQRT,6300,0.005651947,1114660.112,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,secp256k1,256,ISQRT,6300,0.002697499,2335496.289,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,secp256k1,256,ISQRT,6300,0.003576305,1761594.835,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,secp256k1,256,ISQRT,6300,0.000784589,8029682.011,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,secp256k1,256,ISQRT,6300,0.001023879,6153069.171,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,secp256k1,256,ISQRT,6300,0.000658046,9573792.184,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,secp256k1,256,ISQRT,6300,0.000863515,7295763.435,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,secp256k1,256,ISQRT,6300,0.000711421,8855521.522,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,secp256k1,256,ISQRT,6300,0.000953374,6608106.740,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,secp256k1,256,ISQRT,6300,0.000662418,9510609.942,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,secp256k1,256,ISQRT,6300,0.000951229,6623013.207,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,secp256k1,256,ISQRT,6300,0.057938054,108736.824,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,secp256k1,256,ISQRT,6300,0.061079914,103143.564,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,secp256k1,256,MODMUL_R2,50000,0.007093824,7048384.180,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,secp256k1,256,MODMUL_R2,50000,0.000961339,52010783.643,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,secp256k1,256,MODMUL_R2,50000,0.003501805,14278352.589,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,secp256k1,256,MODMUL_R2,50000,0.000111616,447964449.541,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,secp256k1,256,MODMUL_R2,50000,0.000942498,53050485.375,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,secp256k1,256,MODMUL_R2,50000,0.001997348,25033195.004,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,secp256k1,256,MODMUL_R2,50000,0.000249038,200772960.561,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,secp256k1,256,MODMUL_R2,50000,0.001229754,40658547.543,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,secp256k1,256,MODMUL_R2,50000,0.000272516,183475357.128,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,secp256k1,256,MODMUL_R2,50000,0.001257416,39764092.020,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,secp256k1,256,MODMUL_R2,50000,0.000244062,204865646.035,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,secp256k1,256,MODMUL_R2,50000,0.001271196,39333046.531,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,secp256k1,256,MODMUL_R2,50000,0.000250222,199822428.668,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,secp256k1,256,MODMUL_R2,50000,0.001184830,42200131.111,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,secp256k1,256,MODMUL_R2,50000,0.000207230,241277284.821,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,secp256k1,256,MODMUL_R2,50000,0.001254758,39848325.597,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,secp256k1,256,MODMUL_R2,50000,0.035808062,1396333.590,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,secp256k1,256,MODMUL_R2,50000,0.038969450,1283056.353,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,rsa256(composite),256,ADD,50000,0.001283498,38956025.849,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,rsa256(composite),256,ADD,50000,0.000187572,266564175.488,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,rsa256(composite),256,ADD,50000,0.000198359,252068639.253,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,rsa256(composite),256,ADD,50000,0.000052224,957414215.686,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,rsa256(composite),256,ADD,50000,0.000218578,228751624.229,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,rsa256(composite),256,ADD,50000,0.001199767,41674758.741,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,rsa256(composite),256,ADD,50000,0.000113197,441709102.876,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,rsa256(composite),256,ADD,50000,0.001129419,44270564.952,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32,rsa256(composite),256,ADD,50000,0.000069460,719839789.762,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32,rsa256(composite),256,ADD,50000,0.001110530,45023566.489,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,rsa256(composite),256,ADD,50000,0.000069719,717166593.641,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,rsa256(composite),256,ADD,50000,0.001048429,47690404.373,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,rsa256(composite),256,ADD,50000,0.000070179,712465047.642,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,rsa256(composite),256,ADD,50000,0.001084177,46117938.436,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,rsa256(composite),256,ADD,50000,0.000046168,1083012410.232,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,rsa256(composite),256,ADD,50000,0.000899313,55597994.263,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,rsa256(composite),256,ADD,50000,0.000042565,1174669420.620,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,rsa256(composite),256,ADD,50000,0.001022603,48894812.807,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,rsa256(composite),256,ADD,50000,0.000394657,126692210.685,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,rsa256(composite),256,ADD,50000,0.001496468,33412013.293,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACT,50000,0.001049764,47629732.358,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACT,50000,0.000165094,302858334.274,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACT,50000,0.000194704,256799854.589,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,rsa256(composite),256,SUBTRACT,50000,0.000052800,946969696.970,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,rsa256(composite),256,SUBTRACT,50000,0.000220066,227204631.518,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,rsa256(composite),256,SUBTRACT,50000,0.001209356,41344322.529,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,rsa256(composite),256,SUBTRACT,50000,0.000113137,441941810.998,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,rsa256(composite),256,SUBTRACT,50000,0.001079535,46316233.214,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32,rsa256(composite),256,SUBTRACT,50000,0.000069775,716592247.731,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32,rsa256(composite),256,SUBTRACT,50000,0.001053479,47461810.312,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,rsa256(composite),256,SUBTRACT,50000,0.000069266,721852948.611,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,rsa256(composite),256,SUBTRACT,50000,0.001045492,47824394.346,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,rsa256(composite),256,SUBTRACT,50000,0.000070324,710993129.387,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,rsa256(composite),256,SUBTRACT,50000,0.001092426,45769670.516,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,rsa256(composite),256,SUBTRACT,50000,0.000046074,1085201552.393,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,rsa256(composite),256,SUBTRACT,50000,0.001064289,46979704.788,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,rsa256(composite),256,SUBTRACT,50000,0.000042133,1186717312.113,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,rsa256(composite),256,SUBTRACT,50000,0.001022518,48898909.939,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,rsa256(composite),256,SUBTRACT,50000,0.000602491,82988763.989,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,rsa256(composite),256,SUBTRACT,50000,0.001687504,29629553.231,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,rsa256(composite),256,ADDMOD,50000,0.003330626,15012194.735,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,rsa256(composite),256,ADDMOD,50000,0.000537155,93082968.136,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,rsa256(composite),256,ADDMOD,50000,0.001587067,31504660.055,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,rsa256(composite),256,ADDMOD,50000,0.000052224,957414215.686,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,rsa256(composite),256,ADDMOD,50000,0.000283273,176508213.386,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,rsa256(composite),256,ADDMOD,50000,0.001272520,39292111.731,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,rsa256(composite),256,ADDMOD,50000,0.000146557,341165013.599,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,rsa256(composite),256,ADDMOD,50000,0.001099158,45489362.250,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32,rsa256(composite),256,ADDMOD,50000,0.000085266,586398095.113,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32,rsa256(composite),256,ADDMOD,50000,0.001124492,44464526.004,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,rsa256(composite),256,ADDMOD,50000,0.000073535,679944922.617,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,rsa256(composite),256,ADDMOD,50000,0.001069991,46729368.593,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,rsa256(composite),256,ADDMOD,50000,0.000073390,681290972.310,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,rsa256(composite),256,ADDMOD,50000,0.001089614,45887814.668,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,rsa256(composite),256,ADDMOD,50000,0.000039708,1259196247.303,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,rsa256(composite),256,ADDMOD,50000,0.001062410,47062812.140,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,rsa256(composite),256,ADDMOD,50000,0.000039257,1273654659.328,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,rsa256(composite),256,ADDMOD,50000,0.001018943,49070445.542,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,rsa256(composite),256,ADDMOD,50000,0.000661129,75628190.601,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,rsa256(composite),256,ADDMOD,50000,0.001693038,29532705.131,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACTMOD,50000,0.003039259,16451375.935,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACTMOD,50000,0.000490714,101892372.746,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACTMOD,50000,0.001741877,28704670.309,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,rsa256(composite),256,SUBTRACTMOD,50000,0.000052800,946969696.970,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,rsa256(composite),256,SUBTRACTMOD,50000,0.000303712,164629788.904,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,rsa256(composite),256,SUBTRACTMOD,50000,0.001289124,38786038.293,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,rsa256(composite),256,SUBTRACTMOD,50000,0.000157541,317378375.247,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,rsa256(composite),256,SUBTRACTMOD,50000,0.001140185,43852544.222,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32,rsa256(composite),256,SUBTRACTMOD,50000,0.000089664,557637326.021,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32,rsa256(composite),256,SUBTRACTMOD,50000,0.001064492,46970744.451,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,rsa256(composite),256,SUBTRACTMOD,50000,0.000072185,692665159.725,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,rsa256(composite),256,SUBTRACTMOD,50000,0.001045194,47838030.845,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,rsa256(composite),256,SUBTRACTMOD,50000,0.000073779,677696177.733,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,rsa256(composite),256,SUBTRACTMOD,50000,0.001099247,45485662.386,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,rsa256(composite),256,SUBTRACTMOD,50000,0.000039384,1269558531.971,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,rsa256(composite),256,SUBTRACTMOD,50000,0.000998126,50093858.714,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,rsa256(composite),256,SUBTRACTMOD,50000,0.000039048,1280459149.017,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,rsa256(composite),256,SUBTRACTMOD,50000,0.001036128,48256584.698,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,rsa256(composite),256,SUBTRACTMOD,50000,0.000664702,75221713.898,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,rsa256(composite),256,SUBTRACTMOD,50000,0.001737565,28775905.185,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001584690,31551911.083,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000255490,195702588.853,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000506146,98785749.351,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.008206690,6092590.514,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.009438662,5297361.083,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.002201585,22710909.526,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.003367595,14847389.439,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000629243,79460620.090,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001824139,27410196.594,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000117162,426758646.126,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001286019,38879685.296,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000120247,415811540.189,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001282159,38996716.225,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000053395,936424530.803,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001225067,40814084.169,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000056004,892790953.537,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001187393,42109041.738,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.008710414,5740255.172,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.011504056,4346293.194,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001580894,31627674.094,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000255192,195931138.280,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000515964,96906007.812,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000051200,976562500.000,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.004577631,10922680.000,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.005788656,8637583.649,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001080954,46255418.182,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.002259845,22125412.097,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000289356,172797322.124,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001424052,35111081.812,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000285156,175342575.706,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001424057,35110944.038,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000288561,173273596.695,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001467394,34074018.185,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000183497,272484576.811,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001356859,36849819.207,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000157259,317946008.433,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001289163,38784861.453,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.009406500,5315473.589,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.010790978,4633500.479,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.012295127,4066651.800,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.001858870,26898067.984,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000448216,111553419.716,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000064512,775049603.175,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000920692,54306950.736,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.001954015,25588335.396,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000238253,209861119.059,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.001256829,39782655.334,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000079939,625475816.110,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.001060572,47144392.420,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000085108,587488960.868,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.001067594,46834296.880,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000083957,595543896.703,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.001095653,45634903.710,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000065833,759493707.560,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.001085436,46064439.826,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000065120,767814009.897,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.001057066,47300735.144,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.036582025,1366791.486,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.039777344,1256996.945,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,rsa256(composite),256,COMPARE,50000,0.000542060,92240789.510,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,rsa256(composite),256,COMPARE,50000,0.000082063,609291272.669,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,rsa256(composite),256,COMPARE,50000,0.000130860,382087333.286,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,rsa256(composite),256,COMPARE,50000,0.000052224,957414215.686,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,rsa256(composite),256,COMPARE,50000,0.000166604,300112310.358,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,rsa256(composite),256,COMPARE,50000,0.001153866,43332594.430,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,rsa256(composite),256,COMPARE,50000,0.000088099,567540817.794,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,rsa256(composite),256,COMPARE,50000,0.001104882,45253701.834,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,rsa256(composite),256,COMPARE,50000,0.000052383,954505052.804,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,rsa256(composite),256,COMPARE,50000,0.001030713,48510093.971,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,rsa256(composite),256,COMPARE,50000,0.000052404,954131854.695,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,rsa256(composite),256,COMPARE,50000,0.001071388,46668438.118,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,rsa256(composite),256,COMPARE,50000,0.000020873,2395461859.718,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,rsa256(composite),256,COMPARE,50000,0.001051681,47542928.392,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,rsa256(composite),256,COMPARE,50000,0.000020310,2461807190.022,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,rsa256(composite),256,COMPARE,50000,0.001005542,49724450.306,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,rsa256(composite),256,COMPARE,50000,0.000272585,183428968.929,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,rsa256(composite),256,COMPARE,50000,0.001395328,35833869.213,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,rsa256(composite),256,REDUCE,6250,0.000307573,20320378.860,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,rsa256(composite),256,REDUCE,6250,0.000054555,114563255.830,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,rsa256(composite),256,REDUCE,6250,0.000215190,29044163.803,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,rsa256(composite),256,REDUCE,50000,0.000055296,904224537.037,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,rsa256(composite),256,REDUCE,6300,0.000188544,33413881.546,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,rsa256(composite),256,REDUCE,6300,0.000427542,14735385.652,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,rsa256(composite),256,REDUCE,6300,0.000106283,59275968.202,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,rsa256(composite),256,REDUCE,6300,0.000344066,18310443.136,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,rsa256(composite),256,REDUCE,6300,0.000049561,127115406.855,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,rsa256(composite),256,REDUCE,6300,0.000290493,21687303.203,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,rsa256(composite),256,REDUCE,6300,0.000050597,124513574.790,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,rsa256(composite),256,REDUCE,6300,0.000251556,25044143.748,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,rsa256(composite),256,REDUCE,6300,0.000047639,132244555.271,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,rsa256(composite),256,REDUCE,6300,0.000250380,25161705.268,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,rsa256(composite),256,REDUCE,6300,0.000046093,136680140.047,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,rsa256(composite),256,REDUCE,6300,0.000287989,21875823.803,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,rsa256(composite),256,REDUCE,6300,0.005254842,1198894.199,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,rsa256(composite),256,REDUCE,6300,0.005542010,1136771.663,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,rsa256(composite),256,MODMUL,3125,0.000440035,7101707.571,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,rsa256(composite),256,MODMUL,3125,0.000080070,39028580.734,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,rsa256(composite),256,MODMUL,3125,0.000261085,11969277.087,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,rsa256(composite),256,MODMUL,50000,0.000218112,229240023.474,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,rsa256(composite),256,MODMUL,6300,0.000547308,11510876.024,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,rsa256(composite),256,MODMUL,6300,0.000794699,7927525.649,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,rsa256(composite),256,MODMUL,6300,0.000265831,23699256.190,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,rsa256(composite),256,MODMUL,6300,0.000493709,12760553.334,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,rsa256(composite),256,MODMUL,6300,0.000119988,52505305.126,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,rsa256(composite),256,MODMUL,6300,0.000412831,15260478.826,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,rsa256(composite),256,MODMUL,6300,0.000106739,59022541.586,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,rsa256(composite),256,MODMUL,6300,0.000331772,18988916.093,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,rsa256(composite),256,MODMUL,6300,0.000119912,52538744.359,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,rsa256(composite),256,MODMUL,6300,0.000379803,16587561.588,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,rsa256(composite),256,MODMUL,6300,0.000102282,61594674.126,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,rsa256(composite),256,MODMUL,6300,0.000367824,17127757.303,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,rsa256(composite),256,MODMUL,6300,0.018321147,343864.930,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,rsa256(composite),256,MODMUL,6300,0.018587293,338941.231,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,rsa256(composite),256,MODEXP,781,0.010933070,71434.651,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,rsa256(composite),256,MODEXP,781,0.001546755,504927.918,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,rsa256(composite),256,MODEXP,781,0.002854526,273600.587,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,rsa256(composite),256,MODEXP,50000,0.055681024,897971.991,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,rsa256(composite),256,MODEXP,6300,0.024287036,259397.651,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,rsa256(composite),256,MODEXP,6300,0.024589261,256209.407,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,rsa256(composite),256,MODEXP,6300,0.005880192,1071393.601,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,rsa256(composite),256,MODEXP,6300,0.006151387,1024159.189,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,rsa256(composite),256,MODEXP,6300,0.003827265,1646084.285,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,rsa256(composite),256,MODEXP,6300,0.004099967,1536597.749,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,rsa256(composite),256,MODEXP,6300,0.003157940,1994971.556,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,rsa256(composite),256,MODEXP,6300,0.003488837,1805759.439,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,rsa256(composite),256,MODEXP,6300,0.003865983,1629598.333,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,rsa256(composite),256,MODEXP,6300,0.004248193,1482983.634,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,rsa256(composite),256,MODEXP,6300,0.003139116,2006934.519,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,rsa256(composite),256,MODEXP,6300,0.003505219,1797320.055,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,rsa256(composite),256,MODEXP,6300,0.634478120,9929.420,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,rsa256(composite),256,MODEXP,6300,0.622299096,10123.749,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,rsa256(composite),256,EXPONENTIATION,781,0.004026333,193973.033,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,rsa256(composite),256,EXPONENTIATION,781,0.000576690,1354280.840,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,rsa256(composite),256,EXPONENTIATION,781,0.010322407,75660.649,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,rsa256(composite),256,EXPONENTIATION,6300,0.120584909,52245.344,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,rsa256(composite),256,EXPONENTIATION,6300,0.120880015,52117.796,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,rsa256(composite),256,EXPONENTIATION,6300,0.021970127,286753.011,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,rsa256(composite),256,EXPONENTIATION,6300,0.022250138,283144.308,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,rsa256(composite),256,EXPONENTIATION,6300,0.000654148,9630848.979,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,rsa256(composite),256,EXPONENTIATION,6300,0.000942463,6684612.161,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,rsa256(composite),256,EXPONENTIATION,6300,0.000620825,10147783.685,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,rsa256(composite),256,EXPONENTIATION,6300,0.000849161,7419085.186,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,rsa256(composite),256,EXPONENTIATION,6300,0.000654850,9620521.505,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,rsa256(composite),256,EXPONENTIATION,6300,0.000931129,6765980.815,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,rsa256(composite),256,EXPONENTIATION,6300,0.000614256,10256315.637,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,rsa256(composite),256,EXPONENTIATION,6300,0.000879040,7166909.454,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,rsa256(composite),256,EXPONENTIATION,6300,0.255430369,24664.256,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,rsa256(composite),256,EXPONENTIATION,6300,0.263561772,23903.315,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,rsa256(composite),256,DIVIDE,6250,0.000459699,13595854.116,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,rsa256(composite),256,DIVIDE,6250,0.000054242,115224174.994,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,rsa256(composite),256,DIVIDE,6250,0.000163402,38249130.246,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,rsa256(composite),256,DIVIDE,50000,0.000071200,702247191.011,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,rsa256(composite),256,DIVIDE,6300,0.000357004,17646866.906,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,rsa256(composite),256,DIVIDE,6300,0.000627605,10038157.034,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,rsa256(composite),256,DIVIDE,6300,0.000226492,27815544.344,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,rsa256(composite),256,DIVIDE,6300,0.000524718,12006441.984,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,rsa256(composite),256,DIVIDE,6300,0.000111800,56350783.806,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,rsa256(composite),256,DIVIDE,6300,0.000447271,14085407.912,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,rsa256(composite),256,DIVIDE,6300,0.000114907,54827147.765,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,rsa256(composite),256,DIVIDE,6300,0.000394508,15969248.091,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,rsa256(composite),256,DIVIDE,6300,0.000108741,57935709.928,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,rsa256(composite),256,DIVIDE,6300,0.000409091,15400000.663,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,rsa256(composite),256,DIVIDE,6300,0.000105457,59739773.312,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,rsa256(composite),256,DIVIDE,6300,0.000396289,15897491.707,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,rsa256(composite),256,DIVIDE,6300,0.008403827,749658.503,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,rsa256(composite),256,DIVIDE,6300,0.008212000,767169.983,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,rsa256(composite),256,ISQRT,1562,0.000224756,6949756.471,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,rsa256(composite),256,ISQRT,1562,0.000032576,47949703.502,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,rsa256(composite),256,ISQRT,6300,0.005350396,1177482.930,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,rsa256(composite),256,ISQRT,6300,0.005633309,1118347.930,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,rsa256(composite),256,ISQRT,6300,0.002727805,2309549.502,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,rsa256(composite),256,ISQRT,6300,0.002992010,2105608.116,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,rsa256(composite),256,ISQRT,6300,0.000820639,7676948.060,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,rsa256(composite),256,ISQRT,6300,0.001120728,5621348.586,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,rsa256(composite),256,ISQRT,6300,0.000698013,9025616.878,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,rsa256(composite),256,ISQRT,6300,0.000912562,6903640.424,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,rsa256(composite),256,ISQRT,6300,0.000751086,8387859.084,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,rsa256(composite),256,ISQRT,6300,0.001016261,6199194.550,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,rsa256(composite),256,ISQRT,6300,0.000701465,8981207.303,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,rsa256(composite),256,ISQRT,6300,0.000994558,6334475.288,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,rsa256(composite),256,ISQRT,6300,0.070880299,88882.243,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,rsa256(composite),256,ISQRT,6300,0.069345431,90849.532,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,rsa256(composite),256,MODMUL_R2,50000,0.007111194,7031168.426,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,rsa256(composite),256,MODMUL_R2,50000,0.000973275,51372945.497,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,rsa256(composite),256,MODMUL_R2,50000,0.003510859,14241529.478,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,rsa256(composite),256,MODMUL_R2,50000,0.000112512,444397042.093,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,rsa256(composite),256,MODMUL_R2,50000,0.000949344,52667966.057,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,rsa256(composite),256,MODMUL_R2,50000,0.001863815,26826698.433,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,rsa256(composite),256,MODMUL_R2,50000,0.000247136,202317949.955,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,rsa256(composite),256,MODMUL_R2,50000,0.001259618,39694589.591,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,rsa256(composite),256,MODMUL_R2,50000,0.000273937,182523479.455,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,rsa256(composite),256,MODMUL_R2,50000,0.001302581,38385320.555,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,rsa256(composite),256,MODMUL_R2,50000,0.000247050,202388117.677,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,rsa256(composite),256,MODMUL_R2,50000,0.001265837,39499560.177,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,rsa256(composite),256,MODMUL_R2,50000,0.000251178,199062258.806,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,rsa256(composite),256,MODMUL_R2,50000,0.001281230,39025006.142,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,rsa256(composite),256,MODMUL_R2,50000,0.000202851,246485887.700,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,rsa256(composite),256,MODMUL_R2,50000,0.001249816,40005880.266,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,rsa256(composite),256,MODMUL_R2,50000,0.035996651,1389018.095,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,rsa256(composite),256,MODMUL_R2,50000,0.039309386,1271960.844,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,brainpoolP512r1,512,ADD,25000,0.000687648,36355830.521,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,brainpoolP512r1,512,ADD,25000,0.000106195,235416006.876,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,brainpoolP512r1,512,ADD,25000,0.000129825,192567652.334,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,brainpoolP512r1,512,ADD,50000,0.000088064,567768895.349,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,brainpoolP512r1,512,ADD,25000,0.000243731,102572163.972,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,brainpoolP512r1,512,ADD,25000,0.001262423,19803193.470,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,brainpoolP512r1,512,ADD,25000,0.000131898,189540936.564,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,brainpoolP512r1,512,ADD,25000,0.001151511,21710596.186,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32,brainpoolP512r1,512,ADD,25000,0.000076218,328008328.649,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32,brainpoolP512r1,512,ADD,25000,0.001068007,23408081.872,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,brainpoolP512r1,512,ADD,25000,0.000072448,345077074.174,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,brainpoolP512r1,512,ADD,25000,0.001041422,24005646.126,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,brainpoolP512r1,512,ADD,25000,0.000073073,342121607.912,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,brainpoolP512r1,512,ADD,25000,0.001086863,23001984.216,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,brainpoolP512r1,512,ADD,25000,0.000077154,324025223.311,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,brainpoolP512r1,512,ADD,25000,0.001059830,23588684.551,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,brainpoolP512r1,512,ADD,25000,0.000052100,479846011.941,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,brainpoolP512r1,512,ADD,25000,0.001028599,24304898.909,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,brainpoolP512r1,512,ADD,25000,0.000467356,53492378.941,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,brainpoolP512r1,512,ADD,25000,0.001555072,16076424.558,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACT,25000,0.000581304,43006785.972,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACT,25000,0.000098379,254118423.992,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACT,25000,0.000128230,194962055.692,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,brainpoolP512r1,512,SUBTRACT,50000,0.000088064,567768895.349,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,brainpoolP512r1,512,SUBTRACT,25000,0.000243319,102745694.360,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,brainpoolP512r1,512,SUBTRACT,25000,0.001236787,20213665.464,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,brainpoolP512r1,512,SUBTRACT,25000,0.000129934,192404782.247,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,brainpoolP512r1,512,SUBTRACT,25000,0.001150332,21732848.862,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32,brainpoolP512r1,512,SUBTRACT,25000,0.000075877,329481853.888,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32,brainpoolP512r1,512,SUBTRACT,25000,0.001067098,23428021.240,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,brainpoolP512r1,512,SUBTRACT,25000,0.000073239,341347222.787,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,brainpoolP512r1,512,SUBTRACT,25000,0.001047229,23872515.550,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,brainpoolP512r1,512,SUBTRACT,25000,0.000073329,340931030.278,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,brainpoolP512r1,512,SUBTRACT,25000,0.001083877,23065349.148,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,brainpoolP512r1,512,SUBTRACT,25000,0.000077400,322995928.190,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,brainpoolP512r1,512,SUBTRACT,25000,0.001049060,23830849.594,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,brainpoolP512r1,512,SUBTRACT,25000,0.000052234,478614014.192,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,brainpoolP512r1,512,SUBTRACT,25000,0.001046335,23892914.070,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,brainpoolP512r1,512,SUBTRACT,25000,0.000575937,43407500.550,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,brainpoolP512r1,512,SUBTRACT,25000,0.001571791,15905420.046,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,brainpoolP512r1,512,ADDMOD,25000,0.001855368,13474417.349,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,brainpoolP512r1,512,ADDMOD,25000,0.000298876,83646641.489,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,brainpoolP512r1,512,ADDMOD,25000,0.000972854,25697587.388,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,brainpoolP512r1,512,ADDMOD,50000,0.000087040,574448529.412,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,brainpoolP512r1,512,ADDMOD,25000,0.000330094,75735945.559,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,brainpoolP512r1,512,ADDMOD,25000,0.001356432,18430703.529,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,brainpoolP512r1,512,ADDMOD,25000,0.000162249,154083745.279,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,brainpoolP512r1,512,ADDMOD,25000,0.001171971,21331590.583,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32,brainpoolP512r1,512,ADDMOD,25000,0.000092421,270502091.983,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32,brainpoolP512r1,512,ADDMOD,25000,0.001078231,23186120.569,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,brainpoolP512r1,512,ADDMOD,25000,0.000096003,260409631.167,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,brainpoolP512r1,512,ADDMOD,25000,0.001083042,23083120.591,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,brainpoolP512r1,512,ADDMOD,25000,0.000095539,261673804.882,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,brainpoolP512r1,512,ADDMOD,25000,0.001110367,22515068.677,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,brainpoolP512r1,512,ADDMOD,25000,0.000049716,502857622.420,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,brainpoolP512r1,512,ADDMOD,25000,0.001030985,24248649.156,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,brainpoolP512r1,512,ADDMOD,25000,0.000043459,575251705.812,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,brainpoolP512r1,512,ADDMOD,25000,0.001030199,24267150.740,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,brainpoolP512r1,512,ADDMOD,25000,0.000555824,44978226.979,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,brainpoolP512r1,512,ADDMOD,25000,0.001689719,14795359.148,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACTMOD,25000,0.001665439,15011053.028,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000267992,93286437.721,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACTMOD,25000,0.001040546,24025842.807,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000087712,570047427.946,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000364279,68628645.351,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,brainpoolP512r1,512,SUBTRACTMOD,25000,0.001334798,18729431.465,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000176154,141921211.352,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,brainpoolP512r1,512,SUBTRACTMOD,25000,0.001193821,20941155.325,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000100972,247593071.261,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32,brainpoolP512r1,512,SUBTRACTMOD,25000,0.001074942,23257072.456,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000093903,266231062.800,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,brainpoolP512r1,512,SUBTRACTMOD,25000,0.001085177,23037715.070,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000095371,262133760.400,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,brainpoolP512r1,512,SUBTRACTMOD,25000,0.001116062,22400197.603,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000051072,489506283.964,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,brainpoolP512r1,512,SUBTRACTMOD,25000,0.001083085,23082207.551,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000045149,553726341.846,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,brainpoolP512r1,512,SUBTRACTMOD,25000,0.001037048,24106883.923,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000550617,45403649.403,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,brainpoolP512r1,512,SUBTRACTMOD,25000,0.001706870,14646690.143,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.001798162,13903083.653,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.000293612,85146244.417,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.000310177,80599146.074,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.022353178,1118409.209,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.023550430,1061551.728,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.005723687,4367814.048,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.006903825,3621180.840,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.004690422,5330011.107,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.005890258,4244296.507,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.000213699,116986749.645,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.001344528,18593885.055,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.000213452,117122524.346,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.001389729,17989120.566,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.000147061,169996995.681,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.001291193,19361937.900,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.000115652,216166416.492,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.001246497,20056204.947,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.022588847,1106740.852,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.023897342,1046141.433,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.001793329,13940556.676,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000293639,85138682.872,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000310130,80611248.048,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000104448,478707107.843,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.008071313,3097389.605,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.009205364,2715807.961,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.002051033,12188979.421,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.003221482,7760403.349,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000592358,42204178.354,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.001788281,13979906.591,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000548052,45616136.871,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.001738107,14383465.700,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000556938,44888271.729,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.001730153,14449586.325,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000726027,34433961.461,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.001913117,13067677.931,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000399565,62567991.982,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.001445889,17290394.921,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.016714377,1495718.297,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.018213941,1372575.013,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.013713608,1823006.811,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.004366063,5725982.973,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000473095,52843497.604,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000160640,311254980.080,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.002373856,10531390.541,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.003396936,7359573.969,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000513330,48701605.272,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.001535388,16282532.315,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000123871,201822065.170,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.001113502,22451682.232,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000161164,155121963.848,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.001136366,21999948.859,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000412261,60641179.049,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.001449058,17252589.539,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000161761,154548595.774,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.001146747,21800802.073,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000321416,77780775.271,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.001379622,18120903.657,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.062193183,401973.314,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.062907018,397411.940,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,brainpoolP512r1,512,COMPARE,25000,0.000291983,85621521.208,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,brainpoolP512r1,512,COMPARE,25000,0.000049209,508034853.704,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,brainpoolP512r1,512,COMPARE,25000,0.000090873,275109615.266,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,brainpoolP512r1,512,COMPARE,50000,0.000088064,567768895.349,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,brainpoolP512r1,512,COMPARE,25000,0.000177437,140894729.218,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,brainpoolP512r1,512,COMPARE,25000,0.001194561,20928192.102,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,brainpoolP512r1,512,COMPARE,25000,0.000094958,263275260.887,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,brainpoolP512r1,512,COMPARE,25000,0.001110127,22519941.946,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,brainpoolP512r1,512,COMPARE,25000,0.000054179,461435445.388,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,brainpoolP512r1,512,COMPARE,25000,0.001069142,23383246.079,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,brainpoolP512r1,512,COMPARE,25000,0.000053927,463587068.251,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,brainpoolP512r1,512,COMPARE,25000,0.001068743,23391967.263,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,brainpoolP512r1,512,COMPARE,25000,0.000021517,1161857063.712,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,brainpoolP512r1,512,COMPARE,25000,0.001000429,24989290.282,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,brainpoolP512r1,512,COMPARE,25000,0.000018911,1321951423.225,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,brainpoolP512r1,512,COMPARE,25000,0.001001505,24962427.024,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,brainpoolP512r1,512,COMPARE,25000,0.000229653,108859902.347,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,brainpoolP512r1,512,COMPARE,25000,0.001339374,18665434.707,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,brainpoolP512r1,512,REDUCE,3125,0.000162832,19191507.664,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,brainpoolP512r1,512,REDUCE,3125,0.000036137,86476037.318,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,brainpoolP512r1,512,REDUCE,3125,0.000152100,20545710.157,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,brainpoolP512r1,512,REDUCE,50000,0.000113536,440388951.522,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,brainpoolP512r1,512,REDUCE,6300,0.000611119,10308958.297,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,brainpoolP512r1,512,REDUCE,6300,0.001019152,6181610.529,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,brainpoolP512r1,512,REDUCE,6300,0.000350693,17964418.095,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,brainpoolP512r1,512,REDUCE,6300,0.000754774,8346873.632,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,brainpoolP512r1,512,REDUCE,6300,0.000145599,43269454.836,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,brainpoolP512r1,512,REDUCE,6300,0.000609294,10339843.068,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,brainpoolP512r1,512,REDUCE,6300,0.000152340,41354821.006,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,brainpoolP512r1,512,REDUCE,6300,0.000573067,10993479.052,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,brainpoolP512r1,512,REDUCE,6300,0.000178440,35306075.696,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,brainpoolP512r1,512,REDUCE,6300,0.000614125,10258493.159,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,brainpoolP512r1,512,REDUCE,6300,0.000157161,40086361.429,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,brainpoolP512r1,512,REDUCE,6300,0.000596773,10556780.004,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,brainpoolP512r1,512,REDUCE,6300,0.018524382,340092.312,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,brainpoolP512r1,512,REDUCE,6300,0.018960007,332278.366,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL,1562,0.000429239,3638997.269,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL,1562,0.000075104,20797905.919,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL,1562,0.000223843,6978093.318,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,brainpoolP512r1,512,MODMUL,50000,0.000567008,88182177.324,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,brainpoolP512r1,512,MODMUL,6300,0.001986133,3171993.118,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,brainpoolP512r1,512,MODMUL,6300,0.002356235,2673757.145,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,brainpoolP512r1,512,MODMUL,6300,0.001029518,6119371.444,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,brainpoolP512r1,512,MODMUL,6300,0.001435656,4388237.544,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,brainpoolP512r1,512,MODMUL,6300,0.000472115,13344196.420,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,brainpoolP512r1,512,MODMUL,6300,0.000905452,6957848.610,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,brainpoolP512r1,512,MODMUL,6300,0.000555791,11335196.944,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,brainpoolP512r1,512,MODMUL,6300,0.000977024,6448149.880,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,brainpoolP512r1,512,MODMUL,6300,0.000580121,10859806.536,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,brainpoolP512r1,512,MODMUL,6300,0.000977490,6445078.081,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,brainpoolP512r1,512,MODMUL,6300,0.000516033,12208526.278,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,brainpoolP512r1,512,MODMUL,6300,0.000956835,6584205.759,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,brainpoolP512r1,512,MODMUL,6300,0.050821500,123963.284,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,brainpoolP512r1,512,MODMUL,6300,0.055556310,113398.459,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,brainpoolP512r1,512,MODEXP,390,0.029694298,13133.835,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,brainpoolP512r1,512,MODEXP,390,0.004656041,83762.148,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,brainpoolP512r1,512,MODEXP,390,0.004389102,88856.452,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,brainpoolP512r1,512,MODEXP,50000,0.282700807,176865.431,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,brainpoolP512r1,512,MODEXP,6300,0.428206285,14712.535,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,brainpoolP512r1,512,MODEXP,6300,0.428977944,14686.070,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,brainpoolP512r1,512,MODEXP,6300,0.040689819,154829.887,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,brainpoolP512r1,512,MODEXP,6300,0.041100567,153282.556,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,brainpoolP512r1,512,MODEXP,6300,0.037343079,168705.959,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,brainpoolP512r1,512,MODEXP,6300,0.037865946,166376.406,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,brainpoolP512r1,512,MODEXP,6300,0.057152856,110230.712,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,brainpoolP512r1,512,MODEXP,6300,0.057564143,109443.130,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,brainpoolP512r1,512,MODEXP,6300,0.038099898,165354.773,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,brainpoolP512r1,512,MODEXP,6300,0.038655063,162979.943,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,brainpoolP512r1,512,MODEXP,6300,0.056417322,111667.832,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,brainpoolP512r1,512,MODEXP,6300,0.056847498,110822.820,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,brainpoolP512r1,512,MODEXP,6300,5.299319079,1188.832,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,brainpoolP512r1,512,MODEXP,6300,5.242849780,1201.637,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,brainpoolP512r1,512,EXPONENTIATION,390,0.007298751,53433.801,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,brainpoolP512r1,512,EXPONENTIATION,390,0.001019079,382698.406,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,brainpoolP512r1,512,EXPONENTIATION,390,0.012200475,31965.969,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,brainpoolP512r1,512,EXPONENTIATION,6300,0.975597324,6457.582,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,brainpoolP512r1,512,EXPONENTIATION,6300,0.977079289,6447.788,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,brainpoolP512r1,512,EXPONENTIATION,6300,0.241951499,26038.276,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,brainpoolP512r1,512,EXPONENTIATION,6300,0.242427506,25987.150,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,brainpoolP512r1,512,EXPONENTIATION,6300,0.024219250,260123.659,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,brainpoolP512r1,512,EXPONENTIATION,6300,0.024630085,255784.748,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,brainpoolP512r1,512,EXPONENTIATION,6300,0.023451528,268639.217,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,brainpoolP512r1,512,EXPONENTIATION,6300,0.023741169,265361.829,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,brainpoolP512r1,512,EXPONENTIATION,6300,0.026668375,236234.863,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,brainpoolP512r1,512,EXPONENTIATION,6300,0.027130406,232211.782,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,brainpoolP512r1,512,EXPONENTIATION,6300,0.023533210,267706.784,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,brainpoolP512r1,512,EXPONENTIATION,6300,0.023982635,262690.070,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,brainpoolP512r1,512,EXPONENTIATION,6300,2.321856018,2713.347,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,brainpoolP512r1,512,EXPONENTIATION,6300,2.352634456,2677.849,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,brainpoolP512r1,512,DIVIDE,3125,0.000285666,10939332.055,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,brainpoolP512r1,512,DIVIDE,3125,0.000041269,75723126.918,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,brainpoolP512r1,512,DIVIDE,3125,0.000112150,27864500.913,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,brainpoolP512r1,512,DIVIDE,50000,0.000158560,315338042.381,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,brainpoolP512r1,512,DIVIDE,6300,0.001569368,4014354.963,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,brainpoolP512r1,512,DIVIDE,6300,0.002055746,3064581.574,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,brainpoolP512r1,512,DIVIDE,6300,0.001992302,3162171.172,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,brainpoolP512r1,512,DIVIDE,6300,0.002480064,2540257.448,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,brainpoolP512r1,512,DIVIDE,6300,0.000388000,16237106.687,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,brainpoolP512r1,512,DIVIDE,6300,0.000854639,7371530.576,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,brainpoolP512r1,512,DIVIDE,6300,0.000345463,18236399.321,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,brainpoolP512r1,512,DIVIDE,6300,0.000788962,7985170.717,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,brainpoolP512r1,512,DIVIDE,6300,0.000307949,20457915.705,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,brainpoolP512r1,512,DIVIDE,6300,0.000805693,7819357.964,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,brainpoolP512r1,512,DIVIDE,6300,0.000336036,18747986.484,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,brainpoolP512r1,512,DIVIDE,6300,0.000850350,7408716.980,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,brainpoolP512r1,512,DIVIDE,6300,0.029368442,214515.979,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,brainpoolP512r1,512,DIVIDE,6300,0.033425631,188478.117,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,brainpoolP512r1,512,ISQRT,781,0.000214463,3641652.110,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,brainpoolP512r1,512,ISQRT,781,0.000028413,27487621.756,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,brainpoolP512r1,512,ISQRT,6300,0.038189756,164965.705,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,brainpoolP512r1,512,ISQRT,6300,0.038628567,163091.735,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,brainpoolP512r1,512,ISQRT,6300,0.023045151,273376.387,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,brainpoolP512r1,512,ISQRT,6300,0.023463955,268496.933,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,brainpoolP512r1,512,ISQRT,6300,0.004559811,1381636.151,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,brainpoolP512r1,512,ISQRT,6300,0.004967744,1268181.411,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,brainpoolP512r1,512,ISQRT,6300,0.003995966,1576589.925,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,brainpoolP512r1,512,ISQRT,6300,0.004413588,1427410.089,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,brainpoolP512r1,512,ISQRT,6300,0.005048292,1247946.890,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,brainpoolP512r1,512,ISQRT,6300,0.005488815,1147788.775,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,brainpoolP512r1,512,ISQRT,6300,0.003693579,1705662.828,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,brainpoolP512r1,512,ISQRT,6300,0.004174035,1509330.884,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,brainpoolP512r1,512,ISQRT,6300,0.231606105,27201.355,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,brainpoolP512r1,512,ISQRT,6300,0.232692676,27074.337,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL_R2,25000,0.006931501,3606722.666,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL_R2,25000,0.000814661,30687596.715,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL_R2,25000,0.002844620,8788518.551,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,brainpoolP512r1,512,MODMUL_R2,50000,0.000303104,164959881.757,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,brainpoolP512r1,512,MODMUL_R2,25000,0.002322750,10763103.310,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,brainpoolP512r1,512,MODMUL_R2,25000,0.003371656,7414754.185,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,brainpoolP512r1,512,MODMUL_R2,25000,0.000458883,54480103.588,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,brainpoolP512r1,512,MODMUL_R2,25000,0.001496691,16703511.754,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,brainpoolP512r1,512,MODMUL_R2,25000,0.000277998,89928728.501,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,brainpoolP512r1,512,MODMUL_R2,25000,0.001319854,18941494.882,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,brainpoolP512r1,512,MODMUL_R2,25000,0.000785137,31841593.103,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,brainpoolP512r1,512,MODMUL_R2,25000,0.001799472,13892966.657,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,brainpoolP512r1,512,MODMUL_R2,25000,0.000307953,81181223.122,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,brainpoolP512r1,512,MODMUL_R2,25000,0.001285525,19447306.933,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,brainpoolP512r1,512,MODMUL_R2,25000,0.000745792,33521413.408,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,brainpoolP512r1,512,MODMUL_R2,25000,0.001768779,14134043.803,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,brainpoolP512r1,512,MODMUL_R2,25000,0.049854765,501456.584,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,brainpoolP512r1,512,MODMUL_R2,25000,0.085664747,291835.332,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,p1024,1024,ADD,12500,0.000427425,29244903.082,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,p1024,1024,ADD,12500,0.000064261,194518446.377,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,p1024,1024,ADD,12500,0.000099016,126241772.795,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,p1024,1024,ADD,50000,0.000168480,296771130.104,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p1024,1024,ADD,12500,0.000366215,34132986.115,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p1024,1024,ADD,12500,0.001412915,8846957.826,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p1024,1024,ADD,12500,0.000187885,66530052.543,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p1024,1024,ADD,12500,0.001208495,10343440.720,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32,p1024,1024,ADD,12500,0.000130173,96026191.226,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32,p1024,1024,ADD,12500,0.001102008,11342931.127,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p1024,1024,ADD,12500,0.000101333,123355079.683,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p1024,1024,ADD,12500,0.001084283,11528355.665,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p1024,1024,ADD,12500,0.000100302,124624160.152,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p1024,1024,ADD,12500,0.001109652,11264790.680,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p1024,1024,ADD,12500,0.000076681,163012203.653,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p1024,1024,ADD,12500,0.001053201,11868579.301,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p1024,1024,ADD,12500,0.000076823,162711822.326,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p1024,1024,ADD,12500,0.001059074,11802764.050,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,p1024,1024,ADD,12500,0.000418758,29850174.585,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,p1024,1024,ADD,12500,0.001540454,8114490.372,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,p1024,1024,SUBTRACT,12500,0.000344139,36322575.477,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,p1024,1024,SUBTRACT,12500,0.000063447,197013956.492,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,p1024,1024,SUBTRACT,12500,0.000098981,126286910.049,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,p1024,1024,SUBTRACT,50000,0.000167776,298016402.823,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p1024,1024,SUBTRACT,12500,0.000364589,34285221.496,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p1024,1024,SUBTRACT,12500,0.001394048,8966690.539,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p1024,1024,SUBTRACT,12500,0.000187738,66582198.807,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p1024,1024,SUBTRACT,12500,0.001199089,10424580.742,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32,p1024,1024,SUBTRACT,12500,0.000129763,96329434.731,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32,p1024,1024,SUBTRACT,12500,0.001102066,11342336.821,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p1024,1024,SUBTRACT,12500,0.000101073,123673339.108,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p1024,1024,SUBTRACT,12500,0.000975264,12817039.254,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p1024,1024,SUBTRACT,12500,0.000100192,124760855.178,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p1024,1024,SUBTRACT,12500,0.001054017,11859392.688,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p1024,1024,SUBTRACT,12500,0.000077285,161739284.681,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p1024,1024,SUBTRACT,12500,0.001060303,11789079.589,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p1024,1024,SUBTRACT,12500,0.000077121,162083045.116,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p1024,1024,SUBTRACT,12500,0.001045940,11950967.259,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,p1024,1024,SUBTRACT,12500,0.000455253,27457270.510,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,p1024,1024,SUBTRACT,12500,0.001392411,8977234.013,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,p1024,1024,ADDMOD,12500,0.001353014,9238631.376,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,p1024,1024,ADDMOD,12500,0.000206018,60674349.261,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,p1024,1024,ADDMOD,12500,0.000705549,17716687.655,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,p1024,1024,ADDMOD,50000,0.000167552,298414820.474,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p1024,1024,ADDMOD,12500,0.000518285,24118016.762,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p1024,1024,ADDMOD,12500,0.001558928,8018331.475,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p1024,1024,ADDMOD,12500,0.000263648,47411681.091,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p1024,1024,ADDMOD,12500,0.001284771,9729362.847,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32,p1024,1024,ADDMOD,12500,0.000158217,79005514.351,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32,p1024,1024,ADDMOD,12500,0.001201313,10405281.650,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p1024,1024,ADDMOD,12500,0.000126835,98553270.479,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p1024,1024,ADDMOD,12500,0.001113106,11229842.334,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p1024,1024,ADDMOD,12500,0.000127751,97846300.994,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p1024,1024,ADDMOD,12500,0.001137877,10985371.252,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p1024,1024,ADDMOD,12500,0.000049809,250958692.644,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p1024,1024,ADDMOD,12500,0.001087768,11491420.929,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p1024,1024,ADDMOD,12500,0.000046853,266792017.174,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p1024,1024,ADDMOD,12500,0.001030510,12129912.824,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,p1024,1024,ADDMOD,12500,0.000611885,20428691.367,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,p1024,1024,ADDMOD,12500,0.001628246,7676971.898,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,p1024,1024,SUBTRACTMOD,12500,0.000962127,12992047.826,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,p1024,1024,SUBTRACTMOD,12500,0.000156172,80040150.756,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,p1024,1024,SUBTRACTMOD,12500,0.000677625,18446788.731,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,p1024,1024,SUBTRACTMOD,50000,0.000166912,299559049.080,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p1024,1024,SUBTRACTMOD,12500,0.000519700,24052321.764,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p1024,1024,SUBTRACTMOD,12500,0.001564607,7989226.564,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p1024,1024,SUBTRACTMOD,12500,0.000263410,47454594.568,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p1024,1024,SUBTRACTMOD,12500,0.001233980,10129823.333,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32,p1024,1024,SUBTRACTMOD,12500,0.000159875,78186298.816,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32,p1024,1024,SUBTRACTMOD,12500,0.001141604,10949505.787,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p1024,1024,SUBTRACTMOD,12500,0.000126721,98641635.677,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p1024,1024,SUBTRACTMOD,12500,0.001095211,11413324.444,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p1024,1024,SUBTRACTMOD,12500,0.000124643,100286719.369,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p1024,1024,SUBTRACTMOD,12500,0.001138445,10979889.332,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p1024,1024,SUBTRACTMOD,12500,0.000050601,247032555.400,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p1024,1024,SUBTRACTMOD,12500,0.001027456,12165976.386,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p1024,1024,SUBTRACTMOD,12500,0.000048093,259910395.043,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p1024,1024,SUBTRACTMOD,12500,0.001025517,12188979.421,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,p1024,1024,SUBTRACTMOD,12500,0.000565600,22100437.999,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,p1024,1024,SUBTRACTMOD,12500,0.001462992,8544132.573,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.003011836,4150292.739,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.000417121,29967341.252,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.000587998,21258577.226,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.055363819,225779.222,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.056536276,221096.982,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.013933213,897136.907,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.015113557,827072.026,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.012733333,981675.433,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.013979321,894177.885,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.000311464,40133039.900,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.001475513,8471629.974,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.000303470,41190287.496,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.001471939,8492202.252,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.000225903,55333369.613,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.001421170,8795567.934,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.000206325,60583970.389,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.001395863,8955036.443,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.038889652,321422.263,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.051728645,241645.609,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.003011039,4151391.584,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.000428231,29189823.623,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.000587484,21277179.989,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000348000,143678160.920,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.022923015,545303.485,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.024044229,519875.274,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.005816353,2149112.847,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.006993899,1787271.990,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.001548566,8071984.411,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.002731225,4576701.732,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.001696786,7366869.275,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.002879454,4341101.045,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.001670042,7484841.462,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.002837576,4405168.398,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.001429597,8743723.396,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.002627010,4758262.152,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.001379110,9063817.051,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.002568100,4867412.373,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.026272988,475773.826,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.025739979,485625.882,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.020108018,621642.568,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.002594059,4818702.425,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.000686117,18218479.949,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000572416,87349060.823,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.007100012,1760560.367,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.008064523,1549998.614,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.001257308,9941877.181,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.002278974,5484923.683,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.000251455,49710637.857,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.001247259,10021977.339,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.000388112,32207204.630,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.001436708,8700442.743,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.000346469,36078288.685,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.001359146,9196951.017,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.000306299,40809802.789,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.001342688,9309684.955,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.000264656,47231158.594,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.001262370,9902005.814,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.084568102,147809.868,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.108418245,115294.248,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,p1024,1024,COMPARE,12500,0.000164997,75759030.051,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,p1024,1024,COMPARE,12500,0.000029767,419929065.766,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,p1024,1024,COMPARE,12500,0.000070773,176620865.354,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,p1024,1024,COMPARE,50000,0.000166912,299559049.080,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p1024,1024,COMPARE,12500,0.000238299,52455027.514,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p1024,1024,COMPARE,12500,0.001277609,9783902.576,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p1024,1024,COMPARE,12500,0.000122890,101717085.001,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p1024,1024,COMPARE,12500,0.001147622,10892085.859,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p1024,1024,COMPARE,12500,0.000066746,187277066.473,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p1024,1024,COMPARE,12500,0.001051025,11893146.660,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p1024,1024,COMPARE,12500,0.000066884,186891121.756,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p1024,1024,COMPARE,12500,0.001079084,11583895.161,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p1024,1024,COMPARE,12500,0.000022823,547693332.245,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p1024,1024,COMPARE,12500,0.000987209,12661954.249,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p1024,1024,COMPARE,12500,0.000023123,540590172.386,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p1024,1024,COMPARE,12500,0.001006398,12420529.108,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,p1024,1024,COMPARE,12500,0.000157373,79429113.849,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,p1024,1024,COMPARE,12500,0.001161572,10761282.460,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,p1024,1024,REDUCE,1562,0.000060063,26006089.578,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,p1024,1024,REDUCE,1562,0.000020616,75767289.894,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,p1024,1024,REDUCE,1562,0.000116425,13416404.520,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,p1024,1024,REDUCE,50000,0.000222208,225014400.922,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p1024,1024,REDUCE,6300,0.003852244,1635410.228,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p1024,1024,REDUCE,6300,0.004558673,1381981.077,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p1024,1024,REDUCE,6300,0.001263846,4984785.675,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p1024,1024,REDUCE,6300,0.001939548,3248179.172,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p1024,1024,REDUCE,6300,0.000440806,14291996.593,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p1024,1024,REDUCE,6300,0.001149956,5478469.828,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p1024,1024,REDUCE,6300,0.000468872,13436488.955,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p1024,1024,REDUCE,6300,0.001154261,5458039.025,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p1024,1024,REDUCE,6300,0.000429669,14662435.985,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p1024,1024,REDUCE,6300,0.001141785,5517677.659,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p1024,1024,REDUCE,6300,0.000417108,15104012.118,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p1024,1024,REDUCE,6300,0.001129851,5575957.733,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,p1024,1024,REDUCE,6300,0.057399796,109756.487,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,p1024,1024,REDUCE,6300,0.057387333,109780.323,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,p1024,1024,MODMUL,781,0.000558734,1397803.040,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,p1024,1024,MODMUL,781,0.000088055,8869488.139,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,p1024,1024,MODMUL,781,0.000259785,3006332.372,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,p1024,1024,MODMUL,50000,0.002072384,24126802.755,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p1024,1024,MODMUL,6300,0.013885718,453703.588,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p1024,1024,MODMUL,6300,0.014597600,431577.784,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p1024,1024,MODMUL,6300,0.004487483,1403905.093,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p1024,1024,MODMUL,6300,0.005166540,1219384.752,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p1024,1024,MODMUL,6300,0.002073813,3037882.057,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p1024,1024,MODMUL,6300,0.002782835,2263878.632,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p1024,1024,MODMUL,6300,0.001976632,3187240.383,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p1024,1024,MODMUL,6300,0.002664430,2364483.151,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p1024,1024,MODMUL,6300,0.002080269,3028454.263,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p1024,1024,MODMUL,6300,0.002780424,2265841.118,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p1024,1024,MODMUL,6300,0.001951741,3227887.216,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p1024,1024,MODMUL,6300,0.002653698,2374046.019,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,p1024,1024,MODMUL,6300,0.149659688,42095.504,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,p1024,1024,MODMUL,6300,0.134218561,46938.367,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,p1024,1024,MODEXP,195,0.095151147,2049.371,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,p1024,1024,MODEXP,195,0.011156194,17479.079,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,p1024,1024,MODEXP,195,0.006176598,31570.775,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,p1024,1024,MODEXP,50000,1.843633175,27120.362,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p1024,1024,MODEXP,6300,3.609580567,1745.355,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p1024,1024,MODEXP,6300,3.615185320,1742.649,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p1024,1024,MODEXP,6300,0.489240279,12877.108,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p1024,1024,MODEXP,6300,0.490118438,12854.036,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p1024,1024,MODEXP,6300,0.303508887,20757.218,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p1024,1024,MODEXP,6300,0.306604631,20547.635,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p1024,1024,MODEXP,6300,0.279276066,22558.324,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p1024,1024,MODEXP,6300,0.280196240,22484.242,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p1024,1024,MODEXP,6300,0.309697928,20342.403,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p1024,1024,MODEXP,6300,0.313023757,20126.268,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p1024,1024,MODEXP,6300,0.267913615,23515.042,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p1024,1024,MODEXP,6300,0.269130569,23408.712,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,p1024,1024,MODEXP,6300,40.614216564,155.118,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,p1024,1024,MODEXP,6300,40.783727687,154.473,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,p1024,1024,EXPONENTIATION,195,0.013919197,14009.429,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,p1024,1024,EXPONENTIATION,195,0.001873951,104058.195,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,p1024,1024,EXPONENTIATION,195,0.023517508,8291.695,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p1024,1024,EXPONENTIATION,6300,7.977891285,789.682,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p1024,1024,EXPONENTIATION,6300,7.980658535,789.409,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p1024,1024,EXPONENTIATION,6300,1.936335890,3253.568,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p1024,1024,EXPONENTIATION,6300,1.937755270,3251.185,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p1024,1024,EXPONENTIATION,6300,0.227861211,27648.409,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p1024,1024,EXPONENTIATION,6300,0.229794564,27415.792,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p1024,1024,EXPONENTIATION,6300,0.221412042,28453.737,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p1024,1024,EXPONENTIATION,6300,0.219748324,28669.161,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p1024,1024,EXPONENTIATION,6300,0.236356409,26654.661,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p1024,1024,EXPONENTIATION,6300,0.242622666,25966.247,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p1024,1024,EXPONENTIATION,6300,0.211422933,29798.092,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p1024,1024,EXPONENTIATION,6300,0.212219853,29686.195,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,p1024,1024,EXPONENTIATION,6300,16.748579642,376.151,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,p1024,1024,EXPONENTIATION,6300,16.719847590,376.798,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,p1024,1024,DIVIDE,1562,0.000195082,8006878.039,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,p1024,1024,DIVIDE,1562,0.000024751,63108997.934,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,p1024,1024,DIVIDE,1562,0.000105314,14831842.316,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,p1024,1024,DIVIDE,50000,0.000296960,168372844.828,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p1024,1024,DIVIDE,6300,0.077064482,81749.722,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p1024,1024,DIVIDE,6300,0.079895975,78852.533,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p1024,1024,DIVIDE,6300,0.024495738,257187.601,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p1024,1024,DIVIDE,6300,0.025199832,250001.663,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p1024,1024,DIVIDE,6300,0.002713034,2322123.558,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p1024,1024,DIVIDE,6300,0.003481973,1809319.046,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p1024,1024,DIVIDE,6300,0.001163714,5413703.007,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p1024,1024,DIVIDE,6300,0.001926707,3269827.566,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p1024,1024,DIVIDE,6300,0.002574710,2446877.346,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p1024,1024,DIVIDE,6300,0.003316117,1899812.083,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p1024,1024,DIVIDE,6300,0.001152657,5465633.012,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p1024,1024,DIVIDE,6300,0.001942370,3243460.164,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,p1024,1024,DIVIDE,6300,0.078524563,80229.672,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,p1024,1024,DIVIDE,6300,0.085149428,73987.579,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,p1024,1024,ISQRT,390,0.000278419,1400767.056,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,p1024,1024,ISQRT,390,0.000032010,12183861.256,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p1024,1024,ISQRT,6300,1.028706882,6124.193,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p1024,1024,ISQRT,6300,1.030220335,6115.197,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p1024,1024,ISQRT,6300,0.244301297,25787.829,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p1024,1024,ISQRT,6300,0.245100787,25703.712,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p1024,1024,ISQRT,6300,0.027895868,225839.900,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p1024,1024,ISQRT,6300,0.028585322,220392.828,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p1024,1024,ISQRT,6300,0.022255145,283080.609,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p1024,1024,ISQRT,6300,0.022933997,274701.350,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p1024,1024,ISQRT,6300,0.022118570,284828.536,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p1024,1024,ISQRT,6300,0.022811595,276175.336,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p1024,1024,ISQRT,6300,0.021239685,296614.570,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p1024,1024,ISQRT,6300,0.021954399,286958.441,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,p1024,1024,ISQRT,6300,0.985089395,6395.359,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,p1024,1024,ISQRT,6300,0.979181793,6433.943,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,p1024,1024,MODMUL_R2,12500,0.009116340,1371164.241,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,p1024,1024,MODMUL_R2,12500,0.000971010,12873194.490,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,p1024,1024,MODMUL_R2,12500,0.003363425,3716449.847,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,p1024,1024,MODMUL_R2,50000,0.001089184,45905925.904,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p1024,1024,MODMUL_R2,12500,0.007540008,1657823.163,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p1024,1024,MODMUL_R2,12500,0.008585740,1455902.482,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p1024,1024,MODMUL_R2,12500,0.000939127,13310233.047,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p1024,1024,MODMUL_R2,12500,0.001984125,6300006.290,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p1024,1024,MODMUL_R2,12500,0.000666833,18745335.654,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p1024,1024,MODMUL_R2,12500,0.001683006,7427186.553,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p1024,1024,MODMUL_R2,12500,0.000574360,21763366.488,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p1024,1024,MODMUL_R2,12500,0.001587214,7875434.823,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p1024,1024,MODMUL_R2,12500,0.000600012,20832916.006,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p1024,1024,MODMUL_R2,12500,0.001594692,7838501.934,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p1024,1024,MODMUL_R2,12500,0.000503786,24812127.128,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p1024,1024,MODMUL_R2,12500,0.001553304,8047359.711,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,p1024,1024,MODMUL_R2,12500,0.080480924,155316.308,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,p1024,1024,MODMUL_R2,12500,0.089085247,140315.039,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,p2048,2048,ADD,6250,0.000291863,21414123.158,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,p2048,2048,ADD,6250,0.000042338,147621786.186,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,p2048,2048,ADD,6250,0.000101170,61777468.471,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,p2048,2048,ADD,50000,0.000330720,151185292.695,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p2048,2048,ADD,6300,0.000403600,15609521.581,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p2048,2048,ADD,6300,0.001440113,4374655.464,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p2048,2048,ADD,6300,0.000204321,30833835.448,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p2048,2048,ADD,6300,0.001211241,5201277.520,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32,p2048,2048,ADD,6300,0.000106921,58921776.660,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32,p2048,2048,ADD,6300,0.001090270,5778387.631,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p2048,2048,ADD,6300,0.000107205,58766167.068,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p2048,2048,ADD,6300,0.001127550,5587333.498,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p2048,2048,ADD,6300,0.000108344,58147862.974,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p2048,2048,ADD,6300,0.001111642,5667294.583,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p2048,2048,ADD,6300,0.000076527,82324126.703,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p2048,2048,ADD,6300,0.001098394,5735644.715,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p2048,2048,ADD,6300,0.000076443,82414394.386,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p2048,2048,ADD,6300,0.001041630,6048211.254,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,p2048,2048,ADD,6300,0.000440631,14297675.644,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,p2048,2048,ADD,6300,0.001411451,4463491.740,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,p2048,2048,SUBTRACT,6250,0.000257967,24227901.368,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,p2048,2048,SUBTRACT,6250,0.000043293,144363601.945,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,p2048,2048,SUBTRACT,6250,0.000098227,63628391.012,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,p2048,2048,SUBTRACT,50000,0.000329728,151640139.752,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p2048,2048,SUBTRACT,6300,0.000402797,15640632.350,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p2048,2048,SUBTRACT,6300,0.001418468,4441412.656,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p2048,2048,SUBTRACT,6300,0.000203107,31018201.662,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p2048,2048,SUBTRACT,6300,0.001225239,5141855.359,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32,p2048,2048,SUBTRACT,6300,0.000106957,58902280.408,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32,p2048,2048,SUBTRACT,6300,0.001102380,5714905.389,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p2048,2048,SUBTRACT,6300,0.000107067,58841821.569,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p2048,2048,SUBTRACT,6300,0.001126043,5594810.510,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p2048,2048,SUBTRACT,6300,0.000108574,58025162.903,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p2048,2048,SUBTRACT,6300,0.001115449,5647950.989,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p2048,2048,SUBTRACT,6300,0.000076436,82422427.761,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p2048,2048,SUBTRACT,6300,0.001046285,6021303.756,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p2048,2048,SUBTRACT,6300,0.000076657,82184102.675,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p2048,2048,SUBTRACT,6300,0.001039172,6062521.389,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,p2048,2048,SUBTRACT,6300,0.000527959,11932737.614,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,p2048,2048,SUBTRACT,6300,0.001646783,3825640.187,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,p2048,2048,ADDMOD,6250,0.000795593,7855770.712,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,p2048,2048,ADDMOD,6250,0.000124861,50055840.320,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,p2048,2048,ADDMOD,6250,0.000515403,12126427.784,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,p2048,2048,ADDMOD,50000,0.000331648,150762253.956,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p2048,2048,ADDMOD,6300,0.000544207,11576473.704,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p2048,2048,ADDMOD,6300,0.001571637,4008560.129,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p2048,2048,ADDMOD,6300,0.000269752,23354785.499,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p2048,2048,ADDMOD,6300,0.001233865,5105908.804,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32,p2048,2048,ADDMOD,6300,0.000141134,44638275.140,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32,p2048,2048,ADDMOD,6300,0.001142250,5515428.267,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p2048,2048,ADDMOD,6300,0.000134587,46809770.062,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p2048,2048,ADDMOD,6300,0.001153424,5461996.553,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p2048,2048,ADDMOD,6300,0.000132956,47384235.719,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p2048,2048,ADDMOD,6300,0.001138594,5533139.988,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p2048,2048,ADDMOD,6300,0.000070296,89620740.477,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p2048,2048,ADDMOD,6300,0.001026509,6137304.180,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p2048,2048,ADDMOD,6300,0.000054844,114871849.803,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p2048,2048,ADDMOD,6300,0.001021205,6169185.408,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,p2048,2048,ADDMOD,6300,0.000618540,10185279.110,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,p2048,2048,ADDMOD,6300,0.001691097,3725391.585,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,p2048,2048,SUBTRACTMOD,6250,0.000631955,9889952.193,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,p2048,2048,SUBTRACTMOD,6250,0.000101943,61309029.783,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,p2048,2048,SUBTRACTMOD,6250,0.000504183,12396301.181,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,p2048,2048,SUBTRACTMOD,50000,0.000331776,150704089.506,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p2048,2048,SUBTRACTMOD,6300,0.000557162,11307306.487,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p2048,2048,SUBTRACTMOD,6300,0.001576897,3995188.634,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p2048,2048,SUBTRACTMOD,6300,0.000276854,22755654.762,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p2048,2048,SUBTRACTMOD,6300,0.001256736,5012986.039,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32,p2048,2048,SUBTRACTMOD,6300,0.000142720,44142501.443,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32,p2048,2048,SUBTRACTMOD,6300,0.001134545,5552888.744,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p2048,2048,SUBTRACTMOD,6300,0.000129905,48497128.640,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p2048,2048,SUBTRACTMOD,6300,0.001138708,5532587.886,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p2048,2048,SUBTRACTMOD,6300,0.000132324,47610348.188,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p2048,2048,SUBTRACTMOD,6300,0.001066463,5907377.563,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p2048,2048,SUBTRACTMOD,6300,0.000053478,117804560.816,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p2048,2048,SUBTRACTMOD,6300,0.001022292,6162620.997,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p2048,2048,SUBTRACTMOD,6300,0.000056434,111633993.848,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p2048,2048,SUBTRACTMOD,6300,0.001020433,6173847.422,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,p2048,2048,SUBTRACTMOD,6300,0.000630509,9991925.441,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,p2048,2048,SUBTRACTMOD,6300,0.001590658,3960625.151,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.004975693,1256106.360,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.000659486,9477073.595,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.000898387,6956910.586,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,6300,0.123433614,51039.581,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,6300,0.124982098,50407.219,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p2048,2048,MULTIPLYOPERANDSCANNING,6300,0.033460351,188282.545,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p2048,2048,MULTIPLYOPERANDSCANNING,6300,0.034587178,182148.423,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32,p2048,2048,MULTIPLYOPERANDSCANNING,6300,0.009532450,660900.414,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32,p2048,2048,MULTIPLYOPERANDSCANNING,6300,0.010744391,586352.451,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p2048,2048,MULTIPLYOPERANDSCANNING,6300,0.000708690,8889642.541,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p2048,2048,MULTIPLYOPERANDSCANNING,6300,0.001880670,3349870.053,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p2048,2048,MULTIPLYOPERANDSCANNING,6300,0.000546880,11519893.277,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p2048,2048,MULTIPLYOPERANDSCANNING,6300,0.001751943,3596008.307,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p2048,2048,MULTIPLYOPERANDSCANNING,6300,0.000630302,9995203.021,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p2048,2048,MULTIPLYOPERANDSCANNING,6300,0.001803976,3492286.795,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p2048,2048,MULTIPLYOPERANDSCANNING,6300,0.000498954,12626409.228,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p2048,2048,MULTIPLYOPERANDSCANNING,6300,0.001689900,3728031.878,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,6300,0.066953955,94094.516,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,6300,0.070808182,88972.769,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.004970046,1257533.693,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.000662129,9439242.935,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.000890264,7020387.190,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.001241664,40268542.859,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,6300,0.047258060,133310.592,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,6300,0.048606636,129611.932,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p2048,2048,MULTIPLYPRODUCTSCANNING,6300,0.011841785,532014.381,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p2048,2048,MULTIPLYPRODUCTSCANNING,6300,0.013091199,481239.339,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32,p2048,2048,MULTIPLYPRODUCTSCANNING,6300,0.003053084,2063487.270,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32,p2048,2048,MULTIPLYPRODUCTSCANNING,6300,0.004243122,1484755.655,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p2048,2048,MULTIPLYPRODUCTSCANNING,6300,0.003371188,1868777.185,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p2048,2048,MULTIPLYPRODUCTSCANNING,6300,0.004567293,1379372.722,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p2048,2048,MULTIPLYPRODUCTSCANNING,6300,0.003357930,1876155.785,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p2048,2048,MULTIPLYPRODUCTSCANNING,6300,0.004527258,1391570.867,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p2048,2048,MULTIPLYPRODUCTSCANNING,6300,0.003449742,1826223.658,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p2048,2048,MULTIPLYPRODUCTSCANNING,6300,0.004574528,1377191.275,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p2048,2048,MULTIPLYPRODUCTSCANNING,6300,0.003213130,1960705.020,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p2048,2048,MULTIPLYPRODUCTSCANNING,6300,0.004416505,1426467.349,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,6300,0.051025720,123467.145,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,6300,0.059116665,106568.934,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.031967930,195508.433,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.004206050,1485954.689,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.001210442,5163403.919,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.002159136,23157411.113,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,6300,0.081140833,77642.782,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,6300,0.082269894,76577.223,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p2048,2048,MONTGOMERYMULTIPLICATION,6300,0.003207689,1964030.710,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p2048,2048,MONTGOMERYMULTIPLICATION,6300,0.004237233,1486819.441,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32,p2048,2048,MONTGOMERYMULTIPLICATION,6300,0.000703862,8950618.963,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32,p2048,2048,MONTGOMERYMULTIPLICATION,6300,0.001726776,3648416.865,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p2048,2048,MONTGOMERYMULTIPLICATION,6300,0.000760484,8284192.655,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p2048,2048,MONTGOMERYMULTIPLICATION,6300,0.001776490,3546318.145,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p2048,2048,MONTGOMERYMULTIPLICATION,6300,0.000674786,9336292.646,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p2048,2048,MONTGOMERYMULTIPLICATION,6300,0.001616225,3897972.979,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p2048,2048,MONTGOMERYMULTIPLICATION,6300,0.000745159,8454575.496,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p2048,2048,MONTGOMERYMULTIPLICATION,6300,0.001764379,3570660.512,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p2048,2048,MONTGOMERYMULTIPLICATION,6300,0.000628782,10019363.832,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p2048,2048,MONTGOMERYMULTIPLICATION,6300,0.001647286,3824472.224,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,6300,0.188204734,33474.185,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,6300,0.207588227,30348.542,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,p2048,2048,COMPARE,6250,0.000089739,69646793.142,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,p2048,2048,COMPARE,6250,0.000022223,281237381.611,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,p2048,2048,COMPARE,6250,0.000076560,81634995.013,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,p2048,2048,COMPARE,50000,0.000331328,150907861.696,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p2048,2048,COMPARE,6300,0.000248294,25373113.274,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p2048,2048,COMPARE,6300,0.001205234,5227201.313,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p2048,2048,COMPARE,6300,0.000124872,50451771.265,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p2048,2048,COMPARE,6300,0.001107570,5688129.174,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p2048,2048,COMPARE,6300,0.000068739,91650952.352,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p2048,2048,COMPARE,6300,0.001090730,5775950.288,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p2048,2048,COMPARE,6300,0.000068137,92460204.631,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p2048,2048,COMPARE,6300,0.001024982,6146449.637,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p2048,2048,COMPARE,6300,0.000034027,185148168.688,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p2048,2048,COMPARE,6300,0.001034658,6088966.162,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p2048,2048,COMPARE,6300,0.000033241,189526322.179,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p2048,2048,COMPARE,6300,0.000991294,6355328.470,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,p2048,2048,COMPARE,6300,0.000197602,31882198.060,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,p2048,2048,COMPARE,6300,0.001152704,5465412.215,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,p2048,2048,REDUCE,781,0.000044640,17495459.496,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,p2048,2048,REDUCE,781,0.000018530,42148791.945,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,p2048,2048,REDUCE,781,0.000113014,6910640.180,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,p2048,2048,REDUCE,50000,0.000378688,132034814.940,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p2048,2048,REDUCE,6300,0.390049620,16151.791,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p2048,2048,REDUCE,6300,0.390747877,16122.928,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p2048,2048,REDUCE,6300,0.007483032,841904.770,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p2048,2048,REDUCE,6300,0.008572921,734872.033,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p2048,2048,REDUCE,6300,0.001468740,4289389.728,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p2048,2048,REDUCE,6300,0.002499267,2520738.618,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p2048,2048,REDUCE,6300,0.001398182,4505852.643,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p2048,2048,REDUCE,6300,0.002422232,2600906.738,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p2048,2048,REDUCE,6300,0.001289878,4884182.692,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p2048,2048,REDUCE,6300,0.002321279,2714021.399,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p2048,2048,REDUCE,6300,0.001362454,4624008.828,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p2048,2048,REDUCE,6300,0.002398925,2626176.415,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,p2048,2048,REDUCE,6300,0.129708856,48570.315,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,p2048,2048,REDUCE,6300,0.153340165,41085.126,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,p2048,2048,MODMUL,390,0.000824545,472988.374,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,p2048,2048,MODMUL,390,0.000130929,2978712.452,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,p2048,2048,MODMUL,390,0.000381326,1022746.129,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,p2048,2048,MODMUL,50000,0.007330752,6820582.663,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p2048,2048,MODMUL,6300,0.973220401,6473.354,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p2048,2048,MODMUL,6300,0.974712349,6463.445,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p2048,2048,MODMUL,6300,0.026211048,240356.666,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p2048,2048,MODMUL,6300,0.027221460,231435.052,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p2048,2048,MODMUL,6300,0.010947922,575451.657,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p2048,2048,MODMUL,6300,0.011976715,526020.692,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p2048,2048,MODMUL,6300,0.007138656,882519.006,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p2048,2048,MODMUL,6300,0.008197835,768495.611,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p2048,2048,MODMUL,6300,0.011911629,528894.925,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p2048,2048,MODMUL,6300,0.012911849,487923.932,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p2048,2048,MODMUL,6300,0.006945560,907054.323,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p2048,2048,MODMUL,6300,0.007931793,794271.852,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,p2048,2048,MODMUL,6300,0.544292999,11574.648,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,p2048,2048,MODMUL,6300,0.576493623,10928.135,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,p2048,2048,MODEXP,97,0.331504490,292.605,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,p2048,2048,MODEXP,97,0.030087037,3223.980,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,p2048,2048,MODEXP,97,0.014717605,6590.746,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,p2048,2048,MODEXP,50000,12.551837921,3983.480,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p2048,2048,MODEXP,6300,107.140816340,58.801,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p2048,2048,MODEXP,6300,107.151448937,58.795,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p2048,2048,MODEXP,6300,16.450120445,382.976,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p2048,2048,MODEXP,6300,16.400271410,384.140,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p2048,2048,MODEXP,6300,5.276535261,1193.965,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p2048,2048,MODEXP,6300,5.175839795,1217.194,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p2048,2048,MODEXP,6300,3.066316526,2054.582,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p2048,2048,MODEXP,6300,3.024475630,2083.006,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p2048,2048,MODEXP,6300,5.219526570,1207.006,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p2048,2048,MODEXP,6300,5.497608488,1145.953,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p2048,2048,MODEXP,6300,3.008573443,2094.016,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p2048,2048,MODEXP,6300,3.023180665,2083.898,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,p2048,2048,MODEXP,6300,0.000000000,inf,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,p2048,2048,MODEXP,6300,0.000000000,inf,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,p2048,2048,EXPONENTIATION,97,0.039640121,2447.016,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,p2048,2048,EXPONENTIATION,97,0.004751116,20416.256,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,p2048,2048,EXPONENTIATION,97,0.054580262,1777.199,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p2048,2048,EXPONENTIATION,6300,68.612949066,91.819,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p2048,2048,EXPONENTIATION,6300,68.581670670,91.861,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p2048,2048,EXPONENTIATION,6300,15.950623930,394.969,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p2048,2048,EXPONENTIATION,6300,15.956176817,394.831,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p2048,2048,EXPONENTIATION,6300,1.912643425,3293.871,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p2048,2048,EXPONENTIATION,6300,1.912305087,3294.453,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p2048,2048,EXPONENTIATION,6300,1.662305288,3789.918,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p2048,2048,EXPONENTIATION,6300,1.668777268,3775.219,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p2048,2048,EXPONENTIATION,6300,1.900033200,3315.732,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p2048,2048,EXPONENTIATION,6300,1.903256398,3310.116,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p2048,2048,EXPONENTIATION,6300,1.688371986,3731.405,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p2048,2048,EXPONENTIATION,6300,1.679140467,3751.920,0
opencl-kernel,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,p2048,2048,EXPONENTIATION,6300,0.000000000,inf,0
opencl-e2e,cpu-haswell-Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,CPU,w8,p2048,2048,EXPONENTIATION,6300,0.000000000,inf,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,p2048,2048,DIVIDE,781,0.000136595,5717623.234,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,p2048,2048,DIVIDE,781,0.000024687,31635444.566,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,p2048,2048,DIVIDE,781,0.000123322,6333013.870,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,p2048,2048,DIVIDE,50000,0.000476736,104879849.644,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p2048,2048,DIVIDE,6300,1.110800717,5671.584,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p2048,2048,DIVIDE,6300,1.115655683,5646.904,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p2048,2048,DIVIDE,6300,0.385384532,16347.309,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p2048,2048,DIVIDE,6300,0.388686214,16208.447,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p2048,2048,DIVIDE,6300,0.049040994,128463.955,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p2048,2048,DIVIDE,6300,0.050236849,125405.954,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p2048,2048,DIVIDE,6300,0.004802009,1311950.869,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p2048,2048,DIVIDE,6300,0.005987883,1052124.834,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p2048,2048,DIVIDE,6300,0.049089419,128337.230,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p2048,2048,DIVIDE,6300,0.050315039,125211.072,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p2048,2048,DIVIDE,6300,0.004775159,1319327.779,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p2048,2048,DIVIDE,6300,0.005961876,1056714.299,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,p2048,2048,ISQRT,195,0.000251822,774355.956,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,p2048,2048,ISQRT,195,0.000035511,5491205.237,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p2048,2048,ISQRT,6300,16.203781903,388.798,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p2048,2048,ISQRT,6300,16.196858607,388.964,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p2048,2048,ISQRT,6300,9.255594144,680.669,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p2048,2048,ISQRT,6300,9.256136701,680.630,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p2048,2048,ISQRT,6300,0.318232970,19796.817,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p2048,2048,ISQRT,6300,0.319179272,19738.124,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p2048,2048,ISQRT,6300,0.088166190,71455.963,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p2048,2048,ISQRT,6300,0.089252943,70585.908,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p2048,2048,ISQRT,6300,0.318120372,19803.824,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p2048,2048,ISQRT,6300,0.319115022,19742.098,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p2048,2048,ISQRT,6300,0.087963967,71620.235,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p2048,2048,ISQRT,6300,0.089047985,70748.372,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-1t,p2048,2048,MODMUL_R2,6250,0.013651818,457814.482,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,gmp-nt,p2048,2048,MODMUL_R2,6250,0.001809094,3454767.773,0
library,Intel(R) Xeon(R) CPU E5-2697A v4 @ 2.60GHz,host-cpu,openssl-nt,p2048,2048,MODMUL_R2,6250,0.004906500,1273820.505,0
library,NVIDIA GeForce GTX 1060 3GB,gpu,cgbn,p2048,2048,MODMUL_R2,50000,0.004216832,11857242.593,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p2048,2048,MODMUL_R2,6300,0.043169428,145936.612,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w8,p2048,2048,MODMUL_R2,6300,0.044410888,141857.105,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p2048,2048,MODMUL_R2,6300,0.006070483,1037808.607,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w16,p2048,2048,MODMUL_R2,6300,0.007110897,885964.120,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p2048,2048,MODMUL_R2,6300,0.001712343,3679170.011,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-opt,p2048,2048,MODMUL_R2,6300,0.002785655,2261586.805,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p2048,2048,MODMUL_R2,6300,0.001346827,4677662.454,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-o64,p2048,2048,MODMUL_R2,6300,0.002179878,2890070.217,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p2048,2048,MODMUL_R2,6300,0.001895444,3323758.679,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il,p2048,2048,MODMUL_R2,6300,0.002914924,2161291.211,0
opencl-kernel,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p2048,2048,MODMUL_R2,6300,0.001387421,4540798.943,0
opencl-e2e,NVIDIA GeForce GTX 1060 3GB,GPU,w32-il64,p2048,2048,MODMUL_R2,6300,0.002365118,2663714.997,0
```
