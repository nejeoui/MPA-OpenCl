# MPA-OpenCL benchmark report - NVIDIA RTX A2000


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

### Device 0 - NVIDIA RTX A2000 (GPU)

| Property | Value |
|---|---|
| Model | NVIDIA RTX A2000 |
| Type | GPU |
| Vendor | NVIDIA Corporation |
| Device memory | 5.66 GiB |
| Max single allocation | 1.42 GiB |
| Local memory | 48 KiB |
| Global cache | 728 KiB |
| Compute units | 26 |
| Max clock | 1200 MHz |
| Max work-group size | 1024 |
| OpenCL version | OpenCL 3.0 CUDA |
| Driver | 590.48.01 |

### Device 1 - cpu-haswell-AMD EPYC 7282 16-Core Processor (CPU)

| Property | Value |
|---|---|
| Model | cpu-haswell-AMD EPYC 7282 16-Core Processor |
| Type | CPU |
| Vendor | AuthenticAMD |
| Device memory | 249.57 GiB |
| Max single allocation | 64.00 GiB |
| Local memory | 512 KiB |
| Global cache | 16384 KiB |
| Compute units | 32 |
| Max clock | 2800 MHz |
| Max work-group size | 4096 |
| OpenCL version | OpenCL 3.0 PoCL HSTR: cpu-x86_64-pc-linux-gnu-haswell |
| Driver | 5.0+debian |

### Host

| Property | Value |
|---|---|
| CPU | AMD EPYC 7282 16-Core Processor |
| Logical cores | 32 |
| OpenMP threads used | 32 |
| RAM | 251.6 GB |
| OS | Ubuntu 24.04.4 LTS |
| Kernel | 5.15.0-181-generic |
| Arch | x86_64 |
| GMP | 6.3.0 |
| OpenSSL | OpenSSL 3.0.13 30 Jan 2024 |

## 2. Method

- Workload auto-sized from the device and host: --min-items from 700 x compute units, --items from ten times that capped by host RAM. Either flag, given explicitly, overrides its half.
- Base workload 50000 items, scaled down per operator by its cost weight and by modulus size. Device rows honour --min-items (18200) so the GPU is not left idle; the CPU libraries keep the smaller count because a full-width MODEXP there costs minutes. Both counts appear in every row as dev/cpu, and throughput is per-second so they remain comparable.
- 5 timed repetitions, **minimum** reported. Two untimed warm-up launches precede them.
- `kernel` times `clEnqueueNDRangeKernel` + `clFinish` only. `e2e` adds the host->device operand writes and the device->host result read.
- Every OpenCL device runs the same kernels on the same operands, so GPU and CPU-OpenCL columns are directly comparable.
- CPU library baselines (GMP, OpenSSL) run those same operands, with temporaries preallocated outside the timed region, so the figure is the arithmetic and not marshalling. The generator is reseeded per modulus and operation so every backend sees identical inputs.
- OpenSSL rows time the nearest BN primitive, which is not always semantically identical (its Montgomery routine expects Montgomery-domain inputs); they measure comparable work, not identical results. Correctness is judged against GMP only.
- Every device cell is checked word-for-word against GMP before it is timed. A cell that mismatches is reported and excluded from the speedup tables.
- Total wall time 6317.4 s.

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

### Device 0 - NVIDIA RTX A2000 (GPU)

#### secp256k1 (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | 325.32 M | 806.20 M | 1.15 G | 1.34 G | 1.28 G | 1.78 G | 1.83 G | 56.08 M |
| SUBTRACT | 50000 / 50000 | 324.28 M | 874.90 M | 1.18 G | 1.34 G | 1.35 G | 1.77 G | 1.83 G | 71.80 M |
| ADDMOD | 50000 / 50000 | 212.89 M | 615.18 M | 996.68 M | 1.77 G | 1.79 G | 1.93 G | 1.97 G | 18.34 M |
| SUBTRACTMOD | 50000 / 50000 | 213.29 M | 613.60 M | 993.87 M | 1.77 G | 1.78 G | 1.93 G | 1.96 G | 19.63 M |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | 8.46 M | 45.03 M | 162.07 M | 782.55 M | 743.27 M | 1.44 G | 1.46 G | 44.71 M |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | 47.96 M | 222.52 M | 511.99 M | 641.02 M | 663.43 M | 937.32 M | 986.42 M | 44.72 M |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | 101.68 M | 640.52 M | 1.72 G | 1.30 G | 1.28 G | 1.28 G | 1.34 G | 5.48 M |
| COMPARE | 50000 / 50000 | 346.20 M | 941.40 M | - | 1.75 G | 1.79 G | 4.12 G | 4.20 G | 129.12 M |
| REDUCE | 18200 / 6250 | 42.47 M | 106.48 M | - | 311.42 M | 325.01 M | 300.69 M | 331.66 M | 44.56 M |
| MODMUL | 18200 / 3125 | 14.65 M | 43.02 M | - | 98.20 M | 122.18 M | 99.25 M | 122.97 M | 9.30 M |
| MODEXP | 18200 / 781 | 449.68 k | 2.45 M | - | 4.54 M | 7.90 M | 4.54 M | 8.07 M | 93.28 k |
| EXPONENTIATION | 18200 / 781 | 344.81 k | 1.28 M | - | 22.46 M | 32.61 M | 22.53 M | 32.52 M | 280.12 k |
| DIVIDE | 18200 / 6250 | 33.17 M | 49.11 M | - | 145.16 M | 145.58 M | 147.99 M | 149.11 M | 22.02 M |
| ISQRT | 18200 / 1562 | 2.45 M | 3.42 M | - | 15.21 M | 17.80 M | 17.23 M | 18.90 M | 9.63 M |
| MODMUL_R2 | 50000 / 50000 | 136.45 M | 655.33 M | - | 892.08 M | 884.03 M | 908.50 M | 952.15 M | 9.58 M |

#### rsa256(composite) (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | 496.41 M | 880.92 M | 1.19 G | 1.35 G | 1.34 G | 1.76 G | 1.78 G | 55.19 M |
| SUBTRACT | 50000 / 50000 | 491.56 M | 865.93 M | 1.18 G | 1.07 G | 1.36 G | 1.78 G | 1.77 G | 72.04 M |
| ADDMOD | 50000 / 50000 | 356.45 M | 662.73 M | 1.07 G | 1.47 G | 1.79 G | 1.93 G | 2.01 G | 20.33 M |
| SUBTRACTMOD | 50000 / 50000 | 324.35 M | 611.87 M | 985.84 M | 1.52 G | 1.84 G | 1.97 G | 1.97 G | 19.42 M |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | 13.14 M | 45.04 M | 161.82 M | 539.03 M | 746.84 M | 1.45 G | 1.44 G | 44.57 M |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | 74.25 M | 222.26 M | 517.03 M | 462.67 M | 666.19 M | 941.22 M | 989.17 M | 44.54 M |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | 157.73 M | 636.44 M | 1.73 G | 1.17 G | 1.26 G | 1.28 G | 1.35 G | 5.48 M |
| COMPARE | 50000 / 50000 | 532.59 M | 959.88 M | - | 1.46 G | 1.78 G | 4.17 G | 4.19 G | 129.13 M |
| REDUCE | 18200 / 6250 | 64.88 M | 106.19 M | - | 212.36 M | 323.04 M | 298.77 M | 328.90 M | 27.27 M |
| MODMUL | 18200 / 3125 | 22.67 M | 42.88 M | - | 66.68 M | 122.17 M | 99.00 M | 122.74 M | 9.61 M |
| MODEXP | 18200 / 781 | 453.62 k | 2.69 M | - | 3.04 M | 5.25 M | 3.04 M | 5.42 M | 98.95 k |
| EXPONENTIATION | 18200 / 781 | 341.95 k | 1.28 M | - | 15.05 M | 21.74 M | 15.11 M | 21.73 M | 281.10 k |
| DIVIDE | 18200 / 6250 | 32.65 M | 50.12 M | - | 100.11 M | 99.40 M | 99.05 M | 102.84 M | 21.98 M |
| ISQRT | 18200 / 1562 | 2.43 M | 3.43 M | - | 9.97 M | 11.63 M | 11.34 M | 12.43 M | 10.68 M |
| MODMUL_R2 | 50000 / 50000 | 135.59 M | 643.84 M | - | 757.83 M | 840.54 M | 743.94 M | 887.65 M | 9.55 M |

#### brainpoolP512r1 (512-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 25000 / 25000 | 207.23 M | 359.69 M | 316.52 M | 438.74 M | 430.88 M | 758.66 M | 553.02 M | 51.24 M |
| SUBTRACT | 25000 / 25000 | 207.76 M | 364.10 M | 312.56 M | 445.09 M | 446.44 M | 753.38 M | 554.74 M | 65.22 M |
| ADDMOD | 25000 / 25000 | 153.58 M | 276.90 M | 287.33 M | 372.19 M | 370.31 M | 716.59 M | 740.63 M | 19.02 M |
| SUBTRACTMOD | 25000 / 25000 | 133.12 M | 249.72 M | 266.58 M | 374.82 M | 372.41 M | 713.11 M | 737.78 M | 18.59 M |
| MULTIPLYOPERANDSCANNING | 25000 / 25000 | 2.22 M | 8.77 M | 12.73 M | 126.24 M | 125.69 M | 253.34 M | 253.37 M | 18.59 M |
| MULTIPLYPRODUCTSCANNING | 25000 / 25000 | 10.43 M | 37.59 M | 74.31 M | 88.16 M | 87.83 M | 191.17 M | 136.60 M | 18.60 M |
| MONTGOMERYMULTIPLICATION | 25000 / 25000 | 42.44 M | 136.71 M | 515.85 M | 277.18 M | 383.17 M | 314.38 M | 428.14 M | 2.35 M |
| COMPARE | 25000 / 25000 | 209.87 M | 373.26 M | - | 488.86 M | 487.62 M | 1.54 G | 1.56 G | 129.33 M |
| REDUCE | 18200 / 3125 | 21.20 M | 26.44 M | - | 72.50 M | 77.71 M | 77.07 M | 61.87 M | 26.51 M |
| MODMUL | 18200 / 1562 | 7.16 M | 10.17 M | - | 19.03 M | 23.32 M | 18.73 M | 18.22 M | 4.82 M |
| MODEXP | 18200 / 390 | 37.15 k | 305.92 k | - | 405.30 k | 859.49 k | 412.15 k | 788.17 k | 18.35 k |
| EXPONENTIATION | 18200 / 390 | 42.77 k | 156.03 k | - | 501.50 k | 562.48 k | 506.40 k | 560.36 k | 80.84 k |
| DIVIDE | 18200 / 3125 | 11.13 M | 13.63 M | - | 38.82 M | 38.31 M | 38.46 M | 45.05 M | 18.63 M |
| ISQRT | 18200 / 781 | 483.50 k | 460.17 k | - | 2.48 M | 2.79 M | 2.67 M | 2.09 M | 6.25 M |
| MODMUL_R2 | 25000 / 25000 | 28.85 M | 144.77 M | - | 213.96 M | 320.18 M | 224.22 M | 331.98 M | 4.79 M |

#### p1024 (1024-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 18200 / 12500 | 64.41 M | 120.18 M | 126.70 M | 239.20 M | 234.12 M | 285.25 M | 283.25 M | 43.14 M |
| SUBTRACT | 18200 / 12500 | 64.59 M | 119.99 M | 126.67 M | 238.76 M | 231.78 M | 286.47 M | 283.82 M | 54.10 M |
| ADDMOD | 18200 / 12500 | 43.97 M | 83.46 M | 100.04 M | 184.34 M | 180.98 M | 421.08 M | 422.73 M | 14.58 M |
| SUBTRACTMOD | 18200 / 12500 | 44.30 M | 84.85 M | 98.93 M | 182.43 M | 183.52 M | 419.03 M | 414.25 M | 16.29 M |
| MULTIPLYOPERANDSCANNING | 18200 / 12500 | 437.41 k | 1.71 M | 3.07 M | 68.37 M | 69.51 M | 99.17 M | 110.31 M | 5.30 M |
| MULTIPLYPRODUCTSCANNING | 18200 / 12500 | 1.25 M | 4.84 M | 12.43 M | 18.98 M | 18.96 M | 38.06 M | 40.73 M | 5.30 M |
| MONTGOMERYMULTIPLICATION | 18200 / 12500 | 6.98 M | 31.71 M | 90.14 M | 80.64 M | 104.58 M | 100.93 M | 145.48 M | 776.05 k |
| COMPARE | 18200 / 12500 | 74.97 M | 140.74 M | - | 264.18 M | 254.59 M | 1.20 G | 1.19 G | 134.01 M |
| REDUCE | 18200 / 1562 | 3.89 M | 7.13 M | - | 27.53 M | 26.87 M | 32.58 M | 29.96 M | 33.79 M |
| MODMUL | 18200 / 781 | 1.01 M | 2.63 M | - | 6.01 M | 8.09 M | 6.17 M | 8.22 M | 1.66 M |
| MODEXP | 18200 / 195 | 3.97 k | 24.46 k | - | 50.12 k | 82.59 k | 50.09 k | 81.06 k | 2.78 k |
| EXPONENTIATION | 18200 / 195 | 4.95 k | 19.70 k | - | 60.48 k | 68.31 k | 60.23 k | 68.37 k | 18.99 k |
| DIVIDE | 18200 / 1562 | 342.31 k | 2.47 M | - | 10.83 M | 11.05 M | 10.85 M | 10.96 M | 17.39 M |
| ISQRT | 18200 / 390 | 19.75 k | 63.83 k | - | 328.02 k | 347.46 k | 326.74 k | 353.61 k | 2.77 M |
| MODMUL_R2 | 18200 / 12500 | 3.64 M | 27.49 M | - | 49.48 M | 63.27 M | 54.95 M | 72.41 M | 1.74 M |

#### p2048 (2048-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 18200 / 6250 | 31.87 M | 63.22 M | 88.60 M | 131.29 M | 130.95 M | 154.91 M | 158.05 M | 31.74 M |
| SUBTRACT | 18200 / 6250 | 31.99 M | 62.92 M | 88.79 M | 132.22 M | 131.11 M | 154.69 M | 157.80 M | 38.25 M |
| ADDMOD | 18200 / 6250 | 24.62 M | 47.27 M | 70.16 M | 97.82 M | 95.71 M | 209.63 M | 207.58 M | 11.25 M |
| SUBTRACTMOD | 18200 / 6250 | 23.03 M | 43.48 M | 62.30 M | 94.47 M | 94.66 M | 196.93 M | 209.56 M | 13.32 M |
| MULTIPLYOPERANDSCANNING | 18200 / 6250 | 99.91 k | 361.72 k | 1.27 M | 23.31 M | 24.22 M | 26.63 M | 31.95 M | 1.55 M |
| MULTIPLYPRODUCTSCANNING | 18200 / 6250 | 313.07 k | 1.20 M | 4.87 M | 4.86 M | 4.88 M | 8.37 M | 8.97 M | 1.55 M |
| MONTGOMERYMULTIPLICATION | 18200 / 6250 | 161.18 k | 6.52 M | 37.91 M | 25.12 M | 30.19 M | 26.56 M | 31.86 M | 238.13 k |
| COMPARE | 18200 / 6250 | 37.28 M | 74.00 M | - | 141.66 M | 94.58 M | 349.38 M | 359.27 M | 128.33 M |
| REDUCE | 18200 / 781 | 26.83 k | 1.51 M | - | 5.87 M | 5.20 M | 6.49 M | 5.47 M | 24.42 M |
| MODMUL | 18200 / 390 | 18.67 k | 482.31 k | - | 1.26 M | 1.34 M | 1.46 M | 1.30 M | 565.83 k |
| MODEXP | 18200 / 97 | 112.8 | 958.3 | - | 6.12 k | 4.03 k | 6.58 k | 4.03 k | 380.0 |
| EXPONENTIATION | 18200 / 97 | 388.5 | 2.08 k | - | 7.34 k | 8.07 k | 7.30 k | 8.01 k | 3.18 k |
| DIVIDE | 18200 / 781 | 10.35 k | 53.76 k | - | 987.69 k | 1.09 M | 1.02 M | 1.14 M | 13.60 M |
| ISQRT | 18200 / 195 | 783.8 | 2.05 k | - | 51.93 k | 159.70 k | 51.88 k | 160.02 k | 1.78 M |
| MODMUL_R2 | 18200 / 6250 | 244.74 k | 5.63 M | - | 13.54 M | 17.58 M | 14.59 M | 18.17 M | 553.69 k |

### Device 1 - cpu-haswell-AMD EPYC 7282 16-Core Processor (CPU)

#### secp256k1 (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | 122.72 M | - | - | - | - | - | - | 56.08 M |
| SUBTRACT | 50000 / 50000 | 96.44 M | - | - | - | - | - | - | 71.80 M |
| ADDMOD | 50000 / 50000 | 91.66 M | - | - | - | - | - | - | 18.34 M |
| SUBTRACTMOD | 50000 / 50000 | 86.16 M | - | - | - | - | - | - | 19.63 M |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | 5.82 M | - | - | - | - | - | - | 44.71 M |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | 8.87 M | - | - | - | - | - | - | 44.72 M |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | 2.12 M | - | - | - | - | - | - | 5.48 M |
| COMPARE | 50000 / 50000 | 139.68 M | - | - | - | - | - | - | 129.12 M |
| REDUCE | 18200 / 6250 | 2.36 M | - | - | - | - | - | - | 44.56 M |
| MODMUL | 18200 / 3125 | 574.60 k | - | - | - | - | - | - | 9.30 M |
| MODEXP | 18200 / 781 | 12.72 k | - | - | - | - | - | - | 93.28 k |
| EXPONENTIATION | 18200 / 781 | 30.26 k | - | - | - | - | - | - | 280.12 k |
| DIVIDE | 18200 / 6250 | 1.19 M | - | - | - | - | - | - | 22.02 M |
| ISQRT | 18200 / 1562 | 168.56 k | - | - | - | - | - | - | 9.63 M |
| MODMUL_R2 | 50000 / 50000 | 1.46 M | - | - | - | - | - | - | 9.58 M |

#### rsa256(composite) (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | 127.83 M | - | - | - | - | - | - | 55.19 M |
| SUBTRACT | 50000 / 50000 | 93.68 M | - | - | - | - | - | - | 72.04 M |
| ADDMOD | 50000 / 50000 | 97.23 M | - | - | - | - | - | - | 20.33 M |
| SUBTRACTMOD | 50000 / 50000 | 104.31 M | - | - | - | - | - | - | 19.42 M |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | 5.96 M | - | - | - | - | - | - | 44.57 M |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | 8.93 M | - | - | - | - | - | - | 44.54 M |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | 1.73 M | - | - | - | - | - | - | 5.48 M |
| COMPARE | 50000 / 50000 | 160.32 M | - | - | - | - | - | - | 129.13 M |
| REDUCE | 18200 / 6250 | 2.33 M | - | - | - | - | - | - | 27.27 M |
| MODMUL | 18200 / 3125 | 578.91 k | - | - | - | - | - | - | 9.61 M |
| MODEXP | 18200 / 781 | 12.57 k | - | - | - | - | - | - | 98.95 k |
| EXPONENTIATION | 18200 / 781 | 30.40 k | - | - | - | - | - | - | 281.10 k |
| DIVIDE | 18200 / 6250 | 1.40 M | - | - | - | - | - | - | 21.98 M |
| ISQRT | 18200 / 1562 | 208.71 k | - | - | - | - | - | - | 10.68 M |
| MODMUL_R2 | 50000 / 50000 | 1.48 M | - | - | - | - | - | - | 9.55 M |

#### brainpoolP512r1 (512-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 25000 / 25000 | 49.07 M | - | - | - | - | - | - | 51.24 M |
| SUBTRACT | 25000 / 25000 | 55.65 M | - | - | - | - | - | - | 65.22 M |
| ADDMOD | 25000 / 25000 | 50.97 M | - | - | - | - | - | - | 19.02 M |
| SUBTRACTMOD | 25000 / 25000 | 44.31 M | - | - | - | - | - | - | 18.59 M |
| MULTIPLYOPERANDSCANNING | 25000 / 25000 | 1.56 M | - | - | - | - | - | - | 18.59 M |
| MULTIPLYPRODUCTSCANNING | 25000 / 25000 | 1.74 M | - | - | - | - | - | - | 18.60 M |
| MONTGOMERYMULTIPLICATION | 25000 / 25000 | 558.23 k | - | - | - | - | - | - | 2.35 M |
| COMPARE | 25000 / 25000 | 76.12 M | - | - | - | - | - | - | 129.33 M |
| REDUCE | 18200 / 3125 | 624.78 k | - | - | - | - | - | - | 26.51 M |
| MODMUL | 18200 / 1562 | 256.52 k | - | - | - | - | - | - | 4.82 M |
| MODEXP | 18200 / 390 | 1.49 k | - | - | - | - | - | - | 18.35 k |
| EXPONENTIATION | 18200 / 390 | 2.76 k | - | - | - | - | - | - | 80.84 k |
| DIVIDE | 18200 / 3125 | 378.86 k | - | - | - | - | - | - | 18.63 M |
| ISQRT | 18200 / 781 | 36.09 k | - | - | - | - | - | - | 6.25 M |
| MODMUL_R2 | 25000 / 25000 | 490.80 k | - | - | - | - | - | - | 4.79 M |

#### p1024 (1024-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 18200 / 12500 | 30.07 M | - | - | - | - | - | - | 43.14 M |
| SUBTRACT | 18200 / 12500 | 25.59 M | - | - | - | - | - | - | 54.10 M |
| ADDMOD | 18200 / 12500 | 27.11 M | - | - | - | - | - | - | 14.58 M |
| SUBTRACTMOD | 18200 / 12500 | 24.67 M | - | - | - | - | - | - | 16.29 M |
| MULTIPLYOPERANDSCANNING | 18200 / 12500 | 371.27 k | - | - | - | - | - | - | 5.30 M |
| MULTIPLYPRODUCTSCANNING | 18200 / 12500 | 541.22 k | - | - | - | - | - | - | 5.30 M |
| MONTGOMERYMULTIPLICATION | 18200 / 12500 | 191.49 k | - | - | - | - | - | - | 776.05 k |
| COMPARE | 18200 / 12500 | 48.93 M | - | - | - | - | - | - | 134.01 M |
| REDUCE | 18200 / 1562 | 300.71 k | - | - | - | - | - | - | 33.79 M |
| MODMUL | 18200 / 781 | 70.43 k | - | - | - | - | - | - | 1.66 M |
| MODEXP | 18200 / 195 | over budget | - | - | - | - | - | - | 2.78 k |
| EXPONENTIATION | 18200 / 195 | over budget | - | - | - | - | - | - | 18.99 k |
| DIVIDE | 18200 / 1562 | 175.79 k | - | - | - | - | - | - | 17.39 M |
| ISQRT | 18200 / 390 | 8.09 k | - | - | - | - | - | - | 2.77 M |
| MODMUL_R2 | 18200 / 12500 | 211.65 k | - | - | - | - | - | - | 1.74 M |

#### p2048 (2048-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 18200 / 6250 | 15.97 M | - | - | - | - | - | - | 31.74 M |
| SUBTRACT | 18200 / 6250 | 15.39 M | - | - | - | - | - | - | 38.25 M |
| ADDMOD | 18200 / 6250 | 13.61 M | - | - | - | - | - | - | 11.25 M |
| SUBTRACTMOD | 18200 / 6250 | 12.35 M | - | - | - | - | - | - | 13.32 M |
| MULTIPLYOPERANDSCANNING | 18200 / 6250 | 142.22 k | - | - | - | - | - | - | 1.55 M |
| MULTIPLYPRODUCTSCANNING | 18200 / 6250 | 188.09 k | - | - | - | - | - | - | 1.55 M |
| MONTGOMERYMULTIPLICATION | 18200 / 6250 | 37.95 k | - | - | - | - | - | - | 238.13 k |
| COMPARE | 18200 / 6250 | 46.88 M | - | - | - | - | - | - | 128.33 M |
| REDUCE | 18200 / 781 | 81.09 k | - | - | - | - | - | - | 24.42 M |
| MODMUL | 18200 / 390 | 15.91 k | - | - | - | - | - | - | 565.83 k |
| MODEXP | 18200 / 97 | over budget | - | - | - | - | - | - | 380.0 |
| EXPONENTIATION | 18200 / 97 | - | - | - | - | - | - | - | 3.18 k |
| DIVIDE | 18200 / 781 | - | - | - | - | - | - | - | 13.60 M |
| ISQRT | 18200 / 195 | - | - | - | - | - | - | - | 1.78 M |
| MODMUL_R2 | 18200 / 6250 | - | - | - | - | - | - | - | 553.69 k |

## 5. Head to head

Best OpenCL GPU result against best OpenCL CPU result and the CPU libraries.
Ratios above 1.00x mean the GPU is faster than that baseline.

### secp256k1 (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w32-il64 | 1.83 G | w8 | 122.72 M | 56.08 M | 14.88x |
| SUBTRACT | w32-il64 | 1.83 G | w8 | 96.44 M | 71.80 M | 19.02x |
| ADDMOD | w32-il64 | 1.97 G | w8 | 91.66 M | 18.34 M | 21.47x |
| SUBTRACTMOD | w32-il64 | 1.96 G | w8 | 86.16 M | 19.63 M | 22.70x |
| MULTIPLYOPERANDSCANNING | w32-il64 | 1.46 G | w8 | 5.82 M | 44.71 M | 250.31x |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 986.42 M | w8 | 8.87 M | 44.72 M | 111.26x |
| MONTGOMERYMULTIPLICATION | w32 | 1.72 G | w8 | 2.12 M | 5.48 M | 812.15x |
| COMPARE | w32-il64 | 4.20 G | w8 | 139.68 M | 129.12 M | 30.05x |
| REDUCE | w32-il64 | 113.89 M | w8 | 811.86 k | 44.56 M | 140.29x |
| MODMUL | w32-il64 | 21.11 M | w8 | 98.66 k | 9.30 M | 214.01x |
| MODEXP | w32-il64 | 346.45 k | w8 | 545.7 | 93.28 k | 634.86x |
| EXPONENTIATION | w32-o64 | 1.40 M | w8 | 1.30 k | 280.12 k | 1077.70x |
| DIVIDE | w32-il64 | 51.21 M | w8 | 407.43 k | 22.02 M | 125.68x |
| ISQRT | w32-il64 | 1.62 M | w8 | 14.47 k | 9.63 M | 112.14x |
| MODMUL_R2 | w32-il64 | 952.15 M | w8 | 1.46 M | 9.58 M | 652.35x |

### rsa256(composite) (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w32-il64 | 1.78 G | w8 | 127.83 M | 55.19 M | 13.94x |
| SUBTRACT | w32-il | 1.78 G | w8 | 93.68 M | 72.04 M | 19.04x |
| ADDMOD | w32-il64 | 2.01 G | w8 | 97.23 M | 20.33 M | 20.66x |
| SUBTRACTMOD | w32-il | 1.97 G | w8 | 104.31 M | 19.42 M | 18.92x |
| MULTIPLYOPERANDSCANNING | w32-il | 1.45 G | w8 | 5.96 M | 44.57 M | 242.62x |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 989.17 M | w8 | 8.93 M | 44.54 M | 110.75x |
| MONTGOMERYMULTIPLICATION | w32 | 1.73 G | w8 | 1.73 M | 5.48 M | 1000.33x |
| COMPARE | w32-il64 | 4.19 G | w8 | 160.32 M | 129.13 M | 26.14x |
| REDUCE | w32-il64 | 112.95 M | w8 | 801.21 k | 27.27 M | 140.97x |
| MODMUL | w32-il64 | 21.07 M | w8 | 99.40 k | 9.61 M | 212.01x |
| MODEXP | w32-il64 | 232.62 k | w8 | 539.3 | 98.95 k | 431.34x |
| EXPONENTIATION | w32-o64 | 932.77 k | w8 | 1.30 k | 281.10 k | 715.03x |
| DIVIDE | w32-il64 | 35.32 M | w8 | 479.47 k | 21.98 M | 73.66x |
| ISQRT | w32-il64 | 1.07 M | w8 | 17.91 k | 10.68 M | 59.58x |
| MODMUL_R2 | w32-il64 | 887.65 M | w8 | 1.48 M | 9.55 M | 598.96x |

### brainpoolP512r1 (512-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w32-il | 758.66 M | w8 | 49.07 M | 51.24 M | 15.46x |
| SUBTRACT | w32-il | 753.38 M | w8 | 55.65 M | 65.22 M | 13.54x |
| ADDMOD | w32-il64 | 740.63 M | w8 | 50.97 M | 19.02 M | 14.53x |
| SUBTRACTMOD | w32-il64 | 737.78 M | w8 | 44.31 M | 18.59 M | 16.65x |
| MULTIPLYOPERANDSCANNING | w32-il64 | 253.37 M | w8 | 1.56 M | 18.59 M | 162.02x |
| MULTIPLYPRODUCTSCANNING | w32-il | 191.17 M | w8 | 1.74 M | 18.60 M | 109.68x |
| MONTGOMERYMULTIPLICATION | w32 | 515.85 M | w8 | 558.23 k | 2.35 M | 924.09x |
| COMPARE | w32-il64 | 1.56 G | w8 | 76.12 M | 129.33 M | 20.50x |
| REDUCE | w32-o64 | 13.34 M | w8 | 107.28 k | 26.51 M | 124.38x |
| MODMUL | w32-o64 | 2.00 M | w8 | 22.02 k | 4.82 M | 90.91x |
| MODEXP | w32-o64 | 18.42 k | w8 | 31.9 | 18.35 k | 576.52x |
| EXPONENTIATION | w32-o64 | 12.05 k | w8 | 59.1 | 80.84 k | 203.92x |
| DIVIDE | w32-il64 | 7.74 M | w8 | 65.05 k | 18.63 M | 118.92x |
| ISQRT | w32-o64 | 119.75 k | w8 | 1.55 k | 6.25 M | 77.33x |
| MODMUL_R2 | w32-il64 | 331.98 M | w8 | 490.80 k | 4.79 M | 676.41x |

### p1024 (1024-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w32-il | 195.92 M | w8 | 20.65 M | 43.14 M | 9.49x |
| SUBTRACT | w32-il | 196.75 M | w8 | 17.57 M | 54.10 M | 11.19x |
| ADDMOD | w32-il64 | 290.34 M | w8 | 18.62 M | 14.58 M | 15.59x |
| SUBTRACTMOD | w32-il | 287.79 M | w8 | 16.94 M | 16.29 M | 16.99x |
| MULTIPLYOPERANDSCANNING | w32-il64 | 75.76 M | w8 | 254.99 k | 5.30 M | 297.12x |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 27.97 M | w8 | 371.71 k | 5.30 M | 75.26x |
| MONTGOMERYMULTIPLICATION | w32-il64 | 99.92 M | w8 | 131.52 k | 776.05 k | 759.74x |
| COMPARE | w32-il | 821.26 M | w8 | 33.61 M | 134.01 M | 24.44x |
| REDUCE | w32-il | 2.80 M | w8 | 25.81 k | 33.79 M | 108.36x |
| MODMUL | w32-il64 | 352.85 k | w8 | 3.02 k | 1.66 M | 116.75x |
| MODEXP | w32-o64 | 884.9 | none | n/a | 2.78 k | n/a |
| EXPONENTIATION | w32-il64 | 732.6 | none | n/a | 18.99 k | n/a |
| DIVIDE | w32-o64 | 948.23 k | w8 | 15.09 k | 17.39 M | 62.85x |
| ISQRT | w32-il64 | 7.58 k | w8 | 173.3 | 2.77 M | 43.72x |
| MODMUL_R2 | w32-il64 | 49.73 M | w8 | 145.37 k | 1.74 M | 342.13x |

### p2048 (2048-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w32-il64 | 54.28 M | w8 | 5.49 M | 31.74 M | 9.89x |
| SUBTRACT | w32-il64 | 54.19 M | w8 | 5.29 M | 38.25 M | 10.25x |
| ADDMOD | w32-il | 71.99 M | w8 | 4.68 M | 11.25 M | 15.40x |
| SUBTRACTMOD | w32-il64 | 71.96 M | w8 | 4.24 M | 13.32 M | 16.96x |
| MULTIPLYOPERANDSCANNING | w32-il64 | 10.97 M | w8 | 48.84 k | 1.55 M | 224.65x |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 3.08 M | w8 | 64.59 k | 1.55 M | 47.68x |
| MONTGOMERYMULTIPLICATION | w32 | 13.02 M | w8 | 13.03 k | 238.13 k | 998.95x |
| COMPARE | w32-il64 | 123.38 M | w8 | 16.10 M | 128.33 M | 7.66x |
| REDUCE | w32-il | 278.36 k | w8 | 3.48 k | 24.42 M | 79.99x |
| MODMUL | w32-il | 31.32 k | w8 | 340.9 | 565.83 k | 91.88x |
| MODEXP | w32-il | 35.1 | none | n/a | 380.0 | n/a |
| EXPONENTIATION | w32-o64 | 43.0 | none | n/a | 3.18 k | n/a |
| DIVIDE | w32-il64 | 48.74 k | none | n/a | 13.60 M | n/a |
| ISQRT | w32-il64 | 1.71 k | none | n/a | 1.78 M | n/a |
| MODMUL_R2 | w32-il64 | 6.24 M | none | n/a | 553.69 k | n/a |

## 6. Raw data

Also written to `NVIDIA_RTX_A2000_Report.csv` for analysis.

```csv
kind,device,device_type,kernel,modulus,bits,operation,items,seconds,ops_per_sec,mismatches
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,ADD,50000,0.000891583,56080025.905,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,ADD,50000,0.000123928,403459093.539,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,ADD,50000,0.000121223,412463631.475,0
library,NVIDIA RTX A2000,gpu,cgbn,secp256k1,256,ADD,50000,0.000034816,1436121323.529,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,secp256k1,256,ADD,50000,0.000153695,325319133.001,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,secp256k1,256,ADD,50000,0.000748379,66811053.728,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,secp256k1,256,ADD,50000,0.000062020,806197215.924,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,secp256k1,256,ADD,50000,0.000688154,72658128.569,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32,secp256k1,256,ADD,50000,0.000043305,1154610761.753,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32,secp256k1,256,ADD,50000,0.000674074,74175878.063,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,secp256k1,256,ADD,50000,0.000037221,1343319101.236,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,secp256k1,256,ADD,50000,0.000653904,76463831.460,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,secp256k1,256,ADD,50000,0.000039197,1275621716.920,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,secp256k1,256,ADD,50000,0.000658400,75941635.641,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,secp256k1,256,ADD,50000,0.000028074,1781020806.794,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,secp256k1,256,ADD,50000,0.000650737,76836015.417,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,secp256k1,256,ADD,50000,0.000027384,1825905220.556,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,secp256k1,256,ADD,50000,0.000698077,71625381.494,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,secp256k1,256,ADD,50000,0.000407446,122715595.257,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,secp256k1,256,ADD,50000,0.001751798,28542100.257,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,SUBTRACT,50000,0.000696379,71800007.222,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,SUBTRACT,50000,0.000104280,479477460.034,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,SUBTRACT,50000,0.000118237,422879511.012,0
library,NVIDIA RTX A2000,gpu,cgbn,secp256k1,256,SUBTRACT,50000,0.000034816,1436121323.529,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,secp256k1,256,SUBTRACT,50000,0.000154188,324279655.468,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,secp256k1,256,SUBTRACT,50000,0.000771482,64810364.426,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,secp256k1,256,SUBTRACT,50000,0.000057150,874895560.915,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,secp256k1,256,SUBTRACT,50000,0.000682062,73307168.352,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32,secp256k1,256,SUBTRACT,50000,0.000042542,1175312314.193,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32,secp256k1,256,SUBTRACT,50000,0.000661929,75536786.537,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,secp256k1,256,SUBTRACT,50000,0.000037341,1338997161.741,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,secp256k1,256,SUBTRACT,50000,0.000663934,75308659.071,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,secp256k1,256,SUBTRACT,50000,0.000036920,1354264086.976,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,secp256k1,256,SUBTRACT,50000,0.000634726,78774131.222,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,secp256k1,256,SUBTRACT,50000,0.000028265,1768990451.086,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,secp256k1,256,SUBTRACT,50000,0.000646798,77303894.922,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,secp256k1,256,SUBTRACT,50000,0.000027253,1834640713.529,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,secp256k1,256,SUBTRACT,50000,0.000691684,72287348.961,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,secp256k1,256,SUBTRACT,50000,0.000518450,96441393.186,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,secp256k1,256,SUBTRACT,50000,0.001742709,28690956.275,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,ADDMOD,50000,0.002725989,18341969.410,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,ADDMOD,50000,0.000319895,156301122.902,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,ADDMOD,50000,0.000831327,60144775.613,0
library,NVIDIA RTX A2000,gpu,cgbn,secp256k1,256,ADDMOD,50000,0.000033792,1479640151.515,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,secp256k1,256,ADDMOD,50000,0.000234861,212891947.022,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,secp256k1,256,ADDMOD,50000,0.000836347,59783781.448,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,secp256k1,256,ADDMOD,50000,0.000081277,615183811.161,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,secp256k1,256,ADDMOD,50000,0.000728320,68651111.213,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32,secp256k1,256,ADDMOD,50000,0.000050167,996678632.161,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32,secp256k1,256,ADDMOD,50000,0.000661858,75544864.599,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,secp256k1,256,ADDMOD,50000,0.000028255,1769573525.825,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,secp256k1,256,ADDMOD,50000,0.000652540,76623599.463,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,secp256k1,256,ADDMOD,50000,0.000027994,1786116548.007,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,secp256k1,256,ADDMOD,50000,0.000623103,80243525.466,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,secp256k1,256,ADDMOD,50000,0.000025890,1931259800.712,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,secp256k1,256,ADDMOD,50000,0.000649364,76998448.474,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,secp256k1,256,ADDMOD,50000,0.000025408,1967857605.747,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,secp256k1,256,ADDMOD,50000,0.000660804,75665390.052,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,secp256k1,256,ADDMOD,50000,0.000545522,91655298.677,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,secp256k1,256,ADDMOD,50000,0.001686271,29651219.582,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,SUBTRACTMOD,50000,0.002547456,19627424.753,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,SUBTRACTMOD,50000,0.000285639,175045862.610,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,SUBTRACTMOD,50000,0.000885782,56447308.120,0
library,NVIDIA RTX A2000,gpu,cgbn,secp256k1,256,SUBTRACTMOD,50000,0.000033792,1479640151.515,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.000234421,213291159.598,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.000826238,60515271.325,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,secp256k1,256,SUBTRACTMOD,50000,0.000081486,613601819.533,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,secp256k1,256,SUBTRACTMOD,50000,0.000709545,70467707.267,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32,secp256k1,256,SUBTRACTMOD,50000,0.000050308,993874101.226,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32,secp256k1,256,SUBTRACTMOD,50000,0.000647993,77161347.961,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,secp256k1,256,SUBTRACTMOD,50000,0.000028204,1772787320.037,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,secp256k1,256,SUBTRACTMOD,50000,0.000690344,72427681.312,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,secp256k1,256,SUBTRACTMOD,50000,0.000028124,1777835989.138,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,secp256k1,256,SUBTRACTMOD,50000,0.000631548,79170487.819,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,secp256k1,256,SUBTRACTMOD,50000,0.000025950,1926754636.807,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,secp256k1,256,SUBTRACTMOD,50000,0.000677407,73810893.093,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,secp256k1,256,SUBTRACTMOD,50000,0.000025569,1955457701.694,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,secp256k1,256,SUBTRACTMOD,50000,0.000692986,72151534.765,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.000580328,86158234.182,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.001805922,27686688.319,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001118229,44713569.625,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000148767,336096779.081,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000281329,177727687.495,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.005913565,8455136.715,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.006803735,7348904.510,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001110305,45032668.026,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001889404,26463372.787,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000308516,162065920.240,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001083786,46134583.134,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000063893,782553621.456,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000848118,58954045.637,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000067270,743269388.490,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000835475,59846225.326,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000034837,1435253467.358,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000799512,62538109.893,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000034336,1456197548.009,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000855180,58467193.542,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.008594811,5817463.609,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.010244683,4880580.320,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001118188,44715208.242,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000152153,328616669.829,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000282363,177077003.150,0
library,NVIDIA RTX A2000,gpu,cgbn,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000034816,1436121323.529,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001042513,47961024.487,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001819286,27483301.902,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000224702,222516873.902,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001210008,41322047.111,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000097658,511988281.518,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000908787,55018370.716,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000078001,641016933.125,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000885772,56447960.971,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000075365,663434267.143,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000827468,60425297.388,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000053343,937323728.547,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000851423,58725247.646,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000050688,986423606.365,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000878465,56917440.538,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.005639633,8865824.877,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.007343669,6808585.344,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.009123294,5480476.869,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000993408,50331773.829,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000302953,165042257.903,0
library,NVIDIA RTX A2000,gpu,cgbn,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000047104,1061480978.261,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000491725,101682793.104,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.001107738,45137012.590,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000078062,640519831.063,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000894129,55920325.312,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000028996,1724387846.085,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000626170,79850510.080,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000038505,1298546130.031,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000695063,71935945.803,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000039104,1278629398.876,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000641269,77970438.468,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000038994,1282232892.286,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000668019,74848164.176,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000037212,1343655300.831,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000642358,77838175.394,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.023549044,2123228.480,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.033526321,1491365.548,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,COMPARE,50000,0.000387234,129120999.353,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,COMPARE,50000,0.000059566,839411663.904,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,COMPARE,50000,0.000090674,551428627.773,0
library,NVIDIA RTX A2000,gpu,cgbn,secp256k1,256,COMPARE,50000,0.000034816,1436121323.529,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,secp256k1,256,COMPARE,50000,0.000144427,346196348.911,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,secp256k1,256,COMPARE,50000,0.000738890,67669079.391,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,secp256k1,256,COMPARE,50000,0.000053112,941399835.172,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,secp256k1,256,COMPARE,50000,0.000857639,58299578.991,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,secp256k1,256,COMPARE,50000,0.000028616,1747285399.987,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,secp256k1,256,COMPARE,50000,0.000690054,72458081.896,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,secp256k1,256,COMPARE,50000,0.000027934,1789927692.205,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,secp256k1,256,COMPARE,50000,0.000639915,78135320.417,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,secp256k1,256,COMPARE,50000,0.000012133,4120900460.547,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,secp256k1,256,COMPARE,50000,0.000651186,76782938.362,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,secp256k1,256,COMPARE,50000,0.000011913,4197255195.059,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,secp256k1,256,COMPARE,50000,0.000681355,73383221.136,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,secp256k1,256,COMPARE,50000,0.000357960,139680274.328,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,secp256k1,256,COMPARE,50000,0.001557102,32110925.484,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,REDUCE,6250,0.000140260,44560111.020,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,REDUCE,6250,0.000029347,212969642.347,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,REDUCE,6250,0.000112335,55637058.838,0
library,NVIDIA RTX A2000,gpu,cgbn,secp256k1,256,REDUCE,50000,0.000073728,678168402.778,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,secp256k1,256,REDUCE,18200,0.000428572,42466580.606,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,secp256k1,256,REDUCE,18200,0.000766222,23752895.800,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,secp256k1,256,REDUCE,18200,0.000170928,106477315.779,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,secp256k1,256,REDUCE,18200,0.000593311,30675304.007,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,secp256k1,256,REDUCE,18200,0.000058442,311417981.846,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,secp256k1,256,REDUCE,18200,0.000472323,38532950.011,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,secp256k1,256,REDUCE,18200,0.000055998,325013740.862,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,secp256k1,256,REDUCE,18200,0.000377129,48259370.416,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,secp256k1,256,REDUCE,18200,0.000060528,300689344.629,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,secp256k1,256,REDUCE,18200,0.000473743,38417429.487,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,secp256k1,256,REDUCE,18200,0.000054875,331660520.634,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,secp256k1,256,REDUCE,18200,0.000431602,42168490.096,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,secp256k1,256,REDUCE,18200,0.007698383,2364132.776,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,secp256k1,256,REDUCE,18200,0.008369719,2174505.584,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,MODMUL,3125,0.000335896,9303464.147,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,MODMUL,3125,0.000048113,64951185.614,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,MODMUL,3125,0.000129539,24124085.671,0
library,NVIDIA RTX A2000,gpu,cgbn,secp256k1,256,MODMUL,50000,0.000261120,191482843.137,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,secp256k1,256,MODMUL,18200,0.001242588,14646846.730,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,secp256k1,256,MODMUL,18200,0.001603422,11350725.373,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,secp256k1,256,MODMUL,18200,0.000423034,43022574.950,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,secp256k1,256,MODMUL,18200,0.000848302,21454631.404,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,secp256k1,256,MODMUL,18200,0.000185329,98203487.491,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,secp256k1,256,MODMUL,18200,0.000597236,30473728.509,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,secp256k1,256,MODMUL,18200,0.000148957,122183187.530,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,secp256k1,256,MODMUL,18200,0.000476060,38230514.994,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,secp256k1,256,MODMUL,18200,0.000183374,99250879.636,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,secp256k1,256,MODMUL,18200,0.000600019,30332396.136,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,secp256k1,256,MODMUL,18200,0.000148005,122968941.390,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,secp256k1,256,MODMUL,18200,0.000532467,34180527.234,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,secp256k1,256,MODMUL,18200,0.031674099,574601.977,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,secp256k1,256,MODMUL,18200,0.033415066,544664.492,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,MODEXP,781,0.008372962,93276.435,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,MODEXP,781,0.000981094,796049.868,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,MODEXP,781,0.001910752,408739.663,0
library,NVIDIA RTX A2000,gpu,cgbn,secp256k1,256,MODEXP,50000,0.062068928,805556.042,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,secp256k1,256,MODEXP,18200,0.040473090,449681.499,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,secp256k1,256,MODEXP,18200,0.040752018,446603.652,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,secp256k1,256,MODEXP,18200,0.007417054,2453804.463,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,secp256k1,256,MODEXP,18200,0.007984958,2279285.724,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,secp256k1,256,MODEXP,18200,0.004009158,4539606.174,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,secp256k1,256,MODEXP,18200,0.005020886,3624858.438,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,secp256k1,256,MODEXP,18200,0.002303498,7901027.507,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,secp256k1,256,MODEXP,18200,0.002699994,6740756.842,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,secp256k1,256,MODEXP,18200,0.004006061,4543116.302,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,secp256k1,256,MODEXP,18200,0.004485295,4057703.980,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,secp256k1,256,MODEXP,18200,0.002254268,8073573.662,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,secp256k1,256,MODEXP,18200,0.002643090,6885880.458,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,secp256k1,256,MODEXP,18200,1.431145079,12717.089,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,secp256k1,256,MODEXP,18200,1.430145247,12725.980,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,EXPONENTIATION,781,0.002788108,280118.290,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,EXPONENTIATION,781,0.000313583,2490570.034,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,EXPONENTIATION,781,0.005023023,155484.055,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,secp256k1,256,EXPONENTIATION,18200,0.052782248,344812.901,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,secp256k1,256,EXPONENTIATION,18200,0.053223775,341952.447,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,secp256k1,256,EXPONENTIATION,18200,0.014241720,1277935.564,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,secp256k1,256,EXPONENTIATION,18200,0.014678329,1239923.135,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,secp256k1,256,EXPONENTIATION,18200,0.000810185,22464018.554,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,secp256k1,256,EXPONENTIATION,18200,0.001198166,15189881.235,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,secp256k1,256,EXPONENTIATION,18200,0.000558178,32606089.703,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,secp256k1,256,EXPONENTIATION,18200,0.000915829,19872703.492,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,secp256k1,256,EXPONENTIATION,18200,0.000807929,22526736.034,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,secp256k1,256,EXPONENTIATION,18200,0.001189255,15303696.531,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,secp256k1,256,EXPONENTIATION,18200,0.000559620,32522090.171,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,secp256k1,256,EXPONENTIATION,18200,0.000951244,19132830.813,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,secp256k1,256,EXPONENTIATION,18200,0.601546217,30255.364,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,secp256k1,256,EXPONENTIATION,18200,0.611505903,29762.591,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,DIVIDE,6250,0.000283886,22015899.219,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,DIVIDE,6250,0.000046410,134670219.939,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,DIVIDE,6250,0.000117416,53229741.263,0
library,NVIDIA RTX A2000,gpu,cgbn,secp256k1,256,DIVIDE,50000,0.000098304,508626302.083,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,secp256k1,256,DIVIDE,18200,0.000548764,33165432.906,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,secp256k1,256,DIVIDE,18200,0.000973211,18700986.909,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,secp256k1,256,DIVIDE,18200,0.000370592,49110628.259,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,secp256k1,256,DIVIDE,18200,0.000892648,20388786.327,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,secp256k1,256,DIVIDE,18200,0.000125382,145156291.387,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,secp256k1,256,DIVIDE,18200,0.000611512,29762294.583,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,secp256k1,256,DIVIDE,18200,0.000125021,145575843.242,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,secp256k1,256,DIVIDE,18200,0.000618976,29403420.296,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,secp256k1,256,DIVIDE,18200,0.000122977,147994647.296,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,secp256k1,256,DIVIDE,18200,0.000616469,29522956.687,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,secp256k1,256,DIVIDE,18200,0.000122056,149111464.453,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,secp256k1,256,DIVIDE,18200,0.000528189,34457339.729,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,secp256k1,256,DIVIDE,18200,0.015340055,1186436.398,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,secp256k1,256,DIVIDE,18200,0.015258919,1192745.008,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,ISQRT,1562,0.000162222,9628755.226,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,ISQRT,1562,0.000027784,56219110.686,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,secp256k1,256,ISQRT,18200,0.007421643,2452287.026,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,secp256k1,256,ISQRT,18200,0.007970347,2283463.925,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,secp256k1,256,ISQRT,18200,0.005317517,3422650.324,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,secp256k1,256,ISQRT,18200,0.006238218,2917499.652,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,secp256k1,256,ISQRT,18200,0.001196383,15212525.239,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,secp256k1,256,ISQRT,18200,0.001587238,11466458.171,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,secp256k1,256,ISQRT,18200,0.001022643,17797014.901,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,secp256k1,256,ISQRT,18200,0.001406335,12941437.482,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,secp256k1,256,ISQRT,18200,0.001056129,17232741.507,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,secp256k1,256,ISQRT,18200,0.001443545,12607847.894,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,secp256k1,256,ISQRT,18200,0.000962848,18902259.706,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,secp256k1,256,ISQRT,18200,0.001352037,13461167.053,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,secp256k1,256,ISQRT,18200,0.107973239,168560.285,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,secp256k1,256,ISQRT,18200,0.108669010,167481.051,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,MODMUL_R2,50000,0.005217688,9582788.102,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,MODMUL_R2,50000,0.000287903,173669318.807,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,MODMUL_R2,50000,0.001225735,40791849.394,0
library,NVIDIA RTX A2000,gpu,cgbn,secp256k1,256,MODMUL_R2,50000,0.000053248,939002403.846,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,secp256k1,256,MODMUL_R2,50000,0.000366434,136450050.832,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,secp256k1,256,MODMUL_R2,50000,0.000966117,51753577.543,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,secp256k1,256,MODMUL_R2,50000,0.000076298,655328001.562,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,secp256k1,256,MODMUL_R2,50000,0.000927594,53902864.366,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,secp256k1,256,MODMUL_R2,50000,0.000056049,892078880.728,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,secp256k1,256,MODMUL_R2,50000,0.000675065,74066893.106,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,secp256k1,256,MODMUL_R2,50000,0.000056559,884029165.157,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,secp256k1,256,MODMUL_R2,50000,0.000717577,69678972.916,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,secp256k1,256,MODMUL_R2,50000,0.000055036,908503252.445,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,secp256k1,256,MODMUL_R2,50000,0.000667769,74876244.855,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,secp256k1,256,MODMUL_R2,50000,0.000052513,952152012.060,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,secp256k1,256,MODMUL_R2,50000,0.000709097,70512224.664,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,secp256k1,256,MODMUL_R2,50000,0.034256617,1459572.046,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,secp256k1,256,MODMUL_R2,50000,0.034493994,1449527.708,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,ADD,50000,0.000906001,55187586.489,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,ADD,50000,0.000102026,490069294.386,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,ADD,50000,0.000142124,351806579.120,0
library,NVIDIA RTX A2000,gpu,cgbn,rsa256(composite),256,ADD,50000,0.000033792,1479640151.515,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,rsa256(composite),256,ADD,50000,0.000100723,496408643.471,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,rsa256(composite),256,ADD,50000,0.000685517,72937678.584,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,rsa256(composite),256,ADD,50000,0.000056759,880924967.183,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,rsa256(composite),256,ADD,50000,0.000707351,70686297.517,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32,rsa256(composite),256,ADD,50000,0.000042011,1190163630.318,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32,rsa256(composite),256,ADD,50000,0.000641820,77903459.484,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,rsa256(composite),256,ADD,50000,0.000037093,1347973566.335,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,rsa256(composite),256,ADD,50000,0.000707437,70677643.249,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,rsa256(composite),256,ADD,50000,0.000037291,1340802956.969,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,rsa256(composite),256,ADD,50000,0.000666577,75010047.336,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,rsa256(composite),256,ADD,50000,0.000028485,1755340565.637,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,rsa256(composite),256,ADD,50000,0.000656417,76171134.272,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,rsa256(composite),256,ADD,50000,0.000028054,1782262430.701,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,rsa256(composite),256,ADD,50000,0.000697455,71689175.002,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,rsa256(composite),256,ADD,50000,0.000391133,127833712.403,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,rsa256(composite),256,ADD,50000,0.001636566,30551783.612,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACT,50000,0.000694034,72042612.246,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACT,50000,0.000138114,362019239.509,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACT,50000,0.000153816,325063067.711,0
library,NVIDIA RTX A2000,gpu,cgbn,rsa256(composite),256,SUBTRACT,50000,0.000033792,1479640151.515,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,rsa256(composite),256,SUBTRACT,50000,0.000101716,491563503.850,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,rsa256(composite),256,SUBTRACT,50000,0.000701458,71280085.504,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,rsa256(composite),256,SUBTRACT,50000,0.000057741,865934792.497,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,rsa256(composite),256,SUBTRACT,50000,0.000674788,74097356.145,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32,rsa256(composite),256,SUBTRACT,50000,0.000042392,1179469467.024,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32,rsa256(composite),256,SUBTRACT,50000,0.000652932,76577696.126,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,rsa256(composite),256,SUBTRACT,50000,0.000046891,1066299056.585,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,rsa256(composite),256,SUBTRACT,50000,0.000707488,70672619.162,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,rsa256(composite),256,SUBTRACT,50000,0.000036761,1360131009.323,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,rsa256(composite),256,SUBTRACT,50000,0.000655716,76252490.803,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,rsa256(composite),256,SUBTRACT,50000,0.000028035,1783505786.991,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,rsa256(composite),256,SUBTRACT,50000,0.000656376,76175889.713,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,rsa256(composite),256,SUBTRACT,50000,0.000028234,1770916057.527,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,rsa256(composite),256,SUBTRACT,50000,0.000689921,72472068.866,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,rsa256(composite),256,SUBTRACT,50000,0.000533759,93675294.660,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,rsa256(composite),256,SUBTRACT,50000,0.001764262,28340460.377,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,ADDMOD,50000,0.002459385,20330282.250,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,ADDMOD,50000,0.000299013,167216687.068,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,ADDMOD,50000,0.000775659,64461278.153,0
library,NVIDIA RTX A2000,gpu,cgbn,rsa256(composite),256,ADDMOD,50000,0.000033792,1479640151.515,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,rsa256(composite),256,ADDMOD,50000,0.000140270,356454852.803,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,rsa256(composite),256,ADDMOD,50000,0.000741686,67413996.602,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,rsa256(composite),256,ADDMOD,50000,0.000075446,662729958.400,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,rsa256(composite),256,ADDMOD,50000,0.000727749,68704966.247,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32,rsa256(composite),256,ADDMOD,50000,0.000046761,1069272265.928,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32,rsa256(composite),256,ADDMOD,50000,0.000687168,72762412.498,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,rsa256(composite),256,ADDMOD,50000,0.000033935,1473422378.352,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,rsa256(composite),256,ADDMOD,50000,0.000699081,71522423.341,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,rsa256(composite),256,ADDMOD,50000,0.000027894,1792497452.506,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,rsa256(composite),256,ADDMOD,50000,0.000657910,75998288.849,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,rsa256(composite),256,ADDMOD,50000,0.000025960,1926063399.584,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,rsa256(composite),256,ADDMOD,50000,0.000692315,72221418.194,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,rsa256(composite),256,ADDMOD,50000,0.000024888,2009021861.318,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,rsa256(composite),256,ADDMOD,50000,0.000692866,72164045.599,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,rsa256(composite),256,ADDMOD,50000,0.000514242,97230512.005,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,rsa256(composite),256,ADDMOD,50000,0.001686832,29641364.339,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACTMOD,50000,0.002574317,19422626.252,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACTMOD,50000,0.000287313,174026227.553,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACTMOD,50000,0.000887857,56315386.608,0
library,NVIDIA RTX A2000,gpu,cgbn,rsa256(composite),256,SUBTRACTMOD,50000,0.000033792,1479640151.515,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,rsa256(composite),256,SUBTRACTMOD,50000,0.000154155,324348224.718,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,rsa256(composite),256,SUBTRACTMOD,50000,0.000766313,65247515.502,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,rsa256(composite),256,SUBTRACTMOD,50000,0.000081717,611867513.078,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,rsa256(composite),256,SUBTRACTMOD,50000,0.000731848,68320180.806,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32,rsa256(composite),256,SUBTRACTMOD,50000,0.000050718,985843975.174,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32,rsa256(composite),256,SUBTRACTMOD,50000,0.000687649,72711562.564,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,rsa256(composite),256,SUBTRACTMOD,50000,0.000032824,1523297332.879,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,rsa256(composite),256,SUBTRACTMOD,50000,0.000642111,77868206.033,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,rsa256(composite),256,SUBTRACTMOD,50000,0.000027242,1835393360.911,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,rsa256(composite),256,SUBTRACTMOD,50000,0.000675115,74061477.806,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,rsa256(composite),256,SUBTRACTMOD,50000,0.000025339,1973209761.835,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,rsa256(composite),256,SUBTRACTMOD,50000,0.000679171,73619193.313,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,rsa256(composite),256,SUBTRACTMOD,50000,0.000025408,1967857605.747,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,rsa256(composite),256,SUBTRACTMOD,50000,0.000677566,73793544.468,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,rsa256(composite),256,SUBTRACTMOD,50000,0.000479326,104313215.407,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,rsa256(composite),256,SUBTRACTMOD,50000,0.001762999,28360761.244,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001121895,44567435.204,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000146391,341551354.446,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000281331,177726510.790,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.003805463,13139004.989,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.004770087,10481988.941,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001110185,45037541.315,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001926605,25952392.351,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000308987,161819234.896,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001109695,45057423.170,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000092759,539032432.053,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000882707,56643962.697,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000066949,746836535.626,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000877354,56989519.899,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000034556,1446935403.191,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000816847,61210958.598,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000034707,1440645392.583,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000856191,58398189.982,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.008384093,5963674.511,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.009637992,5187802.660,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001122466,44544767.635,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000145630,343335899.059,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000281322,177732394.469,0
library,NVIDIA RTX A2000,gpu,cgbn,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000033792,1479640151.515,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000673424,74247378.175,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001445598,34587763.038,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000224963,222258938.862,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001010993,49456302.358,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000096707,517027400.374,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000833050,60020381.875,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000108069,462668188.007,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000880071,56813600.446,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000075053,666192127.860,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000859701,58159750.492,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000053123,941218288.920,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000868496,57570818.624,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000050548,989167963.151,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000876150,57067847.557,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.005598053,8931676.245,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.007352865,6800070.170,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.009121890,5481320.101,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.001051491,47551518.737,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000302632,165216976.202,0
library,NVIDIA RTX A2000,gpu,cgbn,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000033984,1471280602.637,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000316989,157734340.882,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000899320,55597591.227,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000078562,636442311.659,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000732759,68235257.513,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000028877,1731506521.318,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000679514,73581961.095,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000042583,1174181291.691,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000694272,72017872.239,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000039697,1259550750.751,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000691446,72312274.492,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000038915,1284841239.679,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000694691,71974425.105,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000036932,1353854273.106,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000691082,72350280.037,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.028886212,1730929.619,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.031219300,1601573.392,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,COMPARE,50000,0.000387214,129127521.118,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,COMPARE,50000,0.000059504,840278770.425,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,COMPARE,50000,0.000091797,544681646.815,0
library,NVIDIA RTX A2000,gpu,cgbn,rsa256(composite),256,COMPARE,50000,0.000033792,1479640151.515,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,rsa256(composite),256,COMPARE,50000,0.000093881,532588897.266,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,rsa256(composite),256,COMPARE,50000,0.000701848,71240454.138,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,rsa256(composite),256,COMPARE,50000,0.000052090,959880767.374,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,rsa256(composite),256,COMPARE,50000,0.000665891,75087365.926,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,rsa256(composite),256,COMPARE,50000,0.000034236,1460436092.598,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,rsa256(composite),256,COMPARE,50000,0.000691446,72312274.492,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,rsa256(composite),256,COMPARE,50000,0.000028154,1775954058.882,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,rsa256(composite),256,COMPARE,50000,0.000668650,74777481.691,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,rsa256(composite),256,COMPARE,50000,0.000011994,4168899767.045,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,rsa256(composite),256,COMPARE,50000,0.000642891,77773676.624,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,rsa256(composite),256,COMPARE,50000,0.000011933,4190048482.010,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,rsa256(composite),256,COMPARE,50000,0.000619655,80689999.549,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,rsa256(composite),256,COMPARE,50000,0.000311871,160322667.073,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,rsa256(composite),256,COMPARE,50000,0.001704366,29336422.795,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,REDUCE,6250,0.000229200,27268719.474,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,REDUCE,6250,0.000036609,170721371.696,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,REDUCE,6250,0.000112866,55375375.653,0
library,NVIDIA RTX A2000,gpu,cgbn,rsa256(composite),256,REDUCE,50000,0.000048128,1038896276.596,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,rsa256(composite),256,REDUCE,18200,0.000280519,64879737.046,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,rsa256(composite),256,REDUCE,18200,0.000588831,30908672.947,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,rsa256(composite),256,REDUCE,18200,0.000171399,106185142.182,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,rsa256(composite),256,REDUCE,18200,0.000472960,38481050.409,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,rsa256(composite),256,REDUCE,18200,0.000085705,212356437.890,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,rsa256(composite),256,REDUCE,18200,0.000490178,37129342.510,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,rsa256(composite),256,REDUCE,18200,0.000056339,323041974.358,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,rsa256(composite),256,REDUCE,18200,0.000398078,45719682.843,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,rsa256(composite),256,REDUCE,18200,0.000060917,298767771.970,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,rsa256(composite),256,REDUCE,18200,0.000472653,38506072.225,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,rsa256(composite),256,REDUCE,18200,0.000055336,328897473.733,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,rsa256(composite),256,REDUCE,18200,0.000385583,47201204.774,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,rsa256(composite),256,REDUCE,18200,0.007800700,2333123.908,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,rsa256(composite),256,REDUCE,18200,0.008699493,2092076.048,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MODMUL,3125,0.000325095,9612577.492,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MODMUL,3125,0.000050988,61288872.653,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MODMUL,3125,0.000129549,24122177.970,0
library,NVIDIA RTX A2000,gpu,cgbn,rsa256(composite),256,MODMUL,50000,0.000167840,297902764.538,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,rsa256(composite),256,MODMUL,18200,0.000802672,22674281.030,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,rsa256(composite),256,MODMUL,18200,0.001125804,16166222.182,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,rsa256(composite),256,MODMUL,18200,0.000424477,42876265.025,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,rsa256(composite),256,MODMUL,18200,0.000734272,24786441.216,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,rsa256(composite),256,MODMUL,18200,0.000272929,66684073.627,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,rsa256(composite),256,MODMUL,18200,0.000689502,26395867.339,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,rsa256(composite),256,MODMUL,18200,0.000148968,122174021.099,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,rsa256(composite),256,MODMUL,18200,0.000538229,33814602.387,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,rsa256(composite),256,MODMUL,18200,0.000183835,99001986.903,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,rsa256(composite),256,MODMUL,18200,0.000558938,32561756.874,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,rsa256(composite),256,MODMUL,18200,0.000148286,122735701.929,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,rsa256(composite),256,MODMUL,18200,0.000538941,33769959.022,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,rsa256(composite),256,MODMUL,18200,0.031438502,578907.998,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,rsa256(composite),256,MODMUL,18200,0.030399920,598685.778,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MODEXP,781,0.007893038,98947.959,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MODEXP,781,0.000954885,817899.507,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MODEXP,781,0.001919789,406815.490,0
library,NVIDIA RTX A2000,gpu,cgbn,rsa256(composite),256,MODEXP,50000,0.061051901,818975.317,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,rsa256(composite),256,MODEXP,18200,0.040122085,453615.506,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,rsa256(composite),256,MODEXP,18200,0.040536704,448975.827,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,rsa256(composite),256,MODEXP,18200,0.006777455,2685373.748,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,rsa256(composite),256,MODEXP,18200,0.007266195,2504749.887,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,rsa256(composite),256,MODEXP,18200,0.005978416,3044284.782,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,rsa256(composite),256,MODEXP,18200,0.006367087,2858449.862,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,rsa256(composite),256,MODEXP,18200,0.003469229,5246122.796,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,rsa256(composite),256,MODEXP,18200,0.003887076,4682182.783,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,rsa256(composite),256,MODEXP,18200,0.005984558,3041160.372,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,rsa256(composite),256,MODEXP,18200,0.006345014,2868393.946,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,rsa256(composite),256,MODEXP,18200,0.003357436,5420803.939,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,rsa256(composite),256,MODEXP,18200,0.003725015,4885886.855,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,rsa256(composite),256,MODEXP,18200,1.448179659,12567.501,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,rsa256(composite),256,MODEXP,18200,1.437216697,12663.365,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,EXPONENTIATION,781,0.002778359,281101.207,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,EXPONENTIATION,781,0.000381353,2047969.631,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,EXPONENTIATION,781,0.005662932,137914.422,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,rsa256(composite),256,EXPONENTIATION,18200,0.053224236,341949.485,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,rsa256(composite),256,EXPONENTIATION,18200,0.053521077,340052.947,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,rsa256(composite),256,EXPONENTIATION,18200,0.014179890,1283507.837,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,rsa256(composite),256,EXPONENTIATION,18200,0.014502491,1254956.832,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,rsa256(composite),256,EXPONENTIATION,18200,0.001208966,15054190.953,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,rsa256(composite),256,EXPONENTIATION,18200,0.001526061,11926132.569,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,rsa256(composite),256,EXPONENTIATION,18200,0.000837287,21736873.293,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,rsa256(composite),256,EXPONENTIATION,18200,0.001235075,14735943.690,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,rsa256(composite),256,EXPONENTIATION,18200,0.001204816,15106044.935,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,rsa256(composite),256,EXPONENTIATION,18200,0.001518480,11985666.094,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,rsa256(composite),256,EXPONENTIATION,18200,0.000837435,21733029.649,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,rsa256(composite),256,EXPONENTIATION,18200,0.001197730,15195408.896,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,rsa256(composite),256,EXPONENTIATION,18200,0.598681115,30400.157,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,rsa256(composite),256,EXPONENTIATION,18200,0.605586187,30053.526,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,DIVIDE,6250,0.000284397,21976318.487,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,DIVIDE,6250,0.000046618,134067572.319,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,DIVIDE,6250,0.000119410,52340883.672,0
library,NVIDIA RTX A2000,gpu,cgbn,rsa256(composite),256,DIVIDE,50000,0.000063488,787550403.226,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,rsa256(composite),256,DIVIDE,18200,0.000557442,32649125.129,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,rsa256(composite),256,DIVIDE,18200,0.000957540,19007034.178,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,rsa256(composite),256,DIVIDE,18200,0.000363159,50115790.410,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,rsa256(composite),256,DIVIDE,18200,0.000827441,21995526.178,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,rsa256(composite),256,DIVIDE,18200,0.000181801,100109633.347,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,rsa256(composite),256,DIVIDE,18200,0.000586443,31034577.754,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,rsa256(composite),256,DIVIDE,18200,0.000183105,99396772.225,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,rsa256(composite),256,DIVIDE,18200,0.000653934,27831566.200,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,rsa256(composite),256,DIVIDE,18200,0.000183744,99050661.433,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,rsa256(composite),256,DIVIDE,18200,0.000562204,32372587.966,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,rsa256(composite),256,DIVIDE,18200,0.000176971,102841797.469,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,rsa256(composite),256,DIVIDE,18200,0.000606600,30003287.412,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,rsa256(composite),256,DIVIDE,18200,0.013035130,1396226.991,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,rsa256(composite),256,DIVIDE,18200,0.012898487,1411018.203,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,ISQRT,1562,0.000146301,10676584.945,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,ISQRT,1562,0.000023886,65394967.407,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,rsa256(composite),256,ISQRT,18200,0.007485727,2431293.630,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,rsa256(composite),256,ISQRT,18200,0.007991639,2277380.157,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,rsa256(composite),256,ISQRT,18200,0.005303099,3431955.658,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,rsa256(composite),256,ISQRT,18200,0.005813540,3130622.672,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,rsa256(composite),256,ISQRT,18200,0.001824937,9972947.943,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,rsa256(composite),256,ISQRT,18200,0.002199962,8272870.905,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,rsa256(composite),256,ISQRT,18200,0.001565033,11629150.441,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,rsa256(composite),256,ISQRT,18200,0.001964593,9264003.874,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,rsa256(composite),256,ISQRT,18200,0.001604606,11342345.406,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,rsa256(composite),256,ISQRT,18200,0.001935505,9403229.184,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,rsa256(composite),256,ISQRT,18200,0.001463713,12434130.498,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,rsa256(composite),256,ISQRT,18200,0.001830330,9943561.446,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,rsa256(composite),256,ISQRT,18200,0.087204248,208705.430,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,rsa256(composite),256,ISQRT,18200,0.108105052,168354.759,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MODMUL_R2,50000,0.005235582,9550037.275,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MODMUL_R2,50000,0.000310888,160829839.461,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MODMUL_R2,50000,0.000708440,70577575.386,0
library,NVIDIA RTX A2000,gpu,cgbn,rsa256(composite),256,MODMUL_R2,50000,0.000053248,939002403.846,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,rsa256(composite),256,MODMUL_R2,50000,0.000368749,135593664.713,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,rsa256(composite),256,MODMUL_R2,50000,0.000962871,51928029.210,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,rsa256(composite),256,MODMUL_R2,50000,0.000077659,643838188.665,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,rsa256(composite),256,MODMUL_R2,50000,0.000700608,71366595.061,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,rsa256(composite),256,MODMUL_R2,50000,0.000065978,757831983.400,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,rsa256(composite),256,MODMUL_R2,50000,0.000729140,68573946.364,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,rsa256(composite),256,MODMUL_R2,50000,0.000059485,840541883.768,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,rsa256(composite),256,MODMUL_R2,50000,0.000680334,73493320.630,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,rsa256(composite),256,MODMUL_R2,50000,0.000067210,743938852.091,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,rsa256(composite),256,MODMUL_R2,50000,0.000677257,73827234.636,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,rsa256(composite),256,MODMUL_R2,50000,0.000056328,887654032.605,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,rsa256(composite),256,MODMUL_R2,50000,0.000700121,71416245.582,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,rsa256(composite),256,MODMUL_R2,50000,0.033738444,1481988.931,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,rsa256(composite),256,MODMUL_R2,50000,0.034126762,1465125.802,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,ADD,25000,0.000487869,51243298.304,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,ADD,25000,0.000070516,354529367.637,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,ADD,25000,0.000081677,306083758.267,0
library,NVIDIA RTX A2000,gpu,cgbn,brainpoolP512r1,512,ADD,50000,0.000057344,871930803.571,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,brainpoolP512r1,512,ADD,25000,0.000120642,207225258.997,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,brainpoolP512r1,512,ADD,25000,0.000702230,35600858.340,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,brainpoolP512r1,512,ADD,25000,0.000069504,359693223.814,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,brainpoolP512r1,512,ADD,25000,0.000732799,34115764.256,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32,brainpoolP512r1,512,ADD,25000,0.000078984,316521384.775,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32,brainpoolP512r1,512,ADD,25000,0.000731985,34153701.383,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,brainpoolP512r1,512,ADD,25000,0.000056981,438741898.893,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,brainpoolP512r1,512,ADD,25000,0.000711365,35143718.825,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,brainpoolP512r1,512,ADD,25000,0.000058021,430875531.300,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,brainpoolP512r1,512,ADD,25000,0.000679562,36788409.070,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,brainpoolP512r1,512,ADD,25000,0.000032953,758656575.191,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,brainpoolP512r1,512,ADD,25000,0.000637240,39231700.405,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,brainpoolP512r1,512,ADD,25000,0.000045206,553019068.809,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,brainpoolP512r1,512,ADD,25000,0.000668800,37380358.409,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,brainpoolP512r1,512,ADD,25000,0.000509463,49071252.870,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,brainpoolP512r1,512,ADD,25000,0.001698584,14718139.136,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACT,25000,0.000383307,65221833.308,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACT,25000,0.000078429,318757740.491,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACT,25000,0.000096787,258299772.911,0
library,NVIDIA RTX A2000,gpu,cgbn,brainpoolP512r1,512,SUBTRACT,50000,0.000057344,871930803.571,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,brainpoolP512r1,512,SUBTRACT,25000,0.000120331,207760948.577,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,brainpoolP512r1,512,SUBTRACT,25000,0.000692500,36101093.380,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,brainpoolP512r1,512,SUBTRACT,25000,0.000068662,364103704.307,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,brainpoolP512r1,512,SUBTRACT,25000,0.000689856,36239458.802,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32,brainpoolP512r1,512,SUBTRACT,25000,0.000079985,312559477.429,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32,brainpoolP512r1,512,SUBTRACT,25000,0.000746162,33504803.011,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,brainpoolP512r1,512,SUBTRACT,25000,0.000056169,445085400.673,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,brainpoolP512r1,512,SUBTRACT,25000,0.000718619,34788961.796,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,brainpoolP512r1,512,SUBTRACT,25000,0.000055999,446440021.288,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,brainpoolP512r1,512,SUBTRACT,25000,0.000700474,35690129.340,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,brainpoolP512r1,512,SUBTRACT,25000,0.000033184,753376149.982,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,brainpoolP512r1,512,SUBTRACT,25000,0.000637290,39228604.456,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,brainpoolP512r1,512,SUBTRACT,25000,0.000045066,554744789.105,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,brainpoolP512r1,512,SUBTRACT,25000,0.000687786,36348495.406,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,brainpoolP512r1,512,SUBTRACT,25000,0.000449247,55648708.163,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,brainpoolP512r1,512,SUBTRACT,25000,0.001800471,13885255.692,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,ADDMOD,25000,0.001314485,19018848.846,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,ADDMOD,25000,0.000170858,146320639.714,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,ADDMOD,25000,0.000472168,52947206.930,0
library,NVIDIA RTX A2000,gpu,cgbn,brainpoolP512r1,512,ADDMOD,50000,0.000057344,871930803.571,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,brainpoolP512r1,512,ADDMOD,25000,0.000162783,153578616.259,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,brainpoolP512r1,512,ADDMOD,25000,0.000766923,32597808.561,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,brainpoolP512r1,512,ADDMOD,25000,0.000090284,276903154.463,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,brainpoolP512r1,512,ADDMOD,25000,0.000748710,33390775.462,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32,brainpoolP512r1,512,ADDMOD,25000,0.000087009,287327220.765,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32,brainpoolP512r1,512,ADDMOD,25000,0.000741423,33718938.270,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,brainpoolP512r1,512,ADDMOD,25000,0.000067170,372191195.596,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,brainpoolP512r1,512,ADDMOD,25000,0.000692308,36111097.718,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,brainpoolP512r1,512,ADDMOD,25000,0.000067512,370306878.190,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,brainpoolP512r1,512,ADDMOD,25000,0.000724188,34521420.130,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,brainpoolP512r1,512,ADDMOD,25000,0.000034887,716592247.731,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,brainpoolP512r1,512,ADDMOD,25000,0.000633102,39488113.353,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,brainpoolP512r1,512,ADDMOD,25000,0.000033755,740634190.487,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,brainpoolP512r1,512,ADDMOD,25000,0.000663530,37677266.934,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,brainpoolP512r1,512,ADDMOD,25000,0.000490506,50967757.468,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,brainpoolP512r1,512,ADDMOD,25000,0.001787246,13988000.048,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACTMOD,25000,0.001344574,18593253.978,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000161780,154530801.911,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000534237,46795678.770,0
library,NVIDIA RTX A2000,gpu,cgbn,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000057344,871930803.571,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000187800,133120152.344,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000784678,31860186.911,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000100112,249721338.865,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000710937,35164850.274,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000093781,266577411.442,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000720252,34710060.062,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000066699,374820860.969,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000698540,35788912.769,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000067130,372413229.745,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000726722,34401041.635,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000035058,713108562.017,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000636348,39286706.192,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000033885,737784344.767,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000666495,37509635.517,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000564198,44310684.496,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,brainpoolP512r1,512,SUBTRACTMOD,25000,0.001666051,15005540.025,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.001344834,18589661.524,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.000173564,144039029.206,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.000198781,125766236.882,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.011251373,2221951.003,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.012421061,2012710.462,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.002849831,8772449.259,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.003730885,6700823.438,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.001963917,12729661.125,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.002848306,8777147.640,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.000198033,126241772.795,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.000983421,25421468.644,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.000198904,125688505.984,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.001019679,24517519.491,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.000098680,253343767.755,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.000865971,28869335.760,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.000098670,253370071.546,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.000883123,28308630.461,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.015986049,1563863.629,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.019438074,1286135.655,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.001344253,18597685.296,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000175016,142844082.119,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000199273,125455888.731,0
library,NVIDIA RTX A2000,gpu,cgbn,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000057344,871930803.571,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.002396095,10433644.760,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.003236510,7724370.104,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000665128,37586737.163,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.001465809,17055432.747,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000336441,74307186.970,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.001147607,21784454.576,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000283588,88156143.186,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.001109454,22533609.678,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000284640,87830205.150,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.001099593,22735684.805,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000130772,191172920.272,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000890137,28085573.936,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000183013,136602117.969,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000969401,25789131.182,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.014342629,1743055.572,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.016139684,1548977.020,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.010628715,2352118.730,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.001127438,22174174.464,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000344573,72553652.880,0
library,NVIDIA RTX A2000,gpu,cgbn,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000057344,871930803.571,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000589052,42441059.305,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.001186019,21078923.703,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000182871,136708557.925,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000793657,31499743.130,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000048463,515854980.110,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000701556,35635074.938,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000090195,277177638.725,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000741923,33696208.799,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000065246,383167215.268,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000694351,36004830.770,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000079523,314375087.543,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000680864,36718060.054,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000058392,428140380.873,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000649954,38464249.165,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.044784497,558228.883,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.043345714,576758.293,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,COMPARE,25000,0.000193311,129325350.009,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,COMPARE,25000,0.000038493,649461569.728,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,COMPARE,25000,0.000062410,400578188.981,0
library,NVIDIA RTX A2000,gpu,cgbn,brainpoolP512r1,512,COMPARE,50000,0.000057344,871930803.571,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,brainpoolP512r1,512,COMPARE,25000,0.000119119,209874245.327,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,brainpoolP512r1,512,COMPARE,25000,0.000724472,34507884.859,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,brainpoolP512r1,512,COMPARE,25000,0.000066978,373257304.949,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,brainpoolP512r1,512,COMPARE,25000,0.000660120,37871926.275,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,brainpoolP512r1,512,COMPARE,25000,0.000051139,488864425.423,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,brainpoolP512r1,512,COMPARE,25000,0.000671769,37215197.988,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,brainpoolP512r1,512,COMPARE,25000,0.000051269,487621173.479,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,brainpoolP512r1,512,COMPARE,25000,0.000639064,39119698.043,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,brainpoolP512r1,512,COMPARE,25000,0.000016261,1537431019.473,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,brainpoolP512r1,512,COMPARE,25000,0.000626769,39887108.834,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,brainpoolP512r1,512,COMPARE,25000,0.000016020,1560580524.388,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,brainpoolP512r1,512,COMPARE,25000,0.000656836,38061266.575,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,brainpoolP512r1,512,COMPARE,25000,0.000328435,76118647.740,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,brainpoolP512r1,512,COMPARE,25000,0.001664449,15019989.839,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,REDUCE,3125,0.000117877,26510572.806,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,REDUCE,3125,0.000027202,114880964.119,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,REDUCE,3125,0.000079943,39090416.832,0
library,NVIDIA RTX A2000,gpu,cgbn,brainpoolP512r1,512,REDUCE,50000,0.000084992,588290662.651,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,brainpoolP512r1,512,REDUCE,18200,0.000858350,21203478.556,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,brainpoolP512r1,512,REDUCE,18200,0.001373310,13252655.805,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,brainpoolP512r1,512,REDUCE,18200,0.000688444,26436431.813,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,brainpoolP512r1,512,REDUCE,18200,0.001427766,12747187.438,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,brainpoolP512r1,512,REDUCE,18200,0.000251025,72502749.899,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,brainpoolP512r1,512,REDUCE,18200,0.000858790,21192602.253,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,brainpoolP512r1,512,REDUCE,18200,0.000234202,77710534.319,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,brainpoolP512r1,512,REDUCE,18200,0.000791057,23007191.266,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,brainpoolP512r1,512,REDUCE,18200,0.000236146,77070914.958,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,brainpoolP512r1,512,REDUCE,18200,0.000766930,23730974.088,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,brainpoolP512r1,512,REDUCE,18200,0.000294157,61871657.649,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,brainpoolP512r1,512,REDUCE,18200,0.000850741,21393119.465,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,brainpoolP512r1,512,REDUCE,18200,0.029130058,624784.193,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,brainpoolP512r1,512,REDUCE,18200,0.030905111,588899.362,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL,1562,0.000324083,4819745.645,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL,1562,0.000047041,33205003.546,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL,1562,0.000119559,13064628.350,0
library,NVIDIA RTX A2000,gpu,cgbn,brainpoolP512r1,512,MODMUL,50000,0.000542432,92177452.658,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,brainpoolP512r1,512,MODMUL,18200,0.002540883,7162865.159,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,brainpoolP512r1,512,MODMUL,18200,0.003172829,5736206.098,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,brainpoolP512r1,512,MODMUL,18200,0.001789983,10167692.274,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,brainpoolP512r1,512,MODMUL,18200,0.002544784,7151884.047,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,brainpoolP512r1,512,MODMUL,18200,0.000956539,19026928.092,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,brainpoolP512r1,512,MODMUL,18200,0.001574563,11558763.663,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,brainpoolP512r1,512,MODMUL,18200,0.000780467,23319373.281,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,brainpoolP512r1,512,MODMUL,18200,0.001388652,13106233.038,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,brainpoolP512r1,512,MODMUL,18200,0.000971916,18725895.586,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,brainpoolP512r1,512,MODMUL,18200,0.001500184,12131847.454,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,brainpoolP512r1,512,MODMUL,18200,0.000998996,18218287.550,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,brainpoolP512r1,512,MODMUL,18200,0.001597661,11391656.070,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,brainpoolP512r1,512,MODMUL,18200,0.070949179,256521.643,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,brainpoolP512r1,512,MODMUL,18200,0.071263956,255388.573,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MODEXP,390,0.021250448,18352.554,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MODEXP,390,0.002155202,180957.479,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MODEXP,390,0.002153649,181088.005,0
library,NVIDIA RTX A2000,gpu,cgbn,brainpoolP512r1,512,MODEXP,50000,0.171768829,291088.903,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,brainpoolP512r1,512,MODEXP,18200,0.489869938,37152.719,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,brainpoolP512r1,512,MODEXP,18200,0.490906550,37074.266,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,brainpoolP512r1,512,MODEXP,18200,0.059492115,305922.895,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,brainpoolP512r1,512,MODEXP,18200,0.060323814,301705.062,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,brainpoolP512r1,512,MODEXP,18200,0.044904671,405303.050,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,brainpoolP512r1,512,MODEXP,18200,0.045462669,400328.455,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,brainpoolP512r1,512,MODEXP,18200,0.021175249,859494.006,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,brainpoolP512r1,512,MODEXP,18200,0.021844693,833154.287,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,brainpoolP512r1,512,MODEXP,18200,0.044159213,412145.025,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,brainpoolP512r1,512,MODEXP,18200,0.044772957,406495.372,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,brainpoolP512r1,512,MODEXP,18200,0.023091368,788173.299,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,brainpoolP512r1,512,MODEXP,18200,0.023565574,772313.048,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,brainpoolP512r1,512,MODEXP,18200,12.207959394,1490.831,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,brainpoolP512r1,512,MODEXP,18200,12.584155072,1446.263,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,EXPONENTIATION,390,0.004824261,80841.394,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,EXPONENTIATION,390,0.000594102,656453.102,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,EXPONENTIATION,390,0.007098181,54943.653,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,brainpoolP512r1,512,EXPONENTIATION,18200,0.425519317,42771.266,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,brainpoolP512r1,512,EXPONENTIATION,18200,0.425869111,42736.135,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,brainpoolP512r1,512,EXPONENTIATION,18200,0.116647294,156025.909,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,brainpoolP512r1,512,EXPONENTIATION,18200,0.117261156,155209.112,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,brainpoolP512r1,512,EXPONENTIATION,18200,0.036290882,501503.378,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,brainpoolP512r1,512,EXPONENTIATION,18200,0.036575342,497603.008,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,brainpoolP512r1,512,EXPONENTIATION,18200,0.032356853,562477.451,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,brainpoolP512r1,512,EXPONENTIATION,18200,0.033254475,547294.757,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,brainpoolP512r1,512,EXPONENTIATION,18200,0.035940032,506399.098,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,brainpoolP512r1,512,EXPONENTIATION,18200,0.036420389,499720.087,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,brainpoolP512r1,512,EXPONENTIATION,18200,0.032478904,560363.730,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,brainpoolP512r1,512,EXPONENTIATION,18200,0.032853216,553979.252,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,brainpoolP512r1,512,EXPONENTIATION,18200,6.598207767,2758.325,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,brainpoolP512r1,512,EXPONENTIATION,18200,6.768749940,2688.827,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,DIVIDE,3125,0.000167753,18628621.553,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,DIVIDE,3125,0.000025028,124858346.357,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,DIVIDE,3125,0.000065655,47596964.410,0
library,NVIDIA RTX A2000,gpu,cgbn,brainpoolP512r1,512,DIVIDE,50000,0.000115648,432346430.548,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,brainpoolP512r1,512,DIVIDE,18200,0.001635524,11127935.402,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,brainpoolP512r1,512,DIVIDE,18200,0.002344975,7761275.952,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,brainpoolP512r1,512,DIVIDE,18200,0.001335509,13627761.659,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,brainpoolP512r1,512,DIVIDE,18200,0.002200472,8270952.144,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,brainpoolP512r1,512,DIVIDE,18200,0.000468796,38822846.988,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,brainpoolP512r1,512,DIVIDE,18200,0.001176143,15474307.965,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,brainpoolP512r1,512,DIVIDE,18200,0.000475108,38307104.472,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,brainpoolP512r1,512,DIVIDE,18200,0.001151082,15811209.215,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,brainpoolP512r1,512,DIVIDE,18200,0.000473213,38460450.802,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,brainpoolP512r1,512,DIVIDE,18200,0.001179055,15436086.699,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,brainpoolP512r1,512,DIVIDE,18200,0.000403958,45054147.226,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,brainpoolP512r1,512,DIVIDE,18200,0.001022482,17799819.285,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,brainpoolP512r1,512,DIVIDE,18200,0.048038690,378861.286,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,brainpoolP512r1,512,DIVIDE,18200,0.037372726,486986.151,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,ISQRT,781,0.000124880,6253998.199,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,ISQRT,781,0.000023324,33484761.402,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,brainpoolP512r1,512,ISQRT,18200,0.037642193,483499.991,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,brainpoolP512r1,512,ISQRT,18200,0.038328501,474842.463,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,brainpoolP512r1,512,ISQRT,18200,0.039550192,460174.760,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,brainpoolP512r1,512,ISQRT,18200,0.040412190,450359.157,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,brainpoolP512r1,512,ISQRT,18200,0.007340278,2479470.277,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,brainpoolP512r1,512,ISQRT,18200,0.008054587,2259582.008,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,brainpoolP512r1,512,ISQRT,18200,0.006521723,2790673.469,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,brainpoolP512r1,512,ISQRT,18200,0.007883373,2308656.543,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,brainpoolP512r1,512,ISQRT,18200,0.006819339,2668880.552,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,brainpoolP512r1,512,ISQRT,18200,0.008015756,2270528.054,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,brainpoolP512r1,512,ISQRT,18200,0.008694560,2093262.850,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,brainpoolP512r1,512,ISQRT,18200,0.009540642,1907628.502,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,brainpoolP512r1,512,ISQRT,18200,0.504350289,36086.031,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,brainpoolP512r1,512,ISQRT,18200,0.496879919,36628.568,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL_R2,25000,0.005219961,4789308.193,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL_R2,25000,0.000268035,93271527.450,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL_R2,25000,0.000650842,38411795.614,0
library,NVIDIA RTX A2000,gpu,cgbn,brainpoolP512r1,512,MODMUL_R2,50000,0.000084992,588290662.651,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,brainpoolP512r1,512,MODMUL_R2,25000,0.000866475,28852548.556,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,brainpoolP512r1,512,MODMUL_R2,25000,0.001520833,16438359.761,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,brainpoolP512r1,512,MODMUL_R2,25000,0.000172682,144774699.055,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,brainpoolP512r1,512,MODMUL_R2,25000,0.000829946,30122432.637,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,brainpoolP512r1,512,MODMUL_R2,25000,0.000116846,213957577.593,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,brainpoolP512r1,512,MODMUL_R2,25000,0.000770770,32435096.809,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,brainpoolP512r1,512,MODMUL_R2,25000,0.000078081,320179696.800,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,brainpoolP512r1,512,MODMUL_R2,25000,0.000749066,33374875.172,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,brainpoolP512r1,512,MODMUL_R2,25000,0.000111496,224223137.707,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,brainpoolP512r1,512,MODMUL_R2,25000,0.000772030,32382157.424,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,brainpoolP512r1,512,MODMUL_R2,25000,0.000075306,331979688.099,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,brainpoolP512r1,512,MODMUL_R2,25000,0.000687757,36350070.483,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,brainpoolP512r1,512,MODMUL_R2,25000,0.050937327,490799.217,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,brainpoolP512r1,512,MODMUL_R2,25000,0.056471383,442702.101,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p1024,1024,ADD,12500,0.000289727,43144062.850,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p1024,1024,ADD,12500,0.000044675,279800971.461,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p1024,1024,ADD,12500,0.000056158,222586988.176,0
library,NVIDIA RTX A2000,gpu,cgbn,p1024,1024,ADD,50000,0.000105472,474059466.019,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,p1024,1024,ADD,18200,0.000282583,64405895.428,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,p1024,1024,ADD,18200,0.001168636,15573715.922,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,p1024,1024,ADD,18200,0.000151441,120179211.335,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,p1024,1024,ADD,18200,0.000974875,18669061.218,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32,p1024,1024,ADD,18200,0.000143648,126698486.115,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32,p1024,1024,ADD,18200,0.000963993,18879797.887,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,p1024,1024,ADD,18200,0.000076086,239202188.536,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,p1024,1024,ADD,18200,0.000918284,19819575.250,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,p1024,1024,ADD,18200,0.000077739,234115645.927,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,p1024,1024,ADD,18200,0.000906001,20088281.482,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,p1024,1024,ADD,18200,0.000063803,285252834.659,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,p1024,1024,ADD,18200,0.000953631,19084959.170,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,p1024,1024,ADD,18200,0.000064254,283251698.701,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,p1024,1024,ADD,18200,0.000945855,19241851.348,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,p1024,1024,ADD,18200,0.000605256,30069906.087,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,p1024,1024,ADD,18200,0.002261508,8047726.534,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p1024,1024,SUBTRACT,12500,0.000231064,54097583.665,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p1024,1024,SUBTRACT,12500,0.000038884,321471888.098,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p1024,1024,SUBTRACT,12500,0.000056979,219378120.658,0
library,NVIDIA RTX A2000,gpu,cgbn,p1024,1024,SUBTRACT,50000,0.000106208,470774329.617,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,p1024,1024,SUBTRACT,18200,0.000281792,64586614.746,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,p1024,1024,SUBTRACT,18200,0.001158596,15708667.854,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,p1024,1024,SUBTRACT,18200,0.000151683,119987359.069,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,p1024,1024,SUBTRACT,18200,0.000992038,18346066.109,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32,p1024,1024,SUBTRACT,18200,0.000143677,126673026.841,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32,p1024,1024,SUBTRACT,18200,0.000973001,18705014.393,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,p1024,1024,SUBTRACT,18200,0.000076227,238760888.437,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,p1024,1024,SUBTRACT,18200,0.000904237,20127468.359,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,p1024,1024,SUBTRACT,18200,0.000078522,231783153.013,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,p1024,1024,SUBTRACT,18200,0.000874048,20822652.692,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,p1024,1024,SUBTRACT,18200,0.000063532,286469665.872,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,p1024,1024,SUBTRACT,18200,0.000901069,20198220.181,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,p1024,1024,SUBTRACT,18200,0.000064124,283823525.435,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,p1024,1024,SUBTRACT,18200,0.000944492,19269609.647,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,p1024,1024,SUBTRACT,18200,0.000711241,25589082.986,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,p1024,1024,SUBTRACT,18200,0.002406718,7562164.914,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p1024,1024,ADDMOD,12500,0.000857318,14580357.139,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p1024,1024,ADDMOD,12500,0.000109660,113989203.880,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p1024,1024,ADDMOD,12500,0.000358540,34863649.186,0
library,NVIDIA RTX A2000,gpu,cgbn,p1024,1024,ADDMOD,50000,0.000105472,474059466.019,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,p1024,1024,ADDMOD,18200,0.000413896,43972442.991,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,p1024,1024,ADDMOD,18200,0.001297534,14026612.669,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,p1024,1024,ADDMOD,18200,0.000218079,83455832.512,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,p1024,1024,ADDMOD,18200,0.001041634,17472547.800,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32,p1024,1024,ADDMOD,18200,0.000181932,100037375.334,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32,p1024,1024,ADDMOD,18200,0.001006526,18081997.945,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,p1024,1024,ADDMOD,18200,0.000098731,184338576.735,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,p1024,1024,ADDMOD,18200,0.000958121,18995505.525,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,p1024,1024,ADDMOD,18200,0.000100565,180977220.037,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,p1024,1024,ADDMOD,18200,0.000892163,20399853.852,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,p1024,1024,ADDMOD,18200,0.000043223,421075225.098,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,p1024,1024,ADDMOD,18200,0.000891310,20419378.976,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,p1024,1024,ADDMOD,18200,0.000043053,422733001.575,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,p1024,1024,ADDMOD,18200,0.000927079,19631562.226,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,p1024,1024,ADDMOD,18200,0.000671344,27109802.590,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,p1024,1024,ADDMOD,18200,0.002363134,7701636.322,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p1024,1024,SUBTRACTMOD,12500,0.000767505,16286543.171,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p1024,1024,SUBTRACTMOD,12500,0.000096684,129286732.040,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p1024,1024,SUBTRACTMOD,12500,0.000378708,33006929.556,0
library,NVIDIA RTX A2000,gpu,cgbn,p1024,1024,SUBTRACTMOD,50000,0.000106304,470349187.237,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,p1024,1024,SUBTRACTMOD,18200,0.000410870,44296279.868,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,p1024,1024,SUBTRACTMOD,18200,0.001272035,14307784.309,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,p1024,1024,SUBTRACTMOD,18200,0.000214492,84851661.015,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,p1024,1024,SUBTRACTMOD,18200,0.001070297,17004620.684,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32,p1024,1024,SUBTRACTMOD,18200,0.000183976,98925809.554,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32,p1024,1024,SUBTRACTMOD,18200,0.001005884,18093532.949,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,p1024,1024,SUBTRACTMOD,18200,0.000099763,182431863.301,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,p1024,1024,SUBTRACTMOD,18200,0.000922964,19719079.939,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,p1024,1024,SUBTRACTMOD,18200,0.000099172,183519755.804,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,p1024,1024,SUBTRACTMOD,18200,0.000907663,20051489.229,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,p1024,1024,SUBTRACTMOD,18200,0.000043434,419025691.978,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,p1024,1024,SUBTRACTMOD,18200,0.000902763,20160338.126,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,p1024,1024,SUBTRACTMOD,18200,0.000043935,414246978.205,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,p1024,1024,SUBTRACTMOD,18200,0.000898393,20258397.620,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,p1024,1024,SUBTRACTMOD,18200,0.000737851,24666210.416,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,p1024,1024,SUBTRACTMOD,18200,0.002621451,6942718.850,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.002360477,5295539.590,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.000290568,43019191.974,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.000361033,34622894.525,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,18200,0.041608634,437409.212,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,18200,0.042737816,425852.360,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,p1024,1024,MULTIPLYOPERANDSCANNING,18200,0.010666715,1706242.272,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,p1024,1024,MULTIPLYOPERANDSCANNING,18200,0.012168582,1495655.007,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32,p1024,1024,MULTIPLYOPERANDSCANNING,18200,0.005922798,3072871.942,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32,p1024,1024,MULTIPLYOPERANDSCANNING,18200,0.006950805,2618401.767,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,p1024,1024,MULTIPLYOPERANDSCANNING,18200,0.000266214,68366076.709,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,p1024,1024,MULTIPLYOPERANDSCANNING,18200,0.001351883,13462706.455,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,p1024,1024,MULTIPLYOPERANDSCANNING,18200,0.000261845,69506751.449,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,p1024,1024,MULTIPLYOPERANDSCANNING,18200,0.001297033,14032021.160,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,p1024,1024,MULTIPLYOPERANDSCANNING,18200,0.000183514,99174822.234,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,p1024,1024,MULTIPLYOPERANDSCANNING,18200,0.001264380,14394403.431,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,p1024,1024,MULTIPLYOPERANDSCANNING,18200,0.000164988,110310751.584,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,p1024,1024,MULTIPLYOPERANDSCANNING,18200,0.001268036,14352907.838,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,18200,0.049021275,371267.375,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,18200,0.050119992,363128.547,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.002359564,5297587.942,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.000290799,42985023.844,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.000359641,34756845.055,0
library,NVIDIA RTX A2000,gpu,cgbn,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000111584,448092916.547,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,18200,0.014538400,1251857.161,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,18200,0.015686988,1160197.236,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,p1024,1024,MULTIPLYPRODUCTSCANNING,18200,0.003759119,4841560.301,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,p1024,1024,MULTIPLYPRODUCTSCANNING,18200,0.004863143,3742435.730,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32,p1024,1024,MULTIPLYPRODUCTSCANNING,18200,0.001463910,12432461.392,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32,p1024,1024,MULTIPLYPRODUCTSCANNING,18200,0.002484693,7324848.017,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,p1024,1024,MULTIPLYPRODUCTSCANNING,18200,0.000959015,18977814.850,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,p1024,1024,MULTIPLYPRODUCTSCANNING,18200,0.002046466,8893380.649,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,p1024,1024,MULTIPLYPRODUCTSCANNING,18200,0.000960104,18956276.430,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,p1024,1024,MULTIPLYPRODUCTSCANNING,18200,0.002027645,8975931.197,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,p1024,1024,MULTIPLYPRODUCTSCANNING,18200,0.000478134,38064678.057,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,p1024,1024,MULTIPLYPRODUCTSCANNING,18200,0.001580930,11512208.873,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,p1024,1024,MULTIPLYPRODUCTSCANNING,18200,0.000446851,40729426.862,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,p1024,1024,MULTIPLYPRODUCTSCANNING,18200,0.001550101,11741172.355,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,18200,0.033628008,541215.520,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,18200,0.035178135,517366.822,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.016107164,776052.181,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.001734363,7207257.426,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.000570095,21926164.653,0
library,NVIDIA RTX A2000,gpu,cgbn,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000164864,303280279.503,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,18200,0.002607001,6981201.794,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,18200,0.003553198,5122146.708,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,p1024,1024,MONTGOMERYMULTIPLICATION,18200,0.000573885,31713687.670,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,p1024,1024,MONTGOMERYMULTIPLICATION,18200,0.001431663,12712492.549,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32,p1024,1024,MONTGOMERYMULTIPLICATION,18200,0.000201900,90143832.670,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32,p1024,1024,MONTGOMERYMULTIPLICATION,18200,0.001003750,18132010.968,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,p1024,1024,MONTGOMERYMULTIPLICATION,18200,0.000225707,80635529.446,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,p1024,1024,MONTGOMERYMULTIPLICATION,18200,0.001085578,16765255.211,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,p1024,1024,MONTGOMERYMULTIPLICATION,18200,0.000174026,104582070.956,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,p1024,1024,MONTGOMERYMULTIPLICATION,18200,0.000988919,18403928.634,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,p1024,1024,MONTGOMERYMULTIPLICATION,18200,0.000180318,100932785.150,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,p1024,1024,MONTGOMERYMULTIPLICATION,18200,0.001032492,17627252.732,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,p1024,1024,MONTGOMERYMULTIPLICATION,18200,0.000125101,145482640.716,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,p1024,1024,MONTGOMERYMULTIPLICATION,18200,0.001005499,18100471.077,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,18200,0.095043669,191490.923,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,18200,0.113817541,159905.054,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p1024,1024,COMPARE,12500,0.000093278,134007336.482,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p1024,1024,COMPARE,12500,0.000026480,472049126.016,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p1024,1024,COMPARE,12500,0.000048263,258997584.038,0
library,NVIDIA RTX A2000,gpu,cgbn,p1024,1024,COMPARE,50000,0.000105440,474203338.392,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,p1024,1024,COMPARE,18200,0.000242757,74972190.367,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,p1024,1024,COMPARE,18200,0.001116947,16294412.882,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,p1024,1024,COMPARE,18200,0.000129319,140737468.559,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,p1024,1024,COMPARE,18200,0.000950057,19156744.132,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,p1024,1024,COMPARE,18200,0.000068893,264178838.181,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,p1024,1024,COMPARE,18200,0.000931441,19539619.665,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,p1024,1024,COMPARE,18200,0.000071488,254587040.083,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,p1024,1024,COMPARE,18200,0.000875401,20790464.646,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,p1024,1024,COMPARE,18200,0.000015221,1195747488.025,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,p1024,1024,COMPARE,18200,0.000874998,20800046.403,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,p1024,1024,COMPARE,18200,0.000015261,1192609617.771,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,p1024,1024,COMPARE,18200,0.000859506,21174943.435,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,p1024,1024,COMPARE,18200,0.000371926,48934404.391,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,p1024,1024,COMPARE,18200,0.002097072,8678767.018,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p1024,1024,REDUCE,1562,0.000046229,33788322.033,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p1024,1024,REDUCE,1562,0.000018585,84044133.548,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p1024,1024,REDUCE,1562,0.000065235,23944047.185,0
library,NVIDIA RTX A2000,gpu,cgbn,p1024,1024,REDUCE,50000,0.000128000,390625000.000,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,p1024,1024,REDUCE,18200,0.004676146,3892093.685,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,p1024,1024,REDUCE,18200,0.005808252,3133472.908,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,p1024,1024,REDUCE,18200,0.002553101,7128586.826,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,p1024,1024,REDUCE,18200,0.003442080,5287500.950,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,p1024,1024,REDUCE,18200,0.000660978,27534945.044,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,p1024,1024,REDUCE,18200,0.001500510,12129211.994,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,p1024,1024,REDUCE,18200,0.000677358,26869085.995,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,p1024,1024,REDUCE,18200,0.001517721,11991660.247,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,p1024,1024,REDUCE,18200,0.000558558,32583908.209,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,p1024,1024,REDUCE,18200,0.001404139,12961677.751,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,p1024,1024,REDUCE,18200,0.000607432,29962208.129,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,p1024,1024,REDUCE,18200,0.001480335,12294511.347,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,p1024,1024,REDUCE,18200,0.060523441,300709.935,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,p1024,1024,REDUCE,18200,0.069213729,262953.610,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p1024,1024,MODMUL,781,0.000471376,1656851.843,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p1024,1024,MODMUL,781,0.000060496,12909961.429,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p1024,1024,MODMUL,781,0.000155879,5010290.514,0
library,NVIDIA RTX A2000,gpu,cgbn,p1024,1024,MODMUL,50000,0.001363872,36660331.761,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,p1024,1024,MODMUL,18200,0.018083404,1006447.690,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,p1024,1024,MODMUL,18200,0.018993571,958218.948,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,p1024,1024,MODMUL,18200,0.006923164,2628855.732,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,p1024,1024,MODMUL,18200,0.008387952,2169778.720,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,p1024,1024,MODMUL,18200,0.003025989,6014562.416,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,p1024,1024,MODMUL,18200,0.004064787,4477479.188,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,p1024,1024,MODMUL,18200,0.002249143,8091970.720,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,p1024,1024,MODMUL,18200,0.003131577,5811768.666,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,p1024,1024,MODMUL,18200,0.002949852,6169800.023,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,p1024,1024,MODMUL,18200,0.003963378,4592042.146,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,p1024,1024,MODMUL,18200,0.002213391,8222679.025,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,p1024,1024,MODMUL,18200,0.003128033,5818352.695,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,p1024,1024,MODMUL,18200,0.258423735,70426.967,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,p1024,1024,MODMUL,18200,0.293168947,62080.245,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p1024,1024,MODEXP,195,0.070172154,2778.880,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p1024,1024,MODEXP,195,0.006987307,27907.747,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p1024,1024,MODEXP,195,0.005855894,33299.785,0
library,NVIDIA RTX A2000,gpu,cgbn,p1024,1024,MODEXP,50000,0.923003912,54170.951,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,p1024,1024,MODEXP,18200,4.579885645,3973.898,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,p1024,1024,MODEXP,18200,4.594887288,3960.924,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,p1024,1024,MODEXP,18200,0.744081890,24459.673,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,p1024,1024,MODEXP,18200,0.744883521,24433.350,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,p1024,1024,MODEXP,18200,0.363093211,50124.870,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,p1024,1024,MODEXP,18200,0.364676332,49907.270,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,p1024,1024,MODEXP,18200,0.220361636,82591.509,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,p1024,1024,MODEXP,18200,0.221022392,82344.598,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,p1024,1024,MODEXP,18200,0.363374891,50086.014,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,p1024,1024,MODEXP,18200,0.365363297,49813.433,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,p1024,1024,MODEXP,18200,0.224531260,81057.755,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,p1024,1024,MODEXP,18200,0.225935308,80554.032,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,p1024,1024,MODEXP,18200,0.000000000,inf,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,p1024,1024,MODEXP,18200,0.000000000,inf,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p1024,1024,EXPONENTIATION,195,0.010269205,18988.812,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p1024,1024,EXPONENTIATION,195,0.000856937,227554.675,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p1024,1024,EXPONENTIATION,195,0.013943926,13984.583,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,p1024,1024,EXPONENTIATION,18200,3.674644000,4952.861,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,p1024,1024,EXPONENTIATION,18200,3.679426773,4946.423,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,p1024,1024,EXPONENTIATION,18200,0.924042178,19696.071,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,p1024,1024,EXPONENTIATION,18200,0.931364002,19541.232,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,p1024,1024,EXPONENTIATION,18200,0.300942768,60476.615,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,p1024,1024,EXPONENTIATION,18200,0.302562956,60152.770,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,p1024,1024,EXPONENTIATION,18200,0.266424168,68312.121,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,p1024,1024,EXPONENTIATION,18200,0.267583164,68016.237,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,p1024,1024,EXPONENTIATION,18200,0.302169647,60231.066,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,p1024,1024,EXPONENTIATION,18200,0.302937767,60078.346,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,p1024,1024,EXPONENTIATION,18200,0.266180514,68374.652,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,p1024,1024,EXPONENTIATION,18200,0.267030001,68157.136,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,p1024,1024,EXPONENTIATION,18200,0.000000000,inf,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,p1024,1024,EXPONENTIATION,18200,0.000000000,inf,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p1024,1024,DIVIDE,1562,0.000089843,17385917.912,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p1024,1024,DIVIDE,1562,0.000018074,86421638.022,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p1024,1024,DIVIDE,1562,0.000037402,41762567.955,0
library,NVIDIA RTX A2000,gpu,cgbn,p1024,1024,DIVIDE,50000,0.000164864,303280279.503,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,p1024,1024,DIVIDE,18200,0.053167564,342313.973,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,p1024,1024,DIVIDE,18200,0.053249242,341788.904,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,p1024,1024,DIVIDE,18200,0.007358569,2473307.064,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,p1024,1024,DIVIDE,18200,0.008775426,2073973.371,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,p1024,1024,DIVIDE,18200,0.001680348,10831089.309,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,p1024,1024,DIVIDE,18200,0.002761696,6590152.187,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,p1024,1024,DIVIDE,18200,0.001647281,11048512.793,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,p1024,1024,DIVIDE,18200,0.002715805,6701513.027,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,p1024,1024,DIVIDE,18200,0.001676675,10854817.264,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,p1024,1024,DIVIDE,18200,0.002771087,6567819.798,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,p1024,1024,DIVIDE,18200,0.001660593,10959941.178,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,p1024,1024,DIVIDE,18200,0.002747861,6623334.201,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,p1024,1024,DIVIDE,18200,0.103533342,175788.780,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,p1024,1024,DIVIDE,18200,0.129260127,140801.348,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p1024,1024,ISQRT,390,0.000140870,2768510.170,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p1024,1024,ISQRT,390,0.000017965,21708621.636,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,p1024,1024,ISQRT,18200,0.921527977,19749.807,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,p1024,1024,ISQRT,18200,0.926639414,19640.865,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,p1024,1024,ISQRT,18200,0.285151768,63825.661,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,p1024,1024,ISQRT,18200,0.284533334,63964.386,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,p1024,1024,ISQRT,18200,0.055484316,328020.623,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,p1024,1024,ISQRT,18200,0.056265156,323468.401,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,p1024,1024,ISQRT,18200,0.052379630,347463.320,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,p1024,1024,ISQRT,18200,0.053277423,341608.116,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,p1024,1024,ISQRT,18200,0.055701789,326739.952,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,p1024,1024,ISQRT,18200,0.056535749,321920.208,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,p1024,1024,ISQRT,18200,0.051468747,353612.651,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,p1024,1024,ISQRT,18200,0.052359846,347594.602,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,p1024,1024,ISQRT,18200,2.250039481,8088.747,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,p1024,1024,ISQRT,18200,2.369217576,7681.861,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p1024,1024,MODMUL_R2,12500,0.007186130,1739462.067,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p1024,1024,MODMUL_R2,12500,0.000375831,33259585.772,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p1024,1024,MODMUL_R2,12500,0.000877758,14240834.390,0
library,NVIDIA RTX A2000,gpu,cgbn,p1024,1024,MODMUL_R2,50000,0.000313152,159666871.040,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,p1024,1024,MODMUL_R2,18200,0.005000762,3639445.560,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,p1024,1024,MODMUL_R2,18200,0.005917665,3075537.594,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,p1024,1024,MODMUL_R2,18200,0.000661978,27493378.799,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,p1024,1024,MODMUL_R2,18200,0.001475120,12337979.590,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,p1024,1024,MODMUL_R2,18200,0.000367860,49475302.217,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,p1024,1024,MODMUL_R2,18200,0.001200580,15159339.171,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,p1024,1024,MODMUL_R2,18200,0.000287666,63267820.722,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,p1024,1024,MODMUL_R2,18200,0.001124211,16189136.714,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,p1024,1024,MODMUL_R2,18200,0.000331229,54946946.180,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,p1024,1024,MODMUL_R2,18200,0.001215525,14972954.425,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,p1024,1024,MODMUL_R2,18200,0.000251335,72413286.435,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,p1024,1024,MODMUL_R2,18200,0.001129238,16117063.939,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,p1024,1024,MODMUL_R2,18200,0.085990354,211651.645,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,p1024,1024,MODMUL_R2,18200,0.114373728,159127.453,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p2048,2048,ADD,6250,0.000196919,31738963.299,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p2048,2048,ADD,6250,0.000032833,190358155.103,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p2048,2048,ADD,6250,0.000068351,91440182.039,0
library,NVIDIA RTX A2000,gpu,cgbn,p2048,2048,ADD,50000,0.000206784,241798204.890,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,p2048,2048,ADD,18200,0.000571108,31867854.266,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,p2048,2048,ADD,18200,0.002092673,8697012.980,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,p2048,2048,ADD,18200,0.000287870,63222994.639,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,p2048,2048,ADD,18200,0.001747872,10412663.291,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32,p2048,2048,ADD,18200,0.000205427,88595772.871,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32,p2048,2048,ADD,18200,0.001633557,11141328.114,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,p2048,2048,ADD,18200,0.000138628,131286327.917,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,p2048,2048,ADD,18200,0.001669787,10899594.626,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,p2048,2048,ADD,18200,0.000138989,130945879.714,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,p2048,2048,ADD,18200,0.001650107,11029587.156,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,p2048,2048,ADD,18200,0.000117487,154910394.660,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,p2048,2048,ADD,18200,0.001619044,11241203.264,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,p2048,2048,ADD,18200,0.000115152,158052628.914,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,p2048,2048,ADD,18200,0.001771267,10275133.222,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,p2048,2048,ADD,18200,0.001139346,15974082.261,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,p2048,2048,ADD,18200,0.003773071,4823657.074,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p2048,2048,SUBTRACT,6250,0.000163394,38251092.378,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p2048,2048,SUBTRACT,6250,0.000030640,203984510.167,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p2048,2048,SUBTRACT,6250,0.000056199,111212342.774,0
library,NVIDIA RTX A2000,gpu,cgbn,p2048,2048,SUBTRACT,50000,0.000206848,241723391.089,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,p2048,2048,SUBTRACT,18200,0.000568994,31986259.494,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,p2048,2048,SUBTRACT,18200,0.002085850,8725461.139,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,p2048,2048,SUBTRACT,18200,0.000289273,62916450.518,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,p2048,2048,SUBTRACT,18200,0.001785955,10190624.068,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32,p2048,2048,SUBTRACT,18200,0.000204986,88786568.031,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32,p2048,2048,SUBTRACT,18200,0.001643256,11075569.997,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,p2048,2048,SUBTRACT,18200,0.000137646,132223478.286,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,p2048,2048,SUBTRACT,18200,0.001673754,10873758.442,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,p2048,2048,SUBTRACT,18200,0.000138818,131106646.518,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,p2048,2048,SUBTRACT,18200,0.001632342,11149623.496,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,p2048,2048,SUBTRACT,18200,0.000117658,154686000.576,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,p2048,2048,SUBTRACT,18200,0.001621891,11221470.553,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,p2048,2048,SUBTRACT,18200,0.000115333,157803753.265,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,p2048,2048,SUBTRACT,18200,0.001767328,10298031.779,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,p2048,2048,SUBTRACT,18200,0.001182588,15389977.584,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,p2048,2048,SUBTRACT,18200,0.004120731,4416692.307,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p2048,2048,ADDMOD,6250,0.000555397,11253211.437,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p2048,2048,ADDMOD,6250,0.000084573,73900301.729,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p2048,2048,ADDMOD,6250,0.000280890,22250728.275,0
library,NVIDIA RTX A2000,gpu,cgbn,p2048,2048,ADDMOD,50000,0.000206848,241723391.089,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,p2048,2048,ADDMOD,18200,0.000739100,24624530.397,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,p2048,2048,ADDMOD,18200,0.002269122,8020724.076,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,p2048,2048,ADDMOD,18200,0.000384989,47274054.069,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,p2048,2048,ADDMOD,18200,0.001876862,9697037.956,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32,p2048,2048,ADDMOD,18200,0.000259402,70161314.887,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32,p2048,2048,ADDMOD,18200,0.001680418,10830639.098,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,p2048,2048,ADDMOD,18200,0.000186049,97823492.117,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,p2048,2048,ADDMOD,18200,0.001712610,10627054.527,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,p2048,2048,ADDMOD,18200,0.000190157,95710632.322,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,p2048,2048,ADDMOD,18200,0.001694722,10739223.846,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,p2048,2048,ADDMOD,18200,0.000086818,209634211.508,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,p2048,2048,ADDMOD,18200,0.001569348,11597169.976,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,p2048,2048,ADDMOD,18200,0.000087678,207576703.739,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,p2048,2048,ADDMOD,18200,0.001738353,10469682.009,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,p2048,2048,ADDMOD,18200,0.001336777,13614839.832,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,p2048,2048,ADDMOD,18200,0.004063951,4478400.617,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p2048,2048,SUBTRACTMOD,6250,0.000469232,13319637.817,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p2048,2048,SUBTRACTMOD,6250,0.000068331,91466354.096,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p2048,2048,SUBTRACTMOD,6250,0.000286060,21848533.802,0
library,NVIDIA RTX A2000,gpu,cgbn,p2048,2048,SUBTRACTMOD,50000,0.000205824,242925995.025,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,p2048,2048,SUBTRACTMOD,18200,0.000790379,23026927.292,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,p2048,2048,SUBTRACTMOD,18200,0.002268460,8023065.350,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,p2048,2048,SUBTRACTMOD,18200,0.000418615,43476744.883,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,p2048,2048,SUBTRACTMOD,18200,0.001914534,9506230.775,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32,p2048,2048,SUBTRACTMOD,18200,0.000292135,62299878.209,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32,p2048,2048,SUBTRACTMOD,18200,0.001728681,10528258.761,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,p2048,2048,SUBTRACTMOD,18200,0.000192652,94471092.231,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,p2048,2048,SUBTRACTMOD,18200,0.001711177,10635955.896,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,p2048,2048,SUBTRACTMOD,18200,0.000192261,94662835.980,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,p2048,2048,SUBTRACTMOD,18200,0.001704119,10680004.414,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,p2048,2048,SUBTRACTMOD,18200,0.000092418,196931476.392,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,p2048,2048,SUBTRACTMOD,18200,0.001488262,12229030.748,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,p2048,2048,SUBTRACTMOD,18200,0.000086849,209560026.989,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,p2048,2048,SUBTRACTMOD,18200,0.001750967,10394259.197,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,p2048,2048,SUBTRACTMOD,18200,0.001473349,12352813.279,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,p2048,2048,SUBTRACTMOD,18200,0.004399377,4136949.465,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.004023353,1553430.823,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.000468059,13353004.825,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.000571429,10937495.661,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,18200,0.182168176,99907.681,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,18200,0.186742948,97460.173,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,p2048,2048,MULTIPLYOPERANDSCANNING,18200,0.050315293,361719.047,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,p2048,2048,MULTIPLYOPERANDSCANNING,18200,0.051829498,351151.383,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32,p2048,2048,MULTIPLYOPERANDSCANNING,18200,0.014359866,1267421.274,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32,p2048,2048,MULTIPLYOPERANDSCANNING,18200,0.016332088,1114370.677,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,p2048,2048,MULTIPLYOPERANDSCANNING,18200,0.000780899,23306468.814,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,p2048,2048,MULTIPLYOPERANDSCANNING,18200,0.002739504,6643538.442,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,p2048,2048,MULTIPLYOPERANDSCANNING,18200,0.000751381,24222071.938,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,p2048,2048,MULTIPLYOPERANDSCANNING,18200,0.002713740,6706611.856,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,p2048,2048,MULTIPLYOPERANDSCANNING,18200,0.000683469,26628868.636,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,p2048,2048,MULTIPLYOPERANDSCANNING,18200,0.002613383,6964152.440,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,p2048,2048,MULTIPLYOPERANDSCANNING,18200,0.000569649,31949496.444,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,p2048,2048,MULTIPLYOPERANDSCANNING,18200,0.002687162,6772945.373,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,18200,0.127969082,142221.853,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,18200,0.108610491,167571.289,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.004028861,1551306.774,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.000471805,13246991.988,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.000574174,10885195.656,0
library,NVIDIA RTX A2000,gpu,cgbn,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.000367808,135940490.691,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,18200,0.058133662,313071.626,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,18200,0.060345564,301596.319,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,p2048,2048,MULTIPLYPRODUCTSCANNING,18200,0.015146850,1201569.989,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,p2048,2048,MULTIPLYPRODUCTSCANNING,18200,0.017132191,1062327.639,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32,p2048,2048,MULTIPLYPRODUCTSCANNING,18200,0.003740651,4865463.776,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32,p2048,2048,MULTIPLYPRODUCTSCANNING,18200,0.005688515,3199428.748,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,p2048,2048,MULTIPLYPRODUCTSCANNING,18200,0.003741470,4864398.003,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,p2048,2048,MULTIPLYPRODUCTSCANNING,18200,0.005704383,3190528.903,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,p2048,2048,MULTIPLYPRODUCTSCANNING,18200,0.003731814,4876984.480,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,p2048,2048,MULTIPLYPRODUCTSCANNING,18200,0.005750773,3164792.185,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,p2048,2048,MULTIPLYPRODUCTSCANNING,18200,0.002174797,8368598.944,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,p2048,2048,MULTIPLYPRODUCTSCANNING,18200,0.004112917,4425083.221,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,p2048,2048,MULTIPLYPRODUCTSCANNING,18200,0.002029474,8967841.400,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,p2048,2048,MULTIPLYPRODUCTSCANNING,18200,0.004131359,4405330.011,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,18200,0.096761057,188092.199,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,18200,0.094742900,192098.828,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.026246310,238128.710,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.002834708,2204812.295,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.001074815,5814953.287,0
library,NVIDIA RTX A2000,gpu,cgbn,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000780288,64078904.199,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,18200,0.112913789,161184.919,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,18200,0.114139454,159454.065,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,p2048,2048,MONTGOMERYMULTIPLICATION,18200,0.002790028,6523231.489,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,p2048,2048,MONTGOMERYMULTIPLICATION,18200,0.004467362,4073993.019,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32,p2048,2048,MONTGOMERYMULTIPLICATION,18200,0.000480108,37908140.250,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32,p2048,2048,MONTGOMERYMULTIPLICATION,18200,0.002003932,9082142.819,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,p2048,2048,MONTGOMERYMULTIPLICATION,18200,0.000724570,25118349.715,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,p2048,2048,MONTGOMERYMULTIPLICATION,18200,0.002240661,8122604.556,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,p2048,2048,MONTGOMERYMULTIPLICATION,18200,0.000602804,30192228.727,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,p2048,2048,MONTGOMERYMULTIPLICATION,18200,0.002114292,8608081.563,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,p2048,2048,MONTGOMERYMULTIPLICATION,18200,0.000685313,26557216.490,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,p2048,2048,MONTGOMERYMULTIPLICATION,18200,0.002149107,8468633.974,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,p2048,2048,MONTGOMERYMULTIPLICATION,18200,0.000571302,31857048.628,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,p2048,2048,MONTGOMERYMULTIPLICATION,18200,0.002467979,7374455.729,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,18200,0.479602229,37948.114,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,18200,0.526629282,34559.415,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p2048,2048,COMPARE,6250,0.000048704,128327495.937,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p2048,2048,COMPARE,6250,0.000021491,290816710.002,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p2048,2048,COMPARE,6250,0.000053071,117765840.133,0
library,NVIDIA RTX A2000,gpu,cgbn,p2048,2048,COMPARE,50000,0.000207008,241536558.974,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,p2048,2048,COMPARE,18200,0.000488218,37278435.002,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,p2048,2048,COMPARE,18200,0.001981679,9184129.611,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,p2048,2048,COMPARE,18200,0.000245938,74002367.506,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,p2048,2048,COMPARE,18200,0.001716409,10603534.057,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,p2048,2048,COMPARE,18200,0.000128478,141658701.554,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,p2048,2048,COMPARE,18200,0.001652462,11013866.284,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,p2048,2048,COMPARE,18200,0.000192432,94578995.445,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,p2048,2048,COMPARE,18200,0.001712726,10626332.198,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,p2048,2048,COMPARE,18200,0.000052092,349384105.927,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,p2048,2048,COMPARE,18200,0.001515163,12011908.088,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,p2048,2048,COMPARE,18200,0.000050658,359269426.716,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,p2048,2048,COMPARE,18200,0.001724215,10555526.794,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,p2048,2048,COMPARE,18200,0.000388248,46877266.723,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,p2048,2048,COMPARE,18200,0.002784365,6536499.682,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p2048,2048,REDUCE,781,0.000031981,24420989.678,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p2048,2048,REDUCE,781,0.000018025,43329149.765,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p2048,2048,REDUCE,781,0.000070516,11075497.445,0
library,NVIDIA RTX A2000,gpu,cgbn,p2048,2048,REDUCE,50000,0.000205824,242925995.025,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,p2048,2048,REDUCE,18200,0.678267646,26833.065,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,p2048,2048,REDUCE,18200,0.680593605,26741.362,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,p2048,2048,REDUCE,18200,0.012035996,1512130.813,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,p2048,2048,REDUCE,18200,0.014354226,1267919.271,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,p2048,2048,REDUCE,18200,0.003101534,5868063.599,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,p2048,2048,REDUCE,18200,0.004650136,3913864.425,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,p2048,2048,REDUCE,18200,0.003498455,5202297.069,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,p2048,2048,REDUCE,18200,0.005010906,3632077.947,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,p2048,2048,REDUCE,18200,0.002805673,6486858.121,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,p2048,2048,REDUCE,18200,0.004299527,4233023.739,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,p2048,2048,REDUCE,18200,0.003328599,5467765.809,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,p2048,2048,REDUCE,18200,0.005066341,3592336.315,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,p2048,2048,REDUCE,18200,0.224430951,81093.984,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,p2048,2048,REDUCE,18200,0.254392061,71543.113,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p2048,2048,MODMUL,390,0.000689255,565828.261,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p2048,2048,MODMUL,390,0.000091375,4268132.779,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p2048,2048,MODMUL,390,0.000235774,1654129.054,0
library,NVIDIA RTX A2000,gpu,cgbn,p2048,2048,MODMUL,50000,0.003918848,12758851.581,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,p2048,2048,MODMUL,18200,0.975081678,18665.103,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,p2048,2048,MODMUL,18200,0.977056061,18627.386,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,p2048,2048,MODMUL,18200,0.037735456,482305.028,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,p2048,2048,MODMUL,18200,0.039248398,463713.199,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,p2048,2048,MODMUL,18200,0.014393114,1264493.520,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,p2048,2048,MODMUL,18200,0.015737021,1156508.567,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,p2048,2048,MODMUL,18200,0.013599820,1338252.963,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,p2048,2048,MODMUL,18200,0.012784996,1423543.652,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,p2048,2048,MODMUL,18200,0.012450744,1461759.989,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,p2048,2048,MODMUL,18200,0.013913426,1308089.072,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,p2048,2048,MODMUL,18200,0.013952360,1304438.860,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,p2048,2048,MODMUL,18200,0.015323815,1187693.801,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,p2048,2048,MODMUL,18200,1.143943502,15909.877,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,p2048,2048,MODMUL,18200,1.151127252,15810.589,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p2048,2048,MODEXP,97,0.255280212,379.975,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p2048,2048,MODEXP,97,0.022708601,4271.509,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p2048,2048,MODEXP,97,0.016810416,5770.232,0
library,NVIDIA RTX A2000,gpu,cgbn,p2048,2048,MODEXP,50000,5.265794754,9495.243,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,p2048,2048,MODEXP,18200,161.415891325,112.752,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,p2048,2048,MODEXP,18200,161.383357016,112.775,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,p2048,2048,MODEXP,18200,18.991127291,958.342,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,p2048,2048,MODEXP,18200,19.086419251,953.558,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,p2048,2048,MODEXP,18200,2.973823787,6120.067,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,p2048,2048,MODEXP,18200,2.985122719,6096.902,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,p2048,2048,MODEXP,18200,4.517140778,4029.097,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,p2048,2048,MODEXP,18200,4.547381489,4002.303,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,p2048,2048,MODEXP,18200,2.765455129,6581.195,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,p2048,2048,MODEXP,18200,2.772803010,6563.755,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,p2048,2048,MODEXP,18200,4.519500742,4026.993,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,p2048,2048,MODEXP,18200,4.555589012,3995.093,0
opencl-kernel,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,p2048,2048,MODEXP,18200,0.000000000,inf,0
opencl-e2e,cpu-haswell-AMD EPYC 7282 16-Core Processor,CPU,w8,p2048,2048,MODEXP,18200,0.000000000,inf,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p2048,2048,EXPONENTIATION,97,0.030550611,3175.059,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p2048,2048,EXPONENTIATION,97,0.002298367,42203.874,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p2048,2048,EXPONENTIATION,97,0.031167347,3112.232,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,p2048,2048,EXPONENTIATION,18200,46.846499763,388.503,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,p2048,2048,EXPONENTIATION,18200,46.925725558,387.847,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,p2048,2048,EXPONENTIATION,18200,8.749359471,2080.152,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,p2048,2048,EXPONENTIATION,18200,8.757405303,2078.241,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,p2048,2048,EXPONENTIATION,18200,2.481098787,7335.460,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,p2048,2048,EXPONENTIATION,18200,2.491269387,7305.513,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,p2048,2048,EXPONENTIATION,18200,2.256218618,8066.594,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,p2048,2048,EXPONENTIATION,18200,2.258521968,8058.367,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,p2048,2048,EXPONENTIATION,18200,2.493382022,7299.323,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,p2048,2048,EXPONENTIATION,18200,2.505407706,7264.287,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,p2048,2048,EXPONENTIATION,18200,2.273454266,8005.439,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,p2048,2048,EXPONENTIATION,18200,2.276742992,7993.875,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p2048,2048,DIVIDE,781,0.000057441,13596516.765,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p2048,2048,DIVIDE,781,0.000016212,48175582.498,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p2048,2048,DIVIDE,781,0.000069333,11264438.177,0
library,NVIDIA RTX A2000,gpu,cgbn,p2048,2048,DIVIDE,50000,0.000203776,245367462.312,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,p2048,2048,DIVIDE,18200,1.758353911,10350.590,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,p2048,2048,DIVIDE,18200,1.757527732,10355.455,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,p2048,2048,DIVIDE,18200,0.338556357,53757.667,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,p2048,2048,DIVIDE,18200,0.341257147,53332.216,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,p2048,2048,DIVIDE,18200,0.018426754,987694.328,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,p2048,2048,DIVIDE,18200,0.020695435,879420.998,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,p2048,2048,DIVIDE,18200,0.016695458,1090116.877,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,p2048,2048,DIVIDE,18200,0.018936392,961112.343,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,p2048,2048,DIVIDE,18200,0.017897099,1016924.572,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,p2048,2048,DIVIDE,18200,0.020044643,907973.290,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,p2048,2048,DIVIDE,18200,0.016023595,1135825.021,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,p2048,2048,DIVIDE,18200,0.018653536,975686.314,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p2048,2048,ISQRT,195,0.000109731,1777069.466,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p2048,2048,ISQRT,195,0.000017364,11230404.188,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,p2048,2048,ISQRT,18200,23.219770738,783.815,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,p2048,2048,ISQRT,18200,23.225868860,783.609,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,p2048,2048,ISQRT,18200,8.859234477,2054.354,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,p2048,2048,ISQRT,18200,8.874127184,2050.906,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,p2048,2048,ISQRT,18200,0.350448309,51933.479,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,p2048,2048,ISQRT,18200,0.351582536,51765.939,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,p2048,2048,ISQRT,18200,0.113963326,159700.498,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,p2048,2048,ISQRT,18200,0.115368688,157755.110,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,p2048,2048,ISQRT,18200,0.350810353,51879.883,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,p2048,2048,ISQRT,18200,0.352459587,51637.126,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,p2048,2048,ISQRT,18200,0.113732295,160024.907,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,p2048,2048,ISQRT,18200,0.115181538,158011.434,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p2048,2048,MODMUL_R2,6250,0.011287872,553691.620,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p2048,2048,MODMUL_R2,6250,0.000932993,6698868.630,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p2048,2048,MODMUL_R2,6250,0.002292566,2726202.864,0
library,NVIDIA RTX A2000,gpu,cgbn,p2048,2048,MODMUL_R2,50000,0.001053696,47452016.521,0
opencl-kernel,NVIDIA RTX A2000,GPU,w8,p2048,2048,MODMUL_R2,18200,0.074363585,244743.447,0
opencl-e2e,NVIDIA RTX A2000,GPU,w8,p2048,2048,MODMUL_R2,18200,0.077464138,234947.428,0
opencl-kernel,NVIDIA RTX A2000,GPU,w16,p2048,2048,MODMUL_R2,18200,0.003235064,5625855.540,0
opencl-e2e,NVIDIA RTX A2000,GPU,w16,p2048,2048,MODMUL_R2,18200,0.004796912,3794107.572,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-opt,p2048,2048,MODMUL_R2,18200,0.001343704,13544651.441,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-opt,p2048,2048,MODMUL_R2,18200,0.002814366,6466821.623,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-o64,p2048,2048,MODMUL_R2,18200,0.001035227,17580677.725,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-o64,p2048,2048,MODMUL_R2,18200,0.002784122,6537070.368,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il,p2048,2048,MODMUL_R2,18200,0.001247126,14593557.269,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il,p2048,2048,MODMUL_R2,18200,0.002923810,6224754.850,0
opencl-kernel,NVIDIA RTX A2000,GPU,w32-il64,p2048,2048,MODMUL_R2,18200,0.001001842,18166531.436,0
opencl-e2e,NVIDIA RTX A2000,GPU,w32-il64,p2048,2048,MODMUL_R2,18200,0.002524157,7210327.886,0
```
