# MPA-OpenCL benchmark report - NVIDIA GeForce GTX 1660 SUPER


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

### Device 0 - NVIDIA GeForce GTX 1660 SUPER (GPU)

| Property | Value |
|---|---|
| Model | NVIDIA GeForce GTX 1660 SUPER |
| Type | GPU |
| Vendor | NVIDIA Corporation |
| Device memory | 5.61 GiB |
| Max single allocation | 1.40 GiB |
| Local memory | 48 KiB |
| Global cache | 704 KiB |
| Compute units | 22 |
| Max clock | 1830 MHz |
| Max work-group size | 1024 |
| OpenCL version | OpenCL 3.0 CUDA |
| Driver | 580.173.02 |

### Device 1 - cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz (CPU)

| Property | Value |
|---|---|
| Model | cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz |
| Type | CPU |
| Vendor | GenuineIntel |
| Device memory | 13.56 GiB |
| Max single allocation | 4.00 GiB |
| Local memory | 256 KiB |
| Global cache | 3072 KiB |
| Compute units | 4 |
| Max clock | 3300 MHz |
| Max work-group size | 4096 |
| OpenCL version | OpenCL 3.0 PoCL HSTR: cpu-x86_64-pc-linux-gnu-sandybridge |
| Driver | 5.0+debian |

### Host

| Property | Value |
|---|---|
| CPU | Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz |
| Logical cores | 4 |
| OpenMP threads used | 4 |
| RAM | 15.6 GB |
| OS | Ubuntu 24.04.4 LTS |
| Kernel | 6.8.0-138-generic |
| Arch | x86_64 |
| GMP | 6.3.0 |
| OpenSSL | OpenSSL 3.0.13 30 Jan 2024 |

## 2. Method

- Workload auto-sized from the device and host: --min-items from 700 x compute units, --items from ten times that capped by host RAM. Either flag, given explicitly, overrides its half.
- Base workload 50000 items, scaled down per operator by its cost weight and by modulus size. Device rows honour --min-items (15400) so the GPU is not left idle; the CPU libraries keep the smaller count because a full-width MODEXP there costs minutes. Both counts appear in every row as dev/cpu, and throughput is per-second so they remain comparable.
- 5 timed repetitions, **minimum** reported. Two untimed warm-up launches precede them.
- `kernel` times `clEnqueueNDRangeKernel` + `clFinish` only. `e2e` adds the host->device operand writes and the device->host result read.
- Every OpenCL device runs the same kernels on the same operands, so GPU and CPU-OpenCL columns are directly comparable.
- CPU library baselines (GMP, OpenSSL) run those same operands, with temporaries preallocated outside the timed region, so the figure is the arithmetic and not marshalling. The generator is reseeded per modulus and operation so every backend sees identical inputs.
- OpenSSL rows time the nearest BN primitive, which is not always semantically identical (its Montgomery routine expects Montgomery-domain inputs); they measure comparable work, not identical results. Correctness is judged against GMP only.
- Every device cell is checked word-for-word against GMP before it is timed. A cell that mismatches is reported and excluded from the speedup tables.
- Total wall time 5407.7 s.

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
| [1] CPU | `mpaKernels_8bits.cl` (w8) | 26 | 26 | 0 | 0 |

**All configurations correct** - 511 configurations, 0 problems.

## 4. Throughput per device

Operations per second, higher is better. Kernel-only timings.

### Device 0 - NVIDIA GeForce GTX 1660 SUPER (GPU)

#### secp256k1 (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | 325.22 M | 696.50 M | 445.35 M | 1.07 G | 442.01 M | 1.72 G | 1.73 G | 37.20 M |
| SUBTRACT | 50000 / 50000 | 327.16 M | 690.05 M | 467.75 M | 1.07 G | 457.77 M | 1.71 G | 1.74 G | 40.04 M |
| ADDMOD | 50000 / 50000 | 214.17 M | 446.96 M | 510.44 M | 1.51 G | 1.45 G | 1.87 G | 1.87 G | 12.17 M |
| SUBTRACTMOD | 50000 / 50000 | 214.97 M | 487.09 M | 488.82 M | 1.59 G | 1.47 G | 1.85 G | 1.86 G | 16.09 M |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | 7.60 M | 31.71 M | 46.69 M | 619.77 M | 631.31 M | 1.40 G | 1.41 G | 25.47 M |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | 27.29 M | 112.61 M | 137.51 M | 306.49 M | 261.64 M | 702.71 M | 713.96 M | 25.45 M |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | 108.26 M | 557.10 M | 1.45 G | 1.14 G | 1.08 G | 1.11 G | 1.12 G | 3.82 M |
| COMPARE | 50000 / 50000 | 346.62 M | 742.33 M | - | 1.48 G | 1.23 G | 2.79 G | 2.93 G | 65.51 M |
| REDUCE | 15400 / 6250 | 46.44 M | 92.27 M | - | 244.67 M | 260.49 M | 254.56 M | 274.34 M | 38.07 M |
| MODMUL | 15400 / 3125 | 16.54 M | 37.53 M | - | 85.85 M | 104.18 M | 85.97 M | 104.77 M | 6.98 M |
| MODEXP | 15400 / 781 | 431.73 k | 2.55 M | - | 2.95 M | 6.86 M | 3.64 M | 6.74 M | 68.20 k |
| EXPONENTIATION | 15400 / 781 | 326.89 k | 1.24 M | - | 16.15 M | 23.19 M | 19.76 M | 21.62 M | 213.05 k |
| DIVIDE | 15400 / 6250 | 27.95 M | 49.69 M | - | 104.03 M | 131.00 M | 128.16 M | 132.03 M | 17.02 M |
| ISQRT | 15400 / 1562 | 2.22 M | 3.11 M | - | 10.60 M | 15.46 M | 13.39 M | 15.29 M | 8.65 M |
| MODMUL_R2 | 50000 / 50000 | 111.69 M | 518.58 M | - | 556.09 M | 781.64 M | 653.25 M | 810.02 M | 6.67 M |

#### rsa256(composite) (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | 403.00 M | 699.44 M | 459.74 M | 914.60 M | 452.23 M | 1.71 G | 1.73 G | 36.64 M |
| SUBTRACT | 50000 / 50000 | 399.94 M | 690.67 M | 447.68 M | 929.28 M | 458.57 M | 1.74 G | 1.76 G | 39.96 M |
| ADDMOD | 50000 / 50000 | 290.10 M | 537.05 M | 516.48 M | 1.28 G | 1.47 G | 1.75 G | 1.86 G | 14.72 M |
| SUBTRACTMOD | 50000 / 50000 | 265.27 M | 489.94 M | 469.49 M | 1.30 G | 1.50 G | 1.83 G | 1.85 G | 16.14 M |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | 9.44 M | 31.87 M | 46.11 M | 515.71 M | 626.14 M | 1.37 G | 1.45 G | 25.28 M |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | 32.69 M | 113.32 M | 139.78 M | 249.68 M | 314.09 M | 611.28 M | 637.23 M | 25.40 M |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | 128.60 M | 567.64 M | 1.43 G | 1.09 G | 1.10 G | 1.04 G | 1.10 G | 3.82 M |
| COMPARE | 50000 / 50000 | 409.69 M | 718.78 M | - | 1.24 G | 1.05 G | 2.12 G | 2.82 G | 65.61 M |
| REDUCE | 15400 / 6250 | 52.31 M | 92.18 M | - | 204.07 M | 217.34 M | 202.54 M | 231.22 M | 22.91 M |
| MODMUL | 15400 / 3125 | 18.74 M | 37.73 M | - | 70.43 M | 85.99 M | 69.62 M | 86.03 M | 6.95 M |
| MODEXP | 15400 / 781 | 430.80 k | 2.57 M | - | 2.98 M | 5.53 M | 2.98 M | 5.43 M | 70.62 k |
| EXPONENTIATION | 15400 / 781 | 325.07 k | 1.24 M | - | 16.14 M | 18.87 M | 16.09 M | 18.09 M | 213.87 k |
| DIVIDE | 15400 / 6250 | 28.01 M | 49.29 M | - | 101.85 M | 105.60 M | 104.26 M | 106.31 M | 16.98 M |
| ISQRT | 15400 / 1562 | 2.22 M | 3.11 M | - | 10.50 M | 12.44 M | 10.79 M | 12.33 M | 8.76 M |
| MODMUL_R2 | 50000 / 50000 | 111.80 M | 514.51 M | - | 550.52 M | 731.41 M | 633.33 M | 770.84 M | 6.68 M |

#### brainpoolP512r1 (512-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 25000 / 25000 | 139.89 M | 257.80 M | 94.40 M | 116.33 M | 113.56 M | 512.48 M | 526.19 M | 31.65 M |
| SUBTRACT | 25000 / 25000 | 139.38 M | 247.57 M | 94.34 M | 117.52 M | 114.78 M | 522.99 M | 520.13 M | 33.52 M |
| ADDMOD | 25000 / 25000 | 103.82 M | 203.16 M | 98.59 M | 298.88 M | 293.13 M | 752.04 M | 725.12 M | 13.15 M |
| SUBTRACTMOD | 25000 / 25000 | 91.70 M | 185.86 M | 96.06 M | 294.19 M | 290.10 M | 702.52 M | 720.03 M | 13.94 M |
| MULTIPLYOPERANDSCANNING | 25000 / 25000 | 1.27 M | 4.95 M | 4.43 M | 111.07 M | 111.31 M | 217.37 M | 219.16 M | 11.02 M |
| MULTIPLYPRODUCTSCANNING | 25000 / 25000 | 4.41 M | 16.76 M | 22.33 M | 24.60 M | 24.00 M | 112.36 M | 118.58 M | 11.05 M |
| MONTGOMERYMULTIPLICATION | 25000 / 25000 | 31.96 M | 111.43 M | 345.33 M | 236.95 M | 328.67 M | 268.56 M | 391.08 M | 1.64 M |
| COMPARE | 25000 / 25000 | 161.65 M | 297.14 M | - | 421.07 M | 434.16 M | 1.27 G | 1.26 G | 54.09 M |
| REDUCE | 15400 / 3125 | 17.87 M | 22.80 M | - | 76.74 M | 79.20 M | 68.08 M | 69.15 M | 21.82 M |
| MODMUL | 15400 / 1562 | 6.29 M | 8.86 M | - | 20.06 M | 24.94 M | 15.51 M | 19.38 M | 3.46 M |
| MODEXP | 15400 / 390 | 30.01 k | 329.75 k | - | 471.94 k | 765.12 k | 401.80 k | 740.87 k | 11.60 k |
| EXPONENTIATION | 15400 / 390 | 39.49 k | 154.06 k | - | 535.07 k | 558.32 k | 506.50 k | 537.06 k | 56.17 k |
| DIVIDE | 15400 / 3125 | 9.23 M | 11.77 M | - | 35.17 M | 36.51 M | 42.20 M | 39.83 M | 15.51 M |
| ISQRT | 15400 / 781 | 398.40 k | 430.29 k | - | 2.30 M | 2.54 M | 1.79 M | 1.94 M | 4.98 M |
| MODMUL_R2 | 25000 / 25000 | 21.91 M | 104.99 M | - | 179.45 M | 175.14 M | 178.25 M | 170.59 M | 3.33 M |

#### p1024 (1024-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 15400 / 12500 | 20.69 M | 42.35 M | 41.42 M | 79.55 M | 79.92 M | 258.48 M | 255.55 M | 22.34 M |
| SUBTRACT | 15400 / 12500 | 20.67 M | 40.96 M | 42.71 M | 81.25 M | 79.41 M | 262.90 M | 251.42 M | 23.01 M |
| ADDMOD | 15400 / 12500 | 17.04 M | 37.05 M | 40.13 M | 113.69 M | 112.70 M | 362.24 M | 349.15 M | 8.84 M |
| SUBTRACTMOD | 15400 / 12500 | 17.84 M | 35.60 M | 40.43 M | 113.97 M | 113.20 M | 365.38 M | 348.38 M | 11.25 M |
| MULTIPLYOPERANDSCANNING | 15400 / 12500 | 63.36 k | 228.60 k | 1.09 M | 50.54 M | 51.06 M | 72.35 M | 77.44 M | 3.39 M |
| MULTIPLYPRODUCTSCANNING | 15400 / 12500 | 571.50 k | 2.26 M | 6.90 M | 8.69 M | 8.75 M | 27.06 M | 28.83 M | 3.39 M |
| MONTGOMERYMULTIPLICATION | 15400 / 12500 | 4.98 M | 24.41 M | 87.91 M | 59.43 M | 73.78 M | 81.28 M | 119.87 M | 542.15 k |
| COMPARE | 15400 / 12500 | 60.45 M | 112.45 M | - | 204.42 M | 203.30 M | 664.39 M | 647.22 M | 33.01 M |
| REDUCE | 15400 / 1562 | 3.32 M | 6.11 M | - | 25.51 M | 25.33 M | 27.67 M | 28.26 M | 30.82 M |
| MODMUL | 15400 / 781 | 840.22 k | 2.28 M | - | 5.30 M | 6.85 M | 5.43 M | 7.24 M | 1.25 M |
| MODEXP | 15400 / 195 | 3.59 k | 24.41 k | - | 52.40 k | 85.80 k | 51.97 k | 85.28 k | 1.76 k |
| EXPONENTIATION | 15400 / 195 | 4.19 k | 19.21 k | - | 61.78 k | 60.25 k | 63.05 k | 61.48 k | 13.02 k |
| DIVIDE | 15400 / 1562 | 151.19 k | 1.20 M | - | 8.70 M | 9.32 M | 8.71 M | 9.13 M | 14.12 M |
| ISQRT | 15400 / 390 | 11.64 k | 55.84 k | - | 321.33 k | 337.91 k | 316.04 k | 341.21 k | 2.53 M |
| MODMUL_R2 | 15400 / 12500 | 3.20 M | 21.44 M | - | 36.93 M | 49.83 M | 44.06 M | 61.43 M | 1.21 M |

#### p2048 (2048-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 15400 / 6250 | 8.31 M | 18.75 M | 35.48 M | 33.95 M | 33.69 M | 131.33 M | 141.05 M | 13.81 M |
| SUBTRACT | 15400 / 6250 | 8.43 M | 18.57 M | 35.13 M | 33.93 M | 34.16 M | 140.79 M | 141.81 M | 13.76 M |
| ADDMOD | 15400 / 6250 | 8.15 M | 16.85 M | 34.35 M | 59.95 M | 59.44 M | 150.29 M | 146.38 M | 6.93 M |
| SUBTRACTMOD | 15400 / 6250 | 7.13 M | 16.77 M | 33.20 M | 59.34 M | 59.09 M | 148.79 M | 153.18 M | 8.03 M |
| MULTIPLYOPERANDSCANNING | 15400 / 6250 | 15.08 k | 54.51 k | 287.54 k | 10.54 M | 11.60 M | 12.17 M | 13.20 M | 1.06 M |
| MULTIPLYPRODUCTSCANNING | 15400 / 6250 | 142.59 k | 568.89 k | 2.24 M | 2.23 M | 2.23 M | 4.23 M | 4.40 M | 1.06 M |
| MONTGOMERYMULTIPLICATION | 15400 / 6250 | 148.76 k | 5.80 M | 32.66 M | 18.87 M | 23.27 M | 17.46 M | 23.31 M | 166.44 k |
| COMPARE | 15400 / 6250 | 29.68 M | 59.97 M | - | 112.83 M | 115.73 M | 342.49 M | 347.68 M | 39.00 M |
| REDUCE | 15400 / 781 | 19.13 k | 1.53 M | - | 6.90 M | 7.56 M | 6.55 M | 6.54 M | 20.53 M |
| MODMUL | 15400 / 390 | 11.21 k | 510.89 k | - | 1.28 M | 1.70 M | 1.07 M | 1.39 M | 407.91 k |
| MODEXP | 15400 / 97 | 81.4 | 607.3 | - | 5.76 k | 2.11 k | 5.74 k | 2.55 k | 249.9 |
| EXPONENTIATION | 15400 / 97 | 332.9 | 2.01 k | - | 6.87 k | 7.11 k | 7.00 k | 7.09 k | 2.14 k |
| DIVIDE | 15400 / 781 | 8.20 k | 39.17 k | - | 443.12 k | 478.26 k | 469.17 k | 483.05 k | 11.33 M |
| ISQRT | 15400 / 195 | 618.0 | 1.70 k | - | 130.86 k | 151.43 k | 130.98 k | 149.96 k | 1.14 M |
| MODMUL_R2 | 15400 / 6250 | 206.25 k | 5.47 M | - | 10.86 M | 14.46 M | 11.83 M | 15.98 M | 391.49 k |

### Device 1 - cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz (CPU)

#### secp256k1 (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | 53.21 M | - | - | - | - | - | - | 37.20 M |
| SUBTRACT | 50000 / 50000 | 41.17 M | - | - | - | - | - | - | 40.04 M |
| ADDMOD | 50000 / 50000 | 28.93 M | - | - | - | - | - | - | 12.17 M |
| SUBTRACTMOD | 50000 / 50000 | 33.55 M | - | - | - | - | - | - | 16.09 M |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | 1.49 M | - | - | - | - | - | - | 25.47 M |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | 2.01 M | - | - | - | - | - | - | 25.45 M |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | 433.75 k | - | - | - | - | - | - | 3.82 M |
| COMPARE | 50000 / 50000 | 100.72 M | - | - | - | - | - | - | 65.51 M |
| REDUCE | 15400 / 6250 | 442.31 k | - | - | - | - | - | - | 38.07 M |
| MODMUL | 15400 / 3125 | 122.56 k | - | - | - | - | - | - | 6.98 M |
| MODEXP | 15400 / 781 | 1.87 k | - | - | - | - | - | - | 68.20 k |
| EXPONENTIATION | 15400 / 781 | 3.79 k | - | - | - | - | - | - | 213.05 k |
| DIVIDE | 15400 / 6250 | 279.62 k | - | - | - | - | - | - | 17.02 M |
| ISQRT | 15400 / 1562 | 20.07 k | - | - | - | - | - | - | 8.65 M |
| MODMUL_R2 | 50000 / 50000 | 448.86 k | - | - | - | - | - | - | 6.67 M |

#### rsa256(composite) (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | 51.76 M | - | - | - | - | - | - | 36.64 M |
| SUBTRACT | 50000 / 50000 | 37.23 M | - | - | - | - | - | - | 39.96 M |
| ADDMOD | 50000 / 50000 | 38.44 M | - | - | - | - | - | - | 14.72 M |
| SUBTRACTMOD | 50000 / 50000 | 32.61 M | - | - | - | - | - | - | 16.14 M |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | 1.49 M | - | - | - | - | - | - | 25.28 M |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | 2.01 M | - | - | - | - | - | - | 25.40 M |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | 435.07 k | - | - | - | - | - | - | 3.82 M |
| COMPARE | 50000 / 50000 | 99.53 M | - | - | - | - | - | - | 65.61 M |
| REDUCE | 15400 / 6250 | 445.03 k | - | - | - | - | - | - | 22.91 M |
| MODMUL | 15400 / 3125 | over budget | - | - | - | - | - | - | 6.95 M |
| MODEXP | 15400 / 781 | over budget | - | - | - | - | - | - | 70.62 k |
| EXPONENTIATION | 15400 / 781 | - | - | - | - | - | - | - | 213.87 k |
| DIVIDE | 15400 / 6250 | - | - | - | - | - | - | - | 16.98 M |
| ISQRT | 15400 / 1562 | - | - | - | - | - | - | - | 8.76 M |
| MODMUL_R2 | 50000 / 50000 | - | - | - | - | - | - | - | 6.68 M |

#### brainpoolP512r1 (512-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 25000 / 25000 | - | - | - | - | - | - | - | 31.65 M |
| SUBTRACT | 25000 / 25000 | - | - | - | - | - | - | - | 33.52 M |
| ADDMOD | 25000 / 25000 | - | - | - | - | - | - | - | 13.15 M |
| SUBTRACTMOD | 25000 / 25000 | - | - | - | - | - | - | - | 13.94 M |
| MULTIPLYOPERANDSCANNING | 25000 / 25000 | - | - | - | - | - | - | - | 11.02 M |
| MULTIPLYPRODUCTSCANNING | 25000 / 25000 | - | - | - | - | - | - | - | 11.05 M |
| MONTGOMERYMULTIPLICATION | 25000 / 25000 | - | - | - | - | - | - | - | 1.64 M |
| COMPARE | 25000 / 25000 | - | - | - | - | - | - | - | 54.09 M |
| REDUCE | 15400 / 3125 | - | - | - | - | - | - | - | 21.82 M |
| MODMUL | 15400 / 1562 | - | - | - | - | - | - | - | 3.46 M |
| MODEXP | 15400 / 390 | - | - | - | - | - | - | - | 11.60 k |
| EXPONENTIATION | 15400 / 390 | - | - | - | - | - | - | - | 56.17 k |
| DIVIDE | 15400 / 3125 | - | - | - | - | - | - | - | 15.51 M |
| ISQRT | 15400 / 781 | - | - | - | - | - | - | - | 4.98 M |
| MODMUL_R2 | 25000 / 25000 | - | - | - | - | - | - | - | 3.33 M |

#### p1024 (1024-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 15400 / 12500 | - | - | - | - | - | - | - | 22.34 M |
| SUBTRACT | 15400 / 12500 | - | - | - | - | - | - | - | 23.01 M |
| ADDMOD | 15400 / 12500 | - | - | - | - | - | - | - | 8.84 M |
| SUBTRACTMOD | 15400 / 12500 | - | - | - | - | - | - | - | 11.25 M |
| MULTIPLYOPERANDSCANNING | 15400 / 12500 | - | - | - | - | - | - | - | 3.39 M |
| MULTIPLYPRODUCTSCANNING | 15400 / 12500 | - | - | - | - | - | - | - | 3.39 M |
| MONTGOMERYMULTIPLICATION | 15400 / 12500 | - | - | - | - | - | - | - | 542.15 k |
| COMPARE | 15400 / 12500 | - | - | - | - | - | - | - | 33.01 M |
| REDUCE | 15400 / 1562 | - | - | - | - | - | - | - | 30.82 M |
| MODMUL | 15400 / 781 | - | - | - | - | - | - | - | 1.25 M |
| MODEXP | 15400 / 195 | - | - | - | - | - | - | - | 1.76 k |
| EXPONENTIATION | 15400 / 195 | - | - | - | - | - | - | - | 13.02 k |
| DIVIDE | 15400 / 1562 | - | - | - | - | - | - | - | 14.12 M |
| ISQRT | 15400 / 390 | - | - | - | - | - | - | - | 2.53 M |
| MODMUL_R2 | 15400 / 12500 | - | - | - | - | - | - | - | 1.21 M |

#### p2048 (2048-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 15400 / 6250 | - | - | - | - | - | - | - | 13.81 M |
| SUBTRACT | 15400 / 6250 | - | - | - | - | - | - | - | 13.76 M |
| ADDMOD | 15400 / 6250 | - | - | - | - | - | - | - | 6.93 M |
| SUBTRACTMOD | 15400 / 6250 | - | - | - | - | - | - | - | 8.03 M |
| MULTIPLYOPERANDSCANNING | 15400 / 6250 | - | - | - | - | - | - | - | 1.06 M |
| MULTIPLYPRODUCTSCANNING | 15400 / 6250 | - | - | - | - | - | - | - | 1.06 M |
| MONTGOMERYMULTIPLICATION | 15400 / 6250 | - | - | - | - | - | - | - | 166.44 k |
| COMPARE | 15400 / 6250 | - | - | - | - | - | - | - | 39.00 M |
| REDUCE | 15400 / 781 | - | - | - | - | - | - | - | 20.53 M |
| MODMUL | 15400 / 390 | - | - | - | - | - | - | - | 407.91 k |
| MODEXP | 15400 / 97 | - | - | - | - | - | - | - | 249.9 |
| EXPONENTIATION | 15400 / 97 | - | - | - | - | - | - | - | 2.14 k |
| DIVIDE | 15400 / 781 | - | - | - | - | - | - | - | 11.33 M |
| ISQRT | 15400 / 195 | - | - | - | - | - | - | - | 1.14 M |
| MODMUL_R2 | 15400 / 6250 | - | - | - | - | - | - | - | 391.49 k |

## 5. Head to head

Best OpenCL GPU result against best OpenCL CPU result and the CPU libraries.
Ratios above 1.00x mean the GPU is faster than that baseline.

### secp256k1 (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w32-il64 | 1.73 G | w8 | 53.21 M | 37.20 M | 32.44x |
| SUBTRACT | w32-il64 | 1.74 G | w8 | 41.17 M | 40.04 M | 42.35x |
| ADDMOD | w32-il | 1.87 G | w8 | 28.93 M | 12.17 M | 64.58x |
| SUBTRACTMOD | w32-il64 | 1.86 G | w8 | 33.55 M | 16.09 M | 55.58x |
| MULTIPLYOPERANDSCANNING | w32-il64 | 1.41 G | w8 | 1.49 M | 25.47 M | 945.30x |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 713.96 M | w8 | 2.01 M | 25.45 M | 355.47x |
| MONTGOMERYMULTIPLICATION | w32 | 1.45 G | w8 | 433.75 k | 3.82 M | 3351.75x |
| COMPARE | w32-il64 | 2.93 G | w8 | 100.72 M | 65.51 M | 29.14x |
| REDUCE | w32-il64 | 111.34 M | w8 | 179.51 k | 38.07 M | 620.24x |
| MODMUL | w32-il64 | 21.26 M | w8 | 24.87 k | 6.98 M | 854.80x |
| MODEXP | w32-o64 | 347.73 k | w8 | 94.6 | 68.20 k | 3676.25x |
| EXPONENTIATION | w32-o64 | 1.18 M | w8 | 192.0 | 213.05 k | 6125.93x |
| DIVIDE | w32-il64 | 53.58 M | w8 | 113.48 k | 17.02 M | 472.18x |
| ISQRT | w32-o64 | 1.57 M | w8 | 2.04 k | 8.65 M | 770.23x |
| MODMUL_R2 | w32-il64 | 810.02 M | w8 | 448.86 k | 6.67 M | 1804.63x |

### rsa256(composite) (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w32-il64 | 1.73 G | w8 | 51.76 M | 36.64 M | 33.34x |
| SUBTRACT | w32-il64 | 1.76 G | w8 | 37.23 M | 39.96 M | 47.20x |
| ADDMOD | w32-il64 | 1.86 G | w8 | 38.44 M | 14.72 M | 48.44x |
| SUBTRACTMOD | w32-il64 | 1.85 G | w8 | 32.61 M | 16.14 M | 56.76x |
| MULTIPLYOPERANDSCANNING | w32-il64 | 1.45 G | w8 | 1.49 M | 25.28 M | 970.82x |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 637.23 M | w8 | 2.01 M | 25.40 M | 317.46x |
| MONTGOMERYMULTIPLICATION | w32 | 1.43 G | w8 | 435.07 k | 3.82 M | 3278.83x |
| COMPARE | w32-il64 | 2.82 G | w8 | 99.53 M | 65.61 M | 28.30x |
| REDUCE | w32-il64 | 93.84 M | w8 | 180.61 k | 22.91 M | 519.55x |
| MODMUL | w32-il64 | 17.46 M | none | n/a | 6.95 M | n/a |
| MODEXP | w32-o64 | 280.35 k | none | n/a | 70.62 k | n/a |
| EXPONENTIATION | w32-o64 | 956.91 k | none | n/a | 213.87 k | n/a |
| DIVIDE | w32-il64 | 43.14 M | none | n/a | 16.98 M | n/a |
| ISQRT | w32-o64 | 1.26 M | none | n/a | 8.76 M | n/a |
| MODMUL_R2 | w32-il64 | 770.84 M | none | n/a | 6.68 M | n/a |

### brainpoolP512r1 (512-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w32-il64 | 526.19 M | none | n/a | 31.65 M | n/a |
| SUBTRACT | w32-il | 522.99 M | none | n/a | 33.52 M | n/a |
| ADDMOD | w32-il | 752.04 M | none | n/a | 13.15 M | n/a |
| SUBTRACTMOD | w32-il64 | 720.03 M | none | n/a | 13.94 M | n/a |
| MULTIPLYOPERANDSCANNING | w32-il64 | 219.16 M | none | n/a | 11.02 M | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 118.58 M | none | n/a | 11.05 M | n/a |
| MONTGOMERYMULTIPLICATION | w32-il64 | 391.08 M | none | n/a | 1.64 M | n/a |
| COMPARE | w32-il | 1.27 G | none | n/a | 54.09 M | n/a |
| REDUCE | w32-o64 | 16.07 M | none | n/a | 21.82 M | n/a |
| MODMUL | w32-o64 | 2.53 M | none | n/a | 3.46 M | n/a |
| MODEXP | w32-o64 | 19.38 k | none | n/a | 11.60 k | n/a |
| EXPONENTIATION | w32-o64 | 14.14 k | none | n/a | 56.17 k | n/a |
| DIVIDE | w32-il | 8.56 M | none | n/a | 15.51 M | n/a |
| ISQRT | w32-o64 | 128.64 k | none | n/a | 4.98 M | n/a |
| MODMUL_R2 | w32-opt | 179.45 M | none | n/a | 3.33 M | n/a |

### p1024 (1024-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w32-il | 209.81 M | none | n/a | 22.34 M | n/a |
| SUBTRACT | w32-il | 213.39 M | none | n/a | 23.01 M | n/a |
| ADDMOD | w32-il | 294.03 M | none | n/a | 8.84 M | n/a |
| SUBTRACTMOD | w32-il | 296.57 M | none | n/a | 11.25 M | n/a |
| MULTIPLYOPERANDSCANNING | w32-il64 | 62.86 M | none | n/a | 3.39 M | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 23.40 M | none | n/a | 3.39 M | n/a |
| MONTGOMERYMULTIPLICATION | w32-il64 | 97.30 M | none | n/a | 542.15 k | n/a |
| COMPARE | w32-il | 539.28 M | none | n/a | 33.01 M | n/a |
| REDUCE | w32-il64 | 2.87 M | none | n/a | 30.82 M | n/a |
| MODMUL | w32-il64 | 367.20 k | none | n/a | 1.25 M | n/a |
| MODEXP | w32-o64 | 1.09 k | none | n/a | 1.76 k | n/a |
| EXPONENTIATION | w32-il | 798.3 | none | n/a | 13.02 k | n/a |
| DIVIDE | w32-o64 | 945.03 k | none | n/a | 14.12 M | n/a |
| ISQRT | w32-il64 | 8.64 k | none | n/a | 2.53 M | n/a |
| MODMUL_R2 | w32-il64 | 49.86 M | none | n/a | 1.21 M | n/a |

### p2048 (2048-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w32-il64 | 57.24 M | none | n/a | 13.81 M | n/a |
| SUBTRACT | w32-il64 | 57.55 M | none | n/a | 13.76 M | n/a |
| ADDMOD | w32-il | 60.99 M | none | n/a | 6.93 M | n/a |
| SUBTRACTMOD | w32-il64 | 62.17 M | none | n/a | 8.03 M | n/a |
| MULTIPLYOPERANDSCANNING | w32-il64 | 5.36 M | none | n/a | 1.06 M | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 1.79 M | none | n/a | 1.06 M | n/a |
| MONTGOMERYMULTIPLICATION | w32 | 13.25 M | none | n/a | 166.44 k | n/a |
| COMPARE | w32-il64 | 141.11 M | none | n/a | 39.00 M | n/a |
| REDUCE | w32-o64 | 383.43 k | none | n/a | 20.53 M | n/a |
| MODMUL | w32-o64 | 43.05 k | none | n/a | 407.91 k | n/a |
| MODEXP | w32-opt | 36.3 | none | n/a | 249.9 | n/a |
| EXPONENTIATION | w32-o64 | 44.8 | none | n/a | 2.14 k | n/a |
| DIVIDE | w32-il64 | 24.50 k | none | n/a | 11.33 M | n/a |
| ISQRT | w32-o64 | 1.92 k | none | n/a | 1.14 M | n/a |
| MODMUL_R2 | w32-il64 | 6.48 M | none | n/a | 391.49 k | n/a |

## 6. Raw data

Also written to `NVIDIA_GeForce_GTX_1660_SUPER_Report.csv` for analysis.

```csv
kind,device,device_type,kernel,modulus,bits,operation,items,seconds,ops_per_sec,mismatches
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,secp256k1,256,ADD,50000,0.001344129,37198811.153,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,secp256k1,256,ADD,50000,0.000641851,77899700.986,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,secp256k1,256,ADD,50000,0.000592272,84420667.280,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,secp256k1,256,ADD,50000,0.000063488,787550403.226,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,secp256k1,256,ADD,50000,0.000153741,325222322.461,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,secp256k1,256,ADD,50000,0.001428541,35000744.808,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,secp256k1,256,ADD,50000,0.000071788,696495489.529,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,secp256k1,256,ADD,50000,0.001361065,36735936.282,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,secp256k1,256,ADD,50000,0.000112271,445351233.513,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,secp256k1,256,ADD,50000,0.001373316,36408228.407,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,secp256k1,256,ADD,50000,0.000046835,1067576569.312,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,secp256k1,256,ADD,50000,0.001273895,39249701.705,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,secp256k1,256,ADD,50000,0.000113119,442012308.104,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,secp256k1,256,ADD,50000,0.001346392,37136287.894,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,secp256k1,256,ADD,50000,0.000029053,1720988322.034,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,secp256k1,256,ADD,50000,0.001280997,39032099.225,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,secp256k1,256,ADD,50000,0.000028969,1725988601.557,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,secp256k1,256,ADD,50000,0.001308219,38219899.273,0
opencl-kernel,cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,CPU,w8,secp256k1,256,ADD,50000,0.000939620,53212996.658,0
opencl-e2e,cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,CPU,w8,secp256k1,256,ADD,50000,0.001529155,32697795.733,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,secp256k1,256,SUBTRACT,50000,0.001248737,40040457.447,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,secp256k1,256,SUBTRACT,50000,0.000580193,86178218.636,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,secp256k1,256,SUBTRACT,50000,0.000604737,82680558.698,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,secp256k1,256,SUBTRACT,50000,0.000063488,787550403.226,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,secp256k1,256,SUBTRACT,50000,0.000152831,327158830.233,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,secp256k1,256,SUBTRACT,50000,0.001440892,34700727.145,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,secp256k1,256,SUBTRACT,50000,0.000072459,690045483.260,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,secp256k1,256,SUBTRACT,50000,0.001341716,37265709.447,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,secp256k1,256,SUBTRACT,50000,0.000106895,467748755.038,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,secp256k1,256,SUBTRACT,50000,0.001356432,36861413.385,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,secp256k1,256,SUBTRACT,50000,0.000046740,1069748985.910,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,secp256k1,256,SUBTRACT,50000,0.001265733,39502800.774,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,secp256k1,256,SUBTRACT,50000,0.000109226,457766561.364,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,secp256k1,256,SUBTRACT,50000,0.001347740,37099141.980,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,secp256k1,256,SUBTRACT,50000,0.000029161,1714619405.887,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,secp256k1,256,SUBTRACT,50000,0.001255058,39838793.077,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,secp256k1,256,SUBTRACT,50000,0.000028680,1743370391.297,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,secp256k1,256,SUBTRACT,50000,0.001292171,38694570.358,0
opencl-kernel,cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,CPU,w8,secp256k1,256,SUBTRACT,50000,0.001214462,41170494.171,0
opencl-e2e,cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,CPU,w8,secp256k1,256,SUBTRACT,50000,0.001911446,26158207.113,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,secp256k1,256,ADDMOD,50000,0.004107464,12172961.258,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,secp256k1,256,ADDMOD,50000,0.001973779,25332115.983,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,secp256k1,256,ADDMOD,50000,0.005368908,9312880.868,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,secp256k1,256,ADDMOD,50000,0.000061728,810005184.033,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,secp256k1,256,ADDMOD,50000,0.000233463,214166794.869,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,secp256k1,256,ADDMOD,50000,0.001512533,33057128.818,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,secp256k1,256,ADDMOD,50000,0.000111867,446959434.714,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,secp256k1,256,ADDMOD,50000,0.001380547,36217529.586,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,secp256k1,256,ADDMOD,50000,0.000097954,510443952.471,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,secp256k1,256,ADDMOD,50000,0.001345519,37160382.629,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,secp256k1,256,ADDMOD,50000,0.000033165,1507614404.409,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,secp256k1,256,ADDMOD,50000,0.001257465,39762538.497,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,secp256k1,256,ADDMOD,50000,0.000034492,1449611621.282,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,secp256k1,256,ADDMOD,50000,0.001275434,39202340.804,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,secp256k1,256,ADDMOD,50000,0.000026759,1868522570.793,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,secp256k1,256,ADDMOD,50000,0.001283644,38951609.634,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,secp256k1,256,ADDMOD,50000,0.000026803,1865462958.008,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,secp256k1,256,ADDMOD,50000,0.001292232,38692740.241,0
opencl-kernel,cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,CPU,w8,secp256k1,256,ADDMOD,50000,0.001728133,28932959.135,0
opencl-e2e,cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,CPU,w8,secp256k1,256,ADDMOD,50000,0.002303502,21706081.130,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,secp256k1,256,SUBTRACTMOD,50000,0.003108285,16086040.922,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,secp256k1,256,SUBTRACTMOD,50000,0.001600521,31239829.129,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,secp256k1,256,SUBTRACTMOD,50000,0.005269984,9487694.810,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,secp256k1,256,SUBTRACTMOD,50000,0.000061504,812955254.943,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.000232592,214968707.725,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.001513328,33039765.380,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,secp256k1,256,SUBTRACTMOD,50000,0.000102650,487091942.111,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,secp256k1,256,SUBTRACTMOD,50000,0.001378475,36271967.476,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,secp256k1,256,SUBTRACTMOD,50000,0.000102288,488815463.758,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,secp256k1,256,SUBTRACTMOD,50000,0.001367314,36568047.661,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,secp256k1,256,SUBTRACTMOD,50000,0.000031536,1585490636.859,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,secp256k1,256,SUBTRACTMOD,50000,0.001253327,39893818.436,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,secp256k1,256,SUBTRACTMOD,50000,0.000034007,1470285054.670,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,secp256k1,256,SUBTRACTMOD,50000,0.001246506,40112121.434,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,secp256k1,256,SUBTRACTMOD,50000,0.000027020,1850481385.610,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,secp256k1,256,SUBTRACTMOD,50000,0.001299052,38489605.272,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,secp256k1,256,SUBTRACTMOD,50000,0.000026815,1864628784.531,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,secp256k1,256,SUBTRACTMOD,50000,0.001264906,39528631.934,0
opencl-kernel,cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,CPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.001490307,33550136.010,0
opencl-e2e,cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,CPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.002596737,19254934.158,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001963408,25465925.005,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000920508,54317836.710,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001583023,31585138.120,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.006578204,7600858.849,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.008108207,6166591.436,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001576652,31712767.532,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.003162803,15808761.978,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001070883,46690441.062,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.002630225,19009780.290,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000080675,619770660.887,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001594433,31359108.229,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000079201,631305540.391,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001643904,30415401.032,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000035798,1396728247.621,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001589755,31451385.937,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000035550,1406470565.736,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001597037,31307977.027,0
opencl-kernel,cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,CPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.033605496,1487851.866,0
opencl-e2e,cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,CPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.034579877,1445927.640,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001964927,25446237.059,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000912912,54769794.184,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001575055,31744924.199,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000061536,812532501.300,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001832460,27285726.140,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.003377122,14805506.036,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000444007,112610839.826,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.002030988,24618560.140,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000363606,137511515.549,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001918225,26065764.370,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000163140,306485211.958,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001684417,29683860.648,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000191102,261640329.590,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001699836,29414602.465,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000071153,702710617.801,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001627284,30726044.402,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000070032,713959688.149,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001633510,30608934.256,0
opencl-kernel,cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,CPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.024894038,2008513.044,0
opencl-e2e,cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,CPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.025840102,1934976.882,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.013076159,3823752.815,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.006192304,8074538.875,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.001507837,33160083.479,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000060448,827157226.046,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000461845,108261439.024,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.001711246,29218474.326,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000089750,557102796.431,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.001381013,36205308.303,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000034392,1453826621.308,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.001279249,39085432.998,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000043798,1141604503.724,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.001273537,39260734.360,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000046192,1082439223.158,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.001280892,39035299.049,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000045000,1111111273.920,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.001308530,38210813.698,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000044739,1117590486.748,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.001281739,39009499.008,0
opencl-kernel,cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,CPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.115273366,433751.540,0
opencl-e2e,cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,CPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.115719114,432080.737,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,secp256k1,256,COMPARE,50000,0.000763260,65508481.560,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,secp256k1,256,COMPARE,50000,0.000372067,134384406.659,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,secp256k1,256,COMPARE,50000,0.000525680,95114894.929,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,secp256k1,256,COMPARE,50000,0.000063104,792342799.189,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,secp256k1,256,COMPARE,50000,0.000144251,346617951.154,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,secp256k1,256,COMPARE,50000,0.001389608,35981369.769,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,secp256k1,256,COMPARE,50000,0.000067356,742325179.966,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,secp256k1,256,COMPARE,50000,0.001351758,36988868.686,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,secp256k1,256,COMPARE,50000,0.000033866,1476405619.665,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,secp256k1,256,COMPARE,50000,0.001259160,39709012.507,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,secp256k1,256,COMPARE,50000,0.000040686,1228924231.308,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,secp256k1,256,COMPARE,50000,0.001283820,38946269.119,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,secp256k1,256,COMPARE,50000,0.000017931,2788469022.113,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,secp256k1,256,COMPARE,50000,0.001295168,38605028.327,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,secp256k1,256,COMPARE,50000,0.000017037,2934803340.030,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,secp256k1,256,COMPARE,50000,0.001269026,39400293.573,0
opencl-kernel,cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,CPU,w8,secp256k1,256,COMPARE,50000,0.000496441,100716919.891,0
opencl-e2e,cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,CPU,w8,secp256k1,256,COMPARE,50000,0.001273292,39268290.372,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,secp256k1,256,REDUCE,6250,0.000164155,38073763.989,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,secp256k1,256,REDUCE,6250,0.000083893,74499699.848,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,secp256k1,256,REDUCE,6250,0.000540191,11569982.619,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,secp256k1,256,REDUCE,50000,0.000060960,820209973.753,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,secp256k1,256,REDUCE,15400,0.000331635,46436593.771,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,secp256k1,256,REDUCE,15400,0.000801464,19214837.791,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,secp256k1,256,REDUCE,15400,0.000166897,92272499.618,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,secp256k1,256,REDUCE,15400,0.000661523,23279615.094,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,secp256k1,256,REDUCE,15400,0.000062941,244673616.303,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,secp256k1,256,REDUCE,15400,0.000508527,30283549.066,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,secp256k1,256,REDUCE,15400,0.000059120,260487110.464,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,secp256k1,256,REDUCE,15400,0.000505694,30453201.659,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,secp256k1,256,REDUCE,15400,0.000060497,254558210.995,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,secp256k1,256,REDUCE,15400,0.000571769,26933955.806,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,secp256k1,256,REDUCE,15400,0.000056135,274338635.569,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,secp256k1,256,REDUCE,15400,0.000546931,28157119.806,0
opencl-kernel,cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,CPU,w8,secp256k1,256,REDUCE,15400,0.034816988,442312.815,0
opencl-e2e,cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,CPU,w8,secp256k1,256,REDUCE,15400,0.035210266,437372.441,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,secp256k1,256,MODMUL,3125,0.000447959,6976083.976,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,secp256k1,256,MODMUL,3125,0.000225760,13842130.623,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,secp256k1,256,MODMUL,3125,0.000686290,4553468.813,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,secp256k1,256,MODMUL,50000,0.000119296,419125536.481,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,secp256k1,256,MODMUL,15400,0.000930853,16543966.385,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,secp256k1,256,MODMUL,15400,0.001400476,10996261.488,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,secp256k1,256,MODMUL,15400,0.000410284,37534983.485,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,secp256k1,256,MODMUL,15400,0.000887241,17357178.397,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,secp256k1,256,MODMUL,15400,0.000179390,85846444.442,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,secp256k1,256,MODMUL,15400,0.000624622,24654911.863,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,secp256k1,256,MODMUL,15400,0.000147817,104182904.439,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,secp256k1,256,MODMUL,15400,0.000597602,25769659.868,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,secp256k1,256,MODMUL,15400,0.000179129,85971584.381,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,secp256k1,256,MODMUL,15400,0.000686419,22435277.603,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,secp256k1,256,MODMUL,15400,0.000146994,104766168.106,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,secp256k1,256,MODMUL,15400,0.000654039,23545997.459,0
opencl-kernel,cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,CPU,w8,secp256k1,256,MODMUL,15400,0.125650131,122562.546,0
opencl-e2e,cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,CPU,w8,secp256k1,256,MODMUL,15400,0.126031388,122191.783,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,secp256k1,256,MODEXP,781,0.011451136,68202.840,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,secp256k1,256,MODEXP,781,0.005380383,145156.953,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,secp256k1,256,MODEXP,781,0.008833202,88416.409,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,secp256k1,256,MODEXP,50000,0.206239924,242436.086,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,secp256k1,256,MODEXP,15400,0.035670747,431726.310,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,secp256k1,256,MODEXP,15400,0.036167387,425797.972,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,secp256k1,256,MODEXP,15400,0.006027928,2554775.053,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,secp256k1,256,MODEXP,15400,0.006528835,2358766.900,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,secp256k1,256,MODEXP,15400,0.005216447,2952201.026,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,secp256k1,256,MODEXP,15400,0.005664129,2718864.610,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,secp256k1,256,MODEXP,15400,0.002245972,6856719.345,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,secp256k1,256,MODEXP,15400,0.002697598,5708782.252,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,secp256k1,256,MODEXP,15400,0.004234265,3636994.890,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,secp256k1,256,MODEXP,15400,0.004752886,3240136.614,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,secp256k1,256,MODEXP,15400,0.002283430,6744240.144,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,secp256k1,256,MODEXP,15400,0.002781688,5536206.787,0
opencl-kernel,cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,CPU,w8,secp256k1,256,MODEXP,15400,8.256748403,1865.141,0
opencl-e2e,cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,CPU,w8,secp256k1,256,MODEXP,15400,8.283311897,1859.160,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,secp256k1,256,EXPONENTIATION,781,0.003665852,213047.335,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,secp256k1,256,EXPONENTIATION,781,0.001476106,529094.794,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,secp256k1,256,EXPONENTIATION,781,0.028151537,27742.713,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,secp256k1,256,EXPONENTIATION,15400,0.047110202,326893.100,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,secp256k1,256,EXPONENTIATION,15400,0.047633654,323300.833,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,secp256k1,256,EXPONENTIATION,15400,0.012389386,1242999.456,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,secp256k1,256,EXPONENTIATION,15400,0.012889208,1194798.001,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,secp256k1,256,EXPONENTIATION,15400,0.000953663,16148261.609,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,secp256k1,256,EXPONENTIATION,15400,0.001397838,11017013.407,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,secp256k1,256,EXPONENTIATION,15400,0.000664111,23188894.838,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,secp256k1,256,EXPONENTIATION,15400,0.001121866,13727128.996,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,secp256k1,256,EXPONENTIATION,15400,0.000779390,19759042.687,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,secp256k1,256,EXPONENTIATION,15400,0.001286305,11972276.981,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,secp256k1,256,EXPONENTIATION,15400,0.000712363,21618190.350,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,secp256k1,256,EXPONENTIATION,15400,0.001202641,12805151.983,0
opencl-kernel,cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,CPU,w8,secp256k1,256,EXPONENTIATION,15400,4.068294179,3785.370,0
opencl-e2e,cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,CPU,w8,secp256k1,256,EXPONENTIATION,15400,4.066492999,3787.047,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,secp256k1,256,DIVIDE,6250,0.000367262,17017821.404,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,secp256k1,256,DIVIDE,6250,0.000179758,34768955.061,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,secp256k1,256,DIVIDE,6250,0.000536105,11658163.687,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,secp256k1,256,DIVIDE,50000,0.000069760,716743119.266,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,secp256k1,256,DIVIDE,15400,0.000550984,27949992.091,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,secp256k1,256,DIVIDE,15400,0.001165556,13212578.115,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,secp256k1,256,DIVIDE,15400,0.000309899,49693612.469,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,secp256k1,256,DIVIDE,15400,0.000969679,15881544.658,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,secp256k1,256,DIVIDE,15400,0.000148029,104033704.073,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,secp256k1,256,DIVIDE,15400,0.000743830,20703655.559,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,secp256k1,256,DIVIDE,15400,0.000117554,131003516.320,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,secp256k1,256,DIVIDE,15400,0.000724237,21263759.151,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,secp256k1,256,DIVIDE,15400,0.000120160,128162393.021,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,secp256k1,256,DIVIDE,15400,0.000775261,19864278.806,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,secp256k1,256,DIVIDE,15400,0.000116641,132029120.282,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,secp256k1,256,DIVIDE,15400,0.000775451,19859408.966,0
opencl-kernel,cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,CPU,w8,secp256k1,256,DIVIDE,15400,0.055075086,279618.265,0
opencl-e2e,cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,CPU,w8,secp256k1,256,DIVIDE,15400,0.055508499,277434.992,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,secp256k1,256,ISQRT,1562,0.000180499,8653785.345,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,secp256k1,256,ISQRT,1562,0.000081323,19207337.713,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,secp256k1,256,ISQRT,15400,0.006931349,2221789.710,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,secp256k1,256,ISQRT,15400,0.007400467,2080949.739,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,secp256k1,256,ISQRT,15400,0.004957234,3106571.111,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,secp256k1,256,ISQRT,15400,0.005457946,2821574.279,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,secp256k1,256,ISQRT,15400,0.001452735,10600694.237,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,secp256k1,256,ISQRT,15400,0.001900524,8103028.354,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,secp256k1,256,ISQRT,15400,0.000996234,15458215.456,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,secp256k1,256,ISQRT,15400,0.001451880,10606936.562,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,secp256k1,256,ISQRT,15400,0.001150431,13386287.458,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,secp256k1,256,ISQRT,15400,0.001654954,9305394.992,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,secp256k1,256,ISQRT,15400,0.001006983,15293205.950,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,secp256k1,256,ISQRT,15400,0.001514849,10166030.395,0
opencl-kernel,cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,CPU,w8,secp256k1,256,ISQRT,15400,0.767333961,20069.488,0
opencl-e2e,cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,CPU,w8,secp256k1,256,ISQRT,15400,0.770826110,19978.566,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,secp256k1,256,MODMUL_R2,50000,0.007495067,6671054.472,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,secp256k1,256,MODMUL_R2,50000,0.003685548,13566503.224,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,secp256k1,256,MODMUL_R2,50000,0.010859282,4604355.983,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,secp256k1,256,MODMUL_R2,50000,0.000065824,759601361.206,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,secp256k1,256,MODMUL_R2,50000,0.000447664,111690895.913,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,secp256k1,256,MODMUL_R2,50000,0.001698423,29439074.744,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,secp256k1,256,MODMUL_R2,50000,0.000096417,518580575.116,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,secp256k1,256,MODMUL_R2,50000,0.001364561,36641820.246,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,secp256k1,256,MODMUL_R2,50000,0.000089914,556086480.328,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,secp256k1,256,MODMUL_R2,50000,0.001319413,37895637.848,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,secp256k1,256,MODMUL_R2,50000,0.000063968,781640729.343,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,secp256k1,256,MODMUL_R2,50000,0.001305148,38309831.368,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,secp256k1,256,MODMUL_R2,50000,0.000076540,653252797.221,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,secp256k1,256,MODMUL_R2,50000,0.001351223,37003515.197,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,secp256k1,256,MODMUL_R2,50000,0.000061727,810019669.955,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,secp256k1,256,MODMUL_R2,50000,0.001332174,37532633.068,0
opencl-kernel,cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,CPU,w8,secp256k1,256,MODMUL_R2,50000,0.111394453,448855.384,0
opencl-e2e,cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,CPU,w8,secp256k1,256,MODMUL_R2,50000,0.112224248,445536.512,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,rsa256(composite),256,ADD,50000,0.001364719,36637578.704,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,rsa256(composite),256,ADD,50000,0.000640948,78009454.128,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,rsa256(composite),256,ADD,50000,0.000591640,84510849.564,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,rsa256(composite),256,ADD,50000,0.000065248,766307013.242,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,rsa256(composite),256,ADD,50000,0.000124069,403001780.540,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,rsa256(composite),256,ADD,50000,0.001365429,36618530.448,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,rsa256(composite),256,ADD,50000,0.000071486,699437725.304,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,rsa256(composite),256,ADD,50000,0.001336104,37422236.488,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,rsa256(composite),256,ADD,50000,0.000108756,459744713.538,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,rsa256(composite),256,ADD,50000,0.001384708,36108696.293,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,rsa256(composite),256,ADD,50000,0.000054669,914595370.132,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,rsa256(composite),256,ADD,50000,0.001278809,39098879.096,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,rsa256(composite),256,ADD,50000,0.000110563,452230825.425,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,rsa256(composite),256,ADD,50000,0.001339531,37326499.296,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,rsa256(composite),256,ADD,50000,0.000029172,1713976214.059,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,rsa256(composite),256,ADD,50000,0.001300869,38435847.719,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,rsa256(composite),256,ADD,50000,0.000028980,1725329922.550,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,rsa256(composite),256,ADD,50000,0.001266821,39468877.161,0
opencl-kernel,cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,CPU,w8,rsa256(composite),256,ADD,50000,0.000966069,51756140.753,0
opencl-e2e,cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,CPU,w8,rsa256(composite),256,ADD,50000,0.001592902,31389250.543,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACT,50000,0.001251216,39961128.049,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACT,50000,0.000578800,86385610.522,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACT,50000,0.000604397,82727076.884,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,rsa256(composite),256,SUBTRACT,50000,0.000062624,798415942.770,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,rsa256(composite),256,SUBTRACT,50000,0.000125020,399935869.725,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,rsa256(composite),256,SUBTRACT,50000,0.001374738,36370565.073,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,rsa256(composite),256,SUBTRACT,50000,0.000072394,690665777.290,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,rsa256(composite),256,SUBTRACT,50000,0.001339609,37324322.730,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,rsa256(composite),256,SUBTRACT,50000,0.000111687,447679687.800,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,rsa256(composite),256,SUBTRACT,50000,0.001356927,36847966.586,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,rsa256(composite),256,SUBTRACT,50000,0.000053805,929280520.661,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,rsa256(composite),256,SUBTRACT,50000,0.001278639,39104079.975,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,rsa256(composite),256,SUBTRACT,50000,0.000109035,458568603.346,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,rsa256(composite),256,SUBTRACT,50000,0.001330991,37565992.853,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,rsa256(composite),256,SUBTRACT,50000,0.000028767,1738100214.077,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,rsa256(composite),256,SUBTRACT,50000,0.001308707,38205647.185,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,rsa256(composite),256,SUBTRACT,50000,0.000028458,1756970581.666,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,rsa256(composite),256,SUBTRACT,50000,0.001276612,39166166.448,0
opencl-kernel,cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,CPU,w8,rsa256(composite),256,SUBTRACT,50000,0.001343166,37225480.971,0
opencl-e2e,cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,CPU,w8,rsa256(composite),256,SUBTRACT,50000,0.001951831,25616972.916,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,rsa256(composite),256,ADDMOD,50000,0.003397767,14715547.150,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,rsa256(composite),256,ADDMOD,50000,0.001693177,29530284.730,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,rsa256(composite),256,ADDMOD,50000,0.004863232,10281228.562,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,rsa256(composite),256,ADDMOD,50000,0.000061792,809166235.111,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,rsa256(composite),256,ADDMOD,50000,0.000172357,290095627.324,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,rsa256(composite),256,ADDMOD,50000,0.001417862,35264362.360,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,rsa256(composite),256,ADDMOD,50000,0.000093102,537045451.772,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,rsa256(composite),256,ADDMOD,50000,0.001386275,36067878.467,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,rsa256(composite),256,ADDMOD,50000,0.000096810,516475923.381,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,rsa256(composite),256,ADDMOD,50000,0.001375563,36348753.808,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,rsa256(composite),256,ADDMOD,50000,0.000038948,1283762093.011,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,rsa256(composite),256,ADDMOD,50000,0.001274479,39231714.740,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,rsa256(composite),256,ADDMOD,50000,0.000034010,1470154203.407,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,rsa256(composite),256,ADDMOD,50000,0.001265218,39518884.475,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,rsa256(composite),256,ADDMOD,50000,0.000028527,1752725946.663,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,rsa256(composite),256,ADDMOD,50000,0.001346122,37143735.642,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,rsa256(composite),256,ADDMOD,50000,0.000026853,1861985423.079,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,rsa256(composite),256,ADDMOD,50000,0.001295074,38607832.281,0
opencl-kernel,cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,CPU,w8,rsa256(composite),256,ADDMOD,50000,0.001300640,38442614.683,0
opencl-e2e,cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,CPU,w8,rsa256(composite),256,ADDMOD,50000,0.002078922,24050926.396,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACTMOD,50000,0.003098184,16138486.235,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACTMOD,50000,0.001599623,31257364.884,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACTMOD,50000,0.005275861,9477126.084,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,rsa256(composite),256,SUBTRACTMOD,50000,0.000063200,791139240.506,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,rsa256(composite),256,SUBTRACTMOD,50000,0.000188489,265267498.894,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,rsa256(composite),256,SUBTRACTMOD,50000,0.001442634,34658827.412,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,rsa256(composite),256,SUBTRACTMOD,50000,0.000102054,489936803.489,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,rsa256(composite),256,SUBTRACTMOD,50000,0.001395301,35834559.847,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,rsa256(composite),256,SUBTRACTMOD,50000,0.000106499,469487691.130,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,rsa256(composite),256,SUBTRACTMOD,50000,0.001355871,36876662.004,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,rsa256(composite),256,SUBTRACTMOD,50000,0.000038436,1300862693.882,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,rsa256(composite),256,SUBTRACTMOD,50000,0.001270270,39361711.275,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,rsa256(composite),256,SUBTRACTMOD,50000,0.000033361,1498755027.934,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,rsa256(composite),256,SUBTRACTMOD,50000,0.001267996,39432299.973,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,rsa256(composite),256,SUBTRACTMOD,50000,0.000027392,1825354255.722,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,rsa256(composite),256,SUBTRACTMOD,50000,0.001287400,38837967.273,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,rsa256(composite),256,SUBTRACTMOD,50000,0.000027019,1850545170.019,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,rsa256(composite),256,SUBTRACTMOD,50000,0.001265010,39525379.822,0
opencl-kernel,cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,CPU,w8,rsa256(composite),256,SUBTRACTMOD,50000,0.001533488,32605405.812,0
opencl-e2e,cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,CPU,w8,rsa256(composite),256,SUBTRACTMOD,50000,0.002514290,19886329.413,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001977647,25282571.211,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000909886,54951940.624,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001583721,31571216.946,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.005294275,9444163.735,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.006802553,7350181.657,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001568757,31872367.790,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.003145239,15897043.287,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001084347,46110699.761,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.002630024,19011233.471,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000096954,515708182.112,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001613287,30992625.186,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000079854,626142379.009,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001587046,31505071.409,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000036386,1374156477.781,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001602335,31204460.637,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000034601,1445046529.843,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001591476,31417377.783,0
opencl-kernel,cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,CPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.033591199,1488485.120,0
opencl-e2e,cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,CPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.034511019,1448812.629,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001968151,25404553.867,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000928018,53878264.753,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001600488,31240472.191,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000061632,811266874.351,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001529724,32685635.052,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.003050478,16390873.507,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000441220,113322150.441,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.002025408,24686384.519,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000357696,139783514.402,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001915733,26099670.728,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000200258,249677933.128,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001719706,29074737.307,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000159189,314092076.556,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001683627,29697791.014,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000081796,611277087.253,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001664239,30043759.968,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000078464,637234558.406,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001612355,31010540.458,0
opencl-kernel,cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,CPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.024909470,2007268.722,0
opencl-e2e,cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,CPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.025804954,1937612.445,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.013074673,3824187.415,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.006193971,8072365.817,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.001513420,33037755.065,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000062368,801693175.988,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000388795,128602462.239,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.001623300,30801453.438,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000088084,567640579.143,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.001392387,35909557.525,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000035050,1426534506.455,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.001300960,38433154.669,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000045934,1088518462.123,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.001259962,39683737.111,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000045393,1101491139.533,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.001279623,39074008.063,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000048078,1039976390.540,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.001318848,37911874.796,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000045572,1097163522.292,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.001318877,37911041.536,0
opencl-kernel,cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,CPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.114922972,435074.025,0
opencl-e2e,cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,CPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.115621629,432445.040,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,rsa256(composite),256,COMPARE,50000,0.000762040,65613361.212,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,rsa256(composite),256,COMPARE,50000,0.000373128,134002277.463,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,rsa256(composite),256,COMPARE,50000,0.000529931,94351898.487,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,rsa256(composite),256,COMPARE,50000,0.000063488,787550403.226,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,rsa256(composite),256,COMPARE,50000,0.000122042,409694943.878,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,rsa256(composite),256,COMPARE,50000,0.001366016,36602792.640,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,rsa256(composite),256,COMPARE,50000,0.000069562,718783279.919,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,rsa256(composite),256,COMPARE,50000,0.001354534,36913065.296,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,rsa256(composite),256,COMPARE,50000,0.000040423,1236919311.928,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,rsa256(composite),256,COMPARE,50000,0.001275680,39194781.533,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,rsa256(composite),256,COMPARE,50000,0.000047715,1047887811.959,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,rsa256(composite),256,COMPARE,50000,0.001264896,39528944.805,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,rsa256(composite),256,COMPARE,50000,0.000023592,2119359744.195,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,rsa256(composite),256,COMPARE,50000,0.001398542,35751519.272,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,rsa256(composite),256,COMPARE,50000,0.000017754,2816261193.658,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,rsa256(composite),256,COMPARE,50000,0.001255873,39812942.596,0
opencl-kernel,cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,CPU,w8,rsa256(composite),256,COMPARE,50000,0.000502382,99525847.743,0
opencl-e2e,cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,CPU,w8,rsa256(composite),256,COMPARE,50000,0.001245285,40151453.896,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,rsa256(composite),256,REDUCE,6250,0.000272850,22906365.191,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,rsa256(composite),256,REDUCE,6250,0.000134203,46571238.772,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,rsa256(composite),256,REDUCE,6250,0.000541191,11548603.713,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,rsa256(composite),256,REDUCE,50000,0.000062752,796787353.391,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,rsa256(composite),256,REDUCE,15400,0.000294373,52314576.667,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,rsa256(composite),256,REDUCE,15400,0.000787262,19561467.717,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,rsa256(composite),256,REDUCE,15400,0.000167056,92184663.914,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,rsa256(composite),256,REDUCE,15400,0.000663998,23192842.525,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,rsa256(composite),256,REDUCE,15400,0.000075465,204067926.565,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,rsa256(composite),256,REDUCE,15400,0.000569031,27063552.290,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,rsa256(composite),256,REDUCE,15400,0.000070857,217339146.213,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,rsa256(composite),256,REDUCE,15400,0.000526670,29240326.768,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,rsa256(composite),256,REDUCE,15400,0.000076035,202538204.473,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,rsa256(composite),256,REDUCE,15400,0.000555390,27728267.515,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,rsa256(composite),256,REDUCE,15400,0.000066604,231217345.745,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,rsa256(composite),256,REDUCE,15400,0.000572564,26896558.236,0
opencl-kernel,cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,CPU,w8,rsa256(composite),256,REDUCE,15400,0.034604085,445034.163,0
opencl-e2e,cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,CPU,w8,rsa256(composite),256,REDUCE,15400,0.035085295,438930.326,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,rsa256(composite),256,MODMUL,3125,0.000449514,6951952.298,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,rsa256(composite),256,MODMUL,3125,0.000225868,13835509.888,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,rsa256(composite),256,MODMUL,3125,0.000687102,4548087.658,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,rsa256(composite),256,MODMUL,50000,0.000096544,517898574.743,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,rsa256(composite),256,MODMUL,15400,0.000821970,18735478.388,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,rsa256(composite),256,MODMUL,15400,0.001346046,11440916.678,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,rsa256(composite),256,MODMUL,15400,0.000408137,37732429.530,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,rsa256(composite),256,MODMUL,15400,0.000914311,16843285.651,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,rsa256(composite),256,MODMUL,15400,0.000218660,70428980.164,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,rsa256(composite),256,MODMUL,15400,0.000719706,21397624.403,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,rsa256(composite),256,MODMUL,15400,0.000179094,85988349.445,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,rsa256(composite),256,MODMUL,15400,0.000626297,24588969.605,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,rsa256(composite),256,MODMUL,15400,0.000221212,69616462.249,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,rsa256(composite),256,MODMUL,15400,0.000760979,20237088.820,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,rsa256(composite),256,MODMUL,15400,0.000179017,86025367.497,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,rsa256(composite),256,MODMUL,15400,0.000681497,22597310.174,0
opencl-kernel,cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,CPU,w8,rsa256(composite),256,MODMUL,15400,0.000000000,inf,0
opencl-e2e,cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,CPU,w8,rsa256(composite),256,MODMUL,15400,0.000000000,inf,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,rsa256(composite),256,MODEXP,781,0.011059823,70615.958,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,rsa256(composite),256,MODEXP,781,0.005226825,149421.494,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,rsa256(composite),256,MODEXP,781,0.008818922,88559.577,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,rsa256(composite),256,MODEXP,50000,0.199450627,250688.608,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,rsa256(composite),256,MODEXP,15400,0.035747119,430803.948,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,rsa256(composite),256,MODEXP,15400,0.036302879,424208.780,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,rsa256(composite),256,MODEXP,15400,0.005990098,2570909.525,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,rsa256(composite),256,MODEXP,15400,0.006487881,2373656.379,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,rsa256(composite),256,MODEXP,15400,0.005168229,2979744.062,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,rsa256(composite),256,MODEXP,15400,0.005680823,2710874.815,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,rsa256(composite),256,MODEXP,15400,0.002785834,5527967.531,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,rsa256(composite),256,MODEXP,15400,0.003233907,4762041.808,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,rsa256(composite),256,MODEXP,15400,0.005159431,2984825.253,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,rsa256(composite),256,MODEXP,15400,0.005672623,2714793.521,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,rsa256(composite),256,MODEXP,15400,0.002836207,5429786.951,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,rsa256(composite),256,MODEXP,15400,0.003341375,4608881.026,0
opencl-kernel,cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,CPU,w8,rsa256(composite),256,MODEXP,15400,0.000000000,inf,0
opencl-e2e,cpu-sandybridge-Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,CPU,w8,rsa256(composite),256,MODEXP,15400,0.000000000,inf,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,rsa256(composite),256,EXPONENTIATION,781,0.003651717,213871.996,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,rsa256(composite),256,EXPONENTIATION,781,0.001532941,509478.183,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,rsa256(composite),256,EXPONENTIATION,781,0.028970935,26958.053,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,rsa256(composite),256,EXPONENTIATION,15400,0.047374201,325071.445,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,rsa256(composite),256,EXPONENTIATION,15400,0.048124617,320002.547,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,rsa256(composite),256,EXPONENTIATION,15400,0.012402680,1241667.123,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,rsa256(composite),256,EXPONENTIATION,15400,0.012779403,1205064.116,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,rsa256(composite),256,EXPONENTIATION,15400,0.000954170,16139681.391,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,rsa256(composite),256,EXPONENTIATION,15400,0.001450637,10616025.880,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,rsa256(composite),256,EXPONENTIATION,15400,0.000816169,18868639.362,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,rsa256(composite),256,EXPONENTIATION,15400,0.001270209,12123988.187,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,rsa256(composite),256,EXPONENTIATION,15400,0.000957357,16085953.329,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,rsa256(composite),256,EXPONENTIATION,15400,0.001472478,10458560.917,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,rsa256(composite),256,EXPONENTIATION,15400,0.000851173,18092680.143,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,rsa256(composite),256,EXPONENTIATION,15400,0.001376844,11184999.807,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,rsa256(composite),256,DIVIDE,6250,0.000368023,16982631.543,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,rsa256(composite),256,DIVIDE,6250,0.000178861,34943342.506,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,rsa256(composite),256,DIVIDE,6250,0.000535914,11662319.478,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,rsa256(composite),256,DIVIDE,50000,0.000057824,864692861.096,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,rsa256(composite),256,DIVIDE,15400,0.000549725,28014005.779,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,rsa256(composite),256,DIVIDE,15400,0.001190714,12933415.634,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,rsa256(composite),256,DIVIDE,15400,0.000312452,49287556.672,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,rsa256(composite),256,DIVIDE,15400,0.000962522,15999634.339,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,rsa256(composite),256,DIVIDE,15400,0.000151205,101848488.550,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,rsa256(composite),256,DIVIDE,15400,0.000797900,19300663.711,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,rsa256(composite),256,DIVIDE,15400,0.000145838,105596610.880,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,rsa256(composite),256,DIVIDE,15400,0.000749749,20540207.043,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,rsa256(composite),256,DIVIDE,15400,0.000147702,104263951.490,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,rsa256(composite),256,DIVIDE,15400,0.000806234,19101155.154,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,rsa256(composite),256,DIVIDE,15400,0.000144861,106308764.932,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,rsa256(composite),256,DIVIDE,15400,0.000795141,19367634.682,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,rsa256(composite),256,ISQRT,1562,0.000178307,8760174.735,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,rsa256(composite),256,ISQRT,1562,0.000080000,19524993.390,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,rsa256(composite),256,ISQRT,15400,0.006932671,2221366.067,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,rsa256(composite),256,ISQRT,15400,0.007400573,2080919.951,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,rsa256(composite),256,ISQRT,15400,0.004957219,3106580.522,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,rsa256(composite),256,ISQRT,15400,0.005450360,2825501.489,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,rsa256(composite),256,ISQRT,15400,0.001467203,10496162.413,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,rsa256(composite),256,ISQRT,15400,0.001959583,7858815.232,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,rsa256(composite),256,ISQRT,15400,0.001237601,12443429.610,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,rsa256(composite),256,ISQRT,15400,0.001695386,9083477.060,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,rsa256(composite),256,ISQRT,15400,0.001427759,10786133.567,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,rsa256(composite),256,ISQRT,15400,0.001941991,7930005.697,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,rsa256(composite),256,ISQRT,15400,0.001249474,12325186.363,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,rsa256(composite),256,ISQRT,15400,0.001750607,8796948.497,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,rsa256(composite),256,MODMUL_R2,50000,0.007480983,6683613.601,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,rsa256(composite),256,MODMUL_R2,50000,0.003691824,13543440.723,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,rsa256(composite),256,MODMUL_R2,50000,0.010893771,4589778.843,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,rsa256(composite),256,MODMUL_R2,50000,0.000066144,755926463.474,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,rsa256(composite),256,MODMUL_R2,50000,0.000447233,111798525.185,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,rsa256(composite),256,MODMUL_R2,50000,0.001673062,29885324.343,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,rsa256(composite),256,MODMUL_R2,50000,0.000097180,514509061.897,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,rsa256(composite),256,MODMUL_R2,50000,0.001386030,36074255.376,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,rsa256(composite),256,MODMUL_R2,50000,0.000090824,550515437.109,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,rsa256(composite),256,MODMUL_R2,50000,0.001358020,36818308.929,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,rsa256(composite),256,MODMUL_R2,50000,0.000068361,731410585.184,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,rsa256(composite),256,MODMUL_R2,50000,0.001300798,38437942.574,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,rsa256(composite),256,MODMUL_R2,50000,0.000078948,633328412.553,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,rsa256(composite),256,MODMUL_R2,50000,0.001359650,36774171.612,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,rsa256(composite),256,MODMUL_R2,50000,0.000064864,770842943.548,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,rsa256(composite),256,MODMUL_R2,50000,0.001337366,37386924.766,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,ADD,25000,0.000789830,31652381.073,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,ADD,25000,0.000362921,68885529.438,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,ADD,25000,0.000382849,65299893.695,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,brainpoolP512r1,512,ADD,50000,0.000075168,665176670.924,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,brainpoolP512r1,512,ADD,25000,0.000178711,139890734.440,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,brainpoolP512r1,512,ADD,25000,0.001432385,17453408.780,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,brainpoolP512r1,512,ADD,25000,0.000096975,257798682.366,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,brainpoolP512r1,512,ADD,25000,0.001364189,18325905.330,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,brainpoolP512r1,512,ADD,25000,0.000264844,94395176.073,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,brainpoolP512r1,512,ADD,25000,0.001527302,16368734.237,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,brainpoolP512r1,512,ADD,25000,0.000214901,116332670.525,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,brainpoolP512r1,512,ADD,25000,0.001472095,16982599.310,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,brainpoolP512r1,512,ADD,25000,0.000220150,113558883.069,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,brainpoolP512r1,512,ADD,25000,0.001463920,17077435.991,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,brainpoolP512r1,512,ADD,25000,0.000048782,512483121.458,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,brainpoolP512r1,512,ADD,25000,0.001281487,19508587.414,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,brainpoolP512r1,512,ADD,25000,0.000047511,526193137.752,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,brainpoolP512r1,512,ADD,25000,0.001316969,18982982.173,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACT,25000,0.000745801,33521005.273,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACT,25000,0.000343449,72790997.507,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACT,25000,0.000391044,63931431.878,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,brainpoolP512r1,512,SUBTRACT,50000,0.000076032,657617845.118,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,brainpoolP512r1,512,SUBTRACT,25000,0.000179363,139382186.464,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,brainpoolP512r1,512,SUBTRACT,25000,0.001441687,17340795.932,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,brainpoolP512r1,512,SUBTRACT,25000,0.000100980,247573661.420,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,brainpoolP512r1,512,SUBTRACT,25000,0.001355644,18441420.183,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,brainpoolP512r1,512,SUBTRACT,25000,0.000265010,94336086.255,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,brainpoolP512r1,512,SUBTRACT,25000,0.001529876,16341193.384,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,brainpoolP512r1,512,SUBTRACT,25000,0.000212727,117521519.216,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,brainpoolP512r1,512,SUBTRACT,25000,0.001469577,17011696.991,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,brainpoolP512r1,512,SUBTRACT,25000,0.000217812,114777932.966,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,brainpoolP512r1,512,SUBTRACT,25000,0.001452954,17206325.630,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,brainpoolP512r1,512,SUBTRACT,25000,0.000047802,522990737.818,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,brainpoolP512r1,512,SUBTRACT,25000,0.001287102,19423480.023,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,brainpoolP512r1,512,SUBTRACT,25000,0.000048065,520129251.711,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,brainpoolP512r1,512,SUBTRACT,25000,0.001319234,18950390.527,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,ADDMOD,25000,0.001901269,13149112.797,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,ADDMOD,25000,0.000961101,26011833.482,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,ADDMOD,25000,0.002958204,8451073.649,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,brainpoolP512r1,512,ADDMOD,50000,0.000077280,646997929.607,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,brainpoolP512r1,512,ADDMOD,25000,0.000240806,103818000.607,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,brainpoolP512r1,512,ADDMOD,25000,0.001515129,16500245.204,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,brainpoolP512r1,512,ADDMOD,25000,0.000123058,203156268.140,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,brainpoolP512r1,512,ADDMOD,25000,0.001395453,17915329.282,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,brainpoolP512r1,512,ADDMOD,25000,0.000253567,98593268.703,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,brainpoolP512r1,512,ADDMOD,25000,0.001483038,16857289.442,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,brainpoolP512r1,512,ADDMOD,25000,0.000083646,298878885.699,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,brainpoolP512r1,512,ADDMOD,25000,0.001303867,19173735.167,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,brainpoolP512r1,512,ADDMOD,25000,0.000085285,293134611.988,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,brainpoolP512r1,512,ADDMOD,25000,0.001331182,18780302.527,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,brainpoolP512r1,512,ADDMOD,25000,0.000033243,752035904.691,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,brainpoolP512r1,512,ADDMOD,25000,0.001288311,19405252.696,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,brainpoolP512r1,512,ADDMOD,25000,0.000034477,725121523.527,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,brainpoolP512r1,512,ADDMOD,25000,0.001285507,19447579.910,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACTMOD,25000,0.001792966,13943377.108,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000879805,28415389.019,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACTMOD,25000,0.003152220,7930918.416,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000074304,672911283.376,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000272630,91699367.600,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,brainpoolP512r1,512,SUBTRACTMOD,25000,0.001527960,16361685.433,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000134512,185857069.243,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,brainpoolP512r1,512,SUBTRACTMOD,25000,0.001374165,18192865.394,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000260253,96060382.652,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,brainpoolP512r1,512,SUBTRACTMOD,25000,0.001528482,16356096.354,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000084979,294190350.209,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,brainpoolP512r1,512,SUBTRACTMOD,25000,0.001340202,18653905.254,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000086177,290100917.793,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,brainpoolP512r1,512,SUBTRACTMOD,25000,0.001333325,18750118.072,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000035586,702522113.831,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,brainpoolP512r1,512,SUBTRACTMOD,25000,0.001277114,19575385.891,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000034721,720025632.102,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,brainpoolP512r1,512,SUBTRACTMOD,25000,0.001253357,19946431.215,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.002268733,11019366.359,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.001089545,22945356.172,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.001157719,21594185.383,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.019683347,1270109.189,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.021279425,1174843.775,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.005048626,4951842.356,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.006604319,3785401.660,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.005638094,4434122.658,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.007212917,3466004.128,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.000225082,111070669.827,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.001762316,14185878.288,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.000224589,111314435.117,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.001758103,14219872.004,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.000115010,217372265.416,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.001667585,14991739.632,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.000114073,219157820.141,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.001660498,15055724.468,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.002262730,11048600.413,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.001086822,23002846.569,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.001163177,21492858.242,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000073376,681421718.273,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.005664698,4413297.965,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.007219391,3462895.997,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.001491870,16757491.908,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.003051125,8193699.102,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.001119561,22330180.652,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.002697964,9266246.812,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.001016105,24603755.060,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.002550353,9802564.790,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.001041767,23997686.907,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.002557328,9775828.632,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000222505,112356990.201,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.001775868,14077622.853,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000210831,118578354.689,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.001755006,14244965.307,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.015283882,1635710.086,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.006762248,3696995.438,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.001971678,12679555.223,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000071680,697544642.857,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000782308,31956720.791,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.002044102,12230309.084,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000224356,111430011.976,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.001471390,16990736.461,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000072395,345327890.806,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.001334157,18738424.936,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000105509,236946456.644,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.001357397,18417604.396,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000076064,328670484.936,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.001315485,19004397.797,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000093090,268557313.881,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.001368482,18268416.921,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000063926,391077296.037,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.001306857,19129866.800,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,COMPARE,25000,0.000462219,54086915.199,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,COMPARE,25000,0.000252604,98969130.060,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,COMPARE,25000,0.000354936,70435233.957,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,brainpoolP512r1,512,COMPARE,50000,0.000075296,664045898.853,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,brainpoolP512r1,512,COMPARE,25000,0.000154655,161650043.471,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,brainpoolP512r1,512,COMPARE,25000,0.001414802,17670316.596,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,brainpoolP512r1,512,COMPARE,25000,0.000084136,297137858.903,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,brainpoolP512r1,512,COMPARE,25000,0.001362164,18353149.313,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,brainpoolP512r1,512,COMPARE,25000,0.000059373,421066968.883,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,brainpoolP512r1,512,COMPARE,25000,0.001289062,19393947.381,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,brainpoolP512r1,512,COMPARE,25000,0.000057583,434156092.788,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,brainpoolP512r1,512,COMPARE,25000,0.001303502,19179103.540,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,brainpoolP512r1,512,COMPARE,25000,0.000019648,1272394399.644,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,brainpoolP512r1,512,COMPARE,25000,0.001273592,19629520.244,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,brainpoolP512r1,512,COMPARE,25000,0.000019908,1255772297.364,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,brainpoolP512r1,512,COMPARE,25000,0.001254463,19928846.266,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,REDUCE,3125,0.000143221,21819426.316,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,REDUCE,3125,0.000071835,43502467.515,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,REDUCE,3125,0.000338212,9239766.641,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,brainpoolP512r1,512,REDUCE,50000,0.000071680,697544642.857,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,brainpoolP512r1,512,REDUCE,15400,0.000861893,17867648.339,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,brainpoolP512r1,512,REDUCE,15400,0.001749911,8800447.584,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,brainpoolP512r1,512,REDUCE,15400,0.000675421,22800590.903,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,brainpoolP512r1,512,REDUCE,15400,0.001589987,9685614.017,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,brainpoolP512r1,512,REDUCE,15400,0.000200667,76744032.064,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,brainpoolP512r1,512,REDUCE,15400,0.001114918,13812674.574,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,brainpoolP512r1,512,REDUCE,15400,0.000194433,79204650.024,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,brainpoolP512r1,512,REDUCE,15400,0.001089698,14132355.250,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,brainpoolP512r1,512,REDUCE,15400,0.000226220,68075325.024,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,brainpoolP512r1,512,REDUCE,15400,0.001150203,13388941.637,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,brainpoolP512r1,512,REDUCE,15400,0.000222701,69151014.992,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,brainpoolP512r1,512,REDUCE,15400,0.001144118,13460149.796,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL,1562,0.000451796,3457312.672,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL,1562,0.000214670,7276285.551,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL,1562,0.000583024,2679135.006,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,brainpoolP512r1,512,MODMUL,50000,0.000292864,170727709.790,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,brainpoolP512r1,512,MODMUL,15400,0.002447460,6292237.762,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,brainpoolP512r1,512,MODMUL,15400,0.003351692,4594694.270,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,brainpoolP512r1,512,MODMUL,15400,0.001737513,8863243.051,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,brainpoolP512r1,512,MODMUL,15400,0.002656337,5797457.196,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,brainpoolP512r1,512,MODMUL,15400,0.000767591,20062766.230,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,brainpoolP512r1,512,MODMUL,15400,0.001685384,9137383.732,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,brainpoolP512r1,512,MODMUL,15400,0.000617469,24940524.847,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,brainpoolP512r1,512,MODMUL,15400,0.001513173,10177289.784,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,brainpoolP512r1,512,MODMUL,15400,0.000992717,15512979.401,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,brainpoolP512r1,512,MODMUL,15400,0.001906079,8079413.257,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,brainpoolP512r1,512,MODMUL,15400,0.000794482,19383700.242,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,brainpoolP512r1,512,MODMUL,15400,0.001713248,8988774.697,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,MODEXP,390,0.033632021,11596.092,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,MODEXP,390,0.016593210,23503.590,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,MODEXP,390,0.018072790,21579.402,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,brainpoolP512r1,512,MODEXP,50000,1.466943502,34084.476,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,brainpoolP512r1,512,MODEXP,15400,0.513077256,30014.973,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,brainpoolP512r1,512,MODEXP,15400,0.513658057,29981.035,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,brainpoolP512r1,512,MODEXP,15400,0.046701486,329753.961,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,brainpoolP512r1,512,MODEXP,15400,0.047655997,323149.257,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,brainpoolP512r1,512,MODEXP,15400,0.032631237,471940.430,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,brainpoolP512r1,512,MODEXP,15400,0.033451727,460364.871,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,brainpoolP512r1,512,MODEXP,15400,0.020127451,765124.210,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,brainpoolP512r1,512,MODEXP,15400,0.020330820,757470.678,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,brainpoolP512r1,512,MODEXP,15400,0.038327137,401804.079,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,brainpoolP512r1,512,MODEXP,15400,0.039228494,392571.787,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,brainpoolP512r1,512,MODEXP,15400,0.020786420,740868.318,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,brainpoolP512r1,512,MODEXP,15400,0.021628547,712021.940,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,EXPONENTIATION,390,0.006943389,56168.537,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,EXPONENTIATION,390,0.003515415,110939.960,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,EXPONENTIATION,390,0.043401529,8985.859,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,brainpoolP512r1,512,EXPONENTIATION,15400,0.389989505,39488.242,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,brainpoolP512r1,512,EXPONENTIATION,15400,0.390791042,39407.249,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,brainpoolP512r1,512,EXPONENTIATION,15400,0.099963632,154056.027,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,brainpoolP512r1,512,EXPONENTIATION,15400,0.101530651,151678.334,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,brainpoolP512r1,512,EXPONENTIATION,15400,0.028781254,535070.502,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,brainpoolP512r1,512,EXPONENTIATION,15400,0.029477211,522437.487,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,brainpoolP512r1,512,EXPONENTIATION,15400,0.027582814,558318.669,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,brainpoolP512r1,512,EXPONENTIATION,15400,0.028403087,542194.585,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,brainpoolP512r1,512,EXPONENTIATION,15400,0.030404953,506496.424,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,brainpoolP512r1,512,EXPONENTIATION,15400,0.031399630,490451.639,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,brainpoolP512r1,512,EXPONENTIATION,15400,0.028674621,537060.280,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,brainpoolP512r1,512,EXPONENTIATION,15400,0.029485601,522288.828,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,DIVIDE,3125,0.000201425,15514464.634,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,DIVIDE,3125,0.000101856,30680581.894,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,DIVIDE,3125,0.000335259,9321151.323,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,brainpoolP512r1,512,DIVIDE,50000,0.000143392,348694487.838,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,brainpoolP512r1,512,DIVIDE,15400,0.001669098,9226540.740,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,brainpoolP512r1,512,DIVIDE,15400,0.002783127,5533344.303,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,brainpoolP512r1,512,DIVIDE,15400,0.001308539,11768849.997,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,brainpoolP512r1,512,DIVIDE,15400,0.002409613,6391067.796,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,brainpoolP512r1,512,DIVIDE,15400,0.000437847,35172104.351,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,brainpoolP512r1,512,DIVIDE,15400,0.001543029,9980369.917,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,brainpoolP512r1,512,DIVIDE,15400,0.000421810,36509331.971,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,brainpoolP512r1,512,DIVIDE,15400,0.001503709,10241343.203,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,brainpoolP512r1,512,DIVIDE,15400,0.000364964,42195941.112,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,brainpoolP512r1,512,DIVIDE,15400,0.001486492,10359961.616,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,brainpoolP512r1,512,DIVIDE,15400,0.000386601,39834354.816,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,brainpoolP512r1,512,DIVIDE,15400,0.001507045,10218672.948,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,ISQRT,781,0.000156712,4983663.671,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,ISQRT,781,0.000071526,10919100.325,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,brainpoolP512r1,512,ISQRT,15400,0.038654464,398401.593,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,brainpoolP512r1,512,ISQRT,15400,0.039635694,388538.674,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,brainpoolP512r1,512,ISQRT,15400,0.035789407,430294.920,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,brainpoolP512r1,512,ISQRT,15400,0.036674582,419909.352,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,brainpoolP512r1,512,ISQRT,15400,0.006702363,2297697.087,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,brainpoolP512r1,512,ISQRT,15400,0.007606546,2024571.982,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,brainpoolP512r1,512,ISQRT,15400,0.006071403,2536481.303,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,brainpoolP512r1,512,ISQRT,15400,0.006967186,2210361.545,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,brainpoolP512r1,512,ISQRT,15400,0.008581971,1794459.585,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,brainpoolP512r1,512,ISQRT,15400,0.009494035,1622071.115,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,brainpoolP512r1,512,ISQRT,15400,0.007918806,1944737.617,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,brainpoolP512r1,512,ISQRT,15400,0.008845990,1740901.798,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL_R2,25000,0.007503491,3331782.495,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL_R2,25000,0.003485160,7173271.921,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL_R2,25000,0.009139511,2735376.094,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,brainpoolP512r1,512,MODMUL_R2,50000,0.000098208,509123492.994,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,brainpoolP512r1,512,MODMUL_R2,25000,0.001140946,21911642.890,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,brainpoolP512r1,512,MODMUL_R2,25000,0.002408428,10380214.807,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,brainpoolP512r1,512,MODMUL_R2,25000,0.000238109,104993969.120,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,brainpoolP512r1,512,MODMUL_R2,25000,0.001456316,17166602.833,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,brainpoolP512r1,512,MODMUL_R2,25000,0.000139311,179454658.252,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,brainpoolP512r1,512,MODMUL_R2,25000,0.001396437,17902705.976,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,brainpoolP512r1,512,MODMUL_R2,25000,0.000142746,175136227.908,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,brainpoolP512r1,512,MODMUL_R2,25000,0.001381557,18095525.960,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,brainpoolP512r1,512,MODMUL_R2,25000,0.000140251,178251836.112,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,brainpoolP512r1,512,MODMUL_R2,25000,0.001410755,17721007.802,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,brainpoolP512r1,512,MODMUL_R2,25000,0.000146550,170590230.646,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,brainpoolP512r1,512,MODMUL_R2,25000,0.001418344,17626189.689,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,p1024,1024,ADD,12500,0.000559480,22342173.043,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,p1024,1024,ADD,12500,0.000240080,52065990.551,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,p1024,1024,ADD,12500,0.000276720,45172028.342,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,p1024,1024,ADD,50000,0.000113792,439398200.225,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p1024,1024,ADD,15400,0.000744181,20693890.701,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p1024,1024,ADD,15400,0.002211815,6962607.698,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p1024,1024,ADD,15400,0.000363628,42350984.049,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p1024,1024,ADD,15400,0.001833860,8397587.309,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,p1024,1024,ADD,15400,0.000371772,41423240.015,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,p1024,1024,ADD,15400,0.001832520,8403728.167,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p1024,1024,ADD,15400,0.000193597,79546714.322,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p1024,1024,ADD,15400,0.001680715,9162766.984,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p1024,1024,ADD,15400,0.000192688,79921913.363,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p1024,1024,ADD,15400,0.001625697,9472859.878,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p1024,1024,ADD,15400,0.000059579,258480192.108,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p1024,1024,ADD,15400,0.001505099,10231885.051,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p1024,1024,ADD,15400,0.000060263,255546633.794,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p1024,1024,ADD,15400,0.001535255,10030906.942,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,p1024,1024,SUBTRACT,12500,0.000543231,23010467.651,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,p1024,1024,SUBTRACT,12500,0.000226416,55208130.401,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,p1024,1024,SUBTRACT,12500,0.000280352,44586811.505,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,p1024,1024,SUBTRACT,50000,0.000112640,443892045.455,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p1024,1024,SUBTRACT,15400,0.000744970,20671975.219,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p1024,1024,SUBTRACT,15400,0.002212590,6960168.991,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p1024,1024,SUBTRACT,15400,0.000375982,40959404.940,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p1024,1024,SUBTRACT,15400,0.001820988,8456947.446,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,p1024,1024,SUBTRACT,15400,0.000360532,42714659.124,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,p1024,1024,SUBTRACT,15400,0.001806163,8526362.256,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p1024,1024,SUBTRACT,15400,0.000189537,81250648.124,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p1024,1024,SUBTRACT,15400,0.001679195,9171061.272,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p1024,1024,SUBTRACT,15400,0.000193919,79414626.026,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p1024,1024,SUBTRACT,15400,0.001625520,9473891.756,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p1024,1024,SUBTRACT,15400,0.000058577,262902134.294,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p1024,1024,SUBTRACT,15400,0.001502743,10247926.687,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p1024,1024,SUBTRACT,15400,0.000061252,251420205.828,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p1024,1024,SUBTRACT,15400,0.001518797,10139603.824,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,p1024,1024,ADDMOD,12500,0.001413472,8843471.981,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,p1024,1024,ADDMOD,12500,0.000661236,18903991.884,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,p1024,1024,ADDMOD,12500,0.002128721,5872070.650,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,p1024,1024,ADDMOD,50000,0.000112736,443514050.525,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p1024,1024,ADDMOD,15400,0.000903920,17036905.211,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p1024,1024,ADDMOD,15400,0.002372409,6491292.228,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p1024,1024,ADDMOD,15400,0.000415668,37048797.855,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p1024,1024,ADDMOD,15400,0.001898938,8109795.716,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,p1024,1024,ADDMOD,15400,0.000383770,40128204.202,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,p1024,1024,ADDMOD,15400,0.001846433,8340405.369,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p1024,1024,ADDMOD,15400,0.000135459,113687531.555,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p1024,1024,ADDMOD,15400,0.001606740,9584624.425,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p1024,1024,ADDMOD,15400,0.000136652,112695092.153,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p1024,1024,ADDMOD,15400,0.001572974,9790371.451,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p1024,1024,ADDMOD,15400,0.000042513,362242027.900,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p1024,1024,ADDMOD,15400,0.001490197,10334203.941,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p1024,1024,ADDMOD,15400,0.000044107,349151154.248,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p1024,1024,ADDMOD,15400,0.001516019,10158184.072,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,p1024,1024,SUBTRACTMOD,12500,0.001111532,11245740.322,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,p1024,1024,SUBTRACTMOD,12500,0.000519725,24051179.598,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,p1024,1024,SUBTRACTMOD,12500,0.002049047,6100396.814,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,p1024,1024,SUBTRACTMOD,50000,0.000112640,443892045.455,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p1024,1024,SUBTRACTMOD,15400,0.000863406,17836336.715,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p1024,1024,SUBTRACTMOD,15400,0.002341166,6577918.804,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p1024,1024,SUBTRACTMOD,15400,0.000432570,35601174.438,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p1024,1024,SUBTRACTMOD,15400,0.001879560,8193407.202,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,p1024,1024,SUBTRACTMOD,15400,0.000380907,40429806.510,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,p1024,1024,SUBTRACTMOD,15400,0.001853841,8307077.002,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p1024,1024,SUBTRACTMOD,15400,0.000135125,113968637.164,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p1024,1024,SUBTRACTMOD,15400,0.001579160,9752020.141,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p1024,1024,SUBTRACTMOD,15400,0.000136046,113197030.288,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p1024,1024,SUBTRACTMOD,15400,0.001557678,9886510.663,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p1024,1024,SUBTRACTMOD,15400,0.000042148,365378699.338,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p1024,1024,SUBTRACTMOD,15400,0.001513672,10173935.014,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p1024,1024,SUBTRACTMOD,15400,0.000044205,348376934.243,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p1024,1024,SUBTRACTMOD,15400,0.001515648,10160670.695,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.003687868,3389492.238,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.001665919,7503366.139,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.002137652,5847537.347,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,15400,0.243061193,63358.530,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,15400,0.244863583,62892.161,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p1024,1024,MULTIPLYOPERANDSCANNING,15400,0.067366494,228600.289,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p1024,1024,MULTIPLYOPERANDSCANNING,15400,0.069062177,222987.468,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,p1024,1024,MULTIPLYOPERANDSCANNING,15400,0.014131657,1089751.898,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,p1024,1024,MULTIPLYOPERANDSCANNING,15400,0.015720009,979643.206,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p1024,1024,MULTIPLYOPERANDSCANNING,15400,0.000304723,50537694.477,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p1024,1024,MULTIPLYOPERANDSCANNING,15400,0.002159166,7132383.412,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p1024,1024,MULTIPLYOPERANDSCANNING,15400,0.000301610,51059298.697,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p1024,1024,MULTIPLYOPERANDSCANNING,15400,0.002073259,7427918.707,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p1024,1024,MULTIPLYOPERANDSCANNING,15400,0.000212848,72352116.781,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p1024,1024,MULTIPLYOPERANDSCANNING,15400,0.002018664,7628808.282,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p1024,1024,MULTIPLYOPERANDSCANNING,15400,0.000198866,77439121.857,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p1024,1024,MULTIPLYOPERANDSCANNING,15400,0.002006702,7674283.443,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.003688746,3388685.570,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.001669550,7487047.042,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.002135204,5854241.786,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000124928,400230532.787,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,15400,0.026946441,571504.043,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,15400,0.028717824,536252.331,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p1024,1024,MULTIPLYPRODUCTSCANNING,15400,0.006818079,2258700.714,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p1024,1024,MULTIPLYPRODUCTSCANNING,15400,0.008628727,1784736.047,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,p1024,1024,MULTIPLYPRODUCTSCANNING,15400,0.002231423,6901425.544,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,p1024,1024,MULTIPLYPRODUCTSCANNING,15400,0.004038016,3813754.117,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p1024,1024,MULTIPLYPRODUCTSCANNING,15400,0.001772484,8688372.156,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p1024,1024,MULTIPLYPRODUCTSCANNING,15400,0.003609983,4265948.046,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p1024,1024,MULTIPLYPRODUCTSCANNING,15400,0.001760336,8748329.922,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p1024,1024,MULTIPLYPRODUCTSCANNING,15400,0.003536766,4354260.305,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p1024,1024,MULTIPLYPRODUCTSCANNING,15400,0.000569019,27064122.593,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p1024,1024,MULTIPLYPRODUCTSCANNING,15400,0.002365173,6511151.734,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p1024,1024,MULTIPLYPRODUCTSCANNING,15400,0.000534258,28825026.342,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p1024,1024,MULTIPLYPRODUCTSCANNING,15400,0.002335432,6594069.147,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.023056166,542154.320,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.009996348,1250456.659,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.003312316,3773794.501,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000174816,286015010.068,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,15400,0.003089623,4984426.993,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,15400,0.004540580,3391637.210,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p1024,1024,MONTGOMERYMULTIPLICATION,15400,0.000630875,24410540.851,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p1024,1024,MONTGOMERYMULTIPLICATION,15400,0.002099857,7333832.441,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,p1024,1024,MONTGOMERYMULTIPLICATION,15400,0.000175178,87910615.209,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,p1024,1024,MONTGOMERYMULTIPLICATION,15400,0.001605496,9592051.740,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p1024,1024,MONTGOMERYMULTIPLICATION,15400,0.000259128,59430077.150,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p1024,1024,MONTGOMERYMULTIPLICATION,15400,0.001749782,8801096.324,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p1024,1024,MONTGOMERYMULTIPLICATION,15400,0.000208723,73781994.770,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p1024,1024,MONTGOMERYMULTIPLICATION,15400,0.001636186,9412132.663,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p1024,1024,MONTGOMERYMULTIPLICATION,15400,0.000189480,81275058.946,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p1024,1024,MONTGOMERYMULTIPLICATION,15400,0.001648307,9342919.945,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p1024,1024,MONTGOMERYMULTIPLICATION,15400,0.000128471,119871354.942,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p1024,1024,MONTGOMERYMULTIPLICATION,15400,0.001598721,9632700.014,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,p1024,1024,COMPARE,12500,0.000378726,33005397.526,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,p1024,1024,COMPARE,12500,0.000178834,69897239.632,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,p1024,1024,COMPARE,12500,0.000253310,49346657.867,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,p1024,1024,COMPARE,50000,0.000112896,442885487.528,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p1024,1024,COMPARE,15400,0.000254777,60445037.766,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p1024,1024,COMPARE,15400,0.001722427,8940872.548,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p1024,1024,COMPARE,15400,0.000136945,112453865.916,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p1024,1024,COMPARE,15400,0.001583473,9725457.543,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p1024,1024,COMPARE,15400,0.000075336,204417325.548,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p1024,1024,COMPARE,15400,0.001562581,9855489.521,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p1024,1024,COMPARE,15400,0.000075750,203300187.981,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p1024,1024,COMPARE,15400,0.001504491,10236020.248,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p1024,1024,COMPARE,15400,0.000023179,664394808.377,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p1024,1024,COMPARE,15400,0.001473452,10451647.132,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p1024,1024,COMPARE,15400,0.000023794,647221683.734,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p1024,1024,COMPARE,15400,0.001492825,10316011.675,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,p1024,1024,REDUCE,1562,0.000050681,30820262.120,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,p1024,1024,REDUCE,1562,0.000026914,58036834.939,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,p1024,1024,REDUCE,1562,0.000225267,6933996.734,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,p1024,1024,REDUCE,50000,0.000106816,468094667.466,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p1024,1024,REDUCE,15400,0.004634016,3323251.394,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p1024,1024,REDUCE,15400,0.006097175,2525759.854,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p1024,1024,REDUCE,15400,0.002521482,6107519.436,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p1024,1024,REDUCE,15400,0.003981728,3867667.515,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p1024,1024,REDUCE,15400,0.000603623,25512613.439,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p1024,1024,REDUCE,15400,0.002089800,7369126.360,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p1024,1024,REDUCE,15400,0.000608088,25325282.576,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p1024,1024,REDUCE,15400,0.002036743,7561091.474,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p1024,1024,REDUCE,15400,0.000556603,27667837.168,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p1024,1024,REDUCE,15400,0.001993306,7725858.521,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p1024,1024,REDUCE,15400,0.000544899,28262122.946,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p1024,1024,REDUCE,15400,0.001998537,7705636.722,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,p1024,1024,MODMUL,781,0.000623426,1252755.070,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,p1024,1024,MODMUL,781,0.000280127,2788021.851,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,p1024,1024,MODMUL,781,0.000704777,1108151.818,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,p1024,1024,MODMUL,50000,0.000886784,56383516.166,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p1024,1024,MODMUL,15400,0.018328490,840221.973,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p1024,1024,MODMUL,15400,0.019781291,778513.394,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p1024,1024,MODMUL,15400,0.006746530,2282654.952,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p1024,1024,MODMUL,15400,0.008212678,1875149.639,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p1024,1024,MODMUL,15400,0.002905831,5299688.634,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p1024,1024,MODMUL,15400,0.004396719,3502611.800,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p1024,1024,MODMUL,15400,0.002248669,6848495.933,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p1024,1024,MODMUL,15400,0.003782664,4071204.887,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p1024,1024,MODMUL,15400,0.002836176,5429846.459,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p1024,1024,MODMUL,15400,0.004304894,3577323.941,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p1024,1024,MODMUL,15400,0.002126895,7240601.884,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p1024,1024,MODMUL,15400,0.003593644,4285343.777,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,p1024,1024,MODEXP,195,0.111094517,1755.262,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,p1024,1024,MODEXP,195,0.055187032,3533.439,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,p1024,1024,MODEXP,195,0.050555654,3857.135,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,p1024,1024,MODEXP,50000,3.797446251,13166.743,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p1024,1024,MODEXP,15400,4.295223035,3585.378,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p1024,1024,MODEXP,15400,4.293325385,3586.963,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p1024,1024,MODEXP,15400,0.630977161,24406.589,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p1024,1024,MODEXP,15400,0.632277407,24356.398,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p1024,1024,MODEXP,15400,0.293891055,52400.370,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p1024,1024,MODEXP,15400,0.295335629,52144.064,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p1024,1024,MODEXP,15400,0.179487110,85800.033,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p1024,1024,MODEXP,15400,0.181894623,84664.405,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p1024,1024,MODEXP,15400,0.296352724,51965.104,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p1024,1024,MODEXP,15400,0.299107236,51486.551,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p1024,1024,MODEXP,15400,0.180588786,85276.613,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p1024,1024,MODEXP,15400,0.182093616,84571.883,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,p1024,1024,EXPONENTIATION,195,0.014981719,13015.863,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,p1024,1024,EXPONENTIATION,195,0.007780987,25061.088,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,p1024,1024,EXPONENTIATION,195,0.112339067,1735.816,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p1024,1024,EXPONENTIATION,15400,3.677023927,4188.170,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p1024,1024,EXPONENTIATION,15400,3.682283631,4182.187,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p1024,1024,EXPONENTIATION,15400,0.801796342,19206.872,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p1024,1024,EXPONENTIATION,15400,0.802545299,19188.948,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p1024,1024,EXPONENTIATION,15400,0.249285335,61776.598,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p1024,1024,EXPONENTIATION,15400,0.250725774,61421.687,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p1024,1024,EXPONENTIATION,15400,0.255595205,60251.522,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p1024,1024,EXPONENTIATION,15400,0.257133672,59891.028,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p1024,1024,EXPONENTIATION,15400,0.244267555,63045.622,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p1024,1024,EXPONENTIATION,15400,0.245175309,62812.198,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p1024,1024,EXPONENTIATION,15400,0.250491281,61479.186,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p1024,1024,EXPONENTIATION,15400,0.252728600,60934.932,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,p1024,1024,DIVIDE,1562,0.000110605,14122322.972,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,p1024,1024,DIVIDE,1562,0.000056542,27625486.590,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,p1024,1024,DIVIDE,1562,0.000224640,6953347.242,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,p1024,1024,DIVIDE,50000,0.000207616,240829223.181,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p1024,1024,DIVIDE,15400,0.101860850,151186.643,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p1024,1024,DIVIDE,15400,0.103070311,149412.570,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p1024,1024,DIVIDE,15400,0.012852881,1198174.947,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p1024,1024,DIVIDE,15400,0.014636001,1052199.985,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p1024,1024,DIVIDE,15400,0.001770295,8699115.267,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p1024,1024,DIVIDE,15400,0.003591707,4287654.756,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p1024,1024,DIVIDE,15400,0.001652864,9317161.427,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p1024,1024,DIVIDE,15400,0.003402888,4525567.748,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p1024,1024,DIVIDE,15400,0.001768043,8710195.240,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p1024,1024,DIVIDE,15400,0.003528388,4364599.234,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p1024,1024,DIVIDE,15400,0.001686415,9131797.047,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p1024,1024,DIVIDE,15400,0.003442397,4473626.963,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,p1024,1024,ISQRT,390,0.000153909,2533965.341,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,p1024,1024,ISQRT,390,0.000073310,5319875.392,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p1024,1024,ISQRT,15400,1.322892122,11641.161,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p1024,1024,ISQRT,15400,1.324199225,11629.670,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p1024,1024,ISQRT,15400,0.275790124,55839.563,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p1024,1024,ISQRT,15400,0.276827640,55630.283,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p1024,1024,ISQRT,15400,0.047926203,321327.354,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p1024,1024,ISQRT,15400,0.049389708,311805.852,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p1024,1024,ISQRT,15400,0.045573774,337913.643,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p1024,1024,ISQRT,15400,0.046979896,327799.790,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p1024,1024,ISQRT,15400,0.048727245,316044.956,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p1024,1024,ISQRT,15400,0.050136034,307164.304,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p1024,1024,ISQRT,15400,0.045134000,341206.186,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p1024,1024,ISQRT,15400,0.046541230,330889.407,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,p1024,1024,MODMUL_R2,12500,0.010335906,1209376.317,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,p1024,1024,MODMUL_R2,12500,0.004499293,2778214.236,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,p1024,1024,MODMUL_R2,12500,0.011125525,1123542.485,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,p1024,1024,MODMUL_R2,50000,0.000323808,154412491.353,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p1024,1024,MODMUL_R2,15400,0.004813699,3199202.984,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p1024,1024,MODMUL_R2,15400,0.006256485,2461446.046,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p1024,1024,MODMUL_R2,15400,0.000718296,21439631.630,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p1024,1024,MODMUL_R2,15400,0.002199752,7000789.105,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p1024,1024,MODMUL_R2,15400,0.000417029,36927885.787,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p1024,1024,MODMUL_R2,15400,0.001885442,8167845.780,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p1024,1024,MODMUL_R2,15400,0.000309058,49828834.813,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p1024,1024,MODMUL_R2,15400,0.001725694,8923945.711,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p1024,1024,MODMUL_R2,15400,0.000349508,44061936.433,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p1024,1024,MODMUL_R2,15400,0.001789340,8606525.066,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p1024,1024,MODMUL_R2,15400,0.000250683,61432155.146,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p1024,1024,MODMUL_R2,15400,0.001684994,9139498.576,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,p2048,2048,ADD,6250,0.000452651,13807548.810,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,p2048,2048,ADD,6250,0.000174887,35737359.246,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,p2048,2048,ADD,6250,0.000231046,27050890.677,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,p2048,2048,ADD,50000,0.000212928,234821160.204,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p2048,2048,ADD,15400,0.001852313,8313929.934,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p2048,2048,ADD,15400,0.004320092,3564738.871,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p2048,2048,ADD,15400,0.000821184,18753409.307,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p2048,2048,ADD,15400,0.003296296,4671910.538,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,p2048,2048,ADD,15400,0.000434086,35476844.299,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,p2048,2048,ADD,15400,0.002883681,5340396.565,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p2048,2048,ADD,15400,0.000453546,33954663.640,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p2048,2048,ADD,15400,0.002877755,5351393.539,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p2048,2048,ADD,15400,0.000457151,33686892.992,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p2048,2048,ADD,15400,0.002859023,5386455.348,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p2048,2048,ADD,15400,0.000117259,131333220.865,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p2048,2048,ADD,15400,0.002697768,5708422.584,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p2048,2048,ADD,15400,0.000109180,141051483.466,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p2048,2048,ADD,15400,0.002614882,5889367.183,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,p2048,2048,SUBTRACT,6250,0.000454327,13756612.180,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,p2048,2048,SUBTRACT,6250,0.000172473,36237575.319,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,p2048,2048,SUBTRACT,6250,0.000238040,26256091.802,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,p2048,2048,SUBTRACT,50000,0.000211840,236027190.332,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p2048,2048,SUBTRACT,15400,0.001826722,8430401.291,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p2048,2048,SUBTRACT,15400,0.004292600,3587569.328,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p2048,2048,SUBTRACT,15400,0.000829306,18569744.632,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p2048,2048,SUBTRACT,15400,0.003336678,4615368.903,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,p2048,2048,SUBTRACT,15400,0.000438397,35127982.790,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,p2048,2048,SUBTRACT,15400,0.002867930,5369726.622,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p2048,2048,SUBTRACT,15400,0.000453821,33934081.588,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p2048,2048,SUBTRACT,15400,0.002899267,5311687.593,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p2048,2048,SUBTRACT,15400,0.000450819,34160048.485,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p2048,2048,SUBTRACT,15400,0.002865800,5373717.538,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p2048,2048,SUBTRACT,15400,0.000109380,140793570.567,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p2048,2048,SUBTRACT,15400,0.002546286,6048024.509,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p2048,2048,SUBTRACT,15400,0.000108599,141806140.610,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p2048,2048,SUBTRACT,15400,0.002618835,5880477.481,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,p2048,2048,ADDMOD,6250,0.000902056,6928616.439,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,p2048,2048,ADDMOD,6250,0.000405090,15428672.850,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,p2048,2048,ADDMOD,6250,0.001506658,4148253.845,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,p2048,2048,ADDMOD,50000,0.000213600,234082397.004,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p2048,2048,ADDMOD,15400,0.001888454,8154818.372,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p2048,2048,ADDMOD,15400,0.004431661,3474995.053,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p2048,2048,ADDMOD,15400,0.000913992,16849161.715,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p2048,2048,ADDMOD,15400,0.003363554,4578490.467,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,p2048,2048,ADDMOD,15400,0.000448309,34351306.684,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,p2048,2048,ADDMOD,15400,0.002888738,5331047.562,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p2048,2048,ADDMOD,15400,0.000256888,59948307.511,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p2048,2048,ADDMOD,15400,0.002719913,5661945.748,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p2048,2048,ADDMOD,15400,0.000259101,59436272.063,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p2048,2048,ADDMOD,15400,0.002689491,5725990.724,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p2048,2048,ADDMOD,15400,0.000102471,150286454.526,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p2048,2048,ADDMOD,15400,0.002524520,6100169.727,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p2048,2048,ADDMOD,15400,0.000105206,146379590.751,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p2048,2048,ADDMOD,15400,0.002616939,5884737.814,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,p2048,2048,SUBTRACTMOD,6250,0.000778290,8030425.909,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,p2048,2048,SUBTRACTMOD,6250,0.000346416,18041890.505,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,p2048,2048,SUBTRACTMOD,6250,0.001478744,4226559.584,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,p2048,2048,SUBTRACTMOD,50000,0.000212992,234750600.962,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p2048,2048,SUBTRACTMOD,15400,0.002161072,7126093.025,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p2048,2048,SUBTRACTMOD,15400,0.004593202,3352780.978,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p2048,2048,SUBTRACTMOD,15400,0.000918339,16769406.384,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p2048,2048,SUBTRACTMOD,15400,0.003330811,4623498.541,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,p2048,2048,SUBTRACTMOD,15400,0.000463811,33203181.019,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,p2048,2048,SUBTRACTMOD,15400,0.002895912,5317841.078,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p2048,2048,SUBTRACTMOD,15400,0.000259512,59342152.699,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p2048,2048,SUBTRACTMOD,15400,0.002682623,5740650.129,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p2048,2048,SUBTRACTMOD,15400,0.000260609,59092363.720,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p2048,2048,SUBTRACTMOD,15400,0.002685430,5734649.585,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p2048,2048,SUBTRACTMOD,15400,0.000103501,148790970.206,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p2048,2048,SUBTRACTMOD,15400,0.002513011,6128106.914,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p2048,2048,SUBTRACTMOD,15400,0.000100537,153177482.228,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p2048,2048,SUBTRACTMOD,15400,0.002616623,5885448.645,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.005910978,1057354.625,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.002635518,2371450.343,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.003223524,1938871.820,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,15400,1.021380428,15077.634,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,15400,1.022306149,15063.981,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p2048,2048,MULTIPLYOPERANDSCANNING,15400,0.282521848,54509.059,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p2048,2048,MULTIPLYOPERANDSCANNING,15400,0.285724738,53898.028,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,p2048,2048,MULTIPLYOPERANDSCANNING,15400,0.053558530,287535.897,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,p2048,2048,MULTIPLYOPERANDSCANNING,15400,0.056470149,272710.455,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p2048,2048,MULTIPLYOPERANDSCANNING,15400,0.001461616,10536282.992,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p2048,2048,MULTIPLYOPERANDSCANNING,15400,0.004518725,3408041.002,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p2048,2048,MULTIPLYOPERANDSCANNING,15400,0.001327066,11604547.316,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p2048,2048,MULTIPLYOPERANDSCANNING,15400,0.004393207,3505411.857,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p2048,2048,MULTIPLYOPERANDSCANNING,15400,0.001264981,12174095.952,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p2048,2048,MULTIPLYOPERANDSCANNING,15400,0.004492424,3427993.468,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p2048,2048,MULTIPLYOPERANDSCANNING,15400,0.001166530,13201546.993,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p2048,2048,MULTIPLYOPERANDSCANNING,15400,0.004335890,3551750.668,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.005909441,1057629.640,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.002641617,2365975.108,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.003225273,1937820.465,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.000430784,116067449.116,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,15400,0.108004561,142586.571,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,15400,0.111135342,138569.781,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p2048,2048,MULTIPLYPRODUCTSCANNING,15400,0.027070181,568891.653,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p2048,2048,MULTIPLYPRODUCTSCANNING,15400,0.030152546,510736.307,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,p2048,2048,MULTIPLYPRODUCTSCANNING,15400,0.006865120,2243223.708,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,p2048,2048,MULTIPLYPRODUCTSCANNING,15400,0.009937985,1549609.908,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p2048,2048,MULTIPLYPRODUCTSCANNING,15400,0.006894392,2233699.480,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p2048,2048,MULTIPLYPRODUCTSCANNING,15400,0.009981884,1542794.939,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p2048,2048,MULTIPLYPRODUCTSCANNING,15400,0.006917728,2226164.430,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p2048,2048,MULTIPLYPRODUCTSCANNING,15400,0.009971116,1544461.026,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p2048,2048,MULTIPLYPRODUCTSCANNING,15400,0.003639640,4231187.637,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p2048,2048,MULTIPLYPRODUCTSCANNING,15400,0.006837242,2252370.164,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p2048,2048,MULTIPLYPRODUCTSCANNING,15400,0.003500663,4399166.701,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p2048,2048,MULTIPLYPRODUCTSCANNING,15400,0.006657954,2313022.864,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.037551303,166438.965,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.016542471,377815.382,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.006277034,995693.198,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000767200,65172054.223,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,15400,0.103521916,148760.771,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,15400,0.105851685,145486.584,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p2048,2048,MONTGOMERYMULTIPLICATION,15400,0.002655037,5800295.858,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p2048,2048,MONTGOMERYMULTIPLICATION,15400,0.005101147,3018928.885,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,p2048,2048,MONTGOMERYMULTIPLICATION,15400,0.000471526,32659920.970,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32,p2048,2048,MONTGOMERYMULTIPLICATION,15400,0.002885347,5337312.983,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p2048,2048,MONTGOMERYMULTIPLICATION,15400,0.000816006,18872410.701,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p2048,2048,MONTGOMERYMULTIPLICATION,15400,0.003250928,4737108.892,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p2048,2048,MONTGOMERYMULTIPLICATION,15400,0.000661671,23274409.271,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p2048,2048,MONTGOMERYMULTIPLICATION,15400,0.003132083,4916855.671,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p2048,2048,MONTGOMERYMULTIPLICATION,15400,0.000881793,17464416.698,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p2048,2048,MONTGOMERYMULTIPLICATION,15400,0.003413654,4511295.024,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p2048,2048,MONTGOMERYMULTIPLICATION,15400,0.000660645,23310553.957,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p2048,2048,MONTGOMERYMULTIPLICATION,15400,0.003117173,4940373.892,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,p2048,2048,COMPARE,6250,0.000160239,39004252.413,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,p2048,2048,COMPARE,6250,0.000068431,91332846.561,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,p2048,2048,COMPARE,6250,0.000169491,36875129.609,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,p2048,2048,COMPARE,50000,0.000217760,229610580.456,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p2048,2048,COMPARE,15400,0.000518954,29675076.617,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p2048,2048,COMPARE,15400,0.002966277,5191693.036,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p2048,2048,COMPARE,15400,0.000256781,59973284.477,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p2048,2048,COMPARE,15400,0.002671177,5765248.982,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p2048,2048,COMPARE,15400,0.000136490,112828795.262,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p2048,2048,COMPARE,15400,0.002826774,5447906.513,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p2048,2048,COMPARE,15400,0.000133074,115725080.826,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p2048,2048,COMPARE,15400,0.002556437,6024009.093,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p2048,2048,COMPARE,15400,0.000044965,342488964.848,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p2048,2048,COMPARE,15400,0.002562719,6009242.415,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p2048,2048,COMPARE,15400,0.000044293,347684710.957,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p2048,2048,COMPARE,15400,0.002468407,6238841.687,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,p2048,2048,REDUCE,781,0.000038037,20532660.769,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,p2048,2048,REDUCE,781,0.000020326,38423696.106,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,p2048,2048,REDUCE,781,0.000173816,4493257.451,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,p2048,2048,REDUCE,50000,0.000236000,211864406.780,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p2048,2048,REDUCE,15400,0.805117788,19127.636,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p2048,2048,REDUCE,15400,0.809307843,19028.606,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p2048,2048,REDUCE,15400,0.010040739,1533751.654,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p2048,2048,REDUCE,15400,0.012470766,1234888.060,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p2048,2048,REDUCE,15400,0.002231247,6901969.988,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p2048,2048,REDUCE,15400,0.004719631,3262967.004,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p2048,2048,REDUCE,15400,0.002036904,7560493.825,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p2048,2048,REDUCE,15400,0.004496428,3424940.901,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p2048,2048,REDUCE,15400,0.002351900,6547897.499,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p2048,2048,REDUCE,15400,0.004881920,3154496.553,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p2048,2048,REDUCE,15400,0.002355322,6538384.111,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p2048,2048,REDUCE,15400,0.004837307,3183589.526,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,p2048,2048,MODMUL,390,0.000956102,407906.252,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,p2048,2048,MODMUL,390,0.000438818,888751.126,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,p2048,2048,MODMUL,390,0.001029289,378902.338,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,p2048,2048,MODMUL,50000,0.003799360,13160111.177,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p2048,2048,MODMUL,15400,1.373425348,11212.841,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p2048,2048,MODMUL,15400,1.373488874,11212.322,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p2048,2048,MODMUL,15400,0.030143308,510892.833,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p2048,2048,MODMUL,15400,0.032587492,472573.955,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p2048,2048,MODMUL,15400,0.011993843,1283992.127,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p2048,2048,MODMUL,15400,0.014507403,1061527.000,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p2048,2048,MODMUL,15400,0.009060005,1699778.297,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p2048,2048,MODMUL,15400,0.011508392,1338153.939,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p2048,2048,MODMUL,15400,0.014435883,1066786.145,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p2048,2048,MODMUL,15400,0.014955657,1029710.701,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p2048,2048,MODMUL,15400,0.011091526,1388447.359,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p2048,2048,MODMUL,15400,0.013554512,1136153.034,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,p2048,2048,MODEXP,97,0.388211226,249.864,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,p2048,2048,MODEXP,97,0.189261011,512.520,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,p2048,2048,MODEXP,97,0.169269301,573.051,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,p2048,2048,MODEXP,50000,5.248878002,9525.845,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p2048,2048,MODEXP,15400,189.087005490,81.444,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p2048,2048,MODEXP,15400,189.374759020,81.320,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p2048,2048,MODEXP,15400,25.356944182,607.329,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p2048,2048,MODEXP,15400,25.432302041,605.529,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p2048,2048,MODEXP,15400,2.672018928,5763.432,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p2048,2048,MODEXP,15400,2.675358453,5756.238,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p2048,2048,MODEXP,15400,7.285036055,2113.922,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p2048,2048,MODEXP,15400,7.227530671,2130.742,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p2048,2048,MODEXP,15400,2.683055631,5739.724,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p2048,2048,MODEXP,15400,2.685369383,5734.779,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p2048,2048,MODEXP,15400,6.040691777,2549.377,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p2048,2048,MODEXP,15400,6.095465427,2526.468,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,p2048,2048,EXPONENTIATION,97,0.045352332,2138.810,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,p2048,2048,EXPONENTIATION,97,0.022857602,4243.665,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,p2048,2048,EXPONENTIATION,97,0.267071841,363.198,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p2048,2048,EXPONENTIATION,15400,46.253440989,332.948,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p2048,2048,EXPONENTIATION,15400,46.206609562,333.286,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p2048,2048,EXPONENTIATION,15400,7.674208960,2006.721,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p2048,2048,EXPONENTIATION,15400,7.677863450,2005.766,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p2048,2048,EXPONENTIATION,15400,2.242512440,6867.297,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p2048,2048,EXPONENTIATION,15400,2.246588375,6854.838,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p2048,2048,EXPONENTIATION,15400,2.166167747,7109.329,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p2048,2048,EXPONENTIATION,15400,2.165618911,7111.131,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p2048,2048,EXPONENTIATION,15400,2.201539184,6995.106,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p2048,2048,EXPONENTIATION,15400,2.209595568,6969.601,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p2048,2048,EXPONENTIATION,15400,2.172866523,7087.412,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p2048,2048,EXPONENTIATION,15400,2.176366532,7076.014,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,p2048,2048,DIVIDE,781,0.000068957,11325901.036,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,p2048,2048,DIVIDE,781,0.000034789,22449642.500,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,p2048,2048,DIVIDE,781,0.000176466,4425782.898,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,p2048,2048,DIVIDE,50000,0.000290816,171930017.606,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p2048,2048,DIVIDE,15400,1.878366854,8198.611,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p2048,2048,DIVIDE,15400,1.880863950,8187.727,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p2048,2048,DIVIDE,15400,0.393200873,39165.732,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p2048,2048,DIVIDE,15400,0.395947117,38894.083,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p2048,2048,DIVIDE,15400,0.034753930,443115.353,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p2048,2048,DIVIDE,15400,0.037920648,406111.204,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p2048,2048,DIVIDE,15400,0.032199857,478262.994,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p2048,2048,DIVIDE,15400,0.035175915,437799.557,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p2048,2048,DIVIDE,15400,0.032823787,469171.944,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p2048,2048,DIVIDE,15400,0.035989357,427904.282,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p2048,2048,DIVIDE,15400,0.031881033,483045.829,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p2048,2048,DIVIDE,15400,0.035001157,439985.454,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,p2048,2048,ISQRT,195,0.000170660,1142622.944,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,p2048,2048,ISQRT,195,0.000061684,3161270.756,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p2048,2048,ISQRT,15400,24.920779523,617.958,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p2048,2048,ISQRT,15400,24.897989725,618.524,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p2048,2048,ISQRT,15400,9.068271493,1698.229,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p2048,2048,ISQRT,15400,9.066793493,1698.506,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p2048,2048,ISQRT,15400,0.117686918,130855.666,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p2048,2048,ISQRT,15400,0.120108925,128216.950,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p2048,2048,ISQRT,15400,0.101694710,151433.639,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p2048,2048,ISQRT,15400,0.104231442,147748.124,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p2048,2048,ISQRT,15400,0.117573281,130982.140,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p2048,2048,ISQRT,15400,0.119999499,128333.869,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p2048,2048,ISQRT,15400,0.102692560,149962.178,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p2048,2048,ISQRT,15400,0.105105733,146519.125,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-1t,p2048,2048,MODMUL_R2,6250,0.015964649,391489.972,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,gmp-nt,p2048,2048,MODMUL_R2,6250,0.007046704,886939.489,0
library,Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz,host-cpu,openssl-nt,p2048,2048,MODMUL_R2,6250,0.015841321,394537.805,0
library,NVIDIA GeForce GTX 1660 SUPER,gpu,cgbn,p2048,2048,MODMUL_R2,50000,0.001201152,41626705.030,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p2048,2048,MODMUL_R2,15400,0.074668046,206246.190,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w8,p2048,2048,MODMUL_R2,15400,0.077436590,198872.394,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p2048,2048,MODMUL_R2,15400,0.002816275,5468216.004,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w16,p2048,2048,MODMUL_R2,15400,0.005262405,2926418.608,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p2048,2048,MODMUL_R2,15400,0.001417862,10861423.607,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-opt,p2048,2048,MODMUL_R2,15400,0.003912213,3936391.080,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p2048,2048,MODMUL_R2,15400,0.001065319,14455763.211,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-o64,p2048,2048,MODMUL_R2,15400,0.003542966,4346640.525,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p2048,2048,MODMUL_R2,15400,0.001302042,11827575.049,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il,p2048,2048,MODMUL_R2,15400,0.003837848,4012665.467,0
opencl-kernel,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p2048,2048,MODMUL_R2,15400,0.000963998,15975136.428,0
opencl-e2e,NVIDIA GeForce GTX 1660 SUPER,GPU,w32-il64,p2048,2048,MODMUL_R2,15400,0.003420006,4502916.206,0
```
