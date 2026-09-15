# MPA-OpenCL benchmark report - NVIDIA A40

> **Note.** The multi-threaded GMP and OpenSSL baseline columns have been
> removed from this report: they predate the 2026-09-12 timing fix and were
> understated (see `reports/README.md`). The single-threaded GMP column, the
> OpenCL-on-CPU rows and all MPA measurements are unaffected and were verified
> against GMP before timing.


> **Partial report.** The run was interrupted or hit its time budget.
> Rows that never ran are marked `n/a`.

## 1. System under test

2 OpenCL device(s) exercised with the identical kernels and operands.

### Device 0 - NVIDIA A40 (GPU)

| Property | Value |
|---|---|
| Model | NVIDIA A40 |
| Type | GPU |
| Vendor | NVIDIA Corporation |
| Device memory | 44.43 GiB |
| Max single allocation | 11.11 GiB |
| Local memory | 48 KiB |
| Global cache | 2352 KiB |
| Compute units | 84 |
| Max clock | 1740 MHz |
| Max work-group size | 1024 |
| OpenCL version | OpenCL 3.0 CUDA |
| Driver | 570.86.10 |

### Device 1 - cpu-haswell-AMD Eng Sample: 100-000000020-02_30/20_N (CPU)

| Property | Value |
|---|---|
| Model | cpu-haswell-AMD Eng Sample: 100-000000020-02_30/20_N |
| Type | CPU |
| Vendor | AuthenticAMD |
| Device memory | 249.52 GiB |
| Max single allocation | 64.00 GiB |
| Local memory | 512 KiB |
| Global cache | 16384 KiB |
| Compute units | 128 |
| Max clock | 2000 MHz |
| Max work-group size | 4096 |
| OpenCL version | OpenCL 3.0 PoCL HSTR: cpu-x86_64-pc-linux-gnu-haswell |
| Driver | 5.0+debian |

### Host

| Property | Value |
|---|---|
| CPU | AMD Eng Sample: 100-000000020-02_30/20_N |
| Logical cores | 128 |
| OpenMP threads used | 128 |
| RAM | 251.5 GB |
| OS | Ubuntu 24.04.4 LTS |
| Kernel | 6.8.0-65-generic |
| Arch | x86_64 |
| GMP | 6.3.0 |
| OpenSSL | OpenSSL 3.0.13 30 Jan 2024 |
| CGBN | cgbn_results.tsv loaded |

## 2. Method

- Workload auto-sized from the device and host: --min-items from 700 x compute units, --items from ten times that capped by host RAM. Either flag, given explicitly, overrides its half.
- Base workload 50000 items, scaled down per operator by its cost weight and by modulus size. Device rows honour --min-items (58800) so the GPU is not left idle; the CPU libraries keep the smaller count because a full-width MODEXP there costs minutes. Both counts appear in every row as dev/cpu, and throughput is per-second so they remain comparable.
- 5 timed repetitions, **minimum** reported. Two untimed warm-up launches precede them.
- `kernel` times `clEnqueueNDRangeKernel` + `clFinish` only. `e2e` adds the host->device operand writes and the device->host result read.
- Every OpenCL device runs the same kernels on the same operands, so GPU and CPU-OpenCL columns are directly comparable.
- CPU library baselines (GMP, OpenSSL) run those same operands, with temporaries preallocated outside the timed region, so the figure is the arithmetic and not marshalling. The generator is reseeded per modulus and operation so every backend sees identical inputs.
- OpenSSL rows time the nearest BN primitive, which is not always semantically identical (its Montgomery routine expects Montgomery-domain inputs); they measure comparable work, not identical results. Correctness is judged against GMP only.
- Every device cell is checked word-for-word against GMP before it is timed. A cell that mismatches is reported and excluded from the speedup tables.
- Total wall time 5401.7 s.

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
| [1] CPU | `mpaKernels_8bits.cl` (w8) | 11 | 11 | 0 | 0 |

**All configurations correct** - 496 configurations, 0 problems.

## 4. Throughput per device

Operations per second, higher is better. Kernel-only timings.

### Device 0 - NVIDIA A40 (GPU)

#### secp256k1 (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | 1.05 G | 1.69 G | 2.60 G | 2.63 G | 2.63 G | 4.11 G | 4.12 G | 52.25 M | 1.95 G |
| SUBTRACT | 50000 / 50000 | 1.03 G | 1.74 G | 2.80 G | 2.66 G | 2.72 G | 4.10 G | 4.14 G | 65.35 M | 2.03 G |
| ADDMOD | 50000 / 50000 | 716.00 M | 1.26 G | 2.39 G | 3.38 G | 3.32 G | 3.96 G | 4.12 G | 16.55 M | 2.03 G |
| SUBTRACTMOD | 50000 / 50000 | 718.26 M | 1.23 G | 2.37 G | 3.51 G | 3.28 G | 4.00 G | 4.12 G | 17.95 M | 2.03 G |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | 33.30 M | 110.10 M | 436.89 M | 1.64 G | 1.65 G | 3.23 G | 3.30 G | 41.46 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | 179.15 M | 516.13 M | 1.45 G | 1.41 G | 1.45 G | 1.92 G | 2.05 G | 41.36 M | 2.03 G |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | 359.60 M | 1.19 G | 2.93 G | 2.14 G | 2.14 G | 2.00 G | 2.11 G | 5.09 M | 1.96 G |
| COMPARE | 50000 / 50000 | 1.06 G | 1.79 G | - | 3.40 G | 3.15 G | 4.16 G | 4.00 G | 126.60 M | 2.03 G |
| REDUCE | 50000 / 6250 | 202.01 M | 317.87 M | - | 864.90 M | 882.79 M | 866.10 M | 872.32 M | 41.89 M | 2.03 G |
| MODMUL | 50000 / 3125 | 78.27 M | 132.84 M | - | 292.20 M | 354.49 M | 301.45 M | 352.59 M | 9.01 M | 813.80 M |
| MODEXP | 50000 / 781 | 1.73 M | 9.14 M | - | 14.16 M | 24.55 M | 14.11 M | 24.41 M | 87.15 k | 2.64 M |
| EXPONENTIATION | 50000 / 781 | 1.10 M | 4.14 M | - | 68.52 M | 73.35 M | 69.78 M | 77.60 M | 261.16 k | n/a |
| DIVIDE | 50000 / 6250 | 106.30 M | 162.64 M | - | 389.58 M | 391.75 M | 388.34 M | 394.22 M | 18.86 M | 1.88 G |
| ISQRT | 50000 / 1562 | 8.41 M | 11.05 M | - | 48.41 M | 55.28 M | 54.86 M | 54.83 M | 8.91 M | n/a |
| MODMUL_R2 | 50000 / 50000 | 317.68 M | 1.01 G | - | 1.39 G | 1.41 G | 1.38 G | 1.45 G | 8.60 M | 1.74 G |

#### rsa256(composite) (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | 1.04 G | 1.68 G | 2.73 G | 2.66 G | 2.64 G | 3.96 G | 4.09 G | 51.08 M | 2.03 G |
| SUBTRACT | 50000 / 50000 | 1.04 G | 1.72 G | 2.63 G | 2.67 G | 2.75 G | 4.17 G | 4.09 G | 67.37 M | 2.03 G |
| ADDMOD | 50000 / 50000 | 775.28 M | 1.36 G | 2.52 G | 3.36 G | 3.28 G | 4.04 G | 4.07 G | 19.04 M | 1.95 G |
| SUBTRACTMOD | 50000 / 50000 | 716.31 M | 1.25 G | 2.37 G | 3.57 G | 3.31 G | 3.88 G | 4.18 G | 17.97 M | 1.95 G |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | 33.22 M | 110.23 M | 435.55 M | 1.65 G | 1.67 G | 3.22 G | 3.37 G | 41.72 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | 179.08 M | 516.19 M | 1.46 G | 1.43 G | 1.47 G | 1.91 G | 2.02 G | 41.56 M | 2.03 G |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | 357.59 M | 1.18 G | 2.96 G | 2.12 G | 2.10 G | 2.01 G | 2.09 G | 5.12 M | 2.03 G |
| COMPARE | 50000 / 50000 | 1.07 G | 1.73 G | - | 3.48 G | 3.17 G | 4.17 G | 4.13 G | 125.93 M | 2.12 G |
| REDUCE | 50000 / 6250 | 201.97 M | 315.40 M | - | 852.22 M | 867.77 M | 858.38 M | 870.95 M | 25.18 M | 1.95 G |
| MODMUL | 50000 / 3125 | 77.92 M | 132.83 M | - | 292.87 M | 354.24 M | 302.20 M | 352.16 M | 8.87 M | 813.80 M |
| MODEXP | 50000 / 781 | 1.73 M | 9.15 M | - | 14.17 M | 24.57 M | 14.14 M | 24.50 M | 90.64 k | 2.67 M |
| EXPONENTIATION | 50000 / 781 | 1.10 M | 4.13 M | - | 67.64 M | 73.66 M | 69.79 M | 77.40 M | 258.53 k | n/a |
| DIVIDE | 50000 / 6250 | 105.97 M | 159.19 M | - | 374.18 M | 377.47 M | 375.17 M | 382.73 M | 18.98 M | 1.88 G |
| ISQRT | 50000 / 1562 | 8.41 M | 11.04 M | - | 47.79 M | 54.89 M | 54.71 M | 54.50 M | 10.41 M | n/a |
| MODMUL_R2 | 50000 / 50000 | 317.56 M | 983.14 M | - | 1.42 G | 1.43 G | 1.36 G | 1.43 G | 8.95 M | 1.74 G |

#### brainpoolP512r1 (512-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 25000 | 492.31 M | 845.84 M | 1.22 G | 1.28 G | 1.28 G | 1.61 G | 1.45 G | 47.91 M | 1.44 G |
| SUBTRACT | 50000 / 25000 | 492.55 M | 849.59 M | 1.21 G | 1.28 G | 1.29 G | 1.61 G | 1.48 G | 60.97 M | 1.40 G |
| ADDMOD | 50000 / 25000 | 364.43 M | 659.42 M | 1.06 G | 1.04 G | 1.01 G | 1.40 G | 1.57 G | 17.89 M | 1.44 G |
| SUBTRACTMOD | 50000 / 25000 | 319.54 M | 593.48 M | 998.90 M | 1.02 G | 1.02 G | 1.38 G | 1.54 G | 17.12 M | 1.44 G |
| MULTIPLYOPERANDSCANNING | 50000 / 25000 | 5.76 M | 22.03 M | 15.52 M | 400.68 M | 393.48 M | 535.63 M | 683.16 M | 17.18 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 25000 | 26.07 M | 94.62 M | 175.60 M | 287.75 M | 287.97 M | 559.72 M | 463.20 M | 17.12 M | 1.44 G |
| MONTGOMERYMULTIPLICATION | 50000 / 25000 | 107.49 M | 331.53 M | 1.09 G | 732.17 M | 873.09 M | 763.54 M | 975.07 M | 2.19 M | 1.40 G |
| COMPARE | 50000 / 25000 | 492.99 M | 851.47 M | - | 1.44 G | 1.44 G | 3.52 G | 3.67 G | 125.11 M | 1.44 G |
| REDUCE | 50000 / 3125 | 70.23 M | 88.75 M | - | 318.62 M | 318.01 M | 324.38 M | 292.84 M | 24.77 M | 1.44 G |
| MODMUL | 50000 / 1562 | 24.65 M | 35.13 M | - | 87.67 M | 108.70 M | 86.97 M | 90.02 M | 4.08 M | 259.72 M |
| MODEXP | 50000 / 390 | 118.61 k | 1.30 M | - | 1.59 M | 2.55 M | 1.52 M | 3.43 M | 17.09 k | 1.08 M |
| EXPONENTIATION | 50000 / 390 | 139.10 k | 523.58 k | - | 1.51 M | 1.30 M | 1.51 M | 1.97 M | 74.73 k | n/a |
| DIVIDE | 50000 / 3125 | 36.09 M | 40.27 M | - | 114.74 M | 119.40 M | 119.50 M | 121.75 M | 17.15 M | 1.09 G |
| ISQRT | 50000 / 781 | 1.77 M | 1.74 M | - | 8.39 M | 8.89 M | 8.40 M | 7.55 M | 5.15 M | n/a |
| MODMUL_R2 | 50000 / 25000 | 70.96 M | 352.31 M | - | 466.36 M | 567.68 M | 468.29 M | 627.81 M | 4.46 M | 1.19 G |

#### p1024 (1024-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 12500 | 181.52 M | 326.90 M | 107.47 M | 553.33 M | 551.20 M | 852.35 M | 836.50 M | 39.92 M | 841.86 M |
| SUBTRACT | 50000 / 12500 | 177.67 M | 322.57 M | 107.37 M | 560.86 M | 551.25 M | 845.00 M | 837.06 M | 50.34 M | 841.86 M |
| ADDMOD | 50000 / 12500 | 119.56 M | 229.92 M | 109.57 M | 444.11 M | 441.44 M | 800.65 M | 784.05 M | 13.54 M | 841.86 M |
| SUBTRACTMOD | 50000 / 12500 | 120.46 M | 228.70 M | 108.08 M | 442.97 M | 441.09 M | 805.82 M | 813.71 M | 14.82 M | 841.86 M |
| MULTIPLYOPERANDSCANNING | 50000 / 12500 | 1.20 M | 3.72 M | 2.51 M | 169.92 M | 167.39 M | 218.90 M | 222.20 M | 4.91 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 12500 | 3.46 M | 13.55 M | 26.36 M | 52.12 M | 52.09 M | 121.44 M | 131.25 M | 4.93 M | 831.12 M |
| MONTGOMERYMULTIPLICATION | 50000 / 12500 | 19.64 M | 100.00 M | 350.31 M | 230.49 M | 285.81 M | 269.41 M | 359.24 M | 718.77 k | 856.63 M |
| COMPARE | 50000 / 12500 | 198.82 M | 375.99 M | - | 666.37 M | 658.89 M | 1.59 G | 1.65 G | 108.53 M | 842.32 M |
| REDUCE | 50000 / 1562 | 13.21 M | 24.28 M | - | 95.79 M | 94.60 M | 110.75 M | 100.02 M | 35.00 M | 841.86 M |
| MODMUL | 50000 / 781 | 3.25 M | 9.14 M | - | 20.62 M | 27.51 M | 21.15 M | 27.89 M | 1.56 M | 106.17 M |
| MODEXP | 50000 / 195 | 14.06 k | 96.65 k | - | 212.53 k | 335.75 k | 210.12 k | 333.92 k | 2.59 k | 206.96 k |
| EXPONENTIATION | 50000 / 195 | 17.88 k | 67.84 k | - | 246.59 k | 262.20 k | 243.41 k | 257.73 k | 17.50 k | n/a |
| DIVIDE | 50000 / 1562 | 561.72 k | 5.21 M | - | 27.74 M | 28.54 M | 27.80 M | 28.36 M | 16.23 M | 814.23 M |
| ISQRT | 50000 / 390 | 61.82 k | 245.42 k | - | 1.26 M | 1.34 M | 1.27 M | 1.37 M | 2.57 M | n/a |
| MODMUL_R2 | 50000 / 12500 | 12.16 M | 86.55 M | - | 145.67 M | 173.03 M | 154.84 M | 191.95 M | 1.60 M | 469.50 M |

#### p2048 (2048-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 6250 | 88.49 M | 170.92 M | 256.09 M | 268.61 M | 265.83 M | 467.19 M | 463.80 M | 29.63 M | 447.96 M |
| SUBTRACT | 50000 / 6250 | 82.71 M | 170.70 M | 262.11 M | 268.28 M | 269.19 M | 462.55 M | 462.00 M | 35.15 M | 444.40 M |
| ADDMOD | 50000 / 6250 | 65.37 M | 128.75 M | 225.16 M | 206.40 M | 221.93 M | 369.28 M | 374.66 M | 10.36 M | 447.96 M |
| SUBTRACTMOD | 50000 / 6250 | 61.27 M | 117.93 M | 207.13 M | 224.80 M | 225.93 M | 360.01 M | 383.30 M | 12.32 M | 452.37 M |
| MULTIPLYOPERANDSCANNING | 50000 / 6250 | 233.17 k | 619.32 k | 2.12 M | 31.77 M | 30.82 M | 31.55 M | 30.71 M | 1.43 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 6250 | 874.88 k | 3.42 M | 13.07 M | 11.74 M | 11.79 M | 24.76 M | 26.12 M | 1.45 M | 417.33 M |
| MONTGOMERYMULTIPLICATION | 50000 / 6250 | 453.08 k | 22.25 M | 116.32 M | 74.44 M | 83.08 M | 74.03 M | 83.08 M | 222.41 k | 283.88 M |
| COMPARE | 50000 / 6250 | 98.04 M | 199.72 M | - | 367.84 M | 367.87 M | 836.08 M | 832.72 M | 123.16 M | 447.96 M |
| REDUCE | 50000 / 781 | 58.58 k | 5.63 M | - | 29.30 M | 27.96 M | 30.87 M | 27.37 M | 23.48 M | 443.89 M |
| MODMUL | 50000 / 390 | 33.95 k | 2.06 M | - | 5.13 M | 6.53 M | 5.45 M | 6.80 M | 522.13 k | 58.48 M |
| MODEXP | 50000 / 97 | 315.8 | 1.91 k | - | 24.70 k | 5.08 k | 25.26 k | 5.09 k | 355.0 | 37.22 k |
| EXPONENTIATION | 50000 / 97 | 1.59 k | 8.25 k | - | 27.66 k | 29.42 k | 27.84 k | 28.80 k | 2.97 k | n/a |
| DIVIDE | 50000 / 781 | 21.99 k | 102.04 k | - | 1.13 M | 1.35 M | 1.22 M | 1.31 M | 12.34 M | 443.89 M |
| ISQRT | 50000 / 195 | 1.64 k | 6.22 k | - | 497.04 k | 577.56 k | 492.74 k | 575.79 k | 1.65 M | n/a |
| MODMUL_R2 | 50000 / 6250 | 684.76 k | 18.01 M | - | 42.50 M | 52.29 M | 43.22 M | 53.54 M | 516.80 k | 148.87 M |

### Device 1 - cpu-haswell-AMD Eng Sample: 100-000000020-02_30/20_N (CPU)

#### secp256k1 (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | 67.36 M | - | - | - | - | - | - | 52.25 M | 1.95 G |
| SUBTRACT | 50000 / 50000 | 72.98 M | - | - | - | - | - | - | 65.35 M | 2.03 G |
| ADDMOD | 50000 / 50000 | 65.50 M | - | - | - | - | - | - | 16.55 M | 2.03 G |
| SUBTRACTMOD | 50000 / 50000 | 74.48 M | - | - | - | - | - | - | 17.95 M | 2.03 G |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | 16.47 M | - | - | - | - | - | - | 41.46 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | 31.87 M | - | - | - | - | - | - | 41.36 M | 2.03 G |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | 7.72 M | - | - | - | - | - | - | 5.09 M | 1.96 G |
| COMPARE | 50000 / 50000 | 63.75 M | - | - | - | - | - | - | 126.60 M | 2.03 G |
| REDUCE | 50000 / 6250 | 10.58 M | - | - | - | - | - | - | 41.89 M | 2.03 G |
| MODMUL | 50000 / 3125 | 2.81 M | - | - | - | - | - | - | 9.01 M | 813.80 M |
| MODEXP | 50000 / 781 | over budget | - | - | - | - | - | - | 87.15 k | 2.64 M |
| EXPONENTIATION | 50000 / 781 | - | - | - | - | - | - | - | 261.16 k | n/a |
| DIVIDE | 50000 / 6250 | - | - | - | - | - | - | - | 18.86 M | 1.88 G |
| ISQRT | 50000 / 1562 | - | - | - | - | - | - | - | 8.91 M | n/a |
| MODMUL_R2 | 50000 / 50000 | - | - | - | - | - | - | - | 8.60 M | 1.74 G |

#### rsa256(composite) (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | - | - | - | - | - | - | - | 51.08 M | 2.03 G |
| SUBTRACT | 50000 / 50000 | - | - | - | - | - | - | - | 67.37 M | 2.03 G |
| ADDMOD | 50000 / 50000 | - | - | - | - | - | - | - | 19.04 M | 1.95 G |
| SUBTRACTMOD | 50000 / 50000 | - | - | - | - | - | - | - | 17.97 M | 1.95 G |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | - | - | - | - | - | - | - | 41.72 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | - | - | - | - | - | - | - | 41.56 M | 2.03 G |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | - | - | - | - | - | - | - | 5.12 M | 2.03 G |
| COMPARE | 50000 / 50000 | - | - | - | - | - | - | - | 125.93 M | 2.12 G |
| REDUCE | 50000 / 6250 | - | - | - | - | - | - | - | 25.18 M | 1.95 G |
| MODMUL | 50000 / 3125 | - | - | - | - | - | - | - | 8.87 M | 813.80 M |
| MODEXP | 50000 / 781 | - | - | - | - | - | - | - | 90.64 k | 2.67 M |
| EXPONENTIATION | 50000 / 781 | - | - | - | - | - | - | - | 258.53 k | n/a |
| DIVIDE | 50000 / 6250 | - | - | - | - | - | - | - | 18.98 M | 1.88 G |
| ISQRT | 50000 / 1562 | - | - | - | - | - | - | - | 10.41 M | n/a |
| MODMUL_R2 | 50000 / 50000 | - | - | - | - | - | - | - | 8.95 M | 1.74 G |

#### brainpoolP512r1 (512-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 25000 | - | - | - | - | - | - | - | 47.91 M | 1.44 G |
| SUBTRACT | 50000 / 25000 | - | - | - | - | - | - | - | 60.97 M | 1.40 G |
| ADDMOD | 50000 / 25000 | - | - | - | - | - | - | - | 17.89 M | 1.44 G |
| SUBTRACTMOD | 50000 / 25000 | - | - | - | - | - | - | - | 17.12 M | 1.44 G |
| MULTIPLYOPERANDSCANNING | 50000 / 25000 | - | - | - | - | - | - | - | 17.18 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 25000 | - | - | - | - | - | - | - | 17.12 M | 1.44 G |
| MONTGOMERYMULTIPLICATION | 50000 / 25000 | - | - | - | - | - | - | - | 2.19 M | 1.40 G |
| COMPARE | 50000 / 25000 | - | - | - | - | - | - | - | 125.11 M | 1.44 G |
| REDUCE | 50000 / 3125 | - | - | - | - | - | - | - | 24.77 M | 1.44 G |
| MODMUL | 50000 / 1562 | - | - | - | - | - | - | - | 4.08 M | 259.72 M |
| MODEXP | 50000 / 390 | - | - | - | - | - | - | - | 17.09 k | 1.08 M |
| EXPONENTIATION | 50000 / 390 | - | - | - | - | - | - | - | 74.73 k | n/a |
| DIVIDE | 50000 / 3125 | - | - | - | - | - | - | - | 17.15 M | 1.09 G |
| ISQRT | 50000 / 781 | - | - | - | - | - | - | - | 5.15 M | n/a |
| MODMUL_R2 | 50000 / 25000 | - | - | - | - | - | - | - | 4.46 M | 1.19 G |

#### p1024 (1024-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 12500 | - | - | - | - | - | - | - | 39.92 M | 841.86 M |
| SUBTRACT | 50000 / 12500 | - | - | - | - | - | - | - | 50.34 M | 841.86 M |
| ADDMOD | 50000 / 12500 | - | - | - | - | - | - | - | 13.54 M | 841.86 M |
| SUBTRACTMOD | 50000 / 12500 | - | - | - | - | - | - | - | 14.82 M | 841.86 M |
| MULTIPLYOPERANDSCANNING | 50000 / 12500 | - | - | - | - | - | - | - | 4.91 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 12500 | - | - | - | - | - | - | - | 4.93 M | 831.12 M |
| MONTGOMERYMULTIPLICATION | 50000 / 12500 | - | - | - | - | - | - | - | 718.77 k | 856.63 M |
| COMPARE | 50000 / 12500 | - | - | - | - | - | - | - | 108.53 M | 842.32 M |
| REDUCE | 50000 / 1562 | - | - | - | - | - | - | - | 35.00 M | 841.86 M |
| MODMUL | 50000 / 781 | - | - | - | - | - | - | - | 1.56 M | 106.17 M |
| MODEXP | 50000 / 195 | - | - | - | - | - | - | - | 2.59 k | 206.96 k |
| EXPONENTIATION | 50000 / 195 | - | - | - | - | - | - | - | 17.50 k | n/a |
| DIVIDE | 50000 / 1562 | - | - | - | - | - | - | - | 16.23 M | 814.23 M |
| ISQRT | 50000 / 390 | - | - | - | - | - | - | - | 2.57 M | n/a |
| MODMUL_R2 | 50000 / 12500 | - | - | - | - | - | - | - | 1.60 M | 469.50 M |

#### p2048 (2048-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 6250 | - | - | - | - | - | - | - | 29.63 M | 447.96 M |
| SUBTRACT | 50000 / 6250 | - | - | - | - | - | - | - | 35.15 M | 444.40 M |
| ADDMOD | 50000 / 6250 | - | - | - | - | - | - | - | 10.36 M | 447.96 M |
| SUBTRACTMOD | 50000 / 6250 | - | - | - | - | - | - | - | 12.32 M | 452.37 M |
| MULTIPLYOPERANDSCANNING | 50000 / 6250 | - | - | - | - | - | - | - | 1.43 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 6250 | - | - | - | - | - | - | - | 1.45 M | 417.33 M |
| MONTGOMERYMULTIPLICATION | 50000 / 6250 | - | - | - | - | - | - | - | 222.41 k | 283.88 M |
| COMPARE | 50000 / 6250 | - | - | - | - | - | - | - | 123.16 M | 447.96 M |
| REDUCE | 50000 / 781 | - | - | - | - | - | - | - | 23.48 M | 443.89 M |
| MODMUL | 50000 / 390 | - | - | - | - | - | - | - | 522.13 k | 58.48 M |
| MODEXP | 50000 / 97 | - | - | - | - | - | - | - | 355.0 | 37.22 k |
| EXPONENTIATION | 50000 / 97 | - | - | - | - | - | - | - | 2.97 k | n/a |
| DIVIDE | 50000 / 781 | - | - | - | - | - | - | - | 12.34 M | 443.89 M |
| ISQRT | 50000 / 195 | - | - | - | - | - | - | - | 1.65 M | n/a |
| MODMUL_R2 | 50000 / 6250 | - | - | - | - | - | - | - | 516.80 k | 148.87 M |

## 5. Head to head

Best OpenCL GPU result against best OpenCL CPU result and the CPU libraries.
Ratios above 1.00x mean the GPU is faster than that baseline.

### secp256k1 (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w32-il64 | 4.12 G | w8 | 67.36 M | 52.25 M | 1.95 G | 61.13x | 2.11x |
| SUBTRACT | w32-il64 | 4.14 G | w8 | 72.98 M | 65.35 M | 2.03 G | 56.70x | 2.03x |
| ADDMOD | w32-il64 | 4.12 G | w8 | 65.50 M | 16.55 M | 2.03 G | 62.92x | 2.03x |
| SUBTRACTMOD | w32-il64 | 4.12 G | w8 | 74.48 M | 17.95 M | 2.03 G | 55.38x | 2.03x |
| MULTIPLYOPERANDSCANNING | w32-il64 | 3.30 G | w8 | 16.47 M | 41.46 M | n/a | 200.31x | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 2.05 G | w8 | 31.87 M | 41.36 M | 2.03 G | 64.41x | 1.01x |
| MONTGOMERYMULTIPLICATION | w32 | 2.93 G | w8 | 7.72 M | 5.09 M | 1.96 G | 379.72x | 1.50x |
| COMPARE | w32-il | 4.16 G | w8 | 63.75 M | 126.60 M | 2.03 G | 65.24x | 2.04x |
| REDUCE | w32-o64 | 110.35 M | w8 | 1.32 M | 41.89 M | 2.03 G | 83.44x | 0.05x |
| MODMUL | w32-o64 | 22.16 M | w8 | 175.48 k | 9.01 M | 813.80 M | 126.25x | 0.03x |
| MODEXP | w32-o64 | 383.47 k | none | n/a | 87.15 k | 2.64 M | n/a | 0.15x |
| EXPONENTIATION | w32-il64 | 1.21 M | none | n/a | 261.16 k | n/a | n/a | n/a |
| DIVIDE | w32-il64 | 49.28 M | none | n/a | 18.86 M | 1.88 G | n/a | 0.03x |
| ISQRT | w32-o64 | 1.73 M | none | n/a | 8.91 M | n/a | n/a | n/a |
| MODMUL_R2 | w32-il64 | 1.45 G | none | n/a | 8.60 M | 1.74 G | n/a | 0.83x |

### rsa256(composite) (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w32-il64 | 4.09 G | none | n/a | 51.08 M | 2.03 G | n/a | 2.01x |
| SUBTRACT | w32-il | 4.17 G | none | n/a | 67.37 M | 2.03 G | n/a | 2.05x |
| ADDMOD | w32-il64 | 4.07 G | none | n/a | 19.04 M | 1.95 G | n/a | 2.08x |
| SUBTRACTMOD | w32-il64 | 4.18 G | none | n/a | 17.97 M | 1.95 G | n/a | 2.14x |
| MULTIPLYOPERANDSCANNING | w32-il64 | 3.37 G | none | n/a | 41.72 M | n/a | n/a | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 2.02 G | none | n/a | 41.56 M | 2.03 G | n/a | 0.99x |
| MONTGOMERYMULTIPLICATION | w32 | 2.96 G | none | n/a | 5.12 M | 2.03 G | n/a | 1.45x |
| COMPARE | w32-il | 4.17 G | none | n/a | 125.93 M | 2.12 G | n/a | 1.96x |
| REDUCE | w32-il64 | 108.87 M | none | n/a | 25.18 M | 1.95 G | n/a | 0.06x |
| MODMUL | w32-o64 | 22.14 M | none | n/a | 8.87 M | 813.80 M | n/a | 0.03x |
| MODEXP | w32-o64 | 383.75 k | none | n/a | 90.64 k | 2.67 M | n/a | 0.14x |
| EXPONENTIATION | w32-il64 | 1.21 M | none | n/a | 258.53 k | n/a | n/a | n/a |
| DIVIDE | w32-il64 | 47.84 M | none | n/a | 18.98 M | 1.88 G | n/a | 0.03x |
| ISQRT | w32-o64 | 1.71 M | none | n/a | 10.41 M | n/a | n/a | n/a |
| MODMUL_R2 | w32-il64 | 1.43 G | none | n/a | 8.95 M | 1.74 G | n/a | 0.82x |

### brainpoolP512r1 (512-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w32-il | 803.36 M | none | n/a | 47.91 M | 1.44 G | n/a | 0.56x |
| SUBTRACT | w32-il | 806.48 M | none | n/a | 60.97 M | 1.40 G | n/a | 0.58x |
| ADDMOD | w32-il64 | 784.90 M | none | n/a | 17.89 M | 1.44 G | n/a | 0.55x |
| SUBTRACTMOD | w32-il64 | 771.59 M | none | n/a | 17.12 M | 1.44 G | n/a | 0.54x |
| MULTIPLYOPERANDSCANNING | w32-il64 | 341.58 M | none | n/a | 17.18 M | n/a | n/a | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il | 279.86 M | none | n/a | 17.12 M | 1.44 G | n/a | 0.19x |
| MONTGOMERYMULTIPLICATION | w32 | 546.72 M | none | n/a | 2.19 M | 1.40 G | n/a | 0.39x |
| COMPARE | w32-il64 | 1.83 G | none | n/a | 125.11 M | 1.44 G | n/a | 1.28x |
| REDUCE | w32-il | 20.27 M | none | n/a | 24.77 M | 1.44 G | n/a | 0.01x |
| MODMUL | w32-o64 | 3.40 M | none | n/a | 4.08 M | 259.72 M | n/a | 0.01x |
| MODEXP | w32-il64 | 26.76 k | none | n/a | 17.09 k | 1.08 M | n/a | 0.02x |
| EXPONENTIATION | w32-il64 | 15.37 k | none | n/a | 74.73 k | n/a | n/a | n/a |
| DIVIDE | w32-il64 | 7.61 M | none | n/a | 17.15 M | 1.09 G | n/a | 0.01x |
| ISQRT | w32-o64 | 138.88 k | none | n/a | 5.15 M | n/a | n/a | n/a |
| MODMUL_R2 | w32-il64 | 313.90 M | none | n/a | 4.46 M | 1.19 G | n/a | 0.26x |

### p1024 (1024-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w32-il | 213.09 M | none | n/a | 39.92 M | 841.86 M | n/a | 0.25x |
| SUBTRACT | w32-il | 211.25 M | none | n/a | 50.34 M | 841.86 M | n/a | 0.25x |
| ADDMOD | w32-il | 200.16 M | none | n/a | 13.54 M | 841.86 M | n/a | 0.24x |
| SUBTRACTMOD | w32-il64 | 203.43 M | none | n/a | 14.82 M | 841.86 M | n/a | 0.24x |
| MULTIPLYOPERANDSCANNING | w32-il64 | 55.55 M | none | n/a | 4.91 M | n/a | n/a | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 32.81 M | none | n/a | 4.93 M | 831.12 M | n/a | 0.04x |
| MONTGOMERYMULTIPLICATION | w32-il64 | 89.81 M | none | n/a | 718.77 k | 856.63 M | n/a | 0.10x |
| COMPARE | w32-il64 | 413.68 M | none | n/a | 108.53 M | 842.32 M | n/a | 0.49x |
| REDUCE | w32-il | 3.46 M | none | n/a | 35.00 M | 841.86 M | n/a | 0.00x |
| MODMUL | w32-il64 | 435.67 k | none | n/a | 1.56 M | 106.17 M | n/a | 0.00x |
| MODEXP | w32-o64 | 1.31 k | none | n/a | 2.59 k | 206.96 k | n/a | 0.01x |
| EXPONENTIATION | w32-o64 | 1.02 k | none | n/a | 17.50 k | n/a | n/a | n/a |
| DIVIDE | w32-o64 | 891.70 k | none | n/a | 16.23 M | 814.23 M | n/a | 0.00x |
| ISQRT | w32-il64 | 10.65 k | none | n/a | 2.57 M | n/a | n/a | n/a |
| MODMUL_R2 | w32-il64 | 47.99 M | none | n/a | 1.60 M | 469.50 M | n/a | 0.10x |

### p2048 (2048-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w32-il | 58.40 M | none | n/a | 29.63 M | 447.96 M | n/a | 0.13x |
| SUBTRACT | w32-il | 57.82 M | none | n/a | 35.15 M | 444.40 M | n/a | 0.13x |
| ADDMOD | w32-il64 | 46.83 M | none | n/a | 10.36 M | 447.96 M | n/a | 0.10x |
| SUBTRACTMOD | w32-il64 | 47.91 M | none | n/a | 12.32 M | 452.37 M | n/a | 0.11x |
| MULTIPLYOPERANDSCANNING | w32-opt | 3.97 M | none | n/a | 1.43 M | n/a | n/a | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 3.26 M | none | n/a | 1.45 M | 417.33 M | n/a | 0.01x |
| MONTGOMERYMULTIPLICATION | w32 | 14.54 M | none | n/a | 222.41 k | 283.88 M | n/a | 0.05x |
| COMPARE | w32-il | 104.51 M | none | n/a | 123.16 M | 447.96 M | n/a | 0.23x |
| REDUCE | w32-il | 482.22 k | none | n/a | 23.48 M | 443.89 M | n/a | 0.00x |
| MODMUL | w32-il64 | 53.01 k | none | n/a | 522.13 k | 58.48 M | n/a | 0.00x |
| MODEXP | w32-il | 49.0 | none | n/a | 355.0 | 37.22 k | n/a | 0.00x |
| EXPONENTIATION | w32-o64 | 57.1 | none | n/a | 2.97 k | n/a | n/a | n/a |
| DIVIDE | w32-o64 | 21.03 k | none | n/a | 12.34 M | 443.89 M | n/a | 0.00x |
| ISQRT | w32-o64 | 2.25 k | none | n/a | 1.65 M | n/a | n/a | n/a |
| MODMUL_R2 | w32-il64 | 6.69 M | none | n/a | 516.80 k | 148.87 M | n/a | 0.04x |

## 6. Raw data

Also written to `NVIDIA_A40_Report.csv` for analysis.

```csv
kind,device,device_type,kernel,modulus,bits,operation,items,seconds,ops_per_sec,mismatches
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,secp256k1,256,ADD,50000,0.000956858,52254380.593,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,secp256k1,256,ADD,50000,0.012999140,3846408.382,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,secp256k1,256,ADD,50000,0.008996839,5557507.659,0
library,NVIDIA A40,gpu,cgbn,secp256k1,256,ADD,50000,0.000025600,1953125000.000,0
opencl-kernel,NVIDIA A40,GPU,w8,secp256k1,256,ADD,50000,0.000047810,1045798098.800,0
opencl-e2e,NVIDIA A40,GPU,w8,secp256k1,256,ADD,50000,0.000606262,82472577.342,0
opencl-kernel,NVIDIA A40,GPU,w16,secp256k1,256,ADD,50000,0.000029576,1690559284.567,0
opencl-e2e,NVIDIA A40,GPU,w16,secp256k1,256,ADD,50000,0.000601484,83127798.612,0
opencl-kernel,NVIDIA A40,GPU,w32,secp256k1,256,ADD,50000,0.000019256,2596589823.950,0
opencl-e2e,NVIDIA A40,GPU,w32,secp256k1,256,ADD,50000,0.000575172,86930550.887,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,secp256k1,256,ADD,50000,0.000019026,2627984296.833,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,secp256k1,256,ADD,50000,0.000596805,83779522.919,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,secp256k1,256,ADD,50000,0.000018986,2633527479.643,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,secp256k1,256,ADD,50000,0.000556297,89880099.578,0
opencl-kernel,NVIDIA A40,GPU,w32-il,secp256k1,256,ADD,50000,0.000012173,4107343829.852,0
opencl-e2e,NVIDIA A40,GPU,w32-il,secp256k1,256,ADD,50000,0.000523113,95581695.176,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,secp256k1,256,ADD,50000,0.000012143,4117739776.039,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,secp256k1,256,ADD,50000,0.000479821,104205501.510,0
opencl-kernel,cpu-haswell-AMD Eng Sample: 100-000000020-02_30/20_N,CPU,w8,secp256k1,256,ADD,50000,0.000742229,67364681.378,0
opencl-e2e,cpu-haswell-AMD Eng Sample: 100-000000020-02_30/20_N,CPU,w8,secp256k1,256,ADD,50000,0.002368077,21114180.933,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,secp256k1,256,SUBTRACT,50000,0.000765122,65349015.025,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,secp256k1,256,SUBTRACT,50000,0.007441795,6718808.715,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,secp256k1,256,SUBTRACT,50000,0.011996716,4167807.166,0
library,NVIDIA A40,gpu,cgbn,secp256k1,256,SUBTRACT,50000,0.000024576,2034505208.333,0
opencl-kernel,NVIDIA A40,GPU,w8,secp256k1,256,SUBTRACT,50000,0.000048632,1028133808.265,0
opencl-e2e,NVIDIA A40,GPU,w8,secp256k1,256,SUBTRACT,50000,0.000548081,91227467.940,0
opencl-kernel,NVIDIA A40,GPU,w16,secp256k1,256,SUBTRACT,50000,0.000028785,1736996609.292,0
opencl-e2e,NVIDIA A40,GPU,w16,secp256k1,256,SUBTRACT,50000,0.000609809,81992939.886,0
opencl-kernel,NVIDIA A40,GPU,w32,secp256k1,256,SUBTRACT,50000,0.000017875,2797222487.365,0
opencl-e2e,NVIDIA A40,GPU,w32,secp256k1,256,SUBTRACT,50000,0.000574460,87038223.669,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,secp256k1,256,SUBTRACT,50000,0.000018827,2655804659.906,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,secp256k1,256,SUBTRACT,50000,0.000547530,91319175.514,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,secp256k1,256,SUBTRACT,50000,0.000018375,2721089265.079,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,secp256k1,256,SUBTRACT,50000,0.000543713,91960341.688,0
opencl-kernel,NVIDIA A40,GPU,w32-il,secp256k1,256,SUBTRACT,50000,0.000012194,4100442312.686,0
opencl-e2e,NVIDIA A40,GPU,w32-il,secp256k1,256,SUBTRACT,50000,0.000501993,99602960.236,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,secp256k1,256,SUBTRACT,50000,0.000012083,4138052350.856,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,secp256k1,256,SUBTRACT,50000,0.000514126,97252352.098,0
opencl-kernel,cpu-haswell-AMD Eng Sample: 100-000000020-02_30/20_N,CPU,w8,secp256k1,256,SUBTRACT,50000,0.000685132,72978626.122,0
opencl-e2e,cpu-haswell-AMD Eng Sample: 100-000000020-02_30/20_N,CPU,w8,secp256k1,256,SUBTRACT,50000,0.002438301,20506080.044,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,secp256k1,256,ADDMOD,50000,0.003020577,16553128.029,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,secp256k1,256,ADDMOD,50000,0.009990335,5004837.420,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,secp256k1,256,ADDMOD,50000,0.008145561,6138313.003,0
library,NVIDIA A40,gpu,cgbn,secp256k1,256,ADDMOD,50000,0.000024576,2034505208.333,0
opencl-kernel,NVIDIA A40,GPU,w8,secp256k1,256,ADDMOD,50000,0.000069832,715999722.600,0
opencl-e2e,NVIDIA A40,GPU,w8,secp256k1,256,ADDMOD,50000,0.000599588,83390557.671,0
opencl-kernel,NVIDIA A40,GPU,w16,secp256k1,256,ADDMOD,50000,0.000039635,1261504093.237,0
opencl-e2e,NVIDIA A40,GPU,w16,secp256k1,256,ADDMOD,50000,0.000577637,86559552.942,0
opencl-kernel,NVIDIA A40,GPU,w32,secp256k1,256,ADDMOD,50000,0.000020950,2386623302.956,0
opencl-e2e,NVIDIA A40,GPU,w32,secp256k1,256,ADDMOD,50000,0.000584380,85560798.951,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,secp256k1,256,ADDMOD,50000,0.000014808,3376546616.352,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,secp256k1,256,ADDMOD,50000,0.000570794,87597332.615,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,secp256k1,256,ADDMOD,50000,0.000015069,3318114412.855,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,secp256k1,256,ADDMOD,50000,0.000572028,87408363.915,0
opencl-kernel,NVIDIA A40,GPU,w32-il,secp256k1,256,ADDMOD,50000,0.000012614,3963902185.470,0
opencl-e2e,NVIDIA A40,GPU,w32-il,secp256k1,256,ADDMOD,50000,0.000517373,96642079.474,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,secp256k1,256,ADDMOD,50000,0.000012132,4121216795.885,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,secp256k1,256,ADDMOD,50000,0.000485151,103060680.671,0
opencl-kernel,cpu-haswell-AMD Eng Sample: 100-000000020-02_30/20_N,CPU,w8,secp256k1,256,ADDMOD,50000,0.000763400,65496423.914,0
opencl-e2e,cpu-haswell-AMD Eng Sample: 100-000000020-02_30/20_N,CPU,w8,secp256k1,256,ADDMOD,50000,0.002228272,22438910.820,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,secp256k1,256,SUBTRACTMOD,50000,0.002785110,17952612.834,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,secp256k1,256,SUBTRACTMOD,50000,0.009937775,5031307.147,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,secp256k1,256,SUBTRACTMOD,50000,0.008364136,5977903.518,0
library,NVIDIA A40,gpu,cgbn,secp256k1,256,SUBTRACTMOD,50000,0.000024576,2034505208.333,0
opencl-kernel,NVIDIA A40,GPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.000069613,718260391.192,0
opencl-e2e,NVIDIA A40,GPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.000588418,84973609.429,0
opencl-kernel,NVIDIA A40,GPU,w16,secp256k1,256,SUBTRACTMOD,50000,0.000040648,1230085718.868,0
opencl-e2e,NVIDIA A40,GPU,w16,secp256k1,256,SUBTRACTMOD,50000,0.000592305,84415921.809,0
opencl-kernel,NVIDIA A40,GPU,w32,secp256k1,256,SUBTRACTMOD,50000,0.000021090,2370814360.786,0
opencl-e2e,NVIDIA A40,GPU,w32,secp256k1,256,SUBTRACTMOD,50000,0.000580894,86074248.794,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,secp256k1,256,SUBTRACTMOD,50000,0.000014246,3509648375.498,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,secp256k1,256,SUBTRACTMOD,50000,0.000578950,86363220.048,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,secp256k1,256,SUBTRACTMOD,50000,0.000015238,3281205916.147,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,secp256k1,256,SUBTRACTMOD,50000,0.000568881,87891890.228,0
opencl-kernel,NVIDIA A40,GPU,w32-il,secp256k1,256,SUBTRACTMOD,50000,0.000012514,3995467083.426,0
opencl-e2e,NVIDIA A40,GPU,w32-il,secp256k1,256,SUBTRACTMOD,50000,0.000508947,98241999.129,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,secp256k1,256,SUBTRACTMOD,50000,0.000012123,4124382822.463,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,secp256k1,256,SUBTRACTMOD,50000,0.000501904,99620703.079,0
opencl-kernel,cpu-haswell-AMD Eng Sample: 100-000000020-02_30/20_N,CPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.000671365,74475103.382,0
opencl-e2e,cpu-haswell-AMD Eng Sample: 100-000000020-02_30/20_N,CPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.002343301,21337423.473,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001205930,41461791.688,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.009000014,5555547.171,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.008140251,6142316.711,0
opencl-kernel,NVIDIA A40,GPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001501522,33299544.797,0
opencl-e2e,NVIDIA A40,GPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.002363028,21159292.187,0
opencl-kernel,NVIDIA A40,GPU,w16,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000454132,110100264.960,0
opencl-e2e,NVIDIA A40,GPU,w16,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001118183,44715394.456,0
opencl-kernel,NVIDIA A40,GPU,w32,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000114447,436885334.375,0
opencl-e2e,NVIDIA A40,GPU,w32,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000749654,66697424.264,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000030397,1644875492.509,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000660765,75669869.244,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000030257,1652520659.936,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000704478,70974489.575,0
opencl-kernel,NVIDIA A40,GPU,w32-il,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000015490,3227939586.340,0
opencl-e2e,NVIDIA A40,GPU,w32-il,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000629625,79412312.904,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000015159,3298340677.029,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000652119,76673171.219,0
opencl-kernel,cpu-haswell-AMD Eng Sample: 100-000000020-02_30/20_N,CPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.003036599,16465791.286,0
opencl-e2e,cpu-haswell-AMD Eng Sample: 100-000000020-02_30/20_N,CPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.005090809,9821620.847,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001208926,41359037.466,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000056046,892123352.000,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000107222,466321180.589,0
library,NVIDIA A40,gpu,cgbn,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000024576,2034505208.333,0
opencl-kernel,NVIDIA A40,GPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000279101,179146866.345,0
opencl-e2e,NVIDIA A40,GPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000914357,54683256.076,0
opencl-kernel,NVIDIA A40,GPU,w16,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000096874,516132700.110,0
opencl-e2e,NVIDIA A40,GPU,w16,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000770234,64915373.528,0
opencl-kernel,NVIDIA A40,GPU,w32,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000034596,1445260484.023,0
opencl-e2e,NVIDIA A40,GPU,w32,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000704427,70979650.516,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000035557,1406194274.339,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000674550,74123443.413,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000034366,1454934720.867,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000724074,69053674.423,0
opencl-kernel,NVIDIA A40,GPU,w32-il,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000025989,1923923712.596,0
opencl-e2e,NVIDIA A40,GPU,w32-il,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000628916,79501921.684,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000024358,2052729647.473,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000636789,78518942.240,0
opencl-kernel,cpu-haswell-AMD Eng Sample: 100-000000020-02_30/20_N,CPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001568859,31870286.545,0
opencl-e2e,cpu-haswell-AMD Eng Sample: 100-000000020-02_30/20_N,CPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.003585029,13946887.044,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.009823999,5089576.875,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000285633,175049857.840,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.007683815,6507184.677,0
library,NVIDIA A40,gpu,cgbn,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000025504,1960476787.955,0
opencl-kernel,NVIDIA A40,GPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000139045,359596854.613,0
opencl-e2e,NVIDIA A40,GPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000661387,75598691.561,0
opencl-kernel,NVIDIA A40,GPU,w16,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000042181,1185381007.264,0
opencl-e2e,NVIDIA A40,GPU,w16,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000575994,86806438.480,0
opencl-kernel,NVIDIA A40,GPU,w32,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000017053,2932118580.011,0
opencl-e2e,NVIDIA A40,GPU,w32,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000539805,92626073.694,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000023405,2136289491.067,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000542731,92126666.358,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000023405,2136289491.067,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000596643,83802146.912,0
opencl-kernel,NVIDIA A40,GPU,w32-il,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000024998,2000189679.967,0
opencl-e2e,NVIDIA A40,GPU,w32-il,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000536287,93233624.968,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000023746,2105623845.943,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000531198,94126779.459,0
opencl-kernel,cpu-haswell-AMD Eng Sample: 100-000000020-02_30/20_N,CPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.006475259,7721699.309,0
opencl-e2e,cpu-haswell-AMD Eng Sample: 100-000000020-02_30/20_N,CPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.008368063,5975098.571,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,secp256k1,256,COMPARE,50000,0.000394939,126601687.013,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,secp256k1,256,COMPARE,50000,0.000049152,1017244087.387,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,secp256k1,256,COMPARE,50000,0.000074161,674206846.666,0
library,NVIDIA A40,gpu,cgbn,secp256k1,256,COMPARE,50000,0.000024576,2034505208.333,0
opencl-kernel,NVIDIA A40,GPU,w8,secp256k1,256,COMPARE,50000,0.000047129,1060925839.855,0
opencl-e2e,NVIDIA A40,GPU,w8,secp256k1,256,COMPARE,50000,0.000580613,86115944.744,0
opencl-kernel,NVIDIA A40,GPU,w16,secp256k1,256,COMPARE,50000,0.000027914,1791241532.097,0
opencl-e2e,NVIDIA A40,GPU,w16,secp256k1,256,COMPARE,50000,0.000579040,86349746.196,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,secp256k1,256,COMPARE,50000,0.000014688,3404165316.086,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,secp256k1,256,COMPARE,50000,0.000521170,95937990.214,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,secp256k1,256,COMPARE,50000,0.000015880,3148618333.236,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,secp256k1,256,COMPARE,50000,0.000539646,92653408.809,0
opencl-kernel,NVIDIA A40,GPU,w32-il,secp256k1,256,COMPARE,50000,0.000012023,4158566320.682,0
opencl-e2e,NVIDIA A40,GPU,w32-il,secp256k1,256,COMPARE,50000,0.000526389,94986768.010,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,secp256k1,256,COMPARE,50000,0.000012494,4002019470.742,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,secp256k1,256,COMPARE,50000,0.000534436,93556618.326,0
opencl-kernel,cpu-haswell-AMD Eng Sample: 100-000000020-02_30/20_N,CPU,w8,secp256k1,256,COMPARE,50000,0.000784359,63746325.037,0
opencl-e2e,cpu-haswell-AMD Eng Sample: 100-000000020-02_30/20_N,CPU,w8,secp256k1,256,COMPARE,50000,0.002377185,21033280.613,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,secp256k1,256,REDUCE,6250,0.000149193,41891984.144,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,secp256k1,256,REDUCE,6250,0.000042922,145613435.459,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,secp256k1,256,REDUCE,6250,0.000085632,72986463.941,0
library,NVIDIA A40,gpu,cgbn,secp256k1,256,REDUCE,50000,0.000024576,2034505208.333,0
opencl-kernel,NVIDIA A40,GPU,w8,secp256k1,256,REDUCE,50000,0.000247510,202011917.430,0
opencl-e2e,NVIDIA A40,GPU,w8,secp256k1,256,REDUCE,50000,0.000768569,65055942.417,0
opencl-kernel,NVIDIA A40,GPU,w16,secp256k1,256,REDUCE,50000,0.000157299,317866944.546,0
opencl-e2e,NVIDIA A40,GPU,w16,secp256k1,256,REDUCE,50000,0.000706611,70760270.878,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,secp256k1,256,REDUCE,50000,0.000057810,864902472.895,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,secp256k1,256,REDUCE,50000,0.000607094,82359596.356,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,secp256k1,256,REDUCE,50000,0.000056638,882793573.954,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,secp256k1,256,REDUCE,50000,0.000635437,78686039.323,0
opencl-kernel,NVIDIA A40,GPU,w32-il,secp256k1,256,REDUCE,50000,0.000057730,866102427.928,0
opencl-e2e,NVIDIA A40,GPU,w32-il,secp256k1,256,REDUCE,50000,0.000611692,81740519.093,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,secp256k1,256,REDUCE,50000,0.000057318,872322547.729,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,secp256k1,256,REDUCE,50000,0.000610230,81936378.274,0
opencl-kernel,cpu-haswell-AMD Eng Sample: 100-000000020-02_30/20_N,CPU,w8,secp256k1,256,REDUCE,50000,0.004726127,10579487.648,0
opencl-e2e,cpu-haswell-AMD Eng Sample: 100-000000020-02_30/20_N,CPU,w8,secp256k1,256,REDUCE,50000,0.006461222,7738474.569,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,secp256k1,256,MODMUL,3125,0.000347019,9005266.110,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,secp256k1,256,MODMUL,3125,0.000048372,64603538.767,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,secp256k1,256,MODMUL,3125,0.000088138,35455933.726,0
library,NVIDIA A40,gpu,cgbn,secp256k1,256,MODMUL,50000,0.000061440,813802083.333,0
opencl-kernel,NVIDIA A40,GPU,w8,secp256k1,256,MODMUL,50000,0.000638844,78266428.166,0
opencl-e2e,NVIDIA A40,GPU,w8,secp256k1,256,MODMUL,50000,0.001184409,42215129.703,0
opencl-kernel,NVIDIA A40,GPU,w16,secp256k1,256,MODMUL,50000,0.000376404,132835902.524,0
opencl-e2e,NVIDIA A40,GPU,w16,secp256k1,256,MODMUL,50000,0.000927893,53885551.687,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,secp256k1,256,MODMUL,50000,0.000171115,292201679.611,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,secp256k1,256,MODMUL,50000,0.000713784,70049178.257,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,secp256k1,256,MODMUL,50000,0.000141048,354489572.067,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,secp256k1,256,MODMUL,50000,0.000674231,74158562.332,0
opencl-kernel,NVIDIA A40,GPU,w32-il,secp256k1,256,MODMUL,50000,0.000165865,301450291.977,0
opencl-e2e,NVIDIA A40,GPU,w32-il,secp256k1,256,MODMUL,50000,0.000716962,69738710.784,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,secp256k1,256,MODMUL,50000,0.000141809,352587519.210,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,secp256k1,256,MODMUL,50000,0.000699419,71487947.590,0
opencl-kernel,cpu-haswell-AMD Eng Sample: 100-000000020-02_30/20_N,CPU,w8,secp256k1,256,MODMUL,50000,0.017807802,2807758.053,0
opencl-e2e,cpu-haswell-AMD Eng Sample: 100-000000020-02_30/20_N,CPU,w8,secp256k1,256,MODMUL,50000,0.019411588,2575781.012,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,secp256k1,256,MODEXP,781,0.008961861,87147.078,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,secp256k1,256,MODEXP,781,0.008090096,96537.796,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,secp256k1,256,MODEXP,781,0.000545867,1430752.004,0
library,NVIDIA A40,gpu,cgbn,secp256k1,256,MODEXP,50000,0.018923519,2642214.696,0
opencl-kernel,NVIDIA A40,GPU,w8,secp256k1,256,MODEXP,50000,0.028947090,1727289.353,0
opencl-e2e,NVIDIA A40,GPU,w8,secp256k1,256,MODEXP,50000,0.029595791,1689429.399,0
opencl-kernel,NVIDIA A40,GPU,w16,secp256k1,256,MODEXP,50000,0.005468527,9143230.220,0
opencl-e2e,NVIDIA A40,GPU,w16,secp256k1,256,MODEXP,50000,0.006071392,8235344.381,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,secp256k1,256,MODEXP,50000,0.003530664,14161640.634,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,secp256k1,256,MODEXP,50000,0.004165649,12002930.403,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,secp256k1,256,MODEXP,50000,0.002036688,24549662.007,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,secp256k1,256,MODEXP,50000,0.002687715,18603162.952,0
opencl-kernel,NVIDIA A40,GPU,w32-il,secp256k1,256,MODEXP,50000,0.003542419,14114650.402,0
opencl-e2e,NVIDIA A40,GPU,w32-il,secp256k1,256,MODEXP,50000,0.004144422,12064409.143,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,secp256k1,256,MODEXP,50000,0.002048550,24407505.681,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,secp256k1,256,MODEXP,50000,0.002694117,18558957.075,0
opencl-kernel,cpu-haswell-AMD Eng Sample: 100-000000020-02_30/20_N,CPU,w8,secp256k1,256,MODEXP,50000,0.000000000,inf,0
opencl-e2e,cpu-haswell-AMD Eng Sample: 100-000000020-02_30/20_N,CPU,w8,secp256k1,256,MODEXP,50000,0.000000000,inf,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,secp256k1,256,EXPONENTIATION,781,0.002990500,261160.335,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,secp256k1,256,EXPONENTIATION,781,0.012976376,60186.294,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,secp256k1,256,EXPONENTIATION,781,0.003944451,197999.671,0
opencl-kernel,NVIDIA A40,GPU,w8,secp256k1,256,EXPONENTIATION,50000,0.045389568,1101574.695,0
opencl-e2e,NVIDIA A40,GPU,w8,secp256k1,256,EXPONENTIATION,50000,0.046032369,1086192.193,0
opencl-kernel,NVIDIA A40,GPU,w16,secp256k1,256,EXPONENTIATION,50000,0.012077972,4139767.735,0
opencl-e2e,NVIDIA A40,GPU,w16,secp256k1,256,EXPONENTIATION,50000,0.012698211,3937562.613,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,secp256k1,256,EXPONENTIATION,50000,0.000729766,68515137.191,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,secp256k1,256,EXPONENTIATION,50000,0.001287144,38845702.213,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,secp256k1,256,EXPONENTIATION,50000,0.000681634,73353141.895,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,secp256k1,256,EXPONENTIATION,50000,0.001264301,39547540.428,0
opencl-kernel,NVIDIA A40,GPU,w32-il,secp256k1,256,EXPONENTIATION,50000,0.000716502,69783490.741,0
opencl-e2e,NVIDIA A40,GPU,w32-il,secp256k1,256,EXPONENTIATION,50000,0.001268419,39419137.343,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,secp256k1,256,EXPONENTIATION,50000,0.000644314,77601902.210,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,secp256k1,256,EXPONENTIATION,50000,0.001218564,41031907.265,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,secp256k1,256,DIVIDE,6250,0.000331439,18857160.841,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,secp256k1,256,DIVIDE,6250,0.005409935,1155281.798,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,secp256k1,256,DIVIDE,6250,0.000073720,84780514.427,0
library,NVIDIA A40,gpu,cgbn,secp256k1,256,DIVIDE,50000,0.000026624,1878004807.692,0
opencl-kernel,NVIDIA A40,GPU,w8,secp256k1,256,DIVIDE,50000,0.000470363,106300967.829,0
opencl-e2e,NVIDIA A40,GPU,w8,secp256k1,256,DIVIDE,50000,0.001143342,43731451.370,0
opencl-kernel,NVIDIA A40,GPU,w16,secp256k1,256,DIVIDE,50000,0.000307434,162636407.204,0
opencl-e2e,NVIDIA A40,GPU,w16,secp256k1,256,DIVIDE,50000,0.000983840,50821279.262,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,secp256k1,256,DIVIDE,50000,0.000128344,389578915.593,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,secp256k1,256,DIVIDE,50000,0.000785412,63660834.074,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,secp256k1,256,DIVIDE,50000,0.000127633,391747901.784,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,secp256k1,256,DIVIDE,50000,0.000811311,61628616.329,0
opencl-kernel,NVIDIA A40,GPU,w32-il,secp256k1,256,DIVIDE,50000,0.000128754,388336199.177,0
opencl-e2e,NVIDIA A40,GPU,w32-il,secp256k1,256,DIVIDE,50000,0.000792645,63079947.738,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,secp256k1,256,DIVIDE,50000,0.000126831,394224660.753,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,secp256k1,256,DIVIDE,50000,0.000811202,61636965.354,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,secp256k1,256,ISQRT,1562,0.000175224,8914319.961,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,secp256k1,256,ISQRT,1562,0.000036559,42725378.400,0
opencl-kernel,NVIDIA A40,GPU,w8,secp256k1,256,ISQRT,50000,0.005942857,8413461.892,0
opencl-e2e,NVIDIA A40,GPU,w8,secp256k1,256,ISQRT,50000,0.006577494,7601679.738,0
opencl-kernel,NVIDIA A40,GPU,w16,secp256k1,256,ISQRT,50000,0.004526668,11045651.565,0
opencl-e2e,NVIDIA A40,GPU,w16,secp256k1,256,ISQRT,50000,0.005137909,9731584.973,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,secp256k1,256,ISQRT,50000,0.001032921,48406389.737,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,secp256k1,256,ISQRT,50000,0.001606280,31127823.601,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,secp256k1,256,ISQRT,50000,0.000904518,55278048.720,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,secp256k1,256,ISQRT,50000,0.001469952,34014714.870,0
opencl-kernel,NVIDIA A40,GPU,w32-il,secp256k1,256,ISQRT,50000,0.000911331,54864819.414,0
opencl-e2e,NVIDIA A40,GPU,w32-il,secp256k1,256,ISQRT,50000,0.001464922,34131510.514,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,secp256k1,256,ISQRT,50000,0.000911842,54834055.131,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,secp256k1,256,ISQRT,50000,0.001483909,33694791.867,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,secp256k1,256,MODMUL_R2,50000,0.005816537,8596180.233,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,secp256k1,256,MODMUL_R2,50000,0.000159693,313100859.047,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,secp256k1,256,MODMUL_R2,50000,0.000281996,177307420.630,0
library,NVIDIA A40,gpu,cgbn,secp256k1,256,MODMUL_R2,50000,0.000028672,1743861607.143,0
opencl-kernel,NVIDIA A40,GPU,w8,secp256k1,256,MODMUL_R2,50000,0.000157389,317684494.808,0
opencl-e2e,NVIDIA A40,GPU,w8,secp256k1,256,MODMUL_R2,50000,0.000701621,71263527.663,0
opencl-kernel,NVIDIA A40,GPU,w16,secp256k1,256,MODMUL_R2,50000,0.000049302,1014150349.465,0
opencl-e2e,NVIDIA A40,GPU,w16,secp256k1,256,MODMUL_R2,50000,0.000619587,80698853.560,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,secp256k1,256,MODMUL_R2,50000,0.000035989,1389309608.467,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,secp256k1,256,MODMUL_R2,50000,0.000593217,84286175.941,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,secp256k1,256,MODMUL_R2,50000,0.000035346,1414567786.473,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,secp256k1,256,MODMUL_R2,50000,0.000590642,84653650.522,0
opencl-kernel,NVIDIA A40,GPU,w32-il,secp256k1,256,MODMUL_R2,50000,0.000036268,1378642370.705,0
opencl-e2e,NVIDIA A40,GPU,w32-il,secp256k1,256,MODMUL_R2,50000,0.000542490,92167629.538,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,secp256k1,256,MODMUL_R2,50000,0.000034466,1450688802.421,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,secp256k1,256,MODMUL_R2,50000,0.000587075,85167991.345,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,rsa256(composite),256,ADD,50000,0.000978930,51076179.078,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,rsa256(composite),256,ADD,50000,0.007956642,6284058.240,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,rsa256(composite),256,ADD,50000,0.009994332,5002835.734,0
library,NVIDIA A40,gpu,cgbn,rsa256(composite),256,ADD,50000,0.000024576,2034505208.333,0
opencl-kernel,NVIDIA A40,GPU,w8,rsa256(composite),256,ADD,50000,0.000048081,1039903368.393,0
opencl-e2e,NVIDIA A40,GPU,w8,rsa256(composite),256,ADD,50000,0.000605249,82610649.027,0
opencl-kernel,NVIDIA A40,GPU,w16,rsa256(composite),256,ADD,50000,0.000029676,1684882349.987,0
opencl-e2e,NVIDIA A40,GPU,w16,rsa256(composite),256,ADD,50000,0.000601502,83125224.430,0
opencl-kernel,NVIDIA A40,GPU,w32,rsa256(composite),256,ADD,50000,0.000018305,2731472459.934,0
opencl-e2e,NVIDIA A40,GPU,w32,rsa256(composite),256,ADD,50000,0.000549604,90974562.006,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,rsa256(composite),256,ADD,50000,0.000018816,2657250603.841,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,rsa256(composite),256,ADD,50000,0.000583468,85694501.783,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,rsa256(composite),256,ADD,50000,0.000018966,2636243123.005,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,rsa256(composite),256,ADD,50000,0.000569913,87732749.833,0
opencl-kernel,NVIDIA A40,GPU,w32-il,rsa256(composite),256,ADD,50000,0.000012623,3960977659.731,0
opencl-e2e,NVIDIA A40,GPU,w32-il,rsa256(composite),256,ADD,50000,0.000561967,88973118.023,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,rsa256(composite),256,ADD,50000,0.000012213,4093876101.876,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,rsa256(composite),256,ADD,50000,0.000527492,94788204.594,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACT,50000,0.000742139,67372796.932,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACT,50000,0.000062990,793776760.553,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACT,50000,0.000086244,579749159.864,0
library,NVIDIA A40,gpu,cgbn,rsa256(composite),256,SUBTRACT,50000,0.000024576,2034505208.333,0
opencl-kernel,NVIDIA A40,GPU,w8,rsa256(composite),256,SUBTRACT,50000,0.000048102,1039460419.369,0
opencl-e2e,NVIDIA A40,GPU,w8,rsa256(composite),256,SUBTRACT,50000,0.000570974,87569613.704,0
opencl-kernel,NVIDIA A40,GPU,w16,rsa256(composite),256,SUBTRACT,50000,0.000028995,1724443233.868,0
opencl-e2e,NVIDIA A40,GPU,w16,rsa256(composite),256,SUBTRACT,50000,0.000576375,86749070.415,0
opencl-kernel,NVIDIA A40,GPU,w32,rsa256(composite),256,SUBTRACT,50000,0.000018995,2632236281.624,0
opencl-e2e,NVIDIA A40,GPU,w32,rsa256(composite),256,SUBTRACT,50000,0.000535096,93441169.601,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,rsa256(composite),256,SUBTRACT,50000,0.000018735,2668742416.861,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,rsa256(composite),256,SUBTRACT,50000,0.000539114,92744802.710,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,rsa256(composite),256,SUBTRACT,50000,0.000018215,2745019490.745,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,rsa256(composite),256,SUBTRACT,50000,0.000582346,85859644.420,0
opencl-kernel,NVIDIA A40,GPU,w32-il,rsa256(composite),256,SUBTRACT,50000,0.000011981,4173112413.525,0
opencl-e2e,NVIDIA A40,GPU,w32-il,rsa256(composite),256,SUBTRACT,50000,0.000551827,90608065.888,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,rsa256(composite),256,SUBTRACT,50000,0.000012213,4093876101.876,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,rsa256(composite),256,SUBTRACT,50000,0.000561057,89117559.218,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,rsa256(composite),256,ADDMOD,50000,0.002626529,19036530.987,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,rsa256(composite),256,ADDMOD,50000,0.000111221,449554032.305,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,rsa256(composite),256,ADDMOD,50000,0.008952575,5584985.447,0
library,NVIDIA A40,gpu,cgbn,rsa256(composite),256,ADDMOD,50000,0.000025600,1953125000.000,0
opencl-kernel,NVIDIA A40,GPU,w8,rsa256(composite),256,ADDMOD,50000,0.000064493,775276050.196,0
opencl-e2e,NVIDIA A40,GPU,w8,rsa256(composite),256,ADDMOD,50000,0.000669232,74712546.654,0
opencl-kernel,NVIDIA A40,GPU,w16,rsa256(composite),256,ADDMOD,50000,0.000036799,1358719692.253,0
opencl-e2e,NVIDIA A40,GPU,w16,rsa256(composite),256,ADDMOD,50000,0.000584571,85532854.744,0
opencl-kernel,NVIDIA A40,GPU,w32,rsa256(composite),256,ADDMOD,50000,0.000019847,2519219708.132,0
opencl-e2e,NVIDIA A40,GPU,w32,rsa256(composite),256,ADDMOD,50000,0.000579441,86290067.345,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,rsa256(composite),256,ADDMOD,50000,0.000014878,3360694284.820,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,rsa256(composite),256,ADDMOD,50000,0.000577917,86517565.838,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,rsa256(composite),256,ADDMOD,50000,0.000015249,3278801221.449,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,rsa256(composite),256,ADDMOD,50000,0.000544875,91764178.069,0
opencl-kernel,NVIDIA A40,GPU,w32-il,rsa256(composite),256,ADDMOD,50000,0.000012383,4037837785.800,0
opencl-e2e,NVIDIA A40,GPU,w32-il,rsa256(composite),256,ADDMOD,50000,0.000564411,88587881.828,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,rsa256(composite),256,ADDMOD,50000,0.000012283,4070596042.156,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,rsa256(composite),256,ADDMOD,50000,0.000560484,89208629.230,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACTMOD,50000,0.002781974,17972848.514,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACTMOD,50000,0.000092206,542266463.310,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACTMOD,50000,0.000282126,177226062.787,0
library,NVIDIA A40,gpu,cgbn,rsa256(composite),256,SUBTRACTMOD,50000,0.000025600,1953125000.000,0
opencl-kernel,NVIDIA A40,GPU,w8,rsa256(composite),256,SUBTRACTMOD,50000,0.000069803,716305419.613,0
opencl-e2e,NVIDIA A40,GPU,w8,rsa256(composite),256,SUBTRACTMOD,50000,0.000624426,80073576.823,0
opencl-kernel,NVIDIA A40,GPU,w16,rsa256(composite),256,SUBTRACTMOD,50000,0.000039885,1253609751.086,0
opencl-e2e,NVIDIA A40,GPU,w16,rsa256(composite),256,SUBTRACTMOD,50000,0.000591986,84461473.790,0
opencl-kernel,NVIDIA A40,GPU,w32,rsa256(composite),256,SUBTRACTMOD,50000,0.000021131,2366216721.759,0
opencl-e2e,NVIDIA A40,GPU,w32,rsa256(composite),256,SUBTRACTMOD,50000,0.000592205,84430259.406,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,rsa256(composite),256,SUBTRACTMOD,50000,0.000013997,3572233096.014,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,rsa256(composite),256,SUBTRACTMOD,50000,0.000583979,85619472.987,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,rsa256(composite),256,SUBTRACTMOD,50000,0.000015099,3311564964.224,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,rsa256(composite),256,SUBTRACTMOD,50000,0.000572758,87296934.929,0
opencl-kernel,NVIDIA A40,GPU,w32-il,rsa256(composite),256,SUBTRACTMOD,50000,0.000012884,3880807517.710,0
opencl-e2e,NVIDIA A40,GPU,w32-il,rsa256(composite),256,SUBTRACTMOD,50000,0.000557389,89703941.239,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,rsa256(composite),256,SUBTRACTMOD,50000,0.000011973,4176033851.898,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,rsa256(composite),256,SUBTRACTMOD,50000,0.000538533,92844886.164,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001198416,41721751.740,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.007914152,6317796.275,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000120819,413843512.580,0
opencl-kernel,NVIDIA A40,GPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001505329,33215324.298,0
opencl-e2e,NVIDIA A40,GPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.002188066,22851230.049,0
opencl-kernel,NVIDIA A40,GPU,w16,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000453592,110231153.435,0
opencl-e2e,NVIDIA A40,GPU,w16,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001119056,44680525.010,0
opencl-kernel,NVIDIA A40,GPU,w32,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000114799,435545586.708,0
opencl-e2e,NVIDIA A40,GPU,w32,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000792104,63123038.630,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000030218,1654659779.326,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000712884,70137671.679,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000030028,1665129061.473,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000705319,70889863.481,0
opencl-kernel,NVIDIA A40,GPU,w32-il,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000015539,3217686017.381,0
opencl-e2e,NVIDIA A40,GPU,w32-il,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000631450,79182865.251,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000014838,3369764699.975,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000677356,73816373.301,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001203045,41561196.507,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000066857,747866482.789,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000124057,403041111.069,0
library,NVIDIA A40,gpu,cgbn,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000024576,2034505208.333,0
opencl-kernel,NVIDIA A40,GPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000279211,179076354.903,0
opencl-e2e,NVIDIA A40,GPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000944774,52922728.236,0
opencl-kernel,NVIDIA A40,GPU,w16,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000096864,516187287.394,0
opencl-e2e,NVIDIA A40,GPU,w16,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000760715,65727678.986,0
opencl-kernel,NVIDIA A40,GPU,w32,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000034234,1460515552.654,0
opencl-e2e,NVIDIA A40,GPU,w32,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000699027,71527950.172,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000034947,1430740091.675,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000682686,73240164.359,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000033934,1473462816.994,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000703766,71046341.028,0
opencl-kernel,NVIDIA A40,GPU,w32-il,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000026119,1914319529.328,0
opencl-e2e,NVIDIA A40,GPU,w32-il,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000682666,73242262.625,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000024757,2019602422.601,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000685292,72961567.325,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.009764486,5120597.238,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.007937836,6298945.594,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.009016845,5545176.602,0
library,NVIDIA A40,gpu,cgbn,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000024576,2034505208.333,0
opencl-kernel,NVIDIA A40,GPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000139825,357589726.648,0
opencl-e2e,NVIDIA A40,GPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000720779,69369351.014,0
opencl-kernel,NVIDIA A40,GPU,w16,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000042470,1177296855.401,0
opencl-e2e,NVIDIA A40,GPU,w16,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000600912,83206903.649,0
opencl-kernel,NVIDIA A40,GPU,w32,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000016891,2960086629.542,0
opencl-e2e,NVIDIA A40,GPU,w32,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000499779,100044147.780,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000023616,2117249327.602,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000549885,90928183.675,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000023806,2100351754.626,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000538573,92837982.459,0
opencl-kernel,NVIDIA A40,GPU,w32-il,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000024927,2005869277.041,0
opencl-e2e,NVIDIA A40,GPU,w32-il,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000576946,86663230.297,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000023915,2090700229.760,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000568259,87988113.340,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,rsa256(composite),256,COMPARE,50000,0.000397034,125933802.473,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,rsa256(composite),256,COMPARE,50000,0.000086784,576140659.341,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,rsa256(composite),256,COMPARE,50000,0.001043522,47914667.466,0
library,NVIDIA A40,gpu,cgbn,rsa256(composite),256,COMPARE,50000,0.000023552,2122961956.522,0
opencl-kernel,NVIDIA A40,GPU,w8,rsa256(composite),256,COMPARE,50000,0.000046900,1066108487.231,0
opencl-e2e,NVIDIA A40,GPU,w8,rsa256(composite),256,COMPARE,50000,0.000535877,93304920.107,0
opencl-kernel,NVIDIA A40,GPU,w16,rsa256(composite),256,COMPARE,50000,0.000028825,1734583412.491,0
opencl-e2e,NVIDIA A40,GPU,w16,rsa256(composite),256,COMPARE,50000,0.000574601,87016921.648,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,rsa256(composite),256,COMPARE,50000,0.000014387,3475342516.831,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,rsa256(composite),256,COMPARE,50000,0.000531259,94116053.889,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,rsa256(composite),256,COMPARE,50000,0.000015770,3170559924.408,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,rsa256(composite),256,COMPARE,50000,0.000563460,88737380.312,0
opencl-kernel,NVIDIA A40,GPU,w32-il,rsa256(composite),256,COMPARE,50000,0.000012003,4165665052.762,0
opencl-e2e,NVIDIA A40,GPU,w32-il,rsa256(composite),256,COMPARE,50000,0.000581563,85975141.605,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,rsa256(composite),256,COMPARE,50000,0.000012093,4134546877.166,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,rsa256(composite),256,COMPARE,50000,0.000536680,93165348.734,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,rsa256(composite),256,REDUCE,6250,0.000248232,25178066.835,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,rsa256(composite),256,REDUCE,6250,0.000038173,163728076.510,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,rsa256(composite),256,REDUCE,6250,0.000061477,101664693.228,0
library,NVIDIA A40,gpu,cgbn,rsa256(composite),256,REDUCE,50000,0.000025600,1953125000.000,0
opencl-kernel,NVIDIA A40,GPU,w8,rsa256(composite),256,REDUCE,50000,0.000247560,201971638.922,0
opencl-e2e,NVIDIA A40,GPU,w8,rsa256(composite),256,REDUCE,50000,0.000826920,60465313.285,0
opencl-kernel,NVIDIA A40,GPU,w16,rsa256(composite),256,REDUCE,50000,0.000158531,315396403.499,0
opencl-e2e,NVIDIA A40,GPU,w16,rsa256(composite),256,REDUCE,50000,0.000659432,75822799.198,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,rsa256(composite),256,REDUCE,50000,0.000058671,852216632.538,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,rsa256(composite),256,REDUCE,50000,0.000627672,79659401.921,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,rsa256(composite),256,REDUCE,50000,0.000057619,867768332.579,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,rsa256(composite),256,REDUCE,50000,0.000616340,81124059.294,0
opencl-kernel,NVIDIA A40,GPU,w32-il,rsa256(composite),256,REDUCE,50000,0.000058250,858375428.891,0
opencl-e2e,NVIDIA A40,GPU,w32-il,rsa256(composite),256,REDUCE,50000,0.000608747,82135942.733,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,rsa256(composite),256,REDUCE,50000,0.000057409,870949858.862,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,rsa256(composite),256,REDUCE,50000,0.000616662,81081667.814,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,rsa256(composite),256,MODMUL,3125,0.000352439,8866770.781,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,rsa256(composite),256,MODMUL,3125,0.000041850,74671604.059,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,rsa256(composite),256,MODMUL,3125,0.000077207,40475792.521,0
library,NVIDIA A40,gpu,cgbn,rsa256(composite),256,MODMUL,50000,0.000061440,813802083.333,0
opencl-kernel,NVIDIA A40,GPU,w8,rsa256(composite),256,MODMUL,50000,0.000641649,77924264.948,0
opencl-e2e,NVIDIA A40,GPU,w8,rsa256(composite),256,MODMUL,50000,0.001190601,41995600.118,0
opencl-kernel,NVIDIA A40,GPU,w16,rsa256(composite),256,MODMUL,50000,0.000376425,132828672.161,0
opencl-e2e,NVIDIA A40,GPU,w16,rsa256(composite),256,MODMUL,50000,0.000934405,53509973.189,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,rsa256(composite),256,MODMUL,50000,0.000170724,292869563.699,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,rsa256(composite),256,MODMUL,50000,0.000736398,67898094.476,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,rsa256(composite),256,MODMUL,50000,0.000141148,354239299.005,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,rsa256(composite),256,MODMUL,50000,0.000698256,71606943.676,0
opencl-kernel,NVIDIA A40,GPU,w32-il,rsa256(composite),256,MODMUL,50000,0.000165454,302198593.904,0
opencl-e2e,NVIDIA A40,GPU,w32-il,rsa256(composite),256,MODMUL,50000,0.000717822,69655106.415,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,rsa256(composite),256,MODMUL,50000,0.000141980,352161962.611,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,rsa256(composite),256,MODMUL,50000,0.000660224,75731885.853,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,rsa256(composite),256,MODEXP,781,0.008616195,90643.259,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,rsa256(composite),256,MODEXP,781,0.000243893,3202225.328,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,rsa256(composite),256,MODEXP,781,0.000445317,1753808.628,0
library,NVIDIA A40,gpu,cgbn,rsa256(composite),256,MODEXP,50000,0.018744320,2667474.734,0
opencl-kernel,NVIDIA A40,GPU,w8,rsa256(composite),256,MODEXP,50000,0.028940377,1727690.013,0
opencl-e2e,NVIDIA A40,GPU,w8,rsa256(composite),256,MODEXP,50000,0.029589880,1689766.892,0
opencl-kernel,NVIDIA A40,GPU,w16,rsa256(composite),256,MODEXP,50000,0.005465271,9148677.247,0
opencl-e2e,NVIDIA A40,GPU,w16,rsa256(composite),256,MODEXP,50000,0.006068968,8238632.703,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,rsa256(composite),256,MODEXP,50000,0.003527669,14173664.438,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,rsa256(composite),256,MODEXP,50000,0.004165961,12002031.491,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,rsa256(composite),256,MODEXP,50000,0.002035205,24567546.659,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,rsa256(composite),256,MODEXP,50000,0.002648690,18877257.843,0
opencl-kernel,NVIDIA A40,GPU,w32-il,rsa256(composite),256,MODEXP,50000,0.003537279,14135160.200,0
opencl-e2e,NVIDIA A40,GPU,w32-il,rsa256(composite),256,MODEXP,50000,0.004146135,12059422.819,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,rsa256(composite),256,MODEXP,50000,0.002041217,24495191.583,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,rsa256(composite),256,MODEXP,50000,0.002685631,18617600.741,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,rsa256(composite),256,EXPONENTIATION,781,0.003020948,258528.135,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,rsa256(composite),256,EXPONENTIATION,781,0.000137379,5684986.540,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,rsa256(composite),256,EXPONENTIATION,781,0.001534525,508952.180,0
opencl-kernel,NVIDIA A40,GPU,w8,rsa256(composite),256,EXPONENTIATION,50000,0.045318863,1103293.337,0
opencl-e2e,NVIDIA A40,GPU,w8,rsa256(composite),256,EXPONENTIATION,50000,0.046084588,1084961.410,0
opencl-kernel,NVIDIA A40,GPU,w16,rsa256(composite),256,EXPONENTIATION,50000,0.012097379,4133126.622,0
opencl-e2e,NVIDIA A40,GPU,w16,rsa256(composite),256,EXPONENTIATION,50000,0.014040490,3561129.358,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,rsa256(composite),256,EXPONENTIATION,50000,0.000739184,67642137.703,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,rsa256(composite),256,EXPONENTIATION,50000,0.001293727,38648024.194,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,rsa256(composite),256,EXPONENTIATION,50000,0.000678759,73663840.881,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,rsa256(composite),256,EXPONENTIATION,50000,0.001236368,40441034.696,0
opencl-kernel,NVIDIA A40,GPU,w32-il,rsa256(composite),256,EXPONENTIATION,50000,0.000716461,69787482.029,0
opencl-e2e,NVIDIA A40,GPU,w32-il,rsa256(composite),256,EXPONENTIATION,50000,0.001268809,39407042.858,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,rsa256(composite),256,EXPONENTIATION,50000,0.000645956,77404650.009,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,rsa256(composite),256,EXPONENTIATION,50000,0.001212033,41253018.803,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,rsa256(composite),256,DIVIDE,6250,0.000329235,18983421.826,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,rsa256(composite),256,DIVIDE,6250,0.000045175,138351676.081,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,rsa256(composite),256,DIVIDE,6250,0.000079401,78714535.047,0
library,NVIDIA A40,gpu,cgbn,rsa256(composite),256,DIVIDE,50000,0.000026624,1878004807.692,0
opencl-kernel,NVIDIA A40,GPU,w8,rsa256(composite),256,DIVIDE,50000,0.000471815,105973634.842,0
opencl-e2e,NVIDIA A40,GPU,w8,rsa256(composite),256,DIVIDE,50000,0.001141558,43799774.014,0
opencl-kernel,NVIDIA A40,GPU,w16,rsa256(composite),256,DIVIDE,50000,0.000314097,159186533.792,0
opencl-e2e,NVIDIA A40,GPU,w16,rsa256(composite),256,DIVIDE,50000,0.000983949,50815651.198,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,rsa256(composite),256,DIVIDE,50000,0.000133624,374183437.182,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,rsa256(composite),256,DIVIDE,50000,0.000781594,63971844.667,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,rsa256(composite),256,DIVIDE,50000,0.000132461,377469371.225,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,rsa256(composite),256,DIVIDE,50000,0.000810578,61684342.887,0
opencl-kernel,NVIDIA A40,GPU,w32-il,rsa256(composite),256,DIVIDE,50000,0.000133272,375171846.261,0
opencl-e2e,NVIDIA A40,GPU,w32-il,rsa256(composite),256,DIVIDE,50000,0.000792555,63087137.823,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,rsa256(composite),256,DIVIDE,50000,0.000130639,382732893.714,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,rsa256(composite),256,DIVIDE,50000,0.000814196,61410295.618,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,rsa256(composite),256,ISQRT,1562,0.000149995,10413674.306,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,rsa256(composite),256,ISQRT,1562,0.000040807,38277905.995,0
opencl-kernel,NVIDIA A40,GPU,w8,rsa256(composite),256,ISQRT,50000,0.005944460,8411193.368,0
opencl-e2e,NVIDIA A40,GPU,w8,rsa256(composite),256,ISQRT,50000,0.006584827,7593214.079,0
opencl-kernel,NVIDIA A40,GPU,w16,rsa256(composite),256,ISQRT,50000,0.004528071,11042230.170,0
opencl-e2e,NVIDIA A40,GPU,w16,rsa256(composite),256,ISQRT,50000,0.005122760,9760362.962,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,rsa256(composite),256,ISQRT,50000,0.001046248,47789826.598,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,rsa256(composite),256,ISQRT,50000,0.001606421,31125098.601,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,rsa256(composite),256,ISQRT,50000,0.000910910,54890174.034,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,rsa256(composite),256,ISQRT,50000,0.001424596,35097676.851,0
opencl-kernel,NVIDIA A40,GPU,w32-il,rsa256(composite),256,ISQRT,50000,0.000913914,54709725.387,0
opencl-e2e,NVIDIA A40,GPU,w32-il,rsa256(composite),256,ISQRT,50000,0.001461156,34219487.744,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,rsa256(composite),256,ISQRT,50000,0.000917462,54498187.728,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,rsa256(composite),256,ISQRT,50000,0.001489318,33572414.224,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,rsa256(composite),256,MODMUL_R2,50000,0.005587403,8948701.643,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,rsa256(composite),256,MODMUL_R2,50000,0.000122854,406988630.385,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,rsa256(composite),256,MODMUL_R2,50000,0.000404199,123701549.982,0
library,NVIDIA A40,gpu,cgbn,rsa256(composite),256,MODMUL_R2,50000,0.000028672,1743861607.143,0
opencl-kernel,NVIDIA A40,GPU,w8,rsa256(composite),256,MODMUL_R2,50000,0.000157448,317564230.239,0
opencl-e2e,NVIDIA A40,GPU,w8,rsa256(composite),256,MODMUL_R2,50000,0.000716612,69772789.087,0
opencl-kernel,NVIDIA A40,GPU,w16,rsa256(composite),256,MODMUL_R2,50000,0.000050858,983136009.376,0
opencl-e2e,NVIDIA A40,GPU,w16,rsa256(composite),256,MODMUL_R2,50000,0.000599247,83437991.909,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,rsa256(composite),256,MODMUL_R2,50000,0.000035136,1423041620.060,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,rsa256(composite),256,MODMUL_R2,50000,0.000598517,83539781.500,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,rsa256(composite),256,MODMUL_R2,50000,0.000034957,1430320799.254,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,rsa256(composite),256,MODMUL_R2,50000,0.000533914,93648006.754,0
opencl-kernel,NVIDIA A40,GPU,w32-il,rsa256(composite),256,MODMUL_R2,50000,0.000036680,1363135488.130,0
opencl-e2e,NVIDIA A40,GPU,w32-il,rsa256(composite),256,MODMUL_R2,50000,0.000559422,89377935.156,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,rsa256(composite),256,MODMUL_R2,50000,0.000034876,1433643751.335,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,rsa256(composite),256,MODMUL_R2,50000,0.000586935,85188262.456,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,brainpoolP512r1,512,ADD,25000,0.000521811,47910092.274,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,brainpoolP512r1,512,ADD,25000,0.000051057,489649148.152,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,brainpoolP512r1,512,ADD,25000,0.000079521,314382451.250,0
library,NVIDIA A40,gpu,cgbn,brainpoolP512r1,512,ADD,50000,0.000034816,1436121323.529,0
opencl-kernel,NVIDIA A40,GPU,w8,brainpoolP512r1,512,ADD,50000,0.000101563,492307258.922,0
opencl-e2e,NVIDIA A40,GPU,w8,brainpoolP512r1,512,ADD,50000,0.000953321,52448211.985,0
opencl-kernel,NVIDIA A40,GPU,w16,brainpoolP512r1,512,ADD,50000,0.000059113,845838971.515,0
opencl-e2e,NVIDIA A40,GPU,w16,brainpoolP512r1,512,ADD,50000,0.001002285,49886025.488,0
opencl-kernel,NVIDIA A40,GPU,w32,brainpoolP512r1,512,ADD,50000,0.000041088,1216897665.352,0
opencl-e2e,NVIDIA A40,GPU,w32,brainpoolP512r1,512,ADD,50000,0.000946278,52838608.820,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,brainpoolP512r1,512,ADD,50000,0.000039014,1281590107.660,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,brainpoolP512r1,512,ADD,50000,0.000909216,54992446.912,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,brainpoolP512r1,512,ADD,50000,0.000038965,1283213614.418,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,brainpoolP512r1,512,ADD,50000,0.001025768,48743965.644,0
opencl-kernel,NVIDIA A40,GPU,w32-il,brainpoolP512r1,512,ADD,50000,0.000031119,1606724462.800,0
opencl-e2e,NVIDIA A40,GPU,w32-il,brainpoolP512r1,512,ADD,50000,0.000834934,59884943.291,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,brainpoolP512r1,512,ADD,50000,0.000034445,1451590947.681,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,brainpoolP512r1,512,ADD,50000,0.000891413,56090748.025,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACT,25000,0.000410018,60972860.213,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACT,25000,0.000056147,445262587.291,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACT,25000,0.009998490,2500377.531,0
library,NVIDIA A40,gpu,cgbn,brainpoolP512r1,512,SUBTRACT,50000,0.000035840,1395089285.714,0
opencl-kernel,NVIDIA A40,GPU,w8,brainpoolP512r1,512,SUBTRACT,50000,0.000101512,492551158.737,0
opencl-e2e,NVIDIA A40,GPU,w8,brainpoolP512r1,512,SUBTRACT,50000,0.000961787,51986558.873,0
opencl-kernel,NVIDIA A40,GPU,w16,brainpoolP512r1,512,SUBTRACT,50000,0.000058852,849586833.776,0
opencl-e2e,NVIDIA A40,GPU,w16,brainpoolP512r1,512,SUBTRACT,50000,0.000984089,50808437.554,0
opencl-kernel,NVIDIA A40,GPU,w32,brainpoolP512r1,512,SUBTRACT,50000,0.000041299,1210695724.337,0
opencl-e2e,NVIDIA A40,GPU,w32,brainpoolP512r1,512,SUBTRACT,50000,0.000913215,54751627.085,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,brainpoolP512r1,512,SUBTRACT,50000,0.000039013,1281620701.838,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,brainpoolP512r1,512,SUBTRACT,50000,0.000900229,55541396.255,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,brainpoolP512r1,512,SUBTRACT,50000,0.000038634,1294195000.362,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,brainpoolP512r1,512,SUBTRACT,50000,0.001017803,49125404.398,0
opencl-kernel,NVIDIA A40,GPU,w32-il,brainpoolP512r1,512,SUBTRACT,50000,0.000030999,1612951515.698,0
opencl-e2e,NVIDIA A40,GPU,w32-il,brainpoolP512r1,512,SUBTRACT,50000,0.000833453,59991408.327,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,brainpoolP512r1,512,SUBTRACT,50000,0.000033744,1481758975.491,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,brainpoolP512r1,512,SUBTRACT,50000,0.000876905,57018693.459,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,brainpoolP512r1,512,ADDMOD,25000,0.001397264,17892106.717,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,brainpoolP512r1,512,ADDMOD,25000,0.000080183,311786210.741,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,brainpoolP512r1,512,ADDMOD,25000,0.007983954,3131280.687,0
library,NVIDIA A40,gpu,cgbn,brainpoolP512r1,512,ADDMOD,50000,0.000034816,1436121323.529,0
opencl-kernel,NVIDIA A40,GPU,w8,brainpoolP512r1,512,ADDMOD,50000,0.000137200,364432422.599,0
opencl-e2e,NVIDIA A40,GPU,w8,brainpoolP512r1,512,ADDMOD,50000,0.001029696,48558006.129,0
opencl-kernel,NVIDIA A40,GPU,w16,brainpoolP512r1,512,ADDMOD,50000,0.000075825,659416959.811,0
opencl-e2e,NVIDIA A40,GPU,w16,brainpoolP512r1,512,ADDMOD,50000,0.001099219,45486857.068,0
opencl-kernel,NVIDIA A40,GPU,w32,brainpoolP512r1,512,ADDMOD,50000,0.000046989,1064079977.802,0
opencl-e2e,NVIDIA A40,GPU,w32,brainpoolP512r1,512,ADDMOD,50000,0.000899679,55575375.742,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,brainpoolP512r1,512,ADDMOD,50000,0.000048181,1037752565.044,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,brainpoolP512r1,512,ADDMOD,50000,0.000899908,55561226.946,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,brainpoolP512r1,512,ADDMOD,50000,0.000049494,1010219238.296,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,brainpoolP512r1,512,ADDMOD,50000,0.000995641,50218922.213,0
opencl-kernel,NVIDIA A40,GPU,w32-il,brainpoolP512r1,512,ADDMOD,50000,0.000035688,1401020125.261,0
opencl-e2e,NVIDIA A40,GPU,w32-il,brainpoolP512r1,512,ADDMOD,50000,0.000832271,60076597.840,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,brainpoolP512r1,512,ADDMOD,50000,0.000031851,1569797988.304,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,brainpoolP512r1,512,ADDMOD,50000,0.000895401,55840931.780,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACTMOD,25000,0.001459893,17124544.575,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACTMOD,25000,0.009999171,2500207.292,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACTMOD,25000,0.007907027,3161744.459,0
library,NVIDIA A40,gpu,cgbn,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000034816,1436121323.529,0
opencl-kernel,NVIDIA A40,GPU,w8,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000156477,319535587.087,0
opencl-e2e,NVIDIA A40,GPU,w8,brainpoolP512r1,512,SUBTRACTMOD,50000,0.001027172,48677363.029,0
opencl-kernel,NVIDIA A40,GPU,w16,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000084249,593476721.717,0
opencl-e2e,NVIDIA A40,GPU,w16,brainpoolP512r1,512,SUBTRACTMOD,50000,0.001029085,48586834.056,0
opencl-kernel,NVIDIA A40,GPU,w32,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000050055,998903940.758,0
opencl-e2e,NVIDIA A40,GPU,w32,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000905610,55211423.388,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000049053,1019310636.036,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000913385,54741410.771,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000049043,1019504200.532,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000923574,54137515.328,0
opencl-kernel,NVIDIA A40,GPU,w32-il,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000036208,1380911857.606,0
opencl-e2e,NVIDIA A40,GPU,w32-il,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000870003,57471042.082,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000032401,1543175947.111,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000884460,56531650.805,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.001455254,17179132.543,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.008284206,3017790.878,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.000097875,255426580.267,0
opencl-kernel,NVIDIA A40,GPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.008684524,5757367.712,0
opencl-e2e,NVIDIA A40,GPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.010078472,4961069.453,0
opencl-kernel,NVIDIA A40,GPU,w16,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.002270132,22025148.018,0
opencl-e2e,NVIDIA A40,GPU,w16,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.003453871,14476509.573,0
opencl-kernel,NVIDIA A40,GPU,w32,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.003222099,15517836.854,0
opencl-e2e,NVIDIA A40,GPU,w32,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.004377804,11421251.998,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000124787,400682826.202,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.001242571,40239162.944,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000127072,393476332.241,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.001260964,39652196.314,0
opencl-kernel,NVIDIA A40,GPU,w32-il,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000093347,535633598.388,0
opencl-e2e,NVIDIA A40,GPU,w32-il,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.001236278,40443989.836,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000073189,683163555.850,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.001171065,42696194.194,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.001460545,17116900.898,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000073099,342003919.020,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000091393,273545281.865,0
library,NVIDIA A40,gpu,cgbn,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000034816,1436121323.529,0
opencl-kernel,NVIDIA A40,GPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.001917962,26069338.389,0
opencl-e2e,NVIDIA A40,GPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.003007352,16625920.127,0
opencl-kernel,NVIDIA A40,GPU,w16,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000528453,94615808.340,0
opencl-e2e,NVIDIA A40,GPU,w16,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.001712383,29199079.110,0
opencl-kernel,NVIDIA A40,GPU,w32,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000284730,175604677.376,0
opencl-e2e,NVIDIA A40,GPU,w32,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.001385181,36096360.249,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000173760,287753807.893,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.001283607,38952747.155,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000173629,287969893.742,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.001308585,38209222.807,0
opencl-kernel,NVIDIA A40,GPU,w32-il,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000089331,559718626.327,0
opencl-e2e,NVIDIA A40,GPU,w32-il,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.001301381,38420729.568,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000107945,463199095.811,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.001118053,44720609.080,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.011398079,2193352.079,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000327181,76410288.408,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000139254,179527869.291,0
library,NVIDIA A40,gpu,cgbn,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000035840,1395089285.714,0
opencl-kernel,NVIDIA A40,GPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000465172,107487043.796,0
opencl-e2e,NVIDIA A40,GPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.001368400,36539031.435,0
opencl-kernel,NVIDIA A40,GPU,w16,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000150817,331526631.633,0
opencl-e2e,NVIDIA A40,GPU,w16,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.001086604,46014929.907,0
opencl-kernel,NVIDIA A40,GPU,w32,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000045727,1093445715.799,0
opencl-e2e,NVIDIA A40,GPU,w32,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000929285,53804816.138,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000068290,732169915.173,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000939274,53232591.335,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000057268,873088601.584,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000944754,52923823.813,0
opencl-kernel,NVIDIA A40,GPU,w32-il,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000065484,763544311.863,0
opencl-e2e,NVIDIA A40,GPU,w32-il,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000997985,50100964.372,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000051279,975065223.393,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000898967,55619421.129,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,brainpoolP512r1,512,COMPARE,25000,0.000199820,125112654.564,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,brainpoolP512r1,512,COMPARE,25000,0.000050205,497960294.581,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,brainpoolP512r1,512,COMPARE,25000,0.007782161,3212475.204,0
library,NVIDIA A40,gpu,cgbn,brainpoolP512r1,512,COMPARE,50000,0.000034816,1436121323.529,0
opencl-kernel,NVIDIA A40,GPU,w8,brainpoolP512r1,512,COMPARE,50000,0.000101422,492989882.554,0
opencl-e2e,NVIDIA A40,GPU,w8,brainpoolP512r1,512,COMPARE,50000,0.000995451,50228506.900,0
opencl-kernel,NVIDIA A40,GPU,w16,brainpoolP512r1,512,COMPARE,50000,0.000058722,851473247.478,0
opencl-e2e,NVIDIA A40,GPU,w16,brainpoolP512r1,512,COMPARE,50000,0.001060514,47146959.300,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,brainpoolP512r1,512,COMPARE,50000,0.000034776,1437790337.440,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,brainpoolP512r1,512,COMPARE,50000,0.000892475,56023962.734,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,brainpoolP512r1,512,COMPARE,50000,0.000034676,1441922251.766,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,brainpoolP512r1,512,COMPARE,50000,0.000898356,55657246.408,0
opencl-kernel,NVIDIA A40,GPU,w32-il,brainpoolP512r1,512,COMPARE,50000,0.000014187,3524393829.187,0
opencl-e2e,NVIDIA A40,GPU,w32-il,brainpoolP512r1,512,COMPARE,50000,0.001015148,49253895.801,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,brainpoolP512r1,512,COMPARE,50000,0.000013635,3666900566.901,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,brainpoolP512r1,512,COMPARE,50000,0.000885872,56441551.812,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,brainpoolP512r1,512,REDUCE,3125,0.000126169,24768353.842,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,brainpoolP512r1,512,REDUCE,3125,0.000430447,7259894.113,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,brainpoolP512r1,512,REDUCE,3125,0.000064933,48126722.221,0
library,NVIDIA A40,gpu,cgbn,brainpoolP512r1,512,REDUCE,50000,0.000034816,1436121323.529,0
opencl-kernel,NVIDIA A40,GPU,w8,brainpoolP512r1,512,REDUCE,50000,0.000711951,70229512.383,0
opencl-e2e,NVIDIA A40,GPU,w8,brainpoolP512r1,512,REDUCE,50000,0.001613336,30991690.383,0
opencl-kernel,NVIDIA A40,GPU,w16,brainpoolP512r1,512,REDUCE,50000,0.000563409,88745447.933,0
opencl-e2e,NVIDIA A40,GPU,w16,brainpoolP512r1,512,REDUCE,50000,0.001554972,32154909.654,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,brainpoolP512r1,512,REDUCE,50000,0.000156928,318617751.929,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,brainpoolP512r1,512,REDUCE,50000,0.001049602,47637086.003,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,brainpoolP512r1,512,REDUCE,50000,0.000157230,318006273.989,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,brainpoolP512r1,512,REDUCE,50000,0.001045646,47817323.474,0
opencl-kernel,NVIDIA A40,GPU,w32-il,brainpoolP512r1,512,REDUCE,50000,0.000154141,324377620.417,0
opencl-e2e,NVIDIA A40,GPU,w32-il,brainpoolP512r1,512,REDUCE,50000,0.001142991,43744884.977,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,brainpoolP512r1,512,REDUCE,50000,0.000170744,292836017.127,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,brainpoolP512r1,512,REDUCE,50000,0.001063360,47020768.734,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL,1562,0.000383228,4075901.920,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL,1562,0.000040858,38229917.920,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL,1562,0.000073761,21176574.862,0
library,NVIDIA A40,gpu,cgbn,brainpoolP512r1,512,MODMUL,50000,0.000192512,259724069.149,0
opencl-kernel,NVIDIA A40,GPU,w8,brainpoolP512r1,512,MODMUL,50000,0.002028292,24651282.404,0
opencl-e2e,NVIDIA A40,GPU,w8,brainpoolP512r1,512,MODMUL,50000,0.002934704,17037495.692,0
opencl-kernel,NVIDIA A40,GPU,w16,brainpoolP512r1,512,MODMUL,50000,0.001423294,35129783.216,0
opencl-e2e,NVIDIA A40,GPU,w16,brainpoolP512r1,512,MODMUL,50000,0.002412262,20727430.219,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,brainpoolP512r1,512,MODMUL,50000,0.000570344,87666420.423,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,brainpoolP512r1,512,MODMUL,50000,0.001488787,33584385.073,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,brainpoolP512r1,512,MODMUL,50000,0.000459984,108699446.046,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,brainpoolP512r1,512,MODMUL,50000,0.001367878,36552962.932,0
opencl-kernel,NVIDIA A40,GPU,w32-il,brainpoolP512r1,512,MODMUL,50000,0.000574922,86968290.573,0
opencl-e2e,NVIDIA A40,GPU,w32-il,brainpoolP512r1,512,MODMUL,50000,0.001556917,32114747.922,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,brainpoolP512r1,512,MODMUL,50000,0.000555405,90024332.872,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,brainpoolP512r1,512,MODMUL,50000,0.001449293,34499591.432,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,brainpoolP512r1,512,MODEXP,390,0.022816706,17092.739,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,brainpoolP512r1,512,MODEXP,390,0.000639825,609541.550,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,brainpoolP512r1,512,MODEXP,390,0.000663680,587632.572,0
library,NVIDIA A40,gpu,cgbn,brainpoolP512r1,512,MODEXP,50000,0.046244767,1081203.415,0
opencl-kernel,NVIDIA A40,GPU,w8,brainpoolP512r1,512,MODEXP,50000,0.421547231,118610.671,0
opencl-e2e,NVIDIA A40,GPU,w8,brainpoolP512r1,512,MODEXP,50000,0.422761106,118270.104,0
opencl-kernel,NVIDIA A40,GPU,w16,brainpoolP512r1,512,MODEXP,50000,0.038392247,1302346.289,0
opencl-e2e,NVIDIA A40,GPU,w16,brainpoolP512r1,512,MODEXP,50000,0.039524646,1265033.453,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,brainpoolP512r1,512,MODEXP,50000,0.031361004,1594336.695,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,brainpoolP512r1,512,MODEXP,50000,0.032356905,1545265.231,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,brainpoolP512r1,512,MODEXP,50000,0.019579453,2553697.440,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,brainpoolP512r1,512,MODEXP,50000,0.020524529,2436109.519,0
opencl-kernel,NVIDIA A40,GPU,w32-il,brainpoolP512r1,512,MODEXP,50000,0.032876932,1520823.169,0
opencl-e2e,NVIDIA A40,GPU,w32-il,brainpoolP512r1,512,MODEXP,50000,0.033868354,1476304.407,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,brainpoolP512r1,512,MODEXP,50000,0.014576134,3430264.779,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,brainpoolP512r1,512,MODEXP,50000,0.015620408,3200940.790,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,brainpoolP512r1,512,EXPONENTIATION,390,0.005219072,74725.927,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,brainpoolP512r1,512,EXPONENTIATION,390,0.009998721,39004.988,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,brainpoolP512r1,512,EXPONENTIATION,390,0.013983189,27890.633,0
opencl-kernel,NVIDIA A40,GPU,w8,brainpoolP512r1,512,EXPONENTIATION,50000,0.359463764,139096.079,0
opencl-e2e,NVIDIA A40,GPU,w8,brainpoolP512r1,512,EXPONENTIATION,50000,0.360707256,138616.563,0
opencl-kernel,NVIDIA A40,GPU,w16,brainpoolP512r1,512,EXPONENTIATION,50000,0.095496615,523578.766,0
opencl-e2e,NVIDIA A40,GPU,w16,brainpoolP512r1,512,EXPONENTIATION,50000,0.096725219,516928.269,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,brainpoolP512r1,512,EXPONENTIATION,50000,0.033094938,1510805.070,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,brainpoolP512r1,512,EXPONENTIATION,50000,0.034241394,1460220.918,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,brainpoolP512r1,512,EXPONENTIATION,50000,0.038447828,1300463.578,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,brainpoolP512r1,512,EXPONENTIATION,50000,0.039417319,1268477.950,0
opencl-kernel,NVIDIA A40,GPU,w32-il,brainpoolP512r1,512,EXPONENTIATION,50000,0.033214152,1505382.411,0
opencl-e2e,NVIDIA A40,GPU,w32-il,brainpoolP512r1,512,EXPONENTIATION,50000,0.033990106,1471016.318,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,brainpoolP512r1,512,EXPONENTIATION,50000,0.025381757,1969918.793,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,brainpoolP512r1,512,EXPONENTIATION,50000,0.026400292,1893918.425,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,brainpoolP512r1,512,DIVIDE,3125,0.000182176,17153740.606,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,brainpoolP512r1,512,DIVIDE,3125,0.005111979,611309.202,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,brainpoolP512r1,512,DIVIDE,3125,0.000073180,42703156.180,0
library,NVIDIA A40,gpu,cgbn,brainpoolP512r1,512,DIVIDE,50000,0.000046080,1085069444.444,0
opencl-kernel,NVIDIA A40,GPU,w8,brainpoolP512r1,512,DIVIDE,50000,0.001385591,36085684.922,0
opencl-e2e,NVIDIA A40,GPU,w8,brainpoolP512r1,512,DIVIDE,50000,0.002535566,19719464.927,0
opencl-kernel,NVIDIA A40,GPU,w16,brainpoolP512r1,512,DIVIDE,50000,0.001241478,40274571.483,0
opencl-e2e,NVIDIA A40,GPU,w16,brainpoolP512r1,512,DIVIDE,50000,0.002459582,20328657.956,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,brainpoolP512r1,512,DIVIDE,50000,0.000435778,114737336.214,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,brainpoolP512r1,512,DIVIDE,50000,0.001571935,31807937.663,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,brainpoolP512r1,512,DIVIDE,50000,0.000418775,119395918.658,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,brainpoolP512r1,512,DIVIDE,50000,0.001559162,32068517.212,0
opencl-kernel,NVIDIA A40,GPU,w32-il,brainpoolP512r1,512,DIVIDE,50000,0.000418415,119498500.234,0
opencl-e2e,NVIDIA A40,GPU,w32-il,brainpoolP512r1,512,DIVIDE,50000,0.001530677,32665282.997,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,brainpoolP512r1,512,DIVIDE,50000,0.000410690,121746337.547,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,brainpoolP512r1,512,DIVIDE,50000,0.001528212,32717976.482,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,brainpoolP512r1,512,ISQRT,781,0.000151638,5150426.020,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,brainpoolP512r1,512,ISQRT,781,0.000040858,19114958.960,0
opencl-kernel,NVIDIA A40,GPU,w8,brainpoolP512r1,512,ISQRT,50000,0.028314088,1765905.397,0
opencl-e2e,NVIDIA A40,GPU,w8,brainpoolP512r1,512,ISQRT,50000,0.029279811,1707661.278,0
opencl-kernel,NVIDIA A40,GPU,w16,brainpoolP512r1,512,ISQRT,50000,0.028748908,1739196.489,0
opencl-e2e,NVIDIA A40,GPU,w16,brainpoolP512r1,512,ISQRT,50000,0.029780638,1678943.189,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,brainpoolP512r1,512,ISQRT,50000,0.005962422,8385853.771,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,brainpoolP512r1,512,ISQRT,50000,0.006899702,7246689.404,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,brainpoolP512r1,512,ISQRT,50000,0.005623629,8891055.435,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,brainpoolP512r1,512,ISQRT,50000,0.006562433,7619125.168,0
opencl-kernel,NVIDIA A40,GPU,w32-il,brainpoolP512r1,512,ISQRT,50000,0.005950501,8402653.561,0
opencl-e2e,NVIDIA A40,GPU,w32-il,brainpoolP512r1,512,ISQRT,50000,0.006881249,7266122.619,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,brainpoolP512r1,512,ISQRT,50000,0.006621105,7551609.970,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,brainpoolP512r1,512,ISQRT,50000,0.007549348,6623088.258,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL_R2,25000,0.005602211,4462524.005,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL_R2,25000,0.000170755,146408425.554,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL_R2,25000,0.000357749,69881408.898,0
library,NVIDIA A40,gpu,cgbn,brainpoolP512r1,512,MODMUL_R2,50000,0.000041984,1190929878.049,0
opencl-kernel,NVIDIA A40,GPU,w8,brainpoolP512r1,512,MODMUL_R2,50000,0.000704658,70956385.297,0
opencl-e2e,NVIDIA A40,GPU,w8,brainpoolP512r1,512,MODMUL_R2,50000,0.001616661,30927952.918,0
opencl-kernel,NVIDIA A40,GPU,w16,brainpoolP512r1,512,MODMUL_R2,50000,0.000141920,352312177.708,0
opencl-e2e,NVIDIA A40,GPU,w16,brainpoolP512r1,512,MODMUL_R2,50000,0.001127171,44358828.744,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,brainpoolP512r1,512,MODMUL_R2,50000,0.000107214,466357637.248,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,brainpoolP512r1,512,MODMUL_R2,50000,0.000974100,51329427.917,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,brainpoolP512r1,512,MODMUL_R2,50000,0.000088078,567678842.799,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,brainpoolP512r1,512,MODMUL_R2,50000,0.000966306,51743451.916,0
opencl-kernel,NVIDIA A40,GPU,w32-il,brainpoolP512r1,512,MODMUL_R2,50000,0.000106771,468289861.747,0
opencl-e2e,NVIDIA A40,GPU,w32-il,brainpoolP512r1,512,MODMUL_R2,50000,0.001000881,49955979.038,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,brainpoolP512r1,512,MODMUL_R2,50000,0.000079642,627809053.382,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,brainpoolP512r1,512,MODMUL_R2,50000,0.000926249,53981180.672,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,p1024,1024,ADD,12500,0.000313105,39922701.773,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,p1024,1024,ADD,12500,0.000040647,307528475.850,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,p1024,1024,ADD,12500,0.000064272,194487441.133,0
library,NVIDIA A40,gpu,cgbn,p1024,1024,ADD,50000,0.000059392,841864224.138,0
opencl-kernel,NVIDIA A40,GPU,w8,p1024,1024,ADD,50000,0.000275454,181518806.083,0
opencl-e2e,NVIDIA A40,GPU,w8,p1024,1024,ADD,50000,0.002226729,22454461.775,0
opencl-kernel,NVIDIA A40,GPU,w16,p1024,1024,ADD,50000,0.000152950,326903842.805,0
opencl-e2e,NVIDIA A40,GPU,w16,p1024,1024,ADD,50000,0.002095100,23865214.554,0
opencl-kernel,NVIDIA A40,GPU,w32,p1024,1024,ADD,50000,0.000465253,107468324.652,0
opencl-e2e,NVIDIA A40,GPU,w32,p1024,1024,ADD,50000,0.002228873,22432863.314,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,p1024,1024,ADD,50000,0.000090362,553332555.527,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,p1024,1024,ADD,50000,0.002278116,21947964.358,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,p1024,1024,ADD,50000,0.000090712,551196509.276,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,p1024,1024,ADD,50000,0.001818594,27493773.339,0
opencl-kernel,NVIDIA A40,GPU,w32-il,p1024,1024,ADD,50000,0.000058661,852351932.939,0
opencl-e2e,NVIDIA A40,GPU,w32-il,p1024,1024,ADD,50000,0.001727359,28945916.458,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,p1024,1024,ADD,50000,0.000059773,836495087.331,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,p1024,1024,ADD,50000,0.001767858,28282815.665,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,p1024,1024,SUBTRACT,12500,0.000248291,50344234.059,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,p1024,1024,SUBTRACT,12500,0.010005172,1249353.791,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,p1024,1024,SUBTRACT,12500,0.000060384,207007924.488,0
library,NVIDIA A40,gpu,cgbn,p1024,1024,SUBTRACT,50000,0.000059392,841864224.138,0
opencl-kernel,NVIDIA A40,GPU,w8,p1024,1024,SUBTRACT,50000,0.000281414,177674163.208,0
opencl-e2e,NVIDIA A40,GPU,w8,p1024,1024,SUBTRACT,50000,0.002212723,22596594.712,0
opencl-kernel,NVIDIA A40,GPU,w16,p1024,1024,SUBTRACT,50000,0.000155005,322570920.780,0
opencl-e2e,NVIDIA A40,GPU,w16,p1024,1024,SUBTRACT,50000,0.002015429,24808618.655,0
opencl-kernel,NVIDIA A40,GPU,w32,p1024,1024,SUBTRACT,50000,0.000465694,107366666.733,0
opencl-e2e,NVIDIA A40,GPU,w32,p1024,1024,SUBTRACT,50000,0.002202282,22703725.575,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,p1024,1024,SUBTRACT,50000,0.000089149,560858844.792,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,p1024,1024,SUBTRACT,50000,0.002259551,22128293.849,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,p1024,1024,SUBTRACT,50000,0.000090702,551253105.523,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,p1024,1024,SUBTRACT,50000,0.001821229,27453984.955,0
opencl-kernel,NVIDIA A40,GPU,w32-il,p1024,1024,SUBTRACT,50000,0.000059172,845000254.978,0
opencl-e2e,NVIDIA A40,GPU,w32-il,p1024,1024,SUBTRACT,50000,0.001711160,29219945.312,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,p1024,1024,SUBTRACT,50000,0.000059733,837055898.219,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,p1024,1024,SUBTRACT,50000,0.001716160,29134808.947,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,p1024,1024,ADDMOD,12500,0.000923134,13540823.722,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,p1024,1024,ADDMOD,12500,0.000076365,163688140.885,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,p1024,1024,ADDMOD,12500,0.000144594,86448744.984,0
library,NVIDIA A40,gpu,cgbn,p1024,1024,ADDMOD,50000,0.000059392,841864224.138,0
opencl-kernel,NVIDIA A40,GPU,w8,p1024,1024,ADDMOD,50000,0.000418194,119561837.909,0
opencl-e2e,NVIDIA A40,GPU,w8,p1024,1024,ADDMOD,50000,0.002301271,21727123.143,0
opencl-kernel,NVIDIA A40,GPU,w16,p1024,1024,ADDMOD,50000,0.000217463,229924287.470,0
opencl-e2e,NVIDIA A40,GPU,w16,p1024,1024,ADDMOD,50000,0.002153660,23216290.425,0
opencl-kernel,NVIDIA A40,GPU,w32,p1024,1024,ADDMOD,50000,0.000456347,109565715.848,0
opencl-e2e,NVIDIA A40,GPU,w32,p1024,1024,ADDMOD,50000,0.002211800,22606023.828,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,p1024,1024,ADDMOD,50000,0.000112585,444109715.685,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,p1024,1024,ADDMOD,50000,0.002219033,22532332.943,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,p1024,1024,ADDMOD,50000,0.000113265,441443969.182,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,p1024,1024,ADDMOD,50000,0.001826949,27368026.459,0
opencl-kernel,NVIDIA A40,GPU,w32-il,p1024,1024,ADDMOD,50000,0.000062449,800654564.977,0
opencl-e2e,NVIDIA A40,GPU,w32-il,p1024,1024,ADDMOD,50000,0.001713063,29187490.833,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,p1024,1024,ADDMOD,50000,0.000063771,784050752.110,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,p1024,1024,ADDMOD,50000,0.001648812,30324871.131,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,p1024,1024,SUBTRACTMOD,12500,0.000843601,14817426.171,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,p1024,1024,SUBTRACTMOD,12500,0.000053432,233942912.919,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,p1024,1024,SUBTRACTMOD,12500,0.000152359,82043184.958,0
library,NVIDIA A40,gpu,cgbn,p1024,1024,SUBTRACTMOD,50000,0.000059392,841864224.138,0
opencl-kernel,NVIDIA A40,GPU,w8,p1024,1024,SUBTRACTMOD,50000,0.000415089,120456208.464,0
opencl-e2e,NVIDIA A40,GPU,w8,p1024,1024,SUBTRACTMOD,50000,0.002287074,21861994.920,0
opencl-kernel,NVIDIA A40,GPU,w16,p1024,1024,SUBTRACTMOD,50000,0.000218625,228701926.755,0
opencl-e2e,NVIDIA A40,GPU,w16,p1024,1024,SUBTRACTMOD,50000,0.002125577,23523023.799,0
opencl-kernel,NVIDIA A40,GPU,w32,p1024,1024,SUBTRACTMOD,50000,0.000462599,108084949.025,0
opencl-e2e,NVIDIA A40,GPU,w32,p1024,1024,SUBTRACTMOD,50000,0.002192754,22802372.357,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,p1024,1024,SUBTRACTMOD,50000,0.000112874,442970108.418,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,p1024,1024,SUBTRACTMOD,50000,0.002213263,22591079.806,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,p1024,1024,SUBTRACTMOD,50000,0.000113355,441092160.310,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,p1024,1024,SUBTRACTMOD,50000,0.001800070,27776703.970,0
opencl-kernel,NVIDIA A40,GPU,w32-il,p1024,1024,SUBTRACTMOD,50000,0.000062048,805822094.140,0
opencl-e2e,NVIDIA A40,GPU,w32-il,p1024,1024,SUBTRACTMOD,50000,0.001720637,29059004.070,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,p1024,1024,SUBTRACTMOD,50000,0.000061447,813712013.095,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,p1024,1024,SUBTRACTMOD,50000,0.001641148,30466482.613,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.002543761,4913982.901,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.000278289,44917264.759,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.013968059,894898.863,0
opencl-kernel,NVIDIA A40,GPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.041558023,1203137.129,0
opencl-e2e,NVIDIA A40,GPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.044271655,1129390.794,0
opencl-kernel,NVIDIA A40,GPU,w16,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.013452673,3716733.379,0
opencl-e2e,NVIDIA A40,GPU,w16,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.015872629,3150076.706,0
opencl-kernel,NVIDIA A40,GPU,w32,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.019933823,2508299.525,0
opencl-e2e,NVIDIA A40,GPU,w32,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.022182073,2254072.421,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.000294259,169918442.325,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.003099116,16133635.208,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.000298707,167388213.286,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.002406711,20775242.649,0
opencl-kernel,NVIDIA A40,GPU,w32-il,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.000228414,218900468.894,0
opencl-e2e,NVIDIA A40,GPU,w32-il,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.002365083,21140911.554,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.000225018,222204664.523,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.002218463,22538121.968,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.002536427,4928191.850,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.000098597,126778373.068,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.000096804,129126278.826,0
library,NVIDIA A40,gpu,cgbn,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000060160,831117021.277,0
opencl-kernel,NVIDIA A40,GPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.014430629,3464852.411,0
opencl-e2e,NVIDIA A40,GPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.016915837,2955809.937,0
opencl-kernel,NVIDIA A40,GPU,w16,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.003691301,13545360.671,0
opencl-e2e,NVIDIA A40,GPU,w16,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.006177313,8094134.533,0
opencl-kernel,NVIDIA A40,GPU,w32,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.001896592,26363078.494,0
opencl-e2e,NVIDIA A40,GPU,w32,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.004058477,12319893.817,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000959353,52118481.269,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.003779838,13228080.629,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000959893,52089152.276,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.003061785,16330341.313,0
opencl-kernel,NVIDIA A40,GPU,w32-il,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000411732,121438181.020,0
opencl-e2e,NVIDIA A40,GPU,w32-il,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.002543240,19659962.436,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000380954,131249520.103,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.002458619,20336620.249,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.017390821,718769.971,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.000429446,29107276.725,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.000201533,62024699.391,0
library,NVIDIA A40,gpu,cgbn,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000058368,856633771.930,0
opencl-kernel,NVIDIA A40,GPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.002545225,19644632.506,0
opencl-e2e,NVIDIA A40,GPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.004528491,11041205.981,0
opencl-kernel,NVIDIA A40,GPU,w16,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000500009,99998120.999,0
opencl-e2e,NVIDIA A40,GPU,w16,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.002453200,20381545.900,0
opencl-kernel,NVIDIA A40,GPU,w32,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000142732,350307595.738,0
opencl-e2e,NVIDIA A40,GPU,w32,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.001844263,27111105.096,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000216932,230486934.645,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.002302934,21711430.278,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000174942,285808314.390,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.001839894,27175480.507,0
opencl-kernel,NVIDIA A40,GPU,w32-il,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000185592,269408018.948,0
opencl-e2e,NVIDIA A40,GPU,w32-il,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.001886894,26498573.915,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000139184,359235929.554,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.001716701,29125625.755,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,p1024,1024,COMPARE,12500,0.000115179,108526297.575,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,p1024,1024,COMPARE,12500,0.000045616,274025577.787,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,p1024,1024,COMPARE,12500,0.000067659,184750754.322,0
library,NVIDIA A40,gpu,cgbn,p1024,1024,COMPARE,50000,0.000059360,842318059.299,0
opencl-kernel,NVIDIA A40,GPU,w8,p1024,1024,COMPARE,50000,0.000251489,198816042.424,0
opencl-e2e,NVIDIA A40,GPU,w8,p1024,1024,COMPARE,50000,0.002205148,22674211.575,0
opencl-kernel,NVIDIA A40,GPU,w16,p1024,1024,COMPARE,50000,0.000132982,375991618.343,0
opencl-e2e,NVIDIA A40,GPU,w16,p1024,1024,COMPARE,50000,0.001995530,25056000.418,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,p1024,1024,COMPARE,50000,0.000075033,666374043.641,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,p1024,1024,COMPARE,50000,0.002103725,23767370.489,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,p1024,1024,COMPARE,50000,0.000075885,658890921.810,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,p1024,1024,COMPARE,50000,0.001803366,27725937.421,0
opencl-kernel,NVIDIA A40,GPU,w32-il,p1024,1024,COMPARE,50000,0.000031420,1591341589.353,0
opencl-e2e,NVIDIA A40,GPU,w32-il,p1024,1024,COMPARE,50000,0.001733412,28844844.155,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,p1024,1024,COMPARE,50000,0.000030217,1654710778.240,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,p1024,1024,COMPARE,50000,0.001678388,29790495.934,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,p1024,1024,REDUCE,1562,0.000044625,35002603.078,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,p1024,1024,REDUCE,1562,0.009849366,158588.891,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,p1024,1024,REDUCE,1562,0.008043708,194189.041,0
library,NVIDIA A40,gpu,cgbn,p1024,1024,REDUCE,50000,0.000059392,841864224.138,0
opencl-kernel,NVIDIA A40,GPU,w8,p1024,1024,REDUCE,50000,0.003785229,13209239.440,0
opencl-e2e,NVIDIA A40,GPU,w8,p1024,1024,REDUCE,50000,0.005750973,8694181.304,0
opencl-kernel,NVIDIA A40,GPU,w16,p1024,1024,REDUCE,50000,0.002059091,24282561.020,0
opencl-e2e,NVIDIA A40,GPU,w16,p1024,1024,REDUCE,50000,0.003981583,12557819.906,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,p1024,1024,REDUCE,50000,0.000521971,95790778.322,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,p1024,1024,REDUCE,50000,0.002617040,19105556.427,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,p1024,1024,REDUCE,50000,0.000528534,94601303.592,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,p1024,1024,REDUCE,50000,0.002292464,21810588.682,0
opencl-kernel,NVIDIA A40,GPU,w32-il,p1024,1024,REDUCE,50000,0.000451457,110752579.072,0
opencl-e2e,NVIDIA A40,GPU,w32-il,p1024,1024,REDUCE,50000,0.002148370,23273455.682,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,p1024,1024,REDUCE,50000,0.000499879,100024203.853,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,p1024,1024,REDUCE,50000,0.002191392,22816550.014,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,p1024,1024,MODMUL,781,0.000500941,1559066.717,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,p1024,1024,MODMUL,781,0.000053351,14638951.986,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,p1024,1024,MODMUL,781,0.000094289,8283048.187,0
library,NVIDIA A40,gpu,cgbn,p1024,1024,MODMUL,50000,0.000470944,106169735.680,0
opencl-kernel,NVIDIA A40,GPU,w8,p1024,1024,MODMUL,50000,0.015391163,3248617.490,0
opencl-e2e,NVIDIA A40,GPU,w8,p1024,1024,MODMUL,50000,0.017418623,2870490.973,0
opencl-kernel,NVIDIA A40,GPU,w16,p1024,1024,MODMUL,50000,0.005471664,9137988.757,0
opencl-e2e,NVIDIA A40,GPU,w16,p1024,1024,MODMUL,50000,0.007449010,6712301.017,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,p1024,1024,MODMUL,50000,0.002424825,20620044.891,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,p1024,1024,MODMUL,50000,0.004520865,11059829.980,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,p1024,1024,MODMUL,50000,0.001817523,27509974.738,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,p1024,1024,MODMUL,50000,0.003559420,14047232.867,0
opencl-kernel,NVIDIA A40,GPU,w32-il,p1024,1024,MODMUL,50000,0.002364220,21148623.197,0
opencl-e2e,NVIDIA A40,GPU,w32-il,p1024,1024,MODMUL,50000,0.004072235,12278269.898,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,p1024,1024,MODMUL,50000,0.001792636,27891890.184,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,p1024,1024,MODMUL,50000,0.003504396,14267792.242,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,p1024,1024,MODEXP,195,0.075262667,2590.926,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,p1024,1024,MODEXP,195,0.001581413,123307.467,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,p1024,1024,MODEXP,195,0.001387916,140498.459,0
library,NVIDIA A40,gpu,cgbn,p1024,1024,MODEXP,50000,0.241587207,206964.601,0
opencl-kernel,NVIDIA A40,GPU,w8,p1024,1024,MODEXP,50000,3.555417262,14063.047,0
opencl-e2e,NVIDIA A40,GPU,w8,p1024,1024,MODEXP,50000,3.559102121,14048.487,0
opencl-kernel,NVIDIA A40,GPU,w16,p1024,1024,MODEXP,50000,0.517323680,96651.288,0
opencl-e2e,NVIDIA A40,GPU,w16,p1024,1024,MODEXP,50000,0.519471068,96251.751,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,p1024,1024,MODEXP,50000,0.235262679,212528.397,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,p1024,1024,MODEXP,50000,0.238019865,210066.500,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,p1024,1024,MODEXP,50000,0.148921516,335747.321,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,p1024,1024,MODEXP,50000,0.151130361,330840.209,0
opencl-kernel,NVIDIA A40,GPU,w32-il,p1024,1024,MODEXP,50000,0.237956199,210122.704,0
opencl-e2e,NVIDIA A40,GPU,w32-il,p1024,1024,MODEXP,50000,0.239855788,208458.593,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,p1024,1024,MODEXP,50000,0.149737145,333918.481,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,p1024,1024,MODEXP,50000,0.151919069,329122.607,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,p1024,1024,EXPONENTIATION,195,0.011139928,17504.601,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,p1024,1024,EXPONENTIATION,195,0.000399679,487891.599,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,p1024,1024,EXPONENTIATION,195,0.006217078,31365.215,0
opencl-kernel,NVIDIA A40,GPU,w8,p1024,1024,EXPONENTIATION,50000,2.796067521,17882.258,0
opencl-e2e,NVIDIA A40,GPU,w8,p1024,1024,EXPONENTIATION,50000,2.810551320,17790.104,0
opencl-kernel,NVIDIA A40,GPU,w16,p1024,1024,EXPONENTIATION,50000,0.737079193,67835.316,0
opencl-e2e,NVIDIA A40,GPU,w16,p1024,1024,EXPONENTIATION,50000,0.739072878,67652.327,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,p1024,1024,EXPONENTIATION,50000,0.202767520,246587.816,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,p1024,1024,EXPONENTIATION,50000,0.204829918,244104.965,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,p1024,1024,EXPONENTIATION,50000,0.190693296,262201.143,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,p1024,1024,EXPONENTIATION,50000,0.193558209,258320.225,0
opencl-kernel,NVIDIA A40,GPU,w32-il,p1024,1024,EXPONENTIATION,50000,0.205418763,243405.224,0
opencl-e2e,NVIDIA A40,GPU,w32-il,p1024,1024,EXPONENTIATION,50000,0.207359219,241127.451,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,p1024,1024,EXPONENTIATION,50000,0.194004801,257725.581,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,p1024,1024,EXPONENTIATION,50000,0.196182598,254864.603,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,p1024,1024,DIVIDE,1562,0.000096263,16226318.464,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,p1024,1024,DIVIDE,1562,0.000035388,44138763.332,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,p1024,1024,DIVIDE,1562,0.000054575,28621388.233,0
library,NVIDIA A40,gpu,cgbn,p1024,1024,DIVIDE,50000,0.000061408,814226159.458,0
opencl-kernel,NVIDIA A40,GPU,w8,p1024,1024,DIVIDE,50000,0.089012693,561717.644,0
opencl-e2e,NVIDIA A40,GPU,w8,p1024,1024,DIVIDE,50000,0.091594867,545882.116,0
opencl-kernel,NVIDIA A40,GPU,w16,p1024,1024,DIVIDE,50000,0.009600456,5208086.018,0
opencl-e2e,NVIDIA A40,GPU,w16,p1024,1024,DIVIDE,50000,0.011975006,4175363.183,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,p1024,1024,DIVIDE,50000,0.001802413,27740593.161,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,p1024,1024,DIVIDE,50000,0.004653587,10744399.951,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,p1024,1024,DIVIDE,50000,0.001751716,28543435.638,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,p1024,1024,DIVIDE,50000,0.004092623,12217102.367,0
opencl-kernel,NVIDIA A40,GPU,w32-il,p1024,1024,DIVIDE,50000,0.001798716,27797600.961,0
opencl-e2e,NVIDIA A40,GPU,w32-il,p1024,1024,DIVIDE,50000,0.003910067,12787506.104,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,p1024,1024,DIVIDE,50000,0.001763139,28358514.143,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,p1024,1024,DIVIDE,50000,0.004076833,12264422.202,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,p1024,1024,ISQRT,390,0.000151488,2574461.366,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,p1024,1024,ISQRT,390,0.000036519,10679366.300,0
opencl-kernel,NVIDIA A40,GPU,w8,p1024,1024,ISQRT,50000,0.808797454,61820.175,0
opencl-e2e,NVIDIA A40,GPU,w8,p1024,1024,ISQRT,50000,0.809978216,61730.055,0
opencl-kernel,NVIDIA A40,GPU,w16,p1024,1024,ISQRT,50000,0.203732758,245419.541,0
opencl-e2e,NVIDIA A40,GPU,w16,p1024,1024,ISQRT,50000,0.205735791,243030.149,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,p1024,1024,ISQRT,50000,0.039805602,1256104.597,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,p1024,1024,ISQRT,50000,0.042089869,1187934.315,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,p1024,1024,ISQRT,50000,0.037426479,1335952.543,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,p1024,1024,ISQRT,50000,0.039436847,1267849.837,0
opencl-kernel,NVIDIA A40,GPU,w32-il,p1024,1024,ISQRT,50000,0.039385403,1269505.876,0
opencl-e2e,NVIDIA A40,GPU,w32-il,p1024,1024,ISQRT,50000,0.041000601,1219494.322,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,p1024,1024,ISQRT,50000,0.036620713,1365347.535,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,p1024,1024,ISQRT,50000,0.038581046,1295973.148,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,p1024,1024,MODMUL_R2,12500,0.007795125,1603566.281,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,p1024,1024,MODMUL_R2,12500,0.000160264,77996378.471,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,p1024,1024,MODMUL_R2,12500,0.008055829,1551671.381,0
library,NVIDIA A40,gpu,cgbn,p1024,1024,MODMUL_R2,50000,0.000106496,469501201.923,0
opencl-kernel,NVIDIA A40,GPU,w8,p1024,1024,MODMUL_R2,50000,0.004110618,12163619.676,0
opencl-e2e,NVIDIA A40,GPU,w8,p1024,1024,MODMUL_R2,50000,0.005905726,8466359.607,0
opencl-kernel,NVIDIA A40,GPU,w16,p1024,1024,MODMUL_R2,50000,0.000577728,86545878.249,0
opencl-e2e,NVIDIA A40,GPU,w16,p1024,1024,MODMUL_R2,50000,0.002450715,20402210.657,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,p1024,1024,MODMUL_R2,50000,0.000343251,145665981.664,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,p1024,1024,MODMUL_R2,50000,0.002531828,19748574.119,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,p1024,1024,MODMUL_R2,50000,0.000288969,173028997.222,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,p1024,1024,MODMUL_R2,50000,0.002213924,22584332.464,0
opencl-kernel,NVIDIA A40,GPU,w32-il,p1024,1024,MODMUL_R2,50000,0.000322923,154835642.102,0
opencl-e2e,NVIDIA A40,GPU,w32-il,p1024,1024,MODMUL_R2,50000,0.001869099,26750853.891,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,p1024,1024,MODMUL_R2,50000,0.000260485,191949384.685,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,p1024,1024,MODMUL_R2,50000,0.002173368,23005769.195,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,p2048,2048,ADD,6250,0.000210901,29634786.909,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,p2048,2048,ADD,6250,0.010895594,573626.366,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,p2048,2048,ADD,6250,0.007997590,781485.411,0
library,NVIDIA A40,gpu,cgbn,p2048,2048,ADD,50000,0.000111616,447964449.541,0
opencl-kernel,NVIDIA A40,GPU,w8,p2048,2048,ADD,50000,0.000565064,88485530.209,0
opencl-e2e,NVIDIA A40,GPU,w8,p2048,2048,ADD,50000,0.004083246,12245159.878,0
opencl-kernel,NVIDIA A40,GPU,w16,p2048,2048,ADD,50000,0.000292527,170924651.542,0
opencl-e2e,NVIDIA A40,GPU,w16,p2048,2048,ADD,50000,0.003981953,12556650.942,0
opencl-kernel,NVIDIA A40,GPU,w32,p2048,2048,ADD,50000,0.000195242,256093051.388,0
opencl-e2e,NVIDIA A40,GPU,w32,p2048,2048,ADD,50000,0.003805637,13138403.709,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,p2048,2048,ADD,50000,0.000186143,268610052.534,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,p2048,2048,ADD,50000,0.004106319,12176356.765,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,p2048,2048,ADD,50000,0.000188088,265832951.406,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,p2048,2048,ADD,50000,0.004148058,12053834.367,0
opencl-kernel,NVIDIA A40,GPU,w32-il,p2048,2048,ADD,50000,0.000107024,467185519.858,0
opencl-e2e,NVIDIA A40,GPU,w32-il,p2048,2048,ADD,50000,0.003327420,15026657.374,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,p2048,2048,ADD,50000,0.000107804,463803334.658,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,p2048,2048,ADD,50000,0.004027951,12413260.550,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,p2048,2048,SUBTRACT,6250,0.000177787,35154488.546,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,p2048,2048,SUBTRACT,6250,0.011965927,522316.419,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,p2048,2048,SUBTRACT,6250,0.007606608,821654.029,0
library,NVIDIA A40,gpu,cgbn,p2048,2048,SUBTRACT,50000,0.000112512,444397042.093,0
opencl-kernel,NVIDIA A40,GPU,w8,p2048,2048,SUBTRACT,50000,0.000604558,82705077.341,0
opencl-e2e,NVIDIA A40,GPU,w8,p2048,2048,SUBTRACT,50000,0.004003993,12487533.584,0
opencl-kernel,NVIDIA A40,GPU,w16,p2048,2048,SUBTRACT,50000,0.000292906,170703458.451,0
opencl-e2e,NVIDIA A40,GPU,w16,p2048,2048,SUBTRACT,50000,0.003986652,12541852.132,0
opencl-kernel,NVIDIA A40,GPU,w32,p2048,2048,SUBTRACT,50000,0.000190763,262105605.624,0
opencl-e2e,NVIDIA A40,GPU,w32,p2048,2048,SUBTRACT,50000,0.003813501,13111310.088,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,p2048,2048,SUBTRACT,50000,0.000186374,268277172.468,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,p2048,2048,SUBTRACT,50000,0.004151414,12044088.643,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,p2048,2048,SUBTRACT,50000,0.000185742,269190535.452,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,p2048,2048,SUBTRACT,50000,0.004130935,12103796.456,0
opencl-kernel,NVIDIA A40,GPU,w32-il,p2048,2048,SUBTRACT,50000,0.000108096,462552587.730,0
opencl-e2e,NVIDIA A40,GPU,w32-il,p2048,2048,SUBTRACT,50000,0.003360381,14879263.629,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,p2048,2048,SUBTRACT,50000,0.000108225,461999304.683,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,p2048,2048,SUBTRACT,50000,0.004008563,12473297.055,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,p2048,2048,ADDMOD,6250,0.000603386,10358208.375,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,p2048,2048,ADDMOD,6250,0.009999943,625003.565,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,p2048,2048,ADDMOD,6250,0.000091354,68414904.527,0
library,NVIDIA A40,gpu,cgbn,p2048,2048,ADDMOD,50000,0.000111616,447964449.541,0
opencl-kernel,NVIDIA A40,GPU,w8,p2048,2048,ADDMOD,50000,0.000764873,65370339.813,0
opencl-e2e,NVIDIA A40,GPU,w8,p2048,2048,ADDMOD,50000,0.004232739,11812681.774,0
opencl-kernel,NVIDIA A40,GPU,w16,p2048,2048,ADDMOD,50000,0.000388337,128754049.893,0
opencl-e2e,NVIDIA A40,GPU,w16,p2048,2048,ADDMOD,50000,0.004084338,12241887.446,0
opencl-kernel,NVIDIA A40,GPU,w32,p2048,2048,ADDMOD,50000,0.000222062,225162584.666,0
opencl-e2e,NVIDIA A40,GPU,w32,p2048,2048,ADDMOD,50000,0.003878776,12890664.531,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,p2048,2048,ADDMOD,50000,0.000242250,206398314.585,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,p2048,2048,ADDMOD,50000,0.004206167,11887308.498,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,p2048,2048,ADDMOD,50000,0.000225298,221928185.456,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,p2048,2048,ADDMOD,50000,0.004201189,11901393.574,0
opencl-kernel,NVIDIA A40,GPU,w32-il,p2048,2048,ADDMOD,50000,0.000135398,369282931.862,0
opencl-e2e,NVIDIA A40,GPU,w32-il,p2048,2048,ADDMOD,50000,0.003358969,14885517.872,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,p2048,2048,ADDMOD,50000,0.000133453,374663916.144,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,p2048,2048,ADDMOD,50000,0.004050463,12344266.468,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,p2048,2048,SUBTRACTMOD,6250,0.000507193,12322732.898,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,p2048,2048,SUBTRACTMOD,6250,0.008174785,764546.063,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,p2048,2048,SUBTRACTMOD,6250,0.000128815,48519212.806,0
library,NVIDIA A40,gpu,cgbn,p2048,2048,SUBTRACTMOD,50000,0.000110528,452374059.062,0
opencl-kernel,NVIDIA A40,GPU,w8,p2048,2048,SUBTRACTMOD,50000,0.000816100,61266980.797,0
opencl-e2e,NVIDIA A40,GPU,w8,p2048,2048,SUBTRACTMOD,50000,0.004384328,11404257.084,0
opencl-kernel,NVIDIA A40,GPU,w16,p2048,2048,SUBTRACTMOD,50000,0.000423986,117928292.902,0
opencl-e2e,NVIDIA A40,GPU,w16,p2048,2048,SUBTRACTMOD,50000,0.004123392,12125940.257,0
opencl-kernel,NVIDIA A40,GPU,w32,p2048,2048,SUBTRACTMOD,50000,0.000241389,207134913.905,0
opencl-e2e,NVIDIA A40,GPU,w32,p2048,2048,SUBTRACTMOD,50000,0.003852185,12979647.081,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,p2048,2048,SUBTRACTMOD,50000,0.000222422,224797722.172,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,p2048,2048,SUBTRACTMOD,50000,0.004196239,11915432.731,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,p2048,2048,SUBTRACTMOD,50000,0.000221310,225927244.876,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,p2048,2048,SUBTRACTMOD,50000,0.004162885,12010900.376,0
opencl-kernel,NVIDIA A40,GPU,w32-il,p2048,2048,SUBTRACTMOD,50000,0.000138884,360011608.975,0
opencl-e2e,NVIDIA A40,GPU,w32-il,p2048,2048,SUBTRACTMOD,50000,0.003871514,12914845.733,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,p2048,2048,SUBTRACTMOD,50000,0.000130448,383295788.444,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,p2048,2048,SUBTRACTMOD,50000,0.004022189,12431040.617,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.004373347,1429111.261,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.000111070,56270586.361,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.000204799,30517762.084,0
opencl-kernel,NVIDIA A40,GPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.214432395,233173.724,0
opencl-e2e,NVIDIA A40,GPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.219525489,227763.984,0
opencl-kernel,NVIDIA A40,GPU,w16,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.080733873,619318.734,0
opencl-e2e,NVIDIA A40,GPU,w16,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.085538289,584533.556,0
opencl-kernel,NVIDIA A40,GPU,w32,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.023572945,2121075.662,0
opencl-e2e,NVIDIA A40,GPU,w32,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.026989413,1852578.239,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.001573740,31771457.518,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.006360729,7860733.924,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.001622400,30818535.750,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.006464697,7734314.036,0
opencl-kernel,NVIDIA A40,GPU,w32-il,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.001584860,31548518.068,0
opencl-e2e,NVIDIA A40,GPU,w32-il,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.005699124,8773277.851,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.001627872,30714950.049,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.006360851,7860583.152,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.004309635,1450238.934,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.000158351,39469301.527,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.007601418,822215.056,0
library,NVIDIA A40,gpu,cgbn,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.000119808,417334401.709,0
opencl-kernel,NVIDIA A40,GPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.057150527,874882.573,0
opencl-e2e,NVIDIA A40,GPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.061546448,812394.572,0
opencl-kernel,NVIDIA A40,GPU,w16,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.014625110,3418777.768,0
opencl-e2e,NVIDIA A40,GPU,w16,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.019551297,2557375.131,0
opencl-kernel,NVIDIA A40,GPU,w32,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.003824513,13073558.568,0
opencl-e2e,NVIDIA A40,GPU,w32,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.008559664,5841350.932,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.004259158,11739409.731,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.008695743,5749939.911,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.004241164,11789216.795,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.009121023,5481841.165,0
opencl-kernel,NVIDIA A40,GPU,w32-il,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.002019265,24761488.242,0
opencl-e2e,NVIDIA A40,GPU,w32-il,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.006264948,7980911.847,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.001914507,26116387.036,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.006800456,7352448.160,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.028101865,222405.166,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.007113330,878632.081,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.000219257,28505411.065,0
library,NVIDIA A40,gpu,cgbn,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000176128,283884447.674,0
opencl-kernel,NVIDIA A40,GPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.110356824,453075.745,0
opencl-e2e,NVIDIA A40,GPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.114300343,437444.006,0
opencl-kernel,NVIDIA A40,GPU,w16,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.002247499,22246947.954,0
opencl-e2e,NVIDIA A40,GPU,w16,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.006233781,8020814.791,0
opencl-kernel,NVIDIA A40,GPU,w32,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000429836,116323407.427,0
opencl-e2e,NVIDIA A40,GPU,w32,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.004101991,12189203.581,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000671687,74439477.661,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.004630303,10798429.366,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000601803,83083673.460,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.004546304,10997944.140,0
opencl-kernel,NVIDIA A40,GPU,w32-il,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000675433,74026552.829,0
opencl-e2e,NVIDIA A40,GPU,w32-il,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.003911459,12782954.249,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000601794,83084959.245,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.004520425,11060907.761,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,p2048,2048,COMPARE,6250,0.000050747,123160388.335,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,p2048,2048,COMPARE,6250,0.000036799,169839961.532,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,p2048,2048,COMPARE,6250,0.000066726,93665979.036,0
library,NVIDIA A40,gpu,cgbn,p2048,2048,COMPARE,50000,0.000111616,447964449.541,0
opencl-kernel,NVIDIA A40,GPU,w8,p2048,2048,COMPARE,50000,0.000510018,98035694.303,0
opencl-e2e,NVIDIA A40,GPU,w8,p2048,2048,COMPARE,50000,0.003969008,12597605.926,0
opencl-kernel,NVIDIA A40,GPU,w16,p2048,2048,COMPARE,50000,0.000250345,199724303.773,0
opencl-e2e,NVIDIA A40,GPU,w16,p2048,2048,COMPARE,50000,0.004156475,12029424.100,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,p2048,2048,COMPARE,50000,0.000135927,367843256.984,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,p2048,2048,COMPARE,50000,0.004059860,12315694.150,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,p2048,2048,COMPARE,50000,0.000135918,367868461.913,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,p2048,2048,COMPARE,50000,0.004086281,12236064.501,0
opencl-kernel,NVIDIA A40,GPU,w32-il,p2048,2048,COMPARE,50000,0.000059803,836078227.150,0
opencl-e2e,NVIDIA A40,GPU,w32-il,p2048,2048,COMPARE,50000,0.003289738,15198779.955,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,p2048,2048,COMPARE,50000,0.000060044,832719493.734,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,p2048,2048,COMPARE,50000,0.003968768,12598368.624,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,p2048,2048,REDUCE,781,0.000033263,23479459.193,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,p2048,2048,REDUCE,781,0.000040827,19129348.158,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,p2048,2048,REDUCE,781,0.000078209,9986095.605,0
library,NVIDIA A40,gpu,cgbn,p2048,2048,REDUCE,50000,0.000112640,443892045.455,0
opencl-kernel,NVIDIA A40,GPU,w8,p2048,2048,REDUCE,50000,0.853517662,58581.096,0
opencl-e2e,NVIDIA A40,GPU,w8,p2048,2048,REDUCE,50000,0.857994085,58275.460,0
opencl-kernel,NVIDIA A40,GPU,w16,p2048,2048,REDUCE,50000,0.008874747,5633963.443,0
opencl-e2e,NVIDIA A40,GPU,w16,p2048,2048,REDUCE,50000,0.012913720,3871851.098,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,p2048,2048,REDUCE,50000,0.001706271,29303677.213,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,p2048,2048,REDUCE,50000,0.005712357,8752953.805,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,p2048,2048,REDUCE,50000,0.001788517,27956115.113,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,p2048,2048,REDUCE,50000,0.005734290,8719475.302,0
opencl-kernel,NVIDIA A40,GPU,w32-il,p2048,2048,REDUCE,50000,0.001619596,30871895.795,0
opencl-e2e,NVIDIA A40,GPU,w32-il,p2048,2048,REDUCE,50000,0.005172835,9665880.043,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,p2048,2048,REDUCE,50000,0.001826828,27369840.258,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,p2048,2048,REDUCE,50000,0.005811668,8603382.022,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,p2048,2048,MODMUL,390,0.000746939,522130.759,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,p2048,2048,MODMUL,390,0.000056488,6904181.349,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,p2048,2048,MODMUL,390,0.000116952,3334708.156,0
library,NVIDIA A40,gpu,cgbn,p2048,2048,MODMUL,50000,0.000855040,58476796.407,0
opencl-kernel,NVIDIA A40,GPU,w8,p2048,2048,MODMUL,50000,1.472832921,33948.182,0
opencl-e2e,NVIDIA A40,GPU,w8,p2048,2048,MODMUL,50000,1.471441528,33980.283,0
opencl-kernel,NVIDIA A40,GPU,w16,p2048,2048,MODMUL,50000,0.024287533,2058669.358,0
opencl-e2e,NVIDIA A40,GPU,w16,p2048,2048,MODMUL,50000,0.028413290,1759739.887,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,p2048,2048,MODMUL,50000,0.009752500,5126890.609,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,p2048,2048,MODMUL,50000,0.013866974,3605689.306,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,p2048,2048,MODMUL,50000,0.007658775,6528459.216,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,p2048,2048,MODMUL,50000,0.011725307,4264280.510,0
opencl-kernel,NVIDIA A40,GPU,w32-il,p2048,2048,MODMUL,50000,0.009168663,5453357.775,0
opencl-e2e,NVIDIA A40,GPU,w32-il,p2048,2048,MODMUL,50000,0.012727723,3928432.305,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,p2048,2048,MODMUL,50000,0.007356783,6796448.590,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,p2048,2048,MODMUL,50000,0.011412255,4381255.242,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,p2048,2048,MODEXP,97,0.273221252,355.024,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,p2048,2048,MODEXP,97,0.010536372,9206.205,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,p2048,2048,MODEXP,97,0.008403861,11542.314,0
library,NVIDIA A40,gpu,cgbn,p2048,2048,MODEXP,50000,1.343389630,37219.284,0
opencl-kernel,NVIDIA A40,GPU,w8,p2048,2048,MODEXP,50000,158.312622023,315.831,0
opencl-e2e,NVIDIA A40,GPU,w8,p2048,2048,MODEXP,50000,158.219258104,316.017,0
opencl-kernel,NVIDIA A40,GPU,w16,p2048,2048,MODEXP,50000,26.145272882,1912.392,0
opencl-e2e,NVIDIA A40,GPU,w16,p2048,2048,MODEXP,50000,25.997796724,1923.240,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,p2048,2048,MODEXP,50000,2.024635421,24695.804,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,p2048,2048,MODEXP,50000,2.030544695,24623.935,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,p2048,2048,MODEXP,50000,9.839008094,5081.813,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,p2048,2048,MODEXP,50000,9.842412126,5080.056,0
opencl-kernel,NVIDIA A40,GPU,w32-il,p2048,2048,MODEXP,50000,1.979190901,25262.849,0
opencl-e2e,NVIDIA A40,GPU,w32-il,p2048,2048,MODEXP,50000,1.987903680,25152.124,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,p2048,2048,MODEXP,50000,9.824901272,5089.110,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,p2048,2048,MODEXP,50000,9.815495946,5093.986,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,p2048,2048,EXPONENTIATION,97,0.032706278,2965.791,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,p2048,2048,EXPONENTIATION,97,0.001143782,84806.341,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,p2048,2048,EXPONENTIATION,97,0.015137440,6407.953,0
opencl-kernel,NVIDIA A40,GPU,w8,p2048,2048,EXPONENTIATION,50000,31.488708859,1587.871,0
opencl-e2e,NVIDIA A40,GPU,w8,p2048,2048,EXPONENTIATION,50000,31.483593524,1588.129,0
opencl-kernel,NVIDIA A40,GPU,w16,p2048,2048,EXPONENTIATION,50000,6.062017352,8248.079,0
opencl-e2e,NVIDIA A40,GPU,w16,p2048,2048,EXPONENTIATION,50000,6.061371471,8248.958,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,p2048,2048,EXPONENTIATION,50000,1.807391050,27664.185,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,p2048,2048,EXPONENTIATION,50000,1.812702791,27583.121,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,p2048,2048,EXPONENTIATION,50000,1.699493535,29420.530,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,p2048,2048,EXPONENTIATION,50000,1.703942736,29343.709,0
opencl-kernel,NVIDIA A40,GPU,w32-il,p2048,2048,EXPONENTIATION,50000,1.796247030,27835.815,0
opencl-e2e,NVIDIA A40,GPU,w32-il,p2048,2048,EXPONENTIATION,50000,1.798516359,27800.692,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,p2048,2048,EXPONENTIATION,50000,1.736227867,28798.063,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,p2048,2048,EXPONENTIATION,50000,1.738879154,28754.155,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,p2048,2048,DIVIDE,781,0.000063282,12341678.409,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,p2048,2048,DIVIDE,781,0.000038653,20205584.284,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,p2048,2048,DIVIDE,781,0.000075453,10350819.761,0
library,NVIDIA A40,gpu,cgbn,p2048,2048,DIVIDE,50000,0.000112640,443892045.455,0
opencl-kernel,NVIDIA A40,GPU,w8,p2048,2048,DIVIDE,50000,2.274110543,21986.618,0
opencl-e2e,NVIDIA A40,GPU,w8,p2048,2048,DIVIDE,50000,2.281710173,21913.388,0
opencl-kernel,NVIDIA A40,GPU,w16,p2048,2048,DIVIDE,50000,0.489989736,102042.954,0
opencl-e2e,NVIDIA A40,GPU,w16,p2048,2048,DIVIDE,50000,0.493964094,101221.932,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,p2048,2048,DIVIDE,50000,0.044296537,1128756.419,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,p2048,2048,DIVIDE,50000,0.047810749,1045789.930,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,p2048,2048,DIVIDE,50000,0.037129529,1346637.073,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,p2048,2048,DIVIDE,50000,0.042361886,1180306.271,0
opencl-kernel,NVIDIA A40,GPU,w32-il,p2048,2048,DIVIDE,50000,0.040915521,1222030.148,0
opencl-e2e,NVIDIA A40,GPU,w32-il,p2048,2048,DIVIDE,50000,0.045079519,1109151.141,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,p2048,2048,DIVIDE,50000,0.038190938,1309211.115,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,p2048,2048,DIVIDE,50000,0.042988930,1163090.129,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,p2048,2048,ISQRT,195,0.000117934,1653462.862,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,p2048,2048,ISQRT,195,0.000037692,5173444.744,0
opencl-kernel,NVIDIA A40,GPU,w8,p2048,2048,ISQRT,50000,30.468821829,1641.022,0
opencl-e2e,NVIDIA A40,GPU,w8,p2048,2048,ISQRT,50000,30.456795557,1641.670,0
opencl-kernel,NVIDIA A40,GPU,w16,p2048,2048,ISQRT,50000,8.043744151,6216.011,0
opencl-e2e,NVIDIA A40,GPU,w16,p2048,2048,ISQRT,50000,8.053599242,6208.404,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,p2048,2048,ISQRT,50000,0.100595706,497039.108,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,p2048,2048,ISQRT,50000,0.104123445,480199.247,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,p2048,2048,ISQRT,50000,0.086571050,577560.282,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,p2048,2048,ISQRT,50000,0.090084203,555036.272,0
opencl-kernel,NVIDIA A40,GPU,w32-il,p2048,2048,ISQRT,50000,0.101472706,492743.341,0
opencl-e2e,NVIDIA A40,GPU,w32-il,p2048,2048,ISQRT,50000,0.104933227,476493.492,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,p2048,2048,ISQRT,50000,0.086837031,575791.220,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,p2048,2048,ISQRT,50000,0.090779219,550786.850,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-1t,p2048,2048,MODMUL_R2,6250,0.012093690,516798.420,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,gmp-nt,p2048,2048,MODMUL_R2,6250,0.000276404,22611793.643,0
library,AMD Eng Sample: 100-000000020-02_30/20_N,host-cpu,openssl-nt,p2048,2048,MODMUL_R2,6250,0.000470573,13281677.664,0
library,NVIDIA A40,gpu,cgbn,p2048,2048,MODMUL_R2,50000,0.000335872,148866234.756,0
opencl-kernel,NVIDIA A40,GPU,w8,p2048,2048,MODMUL_R2,50000,0.073018248,684760.334,0
opencl-e2e,NVIDIA A40,GPU,w8,p2048,2048,MODMUL_R2,50000,0.077228664,647428.003,0
opencl-kernel,NVIDIA A40,GPU,w16,p2048,2048,MODMUL_R2,50000,0.002775601,18014115.860,0
opencl-e2e,NVIDIA A40,GPU,w16,p2048,2048,MODMUL_R2,50000,0.006324300,7906013.800,0
opencl-kernel,NVIDIA A40,GPU,w32-opt,p2048,2048,MODMUL_R2,50000,0.001176384,42503118.199,0
opencl-e2e,NVIDIA A40,GPU,w32-opt,p2048,2048,MODMUL_R2,50000,0.004785337,10448583.857,0
opencl-kernel,NVIDIA A40,GPU,w32-o64,p2048,2048,MODMUL_R2,50000,0.000956236,52288326.186,0
opencl-e2e,NVIDIA A40,GPU,w32-o64,p2048,2048,MODMUL_R2,50000,0.004509826,11086901.775,0
opencl-kernel,NVIDIA A40,GPU,w32-il,p2048,2048,MODMUL_R2,50000,0.001156827,43221658.211,0
opencl-e2e,NVIDIA A40,GPU,w32-il,p2048,2048,MODMUL_R2,50000,0.004453408,11227357.674,0
opencl-kernel,NVIDIA A40,GPU,w32-il64,p2048,2048,MODMUL_R2,50000,0.000933885,53539803.201,0
opencl-e2e,NVIDIA A40,GPU,w32-il64,p2048,2048,MODMUL_R2,50000,0.004864049,10279501.403,0
```
