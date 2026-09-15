# MPA-OpenCL benchmark report - NVIDIA A100-SXM4-40GB

> **Note.** The multi-threaded GMP and OpenSSL baseline columns have been
> removed from this report: they predate the 2026-09-12 timing fix and were
> understated (see `reports/README.md`). The single-threaded GMP column, the
> OpenCL-on-CPU rows and all MPA measurements are unaffected and were verified
> against GMP before timing.


> **Partial report.** The run was interrupted or hit its time budget.
> Rows that never ran are marked `n/a`.

## 1. System under test

2 OpenCL device(s) exercised with the identical kernels and operands.

### Device 0 - NVIDIA A100-SXM4-40GB (GPU)

| Property | Value |
|---|---|
| Model | NVIDIA A100-SXM4-40GB |
| Type | GPU |
| Vendor | NVIDIA Corporation |
| Device memory | 39.49 GiB |
| Max single allocation | 9.87 GiB |
| Local memory | 48 KiB |
| Global cache | 3024 KiB |
| Compute units | 108 |
| Max clock | 1410 MHz |
| Max work-group size | 1024 |
| OpenCL version | OpenCL 3.0 CUDA |
| Driver | 580.159.03 |

### Device 1 - cpu-haswell-AMD EPYC 7532 32-Core Processor (CPU)

| Property | Value |
|---|---|
| Model | cpu-haswell-AMD EPYC 7532 32-Core Processor |
| Type | CPU |
| Vendor | AuthenticAMD |
| Device memory | 501.70 GiB |
| Max single allocation | 128.00 GiB |
| Local memory | 512 KiB |
| Global cache | 16384 KiB |
| Compute units | 128 |
| Max clock | 2400 MHz |
| Max work-group size | 4096 |
| OpenCL version | OpenCL 3.0 PoCL HSTR: cpu-x86_64-pc-linux-gnu-haswell |
| Driver | 5.0+debian |

### Host

| Property | Value |
|---|---|
| CPU | AMD EPYC 7532 32-Core Processor |
| Logical cores | 128 |
| OpenMP threads used | 128 |
| RAM | 503.7 GB |
| OS | Ubuntu 24.04.4 LTS |
| Kernel | 6.8.0-124-generic |
| Arch | x86_64 |
| GMP | 6.3.0 |
| OpenSSL | OpenSSL 3.0.13 30 Jan 2024 |
| CGBN | cgbn_results.tsv loaded |

## 2. Method

- Workload auto-sized from the device and host: --min-items from 700 x compute units, --items from ten times that capped by host RAM. Either flag, given explicitly, overrides its half.
- Base workload 50000 items, scaled down per operator by its cost weight and by modulus size. Device rows honour --min-items (75600) so the GPU is not left idle; the CPU libraries keep the smaller count because a full-width MODEXP there costs minutes. Both counts appear in every row as dev/cpu, and throughput is per-second so they remain comparable.
- 5 timed repetitions, **minimum** reported. Two untimed warm-up launches precede them.
- `kernel` times `clEnqueueNDRangeKernel` + `clFinish` only. `e2e` adds the host->device operand writes and the device->host result read.
- Every OpenCL device runs the same kernels on the same operands, so GPU and CPU-OpenCL columns are directly comparable.
- CPU library baselines (GMP, OpenSSL) run those same operands, with temporaries preallocated outside the timed region, so the figure is the arithmetic and not marshalling. The generator is reseeded per modulus and operation so every backend sees identical inputs.
- OpenSSL rows time the nearest BN primitive, which is not always semantically identical (its Montgomery routine expects Montgomery-domain inputs); they measure comparable work, not identical results. Correctness is judged against GMP only.
- Every device cell is checked word-for-word against GMP before it is timed. A cell that mismatches is reported and excluded from the speedup tables.
- Total wall time 3454.5 s.

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

### Device 0 - NVIDIA A100-SXM4-40GB (GPU)

#### secp256k1 (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | 935.83 M | 1.71 G | 2.38 G | 2.03 G | 2.24 G | 2.84 G | 3.08 G | 41.71 M | 4.44 G |
| SUBTRACT | 50000 / 50000 | 941.26 M | 1.70 G | 2.42 G | 2.05 G | 2.40 G | 3.19 G | 3.21 G | 53.61 M | 4.88 G |
| ADDMOD | 50000 / 50000 | 659.56 M | 1.32 G | 2.24 G | 2.51 G | 2.81 G | 3.17 G | 3.41 G | 13.61 M | 3.49 G |
| SUBTRACTMOD | 50000 / 50000 | 669.18 M | 1.28 G | 2.16 G | 2.42 G | 2.87 G | 3.32 G | 3.29 G | 14.87 M | 3.26 G |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | 31.61 M | 140.48 M | 489.86 M | 1.28 G | 1.67 G | 2.81 G | 2.94 G | 33.41 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | 163.19 M | 581.21 M | 1.57 G | 1.27 G | 1.58 G | 1.85 G | 2.02 G | 33.39 M | 4.44 G |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | 331.46 M | 1.16 G | 2.64 G | 2.00 G | 2.71 G | 2.38 G | 2.78 G | 4.06 M | 3.05 G |
| COMPARE | 50000 / 50000 | 967.70 M | 1.78 G | - | 2.66 G | 2.84 G | 3.15 G | 3.18 G | 101.37 M | 4.88 G |
| REDUCE | 50000 / 6250 | 189.87 M | 362.38 M | - | 791.52 M | 947.71 M | 958.79 M | 953.31 M | 33.62 M | 2.22 G |
| MODMUL | 50000 / 3125 | 73.42 M | 158.65 M | - | 286.82 M | 423.77 M | 363.38 M | 416.15 M | 7.00 M | 762.94 M |
| MODEXP | 50000 / 781 | 2.07 M | 10.81 M | - | 12.80 M | 22.47 M | 12.85 M | 22.56 M | 69.66 k | 2.48 M |
| EXPONENTIATION | 50000 / 781 | 1.24 M | 3.99 M | - | 64.79 M | 74.02 M | 64.93 M | 68.03 M | 208.99 k | n/a |
| DIVIDE | 50000 / 6250 | 113.25 M | 105.25 M | - | 283.24 M | 291.89 M | 292.57 M | 301.57 M | 15.14 M | 1.74 G |
| ISQRT | 50000 / 1562 | 9.98 M | 11.80 M | - | 42.87 M | 48.35 M | 42.57 M | 48.42 M | 8.40 M | n/a |
| MODMUL_R2 | 50000 / 50000 | 371.84 M | 1.08 G | - | 1.55 G | 1.84 G | 1.60 G | 1.97 G | 7.16 M | 2.22 G |

#### rsa256(composite) (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | 1.13 G | 1.56 G | 2.43 G | 2.09 G | 2.00 G | 2.79 G | 2.77 G | 41.61 M | 4.88 G |
| SUBTRACT | 50000 / 50000 | 1.16 G | 1.61 G | 2.48 G | 2.05 G | 2.05 G | 2.86 G | 2.81 G | 53.25 M | 4.88 G |
| ADDMOD | 50000 / 50000 | 886.85 M | 1.28 G | 2.25 G | 2.42 G | 2.42 G | 2.97 G | 2.87 G | 15.24 M | 3.76 G |
| SUBTRACTMOD | 50000 / 50000 | 840.07 M | 1.22 G | 2.19 G | 2.44 G | 2.55 G | 2.86 G | 2.92 G | 15.05 M | 3.76 G |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | 40.03 M | 127.47 M | 389.56 M | 1.28 G | 1.37 G | 2.43 G | 2.60 G | 33.36 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | 207.89 M | 531.42 M | 1.29 G | 1.31 G | 1.34 G | 1.60 G | 1.68 G | 33.22 M | 4.88 G |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | 420.74 M | 1.11 G | 2.31 G | 1.89 G | 2.22 G | 2.01 G | 2.47 G | 4.10 M | 3.05 G |
| COMPARE | 50000 / 50000 | 1.11 G | 1.64 G | - | 2.67 G | 2.64 G | 2.68 G | 2.78 G | 87.90 M | 4.88 G |
| REDUCE | 50000 / 6250 | 189.15 M | 330.89 M | - | 791.77 M | 765.48 M | 776.16 M | 769.12 M | 19.83 M | 2.44 G |
| MODMUL | 50000 / 3125 | 73.41 M | 144.18 M | - | 286.69 M | 333.74 M | 286.44 M | 325.02 M | 6.26 M | 827.60 M |
| MODEXP | 50000 / 781 | 2.07 M | 8.42 M | - | 12.83 M | 22.59 M | 12.91 M | 22.62 M | 74.08 k | 2.72 M |
| EXPONENTIATION | 50000 / 781 | 1.23 M | 3.94 M | - | 64.78 M | 74.14 M | 64.83 M | 67.58 M | 210.22 k | n/a |
| DIVIDE | 50000 / 6250 | 111.43 M | 130.32 M | - | 271.11 M | 278.45 M | 278.11 M | 289.36 M | 16.81 M | 1.88 G |
| ISQRT | 50000 / 1562 | 9.99 M | 12.66 M | - | 41.09 M | 44.21 M | 40.90 M | 44.28 M | 8.43 M | n/a |
| MODMUL_R2 | 50000 / 50000 | 373.25 M | 1.15 G | - | 1.50 G | 1.81 G | 1.54 G | 1.98 G | 7.18 M | 2.22 G |

#### brainpoolP512r1 (512-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 25000 | 602.34 M | 1.06 G | 1.37 G | 1.37 G | 1.34 G | 1.95 G | 2.00 G | 38.13 M | 4.44 G |
| SUBTRACT | 50000 / 25000 | 615.24 M | 1.06 G | 1.39 G | 1.34 G | 1.35 G | 1.94 G | 2.00 G | 48.74 M | 4.44 G |
| ADDMOD | 50000 / 25000 | 470.59 M | 817.94 M | 1.18 G | 1.04 G | 1.05 G | 2.11 G | 2.12 G | 14.37 M | 3.49 G |
| SUBTRACTMOD | 50000 / 25000 | 405.85 M | 752.23 M | 1.04 G | 1.03 G | 1.05 G | 2.11 G | 2.14 G | 14.01 M | 3.26 G |
| MULTIPLYOPERANDSCANNING | 50000 / 25000 | 5.64 M | 27.62 M | 82.05 M | 430.86 M | 424.78 M | 844.89 M | 853.39 M | 13.86 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 25000 | 24.44 M | 113.30 M | 305.52 M | 308.53 M | 310.80 M | 478.02 M | 518.90 M | 13.88 M | 3.26 G |
| MONTGOMERYMULTIPLICATION | 50000 / 25000 | 98.92 M | 403.69 M | 955.31 M | 797.97 M | 972.22 M | 930.07 M | 1.35 G | 1.76 M | 2.57 G |
| COMPARE | 50000 / 25000 | 466.90 M | 837.26 M | - | 1.36 G | 1.33 G | 2.50 G | 2.64 G | 100.26 M | 4.44 G |
| REDUCE | 50000 / 3125 | 67.13 M | 83.77 M | - | 315.80 M | 313.58 M | 244.34 M | 258.90 M | 19.66 M | 1.53 G |
| MODMUL | 50000 / 1562 | 23.04 M | 33.04 M | - | 86.00 M | 105.94 M | 68.41 M | 85.23 M | 3.62 M | 269.77 M |
| MODEXP | 50000 / 390 | 143.12 k | 1.54 M | - | 2.06 M | 3.22 M | 1.92 M | 3.17 M | 13.75 k | 1.12 M |
| EXPONENTIATION | 50000 / 390 | 172.31 k | 732.24 k | - | 2.11 M | 2.14 M | 2.50 M | 2.54 M | 60.39 k | n/a |
| DIVIDE | 50000 / 3125 | 32.73 M | 32.42 M | - | 97.47 M | 100.78 M | 104.61 M | 101.06 M | 13.97 M | 1.11 G |
| ISQRT | 50000 / 781 | 2.12 M | 2.09 M | - | 10.50 M | 11.23 M | 8.94 M | 9.38 M | 4.17 M | n/a |
| MODMUL_R2 | 50000 / 25000 | 85.04 M | 412.28 M | - | 660.60 M | 963.97 M | 696.68 M | 1.07 G | 3.58 M | 1.63 G |

#### p1024 (1024-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 12500 | 188.93 M | 367.95 M | 564.79 M | 774.12 M | 769.12 M | 1.69 G | 1.52 G | 32.41 M | 3.05 G |
| SUBTRACT | 50000 / 12500 | 189.39 M | 366.82 M | 564.41 M | 785.19 M | 771.50 M | 1.66 G | 1.67 G | 40.69 M | 3.05 G |
| ADDMOD | 50000 / 12500 | 106.31 M | 264.93 M | 416.85 M | 572.82 M | 568.90 M | 1.83 G | 1.84 G | 10.84 M | 2.71 G |
| SUBTRACTMOD | 50000 / 12500 | 107.53 M | 263.70 M | 418.00 M | 569.94 M | 573.99 M | 1.82 G | 1.78 G | 12.36 M | 2.57 G |
| MULTIPLYOPERANDSCANNING | 50000 / 12500 | 1.52 M | 4.73 M | 17.19 M | 184.33 M | 182.91 M | 305.57 M | 308.30 M | 3.95 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 12500 | 4.21 M | 16.37 M | 49.19 M | 49.03 M | 48.99 M | 132.47 M | 142.10 M | 3.95 M | 1.40 G |
| MONTGOMERYMULTIPLICATION | 50000 / 12500 | 23.79 M | 121.61 M | 325.25 M | 221.03 M | 283.26 M | 329.69 M | 490.49 M | 581.22 k | 957.41 M |
| COMPARE | 50000 / 12500 | 247.01 M | 476.93 M | - | 672.05 M | 652.41 M | 2.10 G | 2.07 G | 98.84 M | 3.05 G |
| REDUCE | 50000 / 1562 | 15.95 M | 29.50 M | - | 85.39 M | 86.18 M | 97.30 M | 102.95 M | 28.11 M | 1.16 G |
| MODMUL | 50000 / 781 | 3.95 M | 8.54 M | - | 19.68 M | 26.08 M | 19.94 M | 26.27 M | 1.28 M | 110.97 M |
| MODEXP | 50000 / 195 | 17.64 k | 113.98 k | - | 263.31 k | 416.73 k | 259.47 k | 425.93 k | 2.05 k | 242.91 k |
| EXPONENTIATION | 50000 / 195 | 22.71 k | 92.62 k | - | 334.93 k | 354.94 k | 331.16 k | 352.46 k | 14.22 k | n/a |
| DIVIDE | 50000 / 1562 | 2.76 M | 4.61 M | - | 22.75 M | 23.23 M | 23.03 M | 23.20 M | 13.08 M | 887.78 M |
| ISQRT | 50000 / 390 | 135.42 k | 295.66 k | - | 1.49 M | 1.60 M | 1.47 M | 1.61 M | 2.11 M | n/a |
| MODMUL_R2 | 50000 / 12500 | 15.39 M | 105.45 M | - | 185.31 M | 247.37 M | 223.89 M | 328.71 M | 1.30 M | 503.38 M |

#### p2048 (2048-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 6250 | 76.32 M | 193.25 M | 334.61 M | 422.66 M | 415.22 M | 511.99 M | 513.67 M | 23.55 M | 1.16 G |
| SUBTRACT | 50000 / 6250 | 76.54 M | 151.51 M | 335.55 M | 421.81 M | 413.23 M | 510.99 M | 513.88 M | 28.79 M | 1.16 G |
| ADDMOD | 50000 / 6250 | 60.31 M | 116.16 M | 273.08 M | 241.17 M | 243.54 M | 734.33 M | 731.97 M | 8.43 M | 1.16 G |
| SUBTRACTMOD | 50000 / 6250 | 56.63 M | 108.91 M | 251.88 M | 246.43 M | 244.48 M | 704.24 M | 738.89 M | 10.07 M | 1.16 G |
| MULTIPLYOPERANDSCANNING | 50000 / 6250 | 384.02 k | 1.55 M | 5.21 M | 64.99 M | 65.97 M | 91.05 M | 91.51 M | 1.16 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 6250 | 1.06 M | 4.16 M | 16.46 M | 12.74 M | 12.75 M | 28.68 M | 31.16 M | 1.16 M | 439.89 M |
| MONTGOMERYMULTIPLICATION | 50000 / 6250 | 727.16 k | 26.45 M | 115.48 M | 73.93 M | 95.06 M | 87.27 M | 105.02 M | 178.03 k | 239.35 M |
| COMPARE | 50000 / 6250 | 130.19 M | 203.32 M | - | 377.85 M | 364.68 M | 1.40 G | 1.44 G | 100.52 M | 1.16 G |
| REDUCE | 50000 / 781 | 140.89 k | 5.21 M | - | 26.08 M | 24.68 M | 31.44 M | 26.08 M | 18.21 M | 1.02 G |
| MODMUL | 50000 / 390 | 90.47 k | 2.46 M | - | 4.85 M | 6.54 M | 5.34 M | 6.52 M | 421.36 k | 47.68 M |
| MODEXP | 50000 / 97 | 871.5 | 6.54 k | - | 27.63 k | 22.29 k | 31.79 k | 22.66 k | 283.9 | 46.35 k |
| EXPONENTIATION | 50000 / 97 | 2.57 k | 11.68 k | - | 39.18 k | 43.05 k | 40.37 k | 43.24 k | 2.38 k | n/a |
| DIVIDE | 50000 / 781 | 43.54 k | 233.07 k | - | 2.40 M | 2.59 M | 2.40 M | 2.55 M | 10.27 M | 651.04 M |
| ISQRT | 50000 / 195 | 3.39 k | 9.83 k | - | 236.92 k | 702.37 k | 236.54 k | 700.63 k | 1.35 M | n/a |
| MODMUL_R2 | 50000 / 6250 | 1.56 M | 22.84 M | - | 53.11 M | 70.98 M | 59.50 M | 74.61 M | 414.50 k | 157.00 M |

### Device 1 - cpu-haswell-AMD EPYC 7532 32-Core Processor (CPU)

#### secp256k1 (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | - | - | - | - | - | - | - | 41.71 M | 4.44 G |
| SUBTRACT | 50000 / 50000 | - | - | - | - | - | - | - | 53.61 M | 4.88 G |
| ADDMOD | 50000 / 50000 | - | - | - | - | - | - | - | 13.61 M | 3.49 G |
| SUBTRACTMOD | 50000 / 50000 | - | - | - | - | - | - | - | 14.87 M | 3.26 G |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | - | - | - | - | - | - | - | 33.41 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | - | - | - | - | - | - | - | 33.39 M | 4.44 G |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | - | - | - | - | - | - | - | 4.06 M | 3.05 G |
| COMPARE | 50000 / 50000 | - | - | - | - | - | - | - | 101.37 M | 4.88 G |
| REDUCE | 50000 / 6250 | - | - | - | - | - | - | - | 33.62 M | 2.22 G |
| MODMUL | 50000 / 3125 | - | - | - | - | - | - | - | 7.00 M | 762.94 M |
| MODEXP | 50000 / 781 | - | - | - | - | - | - | - | 69.66 k | 2.48 M |
| EXPONENTIATION | 50000 / 781 | - | - | - | - | - | - | - | 208.99 k | n/a |
| DIVIDE | 50000 / 6250 | - | - | - | - | - | - | - | 15.14 M | 1.74 G |
| ISQRT | 50000 / 1562 | - | - | - | - | - | - | - | 8.40 M | n/a |
| MODMUL_R2 | 50000 / 50000 | - | - | - | - | - | - | - | 7.16 M | 2.22 G |

#### rsa256(composite) (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | - | - | - | - | - | - | - | 41.61 M | 4.88 G |
| SUBTRACT | 50000 / 50000 | - | - | - | - | - | - | - | 53.25 M | 4.88 G |
| ADDMOD | 50000 / 50000 | - | - | - | - | - | - | - | 15.24 M | 3.76 G |
| SUBTRACTMOD | 50000 / 50000 | - | - | - | - | - | - | - | 15.05 M | 3.76 G |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | - | - | - | - | - | - | - | 33.36 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | - | - | - | - | - | - | - | 33.22 M | 4.88 G |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | - | - | - | - | - | - | - | 4.10 M | 3.05 G |
| COMPARE | 50000 / 50000 | - | - | - | - | - | - | - | 87.90 M | 4.88 G |
| REDUCE | 50000 / 6250 | - | - | - | - | - | - | - | 19.83 M | 2.44 G |
| MODMUL | 50000 / 3125 | - | - | - | - | - | - | - | 6.26 M | 827.60 M |
| MODEXP | 50000 / 781 | - | - | - | - | - | - | - | 74.08 k | 2.72 M |
| EXPONENTIATION | 50000 / 781 | - | - | - | - | - | - | - | 210.22 k | n/a |
| DIVIDE | 50000 / 6250 | - | - | - | - | - | - | - | 16.81 M | 1.88 G |
| ISQRT | 50000 / 1562 | - | - | - | - | - | - | - | 8.43 M | n/a |
| MODMUL_R2 | 50000 / 50000 | - | - | - | - | - | - | - | 7.18 M | 2.22 G |

#### brainpoolP512r1 (512-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 25000 | - | - | - | - | - | - | - | 38.13 M | 4.44 G |
| SUBTRACT | 50000 / 25000 | - | - | - | - | - | - | - | 48.74 M | 4.44 G |
| ADDMOD | 50000 / 25000 | - | - | - | - | - | - | - | 14.37 M | 3.49 G |
| SUBTRACTMOD | 50000 / 25000 | - | - | - | - | - | - | - | 14.01 M | 3.26 G |
| MULTIPLYOPERANDSCANNING | 50000 / 25000 | - | - | - | - | - | - | - | 13.86 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 25000 | - | - | - | - | - | - | - | 13.88 M | 3.26 G |
| MONTGOMERYMULTIPLICATION | 50000 / 25000 | - | - | - | - | - | - | - | 1.76 M | 2.57 G |
| COMPARE | 50000 / 25000 | - | - | - | - | - | - | - | 100.26 M | 4.44 G |
| REDUCE | 50000 / 3125 | - | - | - | - | - | - | - | 19.66 M | 1.53 G |
| MODMUL | 50000 / 1562 | - | - | - | - | - | - | - | 3.62 M | 269.77 M |
| MODEXP | 50000 / 390 | - | - | - | - | - | - | - | 13.75 k | 1.12 M |
| EXPONENTIATION | 50000 / 390 | - | - | - | - | - | - | - | 60.39 k | n/a |
| DIVIDE | 50000 / 3125 | - | - | - | - | - | - | - | 13.97 M | 1.11 G |
| ISQRT | 50000 / 781 | - | - | - | - | - | - | - | 4.17 M | n/a |
| MODMUL_R2 | 50000 / 25000 | - | - | - | - | - | - | - | 3.58 M | 1.63 G |

#### p1024 (1024-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 12500 | - | - | - | - | - | - | - | 32.41 M | 3.05 G |
| SUBTRACT | 50000 / 12500 | - | - | - | - | - | - | - | 40.69 M | 3.05 G |
| ADDMOD | 50000 / 12500 | - | - | - | - | - | - | - | 10.84 M | 2.71 G |
| SUBTRACTMOD | 50000 / 12500 | - | - | - | - | - | - | - | 12.36 M | 2.57 G |
| MULTIPLYOPERANDSCANNING | 50000 / 12500 | - | - | - | - | - | - | - | 3.95 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 12500 | - | - | - | - | - | - | - | 3.95 M | 1.40 G |
| MONTGOMERYMULTIPLICATION | 50000 / 12500 | - | - | - | - | - | - | - | 581.22 k | 957.41 M |
| COMPARE | 50000 / 12500 | - | - | - | - | - | - | - | 98.84 M | 3.05 G |
| REDUCE | 50000 / 1562 | - | - | - | - | - | - | - | 28.11 M | 1.16 G |
| MODMUL | 50000 / 781 | - | - | - | - | - | - | - | 1.28 M | 110.97 M |
| MODEXP | 50000 / 195 | - | - | - | - | - | - | - | 2.05 k | 242.91 k |
| EXPONENTIATION | 50000 / 195 | - | - | - | - | - | - | - | 14.22 k | n/a |
| DIVIDE | 50000 / 1562 | - | - | - | - | - | - | - | 13.08 M | 887.78 M |
| ISQRT | 50000 / 390 | - | - | - | - | - | - | - | 2.11 M | n/a |
| MODMUL_R2 | 50000 / 12500 | - | - | - | - | - | - | - | 1.30 M | 503.38 M |

#### p2048 (2048-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 6250 | - | - | - | - | - | - | - | 23.55 M | 1.16 G |
| SUBTRACT | 50000 / 6250 | - | - | - | - | - | - | - | 28.79 M | 1.16 G |
| ADDMOD | 50000 / 6250 | - | - | - | - | - | - | - | 8.43 M | 1.16 G |
| SUBTRACTMOD | 50000 / 6250 | - | - | - | - | - | - | - | 10.07 M | 1.16 G |
| MULTIPLYOPERANDSCANNING | 50000 / 6250 | - | - | - | - | - | - | - | 1.16 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 6250 | - | - | - | - | - | - | - | 1.16 M | 439.89 M |
| MONTGOMERYMULTIPLICATION | 50000 / 6250 | - | - | - | - | - | - | - | 178.03 k | 239.35 M |
| COMPARE | 50000 / 6250 | - | - | - | - | - | - | - | 100.52 M | 1.16 G |
| REDUCE | 50000 / 781 | - | - | - | - | - | - | - | 18.21 M | 1.02 G |
| MODMUL | 50000 / 390 | - | - | - | - | - | - | - | 421.36 k | 47.68 M |
| MODEXP | 50000 / 97 | - | - | - | - | - | - | - | 283.9 | 46.35 k |
| EXPONENTIATION | 50000 / 97 | - | - | - | - | - | - | - | 2.38 k | n/a |
| DIVIDE | 50000 / 781 | - | - | - | - | - | - | - | 10.27 M | 651.04 M |
| ISQRT | 50000 / 195 | - | - | - | - | - | - | - | 1.35 M | n/a |
| MODMUL_R2 | 50000 / 6250 | - | - | - | - | - | - | - | 414.50 k | 157.00 M |

## 5. Head to head

Best OpenCL GPU result against best OpenCL CPU result and the CPU libraries.
Ratios above 1.00x mean the GPU is faster than that baseline.

### secp256k1 (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w32-il64 | 3.08 G | none | n/a | 41.71 M | 4.44 G | n/a | 0.69x |
| SUBTRACT | w32-il64 | 3.21 G | none | n/a | 53.61 M | 4.88 G | n/a | 0.66x |
| ADDMOD | w32-il64 | 3.41 G | none | n/a | 13.61 M | 3.49 G | n/a | 0.98x |
| SUBTRACTMOD | w32-il | 3.32 G | none | n/a | 14.87 M | 3.26 G | n/a | 1.02x |
| MULTIPLYOPERANDSCANNING | w32-il64 | 2.94 G | none | n/a | 33.41 M | n/a | n/a | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 2.02 G | none | n/a | 33.39 M | 4.44 G | n/a | 0.46x |
| MONTGOMERYMULTIPLICATION | w32-il64 | 2.78 G | none | n/a | 4.06 M | 3.05 G | n/a | 0.91x |
| COMPARE | w32-il64 | 3.18 G | none | n/a | 101.37 M | 4.88 G | n/a | 0.65x |
| REDUCE | w32-il | 119.85 M | none | n/a | 33.62 M | 2.22 G | n/a | 0.05x |
| MODMUL | w32-o64 | 26.49 M | none | n/a | 7.00 M | 762.94 M | n/a | 0.03x |
| MODEXP | w32-il64 | 352.35 k | none | n/a | 69.66 k | 2.48 M | n/a | 0.14x |
| EXPONENTIATION | w32-o64 | 1.16 M | none | n/a | 208.99 k | n/a | n/a | n/a |
| DIVIDE | w32-il64 | 37.70 M | none | n/a | 15.14 M | 1.74 G | n/a | 0.02x |
| ISQRT | w32-il64 | 1.51 M | none | n/a | 8.40 M | n/a | n/a | n/a |
| MODMUL_R2 | w32-il64 | 1.97 G | none | n/a | 7.16 M | 2.22 G | n/a | 0.89x |

### rsa256(composite) (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w32-il | 2.79 G | none | n/a | 41.61 M | 4.88 G | n/a | 0.57x |
| SUBTRACT | w32-il | 2.86 G | none | n/a | 53.25 M | 4.88 G | n/a | 0.58x |
| ADDMOD | w32-il | 2.97 G | none | n/a | 15.24 M | 3.76 G | n/a | 0.79x |
| SUBTRACTMOD | w32-il64 | 2.92 G | none | n/a | 15.05 M | 3.76 G | n/a | 0.78x |
| MULTIPLYOPERANDSCANNING | w32-il64 | 2.60 G | none | n/a | 33.36 M | n/a | n/a | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 1.68 G | none | n/a | 33.22 M | 4.88 G | n/a | 0.34x |
| MONTGOMERYMULTIPLICATION | w32-il64 | 2.47 G | none | n/a | 4.10 M | 3.05 G | n/a | 0.81x |
| COMPARE | w32-il64 | 2.78 G | none | n/a | 87.90 M | 4.88 G | n/a | 0.57x |
| REDUCE | w32-opt | 98.97 M | none | n/a | 19.83 M | 2.44 G | n/a | 0.04x |
| MODMUL | w32-o64 | 20.86 M | none | n/a | 6.26 M | 827.60 M | n/a | 0.03x |
| MODEXP | w32-il64 | 353.37 k | none | n/a | 74.08 k | 2.72 M | n/a | 0.13x |
| EXPONENTIATION | w32-o64 | 1.16 M | none | n/a | 210.22 k | n/a | n/a | n/a |
| DIVIDE | w32-il64 | 36.17 M | none | n/a | 16.81 M | 1.88 G | n/a | 0.02x |
| ISQRT | w32-il64 | 1.38 M | none | n/a | 8.43 M | n/a | n/a | n/a |
| MODMUL_R2 | w32-il64 | 1.98 G | none | n/a | 7.18 M | 2.22 G | n/a | 0.89x |

### brainpoolP512r1 (512-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w32-il64 | 1.00 G | none | n/a | 38.13 M | 4.44 G | n/a | 0.23x |
| SUBTRACT | w32-il64 | 997.60 M | none | n/a | 48.74 M | 4.44 G | n/a | 0.22x |
| ADDMOD | w32-il64 | 1.06 G | none | n/a | 14.37 M | 3.49 G | n/a | 0.30x |
| SUBTRACTMOD | w32-il64 | 1.07 G | none | n/a | 14.01 M | 3.26 G | n/a | 0.33x |
| MULTIPLYOPERANDSCANNING | w32-il64 | 426.69 M | none | n/a | 13.86 M | n/a | n/a | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 259.45 M | none | n/a | 13.88 M | 3.26 G | n/a | 0.08x |
| MONTGOMERYMULTIPLICATION | w32-il64 | 674.60 M | none | n/a | 1.76 M | 2.57 G | n/a | 0.26x |
| COMPARE | w32-il64 | 1.32 G | none | n/a | 100.26 M | 4.44 G | n/a | 0.30x |
| REDUCE | w32-opt | 19.74 M | none | n/a | 19.66 M | 1.53 G | n/a | 0.01x |
| MODMUL | w32-o64 | 3.31 M | none | n/a | 3.62 M | 269.77 M | n/a | 0.01x |
| MODEXP | w32-o64 | 25.13 k | none | n/a | 13.75 k | 1.12 M | n/a | 0.02x |
| EXPONENTIATION | w32-il64 | 19.84 k | none | n/a | 60.39 k | n/a | n/a | n/a |
| DIVIDE | w32-il | 6.54 M | none | n/a | 13.97 M | 1.11 G | n/a | 0.01x |
| ISQRT | w32-o64 | 175.47 k | none | n/a | 4.17 M | n/a | n/a | n/a |
| MODMUL_R2 | w32-il64 | 534.77 M | none | n/a | 3.58 M | 1.63 G | n/a | 0.33x |

### p1024 (1024-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w32-il | 421.87 M | none | n/a | 32.41 M | 3.05 G | n/a | 0.14x |
| SUBTRACT | w32-il64 | 416.66 M | none | n/a | 40.69 M | 3.05 G | n/a | 0.14x |
| ADDMOD | w32-il64 | 460.42 M | none | n/a | 10.84 M | 2.71 G | n/a | 0.17x |
| SUBTRACTMOD | w32-il | 455.22 M | none | n/a | 12.36 M | 2.57 G | n/a | 0.18x |
| MULTIPLYOPERANDSCANNING | w32-il64 | 77.08 M | none | n/a | 3.95 M | n/a | n/a | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 35.53 M | none | n/a | 3.95 M | 1.40 G | n/a | 0.03x |
| MONTGOMERYMULTIPLICATION | w32-il64 | 122.62 M | none | n/a | 581.22 k | 957.41 M | n/a | 0.13x |
| COMPARE | w32-il | 525.44 M | none | n/a | 98.84 M | 3.05 G | n/a | 0.17x |
| REDUCE | w32-il64 | 3.22 M | none | n/a | 28.11 M | 1.16 G | n/a | 0.00x |
| MODMUL | w32-il64 | 410.32 k | none | n/a | 1.28 M | 110.97 M | n/a | 0.00x |
| MODEXP | w32-il64 | 1.66 k | none | n/a | 2.05 k | 242.91 k | n/a | 0.01x |
| EXPONENTIATION | w32-o64 | 1.38 k | none | n/a | 14.22 k | n/a | n/a | n/a |
| DIVIDE | w32-o64 | 725.56 k | none | n/a | 13.08 M | 887.78 M | n/a | 0.00x |
| ISQRT | w32-il64 | 12.59 k | none | n/a | 2.11 M | n/a | n/a | n/a |
| MODMUL_R2 | w32-il64 | 82.18 M | none | n/a | 1.30 M | 503.38 M | n/a | 0.16x |

### p2048 (2048-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w32-il64 | 64.21 M | none | n/a | 23.55 M | 1.16 G | n/a | 0.06x |
| SUBTRACT | w32-il64 | 64.24 M | none | n/a | 28.79 M | 1.16 G | n/a | 0.06x |
| ADDMOD | w32-il | 91.79 M | none | n/a | 8.43 M | 1.16 G | n/a | 0.08x |
| SUBTRACTMOD | w32-il64 | 92.36 M | none | n/a | 10.07 M | 1.16 G | n/a | 0.08x |
| MULTIPLYOPERANDSCANNING | w32-il64 | 11.44 M | none | n/a | 1.16 M | n/a | n/a | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 3.89 M | none | n/a | 1.16 M | 439.89 M | n/a | 0.01x |
| MONTGOMERYMULTIPLICATION | w32 | 14.44 M | none | n/a | 178.03 k | 239.35 M | n/a | 0.06x |
| COMPARE | w32-il64 | 180.49 M | none | n/a | 100.52 M | 1.16 G | n/a | 0.16x |
| REDUCE | w32-il | 491.13 k | none | n/a | 18.21 M | 1.02 G | n/a | 0.00x |
| MODMUL | w32-o64 | 51.04 k | none | n/a | 421.36 k | 47.68 M | n/a | 0.00x |
| MODEXP | w32-il | 61.7 | none | n/a | 283.9 | 46.35 k | n/a | 0.00x |
| EXPONENTIATION | w32-il64 | 83.9 | none | n/a | 2.38 k | n/a | n/a | n/a |
| DIVIDE | w32-o64 | 40.51 k | none | n/a | 10.27 M | 651.04 M | n/a | 0.00x |
| ISQRT | w32-o64 | 2.74 k | none | n/a | 1.35 M | n/a | n/a | n/a |
| MODMUL_R2 | w32-il64 | 9.33 M | none | n/a | 414.50 k | 157.00 M | n/a | 0.06x |

## 6. Raw data

Also written to `NVIDIA_A100-SXM4-40GB_Report.csv` for analysis.

```csv
kind,device,device_type,kernel,modulus,bits,operation,items,seconds,ops_per_sec,mismatches
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,secp256k1,256,ADD,50000,0.001198843,41706891.021,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,secp256k1,256,ADD,50000,0.018179942,2750283.814,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,secp256k1,256,ADD,50000,0.008967528,5575672.508,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,secp256k1,256,ADD,50000,0.000011264,4438920454.545,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,secp256k1,256,ADD,50000,0.000053429,935828742.254,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,secp256k1,256,ADD,50000,0.000706230,70798435.994,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,secp256k1,256,ADD,50000,0.000029199,1712397652.462,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,secp256k1,256,ADD,50000,0.000576022,86802227.974,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32,secp256k1,256,ADD,50000,0.000021040,2376428798.442,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32,secp256k1,256,ADD,50000,0.000568652,87927229.192,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,secp256k1,256,ADD,50000,0.000024610,2031678001.892,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,secp256k1,256,ADD,50000,0.000582532,85832190.814,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,secp256k1,256,ADD,50000,0.000022319,2240229134.154,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,secp256k1,256,ADD,50000,0.000824378,60651797.792,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,secp256k1,256,ADD,50000,0.000017630,2836085113.576,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,secp256k1,256,ADD,50000,0.000684190,73079107.349,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,secp256k1,256,ADD,50000,0.000016240,3078829602.867,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,secp256k1,256,ADD,50000,0.000661420,75594912.661,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,secp256k1,256,SUBTRACT,50000,0.000932676,53609170.271,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,secp256k1,256,SUBTRACT,50000,0.009007187,5551122.617,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,secp256k1,256,SUBTRACT,50000,0.008999427,5555909.375,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,secp256k1,256,SUBTRACT,50000,0.000010240,4882812500.000,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,secp256k1,256,SUBTRACT,50000,0.000053120,941259543.283,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,secp256k1,256,SUBTRACT,50000,0.000604351,82733371.448,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,secp256k1,256,SUBTRACT,50000,0.000029439,1698447972.920,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,secp256k1,256,SUBTRACT,50000,0.000567032,88178441.182,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32,secp256k1,256,SUBTRACT,50000,0.000020689,2416704532.973,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32,secp256k1,256,SUBTRACT,50000,0.000594121,84157949.626,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,secp256k1,256,SUBTRACT,50000,0.000024400,2049164724.518,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,secp256k1,256,SUBTRACT,50000,0.000600291,83292943.143,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,secp256k1,256,SUBTRACT,50000,0.000020810,2402698256.842,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,secp256k1,256,SUBTRACT,50000,0.000824548,60639295.466,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,secp256k1,256,SUBTRACT,50000,0.000015690,3186744892.266,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,secp256k1,256,SUBTRACT,50000,0.000645441,77466414.059,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,secp256k1,256,SUBTRACT,50000,0.000015600,3205103799.887,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,secp256k1,256,SUBTRACT,50000,0.000654101,76440805.285,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,secp256k1,256,ADDMOD,50000,0.003672696,13613978.350,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,secp256k1,256,ADDMOD,50000,0.008990258,5561575.377,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,secp256k1,256,ADDMOD,50000,0.008998107,5556724.218,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,secp256k1,256,ADDMOD,50000,0.000014336,3487723214.286,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,secp256k1,256,ADDMOD,50000,0.000075808,659562780.412,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,secp256k1,256,ADDMOD,50000,0.000626281,79836379.643,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,secp256k1,256,ADDMOD,50000,0.000037869,1320342122.155,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,secp256k1,256,ADDMOD,50000,0.000573622,87165476.769,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32,secp256k1,256,ADDMOD,50000,0.000022280,2244162153.576,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32,secp256k1,256,ADDMOD,50000,0.000562692,88858622.393,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,secp256k1,256,ADDMOD,50000,0.000019900,2512558380.718,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,secp256k1,256,ADDMOD,50000,0.000599601,83388809.078,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,secp256k1,256,ADDMOD,50000,0.000017810,2807461758.092,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,secp256k1,256,ADDMOD,50000,0.000793978,62974028.638,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,secp256k1,256,ADDMOD,50000,0.000015750,3174590734.116,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,secp256k1,256,ADDMOD,50000,0.000690589,72401992.680,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,secp256k1,256,ADDMOD,50000,0.000014679,3406217124.005,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,secp256k1,256,ADDMOD,50000,0.000653061,76562516.418,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,secp256k1,256,SUBTRACTMOD,50000,0.003361580,14873958.247,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,secp256k1,256,SUBTRACTMOD,50000,0.008940629,5592447.568,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,secp256k1,256,SUBTRACTMOD,50000,0.008999588,5555809.908,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,secp256k1,256,SUBTRACTMOD,50000,0.000015360,3255208333.333,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.000074719,669177302.330,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.000619851,80664539.936,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,secp256k1,256,SUBTRACTMOD,50000,0.000039050,1280413341.442,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,secp256k1,256,SUBTRACTMOD,50000,0.000574702,87001621.666,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32,secp256k1,256,SUBTRACTMOD,50000,0.000023200,2155199261.356,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32,secp256k1,256,SUBTRACTMOD,50000,0.000583151,85741102.142,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,secp256k1,256,SUBTRACTMOD,50000,0.000020700,2415508467.561,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,secp256k1,256,SUBTRACTMOD,50000,0.000565971,88343710.374,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,secp256k1,256,SUBTRACTMOD,50000,0.000017440,2866981266.688,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,secp256k1,256,SUBTRACTMOD,50000,0.000791678,63157011.379,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,secp256k1,256,SUBTRACTMOD,50000,0.000015080,3315655335.968,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,secp256k1,256,SUBTRACTMOD,50000,0.000651000,76804907.540,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,secp256k1,256,SUBTRACTMOD,50000,0.000015180,3293787613.117,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,secp256k1,256,SUBTRACTMOD,50000,0.000655241,76307819.310,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001496708,33406649.335,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.008996448,5557749.005,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.008991998,5560499.650,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001581536,31614841.702,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.002527082,19785666.831,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000355925,140479055.084,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001027136,48679040.221,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000102069,489863600.861,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000780149,64090329.643,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000039050,1280413341.442,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000712469,70178516.158,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000029869,1673955200.798,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000920677,54307857.161,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000017820,2805847768.370,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000836638,59762984.682,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000017020,2937734128.591,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000833847,59963032.000,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001497288,33393703.929,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.008999887,5555625.070,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.008980567,5567577.129,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000011264,4438920454.545,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000306385,163193312.633,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000982415,50894991.942,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000086028,581208394.364,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000783779,63793515.024,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000031809,1571889244.463,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000724299,69032275.784,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000039370,1269993996.239,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000714230,70005471.650,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000031720,1576296755.630,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000929256,53806487.796,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000027019,1850545170.019,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000837618,59693080.665,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000024750,2020172384.339,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000846328,59078767.191,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.012311459,4061257.060,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.009009367,5549779.558,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.008998397,5556545.070,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000016384,3051757812.500,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000150848,331459086.756,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000705400,70881767.585,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000042939,1164440060.296,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000587942,85042390.757,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000018910,2644098165.431,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000546702,91457472.612,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000024980,2001606561.778,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000582391,85852985.308,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000018460,2708596498.663,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000793169,63038247.937,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000021000,2380960649.267,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000690000,72463754.282,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000018010,2776248381.425,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000667361,74921960.134,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,secp256k1,256,COMPARE,50000,0.000493233,101371949.283,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,secp256k1,256,COMPARE,50000,0.009003087,5553650.686,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,secp256k1,256,COMPARE,50000,0.006079170,8224807.111,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,secp256k1,256,COMPARE,50000,0.000010240,4882812500.000,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,secp256k1,256,COMPARE,50000,0.000051669,967701133.762,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,secp256k1,256,COMPARE,50000,0.000629251,79459561.639,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,secp256k1,256,COMPARE,50000,0.000028030,1783802079.942,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,secp256k1,256,COMPARE,50000,0.000564002,88652173.202,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,secp256k1,256,COMPARE,50000,0.000018780,2662389843.789,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,secp256k1,256,COMPARE,50000,0.000560122,89266329.025,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,secp256k1,256,COMPARE,50000,0.000017610,2839309897.665,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,secp256k1,256,COMPARE,50000,0.000797218,62718088.239,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,secp256k1,256,COMPARE,50000,0.000015850,3154538527.528,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,secp256k1,256,COMPARE,50000,0.000645931,77407607.523,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,secp256k1,256,COMPARE,50000,0.000015729,3178914124.996,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,secp256k1,256,COMPARE,50000,0.000606831,82395241.106,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,secp256k1,256,REDUCE,6250,0.000185887,33622601.825,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,secp256k1,256,REDUCE,6250,0.008963938,697238.174,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,secp256k1,256,REDUCE,6250,0.008998818,694535.678,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,secp256k1,256,REDUCE,50000,0.000022528,2219460227.273,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,secp256k1,256,REDUCE,50000,0.000263336,189871412.353,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,secp256k1,256,REDUCE,50000,0.000844078,59236254.684,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,secp256k1,256,REDUCE,50000,0.000137978,362377220.769,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,secp256k1,256,REDUCE,50000,0.000701720,71253502.082,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,secp256k1,256,REDUCE,50000,0.000063169,791523957.097,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,secp256k1,256,REDUCE,50000,0.000623981,80130644.655,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,secp256k1,256,REDUCE,50000,0.000052759,947706355.749,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,secp256k1,256,REDUCE,50000,0.000837198,59723028.986,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,secp256k1,256,REDUCE,50000,0.000052149,958792224.236,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,secp256k1,256,REDUCE,50000,0.000713779,70049680.949,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,secp256k1,256,REDUCE,50000,0.000052449,953310152.442,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,secp256k1,256,REDUCE,50000,0.000700160,71412255.800,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,secp256k1,256,MODMUL,3125,0.000446414,7000234.075,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,secp256k1,256,MODMUL,3125,0.008032971,389021.715,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,secp256k1,256,MODMUL,3125,0.009000177,347215.410,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,secp256k1,256,MODMUL,50000,0.000065536,762939453.125,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,secp256k1,256,MODMUL,50000,0.000681009,73420503.468,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,secp256k1,256,MODMUL,50000,0.001271192,39333161.799,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,secp256k1,256,MODMUL,50000,0.000315155,158651906.423,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,secp256k1,256,MODMUL,50000,0.000890738,56133237.351,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,secp256k1,256,MODMUL,50000,0.000174327,286817595.709,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,secp256k1,256,MODMUL,50000,0.000721560,69294275.523,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,secp256k1,256,MODMUL,50000,0.000117988,423772411.861,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,secp256k1,256,MODMUL,50000,0.000893417,55964890.214,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,secp256k1,256,MODMUL,50000,0.000137598,363376704.457,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,secp256k1,256,MODMUL,50000,0.000786879,63542162.965,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,secp256k1,256,MODMUL,50000,0.000120148,416153193.600,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,secp256k1,256,MODMUL,50000,0.000760799,65720357.154,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,secp256k1,256,MODEXP,781,0.011211044,69663.452,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,secp256k1,256,MODEXP,781,0.008995248,86823.621,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,secp256k1,256,MODEXP,781,0.007140695,109373.105,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,secp256k1,256,MODEXP,50000,0.020189185,2476573.472,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,secp256k1,256,MODEXP,50000,0.024204622,2065721.165,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,secp256k1,256,MODEXP,50000,0.024854103,2011740.266,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,secp256k1,256,MODEXP,50000,0.004625695,10809185.829,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,secp256k1,256,MODEXP,50000,0.005292396,9447516.668,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,secp256k1,256,MODEXP,50000,0.003905313,12803070.805,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,secp256k1,256,MODEXP,50000,0.004535163,11024962.492,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,secp256k1,256,MODEXP,50000,0.002225147,22470424.594,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,secp256k1,256,MODEXP,50000,0.003097834,16140309.307,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,secp256k1,256,MODEXP,50000,0.003890832,12850722.052,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,secp256k1,256,MODEXP,50000,0.004675720,10693540.505,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,secp256k1,256,MODEXP,50000,0.002216518,22557904.644,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,secp256k1,256,MODEXP,50000,0.002937147,17023322.708,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,secp256k1,256,EXPONENTIATION,781,0.003736985,208992.013,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,secp256k1,256,EXPONENTIATION,781,0.009002107,86757.469,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,secp256k1,256,EXPONENTIATION,781,0.009000207,86775.779,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,secp256k1,256,EXPONENTIATION,50000,0.040388513,1237975.761,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,secp256k1,256,EXPONENTIATION,50000,0.041187892,1213949.000,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,secp256k1,256,EXPONENTIATION,50000,0.012527453,3991234.361,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,secp256k1,256,EXPONENTIATION,50000,0.013194373,3789494.239,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,secp256k1,256,EXPONENTIATION,50000,0.000771769,64786236.915,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,secp256k1,256,EXPONENTIATION,50000,0.001353000,36954919.783,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,secp256k1,256,EXPONENTIATION,50000,0.000675490,74020326.981,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,secp256k1,256,EXPONENTIATION,50000,0.001483808,33697075.935,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,secp256k1,256,EXPONENTIATION,50000,0.000770048,64930997.190,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,secp256k1,256,EXPONENTIATION,50000,0.001444569,34612403.375,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,secp256k1,256,EXPONENTIATION,50000,0.000734949,68031929.727,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,secp256k1,256,EXPONENTIATION,50000,0.001394020,35867492.998,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,secp256k1,256,DIVIDE,6250,0.000412854,15138527.200,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,secp256k1,256,DIVIDE,6250,0.005340191,1170370.108,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,secp256k1,256,DIVIDE,6250,0.006141880,1017603.720,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,secp256k1,256,DIVIDE,50000,0.000028672,1743861607.143,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,secp256k1,256,DIVIDE,50000,0.000441514,113246584.014,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,secp256k1,256,DIVIDE,50000,0.001136953,43977191.238,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,secp256k1,256,DIVIDE,50000,0.000475073,105246931.422,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,secp256k1,256,DIVIDE,50000,0.001178274,42434953.943,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,secp256k1,256,DIVIDE,50000,0.000176527,283242674.939,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,secp256k1,256,DIVIDE,50000,0.000844707,59192104.940,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,secp256k1,256,DIVIDE,50000,0.000171297,291890301.149,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,secp256k1,256,DIVIDE,50000,0.001099094,45492021.899,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,secp256k1,256,DIVIDE,50000,0.000170898,292571906.889,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,secp256k1,256,DIVIDE,50000,0.000995285,50236872.999,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,secp256k1,256,DIVIDE,50000,0.000165798,301571363.492,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,secp256k1,256,DIVIDE,50000,0.000995335,50234334.675,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,secp256k1,256,ISQRT,1562,0.000185968,8399299.530,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,secp256k1,256,ISQRT,1562,0.006391545,244385.363,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,secp256k1,256,ISQRT,50000,0.005009596,9980844.292,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,secp256k1,256,ISQRT,50000,0.005653576,8843960.009,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,secp256k1,256,ISQRT,50000,0.004237470,11799492.961,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,secp256k1,256,ISQRT,50000,0.004823322,10366299.534,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,secp256k1,256,ISQRT,50000,0.001166333,42869404.287,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,secp256k1,256,ISQRT,50000,0.001755084,28488658.814,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,secp256k1,256,ISQRT,50000,0.001034044,48353832.692,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,secp256k1,256,ISQRT,50000,0.001869272,26748381.546,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,secp256k1,256,ISQRT,50000,0.001174502,42571231.738,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,secp256k1,256,ISQRT,50000,0.001844163,27112570.077,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,secp256k1,256,ISQRT,50000,0.001032735,48415120.332,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,secp256k1,256,ISQRT,50000,0.001725855,28971150.693,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,secp256k1,256,MODMUL_R2,50000,0.006979547,7163788.735,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,secp256k1,256,MODMUL_R2,50000,0.006413805,7795684.668,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,secp256k1,256,MODMUL_R2,50000,0.008967888,5575448.710,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,secp256k1,256,MODMUL_R2,50000,0.000022528,2219460227.273,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,secp256k1,256,MODMUL_R2,50000,0.000134468,371835460.993,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,secp256k1,256,MODMUL_R2,50000,0.000673740,74212636.737,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,secp256k1,256,MODMUL_R2,50000,0.000046339,1079007380.014,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,secp256k1,256,MODMUL_R2,50000,0.000631071,79230367.432,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,secp256k1,256,MODMUL_R2,50000,0.000032240,1550865637.322,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,secp256k1,256,MODMUL_R2,50000,0.000589642,84797252.982,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,secp256k1,256,MODMUL_R2,50000,0.000027130,1842985571.823,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,secp256k1,256,MODMUL_R2,50000,0.000771078,64844298.459,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,secp256k1,256,MODMUL_R2,50000,0.000031170,1604108077.745,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,secp256k1,256,MODMUL_R2,50000,0.000708319,70589685.511,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,secp256k1,256,MODMUL_R2,50000,0.000025319,1974806562.201,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,secp256k1,256,MODMUL_R2,50000,0.000681770,73338562.341,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,ADD,50000,0.001201722,41606966.231,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,ADD,50000,0.009001527,5554613.134,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,ADD,50000,0.007999782,6250170.329,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,rsa256(composite),256,ADD,50000,0.000010240,4882812500.000,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,rsa256(composite),256,ADD,50000,0.000044210,1130968847.693,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,rsa256(composite),256,ADD,50000,0.000573982,87110742.753,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,rsa256(composite),256,ADD,50000,0.000031990,1562988477.103,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,rsa256(composite),256,ADD,50000,0.000516392,96825700.149,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32,rsa256(composite),256,ADD,50000,0.000020560,2431866065.726,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32,rsa256(composite),256,ADD,50000,0.000511362,97778049.105,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,rsa256(composite),256,ADD,50000,0.000023949,2087773330.741,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,rsa256(composite),256,ADD,50000,0.000606511,82438701.148,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,rsa256(composite),256,ADD,50000,0.000025009,1999295840.316,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,rsa256(composite),256,ADD,50000,0.000570652,87619062.789,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,rsa256(composite),256,ADD,50000,0.000017910,2791695242.057,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,rsa256(composite),256,ADD,50000,0.000659601,75803368.201,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,rsa256(composite),256,ADD,50000,0.000018020,2774670070.805,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,rsa256(composite),256,ADD,50000,0.000682820,73225729.552,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACT,50000,0.000938896,53254029.433,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACT,50000,0.006423995,7783318.758,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACT,50000,0.008966588,5576257.132,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,rsa256(composite),256,SUBTRACT,50000,0.000010240,4882812500.000,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,rsa256(composite),256,SUBTRACT,50000,0.000043020,1162246927.532,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,rsa256(composite),256,SUBTRACT,50000,0.000581831,85935645.137,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,rsa256(composite),256,SUBTRACT,50000,0.000030999,1612975745.467,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,rsa256(composite),256,SUBTRACT,50000,0.000541373,92357737.810,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32,rsa256(composite),256,SUBTRACT,50000,0.000020149,2481549894.844,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32,rsa256(composite),256,SUBTRACT,50000,0.000529842,94367757.489,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,rsa256(composite),256,SUBTRACT,50000,0.000024389,2050103721.241,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,rsa256(composite),256,SUBTRACT,50000,0.000600211,83304057.976,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,rsa256(composite),256,SUBTRACT,50000,0.000024340,2054261271.499,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,rsa256(composite),256,SUBTRACT,50000,0.000572102,87396980.579,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,rsa256(composite),256,SUBTRACT,50000,0.000017510,2855544449.763,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,rsa256(composite),256,SUBTRACT,50000,0.000652990,76570869.973,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,rsa256(composite),256,SUBTRACT,50000,0.000017790,2810548172.966,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,rsa256(composite),256,SUBTRACT,50000,0.000674510,74127895.428,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,ADDMOD,50000,0.003280782,15240268.282,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,ADDMOD,50000,0.006409516,7800900.923,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,ADDMOD,50000,0.009000557,5555211.744,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,rsa256(composite),256,ADDMOD,50000,0.000013312,3756009615.385,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,rsa256(composite),256,ADDMOD,50000,0.000056379,886854892.503,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,rsa256(composite),256,ADDMOD,50000,0.000621701,80424496.571,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,rsa256(composite),256,ADDMOD,50000,0.000038939,1284057621.889,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,rsa256(composite),256,ADDMOD,50000,0.000561342,89072241.684,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32,rsa256(composite),256,ADDMOD,50000,0.000022250,2247168021.431,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32,rsa256(composite),256,ADDMOD,50000,0.000539112,92745123.146,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,rsa256(composite),256,ADDMOD,50000,0.000020650,2421282244.171,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,rsa256(composite),256,ADDMOD,50000,0.000601831,83079816.345,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,rsa256(composite),256,ADDMOD,50000,0.000020619,2424945965.356,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,rsa256(composite),256,ADDMOD,50000,0.000561632,89026232.074,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,rsa256(composite),256,ADDMOD,50000,0.000016860,2965564183.721,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,rsa256(composite),256,ADDMOD,50000,0.000660131,75742570.223,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,rsa256(composite),256,ADDMOD,50000,0.000017449,2865527538.630,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,rsa256(composite),256,ADDMOD,50000,0.000667310,74927658.862,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACTMOD,50000,0.003321951,15051396.955,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACTMOD,50000,0.008986377,5563977.466,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACTMOD,50000,0.006314407,7918400.394,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,rsa256(composite),256,SUBTRACTMOD,50000,0.000013312,3756009615.385,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,rsa256(composite),256,SUBTRACTMOD,50000,0.000059519,840068398.323,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,rsa256(composite),256,SUBTRACTMOD,50000,0.000634170,78843253.206,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,rsa256(composite),256,SUBTRACTMOD,50000,0.000040820,1224893707.506,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,rsa256(composite),256,SUBTRACTMOD,50000,0.000561492,89048529.351,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32,rsa256(composite),256,SUBTRACTMOD,50000,0.000022859,2187336926.806,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32,rsa256(composite),256,SUBTRACTMOD,50000,0.000555662,89982688.330,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,rsa256(composite),256,SUBTRACTMOD,50000,0.000020500,2439047370.693,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,rsa256(composite),256,SUBTRACTMOD,50000,0.000536982,93112995.942,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,rsa256(composite),256,SUBTRACTMOD,50000,0.000019610,2549728875.380,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,rsa256(composite),256,SUBTRACTMOD,50000,0.000546682,91460822.442,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,rsa256(composite),256,SUBTRACTMOD,50000,0.000017459,2863846115.275,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,rsa256(composite),256,SUBTRACTMOD,50000,0.000668430,74802174.104,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,rsa256(composite),256,SUBTRACTMOD,50000,0.000017109,2922462164.884,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,rsa256(composite),256,SUBTRACTMOD,50000,0.000674690,74108146.885,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001498758,33360959.287,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.008986898,5563654.857,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.011000278,4545339.697,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001248992,40032280.487,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001936371,25821498.416,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000392234,127474940.996,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001029075,48587317.742,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000128349,389563367.885,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000785598,63645777.889,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000038989,1282401347.203,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000732199,68287462.931,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000036469,1371020115.684,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000685180,72973517.564,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000020590,2428346165.502,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000839058,59590647.236,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000019200,2604146837.408,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000848177,58949967.443,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001504947,33223762.125,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.008960318,5580159.209,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.006326167,7903680.173,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000010240,4882812500.000,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000240517,207885643.923,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000919366,54385317.094,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000094088,531418557.599,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000718410,69598128.047,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000038699,1292030351.964,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000683930,73106871.588,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000038269,1306541364.289,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000735150,68013313.545,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000037440,1335483170.607,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000725879,68882016.252,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000031159,1604683430.723,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000846598,59059887.264,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000029799,1677905120.873,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000848298,58941553.906,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.012195280,4099946.925,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.012124571,4123857.242,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.008996548,5557687.156,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000016384,3051757812.500,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000118839,420736985.561,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000680800,73443001.572,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000045199,1106220456.606,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000557283,89720956.218,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000021670,2307335877.600,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000674301,74150880.425,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000026449,1890423817.321,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000589701,84788748.970,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000022529,2219346073.872,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000596591,83809472.904,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000024920,2006431512.660,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000710499,70373075.106,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000020280,2465481444.743,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000670370,74585656.522,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,COMPARE,50000,0.000568811,87902611.270,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,COMPARE,50000,0.008994778,5558780.789,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,COMPARE,50000,0.009001428,5554674.053,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,rsa256(composite),256,COMPARE,50000,0.000010240,4882812500.000,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,rsa256(composite),256,COMPARE,50000,0.000044979,1111626039.424,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,rsa256(composite),256,COMPARE,50000,0.000579911,86220084.635,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,rsa256(composite),256,COMPARE,50000,0.000030500,1639350551.162,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,rsa256(composite),256,COMPARE,50000,0.000522203,95748239.646,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,rsa256(composite),256,COMPARE,50000,0.000018750,2666687753.632,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,rsa256(composite),256,COMPARE,50000,0.000593271,84278501.774,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,rsa256(composite),256,COMPARE,50000,0.000018960,2637149582.474,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,rsa256(composite),256,COMPARE,50000,0.000568952,87880884.095,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,rsa256(composite),256,COMPARE,50000,0.000018630,2683817796.441,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,rsa256(composite),256,COMPARE,50000,0.000661960,75533279.495,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,rsa256(composite),256,COMPARE,50000,0.000017969,2782579620.607,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,rsa256(composite),256,COMPARE,50000,0.000689980,72465906.109,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,REDUCE,6250,0.000315135,19832748.379,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,REDUCE,6250,0.008976468,696264.954,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,REDUCE,6250,0.008998437,694565.079,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,rsa256(composite),256,REDUCE,50000,0.000020480,2441406250.000,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,rsa256(composite),256,REDUCE,50000,0.000264336,189152944.907,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,rsa256(composite),256,REDUCE,50000,0.000804098,62181486.425,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,rsa256(composite),256,REDUCE,50000,0.000151108,330889123.642,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,rsa256(composite),256,REDUCE,50000,0.000666331,75037777.677,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,rsa256(composite),256,REDUCE,50000,0.000063149,791774934.372,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,rsa256(composite),256,REDUCE,50000,0.000613871,81450333.238,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,rsa256(composite),256,REDUCE,50000,0.000065319,765476701.528,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,rsa256(composite),256,REDUCE,50000,0.000644611,77566136.624,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,rsa256(composite),256,REDUCE,50000,0.000064420,776155893.047,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,rsa256(composite),256,REDUCE,50000,0.000723629,69096155.811,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,rsa256(composite),256,REDUCE,50000,0.000065009,769122977.522,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,rsa256(composite),256,REDUCE,50000,0.000730000,68493153.439,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MODMUL,3125,0.000499413,6257341.755,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MODMUL,3125,0.008041941,388587.775,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MODMUL,3125,0.008001572,390548.235,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,rsa256(composite),256,MODMUL,50000,0.000060416,827595338.983,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,rsa256(composite),256,MODMUL,50000,0.000681130,73407402.667,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,rsa256(composite),256,MODMUL,50000,0.001261691,39629366.147,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,rsa256(composite),256,MODMUL,50000,0.000346795,144177316.825,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,rsa256(composite),256,MODMUL,50000,0.000920367,54326156.947,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,rsa256(composite),256,MODMUL,50000,0.000174407,286685879.061,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,rsa256(composite),256,MODMUL,50000,0.000745399,67078175.735,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,rsa256(composite),256,MODMUL,50000,0.000149818,333739000.594,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,rsa256(composite),256,MODMUL,50000,0.000724079,69053274.742,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,rsa256(composite),256,MODMUL,50000,0.000174557,286438854.076,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,rsa256(composite),256,MODMUL,50000,0.000808728,61825489.512,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,rsa256(composite),256,MODMUL,50000,0.000153838,325017805.814,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,rsa256(composite),256,MODMUL,50000,0.000820098,60968324.847,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MODEXP,781,0.010543184,74076.295,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MODEXP,781,0.008997778,86799.208,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MODEXP,781,0.009996243,78129.355,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,rsa256(composite),256,MODEXP,50000,0.018352127,2724479.838,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,rsa256(composite),256,MODEXP,50000,0.024195053,2066538.135,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,rsa256(composite),256,MODEXP,50000,0.024842282,2012697.522,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,rsa256(composite),256,MODEXP,50000,0.005940066,8417414.629,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,rsa256(composite),256,MODEXP,50000,0.006560067,7621873.705,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,rsa256(composite),256,MODEXP,50000,0.003898163,12826554.916,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,rsa256(composite),256,MODEXP,50000,0.004527634,11043295.437,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,rsa256(composite),256,MODEXP,50000,0.002213677,22586850.367,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,rsa256(composite),256,MODEXP,50000,0.002886738,17320589.971,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,rsa256(composite),256,MODEXP,50000,0.003872213,12912512.973,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,rsa256(composite),256,MODEXP,50000,0.004632281,10793818.098,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,rsa256(composite),256,MODEXP,50000,0.002210168,22622712.934,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,rsa256(composite),256,MODEXP,50000,0.002967387,16849841.449,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,EXPONENTIATION,781,0.003715126,210221.681,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,EXPONENTIATION,781,0.006130989,127385.643,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,EXPONENTIATION,781,0.006003181,130097.694,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,rsa256(composite),256,EXPONENTIATION,50000,0.040556661,1232843.098,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,rsa256(composite),256,EXPONENTIATION,50000,0.040994424,1219678.060,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,rsa256(composite),256,EXPONENTIATION,50000,0.012706040,3935136.235,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,rsa256(composite),256,EXPONENTIATION,50000,0.013332942,3750110.056,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,rsa256(composite),256,EXPONENTIATION,50000,0.000771838,64780452.112,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,rsa256(composite),256,EXPONENTIATION,50000,0.001351320,37000866.111,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,rsa256(composite),256,EXPONENTIATION,50000,0.000674400,74139974.853,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,rsa256(composite),256,EXPONENTIATION,50000,0.001262211,39613035.306,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,rsa256(composite),256,EXPONENTIATION,50000,0.000771258,64829147.044,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,rsa256(composite),256,EXPONENTIATION,50000,0.001450228,34477336.425,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,rsa256(composite),256,EXPONENTIATION,50000,0.000739909,67575855.396,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,rsa256(composite),256,EXPONENTIATION,50000,0.001364700,36638088.137,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,DIVIDE,6250,0.000371905,16805364.986,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,DIVIDE,6250,0.006001182,1041461.516,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,DIVIDE,6250,0.007003677,892388.374,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,rsa256(composite),256,DIVIDE,50000,0.000026624,1878004807.692,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,rsa256(composite),256,DIVIDE,50000,0.000448694,111434522.081,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,rsa256(composite),256,DIVIDE,50000,0.001156423,43236765.080,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,rsa256(composite),256,DIVIDE,50000,0.000383675,130318595.619,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,rsa256(composite),256,DIVIDE,50000,0.001069485,46751485.021,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,rsa256(composite),256,DIVIDE,50000,0.000184427,271109955.713,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,rsa256(composite),256,DIVIDE,50000,0.000877387,56987402.650,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,rsa256(composite),256,DIVIDE,50000,0.000179567,278447740.634,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,rsa256(composite),256,DIVIDE,50000,0.000877157,57002317.484,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,rsa256(composite),256,DIVIDE,50000,0.000179787,278106613.484,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,rsa256(composite),256,DIVIDE,50000,0.001023925,48831705.855,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,rsa256(composite),256,DIVIDE,50000,0.000172797,289356666.370,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,rsa256(composite),256,DIVIDE,50000,0.000998756,50062281.636,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,ISQRT,1562,0.000185207,8433806.582,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,ISQRT,1562,0.008998037,173593.421,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,rsa256(composite),256,ISQRT,50000,0.005006986,9986047.087,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,rsa256(composite),256,ISQRT,50000,0.005693816,8781456.759,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,rsa256(composite),256,ISQRT,50000,0.003950344,12657125.771,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,rsa256(composite),256,ISQRT,50000,0.004572695,10934470.520,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,rsa256(composite),256,ISQRT,50000,0.001216752,41093008.643,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,rsa256(composite),256,ISQRT,50000,0.001814003,27563348.669,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,rsa256(composite),256,ISQRT,50000,0.001131053,44206590.095,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,rsa256(composite),256,ISQRT,50000,0.001727404,28945167.367,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,rsa256(composite),256,ISQRT,50000,0.001222481,40900430.510,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,rsa256(composite),256,ISQRT,50000,0.001941091,25758711.045,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,rsa256(composite),256,ISQRT,50000,0.001129243,44277447.340,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,rsa256(composite),256,ISQRT,50000,0.001794644,27860683.431,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MODMUL_R2,50000,0.006965967,7177753.987,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MODMUL_R2,50000,0.009103836,5492190.309,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MODMUL_R2,50000,0.012991268,3848739.097,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,rsa256(composite),256,MODMUL_R2,50000,0.000022528,2219460227.273,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,rsa256(composite),256,MODMUL_R2,50000,0.000133958,373250817.421,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,rsa256(composite),256,MODMUL_R2,50000,0.000715340,69896830.126,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,rsa256(composite),256,MODMUL_R2,50000,0.000043360,1153135181.228,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,rsa256(composite),256,MODMUL_R2,50000,0.000607192,82346269.101,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,rsa256(composite),256,MODMUL_R2,50000,0.000033230,1504661963.818,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,rsa256(composite),256,MODMUL_R2,50000,0.000598872,83490348.792,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,rsa256(composite),256,MODMUL_R2,50000,0.000027609,1810999871.816,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,rsa256(composite),256,MODMUL_R2,50000,0.000503723,99260894.654,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,rsa256(composite),256,MODMUL_R2,50000,0.000032520,1537519078.985,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,rsa256(composite),256,MODMUL_R2,50000,0.000710100,70412624.334,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,rsa256(composite),256,MODMUL_R2,50000,0.000025210,1983342243.895,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,rsa256(composite),256,MODMUL_R2,50000,0.000666450,75024355.519,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,ADD,25000,0.000655670,38128926.075,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,ADD,25000,0.010262609,2436027.721,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,ADD,25000,0.006004021,4163876.150,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,brainpoolP512r1,512,ADD,50000,0.000011264,4438920454.545,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,brainpoolP512r1,512,ADD,50000,0.000083009,602342533.701,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,brainpoolP512r1,512,ADD,50000,0.001092814,45753444.009,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,brainpoolP512r1,512,ADD,50000,0.000047329,1056437379.720,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,brainpoolP512r1,512,ADD,50000,0.001064985,46949015.451,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32,brainpoolP512r1,512,ADD,50000,0.000036619,1365423617.081,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32,brainpoolP512r1,512,ADD,50000,0.001329500,37608116.900,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,brainpoolP512r1,512,ADD,50000,0.000036619,1365423617.081,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,brainpoolP512r1,512,ADD,50000,0.001162923,42995110.192,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,brainpoolP512r1,512,ADD,50000,0.000037379,1337646004.161,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,brainpoolP512r1,512,ADD,50000,0.000939876,53198489.672,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,brainpoolP512r1,512,ADD,50000,0.000025650,1949316166.512,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,brainpoolP512r1,512,ADD,50000,0.001113273,44912624.109,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,brainpoolP512r1,512,ADD,50000,0.000024960,2003212298.278,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,brainpoolP512r1,512,ADD,50000,0.001082465,46190887.539,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACT,25000,0.000512903,48742151.217,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACT,25000,0.008986557,2781933.089,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACT,25000,0.009003817,2776600.159,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,brainpoolP512r1,512,SUBTRACT,50000,0.000011264,4438920454.545,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,brainpoolP512r1,512,SUBTRACT,50000,0.000081269,615243735.210,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,brainpoolP512r1,512,SUBTRACT,50000,0.001059874,47175420.705,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,brainpoolP512r1,512,SUBTRACT,50000,0.000046989,1064079977.802,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,brainpoolP512r1,512,SUBTRACT,50000,0.001053475,47461978.146,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32,brainpoolP512r1,512,SUBTRACT,50000,0.000036089,1385455444.446,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32,brainpoolP512r1,512,SUBTRACT,50000,0.001420579,35196918.185,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,brainpoolP512r1,512,SUBTRACT,50000,0.000037219,1343403135.361,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,brainpoolP512r1,512,SUBTRACT,50000,0.001159863,43108536.539,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,brainpoolP512r1,512,SUBTRACT,50000,0.000037070,1348803276.094,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,brainpoolP512r1,512,SUBTRACT,50000,0.000919146,54398322.062,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,brainpoolP512r1,512,SUBTRACT,50000,0.000025709,1944832139.105,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,brainpoolP512r1,512,SUBTRACT,50000,0.001103693,45302465.570,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,brainpoolP512r1,512,SUBTRACT,50000,0.000025060,1995209276.052,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,brainpoolP512r1,512,SUBTRACT,50000,0.001109884,45049748.055,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,ADDMOD,25000,0.001740315,14365215.511,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,ADDMOD,25000,0.006002711,4164784.786,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,ADDMOD,25000,0.008994687,2779418.453,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,brainpoolP512r1,512,ADDMOD,50000,0.000014336,3487723214.286,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,brainpoolP512r1,512,ADDMOD,50000,0.000106249,470592644.017,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,brainpoolP512r1,512,ADDMOD,50000,0.001102163,45365341.080,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,brainpoolP512r1,512,ADDMOD,50000,0.000061129,817939442.692,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,brainpoolP512r1,512,ADDMOD,50000,0.001087505,45976784.448,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32,brainpoolP512r1,512,ADDMOD,50000,0.000042510,1176187779.603,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32,brainpoolP512r1,512,ADDMOD,50000,0.001430199,34960169.490,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,brainpoolP512r1,512,ADDMOD,50000,0.000047869,1044516259.071,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,brainpoolP512r1,512,ADDMOD,50000,0.001141213,43813035.105,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,brainpoolP512r1,512,ADDMOD,50000,0.000047420,1054404053.656,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,brainpoolP512r1,512,ADDMOD,50000,0.000991205,50443664.462,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,brainpoolP512r1,512,ADDMOD,50000,0.000023740,2106160773.622,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,brainpoolP512r1,512,ADDMOD,50000,0.001072874,46603802.552,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,brainpoolP512r1,512,ADDMOD,50000,0.000023600,2118669739.542,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,brainpoolP512r1,512,ADDMOD,50000,0.001109894,45049351.138,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACTMOD,25000,0.001785054,14005179.562,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACTMOD,25000,0.008994907,2779350.538,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000154498,161814845.335,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000015360,3255208333.333,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000123198,405850269.498,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,brainpoolP512r1,512,SUBTRACTMOD,50000,0.001168113,42804087.697,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000066469,752230840.473,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,brainpoolP512r1,512,SUBTRACTMOD,50000,0.001087384,45981903.612,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000047959,1042558887.670,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32,brainpoolP512r1,512,SUBTRACTMOD,50000,0.001385220,36095353.102,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000048650,1027750011.007,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000907866,55074189.743,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000047589,1050658849.086,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000947116,52791821.102,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000023710,2108849524.707,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,brainpoolP512r1,512,SUBTRACTMOD,50000,0.001105644,45222501.625,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000023340,2142256542.037,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,brainpoolP512r1,512,SUBTRACTMOD,50000,0.001075734,46479895.382,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.001803334,13863212.132,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.000093129,268444851.570,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.000101559,246161529.961,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.008867219,5638746.431,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.010494985,4764180.156,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.001810244,27620588.717,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.003112476,16064381.655,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000609401,82047824.140,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.002241997,22301545.865,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000116048,430856512.754,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.001213543,41201667.496,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000117709,424776611.889,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.001308610,38208488.594,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000059179,844893870.293,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.001351340,37000305.102,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000058590,853388404.161,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.001451919,34437186.222,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.001800744,13883151.576,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000086719,288286891.338,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000122698,203752262.687,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000015360,3255208333.333,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.002046060,24437213.215,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.003432470,14566770.092,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000441314,113297966.583,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.001725855,28971142.876,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000163658,305515054.360,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.001735834,28804598.684,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000162057,308533202.304,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.001504728,33228594.479,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000160877,310797100.845,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.001355820,36878055.201,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000104598,478021664.841,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.001530977,32658884.562,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000096358,518898850.319,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.001529778,32684483.405,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.014240160,1755598.237,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000377794,66173583.120,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.005750115,4347739.110,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000019456,2569901315.789,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000505473,98917251.945,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.001556617,32120944.511,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000123858,403688139.468,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.001199743,41675583.686,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000052339,955311817.931,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.001391899,35922138.744,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000062659,797971019.404,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.001201652,41609384.746,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000051429,972222364.680,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.001004486,49776707.708,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000053759,930073387.788,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.001090384,45855401.348,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000037059,1349193084.037,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.001124314,44471560.930,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,COMPARE,25000,0.000249346,100262371.841,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,COMPARE,25000,0.008994767,2779393.848,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,COMPARE,25000,0.008996797,2778766.630,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,brainpoolP512r1,512,COMPARE,50000,0.000011264,4438920454.545,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,brainpoolP512r1,512,COMPARE,50000,0.000107089,466901111.440,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,brainpoolP512r1,512,COMPARE,50000,0.001204003,41528131.786,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,brainpoolP512r1,512,COMPARE,50000,0.000059719,837258235.409,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,brainpoolP512r1,512,COMPARE,50000,0.001138633,43912300.854,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,brainpoolP512r1,512,COMPARE,50000,0.000036889,1355409464.901,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,brainpoolP512r1,512,COMPARE,50000,0.001155473,43272328.752,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,brainpoolP512r1,512,COMPARE,50000,0.000037630,1328740392.778,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,brainpoolP512r1,512,COMPARE,50000,0.000990686,50470077.961,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,brainpoolP512r1,512,COMPARE,50000,0.000019969,2503886910.897,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,brainpoolP512r1,512,COMPARE,50000,0.001028265,48625603.506,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,brainpoolP512r1,512,COMPARE,50000,0.000018950,2638575278.911,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,brainpoolP512r1,512,COMPARE,50000,0.001054134,47432310.937,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,REDUCE,3125,0.000158928,19663008.248,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,REDUCE,3125,0.008994767,347424.231,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,REDUCE,3125,0.008962048,348692.630,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,brainpoolP512r1,512,REDUCE,50000,0.000032768,1525878906.250,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,brainpoolP512r1,512,REDUCE,50000,0.000744809,67131311.032,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,brainpoolP512r1,512,REDUCE,50000,0.001812053,27593013.229,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,brainpoolP512r1,512,REDUCE,50000,0.000596862,83771483.229,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,brainpoolP512r1,512,REDUCE,50000,0.001681576,29734011.122,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,brainpoolP512r1,512,REDUCE,50000,0.000158328,315799917.061,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,brainpoolP512r1,512,REDUCE,50000,0.001276301,39175719.224,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,brainpoolP512r1,512,REDUCE,50000,0.000159447,313583664.029,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,brainpoolP512r1,512,REDUCE,50000,0.001109554,45063152.854,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,brainpoolP512r1,512,REDUCE,50000,0.000204637,244335430.785,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,brainpoolP512r1,512,REDUCE,50000,0.001284851,38915025.391,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,brainpoolP512r1,512,REDUCE,50000,0.000193127,258897041.740,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,brainpoolP512r1,512,REDUCE,50000,0.001280681,39041735.761,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL,1562,0.000431473,3620157.826,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL,1562,0.009012587,173313.169,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL,1562,0.005112795,305508.035,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,brainpoolP512r1,512,MODMUL,50000,0.000185344,269768646.409,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,brainpoolP512r1,512,MODMUL,50000,0.002169818,23043405.716,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,brainpoolP512r1,512,MODMUL,50000,0.003218212,15536578.905,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,brainpoolP512r1,512,MODMUL,50000,0.001513329,33039734.880,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,brainpoolP512r1,512,MODMUL,50000,0.002531225,19753282.593,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,brainpoolP512r1,512,MODMUL,50000,0.000581392,86000482.487,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,brainpoolP512r1,512,MODMUL,50000,0.001440839,34702005.638,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,brainpoolP512r1,512,MODMUL,50000,0.000471973,105938294.753,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,brainpoolP512r1,512,MODMUL,50000,0.001419879,35214267.503,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,brainpoolP512r1,512,MODMUL,50000,0.000730939,68405185.009,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,brainpoolP512r1,512,MODMUL,50000,0.001858702,26900493.946,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,brainpoolP512r1,512,MODMUL,50000,0.000586661,85228089.487,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,brainpoolP512r1,512,MODMUL,50000,0.001667806,29979506.956,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MODEXP,390,0.028372601,13745.656,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MODEXP,390,0.006965557,55989.779,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MODEXP,390,0.007865984,49580.574,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,brainpoolP512r1,512,MODEXP,50000,0.044518400,1123131.110,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,brainpoolP512r1,512,MODEXP,50000,0.349355672,143120.619,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,brainpoolP512r1,512,MODEXP,50000,0.350430696,142681.565,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,brainpoolP512r1,512,MODEXP,50000,0.032432582,1541659.546,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,brainpoolP512r1,512,MODEXP,50000,0.033464846,1494105.178,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,brainpoolP512r1,512,MODEXP,50000,0.024230296,2063532.382,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,brainpoolP512r1,512,MODEXP,50000,0.025319489,1974763.415,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,brainpoolP512r1,512,MODEXP,50000,0.015517601,3222147.490,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,brainpoolP512r1,512,MODEXP,50000,0.013505851,3702099.209,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,brainpoolP512r1,512,MODEXP,50000,0.026059783,1918665.243,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,brainpoolP512r1,512,MODEXP,50000,0.027236936,1835742.454,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,brainpoolP512r1,512,MODEXP,50000,0.015793600,3165839.303,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,brainpoolP512r1,512,MODEXP,50000,0.013713241,3646111.060,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,EXPONENTIATION,390,0.006457645,60393.534,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,EXPONENTIATION,390,0.006209758,62804.383,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,EXPONENTIATION,390,0.011996803,32508.661,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,brainpoolP512r1,512,EXPONENTIATION,50000,0.290177526,172308.313,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,brainpoolP512r1,512,EXPONENTIATION,50000,0.290976285,171835.310,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,brainpoolP512r1,512,EXPONENTIATION,50000,0.068283654,732239.668,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,brainpoolP512r1,512,EXPONENTIATION,50000,0.069548776,718919.914,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,brainpoolP512r1,512,EXPONENTIATION,50000,0.023688453,2110733.042,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,brainpoolP512r1,512,EXPONENTIATION,50000,0.024839796,2012899.006,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,brainpoolP512r1,512,EXPONENTIATION,50000,0.023311245,2144887.605,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,brainpoolP512r1,512,EXPONENTIATION,50000,0.024312301,2056572.097,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,brainpoolP512r1,512,EXPONENTIATION,50000,0.019970024,2503752.624,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,brainpoolP512r1,512,EXPONENTIATION,50000,0.021186816,2359958.188,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,brainpoolP512r1,512,EXPONENTIATION,50000,0.019659574,2543290.085,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,brainpoolP512r1,512,EXPONENTIATION,50000,0.020869406,2395851.619,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,DIVIDE,3125,0.000223636,13973618.960,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,DIVIDE,3125,0.004621431,676197.515,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,DIVIDE,3125,0.009003387,347091.607,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,brainpoolP512r1,512,DIVIDE,50000,0.000045056,1109730113.636,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,brainpoolP512r1,512,DIVIDE,50000,0.001527687,32729215.959,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,brainpoolP512r1,512,DIVIDE,50000,0.002904257,17216108.404,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,brainpoolP512r1,512,DIVIDE,50000,0.001542438,32416217.561,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,brainpoolP512r1,512,DIVIDE,50000,0.002774621,18020479.880,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,brainpoolP512r1,512,DIVIDE,50000,0.000512992,97467400.843,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,brainpoolP512r1,512,DIVIDE,50000,0.001758215,28437932.505,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,brainpoolP512r1,512,DIVIDE,50000,0.000496122,100781651.452,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,brainpoolP512r1,512,DIVIDE,50000,0.001646326,30370657.041,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,brainpoolP512r1,512,DIVIDE,50000,0.000477983,104606197.800,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,brainpoolP512r1,512,DIVIDE,50000,0.001871193,26720916.633,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,brainpoolP512r1,512,DIVIDE,50000,0.000494763,101058434.479,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,brainpoolP512r1,512,DIVIDE,50000,0.001972551,25347887.180,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,ISQRT,781,0.000187127,4173637.879,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,ISQRT,781,0.005081155,153705.215,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,brainpoolP512r1,512,ISQRT,50000,0.023561971,2122063.576,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,brainpoolP512r1,512,ISQRT,50000,0.024672445,2026552.301,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,brainpoolP512r1,512,ISQRT,50000,0.023904952,2091616.855,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,brainpoolP512r1,512,ISQRT,50000,0.024959367,2003255.951,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,brainpoolP512r1,512,ISQRT,50000,0.004760821,10502390.494,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,brainpoolP512r1,512,ISQRT,50000,0.005826035,8582166.015,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,brainpoolP512r1,512,ISQRT,50000,0.004450875,11233746.505,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,brainpoolP512r1,512,ISQRT,50000,0.005427310,9212667.426,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,brainpoolP512r1,512,ISQRT,50000,0.005595557,8935661.036,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,brainpoolP512r1,512,ISQRT,50000,0.006756149,7400665.931,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,brainpoolP512r1,512,ISQRT,50000,0.005329402,9381915.357,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,brainpoolP512r1,512,ISQRT,50000,0.006537245,7648481.750,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL_R2,25000,0.006978627,3582366.409,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL_R2,25000,0.008991317,2780460.337,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL_R2,25000,0.006186049,4041351.780,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,brainpoolP512r1,512,MODMUL_R2,50000,0.000030720,1627604166.667,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,brainpoolP512r1,512,MODMUL_R2,50000,0.000587941,85042525.467,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,brainpoolP512r1,512,MODMUL_R2,50000,0.001626276,30745089.811,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,brainpoolP512r1,512,MODMUL_R2,50000,0.000121278,412275171.151,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,brainpoolP512r1,512,MODMUL_R2,50000,0.001017996,49116123.705,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,brainpoolP512r1,512,MODMUL_R2,50000,0.000075689,660597525.547,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,brainpoolP512r1,512,MODMUL_R2,50000,0.001097284,45567063.303,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,brainpoolP512r1,512,MODMUL_R2,50000,0.000051869,963965439.724,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,brainpoolP512r1,512,MODMUL_R2,50000,0.000927087,53932375.644,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,brainpoolP512r1,512,MODMUL_R2,50000,0.000071769,696678512.617,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,brainpoolP512r1,512,MODMUL_R2,50000,0.001125613,44420231.372,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,brainpoolP512r1,512,MODMUL_R2,50000,0.000046749,1069538537.547,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,brainpoolP512r1,512,MODMUL_R2,50000,0.001225072,40813929.031,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,p1024,1024,ADD,12500,0.000385625,32414886.659,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,p1024,1024,ADD,12500,0.008996007,1389505.361,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,p1024,1024,ADD,12500,0.008012172,1560126.209,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,p1024,1024,ADD,50000,0.000016384,3051757812.500,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,p1024,1024,ADD,50000,0.000264646,188931615.296,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,p1024,1024,ADD,50000,0.002210297,22621392.722,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,p1024,1024,ADD,50000,0.000135888,367950401.793,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,p1024,1024,ADD,50000,0.001968592,25398864.535,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32,p1024,1024,ADD,50000,0.000088528,564794344.386,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32,p1024,1024,ADD,50000,0.001864823,26812195.440,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p1024,1024,ADD,50000,0.000064589,774124627.985,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p1024,1024,ADD,50000,0.001877282,26634253.635,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p1024,1024,ADD,50000,0.000065009,769122977.522,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p1024,1024,ADD,50000,0.001824523,27404425.118,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,p1024,1024,ADD,50000,0.000029630,1687477328.304,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,p1024,1024,ADD,50000,0.002187408,22858103.773,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p1024,1024,ADD,50000,0.000032829,1523059651.910,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p1024,1024,ADD,50000,0.002346326,21309910.558,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,p1024,1024,SUBTRACT,12500,0.000307175,40693433.669,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,p1024,1024,SUBTRACT,12500,0.008985967,1391057.871,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,p1024,1024,SUBTRACT,12500,0.006110700,2045592.053,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,p1024,1024,SUBTRACT,50000,0.000016384,3051757812.500,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,p1024,1024,SUBTRACT,50000,0.000264006,189389490.765,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,p1024,1024,SUBTRACT,50000,0.002239866,22322762.135,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,p1024,1024,SUBTRACT,50000,0.000136308,366816579.723,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,p1024,1024,SUBTRACT,50000,0.001658247,30152326.191,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32,p1024,1024,SUBTRACT,50000,0.000088588,564408397.725,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32,p1024,1024,SUBTRACT,50000,0.001906952,26219847.532,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p1024,1024,SUBTRACT,50000,0.000063679,785191719.135,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p1024,1024,SUBTRACT,50000,0.001834413,27256675.023,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p1024,1024,SUBTRACT,50000,0.000064809,771499269.989,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p1024,1024,SUBTRACT,50000,0.001790283,27928548.773,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,p1024,1024,SUBTRACT,50000,0.000030209,1655144395.974,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,p1024,1024,SUBTRACT,50000,0.002188007,22851847.688,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p1024,1024,SUBTRACT,50000,0.000030000,1666653975.941,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p1024,1024,SUBTRACT,50000,0.002096939,23844275.520,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,p1024,1024,ADDMOD,12500,0.001153213,10839277.095,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,p1024,1024,ADDMOD,12500,0.006989526,1788390.204,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,p1024,1024,ADDMOD,12500,0.004090409,3055928.931,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,p1024,1024,ADDMOD,50000,0.000018432,2712673611.111,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,p1024,1024,ADDMOD,50000,0.000470313,106312124.278,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,p1024,1024,ADDMOD,50000,0.002483503,20132852.979,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,p1024,1024,ADDMOD,50000,0.000188727,264933041.195,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,p1024,1024,ADDMOD,50000,0.001960012,25510046.200,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32,p1024,1024,ADDMOD,50000,0.000119948,416847897.013,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32,p1024,1024,ADDMOD,50000,0.001942581,25738952.039,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p1024,1024,ADDMOD,50000,0.000087288,572818112.660,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p1024,1024,ADDMOD,50000,0.001846643,27076163.691,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p1024,1024,ADDMOD,50000,0.000087889,568899980.926,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p1024,1024,ADDMOD,50000,0.001805753,27689280.047,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,p1024,1024,ADDMOD,50000,0.000027360,1827490126.798,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,p1024,1024,ADDMOD,50000,0.002209957,22624877.082,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p1024,1024,ADDMOD,50000,0.000027149,1841689520.085,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p1024,1024,ADDMOD,50000,0.002111169,23683560.605,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,p1024,1024,SUBTRACTMOD,12500,0.001011185,12361735.443,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,p1024,1024,SUBTRACTMOD,12500,0.000070199,178065456.282,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,p1024,1024,SUBTRACTMOD,12500,0.000131428,95109271.220,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,p1024,1024,SUBTRACTMOD,50000,0.000019456,2569901315.789,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,p1024,1024,SUBTRACTMOD,50000,0.000464993,107528485.777,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,p1024,1024,SUBTRACTMOD,50000,0.002505043,19959736.773,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,p1024,1024,SUBTRACTMOD,50000,0.000189607,263703300.277,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,p1024,1024,SUBTRACTMOD,50000,0.001937023,25812807.914,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32,p1024,1024,SUBTRACTMOD,50000,0.000119618,417996801.595,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32,p1024,1024,SUBTRACTMOD,50000,0.001809513,27631748.123,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p1024,1024,SUBTRACTMOD,50000,0.000087729,569938758.785,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p1024,1024,SUBTRACTMOD,50000,0.001862542,26845035.499,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p1024,1024,SUBTRACTMOD,50000,0.000087109,573993972.149,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p1024,1024,SUBTRACTMOD,50000,0.001867632,26771870.609,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,p1024,1024,SUBTRACTMOD,50000,0.000027459,1820888997.422,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,p1024,1024,SUBTRACTMOD,50000,0.002088769,23937540.190,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p1024,1024,SUBTRACTMOD,50000,0.000028140,1776835717.359,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p1024,1024,SUBTRACTMOD,50000,0.002223707,22484978.643,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.003164013,3950679.195,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.000110058,113576360.282,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.000280126,44622777.218,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.032942392,1517801.139,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.036287533,1377883.689,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.010571310,4729782.856,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.012806979,3904121.336,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.002907878,17194667.603,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.005100985,9802028.321,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.000271246,184334626.161,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.002566462,19482073.616,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.000273356,182911517.679,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.002457133,20348918.604,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.000163627,305573307.986,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.003120284,16024182.096,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.000162177,308304641.474,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.003046055,16414672.314,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.003162194,3952951.599,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.000130048,96118340.280,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.008072931,1548384.376,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000035840,1395089285.714,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.011888214,4205846.255,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.014902070,3355238.530,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.003054767,16367860.914,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.005157587,9694455.731,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.001016375,49194456.642,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.003254902,15361447.402,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.001019755,49031389.187,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.003246152,15402853.111,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.001020705,48985756.681,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.003254472,15363476.128,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000377434,132473548.147,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.003188922,15679279.362,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000351855,142103965.177,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.003247253,15397631.522,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.021506542,581218.500,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.006224409,2008222.720,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.008995257,1389621.242,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000052224,957414215.686,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.002102019,23786656.908,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.004331046,11544554.622,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000411144,121611895.008,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.002204559,22680274.949,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000153728,325250152.668,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.001919482,26048689.380,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000226216,221027763.563,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.002047650,24418234.938,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000176518,283257619.061,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.001988791,25140903.080,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000151658,329689245.063,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.002314366,21604185.057,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000101938,490494643.461,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.002319866,21552970.308,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,p1024,1024,COMPARE,12500,0.000126468,98839218.228,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,p1024,1024,COMPARE,12500,0.006259947,1996821.991,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,p1024,1024,COMPARE,12500,0.006326457,1975829.435,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,p1024,1024,COMPARE,50000,0.000016384,3051757812.500,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,p1024,1024,COMPARE,50000,0.000202417,247014369.847,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,p1024,1024,COMPARE,50000,0.001969611,25385725.867,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,p1024,1024,COMPARE,50000,0.000104838,476926073.786,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,p1024,1024,COMPARE,50000,0.001895163,26382958.506,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p1024,1024,COMPARE,50000,0.000074399,672050512.296,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p1024,1024,COMPARE,50000,0.001961811,25486655.193,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p1024,1024,COMPARE,50000,0.000076639,652409344.943,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p1024,1024,COMPARE,50000,0.001796133,27837583.967,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,p1024,1024,COMPARE,50000,0.000023790,2101749577.200,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,p1024,1024,COMPARE,50000,0.002213047,22593285.445,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p1024,1024,COMPARE,50000,0.000024160,2069544598.423,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p1024,1024,COMPARE,50000,0.002230897,22412513.281,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,p1024,1024,REDUCE,1562,0.000055570,28108613.144,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,p1024,1024,REDUCE,1562,0.009000108,173553.471,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,p1024,1024,REDUCE,1562,0.004314256,362055.482,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,p1024,1024,REDUCE,50000,0.000043008,1162574404.762,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,p1024,1024,REDUCE,50000,0.003134583,15951084.415,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,p1024,1024,REDUCE,50000,0.005971601,8372963.641,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,p1024,1024,REDUCE,50000,0.001694656,29504513.407,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,p1024,1024,REDUCE,50000,0.003461471,14444726.676,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p1024,1024,REDUCE,50000,0.000585532,85392456.299,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p1024,1024,REDUCE,50000,0.002438024,20508414.368,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p1024,1024,REDUCE,50000,0.000580211,86175521.230,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p1024,1024,REDUCE,50000,0.002365895,21133650.598,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,p1024,1024,REDUCE,50000,0.000513852,97304261.665,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,p1024,1024,REDUCE,50000,0.002701750,18506523.377,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p1024,1024,REDUCE,50000,0.000485652,102954352.066,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p1024,1024,REDUCE,50000,0.002690761,18582104.523,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,p1024,1024,MODMUL,781,0.000612431,1275246.242,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,p1024,1024,MODMUL,781,0.006000401,130157.968,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,p1024,1024,MODMUL,781,0.000091608,8525442.512,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,p1024,1024,MODMUL,50000,0.000450560,110973011.364,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,p1024,1024,MODMUL,50000,0.012661112,3949100.203,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,p1024,1024,MODMUL,50000,0.014723352,3395965.774,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,p1024,1024,MODMUL,50000,0.005854547,8540370.382,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,p1024,1024,MODMUL,50000,0.007627322,6555380.884,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p1024,1024,MODMUL,50000,0.002540353,19682302.291,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p1024,1024,MODMUL,50000,0.004432965,11279132.411,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p1024,1024,MODMUL,50000,0.001917532,26075188.034,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p1024,1024,MODMUL,50000,0.003721755,13434522.162,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,p1024,1024,MODMUL,50000,0.002507683,19938721.487,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,p1024,1024,MODMUL,50000,0.004753000,10519671.558,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p1024,1024,MODMUL,50000,0.001903412,26268617.503,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p1024,1024,MODMUL,50000,0.004099461,12196726.016,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,p1024,1024,MODEXP,195,0.095139435,2049.623,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,p1024,1024,MODEXP,195,0.006309427,30906.135,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,p1024,1024,MODEXP,195,0.001583636,123134.372,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,p1024,1024,MODEXP,50000,0.205836281,242911.501,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,p1024,1024,MODEXP,50000,2.834289914,17641.103,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,p1024,1024,MODEXP,50000,2.837901030,17618.655,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,p1024,1024,MODEXP,50000,0.438683223,113977.461,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,p1024,1024,MODEXP,50000,0.440496516,113508.276,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p1024,1024,MODEXP,50000,0.189887638,263313.613,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p1024,1024,MODEXP,50000,0.191853209,260615.917,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p1024,1024,MODEXP,50000,0.119981537,416730.785,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p1024,1024,MODEXP,50000,0.122262984,408954.520,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,p1024,1024,MODEXP,50000,0.192703021,259466.612,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,p1024,1024,MODEXP,50000,0.195309293,256004.204,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p1024,1024,MODEXP,50000,0.117389479,425932.548,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p1024,1024,MODEXP,50000,0.119991030,416697.815,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,p1024,1024,EXPONENTIATION,195,0.013709738,14223.467,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,p1024,1024,EXPONENTIATION,195,0.000498112,391477.993,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,p1024,1024,EXPONENTIATION,195,0.007714696,25276.433,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,p1024,1024,EXPONENTIATION,50000,2.202006485,22706.563,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,p1024,1024,EXPONENTIATION,50000,2.213646144,22587.169,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,p1024,1024,EXPONENTIATION,50000,0.539841480,92619.782,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,p1024,1024,EXPONENTIATION,50000,0.540890352,92440.177,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p1024,1024,EXPONENTIATION,50000,0.149283522,334933.148,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p1024,1024,EXPONENTIATION,50000,0.151052116,331011.583,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p1024,1024,EXPONENTIATION,50000,0.140870859,354935.012,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p1024,1024,EXPONENTIATION,50000,0.143957644,347324.383,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,p1024,1024,EXPONENTIATION,50000,0.150984880,331158.987,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,p1024,1024,EXPONENTIATION,50000,0.153596092,325529.116,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p1024,1024,EXPONENTIATION,50000,0.141859841,352460.568,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p1024,1024,EXPONENTIATION,50000,0.144607761,345762.907,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,p1024,1024,DIVIDE,1562,0.000119389,13083278.565,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,p1024,1024,DIVIDE,1562,0.000057879,26987380.390,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,p1024,1024,DIVIDE,1562,0.000087299,17892546.970,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,p1024,1024,DIVIDE,50000,0.000056320,887784090.909,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,p1024,1024,DIVIDE,50000,0.018138060,2756634.410,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,p1024,1024,DIVIDE,50000,0.021069766,2373068.593,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,p1024,1024,DIVIDE,50000,0.010851206,4607782.657,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,p1024,1024,DIVIDE,50000,0.012999985,3846158.313,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p1024,1024,DIVIDE,50000,0.002198088,22747042.152,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p1024,1024,DIVIDE,50000,0.004424675,11300265.241,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p1024,1024,DIVIDE,50000,0.002152818,23225374.783,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p1024,1024,DIVIDE,50000,0.004989286,10021474.108,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,p1024,1024,DIVIDE,50000,0.002171078,23030036.415,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,p1024,1024,DIVIDE,50000,0.004988086,10023885.025,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p1024,1024,DIVIDE,50000,0.002154919,23202729.850,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p1024,1024,DIVIDE,50000,0.005241254,9539701.979,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,p1024,1024,ISQRT,390,0.000184757,2110878.844,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,p1024,1024,ISQRT,390,0.000058549,6661088.360,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,p1024,1024,ISQRT,50000,0.369222595,135419.665,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,p1024,1024,ISQRT,50000,0.372585894,134197.244,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,p1024,1024,ISQRT,50000,0.169112316,295661.494,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,p1024,1024,ISQRT,50000,0.171169667,292107.830,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p1024,1024,ISQRT,50000,0.033554368,1490118.953,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p1024,1024,ISQRT,50000,0.035392382,1412733.407,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p1024,1024,ISQRT,50000,0.031235979,1600718.201,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p1024,1024,ISQRT,50000,0.033654603,1485680.866,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,p1024,1024,ISQRT,50000,0.033958906,1472367.795,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,p1024,1024,ISQRT,50000,0.036166803,1382483.251,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p1024,1024,ISQRT,50000,0.030976138,1614145.722,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p1024,1024,ISQRT,50000,0.033627600,1486873.880,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,p1024,1024,MODMUL_R2,12500,0.009624807,1298727.319,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,p1024,1024,MODMUL_R2,12500,0.000171527,72874711.146,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,p1024,1024,MODMUL_R2,12500,0.008958157,1395376.280,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,p1024,1024,MODMUL_R2,50000,0.000099328,503382731.959,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,p1024,1024,MODMUL_R2,50000,0.003247822,15394933.763,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,p1024,1024,MODMUL_R2,50000,0.005501018,9089226.913,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,p1024,1024,MODMUL_R2,50000,0.000474173,105446721.773,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,p1024,1024,MODMUL_R2,50000,0.002269458,22031691.899,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p1024,1024,MODMUL_R2,50000,0.000269816,185311294.971,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p1024,1024,MODMUL_R2,50000,0.002091959,23901040.680,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p1024,1024,MODMUL_R2,50000,0.000202127,247369471.783,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p1024,1024,MODMUL_R2,50000,0.002446853,20434410.306,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,p1024,1024,MODMUL_R2,50000,0.000223327,223886983.701,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,p1024,1024,MODMUL_R2,50000,0.002462343,20305864.416,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p1024,1024,MODMUL_R2,50000,0.000152108,328714254.138,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p1024,1024,MODMUL_R2,50000,0.002556902,19554915.337,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,p2048,2048,ADD,6250,0.000265346,23554160.042,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,p2048,2048,ADD,6250,0.000064979,96185155.617,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,p2048,2048,ADD,6250,0.006205029,1007247.532,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,p2048,2048,ADD,50000,0.000043008,1162574404.762,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,p2048,2048,ADD,50000,0.000655140,76319534.751,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,p2048,2048,ADD,50000,0.004569662,10941728.765,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,p2048,2048,ADD,50000,0.000258736,193247285.339,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,p2048,2048,ADD,50000,0.004318648,11577697.416,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32,p2048,2048,ADD,50000,0.000149428,334609504.696,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32,p2048,2048,ADD,50000,0.003738835,13373149.699,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p2048,2048,ADD,50000,0.000118298,422661448.652,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p2048,2048,ADD,50000,0.003745805,13348265.570,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p2048,2048,ADD,50000,0.000120418,415219811.598,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p2048,2048,ADD,50000,0.004760050,10504090.870,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,p2048,2048,ADD,50000,0.000097658,511990722.824,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,p2048,2048,ADD,50000,0.004228667,11824057.391,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p2048,2048,ADD,50000,0.000097339,513668505.602,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p2048,2048,ADD,50000,0.005487390,9111799.418,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,p2048,2048,SUBTRACT,6250,0.000217107,28787631.962,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,p2048,2048,SUBTRACT,6250,0.008995437,694796.701,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,p2048,2048,SUBTRACT,6250,0.000084549,73921466.338,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,p2048,2048,SUBTRACT,50000,0.000043008,1162574404.762,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,p2048,2048,SUBTRACT,50000,0.000653281,76536757.454,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,p2048,2048,SUBTRACT,50000,0.004229897,11820619.633,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,p2048,2048,SUBTRACT,50000,0.000330015,151508439.243,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,p2048,2048,SUBTRACT,50000,0.004378338,11419858.719,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32,p2048,2048,SUBTRACT,50000,0.000149007,335554806.088,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32,p2048,2048,SUBTRACT,50000,0.003734495,13388691.018,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p2048,2048,SUBTRACT,50000,0.000118538,421806356.115,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p2048,2048,SUBTRACT,50000,0.003735005,13386861.538,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p2048,2048,SUBTRACT,50000,0.000120999,413225559.951,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p2048,2048,SUBTRACT,50000,0.005354470,9337992.181,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,p2048,2048,SUBTRACT,50000,0.000097849,510989303.764,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,p2048,2048,SUBTRACT,50000,0.004221617,11843802.294,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p2048,2048,SUBTRACT,50000,0.000097298,513884842.974,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p2048,2048,SUBTRACT,50000,0.005586119,8950758.276,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,p2048,2048,ADDMOD,6250,0.000741279,8431376.156,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,p2048,2048,ADDMOD,6250,0.000069709,89658399.856,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,p2048,2048,ADDMOD,6250,0.008993537,694943.478,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,p2048,2048,ADDMOD,50000,0.000043008,1162574404.762,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,p2048,2048,ADDMOD,50000,0.000829038,60310851.037,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,p2048,2048,ADDMOD,50000,0.004998866,10002268.508,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,p2048,2048,ADDMOD,50000,0.000430434,116161824.417,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,p2048,2048,ADDMOD,50000,0.004584415,10906517.246,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32,p2048,2048,ADDMOD,50000,0.000183098,273077778.230,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32,p2048,2048,ADDMOD,50000,0.003825264,13070993.094,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p2048,2048,ADDMOD,50000,0.000207327,241165108.449,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p2048,2048,ADDMOD,50000,0.003762685,13288382.357,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p2048,2048,ADDMOD,50000,0.000205307,243537961.365,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p2048,2048,ADDMOD,50000,0.005941262,8415720.427,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,p2048,2048,ADDMOD,50000,0.000068089,734333076.187,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,p2048,2048,ADDMOD,50000,0.004247507,11771612.038,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p2048,2048,ADDMOD,50000,0.000068309,731965277.143,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p2048,2048,ADDMOD,50000,0.005557099,8997499.736,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,p2048,2048,SUBTRACTMOD,6250,0.000620950,10065222.426,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,p2048,2048,SUBTRACTMOD,6250,0.009380902,666247.254,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,p2048,2048,SUBTRACTMOD,6250,0.008994408,694876.198,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,p2048,2048,SUBTRACTMOD,50000,0.000043008,1162574404.762,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,p2048,2048,SUBTRACTMOD,50000,0.000882977,56626606.665,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,p2048,2048,SUBTRACTMOD,50000,0.005037975,9924622.429,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,p2048,2048,SUBTRACTMOD,50000,0.000459114,108905392.215,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,p2048,2048,SUBTRACTMOD,50000,0.004559025,10967256.848,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32,p2048,2048,SUBTRACTMOD,50000,0.000198508,251878830.757,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32,p2048,2048,SUBTRACTMOD,50000,0.003821574,13083613.720,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p2048,2048,SUBTRACTMOD,50000,0.000202897,246430449.052,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p2048,2048,SUBTRACTMOD,50000,0.003833954,13041367.521,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p2048,2048,SUBTRACTMOD,50000,0.000204517,244478405.818,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p2048,2048,SUBTRACTMOD,50000,0.005764405,8673922.448,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,p2048,2048,SUBTRACTMOD,50000,0.000070999,704236155.546,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,p2048,2048,SUBTRACTMOD,50000,0.004222647,11840914.506,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p2048,2048,SUBTRACTMOD,50000,0.000067669,738886053.441,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p2048,2048,SUBTRACTMOD,50000,0.005565999,8983113.242,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.005386540,1160299.619,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.009008447,693793.268,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.008965898,697085.792,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.130201621,384019.797,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.137006589,364945.951,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.032200592,1552766.478,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.037309619,1340136.972,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.009589299,5214145.678,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.014035045,3562510.800,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.000769359,64989161.269,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.005331041,9379030.712,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.000757939,65968353.626,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.007398341,6758272.024,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.000549142,91051089.399,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.005678236,8805551.554,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.000546373,91512581.233,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.007180806,6963006.869,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.005387100,1160178.962,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.008995168,694817.491,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.000192778,32420680.793,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.000113664,439893018.018,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.047016330,1063460.290,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.051658081,967902.778,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.012026089,4157627.553,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.017191846,2908355.521,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.003037986,16458272.677,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.007573119,6602299.616,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.003923683,12743129.475,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.008507035,5877488.312,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.003921273,12750960.696,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.010372027,4820658.571,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.001743144,28683805.327,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.007067035,7075102.973,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.001604776,31156989.302,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.008522166,5867052.769,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.035107152,178026.403,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.000481623,12976961.491,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.000268686,23261385.308,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000208896,239353553.922,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.068760716,727159.387,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.072503680,689620.169,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.001890133,26453169.563,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.006031294,8290095.116,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000432974,115480323.724,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.004641332,10772770.163,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000676360,73925130.674,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.004372436,11435272.923,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000525992,95058498.405,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.005686466,8792807.163,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000572941,87269051.131,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.004618091,10826984.469,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000476083,105023648.131,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.005902163,8471470.898,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,p2048,2048,COMPARE,6250,0.000062179,100516541.849,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,p2048,2048,COMPARE,6250,0.000061679,101330812.729,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,p2048,2048,COMPARE,6250,0.000093489,66852817.708,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,p2048,2048,COMPARE,50000,0.000043008,1162574404.762,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,p2048,2048,COMPARE,50000,0.000384055,130189659.922,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,p2048,2048,COMPARE,50000,0.003930171,12722092.688,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,p2048,2048,COMPARE,50000,0.000245917,203320531.033,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,p2048,2048,COMPARE,50000,0.004377838,11421162.099,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p2048,2048,COMPARE,50000,0.000132328,377849268.753,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p2048,2048,COMPARE,50000,0.004225048,11834184.417,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p2048,2048,COMPARE,50000,0.000137108,364676254.683,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p2048,2048,COMPARE,50000,0.005355921,9335463.189,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,p2048,2048,COMPARE,50000,0.000035759,1398246984.061,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,p2048,2048,COMPARE,50000,0.004107039,12174221.022,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p2048,2048,COMPARE,50000,0.000034629,1443880621.260,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p2048,2048,COMPARE,50000,0.005617158,8901298.451,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,p2048,2048,REDUCE,781,0.000042879,18214034.546,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,p2048,2048,REDUCE,781,0.000059249,13181683.872,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,p2048,2048,REDUCE,781,0.000099388,7858095.379,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,p2048,2048,REDUCE,50000,0.000049152,1017252604.167,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,p2048,2048,REDUCE,50000,0.354894312,140887.014,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,p2048,2048,REDUCE,50000,0.358742165,139375.866,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,p2048,2048,REDUCE,50000,0.009596984,5209970.185,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,p2048,2048,REDUCE,50000,0.013855613,3608645.746,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p2048,2048,REDUCE,50000,0.001917502,26075593.302,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p2048,2048,REDUCE,50000,0.006027481,8295339.865,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p2048,2048,REDUCE,50000,0.002026150,24677343.693,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p2048,2048,REDUCE,50000,0.007277892,6870121.211,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,p2048,2048,REDUCE,50000,0.001590197,31442636.452,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,p2048,2048,REDUCE,50000,0.005255572,9513713.089,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p2048,2048,REDUCE,50000,0.001917542,26075048.726,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p2048,2048,REDUCE,50000,0.006971629,7171925.088,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,p2048,2048,MODMUL,390,0.000925566,421363.761,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,p2048,2048,MODMUL,390,0.000079589,4900176.828,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,p2048,2048,MODMUL,390,0.000143178,2723885.826,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,p2048,2048,MODMUL,50000,0.001048576,47683715.820,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,p2048,2048,MODMUL,50000,0.552690513,90466.543,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,p2048,2048,MODMUL,50000,0.554997438,90090.506,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,p2048,2048,MODMUL,50000,0.020348310,2457206.519,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,p2048,2048,MODMUL,50000,0.024619749,2030889.932,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p2048,2048,MODMUL,50000,0.010315738,4846962.889,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p2048,2048,MODMUL,50000,0.014631016,3417397.629,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p2048,2048,MODMUL,50000,0.007641517,6543203.320,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p2048,2048,MODMUL,50000,0.011767986,4248815.305,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,p2048,2048,MODMUL,50000,0.009360041,5341856.941,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,p2048,2048,MODMUL,50000,0.013123816,3809867.396,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p2048,2048,MODMUL,50000,0.007669569,6519271.167,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p2048,2048,MODMUL,50000,0.012931121,3866640.782,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,p2048,2048,MODEXP,97,0.341614634,283.946,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,p2048,2048,MODEXP,97,0.007913153,12258.073,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,p2048,2048,MODEXP,97,0.006675472,14530.808,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,p2048,2048,MODEXP,50000,1.078816772,46347.073,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,p2048,2048,MODEXP,50000,57.369559735,871.542,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,p2048,2048,MODEXP,50000,57.380571011,871.375,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,p2048,2048,MODEXP,50000,7.640384655,6544.173,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,p2048,2048,MODEXP,50000,7.648386335,6537.327,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p2048,2048,MODEXP,50000,1.809708254,27628.763,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p2048,2048,MODEXP,50000,1.916543508,26088.633,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p2048,2048,MODEXP,50000,2.242781403,22293.746,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p2048,2048,MODEXP,50000,2.248492301,22237.123,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,p2048,2048,MODEXP,50000,1.572670872,31793.048,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,p2048,2048,MODEXP,50000,1.576615076,31713.511,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p2048,2048,MODEXP,50000,2.206071228,22664.726,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p2048,2048,MODEXP,50000,2.203835516,22687.719,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,p2048,2048,EXPONENTIATION,97,0.040766668,2379.395,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,p2048,2048,EXPONENTIATION,97,0.001425769,68033.453,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,p2048,2048,EXPONENTIATION,97,0.018467157,5252.568,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,p2048,2048,EXPONENTIATION,50000,19.449075830,2570.816,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,p2048,2048,EXPONENTIATION,50000,19.453724421,2570.202,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,p2048,2048,EXPONENTIATION,50000,4.279596885,11683.343,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,p2048,2048,EXPONENTIATION,50000,4.278319772,11686.831,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p2048,2048,EXPONENTIATION,50000,1.276311816,39175.380,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p2048,2048,EXPONENTIATION,50000,1.282367009,38990.398,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p2048,2048,EXPONENTIATION,50000,1.161505458,43047.581,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p2048,2048,EXPONENTIATION,50000,1.166492576,42863.539,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,p2048,2048,EXPONENTIATION,50000,1.238576303,40368.930,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,p2048,2048,EXPONENTIATION,50000,1.243184456,40219.293,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p2048,2048,EXPONENTIATION,50000,1.156261839,43242.800,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p2048,2048,EXPONENTIATION,50000,1.160681463,43078.141,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,p2048,2048,DIVIDE,781,0.000076069,10266990.267,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,p2048,2048,DIVIDE,781,0.000316375,2468589.996,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,p2048,2048,DIVIDE,781,0.000088249,8849925.225,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,p2048,2048,DIVIDE,50000,0.000076800,651041666.667,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,p2048,2048,DIVIDE,50000,1.148441415,43537.267,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,p2048,2048,DIVIDE,50000,1.152903026,43368.782,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,p2048,2048,DIVIDE,50000,0.214525791,233072.209,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,p2048,2048,DIVIDE,50000,0.219166424,228137.135,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p2048,2048,DIVIDE,50000,0.020790303,2404967.340,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p2048,2048,DIVIDE,50000,0.027385046,1825814.001,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p2048,2048,DIVIDE,50000,0.019280205,2593333.419,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p2048,2048,DIVIDE,50000,0.026516626,1885609.420,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,p2048,2048,DIVIDE,50000,0.020869946,2395789.661,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,p2048,2048,DIVIDE,50000,0.027363542,1827248.826,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p2048,2048,DIVIDE,50000,0.019626631,2547558.986,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p2048,2048,DIVIDE,50000,0.025296198,1976581.640,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,p2048,2048,ISQRT,195,0.000144858,1346146.687,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,p2048,2048,ISQRT,195,0.000057299,3403191.504,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,p2048,2048,ISQRT,50000,14.761405093,3387.211,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,p2048,2048,ISQRT,50000,14.772018182,3384.778,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,p2048,2048,ISQRT,50000,5.087027499,9828.923,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,p2048,2048,ISQRT,50000,5.069623789,9862.665,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p2048,2048,ISQRT,50000,0.211038996,236923.038,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p2048,2048,ISQRT,50000,0.217292493,230104.590,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p2048,2048,ISQRT,50000,0.071187334,702372.137,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p2048,2048,ISQRT,50000,0.076897059,650219.929,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,p2048,2048,ISQRT,50000,0.211384083,236536.258,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,p2048,2048,ISQRT,50000,0.216915363,230504.651,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p2048,2048,ISQRT,50000,0.071364751,700625.999,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p2048,2048,ISQRT,50000,0.076078552,657215.450,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-1t,p2048,2048,MODMUL_R2,6250,0.015078367,414501.127,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,gmp-nt,p2048,2048,MODMUL_R2,6250,0.000236297,26449761.550,0
library,AMD EPYC 7532 32-Core Processor,host-cpu,openssl-nt,p2048,2048,MODMUL_R2,6250,0.008991148,695128.153,0
library,NVIDIA A100-SXM4-40GB,gpu,cgbn,p2048,2048,MODMUL_R2,50000,0.000318464,157003617.363,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w8,p2048,2048,MODMUL_R2,50000,0.031992228,1562879.596,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w8,p2048,2048,MODMUL_R2,50000,0.035720976,1399737.787,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w16,p2048,2048,MODMUL_R2,50000,0.002189538,22835867.909,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w16,p2048,2048,MODMUL_R2,50000,0.006289917,7949230.435,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p2048,2048,MODMUL_R2,50000,0.000941366,53114306.519,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-opt,p2048,2048,MODMUL_R2,50000,0.005961782,8386753.741,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p2048,2048,MODMUL_R2,50000,0.000704429,70979462.833,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-o64,p2048,2048,MODMUL_R2,50000,0.005960591,8388429.746,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il,p2048,2048,MODMUL_R2,50000,0.000840378,59497036.291,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il,p2048,2048,MODMUL_R2,50000,0.006396097,7817267.721,0
opencl-kernel,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p2048,2048,MODMUL_R2,50000,0.000670150,74610170.547,0
opencl-e2e,NVIDIA A100-SXM4-40GB,GPU,w32-il64,p2048,2048,MODMUL_R2,50000,0.004988696,10022659.307,0
```
