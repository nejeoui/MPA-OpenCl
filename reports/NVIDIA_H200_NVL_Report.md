# MPA-OpenCL benchmark report - NVIDIA H200 NVL

> **Note.** The multi-threaded GMP and OpenSSL baseline columns have been
> removed from this report: they predate the 2026-09-12 timing fix and were
> understated (see `reports/README.md`). The single-threaded GMP column, the
> OpenCL-on-CPU rows and all MPA measurements are unaffected and were verified
> against GMP before timing.


> **Partial report.** The run was interrupted or hit its time budget.
> Rows that never ran are marked `n/a`.

## 1. System under test

2 OpenCL device(s) exercised with the identical kernels and operands.

### Device 0 - NVIDIA H200 NVL (GPU)

| Property | Value |
|---|---|
| Model | NVIDIA H200 NVL |
| Type | GPU |
| Vendor | NVIDIA Corporation |
| Device memory | 139.80 GiB |
| Max single allocation | 34.95 GiB |
| Local memory | 48 KiB |
| Global cache | 4224 KiB |
| Compute units | 132 |
| Max clock | 1785 MHz |
| Max work-group size | 1024 |
| OpenCL version | OpenCL 3.0 CUDA |
| Driver | 580.159.03 |

### Device 1 - cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor (CPU)

| Property | Value |
|---|---|
| Model | cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor |
| Type | CPU |
| Vendor | AuthenticAMD |
| Device memory | 1498.55 GiB |
| Max single allocation | 512.00 GiB |
| Local memory | 1024 KiB |
| Global cache | 32768 KiB |
| Compute units | 384 |
| Max clock | 4510 MHz |
| Max work-group size | 4096 |
| OpenCL version | OpenCL 3.0 PoCL HSTR: cpu-x86_64-pc-linux-gnu-skylake-avx512 |
| Driver | 5.0+debian |

### Host

| Property | Value |
|---|---|
| CPU | AMD EPYC 9655 96-Core Processor |
| Logical cores | 384 |
| OpenMP threads used | 384 |
| RAM | 1500.5 GB |
| OS | Ubuntu 24.04.4 LTS |
| Kernel | 6.17.0-40-generic |
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
- Total wall time 5521.5 s.

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
| [1] CPU | `mpaKernels_8bits.cl` (w8) | 75 | 75 | 0 | 0 |
| [1] CPU | `mpaKernel_16bits.cl` (w16) | 75 | 75 | 0 | 0 |
| [1] CPU | `mpaKernel_32bits.cl` (w32) | 35 | 35 | 0 | 0 |
| [1] CPU | `mpaKernel_32bits_opt.cl` (w32-opt) | 71 | 71 | 0 | 0 |

**All configurations correct** - 741 configurations, 0 problems.

## 4. Throughput per device

Operations per second, higher is better. Kernel-only timings.

### Device 0 - NVIDIA H200 NVL (GPU)

#### secp256k1 (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | 1.85 G | 2.84 G | 3.83 G | 3.99 G | 3.87 G | 5.28 G | 5.27 G | 95.19 M | 6.18 G |
| SUBTRACT | 50000 / 50000 | 1.86 G | 2.81 G | 3.95 G | 4.00 G | 3.90 G | 5.36 G | 5.39 G | 151.93 M | 6.30 G |
| ADDMOD | 50000 / 50000 | 1.34 G | 2.17 G | 3.56 G | 4.65 G | 4.73 G | 5.54 G | 5.85 G | 38.67 M | 5.33 G |
| SUBTRACTMOD | 50000 / 50000 | 1.28 G | 2.09 G | 3.40 G | 4.56 G | 4.73 G | 5.71 G | 5.79 G | 45.51 M | 4.98 G |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | 58.75 M | 205.34 M | 735.61 M | 2.51 G | 2.67 G | 4.60 G | 5.22 G | 93.92 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | 268.89 M | 819.93 M | 2.22 G | 2.13 G | 2.33 G | 2.74 G | 2.90 G | 94.87 M | 6.48 G |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | 605.84 M | 1.87 G | 4.02 G | 3.50 G | 4.19 G | 3.71 G | 4.66 G | 9.50 M | 4.57 G |
| COMPARE | 50000 / 50000 | 1.82 G | 2.79 G | - | 4.92 G | 4.84 G | 5.30 G | 5.54 G | 287.54 M | 6.25 G |
| REDUCE | 50000 / 6250 | 331.36 M | 493.68 M | - | 1.36 G | 1.36 G | 1.39 G | 1.45 G | 100.82 M | 3.65 G |
| MODMUL | 50000 / 3125 | 122.25 M | 208.38 M | - | 444.48 M | 568.24 M | 450.17 M | 551.63 M | 17.66 M | 1.40 G |
| MODEXP | 50000 / 781 | 2.68 M | 14.18 M | - | 17.73 M | 37.32 M | 18.02 M | 37.31 M | 156.78 k | 4.86 M |
| EXPONENTIATION | 50000 / 781 | 1.83 M | 7.08 M | - | 97.92 M | 142.48 M | 98.55 M | 136.53 M | 478.24 k | n/a |
| DIVIDE | 50000 / 6250 | 140.26 M | 166.77 M | - | 460.06 M | 471.97 M | 462.36 M | 484.77 M | 49.65 M | 2.97 G |
| ISQRT | 50000 / 1562 | 12.31 M | 16.60 M | - | 68.35 M | 74.92 M | 67.27 M | 75.81 M | 21.14 M | n/a |
| MODMUL_R2 | 50000 / 50000 | 506.53 M | 1.73 G | - | 2.65 G | 3.61 G | 2.68 G | 3.77 G | 17.37 M | 3.16 G |

#### rsa256(composite) (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | 1.82 G | 2.83 G | 3.78 G | 3.96 G | 3.94 G | 5.17 G | 5.20 G | 92.83 M | 6.30 G |
| SUBTRACT | 50000 / 50000 | 1.85 G | 2.83 G | 3.86 G | 4.00 G | 3.91 G | 5.37 G | 5.39 G | 150.84 M | 6.33 G |
| ADDMOD | 50000 / 50000 | 1.46 G | 2.30 G | 3.73 G | 4.62 G | 4.87 G | 5.76 G | 5.82 G | 44.69 M | 4.94 G |
| SUBTRACTMOD | 50000 / 50000 | 1.30 G | 2.15 G | 3.66 G | 4.82 G | 4.67 G | 5.66 G | 5.81 G | 44.83 M | 4.96 G |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | 59.47 M | 209.64 M | 730.03 M | 2.59 G | 2.74 G | 4.60 G | 5.26 G | 90.73 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | 269.77 M | 827.81 M | 2.23 G | 2.20 G | 2.38 G | 2.75 G | 2.87 G | 92.20 M | 6.35 G |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | 612.36 M | 1.85 G | 4.03 G | 3.52 G | 4.10 G | 3.70 G | 4.64 G | 9.48 M | 4.62 G |
| COMPARE | 50000 / 50000 | 1.91 G | 2.76 G | - | 4.90 G | 4.81 G | 5.41 G | 5.51 G | 247.87 M | 6.30 G |
| REDUCE | 50000 / 6250 | 329.90 M | 491.84 M | - | 1.39 G | 1.38 G | 1.38 G | 1.45 G | 62.67 M | 3.69 G |
| MODMUL | 50000 / 3125 | 121.89 M | 208.62 M | - | 443.93 M | 568.50 M | 451.51 M | 550.65 M | 17.56 M | 1.40 G |
| MODEXP | 50000 / 781 | 2.68 M | 14.17 M | - | 17.87 M | 37.40 M | 18.02 M | 37.57 M | 171.61 k | 4.95 M |
| EXPONENTIATION | 50000 / 781 | 1.83 M | 7.10 M | - | 97.68 M | 143.25 M | 98.77 M | 136.75 M | 479.62 k | n/a |
| DIVIDE | 50000 / 6250 | 139.41 M | 163.53 M | - | 439.90 M | 455.37 M | 445.71 M | 467.77 M | 48.52 M | 2.96 G |
| ISQRT | 50000 / 1562 | 12.32 M | 16.18 M | - | 65.85 M | 69.71 M | 64.62 M | 70.38 M | 21.19 M | n/a |
| MODMUL_R2 | 50000 / 50000 | 504.59 M | 1.71 G | - | 2.55 G | 3.50 G | 2.66 G | 3.74 G | 17.48 M | 3.12 G |

#### brainpoolP512r1 (512-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 25000 | 843.30 M | 1.52 G | 2.54 G | 2.62 G | 2.51 G | 4.45 G | 3.20 G | 84.09 M | 5.92 G |
| SUBTRACT | 50000 / 25000 | 835.70 M | 1.52 G | 2.53 G | 2.63 G | 2.51 G | 4.37 G | 3.22 G | 131.33 M | 5.87 G |
| ADDMOD | 50000 / 25000 | 644.00 M | 1.24 G | 2.27 G | 2.09 G | 2.06 G | 4.22 G | 4.17 G | 38.32 M | 4.75 G |
| SUBTRACTMOD | 50000 / 25000 | 578.23 M | 1.13 G | 2.08 G | 2.12 G | 2.09 G | 4.26 G | 4.26 G | 40.93 M | 4.62 G |
| MULTIPLYOPERANDSCANNING | 50000 / 25000 | 9.46 M | 35.57 M | 141.19 M | 834.72 M | 856.46 M | 1.63 G | 1.67 G | 36.05 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 25000 | 39.65 M | 144.03 M | 510.36 M | 507.51 M | 510.62 M | 1.08 G | 830.98 M | 36.28 M | 4.73 G |
| MONTGOMERYMULTIPLICATION | 50000 / 25000 | 173.09 M | 583.23 M | 1.76 G | 1.28 G | 1.83 G | 1.54 G | 2.43 G | 4.15 M | 3.43 G |
| COMPARE | 50000 / 25000 | 888.26 M | 1.66 G | - | 2.64 G | 2.65 G | 4.52 G | 4.59 G | 232.13 M | 5.99 G |
| REDUCE | 50000 / 3125 | 110.43 M | 140.17 M | - | 494.27 M | 495.34 M | 485.44 M | 413.53 M | 60.48 M | 2.58 G |
| MODMUL | 50000 / 1562 | 39.29 M | 54.32 M | - | 137.10 M | 169.73 M | 139.41 M | 137.84 M | 8.62 M | 525.74 M |
| MODEXP | 50000 / 390 | 200.74 k | 2.02 M | - | 2.32 M | 5.19 M | 2.27 M | 5.24 M | 33.69 k | 2.06 M |
| EXPONENTIATION | 50000 / 390 | 240.72 k | 963.66 k | - | 3.05 M | 3.42 M | 3.19 M | 3.57 M | 149.29 k | n/a |
| DIVIDE | 50000 / 3125 | 39.61 M | 40.23 M | - | 120.83 M | 129.28 M | 126.33 M | 141.35 M | 44.68 M | 1.81 G |
| ISQRT | 50000 / 781 | 2.62 M | 2.64 M | - | 13.46 M | 14.83 M | 14.17 M | 11.73 M | 12.81 M | n/a |
| MODMUL_R2 | 50000 / 25000 | 114.06 M | 587.47 M | - | 833.47 M | 1.42 G | 912.57 M | 1.57 G | 8.47 M | 2.17 G |

#### p1024 (1024-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 12500 | 239.18 M | 466.11 M | 1.14 G | 1.14 G | 1.14 G | 2.18 G | 2.20 G | 66.62 M | 4.54 G |
| SUBTRACT | 50000 / 12500 | 239.43 M | 467.94 M | 1.14 G | 1.14 G | 1.13 G | 2.20 G | 2.23 G | 95.89 M | 4.65 G |
| ADDMOD | 50000 / 12500 | 176.12 M | 347.44 M | 914.75 M | 895.08 M | 917.60 M | 3.10 G | 3.16 G | 26.21 M | 4.08 G |
| SUBTRACTMOD | 50000 / 12500 | 175.13 M | 349.09 M | 907.77 M | 890.47 M | 896.70 M | 3.05 G | 3.15 G | 32.20 M | 3.92 G |
| MULTIPLYOPERANDSCANNING | 50000 / 12500 | 1.98 M | 7.71 M | 35.55 M | 351.37 M | 348.58 M | 549.81 M | 607.68 M | 10.20 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 12500 | 5.33 M | 20.74 M | 80.91 M | 80.73 M | 80.87 M | 231.45 M | 248.31 M | 10.18 M | 2.02 G |
| MONTGOMERYMULTIPLICATION | 50000 / 12500 | 24.63 M | 175.92 M | 612.29 M | 364.70 M | 505.51 M | 497.86 M | 885.74 M | 1.33 M | 1.34 G |
| COMPARE | 50000 / 12500 | 388.08 M | 771.23 M | - | 1.37 G | 1.38 G | 3.48 G | 3.49 G | 204.23 M | 4.54 G |
| REDUCE | 50000 / 1562 | 20.62 M | 38.38 M | - | 147.04 M | 140.53 M | 145.15 M | 148.22 M | 80.94 M | 1.84 G |
| MODMUL | 50000 / 781 | 5.38 M | 14.83 M | - | 30.99 M | 41.13 M | 33.19 M | 42.37 M | 3.14 M | 213.54 M |
| MODEXP | 50000 / 195 | 17.11 k | 288.90 k | - | 289.01 k | 552.78 k | 288.33 k | 550.36 k | 5.35 k | 398.62 k |
| EXPONENTIATION | 50000 / 195 | 28.36 k | 123.08 k | - | 400.43 k | 462.25 k | 404.83 k | 460.80 k | 36.33 k | n/a |
| DIVIDE | 50000 / 1562 | 4.09 M | 6.49 M | - | 30.01 M | 29.64 M | 30.40 M | 29.46 M | 39.16 M | 1.35 G |
| ISQRT | 50000 / 390 | 235.12 k | 368.43 k | - | 1.98 M | 2.01 M | 1.97 M | 2.07 M | 6.38 M | n/a |
| MODMUL_R2 | 50000 / 12500 | 16.71 M | 187.40 M | - | 232.49 M | 342.49 M | 278.68 M | 464.86 M | 3.09 M | 746.18 M |

#### p2048 (2048-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 6250 | 120.98 M | 244.38 M | 623.91 M | 607.75 M | 587.68 M | 888.40 M | 742.72 M | 46.86 M | 2.19 G |
| SUBTRACT | 50000 / 6250 | 120.94 M | 245.29 M | 625.39 M | 611.85 M | 584.32 M | 895.08 M | 740.84 M | 57.37 M | 2.18 G |
| ADDMOD | 50000 / 6250 | 96.76 M | 189.44 M | 492.61 M | 453.68 M | 448.19 M | 1.45 G | 1.38 G | 20.29 M | 2.16 G |
| SUBTRACTMOD | 50000 / 6250 | 92.71 M | 181.94 M | 446.95 M | 455.95 M | 447.58 M | 1.31 G | 1.46 G | 22.69 M | 2.16 G |
| MULTIPLYOPERANDSCANNING | 50000 / 6250 | 497.63 k | 1.97 M | 8.74 M | 120.96 M | 125.87 M | 137.69 M | 176.17 M | 3.09 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 6250 | 1.35 M | 5.29 M | 20.93 M | 20.63 M | 20.68 M | 47.85 M | 50.97 M | 3.08 M | 702.25 M |
| MONTGOMERYMULTIPLICATION | 50000 / 6250 | 1.10 M | 31.98 M | 205.97 M | 115.90 M | 160.75 M | 130.14 M | 189.45 M | 425.91 k | 456.87 M |
| COMPARE | 50000 / 6250 | 195.25 M | 397.64 M | - | 744.49 M | 710.53 M | 2.24 G | 2.04 G | 188.28 M | 2.16 G |
| REDUCE | 50000 / 781 | 352.92 k | 9.40 M | - | 41.47 M | 46.42 M | 42.94 M | 44.78 M | 52.79 M | 1.61 G |
| MODMUL | 50000 / 390 | 186.40 k | 3.01 M | - | 8.47 M | 12.00 M | 8.36 M | 12.65 M | 997.46 k | 104.41 M |
| MODEXP | 50000 / 97 | 1.15 k | 7.93 k | - | 37.19 k | 31.78 k | 36.28 k | 31.91 k | 727.2 | 71.30 k |
| EXPONENTIATION | 50000 / 97 | 3.31 k | 14.36 k | - | 48.77 k | 55.89 k | 47.63 k | 55.83 k | 6.34 k | n/a |
| DIVIDE | 50000 / 781 | 81.36 k | 318.79 k | - | 3.87 M | 4.37 M | 3.89 M | 4.22 M | 32.02 M | 1.36 G |
| ISQRT | 50000 / 195 | 5.44 k | 11.48 k | - | 296.90 k | 864.38 k | 295.65 k | 859.95 k | 3.87 M | n/a |
| MODMUL_R2 | 50000 / 6250 | 2.10 M | 37.55 M | - | 64.66 M | 94.69 M | 72.04 M | 106.56 M | 965.60 k | 243.46 M |

### Device 1 - cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor (CPU)

#### secp256k1 (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | 44.59 M | 26.38 M | 25.64 M | 22.51 M | - | - | - | 95.19 M | 6.18 G |
| SUBTRACT | 50000 / 50000 | 46.90 M | 24.87 M | 26.60 M | 20.14 M | - | - | - | 151.93 M | 6.30 G |
| ADDMOD | 50000 / 50000 | 40.97 M | 26.54 M | 21.25 M | 24.90 M | - | - | - | 38.67 M | 5.33 G |
| SUBTRACTMOD | 50000 / 50000 | 51.82 M | 26.94 M | 38.26 M | 20.36 M | - | - | - | 45.51 M | 4.98 G |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | 39.85 M | 24.45 M | 39.98 M | 21.31 M | - | - | - | 93.92 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | 37.07 M | 19.53 M | 20.30 M | 24.65 M | - | - | - | 94.87 M | 6.48 G |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | 24.89 M | 29.11 M | 23.23 M | 24.77 M | - | - | - | 9.50 M | 4.57 G |
| COMPARE | 50000 / 50000 | 52.71 M | 26.46 M | - | 20.34 M | - | - | - | 287.54 M | 6.25 G |
| REDUCE | 50000 / 6250 | 31.05 M | 26.24 M | - | 28.13 M | - | - | - | 100.82 M | 3.65 G |
| MODMUL | 50000 / 3125 | 9.96 M | 13.10 M | - | 23.07 M | - | - | - | 17.66 M | 1.40 G |
| MODEXP | 50000 / 781 | 56.34 k | 248.06 k | - | 475.86 k | - | - | - | 156.78 k | 4.86 M |
| EXPONENTIATION | 50000 / 781 | 158.05 k | 516.48 k | - | 3.63 M | - | - | - | 478.24 k | n/a |
| DIVIDE | 50000 / 6250 | 18.42 M | 18.21 M | - | 25.46 M | - | - | - | 49.65 M | 2.97 G |
| ISQRT | 50000 / 1562 | 2.80 M | 2.76 M | - | 5.35 M | - | - | - | 21.14 M | n/a |
| MODMUL_R2 | 50000 / 50000 | 21.11 M | 25.21 M | - | 24.53 M | - | - | - | 17.37 M | 3.16 G |

#### rsa256(composite) (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | 21.91 M | 22.33 M | 24.89 M | 24.94 M | - | - | - | 92.83 M | 6.30 G |
| SUBTRACT | 50000 / 50000 | 26.11 M | 23.52 M | 21.31 M | 18.74 M | - | - | - | 150.84 M | 6.33 G |
| ADDMOD | 50000 / 50000 | 26.51 M | 23.57 M | 26.62 M | 19.35 M | - | - | - | 44.69 M | 4.94 G |
| SUBTRACTMOD | 50000 / 50000 | 23.96 M | 23.16 M | 21.67 M | 21.93 M | - | - | - | 44.83 M | 4.96 G |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | 24.33 M | 25.07 M | 31.51 M | 24.22 M | - | - | - | 90.73 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | 21.55 M | 25.61 M | 22.75 M | 17.42 M | - | - | - | 92.20 M | 6.35 G |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | 21.24 M | 26.81 M | 23.57 M | 24.19 M | - | - | - | 9.48 M | 4.62 G |
| COMPARE | 50000 / 50000 | 23.89 M | 22.24 M | - | 19.06 M | - | - | - | 247.87 M | 6.30 G |
| REDUCE | 50000 / 6250 | 24.55 M | 26.08 M | - | 22.02 M | - | - | - | 62.67 M | 3.69 G |
| MODMUL | 50000 / 3125 | 9.67 M | 12.70 M | - | 23.94 M | - | - | - | 17.56 M | 1.40 G |
| MODEXP | 50000 / 781 | 61.29 k | 249.95 k | - | 480.18 k | - | - | - | 171.61 k | 4.95 M |
| EXPONENTIATION | 50000 / 781 | 155.64 k | 523.16 k | - | 3.66 M | - | - | - | 479.62 k | n/a |
| DIVIDE | 50000 / 6250 | 17.34 M | 17.88 M | - | 23.32 M | - | - | - | 48.52 M | 2.96 G |
| ISQRT | 50000 / 1562 | 2.72 M | 3.52 M | - | 5.22 M | - | - | - | 21.19 M | n/a |
| MODMUL_R2 | 50000 / 50000 | 20.37 M | 28.80 M | - | 23.03 M | - | - | - | 17.48 M | 3.12 G |

#### brainpoolP512r1 (512-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 25000 | 26.32 M | 21.83 M | 22.23 M | 25.68 M | - | - | - | 84.09 M | 5.92 G |
| SUBTRACT | 50000 / 25000 | 24.86 M | 37.43 M | 23.93 M | 24.35 M | - | - | - | 131.33 M | 5.87 G |
| ADDMOD | 50000 / 25000 | 24.73 M | 21.92 M | 24.63 M | 22.20 M | - | - | - | 38.32 M | 4.75 G |
| SUBTRACTMOD | 50000 / 25000 | 24.10 M | 21.33 M | 24.09 M | 23.89 M | - | - | - | 40.93 M | 4.62 G |
| MULTIPLYOPERANDSCANNING | 50000 / 25000 | 16.60 M | 21.99 M | 22.99 M | 22.77 M | - | - | - | 36.05 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 25000 | 30.72 M | 24.06 M | 21.97 M | 22.55 M | - | - | - | 36.28 M | 4.73 G |
| MONTGOMERYMULTIPLICATION | 50000 / 25000 | 7.89 M | 19.15 M | 24.39 M | 24.69 M | - | - | - | 4.15 M | 3.43 G |
| COMPARE | 50000 / 25000 | 23.61 M | 24.06 M | - | 24.79 M | - | - | - | 232.13 M | 5.99 G |
| REDUCE | 50000 / 3125 | 12.39 M | 12.34 M | - | 13.53 M | - | - | - | 60.48 M | 2.58 G |
| MODMUL | 50000 / 1562 | 3.05 M | 5.02 M | - | 7.27 M | - | - | - | 8.62 M | 525.74 M |
| MODEXP | 50000 / 390 | 7.33 k | 36.29 k | - | 66.86 k | - | - | - | 33.69 k | 2.06 M |
| EXPONENTIATION | 50000 / 390 | 15.44 k | 79.85 k | - | 131.75 k | - | - | - | 149.29 k | n/a |
| DIVIDE | 50000 / 3125 | 6.92 M | 8.45 M | - | 10.45 M | - | - | - | 44.68 M | 1.81 G |
| ISQRT | 50000 / 781 | 172.39 k | 401.16 k | - | 472.42 k | - | - | - | 12.81 M | n/a |
| MODMUL_R2 | 50000 / 25000 | 8.03 M | 19.17 M | - | 26.11 M | - | - | - | 8.47 M | 2.17 G |

#### p1024 (1024-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 12500 | 24.48 M | 25.05 M | 23.73 M | 22.34 M | - | - | - | 66.62 M | 4.54 G |
| SUBTRACT | 50000 / 12500 | 24.71 M | 25.29 M | 24.42 M | 24.78 M | - | - | - | 95.89 M | 4.65 G |
| ADDMOD | 50000 / 12500 | 25.26 M | 25.43 M | 22.00 M | 23.24 M | - | - | - | 26.21 M | 4.08 G |
| SUBTRACTMOD | 50000 / 12500 | 25.06 M | 25.11 M | 23.98 M | 24.31 M | - | - | - | 32.20 M | 3.92 G |
| MULTIPLYOPERANDSCANNING | 50000 / 12500 | 6.61 M | 19.14 M | 25.00 M | 29.66 M | - | - | - | 10.20 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 12500 | 8.87 M | 24.33 M | 23.64 M | 24.53 M | - | - | - | 10.18 M | 2.02 G |
| MONTGOMERYMULTIPLICATION | 50000 / 12500 | 2.75 M | 7.97 M | 18.39 M | 21.49 M | - | - | - | 1.33 M | 1.34 G |
| COMPARE | 50000 / 12500 | 21.75 M | 24.38 M | - | 21.79 M | - | - | - | 204.23 M | 4.54 G |
| REDUCE | 50000 / 1562 | 4.48 M | 4.91 M | - | 5.81 M | - | - | - | 80.94 M | 1.84 G |
| MODMUL | 50000 / 781 | 259.56 k | 512.95 k | - | 2.25 M | - | - | - | 3.14 M | 213.54 M |
| MODEXP | 50000 / 195 | 895.8 | 4.60 k | - | 9.18 k | - | - | - | 5.35 k | 398.62 k |
| EXPONENTIATION | 50000 / 195 | 2.21 k | 9.22 k | - | over budget | - | - | - | 36.33 k | n/a |
| DIVIDE | 50000 / 1562 | 562.00 k | 2.37 M | - | 3.17 M | - | - | - | 39.16 M | 1.35 G |
| ISQRT | 50000 / 390 | 42.06 k | 60.31 k | - | 83.93 k | - | - | - | 6.38 M | n/a |
| MODMUL_R2 | 50000 / 12500 | 2.72 M | 9.13 M | - | 13.48 M | - | - | - | 3.09 M | 746.18 M |

#### p2048 (2048-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 6250 | 23.48 M | 24.09 M | 24.20 M | over budget | - | - | - | 46.86 M | 2.19 G |
| SUBTRACT | 50000 / 6250 | 23.10 M | 22.60 M | 24.43 M | 24.00 M | - | - | - | 57.37 M | 2.18 G |
| ADDMOD | 50000 / 6250 | 22.84 M | 22.14 M | 22.22 M | 21.78 M | - | - | - | 20.29 M | 2.16 G |
| SUBTRACTMOD | 50000 / 6250 | 23.02 M | 24.97 M | 24.94 M | 22.08 M | - | - | - | 22.69 M | 2.16 G |
| MULTIPLYOPERANDSCANNING | 50000 / 6250 | 532.00 k | 8.26 M | 20.07 M | 13.03 M | - | - | - | 3.09 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 6250 | 2.86 M | 11.11 M | 23.99 M | 11.42 M | - | - | - | 3.08 M | 702.25 M |
| MONTGOMERYMULTIPLICATION | 50000 / 6250 | 172.30 k | 2.68 M | 7.64 M | 8.19 M | - | - | - | 425.91 k | 456.87 M |
| COMPARE | 50000 / 6250 | 23.94 M | 23.69 M | - | 20.32 M | - | - | - | 188.28 M | 2.16 G |
| REDUCE | 50000 / 781 | 264.76 k | 528.38 k | - | 504.28 k | - | - | - | 52.79 M | 1.61 G |
| MODMUL | 50000 / 390 | 71.34 k | 123.71 k | - | 170.64 k | - | - | - | 997.46 k | 104.41 M |
| MODEXP | 50000 / 97 | over budget | over budget | - | over budget | - | - | - | 727.2 | 71.30 k |
| EXPONENTIATION | 50000 / 97 | over budget | 1.10 k | - | - | - | - | - | 6.34 k | n/a |
| DIVIDE | 50000 / 781 | 227.31 k | 248.12 k | - | - | - | - | - | 32.02 M | 1.36 G |
| ISQRT | 50000 / 195 | 9.68 k | 14.32 k | - | - | - | - | - | 3.87 M | n/a |
| MODMUL_R2 | 50000 / 6250 | 167.59 k | 2.87 M | - | - | - | - | - | 965.60 k | 243.46 M |

## 5. Head to head

Best OpenCL GPU result against best OpenCL CPU result and the CPU libraries.
Ratios above 1.00x mean the GPU is faster than that baseline.

### secp256k1 (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w32-il | 5.28 G | w8 | 44.59 M | 95.19 M | 6.18 G | 118.40x | 0.85x |
| SUBTRACT | w32-il64 | 5.39 G | w8 | 46.90 M | 151.93 M | 6.30 G | 115.01x | 0.86x |
| ADDMOD | w32-il64 | 5.85 G | w8 | 40.97 M | 38.67 M | 5.33 G | 142.74x | 1.10x |
| SUBTRACTMOD | w32-il64 | 5.79 G | w8 | 51.82 M | 45.51 M | 4.98 G | 111.68x | 1.16x |
| MULTIPLYOPERANDSCANNING | w32-il64 | 5.22 G | w32 | 39.98 M | 93.92 M | n/a | 130.56x | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 2.90 G | w8 | 37.07 M | 94.87 M | 6.48 G | 78.23x | 0.45x |
| MONTGOMERYMULTIPLICATION | w32-il64 | 4.66 G | w16 | 29.11 M | 9.50 M | 4.57 G | 160.20x | 1.02x |
| COMPARE | w32-il64 | 5.54 G | w8 | 52.71 M | 287.54 M | 6.25 G | 105.04x | 0.89x |
| REDUCE | w32-il64 | 180.64 M | w8 | 3.88 M | 100.82 M | 3.65 G | 46.55x | 0.05x |
| MODMUL | w32-o64 | 35.51 M | w32-opt | 1.44 M | 17.66 M | 1.40 G | 24.63x | 0.03x |
| MODEXP | w32-o64 | 582.86 k | w32-opt | 7.43 k | 156.78 k | 4.86 M | 78.42x | 0.12x |
| EXPONENTIATION | w32-o64 | 2.23 M | w32-opt | 56.66 k | 478.24 k | n/a | 39.28x | n/a |
| DIVIDE | w32-il64 | 60.60 M | w32-opt | 3.18 M | 49.65 M | 2.97 G | 19.04x | 0.02x |
| ISQRT | w32-il64 | 2.37 M | w32-opt | 167.29 k | 21.14 M | n/a | 14.16x | n/a |
| MODMUL_R2 | w32-il64 | 3.77 G | w16 | 25.21 M | 17.37 M | 3.16 G | 149.47x | 1.19x |

### rsa256(composite) (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w32-il64 | 5.20 G | w32-opt | 24.94 M | 92.83 M | 6.30 G | 208.59x | 0.83x |
| SUBTRACT | w32-il64 | 5.39 G | w8 | 26.11 M | 150.84 M | 6.33 G | 206.37x | 0.85x |
| ADDMOD | w32-il64 | 5.82 G | w32 | 26.62 M | 44.69 M | 4.94 G | 218.69x | 1.18x |
| SUBTRACTMOD | w32-il64 | 5.81 G | w8 | 23.96 M | 44.83 M | 4.96 G | 242.66x | 1.17x |
| MULTIPLYOPERANDSCANNING | w32-il64 | 5.26 G | w32 | 31.51 M | 90.73 M | n/a | 166.86x | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 2.87 G | w16 | 25.61 M | 92.20 M | 6.35 G | 112.22x | 0.45x |
| MONTGOMERYMULTIPLICATION | w32-il64 | 4.64 G | w16 | 26.81 M | 9.48 M | 4.62 G | 173.16x | 1.00x |
| COMPARE | w32-il64 | 5.51 G | w8 | 23.89 M | 247.87 M | 6.30 G | 230.79x | 0.87x |
| REDUCE | w32-il64 | 180.90 M | w16 | 3.26 M | 62.67 M | 3.69 G | 55.49x | 0.05x |
| MODMUL | w32-o64 | 35.53 M | w32-opt | 1.50 M | 17.56 M | 1.40 G | 23.75x | 0.03x |
| MODEXP | w32-il64 | 586.86 k | w32-opt | 7.50 k | 171.61 k | 4.95 M | 78.24x | 0.12x |
| EXPONENTIATION | w32-o64 | 2.24 M | w32-opt | 57.22 k | 479.62 k | n/a | 39.10x | n/a |
| DIVIDE | w32-il64 | 58.47 M | w32-opt | 2.91 M | 48.52 M | 2.96 G | 20.06x | 0.02x |
| ISQRT | w32-il64 | 2.20 M | w32-opt | 163.16 k | 21.19 M | n/a | 13.48x | n/a |
| MODMUL_R2 | w32-il64 | 3.74 G | w16 | 28.80 M | 17.48 M | 3.12 G | 129.76x | 1.20x |

### brainpoolP512r1 (512-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w32-il | 2.23 G | w8 | 13.16 M | 84.09 M | 5.92 G | 169.18x | 0.38x |
| SUBTRACT | w32-il | 2.18 G | w16 | 18.71 M | 131.33 M | 5.87 G | 116.67x | 0.37x |
| ADDMOD | w32-il | 2.11 G | w8 | 12.36 M | 38.32 M | 4.75 G | 170.62x | 0.44x |
| SUBTRACTMOD | w32-il | 2.13 G | w8 | 12.05 M | 40.93 M | 4.62 G | 176.71x | 0.46x |
| MULTIPLYOPERANDSCANNING | w32-il64 | 835.00 M | w32 | 11.49 M | 36.05 M | n/a | 72.65x | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il | 538.56 M | w8 | 15.36 M | 36.28 M | 4.73 G | 35.07x | 0.11x |
| MONTGOMERYMULTIPLICATION | w32-il64 | 1.21 G | w32-opt | 12.34 M | 4.15 M | 3.43 G | 98.31x | 0.35x |
| COMPARE | w32-il64 | 2.29 G | w32-opt | 12.39 M | 232.13 M | 5.99 G | 185.05x | 0.38x |
| REDUCE | w32-o64 | 30.96 M | w32-opt | 845.34 k | 60.48 M | 2.58 G | 36.62x | 0.01x |
| MODMUL | w32-o64 | 5.30 M | w32-opt | 227.13 k | 8.62 M | 525.74 M | 23.34x | 0.01x |
| MODEXP | w32-il64 | 40.88 k | w32-opt | 521.5 | 33.69 k | 2.06 M | 78.39x | 0.02x |
| EXPONENTIATION | w32-il64 | 27.88 k | w32-opt | 1.03 k | 149.29 k | n/a | 27.13x | n/a |
| DIVIDE | w32-il64 | 8.83 M | w32-opt | 652.82 k | 44.68 M | 1.81 G | 13.53x | 0.00x |
| ISQRT | w32-o64 | 231.60 k | w32-opt | 7.38 k | 12.81 M | n/a | 31.38x | n/a |
| MODMUL_R2 | w32-il64 | 786.91 M | w32-opt | 13.06 M | 8.47 M | 2.17 G | 60.27x | 0.36x |

### p1024 (1024-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w32-il64 | 549.69 M | w16 | 6.26 M | 66.62 M | 4.54 G | 87.79x | 0.12x |
| SUBTRACT | w32-il64 | 556.30 M | w16 | 6.32 M | 95.89 M | 4.65 G | 87.98x | 0.12x |
| ADDMOD | w32-il64 | 790.64 M | w16 | 6.36 M | 26.21 M | 4.08 G | 124.38x | 0.19x |
| SUBTRACTMOD | w32-il64 | 786.66 M | w16 | 6.28 M | 32.20 M | 3.92 G | 125.33x | 0.20x |
| MULTIPLYOPERANDSCANNING | w32-il64 | 151.92 M | w32-opt | 7.41 M | 10.20 M | n/a | 20.49x | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 62.08 M | w32-opt | 6.13 M | 10.18 M | 2.02 G | 10.12x | 0.03x |
| MONTGOMERYMULTIPLICATION | w32-il64 | 221.43 M | w32-opt | 5.37 M | 1.33 M | 1.34 G | 41.21x | 0.17x |
| COMPARE | w32-il64 | 872.91 M | w16 | 6.09 M | 204.23 M | 4.54 G | 143.23x | 0.19x |
| REDUCE | w32-il64 | 4.63 M | w32-opt | 181.58 k | 80.94 M | 1.84 G | 25.50x | 0.00x |
| MODMUL | w32-il64 | 661.77 k | w32-opt | 35.19 k | 3.14 M | 213.54 M | 18.81x | 0.00x |
| MODEXP | w32-o64 | 2.16 k | w32-opt | 35.8 | 5.35 k | 398.62 k | 60.24x | 0.01x |
| EXPONENTIATION | w32-o64 | 1.80 k | w16 | 36.0 | 36.33 k | n/a | 50.14x | n/a |
| DIVIDE | w32-il | 949.70 k | w32-opt | 98.91 k | 39.16 M | 1.35 G | 9.60x | 0.00x |
| ISQRT | w32-il64 | 16.11 k | w32-opt | 654.6 | 6.38 M | n/a | 24.61x | n/a |
| MODMUL_R2 | w32-il64 | 116.21 M | w32-opt | 3.37 M | 3.09 M | 746.18 M | 34.48x | 0.16x |

### p2048 (2048-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w32-il | 111.05 M | w32 | 3.03 M | 46.86 M | 2.19 G | 36.70x | 0.05x |
| SUBTRACT | w32-il | 111.88 M | w32 | 3.05 M | 57.37 M | 2.18 G | 36.64x | 0.05x |
| ADDMOD | w32-il | 181.74 M | w8 | 2.86 M | 20.29 M | 2.16 G | 63.65x | 0.08x |
| SUBTRACTMOD | w32-il64 | 182.11 M | w16 | 3.12 M | 22.69 M | 2.16 G | 58.35x | 0.08x |
| MULTIPLYOPERANDSCANNING | w32-il64 | 22.02 M | w32 | 2.51 M | 3.09 M | n/a | 8.78x | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 6.37 M | w32 | 3.00 M | 3.08 M | 702.25 M | 2.12x | 0.01x |
| MONTGOMERYMULTIPLICATION | w32 | 25.75 M | w32-opt | 1.02 M | 425.91 k | 456.87 M | 25.16x | 0.06x |
| COMPARE | w32-il | 280.14 M | w8 | 2.99 M | 188.28 M | 2.16 G | 93.63x | 0.13x |
| REDUCE | w32-o64 | 725.02 k | w16 | 8.25 k | 52.79 M | 1.61 G | 87.85x | 0.00x |
| MODMUL | w32-il64 | 98.63 k | w32-opt | 1.33 k | 997.46 k | 104.41 M | 74.11x | 0.00x |
| MODEXP | w32-opt | 72.1 | none | n/a | 727.2 | 71.30 k | n/a | 0.00x |
| EXPONENTIATION | w32-o64 | 108.4 | w16 | 2.1 | 6.34 k | n/a | 51.03x | n/a |
| DIVIDE | w32-o64 | 68.30 k | w16 | 3.88 k | 32.02 M | 1.36 G | 17.62x | 0.00x |
| ISQRT | w32-o64 | 3.37 k | w16 | 55.8 | 3.87 M | n/a | 60.38x | n/a |
| MODMUL_R2 | w32-il64 | 13.32 M | w16 | 358.18 k | 965.60 k | 243.46 M | 37.19x | 0.05x |

## 6. Raw data

Also written to `NVIDIA_H200_NVL_Report.csv` for analysis.

```csv
kind,device,device_type,kernel,modulus,bits,operation,items,seconds,ops_per_sec,mismatches
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,secp256k1,256,ADD,50000,0.000525289,95185653.518,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,secp256k1,256,ADD,50000,0.002177506,22962046.621,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,secp256k1,256,ADD,50000,0.001494164,33463517.766,0
library,NVIDIA H200 NVL,gpu,cgbn,secp256k1,256,ADD,50000,0.000008096,6175889328.063,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,secp256k1,256,ADD,50000,0.000027000,1851853720.120,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,secp256k1,256,ADD,50000,0.000301491,165842427.874,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,secp256k1,256,ADD,50000,0.000017590,2842523210.917,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,secp256k1,256,ADD,50000,0.000300611,166327965.462,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32,secp256k1,256,ADD,50000,0.000013040,3834347170.418,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32,secp256k1,256,ADD,50000,0.000297191,168241922.994,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,secp256k1,256,ADD,50000,0.000012530,3990418552.103,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,secp256k1,256,ADD,50000,0.000285610,175063842.581,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,secp256k1,256,ADD,50000,0.000012930,3866970951.129,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,secp256k1,256,ADD,50000,0.000280280,178393025.714,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,secp256k1,256,ADD,50000,0.000009470,5279875219.433,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,secp256k1,256,ADD,50000,0.000284901,175499556.283,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,secp256k1,256,ADD,50000,0.000009480,5274299164.947,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,secp256k1,256,ADD,50000,0.000288411,173363680.456,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,secp256k1,256,ADD,50000,0.001121264,44592528.603,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,secp256k1,256,ADD,50000,0.003147717,15884528.459,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,secp256k1,256,ADD,50000,0.001895432,26379210.493,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,secp256k1,256,ADD,50000,0.004378309,11419934.023,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,secp256k1,256,ADD,50000,0.001950289,25637225.989,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,secp256k1,256,ADD,50000,0.004314240,11589526.906,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,secp256k1,256,ADD,50000,0.002221231,22510040.608,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,secp256k1,256,ADD,50000,0.004128179,12111877.762,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,secp256k1,256,SUBTRACT,50000,0.000329095,151931905.688,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,secp256k1,256,SUBTRACT,50000,0.000770650,64880279.127,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,secp256k1,256,SUBTRACT,50000,0.000728648,68620220.058,0
library,NVIDIA H200 NVL,gpu,cgbn,secp256k1,256,SUBTRACT,50000,0.000007936,6300403225.806,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,secp256k1,256,SUBTRACT,50000,0.000026950,1855285463.868,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,secp256k1,256,SUBTRACT,50000,0.000294611,169715267.220,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,secp256k1,256,SUBTRACT,50000,0.000017780,2812149163.551,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,secp256k1,256,SUBTRACT,50000,0.000299901,166721683.145,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32,secp256k1,256,SUBTRACT,50000,0.000012650,3952557260.521,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32,secp256k1,256,SUBTRACT,50000,0.000290281,172246878.217,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,secp256k1,256,SUBTRACT,50000,0.000012500,3999969542.258,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,secp256k1,256,SUBTRACT,50000,0.000292311,171050653.124,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,secp256k1,256,SUBTRACT,50000,0.000012810,3903202828.140,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,secp256k1,256,SUBTRACT,50000,0.000290860,171903965.829,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,secp256k1,256,SUBTRACT,50000,0.000009320,5364752614.946,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,secp256k1,256,SUBTRACT,50000,0.000292491,170945400.777,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,secp256k1,256,SUBTRACT,50000,0.000009270,5393790244.638,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,secp256k1,256,SUBTRACT,50000,0.000284171,175950417.901,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,secp256k1,256,SUBTRACT,50000,0.001066092,46900270.331,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,secp256k1,256,SUBTRACT,50000,0.003159477,15825403.883,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,secp256k1,256,SUBTRACT,50000,0.002010783,24865934.516,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,secp256k1,256,SUBTRACT,50000,0.004437191,11268390.431,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,secp256k1,256,SUBTRACT,50000,0.001879899,26597172.665,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,secp256k1,256,SUBTRACT,50000,0.004498211,11115530.230,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,secp256k1,256,SUBTRACT,50000,0.002483001,20136923.193,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,secp256k1,256,SUBTRACT,50000,0.005240045,9541902.749,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,secp256k1,256,ADDMOD,50000,0.001293050,38668255.407,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,secp256k1,256,ADDMOD,50000,0.002560684,19526032.655,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,secp256k1,256,ADDMOD,50000,0.000956713,52262275.666,0
library,NVIDIA H200 NVL,gpu,cgbn,secp256k1,256,ADDMOD,50000,0.000009376,5332764505.119,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,secp256k1,256,ADDMOD,50000,0.000037440,1335470713.013,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,secp256k1,256,ADDMOD,50000,0.000310481,161040479.248,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,secp256k1,256,ADDMOD,50000,0.000023000,2173918496.923,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,secp256k1,256,ADDMOD,50000,0.000290401,172175687.726,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32,secp256k1,256,ADDMOD,50000,0.000014030,3563785437.739,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32,secp256k1,256,ADDMOD,50000,0.000293031,170630419.484,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,secp256k1,256,ADDMOD,50000,0.000010750,4651152558.966,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,secp256k1,256,ADDMOD,50000,0.000276561,180791933.461,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,secp256k1,256,ADDMOD,50000,0.000010560,4734833310.550,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,secp256k1,256,ADDMOD,50000,0.000278001,179855498.763,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,secp256k1,256,ADDMOD,50000,0.000009020,5543180733.590,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,secp256k1,256,ADDMOD,50000,0.000284861,175524228.699,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,secp256k1,256,ADDMOD,50000,0.000008550,5847948499.537,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,secp256k1,256,ADDMOD,50000,0.000284131,175975145.155,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,secp256k1,256,ADDMOD,50000,0.001220466,40967958.577,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,secp256k1,256,ADDMOD,50000,0.003285760,15217179.309,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,secp256k1,256,ADDMOD,50000,0.001883923,26540363.417,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,secp256k1,256,ADDMOD,50000,0.004163349,12009562.358,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,secp256k1,256,ADDMOD,50000,0.002353031,21249188.164,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,secp256k1,256,ADDMOD,50000,0.005363325,9322575.071,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,secp256k1,256,ADDMOD,50000,0.002007670,24904491.250,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,secp256k1,256,ADDMOD,50000,0.004183190,11952600.595,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,secp256k1,256,SUBTRACTMOD,50000,0.001098547,45514680.064,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,secp256k1,256,SUBTRACTMOD,50000,0.000052107,959569809.179,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,secp256k1,256,SUBTRACTMOD,50000,0.001057420,47284900.827,0
library,NVIDIA H200 NVL,gpu,cgbn,secp256k1,256,SUBTRACTMOD,50000,0.000010048,4976114649.682,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.000038930,1284357127.564,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.000301461,165858951.084,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,secp256k1,256,SUBTRACTMOD,50000,0.000023960,2086819796.515,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,secp256k1,256,SUBTRACTMOD,50000,0.000291301,171643731.879,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32,secp256k1,256,SUBTRACTMOD,50000,0.000014720,3396761620.651,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32,secp256k1,256,SUBTRACTMOD,50000,0.000292530,170922610.906,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,secp256k1,256,SUBTRACTMOD,50000,0.000010970,4557913315.151,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,secp256k1,256,SUBTRACTMOD,50000,0.000282620,176916023.194,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,secp256k1,256,SUBTRACTMOD,50000,0.000010560,4734885508.604,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,secp256k1,256,SUBTRACTMOD,50000,0.000285910,174880150.036,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,secp256k1,256,SUBTRACTMOD,50000,0.000008760,5707749436.530,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,secp256k1,256,SUBTRACTMOD,50000,0.000293721,170229583.630,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,secp256k1,256,SUBTRACTMOD,50000,0.000008640,5787039756.390,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,secp256k1,256,SUBTRACTMOD,50000,0.000283701,176241890.903,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.000964881,51819859.717,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.003471874,14401444.568,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,secp256k1,256,SUBTRACTMOD,50000,0.001856193,26936854.865,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,secp256k1,256,SUBTRACTMOD,50000,0.004214858,11862795.700,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,secp256k1,256,SUBTRACTMOD,50000,0.001306746,38262981.882,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,secp256k1,256,SUBTRACTMOD,50000,0.005513446,9068738.597,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,secp256k1,256,SUBTRACTMOD,50000,0.002455741,20360453.867,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,secp256k1,256,SUBTRACTMOD,50000,0.005131384,9743959.898,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000532356,93922024.358,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001929416,25914583.892,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001707000,29291153.769,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000851102,58747361.252,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001186724,42132798.058,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000243501,205338030.853,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000553362,90356752.972,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000067971,735607132.446,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000406811,122907228.688,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000019930,2508771887.522,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000348011,143673606.406,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000018750,2666671196.627,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000359601,139042997.991,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000010860,4604036249.424,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000341341,146481078.082,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000009580,5219243047.235,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000331001,151056942.254,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001254857,39845175.923,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.003540065,14124034.353,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.002044794,24452340.522,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.004325190,11560185.857,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001250726,39976779.746,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.004331200,11544144.718,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.002345811,21314589.136,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.004367571,11448010.566,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000527015,94873888.250,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000780414,64068541.420,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000810130,61718460.455,0
library,NVIDIA H200 NVL,gpu,cgbn,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000007712,6483402489.627,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000185950,268889374.876,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000507432,98535363.246,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000060981,819927207.473,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000383931,130231741.371,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000022540,2218280055.987,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000348461,143488041.993,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000023510,2126747856.400,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000349871,142909815.600,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000021470,2328829609.708,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000342601,145942323.834,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000018230,2742723136.754,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000348441,143496335.587,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000017240,2900221685.315,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000362251,138025861.519,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001348739,37071664.051,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.003590566,13925381.164,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.002559807,19532722.316,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.004454301,11225105.580,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.002463511,20296235.121,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.004683692,10675339.009,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.002028159,24652899.717,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.004804493,10406925.507,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.005263825,9498796.123,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.004984315,10031468.679,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.004641294,10772857.710,0
library,NVIDIA H200 NVL,gpu,cgbn,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000010944,4568713450.292,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000082530,605840558.478,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000359511,139077801.783,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000026710,1871959316.065,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000291611,171461256.088,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000012450,4016090005.984,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000285671,175026531.127,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000014281,3501151268.820,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000283821,176167360.513,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000011930,4191111551.748,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000284181,175944219.154,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000013460,3714694818.415,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000283210,176547465.569,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000010720,4664184110.160,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000276551,180798478.503,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.002009072,24887112.162,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.004146538,12058252.037,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.001717351,29114606.852,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.004239999,11792455.637,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.002152471,23229116.803,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.004094239,12212281.771,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.002018400,24772096.682,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.004336361,11530405.444,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,secp256k1,256,COMPARE,50000,0.000173890,287537276.094,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,secp256k1,256,COMPARE,50000,0.000677072,73847390.286,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,secp256k1,256,COMPARE,50000,0.000664865,75203199.823,0
library,NVIDIA H200 NVL,gpu,cgbn,secp256k1,256,COMPARE,50000,0.000008000,6250000000.000,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,secp256k1,256,COMPARE,50000,0.000027421,1823416882.549,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,secp256k1,256,COMPARE,50000,0.000304681,164106056.343,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,secp256k1,256,COMPARE,50000,0.000017920,2790189952.641,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,secp256k1,256,COMPARE,50000,0.000283791,176186005.269,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,secp256k1,256,COMPARE,50000,0.000010170,4916456571.160,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,secp256k1,256,COMPARE,50000,0.000284511,175740192.886,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,secp256k1,256,COMPARE,50000,0.000010330,4840272382.627,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,secp256k1,256,COMPARE,50000,0.000282940,176715918.296,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,secp256k1,256,COMPARE,50000,0.000009430,5302232381.517,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,secp256k1,256,COMPARE,50000,0.000285610,175063913.937,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,secp256k1,256,COMPARE,50000,0.000009030,5537034983.498,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,secp256k1,256,COMPARE,50000,0.000277500,180180244.376,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,secp256k1,256,COMPARE,50000,0.000948510,52714262.661,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,secp256k1,256,COMPARE,50000,0.003861812,12947289.828,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,secp256k1,256,COMPARE,50000,0.001889552,26461297.276,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,secp256k1,256,COMPARE,50000,0.004151419,12044074.458,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,secp256k1,256,COMPARE,50000,0.002458631,20336521.067,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,secp256k1,256,COMPARE,50000,0.004965193,10070102.033,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,secp256k1,256,REDUCE,6250,0.000061989,100823892.420,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,secp256k1,256,REDUCE,6250,0.000094860,65886801.104,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,secp256k1,256,REDUCE,6250,0.000119416,52338137.559,0
library,NVIDIA H200 NVL,gpu,cgbn,secp256k1,256,REDUCE,50000,0.000013696,3650700934.579,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,secp256k1,256,REDUCE,50000,0.000150891,331364979.200,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,secp256k1,256,REDUCE,50000,0.000416801,119961335.455,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,secp256k1,256,REDUCE,50000,0.000101280,493681211.235,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,secp256k1,256,REDUCE,50000,0.000364631,137124904.531,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,secp256k1,256,REDUCE,50000,0.000036870,1356119887.594,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,secp256k1,256,REDUCE,50000,0.000315601,158427885.291,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,secp256k1,256,REDUCE,50000,0.000036710,1362028850.496,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,secp256k1,256,REDUCE,50000,0.000318400,157035185.875,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,secp256k1,256,REDUCE,50000,0.000035951,1390785223.564,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,secp256k1,256,REDUCE,50000,0.000317840,157311901.734,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,secp256k1,256,REDUCE,50000,0.000034600,1445085425.891,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,secp256k1,256,REDUCE,50000,0.000307841,162421513.167,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,secp256k1,256,REDUCE,50000,0.001610474,31046760.869,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,secp256k1,256,REDUCE,50000,0.003617067,13823354.644,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,secp256k1,256,REDUCE,50000,0.001905803,26235660.429,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,secp256k1,256,REDUCE,50000,0.003928596,12727193.380,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,secp256k1,256,REDUCE,50000,0.001777149,28134952.216,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,secp256k1,256,REDUCE,50000,0.003781918,13220804.996,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,secp256k1,256,MODMUL,3125,0.000177000,17655358.594,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,secp256k1,256,MODMUL,3125,0.000023086,135363893.568,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,secp256k1,256,MODMUL,3125,0.000157966,19782696.499,0
library,NVIDIA H200 NVL,gpu,cgbn,secp256k1,256,MODMUL,50000,0.000035680,1401345291.480,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,secp256k1,256,MODMUL,50000,0.000409001,122249084.221,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,secp256k1,256,MODMUL,50000,0.000686722,72809667.570,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,secp256k1,256,MODMUL,50000,0.000239951,208375815.607,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,secp256k1,256,MODMUL,50000,0.000519541,96238799.881,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,secp256k1,256,MODMUL,50000,0.000112490,444483835.568,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,secp256k1,256,MODMUL,50000,0.000382751,130633233.713,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,secp256k1,256,MODMUL,50000,0.000087991,568239884.525,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,secp256k1,256,MODMUL,50000,0.000365371,136847204.265,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,secp256k1,256,MODMUL,50000,0.000111070,450166578.205,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,secp256k1,256,MODMUL,50000,0.000390841,127929284.085,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,secp256k1,256,MODMUL,50000,0.000090640,551633308.887,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,secp256k1,256,MODMUL,50000,0.000357661,139797163.870,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,secp256k1,256,MODMUL,50000,0.005020137,9959887.427,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,secp256k1,256,MODMUL,50000,0.006767804,7387920.769,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,secp256k1,256,MODMUL,50000,0.003817956,13096012.357,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,secp256k1,256,MODMUL,50000,0.006255813,7992566.347,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,secp256k1,256,MODMUL,50000,0.002167040,23072947.255,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,secp256k1,256,MODMUL,50000,0.004306920,11609224.254,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,secp256k1,256,MODEXP,781,0.004981375,156784.022,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,secp256k1,256,MODEXP,781,0.003022429,258401.432,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,secp256k1,256,MODEXP,781,0.002907168,268646.325,0
library,NVIDIA H200 NVL,gpu,cgbn,secp256k1,256,MODEXP,50000,0.010293600,4857387.114,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,secp256k1,256,MODEXP,50000,0.018673935,2677528.870,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,secp256k1,256,MODEXP,50000,0.018956015,2637685.188,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,secp256k1,256,MODEXP,50000,0.003527112,14175903.863,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,secp256k1,256,MODEXP,50000,0.003789774,13193399.183,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,secp256k1,256,MODEXP,50000,0.002819847,17731458.692,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,secp256k1,256,MODEXP,50000,0.003110948,16072271.121,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,secp256k1,256,MODEXP,50000,0.001339933,37315298.042,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,secp256k1,256,MODEXP,50000,0.001618704,30888907.511,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,secp256k1,256,MODEXP,50000,0.002774687,18020050.431,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,secp256k1,256,MODEXP,50000,0.003055168,16365711.687,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,secp256k1,256,MODEXP,50000,0.001340163,37308892.936,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,secp256k1,256,MODEXP,50000,0.001611314,31030576.990,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,secp256k1,256,MODEXP,50000,0.887418094,56343.228,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,secp256k1,256,MODEXP,50000,0.887733641,56323.201,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,secp256k1,256,MODEXP,50000,0.201563744,248060.485,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,secp256k1,256,MODEXP,50000,0.199680451,250400.076,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,secp256k1,256,MODEXP,50000,0.105073005,475859.617,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,secp256k1,256,MODEXP,50000,0.104616612,477935.569,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,secp256k1,256,EXPONENTIATION,781,0.001633088,478235.102,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,secp256k1,256,EXPONENTIATION,781,0.001065700,732851.884,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,secp256k1,256,EXPONENTIATION,781,0.001390357,561726.109,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,secp256k1,256,EXPONENTIATION,50000,0.027275740,1833130.838,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,secp256k1,256,EXPONENTIATION,50000,0.027543141,1815333.993,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,secp256k1,256,EXPONENTIATION,50000,0.007062775,7079370.394,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,secp256k1,256,EXPONENTIATION,50000,0.007310786,6839209.962,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,secp256k1,256,EXPONENTIATION,50000,0.000510601,97923801.160,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,secp256k1,256,EXPONENTIATION,50000,0.000779012,64183864.875,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,secp256k1,256,EXPONENTIATION,50000,0.000350931,142478166.871,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,secp256k1,256,EXPONENTIATION,50000,0.000629072,79482162.753,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,secp256k1,256,EXPONENTIATION,50000,0.000507332,98554785.657,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,secp256k1,256,EXPONENTIATION,50000,0.000792012,63130359.173,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,secp256k1,256,EXPONENTIATION,50000,0.000366221,136529600.796,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,secp256k1,256,EXPONENTIATION,50000,0.000636061,78608818.238,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,secp256k1,256,EXPONENTIATION,50000,0.316363839,158045.876,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,secp256k1,256,EXPONENTIATION,50000,0.317888701,157287.755,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,secp256k1,256,EXPONENTIATION,50000,0.096809992,516475.613,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,secp256k1,256,EXPONENTIATION,50000,0.100304604,498481.605,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,secp256k1,256,EXPONENTIATION,50000,0.013784425,3627282.245,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,secp256k1,256,EXPONENTIATION,50000,0.017124801,2919741.948,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,secp256k1,256,DIVIDE,6250,0.000125886,49647975.212,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,secp256k1,256,DIVIDE,6250,0.000753255,8297320.747,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,secp256k1,256,DIVIDE,6250,0.000364350,17153831.041,0
library,NVIDIA H200 NVL,gpu,cgbn,secp256k1,256,DIVIDE,50000,0.000016832,2970532319.392,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,secp256k1,256,DIVIDE,50000,0.000356471,140263890.672,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,secp256k1,256,DIVIDE,50000,0.000687072,72772583.623,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,secp256k1,256,DIVIDE,50000,0.000299821,166766156.191,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,secp256k1,256,DIVIDE,50000,0.000622182,80362347.998,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,secp256k1,256,DIVIDE,50000,0.000108681,460062352.152,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,secp256k1,256,DIVIDE,50000,0.000430941,116025163.632,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,secp256k1,256,DIVIDE,50000,0.000105940,471965611.192,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,secp256k1,256,DIVIDE,50000,0.000438331,114069036.753,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,secp256k1,256,DIVIDE,50000,0.000108140,462363367.194,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,secp256k1,256,DIVIDE,50000,0.000429651,116373521.368,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,secp256k1,256,DIVIDE,50000,0.000103141,484772967.183,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,secp256k1,256,DIVIDE,50000,0.000404821,123511375.959,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,secp256k1,256,DIVIDE,50000,0.002713948,18423344.730,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,secp256k1,256,DIVIDE,50000,0.005323413,9392470.680,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,secp256k1,256,DIVIDE,50000,0.002745538,18211366.829,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,secp256k1,256,DIVIDE,50000,0.004844433,10321125.450,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,secp256k1,256,DIVIDE,50000,0.001963849,25460205.136,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,secp256k1,256,DIVIDE,50000,0.004441991,11256213.460,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,secp256k1,256,ISQRT,1562,0.000073892,21138995.260,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,secp256k1,256,ISQRT,1562,0.000083011,18816670.012,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,secp256k1,256,ISQRT,50000,0.004060782,12312899.258,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,secp256k1,256,ISQRT,50000,0.004342243,11514786.361,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,secp256k1,256,ISQRT,50000,0.003012921,16595191.384,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,secp256k1,256,ISQRT,50000,0.003291042,15192756.859,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,secp256k1,256,ISQRT,50000,0.000731491,68353539.228,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,secp256k1,256,ISQRT,50000,0.001005263,49738230.114,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,secp256k1,256,ISQRT,50000,0.000667412,74916249.205,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,secp256k1,256,ISQRT,50000,0.000939163,53238899.524,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,secp256k1,256,ISQRT,50000,0.000743322,67265603.366,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,secp256k1,256,ISQRT,50000,0.001030303,48529415.275,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,secp256k1,256,ISQRT,50000,0.000659562,75807877.125,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,secp256k1,256,ISQRT,50000,0.000923543,54139330.563,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,secp256k1,256,ISQRT,50000,0.017857919,2799878.309,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,secp256k1,256,ISQRT,50000,0.021534908,2321811.639,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,secp256k1,256,ISQRT,50000,0.018108843,2761081.974,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,secp256k1,256,ISQRT,50000,0.018003102,2777299.168,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,secp256k1,256,ISQRT,50000,0.009337304,5354864.786,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,secp256k1,256,ISQRT,50000,0.011108442,4501081.259,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,secp256k1,256,MODMUL_R2,50000,0.002878323,17371223.284,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,secp256k1,256,MODMUL_R2,50000,0.002101696,23790310.529,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,secp256k1,256,MODMUL_R2,50000,0.002155502,23196457.866,0
library,NVIDIA H200 NVL,gpu,cgbn,secp256k1,256,MODMUL_R2,50000,0.000015808,3162955465.587,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,secp256k1,256,MODMUL_R2,50000,0.000098710,506534557.360,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,secp256k1,256,MODMUL_R2,50000,0.000369211,135423932.035,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,secp256k1,256,MODMUL_R2,50000,0.000028890,1730697158.331,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,secp256k1,256,MODMUL_R2,50000,0.000290671,172015775.749,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,secp256k1,256,MODMUL_R2,50000,0.000018840,2653933843.321,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,secp256k1,256,MODMUL_R2,50000,0.000300640,166311863.876,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,secp256k1,256,MODMUL_R2,50000,0.000013850,3610126331.008,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,secp256k1,256,MODMUL_R2,50000,0.000284241,175907035.844,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,secp256k1,256,MODMUL_R2,50000,0.000018630,2683851337.874,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,secp256k1,256,MODMUL_R2,50000,0.000284581,175696914.420,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,secp256k1,256,MODMUL_R2,50000,0.000013270,3767911794.224,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,secp256k1,256,MODMUL_R2,50000,0.000278781,179352293.013,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,secp256k1,256,MODMUL_R2,50000,0.002368751,21108170.669,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,secp256k1,256,MODMUL_R2,50000,0.004721891,10588977.976,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,secp256k1,256,MODMUL_R2,50000,0.001983484,25208168.800,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,secp256k1,256,MODMUL_R2,50000,0.005278406,9472556.738,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,secp256k1,256,MODMUL_R2,50000,0.002038449,24528453.571,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,secp256k1,256,MODMUL_R2,50000,0.005490075,9107343.737,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,ADD,50000,0.000538620,92829891.980,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,ADD,50000,0.000679132,73623390.605,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,ADD,50000,0.000788186,63436785.995,0
library,NVIDIA H200 NVL,gpu,cgbn,rsa256(composite),256,ADD,50000,0.000007936,6300403225.806,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,rsa256(composite),256,ADD,50000,0.000027471,1820101916.313,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,rsa256(composite),256,ADD,50000,0.000295800,169033173.889,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,rsa256(composite),256,ADD,50000,0.000017650,2832867646.360,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,rsa256(composite),256,ADD,50000,0.000283131,176596682.424,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32,rsa256(composite),256,ADD,50000,0.000013240,3776426212.730,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32,rsa256(composite),256,ADD,50000,0.000297111,168287276.836,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,rsa256(composite),256,ADD,50000,0.000012640,3955687941.277,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,rsa256(composite),256,ADD,50000,0.000283090,176622317.993,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,rsa256(composite),256,ADD,50000,0.000012690,3940120081.463,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,rsa256(composite),256,ADD,50000,0.000309011,161806554.590,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,rsa256(composite),256,ADD,50000,0.000009680,5165324469.032,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,rsa256(composite),256,ADD,50000,0.000284691,175629020.269,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,rsa256(composite),256,ADD,50000,0.000009610,5202930739.318,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,rsa256(composite),256,ADD,50000,0.000277430,180225608.764,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,rsa256(composite),256,ADD,50000,0.002282489,21905910.851,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,rsa256(composite),256,ADD,50000,0.004831232,10349327.265,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,rsa256(composite),256,ADD,50000,0.002239055,22330849.396,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,rsa256(composite),256,ADD,50000,0.004177329,11969370.818,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,rsa256(composite),256,ADD,50000,0.002008700,24891720.458,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,rsa256(composite),256,ADD,50000,0.004478521,11164400.166,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,rsa256(composite),256,ADD,50000,0.002004580,24942880.925,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,rsa256(composite),256,ADD,50000,0.004705872,10625023.553,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACT,50000,0.000331467,150844581.414,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACT,50000,0.000779218,64166861.524,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACT,50000,0.000764650,65389376.775,0
library,NVIDIA H200 NVL,gpu,cgbn,rsa256(composite),256,SUBTRACT,50000,0.000007904,6325910931.174,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,rsa256(composite),256,SUBTRACT,50000,0.000027031,1849724281.746,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,rsa256(composite),256,SUBTRACT,50000,0.000301900,165617706.538,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,rsa256(composite),256,SUBTRACT,50000,0.000017650,2832867646.360,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,rsa256(composite),256,SUBTRACT,50000,0.000284781,175573522.575,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32,rsa256(composite),256,SUBTRACT,50000,0.000012940,3863979070.479,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32,rsa256(composite),256,SUBTRACT,50000,0.000289670,172610166.756,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,rsa256(composite),256,SUBTRACT,50000,0.000012490,4003175811.127,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,rsa256(composite),256,SUBTRACT,50000,0.000282620,176916023.194,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,rsa256(composite),256,SUBTRACT,50000,0.000012780,3912376042.777,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,rsa256(composite),256,SUBTRACT,50000,0.000308901,161864119.633,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,rsa256(composite),256,SUBTRACT,50000,0.000009310,5370588826.089,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,rsa256(composite),256,SUBTRACT,50000,0.000292410,170992768.722,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,rsa256(composite),256,SUBTRACT,50000,0.000009280,5387903526.312,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,rsa256(composite),256,SUBTRACT,50000,0.000284101,175993749.244,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,rsa256(composite),256,SUBTRACT,50000,0.001915151,26107601.680,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,rsa256(composite),256,SUBTRACT,50000,0.004363043,11459891.637,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,rsa256(composite),256,SUBTRACT,50000,0.002126065,23517625.656,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,rsa256(composite),256,SUBTRACT,50000,0.004084268,12242095.760,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,rsa256(composite),256,SUBTRACT,50000,0.002346111,21311864.651,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,rsa256(composite),256,SUBTRACT,50000,0.005459236,9158790.749,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,rsa256(composite),256,SUBTRACT,50000,0.002667922,18741177.980,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,rsa256(composite),256,SUBTRACT,50000,0.005341465,9360727.810,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,ADDMOD,50000,0.001118737,44693270.672,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,ADDMOD,50000,0.001201524,41613828.828,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,ADDMOD,50000,0.000915666,54605061.438,0
library,NVIDIA H200 NVL,gpu,cgbn,rsa256(composite),256,ADDMOD,50000,0.000010112,4944620253.165,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,rsa256(composite),256,ADDMOD,50000,0.000034250,1459855303.615,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,rsa256(composite),256,ADDMOD,50000,0.000302401,165343356.515,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,rsa256(composite),256,ADDMOD,50000,0.000021750,2298839222.403,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,rsa256(composite),256,ADDMOD,50000,0.000291901,171290917.046,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32,rsa256(composite),256,ADDMOD,50000,0.000013400,3731347288.128,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32,rsa256(composite),256,ADDMOD,50000,0.000276791,180641680.136,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,rsa256(composite),256,ADDMOD,50000,0.000010830,4616804755.506,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,rsa256(composite),256,ADDMOD,50000,0.000276721,180687353.219,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,rsa256(composite),256,ADDMOD,50000,0.000010270,4868528657.092,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,rsa256(composite),256,ADDMOD,50000,0.000291741,171384900.137,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,rsa256(composite),256,ADDMOD,50000,0.000008680,5760340252.947,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,rsa256(composite),256,ADDMOD,50000,0.000292300,171057124.992,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,rsa256(composite),256,ADDMOD,50000,0.000008590,5820764221.340,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,rsa256(composite),256,ADDMOD,50000,0.000283721,176229452.741,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,rsa256(composite),256,ADDMOD,50000,0.001886400,26505514.152,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,rsa256(composite),256,ADDMOD,50000,0.004297941,11633477.541,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,rsa256(composite),256,ADDMOD,50000,0.002120965,23574174.474,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,rsa256(composite),256,ADDMOD,50000,0.004216809,11857307.430,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,rsa256(composite),256,ADDMOD,50000,0.001878549,26616287.420,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,rsa256(composite),256,ADDMOD,50000,0.004217930,11854155.888,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,rsa256(composite),256,ADDMOD,50000,0.002583452,19353949.416,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,rsa256(composite),256,ADDMOD,50000,0.005467636,9144720.068,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACTMOD,50000,0.001115213,44834485.592,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACTMOD,50000,0.001908240,26202148.129,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACTMOD,50000,0.001956430,25556747.396,0
library,NVIDIA H200 NVL,gpu,cgbn,rsa256(composite),256,SUBTRACTMOD,50000,0.000010080,4960317460.317,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,rsa256(composite),256,SUBTRACTMOD,50000,0.000038460,1300055482.371,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,rsa256(composite),256,SUBTRACTMOD,50000,0.000306761,162993335.476,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,rsa256(composite),256,SUBTRACTMOD,50000,0.000023280,2147773597.436,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,rsa256(composite),256,SUBTRACTMOD,50000,0.000284121,175981346.081,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32,rsa256(composite),256,SUBTRACTMOD,50000,0.000013660,3660337909.288,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32,rsa256(composite),256,SUBTRACTMOD,50000,0.000293621,170287560.032,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,rsa256(composite),256,SUBTRACTMOD,50000,0.000010370,4821580295.920,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,rsa256(composite),256,SUBTRACTMOD,50000,0.000292951,170677002.579,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,rsa256(composite),256,SUBTRACTMOD,50000,0.000010710,4668544203.135,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,rsa256(composite),256,SUBTRACTMOD,50000,0.000295560,169170392.947,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,rsa256(composite),256,SUBTRACTMOD,50000,0.000008840,5656110220.583,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,rsa256(composite),256,SUBTRACTMOD,50000,0.000291311,171637901.454,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,rsa256(composite),256,SUBTRACTMOD,50000,0.000008600,5813909218.399,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,rsa256(composite),256,SUBTRACTMOD,50000,0.000276811,180628689.185,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,rsa256(composite),256,SUBTRACTMOD,50000,0.002086854,23959510.830,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,rsa256(composite),256,SUBTRACTMOD,50000,0.004242361,11785889.820,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,rsa256(composite),256,SUBTRACTMOD,50000,0.002158754,23161508.483,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,rsa256(composite),256,SUBTRACTMOD,50000,0.004524541,11050844.502,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,rsa256(composite),256,SUBTRACTMOD,50000,0.002307711,21666490.673,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,rsa256(composite),256,SUBTRACTMOD,50000,0.005181305,9650078.514,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,rsa256(composite),256,SUBTRACTMOD,50000,0.002280371,21926256.440,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,rsa256(composite),256,SUBTRACTMOD,50000,0.004856813,10294816.785,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000551066,90733265.912,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000691224,72335453.124,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000819450,61016505.154,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000840722,59472691.248,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001172324,42650327.372,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000238501,209642725.855,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000553322,90363292.573,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000068490,730034351.140,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000398671,125416690.154,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000019290,2592014059.143,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000340811,146708889.636,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000018230,2742723136.754,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000340161,146989206.769,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000010860,4604036249.424,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000340320,146920522.268,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000009510,5257577084.379,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000331881,150656414.329,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.002055314,24327182.857,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.004199389,11906493.962,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001994374,25070523.692,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.004462871,11203550.355,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001586867,31508626.137,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.004437481,11267654.039,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.002064069,24223995.982,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.004229710,11821141.481,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000542312,92197920.626,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000669990,74627980.524,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000770308,64909069.273,0
library,NVIDIA H200 NVL,gpu,cgbn,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000007872,6351626016.260,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000185341,269773024.376,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000506422,98731884.366,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000060400,827814648.287,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000372361,134278321.092,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000022470,2225187184.482,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000349631,143007886.504,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000022730,2199738434.512,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000349401,143102039.228,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000021030,2377546981.683,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000349231,143171685.192,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000018170,2751790629.104,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000349471,143073389.579,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000017400,2873560563.343,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000343301,145644741.320,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.002319899,21552663.147,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.004488115,11140534.585,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001952613,25606712.575,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.004343530,11511374.574,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.002198000,22747951.760,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.004569761,10941491.276,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.002869994,17421639.025,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.004365841,11452547.365,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.005276665,9475681.919,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.003415890,14637473.526,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.003536620,14137792.803,0
library,NVIDIA H200 NVL,gpu,cgbn,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000010816,4622781065.089,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000081651,612362152.737,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000358331,139535785.630,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000027010,1851167298.525,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000298211,167666516.995,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000012410,4029012200.636,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000276691,180707043.084,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000014190,3523613142.890,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000286021,174812387.475,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000012200,4098368556.352,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000283800,176180368.066,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000013530,3695485618.902,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000308121,162273865.057,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000010770,4642555420.319,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000268230,186407200.451,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.002354330,21237464.280,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.004669379,10708062.088,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.001864952,26810340.993,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.004114208,12153007.203,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.002121300,23570452.405,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.004211040,11873551.481,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.002066760,24192454.758,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.004180739,11959607.950,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,COMPARE,50000,0.000201716,247873664.228,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,COMPARE,50000,0.000677127,73841391.312,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,COMPARE,50000,0.000565700,88386070.332,0
library,NVIDIA H200 NVL,gpu,cgbn,rsa256(composite),256,COMPARE,50000,0.000007936,6300403225.806,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,rsa256(composite),256,COMPARE,50000,0.000026190,1909120418.187,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,rsa256(composite),256,COMPARE,50000,0.000294691,169669207.550,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,rsa256(composite),256,COMPARE,50000,0.000018140,2756346895.476,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,rsa256(composite),256,COMPARE,50000,0.000288751,173159587.671,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,rsa256(composite),256,COMPARE,50000,0.000010200,4901979405.823,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,rsa256(composite),256,COMPARE,50000,0.000283541,176341387.043,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,rsa256(composite),256,COMPARE,50000,0.000010400,4807709515.867,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,rsa256(composite),256,COMPARE,50000,0.000285180,175327902.250,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,rsa256(composite),256,COMPARE,50000,0.000009240,5411255113.329,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,rsa256(composite),256,COMPARE,50000,0.000291921,171279236.173,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,rsa256(composite),256,COMPARE,50000,0.000009070,5512658412.804,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,rsa256(composite),256,COMPARE,50000,0.000284240,175907684.255,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,rsa256(composite),256,COMPARE,50000,0.002093244,23886369.733,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,rsa256(composite),256,COMPARE,50000,0.004130798,12104198.628,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,rsa256(composite),256,COMPARE,50000,0.002248115,22240854.892,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,rsa256(composite),256,COMPARE,50000,0.004697552,10643841.840,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,rsa256(composite),256,COMPARE,50000,0.002623932,19055372.076,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,rsa256(composite),256,COMPARE,50000,0.005553546,9003256.827,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,REDUCE,6250,0.000099735,62666086.519,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,REDUCE,6250,0.000125397,49841690.168,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,REDUCE,6250,0.000199668,31301976.953,0
library,NVIDIA H200 NVL,gpu,cgbn,rsa256(composite),256,REDUCE,50000,0.000013536,3693853427.896,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,rsa256(composite),256,REDUCE,50000,0.000151560,329902472.252,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,rsa256(composite),256,REDUCE,50000,0.000432052,115726823.769,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,rsa256(composite),256,REDUCE,50000,0.000101660,491835952.967,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,rsa256(composite),256,REDUCE,50000,0.000373012,134043931.665,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,rsa256(composite),256,REDUCE,50000,0.000035910,1392367787.593,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,rsa256(composite),256,REDUCE,50000,0.000315971,158242383.344,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,rsa256(composite),256,REDUCE,50000,0.000036310,1377029023.953,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,rsa256(composite),256,REDUCE,50000,0.000317391,157534393.596,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,rsa256(composite),256,REDUCE,50000,0.000036310,1377033438.923,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,rsa256(composite),256,REDUCE,50000,0.000316071,158192317.504,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,rsa256(composite),256,REDUCE,50000,0.000034550,1447179173.939,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,rsa256(composite),256,REDUCE,50000,0.000300910,166162654.011,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,rsa256(composite),256,REDUCE,50000,0.002036993,24545986.075,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,rsa256(composite),256,REDUCE,50000,0.003845452,13002372.735,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,rsa256(composite),256,REDUCE,50000,0.001917073,26081428.333,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,rsa256(composite),256,REDUCE,50000,0.003974437,12580398.132,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,rsa256(composite),256,REDUCE,50000,0.002270321,22023317.287,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,rsa256(composite),256,REDUCE,50000,0.004786992,10444972.586,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MODMUL,3125,0.000177918,17564225.970,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MODMUL,3125,0.000232903,13417617.532,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MODMUL,3125,0.000249878,12506090.651,0
library,NVIDIA H200 NVL,gpu,cgbn,rsa256(composite),256,MODMUL,50000,0.000035776,1397584973.166,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,rsa256(composite),256,MODMUL,50000,0.000410201,121891454.649,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,rsa256(composite),256,MODMUL,50000,0.000679242,73611471.343,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,rsa256(composite),256,MODMUL,50000,0.000239671,208619235.957,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,rsa256(composite),256,MODMUL,50000,0.000507442,98533419.166,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,rsa256(composite),256,MODMUL,50000,0.000112631,443927480.199,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,rsa256(composite),256,MODMUL,50000,0.000383071,130524099.916,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,rsa256(composite),256,MODMUL,50000,0.000087951,568498622.224,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,rsa256(composite),256,MODMUL,50000,0.000357461,139875426.745,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,rsa256(composite),256,MODMUL,50000,0.000110740,451507731.511,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,rsa256(composite),256,MODMUL,50000,0.000394231,126829175.480,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,rsa256(composite),256,MODMUL,50000,0.000090801,550654482.003,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,rsa256(composite),256,MODMUL,50000,0.000366391,136466222.240,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,rsa256(composite),256,MODMUL,50000,0.005168810,9673406.412,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,rsa256(composite),256,MODMUL,50000,0.007114751,7027652.799,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,rsa256(composite),256,MODMUL,50000,0.003935797,12703907.163,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,rsa256(composite),256,MODMUL,50000,0.006277983,7964341.464,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,rsa256(composite),256,MODMUL,50000,0.002088660,23938793.006,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,rsa256(composite),256,MODMUL,50000,0.004556252,10973931.989,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MODEXP,781,0.004550953,171612.408,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MODEXP,781,0.001764970,442500.450,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MODEXP,781,0.001188778,656976.892,0
library,NVIDIA H200 NVL,gpu,cgbn,rsa256(composite),256,MODEXP,50000,0.010100096,4950447.996,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,rsa256(composite),256,MODEXP,50000,0.018682055,2676365.108,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,rsa256(composite),256,MODEXP,50000,0.018970086,2635728.693,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,rsa256(composite),256,MODEXP,50000,0.003529022,14168231.384,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,rsa256(composite),256,MODEXP,50000,0.003804624,13141903.239,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,rsa256(composite),256,MODEXP,50000,0.002797657,17872097.706,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,rsa256(composite),256,MODEXP,50000,0.003076948,16249868.415,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,rsa256(composite),256,MODEXP,50000,0.001336833,37401829.647,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,rsa256(composite),256,MODEXP,50000,0.001611985,31017657.711,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,rsa256(composite),256,MODEXP,50000,0.002774537,18021025.036,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,rsa256(composite),256,MODEXP,50000,0.003060918,16334968.425,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,rsa256(composite),256,MODEXP,50000,0.001330803,37571300.032,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,rsa256(composite),256,MODEXP,50000,0.001604184,31168493.558,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,rsa256(composite),256,MODEXP,50000,0.815781681,61290.908,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,rsa256(composite),256,MODEXP,50000,0.816662769,61224.782,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,rsa256(composite),256,MODEXP,50000,0.200037124,249953.604,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,rsa256(composite),256,MODEXP,50000,0.200709638,249116.089,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,rsa256(composite),256,MODEXP,50000,0.104126710,480184.191,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,rsa256(composite),256,MODEXP,50000,0.102478252,487908.400,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,EXPONENTIATION,781,0.001628358,479624.250,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,EXPONENTIATION,781,0.001238177,630766.044,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,EXPONENTIATION,781,0.002612451,298952.974,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,rsa256(composite),256,EXPONENTIATION,50000,0.027281480,1832745.143,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,rsa256(composite),256,EXPONENTIATION,50000,0.027482471,1819341.502,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,rsa256(composite),256,EXPONENTIATION,50000,0.007044995,7097237.146,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,rsa256(composite),256,EXPONENTIATION,50000,0.007334865,6816758.001,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,rsa256(composite),256,EXPONENTIATION,50000,0.000511881,97678937.372,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,rsa256(composite),256,EXPONENTIATION,50000,0.000787062,63527397.763,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,rsa256(composite),256,EXPONENTIATION,50000,0.000349050,143245985.280,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,rsa256(composite),256,EXPONENTIATION,50000,0.000623722,80163921.929,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,rsa256(composite),256,EXPONENTIATION,50000,0.000506241,98767189.727,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,rsa256(composite),256,EXPONENTIATION,50000,0.000785302,63669771.156,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,rsa256(composite),256,EXPONENTIATION,50000,0.000365641,136746164.583,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,rsa256(composite),256,EXPONENTIATION,50000,0.000629642,79410198.603,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,rsa256(composite),256,EXPONENTIATION,50000,0.321249023,155642.497,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,rsa256(composite),256,EXPONENTIATION,50000,0.321419686,155559.856,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,rsa256(composite),256,EXPONENTIATION,50000,0.095573861,523155.594,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,rsa256(composite),256,EXPONENTIATION,50000,0.099352527,503258.463,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,rsa256(composite),256,EXPONENTIATION,50000,0.013648904,3663297.788,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,rsa256(composite),256,EXPONENTIATION,50000,0.016915870,2955804.221,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,DIVIDE,6250,0.000128825,48515285.916,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,DIVIDE,6250,0.000256499,24366578.569,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,DIVIDE,6250,0.000240625,25974019.046,0
library,NVIDIA H200 NVL,gpu,cgbn,rsa256(composite),256,DIVIDE,50000,0.000016896,2959280303.030,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,rsa256(composite),256,DIVIDE,50000,0.000358662,139407023.766,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,rsa256(composite),256,DIVIDE,50000,0.000689502,72516106.819,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,rsa256(composite),256,DIVIDE,50000,0.000305751,163531768.468,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,rsa256(composite),256,DIVIDE,50000,0.000634362,78819350.472,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,rsa256(composite),256,DIVIDE,50000,0.000113661,439904428.473,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,rsa256(composite),256,DIVIDE,50000,0.000433991,115209771.794,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,rsa256(composite),256,DIVIDE,50000,0.000109800,455373318.419,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,rsa256(composite),256,DIVIDE,50000,0.000432531,115598650.702,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,rsa256(composite),256,DIVIDE,50000,0.000112180,445712184.587,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,rsa256(composite),256,DIVIDE,50000,0.000431261,115939064.879,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,rsa256(composite),256,DIVIDE,50000,0.000106890,467770660.591,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,rsa256(composite),256,DIVIDE,50000,0.000420881,118798433.785,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,rsa256(composite),256,DIVIDE,50000,0.002884181,17335943.295,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,rsa256(composite),256,DIVIDE,50000,0.005269342,9488850.606,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,rsa256(composite),256,DIVIDE,50000,0.002796609,17878795.645,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,rsa256(composite),256,DIVIDE,50000,0.004950494,10100002.298,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,rsa256(composite),256,DIVIDE,50000,0.002144210,23318611.075,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,rsa256(composite),256,DIVIDE,50000,0.004889724,10225526.103,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,ISQRT,1562,0.000073698,21194688.694,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,ISQRT,1562,0.000081405,19188050.023,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,rsa256(composite),256,ISQRT,50000,0.004058071,12321124.797,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,rsa256(composite),256,ISQRT,50000,0.004335122,11533700.914,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,rsa256(composite),256,ISQRT,50000,0.003090801,16177035.960,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,rsa256(composite),256,ISQRT,50000,0.003365762,14855476.715,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,rsa256(composite),256,ISQRT,50000,0.000759312,65849078.528,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,rsa256(composite),256,ISQRT,50000,0.001030383,48525642.987,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,rsa256(composite),256,ISQRT,50000,0.000717252,69710503.562,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,rsa256(composite),256,ISQRT,50000,0.000996753,50162879.303,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,rsa256(composite),256,ISQRT,50000,0.000773751,64620266.845,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,rsa256(composite),256,ISQRT,50000,0.001046083,47797352.107,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,rsa256(composite),256,ISQRT,50000,0.000710421,70380801.479,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,rsa256(composite),256,ISQRT,50000,0.000973393,51366715.407,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,rsa256(composite),256,ISQRT,50000,0.018394172,2718252.287,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,rsa256(composite),256,ISQRT,50000,0.019875223,2515695.038,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,rsa256(composite),256,ISQRT,50000,0.014185157,3524811.176,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,rsa256(composite),256,ISQRT,50000,0.017360059,2880174.549,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,rsa256(composite),256,ISQRT,50000,0.009573555,5222720.249,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,rsa256(composite),256,ISQRT,50000,0.012227137,4089264.726,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MODMUL_R2,50000,0.002859728,17484179.937,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MODMUL_R2,50000,0.004057782,12322002.502,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MODMUL_R2,50000,0.003441540,14528379.489,0
library,NVIDIA H200 NVL,gpu,cgbn,rsa256(composite),256,MODMUL_R2,50000,0.000016000,3125000000.000,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,rsa256(composite),256,MODMUL_R2,50000,0.000099090,504591559.217,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,rsa256(composite),256,MODMUL_R2,50000,0.000373421,133897128.064,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,rsa256(composite),256,MODMUL_R2,50000,0.000029320,1705319802.904,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,rsa256(composite),256,MODMUL_R2,50000,0.000300111,166605078.410,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,rsa256(composite),256,MODMUL_R2,50000,0.000019600,2551016135.374,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,rsa256(composite),256,MODMUL_R2,50000,0.000303771,164597675.227,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,rsa256(composite),256,MODMUL_R2,50000,0.000014290,3498926522.798,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,rsa256(composite),256,MODMUL_R2,50000,0.000291151,171732197.133,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,rsa256(composite),256,MODMUL_R2,50000,0.000018790,2660987761.222,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,rsa256(composite),256,MODMUL_R2,50000,0.000316391,158032307.863,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,rsa256(composite),256,MODMUL_R2,50000,0.000013380,3736931339.128,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,rsa256(composite),256,MODMUL_R2,50000,0.000284061,176018488.680,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,rsa256(composite),256,MODMUL_R2,50000,0.002455002,20366581.803,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,rsa256(composite),256,MODMUL_R2,50000,0.004699180,10640154.196,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,rsa256(composite),256,MODMUL_R2,50000,0.001736142,28799486.098,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,rsa256(composite),256,MODMUL_R2,50000,0.004266909,11718084.520,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,rsa256(composite),256,MODMUL_R2,50000,0.002170990,23030967.562,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,rsa256(composite),256,MODMUL_R2,50000,0.004446531,11244720.723,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,ADD,25000,0.000297286,84094143.200,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,ADD,25000,0.000492165,50795946.263,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,ADD,25000,0.000469378,53262018.680,0
library,NVIDIA H200 NVL,gpu,cgbn,brainpoolP512r1,512,ADD,50000,0.000008448,5918560606.061,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,brainpoolP512r1,512,ADD,50000,0.000059291,843297996.882,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,brainpoolP512r1,512,ADD,50000,0.000487281,102610220.725,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,brainpoolP512r1,512,ADD,50000,0.000032840,1522530538.048,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,brainpoolP512r1,512,ADD,50000,0.000448941,111373204.362,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32,brainpoolP512r1,512,ADD,50000,0.000019650,2544532499.171,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32,brainpoolP512r1,512,ADD,50000,0.000442101,113096318.930,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,brainpoolP512r1,512,ADD,50000,0.000019100,2617797069.508,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,brainpoolP512r1,512,ADD,50000,0.000498392,100322653.565,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,brainpoolP512r1,512,ADD,50000,0.000019930,2508786541.820,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,brainpoolP512r1,512,ADD,50000,0.000495331,100942599.441,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,brainpoolP512r1,512,ADD,50000,0.000011230,4452358156.844,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,brainpoolP512r1,512,ADD,50000,0.000473381,105623161.627,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,brainpoolP512r1,512,ADD,50000,0.000015640,3196944677.177,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,brainpoolP512r1,512,ADD,50000,0.000422531,118334501.424,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,brainpoolP512r1,512,ADD,50000,0.001899910,26317035.475,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,brainpoolP512r1,512,ADD,50000,0.004906275,10191030.956,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,brainpoolP512r1,512,ADD,50000,0.002290736,21827045.323,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,brainpoolP512r1,512,ADD,50000,0.005101865,9800337.641,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,brainpoolP512r1,512,ADD,50000,0.002249431,22227843.580,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,brainpoolP512r1,512,ADD,50000,0.005376645,9299479.684,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,brainpoolP512r1,512,ADD,50000,0.001947259,25677117.598,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,brainpoolP512r1,512,ADD,50000,0.005113624,9777801.571,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACT,25000,0.000190359,131331071.058,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACT,25000,0.000334517,74734642.079,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACT,25000,0.000396494,63052623.695,0
library,NVIDIA H200 NVL,gpu,cgbn,brainpoolP512r1,512,SUBTRACT,50000,0.000008512,5874060150.376,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,brainpoolP512r1,512,SUBTRACT,50000,0.000059830,835700806.326,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,brainpoolP512r1,512,SUBTRACT,50000,0.000477691,104670185.171,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,brainpoolP512r1,512,SUBTRACT,50000,0.000032970,1516530947.354,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,brainpoolP512r1,512,SUBTRACT,50000,0.000418671,119425499.081,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32,brainpoolP512r1,512,SUBTRACT,50000,0.000019760,2530351068.117,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32,brainpoolP512r1,512,SUBTRACT,50000,0.000431901,115767281.357,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,brainpoolP512r1,512,SUBTRACT,50000,0.000019040,2626040057.963,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,brainpoolP512r1,512,SUBTRACT,50000,0.000506981,98623016.927,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,brainpoolP512r1,512,SUBTRACT,50000,0.000019930,2508786541.820,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,brainpoolP512r1,512,SUBTRACT,50000,0.000473081,105690142.098,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,brainpoolP512r1,512,SUBTRACT,50000,0.000011450,4366801175.334,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,brainpoolP512r1,512,SUBTRACT,50000,0.000477122,104794994.201,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,brainpoolP512r1,512,SUBTRACT,50000,0.000015550,3215421635.947,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,brainpoolP512r1,512,SUBTRACT,50000,0.000442361,113029857.373,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,brainpoolP512r1,512,SUBTRACT,50000,0.002010972,24863598.222,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,brainpoolP512r1,512,SUBTRACT,50000,0.005157040,9695484.294,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,brainpoolP512r1,512,SUBTRACT,50000,0.001335879,37428540.320,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,brainpoolP512r1,512,SUBTRACT,50000,0.004657721,10734863.611,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,brainpoolP512r1,512,SUBTRACT,50000,0.002089000,23934897.564,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,brainpoolP512r1,512,SUBTRACT,50000,0.005262074,9501956.934,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,brainpoolP512r1,512,SUBTRACT,50000,0.002053319,24350819.203,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,brainpoolP512r1,512,SUBTRACT,50000,0.005346355,9352166.052,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,ADDMOD,25000,0.000652418,38318962.958,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,ADDMOD,25000,0.000228411,109451982.349,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,ADDMOD,25000,0.000686635,36409429.020,0
library,NVIDIA H200 NVL,gpu,cgbn,brainpoolP512r1,512,ADDMOD,50000,0.000010528,4749240121.581,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,brainpoolP512r1,512,ADDMOD,50000,0.000077640,643997477.370,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,brainpoolP512r1,512,ADDMOD,50000,0.000495152,100979100.247,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,brainpoolP512r1,512,ADDMOD,50000,0.000040471,1235449855.599,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,brainpoolP512r1,512,ADDMOD,50000,0.000424812,117699132.856,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32,brainpoolP512r1,512,ADDMOD,50000,0.000022010,2271700215.800,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32,brainpoolP512r1,512,ADDMOD,50000,0.000444891,112387096.450,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,brainpoolP512r1,512,ADDMOD,50000,0.000023920,2090293225.354,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,brainpoolP512r1,512,ADDMOD,50000,0.000498361,100328887.278,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,brainpoolP512r1,512,ADDMOD,50000,0.000024250,2061854817.432,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,brainpoolP512r1,512,ADDMOD,50000,0.000515761,96944132.804,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,brainpoolP512r1,512,ADDMOD,50000,0.000011851,4219024848.723,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,brainpoolP512r1,512,ADDMOD,50000,0.000482161,103699812.856,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,brainpoolP512r1,512,ADDMOD,50000,0.000012000,4166634939.853,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,brainpoolP512r1,512,ADDMOD,50000,0.000438551,114011807.307,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,brainpoolP512r1,512,ADDMOD,50000,0.002022073,24727099.098,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,brainpoolP512r1,512,ADDMOD,50000,0.005058507,9884339.451,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,brainpoolP512r1,512,ADDMOD,50000,0.002280906,21921113.076,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,brainpoolP512r1,512,ADDMOD,50000,0.005121185,9763365.356,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,brainpoolP512r1,512,ADDMOD,50000,0.002030030,24630179.178,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,brainpoolP512r1,512,ADDMOD,50000,0.005132284,9742251.179,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,brainpoolP512r1,512,ADDMOD,50000,0.002251961,22202870.835,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,brainpoolP512r1,512,ADDMOD,50000,0.005369055,9312625.833,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000610865,40925561.532,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000548765,45556842.460,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000587765,42534006.545,0
library,NVIDIA H200 NVL,gpu,cgbn,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000010816,4622781065.089,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000086471,578227942.826,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000504192,99168577.614,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000044261,1129662964.921,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000449491,111236940.787,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000024050,2079001343.737,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000484891,103115957.382,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000023560,2122239607.864,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000505582,98895932.963,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000023870,2094687060.637,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000496641,100676336.732,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000011740,4258978924.091,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000481671,103805304.141,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000011750,4255349987.615,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000416171,120142910.899,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,brainpoolP512r1,512,SUBTRACTMOD,50000,0.002074594,24101100.910,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,brainpoolP512r1,512,SUBTRACTMOD,50000,0.004955456,10089888.781,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,brainpoolP512r1,512,SUBTRACTMOD,50000,0.002343766,21333187.316,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,brainpoolP512r1,512,SUBTRACTMOD,50000,0.005118064,9769319.023,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,brainpoolP512r1,512,SUBTRACTMOD,50000,0.002075429,24091405.204,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,brainpoolP512r1,512,SUBTRACTMOD,50000,0.004566592,10949084.295,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,brainpoolP512r1,512,SUBTRACTMOD,50000,0.002092730,23892236.216,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,brainpoolP512r1,512,SUBTRACTMOD,50000,0.005053764,9893615.973,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.000693432,36052560.887,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.001267874,19718047.809,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.001522975,16415243.634,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.005285896,9459134.249,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.005805357,8612734.740,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.001405595,35572124.593,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.001912587,26142601.986,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000354130,141191066.532,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000865393,57777212.279,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000059900,834723051.747,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000601481,83128152.575,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000058380,856458342.506,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000640292,78089377.400,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000030630,1632384666.489,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000609801,81993957.333,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000029940,1670004353.320,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000548861,91097747.800,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.003011974,16600409.320,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.006358465,7863533.023,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.002273985,21987831.273,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.004995534,10008940.058,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.002175240,22985969.760,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.005426366,9214269.737,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.002195530,22773544.543,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.005417086,9230054.766,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000689090,36279713.024,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000662427,37740008.255,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000670844,37266505.957,0
library,NVIDIA H200 NVL,gpu,cgbn,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000010560,4734848484.848,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.001261084,39648429.721,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.001788066,27963173.594,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000347151,144029561.865,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000860903,58078556.017,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000097970,510360248.872,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000696582,71779056.132,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000098520,507511378.740,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000694092,72036558.555,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000097921,510615691.417,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000697622,71672055.812,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000046420,1077121292.054,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000630612,79288054.236,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000060170,830977871.249,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000565061,88486022.419,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.001627864,30715095.022,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.005545238,9016745.486,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.002077744,24064562.412,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.005495827,9097811.745,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.002275681,21971444.966,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.005403365,9253492.948,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.002217271,22550244.091,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.005630466,8880259.583,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.006020138,4152728.697,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.004810555,5196905.530,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.003680221,6793070.237,0
library,NVIDIA H200 NVL,gpu,cgbn,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000014560,3434065934.066,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000288871,173087640.884,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000752502,66445015.998,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000085730,583226481.807,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000481352,103874092.903,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000028451,1757409119.775,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000508731,98283769.164,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000039090,1279097773.847,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000515131,97062701.989,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000027300,1831503505.682,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000514111,97255258.979,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000032530,1537040377.052,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000522351,95721074.821,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000020600,2427193417.424,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000436561,114531544.974,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.006334085,7893799.951,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.008575862,5830317.671,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.002610298,19154901.161,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.005522558,9053775.443,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.002050360,24385960.507,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.005001963,9996075.541,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.002025189,24689053.779,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.005364546,9320453.264,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,COMPARE,25000,0.000107699,232127541.012,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,COMPARE,25000,0.000339265,73688673.410,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,COMPARE,25000,0.000401035,62338647.389,0
library,NVIDIA H200 NVL,gpu,cgbn,brainpoolP512r1,512,COMPARE,50000,0.000008352,5986590038.314,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,brainpoolP512r1,512,COMPARE,50000,0.000056290,888256170.480,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,brainpoolP512r1,512,COMPARE,50000,0.000504232,99160701.495,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,brainpoolP512r1,512,COMPARE,50000,0.000030200,1655629296.574,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,brainpoolP512r1,512,COMPARE,50000,0.000414672,120577205.541,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,brainpoolP512r1,512,COMPARE,50000,0.000018920,2642715277.410,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,brainpoolP512r1,512,COMPARE,50000,0.000489962,102048718.976,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,brainpoolP512r1,512,COMPARE,50000,0.000018860,2651116190.758,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,brainpoolP512r1,512,COMPARE,50000,0.000501871,99627196.506,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,brainpoolP512r1,512,COMPARE,50000,0.000011050,4524876258.705,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,brainpoolP512r1,512,COMPARE,50000,0.000485661,102952476.487,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,brainpoolP512r1,512,COMPARE,50000,0.000010900,4587121034.700,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,brainpoolP512r1,512,COMPARE,50000,0.000415811,120246915.300,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,brainpoolP512r1,512,COMPARE,50000,0.002117835,23609015.975,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,brainpoolP512r1,512,COMPARE,50000,0.005056918,9887445.245,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,brainpoolP512r1,512,COMPARE,50000,0.002078194,24059350.899,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,brainpoolP512r1,512,COMPARE,50000,0.005057454,9886397.399,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,brainpoolP512r1,512,COMPARE,50000,0.002017089,24788196.642,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,brainpoolP512r1,512,COMPARE,50000,0.005307696,9420283.193,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,REDUCE,3125,0.000051670,60479789.028,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,REDUCE,3125,0.000058146,53744469.908,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,REDUCE,3125,0.000063338,49338376.106,0
library,NVIDIA H200 NVL,gpu,cgbn,brainpoolP512r1,512,REDUCE,50000,0.000019360,2582644628.099,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,brainpoolP512r1,512,REDUCE,50000,0.000452781,110428666.935,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,brainpoolP512r1,512,REDUCE,50000,0.000874803,57155721.930,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,brainpoolP512r1,512,REDUCE,50000,0.000356721,140165566.198,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,brainpoolP512r1,512,REDUCE,50000,0.000763002,65530630.440,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,brainpoolP512r1,512,REDUCE,50000,0.000101160,494266954.061,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,brainpoolP512r1,512,REDUCE,50000,0.000555442,90018389.377,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,brainpoolP512r1,512,REDUCE,50000,0.000100940,495343767.271,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,brainpoolP512r1,512,REDUCE,50000,0.000588331,84986169.528,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,brainpoolP512r1,512,REDUCE,50000,0.000103000,485437037.489,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,brainpoolP512r1,512,REDUCE,50000,0.000587962,85039511.420,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,brainpoolP512r1,512,REDUCE,50000,0.000120911,413527536.715,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,brainpoolP512r1,512,REDUCE,50000,0.000526682,94933943.509,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,brainpoolP512r1,512,REDUCE,50000,0.004034755,12392326.013,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,brainpoolP512r1,512,REDUCE,50000,0.006699752,7462962.738,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,brainpoolP512r1,512,REDUCE,50000,0.004052838,12337034.068,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,brainpoolP512r1,512,REDUCE,50000,0.007596532,6581950.817,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,brainpoolP512r1,512,REDUCE,50000,0.003696728,13525474.302,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,brainpoolP512r1,512,REDUCE,50000,0.006173009,8099777.527,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL,1562,0.000181114,8624378.549,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL,1562,0.000187163,8345664.162,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL,1562,0.000080431,19420432.591,0
library,NVIDIA H200 NVL,gpu,cgbn,brainpoolP512r1,512,MODMUL,50000,0.000095104,525740242.261,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,brainpoolP512r1,512,MODMUL,50000,0.001272504,39292607.792,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,brainpoolP512r1,512,MODMUL,50000,0.001692335,29544977.711,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,brainpoolP512r1,512,MODMUL,50000,0.000920463,54320488.466,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,brainpoolP512r1,512,MODMUL,50000,0.001321065,37848251.027,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,brainpoolP512r1,512,MODMUL,50000,0.000364701,137098641.740,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,brainpoolP512r1,512,MODMUL,50000,0.000849722,58842777.696,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,brainpoolP512r1,512,MODMUL,50000,0.000294591,169726869.885,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,brainpoolP512r1,512,MODMUL,50000,0.000747082,66927052.479,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,brainpoolP512r1,512,MODMUL,50000,0.000358661,139407431.008,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,brainpoolP512r1,512,MODMUL,50000,0.000825752,60550864.749,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,brainpoolP512r1,512,MODMUL,50000,0.000362731,137843175.828,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,brainpoolP512r1,512,MODMUL,50000,0.000769362,64988905.591,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,brainpoolP512r1,512,MODMUL,50000,0.016387489,3051108.088,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,brainpoolP512r1,512,MODMUL,50000,0.016862519,2965156.042,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,brainpoolP512r1,512,MODMUL,50000,0.009964728,5017698.411,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,brainpoolP512r1,512,MODMUL,50000,0.012892638,3878182.232,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,brainpoolP512r1,512,MODMUL,50000,0.006877082,7270525.440,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,brainpoolP512r1,512,MODMUL,50000,0.010711480,4667889.010,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MODEXP,390,0.011576074,33690.179,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MODEXP,390,0.003608020,108092.530,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MODEXP,390,0.004452383,87593.542,0
library,NVIDIA H200 NVL,gpu,cgbn,brainpoolP512r1,512,MODEXP,50000,0.024278849,2059405.699,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,brainpoolP512r1,512,MODEXP,50000,0.249073699,200743.797,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,brainpoolP512r1,512,MODEXP,50000,0.249679961,200256.359,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,brainpoolP512r1,512,MODEXP,50000,0.024779427,2017802.915,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,brainpoolP512r1,512,MODEXP,50000,0.025181180,1985609.889,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,brainpoolP512r1,512,MODEXP,50000,0.021583965,2316534.515,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,brainpoolP512r1,512,MODEXP,50000,0.022078977,2264597.667,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,brainpoolP512r1,512,MODEXP,50000,0.009630454,5191863.226,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,brainpoolP512r1,512,MODEXP,50000,0.010172956,4914992.296,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,brainpoolP512r1,512,MODEXP,50000,0.022057756,2266776.370,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,brainpoolP512r1,512,MODEXP,50000,0.022505038,2221724.747,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,brainpoolP512r1,512,MODEXP,50000,0.009539804,5241197.788,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,brainpoolP512r1,512,MODEXP,50000,0.009971206,5014438.595,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,brainpoolP512r1,512,MODEXP,50000,6.823305665,7327.826,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,brainpoolP512r1,512,MODEXP,50000,6.833213278,7317.202,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,brainpoolP512r1,512,MODEXP,50000,1.377823224,36289.126,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,brainpoolP512r1,512,MODEXP,50000,1.342330011,37248.664,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,brainpoolP512r1,512,MODEXP,50000,0.747821899,66860.840,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,brainpoolP512r1,512,MODEXP,50000,0.774140591,64587.751,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,EXPONENTIATION,390,0.002612278,149294.984,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,EXPONENTIATION,390,0.001847835,211057.809,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,EXPONENTIATION,390,0.002475807,157524.397,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,brainpoolP512r1,512,EXPONENTIATION,50000,0.207710698,240719.426,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,brainpoolP512r1,512,EXPONENTIATION,50000,0.207973059,240415.755,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,brainpoolP512r1,512,EXPONENTIATION,50000,0.051885703,963656.597,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,brainpoolP512r1,512,EXPONENTIATION,50000,0.052180735,958208.045,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,brainpoolP512r1,512,EXPONENTIATION,50000,0.016382332,3052068.541,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,brainpoolP512r1,512,EXPONENTIATION,50000,0.016778233,2980051.595,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,brainpoolP512r1,512,EXPONENTIATION,50000,0.014617137,3420642.490,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,brainpoolP512r1,512,EXPONENTIATION,50000,0.015082889,3315014.782,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,brainpoolP512r1,512,EXPONENTIATION,50000,0.015678990,3188980.913,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,brainpoolP512r1,512,EXPONENTIATION,50000,0.016201441,3086145.249,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,brainpoolP512r1,512,EXPONENTIATION,50000,0.013987765,3574552.458,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,brainpoolP512r1,512,EXPONENTIATION,50000,0.014418857,3467681.232,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,brainpoolP512r1,512,EXPONENTIATION,50000,3.238132208,15441.000,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,brainpoolP512r1,512,EXPONENTIATION,50000,3.173989835,15753.044,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,brainpoolP512r1,512,EXPONENTIATION,50000,0.626164928,79851.167,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,brainpoolP512r1,512,EXPONENTIATION,50000,0.647116352,77265.858,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,brainpoolP512r1,512,EXPONENTIATION,50000,0.379513545,131747.603,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,brainpoolP512r1,512,EXPONENTIATION,50000,0.310690512,160931.854,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,DIVIDE,3125,0.000069938,44682528.710,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,DIVIDE,3125,0.000176398,17715664.922,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,DIVIDE,3125,0.000107757,29000492.177,0
library,NVIDIA H200 NVL,gpu,cgbn,brainpoolP512r1,512,DIVIDE,50000,0.000027680,1806358381.503,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,brainpoolP512r1,512,DIVIDE,50000,0.001262403,39607004.195,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,brainpoolP512r1,512,DIVIDE,50000,0.001762126,28374815.577,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,brainpoolP512r1,512,DIVIDE,50000,0.001242774,40232577.882,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,brainpoolP512r1,512,DIVIDE,50000,0.001744396,28663217.418,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,brainpoolP512r1,512,DIVIDE,50000,0.000413811,120828121.912,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,brainpoolP512r1,512,DIVIDE,50000,0.000968882,51605871.555,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,brainpoolP512r1,512,DIVIDE,50000,0.000386771,129275446.869,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,brainpoolP512r1,512,DIVIDE,50000,0.000943692,52983385.739,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,brainpoolP512r1,512,DIVIDE,50000,0.000395801,126326097.229,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,brainpoolP512r1,512,DIVIDE,50000,0.000950372,52610979.646,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,brainpoolP512r1,512,DIVIDE,50000,0.000353720,141354728.228,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,brainpoolP512r1,512,DIVIDE,50000,0.000846312,59079872.414,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,brainpoolP512r1,512,DIVIDE,50000,0.007221616,6923658.094,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,brainpoolP512r1,512,DIVIDE,50000,0.009530672,5246219.790,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,brainpoolP512r1,512,DIVIDE,50000,0.005914600,8453657.075,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,brainpoolP512r1,512,DIVIDE,50000,0.010290540,4858831.483,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,brainpoolP512r1,512,DIVIDE,50000,0.004786922,10445125.250,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,brainpoolP512r1,512,DIVIDE,50000,0.008716751,5736082.167,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,ISQRT,781,0.000060964,12810817.018,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,ISQRT,781,0.000084518,9240623.572,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,brainpoolP512r1,512,ISQRT,50000,0.019057135,2623689.224,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,brainpoolP512r1,512,ISQRT,50000,0.019491957,2565160.590,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,brainpoolP512r1,512,ISQRT,50000,0.018954607,2637881.112,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,brainpoolP512r1,512,ISQRT,50000,0.019402059,2577046.087,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,brainpoolP512r1,512,ISQRT,50000,0.003715099,13458591.722,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,brainpoolP512r1,512,ISQRT,50000,0.004201621,11900169.194,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,brainpoolP512r1,512,ISQRT,50000,0.003372259,14826856.540,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,brainpoolP512r1,512,ISQRT,50000,0.003865310,12935573.115,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,brainpoolP512r1,512,ISQRT,50000,0.003528699,14169528.018,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,brainpoolP512r1,512,ISQRT,50000,0.003992690,12522885.554,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,brainpoolP512r1,512,ISQRT,50000,0.004262261,11730863.098,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,brainpoolP512r1,512,ISQRT,50000,0.004664192,10719970.149,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,brainpoolP512r1,512,ISQRT,50000,0.290033770,172393.718,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,brainpoolP512r1,512,ISQRT,50000,0.227770981,219518.745,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,brainpoolP512r1,512,ISQRT,50000,0.124638449,401160.320,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,brainpoolP512r1,512,ISQRT,50000,0.184348447,271225.502,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,brainpoolP512r1,512,ISQRT,50000,0.105837338,472423.069,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,brainpoolP512r1,512,ISQRT,50000,0.111572485,448139.162,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL_R2,25000,0.002952759,8466657.685,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL_R2,25000,0.003388760,7377329.697,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL_R2,25000,0.003416880,7316616.256,0
library,NVIDIA H200 NVL,gpu,cgbn,brainpoolP512r1,512,MODMUL_R2,50000,0.000023040,2170138888.889,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,brainpoolP512r1,512,MODMUL_R2,50000,0.000438381,114056041.537,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,brainpoolP512r1,512,MODMUL_R2,50000,0.000875253,57126339.633,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,brainpoolP512r1,512,MODMUL_R2,50000,0.000085110,587474496.437,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,brainpoolP512r1,512,MODMUL_R2,50000,0.000498052,100391126.245,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,brainpoolP512r1,512,MODMUL_R2,50000,0.000059990,833472530.322,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,brainpoolP512r1,512,MODMUL_R2,50000,0.000555091,90075309.320,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,brainpoolP512r1,512,MODMUL_R2,50000,0.000035330,1415229666.339,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,brainpoolP512r1,512,MODMUL_R2,50000,0.000531292,94110217.730,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,brainpoolP512r1,512,MODMUL_R2,50000,0.000054790,912574349.560,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,brainpoolP512r1,512,MODMUL_R2,50000,0.000508622,98304847.507,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,brainpoolP512r1,512,MODMUL_R2,50000,0.000031770,1573813052.304,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,brainpoolP512r1,512,MODMUL_R2,50000,0.000451001,110864501.794,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,brainpoolP512r1,512,MODMUL_R2,50000,0.006227434,8028989.258,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,brainpoolP512r1,512,MODMUL_R2,50000,0.010019912,4990063.767,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,brainpoolP512r1,512,MODMUL_R2,50000,0.002607928,19172308.387,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,brainpoolP512r1,512,MODMUL_R2,50000,0.005682278,8799287.955,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,brainpoolP512r1,512,MODMUL_R2,50000,0.001914629,26114721.269,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,brainpoolP512r1,512,MODMUL_R2,50000,0.004907083,10189352.821,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,p1024,1024,ADD,12500,0.000187629,66620995.673,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,p1024,1024,ADD,12500,0.000213573,58528087.864,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,p1024,1024,ADD,12500,0.000303752,41152046.231,0
library,NVIDIA H200 NVL,gpu,cgbn,p1024,1024,ADD,50000,0.000011008,4542151162.791,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,p1024,1024,ADD,50000,0.000209051,239176141.678,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,p1024,1024,ADD,50000,0.000877582,56974732.688,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,p1024,1024,ADD,50000,0.000107270,466113182.910,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,p1024,1024,ADD,50000,0.000792483,63092837.315,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32,p1024,1024,ADD,50000,0.000044050,1135072636.216,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32,p1024,1024,ADD,50000,0.000934693,53493504.892,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,p1024,1024,ADD,50000,0.000043800,1141552921.414,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,p1024,1024,ADD,50000,0.000933932,53537086.974,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,p1024,1024,ADD,50000,0.000043760,1142594572.967,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,p1024,1024,ADD,50000,0.000899283,55599851.180,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,p1024,1024,ADD,50000,0.000022970,2176750052.455,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,p1024,1024,ADD,50000,0.000878452,56918307.969,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,p1024,1024,ADD,50000,0.000022740,2198769957.253,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,p1024,1024,ADD,50000,0.000704192,71003365.304,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p1024,1024,ADD,50000,0.002042425,24480703.308,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p1024,1024,ADD,50000,0.007079676,7062470.180,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p1024,1024,ADD,50000,0.001996383,25045294.985,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p1024,1024,ADD,50000,0.007011238,7131408.156,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,p1024,1024,ADD,50000,0.002107480,23725018.020,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,p1024,1024,ADD,50000,0.006993953,7149032.982,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,p1024,1024,ADD,50000,0.002238500,22336385.490,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,p1024,1024,ADD,50000,0.006203570,8059875.131,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,p1024,1024,SUBTRACT,12500,0.000130361,95887372.094,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,p1024,1024,SUBTRACT,12500,0.000369475,33831813.464,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,p1024,1024,SUBTRACT,12500,0.000612189,20418530.931,0
library,NVIDIA H200 NVL,gpu,cgbn,p1024,1024,SUBTRACT,50000,0.000010752,4650297619.048,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,p1024,1024,SUBTRACT,50000,0.000208830,239429206.237,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,p1024,1024,SUBTRACT,50000,0.000910633,54906860.840,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,p1024,1024,SUBTRACT,50000,0.000106851,467941390.412,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,p1024,1024,SUBTRACT,50000,0.000768562,65056563.227,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32,p1024,1024,SUBTRACT,50000,0.000044000,1136364003.313,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32,p1024,1024,SUBTRACT,50000,0.000929882,53770267.082,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,p1024,1024,SUBTRACT,50000,0.000043980,1136881373.268,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,p1024,1024,SUBTRACT,50000,0.000917422,54500545.911,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,p1024,1024,SUBTRACT,50000,0.000044070,1134559907.861,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,p1024,1024,SUBTRACT,50000,0.000915172,54634536.006,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,p1024,1024,SUBTRACT,50000,0.000022700,2202648992.005,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,p1024,1024,SUBTRACT,50000,0.000916082,54580270.029,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,p1024,1024,SUBTRACT,50000,0.000022470,2225187184.482,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,p1024,1024,SUBTRACT,50000,0.000710972,70326257.306,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p1024,1024,SUBTRACT,50000,0.002023865,24705205.101,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p1024,1024,SUBTRACT,50000,0.006951656,7192530.737,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p1024,1024,SUBTRACT,50000,0.001976933,25291702.086,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p1024,1024,SUBTRACT,50000,0.006890047,7256844.555,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,p1024,1024,SUBTRACT,50000,0.002047910,24415135.362,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,p1024,1024,SUBTRACT,50000,0.007451885,6709711.676,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,p1024,1024,SUBTRACT,50000,0.002017430,24784008.440,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,p1024,1024,SUBTRACT,50000,0.006668642,7497778.486,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,p1024,1024,ADDMOD,12500,0.000476964,26207435.451,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,p1024,1024,ADDMOD,12500,0.000468560,26677472.498,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,p1024,1024,ADDMOD,12500,0.000713806,17511766.900,0
library,NVIDIA H200 NVL,gpu,cgbn,p1024,1024,ADDMOD,50000,0.000012256,4079634464.752,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,p1024,1024,ADDMOD,50000,0.000283900,176118310.416,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,p1024,1024,ADDMOD,50000,0.000969643,51565369.611,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,p1024,1024,ADDMOD,50000,0.000143911,347436979.225,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,p1024,1024,ADDMOD,50000,0.000802333,62318264.280,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32,p1024,1024,ADDMOD,50000,0.000054660,914745359.362,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32,p1024,1024,ADDMOD,50000,0.000940752,53148971.154,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,p1024,1024,ADDMOD,50000,0.000055861,895079483.663,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,p1024,1024,ADDMOD,50000,0.000934722,53491845.962,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,p1024,1024,ADDMOD,50000,0.000054490,917598649.763,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,p1024,1024,ADDMOD,50000,0.000898042,55676683.568,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,p1024,1024,ADDMOD,50000,0.000016150,3095962816.447,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,p1024,1024,ADDMOD,50000,0.000916662,54545729.621,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,p1024,1024,ADDMOD,50000,0.000015810,3162575509.182,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,p1024,1024,ADDMOD,50000,0.000668872,74752727.512,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p1024,1024,ADDMOD,50000,0.001979475,25259222.442,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p1024,1024,ADDMOD,50000,0.006866835,7281374.837,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p1024,1024,ADDMOD,50000,0.001966423,25426879.086,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p1024,1024,ADDMOD,50000,0.006659206,7508402.650,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,p1024,1024,ADDMOD,50000,0.002272501,22002190.482,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,p1024,1024,ADDMOD,50000,0.006819432,7331988.919,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,p1024,1024,ADDMOD,50000,0.002151650,23237979.859,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,p1024,1024,ADDMOD,50000,0.006612961,7560909.613,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,p1024,1024,SUBTRACTMOD,12500,0.000388245,32196192.040,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,p1024,1024,SUBTRACTMOD,12500,0.001276964,9788845.502,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,p1024,1024,SUBTRACTMOD,12500,0.001468831,8510168.892,0
library,NVIDIA H200 NVL,gpu,cgbn,p1024,1024,SUBTRACTMOD,50000,0.000012768,3916040100.251,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,p1024,1024,SUBTRACTMOD,50000,0.000285501,175130729.092,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,p1024,1024,SUBTRACTMOD,50000,0.000970133,51539325.116,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,p1024,1024,SUBTRACTMOD,50000,0.000143231,349086418.525,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,p1024,1024,SUBTRACTMOD,50000,0.000822013,60826287.806,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32,p1024,1024,SUBTRACTMOD,50000,0.000055080,907771661.668,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32,p1024,1024,SUBTRACTMOD,50000,0.001001073,49946405.134,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,p1024,1024,SUBTRACTMOD,50000,0.000056150,890471631.369,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,p1024,1024,SUBTRACTMOD,50000,0.000942602,53044660.684,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,p1024,1024,SUBTRACTMOD,50000,0.000055760,896699680.777,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,p1024,1024,SUBTRACTMOD,50000,0.000914972,54646478.437,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,p1024,1024,SUBTRACTMOD,50000,0.000016400,3048778914.641,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,p1024,1024,SUBTRACTMOD,50000,0.000702152,71209649.702,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,p1024,1024,SUBTRACTMOD,50000,0.000015890,3146634501.150,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,p1024,1024,SUBTRACTMOD,50000,0.000693612,72086420.040,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p1024,1024,SUBTRACTMOD,50000,0.001995464,25056829.240,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p1024,1024,SUBTRACTMOD,50000,0.006744445,7413508.375,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p1024,1024,SUBTRACTMOD,50000,0.001991524,25106400.778,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p1024,1024,SUBTRACTMOD,50000,0.006992847,7150163.629,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,p1024,1024,SUBTRACTMOD,50000,0.002085400,23976215.056,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,p1024,1024,SUBTRACTMOD,50000,0.008079048,6188847.983,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,p1024,1024,SUBTRACTMOD,50000,0.002056470,24313507.810,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,p1024,1024,SUBTRACTMOD,50000,0.006640931,7529064.815,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.001225684,10198387.080,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.002021836,6182499.469,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.002008271,6224259.586,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.025213584,1983058.025,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.026212336,1907498.827,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.006482683,7712855.920,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.007428036,6731254.334,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.001406323,35553708.702,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.002472647,20221244.681,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.000142300,351370255.934,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.001225773,40790586.430,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.000143440,348577864.204,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.001212763,41228170.103,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.000090940,549812818.481,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.000968052,51650120.215,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.000082280,607680932.681,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.000969782,51557978.691,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.007568237,6606558.417,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.012806429,3904288.991,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.002611847,19143540.860,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.007930714,6304602.543,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.002000140,24998249.795,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.007404605,6752554.727,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.001685978,29656379.008,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.006818122,7333397.674,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.001227318,10184805.464,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.001995996,6262537.632,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.001752650,7132057.342,0
library,NVIDIA H200 NVL,gpu,cgbn,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000024736,2021345407.503,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.009386027,5327067.566,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.010269830,4868629.761,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.002410469,20742850.332,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.003267502,15302209.161,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000617952,80912422.064,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.001684444,29683384.697,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000619382,80725624.501,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.001710724,29227391.971,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000618281,80869383.484,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.001690005,29585712.121,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000216031,231448284.306,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.001087923,45959136.997,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000201360,248311514.862,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.001086013,46039967.670,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.005639723,8865683.594,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.010694284,4675394.822,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.002054684,24334642.443,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.007089728,7052456.714,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.002115360,23636638.731,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.007008833,7133855.241,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.002038680,24525673.271,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.006920853,7224542.937,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.009387907,1331500.196,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.002927479,4269885.476,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.003592211,3479750.997,0
library,NVIDIA H200 NVL,gpu,cgbn,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000037312,1340051457.976,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.002030186,24628285.219,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.002723678,18357529.679,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000284222,175918852.091,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000961043,52026811.589,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000081660,612294932.526,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000925212,54041672.121,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000137100,364697311.324,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.001037763,48180559.181,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000098910,505510315.785,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000989933,50508467.095,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000100430,497859281.060,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000780822,64035080.100,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000056450,885739241.331,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000734742,68051095.211,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.018177711,2750621.361,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.023407693,2136049.886,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.006270043,7974427.008,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.011549899,4329042.183,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.002718413,18393084.387,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.007208974,6935799.756,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.002326171,21494551.282,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.006997373,7145538.807,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,p1024,1024,COMPARE,12500,0.000061206,204226706.970,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,p1024,1024,COMPARE,12500,0.000749411,16679768.687,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,p1024,1024,COMPARE,12500,0.000446285,28009000.094,0
library,NVIDIA H200 NVL,gpu,cgbn,p1024,1024,COMPARE,50000,0.000011008,4542151162.791,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,p1024,1024,COMPARE,50000,0.000128840,388078297.177,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,p1024,1024,COMPARE,50000,0.000793692,62996723.353,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,p1024,1024,COMPARE,50000,0.000064831,771234666.499,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,p1024,1024,COMPARE,50000,0.000736733,67867195.071,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,p1024,1024,COMPARE,50000,0.000036470,1370989480.808,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,p1024,1024,COMPARE,50000,0.000884493,56529559.999,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,p1024,1024,COMPARE,50000,0.000036310,1377033438.923,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,p1024,1024,COMPARE,50000,0.000885892,56440290.921,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,p1024,1024,COMPARE,50000,0.000014370,3479481270.608,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,p1024,1024,COMPARE,50000,0.000710282,70394574.774,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,p1024,1024,COMPARE,50000,0.000014320,3491644618.599,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,p1024,1024,COMPARE,50000,0.000702491,71175286.029,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p1024,1024,COMPARE,50000,0.002298895,21749579.341,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p1024,1024,COMPARE,50000,0.006404275,7807284.887,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p1024,1024,COMPARE,50000,0.002051034,24377947.773,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p1024,1024,COMPARE,50000,0.005950470,8402697.617,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,p1024,1024,COMPARE,50000,0.002294521,21791040.826,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,p1024,1024,COMPARE,50000,0.006021288,8303871.211,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,p1024,1024,REDUCE,1562,0.000019299,80938442.671,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,p1024,1024,REDUCE,1562,0.000461566,3384133.883,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,p1024,1024,REDUCE,1562,0.000172912,9033499.175,0
library,NVIDIA H200 NVL,gpu,cgbn,p1024,1024,REDUCE,50000,0.000027104,1844746162.928,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,p1024,1024,REDUCE,50000,0.002424328,20624271.913,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,p1024,1024,REDUCE,50000,0.003113329,16059979.210,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,p1024,1024,REDUCE,50000,0.001302744,38380525.171,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,p1024,1024,REDUCE,50000,0.001966997,25419458.563,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,p1024,1024,REDUCE,50000,0.000340041,147041089.516,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,p1024,1024,REDUCE,50000,0.001211753,41262534.573,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,p1024,1024,REDUCE,50000,0.000355791,140531914.763,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,p1024,1024,REDUCE,50000,0.001219793,40990561.849,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,p1024,1024,REDUCE,50000,0.000344461,145154287.035,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,p1024,1024,REDUCE,50000,0.001006042,49699714.193,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,p1024,1024,REDUCE,50000,0.000337331,148222380.680,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,p1024,1024,REDUCE,50000,0.001023883,48833704.631,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p1024,1024,REDUCE,50000,0.011172565,4475248.101,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p1024,1024,REDUCE,50000,0.014755093,3388660.438,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p1024,1024,REDUCE,50000,0.010185640,4908871.671,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p1024,1024,REDUCE,50000,0.015158993,3298372.112,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,p1024,1024,REDUCE,50000,0.008602490,5812270.667,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,p1024,1024,REDUCE,50000,0.012896151,3877125.804,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,p1024,1024,MODMUL,781,0.000248462,3143342.449,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,p1024,1024,MODMUL,781,0.000489846,1594378.955,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,p1024,1024,MODMUL,781,0.000637815,1224493.696,0
library,NVIDIA H200 NVL,gpu,cgbn,p1024,1024,MODMUL,50000,0.000234144,213543802.105,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,p1024,1024,MODMUL,50000,0.009286367,5384236.923,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,p1024,1024,MODMUL,50000,0.009980649,5009694.265,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,p1024,1024,MODMUL,50000,0.003372512,14825744.385,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,p1024,1024,MODMUL,50000,0.004062855,12306616.795,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,p1024,1024,MODMUL,50000,0.001613584,30986921.082,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,p1024,1024,MODMUL,50000,0.002498447,20012431.081,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,p1024,1024,MODMUL,50000,0.001215673,41129483.552,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,p1024,1024,MODMUL,50000,0.002070776,24145538.247,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,p1024,1024,MODMUL,50000,0.001506363,33192529.700,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,p1024,1024,MODMUL,50000,0.002188095,22850927.321,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,p1024,1024,MODMUL,50000,0.001180174,42366636.025,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,p1024,1024,MODMUL,50000,0.001867825,26769104.068,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p1024,1024,MODMUL,50000,0.192630901,259563.755,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p1024,1024,MODMUL,50000,0.115864999,431536.706,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p1024,1024,MODMUL,50000,0.097474644,512953.913,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p1024,1024,MODMUL,50000,0.102424008,488166.797,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,p1024,1024,MODMUL,50000,0.022193815,2252879.915,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,p1024,1024,MODMUL,50000,0.095396009,524130.942,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,p1024,1024,MODEXP,195,0.036424807,5353.494,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,p1024,1024,MODEXP,195,0.006023298,32374.291,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,p1024,1024,MODEXP,195,0.003859362,50526.486,0
library,NVIDIA H200 NVL,gpu,cgbn,p1024,1024,MODEXP,50000,0.125432357,398621.227,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,p1024,1024,MODEXP,50000,2.922903054,17106.281,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,p1024,1024,MODEXP,50000,2.920508327,17120.307,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,p1024,1024,MODEXP,50000,0.173068823,288902.410,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,p1024,1024,MODEXP,50000,0.173856436,287593.610,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,p1024,1024,MODEXP,50000,0.173005121,289008.786,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,p1024,1024,MODEXP,50000,0.173662964,287914.008,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,p1024,1024,MODEXP,50000,0.090451741,552780.957,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,p1024,1024,MODEXP,50000,0.091200453,548242.891,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,p1024,1024,MODEXP,50000,0.173410733,288332.787,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,p1024,1024,MODEXP,50000,0.174067325,287245.179,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,p1024,1024,MODEXP,50000,0.090849922,550358.206,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,p1024,1024,MODEXP,50000,0.091626424,545694.111,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p1024,1024,MODEXP,50000,55.813199363,895.845,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p1024,1024,MODEXP,50000,55.551730201,900.062,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p1024,1024,MODEXP,50000,10.868920803,4600.273,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p1024,1024,MODEXP,50000,10.859768810,4604.150,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,p1024,1024,MODEXP,50000,5.448534024,9176.780,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,p1024,1024,MODEXP,50000,5.447019457,9179.332,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,p1024,1024,EXPONENTIATION,195,0.005367125,36332.301,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,p1024,1024,EXPONENTIATION,195,0.004191792,46519.484,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,p1024,1024,EXPONENTIATION,195,0.006554829,29749.060,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,p1024,1024,EXPONENTIATION,50000,1.763092270,28359.264,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,p1024,1024,EXPONENTIATION,50000,1.762478728,28369.137,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,p1024,1024,EXPONENTIATION,50000,0.406242139,123079.305,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,p1024,1024,EXPONENTIATION,50000,0.407658374,122651.718,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,p1024,1024,EXPONENTIATION,50000,0.124865529,400430.771,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,p1024,1024,EXPONENTIATION,50000,0.125735111,397661.398,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,p1024,1024,EXPONENTIATION,50000,0.108166336,462251.028,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,p1024,1024,EXPONENTIATION,50000,0.108818448,459480.914,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,p1024,1024,EXPONENTIATION,50000,0.123508185,404831.469,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,p1024,1024,EXPONENTIATION,50000,0.124314708,402205.023,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,p1024,1024,EXPONENTIATION,50000,0.108507837,460796.210,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,p1024,1024,EXPONENTIATION,50000,0.109291179,457493.463,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p1024,1024,EXPONENTIATION,50000,22.644258261,2208.065,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p1024,1024,EXPONENTIATION,50000,22.766789344,2196.181,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p1024,1024,EXPONENTIATION,50000,5.423593581,9218.980,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p1024,1024,EXPONENTIATION,50000,5.366645653,9316.807,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,p1024,1024,DIVIDE,1562,0.000039883,39164196.719,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,p1024,1024,DIVIDE,1562,0.000069288,22543649.968,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,p1024,1024,DIVIDE,1562,0.000050351,31022232.637,0
library,NVIDIA H200 NVL,gpu,cgbn,p1024,1024,DIVIDE,50000,0.000037088,1348144952.545,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,p1024,1024,DIVIDE,50000,0.012231326,4087864.239,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,p1024,1024,DIVIDE,50000,0.013662740,3659588.046,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,p1024,1024,DIVIDE,50000,0.007701557,6492193.710,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,p1024,1024,DIVIDE,50000,0.008562430,5839463.769,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,p1024,1024,DIVIDE,50000,0.001666184,30008689.632,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,p1024,1024,DIVIDE,50000,0.002676666,18679955.542,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,p1024,1024,DIVIDE,50000,0.001687044,29637637.588,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,p1024,1024,DIVIDE,50000,0.002517837,19858314.108,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,p1024,1024,DIVIDE,50000,0.001644724,30400236.465,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,p1024,1024,DIVIDE,50000,0.002487576,20099888.530,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,p1024,1024,DIVIDE,50000,0.001697045,29462979.406,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,p1024,1024,DIVIDE,50000,0.002539916,19685690.672,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p1024,1024,DIVIDE,50000,0.088967612,562002.271,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p1024,1024,DIVIDE,50000,0.024121809,2072813.035,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p1024,1024,DIVIDE,50000,0.021139584,2365231.026,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p1024,1024,DIVIDE,50000,0.025441032,1965329.078,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,p1024,1024,DIVIDE,50000,0.015792874,3165984.854,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,p1024,1024,DIVIDE,50000,0.022111774,2261238.733,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,p1024,1024,ISQRT,390,0.000061166,6376056.639,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,p1024,1024,ISQRT,390,0.000122820,3175382.608,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,p1024,1024,ISQRT,50000,0.212656222,235121.265,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,p1024,1024,ISQRT,50000,0.212712713,235058.823,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,p1024,1024,ISQRT,50000,0.135710571,368431.137,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,p1024,1024,ISQRT,50000,0.136551724,366161.616,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,p1024,1024,ISQRT,50000,0.025277264,1978062.187,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,p1024,1024,ISQRT,50000,0.026104847,1915353.117,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,p1024,1024,ISQRT,50000,0.024895933,2008360.161,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,p1024,1024,ISQRT,50000,0.025572396,1955233.291,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,p1024,1024,ISQRT,50000,0.025347004,1972619.718,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,p1024,1024,ISQRT,50000,0.026040277,1920102.456,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,p1024,1024,ISQRT,50000,0.024202891,2065868.904,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,p1024,1024,ISQRT,50000,0.024882464,2009447.289,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p1024,1024,ISQRT,50000,1.188672851,42063.718,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p1024,1024,ISQRT,50000,1.196726237,41780.650,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p1024,1024,ISQRT,50000,0.829063692,60308.997,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p1024,1024,ISQRT,50000,0.830627552,60195.451,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,p1024,1024,ISQRT,50000,0.595749103,83927.948,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,p1024,1024,ISQRT,50000,0.596010624,83891.122,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,p1024,1024,MODMUL_R2,12500,0.004048342,3087683.796,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,p1024,1024,MODMUL_R2,12500,0.003351539,3729629.854,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,p1024,1024,MODMUL_R2,12500,0.006011528,2079338.244,0
library,NVIDIA H200 NVL,gpu,cgbn,p1024,1024,MODMUL_R2,50000,0.000067008,746179560.649,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,p1024,1024,MODMUL_R2,50000,0.002992239,16709894.706,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,p1024,1024,MODMUL_R2,50000,0.003717521,13449823.247,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,p1024,1024,MODMUL_R2,50000,0.000266811,187398573.666,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,p1024,1024,MODMUL_R2,50000,0.000942733,53037285.014,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,p1024,1024,MODMUL_R2,50000,0.000215060,232493299.649,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,p1024,1024,MODMUL_R2,50000,0.001085923,46043778.007,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,p1024,1024,MODMUL_R2,50000,0.000145990,342489088.988,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,p1024,1024,MODMUL_R2,50000,0.000820522,60936812.315,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,p1024,1024,MODMUL_R2,50000,0.000179420,278675744.579,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,p1024,1024,MODMUL_R2,50000,0.000861692,58025372.696,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,p1024,1024,MODMUL_R2,50000,0.000107560,464857007.597,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,p1024,1024,MODMUL_R2,50000,0.000785962,63616308.769,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p1024,1024,MODMUL_R2,50000,0.018400763,2717278.619,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p1024,1024,MODMUL_R2,50000,0.091979086,543601.835,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p1024,1024,MODMUL_R2,50000,0.005478588,9126439.260,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p1024,1024,MODMUL_R2,50000,0.010562092,4733910.691,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,p1024,1024,MODMUL_R2,50000,0.003708197,13483641.647,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,p1024,1024,MODMUL_R2,50000,0.008659461,5774031.401,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,p2048,2048,ADD,6250,0.000133368,46862676.452,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,p2048,2048,ADD,6250,0.001154583,5413207.963,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,p2048,2048,ADD,6250,0.000697160,8964943.485,0
library,NVIDIA H200 NVL,gpu,cgbn,p2048,2048,ADD,50000,0.000022784,2194522471.910,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,p2048,2048,ADD,50000,0.000413281,120983051.709,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,p2048,2048,ADD,50000,0.001770966,28233179.149,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,p2048,2048,ADD,50000,0.000204601,244378111.065,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,p2048,2048,ADD,50000,0.001472505,33955741.110,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32,p2048,2048,ADD,50000,0.000080140,623907579.046,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32,p2048,2048,ADD,50000,0.001759405,28418700.004,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,p2048,2048,ADD,50000,0.000082270,607754883.443,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,p2048,2048,ADD,50000,0.001727695,28940293.369,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,p2048,2048,ADD,50000,0.000085080,587682692.168,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,p2048,2048,ADD,50000,0.001354133,36923994.589,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,p2048,2048,ADD,50000,0.000056281,888399482.056,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,p2048,2048,ADD,50000,0.001283103,38968032.415,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,p2048,2048,ADD,50000,0.000067320,742721840.512,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,p2048,2048,ADD,50000,0.001300463,38447845.490,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p2048,2048,ADD,50000,0.002129180,23483218.808,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p2048,2048,ADD,50000,0.009863174,5069362.085,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p2048,2048,ADD,50000,0.002075794,24087168.149,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p2048,2048,ADD,50000,0.010122890,4939300.930,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,p2048,2048,ADD,50000,0.002065760,24204166.011,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,p2048,2048,ADD,50000,0.009306634,5372511.649,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,p2048,2048,SUBTRACT,6250,0.000108944,57368717.149,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,p2048,2048,SUBTRACT,6250,0.000727438,8591795.369,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,p2048,2048,SUBTRACT,6250,0.000806124,7753145.870,0
library,NVIDIA H200 NVL,gpu,cgbn,p2048,2048,SUBTRACT,50000,0.000022976,2176183844.011,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,p2048,2048,SUBTRACT,50000,0.000413442,120935938.549,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,p2048,2048,SUBTRACT,50000,0.001775265,28164807.782,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,p2048,2048,SUBTRACT,50000,0.000203841,245289200.511,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,p2048,2048,SUBTRACT,50000,0.001462205,34194929.749,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32,p2048,2048,SUBTRACT,50000,0.000079950,625391115.738,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32,p2048,2048,SUBTRACT,50000,0.001806574,27676697.216,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,p2048,2048,SUBTRACT,50000,0.000081720,611845721.962,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,p2048,2048,SUBTRACT,50000,0.001795074,27854005.352,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,p2048,2048,SUBTRACT,50000,0.000085570,584316697.436,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,p2048,2048,SUBTRACT,50000,0.001333214,37503356.738,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,p2048,2048,SUBTRACT,50000,0.000055861,895077618.304,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,p2048,2048,SUBTRACT,50000,0.001288083,38817373.361,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,p2048,2048,SUBTRACT,50000,0.000067491,740838593.586,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,p2048,2048,SUBTRACT,50000,0.001292353,38689122.347,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p2048,2048,SUBTRACT,50000,0.002164590,23099062.042,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p2048,2048,SUBTRACT,50000,0.009894345,5053391.572,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p2048,2048,SUBTRACT,50000,0.002212745,22596367.644,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p2048,2048,SUBTRACT,50000,0.010277460,4865015.295,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,p2048,2048,SUBTRACT,50000,0.002046760,24428852.748,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,p2048,2048,SUBTRACT,50000,0.010037907,4981118.051,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,p2048,2048,SUBTRACT,50000,0.002083250,24000959.017,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,p2048,2048,SUBTRACT,50000,0.009828076,5087465.746,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,p2048,2048,ADDMOD,6250,0.000308109,20285017.430,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,p2048,2048,ADDMOD,6250,0.001383234,4518396.775,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,p2048,2048,ADDMOD,6250,0.001563660,3997033.880,0
library,NVIDIA H200 NVL,gpu,cgbn,p2048,2048,ADDMOD,50000,0.000023200,2155172413.793,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,p2048,2048,ADDMOD,50000,0.000516742,96760084.888,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,p2048,2048,ADDMOD,50000,0.002671018,18719454.460,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,p2048,2048,ADDMOD,50000,0.000263931,189443455.263,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,p2048,2048,ADDMOD,50000,0.001517545,32947952.854,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32,p2048,2048,ADDMOD,50000,0.000101500,492611041.453,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32,p2048,2048,ADDMOD,50000,0.001773454,28193570.032,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,p2048,2048,ADDMOD,50000,0.000110210,453679191.547,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,p2048,2048,ADDMOD,50000,0.001780244,28086036.759,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,p2048,2048,ADDMOD,50000,0.000111560,448189363.378,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,p2048,2048,ADDMOD,50000,0.001357113,36842915.504,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,p2048,2048,ADDMOD,50000,0.000034390,1453915207.155,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,p2048,2048,ADDMOD,50000,0.001263513,39572208.692,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,p2048,2048,ADDMOD,50000,0.000036340,1375895314.552,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,p2048,2048,ADDMOD,50000,0.001264674,39535883.815,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p2048,2048,ADDMOD,50000,0.002188869,22842846.592,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p2048,2048,ADDMOD,50000,0.010831329,4616238.673,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p2048,2048,ADDMOD,50000,0.002258715,22136480.383,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p2048,2048,ADDMOD,50000,0.010331661,4839492.852,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,p2048,2048,ADDMOD,50000,0.002250570,22216594.127,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,p2048,2048,ADDMOD,50000,0.009636585,5188560.108,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,p2048,2048,ADDMOD,50000,0.002295721,21779650.309,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,p2048,2048,ADDMOD,50000,0.009872346,5064652.337,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,p2048,2048,SUBTRACTMOD,6250,0.000275392,22694881.539,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,p2048,2048,SUBTRACTMOD,6250,0.000994825,6282514.093,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,p2048,2048,SUBTRACTMOD,6250,0.000853044,7326701.458,0
library,NVIDIA H200 NVL,gpu,cgbn,p2048,2048,SUBTRACTMOD,50000,0.000023136,2161134163.209,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,p2048,2048,SUBTRACTMOD,50000,0.000539311,92710889.112,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,p2048,2048,SUBTRACTMOD,50000,0.001915536,26102354.560,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,p2048,2048,SUBTRACTMOD,50000,0.000274821,181936560.886,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,p2048,2048,SUBTRACTMOD,50000,0.001506675,33185656.393,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32,p2048,2048,SUBTRACTMOD,50000,0.000111870,446947341.606,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32,p2048,2048,SUBTRACTMOD,50000,0.001852845,26985526.575,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,p2048,2048,SUBTRACTMOD,50000,0.000109660,455954879.338,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,p2048,2048,SUBTRACTMOD,50000,0.001863215,26835335.664,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,p2048,2048,SUBTRACTMOD,50000,0.000111711,447583582.155,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,p2048,2048,SUBTRACTMOD,50000,0.001363463,36671328.779,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,p2048,2048,SUBTRACTMOD,50000,0.000038240,1307531774.025,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,p2048,2048,SUBTRACTMOD,50000,0.001274224,39239567.892,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,p2048,2048,SUBTRACTMOD,50000,0.000034320,1456879200.559,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,p2048,2048,SUBTRACTMOD,50000,0.001254883,39844351.619,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p2048,2048,SUBTRACTMOD,50000,0.002172219,23017937.117,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p2048,2048,SUBTRACTMOD,50000,0.009788155,5108214.994,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p2048,2048,SUBTRACTMOD,50000,0.002002544,24968240.366,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p2048,2048,SUBTRACTMOD,50000,0.010026888,4986592.056,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,p2048,2048,SUBTRACTMOD,50000,0.002004640,24942134.945,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,p2048,2048,SUBTRACTMOD,50000,0.010276909,4865276.132,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,p2048,2048,SUBTRACTMOD,50000,0.002264211,22082746.899,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,p2048,2048,SUBTRACTMOD,50000,0.009767056,5119249.866,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.002022646,3090011.857,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.002016936,3098759.670,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.001991946,3137635.230,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.100475715,497632.687,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.102517000,487723.987,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.025420410,1966923.431,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.027236027,1835803.727,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.005723824,8735418.799,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.007891190,6336179.963,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.000413361,120959643.838,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.002489836,20081644.293,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.000397231,125871318.851,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.002059825,24273906.793,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.000363131,137691380.326,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.001991955,25100968.881,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.000283821,176167360.513,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.001921885,26016123.917,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.093984855,532000.608,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.091785445,544748.680,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.006050841,8263314.048,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.014637009,3415998.444,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.002491712,20066524.638,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.010846861,4609628.513,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.003838668,13025351.518,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.012513029,3995835.048,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.002027271,3082962.231,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.001959771,3189148.917,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.003001258,2082459.758,0
library,NVIDIA H200 NVL,gpu,cgbn,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.000071200,702247191.011,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.037143579,1346127.685,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.038995234,1282207.972,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.009443093,5294875.327,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.011179140,4472615.980,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.002388786,20931134.561,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.004503251,11103089.901,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.002423927,20627684.302,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.004493611,11126909.014,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.002417857,20679468.949,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.004068120,12290689.807,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.001044862,47853210.590,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.002667687,18742828.396,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.000981022,50967261.515,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.002602236,19214245.008,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.017453929,2864684.498,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.025603795,1952835.502,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.004501971,11106246.410,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.013113029,3813001.571,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.002083970,23992667.830,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.010862541,4602974.556,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.004376580,11424445.552,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.012912211,3872303.509,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.014674473,425909.674,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.004387093,1424633.575,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.004108932,1521076.528,0
library,NVIDIA H200 NVL,gpu,cgbn,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000109440,456871345.029,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.045409423,1101093.050,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.047223798,1058788.198,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.001563266,31984321.650,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.002858780,17489978.394,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000242750,205973284.942,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.001893975,26399502.665,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000431391,115904116.908,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.002231605,22405398.788,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000311041,160750562.856,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.001585714,31531536.152,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000384211,130136800.604,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.001624104,30786206.098,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000263921,189450641.711,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.001506793,33183059.127,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.290194121,172298.459,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.226453943,220795.449,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.018650327,2680918.142,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.092472951,540698.652,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.006544281,7640258.824,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.014286947,3499697.990,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.006108388,8185465.529,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.014363407,3481068.235,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,p2048,2048,COMPARE,6250,0.000033196,188277119.553,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,p2048,2048,COMPARE,6250,0.000254219,24585053.490,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,p2048,2048,COMPARE,6250,0.000485103,12883855.139,0
library,NVIDIA H200 NVL,gpu,cgbn,p2048,2048,COMPARE,50000,0.000023168,2158149171.271,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,p2048,2048,COMPARE,50000,0.000256081,195250725.024,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,p2048,2048,COMPARE,50000,0.002106576,23735198.800,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,p2048,2048,COMPARE,50000,0.000125741,397642760.949,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,p2048,2048,COMPARE,50000,0.001294934,38612004.262,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,p2048,2048,COMPARE,50000,0.000067160,744489487.104,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,p2048,2048,COMPARE,50000,0.001733184,28848639.647,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,p2048,2048,COMPARE,50000,0.000070370,710529699.540,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,p2048,2048,COMPARE,50000,0.001258724,39722766.239,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,p2048,2048,COMPARE,50000,0.000022310,2241140927.354,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,p2048,2048,COMPARE,50000,0.001270623,39350773.228,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,p2048,2048,COMPARE,50000,0.000024560,2035828626.955,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,p2048,2048,COMPARE,50000,0.001270973,39339938.549,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p2048,2048,COMPARE,50000,0.002088839,23936741.071,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p2048,2048,COMPARE,50000,0.006741210,7417066.072,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p2048,2048,COMPARE,50000,0.002110645,23689441.555,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p2048,2048,COMPARE,50000,0.006875257,7272455.407,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,p2048,2048,COMPARE,50000,0.002460971,20317183.672,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,p2048,2048,COMPARE,50000,0.007183164,6960720.922,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,p2048,2048,REDUCE,781,0.000014795,52789526.818,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,p2048,2048,REDUCE,781,0.000141601,5515497.705,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,p2048,2048,REDUCE,781,0.000114307,6832465.496,0
library,NVIDIA H200 NVL,gpu,cgbn,p2048,2048,REDUCE,50000,0.000031008,1612487100.103,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,p2048,2048,REDUCE,50000,0.141675795,352918.436,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,p2048,2048,REDUCE,50000,0.143903231,347455.715,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,p2048,2048,REDUCE,50000,0.005321868,9395197.331,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,p2048,2048,REDUCE,50000,0.006563723,7617628.017,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,p2048,2048,REDUCE,50000,0.001205673,41470615.176,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,p2048,2048,REDUCE,50000,0.002873377,17401127.302,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,p2048,2048,REDUCE,50000,0.001077213,46416076.620,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,p2048,2048,REDUCE,50000,0.002326055,21495622.747,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,p2048,2048,REDUCE,50000,0.001164303,42944146.048,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,p2048,2048,REDUCE,50000,0.002411506,20733931.219,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,p2048,2048,REDUCE,50000,0.001116453,44784689.797,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,p2048,2048,REDUCE,50000,0.002351406,21263873.340,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p2048,2048,REDUCE,50000,0.188850234,264760.064,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p2048,2048,REDUCE,50000,0.109257903,457632.799,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p2048,2048,REDUCE,50000,0.094629305,528377.547,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p2048,2048,REDUCE,50000,0.098645342,506866.305,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,p2048,2048,REDUCE,50000,0.099150347,504284.670,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,p2048,2048,REDUCE,50000,0.104497812,478478.918,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,p2048,2048,MODMUL,390,0.000390994,997456.606,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,p2048,2048,MODMUL,390,0.000441237,883878.958,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,p2048,2048,MODMUL,390,0.001141354,341699.288,0
library,NVIDIA H200 NVL,gpu,cgbn,p2048,2048,MODMUL,50000,0.000478880,104410290.678,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,p2048,2048,MODMUL,50000,0.268245715,186396.267,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,p2048,2048,MODMUL,50000,0.267945234,186605.297,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,p2048,2048,MODMUL,50000,0.016600459,3011964.910,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,p2048,2048,MODMUL,50000,0.017850434,2801052.335,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,p2048,2048,MODMUL,50000,0.005905105,8467250.066,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,p2048,2048,MODMUL,50000,0.007610540,6569836.024,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,p2048,2048,MODMUL,50000,0.004167941,11996330.899,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,p2048,2048,MODMUL,50000,0.005414503,9234457.824,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,p2048,2048,MODMUL,50000,0.005982285,8358010.499,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,p2048,2048,MODMUL,50000,0.007250079,6896476.510,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,p2048,2048,MODMUL,50000,0.003954030,12645326.442,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,p2048,2048,MODMUL,50000,0.005207773,9601033.012,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p2048,2048,MODMUL,50000,0.700911957,71335.636,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p2048,2048,MODMUL,50000,0.707853118,70636.123,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p2048,2048,MODMUL,50000,0.404158345,123713.887,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p2048,2048,MODMUL,50000,0.401489247,124536.337,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,p2048,2048,MODMUL,50000,0.293022299,170635.478,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w32-opt,p2048,2048,MODMUL,50000,0.297988232,167791.861,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,p2048,2048,MODEXP,97,0.133397421,727.150,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,p2048,2048,MODEXP,97,0.005736656,16908.805,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,p2048,2048,MODEXP,97,0.006636439,14616.272,0
library,NVIDIA H200 NVL,gpu,cgbn,p2048,2048,MODEXP,50000,0.701264918,71299.731,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,p2048,2048,MODEXP,50000,43.443485951,1150.921,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,p2048,2048,MODEXP,50000,43.445529175,1150.866,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,p2048,2048,MODEXP,50000,6.308410606,7925.927,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,p2048,2048,MODEXP,50000,6.310054951,7923.861,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,p2048,2048,MODEXP,50000,1.344470883,37189.351,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,p2048,2048,MODEXP,50000,1.346453988,37134.578,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,p2048,2048,MODEXP,50000,1.573074507,31784.890,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,p2048,2048,MODEXP,50000,1.571971034,31807.202,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,p2048,2048,MODEXP,50000,1.378118909,36281.339,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,p2048,2048,MODEXP,50000,1.379204862,36252.772,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,p2048,2048,MODEXP,50000,1.566792763,31912.325,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,p2048,2048,MODEXP,50000,1.574413336,31757.861,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,p2048,2048,EXPONENTIATION,97,0.015306585,6337.142,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,p2048,2048,EXPONENTIATION,97,0.004098502,23667.184,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,p2048,2048,EXPONENTIATION,97,0.010775212,9002.143,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,p2048,2048,EXPONENTIATION,50000,15.110040283,3309.058,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,p2048,2048,EXPONENTIATION,50000,15.104804584,3310.205,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,p2048,2048,EXPONENTIATION,50000,3.482968537,14355.570,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,p2048,2048,EXPONENTIATION,50000,3.481442992,14361.861,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,p2048,2048,EXPONENTIATION,50000,1.025173357,48772.239,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,p2048,2048,EXPONENTIATION,50000,1.027237833,48674.220,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,p2048,2048,EXPONENTIATION,50000,0.894597114,55891.081,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,p2048,2048,EXPONENTIATION,50000,0.896082648,55798.425,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,p2048,2048,EXPONENTIATION,50000,1.049814841,47627.446,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,p2048,2048,EXPONENTIATION,50000,1.051185144,47565.360,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,p2048,2048,EXPONENTIATION,50000,0.895522047,55833.355,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,p2048,2048,EXPONENTIATION,50000,0.897752313,55694.649,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p2048,2048,EXPONENTIATION,50000,45.653449844,1095.207,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p2048,2048,EXPONENTIATION,50000,45.734481762,1093.267,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,p2048,2048,DIVIDE,781,0.000024391,32020065.993,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,p2048,2048,DIVIDE,781,0.000046542,16780614.230,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,p2048,2048,DIVIDE,781,0.000048234,16191812.056,0
library,NVIDIA H200 NVL,gpu,cgbn,p2048,2048,DIVIDE,50000,0.000036704,1362249346.120,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,p2048,2048,DIVIDE,50000,0.614563017,81358.622,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,p2048,2048,DIVIDE,50000,0.616451314,81109.406,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,p2048,2048,DIVIDE,50000,0.156840686,318794.831,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,p2048,2048,DIVIDE,50000,0.157984779,316486.185,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,p2048,2048,DIVIDE,50000,0.012916873,3870905.914,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,p2048,2048,DIVIDE,50000,0.015065398,3318863.542,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,p2048,2048,DIVIDE,50000,0.011435659,4372288.463,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,p2048,2048,DIVIDE,50000,0.013496384,3704696.007,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,p2048,2048,DIVIDE,50000,0.012868233,3885537.381,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,p2048,2048,DIVIDE,50000,0.014536887,3439525.932,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,p2048,2048,DIVIDE,50000,0.011860722,4215594.978,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,p2048,2048,DIVIDE,50000,0.013453066,3716624.894,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p2048,2048,DIVIDE,50000,0.219959430,227314.646,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p2048,2048,DIVIDE,50000,0.215559289,231954.745,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p2048,2048,DIVIDE,50000,0.201515688,248119.640,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p2048,2048,DIVIDE,50000,0.203321207,245916.305,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,p2048,2048,ISQRT,195,0.000050395,3869453.425,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,p2048,2048,ISQRT,195,0.000066625,2926807.968,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,p2048,2048,ISQRT,50000,9.192383292,5439.286,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,p2048,2048,ISQRT,50000,9.187155403,5442.381,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,p2048,2048,ISQRT,50000,4.357044334,11475.669,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,p2048,2048,ISQRT,50000,4.442244185,11255.572,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,p2048,2048,ISQRT,50000,0.168405620,296902.206,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,p2048,2048,ISQRT,50000,0.170067305,294001.249,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,p2048,2048,ISQRT,50000,0.057844988,864379.124,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,p2048,2048,ISQRT,50000,0.059617002,838686.923,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,p2048,2048,ISQRT,50000,0.169120952,295646.396,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,p2048,2048,ISQRT,50000,0.170462775,293319.172,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,p2048,2048,ISQRT,50000,0.058143027,859948.348,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,p2048,2048,ISQRT,50000,0.059429604,841331.536,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p2048,2048,ISQRT,50000,5.163005145,9684.282,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p2048,2048,ISQRT,50000,5.178885793,9654.586,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p2048,2048,ISQRT,50000,3.492400971,14316.798,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p2048,2048,ISQRT,50000,3.505534004,14263.162,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-1t,p2048,2048,MODMUL_R2,6250,0.006472639,965603.067,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,gmp-nt,p2048,2048,MODMUL_R2,6250,0.003198609,1953974.304,0
library,AMD EPYC 9655 96-Core Processor,host-cpu,openssl-nt,p2048,2048,MODMUL_R2,6250,0.006998851,893003.713,0
library,NVIDIA H200 NVL,gpu,cgbn,p2048,2048,MODMUL_R2,50000,0.000205376,243455905.266,0
opencl-kernel,NVIDIA H200 NVL,GPU,w8,p2048,2048,MODMUL_R2,50000,0.023804104,2100478.140,0
opencl-e2e,NVIDIA H200 NVL,GPU,w8,p2048,2048,MODMUL_R2,50000,0.025112669,1991026.921,0
opencl-kernel,NVIDIA H200 NVL,GPU,w16,p2048,2048,MODMUL_R2,50000,0.001331584,37549266.158,0
opencl-e2e,NVIDIA H200 NVL,GPU,w16,p2048,2048,MODMUL_R2,50000,0.002593436,19279442.069,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-opt,p2048,2048,MODMUL_R2,50000,0.000773242,64662811.298,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-opt,p2048,2048,MODMUL_R2,50000,0.002534716,19726076.409,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-o64,p2048,2048,MODMUL_R2,50000,0.000528042,94689442.631,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-o64,p2048,2048,MODMUL_R2,50000,0.002228176,22439879.193,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il,p2048,2048,MODMUL_R2,50000,0.000694062,72039675.898,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il,p2048,2048,MODMUL_R2,50000,0.001976595,25296026.403,0
opencl-kernel,NVIDIA H200 NVL,GPU,w32-il64,p2048,2048,MODMUL_R2,50000,0.000469240,106555278.446,0
opencl-e2e,NVIDIA H200 NVL,GPU,w32-il64,p2048,2048,MODMUL_R2,50000,0.001694066,29514788.953,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p2048,2048,MODMUL_R2,50000,0.298339304,167594.411,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w8,p2048,2048,MODMUL_R2,50000,0.305189473,163832.650,0
opencl-kernel,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p2048,2048,MODMUL_R2,50000,0.017449333,2865439.043,0
opencl-e2e,cpu-skylake-avx512-AMD EPYC 9655 96-Core Processor,CPU,w16,p2048,2048,MODMUL_R2,50000,0.092075353,543033.487,0
```
