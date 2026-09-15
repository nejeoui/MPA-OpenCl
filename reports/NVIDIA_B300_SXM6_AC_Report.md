# MPA-OpenCL benchmark report - NVIDIA B300 SXM6 AC

> **Note.** The multi-threaded GMP and OpenSSL baseline columns have been
> removed from this report: they predate the 2026-09-12 timing fix and were
> understated (see `reports/README.md`). The single-threaded GMP column, the
> OpenCL-on-CPU rows and all MPA measurements are unaffected and were verified
> against GMP before timing.


> **CGBN column rebuilt offline.** The run that produced this report loaded the
> wrong `cgbn_results` file. Its CGBN column, and the `GPU vs CGBN` ratios, were
> regenerated from `cgbn_results_NVIDIA_B300_SXM6_AC.tsv` (device `NVIDIA B300 SXM6 AC`, 60 rows)
> by `rebuild_cgbn_column.py`, which also recomputed section 5 from the CSV:
> as written, that section did not normalise the GPU time to the CPU item
> count, so its figures were out by the ratio of the two. Section 4, and every
> measured value in the CSV, are untouched.

> **Partial report.** The run was interrupted or hit its time budget.
> Rows that never ran are marked `n/a`.

## 1. System under test

2 OpenCL device(s) exercised with the identical kernels and operands.

### Device 0 - NVIDIA B300 SXM6 AC (GPU)

| Property | Value |
|---|---|
| Model | NVIDIA B300 SXM6 AC |
| Type | GPU |
| Vendor | NVIDIA Corporation |
| Device memory | 267.68 GiB |
| Max single allocation | 66.92 GiB |
| Local memory | 48 KiB |
| Global cache | 4736 KiB |
| Compute units | 148 |
| Max clock | 2032 MHz |
| Max work-group size | 1024 |
| OpenCL version | OpenCL 3.0 CUDA |
| Driver | 595.84 |

### Device 1 - cpu-skylake-avx512-Intel(R) Xeon(R) 6787P (CPU)

| Property | Value |
|---|---|
| Model | cpu-skylake-avx512-Intel(R) Xeon(R) 6787P |
| Type | CPU |
| Vendor | GenuineIntel |
| Device memory | 2013.37 GiB |
| Max single allocation | 512.00 GiB |
| Local memory | 2048 KiB |
| Global cache | 344064 KiB |
| Compute units | 344 |
| Max clock | 2000 MHz |
| Max work-group size | 4096 |
| OpenCL version | OpenCL 3.0 PoCL HSTR: cpu-x86_64-pc-linux-gnu-skylake-avx512 |
| Driver | 5.0+debian |

### Host

| Property | Value |
|---|---|
| CPU | Intel(R) Xeon(R) 6787P |
| Logical cores | 344 |
| OpenMP threads used | 344 |
| RAM | 2015.4 GB |
| OS | Ubuntu 24.04.4 LTS |
| Kernel | 7.0.0-31-generic |
| Arch | x86_64 |
| GMP | 6.3.0 |
| OpenSSL | OpenSSL 3.0.13 30 Jan 2024 |
| CGBN | 60 rows from `cgbn_results_NVIDIA_B300_SXM6_AC.tsv` |

## 2. Method

- Workload auto-sized from the device and host: --min-items from 700 x compute units, --items from ten times that capped by host RAM. Either flag, given explicitly, overrides its half.
- Base workload 50000 items, scaled down per operator by its cost weight and by modulus size. Device rows honour --min-items (103600) so the GPU is not left idle; the CPU libraries keep the smaller count because a full-width MODEXP there costs minutes. Both counts appear in every row as dev/cpu, and throughput is per-second so they remain comparable.
- 5 timed repetitions, **minimum** reported. Two untimed warm-up launches precede them.
- `kernel` times `clEnqueueNDRangeKernel` + `clFinish` only. `e2e` adds the host->device operand writes and the device->host result read.
- Every OpenCL device runs the same kernels on the same operands, so GPU and CPU-OpenCL columns are directly comparable.
- CPU library baselines (GMP, OpenSSL) run those same operands, with every temporary - including each thread's GMP context, BN_CTX and Montgomery context - allocated outside the timed region, so the figure is the arithmetic and not marshalling. The generator is reseeded per modulus and operation so every backend sees identical inputs.
- Cost weighting drives the wide cells down to a few hundred items, which is tens of microseconds of work - the same order as the cost of entering an OpenMP region. Each baseline pass is therefore repeated until the timed interval reaches 5 ms and the per-pass time is reported; the multi-threaded loop enters one parallel region per interval and partitions the range itself. Without this the multi-threaded GMP figure came out up to 9x slower than the single-threaded one at 2048 bits.
- OpenSSL rows time the nearest BN primitive, which is not always semantically identical (its Montgomery routine expects Montgomery-domain inputs); they measure comparable work, not identical results. Correctness is judged against GMP only.
- Every device cell is checked word-for-word against GMP before it is timed. A cell that mismatches is reported and excluded from the speedup tables.
- Total wall time 2331.6 s.

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

**All configurations correct** - 485 configurations, 0 problems.

## 4. Throughput per device

Operations per second, higher is better. Kernel-only timings.

### Device 0 - NVIDIA B300 SXM6 AC (GPU)

#### secp256k1 (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | 2.17 G | 2.84 G | 3.49 G | 3.74 G | 3.73 G | 4.02 G | 3.92 G | 52.74 M | 6.40 G |
| SUBTRACT | 50000 / 50000 | 2.18 G | 2.86 G | 3.70 G | 3.65 G | 3.68 G | 4.06 G | 4.11 G | 31.29 M | 6.57 G |
| ADDMOD | 50000 / 50000 | 1.71 G | 2.35 G | 3.40 G | 3.56 G | 3.63 G | 4.24 G | 4.40 G | 10.35 M | 5.16 G |
| SUBTRACTMOD | 50000 / 50000 | 1.70 G | 2.39 G | 3.41 G | 3.60 G | 3.67 G | 4.08 G | 4.37 G | 12.04 M | 5.14 G |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | 29.22 M | 81.67 M | 611.53 M | 2.82 G | 3.04 G | 3.93 G | 3.84 G | 20.13 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | 299.10 M | 852.81 M | 2.15 G | 2.13 G | 2.16 G | 2.46 G | 2.51 G | 19.66 M | 6.40 G |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | 682.76 M | 1.80 G | 3.36 G | 3.02 G | 3.48 G | 3.13 G | 3.87 G | 4.21 M | 4.82 G |
| COMPARE | 50000 / 50000 | 2.16 G | 2.82 G | - | 3.97 G | 3.91 G | 4.02 G | 3.96 G | 57.12 M | 6.40 G |
| REDUCE | 50000 / 6250 | 370.91 M | 540.53 M | - | 1.36 G | 1.40 G | 1.39 G | 1.44 G | 10.58 M | 3.67 G |
| MODMUL | 50000 / 3125 | 137.89 M | 230.93 M | - | 477.24 M | 590.91 M | 488.86 M | 591.78 M | 9.22 M | 1.60 G |
| MODEXP | 50000 / 781 | 3.01 M | 16.20 M | - | 20.16 M | 42.09 M | 20.32 M | 42.20 M | 89.45 k | 5.69 M |
| EXPONENTIATION | 50000 / 781 | 2.21 M | 10.95 M | - | 110.65 M | 195.46 M | 109.70 M | 192.27 M | 161.05 k | n/a |
| DIVIDE | 50000 / 6250 | 152.02 M | 178.44 M | - | 466.07 M | 476.65 M | 474.16 M | 491.96 M | 4.45 M | 3.18 G |
| ISQRT | 50000 / 1562 | 13.95 M | 18.72 M | - | 75.94 M | 82.50 M | 74.24 M | 83.46 M | 2.40 M | n/a |
| MODMUL_R2 | 50000 / 50000 | 573.47 M | 1.72 G | - | 2.40 G | 3.07 G | 2.45 G | 3.16 G | 5.32 M | 3.59 G |

#### rsa256(composite) (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | 2.21 G | 2.80 G | 3.72 G | 3.62 G | 3.81 G | 3.89 G | 3.87 G | 18.97 M | 6.33 G |
| SUBTRACT | 50000 / 50000 | 2.25 G | 2.85 G | 3.56 G | 3.68 G | 3.59 G | 4.05 G | 4.09 G | 32.44 M | 6.59 G |
| ADDMOD | 50000 / 50000 | 1.81 G | 2.51 G | 3.57 G | 3.82 G | 3.71 G | 4.23 G | 4.21 G | 12.16 M | 5.28 G |
| SUBTRACTMOD | 50000 / 50000 | 1.63 G | 2.34 G | 3.36 G | 3.67 G | 3.65 G | 4.09 G | 4.40 G | 12.07 M | 5.21 G |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | 29.70 M | 80.96 M | 616.99 M | 2.84 G | 2.99 G | 3.57 G | 3.89 G | 21.12 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | 298.71 M | 864.38 M | 2.06 G | 2.10 G | 2.19 G | 2.44 G | 2.53 G | 21.06 M | 6.18 G |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | 682.43 M | 1.83 G | 3.29 G | 2.93 G | 3.34 G | 3.15 G | 3.79 G | 4.21 M | 4.96 G |
| COMPARE | 50000 / 50000 | 2.18 G | 2.86 G | - | 3.82 G | 3.86 G | 4.02 G | 4.17 G | 57.61 M | 6.38 G |
| REDUCE | 50000 / 6250 | 370.23 M | 535.66 M | - | 1.36 G | 1.38 G | 1.38 G | 1.43 G | 7.22 M | 3.78 G |
| MODMUL | 50000 / 3125 | 137.59 M | 230.15 M | - | 477.87 M | 588.33 M | 486.53 M | 592.75 M | 2.02 M | 1.58 G |
| MODEXP | 50000 / 781 | 3.01 M | 16.23 M | - | 20.20 M | 42.15 M | 20.36 M | 42.31 M | 44.66 k | 5.82 M |
| EXPONENTIATION | 50000 / 781 | 2.21 M | 10.92 M | - | 110.61 M | 195.61 M | 109.75 M | 192.34 M | 60.87 k | n/a |
| DIVIDE | 50000 / 6250 | 145.88 M | 172.02 M | - | 433.42 M | 443.88 M | 438.63 M | 454.18 M | 5.24 M | 3.18 G |
| ISQRT | 50000 / 1562 | 13.95 M | 18.12 M | - | 70.51 M | 73.93 M | 69.33 M | 74.67 M | 2.38 M | n/a |
| MODMUL_R2 | 50000 / 50000 | 570.56 M | 1.75 G | - | 2.35 G | 2.95 G | 2.41 G | 3.11 G | 5.72 M | 3.46 G |

#### brainpoolP512r1 (512-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 25000 | 879.57 M | 1.53 G | 2.60 G | 2.72 G | 2.68 G | 3.46 G | 2.70 G | 17.38 M | 5.94 G |
| SUBTRACT | 50000 / 25000 | 908.64 M | 1.54 G | 2.59 G | 2.62 G | 2.62 G | 3.49 G | 2.75 G | 24.84 M | 5.77 G |
| ADDMOD | 50000 / 25000 | 712.15 M | 1.33 G | 2.41 G | 2.47 G | 2.48 G | 3.61 G | 3.52 G | 9.79 M | 4.93 G |
| SUBTRACTMOD | 50000 / 25000 | 668.44 M | 1.25 G | 2.33 G | 2.47 G | 2.47 G | 3.62 G | 3.47 G | 9.89 M | 4.66 G |
| MULTIPLYOPERANDSCANNING | 50000 / 25000 | 10.59 M | 25.86 M | 57.25 M | 1.37 G | 1.43 G | 1.88 G | 1.73 G | 7.83 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 25000 | 45.13 M | 161.88 M | 553.42 M | 552.71 M | 564.68 M | 1.09 G | 856.81 M | 8.09 M | 4.84 G |
| MONTGOMERYMULTIPLICATION | 50000 / 25000 | 205.02 M | 666.86 M | 1.91 G | 1.37 G | 1.87 G | 1.51 G | 2.35 G | 1.91 M | 3.65 G |
| COMPARE | 50000 / 25000 | 1.23 G | 2.04 G | - | 3.00 G | 3.04 G | 3.46 G | 3.65 G | 48.88 M | 5.90 G |
| REDUCE | 50000 / 3125 | 129.68 M | 157.63 M | - | 542.52 M | 554.32 M | 549.87 M | 493.98 M | 5.97 M | 2.86 G |
| MODMUL | 50000 / 1562 | 44.39 M | 62.91 M | - | 155.90 M | 183.25 M | 157.26 M | 148.95 M | 928.33 k | 629.28 M |
| MODEXP | 50000 / 390 | 216.72 k | 2.32 M | - | 2.62 M | 6.08 M | 2.60 M | 5.97 M | 19.18 k | 2.54 M |
| EXPONENTIATION | 50000 / 390 | 273.99 k | 1.09 M | - | 4.40 M | 5.96 M | 4.46 M | 4.34 M | 24.89 k | n/a |
| DIVIDE | 50000 / 3125 | 44.66 M | 43.92 M | - | 129.99 M | 140.37 M | 135.00 M | 148.40 M | 4.76 M | 2.09 G |
| ISQRT | 50000 / 781 | 2.94 M | 2.95 M | - | 15.14 M | 16.65 M | 15.88 M | 12.63 M | 1.17 M | n/a |
| MODMUL_R2 | 50000 / 25000 | 128.67 M | 694.78 M | - | 911.63 M | 1.59 G | 979.64 M | 1.62 G | 2.64 M | 2.50 G |

#### p1024 (1024-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 12500 | 269.83 M | 513.57 M | 1.25 G | 1.21 G | 1.19 G | 1.84 G | 1.84 G | 12.94 M | 4.87 G |
| SUBTRACT | 50000 / 12500 | 271.66 M | 510.53 M | 1.23 G | 1.21 G | 1.21 G | 1.79 G | 1.84 G | 21.58 M | 5.04 G |
| ADDMOD | 50000 / 12500 | 206.43 M | 394.46 M | 1.08 G | 1.25 G | 1.26 G | 3.02 G | 2.99 G | 5.86 M | 4.29 G |
| SUBTRACTMOD | 50000 / 12500 | 198.02 M | 394.86 M | 1.09 G | 1.29 G | 1.26 G | 2.89 G | 2.91 G | 7.04 M | 4.28 G |
| MULTIPLYOPERANDSCANNING | 50000 / 12500 | 2.26 M | 7.35 M | 11.45 M | 555.73 M | 630.50 M | 709.58 M | 816.45 M | 2.39 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 12500 | 6.06 M | 23.55 M | 91.15 M | 91.43 M | 91.39 M | 259.43 M | 277.36 M | 2.22 M | 2.41 G |
| MONTGOMERYMULTIPLICATION | 50000 / 12500 | 27.99 M | 202.49 M | 781.41 M | 453.84 M | 658.45 M | 547.95 M | 912.99 M | 752.63 k | 1.65 G |
| COMPARE | 50000 / 12500 | 579.84 M | 1.07 G | - | 2.33 G | 2.25 G | 3.04 G | 2.90 G | 34.82 M | 4.65 G |
| REDUCE | 50000 / 1562 | 23.69 M | 44.57 M | - | 161.34 M | 151.52 M | 165.88 M | 158.79 M | 9.14 M | 2.18 G |
| MODMUL | 50000 / 781 | 6.28 M | 16.87 M | - | 35.64 M | 48.48 M | 34.95 M | 47.59 M | 299.39 k | 267.37 M |
| MODEXP | 50000 / 195 | 19.20 k | 332.15 k | - | 330.02 k | 633.84 k | 326.50 k | 632.64 k | 4.07 k | 506.63 k |
| EXPONENTIATION | 50000 / 195 | 31.80 k | 139.04 k | - | 470.11 k | 537.27 k | 470.80 k | 534.46 k | 11.50 k | n/a |
| DIVIDE | 50000 / 1562 | 5.04 M | 7.04 M | - | 33.22 M | 32.52 M | 33.21 M | 32.19 M | 4.47 M | 1.66 G |
| ISQRT | 50000 / 390 | 263.53 k | 416.67 k | - | 2.19 M | 2.28 M | 2.20 M | 2.25 M | 761.52 k | n/a |
| MODMUL_R2 | 50000 / 12500 | 20.72 M | 225.19 M | - | 278.87 M | 422.88 M | 313.01 M | 518.79 M | 1.09 M | 919.66 M |

#### p2048 (2048-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 6250 | 140.28 M | 274.63 M | 692.74 M | 699.38 M | 676.13 M | 1.22 G | 1.15 G | 10.42 M | 3.68 G |
| SUBTRACT | 50000 / 6250 | 140.74 M | 272.63 M | 692.20 M | 686.69 M | 673.01 M | 1.22 G | 1.16 G | 12.10 M | 4.01 G |
| ADDMOD | 50000 / 6250 | 110.64 M | 217.88 M | 568.85 M | 724.63 M | 711.24 M | 2.16 G | 2.08 G | 4.59 M | 4.12 G |
| SUBTRACTMOD | 50000 / 6250 | 108.67 M | 211.25 M | 545.90 M | 731.31 M | 706.59 M | 2.14 G | 1.92 G | 5.28 M | 3.77 G |
| MULTIPLYOPERANDSCANNING | 50000 / 6250 | 564.49 k | 2.15 M | 7.15 M | 201.50 M | 225.92 M | 234.39 M | 267.20 M | 880.32 k | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 6250 | 1.53 M | 6.04 M | 23.79 M | 23.86 M | 23.86 M | 60.35 M | 64.87 M | 956.26 k | 879.29 M |
| MONTGOMERYMULTIPLICATION | 50000 / 6250 | 1.29 M | 39.27 M | 255.11 M | 145.06 M | 190.45 M | 163.25 M | 241.03 M | 316.70 k | 576.78 M |
| COMPARE | 50000 / 6250 | 315.72 M | 607.56 M | - | 1.54 G | 1.46 G | 2.05 G | 2.06 G | 45.00 M | 4.05 G |
| REDUCE | 50000 / 781 | 473.08 k | 10.69 M | - | 52.28 M | 62.96 M | 57.50 M | 58.34 M | 5.65 M | 2.05 G |
| MODMUL | 50000 / 390 | 286.80 k | 3.34 M | - | 9.67 M | 14.13 M | 9.58 M | 14.25 M | 138.07 k | 133.64 M |
| MODEXP | 50000 / 97 | 898.5 | 9.72 k | - | 30.34 k | 50.46 k | 30.03 k | 50.52 k | 587.5 | 90.79 k |
| EXPONENTIATION | 50000 / 97 | 3.80 k | 16.10 k | - | 56.17 k | 63.82 k | 54.96 k | 63.31 k | 4.59 k | n/a |
| DIVIDE | 50000 / 781 | 121.88 k | 342.31 k | - | 3.49 M | 3.80 M | 3.49 M | 3.82 M | 2.88 M | 1.69 G |
| ISQRT | 50000 / 195 | 7.22 k | 11.19 k | - | 332.68 k | 979.73 k | 337.03 k | 979.41 k | 361.86 k | n/a |
| MODMUL_R2 | 50000 / 6250 | 1.73 M | 45.78 M | - | 80.87 M | 118.18 M | 87.45 M | 131.21 M | 446.29 k | 306.37 M |

### Device 1 - cpu-skylake-avx512-Intel(R) Xeon(R) 6787P (CPU)

#### secp256k1 (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | - | - | - | - | - | - | - | 52.74 M | 6.40 G |
| SUBTRACT | 50000 / 50000 | - | - | - | - | - | - | - | 31.29 M | 6.57 G |
| ADDMOD | 50000 / 50000 | - | - | - | - | - | - | - | 10.35 M | 5.16 G |
| SUBTRACTMOD | 50000 / 50000 | - | - | - | - | - | - | - | 12.04 M | 5.14 G |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | - | - | - | - | - | - | - | 20.13 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | - | - | - | - | - | - | - | 19.66 M | 6.40 G |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | - | - | - | - | - | - | - | 4.21 M | 4.82 G |
| COMPARE | 50000 / 50000 | - | - | - | - | - | - | - | 57.12 M | 6.40 G |
| REDUCE | 50000 / 6250 | - | - | - | - | - | - | - | 10.58 M | 3.67 G |
| MODMUL | 50000 / 3125 | - | - | - | - | - | - | - | 9.22 M | 1.60 G |
| MODEXP | 50000 / 781 | - | - | - | - | - | - | - | 89.45 k | 5.69 M |
| EXPONENTIATION | 50000 / 781 | - | - | - | - | - | - | - | 161.05 k | n/a |
| DIVIDE | 50000 / 6250 | - | - | - | - | - | - | - | 4.45 M | 3.18 G |
| ISQRT | 50000 / 1562 | - | - | - | - | - | - | - | 2.40 M | n/a |
| MODMUL_R2 | 50000 / 50000 | - | - | - | - | - | - | - | 5.32 M | 3.59 G |

#### rsa256(composite) (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | - | - | - | - | - | - | - | 18.97 M | 6.33 G |
| SUBTRACT | 50000 / 50000 | - | - | - | - | - | - | - | 32.44 M | 6.59 G |
| ADDMOD | 50000 / 50000 | - | - | - | - | - | - | - | 12.16 M | 5.28 G |
| SUBTRACTMOD | 50000 / 50000 | - | - | - | - | - | - | - | 12.07 M | 5.21 G |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | - | - | - | - | - | - | - | 21.12 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | - | - | - | - | - | - | - | 21.06 M | 6.18 G |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | - | - | - | - | - | - | - | 4.21 M | 4.96 G |
| COMPARE | 50000 / 50000 | - | - | - | - | - | - | - | 57.61 M | 6.38 G |
| REDUCE | 50000 / 6250 | - | - | - | - | - | - | - | 7.22 M | 3.78 G |
| MODMUL | 50000 / 3125 | - | - | - | - | - | - | - | 2.02 M | 1.58 G |
| MODEXP | 50000 / 781 | - | - | - | - | - | - | - | 44.66 k | 5.82 M |
| EXPONENTIATION | 50000 / 781 | - | - | - | - | - | - | - | 60.87 k | n/a |
| DIVIDE | 50000 / 6250 | - | - | - | - | - | - | - | 5.24 M | 3.18 G |
| ISQRT | 50000 / 1562 | - | - | - | - | - | - | - | 2.38 M | n/a |
| MODMUL_R2 | 50000 / 50000 | - | - | - | - | - | - | - | 5.72 M | 3.46 G |

#### brainpoolP512r1 (512-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 25000 | - | - | - | - | - | - | - | 17.38 M | 5.94 G |
| SUBTRACT | 50000 / 25000 | - | - | - | - | - | - | - | 24.84 M | 5.77 G |
| ADDMOD | 50000 / 25000 | - | - | - | - | - | - | - | 9.79 M | 4.93 G |
| SUBTRACTMOD | 50000 / 25000 | - | - | - | - | - | - | - | 9.89 M | 4.66 G |
| MULTIPLYOPERANDSCANNING | 50000 / 25000 | - | - | - | - | - | - | - | 7.83 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 25000 | - | - | - | - | - | - | - | 8.09 M | 4.84 G |
| MONTGOMERYMULTIPLICATION | 50000 / 25000 | - | - | - | - | - | - | - | 1.91 M | 3.65 G |
| COMPARE | 50000 / 25000 | - | - | - | - | - | - | - | 48.88 M | 5.90 G |
| REDUCE | 50000 / 3125 | - | - | - | - | - | - | - | 5.97 M | 2.86 G |
| MODMUL | 50000 / 1562 | - | - | - | - | - | - | - | 928.33 k | 629.28 M |
| MODEXP | 50000 / 390 | - | - | - | - | - | - | - | 19.18 k | 2.54 M |
| EXPONENTIATION | 50000 / 390 | - | - | - | - | - | - | - | 24.89 k | n/a |
| DIVIDE | 50000 / 3125 | - | - | - | - | - | - | - | 4.76 M | 2.09 G |
| ISQRT | 50000 / 781 | - | - | - | - | - | - | - | 1.17 M | n/a |
| MODMUL_R2 | 50000 / 25000 | - | - | - | - | - | - | - | 2.64 M | 2.50 G |

#### p1024 (1024-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 12500 | - | - | - | - | - | - | - | 12.94 M | 4.87 G |
| SUBTRACT | 50000 / 12500 | - | - | - | - | - | - | - | 21.58 M | 5.04 G |
| ADDMOD | 50000 / 12500 | - | - | - | - | - | - | - | 5.86 M | 4.29 G |
| SUBTRACTMOD | 50000 / 12500 | - | - | - | - | - | - | - | 7.04 M | 4.28 G |
| MULTIPLYOPERANDSCANNING | 50000 / 12500 | - | - | - | - | - | - | - | 2.39 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 12500 | - | - | - | - | - | - | - | 2.22 M | 2.41 G |
| MONTGOMERYMULTIPLICATION | 50000 / 12500 | - | - | - | - | - | - | - | 752.63 k | 1.65 G |
| COMPARE | 50000 / 12500 | - | - | - | - | - | - | - | 34.82 M | 4.65 G |
| REDUCE | 50000 / 1562 | - | - | - | - | - | - | - | 9.14 M | 2.18 G |
| MODMUL | 50000 / 781 | - | - | - | - | - | - | - | 299.39 k | 267.37 M |
| MODEXP | 50000 / 195 | - | - | - | - | - | - | - | 4.07 k | 506.63 k |
| EXPONENTIATION | 50000 / 195 | - | - | - | - | - | - | - | 11.50 k | n/a |
| DIVIDE | 50000 / 1562 | - | - | - | - | - | - | - | 4.47 M | 1.66 G |
| ISQRT | 50000 / 390 | - | - | - | - | - | - | - | 761.52 k | n/a |
| MODMUL_R2 | 50000 / 12500 | - | - | - | - | - | - | - | 1.09 M | 919.66 M |

#### p2048 (2048-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 6250 | - | - | - | - | - | - | - | 10.42 M | 3.68 G |
| SUBTRACT | 50000 / 6250 | - | - | - | - | - | - | - | 12.10 M | 4.01 G |
| ADDMOD | 50000 / 6250 | - | - | - | - | - | - | - | 4.59 M | 4.12 G |
| SUBTRACTMOD | 50000 / 6250 | - | - | - | - | - | - | - | 5.28 M | 3.77 G |
| MULTIPLYOPERANDSCANNING | 50000 / 6250 | - | - | - | - | - | - | - | 880.32 k | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 6250 | - | - | - | - | - | - | - | 956.26 k | 879.29 M |
| MONTGOMERYMULTIPLICATION | 50000 / 6250 | - | - | - | - | - | - | - | 316.70 k | 576.78 M |
| COMPARE | 50000 / 6250 | - | - | - | - | - | - | - | 45.00 M | 4.05 G |
| REDUCE | 50000 / 781 | - | - | - | - | - | - | - | 5.65 M | 2.05 G |
| MODMUL | 50000 / 390 | - | - | - | - | - | - | - | 138.07 k | 133.64 M |
| MODEXP | 50000 / 97 | - | - | - | - | - | - | - | 587.5 | 90.79 k |
| EXPONENTIATION | 50000 / 97 | - | - | - | - | - | - | - | 4.59 k | n/a |
| DIVIDE | 50000 / 781 | - | - | - | - | - | - | - | 2.88 M | 1.69 G |
| ISQRT | 50000 / 195 | - | - | - | - | - | - | - | 361.86 k | n/a |
| MODMUL_R2 | 50000 / 6250 | - | - | - | - | - | - | - | 446.29 k | 306.37 M |

## 5. Head to head

Best OpenCL GPU result against best OpenCL CPU result and the CPU libraries.
Ratios above 1.00x mean the GPU is faster than that baseline.

### secp256k1 (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w32-il | 4.02 G | none | n/a | 52.74 M | 6.40 G | n/a | 0.63x |
| SUBTRACT | w32-il64 | 4.11 G | none | n/a | 31.29 M | 6.57 G | n/a | 0.63x |
| ADDMOD | w32-il64 | 4.40 G | none | n/a | 10.35 M | 5.16 G | n/a | 0.85x |
| SUBTRACTMOD | w32-il64 | 4.37 G | none | n/a | 12.04 M | 5.14 G | n/a | 0.85x |
| MULTIPLYOPERANDSCANNING | w32-il | 3.93 G | none | n/a | 20.13 M | n/a | n/a | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 2.51 G | none | n/a | 19.66 M | 6.40 G | n/a | 0.39x |
| MONTGOMERYMULTIPLICATION | w32-il64 | 3.87 G | none | n/a | 4.21 M | 4.82 G | n/a | 0.80x |
| COMPARE | w32-il | 4.02 G | none | n/a | 57.12 M | 6.40 G | n/a | 0.63x |
| REDUCE | w32-il64 | 1.44 G | none | n/a | 10.58 M | 3.67 G | n/a | 0.39x |
| MODMUL | w32-il64 | 591.78 M | none | n/a | 9.22 M | 1.60 G | n/a | 0.37x |
| MODEXP | w32-il64 | 42.20 M | none | n/a | 89.45 k | 5.69 M | n/a | 7.41x |
| EXPONENTIATION | w32-o64 | 195.46 M | none | n/a | 161.05 k | n/a | n/a | n/a |
| DIVIDE | w32-il64 | 491.96 M | none | n/a | 4.45 M | 3.18 G | n/a | 0.15x |
| ISQRT | w32-il64 | 83.46 M | none | n/a | 2.40 M | n/a | n/a | n/a |
| MODMUL_R2 | w32-il64 | 3.16 G | none | n/a | 5.32 M | 3.59 G | n/a | 0.88x |

### rsa256(composite) (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w32-il | 3.89 G | none | n/a | 18.97 M | 6.33 G | n/a | 0.61x |
| SUBTRACT | w32-il64 | 4.09 G | none | n/a | 32.44 M | 6.59 G | n/a | 0.62x |
| ADDMOD | w32-il | 4.23 G | none | n/a | 12.16 M | 5.28 G | n/a | 0.80x |
| SUBTRACTMOD | w32-il64 | 4.40 G | none | n/a | 12.07 M | 5.21 G | n/a | 0.85x |
| MULTIPLYOPERANDSCANNING | w32-il64 | 3.89 G | none | n/a | 21.12 M | n/a | n/a | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 2.53 G | none | n/a | 21.06 M | 6.18 G | n/a | 0.41x |
| MONTGOMERYMULTIPLICATION | w32-il64 | 3.79 G | none | n/a | 4.21 M | 4.96 G | n/a | 0.76x |
| COMPARE | w32-il64 | 4.17 G | none | n/a | 57.61 M | 6.38 G | n/a | 0.65x |
| REDUCE | w32-il64 | 1.43 G | none | n/a | 7.22 M | 3.78 G | n/a | 0.38x |
| MODMUL | w32-il64 | 592.75 M | none | n/a | 2.02 M | 1.58 G | n/a | 0.37x |
| MODEXP | w32-il64 | 42.31 M | none | n/a | 44.66 k | 5.82 M | n/a | 7.27x |
| EXPONENTIATION | w32-o64 | 195.61 M | none | n/a | 60.87 k | n/a | n/a | n/a |
| DIVIDE | w32-il64 | 454.18 M | none | n/a | 5.24 M | 3.18 G | n/a | 0.14x |
| ISQRT | w32-il64 | 74.67 M | none | n/a | 2.38 M | n/a | n/a | n/a |
| MODMUL_R2 | w32-il64 | 3.11 G | none | n/a | 5.72 M | 3.46 G | n/a | 0.90x |

### brainpoolP512r1 (512-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w32-il | 3.46 G | none | n/a | 17.38 M | 5.94 G | n/a | 0.58x |
| SUBTRACT | w32-il | 3.49 G | none | n/a | 24.84 M | 5.77 G | n/a | 0.61x |
| ADDMOD | w32-il | 3.61 G | none | n/a | 9.79 M | 4.93 G | n/a | 0.73x |
| SUBTRACTMOD | w32-il | 3.62 G | none | n/a | 9.89 M | 4.66 G | n/a | 0.78x |
| MULTIPLYOPERANDSCANNING | w32-il | 1.88 G | none | n/a | 7.83 M | n/a | n/a | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il | 1.09 G | none | n/a | 8.09 M | 4.84 G | n/a | 0.23x |
| MONTGOMERYMULTIPLICATION | w32-il64 | 2.35 G | none | n/a | 1.91 M | 3.65 G | n/a | 0.65x |
| COMPARE | w32-il64 | 3.65 G | none | n/a | 48.88 M | 5.90 G | n/a | 0.62x |
| REDUCE | w32-o64 | 554.32 M | none | n/a | 5.97 M | 2.86 G | n/a | 0.19x |
| MODMUL | w32-o64 | 183.25 M | none | n/a | 928.33 k | 629.28 M | n/a | 0.29x |
| MODEXP | w32-o64 | 6.08 M | none | n/a | 19.18 k | 2.54 M | n/a | 2.39x |
| EXPONENTIATION | w32-o64 | 5.96 M | none | n/a | 24.89 k | n/a | n/a | n/a |
| DIVIDE | w32-il64 | 148.40 M | none | n/a | 4.76 M | 2.09 G | n/a | 0.07x |
| ISQRT | w32-o64 | 16.65 M | none | n/a | 1.17 M | n/a | n/a | n/a |
| MODMUL_R2 | w32-il64 | 1.62 G | none | n/a | 2.64 M | 2.50 G | n/a | 0.65x |

### p1024 (1024-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w32-il | 1.84 G | none | n/a | 12.94 M | 4.87 G | n/a | 0.38x |
| SUBTRACT | w32-il64 | 1.84 G | none | n/a | 21.58 M | 5.04 G | n/a | 0.37x |
| ADDMOD | w32-il | 3.02 G | none | n/a | 5.86 M | 4.29 G | n/a | 0.70x |
| SUBTRACTMOD | w32-il64 | 2.91 G | none | n/a | 7.04 M | 4.28 G | n/a | 0.68x |
| MULTIPLYOPERANDSCANNING | w32-il64 | 816.45 M | none | n/a | 2.39 M | n/a | n/a | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 277.36 M | none | n/a | 2.22 M | 2.41 G | n/a | 0.12x |
| MONTGOMERYMULTIPLICATION | w32-il64 | 912.99 M | none | n/a | 752.63 k | 1.65 G | n/a | 0.55x |
| COMPARE | w32-il | 3.04 G | none | n/a | 34.82 M | 4.65 G | n/a | 0.65x |
| REDUCE | w32-il | 165.88 M | none | n/a | 9.14 M | 2.18 G | n/a | 0.08x |
| MODMUL | w32-o64 | 48.48 M | none | n/a | 299.39 k | 267.37 M | n/a | 0.18x |
| MODEXP | w32-o64 | 633.84 k | none | n/a | 4.07 k | 506.63 k | n/a | 1.25x |
| EXPONENTIATION | w32-o64 | 537.27 k | none | n/a | 11.50 k | n/a | n/a | n/a |
| DIVIDE | w32-opt | 33.22 M | none | n/a | 4.47 M | 1.66 G | n/a | 0.02x |
| ISQRT | w32-o64 | 2.28 M | none | n/a | 761.52 k | n/a | n/a | n/a |
| MODMUL_R2 | w32-il64 | 518.79 M | none | n/a | 1.09 M | 919.66 M | n/a | 0.56x |

### p2048 (2048-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w32-il | 1.22 G | none | n/a | 10.42 M | 3.68 G | n/a | 0.33x |
| SUBTRACT | w32-il | 1.22 G | none | n/a | 12.10 M | 4.01 G | n/a | 0.31x |
| ADDMOD | w32-il | 2.16 G | none | n/a | 4.59 M | 4.12 G | n/a | 0.52x |
| SUBTRACTMOD | w32-il | 2.14 G | none | n/a | 5.28 M | 3.77 G | n/a | 0.57x |
| MULTIPLYOPERANDSCANNING | w32-il64 | 267.20 M | none | n/a | 880.32 k | n/a | n/a | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 64.87 M | none | n/a | 956.26 k | 879.29 M | n/a | 0.07x |
| MONTGOMERYMULTIPLICATION | w32 | 255.11 M | none | n/a | 316.70 k | 576.78 M | n/a | 0.44x |
| COMPARE | w32-il64 | 2.06 G | none | n/a | 45.00 M | 4.05 G | n/a | 0.51x |
| REDUCE | w32-o64 | 62.96 M | none | n/a | 5.65 M | 2.05 G | n/a | 0.03x |
| MODMUL | w32-il64 | 14.25 M | none | n/a | 138.07 k | 133.64 M | n/a | 0.11x |
| MODEXP | w32-il64 | 50.52 k | none | n/a | 587.5 | 90.79 k | n/a | 0.56x |
| EXPONENTIATION | w32-o64 | 63.82 k | none | n/a | 4.59 k | n/a | n/a | n/a |
| DIVIDE | w32-il64 | 3.82 M | none | n/a | 2.88 M | 1.69 G | n/a | 0.00x |
| ISQRT | w32-o64 | 979.73 k | none | n/a | 361.86 k | n/a | n/a | n/a |
| MODMUL_R2 | w32-il64 | 131.21 M | none | n/a | 446.29 k | 306.37 M | n/a | 0.43x |

## 6. Raw data

Also written to `NVIDIA_B300_SXM6_AC_Report.csv` for analysis.

```csv
kind,device,device_type,kernel,modulus,bits,operation,items,seconds,ops_per_sec,mismatches
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,secp256k1,256,ADD,50000,0.000948051,52739777.665,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,secp256k1,256,ADD,50000,0.100014896,499925.531,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,secp256k1,256,ADD,50000,0.004245878,11776127.358,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,secp256k1,256,ADD,50000,0.000034816,1436121323.529,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,secp256k1,256,ADD,50000,0.000023070,2167320044.709,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,secp256k1,256,ADD,50000,0.001228570,40697721.517,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,secp256k1,256,ADD,50000,0.000017629,2836234941.096,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,secp256k1,256,ADD,50000,0.000972261,51426519.607,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32,secp256k1,256,ADD,50000,0.000014345,3485524046.647,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32,secp256k1,256,ADD,50000,0.000963612,51888104.667,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,secp256k1,256,ADD,50000,0.000013365,3741114064.344,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,secp256k1,256,ADD,50000,0.001006137,49695021.771,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,secp256k1,256,ADD,50000,0.000013415,3727170189.354,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,secp256k1,256,ADD,50000,0.000968440,51629426.232,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,secp256k1,256,ADD,50000,0.000012428,4023162441.455,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,secp256k1,256,ADD,50000,0.000975740,51243157.686,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,secp256k1,256,ADD,50000,0.000012754,3920339640.001,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,secp256k1,256,ADD,50000,0.001003447,49828242.608,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,secp256k1,256,SUBTRACT,50000,0.001598189,31285420.805,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,secp256k1,256,SUBTRACT,50000,0.004264250,11725391.362,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,secp256k1,256,SUBTRACT,50000,0.003503977,14269500.053,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,secp256k1,256,SUBTRACT,50000,0.000033760,1481042654.028,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,secp256k1,256,SUBTRACT,50000,0.000022949,2178748691.726,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,secp256k1,256,SUBTRACT,50000,0.000995782,50211794.799,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,secp256k1,256,SUBTRACT,50000,0.000017453,2864839445.037,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,secp256k1,256,SUBTRACT,50000,0.000983376,50845251.577,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32,secp256k1,256,SUBTRACT,50000,0.000013512,3700404761.023,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32,secp256k1,256,SUBTRACT,50000,0.000955803,52312036.622,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,secp256k1,256,SUBTRACT,50000,0.000013688,3652835366.860,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,secp256k1,256,SUBTRACT,50000,0.000947217,52786215.257,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,secp256k1,256,SUBTRACT,50000,0.000013581,3681610917.195,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,secp256k1,256,SUBTRACT,50000,0.000968900,51604913.570,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,secp256k1,256,SUBTRACT,50000,0.000012321,4058115391.193,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,secp256k1,256,SUBTRACT,50000,0.000951657,52539937.489,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,secp256k1,256,SUBTRACT,50000,0.000012155,4113539630.593,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,secp256k1,256,SUBTRACT,50000,0.000999809,50009551.326,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,secp256k1,256,ADDMOD,50000,0.004831800,10348110.425,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,secp256k1,256,ADDMOD,50000,0.004615775,10832417.048,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,secp256k1,256,ADDMOD,50000,0.007746756,6454314.587,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,secp256k1,256,ADDMOD,50000,0.000033248,1503849855.630,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,secp256k1,256,ADDMOD,50000,0.000029299,1706542842.101,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,secp256k1,256,ADDMOD,50000,0.000993203,50342174.627,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,secp256k1,256,ADDMOD,50000,0.000021278,2349846559.049,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,secp256k1,256,ADDMOD,50000,0.000973545,51358693.473,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32,secp256k1,256,ADDMOD,50000,0.000014723,3396049874.476,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32,secp256k1,256,ADDMOD,50000,0.000975546,51253348.376,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,secp256k1,256,ADDMOD,50000,0.000014040,3561244161.422,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,secp256k1,256,ADDMOD,50000,0.000962093,51970027.477,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,secp256k1,256,ADDMOD,50000,0.000013773,3630296339.247,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,secp256k1,256,ADDMOD,50000,0.000980986,50969127.438,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,secp256k1,256,ADDMOD,50000,0.000011796,4238719488.389,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,secp256k1,256,ADDMOD,50000,0.000971692,51456635.785,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,secp256k1,256,ADDMOD,50000,0.000011367,4398710892.397,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,secp256k1,256,ADDMOD,50000,0.001003160,49842496.465,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,secp256k1,256,SUBTRACTMOD,50000,0.004153164,12039014.148,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,secp256k1,256,SUBTRACTMOD,50000,0.007055336,7086834.733,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,secp256k1,256,SUBTRACTMOD,50000,0.005216695,9584612.445,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,secp256k1,256,SUBTRACTMOD,50000,0.000033216,1505298651.252,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.000029493,1695317137.636,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.001217692,41061286.957,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,secp256k1,256,SUBTRACTMOD,50000,0.000020924,2389597685.494,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,secp256k1,256,SUBTRACTMOD,50000,0.001010087,49500685.218,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32,secp256k1,256,SUBTRACTMOD,50000,0.000014664,3409719000.496,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32,secp256k1,256,SUBTRACTMOD,50000,0.000981878,50922825.201,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,secp256k1,256,SUBTRACTMOD,50000,0.000013892,3599204981.103,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,secp256k1,256,SUBTRACTMOD,50000,0.000975792,51240428.015,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,secp256k1,256,SUBTRACTMOD,50000,0.000013631,3668106274.260,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,secp256k1,256,SUBTRACTMOD,50000,0.000987908,50612000.922,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,secp256k1,256,SUBTRACTMOD,50000,0.000012259,4078636420.269,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,secp256k1,256,SUBTRACTMOD,50000,0.000981488,50943056.213,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,secp256k1,256,SUBTRACTMOD,50000,0.000011429,4374852095.259,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,secp256k1,256,SUBTRACTMOD,50000,0.000948291,52726432.084,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.002483709,20131182.704,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.006701134,7461423.726,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.004900290,10203477.773,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001710903,29224333.311,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.003073096,16270236.696,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000612184,81674790.870,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001897561,26349614.662,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000081762,611531231.188,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001320517,37863958.355,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000017754,2816270426.968,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001337166,37392515.145,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000016460,3037663277.236,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001114762,44852624.026,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000012716,3932058002.646,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001302816,38378405.705,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000013027,3838184909.876,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001255707,39818205.987,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.002543557,19657510.900,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.007106125,7036183.626,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.004891389,10222045.209,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000033632,1486679352.997,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000167168,299100240.849,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001473216,33939354.661,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000058630,852806349.137,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001256624,39789149.217,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000023307,2145274015.764,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001325993,37707590.170,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000023444,2132741075.317,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001300372,38450535.436,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000023176,2157407723.528,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001250705,39977451.392,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000020290,2464271970.486,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001299174,38485990.778,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000019885,2514455916.774,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001188542,42068349.453,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.011886971,4206286.021,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.005561863,8989793.599,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.005019724,9960706.964,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000033632,1486679352.997,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000073232,682761731.987,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.001218417,41036852.913,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000027776,1800115382.616,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000987282,50644090.385,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000014893,3357279212.069,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000968531,51624573.340,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000016572,3017145734.517,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000967812,51662927.994,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000014369,3479706791.759,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000956018,52300271.062,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000015962,3132440365.393,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000957988,52192722.225,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000012918,3870577793.799,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000965499,51786691.177,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,secp256k1,256,COMPARE,50000,0.000875332,57121215.081,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,secp256k1,256,COMPARE,50000,0.002204377,22682150.654,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,secp256k1,256,COMPARE,50000,0.002212924,22594540.273,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,secp256k1,256,COMPARE,50000,0.000034752,1438766114.180,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,secp256k1,256,COMPARE,50000,0.000023132,2161506420.906,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,secp256k1,256,COMPARE,50000,0.000997942,50103112.244,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,secp256k1,256,COMPARE,50000,0.000017758,2815633470.565,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,secp256k1,256,COMPARE,50000,0.001008634,49571993.030,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,secp256k1,256,COMPARE,50000,0.000012581,3974245670.399,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,secp256k1,256,COMPARE,50000,0.000958816,52147647.030,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,secp256k1,256,COMPARE,50000,0.000012795,3907766274.674,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,secp256k1,256,COMPARE,50000,0.000967825,51662235.100,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,secp256k1,256,COMPARE,50000,0.000012432,4021881539.470,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,secp256k1,256,COMPARE,50000,0.000956713,52262277.786,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,secp256k1,256,COMPARE,50000,0.000012614,3963865602.245,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,secp256k1,256,COMPARE,50000,0.000986211,50699089.652,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,secp256k1,256,REDUCE,6250,0.000590853,10577927.163,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,secp256k1,256,REDUCE,6250,0.001356320,4608055.399,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,secp256k1,256,REDUCE,6250,0.022458600,278289.825,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,secp256k1,256,REDUCE,50000,0.000034656,1442751615.882,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,secp256k1,256,REDUCE,50000,0.000134805,370906006.003,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,secp256k1,256,REDUCE,50000,0.001289013,38789369.555,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,secp256k1,256,REDUCE,50000,0.000092501,540534723.926,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,secp256k1,256,REDUCE,50000,0.001263018,39587717.704,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,secp256k1,256,REDUCE,50000,0.000036667,1363624537.174,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,secp256k1,256,REDUCE,50000,0.000960931,52032871.893,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,secp256k1,256,REDUCE,50000,0.000035765,1398014867.570,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,secp256k1,256,REDUCE,50000,0.000974273,51320316.864,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,secp256k1,256,REDUCE,50000,0.000035931,1391555766.657,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,secp256k1,256,REDUCE,50000,0.000982955,50867029.446,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,secp256k1,256,REDUCE,50000,0.000034780,1437609866.196,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,secp256k1,256,REDUCE,50000,0.000991399,50433781.332,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,secp256k1,256,MODMUL,3125,0.000338840,9222633.834,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,secp256k1,256,MODMUL,3125,0.000633431,4933451.286,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,secp256k1,256,MODMUL,3125,0.009411542,332039.111,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,secp256k1,256,MODMUL,50000,0.000047200,1059322033.898,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,secp256k1,256,MODMUL,50000,0.000362612,137888425.615,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,secp256k1,256,MODMUL,50000,0.001515444,32993630.498,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,secp256k1,256,MODMUL,50000,0.000216515,230930905.753,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,secp256k1,256,MODMUL,50000,0.001515070,33001776.009,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,secp256k1,256,MODMUL,50000,0.000104768,477244836.910,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,secp256k1,256,MODMUL,50000,0.001054334,47423302.807,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,secp256k1,256,MODMUL,50000,0.000084615,590912033.372,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,secp256k1,256,MODMUL,50000,0.001090430,45853470.048,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,secp256k1,256,MODMUL,50000,0.000102279,488858861.104,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,secp256k1,256,MODMUL,50000,0.001260767,39658398.660,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,secp256k1,256,MODMUL,50000,0.000084491,591779139.454,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,secp256k1,256,MODMUL,50000,0.001102583,45348059.257,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,secp256k1,256,MODEXP,781,0.008731146,89449.884,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,secp256k1,256,MODEXP,781,0.010534544,74137.048,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,secp256k1,256,MODEXP,781,0.010164318,76837.423,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,secp256k1,256,MODEXP,50000,0.076185152,656295.862,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,secp256k1,256,MODEXP,50000,0.016602785,3011542.947,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,secp256k1,256,MODEXP,50000,0.017947964,2785831.312,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,secp256k1,256,MODEXP,50000,0.003086865,16197663.291,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,secp256k1,256,MODEXP,50000,0.004414423,11326508.533,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,secp256k1,256,MODEXP,50000,0.002480722,20155422.703,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,secp256k1,256,MODEXP,50000,0.003648719,13703439.505,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,secp256k1,256,MODEXP,50000,0.001187896,42091226.603,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,secp256k1,256,MODEXP,50000,0.002328912,21469252.327,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,secp256k1,256,MODEXP,50000,0.002460207,20323493.309,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,secp256k1,256,MODEXP,50000,0.003623155,13800127.222,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,secp256k1,256,MODEXP,50000,0.001184839,42199826.355,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,secp256k1,256,MODEXP,50000,0.002326507,21491446.682,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,secp256k1,256,EXPONENTIATION,781,0.004849418,161050.254,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,secp256k1,256,EXPONENTIATION,781,0.003386065,230651.212,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,secp256k1,256,EXPONENTIATION,781,0.098306771,7944.519,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,secp256k1,256,EXPONENTIATION,50000,0.022657113,2206812.489,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,secp256k1,256,EXPONENTIATION,50000,0.024059744,2078160.103,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,secp256k1,256,EXPONENTIATION,50000,0.004565816,10950945.106,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,secp256k1,256,EXPONENTIATION,50000,0.005937645,8420847.015,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,secp256k1,256,EXPONENTIATION,50000,0.000451877,110649590.072,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,secp256k1,256,EXPONENTIATION,50000,0.001600425,31241701.581,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,secp256k1,256,EXPONENTIATION,50000,0.000255805,195461361.546,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,secp256k1,256,EXPONENTIATION,50000,0.001393686,35876087.111,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,secp256k1,256,EXPONENTIATION,50000,0.000455800,109697240.870,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,secp256k1,256,EXPONENTIATION,50000,0.001584861,31548507.640,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,secp256k1,256,EXPONENTIATION,50000,0.000260048,192272179.615,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,secp256k1,256,EXPONENTIATION,50000,0.001240762,40297818.591,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,secp256k1,256,DIVIDE,6250,0.001404085,4451297.579,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,secp256k1,256,DIVIDE,6250,0.098207776,63640.582,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,secp256k1,256,DIVIDE,6250,0.006382303,979270.332,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,secp256k1,256,DIVIDE,50000,0.000032928,1518464528.669,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,secp256k1,256,DIVIDE,50000,0.000328905,152019603.635,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,secp256k1,256,DIVIDE,50000,0.001631468,30647245.104,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,secp256k1,256,DIVIDE,50000,0.000280199,178444611.576,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,secp256k1,256,DIVIDE,50000,0.001722309,29030794.614,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,secp256k1,256,DIVIDE,50000,0.000107279,466074488.494,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,secp256k1,256,DIVIDE,50000,0.001396412,35806050.743,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,secp256k1,256,DIVIDE,50000,0.000104898,476653489.252,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,secp256k1,256,DIVIDE,50000,0.001332052,37536070.699,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,secp256k1,256,DIVIDE,50000,0.000105450,474158422.669,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,secp256k1,256,DIVIDE,50000,0.001382448,36167726.866,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,secp256k1,256,DIVIDE,50000,0.000101634,491961302.109,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,secp256k1,256,DIVIDE,50000,0.001392953,35894964.567,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,secp256k1,256,ISQRT,1562,0.000651248,2398471.095,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,secp256k1,256,ISQRT,1562,0.001075869,1451849.879,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,secp256k1,256,ISQRT,50000,0.003584090,13950542.614,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,secp256k1,256,ISQRT,50000,0.004730263,10570236.863,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,secp256k1,256,ISQRT,50000,0.002671525,18715901.975,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,secp256k1,256,ISQRT,50000,0.004013701,12457330.551,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,secp256k1,256,ISQRT,50000,0.000658389,75942944.860,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,secp256k1,256,ISQRT,50000,0.001796520,27831585.137,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,secp256k1,256,ISQRT,50000,0.000606069,82498858.472,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,secp256k1,256,ISQRT,50000,0.001743268,28681763.402,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,secp256k1,256,ISQRT,50000,0.000673450,74244560.957,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,secp256k1,256,ISQRT,50000,0.001814400,27557318.912,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,secp256k1,256,ISQRT,50000,0.000599105,83457828.865,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,secp256k1,256,ISQRT,50000,0.001722129,29033829.577,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,secp256k1,256,MODMUL_R2,50000,0.009401991,5318022.516,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,secp256k1,256,MODMUL_R2,50000,0.005186232,9640910.815,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,secp256k1,256,MODMUL_R2,50000,0.008598911,5814689.794,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,secp256k1,256,MODMUL_R2,50000,0.000036864,1356336805.556,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,secp256k1,256,MODMUL_R2,50000,0.000087188,573473197.234,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,secp256k1,256,MODMUL_R2,50000,0.001246390,40115854.873,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,secp256k1,256,MODMUL_R2,50000,0.000029026,1722593132.436,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,secp256k1,256,MODMUL_R2,50000,0.001323637,37774706.295,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,secp256k1,256,MODMUL_R2,50000,0.000020809,2402805791.393,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,secp256k1,256,MODMUL_R2,50000,0.000947478,52771674.161,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,secp256k1,256,MODMUL_R2,50000,0.000016303,3066913708.745,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,secp256k1,256,MODMUL_R2,50000,0.000945951,52856858.611,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,secp256k1,256,MODMUL_R2,50000,0.000020375,2453993581.287,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,secp256k1,256,MODMUL_R2,50000,0.000962580,51943736.052,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,secp256k1,256,MODMUL_R2,50000,0.000015798,3164952614.514,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,secp256k1,256,MODMUL_R2,50000,0.000965730,51774305.657,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,rsa256(composite),256,ADD,50000,0.002635504,18971703.306,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,rsa256(composite),256,ADD,50000,0.004414158,11327188.560,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,rsa256(composite),256,ADD,50000,0.004512190,11081093.650,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,rsa256(composite),256,ADD,50000,0.000032832,1522904483.431,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,rsa256(composite),256,ADD,50000,0.000022578,2214545625.917,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,rsa256(composite),256,ADD,50000,0.000981338,50950843.047,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,rsa256(composite),256,ADD,50000,0.000017882,2796102545.807,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,rsa256(composite),256,ADD,50000,0.001127636,44340549.535,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32,rsa256(composite),256,ADD,50000,0.000013454,3716366236.621,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32,rsa256(composite),256,ADD,50000,0.000983254,50851560.512,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,rsa256(composite),256,ADD,50000,0.000013813,3619786600.367,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,rsa256(composite),256,ADD,50000,0.000955247,52342484.738,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,rsa256(composite),256,ADD,50000,0.000013130,3808067753.090,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,rsa256(composite),256,ADD,50000,0.000984976,50762657.665,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,rsa256(composite),256,ADD,50000,0.000012854,3889840416.610,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,rsa256(composite),256,ADD,50000,0.000976876,51183568.507,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,rsa256(composite),256,ADD,50000,0.000012933,3866083340.609,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,rsa256(composite),256,ADD,50000,0.000973476,51362335.587,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACT,50000,0.001541091,32444557.549,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACT,50000,0.002257106,22152265.211,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACT,50000,0.002144202,23318698.432,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,rsa256(composite),256,SUBTRACT,50000,0.000032896,1519941634.241,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,rsa256(composite),256,SUBTRACT,50000,0.000022228,2249410040.458,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,rsa256(composite),256,SUBTRACT,50000,0.001002081,49896164.605,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,rsa256(composite),256,SUBTRACT,50000,0.000017523,2853391240.454,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,rsa256(composite),256,SUBTRACT,50000,0.001124618,44459541.209,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32,rsa256(composite),256,SUBTRACT,50000,0.000014027,3564554444.730,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32,rsa256(composite),256,SUBTRACT,50000,0.001054433,47418849.760,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,rsa256(composite),256,SUBTRACT,50000,0.000013581,3681610917.195,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,rsa256(composite),256,SUBTRACT,50000,0.000965976,51761121.341,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,rsa256(composite),256,SUBTRACT,50000,0.000013938,3587315503.251,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,rsa256(composite),256,SUBTRACT,50000,0.000976541,51201129.175,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,rsa256(composite),256,SUBTRACT,50000,0.000012347,4049563733.736,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,rsa256(composite),256,SUBTRACT,50000,0.000993231,50340755.547,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,rsa256(composite),256,SUBTRACT,50000,0.000012219,4091984409.373,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,rsa256(composite),256,SUBTRACT,50000,0.000976559,51200183.106,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,rsa256(composite),256,ADDMOD,50000,0.004111469,12161103.588,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,rsa256(composite),256,ADDMOD,50000,0.007104907,7037389.783,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,rsa256(composite),256,ADDMOD,50000,0.004170797,11988116.563,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,rsa256(composite),256,ADDMOD,50000,0.000033472,1493785850.860,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,rsa256(composite),256,ADDMOD,50000,0.000027685,1806030926.045,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,rsa256(composite),256,ADDMOD,50000,0.000993944,50304644.624,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,rsa256(composite),256,ADDMOD,50000,0.000019925,2509409507.784,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,rsa256(composite),256,ADDMOD,50000,0.001164236,42946619.455,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32,rsa256(composite),256,ADDMOD,50000,0.000014002,3570926282.883,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32,rsa256(composite),256,ADDMOD,50000,0.000954014,52410131.450,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,rsa256(composite),256,ADDMOD,50000,0.000013102,3816205264.537,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,rsa256(composite),256,ADDMOD,50000,0.000952765,52478838.111,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,rsa256(composite),256,ADDMOD,50000,0.000013470,3711949886.999,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,rsa256(composite),256,ADDMOD,50000,0.001009505,49529224.264,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,rsa256(composite),256,ADDMOD,50000,0.000011807,4234770039.883,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,rsa256(composite),256,ADDMOD,50000,0.001013731,49322749.238,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,rsa256(composite),256,ADDMOD,50000,0.000011883,4207699606.168,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,rsa256(composite),256,ADDMOD,50000,0.000965506,51786316.528,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACTMOD,50000,0.004142191,12070906.523,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACTMOD,50000,0.004198897,11907889.173,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACTMOD,50000,0.004264895,11723618.079,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,rsa256(composite),256,SUBTRACTMOD,50000,0.000032960,1516990291.262,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,rsa256(composite),256,SUBTRACTMOD,50000,0.000030670,1630259399.589,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,rsa256(composite),256,SUBTRACTMOD,50000,0.001044797,47856188.507,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,rsa256(composite),256,SUBTRACTMOD,50000,0.000021371,2339619010.162,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,rsa256(composite),256,SUBTRACTMOD,50000,0.001177245,42472042.983,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32,rsa256(composite),256,SUBTRACTMOD,50000,0.000014902,3355259709.235,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32,rsa256(composite),256,SUBTRACTMOD,50000,0.000954382,52389923.249,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,rsa256(composite),256,SUBTRACTMOD,50000,0.000013629,3668638916.910,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,rsa256(composite),256,SUBTRACTMOD,50000,0.001001549,49922669.208,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,rsa256(composite),256,SUBTRACTMOD,50000,0.000013715,3645641999.474,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,rsa256(composite),256,SUBTRACTMOD,50000,0.000976198,51219117.243,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,rsa256(composite),256,SUBTRACTMOD,50000,0.000012223,4090639836.183,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,rsa256(composite),256,SUBTRACTMOD,50000,0.000977411,51155553.144,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,rsa256(composite),256,SUBTRACTMOD,50000,0.000011358,4402182449.559,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,rsa256(composite),256,SUBTRACTMOD,50000,0.000942179,53068472.022,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.002367109,21122812.661,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.004976052,10048126.555,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.006997268,7145646.039,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001683476,29700453.577,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.002986964,16739404.888,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000617625,80955269.960,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001798004,27808614.903,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000081038,616994602.298,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001355420,36888933.669,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000017601,2840746534.031,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001188011,42087153.565,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000016737,2987398089.303,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001267253,39455420.075,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000013987,3574745455.378,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001283688,38950274.368,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000012870,3885002416.047,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001156636,43228812.208,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.002374271,21059095.710,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.004260659,11735273.697,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.004353060,11486172.845,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000032768,1525878906.250,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000167386,298710721.049,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001482468,33727539.916,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000057845,864379407.287,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001491879,33514782.441,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000024319,2056006766.923,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001138415,43920716.072,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000023803,2100572610.183,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001205520,41475877.427,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000022864,2186846892.057,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001306851,38259907.418,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000020467,2442959613.219,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001150130,43473346.369,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000019790,2526525660.166,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001315040,38021657.425,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.011884761,4207068.200,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.010811501,4624704.744,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.008365711,5976778.288,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000033184,1506750241.080,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000073268,682426516.802,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.001221518,40932674.492,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000027315,1830496558.896,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.001331530,37550787.792,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000015211,3287094741.756,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.001001215,49939325.827,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000017081,2927232531.718,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000961785,51986668.991,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000014948,3344924978.875,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000966550,51730379.882,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000015858,3152986951.894,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000982776,50876293.628,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000013201,3787582715.440,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000987508,50632501.976,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,rsa256(composite),256,COMPARE,50000,0.000867838,57614440.770,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,rsa256(composite),256,COMPARE,50000,0.002588451,19316571.988,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,rsa256(composite),256,COMPARE,50000,0.003897622,12828336.461,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,rsa256(composite),256,COMPARE,50000,0.000033344,1499520153.551,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,rsa256(composite),256,COMPARE,50000,0.000022935,2180075780.925,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,rsa256(composite),256,COMPARE,50000,0.000988687,50572120.408,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,rsa256(composite),256,COMPARE,50000,0.000017489,2858937356.511,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,rsa256(composite),256,COMPARE,50000,0.001324829,37740719.598,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,rsa256(composite),256,COMPARE,50000,0.000013078,3823220160.407,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,rsa256(composite),256,COMPARE,50000,0.000964100,51861842.814,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,rsa256(composite),256,COMPARE,50000,0.000012950,3861009170.304,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,rsa256(composite),256,COMPARE,50000,0.000981827,50925469.812,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,rsa256(composite),256,COMPARE,50000,0.000012438,4019942902.337,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,rsa256(composite),256,COMPARE,50000,0.001007157,49644694.389,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,rsa256(composite),256,COMPARE,50000,0.000011985,4171876091.908,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,rsa256(composite),256,COMPARE,50000,0.000974757,51294837.819,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,rsa256(composite),256,REDUCE,6250,0.000866093,7216315.200,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,rsa256(composite),256,REDUCE,6250,0.001922451,3251058.149,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,rsa256(composite),256,REDUCE,6250,0.001190493,5249925.872,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,rsa256(composite),256,REDUCE,50000,0.000033120,1509661835.749,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,rsa256(composite),256,REDUCE,50000,0.000135051,370230587.400,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,rsa256(composite),256,REDUCE,50000,0.001298631,38502083.597,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,rsa256(composite),256,REDUCE,50000,0.000093343,535658983.550,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,rsa256(composite),256,REDUCE,50000,0.001408989,35486438.467,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,rsa256(composite),256,REDUCE,50000,0.000036688,1362843524.381,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,rsa256(composite),256,REDUCE,50000,0.000972669,51404949.203,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,rsa256(composite),256,REDUCE,50000,0.000036325,1376461935.811,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,rsa256(composite),256,REDUCE,50000,0.000976860,51184407.216,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,rsa256(composite),256,REDUCE,50000,0.000036349,1375553803.459,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,rsa256(composite),256,REDUCE,50000,0.001003718,49814788.436,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,rsa256(composite),256,REDUCE,50000,0.000035007,1428285012.703,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,rsa256(composite),256,REDUCE,50000,0.000967843,51661271.914,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,rsa256(composite),256,MODMUL,3125,0.001545463,2022047.775,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,rsa256(composite),256,MODMUL,3125,0.004844906,645007.418,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,rsa256(composite),256,MODMUL,3125,0.003000507,1041490.650,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,rsa256(composite),256,MODMUL,50000,0.000038912,1284950657.895,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,rsa256(composite),256,MODMUL,50000,0.000363392,137592462.631,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,rsa256(composite),256,MODMUL,50000,0.001517961,32938921.924,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,rsa256(composite),256,MODMUL,50000,0.000217250,230149571.353,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,rsa256(composite),256,MODMUL,50000,0.001516719,32965895.474,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,rsa256(composite),256,MODMUL,50000,0.000104631,477869818.814,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,rsa256(composite),256,MODMUL,50000,0.001247092,40093273.777,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,rsa256(composite),256,MODMUL,50000,0.000084986,588332341.950,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,rsa256(composite),256,MODMUL,50000,0.001234287,40509217.436,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,rsa256(composite),256,MODMUL,50000,0.000102768,486532716.714,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,rsa256(composite),256,MODMUL,50000,0.001240140,40318028.553,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,rsa256(composite),256,MODMUL,50000,0.000084352,592754306.093,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,rsa256(composite),256,MODMUL,50000,0.001234006,40518440.899,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,rsa256(composite),256,MODEXP,781,0.017487557,44660.326,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,rsa256(composite),256,MODEXP,781,0.004725928,165258.548,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,rsa256(composite),256,MODEXP,781,0.007909102,98746.988,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,rsa256(composite),256,MODEXP,50000,0.073687039,678545.382,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,rsa256(composite),256,MODEXP,50000,0.016610807,3010088.545,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,rsa256(composite),256,MODEXP,50000,0.018024860,2773946.648,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,rsa256(composite),256,MODEXP,50000,0.003081325,16226785.682,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,rsa256(composite),256,MODEXP,50000,0.004378700,11418914.319,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,rsa256(composite),256,MODEXP,50000,0.002475838,20195182.380,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,rsa256(composite),256,MODEXP,50000,0.003639715,13737339.449,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,rsa256(composite),256,MODEXP,50000,0.001186152,42153113.995,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,rsa256(composite),256,MODEXP,50000,0.002333620,21425938.738,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,rsa256(composite),256,MODEXP,50000,0.002455689,20360884.353,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,rsa256(composite),256,MODEXP,50000,0.003593832,13912725.943,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,rsa256(composite),256,MODEXP,50000,0.001181876,42305622.736,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,rsa256(composite),256,MODEXP,50000,0.002339372,21373257.403,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,rsa256(composite),256,EXPONENTIATION,781,0.012829589,60874.904,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,rsa256(composite),256,EXPONENTIATION,781,0.006844499,114106.234,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,rsa256(composite),256,EXPONENTIATION,781,0.006985734,111799.276,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,rsa256(composite),256,EXPONENTIATION,50000,0.022650227,2207483.394,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,rsa256(composite),256,EXPONENTIATION,50000,0.024063443,2077840.647,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,rsa256(composite),256,EXPONENTIATION,50000,0.004577087,10923978.491,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,rsa256(composite),256,EXPONENTIATION,50000,0.005901754,8472057.685,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,rsa256(composite),256,EXPONENTIATION,50000,0.000452022,110614082.709,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,rsa256(composite),256,EXPONENTIATION,50000,0.001576825,31709289.501,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,rsa256(composite),256,EXPONENTIATION,50000,0.000255614,195607442.977,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,rsa256(composite),256,EXPONENTIATION,50000,0.001374788,36369243.833,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,rsa256(composite),256,EXPONENTIATION,50000,0.000455569,109752856.049,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,rsa256(composite),256,EXPONENTIATION,50000,0.001604701,31158452.865,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,rsa256(composite),256,EXPONENTIATION,50000,0.000259962,192335809.506,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,rsa256(composite),256,EXPONENTIATION,50000,0.001406717,35543751.928,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,rsa256(composite),256,DIVIDE,6250,0.001191843,5243977.787,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,rsa256(composite),256,DIVIDE,6250,0.033496909,186584.378,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,rsa256(composite),256,DIVIDE,6250,0.022365125,279452.947,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,rsa256(composite),256,DIVIDE,50000,0.000033728,1482447817.837,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,rsa256(composite),256,DIVIDE,50000,0.000342742,145882318.640,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,rsa256(composite),256,DIVIDE,50000,0.001664538,30038365.096,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,rsa256(composite),256,DIVIDE,50000,0.000290671,172015810.196,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,rsa256(composite),256,DIVIDE,50000,0.001720341,29064004.677,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,rsa256(composite),256,DIVIDE,50000,0.000115362,433418298.364,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,rsa256(composite),256,DIVIDE,50000,0.001395597,35826961.349,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,rsa256(composite),256,DIVIDE,50000,0.000112642,443884123.761,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,rsa256(composite),256,DIVIDE,50000,0.001364933,36631835.277,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,rsa256(composite),256,DIVIDE,50000,0.000113992,438627417.231,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,rsa256(composite),256,DIVIDE,50000,0.001388683,36005337.851,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,rsa256(composite),256,DIVIDE,50000,0.000110088,454182214.216,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,rsa256(composite),256,DIVIDE,50000,0.001400988,35689100.250,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,rsa256(composite),256,ISQRT,1562,0.000656877,2377917.336,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,rsa256(composite),256,ISQRT,1562,0.001166341,1339230.978,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,rsa256(composite),256,ISQRT,50000,0.003584031,13950772.128,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,rsa256(composite),256,ISQRT,50000,0.004747228,10532462.316,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,rsa256(composite),256,ISQRT,50000,0.002758739,18124222.748,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,rsa256(composite),256,ISQRT,50000,0.003829177,13057636.133,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,rsa256(composite),256,ISQRT,50000,0.000709161,70505852.498,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,rsa256(composite),256,ISQRT,50000,0.001840649,27164331.740,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,rsa256(composite),256,ISQRT,50000,0.000676276,73934312.194,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,rsa256(composite),256,ISQRT,50000,0.001803974,27716586.012,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,rsa256(composite),256,ISQRT,50000,0.000721142,69334467.920,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,rsa256(composite),256,ISQRT,50000,0.001869244,26748781.355,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,rsa256(composite),256,ISQRT,50000,0.000669627,74668436.516,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,rsa256(composite),256,ISQRT,50000,0.001821863,27444434.635,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,rsa256(composite),256,MODMUL_R2,50000,0.008748780,5715082.557,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,rsa256(composite),256,MODMUL_R2,50000,0.004477852,11166067.956,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,rsa256(composite),256,MODMUL_R2,50000,0.004205976,11887847.274,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,rsa256(composite),256,MODMUL_R2,50000,0.000036864,1356336805.556,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,rsa256(composite),256,MODMUL_R2,50000,0.000087633,570561498.695,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,rsa256(composite),256,MODMUL_R2,50000,0.001251777,39943216.923,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,rsa256(composite),256,MODMUL_R2,50000,0.000028493,1754817017.973,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,rsa256(composite),256,MODMUL_R2,50000,0.000979241,51059954.425,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,rsa256(composite),256,MODMUL_R2,50000,0.000021311,2346207416.148,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,rsa256(composite),256,MODMUL_R2,50000,0.000993607,50321707.497,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,rsa256(composite),256,MODMUL_R2,50000,0.000016976,2945339228.171,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,rsa256(composite),256,MODMUL_R2,50000,0.000966565,51729579.259,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,rsa256(composite),256,MODMUL_R2,50000,0.000020763,2408134058.115,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,rsa256(composite),256,MODMUL_R2,50000,0.000995154,50243478.546,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,rsa256(composite),256,MODMUL_R2,50000,0.000016087,3108094709.686,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,rsa256(composite),256,MODMUL_R2,50000,0.000953305,52449111.866,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,brainpoolP512r1,512,ADD,25000,0.001438384,17380615.994,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,brainpoolP512r1,512,ADD,25000,0.003495658,7151728.380,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,brainpoolP512r1,512,ADD,25000,0.002472666,10110542.608,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,brainpoolP512r1,512,ADD,50000,0.000040960,1220703125.000,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,brainpoolP512r1,512,ADD,50000,0.000056846,879569222.010,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,brainpoolP512r1,512,ADD,50000,0.002065428,24208057.501,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,brainpoolP512r1,512,ADD,50000,0.000032711,1528536987.140,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,brainpoolP512r1,512,ADD,50000,0.002014843,24815829.442,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32,brainpoolP512r1,512,ADD,50000,0.000019224,2600913988.119,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32,brainpoolP512r1,512,ADD,50000,0.001995676,25054166.827,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,brainpoolP512r1,512,ADD,50000,0.000018354,2724204564.901,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,brainpoolP512r1,512,ADD,50000,0.001974698,25320326.954,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,brainpoolP512r1,512,ADD,50000,0.000018674,2677518520.527,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,brainpoolP512r1,512,ADD,50000,0.002072394,24126686.211,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,brainpoolP512r1,512,ADD,50000,0.000014436,3463557609.603,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,brainpoolP512r1,512,ADD,50000,0.002009359,24883557.223,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,brainpoolP512r1,512,ADD,50000,0.000018539,2697013667.904,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,brainpoolP512r1,512,ADD,50000,0.002040754,24500748.782,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACT,25000,0.001006493,24838734.494,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACT,25000,0.002445040,10224783.631,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACT,25000,0.002996073,8344255.998,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,brainpoolP512r1,512,SUBTRACT,50000,0.000041344,1209365325.077,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,brainpoolP512r1,512,SUBTRACT,50000,0.000055027,908644521.805,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,brainpoolP512r1,512,SUBTRACT,50000,0.002094201,23875453.657,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,brainpoolP512r1,512,SUBTRACT,50000,0.000032551,1536050905.007,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,brainpoolP512r1,512,SUBTRACT,50000,0.002096792,23845950.853,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32,brainpoolP512r1,512,SUBTRACT,50000,0.000019310,2589334090.546,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32,brainpoolP512r1,512,SUBTRACT,50000,0.001994228,25072358.947,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,brainpoolP512r1,512,SUBTRACT,50000,0.000019081,2620408406.115,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,brainpoolP512r1,512,SUBTRACT,50000,0.002115394,23636258.903,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,brainpoolP512r1,512,SUBTRACT,50000,0.000019057,2623713970.849,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,brainpoolP512r1,512,SUBTRACT,50000,0.001975303,25312572.403,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,brainpoolP512r1,512,SUBTRACT,50000,0.000014326,3490155003.068,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,brainpoolP512r1,512,SUBTRACT,50000,0.002065115,24211726.401,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,brainpoolP512r1,512,SUBTRACT,50000,0.000018178,2750574644.568,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,brainpoolP512r1,512,SUBTRACT,50000,0.002078584,24054837.472,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,brainpoolP512r1,512,ADDMOD,25000,0.002553141,9791860.387,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,brainpoolP512r1,512,ADDMOD,25000,0.004311248,5798784.919,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,brainpoolP512r1,512,ADDMOD,25000,0.004280404,5840570.171,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,brainpoolP512r1,512,ADDMOD,50000,0.000039872,1254012841.091,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,brainpoolP512r1,512,ADDMOD,50000,0.000070210,712149629.829,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,brainpoolP512r1,512,ADDMOD,50000,0.002088164,23944479.693,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,brainpoolP512r1,512,ADDMOD,50000,0.000037709,1325942768.231,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,brainpoolP512r1,512,ADDMOD,50000,0.002132508,23446571.195,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32,brainpoolP512r1,512,ADDMOD,50000,0.000020724,2412665774.994,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32,brainpoolP512r1,512,ADDMOD,50000,0.002052606,24359278.304,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,brainpoolP512r1,512,ADDMOD,50000,0.000020204,2474758022.714,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,brainpoolP512r1,512,ADDMOD,50000,0.002103587,23768923.879,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,brainpoolP512r1,512,ADDMOD,50000,0.000020139,2482747689.062,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,brainpoolP512r1,512,ADDMOD,50000,0.002068478,24172362.787,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,brainpoolP512r1,512,ADDMOD,50000,0.000013832,3614805494.210,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,brainpoolP512r1,512,ADDMOD,50000,0.002007686,24904292.689,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,brainpoolP512r1,512,ADDMOD,50000,0.000014195,3522370541.440,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,brainpoolP512r1,512,ADDMOD,50000,0.001980797,25242364.074,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACTMOD,25000,0.002527615,9890746.822,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACTMOD,25000,0.002463147,10149617.600,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACTMOD,25000,0.002428858,10292901.039,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000040960,1220703125.000,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000074801,668440466.869,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,brainpoolP512r1,512,SUBTRACTMOD,50000,0.002226764,22454107.838,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000039954,1251440054.545,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,brainpoolP512r1,512,SUBTRACTMOD,50000,0.002123605,23544868.342,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000021425,2333720367.638,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32,brainpoolP512r1,512,SUBTRACTMOD,50000,0.002008689,24891856.786,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000020215,2473411229.200,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,brainpoolP512r1,512,SUBTRACTMOD,50000,0.002058029,24295090.106,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000020216,2473290160.895,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,brainpoolP512r1,512,SUBTRACTMOD,50000,0.002050609,24383000.622,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000013812,3620045932.369,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,brainpoolP512r1,512,SUBTRACTMOD,50000,0.002008914,24889069.938,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000014416,3468368397.634,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,brainpoolP512r1,512,SUBTRACTMOD,50000,0.002047802,24416422.708,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.003193058,7829485.052,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.004270038,5854748.774,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.005075007,4926101.601,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.004721592,10589648.563,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.007218023,6927104.550,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.001933578,25858796.566,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.004321071,11571205.460,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000873376,57249109.313,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.003240393,15430227.370,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000036552,1367913286.058,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.002374779,21054591.053,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000034985,1429183278.262,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.002394516,20881046.849,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000026529,1884733017.965,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.002267071,22054889.446,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000028981,1725271012.615,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.002262383,22100590.386,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.003091270,8087290.864,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.007240683,3452712.976,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.004291896,5824931.439,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000040928,1221657544.957,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.001107995,45126557.801,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.003500406,14284057.138,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000308863,161884100.128,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.002727628,18330945.118,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000090348,553415617.778,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.002500628,19994977.223,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000090463,552711982.874,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.002446237,20439557.074,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000088546,564678133.774,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.002348876,21286777.381,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000045756,1092752994.217,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.002265816,22067104.887,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000058356,856809451.907,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.002458378,20338613.723,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.013057127,1914663.155,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.004495396,5561245.343,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.004501586,5553598.264,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000040928,1221657544.957,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000243875,205023045.386,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.002363255,21157259.665,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000074978,666862399.543,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.002432778,20552635.610,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000026212,1907522132.533,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.002069487,24160576.951,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000036600,1366120679.471,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.001998327,25020930.334,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000026684,1873780526.495,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.002101727,23789959.348,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000033116,1509845638.112,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.002025938,24679925.908,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000021234,2354716346.900,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.001984220,25198818.703,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,brainpoolP512r1,512,COMPARE,25000,0.000511460,48879708.712,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,brainpoolP512r1,512,COMPARE,25000,0.001479623,16896199.565,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,brainpoolP512r1,512,COMPARE,25000,0.001524377,16400138.798,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,brainpoolP512r1,512,COMPARE,50000,0.000040960,1220703125.000,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,brainpoolP512r1,512,COMPARE,50000,0.000040685,1228954120.891,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,brainpoolP512r1,512,COMPARE,50000,0.002042491,24479912.161,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,brainpoolP512r1,512,COMPARE,50000,0.000024555,2036243657.017,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,brainpoolP512r1,512,COMPARE,50000,0.002314741,21600688.572,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,brainpoolP512r1,512,COMPARE,50000,0.000016651,3002822672.008,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,brainpoolP512r1,512,COMPARE,50000,0.001972408,25349724.368,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,brainpoolP512r1,512,COMPARE,50000,0.000016438,3041729228.089,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,brainpoolP512r1,512,COMPARE,50000,0.002053963,24343184.129,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,brainpoolP512r1,512,COMPARE,50000,0.000014462,3457340199.231,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,brainpoolP512r1,512,COMPARE,50000,0.002037994,24533928.474,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,brainpoolP512r1,512,COMPARE,50000,0.000013682,3654436022.208,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,brainpoolP512r1,512,COMPARE,50000,0.001963690,25462266.951,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,brainpoolP512r1,512,REDUCE,3125,0.000523642,5967816.136,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,brainpoolP512r1,512,REDUCE,3125,0.000738527,4231395.725,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,brainpoolP512r1,512,REDUCE,3125,0.000746929,4183797.113,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,brainpoolP512r1,512,REDUCE,50000,0.000040864,1223570869.225,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,brainpoolP512r1,512,REDUCE,50000,0.000385554,129683525.885,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,brainpoolP512r1,512,REDUCE,50000,0.002498870,20009043.952,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,brainpoolP512r1,512,REDUCE,50000,0.000317189,157634708.821,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,brainpoolP512r1,512,REDUCE,50000,0.002647347,18886832.802,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,brainpoolP512r1,512,REDUCE,50000,0.000092162,542522984.265,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,brainpoolP512r1,512,REDUCE,50000,0.002046628,24430428.498,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,brainpoolP512r1,512,REDUCE,50000,0.000090200,554323796.740,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,brainpoolP512r1,512,REDUCE,50000,0.002080879,24028307.273,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,brainpoolP512r1,512,REDUCE,50000,0.000090931,549867370.980,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,brainpoolP512r1,512,REDUCE,50000,0.002082562,24008888.868,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,brainpoolP512r1,512,REDUCE,50000,0.000101218,493983283.340,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,brainpoolP512r1,512,REDUCE,50000,0.002061352,24255925.033,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL,1562,0.001682587,928332.390,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL,1562,0.002327962,670973.157,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL,1562,0.036503608,42790.291,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,brainpoolP512r1,512,MODMUL,50000,0.000113568,440264863.342,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,brainpoolP512r1,512,MODMUL,50000,0.001126403,44389086.754,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,brainpoolP512r1,512,MODMUL,50000,0.003233174,15464679.424,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,brainpoolP512r1,512,MODMUL,50000,0.000794752,62912709.174,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,brainpoolP512r1,512,MODMUL,50000,0.003110080,16076756.733,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,brainpoolP512r1,512,MODMUL,50000,0.000320710,155904083.448,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,brainpoolP512r1,512,MODMUL,50000,0.002520159,19840018.004,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,brainpoolP512r1,512,MODMUL,50000,0.000272845,183254205.428,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,brainpoolP512r1,512,MODMUL,50000,0.002391872,20904128.682,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,brainpoolP512r1,512,MODMUL,50000,0.000317945,157259889.269,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,brainpoolP512r1,512,MODMUL,50000,0.002556268,19559764.398,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,brainpoolP512r1,512,MODMUL,50000,0.000335678,148952273.039,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,brainpoolP512r1,512,MODMUL,50000,0.002407437,20768975.334,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,brainpoolP512r1,512,MODEXP,390,0.020333707,19179.975,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,brainpoolP512r1,512,MODEXP,390,0.008758765,44526.825,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,brainpoolP512r1,512,MODEXP,390,0.103374571,3772.688,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,brainpoolP512r1,512,MODEXP,50000,0.674238443,74157.741,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,brainpoolP512r1,512,MODEXP,50000,0.230710321,216721.991,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,brainpoolP512r1,512,MODEXP,50000,0.233135552,214467.504,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,brainpoolP512r1,512,MODEXP,50000,0.021595255,2315323.441,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,brainpoolP512r1,512,MODEXP,50000,0.024147092,2070642.708,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,brainpoolP512r1,512,MODEXP,50000,0.019098882,2617954.289,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,brainpoolP512r1,512,MODEXP,50000,0.021554249,2319728.237,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,brainpoolP512r1,512,MODEXP,50000,0.008219358,6083200.138,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,brainpoolP512r1,512,MODEXP,50000,0.010411561,4802353.871,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,brainpoolP512r1,512,MODEXP,50000,0.019197142,2604554.368,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,brainpoolP512r1,512,MODEXP,50000,0.021654253,2309015.230,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,brainpoolP512r1,512,MODEXP,50000,0.008376452,5969114.352,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,brainpoolP512r1,512,MODEXP,50000,0.010625510,4705656.477,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,brainpoolP512r1,512,EXPONENTIATION,390,0.015669152,24889.669,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,brainpoolP512r1,512,EXPONENTIATION,390,0.008209123,47508.120,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,brainpoolP512r1,512,EXPONENTIATION,390,0.197575054,1973.933,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,brainpoolP512r1,512,EXPONENTIATION,50000,0.182487641,273991.158,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,brainpoolP512r1,512,EXPONENTIATION,50000,0.185169347,270023.094,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,brainpoolP512r1,512,EXPONENTIATION,50000,0.045847444,1090573.336,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,brainpoolP512r1,512,EXPONENTIATION,50000,0.048271423,1035809.531,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,brainpoolP512r1,512,EXPONENTIATION,50000,0.011367081,4398666.632,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,brainpoolP512r1,512,EXPONENTIATION,50000,0.013955796,3582740.814,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,brainpoolP512r1,512,EXPONENTIATION,50000,0.008392205,5957909.753,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,brainpoolP512r1,512,EXPONENTIATION,50000,0.010738626,4656089.147,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,brainpoolP512r1,512,EXPONENTIATION,50000,0.011199641,4464428.799,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,brainpoolP512r1,512,EXPONENTIATION,50000,0.013653368,3662100.079,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,brainpoolP512r1,512,EXPONENTIATION,50000,0.011525972,4338028.935,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,brainpoolP512r1,512,EXPONENTIATION,50000,0.013918972,3592219.312,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,brainpoolP512r1,512,DIVIDE,3125,0.000655860,4764738.774,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,brainpoolP512r1,512,DIVIDE,3125,0.033081184,94464.574,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,brainpoolP512r1,512,DIVIDE,3125,0.030376306,102876.236,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,brainpoolP512r1,512,DIVIDE,50000,0.000058624,852893013.100,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,brainpoolP512r1,512,DIVIDE,50000,0.001119479,44663635.506,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,brainpoolP512r1,512,DIVIDE,50000,0.003645556,13715329.060,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,brainpoolP512r1,512,DIVIDE,50000,0.001138540,43915892.869,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,brainpoolP512r1,512,DIVIDE,50000,0.003485225,14346276.204,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,brainpoolP512r1,512,DIVIDE,50000,0.000384634,129993704.443,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,brainpoolP512r1,512,DIVIDE,50000,0.002626638,19035740.844,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,brainpoolP512r1,512,DIVIDE,50000,0.000356190,140374510.148,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,brainpoolP512r1,512,DIVIDE,50000,0.002700243,18516852.203,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,brainpoolP512r1,512,DIVIDE,50000,0.000370358,135004508.965,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,brainpoolP512r1,512,DIVIDE,50000,0.002740561,18244439.916,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,brainpoolP512r1,512,DIVIDE,50000,0.000336923,148401862.034,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,brainpoolP512r1,512,DIVIDE,50000,0.002593660,19277777.136,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,brainpoolP512r1,512,ISQRT,781,0.000669368,1166772.247,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,brainpoolP512r1,512,ISQRT,781,0.001075159,726404.187,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,brainpoolP512r1,512,ISQRT,50000,0.017026900,2936529.840,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,brainpoolP512r1,512,ISQRT,50000,0.019702931,2537693.504,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,brainpoolP512r1,512,ISQRT,50000,0.016971721,2946077.179,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,brainpoolP512r1,512,ISQRT,50000,0.018551477,2695203.191,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,brainpoolP512r1,512,ISQRT,50000,0.003301640,15143989.126,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,brainpoolP512r1,512,ISQRT,50000,0.005467785,9144470.655,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,brainpoolP512r1,512,ISQRT,50000,0.003003335,16648159.597,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,brainpoolP512r1,512,ISQRT,50000,0.005130399,9745830.650,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,brainpoolP512r1,512,ISQRT,50000,0.003147825,15883983.301,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,brainpoolP512r1,512,ISQRT,50000,0.005304816,9425397.593,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,brainpoolP512r1,512,ISQRT,50000,0.003959517,12627802.853,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,brainpoolP512r1,512,ISQRT,50000,0.006072535,8233793.623,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL_R2,25000,0.009460600,2642538.535,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL_R2,25000,0.004963171,5037102.260,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL_R2,25000,0.004403400,5677431.073,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,brainpoolP512r1,512,MODMUL_R2,50000,0.000048992,1020574787.720,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,brainpoolP512r1,512,MODMUL_R2,50000,0.000388593,128669325.733,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,brainpoolP512r1,512,MODMUL_R2,50000,0.002785433,17950530.327,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,brainpoolP512r1,512,MODMUL_R2,50000,0.000071965,694781784.446,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,brainpoolP512r1,512,MODMUL_R2,50000,0.002194841,22780693.006,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,brainpoolP512r1,512,MODMUL_R2,50000,0.000054847,911627165.211,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,brainpoolP512r1,512,MODMUL_R2,50000,0.002165305,23091434.934,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,brainpoolP512r1,512,MODMUL_R2,50000,0.000031486,1588008428.525,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,brainpoolP512r1,512,MODMUL_R2,50000,0.002064147,24223081.309,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,brainpoolP512r1,512,MODMUL_R2,50000,0.000051039,979642284.374,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,brainpoolP512r1,512,MODMUL_R2,50000,0.002075162,24094504.223,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,brainpoolP512r1,512,MODMUL_R2,50000,0.000030820,1622321843.908,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,brainpoolP512r1,512,MODMUL_R2,50000,0.001964232,25455242.148,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,p1024,1024,ADD,12500,0.000965835,12942173.595,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,p1024,1024,ADD,12500,0.001566950,7977279.031,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,p1024,1024,ADD,12500,0.001365127,9156654.986,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,p1024,1024,ADD,50000,0.000063488,787550403.226,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,p1024,1024,ADD,50000,0.000185302,269829801.499,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,p1024,1024,ADD,50000,0.004880508,10244835.168,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,p1024,1024,ADD,50000,0.000097358,513568695.388,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,p1024,1024,ADD,50000,0.003218691,15534265.418,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32,p1024,1024,ADD,50000,0.000040053,1248346486.655,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32,p1024,1024,ADD,50000,0.004140355,12076259.208,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,p1024,1024,ADD,50000,0.000041466,1205807938.460,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,p1024,1024,ADD,50000,0.004141926,12071678.706,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,p1024,1024,ADD,50000,0.000042114,1187253662.616,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,p1024,1024,ADD,50000,0.004124197,12123572.168,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,p1024,1024,ADD,50000,0.000027140,1842301756.312,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,p1024,1024,ADD,50000,0.004033867,12395054.240,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,p1024,1024,ADD,50000,0.000027159,1841010611.503,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,p1024,1024,ADD,50000,0.004123220,12126444.903,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,p1024,1024,SUBTRACT,12500,0.000579221,21580709.282,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,p1024,1024,SUBTRACT,12500,0.001347818,9274250.966,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,p1024,1024,SUBTRACT,12500,0.001173955,10647772.118,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,p1024,1024,SUBTRACT,50000,0.000063488,787550403.226,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,p1024,1024,SUBTRACT,50000,0.000184051,271663841.173,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,p1024,1024,SUBTRACT,50000,0.004884619,10236212.827,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,p1024,1024,SUBTRACT,50000,0.000097937,510532234.918,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,p1024,1024,SUBTRACT,50000,0.004178659,11965561.198,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32,p1024,1024,SUBTRACT,50000,0.000040745,1227144031.520,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32,p1024,1024,SUBTRACT,50000,0.004152721,12040298.319,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,p1024,1024,SUBTRACT,50000,0.000041229,1212738309.410,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,p1024,1024,SUBTRACT,50000,0.004135944,12089138.516,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,p1024,1024,SUBTRACT,50000,0.000041396,1207847637.367,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,p1024,1024,SUBTRACT,50000,0.004143553,12066938.602,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,p1024,1024,SUBTRACT,50000,0.000027982,1786863353.213,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,p1024,1024,SUBTRACT,50000,0.004134535,12093258.275,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,p1024,1024,SUBTRACT,50000,0.000027164,1840671344.879,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,p1024,1024,SUBTRACT,50000,0.004133251,12097015.070,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,p1024,1024,ADDMOD,12500,0.002131441,5864578.433,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,p1024,1024,ADDMOD,12500,0.002211584,5652057.501,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,p1024,1024,ADDMOD,12500,0.002463627,5073821.051,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,p1024,1024,ADDMOD,50000,0.000062976,793953252.033,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,p1024,1024,ADDMOD,50000,0.000242215,206428174.032,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,p1024,1024,ADDMOD,50000,0.004820976,10371343.896,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,p1024,1024,ADDMOD,50000,0.000126755,394461814.612,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,p1024,1024,ADDMOD,50000,0.004215630,11860623.430,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32,p1024,1024,ADDMOD,50000,0.000046424,1077029456.416,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32,p1024,1024,ADDMOD,50000,0.004143259,12067794.876,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,p1024,1024,ADDMOD,50000,0.000040022,1249312377.304,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,p1024,1024,ADDMOD,50000,0.004151323,12044352.934,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,p1024,1024,ADDMOD,50000,0.000039766,1257356777.629,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,p1024,1024,ADDMOD,50000,0.004126221,12117625.305,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,p1024,1024,ADDMOD,50000,0.000016570,3017495491.989,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,p1024,1024,ADDMOD,50000,0.004116199,12147129.002,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,p1024,1024,ADDMOD,50000,0.000016738,2987211089.280,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,p1024,1024,ADDMOD,50000,0.004121398,12131805.675,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,p1024,1024,SUBTRACTMOD,12500,0.001776122,7037804.862,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,p1024,1024,SUBTRACTMOD,12500,0.002489763,5020557.200,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,p1024,1024,SUBTRACTMOD,12500,0.002209788,5656651.266,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,p1024,1024,SUBTRACTMOD,50000,0.000062336,802104722.793,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,p1024,1024,SUBTRACTMOD,50000,0.000252497,198022174.318,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,p1024,1024,SUBTRACTMOD,50000,0.004961634,10077325.340,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,p1024,1024,SUBTRACTMOD,50000,0.000126627,394860367.853,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,p1024,1024,SUBTRACTMOD,50000,0.004219096,11850879.962,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32,p1024,1024,SUBTRACTMOD,50000,0.000046031,1086225267.891,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32,p1024,1024,SUBTRACTMOD,50000,0.004157625,12026096.565,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,p1024,1024,SUBTRACTMOD,50000,0.000038904,1285214177.861,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,p1024,1024,SUBTRACTMOD,50000,0.004154060,12036417.449,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,p1024,1024,SUBTRACTMOD,50000,0.000039836,1255144758.014,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,p1024,1024,SUBTRACTMOD,50000,0.004038254,12381588.757,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,p1024,1024,SUBTRACTMOD,50000,0.000017324,2886169997.043,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,p1024,1024,SUBTRACTMOD,50000,0.003995044,12515506.586,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,p1024,1024,SUBTRACTMOD,50000,0.000017211,2905125977.232,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,p1024,1024,SUBTRACTMOD,50000,0.004046181,12357331.567,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.005229922,2390093.002,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.007763183,1610164.289,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.004643241,2692085.118,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.022131642,2259208.786,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.028038919,1783235.650,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.006800665,7352222.176,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.011778069,4245178.065,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.004367856,11447263.783,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.009333646,5356963.375,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.000089971,555734629.324,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.004717633,10598535.298,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.000079302,630501115.827,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.004732973,10564184.498,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.000070464,709582375.128,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.004734349,10561114.039,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.000061241,816446784.373,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.004703222,10631009.937,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.005635142,2218222.726,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.008034778,1555736.823,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.004554742,2744392.538,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000063488,787550403.226,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.008247428,6062496.096,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.013035897,3835562.680,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.002122825,23553520.002,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.006836291,7313907.503,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000548521,91154213.113,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.005197255,9620463.004,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000546851,91432580.540,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.005205606,9605029.666,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000547098,91391305.354,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.005078008,9846380.724,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000192729,259431637.793,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.004853652,10301521.464,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000180274,277355557.421,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.004825947,10360660.838,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.016608434,752629.657,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.005077523,2461830.304,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.005311480,2353393.043,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000068032,734948259.643,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.001786211,27992213.973,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.005925429,8438207.603,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000246923,202492267.212,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.004320128,11573731.133,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000063987,781408929.416,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.003595410,13906619.933,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000110171,453839788.157,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.004053498,12335025.095,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000075936,658449496.156,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.003989506,12532880.045,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000091249,547951169.588,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.004178578,11965793.050,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000054765,912991423.980,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.004129712,12107381.810,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,p1024,1024,COMPARE,12500,0.000358965,34822312.394,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,p1024,1024,COMPARE,12500,0.001861432,6715258.295,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,p1024,1024,COMPARE,12500,0.001030689,12127809.608,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,p1024,1024,COMPARE,50000,0.000063488,787550403.226,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,p1024,1024,COMPARE,50000,0.000086230,579844648.513,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,p1024,1024,COMPARE,50000,0.004171810,11985205.472,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,p1024,1024,COMPARE,50000,0.000046792,1068559311.340,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,p1024,1024,COMPARE,50000,0.004135608,12090120.631,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,p1024,1024,COMPARE,50000,0.000021441,2331978094.985,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,p1024,1024,COMPARE,50000,0.004118790,12139487.464,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,p1024,1024,COMPARE,50000,0.000022176,2254688446.930,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,p1024,1024,COMPARE,50000,0.003973284,12584048.743,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,p1024,1024,COMPARE,50000,0.000016425,3044143820.766,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,p1024,1024,COMPARE,50000,0.003970031,12594359.953,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,p1024,1024,COMPARE,50000,0.000017237,2900730960.051,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,p1024,1024,COMPARE,50000,0.004105592,12178511.733,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,p1024,1024,REDUCE,1562,0.000170896,9140083.866,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,p1024,1024,REDUCE,1562,0.000406474,3842805.028,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,p1024,1024,REDUCE,1562,0.000494619,3157986.793,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,p1024,1024,REDUCE,50000,0.000061440,813802083.333,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,p1024,1024,REDUCE,50000,0.002110666,23689205.712,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,p1024,1024,REDUCE,50000,0.006251146,7998533.354,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,p1024,1024,REDUCE,50000,0.001121795,44571421.977,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,p1024,1024,REDUCE,50000,0.005236021,9549236.044,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,p1024,1024,REDUCE,50000,0.000309907,161338715.689,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,p1024,1024,REDUCE,50000,0.004414115,11327298.794,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,p1024,1024,REDUCE,50000,0.000329990,151519743.879,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,p1024,1024,REDUCE,50000,0.004447888,11241290.246,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,p1024,1024,REDUCE,50000,0.000301421,165880955.119,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,p1024,1024,REDUCE,50000,0.004340386,11519712.715,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,p1024,1024,REDUCE,50000,0.000314890,158785607.888,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,p1024,1024,REDUCE,50000,0.004345528,11506081.589,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,p1024,1024,MODMUL,781,0.002608640,299389.720,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,p1024,1024,MODMUL,781,0.002630661,296883.614,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,p1024,1024,MODMUL,781,0.046076976,16949.897,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,p1024,1024,MODMUL,50000,0.000327232,152796792.490,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,p1024,1024,MODMUL,50000,0.007962047,6279792.114,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,p1024,1024,MODMUL,50000,0.012177974,4105773.244,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,p1024,1024,MODMUL,50000,0.002963484,16872032.870,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,p1024,1024,MODMUL,50000,0.007100865,7041395.649,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,p1024,1024,MODMUL,50000,0.001402744,35644421.786,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,p1024,1024,MODMUL,50000,0.005498398,9093557.735,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,p1024,1024,MODMUL,50000,0.001031284,48483250.586,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,p1024,1024,MODMUL,50000,0.005149513,9709656.007,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,p1024,1024,MODMUL,50000,0.001430633,34949563.986,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,p1024,1024,MODMUL,50000,0.005488428,9110076.699,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,p1024,1024,MODMUL,50000,0.001050635,47590266.470,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,p1024,1024,MODMUL,50000,0.005183174,9646598.705,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,p1024,1024,MODEXP,195,0.047929081,4068.511,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,p1024,1024,MODEXP,195,0.010751013,18137.826,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,p1024,1024,MODEXP,195,0.099657757,1956.697,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,p1024,1024,MODEXP,50000,1.544298410,32377.162,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,p1024,1024,MODEXP,50000,2.604771836,19195.539,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,p1024,1024,MODEXP,50000,2.611796451,19143.911,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,p1024,1024,MODEXP,50000,0.150533145,332152.763,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,p1024,1024,MODEXP,50000,0.153681018,325349.224,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,p1024,1024,MODEXP,50000,0.151507549,330016.559,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,p1024,1024,MODEXP,50000,0.156139770,320225.910,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,p1024,1024,MODEXP,50000,0.078884206,633840.442,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,p1024,1024,MODEXP,50000,0.081965343,610013.918,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,p1024,1024,MODEXP,50000,0.153137704,326503.524,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,p1024,1024,MODEXP,50000,0.157798152,316860.491,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,p1024,1024,MODEXP,50000,0.079033555,632642.679,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,p1024,1024,MODEXP,50000,0.083749804,597016.323,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,p1024,1024,EXPONENTIATION,195,0.016962396,11496.017,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,p1024,1024,EXPONENTIATION,195,0.007223065,26996.850,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,p1024,1024,EXPONENTIATION,195,0.207966485,937.651,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,p1024,1024,EXPONENTIATION,50000,1.572323294,31800.076,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,p1024,1024,EXPONENTIATION,50000,1.574135072,31763.475,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,p1024,1024,EXPONENTIATION,50000,0.359606495,139040.870,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,p1024,1024,EXPONENTIATION,50000,0.362923758,137769.983,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,p1024,1024,EXPONENTIATION,50000,0.106358986,470106.024,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,p1024,1024,EXPONENTIATION,50000,0.111131420,449917.764,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,p1024,1024,EXPONENTIATION,50000,0.093062563,537272.974,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,p1024,1024,EXPONENTIATION,50000,0.097724621,511641.790,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,p1024,1024,EXPONENTIATION,50000,0.106203068,470796.192,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,p1024,1024,EXPONENTIATION,50000,0.110514823,452427.997,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,p1024,1024,EXPONENTIATION,50000,0.093552240,534460.746,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,p1024,1024,EXPONENTIATION,50000,0.097809197,511199.371,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,p1024,1024,DIVIDE,1562,0.000349323,4471506.271,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,p1024,1024,DIVIDE,1562,0.000797968,1957471.982,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,p1024,1024,DIVIDE,1562,0.001054262,1481605.380,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,p1024,1024,DIVIDE,50000,0.000082496,606089992.242,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,p1024,1024,DIVIDE,50000,0.009913388,5043684.359,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,p1024,1024,DIVIDE,50000,0.014510686,3445736.477,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,p1024,1024,DIVIDE,50000,0.007101395,7040870.146,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,p1024,1024,DIVIDE,50000,0.011970293,4177007.193,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,p1024,1024,DIVIDE,50000,0.001505231,33217492.442,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,p1024,1024,DIVIDE,50000,0.006190450,8076957.263,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,p1024,1024,DIVIDE,50000,0.001537491,32520514.646,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,p1024,1024,DIVIDE,50000,0.006216012,8043742.516,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,p1024,1024,DIVIDE,50000,0.001505528,33210938.785,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,p1024,1024,DIVIDE,50000,0.006172286,8100726.379,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,p1024,1024,DIVIDE,50000,0.001553081,32194071.508,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,p1024,1024,DIVIDE,50000,0.006241003,8011532.770,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,p1024,1024,ISQRT,390,0.000512133,761520.928,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,p1024,1024,ISQRT,390,0.001271002,306844.522,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,p1024,1024,ISQRT,50000,0.189735191,263525.178,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,p1024,1024,ISQRT,50000,0.194999752,256410.583,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,p1024,1024,ISQRT,50000,0.119997818,416674.243,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,p1024,1024,ISQRT,50000,0.124362010,402052.042,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,p1024,1024,ISQRT,50000,0.022818925,2191163.691,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,p1024,1024,ISQRT,50000,0.027580727,1812860.117,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,p1024,1024,ISQRT,50000,0.021946537,2278263.769,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,p1024,1024,ISQRT,50000,0.026489591,1887533.862,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,p1024,1024,ISQRT,50000,0.022694592,2203168.048,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,p1024,1024,ISQRT,50000,0.026835546,1863200.397,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,p1024,1024,ISQRT,50000,0.022222315,2249990.604,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,p1024,1024,ISQRT,50000,0.026950899,1855225.684,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,p1024,1024,MODMUL_R2,12500,0.011519456,1085120.687,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,p1024,1024,MODMUL_R2,12500,0.004647503,2689616.324,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,p1024,1024,MODMUL_R2,12500,0.010564195,1183242.075,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,p1024,1024,MODMUL_R2,50000,0.000122816,407113079.729,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,p1024,1024,MODMUL_R2,50000,0.002412703,20723644.768,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,p1024,1024,MODMUL_R2,50000,0.006549135,7634596.058,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,p1024,1024,MODMUL_R2,50000,0.000222038,225186667.589,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,p1024,1024,MODMUL_R2,50000,0.004270822,11707348.183,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,p1024,1024,MODMUL_R2,50000,0.000179293,278873155.316,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,p1024,1024,MODMUL_R2,50000,0.004301889,11622801.095,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,p1024,1024,MODMUL_R2,50000,0.000118237,422879511.012,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,p1024,1024,MODMUL_R2,50000,0.004201899,11899381.708,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,p1024,1024,MODMUL_R2,50000,0.000159737,313014490.677,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,p1024,1024,MODMUL_R2,50000,0.004247072,11772816.721,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,p1024,1024,MODMUL_R2,50000,0.000096379,518785091.006,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,p1024,1024,MODMUL_R2,50000,0.004180703,11959711.022,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,p2048,2048,ADD,6250,0.000599530,10424832.757,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,p2048,2048,ADD,6250,0.001451656,4305426.520,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,p2048,2048,ADD,6250,0.000935910,6677989.704,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,p2048,2048,ADD,50000,0.000121280,412269129.288,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,p2048,2048,ADD,50000,0.000356439,140276465.823,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,p2048,2048,ADD,50000,0.008242275,6066286.280,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,p2048,2048,ADD,50000,0.000182063,274630168.562,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,p2048,2048,ADD,50000,0.008043957,6215846.255,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32,p2048,2048,ADD,50000,0.000072177,692741129.943,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32,p2048,2048,ADD,50000,0.008089272,6181025.965,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,p2048,2048,ADD,50000,0.000071492,699379069.762,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,p2048,2048,ADD,50000,0.008119764,6157814.453,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,p2048,2048,ADD,50000,0.000073950,676132163.356,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,p2048,2048,ADD,50000,0.008093060,6178132.874,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,p2048,2048,ADD,50000,0.000040944,1221179391.339,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,p2048,2048,ADD,50000,0.007907201,6323350.061,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,p2048,2048,ADD,50000,0.000043529,1148660321.841,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,p2048,2048,ADD,50000,0.008056739,6205984.831,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,p2048,2048,SUBTRACT,6250,0.000516495,12100790.327,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,p2048,2048,SUBTRACT,6250,0.000889691,7024907.768,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,p2048,2048,SUBTRACT,6250,0.000884607,7065284.423,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,p2048,2048,SUBTRACT,50000,0.000121504,411509086.121,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,p2048,2048,SUBTRACT,50000,0.000355266,140739601.116,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,p2048,2048,SUBTRACT,50000,0.008243823,6065147.190,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,p2048,2048,SUBTRACT,50000,0.000183398,272631077.894,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,p2048,2048,SUBTRACT,50000,0.007009708,7132964.708,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32,p2048,2048,SUBTRACT,50000,0.000072233,692204110.064,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32,p2048,2048,SUBTRACT,50000,0.008098094,6174292.395,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,p2048,2048,SUBTRACT,50000,0.000072813,686690461.885,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,p2048,2048,SUBTRACT,50000,0.008084470,6184697.328,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,p2048,2048,SUBTRACT,50000,0.000074293,673010927.462,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,p2048,2048,SUBTRACT,50000,0.008106908,6167579.547,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,p2048,2048,SUBTRACT,50000,0.000040852,1223930312.667,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,p2048,2048,SUBTRACT,50000,0.008078046,6189615.681,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,p2048,2048,SUBTRACT,50000,0.000043163,1158400638.676,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,p2048,2048,SUBTRACT,50000,0.008030949,6225914.238,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,p2048,2048,ADDMOD,6250,0.001361477,4590601.546,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,p2048,2048,ADDMOD,6250,0.001591285,3927644.168,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,p2048,2048,ADDMOD,6250,0.001469965,4251800.972,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,p2048,2048,ADDMOD,50000,0.000120288,415669060.920,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,p2048,2048,ADDMOD,50000,0.000451899,110644188.405,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,p2048,2048,ADDMOD,50000,0.008358297,5982079.845,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,p2048,2048,ADDMOD,50000,0.000229488,217876305.987,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,p2048,2048,ADDMOD,50000,0.008129203,6150664.467,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32,p2048,2048,ADDMOD,50000,0.000087896,568854018.100,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32,p2048,2048,ADDMOD,50000,0.008092786,6178342.081,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,p2048,2048,ADDMOD,50000,0.000069001,724627273.187,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,p2048,2048,ADDMOD,50000,0.008090788,6179867.790,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,p2048,2048,ADDMOD,50000,0.000070300,711238026.602,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,p2048,2048,ADDMOD,50000,0.003554753,14065675.007,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,p2048,2048,ADDMOD,50000,0.000023128,2161881780.659,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,p2048,2048,ADDMOD,50000,0.008031802,6225253.093,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,p2048,2048,ADDMOD,50000,0.000024054,2078659240.546,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,p2048,2048,ADDMOD,50000,0.008044059,6215767.453,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,p2048,2048,SUBTRACTMOD,6250,0.001184497,5276502.848,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,p2048,2048,SUBTRACTMOD,6250,0.002654238,2354725.074,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,p2048,2048,SUBTRACTMOD,6250,0.001604396,3895546.994,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,p2048,2048,SUBTRACTMOD,50000,0.000118784,420932112.069,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,p2048,2048,SUBTRACTMOD,50000,0.000460113,108668959.468,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,p2048,2048,SUBTRACTMOD,50000,0.008376051,5969400.116,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,p2048,2048,SUBTRACTMOD,50000,0.000236691,211245855.701,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,p2048,2048,SUBTRACTMOD,50000,0.008096170,6175759.671,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32,p2048,2048,SUBTRACTMOD,50000,0.000091592,545899058.873,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32,p2048,2048,SUBTRACTMOD,50000,0.008106797,6167663.997,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,p2048,2048,SUBTRACTMOD,50000,0.000068370,731314690.087,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,p2048,2048,SUBTRACTMOD,50000,0.008090708,6179928.878,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,p2048,2048,SUBTRACTMOD,50000,0.000070762,706593878.333,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,p2048,2048,SUBTRACTMOD,50000,0.003463981,14434259.173,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,p2048,2048,SUBTRACTMOD,50000,0.000023410,2135837991.352,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,p2048,2048,SUBTRACTMOD,50000,0.008026493,6229370.649,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,p2048,2048,SUBTRACTMOD,50000,0.000025979,1924634975.746,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,p2048,2048,SUBTRACTMOD,50000,0.008025073,6230472.889,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.007099654,880324.590,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.004478922,1395425.061,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.008605078,726315.324,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.088575919,564487.510,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.098210866,509108.636,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.023286085,2147205.080,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.032796588,1524548.833,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.006991647,7151390.786,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.016302621,3066991.499,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.000248141,201498337.849,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.009381581,5329592.112,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.000221322,225915182.852,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.004610579,10844624.949,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.000213324,234385280.820,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.009369729,5336333.650,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.000187123,267203979.659,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.009287313,5383688.459,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.006535882,956259.612,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.008220269,760315.752,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.004489267,1392209.477,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.000161216,310142913.855,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.032606749,1533424.873,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.042131234,1186767.992,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.008281767,6037358.964,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.017487801,2859135.915,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.002101849,23788578.443,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.011275469,4434405.354,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.002095628,23859196.979,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.011302469,4423812.175,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.002095398,23861815.624,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.006802591,7350140.588,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.000828569,60345004.586,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.009956858,5021664.463,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.000770734,64873221.293,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.009791617,5106408.885,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.019734492,316704.377,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.008211039,761170.419,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.008079101,773600.926,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000279520,178878076.703,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.038827364,1287751.596,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.047518525,1052221.213,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.001273397,39265052.238,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.009185912,5443117.696,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000195996,255107295.619,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.008189456,6105411.675,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000344677,145063343.548,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.008263309,6050844.792,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000262539,190447881.667,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.005305690,9423844.980,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000306274,163252489.325,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.008320702,6009108.361,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000207440,241033556.167,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.008235220,6071483.228,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,p2048,2048,COMPARE,6250,0.000138900,44996537.151,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,p2048,2048,COMPARE,6250,0.001161380,5381527.573,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,p2048,2048,COMPARE,6250,0.000925389,6753913.351,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,p2048,2048,COMPARE,50000,0.000121536,411400737.230,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,p2048,2048,COMPARE,50000,0.000158370,315716346.807,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,p2048,2048,COMPARE,50000,0.008102740,6170752.082,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,p2048,2048,COMPARE,50000,0.000082296,607563164.558,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,p2048,2048,COMPARE,50000,0.007968025,6275080.720,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,p2048,2048,COMPARE,50000,0.000032551,1536053651.776,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,p2048,2048,COMPARE,50000,0.008062543,6201517.314,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,p2048,2048,COMPARE,50000,0.000034131,1464944180.256,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,p2048,2048,COMPARE,50000,0.004409806,11338367.297,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,p2048,2048,COMPARE,50000,0.000024411,2048255891.344,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,p2048,2048,COMPARE,50000,0.007997165,6252215.608,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,p2048,2048,COMPARE,50000,0.000024215,2064838487.536,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,p2048,2048,COMPARE,50000,0.008006603,6244845.667,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,p2048,2048,REDUCE,781,0.000138342,5645414.063,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,p2048,2048,REDUCE,781,0.000324924,2403642.235,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,p2048,2048,REDUCE,781,0.000482819,1617581.733,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,p2048,2048,REDUCE,50000,0.000118784,420932112.069,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,p2048,2048,REDUCE,50000,0.105691392,473075.423,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,p2048,2048,REDUCE,50000,0.115021665,434700.715,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,p2048,2048,REDUCE,50000,0.004677636,10689160.023,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,p2048,2048,REDUCE,50000,0.012699920,3937032.661,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,p2048,2048,REDUCE,50000,0.000956306,52284519.744,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,p2048,2048,REDUCE,50000,0.009037027,5532793.019,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,p2048,2048,REDUCE,50000,0.000794156,62959918.502,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,p2048,2048,REDUCE,50000,0.004275373,11694886.013,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,p2048,2048,REDUCE,50000,0.000869564,57500079.604,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,p2048,2048,REDUCE,50000,0.008930532,5598770.581,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,p2048,2048,REDUCE,50000,0.000857009,58342442.535,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,p2048,2048,REDUCE,50000,0.008883251,5628569.995,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,p2048,2048,MODMUL,390,0.002824580,138073.649,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,p2048,2048,MODMUL,390,0.003491459,111701.154,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,p2048,2048,MODMUL,390,0.003561674,109499.088,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,p2048,2048,MODMUL,50000,0.001370848,36473773.898,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,p2048,2048,MODMUL,50000,0.174335941,286802.593,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,p2048,2048,MODMUL,50000,0.182784742,273545.808,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,p2048,2048,MODMUL,50000,0.014991731,3335171.895,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,p2048,2048,MODMUL,50000,0.023265833,2149074.135,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,p2048,2048,MODMUL,50000,0.005169577,9671971.292,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,p2048,2048,MODMUL,50000,0.013413397,3727616.514,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,p2048,2048,MODMUL,50000,0.003538353,14130868.154,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,p2048,2048,MODMUL,50000,0.007352334,6800561.526,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,p2048,2048,MODMUL,50000,0.005220460,9577700.037,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,p2048,2048,MODMUL,50000,0.013409331,3728746.786,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,p2048,2048,MODMUL,50000,0.003509234,14248123.678,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,p2048,2048,MODMUL,50000,0.011640560,4295325.999,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,p2048,2048,MODEXP,97,0.165104487,587.507,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,p2048,2048,MODEXP,97,0.285934265,339.239,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,p2048,2048,MODEXP,97,0.303195852,319.925,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,p2048,2048,MODEXP,50000,1.972019672,25354.717,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,p2048,2048,MODEXP,50000,55.649428856,898.482,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,p2048,2048,MODEXP,50000,55.662626701,898.269,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,p2048,2048,MODEXP,50000,5.144748090,9718.649,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,p2048,2048,MODEXP,50000,5.164749828,9681.011,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,p2048,2048,MODEXP,50000,1.647840801,30342.737,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,p2048,2048,MODEXP,50000,1.653014005,30247.778,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,p2048,2048,MODEXP,50000,0.990895486,50459.408,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,p2048,2048,MODEXP,50000,1.001218188,49939.165,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,p2048,2048,MODEXP,50000,1.665073105,30028.712,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,p2048,2048,MODEXP,50000,1.671785845,29908.137,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,p2048,2048,MODEXP,50000,0.989646079,50523.112,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,p2048,2048,MODEXP,50000,0.996141372,50193.679,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,p2048,2048,EXPONENTIATION,97,0.021147438,4586.844,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,p2048,2048,EXPONENTIATION,97,0.011976605,8099.123,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,p2048,2048,EXPONENTIATION,97,0.363137891,267.116,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,p2048,2048,EXPONENTIATION,50000,13.162552446,3798.655,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,p2048,2048,EXPONENTIATION,50000,13.165950792,3797.675,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,p2048,2048,EXPONENTIATION,50000,3.105069145,16102.701,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,p2048,2048,EXPONENTIATION,50000,3.109187242,16081.373,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,p2048,2048,EXPONENTIATION,50000,0.890170092,56169.041,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,p2048,2048,EXPONENTIATION,50000,0.898127281,55671.397,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,p2048,2048,EXPONENTIATION,50000,0.783433054,63821.663,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,p2048,2048,EXPONENTIATION,50000,0.794435597,62937.764,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,p2048,2048,EXPONENTIATION,50000,0.909816235,54956.153,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,p2048,2048,EXPONENTIATION,50000,0.916203283,54573.042,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,p2048,2048,EXPONENTIATION,50000,0.789826683,63305.028,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,p2048,2048,EXPONENTIATION,50000,0.798642854,62606.207,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,p2048,2048,DIVIDE,781,0.000271242,2879343.316,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,p2048,2048,DIVIDE,781,0.001110900,703033.454,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,p2048,2048,DIVIDE,781,0.001118247,698414.695,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,p2048,2048,DIVIDE,50000,0.000116480,429258241.758,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,p2048,2048,DIVIDE,50000,0.410232634,121882.064,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,p2048,2048,DIVIDE,50000,0.418375780,119509.786,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,p2048,2048,DIVIDE,50000,0.146067978,342306.375,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,p2048,2048,DIVIDE,50000,0.154070829,324526.066,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,p2048,2048,DIVIDE,50000,0.014342158,3486225.717,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,p2048,2048,DIVIDE,50000,0.023939506,2088597.815,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,p2048,2048,DIVIDE,50000,0.013149733,3802358.568,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,p2048,2048,DIVIDE,50000,0.022785536,2194374.535,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,p2048,2048,DIVIDE,50000,0.014335178,3487923.202,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,p2048,2048,DIVIDE,50000,0.023776222,2102941.332,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,p2048,2048,DIVIDE,50000,0.013104125,3815592.423,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,p2048,2048,DIVIDE,50000,0.022982820,2175538.076,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,p2048,2048,ISQRT,195,0.000538883,361859.642,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,p2048,2048,ISQRT,195,0.009262638,21052.318,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,p2048,2048,ISQRT,50000,6.923342597,7221.945,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,p2048,2048,ISQRT,50000,6.921819846,7223.534,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,p2048,2048,ISQRT,50000,4.467184167,11192.733,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,p2048,2048,ISQRT,50000,4.496805855,11119.003,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,p2048,2048,ISQRT,50000,0.150296710,332675.280,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,p2048,2048,ISQRT,50000,0.158827739,314806.471,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,p2048,2048,ISQRT,50000,0.051034331,979732.643,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,p2048,2048,ISQRT,50000,0.059597086,838967.194,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,p2048,2048,ISQRT,50000,0.148353562,337032.690,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,p2048,2048,ISQRT,50000,0.156898461,318677.441,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,p2048,2048,ISQRT,50000,0.051051347,979406.087,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,p2048,2048,ISQRT,50000,0.059582156,839177.421,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-1t,p2048,2048,MODMUL_R2,6250,0.014004466,446286.206,0
library,Intel(R) Xeon(R) 6787P,host-cpu,gmp-nt,p2048,2048,MODMUL_R2,6250,0.007179678,870512.583,0
library,Intel(R) Xeon(R) 6787P,host-cpu,openssl-nt,p2048,2048,MODMUL_R2,6250,0.007550830,827723.572,0
library,NVIDIA B300 SXM6 AC,gpu,cgbn,p2048,2048,MODMUL_R2,50000,0.000440032,113628099.775,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w8,p2048,2048,MODMUL_R2,50000,0.028973125,1725737.215,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w8,p2048,2048,MODMUL_R2,50000,0.037962613,1317085.313,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w16,p2048,2048,MODMUL_R2,50000,0.001092144,45781508.788,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w16,p2048,2048,MODMUL_R2,50000,0.008969162,5574656.802,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-opt,p2048,2048,MODMUL_R2,50000,0.000618280,80869512.912,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-opt,p2048,2048,MODMUL_R2,50000,0.008642384,5785440.706,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-o64,p2048,2048,MODMUL_R2,50000,0.000423091,118177886.617,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-o64,p2048,2048,MODMUL_R2,50000,0.008443215,5921914.823,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il,p2048,2048,MODMUL_R2,50000,0.000571763,87448825.423,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il,p2048,2048,MODMUL_R2,50000,0.008611755,5806017.486,0
opencl-kernel,NVIDIA B300 SXM6 AC,GPU,w32-il64,p2048,2048,MODMUL_R2,50000,0.000381079,131206377.626,0
opencl-e2e,NVIDIA B300 SXM6 AC,GPU,w32-il64,p2048,2048,MODMUL_R2,50000,0.008387351,5961357.763,0
```
