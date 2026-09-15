# MPA-OpenCL benchmark report - Tesla V100-PCIE-32GB

> **Note.** The multi-threaded GMP and OpenSSL baseline columns have been
> removed from this report: they predate the 2026-09-12 timing fix and were
> understated (see `reports/README.md`). The single-threaded GMP column, the
> OpenCL-on-CPU rows and all MPA measurements are unaffected and were verified
> against GMP before timing.


> **Partial report.** The run was interrupted or hit its time budget.
> Rows that never ran are marked `n/a`.

## 1. System under test

3 OpenCL device(s) exercised with the identical kernels and operands.

### Device 0 - Tesla V100-PCIE-32GB (GPU)

| Property | Value |
|---|---|
| Model | Tesla V100-PCIE-32GB |
| Type | GPU |
| Vendor | NVIDIA Corporation |
| Device memory | 31.73 GiB |
| Max single allocation | 7.93 GiB |
| Local memory | 48 KiB |
| Global cache | 2560 KiB |
| Compute units | 80 |
| Max clock | 1380 MHz |
| Max work-group size | 1024 |
| OpenCL version | OpenCL 3.0 CUDA |
| Driver | 570.133.07 |

### Device 1 - Tesla V100-PCIE-32GB (GPU)

| Property | Value |
|---|---|
| Model | Tesla V100-PCIE-32GB |
| Type | GPU |
| Vendor | NVIDIA Corporation |
| Device memory | 31.73 GiB |
| Max single allocation | 7.93 GiB |
| Local memory | 48 KiB |
| Global cache | 2560 KiB |
| Compute units | 80 |
| Max clock | 1380 MHz |
| Max work-group size | 1024 |
| OpenCL version | OpenCL 3.0 CUDA |
| Driver | 570.133.07 |

### Device 2 - cpu-haswell-Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz (CPU)

| Property | Value |
|---|---|
| Model | cpu-haswell-Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz |
| Type | CPU |
| Vendor | GenuineIntel |
| Device memory | 155.24 GiB |
| Max single allocation | 64.00 GiB |
| Local memory | 256 KiB |
| Global cache | 46080 KiB |
| Compute units | 72 |
| Max clock | 3000 MHz |
| Max work-group size | 4096 |
| OpenCL version | OpenCL 3.0 PoCL HSTR: cpu-x86_64-pc-linux-gnu-haswell |
| Driver | 5.0+debian |

### Host

| Property | Value |
|---|---|
| CPU | Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz |
| Logical cores | 72 |
| OpenMP threads used | 72 |
| RAM | 157.2 GB |
| OS | Ubuntu 24.04.4 LTS |
| Kernel | 5.15.0-190-generic |
| Arch | x86_64 |
| GMP | 6.3.0 |
| OpenSSL | OpenSSL 3.0.13 30 Jan 2024 |
| CGBN | cgbn_results.tsv loaded |

## 2. Method

- Workload auto-sized from the device and host: --min-items from 700 x compute units, --items from ten times that capped by host RAM. Either flag, given explicitly, overrides its half.
- Base workload 50000 items, scaled down per operator by its cost weight and by modulus size. Device rows honour --min-items (56000) so the GPU is not left idle; the CPU libraries keep the smaller count because a full-width MODEXP there costs minutes. Both counts appear in every row as dev/cpu, and throughput is per-second so they remain comparable.
- 5 timed repetitions, **minimum** reported. Two untimed warm-up launches precede them.
- `kernel` times `clEnqueueNDRangeKernel` + `clFinish` only. `e2e` adds the host->device operand writes and the device->host result read.
- Every OpenCL device runs the same kernels on the same operands, so GPU and CPU-OpenCL columns are directly comparable.
- CPU library baselines (GMP, OpenSSL) run those same operands, with temporaries preallocated outside the timed region, so the figure is the arithmetic and not marshalling. The generator is reseeded per modulus and operation so every backend sees identical inputs.
- OpenSSL rows time the nearest BN primitive, which is not always semantically identical (its Montgomery routine expects Montgomery-domain inputs); they measure comparable work, not identical results. Correctness is judged against GMP only.
- Every device cell is checked word-for-word against GMP before it is timed. A cell that mismatches is reported and excluded from the speedup tables.
- Total wall time 5528.1 s.

## 3. Correctness

| Device | Kernel | Configs run | Passed | Mismatched | Build/launch failed |
|---|---|---|---|---|---|
| [0] GPU | `mpaKernels_8bits.cl` (w8) | 75 | 75 | 0 | 0 |
| [0] GPU | `mpaKernel_16bits.cl` (w16) | 75 | 75 | 0 | 0 |
| [0] GPU | `mpaKernel_32bits.cl` (w32) | 35 | 35 | 0 | 0 |
| [0] GPU | `mpaKernel_32bits_opt.cl` (w32-opt) | 75 | 75 | 0 | 0 |
| [0] GPU | `mpaKernel_32bits_opt.cl` (w32-o64) | 71 | 71 | 0 | 0 |

**All configurations correct** - 331 configurations, 0 problems.

## 4. Throughput per device

Operations per second, higher is better. Kernel-only timings.

### Device 0 - Tesla V100-PCIE-32GB (GPU)

#### secp256k1 (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | 822.02 M | 1.38 G | 2.24 G | 2.11 G | 2.19 G | - | - | 8.77 M | 1.93 G |
| SUBTRACT | 50000 / 50000 | 835.53 M | 1.41 G | 2.32 G | 2.18 G | 2.34 G | - | - | 9.19 M | 2.21 G |
| ADDMOD | 50000 / 50000 | 575.37 M | 999.66 M | 1.76 G | 2.03 G | 2.20 G | - | - | 2.26 M | 2.17 G |
| SUBTRACTMOD | 50000 / 50000 | 565.34 M | 990.77 M | 1.70 G | 2.02 G | 2.18 G | - | - | 4.78 M | 2.24 G |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | 24.63 M | 72.99 M | 327.59 M | 1.14 G | 1.37 G | - | - | 4.16 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | 138.04 M | 428.04 M | 1.17 G | 1.10 G | 1.15 G | - | - | 3.92 M | 2.18 G |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | 275.47 M | 814.31 M | 1.94 G | 1.32 G | 1.64 G | - | - | 1.73 M | 2.26 G |
| COMPARE | 50000 / 50000 | 845.27 M | 1.41 G | - | 2.17 G | 2.26 G | - | - | 11.63 M | 2.25 G |
| REDUCE | 50000 / 6250 | 161.74 M | 239.95 M | - | 613.59 M | 697.49 M | - | - | 14.23 M | 2.24 G |
| MODMUL | 50000 / 3125 | 59.33 M | 103.30 M | - | 216.64 M | 309.68 M | - | - | 2.43 M | 1.19 G |
| MODEXP | 50000 / 781 | 1.32 M | 7.29 M | - | 6.98 M | 19.73 M | - | - | 38.05 k | 628.30 k |
| EXPONENTIATION | 50000 / 781 | 897.89 k | 2.51 M | - | 51.46 M | 63.08 M | - | - | 127.37 k | n/a |
| DIVIDE | 50000 / 6250 | 85.76 M | 116.83 M | - | 303.44 M | 275.37 M | - | - | 5.66 M | 2.00 G |
| ISQRT | 50000 / 1562 | 6.90 M | 8.73 M | - | 27.84 M | 40.36 M | - | - | 4.27 M | n/a |
| MODMUL_R2 | 50000 / 50000 | 241.39 M | 742.20 M | - | 1.07 G | 1.26 G | - | - | 1.66 M | 2.01 G |

#### rsa256(composite) (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | 854.54 M | 1.43 G | 2.35 G | 2.16 G | 2.09 G | - | - | 3.80 M | 2.22 G |
| SUBTRACT | 50000 / 50000 | 841.67 M | 1.39 G | 2.33 G | 2.12 G | 2.16 G | - | - | 4.24 M | 2.12 G |
| ADDMOD | 50000 / 50000 | 615.79 M | 1.07 G | 1.87 G | 1.84 G | 2.01 G | - | - | 2.71 M | 2.21 G |
| SUBTRACTMOD | 50000 / 50000 | 573.78 M | 971.50 M | 1.73 G | 2.04 G | 1.98 G | - | - | 3.02 M | 2.14 G |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | 24.87 M | 73.16 M | 330.68 M | 1.14 G | 1.24 G | - | - | 4.99 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | 139.52 M | 428.37 M | 1.13 G | 1.11 G | 1.03 G | - | - | 4.71 M | 2.26 G |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | 282.10 M | 838.28 M | 1.96 G | 1.34 G | 1.49 G | - | - | 1.37 M | 2.18 G |
| COMPARE | 50000 / 50000 | 844.55 M | 1.37 G | - | 2.16 G | 2.15 G | - | - | 6.68 M | 2.18 G |
| REDUCE | 50000 / 6250 | 162.33 M | 239.45 M | - | 615.66 M | 624.67 M | - | - | 7.41 M | 2.21 G |
| MODMUL | 50000 / 3125 | 60.01 M | 103.64 M | - | 211.46 M | 277.44 M | - | - | 2.77 M | 1.18 G |
| MODEXP | 50000 / 781 | 1.32 M | 7.31 M | - | 6.95 M | 17.68 M | - | - | 62.04 k | 648.94 k |
| EXPONENTIATION | 50000 / 781 | 896.67 k | 2.51 M | - | 51.38 M | 62.93 M | - | - | 166.78 k | n/a |
| DIVIDE | 50000 / 6250 | 85.22 M | 113.87 M | - | 303.83 M | 256.10 M | - | - | 7.31 M | 1.98 G |
| ISQRT | 50000 / 1562 | 6.87 M | 8.43 M | - | 27.80 M | 37.93 M | - | - | 4.25 M | n/a |
| MODMUL_R2 | 50000 / 50000 | 240.05 M | 722.45 M | - | 1.09 G | 1.23 G | - | - | 2.00 M | 1.95 G |

#### brainpoolP512r1 (512-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 25000 | 390.24 M | 670.70 M | 926.72 M | 1.01 G | 912.74 M | - | - | 3.27 M | 1.79 G |
| SUBTRACT | 50000 / 25000 | 387.28 M | 648.38 M | 920.86 M | 1.01 G | 943.75 M | - | - | 3.92 M | 1.80 G |
| ADDMOD | 50000 / 25000 | 274.93 M | 506.00 M | 794.38 M | 862.08 M | 804.09 M | - | - | 2.52 M | 1.86 G |
| SUBTRACTMOD | 50000 / 25000 | 251.71 M | 465.29 M | 734.01 M | 891.28 M | 795.30 M | - | - | 2.94 M | 1.85 G |
| MULTIPLYOPERANDSCANNING | 50000 / 25000 | 4.47 M | 16.97 M | 21.42 M | 331.88 M | 304.74 M | - | - | 2.62 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 25000 | 20.65 M | 74.15 M | 196.69 M | 200.19 M | 181.76 M | - | - | 2.44 M | 1.86 G |
| MONTGOMERYMULTIPLICATION | 50000 / 25000 | 77.80 M | 261.23 M | 847.92 M | 591.17 M | 706.72 M | - | - | 984.19 k | 1.80 G |
| COMPARE | 50000 / 25000 | 392.68 M | 686.37 M | - | 1.01 G | 971.55 M | - | - | 13.15 M | 1.81 G |
| REDUCE | 50000 / 3125 | 56.11 M | 69.84 M | - | 210.17 M | 218.85 M | - | - | 9.07 M | 1.72 G |
| MODMUL | 50000 / 1562 | 19.82 M | 27.84 M | - | 66.83 M | 73.87 M | - | - | 3.35 M | 414.35 M |
| MODEXP | 50000 / 390 | 87.81 k | 1.04 M | - | 992.11 k | 1.58 M | - | - | 13.51 k | 61.86 k |
| EXPONENTIATION | 50000 / 390 | 101.32 k | 407.62 k | - | 1.17 M | 1.28 M | - | - | 30.39 k | n/a |
| DIVIDE | 50000 / 3125 | 27.90 M | 27.34 M | - | 90.98 M | 95.68 M | - | - | 6.73 M | 817.63 M |
| ISQRT | 50000 / 781 | 1.29 M | 1.32 M | - | 5.59 M | 7.17 M | - | - | 2.56 M | n/a |
| MODMUL_R2 | 50000 / 25000 | 54.75 M | 270.37 M | - | 343.59 M | 516.63 M | - | - | 1.28 M | 1.25 G |

#### p1024 (1024-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 12500 | 131.61 M | 249.48 M | 126.30 M | 453.81 M | 474.59 M | - | - | 3.07 M | 1.14 G |
| SUBTRACT | 50000 / 12500 | 137.19 M | 257.85 M | 127.84 M | 451.35 M | 451.39 M | - | - | 4.67 M | 1.15 G |
| ADDMOD | 50000 / 12500 | 90.05 M | 175.50 M | 120.68 M | 355.08 M | 354.18 M | - | - | 2.17 M | 1.15 G |
| SUBTRACTMOD | 50000 / 12500 | 90.74 M | 177.50 M | 117.36 M | 349.46 M | 353.54 M | - | - | 2.41 M | 1.15 G |
| MULTIPLYOPERANDSCANNING | 50000 / 12500 | 964.76 k | 3.80 M | 2.89 M | 134.37 M | 135.90 M | - | - | 1.03 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 12500 | 2.76 M | 10.73 M | 28.53 M | 40.52 M | 39.75 M | - | - | 1.09 M | 1.02 G |
| MONTGOMERYMULTIPLICATION | 50000 / 12500 | 15.25 M | 71.90 M | 262.42 M | 178.04 M | 220.01 M | - | - | 296.93 k | 760.34 M |
| COMPARE | 50000 / 12500 | 153.44 M | 293.50 M | - | 485.46 M | 536.22 M | - | - | 7.45 M | 1.16 G |
| REDUCE | 50000 / 1562 | 10.43 M | 18.82 M | - | 76.27 M | 82.58 M | - | - | 10.89 M | 1.06 G |
| MODMUL | 50000 / 781 | 2.48 M | 7.18 M | - | 16.56 M | 21.96 M | - | - | 714.80 k | 143.47 M |
| MODEXP | 50000 / 195 | 9.89 k | 75.35 k | - | 140.45 k | 265.56 k | - | - | 2.37 k | 24.53 k |
| EXPONENTIATION | 50000 / 195 | 10.49 k | 49.92 k | - | 190.79 k | 202.29 k | - | - | 16.38 k | n/a |
| DIVIDE | 50000 / 1562 | 630.66 k | 4.29 M | - | 15.67 M | 17.40 M | - | - | 15.05 M | 561.85 M |
| ISQRT | 50000 / 390 | 45.54 k | 174.92 k | - | 943.77 k | 993.47 k | - | - | 2.63 M | n/a |
| MODMUL_R2 | 50000 / 12500 | 10.10 M | 66.99 M | - | 112.26 M | 152.79 M | - | - | 1.59 M | 417.67 M |

#### p2048 (2048-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 6250 | 66.65 M | 134.86 M | 264.14 M | 264.44 M | 267.27 M | - | - | 23.97 M | 604.92 M |
| SUBTRACT | 50000 / 6250 | 57.96 M | 135.10 M | 249.08 M | 260.54 M | 269.67 M | - | - | 26.53 M | 603.05 M |
| ADDMOD | 50000 / 6250 | 42.08 M | 95.22 M | 182.26 M | 181.68 M | 182.72 M | - | - | 8.73 M | 582.59 M |
| SUBTRACTMOD | 50000 / 6250 | 42.18 M | 89.78 M | 174.13 M | 183.20 M | 184.45 M | - | - | 11.32 M | 579.99 M |
| MULTIPLYOPERANDSCANNING | 50000 / 6250 | 218.31 k | 729.06 k | 2.78 M | 30.18 M | 35.07 M | - | - | 1.46 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 6250 | 696.99 k | 2.75 M | 10.91 M | 10.51 M | 10.58 M | - | - | 1.46 M | 322.16 M |
| MONTGOMERYMULTIPLICATION | 50000 / 6250 | 412.66 k | 17.08 M | 92.57 M | 55.82 M | 68.27 M | - | - | 226.64 k | 200.04 M |
| COMPARE | 50000 / 6250 | 77.98 M | 155.34 M | - | 296.15 M | 297.16 M | - | - | 61.28 M | 603.98 M |
| REDUCE | 50000 / 781 | 60.77 k | 4.38 M | - | 23.55 M | 26.16 M | - | - | 15.57 M | 593.66 M |
| MODMUL | 50000 / 390 | 42.46 k | 1.48 M | - | 3.97 M | 5.43 M | - | - | 419.99 k | 42.35 M |
| MODEXP | 50000 / 97 | 239.8 | 1.92 k | - | 13.04 k | over budget | - | - | 326.1 | 24.79 k |
| EXPONENTIATION | 50000 / 97 | 1.46 k | 5.26 k | - | 17.87 k | - | - | - | 2.84 k | n/a |
| DIVIDE | 50000 / 781 | 25.87 k | 101.86 k | - | 1.39 M | - | - | - | 9.53 M | 502.90 M |
| ISQRT | 50000 / 195 | 1.87 k | 4.50 k | - | 409.37 k | - | - | - | 1.27 M | n/a |
| MODMUL_R2 | 50000 / 6250 | 579.24 k | 14.80 M | - | 30.90 M | - | - | - | 407.25 k | 115.07 M |

### Device 1 - Tesla V100-PCIE-32GB (GPU)

#### secp256k1 (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | - | - | - | - | - | - | - | 8.77 M | 1.93 G |
| SUBTRACT | 50000 / 50000 | - | - | - | - | - | - | - | 9.19 M | 2.21 G |
| ADDMOD | 50000 / 50000 | - | - | - | - | - | - | - | 2.26 M | 2.17 G |
| SUBTRACTMOD | 50000 / 50000 | - | - | - | - | - | - | - | 4.78 M | 2.24 G |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | - | - | - | - | - | - | - | 4.16 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | - | - | - | - | - | - | - | 3.92 M | 2.18 G |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | - | - | - | - | - | - | - | 1.73 M | 2.26 G |
| COMPARE | 50000 / 50000 | - | - | - | - | - | - | - | 11.63 M | 2.25 G |
| REDUCE | 50000 / 6250 | - | - | - | - | - | - | - | 14.23 M | 2.24 G |
| MODMUL | 50000 / 3125 | - | - | - | - | - | - | - | 2.43 M | 1.19 G |
| MODEXP | 50000 / 781 | - | - | - | - | - | - | - | 38.05 k | 628.30 k |
| EXPONENTIATION | 50000 / 781 | - | - | - | - | - | - | - | 127.37 k | n/a |
| DIVIDE | 50000 / 6250 | - | - | - | - | - | - | - | 5.66 M | 2.00 G |
| ISQRT | 50000 / 1562 | - | - | - | - | - | - | - | 4.27 M | n/a |
| MODMUL_R2 | 50000 / 50000 | - | - | - | - | - | - | - | 1.66 M | 2.01 G |

#### rsa256(composite) (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | - | - | - | - | - | - | - | 3.80 M | 2.22 G |
| SUBTRACT | 50000 / 50000 | - | - | - | - | - | - | - | 4.24 M | 2.12 G |
| ADDMOD | 50000 / 50000 | - | - | - | - | - | - | - | 2.71 M | 2.21 G |
| SUBTRACTMOD | 50000 / 50000 | - | - | - | - | - | - | - | 3.02 M | 2.14 G |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | - | - | - | - | - | - | - | 4.99 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | - | - | - | - | - | - | - | 4.71 M | 2.26 G |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | - | - | - | - | - | - | - | 1.37 M | 2.18 G |
| COMPARE | 50000 / 50000 | - | - | - | - | - | - | - | 6.68 M | 2.18 G |
| REDUCE | 50000 / 6250 | - | - | - | - | - | - | - | 7.41 M | 2.21 G |
| MODMUL | 50000 / 3125 | - | - | - | - | - | - | - | 2.77 M | 1.18 G |
| MODEXP | 50000 / 781 | - | - | - | - | - | - | - | 62.04 k | 648.94 k |
| EXPONENTIATION | 50000 / 781 | - | - | - | - | - | - | - | 166.78 k | n/a |
| DIVIDE | 50000 / 6250 | - | - | - | - | - | - | - | 7.31 M | 1.98 G |
| ISQRT | 50000 / 1562 | - | - | - | - | - | - | - | 4.25 M | n/a |
| MODMUL_R2 | 50000 / 50000 | - | - | - | - | - | - | - | 2.00 M | 1.95 G |

#### brainpoolP512r1 (512-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 25000 | - | - | - | - | - | - | - | 3.27 M | 1.79 G |
| SUBTRACT | 50000 / 25000 | - | - | - | - | - | - | - | 3.92 M | 1.80 G |
| ADDMOD | 50000 / 25000 | - | - | - | - | - | - | - | 2.52 M | 1.86 G |
| SUBTRACTMOD | 50000 / 25000 | - | - | - | - | - | - | - | 2.94 M | 1.85 G |
| MULTIPLYOPERANDSCANNING | 50000 / 25000 | - | - | - | - | - | - | - | 2.62 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 25000 | - | - | - | - | - | - | - | 2.44 M | 1.86 G |
| MONTGOMERYMULTIPLICATION | 50000 / 25000 | - | - | - | - | - | - | - | 984.19 k | 1.80 G |
| COMPARE | 50000 / 25000 | - | - | - | - | - | - | - | 13.15 M | 1.81 G |
| REDUCE | 50000 / 3125 | - | - | - | - | - | - | - | 9.07 M | 1.72 G |
| MODMUL | 50000 / 1562 | - | - | - | - | - | - | - | 3.35 M | 414.35 M |
| MODEXP | 50000 / 390 | - | - | - | - | - | - | - | 13.51 k | 61.86 k |
| EXPONENTIATION | 50000 / 390 | - | - | - | - | - | - | - | 30.39 k | n/a |
| DIVIDE | 50000 / 3125 | - | - | - | - | - | - | - | 6.73 M | 817.63 M |
| ISQRT | 50000 / 781 | - | - | - | - | - | - | - | 2.56 M | n/a |
| MODMUL_R2 | 50000 / 25000 | - | - | - | - | - | - | - | 1.28 M | 1.25 G |

#### p1024 (1024-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 12500 | - | - | - | - | - | - | - | 3.07 M | 1.14 G |
| SUBTRACT | 50000 / 12500 | - | - | - | - | - | - | - | 4.67 M | 1.15 G |
| ADDMOD | 50000 / 12500 | - | - | - | - | - | - | - | 2.17 M | 1.15 G |
| SUBTRACTMOD | 50000 / 12500 | - | - | - | - | - | - | - | 2.41 M | 1.15 G |
| MULTIPLYOPERANDSCANNING | 50000 / 12500 | - | - | - | - | - | - | - | 1.03 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 12500 | - | - | - | - | - | - | - | 1.09 M | 1.02 G |
| MONTGOMERYMULTIPLICATION | 50000 / 12500 | - | - | - | - | - | - | - | 296.93 k | 760.34 M |
| COMPARE | 50000 / 12500 | - | - | - | - | - | - | - | 7.45 M | 1.16 G |
| REDUCE | 50000 / 1562 | - | - | - | - | - | - | - | 10.89 M | 1.06 G |
| MODMUL | 50000 / 781 | - | - | - | - | - | - | - | 714.80 k | 143.47 M |
| MODEXP | 50000 / 195 | - | - | - | - | - | - | - | 2.37 k | 24.53 k |
| EXPONENTIATION | 50000 / 195 | - | - | - | - | - | - | - | 16.38 k | n/a |
| DIVIDE | 50000 / 1562 | - | - | - | - | - | - | - | 15.05 M | 561.85 M |
| ISQRT | 50000 / 390 | - | - | - | - | - | - | - | 2.63 M | n/a |
| MODMUL_R2 | 50000 / 12500 | - | - | - | - | - | - | - | 1.59 M | 417.67 M |

#### p2048 (2048-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 6250 | - | - | - | - | - | - | - | 23.97 M | 604.92 M |
| SUBTRACT | 50000 / 6250 | - | - | - | - | - | - | - | 26.53 M | 603.05 M |
| ADDMOD | 50000 / 6250 | - | - | - | - | - | - | - | 8.73 M | 582.59 M |
| SUBTRACTMOD | 50000 / 6250 | - | - | - | - | - | - | - | 11.32 M | 579.99 M |
| MULTIPLYOPERANDSCANNING | 50000 / 6250 | - | - | - | - | - | - | - | 1.46 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 6250 | - | - | - | - | - | - | - | 1.46 M | 322.16 M |
| MONTGOMERYMULTIPLICATION | 50000 / 6250 | - | - | - | - | - | - | - | 226.64 k | 200.04 M |
| COMPARE | 50000 / 6250 | - | - | - | - | - | - | - | 61.28 M | 603.98 M |
| REDUCE | 50000 / 781 | - | - | - | - | - | - | - | 15.57 M | 593.66 M |
| MODMUL | 50000 / 390 | - | - | - | - | - | - | - | 419.99 k | 42.35 M |
| MODEXP | 50000 / 97 | - | - | - | - | - | - | - | 326.1 | 24.79 k |
| EXPONENTIATION | 50000 / 97 | - | - | - | - | - | - | - | 2.84 k | n/a |
| DIVIDE | 50000 / 781 | - | - | - | - | - | - | - | 9.53 M | 502.90 M |
| ISQRT | 50000 / 195 | - | - | - | - | - | - | - | 1.27 M | n/a |
| MODMUL_R2 | 50000 / 6250 | - | - | - | - | - | - | - | 407.25 k | 115.07 M |

### Device 2 - cpu-haswell-Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz (CPU)

#### secp256k1 (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | - | - | - | - | - | - | - | 8.77 M | 1.93 G |
| SUBTRACT | 50000 / 50000 | - | - | - | - | - | - | - | 9.19 M | 2.21 G |
| ADDMOD | 50000 / 50000 | - | - | - | - | - | - | - | 2.26 M | 2.17 G |
| SUBTRACTMOD | 50000 / 50000 | - | - | - | - | - | - | - | 4.78 M | 2.24 G |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | - | - | - | - | - | - | - | 4.16 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | - | - | - | - | - | - | - | 3.92 M | 2.18 G |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | - | - | - | - | - | - | - | 1.73 M | 2.26 G |
| COMPARE | 50000 / 50000 | - | - | - | - | - | - | - | 11.63 M | 2.25 G |
| REDUCE | 50000 / 6250 | - | - | - | - | - | - | - | 14.23 M | 2.24 G |
| MODMUL | 50000 / 3125 | - | - | - | - | - | - | - | 2.43 M | 1.19 G |
| MODEXP | 50000 / 781 | - | - | - | - | - | - | - | 38.05 k | 628.30 k |
| EXPONENTIATION | 50000 / 781 | - | - | - | - | - | - | - | 127.37 k | n/a |
| DIVIDE | 50000 / 6250 | - | - | - | - | - | - | - | 5.66 M | 2.00 G |
| ISQRT | 50000 / 1562 | - | - | - | - | - | - | - | 4.27 M | n/a |
| MODMUL_R2 | 50000 / 50000 | - | - | - | - | - | - | - | 1.66 M | 2.01 G |

#### rsa256(composite) (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | - | - | - | - | - | - | - | 3.80 M | 2.22 G |
| SUBTRACT | 50000 / 50000 | - | - | - | - | - | - | - | 4.24 M | 2.12 G |
| ADDMOD | 50000 / 50000 | - | - | - | - | - | - | - | 2.71 M | 2.21 G |
| SUBTRACTMOD | 50000 / 50000 | - | - | - | - | - | - | - | 3.02 M | 2.14 G |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | - | - | - | - | - | - | - | 4.99 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | - | - | - | - | - | - | - | 4.71 M | 2.26 G |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | - | - | - | - | - | - | - | 1.37 M | 2.18 G |
| COMPARE | 50000 / 50000 | - | - | - | - | - | - | - | 6.68 M | 2.18 G |
| REDUCE | 50000 / 6250 | - | - | - | - | - | - | - | 7.41 M | 2.21 G |
| MODMUL | 50000 / 3125 | - | - | - | - | - | - | - | 2.77 M | 1.18 G |
| MODEXP | 50000 / 781 | - | - | - | - | - | - | - | 62.04 k | 648.94 k |
| EXPONENTIATION | 50000 / 781 | - | - | - | - | - | - | - | 166.78 k | n/a |
| DIVIDE | 50000 / 6250 | - | - | - | - | - | - | - | 7.31 M | 1.98 G |
| ISQRT | 50000 / 1562 | - | - | - | - | - | - | - | 4.25 M | n/a |
| MODMUL_R2 | 50000 / 50000 | - | - | - | - | - | - | - | 2.00 M | 1.95 G |

#### brainpoolP512r1 (512-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 25000 | - | - | - | - | - | - | - | 3.27 M | 1.79 G |
| SUBTRACT | 50000 / 25000 | - | - | - | - | - | - | - | 3.92 M | 1.80 G |
| ADDMOD | 50000 / 25000 | - | - | - | - | - | - | - | 2.52 M | 1.86 G |
| SUBTRACTMOD | 50000 / 25000 | - | - | - | - | - | - | - | 2.94 M | 1.85 G |
| MULTIPLYOPERANDSCANNING | 50000 / 25000 | - | - | - | - | - | - | - | 2.62 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 25000 | - | - | - | - | - | - | - | 2.44 M | 1.86 G |
| MONTGOMERYMULTIPLICATION | 50000 / 25000 | - | - | - | - | - | - | - | 984.19 k | 1.80 G |
| COMPARE | 50000 / 25000 | - | - | - | - | - | - | - | 13.15 M | 1.81 G |
| REDUCE | 50000 / 3125 | - | - | - | - | - | - | - | 9.07 M | 1.72 G |
| MODMUL | 50000 / 1562 | - | - | - | - | - | - | - | 3.35 M | 414.35 M |
| MODEXP | 50000 / 390 | - | - | - | - | - | - | - | 13.51 k | 61.86 k |
| EXPONENTIATION | 50000 / 390 | - | - | - | - | - | - | - | 30.39 k | n/a |
| DIVIDE | 50000 / 3125 | - | - | - | - | - | - | - | 6.73 M | 817.63 M |
| ISQRT | 50000 / 781 | - | - | - | - | - | - | - | 2.56 M | n/a |
| MODMUL_R2 | 50000 / 25000 | - | - | - | - | - | - | - | 1.28 M | 1.25 G |

#### p1024 (1024-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 12500 | - | - | - | - | - | - | - | 3.07 M | 1.14 G |
| SUBTRACT | 50000 / 12500 | - | - | - | - | - | - | - | 4.67 M | 1.15 G |
| ADDMOD | 50000 / 12500 | - | - | - | - | - | - | - | 2.17 M | 1.15 G |
| SUBTRACTMOD | 50000 / 12500 | - | - | - | - | - | - | - | 2.41 M | 1.15 G |
| MULTIPLYOPERANDSCANNING | 50000 / 12500 | - | - | - | - | - | - | - | 1.03 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 12500 | - | - | - | - | - | - | - | 1.09 M | 1.02 G |
| MONTGOMERYMULTIPLICATION | 50000 / 12500 | - | - | - | - | - | - | - | 296.93 k | 760.34 M |
| COMPARE | 50000 / 12500 | - | - | - | - | - | - | - | 7.45 M | 1.16 G |
| REDUCE | 50000 / 1562 | - | - | - | - | - | - | - | 10.89 M | 1.06 G |
| MODMUL | 50000 / 781 | - | - | - | - | - | - | - | 714.80 k | 143.47 M |
| MODEXP | 50000 / 195 | - | - | - | - | - | - | - | 2.37 k | 24.53 k |
| EXPONENTIATION | 50000 / 195 | - | - | - | - | - | - | - | 16.38 k | n/a |
| DIVIDE | 50000 / 1562 | - | - | - | - | - | - | - | 15.05 M | 561.85 M |
| ISQRT | 50000 / 390 | - | - | - | - | - | - | - | 2.63 M | n/a |
| MODMUL_R2 | 50000 / 12500 | - | - | - | - | - | - | - | 1.59 M | 417.67 M |

#### p2048 (2048-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 6250 | - | - | - | - | - | - | - | 23.97 M | 604.92 M |
| SUBTRACT | 50000 / 6250 | - | - | - | - | - | - | - | 26.53 M | 603.05 M |
| ADDMOD | 50000 / 6250 | - | - | - | - | - | - | - | 8.73 M | 582.59 M |
| SUBTRACTMOD | 50000 / 6250 | - | - | - | - | - | - | - | 11.32 M | 579.99 M |
| MULTIPLYOPERANDSCANNING | 50000 / 6250 | - | - | - | - | - | - | - | 1.46 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 6250 | - | - | - | - | - | - | - | 1.46 M | 322.16 M |
| MONTGOMERYMULTIPLICATION | 50000 / 6250 | - | - | - | - | - | - | - | 226.64 k | 200.04 M |
| COMPARE | 50000 / 6250 | - | - | - | - | - | - | - | 61.28 M | 603.98 M |
| REDUCE | 50000 / 781 | - | - | - | - | - | - | - | 15.57 M | 593.66 M |
| MODMUL | 50000 / 390 | - | - | - | - | - | - | - | 419.99 k | 42.35 M |
| MODEXP | 50000 / 97 | - | - | - | - | - | - | - | 326.1 | 24.79 k |
| EXPONENTIATION | 50000 / 97 | - | - | - | - | - | - | - | 2.84 k | n/a |
| DIVIDE | 50000 / 781 | - | - | - | - | - | - | - | 9.53 M | 502.90 M |
| ISQRT | 50000 / 195 | - | - | - | - | - | - | - | 1.27 M | n/a |
| MODMUL_R2 | 50000 / 6250 | - | - | - | - | - | - | - | 407.25 k | 115.07 M |

## 5. Head to head

Best OpenCL GPU result against best OpenCL CPU result and the CPU libraries.
Ratios above 1.00x mean the GPU is faster than that baseline.

### secp256k1 (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w32 | 2.24 G | none | n/a | 8.77 M | 1.93 G | n/a | 1.16x |
| SUBTRACT | w32-o64 | 2.34 G | none | n/a | 9.19 M | 2.21 G | n/a | 1.06x |
| ADDMOD | w32-o64 | 2.20 G | none | n/a | 2.26 M | 2.17 G | n/a | 1.01x |
| SUBTRACTMOD | w32-o64 | 2.18 G | none | n/a | 4.78 M | 2.24 G | n/a | 0.97x |
| MULTIPLYOPERANDSCANNING | w32-o64 | 1.37 G | none | n/a | 4.16 M | n/a | n/a | n/a |
| MULTIPLYPRODUCTSCANNING | w32 | 1.17 G | none | n/a | 3.92 M | 2.18 G | n/a | 0.54x |
| MONTGOMERYMULTIPLICATION | w32 | 1.94 G | none | n/a | 1.73 M | 2.26 G | n/a | 0.86x |
| COMPARE | w32-o64 | 2.26 G | none | n/a | 11.63 M | 2.25 G | n/a | 1.01x |
| REDUCE | w32-o64 | 87.19 M | none | n/a | 14.23 M | 2.24 G | n/a | 0.04x |
| MODMUL | w32-o64 | 19.35 M | none | n/a | 2.43 M | 1.19 G | n/a | 0.02x |
| MODEXP | w32-o64 | 308.25 k | none | n/a | 38.05 k | 628.30 k | n/a | 0.49x |
| EXPONENTIATION | w32-o64 | 985.35 k | none | n/a | 127.37 k | n/a | n/a | n/a |
| DIVIDE | w32-opt | 37.93 M | none | n/a | 5.66 M | 2.00 G | n/a | 0.02x |
| ISQRT | w32-o64 | 1.26 M | none | n/a | 4.27 M | n/a | n/a | n/a |
| MODMUL_R2 | w32-o64 | 1.26 G | none | n/a | 1.66 M | 2.01 G | n/a | 0.63x |

### rsa256(composite) (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w32 | 2.35 G | none | n/a | 3.80 M | 2.22 G | n/a | 1.06x |
| SUBTRACT | w32 | 2.33 G | none | n/a | 4.24 M | 2.12 G | n/a | 1.10x |
| ADDMOD | w32-o64 | 2.01 G | none | n/a | 2.71 M | 2.21 G | n/a | 0.91x |
| SUBTRACTMOD | w32-opt | 2.04 G | none | n/a | 3.02 M | 2.14 G | n/a | 0.95x |
| MULTIPLYOPERANDSCANNING | w32-o64 | 1.24 G | none | n/a | 4.99 M | n/a | n/a | n/a |
| MULTIPLYPRODUCTSCANNING | w32 | 1.13 G | none | n/a | 4.71 M | 2.26 G | n/a | 0.50x |
| MONTGOMERYMULTIPLICATION | w32 | 1.96 G | none | n/a | 1.37 M | 2.18 G | n/a | 0.90x |
| COMPARE | w32-opt | 2.16 G | none | n/a | 6.68 M | 2.18 G | n/a | 0.99x |
| REDUCE | w32-o64 | 78.08 M | none | n/a | 7.41 M | 2.21 G | n/a | 0.04x |
| MODMUL | w32-o64 | 17.34 M | none | n/a | 2.77 M | 1.18 G | n/a | 0.01x |
| MODEXP | w32-o64 | 276.22 k | none | n/a | 62.04 k | 648.94 k | n/a | 0.43x |
| EXPONENTIATION | w32-o64 | 982.98 k | none | n/a | 166.78 k | n/a | n/a | n/a |
| DIVIDE | w32-opt | 37.98 M | none | n/a | 7.31 M | 1.98 G | n/a | 0.02x |
| ISQRT | w32-o64 | 1.18 M | none | n/a | 4.25 M | n/a | n/a | n/a |
| MODMUL_R2 | w32-o64 | 1.23 G | none | n/a | 2.00 M | 1.95 G | n/a | 0.63x |

### brainpoolP512r1 (512-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w32-opt | 502.93 M | none | n/a | 3.27 M | 1.79 G | n/a | 0.28x |
| SUBTRACT | w32-opt | 506.23 M | none | n/a | 3.92 M | 1.80 G | n/a | 0.28x |
| ADDMOD | w32-opt | 431.04 M | none | n/a | 2.52 M | 1.86 G | n/a | 0.23x |
| SUBTRACTMOD | w32-opt | 445.64 M | none | n/a | 2.94 M | 1.85 G | n/a | 0.24x |
| MULTIPLYOPERANDSCANNING | w32-opt | 165.94 M | none | n/a | 2.62 M | n/a | n/a | n/a |
| MULTIPLYPRODUCTSCANNING | w32-opt | 100.10 M | none | n/a | 2.44 M | 1.86 G | n/a | 0.05x |
| MONTGOMERYMULTIPLICATION | w32 | 423.96 M | none | n/a | 984.19 k | 1.80 G | n/a | 0.23x |
| COMPARE | w32-opt | 505.86 M | none | n/a | 13.15 M | 1.81 G | n/a | 0.28x |
| REDUCE | w32-o64 | 13.68 M | none | n/a | 9.07 M | 1.72 G | n/a | 0.01x |
| MODMUL | w32-o64 | 2.31 M | none | n/a | 3.35 M | 414.35 M | n/a | 0.01x |
| MODEXP | w32-o64 | 12.33 k | none | n/a | 13.51 k | 61.86 k | n/a | 0.20x |
| EXPONENTIATION | w32-o64 | 10.01 k | none | n/a | 30.39 k | n/a | n/a | n/a |
| DIVIDE | w32-o64 | 5.98 M | none | n/a | 6.73 M | 817.63 M | n/a | 0.01x |
| ISQRT | w32-o64 | 111.95 k | none | n/a | 2.56 M | n/a | n/a | n/a |
| MODMUL_R2 | w32-o64 | 258.31 M | none | n/a | 1.28 M | 1.25 G | n/a | 0.21x |

### p1024 (1024-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w32-o64 | 118.65 M | none | n/a | 3.07 M | 1.14 G | n/a | 0.10x |
| SUBTRACT | w32-o64 | 112.85 M | none | n/a | 4.67 M | 1.15 G | n/a | 0.10x |
| ADDMOD | w32-opt | 88.77 M | none | n/a | 2.17 M | 1.15 G | n/a | 0.08x |
| SUBTRACTMOD | w32-o64 | 88.39 M | none | n/a | 2.41 M | 1.15 G | n/a | 0.08x |
| MULTIPLYOPERANDSCANNING | w32-o64 | 33.97 M | none | n/a | 1.03 M | n/a | n/a | n/a |
| MULTIPLYPRODUCTSCANNING | w32-opt | 10.13 M | none | n/a | 1.09 M | 1.02 G | n/a | 0.01x |
| MONTGOMERYMULTIPLICATION | w32 | 65.60 M | none | n/a | 296.93 k | 760.34 M | n/a | 0.09x |
| COMPARE | w32-o64 | 134.05 M | none | n/a | 7.45 M | 1.16 G | n/a | 0.12x |
| REDUCE | w32-o64 | 2.58 M | none | n/a | 10.89 M | 1.06 G | n/a | 0.00x |
| MODMUL | w32-o64 | 342.94 k | none | n/a | 714.80 k | 143.47 M | n/a | 0.00x |
| MODEXP | w32-o64 | 1.04 k | none | n/a | 2.37 k | 24.53 k | n/a | 0.04x |
| EXPONENTIATION | w32-o64 | 788.9 | none | n/a | 16.38 k | n/a | n/a | n/a |
| DIVIDE | w32-o64 | 543.52 k | none | n/a | 15.05 M | 561.85 M | n/a | 0.00x |
| ISQRT | w32-o64 | 7.75 k | none | n/a | 2.63 M | n/a | n/a | n/a |
| MODMUL_R2 | w32-o64 | 38.20 M | none | n/a | 1.59 M | 417.67 M | n/a | 0.09x |

### p2048 (2048-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w32-o64 | 33.41 M | none | n/a | 23.97 M | 604.92 M | n/a | 0.06x |
| SUBTRACT | w32-o64 | 33.71 M | none | n/a | 26.53 M | 603.05 M | n/a | 0.06x |
| ADDMOD | w32-o64 | 22.84 M | none | n/a | 8.73 M | 582.59 M | n/a | 0.04x |
| SUBTRACTMOD | w32-o64 | 23.06 M | none | n/a | 11.32 M | 579.99 M | n/a | 0.04x |
| MULTIPLYOPERANDSCANNING | w32-o64 | 4.38 M | none | n/a | 1.46 M | n/a | n/a | n/a |
| MULTIPLYPRODUCTSCANNING | w32 | 1.36 M | none | n/a | 1.46 M | 322.16 M | n/a | 0.00x |
| MONTGOMERYMULTIPLICATION | w32 | 11.57 M | none | n/a | 226.64 k | 200.04 M | n/a | 0.06x |
| COMPARE | w32-o64 | 37.15 M | none | n/a | 61.28 M | 603.98 M | n/a | 0.06x |
| REDUCE | w32-o64 | 408.59 k | none | n/a | 15.57 M | 593.66 M | n/a | 0.00x |
| MODMUL | w32-o64 | 42.32 k | none | n/a | 419.99 k | 42.35 M | n/a | 0.00x |
| MODEXP | w32-opt | 25.3 | none | n/a | 326.1 | 24.79 k | n/a | 0.00x |
| EXPONENTIATION | w32-opt | 34.7 | none | n/a | 2.84 k | n/a | n/a | n/a |
| DIVIDE | w32-opt | 21.73 k | none | n/a | 9.53 M | 502.90 M | n/a | 0.00x |
| ISQRT | w32-opt | 1.60 k | none | n/a | 1.27 M | n/a | n/a | n/a |
| MODMUL_R2 | w32-opt | 3.86 M | none | n/a | 407.25 k | 115.07 M | n/a | 0.03x |

## 6. Raw data

Also written to `Tesla_V100-PCIE-32GB_Report.csv` for analysis.

```csv
kind,device,device_type,kernel,modulus,bits,operation,items,seconds,ops_per_sec,mismatches
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,secp256k1,256,ADD,50000,0.005703559,8766456.128,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,secp256k1,256,ADD,50000,0.004666208,10715338.882,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,secp256k1,256,ADD,50000,0.004762356,10499005.126,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,secp256k1,256,ADD,50000,0.000025856,1933787128.713,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,secp256k1,256,ADD,50000,0.000060826,822016682.759,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,secp256k1,256,ADD,50000,0.001017086,49160052.569,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,secp256k1,256,ADD,50000,0.000036252,1379235616.984,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,secp256k1,256,ADD,50000,0.001463981,34153448.805,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32,secp256k1,256,ADD,50000,0.000022349,2237230134.860,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32,secp256k1,256,ADD,50000,0.001531616,32645258.987,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,secp256k1,256,ADD,50000,0.000023697,2109968409.676,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,secp256k1,256,ADD,50000,0.001340551,37298094.161,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,secp256k1,256,ADD,50000,0.000022842,2188953369.587,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,secp256k1,256,ADD,50000,0.001435725,34825613.132,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,secp256k1,256,SUBTRACT,50000,0.005437951,9194639.653,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,secp256k1,256,SUBTRACT,50000,0.005126694,9752873.909,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,secp256k1,256,SUBTRACT,50000,0.005073479,9855170.334,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,secp256k1,256,SUBTRACT,50000,0.000022592,2213172804.533,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,secp256k1,256,SUBTRACT,50000,0.000059842,835533353.695,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,secp256k1,256,SUBTRACT,50000,0.001013263,49345529.634,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,secp256k1,256,SUBTRACT,50000,0.000035384,1413069195.186,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,secp256k1,256,SUBTRACT,50000,0.001390731,35952316.703,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32,secp256k1,256,SUBTRACT,50000,0.000021507,2324820993.375,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32,secp256k1,256,SUBTRACT,50000,0.001624495,30778795.420,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,secp256k1,256,SUBTRACT,50000,0.000022937,2179887678.339,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,secp256k1,256,SUBTRACT,50000,0.001332831,37514133.817,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,secp256k1,256,SUBTRACT,50000,0.000021379,2338733587.084,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,secp256k1,256,SUBTRACT,50000,0.001364472,36644211.821,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,secp256k1,256,ADDMOD,50000,0.022078488,2264647.831,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,secp256k1,256,ADDMOD,50000,0.004847938,10313663.333,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,secp256k1,256,ADDMOD,50000,0.005779024,8651979.930,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,secp256k1,256,ADDMOD,50000,0.000023040,2170138888.889,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,secp256k1,256,ADDMOD,50000,0.000086900,575374236.702,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,secp256k1,256,ADDMOD,50000,0.001039366,48106249.156,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,secp256k1,256,ADDMOD,50000,0.000050017,999661880.356,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,secp256k1,256,ADDMOD,50000,0.001396148,35812822.139,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32,secp256k1,256,ADDMOD,50000,0.000028360,1763043251.742,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32,secp256k1,256,ADDMOD,50000,0.001626734,30736434.080,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,secp256k1,256,ADDMOD,50000,0.000024579,2034256908.474,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,secp256k1,256,ADDMOD,50000,0.001384639,36110493.537,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,secp256k1,256,ADDMOD,50000,0.000022729,2199839836.099,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,secp256k1,256,ADDMOD,50000,0.001385298,36093317.744,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,secp256k1,256,SUBTRACTMOD,50000,0.010463644,4778450.034,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,secp256k1,256,SUBTRACTMOD,50000,0.004878091,10249911.317,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,secp256k1,256,SUBTRACTMOD,50000,0.005692881,8782899.191,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,secp256k1,256,SUBTRACTMOD,50000,0.000022336,2238538681.948,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.000088443,565335929.379,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.001085969,46041830.807,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,secp256k1,256,SUBTRACTMOD,50000,0.000050466,990767521.032,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,secp256k1,256,SUBTRACTMOD,50000,0.001413287,35378517.984,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32,secp256k1,256,SUBTRACTMOD,50000,0.000029433,1698770427.326,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32,secp256k1,256,SUBTRACTMOD,50000,0.001343821,37207337.984,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,secp256k1,256,SUBTRACTMOD,50000,0.000024782,2017591130.987,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,secp256k1,256,SUBTRACTMOD,50000,0.001352623,36965215.304,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,secp256k1,256,SUBTRACTMOD,50000,0.000022974,2176375026.476,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,secp256k1,256,SUBTRACTMOD,50000,0.001420317,35203409.191,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.012028825,4156681.969,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.006797946,7355162.868,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.005393819,9269869.809,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.002030121,24629073.980,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.003170425,15770756.506,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000685039,72988547.695,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.002802309,17842429.238,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000152632,327585277.389,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.002218279,22539996.707,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000043835,1140640382.429,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001908696,26195894.673,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000036603,1366005539.123,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.002025349,24687102.509,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.012771380,3915003.698,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.005476250,9130335.503,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.005715016,8748881.916,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000022976,2176183844.011,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000362204,138043739.661,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001427320,35030686.716,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000116812,428037976.177,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.002102367,23782717.314,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000042829,1167434260.583,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001975381,25311572.185,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000045620,1096010007.324,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001932058,25879140.037,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000043538,1148423825.363,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.002032767,24597014.203,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.028883631,1731084.295,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.005242678,9537110.637,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.005145157,9717876.463,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000022112,2261215629.522,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000181509,275468423.093,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.001184213,42222134.909,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000061402,814305975.857,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.001586405,31517803.262,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000025717,1944233480.001,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.001342824,37234959.376,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000038012,1315376484.136,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.001328921,37624507.152,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000030571,1635536264.247,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.001406741,35543145.992,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,secp256k1,256,COMPARE,50000,0.004300209,11627341.847,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,secp256k1,256,COMPARE,50000,0.004917376,10168024.622,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,secp256k1,256,COMPARE,50000,0.007166254,6977145.955,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,secp256k1,256,COMPARE,50000,0.000022208,2251440922.190,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,secp256k1,256,COMPARE,50000,0.000059153,845265502.179,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,secp256k1,256,COMPARE,50000,0.000998883,50055913.242,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,secp256k1,256,COMPARE,50000,0.000035552,1406387666.918,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,secp256k1,256,COMPARE,50000,0.001418352,35252179.730,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,secp256k1,256,COMPARE,50000,0.000023020,2172016575.217,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,secp256k1,256,COMPARE,50000,0.001370617,36479919.168,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,secp256k1,256,COMPARE,50000,0.000022090,2263463518.698,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,secp256k1,256,COMPARE,50000,0.001435926,34820734.247,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,secp256k1,256,REDUCE,6250,0.000439259,14228507.852,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,secp256k1,256,REDUCE,6250,0.004888413,1278533.542,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,secp256k1,256,REDUCE,6250,0.005052318,1237055.946,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,secp256k1,256,REDUCE,50000,0.000022368,2235336194.564,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,secp256k1,256,REDUCE,50000,0.000309131,161743731.055,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,secp256k1,256,REDUCE,50000,0.001243867,40197224.302,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,secp256k1,256,REDUCE,50000,0.000208374,239953165.255,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,secp256k1,256,REDUCE,50000,0.001680193,29758486.021,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,secp256k1,256,REDUCE,50000,0.000081488,613586917.284,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,secp256k1,256,REDUCE,50000,0.001459592,34256147.777,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,secp256k1,256,REDUCE,50000,0.000071686,697486317.472,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,secp256k1,256,REDUCE,50000,0.001499927,33334955.176,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,secp256k1,256,MODMUL,3125,0.001284155,2433506.849,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,secp256k1,256,MODMUL,3125,0.004524139,690739.168,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,secp256k1,256,MODMUL,3125,0.004792895,652006.773,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,secp256k1,256,MODMUL,50000,0.000041856,1194571865.443,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,secp256k1,256,MODMUL,50000,0.000842688,59333939.286,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,secp256k1,256,MODMUL,50000,0.001793825,27873398.876,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,secp256k1,256,MODMUL,50000,0.000484036,103298091.364,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,secp256k1,256,MODMUL,50000,0.002094541,23871578.799,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,secp256k1,256,MODMUL,50000,0.000230794,216643453.409,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,secp256k1,256,MODMUL,50000,0.001621381,30837910.801,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,secp256k1,256,MODMUL,50000,0.000161458,309677982.629,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,secp256k1,256,MODMUL,50000,0.001474099,33919024.414,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,secp256k1,256,MODEXP,781,0.020525776,38049.719,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,secp256k1,256,MODEXP,781,0.005533636,141136.858,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,secp256k1,256,MODEXP,781,0.005978943,130625.095,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,secp256k1,256,MODEXP,50000,0.079579487,628302.618,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,secp256k1,256,MODEXP,50000,0.037811130,1322361.961,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,secp256k1,256,MODEXP,50000,0.038833820,1287537.513,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,secp256k1,256,MODEXP,50000,0.006862134,7286363.133,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,secp256k1,256,MODEXP,50000,0.008575583,5830507.388,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,secp256k1,256,MODEXP,50000,0.007160831,6982429.823,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,secp256k1,256,MODEXP,50000,0.008832452,5660942.169,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,secp256k1,256,MODEXP,50000,0.002533626,19734562.807,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,secp256k1,256,MODEXP,50000,0.004030475,12405485.571,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,secp256k1,256,EXPONENTIATION,781,0.006131917,127366.368,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,secp256k1,256,EXPONENTIATION,781,0.005164757,151217.182,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,secp256k1,256,EXPONENTIATION,781,0.008400294,92972.936,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,secp256k1,256,EXPONENTIATION,50000,0.055685965,897892.315,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,secp256k1,256,EXPONENTIATION,50000,0.056616902,883128.505,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,secp256k1,256,EXPONENTIATION,50000,0.019909517,2511361.774,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,secp256k1,256,EXPONENTIATION,50000,0.021583802,2316552.020,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,secp256k1,256,EXPONENTIATION,50000,0.000971660,51458328.094,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,secp256k1,256,EXPONENTIATION,50000,0.002560829,19524927.822,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,secp256k1,256,EXPONENTIATION,50000,0.000792612,63082569.705,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,secp256k1,256,EXPONENTIATION,50000,0.002155950,23191631.812,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,secp256k1,256,DIVIDE,6250,0.001104825,5657004.494,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,secp256k1,256,DIVIDE,6250,0.005072181,1232211.548,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,secp256k1,256,DIVIDE,6250,0.004963051,1259306.019,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,secp256k1,256,DIVIDE,50000,0.000025056,1995530012.771,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,secp256k1,256,DIVIDE,50000,0.000583005,85762563.208,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,secp256k1,256,DIVIDE,50000,0.001689087,29601790.484,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,secp256k1,256,DIVIDE,50000,0.000427978,116828435.959,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,secp256k1,256,DIVIDE,50000,0.002401897,20816879.453,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,secp256k1,256,DIVIDE,50000,0.000164778,303438617.751,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,secp256k1,256,DIVIDE,50000,0.002054963,24331337.994,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,secp256k1,256,DIVIDE,50000,0.000181573,275371372.443,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,secp256k1,256,DIVIDE,50000,0.002294422,21791980.622,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,secp256k1,256,ISQRT,1562,0.000366234,4265032.881,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,secp256k1,256,ISQRT,1562,0.004957689,315066.151,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,secp256k1,256,ISQRT,50000,0.007250841,6895751.809,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,secp256k1,256,ISQRT,50000,0.008213143,6087803.386,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,secp256k1,256,ISQRT,50000,0.005728113,8728878.080,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,secp256k1,256,ISQRT,50000,0.007487819,6677511.883,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,secp256k1,256,ISQRT,50000,0.001795925,27840804.974,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,secp256k1,256,ISQRT,50000,0.003407417,14673871.812,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,secp256k1,256,ISQRT,50000,0.001238999,40355159.073,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,secp256k1,256,ISQRT,50000,0.003527314,14175091.656,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,secp256k1,256,MODMUL_R2,50000,0.030192364,1656047.867,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,secp256k1,256,MODMUL_R2,50000,0.005708636,8758659.660,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,secp256k1,256,MODMUL_R2,50000,0.007072396,7069739.879,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,secp256k1,256,MODMUL_R2,50000,0.000024832,2013530927.835,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,secp256k1,256,MODMUL_R2,50000,0.000207134,241389633.431,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,secp256k1,256,MODMUL_R2,50000,0.001201690,41608068.636,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,secp256k1,256,MODMUL_R2,50000,0.000067367,742203314.451,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,secp256k1,256,MODMUL_R2,50000,0.001340901,37288363.450,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,secp256k1,256,MODMUL_R2,50000,0.000046629,1072294226.794,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,secp256k1,256,MODMUL_R2,50000,0.001362999,36683813.492,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,secp256k1,256,MODMUL_R2,50000,0.000039588,1263010220.609,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,secp256k1,256,MODMUL_R2,50000,0.002228438,22437240.398,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,rsa256(composite),256,ADD,50000,0.013149819,3802333.709,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,rsa256(composite),256,ADD,50000,0.005625595,8887948.666,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,rsa256(composite),256,ADD,50000,0.005260025,9505658.286,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,rsa256(composite),256,ADD,50000,0.000022528,2219460227.273,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,rsa256(composite),256,ADD,50000,0.000058511,854539598.811,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,rsa256(composite),256,ADD,50000,0.000991533,50426965.799,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,rsa256(composite),256,ADD,50000,0.000034986,1429145233.540,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,rsa256(composite),256,ADD,50000,0.001379365,36248564.044,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32,rsa256(composite),256,ADD,50000,0.000021274,2350290188.353,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32,rsa256(composite),256,ADD,50000,0.001360595,36748628.534,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,rsa256(composite),256,ADD,50000,0.000023145,2160294192.562,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,rsa256(composite),256,ADD,50000,0.001295540,38593944.901,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,rsa256(composite),256,ADD,50000,0.000023932,2089256081.022,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,rsa256(composite),256,ADD,50000,0.002219120,22531453.498,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACT,50000,0.011781581,4243912.582,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACT,50000,0.004909487,10184363.388,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACT,50000,0.004936478,10128678.819,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,rsa256(composite),256,SUBTRACT,50000,0.000023552,2122961956.522,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,rsa256(composite),256,SUBTRACT,50000,0.000059406,841666079.619,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,rsa256(composite),256,SUBTRACT,50000,0.000999672,50016405.925,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,rsa256(composite),256,SUBTRACT,50000,0.000035919,1392024818.744,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,rsa256(composite),256,SUBTRACT,50000,0.001358843,36796007.941,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32,rsa256(composite),256,SUBTRACT,50000,0.000021447,2331319876.892,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32,rsa256(composite),256,SUBTRACT,50000,0.001408124,35508236.659,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,rsa256(composite),256,SUBTRACT,50000,0.000023562,2122050868.343,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,rsa256(composite),256,SUBTRACT,50000,0.001313823,38056875.534,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,rsa256(composite),256,SUBTRACT,50000,0.000023137,2161044200.357,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,rsa256(composite),256,SUBTRACT,50000,0.002226734,22454410.122,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,rsa256(composite),256,ADDMOD,50000,0.018470196,2707063.856,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,rsa256(composite),256,ADDMOD,50000,0.005290615,9450697.193,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,rsa256(composite),256,ADDMOD,50000,0.008106881,6167600.095,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,rsa256(composite),256,ADDMOD,50000,0.000022656,2206920903.955,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,rsa256(composite),256,ADDMOD,50000,0.000081196,615793730.613,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,rsa256(composite),256,ADDMOD,50000,0.001067486,46839021.357,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,rsa256(composite),256,ADDMOD,50000,0.000046615,1072618256.376,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,rsa256(composite),256,ADDMOD,50000,0.001383091,36150911.940,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32,rsa256(composite),256,ADDMOD,50000,0.000026807,1865179417.295,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32,rsa256(composite),256,ADDMOD,50000,0.001384400,36116727.596,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,rsa256(composite),256,ADDMOD,50000,0.000027168,1840395289.923,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,rsa256(composite),256,ADDMOD,50000,0.001356794,36851580.318,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,rsa256(composite),256,ADDMOD,50000,0.000024869,2010535989.102,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,rsa256(composite),256,ADDMOD,50000,0.002197309,22755105.863,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACTMOD,50000,0.016536601,3023595.961,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACTMOD,50000,0.007175369,6968282.733,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACTMOD,50000,0.005503554,9085038.482,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,rsa256(composite),256,SUBTRACTMOD,50000,0.000023360,2140410958.904,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,rsa256(composite),256,SUBTRACTMOD,50000,0.000087141,573782712.510,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,rsa256(composite),256,SUBTRACTMOD,50000,0.001039109,48118146.533,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,rsa256(composite),256,SUBTRACTMOD,50000,0.000051467,971496658.207,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,rsa256(composite),256,SUBTRACTMOD,50000,0.001377327,36302199.355,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32,rsa256(composite),256,SUBTRACTMOD,50000,0.000028919,1728962371.535,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32,rsa256(composite),256,SUBTRACTMOD,50000,0.001357692,36827202.159,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,rsa256(composite),256,SUBTRACTMOD,50000,0.000024549,2036745779.944,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,rsa256(composite),256,SUBTRACTMOD,50000,0.001480285,33777277.589,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,rsa256(composite),256,SUBTRACTMOD,50000,0.000025260,1979420915.195,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,rsa256(composite),256,SUBTRACTMOD,50000,0.002225115,22470749.067,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.010013419,4993299.477,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.004742288,10543433.948,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.005706123,8762517.009,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.002010495,24869497.381,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.003111137,16071295.038,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000683420,73161453.926,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.002847145,17561451.425,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000151203,330681239.029,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.002133484,23435844.860,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000043935,1138044164.399,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001923590,25993065.494,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000040366,1238667274.998,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.002731456,18305255.029,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.010612077,4711612.995,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.004973964,10052344.551,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.005485995,9114116.990,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000022112,2261215629.522,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000358372,139519830.353,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001442790,34655078.652,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000116721,428372252.437,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.002083962,23992760.310,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000044369,1126912366.744,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001990534,25118887.494,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000044974,1111752647.001,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001884628,26530433.425,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000048613,1028530207.408,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.002726450,18338865.690,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.036383080,1374265.180,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.005483836,9117705.208,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.005107087,9790316.820,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000022976,2176183844.011,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000177245,282095396.654,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.001137289,43964197.373,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000059646,838279571.001,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.001537267,32525254.185,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000025488,1961700776.008,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.001377844,36288580.957,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000037253,1342177280.000,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.001352753,36961661.954,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000033478,1493517249.821,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.002222230,22499921.661,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,rsa256(composite),256,COMPARE,50000,0.007485904,6679220.061,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,rsa256(composite),256,COMPARE,50000,0.004693895,10652134.323,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,rsa256(composite),256,COMPARE,50000,0.004888903,10227243.207,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,rsa256(composite),256,COMPARE,50000,0.000022912,2182262569.832,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,rsa256(composite),256,COMPARE,50000,0.000059203,844551626.389,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,rsa256(composite),256,COMPARE,50000,0.000980333,51003079.636,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,rsa256(composite),256,COMPARE,50000,0.000036456,1371519219.298,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,rsa256(composite),256,COMPARE,50000,0.001782549,28049718.560,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,rsa256(composite),256,COMPARE,50000,0.000023191,2156010670.201,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,rsa256(composite),256,COMPARE,50000,0.001319298,37898941.645,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,rsa256(composite),256,COMPARE,50000,0.000023277,2148042138.956,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,rsa256(composite),256,COMPARE,50000,0.002214620,22577236.663,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,rsa256(composite),256,REDUCE,6250,0.000843577,7408926.771,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,rsa256(composite),256,REDUCE,6250,0.004488951,1392307.467,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,rsa256(composite),256,REDUCE,6250,0.004840934,1291073.167,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,rsa256(composite),256,REDUCE,50000,0.000022624,2210042432.815,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,rsa256(composite),256,REDUCE,50000,0.000308016,162329216.968,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,rsa256(composite),256,REDUCE,50000,0.001236756,40428346.973,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,rsa256(composite),256,REDUCE,50000,0.000208809,239453233.842,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,rsa256(composite),256,REDUCE,50000,0.001816902,27519371.485,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,rsa256(composite),256,REDUCE,50000,0.000081214,615657353.040,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,rsa256(composite),256,REDUCE,50000,0.001391504,35932344.744,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,rsa256(composite),256,REDUCE,50000,0.000080042,624671633.438,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,rsa256(composite),256,REDUCE,50000,0.002266517,22060281.042,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,rsa256(composite),256,MODMUL,3125,0.001129248,2767328.334,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,rsa256(composite),256,MODMUL,3125,0.004979243,627605.444,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,rsa256(composite),256,MODMUL,3125,0.005115465,610892.653,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,rsa256(composite),256,MODMUL,50000,0.000042208,1184609552.691,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,rsa256(composite),256,MODMUL,50000,0.000833149,60013274.227,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,rsa256(composite),256,MODMUL,50000,0.001844229,27111600.536,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,rsa256(composite),256,MODMUL,50000,0.000482456,103636380.869,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,rsa256(composite),256,MODMUL,50000,0.002141118,23352285.507,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,rsa256(composite),256,MODMUL,50000,0.000236454,211457660.684,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,rsa256(composite),256,MODMUL,50000,0.001624177,30784822.530,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,rsa256(composite),256,MODMUL,50000,0.000180219,277440121.958,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,rsa256(composite),256,MODMUL,50000,0.002408021,20763938.406,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,rsa256(composite),256,MODEXP,781,0.012588489,62040.806,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,rsa256(composite),256,MODEXP,781,0.005373434,145344.672,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,rsa256(composite),256,MODEXP,781,0.005887809,132646.965,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,rsa256(composite),256,MODEXP,50000,0.077048413,648942.633,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,rsa256(composite),256,MODEXP,50000,0.037919708,1318575.556,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,rsa256(composite),256,MODEXP,50000,0.038914490,1284868.438,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,rsa256(composite),256,MODEXP,50000,0.006844256,7305395.987,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,rsa256(composite),256,MODEXP,50000,0.008605924,5809951.398,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,rsa256(composite),256,MODEXP,50000,0.007190044,6954060.440,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,rsa256(composite),256,MODEXP,50000,0.008848660,5650573.053,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,rsa256(composite),256,MODEXP,50000,0.002827498,17683478.379,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,rsa256(composite),256,MODEXP,50000,0.005085019,9832804.863,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,rsa256(composite),256,EXPONENTIATION,781,0.004682928,166776.001,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,rsa256(composite),256,EXPONENTIATION,781,0.004958010,157522.878,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,rsa256(composite),256,EXPONENTIATION,781,0.008498213,91901.674,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,rsa256(composite),256,EXPONENTIATION,50000,0.055761942,896668.914,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,rsa256(composite),256,EXPONENTIATION,50000,0.056716760,881573.630,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,rsa256(composite),256,EXPONENTIATION,50000,0.019948336,2506474.718,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,rsa256(composite),256,EXPONENTIATION,50000,0.021700888,2304053.174,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,rsa256(composite),256,EXPONENTIATION,50000,0.000973094,51382496.321,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,rsa256(composite),256,EXPONENTIATION,50000,0.002606713,19181244.878,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,rsa256(composite),256,EXPONENTIATION,50000,0.000794521,62931002.568,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,rsa256(composite),256,EXPONENTIATION,50000,0.003005343,16637035.974,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,rsa256(composite),256,DIVIDE,6250,0.000855507,7305609.594,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,rsa256(composite),256,DIVIDE,6250,0.005074800,1231575.631,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,rsa256(composite),256,DIVIDE,6250,0.005188707,1204539.018,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,rsa256(composite),256,DIVIDE,50000,0.000025216,1982868020.305,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,rsa256(composite),256,DIVIDE,50000,0.000586743,85216184.810,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,rsa256(composite),256,DIVIDE,50000,0.001663270,30061264.876,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,rsa256(composite),256,DIVIDE,50000,0.000439094,113870826.178,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,rsa256(composite),256,DIVIDE,50000,0.002572023,19439950.874,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,rsa256(composite),256,DIVIDE,50000,0.000164563,303834877.821,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,rsa256(composite),256,DIVIDE,50000,0.002115092,23639633.559,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,rsa256(composite),256,DIVIDE,50000,0.000195237,256099006.781,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,rsa256(composite),256,DIVIDE,50000,0.002900677,17237355.074,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,rsa256(composite),256,ISQRT,1562,0.000367518,4250132.076,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,rsa256(composite),256,ISQRT,1562,0.005066938,308272.965,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,rsa256(composite),256,ISQRT,50000,0.007276436,6871495.859,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,rsa256(composite),256,ISQRT,50000,0.008278604,6039665.656,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,rsa256(composite),256,ISQRT,50000,0.005929988,8431720.275,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,rsa256(composite),256,ISQRT,50000,0.007804458,6406594.771,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,rsa256(composite),256,ISQRT,50000,0.001798293,27804144.021,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,rsa256(composite),256,ISQRT,50000,0.003397586,14716331.205,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,rsa256(composite),256,ISQRT,50000,0.001318193,37930711.552,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,rsa256(composite),256,ISQRT,50000,0.003579566,13968173.708,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,rsa256(composite),256,MODMUL_R2,50000,0.024991097,2000712.495,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,rsa256(composite),256,MODMUL_R2,50000,0.005684553,8795766.343,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,rsa256(composite),256,MODMUL_R2,50000,0.007370766,6783555.503,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,rsa256(composite),256,MODMUL_R2,50000,0.000025600,1953125000.000,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,rsa256(composite),256,MODMUL_R2,50000,0.000208289,240051067.607,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,rsa256(composite),256,MODMUL_R2,50000,0.001226305,40772891.849,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,rsa256(composite),256,MODMUL_R2,50000,0.000069209,722449128.933,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,rsa256(composite),256,MODMUL_R2,50000,0.001418478,35249049.323,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,rsa256(composite),256,MODMUL_R2,50000,0.000045979,1087451873.495,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,rsa256(composite),256,MODMUL_R2,50000,0.001331667,37546925.675,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,rsa256(composite),256,MODMUL_R2,50000,0.000040641,1230283037.719,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,rsa256(composite),256,MODMUL_R2,50000,0.002177057,22966784.034,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,ADD,25000,0.007636513,3273745.500,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,ADD,25000,0.004596788,5438580.191,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,ADD,25000,0.007491006,3337335.470,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,brainpoolP512r1,512,ADD,50000,0.000027904,1791857798.165,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,brainpoolP512r1,512,ADD,50000,0.000128127,390237828.462,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,brainpoolP512r1,512,ADD,50000,0.001852112,26996207.466,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,brainpoolP512r1,512,ADD,50000,0.000074549,670699843.684,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,brainpoolP512r1,512,ADD,50000,0.003762384,13289446.459,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32,brainpoolP512r1,512,ADD,50000,0.000053954,926716011.919,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32,brainpoolP512r1,512,ADD,50000,0.003248310,15392619.235,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,brainpoolP512r1,512,ADD,50000,0.000049709,1005854208.812,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,brainpoolP512r1,512,ADD,50000,0.003267752,15301038.726,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,brainpoolP512r1,512,ADD,50000,0.000054780,912741133.593,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,brainpoolP512r1,512,ADD,50000,0.004054749,12331219.408,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACT,25000,0.006381351,3917665.725,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACT,25000,0.006006807,4161944.959,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACT,25000,0.006855190,3646871.914,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,brainpoolP512r1,512,SUBTRACT,50000,0.000027744,1802191464.821,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,brainpoolP512r1,512,SUBTRACT,50000,0.000129105,387281671.879,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,brainpoolP512r1,512,SUBTRACT,50000,0.001808490,27647373.991,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,brainpoolP512r1,512,SUBTRACT,50000,0.000077115,648382096.366,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,brainpoolP512r1,512,SUBTRACT,50000,0.003342504,14958845.086,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32,brainpoolP512r1,512,SUBTRACT,50000,0.000054297,920862529.079,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32,brainpoolP512r1,512,SUBTRACT,50000,0.003278607,15250379.029,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,brainpoolP512r1,512,SUBTRACT,50000,0.000049385,1012452982.693,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,brainpoolP512r1,512,SUBTRACT,50000,0.003275917,15262901.829,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,brainpoolP512r1,512,SUBTRACT,50000,0.000052980,943751809.183,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,brainpoolP512r1,512,SUBTRACT,50000,0.004095848,12207484.060,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,ADDMOD,25000,0.009921290,2519833.601,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,ADDMOD,25000,0.005170484,4835137.277,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,ADDMOD,25000,0.005578901,4481169.337,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,brainpoolP512r1,512,ADDMOD,50000,0.000026880,1860119047.619,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,brainpoolP512r1,512,ADDMOD,50000,0.000181862,274933773.487,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,brainpoolP512r1,512,ADDMOD,50000,0.001844586,27106353.557,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,brainpoolP512r1,512,ADDMOD,50000,0.000098814,506001052.773,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,brainpoolP512r1,512,ADDMOD,50000,0.003349174,14929054.459,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32,brainpoolP512r1,512,ADDMOD,50000,0.000062942,794381634.571,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32,brainpoolP512r1,512,ADDMOD,50000,0.003306669,15120957.299,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,brainpoolP512r1,512,ADDMOD,50000,0.000057999,862084895.636,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,brainpoolP512r1,512,ADDMOD,50000,0.003271186,15284976.032,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,brainpoolP512r1,512,ADDMOD,50000,0.000062182,804090181.600,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,brainpoolP512r1,512,ADDMOD,50000,0.004006234,12480549.072,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACTMOD,25000,0.008494392,2943118.223,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACTMOD,25000,0.004731015,5284278.308,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACTMOD,25000,0.005059118,4941572.787,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000027040,1849112426.036,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000198643,251707819.558,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,brainpoolP512r1,512,SUBTRACTMOD,50000,0.001925230,25970923.386,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000107460,465289597.150,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,brainpoolP512r1,512,SUBTRACTMOD,50000,0.003317221,15072857.364,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000068119,734009292.851,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32,brainpoolP512r1,512,SUBTRACTMOD,50000,0.003292788,15184700.917,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000056099,891282854.451,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,brainpoolP512r1,512,SUBTRACTMOD,50000,0.003254108,15365193.711,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000062869,795303929.887,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,brainpoolP512r1,512,SUBTRACTMOD,50000,0.004109986,12165491.876,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.009560021,2615057.004,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.004677164,5345119.371,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.004689995,5330496.126,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.011197037,4465467.091,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.013163803,3798294.468,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.002945901,16972735.951,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.007245647,6900694.975,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.002333952,21422891.324,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.006403347,7808416.570,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000150658,331877334.730,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.004359715,11468639.483,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000164074,304740528.386,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.005075821,9850623.042,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.010244687,2440289.300,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.004681622,5340029.557,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.004849936,5154707.198,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000026848,1862336114.422,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.002421686,20646772.381,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.004432714,11279771.212,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000674284,74152736.740,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.004918836,10165006.640,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000254210,196687802.378,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.004419918,11312427.219,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000249761,200191350.642,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.004463702,11201464.386,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000275094,181756013.271,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.005142117,9723621.546,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.025401706,984185.866,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.005147106,4857098.325,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.004809166,5198406.528,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000027712,1804272517.321,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000642646,77803326.294,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.002388186,20936392.761,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000191405,261226105.293,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.003462309,14441229.753,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000058968,847917954.873,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.003284248,15224184.995,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000084578,591170677.711,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.003336797,14984429.667,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000070749,706724098.406,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.004006788,12478823.380,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,COMPARE,25000,0.001901803,13145420.326,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,COMPARE,25000,0.004510579,5542525.705,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,COMPARE,25000,0.004783276,5226543.479,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,brainpoolP512r1,512,COMPARE,50000,0.000027616,1810544611.819,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,brainpoolP512r1,512,COMPARE,50000,0.000127331,392677338.916,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,brainpoolP512r1,512,COMPARE,50000,0.001791657,27907127.206,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,brainpoolP512r1,512,COMPARE,50000,0.000072847,686370025.138,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,brainpoolP512r1,512,COMPARE,50000,0.003312507,15094307.442,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,brainpoolP512r1,512,COMPARE,50000,0.000049421,1011716042.711,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,brainpoolP512r1,512,COMPARE,50000,0.003286861,15212082.289,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,brainpoolP512r1,512,COMPARE,50000,0.000051464,971553795.762,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,brainpoolP512r1,512,COMPARE,50000,0.003970039,12594334.840,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,REDUCE,3125,0.000344725,9065197.515,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,REDUCE,3125,0.004456657,701198.232,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,REDUCE,3125,0.004493234,695490.147,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,brainpoolP512r1,512,REDUCE,50000,0.000029152,1715148188.804,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,brainpoolP512r1,512,REDUCE,50000,0.000891145,56107593.839,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,brainpoolP512r1,512,REDUCE,50000,0.002569787,19456865.267,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,brainpoolP512r1,512,REDUCE,50000,0.000715938,69838456.420,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,brainpoolP512r1,512,REDUCE,50000,0.003980358,12561684.103,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,brainpoolP512r1,512,REDUCE,50000,0.000237898,210174134.927,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,brainpoolP512r1,512,REDUCE,50000,0.003456892,14463859.527,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,brainpoolP512r1,512,REDUCE,50000,0.000228462,218854736.216,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,brainpoolP512r1,512,REDUCE,50000,0.004167601,11997309.722,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL,1562,0.000466510,3348266.823,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL,1562,0.004516976,345806.574,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL,1562,0.005207484,299952.914,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,brainpoolP512r1,512,MODMUL,50000,0.000120672,414346327.234,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,brainpoolP512r1,512,MODMUL,50000,0.002522489,19821692.005,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,brainpoolP512r1,512,MODMUL,50000,0.004250374,11763670.728,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,brainpoolP512r1,512,MODMUL,50000,0.001795794,27842837.209,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,brainpoolP512r1,512,MODMUL,50000,0.005085793,9831308.559,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,brainpoolP512r1,512,MODMUL,50000,0.000748167,66829994.962,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,brainpoolP512r1,512,MODMUL,50000,0.004027489,12414683.222,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,brainpoolP512r1,512,MODMUL,50000,0.000676883,73868005.829,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,brainpoolP512r1,512,MODMUL,50000,0.004624117,10812875.350,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,MODEXP,390,0.028861927,13512.611,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,MODEXP,390,0.006522103,59796.664,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,MODEXP,390,0.007895099,49397.734,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,brainpoolP512r1,512,MODEXP,50000,0.808304071,61857.909,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,brainpoolP512r1,512,MODEXP,50000,0.569398497,87811.963,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,brainpoolP512r1,512,MODEXP,50000,0.567848307,88051.685,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,brainpoolP512r1,512,MODEXP,50000,0.047921408,1043375.020,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,brainpoolP512r1,512,MODEXP,50000,0.051232763,975937.995,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,brainpoolP512r1,512,MODEXP,50000,0.050397476,992113.175,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,brainpoolP512r1,512,MODEXP,50000,0.053690538,931262.787,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,brainpoolP512r1,512,MODEXP,50000,0.031640771,1580239.624,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,brainpoolP512r1,512,MODEXP,50000,0.035784731,1397243.980,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,EXPONENTIATION,390,0.012835272,30385.020,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,EXPONENTIATION,390,0.004987231,78199.706,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,EXPONENTIATION,390,0.011495298,33926.915,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,brainpoolP512r1,512,EXPONENTIATION,50000,0.493494340,101318.285,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,brainpoolP512r1,512,EXPONENTIATION,50000,0.499426267,100114.878,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,brainpoolP512r1,512,EXPONENTIATION,50000,0.122662150,407623.705,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,brainpoolP512r1,512,EXPONENTIATION,50000,0.125939734,397015.290,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,brainpoolP512r1,512,EXPONENTIATION,50000,0.042878188,1166094.050,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,brainpoolP512r1,512,EXPONENTIATION,50000,0.046292278,1080093.747,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,brainpoolP512r1,512,EXPONENTIATION,50000,0.038943678,1283905.438,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,brainpoolP512r1,512,EXPONENTIATION,50000,0.043093692,1160262.622,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,DIVIDE,3125,0.000464540,6727084.378,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,DIVIDE,3125,0.004896861,638163.921,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,DIVIDE,3125,0.004647418,672416.377,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,brainpoolP512r1,512,DIVIDE,50000,0.000061152,817634746.206,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,brainpoolP512r1,512,DIVIDE,50000,0.001792430,27895092.066,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,brainpoolP512r1,512,DIVIDE,50000,0.006230991,8024405.769,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,brainpoolP512r1,512,DIVIDE,50000,0.001829088,27336028.065,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,brainpoolP512r1,512,DIVIDE,50000,0.006068616,8239110.783,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,brainpoolP512r1,512,DIVIDE,50000,0.000549581,90978416.154,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,brainpoolP512r1,512,DIVIDE,50000,0.004849849,10309599.253,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,brainpoolP512r1,512,DIVIDE,50000,0.000522561,95682605.260,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,brainpoolP512r1,512,DIVIDE,50000,0.005439897,9191350.383,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,ISQRT,781,0.000304811,2562243.540,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,ISQRT,781,0.004708412,165873.335,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,brainpoolP512r1,512,ISQRT,50000,0.038731820,1290928.234,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,brainpoolP512r1,512,ISQRT,50000,0.042577888,1174318.463,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,brainpoolP512r1,512,ISQRT,50000,0.037876269,1320087.785,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,brainpoolP512r1,512,ISQRT,50000,0.041384589,1208179.208,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,brainpoolP512r1,512,ISQRT,50000,0.008939984,5592851.159,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,brainpoolP512r1,512,ISQRT,50000,0.012286563,4069486.298,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,brainpoolP512r1,512,ISQRT,50000,0.006976071,7167358.317,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,brainpoolP512r1,512,ISQRT,50000,0.010925437,4576475.991,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL_R2,25000,0.019552612,1278601.551,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL_R2,25000,0.005236043,4774597.905,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL_R2,25000,0.005964791,4191261.692,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,brainpoolP512r1,512,MODMUL_R2,50000,0.000040000,1250000000.000,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,brainpoolP512r1,512,MODMUL_R2,50000,0.000913215,54751620.106,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,brainpoolP512r1,512,MODMUL_R2,50000,0.004727581,10576233.402,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,brainpoolP512r1,512,MODMUL_R2,50000,0.000184935,270365357.719,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,brainpoolP512r1,512,MODMUL_R2,50000,0.003360564,14878454.901,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,brainpoolP512r1,512,MODMUL_R2,50000,0.000145524,343585839.196,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,brainpoolP512r1,512,MODMUL_R2,50000,0.003398641,14711762.676,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,brainpoolP512r1,512,MODMUL_R2,50000,0.000096782,516625023.125,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,brainpoolP512r1,512,MODMUL_R2,50000,0.004003750,12488292.452,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p1024,1024,ADD,12500,0.004072107,3069663.916,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p1024,1024,ADD,12500,0.004562174,2739921.783,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p1024,1024,ADD,12500,0.004536323,2755535.696,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,p1024,1024,ADD,50000,0.000043776,1142178362.573,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,p1024,1024,ADD,50000,0.000379909,131610463.035,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,p1024,1024,ADD,50000,0.008018534,6235553.788,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,p1024,1024,ADD,50000,0.000200419,249477359.597,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,p1024,1024,ADD,50000,0.006104611,8190529.951,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32,p1024,1024,ADD,50000,0.000395870,126304067.698,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32,p1024,1024,ADD,50000,0.006724310,7435707.168,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,p1024,1024,ADD,50000,0.000110179,453807180.246,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,p1024,1024,ADD,50000,0.005996742,8337860.718,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,p1024,1024,ADD,50000,0.000105355,474585692.486,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,p1024,1024,ADD,50000,0.007727768,6470173.502,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p1024,1024,SUBTRACT,12500,0.002676064,4671039.303,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p1024,1024,SUBTRACT,12500,0.004630746,2699349.107,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p1024,1024,SUBTRACT,12500,0.004550473,2746967.214,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,p1024,1024,SUBTRACT,50000,0.000043296,1154841093.865,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,p1024,1024,SUBTRACT,50000,0.000364466,137187012.268,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,p1024,1024,SUBTRACT,50000,0.007974662,6269858.206,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,p1024,1024,SUBTRACT,50000,0.000193911,257850375.762,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,p1024,1024,SUBTRACT,50000,0.006251990,7997453.636,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32,p1024,1024,SUBTRACT,50000,0.000391118,127838658.824,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32,p1024,1024,SUBTRACT,50000,0.006710598,7450900.785,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,p1024,1024,SUBTRACT,50000,0.000110779,451349255.137,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,p1024,1024,SUBTRACT,50000,0.006106125,8188499.309,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,p1024,1024,SUBTRACT,50000,0.000110770,451385780.241,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,p1024,1024,SUBTRACT,50000,0.007609926,6570366.083,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p1024,1024,ADDMOD,12500,0.005760539,2169935.844,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p1024,1024,ADDMOD,12500,0.004521399,2764631.012,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p1024,1024,ADDMOD,12500,0.004709269,2654339.751,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,p1024,1024,ADDMOD,50000,0.000043648,1145527859.238,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,p1024,1024,ADDMOD,50000,0.000555249,90049691.161,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,p1024,1024,ADDMOD,50000,0.008159092,6128132.895,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,p1024,1024,ADDMOD,50000,0.000284907,175495827.335,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,p1024,1024,ADDMOD,50000,0.006205977,8056749.089,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32,p1024,1024,ADDMOD,50000,0.000414306,120683726.877,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32,p1024,1024,ADDMOD,50000,0.006815925,7335761.587,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,p1024,1024,ADDMOD,50000,0.000140812,355083335.345,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,p1024,1024,ADDMOD,50000,0.006102548,8193298.956,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,p1024,1024,ADDMOD,50000,0.000141171,354180290.768,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,p1024,1024,ADDMOD,50000,0.007677170,6512816.587,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p1024,1024,SUBTRACTMOD,12500,0.005193588,2406813.952,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p1024,1024,SUBTRACTMOD,12500,0.004740972,2636590.126,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p1024,1024,SUBTRACTMOD,12500,0.004843936,2580546.063,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,p1024,1024,SUBTRACTMOD,50000,0.000043488,1149742457.689,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,p1024,1024,SUBTRACTMOD,50000,0.000551036,90738176.979,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,p1024,1024,SUBTRACTMOD,50000,0.008156810,6129847.319,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,p1024,1024,SUBTRACTMOD,50000,0.000281691,177499478.077,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,p1024,1024,SUBTRACTMOD,50000,0.006143962,8138071.277,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32,p1024,1024,SUBTRACTMOD,50000,0.000426026,117363744.192,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32,p1024,1024,SUBTRACTMOD,50000,0.006860877,7287698.019,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,p1024,1024,SUBTRACTMOD,50000,0.000143078,349459639.016,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,p1024,1024,SUBTRACTMOD,50000,0.006146600,8134578.460,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,p1024,1024,SUBTRACTMOD,50000,0.000141426,353541807.645,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,p1024,1024,SUBTRACTMOD,50000,0.007758495,6444548.902,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.012155588,1028333.634,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.004885382,2558653.571,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.005538617,2256881.107,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.051826248,964762.103,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.061404713,814269.745,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.013171139,3796178.920,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.020802468,2403560.968,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.017325163,2885975.733,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.025488732,1961651.139,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.000372099,134372844.653,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.008058435,6204678.661,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.000367928,135896130.142,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.009348038,5348715.954,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.011440041,1092653.426,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.007210393,1733608.688,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.008305130,1505093.839,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000048832,1023918741.809,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.018140263,2756299.622,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.027026520,1850034.707,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.004657931,10734379.874,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.012189724,4101815.589,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.001752685,28527658.064,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.010142956,4929529.397,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.001234063,40516571.791,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.008982936,5566108.875,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.001257757,39753308.213,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.010389743,4812438.564,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.042097692,296928.392,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.010009377,1248828.977,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.009316838,1341656.902,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000065760,760340632.603,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.003278700,15249946.651,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.010882081,4594709.413,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000695404,71900649.053,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.006529896,7657089.864,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000190535,262418847.866,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.006573115,7606743.519,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000280828,178044935.485,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.006275238,7967825.280,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000227264,220008549.203,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.007783056,6424211.772,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p1024,1024,COMPARE,12500,0.001678968,7445049.476,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p1024,1024,COMPARE,12500,0.010625365,1176430.173,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p1024,1024,COMPARE,12500,0.010637735,1175062.171,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,p1024,1024,COMPARE,50000,0.000043232,1156550703.183,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,p1024,1024,COMPARE,50000,0.000325856,153442023.607,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,p1024,1024,COMPARE,50000,0.007916520,6315906.484,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,p1024,1024,COMPARE,50000,0.000170356,293502996.623,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,p1024,1024,COMPARE,50000,0.006053576,8259580.765,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,p1024,1024,COMPARE,50000,0.000102995,485460631.160,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,p1024,1024,COMPARE,50000,0.006059080,8252077.855,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,p1024,1024,COMPARE,50000,0.000093246,536216058.139,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,p1024,1024,COMPARE,50000,0.007601532,6577621.426,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p1024,1024,REDUCE,1562,0.000143417,10891318.250,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p1024,1024,REDUCE,1562,0.010325341,151278.296,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p1024,1024,REDUCE,1562,0.010419275,149914.462,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,p1024,1024,REDUCE,50000,0.000047200,1059322033.898,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,p1024,1024,REDUCE,50000,0.004796114,10425106.673,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,p1024,1024,REDUCE,50000,0.012405995,4030309.544,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,p1024,1024,REDUCE,50000,0.002657057,18817812.205,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,p1024,1024,REDUCE,50000,0.008580346,5827270.828,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,p1024,1024,REDUCE,50000,0.000655586,76267642.616,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,p1024,1024,REDUCE,50000,0.006701420,7461105.320,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,p1024,1024,REDUCE,50000,0.000605451,82583058.058,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,p1024,1024,REDUCE,50000,0.008237098,6070098.964,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p1024,1024,MODMUL,781,0.001092614,714799.538,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p1024,1024,MODMUL,781,0.010616818,73562.530,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p1024,1024,MODMUL,781,0.010548404,74039.637,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,p1024,1024,MODMUL,50000,0.000348512,143467082.912,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,p1024,1024,MODMUL,50000,0.020160405,2480108.904,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,p1024,1024,MODMUL,50000,0.027834360,1796340.924,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,p1024,1024,MODMUL,50000,0.006959731,7184185.719,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,p1024,1024,MODMUL,50000,0.012953785,3859875.700,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,p1024,1024,MODMUL,50000,0.003020061,16555957.272,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,p1024,1024,MODMUL,50000,0.009125812,5478964.515,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,p1024,1024,MODMUL,50000,0.002277349,21955352.421,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,p1024,1024,MODMUL,50000,0.009843001,5079751.579,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p1024,1024,MODEXP,195,0.082149477,2373.722,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p1024,1024,MODEXP,195,0.010758245,18125.633,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p1024,1024,MODEXP,195,0.003434154,56782.544,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,p1024,1024,MODEXP,50000,2.038345814,24529.694,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,p1024,1024,MODEXP,50000,5.054753558,9891.679,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,p1024,1024,MODEXP,50000,5.069735392,9862.448,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,p1024,1024,MODEXP,50000,0.663601702,75346.401,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,p1024,1024,MODEXP,50000,0.667924179,74858.796,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,p1024,1024,MODEXP,50000,0.355998784,140449.918,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,p1024,1024,MODEXP,50000,0.362020591,138113.691,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,p1024,1024,MODEXP,50000,0.188284096,265556.152,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,p1024,1024,MODEXP,50000,0.191521526,261067.260,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p1024,1024,EXPONENTIATION,195,0.011907820,16375.793,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p1024,1024,EXPONENTIATION,195,0.001063585,183342.179,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p1024,1024,EXPONENTIATION,195,0.012579925,15500.887,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,p1024,1024,EXPONENTIATION,50000,4.765562387,10491.941,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,p1024,1024,EXPONENTIATION,50000,4.728163365,10574.931,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,p1024,1024,EXPONENTIATION,50000,1.001696738,49915.307,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,p1024,1024,EXPONENTIATION,50000,1.000041793,49997.910,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,p1024,1024,EXPONENTIATION,50000,0.262069371,190789.178,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,p1024,1024,EXPONENTIATION,50000,0.268012518,186558.450,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,p1024,1024,EXPONENTIATION,50000,0.247173551,202287.016,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,p1024,1024,EXPONENTIATION,50000,0.250778312,199379.283,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p1024,1024,DIVIDE,1562,0.000103768,15052807.290,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p1024,1024,DIVIDE,1562,0.004151660,376235.052,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p1024,1024,DIVIDE,1562,0.000070655,22107417.152,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,p1024,1024,DIVIDE,50000,0.000088992,561848256.023,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,p1024,1024,DIVIDE,50000,0.079282400,630656.993,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,p1024,1024,DIVIDE,50000,0.089298379,559920.578,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,p1024,1024,DIVIDE,50000,0.011668467,4285053.012,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,p1024,1024,DIVIDE,50000,0.018203812,2746677.461,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,p1024,1024,DIVIDE,50000,0.003191012,15669009.968,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,p1024,1024,DIVIDE,50000,0.010912226,4582016.567,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,p1024,1024,DIVIDE,50000,0.002873838,17398336.617,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,p1024,1024,DIVIDE,50000,0.007127721,7014864.826,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p1024,1024,ISQRT,390,0.000148105,2633266.448,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p1024,1024,ISQRT,390,0.000035700,10924375.623,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,p1024,1024,ISQRT,50000,1.097824112,45544.636,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,p1024,1024,ISQRT,50000,1.103120474,45325.965,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,p1024,1024,ISQRT,50000,0.285843492,174920.897,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,p1024,1024,ISQRT,50000,0.292313353,171049.319,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,p1024,1024,ISQRT,50000,0.052979072,943768.891,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,p1024,1024,ISQRT,50000,0.059070910,846440.321,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,p1024,1024,ISQRT,50000,0.050328458,993473.711,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,p1024,1024,ISQRT,50000,0.053504567,934499.666,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p1024,1024,MODMUL_R2,12500,0.007844740,1593424.384,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p1024,1024,MODMUL_R2,12500,0.000379352,32950925.502,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p1024,1024,MODMUL_R2,12500,0.001332220,9382834.800,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,p1024,1024,MODMUL_R2,50000,0.000119712,417669072.441,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,p1024,1024,MODMUL_R2,50000,0.004952588,10095731.881,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,p1024,1024,MODMUL_R2,50000,0.012683259,3942204.453,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,p1024,1024,MODMUL_R2,50000,0.000746411,66987229.844,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,p1024,1024,MODMUL_R2,50000,0.005699837,8772180.681,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,p1024,1024,MODMUL_R2,50000,0.000445409,112256351.445,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,p1024,1024,MODMUL_R2,50000,0.006364630,7855916.270,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,p1024,1024,MODMUL_R2,50000,0.000327245,152790676.105,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,p1024,1024,MODMUL_R2,50000,0.003652361,13689774.958,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p2048,2048,ADD,6250,0.000260763,23968125.532,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p2048,2048,ADD,6250,0.000045364,137774373.451,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p2048,2048,ADD,6250,0.000103221,60549703.805,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,p2048,2048,ADD,50000,0.000082656,604916763.453,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,p2048,2048,ADD,50000,0.000750234,66645873.037,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,p2048,2048,ADD,50000,0.016188448,3088622.207,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,p2048,2048,ADD,50000,0.000370743,134864318.115,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,p2048,2048,ADD,50000,0.010052940,4973669.398,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32,p2048,2048,ADD,50000,0.000189294,264139393.846,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32,p2048,2048,ADD,50000,0.012763481,3917426.600,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,p2048,2048,ADD,50000,0.000189080,264438468.548,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,p2048,2048,ADD,50000,0.011591792,4313396.888,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,p2048,2048,ADD,50000,0.000187074,267273817.172,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,p2048,2048,ADD,50000,0.006324200,7906138.376,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p2048,2048,SUBTRACT,6250,0.000235590,26529142.105,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p2048,2048,SUBTRACT,6250,0.000048192,129689615.089,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p2048,2048,SUBTRACT,6250,0.000104850,59608974.585,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,p2048,2048,SUBTRACT,50000,0.000082912,603049015.824,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,p2048,2048,SUBTRACT,50000,0.000862687,57958451.584,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,p2048,2048,SUBTRACT,50000,0.016286383,3070049.383,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,p2048,2048,SUBTRACT,50000,0.000370104,135097168.348,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,p2048,2048,SUBTRACT,50000,0.011393322,4388535.660,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32,p2048,2048,SUBTRACT,50000,0.000200740,249078477.774,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32,p2048,2048,SUBTRACT,50000,0.012746750,3922568.489,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,p2048,2048,SUBTRACT,50000,0.000191911,260537426.744,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,p2048,2048,SUBTRACT,50000,0.011531628,4335901.239,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,p2048,2048,SUBTRACT,50000,0.000185412,269669700.511,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,p2048,2048,SUBTRACT,50000,0.006581539,7597007.262,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p2048,2048,ADDMOD,6250,0.000716319,8725162.748,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p2048,2048,ADDMOD,6250,0.000083305,75025508.797,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p2048,2048,ADDMOD,6250,0.000286107,21844973.339,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,p2048,2048,ADDMOD,50000,0.000085824,582587621.178,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,p2048,2048,ADDMOD,50000,0.001188320,42076206.674,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,p2048,2048,ADDMOD,50000,0.016685153,2996676.145,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,p2048,2048,ADDMOD,50000,0.000525099,95220140.798,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,p2048,2048,ADDMOD,50000,0.011546893,4330169.146,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32,p2048,2048,ADDMOD,50000,0.000274335,182258894.302,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32,p2048,2048,ADDMOD,50000,0.012807754,3903885.104,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,p2048,2048,ADDMOD,50000,0.000275213,181677438.917,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,p2048,2048,ADDMOD,50000,0.011925813,4192586.268,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,p2048,2048,ADDMOD,50000,0.000273644,182719158.133,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,p2048,2048,ADDMOD,50000,0.006684808,7479646.449,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p2048,2048,SUBTRACTMOD,6250,0.000551959,11323304.112,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p2048,2048,SUBTRACTMOD,6250,0.005653420,1105525.505,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p2048,2048,SUBTRACTMOD,6250,0.000278962,22404483.643,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,p2048,2048,SUBTRACTMOD,50000,0.000086208,579992576.095,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,p2048,2048,SUBTRACTMOD,50000,0.001185376,42180708.106,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,p2048,2048,SUBTRACTMOD,50000,0.016682562,2997141.564,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,p2048,2048,SUBTRACTMOD,50000,0.000556928,89778213.640,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,p2048,2048,SUBTRACTMOD,50000,0.011533045,4335368.465,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32,p2048,2048,SUBTRACTMOD,50000,0.000287140,174131143.683,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32,p2048,2048,SUBTRACTMOD,50000,0.012812104,3902559.654,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,p2048,2048,SUBTRACTMOD,50000,0.000272919,183204568.589,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,p2048,2048,SUBTRACTMOD,50000,0.011529197,4336815.484,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,p2048,2048,SUBTRACTMOD,50000,0.000271074,184451472.617,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,p2048,2048,SUBTRACTMOD,50000,0.006756143,7400672.179,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.004285542,1458391.960,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.000320257,19515577.848,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.000451908,13830247.429,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.229034917,218307.325,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.248019879,201596.744,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.068581678,729057.694,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.081173170,615967.074,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.017964606,2783250.567,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.034423010,1452516.792,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.001656491,30184287.873,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.016192121,3087921.595,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.001425799,35068057.011,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.009612402,5201613.464,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.004285106,1458540.359,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.000320084,19526125.273,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.000452162,13822479.572,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.000155200,322164948.454,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.071736645,696993.845,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.090501924,552474.443,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.018190187,2748734.808,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.031225598,1601250.355,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.004581914,10912470.360,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.020913278,2390825.581,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.004755666,10513774.486,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.019269538,2594769.009,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.004723938,10584389.403,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.012844496,3892717.922,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.027577270,226635.922,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.001884309,3316865.771,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.000587926,10630589.322,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000249952,200038407.374,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.121165896,412657.370,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.136768938,365580.085,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.002927447,17079728.006,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.012344255,4050467.203,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000540105,92574584.451,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.012977474,3852829.905,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000895753,55818956.704,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.012145561,4116730.367,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000732411,68267675.796,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.007143326,6999540.563,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p2048,2048,COMPARE,6250,0.000101983,61284709.898,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p2048,2048,COMPARE,6250,0.000046551,134261362.943,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p2048,2048,COMPARE,6250,0.000098250,63613237.105,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,p2048,2048,COMPARE,50000,0.000082784,603981445.690,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,p2048,2048,COMPARE,50000,0.000641159,77983774.458,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,p2048,2048,COMPARE,50000,0.015767563,3171067.079,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,p2048,2048,COMPARE,50000,0.000321867,155343693.613,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,p2048,2048,COMPARE,50000,0.009896256,5052415.767,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,p2048,2048,COMPARE,50000,0.000168832,296152354.895,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,p2048,2048,COMPARE,50000,0.011551393,4328482.258,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,p2048,2048,COMPARE,50000,0.000168257,297164379.579,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,p2048,2048,COMPARE,50000,0.006647774,7521314.633,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p2048,2048,REDUCE,781,0.000050145,15574824.672,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p2048,2048,REDUCE,781,0.004568340,170959.256,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p2048,2048,REDUCE,781,0.000116254,6718050.134,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,p2048,2048,REDUCE,50000,0.000084224,593655015.198,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,p2048,2048,REDUCE,50000,0.822766769,60770.563,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,p2048,2048,REDUCE,50000,0.838918724,59600.529,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,p2048,2048,REDUCE,50000,0.011406606,4383424.809,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,p2048,2048,REDUCE,50000,0.017944906,2786306.038,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,p2048,2048,REDUCE,50000,0.002123313,23548105.917,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,p2048,2048,REDUCE,50000,0.013351867,3744794.649,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,p2048,2048,REDUCE,50000,0.001911448,26158180.030,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,p2048,2048,REDUCE,50000,0.007467611,6695581.765,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p2048,2048,MODMUL,390,0.000928587,419992.936,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p2048,2048,MODMUL,390,0.000092251,4227596.748,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p2048,2048,MODMUL,390,0.000239104,1631089.576,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,p2048,2048,MODMUL,50000,0.001180576,42352207.736,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,p2048,2048,MODMUL,50000,1.177471322,42463.879,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,p2048,2048,MODMUL,50000,1.199943485,41668.629,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,p2048,2048,MODMUL,50000,0.033788376,1479798.851,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,p2048,2048,MODMUL,50000,0.040365385,1238685.076,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,p2048,2048,MODMUL,50000,0.012583709,3973391.283,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,p2048,2048,MODMUL,50000,0.023828415,2098335.121,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,p2048,2048,MODMUL,50000,0.009216283,5425180.639,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,p2048,2048,MODMUL,50000,0.015690247,3186692.969,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p2048,2048,MODEXP,97,0.297435162,326.121,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p2048,2048,MODEXP,97,0.018483512,5247.920,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p2048,2048,MODEXP,97,0.011098964,8739.554,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,p2048,2048,MODEXP,50000,2.016681194,24793.210,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,p2048,2048,MODEXP,50000,208.522423268,239.782,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,p2048,2048,MODEXP,50000,208.544429070,239.757,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,p2048,2048,MODEXP,50000,26.003051717,1922.851,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,p2048,2048,MODEXP,50000,26.021491674,1921.489,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,p2048,2048,MODEXP,50000,3.833457428,13043.056,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,p2048,2048,MODEXP,50000,3.857114792,12963.057,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-o64,p2048,2048,MODEXP,50000,0.000000000,inf,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-o64,p2048,2048,MODEXP,50000,0.000000000,inf,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p2048,2048,EXPONENTIATION,97,0.034121485,2842.784,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p2048,2048,EXPONENTIATION,97,0.003553854,27294.312,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p2048,2048,EXPONENTIATION,97,0.039862689,2433.353,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,p2048,2048,EXPONENTIATION,50000,34.321634396,1456.807,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,p2048,2048,EXPONENTIATION,50000,34.349362980,1455.631,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,p2048,2048,EXPONENTIATION,50000,9.498624434,5263.920,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,p2048,2048,EXPONENTIATION,50000,9.501773608,5262.175,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,p2048,2048,EXPONENTIATION,50000,2.797280249,17874.505,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,p2048,2048,EXPONENTIATION,50000,2.809839386,17794.611,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p2048,2048,DIVIDE,781,0.000081981,9526603.413,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p2048,2048,DIVIDE,781,0.004160322,187725.853,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p2048,2048,DIVIDE,781,0.004928533,158465.004,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,p2048,2048,DIVIDE,50000,0.000099424,502896684.905,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,p2048,2048,DIVIDE,50000,1.933082038,25865.431,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,p2048,2048,DIVIDE,50000,1.949398668,25648.935,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,p2048,2048,DIVIDE,50000,0.490848695,101864.384,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,p2048,2048,DIVIDE,50000,0.508213361,98383.875,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,p2048,2048,DIVIDE,50000,0.035935688,1391374.504,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,p2048,2048,DIVIDE,50000,0.049306177,1014071.727,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p2048,2048,ISQRT,195,0.000153746,1268325.280,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p2048,2048,ISQRT,195,0.004779904,40795.798,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,p2048,2048,ISQRT,50000,26.806626925,1865.210,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,p2048,2048,ISQRT,50000,26.886296038,1859.683,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,p2048,2048,ISQRT,50000,11.109256371,4500.751,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,p2048,2048,ISQRT,50000,11.119809088,4496.480,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,p2048,2048,ISQRT,50000,0.122139304,409368.633,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,p2048,2048,ISQRT,50000,0.131690178,379679.037,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p2048,2048,MODMUL_R2,6250,0.015346974,407246.405,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p2048,2048,MODMUL_R2,6250,0.004812094,1298810.878,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p2048,2048,MODMUL_R2,6250,0.001650362,3787047.968,0
library,Tesla V100-PCIE-32GB,gpu,cgbn,p2048,2048,MODMUL_R2,50000,0.000434528,115067383.460,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w8,p2048,2048,MODMUL_R2,50000,0.086320458,579236.964,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w8,p2048,2048,MODMUL_R2,50000,0.098060780,509887.847,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w16,p2048,2048,MODMUL_R2,50000,0.003379089,14796887.304,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w16,p2048,2048,MODMUL_R2,50000,0.014874624,3361429.512,0
opencl-kernel,Tesla V100-PCIE-32GB,GPU,w32-opt,p2048,2048,MODMUL_R2,50000,0.001618374,30895208.955,0
opencl-e2e,Tesla V100-PCIE-32GB,GPU,w32-opt,p2048,2048,MODMUL_R2,50000,0.011319055,4417329.896,0
```
