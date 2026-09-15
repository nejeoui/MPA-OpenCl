# MPA-OpenCL benchmark report - Quadro RTX 6000

> **Note.** The multi-threaded GMP and OpenSSL baseline columns have been
> removed from this report: they predate the 2026-09-12 timing fix and were
> understated (see `reports/README.md`). The single-threaded GMP column, the
> OpenCL-on-CPU rows and all MPA measurements are unaffected and were verified
> against GMP before timing.


> **Partial report.** The run was interrupted or hit its time budget.
> Rows that never ran are marked `n/a`.

## 1. System under test

2 OpenCL device(s) exercised with the identical kernels and operands.

### Device 0 - Quadro RTX 6000 (GPU)

| Property | Value |
|---|---|
| Model | Quadro RTX 6000 |
| Type | GPU |
| Vendor | NVIDIA Corporation |
| Device memory | 21.97 GiB |
| Max single allocation | 5.49 GiB |
| Local memory | 48 KiB |
| Global cache | 2304 KiB |
| Compute units | 72 |
| Max clock | 1560 MHz |
| Max work-group size | 1024 |
| OpenCL version | OpenCL 3.0 CUDA |
| Driver | 570.86.10 |

### Device 1 - cpu-haswell-Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz (CPU)

| Property | Value |
|---|---|
| Model | cpu-haswell-Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz |
| Type | CPU |
| Vendor | GenuineIntel |
| Device memory | 375.75 GiB |
| Max single allocation | 128.00 GiB |
| Local memory | 256 KiB |
| Global cache | 56320 KiB |
| Compute units | 88 |
| Max clock | 3600 MHz |
| Max work-group size | 4096 |
| OpenCL version | OpenCL 3.0 PoCL HSTR: cpu-x86_64-pc-linux-gnu-haswell |
| Driver | 5.0+debian |

### Host

| Property | Value |
|---|---|
| CPU | Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz |
| Logical cores | 88 |
| OpenMP threads used | 88 |
| RAM | 377.8 GB |
| OS | Ubuntu 24.04.4 LTS |
| Kernel | 6.8.0-101-generic |
| Arch | x86_64 |
| GMP | 6.3.0 |
| OpenSSL | OpenSSL 3.0.13 30 Jan 2024 |
| CGBN | cgbn_results.tsv loaded |

## 2. Method

- Workload auto-sized from the device and host: --min-items from 700 x compute units, --items from ten times that capped by host RAM. Either flag, given explicitly, overrides its half.
- Base workload 50000 items, scaled down per operator by its cost weight and by modulus size. Device rows honour --min-items (50400) so the GPU is not left idle; the CPU libraries keep the smaller count because a full-width MODEXP there costs minutes. Both counts appear in every row as dev/cpu, and throughput is per-second so they remain comparable.
- 5 timed repetitions, **minimum** reported. Two untimed warm-up launches precede them.
- `kernel` times `clEnqueueNDRangeKernel` + `clFinish` only. `e2e` adds the host->device operand writes and the device->host result read.
- Every OpenCL device runs the same kernels on the same operands, so GPU and CPU-OpenCL columns are directly comparable.
- CPU library baselines (GMP, OpenSSL) run those same operands, with temporaries preallocated outside the timed region, so the figure is the arithmetic and not marshalling. The generator is reseeded per modulus and operation so every backend sees identical inputs.
- OpenSSL rows time the nearest BN primitive, which is not always semantically identical (its Montgomery routine expects Montgomery-domain inputs); they measure comparable work, not identical results. Correctness is judged against GMP only.
- Every device cell is checked word-for-word against GMP before it is timed. A cell that mismatches is reported and excluded from the speedup tables.
- Total wall time 5411.6 s.

## 3. Correctness

| Device | Kernel | Configs run | Passed | Mismatched | Build/launch failed |
|---|---|---|---|---|---|
| [0] GPU | `mpaKernels_8bits.cl` (w8) | 75 | 75 | 0 | 0 |
| [0] GPU | `mpaKernel_16bits.cl` (w16) | 72 | 72 | 0 | 0 |

**All configurations correct** - 147 configurations, 0 problems.

## 4. Throughput per device

Operations per second, higher is better. Kernel-only timings.

### Device 0 - Quadro RTX 6000 (GPU)

#### secp256k1 (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | 802.14 M | 1.52 G | - | - | - | - | - | 47.78 M | 1.44 G |
| SUBTRACT | 50000 / 50000 | 804.30 M | 1.53 G | - | - | - | - | - | 65.74 M | 1.48 G |
| ADDMOD | 50000 / 50000 | 549.54 M | 1.08 G | - | - | - | - | - | 16.18 M | 1.50 G |
| SUBTRACTMOD | 50000 / 50000 | 537.41 M | 1.11 G | - | - | - | - | - | 22.08 M | 1.51 G |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | 20.74 M | 91.38 M | - | - | - | - | - | 43.79 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | 68.73 M | 270.45 M | - | - | - | - | - | 38.60 M | 1.49 G |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | 271.63 M | 1.05 G | - | - | - | - | - | 5.58 M | 1.49 G |
| COMPARE | 50000 / 50000 | 839.28 M | 1.53 G | - | - | - | - | - | 78.01 M | 1.44 G |
| REDUCE | 50000 / 6250 | 149.45 M | 278.45 M | - | - | - | - | - | 27.64 M | 1.44 G |
| MODMUL | 50000 / 3125 | 57.25 M | 117.22 M | - | - | - | - | - | 5.99 M | 1.06 G |
| MODEXP | 50000 / 781 | 1.58 M | 7.84 M | - | - | - | - | - | 57.69 k | 656.30 k |
| EXPONENTIATION | 50000 / 781 | 1.05 M | 2.54 M | - | - | - | - | - | 133.31 k | n/a |
| DIVIDE | 50000 / 6250 | 96.30 M | 149.18 M | - | - | - | - | - | 13.16 M | 1.52 G |
| ISQRT | 50000 / 1562 | 8.09 M | 10.10 M | - | - | - | - | - | 6.52 M | n/a |
| MODMUL_R2 | 50000 / 50000 | 272.56 M | 811.03 M | - | - | - | - | - | 5.86 M | 1.36 G |

#### rsa256(composite) (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | 967.54 M | 1.65 G | - | - | - | - | - | 33.21 M | 1.52 G |
| SUBTRACT | 50000 / 50000 | 975.28 M | 1.65 G | - | - | - | - | - | 40.13 M | 1.52 G |
| ADDMOD | 50000 / 50000 | 710.90 M | 1.17 G | - | - | - | - | - | 12.67 M | 1.49 G |
| SUBTRACTMOD | 50000 / 50000 | 660.36 M | 1.18 G | - | - | - | - | - | 13.94 M | 1.52 G |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | 25.42 M | 96.45 M | - | - | - | - | - | 26.70 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | 84.21 M | 285.74 M | - | - | - | - | - | 42.67 M | 1.53 G |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | 329.86 M | 1.09 G | - | - | - | - | - | 3.44 M | 1.51 G |
| COMPARE | 50000 / 50000 | 999.72 M | 1.60 G | - | - | - | - | - | 114.20 M | 1.50 G |
| REDUCE | 50000 / 6250 | 182.53 M | 294.77 M | - | - | - | - | - | 16.87 M | 1.51 G |
| MODMUL | 50000 / 3125 | 69.98 M | 123.89 M | - | - | - | - | - | 5.97 M | 1.28 G |
| MODEXP | 50000 / 781 | 1.58 M | 7.73 M | - | - | - | - | - | 60.35 k | 678.55 k |
| EXPONENTIATION | 50000 / 781 | 1.05 M | 2.61 M | - | - | - | - | - | 163.95 k | n/a |
| DIVIDE | 50000 / 6250 | 97.41 M | 147.53 M | - | - | - | - | - | 12.99 M | 1.48 G |
| ISQRT | 50000 / 1562 | 8.13 M | 10.09 M | - | - | - | - | - | 6.55 M | n/a |
| MODMUL_R2 | 50000 / 50000 | 274.36 M | 861.86 M | - | - | - | - | - | 5.93 M | 1.36 G |

#### brainpoolP512r1 (512-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 25000 | 373.27 M | 643.16 M | - | - | - | - | - | 49.44 M | 1.22 G |
| SUBTRACT | 50000 / 25000 | 372.02 M | 640.11 M | - | - | - | - | - | 58.57 M | 1.21 G |
| ADDMOD | 50000 / 25000 | 282.75 M | 525.18 M | - | - | - | - | - | 17.56 M | 1.25 G |
| SUBTRACTMOD | 50000 / 25000 | 246.59 M | 487.67 M | - | - | - | - | - | 20.62 M | 1.22 G |
| MULTIPLYOPERANDSCANNING | 50000 / 25000 | 3.62 M | 13.52 M | - | - | - | - | - | 19.26 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 25000 | 12.59 M | 46.34 M | - | - | - | - | - | 17.76 M | 1.22 G |
| MONTGOMERYMULTIPLICATION | 50000 / 25000 | 88.80 M | 293.73 M | - | - | - | - | - | 2.51 M | 1.22 G |
| COMPARE | 50000 / 25000 | 439.45 M | 730.60 M | - | - | - | - | - | 111.31 M | 1.22 G |
| REDUCE | 50000 / 3125 | 63.14 M | 80.46 M | - | - | - | - | - | 26.46 M | 1.22 G |
| MODMUL | 50000 / 1562 | 22.71 M | 26.20 M | - | - | - | - | - | 5.02 M | 440.26 M |
| MODEXP | 50000 / 390 | 108.36 k | 1.13 M | - | - | - | - | - | 18.17 k | 74.16 k |
| EXPONENTIATION | 50000 / 390 | 123.50 k | 310.03 k | - | - | - | - | - | 74.04 k | n/a |
| DIVIDE | 50000 / 3125 | 32.33 M | 34.95 M | - | - | - | - | - | 6.46 M | 852.89 M |
| ISQRT | 50000 / 781 | 1.36 M | 1.19 M | - | - | - | - | - | 2.07 M | n/a |
| MODMUL_R2 | 50000 / 25000 | 61.30 M | 271.64 M | - | - | - | - | - | 3.03 M | 1.02 G |

#### p1024 (1024-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 12500 | 57.51 M | 102.26 M | - | - | - | - | - | 24.24 M | 787.55 M |
| SUBTRACT | 50000 / 12500 | 55.59 M | 97.57 M | - | - | - | - | - | 29.83 M | 787.55 M |
| ADDMOD | 50000 / 12500 | 58.15 M | 109.24 M | - | - | - | - | - | 7.82 M | 793.95 M |
| SUBTRACTMOD | 50000 / 12500 | 50.32 M | 105.02 M | - | - | - | - | - | 10.89 M | 802.10 M |
| MULTIPLYOPERANDSCANNING | 50000 / 12500 | 186.33 k | 713.47 k | - | - | - | - | - | 5.74 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 12500 | 1.62 M | 5.44 M | - | - | - | - | - | 5.75 M | 787.55 M |
| MONTGOMERYMULTIPLICATION | 50000 / 12500 | 16.39 M | 82.20 M | - | - | - | - | - | 850.33 k | 734.95 M |
| COMPARE | 50000 / 12500 | 176.49 M | 324.45 M | - | - | - | - | - | 87.43 M | 787.55 M |
| REDUCE | 50000 / 1562 | 11.97 M | 22.00 M | - | - | - | - | - | 34.40 M | 813.80 M |
| MODMUL | 50000 / 781 | 2.87 M | 7.64 M | - | - | - | - | - | 1.93 M | 152.80 M |
| MODEXP | 50000 / 195 | 12.10 k | 73.92 k | - | - | - | - | - | 2.84 k | 32.38 k |
| EXPONENTIATION | 50000 / 195 | 10.46 k | 44.95 k | - | - | - | - | - | 19.67 k | n/a |
| DIVIDE | 50000 / 1562 | 277.25 k | 3.03 M | - | - | - | - | - | 6.00 M | 606.09 M |
| ISQRT | 50000 / 390 | 25.05 k | 161.82 k | - | - | - | - | - | 1.05 M | n/a |
| MODMUL_R2 | 50000 / 12500 | 11.34 M | 69.71 M | - | - | - | - | - | 1.16 M | 407.11 M |

#### p2048 (2048-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 6250 | 19.21 M | 39.65 M | - | - | - | - | - | 17.80 M | 412.27 M |
| SUBTRACT | 50000 / 6250 | 19.54 M | 39.14 M | - | - | - | - | - | 19.81 M | 411.51 M |
| ADDMOD | 50000 / 6250 | 20.55 M | 39.72 M | - | - | - | - | - | 6.64 M | 415.67 M |
| SUBTRACTMOD | 50000 / 6250 | 19.71 M | 40.50 M | - | - | - | - | - | 13.33 M | 420.93 M |
| MULTIPLYOPERANDSCANNING | 50000 / 6250 | 38.28 k | 151.50 k | - | - | - | - | - | 1.06 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 6250 | 407.75 k | 1.50 M | - | - | - | - | - | 1.04 M | 310.14 M |
| MONTGOMERYMULTIPLICATION | 50000 / 6250 | 361.53 k | 19.83 M | - | - | - | - | - | 269.52 k | 178.88 M |
| COMPARE | 50000 / 6250 | 85.79 M | 143.90 M | - | - | - | - | - | 48.68 M | 411.40 M |
| REDUCE | 50000 / 781 | 39.04 k | 4.41 M | - | - | - | - | - | 14.43 M | 420.93 M |
| MODMUL | 50000 / 390 | 22.51 k | 1.85 M | - | - | - | - | - | 400.09 k | 36.47 M |
| MODEXP | 50000 / 97 | 165.8 | over budget | - | - | - | - | - | 406.8 | 25.35 k |
| EXPONENTIATION | 50000 / 97 | over budget | over budget | - | - | - | - | - | 2.08 k | n/a |
| DIVIDE | 50000 / 781 | 16.45 k | - | - | - | - | - | - | 8.48 M | 429.26 M |
| ISQRT | 50000 / 195 | 1.29 k | - | - | - | - | - | - | 988.81 k | n/a |
| MODMUL_R2 | 50000 / 6250 | 402.60 k | - | - | - | - | - | - | 387.53 k | 113.63 M |

### Device 1 - cpu-haswell-Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz (CPU)

#### secp256k1 (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | - | - | - | - | - | - | - | 47.78 M | 1.44 G |
| SUBTRACT | 50000 / 50000 | - | - | - | - | - | - | - | 65.74 M | 1.48 G |
| ADDMOD | 50000 / 50000 | - | - | - | - | - | - | - | 16.18 M | 1.50 G |
| SUBTRACTMOD | 50000 / 50000 | - | - | - | - | - | - | - | 22.08 M | 1.51 G |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | - | - | - | - | - | - | - | 43.79 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | - | - | - | - | - | - | - | 38.60 M | 1.49 G |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | - | - | - | - | - | - | - | 5.58 M | 1.49 G |
| COMPARE | 50000 / 50000 | - | - | - | - | - | - | - | 78.01 M | 1.44 G |
| REDUCE | 50000 / 6250 | - | - | - | - | - | - | - | 27.64 M | 1.44 G |
| MODMUL | 50000 / 3125 | - | - | - | - | - | - | - | 5.99 M | 1.06 G |
| MODEXP | 50000 / 781 | - | - | - | - | - | - | - | 57.69 k | 656.30 k |
| EXPONENTIATION | 50000 / 781 | - | - | - | - | - | - | - | 133.31 k | n/a |
| DIVIDE | 50000 / 6250 | - | - | - | - | - | - | - | 13.16 M | 1.52 G |
| ISQRT | 50000 / 1562 | - | - | - | - | - | - | - | 6.52 M | n/a |
| MODMUL_R2 | 50000 / 50000 | - | - | - | - | - | - | - | 5.86 M | 1.36 G |

#### rsa256(composite) (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | - | - | - | - | - | - | - | 33.21 M | 1.52 G |
| SUBTRACT | 50000 / 50000 | - | - | - | - | - | - | - | 40.13 M | 1.52 G |
| ADDMOD | 50000 / 50000 | - | - | - | - | - | - | - | 12.67 M | 1.49 G |
| SUBTRACTMOD | 50000 / 50000 | - | - | - | - | - | - | - | 13.94 M | 1.52 G |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | - | - | - | - | - | - | - | 26.70 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | - | - | - | - | - | - | - | 42.67 M | 1.53 G |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | - | - | - | - | - | - | - | 3.44 M | 1.51 G |
| COMPARE | 50000 / 50000 | - | - | - | - | - | - | - | 114.20 M | 1.50 G |
| REDUCE | 50000 / 6250 | - | - | - | - | - | - | - | 16.87 M | 1.51 G |
| MODMUL | 50000 / 3125 | - | - | - | - | - | - | - | 5.97 M | 1.28 G |
| MODEXP | 50000 / 781 | - | - | - | - | - | - | - | 60.35 k | 678.55 k |
| EXPONENTIATION | 50000 / 781 | - | - | - | - | - | - | - | 163.95 k | n/a |
| DIVIDE | 50000 / 6250 | - | - | - | - | - | - | - | 12.99 M | 1.48 G |
| ISQRT | 50000 / 1562 | - | - | - | - | - | - | - | 6.55 M | n/a |
| MODMUL_R2 | 50000 / 50000 | - | - | - | - | - | - | - | 5.93 M | 1.36 G |

#### brainpoolP512r1 (512-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 25000 | - | - | - | - | - | - | - | 49.44 M | 1.22 G |
| SUBTRACT | 50000 / 25000 | - | - | - | - | - | - | - | 58.57 M | 1.21 G |
| ADDMOD | 50000 / 25000 | - | - | - | - | - | - | - | 17.56 M | 1.25 G |
| SUBTRACTMOD | 50000 / 25000 | - | - | - | - | - | - | - | 20.62 M | 1.22 G |
| MULTIPLYOPERANDSCANNING | 50000 / 25000 | - | - | - | - | - | - | - | 19.26 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 25000 | - | - | - | - | - | - | - | 17.76 M | 1.22 G |
| MONTGOMERYMULTIPLICATION | 50000 / 25000 | - | - | - | - | - | - | - | 2.51 M | 1.22 G |
| COMPARE | 50000 / 25000 | - | - | - | - | - | - | - | 111.31 M | 1.22 G |
| REDUCE | 50000 / 3125 | - | - | - | - | - | - | - | 26.46 M | 1.22 G |
| MODMUL | 50000 / 1562 | - | - | - | - | - | - | - | 5.02 M | 440.26 M |
| MODEXP | 50000 / 390 | - | - | - | - | - | - | - | 18.17 k | 74.16 k |
| EXPONENTIATION | 50000 / 390 | - | - | - | - | - | - | - | 74.04 k | n/a |
| DIVIDE | 50000 / 3125 | - | - | - | - | - | - | - | 6.46 M | 852.89 M |
| ISQRT | 50000 / 781 | - | - | - | - | - | - | - | 2.07 M | n/a |
| MODMUL_R2 | 50000 / 25000 | - | - | - | - | - | - | - | 3.03 M | 1.02 G |

#### p1024 (1024-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 12500 | - | - | - | - | - | - | - | 24.24 M | 787.55 M |
| SUBTRACT | 50000 / 12500 | - | - | - | - | - | - | - | 29.83 M | 787.55 M |
| ADDMOD | 50000 / 12500 | - | - | - | - | - | - | - | 7.82 M | 793.95 M |
| SUBTRACTMOD | 50000 / 12500 | - | - | - | - | - | - | - | 10.89 M | 802.10 M |
| MULTIPLYOPERANDSCANNING | 50000 / 12500 | - | - | - | - | - | - | - | 5.74 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 12500 | - | - | - | - | - | - | - | 5.75 M | 787.55 M |
| MONTGOMERYMULTIPLICATION | 50000 / 12500 | - | - | - | - | - | - | - | 850.33 k | 734.95 M |
| COMPARE | 50000 / 12500 | - | - | - | - | - | - | - | 87.43 M | 787.55 M |
| REDUCE | 50000 / 1562 | - | - | - | - | - | - | - | 34.40 M | 813.80 M |
| MODMUL | 50000 / 781 | - | - | - | - | - | - | - | 1.93 M | 152.80 M |
| MODEXP | 50000 / 195 | - | - | - | - | - | - | - | 2.84 k | 32.38 k |
| EXPONENTIATION | 50000 / 195 | - | - | - | - | - | - | - | 19.67 k | n/a |
| DIVIDE | 50000 / 1562 | - | - | - | - | - | - | - | 6.00 M | 606.09 M |
| ISQRT | 50000 / 390 | - | - | - | - | - | - | - | 1.05 M | n/a |
| MODMUL_R2 | 50000 / 12500 | - | - | - | - | - | - | - | 1.16 M | 407.11 M |

#### p2048 (2048-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 6250 | - | - | - | - | - | - | - | 17.80 M | 412.27 M |
| SUBTRACT | 50000 / 6250 | - | - | - | - | - | - | - | 19.81 M | 411.51 M |
| ADDMOD | 50000 / 6250 | - | - | - | - | - | - | - | 6.64 M | 415.67 M |
| SUBTRACTMOD | 50000 / 6250 | - | - | - | - | - | - | - | 13.33 M | 420.93 M |
| MULTIPLYOPERANDSCANNING | 50000 / 6250 | - | - | - | - | - | - | - | 1.06 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 6250 | - | - | - | - | - | - | - | 1.04 M | 310.14 M |
| MONTGOMERYMULTIPLICATION | 50000 / 6250 | - | - | - | - | - | - | - | 269.52 k | 178.88 M |
| COMPARE | 50000 / 6250 | - | - | - | - | - | - | - | 48.68 M | 411.40 M |
| REDUCE | 50000 / 781 | - | - | - | - | - | - | - | 14.43 M | 420.93 M |
| MODMUL | 50000 / 390 | - | - | - | - | - | - | - | 400.09 k | 36.47 M |
| MODEXP | 50000 / 97 | - | - | - | - | - | - | - | 406.8 | 25.35 k |
| EXPONENTIATION | 50000 / 97 | - | - | - | - | - | - | - | 2.08 k | n/a |
| DIVIDE | 50000 / 781 | - | - | - | - | - | - | - | 8.48 M | 429.26 M |
| ISQRT | 50000 / 195 | - | - | - | - | - | - | - | 988.81 k | n/a |
| MODMUL_R2 | 50000 / 6250 | - | - | - | - | - | - | - | 387.53 k | 113.63 M |

## 5. Head to head

Best OpenCL GPU result against best OpenCL CPU result and the CPU libraries.
Ratios above 1.00x mean the GPU is faster than that baseline.

### secp256k1 (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w16 | 1.52 G | none | n/a | 47.78 M | 1.44 G | n/a | 1.06x |
| SUBTRACT | w16 | 1.53 G | none | n/a | 65.74 M | 1.48 G | n/a | 1.04x |
| ADDMOD | w16 | 1.08 G | none | n/a | 16.18 M | 1.50 G | n/a | 0.72x |
| SUBTRACTMOD | w16 | 1.11 G | none | n/a | 22.08 M | 1.51 G | n/a | 0.74x |
| MULTIPLYOPERANDSCANNING | w16 | 91.38 M | none | n/a | 43.79 M | n/a | n/a | n/a |
| MULTIPLYPRODUCTSCANNING | w16 | 270.45 M | none | n/a | 38.60 M | 1.49 G | n/a | 0.18x |
| MONTGOMERYMULTIPLICATION | w16 | 1.05 G | none | n/a | 5.58 M | 1.49 G | n/a | 0.70x |
| COMPARE | w16 | 1.53 G | none | n/a | 78.01 M | 1.44 G | n/a | 1.06x |
| REDUCE | w16 | 34.81 M | none | n/a | 27.64 M | 1.44 G | n/a | 0.02x |
| MODMUL | w16 | 7.33 M | none | n/a | 5.99 M | 1.06 G | n/a | 0.01x |
| MODEXP | w16 | 122.38 k | none | n/a | 57.69 k | 656.30 k | n/a | 0.19x |
| EXPONENTIATION | w16 | 39.71 k | none | n/a | 133.31 k | n/a | n/a | n/a |
| DIVIDE | w16 | 18.65 M | none | n/a | 13.16 M | 1.52 G | n/a | 0.01x |
| ISQRT | w16 | 315.43 k | none | n/a | 6.52 M | n/a | n/a | n/a |
| MODMUL_R2 | w16 | 811.03 M | none | n/a | 5.86 M | 1.36 G | n/a | 0.60x |

### rsa256(composite) (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w16 | 1.65 G | none | n/a | 33.21 M | 1.52 G | n/a | 1.08x |
| SUBTRACT | w16 | 1.65 G | none | n/a | 40.13 M | 1.52 G | n/a | 1.08x |
| ADDMOD | w16 | 1.17 G | none | n/a | 12.67 M | 1.49 G | n/a | 0.78x |
| SUBTRACTMOD | w16 | 1.18 G | none | n/a | 13.94 M | 1.52 G | n/a | 0.78x |
| MULTIPLYOPERANDSCANNING | w16 | 96.45 M | none | n/a | 26.70 M | n/a | n/a | n/a |
| MULTIPLYPRODUCTSCANNING | w16 | 285.74 M | none | n/a | 42.67 M | 1.53 G | n/a | 0.19x |
| MONTGOMERYMULTIPLICATION | w16 | 1.09 G | none | n/a | 3.44 M | 1.51 G | n/a | 0.73x |
| COMPARE | w16 | 1.60 G | none | n/a | 114.20 M | 1.50 G | n/a | 1.06x |
| REDUCE | w16 | 36.85 M | none | n/a | 16.87 M | 1.51 G | n/a | 0.02x |
| MODMUL | w16 | 7.74 M | none | n/a | 5.97 M | 1.28 G | n/a | 0.01x |
| MODEXP | w16 | 120.79 k | none | n/a | 60.35 k | 678.55 k | n/a | 0.18x |
| EXPONENTIATION | w16 | 40.73 k | none | n/a | 163.95 k | n/a | n/a | n/a |
| DIVIDE | w16 | 18.44 M | none | n/a | 12.99 M | 1.48 G | n/a | 0.01x |
| ISQRT | w16 | 315.15 k | none | n/a | 6.55 M | n/a | n/a | n/a |
| MODMUL_R2 | w16 | 861.86 M | none | n/a | 5.93 M | 1.36 G | n/a | 0.64x |

### brainpoolP512r1 (512-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w16 | 321.58 M | none | n/a | 49.44 M | 1.22 G | n/a | 0.26x |
| SUBTRACT | w16 | 320.05 M | none | n/a | 58.57 M | 1.21 G | n/a | 0.26x |
| ADDMOD | w16 | 262.59 M | none | n/a | 17.56 M | 1.25 G | n/a | 0.21x |
| SUBTRACTMOD | w16 | 243.84 M | none | n/a | 20.62 M | 1.22 G | n/a | 0.20x |
| MULTIPLYOPERANDSCANNING | w16 | 6.76 M | none | n/a | 19.26 M | n/a | n/a | n/a |
| MULTIPLYPRODUCTSCANNING | w16 | 23.17 M | none | n/a | 17.76 M | 1.22 G | n/a | 0.02x |
| MONTGOMERYMULTIPLICATION | w16 | 146.86 M | none | n/a | 2.51 M | 1.22 G | n/a | 0.12x |
| COMPARE | w16 | 365.30 M | none | n/a | 111.31 M | 1.22 G | n/a | 0.30x |
| REDUCE | w16 | 5.03 M | none | n/a | 26.46 M | 1.22 G | n/a | 0.00x |
| MODMUL | w16 | 818.64 k | none | n/a | 5.02 M | 440.26 M | n/a | 0.00x |
| MODEXP | w16 | 8.82 k | none | n/a | 18.17 k | 74.16 k | n/a | 0.12x |
| EXPONENTIATION | w16 | 2.42 k | none | n/a | 74.04 k | n/a | n/a | n/a |
| DIVIDE | w16 | 2.18 M | none | n/a | 6.46 M | 852.89 M | n/a | 0.00x |
| ISQRT | w8 | 21.25 k | none | n/a | 2.07 M | n/a | n/a | n/a |
| MODMUL_R2 | w16 | 135.82 M | none | n/a | 3.03 M | 1.02 G | n/a | 0.13x |

### p1024 (1024-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w16 | 25.56 M | none | n/a | 24.24 M | 787.55 M | n/a | 0.03x |
| SUBTRACT | w16 | 24.39 M | none | n/a | 29.83 M | 787.55 M | n/a | 0.03x |
| ADDMOD | w16 | 27.31 M | none | n/a | 7.82 M | 793.95 M | n/a | 0.03x |
| SUBTRACTMOD | w16 | 26.26 M | none | n/a | 10.89 M | 802.10 M | n/a | 0.03x |
| MULTIPLYOPERANDSCANNING | w16 | 178.37 k | none | n/a | 5.74 M | n/a | n/a | n/a |
| MULTIPLYPRODUCTSCANNING | w16 | 1.36 M | none | n/a | 5.75 M | 787.55 M | n/a | 0.00x |
| MONTGOMERYMULTIPLICATION | w16 | 20.55 M | none | n/a | 850.33 k | 734.95 M | n/a | 0.03x |
| COMPARE | w16 | 81.11 M | none | n/a | 87.43 M | 787.55 M | n/a | 0.10x |
| REDUCE | w16 | 687.32 k | none | n/a | 34.40 M | 813.80 M | n/a | 0.00x |
| MODMUL | w16 | 119.37 k | none | n/a | 1.93 M | 152.80 M | n/a | 0.00x |
| MODEXP | w16 | 288.3 | none | n/a | 2.84 k | 32.38 k | n/a | 0.01x |
| EXPONENTIATION | w16 | 175.3 | none | n/a | 19.67 k | n/a | n/a | n/a |
| DIVIDE | w16 | 94.57 k | none | n/a | 6.00 M | 606.09 M | n/a | 0.00x |
| ISQRT | w16 | 1.26 k | none | n/a | 1.05 M | n/a | n/a | n/a |
| MODMUL_R2 | w16 | 17.43 M | none | n/a | 1.16 M | 407.11 M | n/a | 0.04x |

### p2048 (2048-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w16 | 4.96 M | none | n/a | 17.80 M | 412.27 M | n/a | 0.01x |
| SUBTRACT | w16 | 4.89 M | none | n/a | 19.81 M | 411.51 M | n/a | 0.01x |
| ADDMOD | w16 | 4.96 M | none | n/a | 6.64 M | 415.67 M | n/a | 0.01x |
| SUBTRACTMOD | w16 | 5.06 M | none | n/a | 13.33 M | 420.93 M | n/a | 0.01x |
| MULTIPLYOPERANDSCANNING | w16 | 18.94 k | none | n/a | 1.06 M | n/a | n/a | n/a |
| MULTIPLYPRODUCTSCANNING | w16 | 187.50 k | none | n/a | 1.04 M | 310.14 M | n/a | 0.00x |
| MONTGOMERYMULTIPLICATION | w16 | 2.48 M | none | n/a | 269.52 k | 178.88 M | n/a | 0.01x |
| COMPARE | w16 | 17.99 M | none | n/a | 48.68 M | 411.40 M | n/a | 0.04x |
| REDUCE | w16 | 68.90 k | none | n/a | 14.43 M | 420.93 M | n/a | 0.00x |
| MODMUL | w16 | 14.43 k | none | n/a | 400.09 k | 36.47 M | n/a | 0.00x |
| MODEXP | w8 | 0.3 | none | n/a | 406.8 | 25.35 k | n/a | 0.00x |
| EXPONENTIATION | none | n/a | none | n/a | 2.08 k | n/a | n/a | n/a |
| DIVIDE | w8 | 256.9 | none | n/a | 8.48 M | 429.26 M | n/a | 0.00x |
| ISQRT | w8 | 5.0 | none | n/a | 988.81 k | n/a | n/a | n/a |
| MODMUL_R2 | w8 | 50.33 k | none | n/a | 387.53 k | 113.63 M | n/a | 0.00x |

## 6. Raw data

Also written to `Quadro_RTX_6000_Report.csv` for analysis.

```csv
kind,device,device_type,kernel,modulus,bits,operation,items,seconds,ops_per_sec,mismatches
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,secp256k1,256,ADD,50000,0.001046363,47784552.177,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,secp256k1,256,ADD,50000,0.000100622,496909453.731,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,secp256k1,256,ADD,50000,0.000111207,449610505.159,0
library,Quadro RTX 6000,gpu,cgbn,secp256k1,256,ADD,50000,0.000034816,1436121323.529,0
opencl-kernel,Quadro RTX 6000,GPU,w8,secp256k1,256,ADD,50000,0.000062333,802137923.203,0
opencl-e2e,Quadro RTX 6000,GPU,w8,secp256k1,256,ADD,50000,0.001093304,45732943.303,0
opencl-kernel,Quadro RTX 6000,GPU,w16,secp256k1,256,ADD,50000,0.000032822,1523383780.716,0
opencl-e2e,Quadro RTX 6000,GPU,w16,secp256k1,256,ADD,50000,0.001060393,47152342.375,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,secp256k1,256,SUBTRACT,50000,0.000760619,65735968.302,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,secp256k1,256,SUBTRACT,50000,0.000094561,528759737.625,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,secp256k1,256,SUBTRACT,50000,0.000128295,389725973.460,0
library,Quadro RTX 6000,gpu,cgbn,secp256k1,256,SUBTRACT,50000,0.000033760,1481042654.028,0
opencl-kernel,Quadro RTX 6000,GPU,w8,secp256k1,256,SUBTRACT,50000,0.000062166,804300991.760,0
opencl-e2e,Quadro RTX 6000,GPU,w8,secp256k1,256,SUBTRACT,50000,0.001104048,45287905.644,0
opencl-kernel,Quadro RTX 6000,GPU,w16,secp256k1,256,SUBTRACT,50000,0.000032617,1532953320.770,0
opencl-e2e,Quadro RTX 6000,GPU,w16,secp256k1,256,SUBTRACT,50000,0.001057196,47294901.494,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,secp256k1,256,ADDMOD,50000,0.003090981,16176092.194,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,secp256k1,256,ADDMOD,50000,0.000264401,189106972.222,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,secp256k1,256,ADDMOD,50000,0.000703050,71118714.945,0
library,Quadro RTX 6000,gpu,cgbn,secp256k1,256,ADDMOD,50000,0.000033248,1503849855.630,0
opencl-kernel,Quadro RTX 6000,GPU,w8,secp256k1,256,ADDMOD,50000,0.000090985,549543382.398,0
opencl-e2e,Quadro RTX 6000,GPU,w8,secp256k1,256,ADDMOD,50000,0.001137117,43970852.027,0
opencl-kernel,Quadro RTX 6000,GPU,w16,secp256k1,256,ADDMOD,50000,0.000046227,1081615988.396,0
opencl-e2e,Quadro RTX 6000,GPU,w16,secp256k1,256,ADDMOD,50000,0.001080520,46273996.897,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,secp256k1,256,SUBTRACTMOD,50000,0.002263989,22084910.035,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,secp256k1,256,SUBTRACTMOD,50000,0.000228621,218702506.110,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,secp256k1,256,SUBTRACTMOD,50000,0.000703912,71031583.181,0
library,Quadro RTX 6000,gpu,cgbn,secp256k1,256,SUBTRACTMOD,50000,0.000033216,1505298651.252,0
opencl-kernel,Quadro RTX 6000,GPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.000093039,537408320.320,0
opencl-e2e,Quadro RTX 6000,GPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.001139104,43894134.268,0
opencl-kernel,Quadro RTX 6000,GPU,w16,secp256k1,256,SUBTRACTMOD,50000,0.000045167,1106995983.340,0
opencl-e2e,Quadro RTX 6000,GPU,w16,secp256k1,256,SUBTRACTMOD,50000,0.001072101,46637384.204,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001141937,43785235.363,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000141766,352694069.110,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000252780,197800792.867,0
opencl-kernel,Quadro RTX 6000,GPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.002410533,20742301.365,0
opencl-e2e,Quadro RTX 6000,GPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.003715292,13457893.787,0
opencl-kernel,Quadro RTX 6000,GPU,w16,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000547143,91383838.363,0
opencl-e2e,Quadro RTX 6000,GPU,w16,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001842257,27140626.902,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001295263,38602200.495,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000154953,322677552.590,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000253607,197155761.858,0
library,Quadro RTX 6000,gpu,cgbn,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000033632,1486679352.997,0
opencl-kernel,Quadro RTX 6000,GPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000727490,68729505.745,0
opencl-e2e,Quadro RTX 6000,GPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.002061684,24252022.933,0
opencl-kernel,Quadro RTX 6000,GPU,w16,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000184877,270450310.816,0
opencl-e2e,Quadro RTX 6000,GPU,w16,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001456328,34332931.215,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.008963238,5578340.941,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000847330,59008864.669,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000245875,203355572.222,0
library,Quadro RTX 6000,gpu,cgbn,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000033632,1486679352.997,0
opencl-kernel,Quadro RTX 6000,GPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000184076,271627074.121,0
opencl-e2e,Quadro RTX 6000,GPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.001264837,39530796.667,0
opencl-kernel,Quadro RTX 6000,GPU,w16,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000047743,1047266916.355,0
opencl-e2e,Quadro RTX 6000,GPU,w16,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.001055293,47380215.901,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,secp256k1,256,COMPARE,50000,0.000640964,78007484.707,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,secp256k1,256,COMPARE,50000,0.004401250,11360409.462,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,secp256k1,256,COMPARE,50000,0.005461097,9155670.036,0
library,Quadro RTX 6000,gpu,cgbn,secp256k1,256,COMPARE,50000,0.000034752,1438766114.180,0
opencl-kernel,Quadro RTX 6000,GPU,w8,secp256k1,256,COMPARE,50000,0.000059575,839280440.220,0
opencl-e2e,Quadro RTX 6000,GPU,w8,secp256k1,256,COMPARE,50000,0.001115425,44825980.940,0
opencl-kernel,Quadro RTX 6000,GPU,w16,secp256k1,256,COMPARE,50000,0.000032779,1525374792.590,0
opencl-e2e,Quadro RTX 6000,GPU,w16,secp256k1,256,COMPARE,50000,0.001047721,47722622.406,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,secp256k1,256,REDUCE,6250,0.000226086,27644347.045,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,secp256k1,256,REDUCE,6250,0.007887192,792423.963,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,secp256k1,256,REDUCE,6250,0.007863147,794847.134,0
library,Quadro RTX 6000,gpu,cgbn,secp256k1,256,REDUCE,50000,0.000034656,1442751615.882,0
opencl-kernel,Quadro RTX 6000,GPU,w8,secp256k1,256,REDUCE,50000,0.000334563,149448802.730,0
opencl-e2e,Quadro RTX 6000,GPU,w8,secp256k1,256,REDUCE,50000,0.001416050,35309490.211,0
opencl-kernel,Quadro RTX 6000,GPU,w16,secp256k1,256,REDUCE,50000,0.000179565,278451351.099,0
opencl-e2e,Quadro RTX 6000,GPU,w16,secp256k1,256,REDUCE,50000,0.001226222,40775645.997,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,secp256k1,256,MODMUL,3125,0.000521841,5988419.516,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,secp256k1,256,MODMUL,3125,0.007207070,433602.032,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,secp256k1,256,MODMUL,3125,0.007543368,414271.170,0
library,Quadro RTX 6000,gpu,cgbn,secp256k1,256,MODMUL,50000,0.000047200,1059322033.898,0
opencl-kernel,Quadro RTX 6000,GPU,w8,secp256k1,256,MODMUL,50000,0.000873320,57252795.291,0
opencl-e2e,Quadro RTX 6000,GPU,w8,secp256k1,256,MODMUL,50000,0.001977546,25283860.121,0
opencl-kernel,Quadro RTX 6000,GPU,w16,secp256k1,256,MODMUL,50000,0.000426563,117216116.397,0
opencl-e2e,Quadro RTX 6000,GPU,w16,secp256k1,256,MODMUL,50000,0.001515072,33001735.436,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,secp256k1,256,MODEXP,781,0.013538852,57685.835,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,secp256k1,256,MODEXP,781,0.007867822,99265.077,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,secp256k1,256,MODEXP,781,0.005335830,146368.977,0
library,Quadro RTX 6000,gpu,cgbn,secp256k1,256,MODEXP,50000,0.076185152,656295.862,0
opencl-kernel,Quadro RTX 6000,GPU,w8,secp256k1,256,MODEXP,50000,0.031672783,1578642.448,0
opencl-e2e,Quadro RTX 6000,GPU,w8,secp256k1,256,MODEXP,50000,0.032786254,1525029.369,0
opencl-kernel,Quadro RTX 6000,GPU,w16,secp256k1,256,MODEXP,50000,0.006381517,7835127.265,0
opencl-e2e,Quadro RTX 6000,GPU,w16,secp256k1,256,MODEXP,50000,0.008004583,6246421.255,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,secp256k1,256,EXPONENTIATION,781,0.005858392,133313.043,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,secp256k1,256,EXPONENTIATION,781,0.012699474,61498.612,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,secp256k1,256,EXPONENTIATION,781,0.007401038,105525.739,0
opencl-kernel,Quadro RTX 6000,GPU,w8,secp256k1,256,EXPONENTIATION,50000,0.047636244,1049620.947,0
opencl-e2e,Quadro RTX 6000,GPU,w8,secp256k1,256,EXPONENTIATION,50000,0.048808474,1024412.275,0
opencl-kernel,Quadro RTX 6000,GPU,w16,secp256k1,256,EXPONENTIATION,50000,0.019666700,2542368.612,0
opencl-e2e,Quadro RTX 6000,GPU,w16,secp256k1,256,EXPONENTIATION,50000,0.020701161,2415323.509,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,secp256k1,256,DIVIDE,6250,0.000474876,13161336.278,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,secp256k1,256,DIVIDE,6250,0.005077040,1231032.298,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,secp256k1,256,DIVIDE,6250,0.008023659,778946.386,0
library,Quadro RTX 6000,gpu,cgbn,secp256k1,256,DIVIDE,50000,0.000032928,1518464528.669,0
opencl-kernel,Quadro RTX 6000,GPU,w8,secp256k1,256,DIVIDE,50000,0.000519216,96299024.222,0
opencl-e2e,Quadro RTX 6000,GPU,w8,secp256k1,256,DIVIDE,50000,0.001853291,26979034.310,0
opencl-kernel,Quadro RTX 6000,GPU,w16,secp256k1,256,DIVIDE,50000,0.000335164,149180535.734,0
opencl-e2e,Quadro RTX 6000,GPU,w16,secp256k1,256,DIVIDE,50000,0.001619691,30870085.158,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,secp256k1,256,ISQRT,1562,0.000239536,6520935.961,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,secp256k1,256,ISQRT,1562,0.007520288,207704.805,0
opencl-kernel,Quadro RTX 6000,GPU,w8,secp256k1,256,ISQRT,50000,0.006182509,8087332.132,0
opencl-e2e,Quadro RTX 6000,GPU,w8,secp256k1,256,ISQRT,50000,0.007315507,6834796.128,0
opencl-kernel,Quadro RTX 6000,GPU,w16,secp256k1,256,ISQRT,50000,0.004952040,10096849.853,0
opencl-e2e,Quadro RTX 6000,GPU,w16,secp256k1,256,ISQRT,50000,0.006510044,7680440.211,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,secp256k1,256,MODMUL_R2,50000,0.008531127,5860890.275,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,secp256k1,256,MODMUL_R2,50000,0.005645644,8856385.686,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,secp256k1,256,MODMUL_R2,50000,0.007598994,6579818.069,0
library,Quadro RTX 6000,gpu,cgbn,secp256k1,256,MODMUL_R2,50000,0.000036864,1356336805.556,0
opencl-kernel,Quadro RTX 6000,GPU,w8,secp256k1,256,MODMUL_R2,50000,0.000183448,272556510.438,0
opencl-e2e,Quadro RTX 6000,GPU,w8,secp256k1,256,MODMUL_R2,50000,0.001274947,39217314.115,0
opencl-kernel,Quadro RTX 6000,GPU,w16,secp256k1,256,MODMUL_R2,50000,0.000061650,811032255.725,0
opencl-e2e,Quadro RTX 6000,GPU,w16,secp256k1,256,MODMUL_R2,50000,0.001078708,46351742.537,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,rsa256(composite),256,ADD,50000,0.001505516,33211194.307,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,rsa256(composite),256,ADD,50000,0.005862221,8529190.380,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,rsa256(composite),256,ADD,50000,0.006147156,8133843.318,0
library,Quadro RTX 6000,gpu,cgbn,rsa256(composite),256,ADD,50000,0.000032832,1522904483.431,0
opencl-kernel,Quadro RTX 6000,GPU,w8,rsa256(composite),256,ADD,50000,0.000051677,967544175.317,0
opencl-e2e,Quadro RTX 6000,GPU,w8,rsa256(composite),256,ADD,50000,0.001114249,44873264.149,0
opencl-kernel,Quadro RTX 6000,GPU,w16,rsa256(composite),256,ADD,50000,0.000030290,1650691526.258,0
opencl-e2e,Quadro RTX 6000,GPU,w16,rsa256(composite),256,ADD,50000,0.001063814,47000680.404,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACT,50000,0.001245866,40132739.894,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACT,50000,0.006204201,8059055.588,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACT,50000,0.004477825,11166135.305,0
library,Quadro RTX 6000,gpu,cgbn,rsa256(composite),256,SUBTRACT,50000,0.000032896,1519941634.241,0
opencl-kernel,Quadro RTX 6000,GPU,w8,rsa256(composite),256,SUBTRACT,50000,0.000051267,975277779.393,0
opencl-e2e,Quadro RTX 6000,GPU,w8,rsa256(composite),256,SUBTRACT,50000,0.001119548,44660862.852,0
opencl-kernel,Quadro RTX 6000,GPU,w16,rsa256(composite),256,SUBTRACT,50000,0.000030324,1648866437.346,0
opencl-e2e,Quadro RTX 6000,GPU,w16,rsa256(composite),256,SUBTRACT,50000,0.001054855,47399876.748,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,rsa256(composite),256,ADDMOD,50000,0.003946744,12668670.985,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,rsa256(composite),256,ADDMOD,50000,0.004321018,11571346.681,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,rsa256(composite),256,ADDMOD,50000,0.004754407,10516557.904,0
library,Quadro RTX 6000,gpu,cgbn,rsa256(composite),256,ADDMOD,50000,0.000033472,1493785850.860,0
opencl-kernel,Quadro RTX 6000,GPU,w8,rsa256(composite),256,ADDMOD,50000,0.000070333,710898983.051,0
opencl-e2e,Quadro RTX 6000,GPU,w8,rsa256(composite),256,ADDMOD,50000,0.001111910,44967678.473,0
opencl-kernel,Quadro RTX 6000,GPU,w16,rsa256(composite),256,ADDMOD,50000,0.000042737,1169959274.756,0
opencl-e2e,Quadro RTX 6000,GPU,w16,rsa256(composite),256,ADDMOD,50000,0.001057230,47293401.645,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACTMOD,50000,0.003587632,13936771.343,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACTMOD,50000,0.003919262,12757502.407,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACTMOD,50000,0.008105438,6168697.947,0
library,Quadro RTX 6000,gpu,cgbn,rsa256(composite),256,SUBTRACTMOD,50000,0.000032960,1516990291.262,0
opencl-kernel,Quadro RTX 6000,GPU,w8,rsa256(composite),256,SUBTRACTMOD,50000,0.000075717,660357825.338,0
opencl-e2e,Quadro RTX 6000,GPU,w8,rsa256(composite),256,SUBTRACTMOD,50000,0.001118824,44689785.970,0
opencl-kernel,Quadro RTX 6000,GPU,w16,rsa256(composite),256,SUBTRACTMOD,50000,0.000042470,1177296855.401,0
opencl-e2e,Quadro RTX 6000,GPU,w16,rsa256(composite),256,SUBTRACTMOD,50000,0.001063524,47013521.783,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001872571,26701254.614,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.005985957,8352883.760,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.005984955,8354282.342,0
opencl-kernel,Quadro RTX 6000,GPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001966847,25421396.420,0
opencl-e2e,Quadro RTX 6000,GPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.003293602,15180948.729,0
opencl-kernel,Quadro RTX 6000,GPU,w16,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000518383,96453695.047,0
opencl-e2e,Quadro RTX 6000,GPU,w16,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001811361,27603554.259,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001171757,42670980.254,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.005953321,8398673.285,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000278566,179490656.218,0
library,Quadro RTX 6000,gpu,cgbn,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000032768,1525878906.250,0
opencl-kernel,Quadro RTX 6000,GPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000593737,84212403.062,0
opencl-e2e,Quadro RTX 6000,GPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001905486,26240025.024,0
opencl-kernel,Quadro RTX 6000,GPU,w16,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000174986,285736820.480,0
opencl-e2e,Quadro RTX 6000,GPU,w16,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001430674,34948574.315,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.014528222,3441577.275,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000723593,69099624.172,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000244820,204231271.255,0
library,Quadro RTX 6000,gpu,cgbn,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000033184,1506750241.080,0
opencl-kernel,Quadro RTX 6000,GPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000151580,329858386.070,0
opencl-e2e,Quadro RTX 6000,GPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.001195425,41826122.916,0
opencl-kernel,Quadro RTX 6000,GPU,w16,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000045754,1092800260.544,0
opencl-e2e,Quadro RTX 6000,GPU,w16,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.001063600,47010146.126,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,rsa256(composite),256,COMPARE,50000,0.000437824,114201125.694,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,rsa256(composite),256,COMPARE,50000,0.005969189,8376347.122,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,rsa256(composite),256,COMPARE,50000,0.000127921,390866601.628,0
library,Quadro RTX 6000,gpu,cgbn,rsa256(composite),256,COMPARE,50000,0.000033344,1499520153.551,0
opencl-kernel,Quadro RTX 6000,GPU,w8,rsa256(composite),256,COMPARE,50000,0.000050014,999722379.055,0
opencl-e2e,Quadro RTX 6000,GPU,w8,rsa256(composite),256,COMPARE,50000,0.001125423,44427730.232,0
opencl-kernel,Quadro RTX 6000,GPU,w16,rsa256(composite),256,COMPARE,50000,0.000031332,1595835301.112,0
opencl-e2e,Quadro RTX 6000,GPU,w16,rsa256(composite),256,COMPARE,50000,0.001059262,47202671.417,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,rsa256(composite),256,REDUCE,6250,0.000370556,16866524.246,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,rsa256(composite),256,REDUCE,6250,0.000075810,82443321.867,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,rsa256(composite),256,REDUCE,6250,0.000160163,39022680.173,0
library,Quadro RTX 6000,gpu,cgbn,rsa256(composite),256,REDUCE,50000,0.000033120,1509661835.749,0
opencl-kernel,Quadro RTX 6000,GPU,w8,rsa256(composite),256,REDUCE,50000,0.000273930,182528443.885,0
opencl-e2e,Quadro RTX 6000,GPU,w8,rsa256(composite),256,REDUCE,50000,0.001350552,37021903.480,0
opencl-kernel,Quadro RTX 6000,GPU,w16,rsa256(composite),256,REDUCE,50000,0.000169622,294773465.107,0
opencl-e2e,Quadro RTX 6000,GPU,w16,rsa256(composite),256,REDUCE,50000,0.001214797,41159151.105,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,rsa256(composite),256,MODMUL,3125,0.000523442,5970093.338,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,rsa256(composite),256,MODMUL,3125,0.000078943,39585710.915,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,rsa256(composite),256,MODMUL,3125,0.000173893,17970839.135,0
library,Quadro RTX 6000,gpu,cgbn,rsa256(composite),256,MODMUL,50000,0.000038912,1284950657.895,0
opencl-kernel,Quadro RTX 6000,GPU,w8,rsa256(composite),256,MODMUL,50000,0.000714479,69981061.726,0
opencl-e2e,Quadro RTX 6000,GPU,w8,rsa256(composite),256,MODMUL,50000,0.001826253,27378466.070,0
opencl-kernel,Quadro RTX 6000,GPU,w16,rsa256(composite),256,MODMUL,50000,0.000403572,123893668.596,0
opencl-e2e,Quadro RTX 6000,GPU,w16,rsa256(composite),256,MODMUL,50000,0.001490334,33549525.382,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,rsa256(composite),256,MODEXP,781,0.012940142,60354.823,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,rsa256(composite),256,MODEXP,781,0.000606947,1286768.785,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,rsa256(composite),256,MODEXP,781,0.001115950,699852.254,0
library,Quadro RTX 6000,gpu,cgbn,rsa256(composite),256,MODEXP,50000,0.073687039,678545.382,0
opencl-kernel,Quadro RTX 6000,GPU,w8,rsa256(composite),256,MODEXP,50000,0.031585034,1583028.216,0
opencl-e2e,Quadro RTX 6000,GPU,w8,rsa256(composite),256,MODEXP,50000,0.032738388,1527259.089,0
opencl-kernel,Quadro RTX 6000,GPU,w16,rsa256(composite),256,MODEXP,50000,0.006465536,7733311.362,0
opencl-e2e,Quadro RTX 6000,GPU,w16,rsa256(composite),256,MODEXP,50000,0.007772837,6432657.866,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,rsa256(composite),256,EXPONENTIATION,781,0.004763564,163952.869,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,rsa256(composite),256,EXPONENTIATION,781,0.000282556,2764055.626,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,rsa256(composite),256,EXPONENTIATION,781,0.004738379,164824.290,0
opencl-kernel,Quadro RTX 6000,GPU,w8,rsa256(composite),256,EXPONENTIATION,50000,0.047663124,1049029.008,0
opencl-e2e,Quadro RTX 6000,GPU,w8,rsa256(composite),256,EXPONENTIATION,50000,0.048876440,1022987.759,0
opencl-kernel,Quadro RTX 6000,GPU,w16,rsa256(composite),256,EXPONENTIATION,50000,0.019175041,2607556.291,0
opencl-e2e,Quadro RTX 6000,GPU,w16,rsa256(composite),256,EXPONENTIATION,50000,0.018760575,2665163.517,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,rsa256(composite),256,DIVIDE,6250,0.000481050,12992399.966,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,rsa256(composite),256,DIVIDE,6250,0.000106879,58477574.068,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,rsa256(composite),256,DIVIDE,6250,0.000167338,37349516.357,0
library,Quadro RTX 6000,gpu,cgbn,rsa256(composite),256,DIVIDE,50000,0.000033728,1482447817.837,0
opencl-kernel,Quadro RTX 6000,GPU,w8,rsa256(composite),256,DIVIDE,50000,0.000513295,97409926.226,0
opencl-e2e,Quadro RTX 6000,GPU,w8,rsa256(composite),256,DIVIDE,50000,0.001805937,27686459.871,0
opencl-kernel,Quadro RTX 6000,GPU,w16,rsa256(composite),256,DIVIDE,50000,0.000338925,147525242.500,0
opencl-e2e,Quadro RTX 6000,GPU,w16,rsa256(composite),256,DIVIDE,50000,0.001622954,30808013.157,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,rsa256(composite),256,ISQRT,1562,0.000238625,6545826.389,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,rsa256(composite),256,ISQRT,1562,0.000081291,19214819.434,0
opencl-kernel,Quadro RTX 6000,GPU,w8,rsa256(composite),256,ISQRT,50000,0.006151743,8127777.462,0
opencl-e2e,Quadro RTX 6000,GPU,w8,rsa256(composite),256,ISQRT,50000,0.007252129,6894527.306,0
opencl-kernel,Quadro RTX 6000,GPU,w16,rsa256(composite),256,ISQRT,50000,0.004956407,10087951.864,0
opencl-e2e,Quadro RTX 6000,GPU,w16,rsa256(composite),256,ISQRT,50000,0.006527219,7660230.222,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,rsa256(composite),256,MODMUL_R2,50000,0.008429579,5931494.018,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,rsa256(composite),256,MODMUL_R2,50000,0.000472507,105818648.270,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,rsa256(composite),256,MODMUL_R2,50000,0.004830115,10351721.383,0
library,Quadro RTX 6000,gpu,cgbn,rsa256(composite),256,MODMUL_R2,50000,0.000036864,1356336805.556,0
opencl-kernel,Quadro RTX 6000,GPU,w8,rsa256(composite),256,MODMUL_R2,50000,0.000182243,274358863.871,0
opencl-e2e,Quadro RTX 6000,GPU,w8,rsa256(composite),256,MODMUL_R2,50000,0.001236936,40422460.716,0
opencl-kernel,Quadro RTX 6000,GPU,w16,rsa256(composite),256,MODMUL_R2,50000,0.000058014,861861735.054,0
opencl-e2e,Quadro RTX 6000,GPU,w16,rsa256(composite),256,MODMUL_R2,50000,0.001081077,46250158.253,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,brainpoolP512r1,512,ADD,25000,0.000505690,49437448.156,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,brainpoolP512r1,512,ADD,25000,0.000099299,251763666.035,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,brainpoolP512r1,512,ADD,25000,0.000123031,203200096.893,0
library,Quadro RTX 6000,gpu,cgbn,brainpoolP512r1,512,ADD,50000,0.000040960,1220703125.000,0
opencl-kernel,Quadro RTX 6000,GPU,w8,brainpoolP512r1,512,ADD,50000,0.000133950,373272875.935,0
opencl-e2e,Quadro RTX 6000,GPU,w8,brainpoolP512r1,512,ADD,50000,0.002006819,24915046.273,0
opencl-kernel,Quadro RTX 6000,GPU,w16,brainpoolP512r1,512,ADD,50000,0.000077741,643159441.263,0
opencl-e2e,Quadro RTX 6000,GPU,w16,brainpoolP512r1,512,ADD,50000,0.001888383,26477677.058,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACT,25000,0.000426844,58569439.693,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACT,25000,0.000084404,296194836.033,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACT,25000,0.000117518,212733354.466,0
library,Quadro RTX 6000,gpu,cgbn,brainpoolP512r1,512,SUBTRACT,50000,0.000041344,1209365325.077,0
opencl-kernel,Quadro RTX 6000,GPU,w8,brainpoolP512r1,512,SUBTRACT,50000,0.000134403,372015821.057,0
opencl-e2e,Quadro RTX 6000,GPU,w8,brainpoolP512r1,512,SUBTRACT,50000,0.001969757,25383847.452,0
opencl-kernel,Quadro RTX 6000,GPU,w16,brainpoolP512r1,512,SUBTRACT,50000,0.000078112,640107439.908,0
opencl-e2e,Quadro RTX 6000,GPU,w16,brainpoolP512r1,512,SUBTRACT,50000,0.001900418,26310005.008,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,brainpoolP512r1,512,ADDMOD,25000,0.001423445,17563029.864,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,brainpoolP512r1,512,ADDMOD,25000,0.000157058,159176622.391,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,brainpoolP512r1,512,ADDMOD,25000,0.000438647,56993391.848,0
library,Quadro RTX 6000,gpu,cgbn,brainpoolP512r1,512,ADDMOD,50000,0.000039872,1254012841.091,0
opencl-kernel,Quadro RTX 6000,GPU,w8,brainpoolP512r1,512,ADDMOD,50000,0.000176832,282754124.884,0
opencl-e2e,Quadro RTX 6000,GPU,w8,brainpoolP512r1,512,ADDMOD,50000,0.002183501,22899008.239,0
opencl-kernel,Quadro RTX 6000,GPU,w16,brainpoolP512r1,512,ADDMOD,50000,0.000095205,525180396.377,0
opencl-e2e,Quadro RTX 6000,GPU,w16,brainpoolP512r1,512,ADDMOD,50000,0.001970308,25376744.397,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACTMOD,25000,0.001212493,20618682.791,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000155834,160426625.868,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000451455,55376518.011,0
library,Quadro RTX 6000,gpu,cgbn,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000040960,1220703125.000,0
opencl-kernel,Quadro RTX 6000,GPU,w8,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000202766,246590044.002,0
opencl-e2e,Quadro RTX 6000,GPU,w8,brainpoolP512r1,512,SUBTRACTMOD,50000,0.002091330,23908235.852,0
opencl-kernel,Quadro RTX 6000,GPU,w16,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000102527,487674325.994,0
opencl-e2e,Quadro RTX 6000,GPU,w16,brainpoolP512r1,512,SUBTRACTMOD,50000,0.001901995,26288181.478,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.001298344,19255301.036,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.000172157,145216419.622,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.000193452,129230714.719,0
opencl-kernel,Quadro RTX 6000,GPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.013827095,3616088.607,0
opencl-e2e,Quadro RTX 6000,GPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.016105585,3104513.162,0
opencl-kernel,Quadro RTX 6000,GPU,w16,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.003697824,13521466.573,0
opencl-e2e,Quadro RTX 6000,GPU,w16,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.006303489,7932114.794,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.001407461,17762478.478,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000164265,152193275.805,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000186896,133764267.134,0
library,Quadro RTX 6000,gpu,cgbn,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000040928,1221657544.957,0
opencl-kernel,Quadro RTX 6000,GPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.003971621,12589316.815,0
opencl-e2e,Quadro RTX 6000,GPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.006273601,7969904.107,0
opencl-kernel,Quadro RTX 6000,GPU,w16,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.001079097,46335020.782,0
opencl-e2e,Quadro RTX 6000,GPU,w16,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.003404755,14685343.840,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.009967761,2508085.789,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000730343,34230484.060,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000243934,102486792.251,0
library,Quadro RTX 6000,gpu,cgbn,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000040928,1221657544.957,0
opencl-kernel,Quadro RTX 6000,GPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000563039,88803871.920,0
opencl-e2e,Quadro RTX 6000,GPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.002432827,20552222.547,0
opencl-kernel,Quadro RTX 6000,GPU,w16,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000170225,293728409.327,0
opencl-e2e,Quadro RTX 6000,GPU,w16,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.002055308,24327255.887,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,brainpoolP512r1,512,COMPARE,25000,0.000224603,111307338.514,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,brainpoolP512r1,512,COMPARE,25000,0.000068394,365526642.882,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,brainpoolP512r1,512,COMPARE,25000,0.000102485,243939092.347,0
library,Quadro RTX 6000,gpu,cgbn,brainpoolP512r1,512,COMPARE,50000,0.000040960,1220703125.000,0
opencl-kernel,Quadro RTX 6000,GPU,w8,brainpoolP512r1,512,COMPARE,50000,0.000113780,439445782.107,0
opencl-e2e,Quadro RTX 6000,GPU,w8,brainpoolP512r1,512,COMPARE,50000,0.001956776,25552237.715,0
opencl-kernel,Quadro RTX 6000,GPU,w16,brainpoolP512r1,512,COMPARE,50000,0.000068437,730595656.197,0
opencl-e2e,Quadro RTX 6000,GPU,w16,brainpoolP512r1,512,COMPARE,50000,0.001939334,25782047.155,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,brainpoolP512r1,512,REDUCE,3125,0.000118105,26459564.401,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,brainpoolP512r1,512,REDUCE,3125,0.000064574,48393954.079,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,brainpoolP512r1,512,REDUCE,3125,0.000123903,25221310.884,0
library,Quadro RTX 6000,gpu,cgbn,brainpoolP512r1,512,REDUCE,50000,0.000040864,1223570869.225,0
opencl-kernel,Quadro RTX 6000,GPU,w8,brainpoolP512r1,512,REDUCE,50000,0.000791941,63136029.353,0
opencl-e2e,Quadro RTX 6000,GPU,w8,brainpoolP512r1,512,REDUCE,50000,0.002672946,18705951.595,0
opencl-kernel,Quadro RTX 6000,GPU,w16,brainpoolP512r1,512,REDUCE,50000,0.000621401,80463369.423,0
opencl-e2e,Quadro RTX 6000,GPU,w16,brainpoolP512r1,512,REDUCE,50000,0.002505042,19959747.904,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL,1562,0.000311119,5020579.201,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL,1562,0.000078315,19945115.104,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL,1562,0.000155227,10062665.617,0
library,Quadro RTX 6000,gpu,cgbn,brainpoolP512r1,512,MODMUL,50000,0.000113568,440264863.342,0
opencl-kernel,Quadro RTX 6000,GPU,w8,brainpoolP512r1,512,MODMUL,50000,0.002201706,22709660.652,0
opencl-e2e,Quadro RTX 6000,GPU,w8,brainpoolP512r1,512,MODMUL,50000,0.004092166,12218467.568,0
opencl-kernel,Quadro RTX 6000,GPU,w16,brainpoolP512r1,512,MODMUL,50000,0.001908051,26204752.352,0
opencl-e2e,Quadro RTX 6000,GPU,w16,brainpoolP512r1,512,MODMUL,50000,0.003823301,13077704.922,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,brainpoolP512r1,512,MODEXP,390,0.021464465,18169.565,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,brainpoolP512r1,512,MODEXP,390,0.001839718,211989.029,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,brainpoolP512r1,512,MODEXP,390,0.001755847,222115.067,0
library,Quadro RTX 6000,gpu,cgbn,brainpoolP512r1,512,MODEXP,50000,0.674238443,74157.741,0
opencl-kernel,Quadro RTX 6000,GPU,w8,brainpoolP512r1,512,MODEXP,50000,0.461431239,108358.507,0
opencl-e2e,Quadro RTX 6000,GPU,w8,brainpoolP512r1,512,MODEXP,50000,0.463102965,107967.350,0
opencl-kernel,Quadro RTX 6000,GPU,w16,brainpoolP512r1,512,MODEXP,50000,0.044198086,1131270.695,0
opencl-e2e,Quadro RTX 6000,GPU,w16,brainpoolP512r1,512,MODEXP,50000,0.046244949,1081199.157,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,brainpoolP512r1,512,EXPONENTIATION,390,0.005267378,74040.634,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,brainpoolP512r1,512,EXPONENTIATION,390,0.000738392,528174.985,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,brainpoolP512r1,512,EXPONENTIATION,390,0.008676311,44949.979,0
opencl-kernel,Quadro RTX 6000,GPU,w8,brainpoolP512r1,512,EXPONENTIATION,50000,0.404858837,123499.836,0
opencl-e2e,Quadro RTX 6000,GPU,w8,brainpoolP512r1,512,EXPONENTIATION,50000,0.406460391,123013.216,0
opencl-kernel,Quadro RTX 6000,GPU,w16,brainpoolP512r1,512,EXPONENTIATION,50000,0.161273964,310031.445,0
opencl-e2e,Quadro RTX 6000,GPU,w16,brainpoolP512r1,512,EXPONENTIATION,50000,0.169829650,294412.666,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,brainpoolP512r1,512,DIVIDE,3125,0.000483386,6464809.877,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,brainpoolP512r1,512,DIVIDE,3125,0.000082113,38057381.363,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,brainpoolP512r1,512,DIVIDE,3125,0.000151105,20680952.616,0
library,Quadro RTX 6000,gpu,cgbn,brainpoolP512r1,512,DIVIDE,50000,0.000058624,852893013.100,0
opencl-kernel,Quadro RTX 6000,GPU,w8,brainpoolP512r1,512,DIVIDE,50000,0.001546312,32334998.374,0
opencl-e2e,Quadro RTX 6000,GPU,w8,brainpoolP512r1,512,DIVIDE,50000,0.003838906,13024543.752,0
opencl-kernel,Quadro RTX 6000,GPU,w16,brainpoolP512r1,512,DIVIDE,50000,0.001430653,34949074.831,0
opencl-e2e,Quadro RTX 6000,GPU,w16,brainpoolP512r1,512,DIVIDE,50000,0.003705321,13494107.966,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,brainpoolP512r1,512,ISQRT,781,0.000377841,2067005.414,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,brainpoolP512r1,512,ISQRT,781,0.000080002,9762198.372,0
opencl-kernel,Quadro RTX 6000,GPU,w8,brainpoolP512r1,512,ISQRT,50000,0.036745481,1360711.532,0
opencl-e2e,Quadro RTX 6000,GPU,w8,brainpoolP512r1,512,ISQRT,50000,0.038654529,1293509.502,0
opencl-kernel,Quadro RTX 6000,GPU,w16,brainpoolP512r1,512,ISQRT,50000,0.042029435,1189642.455,0
opencl-e2e,Quadro RTX 6000,GPU,w16,brainpoolP512r1,512,ISQRT,50000,0.043414865,1151679.278,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL_R2,25000,0.008243611,3032651.559,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL_R2,25000,0.000412324,60631865.019,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL_R2,25000,0.001170313,21361806.887,0
library,Quadro RTX 6000,gpu,cgbn,brainpoolP512r1,512,MODMUL_R2,50000,0.000048992,1020574787.720,0
opencl-kernel,Quadro RTX 6000,GPU,w8,brainpoolP512r1,512,MODMUL_R2,50000,0.000815690,61297829.741,0
opencl-e2e,Quadro RTX 6000,GPU,w8,brainpoolP512r1,512,MODMUL_R2,50000,0.002704194,18489798.960,0
opencl-kernel,Quadro RTX 6000,GPU,w16,brainpoolP512r1,512,MODMUL_R2,50000,0.000184067,271640817.648,0
opencl-e2e,Quadro RTX 6000,GPU,w16,brainpoolP512r1,512,MODMUL_R2,50000,0.001995681,25054106.175,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,p1024,1024,ADD,12500,0.000515612,24243042.877,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,p1024,1024,ADD,12500,0.000079637,156961440.767,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,p1024,1024,ADD,12500,0.000110546,113074970.092,0
library,Quadro RTX 6000,gpu,cgbn,p1024,1024,ADD,50000,0.000063488,787550403.226,0
opencl-kernel,Quadro RTX 6000,GPU,w8,p1024,1024,ADD,50000,0.000869472,57506192.453,0
opencl-e2e,Quadro RTX 6000,GPU,w8,p1024,1024,ADD,50000,0.004251018,11761888.103,0
opencl-kernel,Quadro RTX 6000,GPU,w16,p1024,1024,ADD,50000,0.000488965,102256841.047,0
opencl-e2e,Quadro RTX 6000,GPU,w16,p1024,1024,ADD,50000,0.003899155,12823290.593,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,p1024,1024,SUBTRACT,12500,0.000419097,29826029.218,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,p1024,1024,SUBTRACT,12500,0.000080889,154532581.113,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,p1024,1024,SUBTRACT,12500,0.000118256,105703225.806,0
library,Quadro RTX 6000,gpu,cgbn,p1024,1024,SUBTRACT,50000,0.000063488,787550403.226,0
opencl-kernel,Quadro RTX 6000,GPU,w8,p1024,1024,SUBTRACT,50000,0.000899417,55591546.396,0
opencl-e2e,Quadro RTX 6000,GPU,w8,p1024,1024,SUBTRACT,50000,0.004298268,11632592.785,0
opencl-kernel,Quadro RTX 6000,GPU,w16,p1024,1024,SUBTRACT,50000,0.000512432,97573863.545,0
opencl-e2e,Quadro RTX 6000,GPU,w16,p1024,1024,SUBTRACT,50000,0.003923045,12745200.219,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,p1024,1024,ADDMOD,12500,0.001597717,7823661.410,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,p1024,1024,ADDMOD,12500,0.000156390,79928614.476,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,p1024,1024,ADDMOD,12500,0.000333803,37447262.135,0
library,Quadro RTX 6000,gpu,cgbn,p1024,1024,ADDMOD,50000,0.000062976,793953252.033,0
opencl-kernel,Quadro RTX 6000,GPU,w8,p1024,1024,ADDMOD,50000,0.000859780,58154395.558,0
opencl-e2e,Quadro RTX 6000,GPU,w8,p1024,1024,ADDMOD,50000,0.004287167,11662714.699,0
opencl-kernel,Quadro RTX 6000,GPU,w16,p1024,1024,ADDMOD,50000,0.000457713,109238667.985,0
opencl-e2e,Quadro RTX 6000,GPU,w16,p1024,1024,ADDMOD,50000,0.003917342,12763756.480,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,p1024,1024,SUBTRACTMOD,12500,0.001147682,10891520.181,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,p1024,1024,SUBTRACTMOD,12500,0.000115562,108166828.922,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,p1024,1024,SUBTRACTMOD,12500,0.000325739,38374236.048,0
library,Quadro RTX 6000,gpu,cgbn,p1024,1024,SUBTRACTMOD,50000,0.000062336,802104722.793,0
opencl-kernel,Quadro RTX 6000,GPU,w8,p1024,1024,SUBTRACTMOD,50000,0.000993598,50322149.694,0
opencl-e2e,Quadro RTX 6000,GPU,w8,p1024,1024,SUBTRACTMOD,50000,0.004409170,11340003.337,0
opencl-kernel,Quadro RTX 6000,GPU,w16,p1024,1024,SUBTRACTMOD,50000,0.000476100,105020052.894,0
opencl-e2e,Quadro RTX 6000,GPU,w16,p1024,1024,SUBTRACTMOD,50000,0.003907835,12794808.015,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.002177916,5739430.409,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.000206444,60548986.773,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.000279343,44747893.926,0
opencl-kernel,Quadro RTX 6000,GPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.268347647,186325.465,0
opencl-e2e,Quadro RTX 6000,GPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.273927748,182529.884,0
opencl-kernel,Quadro RTX 6000,GPU,w16,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.070080148,713468.812,0
opencl-e2e,Quadro RTX 6000,GPU,w16,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.075161457,665234.576,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.002172502,5753735.276,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.000202164,61830972.212,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.000285752,43744207.754,0
library,Quadro RTX 6000,gpu,cgbn,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000063488,787550403.226,0
opencl-kernel,Quadro RTX 6000,GPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.030825313,1622043.536,0
opencl-e2e,Quadro RTX 6000,GPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.035022078,1427670.850,0
opencl-kernel,Quadro RTX 6000,GPU,w16,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.009193357,5438709.664,0
opencl-e2e,Quadro RTX 6000,GPU,w16,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.014338352,3487151.156,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.014700118,850333.284,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.000996159,12548192.719,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.000316788,39458626.245,0
library,Quadro RTX 6000,gpu,cgbn,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000068032,734948259.643,0
opencl-kernel,Quadro RTX 6000,GPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.003051531,16385220.091,0
opencl-e2e,Quadro RTX 6000,GPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.006447753,7754639.559,0
opencl-kernel,Quadro RTX 6000,GPU,w16,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000608236,82204987.995,0
opencl-e2e,Quadro RTX 6000,GPU,w16,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.004050829,12343153.949,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,p1024,1024,COMPARE,12500,0.000142973,87429146.148,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,p1024,1024,COMPARE,12500,0.000065055,192145862.681,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,p1024,1024,COMPARE,12500,0.000098733,126603776.860,0
library,Quadro RTX 6000,gpu,cgbn,p1024,1024,COMPARE,50000,0.000063488,787550403.226,0
opencl-kernel,Quadro RTX 6000,GPU,w8,p1024,1024,COMPARE,50000,0.000283303,176489645.424,0
opencl-e2e,Quadro RTX 6000,GPU,w8,p1024,1024,COMPARE,50000,0.003666213,13638051.845,0
opencl-kernel,Quadro RTX 6000,GPU,w16,p1024,1024,COMPARE,50000,0.000154108,324448191.839,0
opencl-e2e,Quadro RTX 6000,GPU,w16,p1024,1024,COMPARE,50000,0.003589518,13929445.366,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,p1024,1024,REDUCE,1562,0.000045402,34403789.315,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,p1024,1024,REDUCE,1562,0.000069709,22407277.610,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,p1024,1024,REDUCE,1562,0.000113539,13757339.139,0
library,Quadro RTX 6000,gpu,cgbn,p1024,1024,REDUCE,50000,0.000061440,813802083.333,0
opencl-kernel,Quadro RTX 6000,GPU,w8,p1024,1024,REDUCE,50000,0.004176635,11971358.873,0
opencl-e2e,Quadro RTX 6000,GPU,w8,p1024,1024,REDUCE,50000,0.007562557,6611520.210,0
opencl-kernel,Quadro RTX 6000,GPU,w16,p1024,1024,REDUCE,50000,0.002272598,22001247.119,0
opencl-e2e,Quadro RTX 6000,GPU,w16,p1024,1024,REDUCE,50000,0.005905233,8467065.950,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,p1024,1024,MODMUL,781,0.000404814,1929279.871,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,p1024,1024,MODMUL,781,0.000085982,9083342.698,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,p1024,1024,MODMUL,781,0.000180479,4327369.933,0
library,Quadro RTX 6000,gpu,cgbn,p1024,1024,MODMUL,50000,0.000327232,152796792.490,0
opencl-kernel,Quadro RTX 6000,GPU,w8,p1024,1024,MODMUL,50000,0.017451508,2865081.922,0
opencl-e2e,Quadro RTX 6000,GPU,w8,p1024,1024,MODMUL,50000,0.020809913,2402701.053,0
opencl-kernel,Quadro RTX 6000,GPU,w16,p1024,1024,MODMUL,50000,0.006542891,7641881.545,0
opencl-e2e,Quadro RTX 6000,GPU,w16,p1024,1024,MODMUL,50000,0.010741368,4654900.595,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,p1024,1024,MODEXP,195,0.068770196,2835.531,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,p1024,1024,MODEXP,195,0.006813556,28619.417,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,p1024,1024,MODEXP,195,0.005274625,36969.449,0
library,Quadro RTX 6000,gpu,cgbn,p1024,1024,MODEXP,50000,1.544298410,32377.162,0
opencl-kernel,Quadro RTX 6000,GPU,w8,p1024,1024,MODEXP,50000,4.131514436,12102.100,0
opencl-e2e,Quadro RTX 6000,GPU,w8,p1024,1024,MODEXP,50000,4.166463666,12000.585,0
opencl-kernel,Quadro RTX 6000,GPU,w16,p1024,1024,MODEXP,50000,0.676373696,73923.632,0
opencl-e2e,Quadro RTX 6000,GPU,w16,p1024,1024,MODEXP,50000,0.720918896,69355.929,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,p1024,1024,EXPONENTIATION,195,0.009914292,19668.575,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,p1024,1024,EXPONENTIATION,195,0.001372343,142092.764,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,p1024,1024,EXPONENTIATION,195,0.098882837,1972.031,0
opencl-kernel,Quadro RTX 6000,GPU,w8,p1024,1024,EXPONENTIATION,50000,4.778523376,10463.483,0
opencl-e2e,Quadro RTX 6000,GPU,w8,p1024,1024,EXPONENTIATION,50000,4.805269400,10405.244,0
opencl-kernel,Quadro RTX 6000,GPU,w16,p1024,1024,EXPONENTIATION,50000,1.112362232,44949.387,0
opencl-e2e,Quadro RTX 6000,GPU,w16,p1024,1024,EXPONENTIATION,50000,1.114419019,44866.427,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,p1024,1024,DIVIDE,1562,0.000260236,6002250.074,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,p1024,1024,DIVIDE,1562,0.000067990,22973874.433,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,p1024,1024,DIVIDE,1562,0.000134714,11594938.949,0
library,Quadro RTX 6000,gpu,cgbn,p1024,1024,DIVIDE,50000,0.000082496,606089992.242,0
opencl-kernel,Quadro RTX 6000,GPU,w8,p1024,1024,DIVIDE,50000,0.180345442,277245.709,0
opencl-e2e,Quadro RTX 6000,GPU,w8,p1024,1024,DIVIDE,50000,0.185413131,269668.063,0
opencl-kernel,Quadro RTX 6000,GPU,w16,p1024,1024,DIVIDE,50000,0.016517444,3027102.817,0
opencl-e2e,Quadro RTX 6000,GPU,w16,p1024,1024,DIVIDE,50000,0.020769851,2407335.485,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,p1024,1024,ISQRT,390,0.000372875,1045925.568,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,p1024,1024,ISQRT,390,0.000083163,4689564.049,0
opencl-kernel,Quadro RTX 6000,GPU,w8,p1024,1024,ISQRT,50000,1.996044235,25049.545,0
opencl-e2e,Quadro RTX 6000,GPU,w8,p1024,1024,ISQRT,50000,1.975367986,25311.740,0
opencl-kernel,Quadro RTX 6000,GPU,w16,p1024,1024,ISQRT,50000,0.308983272,161821.058,0
opencl-e2e,Quadro RTX 6000,GPU,w16,p1024,1024,ISQRT,50000,0.312931506,159779.374,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,p1024,1024,MODMUL_R2,12500,0.010793075,1158150.021,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,p1024,1024,MODMUL_R2,12500,0.000483979,25827591.462,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,p1024,1024,MODMUL_R2,12500,0.001360195,9189860.445,0
library,Quadro RTX 6000,gpu,cgbn,p1024,1024,MODMUL_R2,50000,0.000122816,407113079.729,0
opencl-kernel,Quadro RTX 6000,GPU,w8,p1024,1024,MODMUL_R2,50000,0.004408233,11342413.502,0
opencl-e2e,Quadro RTX 6000,GPU,w8,p1024,1024,MODMUL_R2,50000,0.008206228,6092933.471,0
opencl-kernel,Quadro RTX 6000,GPU,w16,p1024,1024,MODMUL_R2,50000,0.000717292,69706656.834,0
opencl-e2e,Quadro RTX 6000,GPU,w16,p1024,1024,MODMUL_R2,50000,0.004142627,12069635.657,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,p2048,2048,ADD,6250,0.000351148,17798776.794,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,p2048,2048,ADD,6250,0.000077225,80932059.817,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,p2048,2048,ADD,6250,0.000128489,48642300.890,0
library,Quadro RTX 6000,gpu,cgbn,p2048,2048,ADD,50000,0.000121280,412269129.288,0
opencl-kernel,Quadro RTX 6000,GPU,w8,p2048,2048,ADD,50000,0.002602460,19212592.176,0
opencl-e2e,Quadro RTX 6000,GPU,w8,p2048,2048,ADD,50000,0.009018060,5544429.844,0
opencl-kernel,Quadro RTX 6000,GPU,w16,p2048,2048,ADD,50000,0.001261134,39646866.923,0
opencl-e2e,Quadro RTX 6000,GPU,w16,p2048,2048,ADD,50000,0.007655125,6531571.891,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,p2048,2048,SUBTRACT,6250,0.000315504,19809565.192,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,p2048,2048,SUBTRACT,6250,0.000077313,80840417.279,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,p2048,2048,SUBTRACT,6250,0.000127945,48849078.468,0
library,Quadro RTX 6000,gpu,cgbn,p2048,2048,SUBTRACT,50000,0.000121504,411509086.121,0
opencl-kernel,Quadro RTX 6000,GPU,w8,p2048,2048,SUBTRACT,50000,0.002559174,19537555.042,0
opencl-e2e,Quadro RTX 6000,GPU,w8,p2048,2048,SUBTRACT,50000,0.009494465,5266226.251,0
opencl-kernel,Quadro RTX 6000,GPU,w16,p2048,2048,SUBTRACT,50000,0.001277380,39142629.501,0
opencl-e2e,Quadro RTX 6000,GPU,w16,p2048,2048,SUBTRACT,50000,0.007761912,6441711.436,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,p2048,2048,ADDMOD,6250,0.000940599,6644705.137,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,p2048,2048,ADDMOD,6250,0.000108849,57418857.593,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,p2048,2048,ADDMOD,6250,0.000276932,22568677.065,0
library,Quadro RTX 6000,gpu,cgbn,p2048,2048,ADDMOD,50000,0.000120288,415669060.920,0
opencl-kernel,Quadro RTX 6000,GPU,w8,p2048,2048,ADDMOD,50000,0.002432508,20554913.652,0
opencl-e2e,Quadro RTX 6000,GPU,w8,p2048,2048,ADDMOD,50000,0.009331670,5358097.943,0
opencl-kernel,Quadro RTX 6000,GPU,w16,p2048,2048,ADDMOD,50000,0.001258837,39719199.270,0
opencl-e2e,Quadro RTX 6000,GPU,w16,p2048,2048,ADDMOD,50000,0.007635428,6548421.695,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,p2048,2048,SUBTRACTMOD,6250,0.000469005,13326091.463,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,p2048,2048,SUBTRACTMOD,6250,0.004097404,1525356.070,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,p2048,2048,SUBTRACTMOD,6250,0.000275359,22697678.446,0
library,Quadro RTX 6000,gpu,cgbn,p2048,2048,SUBTRACTMOD,50000,0.000118784,420932112.069,0
opencl-kernel,Quadro RTX 6000,GPU,w8,p2048,2048,SUBTRACTMOD,50000,0.002537094,19707586.272,0
opencl-e2e,Quadro RTX 6000,GPU,w8,p2048,2048,SUBTRACTMOD,50000,0.009479722,5274416.271,0
opencl-kernel,Quadro RTX 6000,GPU,w16,p2048,2048,SUBTRACTMOD,50000,0.001234433,40504436.340,0
opencl-e2e,Quadro RTX 6000,GPU,w16,p2048,2048,SUBTRACTMOD,50000,0.007635003,6548785.939,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.005872170,1064342.572,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.000298291,20952662.604,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.000418162,14946361.454,0
opencl-kernel,Quadro RTX 6000,GPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,50000,1.306237414,38277.881,0
opencl-e2e,Quadro RTX 6000,GPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,50000,1.315309852,38013.857,0
opencl-kernel,Quadro RTX 6000,GPU,w16,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.330040848,151496.399,0
opencl-e2e,Quadro RTX 6000,GPU,w16,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.339932125,147088.187,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.005984377,1044386.054,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.000291601,21433409.985,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.000415802,15031192.661,0
library,Quadro RTX 6000,gpu,cgbn,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.000161216,310142913.855,0
opencl-kernel,Quadro RTX 6000,GPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.122624468,407748.966,0
opencl-e2e,Quadro RTX 6000,GPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.133494264,374547.930,0
opencl-kernel,Quadro RTX 6000,GPU,w16,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.033333981,1499970.837,0
opencl-e2e,Quadro RTX 6000,GPU,w16,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.042469956,1177302.845,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.023189215,269521.845,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.001587085,3938036.289,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.000515727,12118807.132,0
library,Quadro RTX 6000,gpu,cgbn,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000279520,178878076.703,0
opencl-kernel,Quadro RTX 6000,GPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.138301576,361528.781,0
opencl-e2e,Quadro RTX 6000,GPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.148912059,335768.643,0
opencl-kernel,Quadro RTX 6000,GPU,w16,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.002521120,19832455.324,0
opencl-e2e,Quadro RTX 6000,GPU,w16,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.009487102,5270313.422,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,p2048,2048,COMPARE,6250,0.000128394,48678289.884,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,p2048,2048,COMPARE,6250,0.000062708,99668603.339,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,p2048,2048,COMPARE,6250,0.000120837,51722465.086,0
library,Quadro RTX 6000,gpu,cgbn,p2048,2048,COMPARE,50000,0.000121536,411400737.230,0
opencl-kernel,Quadro RTX 6000,GPU,w8,p2048,2048,COMPARE,50000,0.000582818,85790083.669,0
opencl-e2e,Quadro RTX 6000,GPU,w8,p2048,2048,COMPARE,50000,0.007493807,6672176.102,0
opencl-kernel,Quadro RTX 6000,GPU,w16,p2048,2048,COMPARE,50000,0.000347465,143899270.949,0
opencl-e2e,Quadro RTX 6000,GPU,w16,p2048,2048,COMPARE,50000,0.007222198,6923100.190,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,p2048,2048,REDUCE,781,0.000054136,14426650.918,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,p2048,2048,REDUCE,781,0.000062171,12562052.318,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,p2048,2048,REDUCE,781,0.000135023,5784193.437,0
library,Quadro RTX 6000,gpu,cgbn,p2048,2048,REDUCE,50000,0.000118784,420932112.069,0
opencl-kernel,Quadro RTX 6000,GPU,w8,p2048,2048,REDUCE,50000,1.280782944,39038.621,0
opencl-e2e,Quadro RTX 6000,GPU,w8,p2048,2048,REDUCE,50000,1.289765511,38766.737,0
opencl-kernel,Quadro RTX 6000,GPU,w16,p2048,2048,REDUCE,50000,0.011334794,4411196.320,0
opencl-e2e,Quadro RTX 6000,GPU,w16,p2048,2048,REDUCE,50000,0.018293919,2733148.643,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,p2048,2048,MODMUL,390,0.000974791,400085.711,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,p2048,2048,MODMUL,390,0.000096785,4029554.006,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,p2048,2048,MODMUL,390,0.000233533,1670000.524,0
library,Quadro RTX 6000,gpu,cgbn,p2048,2048,MODMUL,50000,0.001370848,36473773.898,0
opencl-kernel,Quadro RTX 6000,GPU,w8,p2048,2048,MODMUL,50000,2.221031133,22512.066,0
opencl-e2e,Quadro RTX 6000,GPU,w8,p2048,2048,MODMUL,50000,2.229554219,22426.008,0
opencl-kernel,Quadro RTX 6000,GPU,w16,p2048,2048,MODMUL,50000,0.027018698,1850570.302,0
opencl-e2e,Quadro RTX 6000,GPU,w16,p2048,2048,MODMUL,50000,0.034395128,1453694.235,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,p2048,2048,MODEXP,97,0.238427982,406.831,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,p2048,2048,MODEXP,97,0.109047588,889.520,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,p2048,2048,MODEXP,97,0.096923640,1000.788,0
library,Quadro RTX 6000,gpu,cgbn,p2048,2048,MODEXP,50000,1.972019672,25354.717,0
opencl-kernel,Quadro RTX 6000,GPU,w8,p2048,2048,MODEXP,50000,301.566856148,165.801,0
opencl-e2e,Quadro RTX 6000,GPU,w8,p2048,2048,MODEXP,50000,301.909952642,165.612,0
opencl-kernel,Quadro RTX 6000,GPU,w16,p2048,2048,MODEXP,50000,0.000000000,inf,0
opencl-e2e,Quadro RTX 6000,GPU,w16,p2048,2048,MODEXP,50000,0.000000000,inf,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,p2048,2048,EXPONENTIATION,97,0.046700343,2077.073,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,p2048,2048,EXPONENTIATION,97,0.003530713,27473.206,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,p2048,2048,EXPONENTIATION,97,0.108162785,896.796,0
opencl-kernel,Quadro RTX 6000,GPU,w8,p2048,2048,EXPONENTIATION,50000,0.000000000,inf,0
opencl-e2e,Quadro RTX 6000,GPU,w8,p2048,2048,EXPONENTIATION,50000,0.000000000,inf,0
opencl-kernel,Quadro RTX 6000,GPU,w16,p2048,2048,EXPONENTIATION,50000,0.000000000,inf,0
opencl-e2e,Quadro RTX 6000,GPU,w16,p2048,2048,EXPONENTIATION,50000,0.000000000,inf,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,p2048,2048,DIVIDE,781,0.000092100,8479880.724,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,p2048,2048,DIVIDE,781,0.000064557,12097757.647,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,p2048,2048,DIVIDE,781,0.000136286,5730595.100,0
library,Quadro RTX 6000,gpu,cgbn,p2048,2048,DIVIDE,50000,0.000116480,429258241.758,0
opencl-kernel,Quadro RTX 6000,GPU,w8,p2048,2048,DIVIDE,50000,3.039860992,16448.121,0
opencl-e2e,Quadro RTX 6000,GPU,w8,p2048,2048,DIVIDE,50000,3.050708240,16389.637,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,p2048,2048,ISQRT,195,0.000197208,988805.930,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,p2048,2048,ISQRT,195,0.000077492,2516400.929,0
opencl-kernel,Quadro RTX 6000,GPU,w8,p2048,2048,ISQRT,50000,38.825405119,1287.817,0
opencl-e2e,Quadro RTX 6000,GPU,w8,p2048,2048,ISQRT,50000,38.823120760,1287.892,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-1t,p2048,2048,MODMUL_R2,6250,0.016127642,387533.398,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,gmp-nt,p2048,2048,MODMUL_R2,6250,0.000696762,8970070.388,0
library,Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz,host-cpu,openssl-nt,p2048,2048,MODMUL_R2,6250,0.001878249,3327568.162,0
library,Quadro RTX 6000,gpu,cgbn,p2048,2048,MODMUL_R2,50000,0.000440032,113628099.775,0
opencl-kernel,Quadro RTX 6000,GPU,w8,p2048,2048,MODMUL_R2,50000,0.124191631,402603.619,0
opencl-e2e,Quadro RTX 6000,GPU,w8,p2048,2048,MODMUL_R2,50000,0.129637642,385690.447,0
```
