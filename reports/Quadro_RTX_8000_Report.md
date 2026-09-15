# MPA-OpenCL benchmark report - Quadro RTX 8000

> **Note.** The multi-threaded GMP and OpenSSL baseline columns have been
> removed from this report: they predate the 2026-09-12 timing fix and were
> understated (see `reports/README.md`). The single-threaded GMP column, the
> OpenCL-on-CPU rows and all MPA measurements are unaffected and were verified
> against GMP before timing.


> **Partial report.** The run was interrupted or hit its time budget.
> Rows that never ran are marked `n/a`.

## 1. System under test

2 OpenCL device(s) exercised with the identical kernels and operands.

### Device 0 - Quadro RTX 8000 (GPU)

| Property | Value |
|---|---|
| Model | Quadro RTX 8000 |
| Type | GPU |
| Vendor | NVIDIA Corporation |
| Device memory | 47.27 GiB |
| Max single allocation | 11.82 GiB |
| Local memory | 48 KiB |
| Global cache | 2304 KiB |
| Compute units | 72 |
| Max clock | 1770 MHz |
| Max work-group size | 1024 |
| OpenCL version | OpenCL 3.0 CUDA |
| Driver | 580.126.20 |

### Device 1 - cpu-haswell-AMD EPYC 7282 16-Core Processor (CPU)

| Property | Value |
|---|---|
| Model | cpu-haswell-AMD EPYC 7282 16-Core Processor |
| Type | CPU |
| Vendor | AuthenticAMD |
| Device memory | 60.66 GiB |
| Max single allocation | 16.00 GiB |
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
| RAM | 62.7 GB |
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
- Total wall time 5400.4 s.

## 3. Correctness

| Device | Kernel | Configs run | Passed | Mismatched | Build/launch failed |
|---|---|---|---|---|---|
| [0] GPU | `mpaKernels_8bits.cl` (w8) | 75 | 75 | 0 | 0 |
| [0] GPU | `mpaKernel_16bits.cl` (w16) | 75 | 75 | 0 | 0 |
| [0] GPU | `mpaKernel_32bits.cl` (w32) | 35 | 35 | 0 | 0 |
| [0] GPU | `mpaKernel_32bits_opt.cl` (w32-opt) | 61 | 61 | 0 | 0 |

**All configurations correct** - 246 configurations, 0 problems.

## 4. Throughput per device

Operations per second, higher is better. Kernel-only timings.

### Device 0 - Quadro RTX 8000 (GPU)

#### secp256k1 (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | 849.15 M | 1.82 G | 2.63 G | 2.74 G | - | - | - | 55.32 M | 1.64 G |
| SUBTRACT | 50000 / 50000 | 850.90 M | 1.93 G | 2.73 G | 2.83 G | - | - | - | 70.65 M | 1.64 G |
| ADDMOD | 50000 / 50000 | 576.54 M | 1.37 G | 2.27 G | 3.07 G | - | - | - | 18.14 M | 1.75 G |
| SUBTRACTMOD | 50000 / 50000 | 577.47 M | 1.37 G | 2.28 G | 3.08 G | - | - | - | 19.62 M | 1.74 G |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | 22.16 M | 114.42 M | 419.44 M | 1.67 G | - | - | - | 43.26 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | 73.25 M | 341.25 M | 1.13 G | 1.12 G | - | - | - | 43.19 M | 1.70 G |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | 285.14 M | 1.35 G | 3.12 G | 2.23 G | - | - | - | 5.41 M | 1.69 G |
| COMPARE | 50000 / 50000 | 884.05 M | 1.97 G | - | 3.26 G | - | - | - | 130.63 M | 1.69 G |
| REDUCE | 50000 / 6250 | 162.93 M | 358.74 M | - | 905.23 M | - | - | - | 44.23 M | 1.70 G |
| MODMUL | 50000 / 3125 | 61.08 M | 148.27 M | - | 235.97 M | - | - | - | 9.52 M | 1.11 G |
| MODEXP | 50000 / 781 | 1.49 M | 9.72 M | - | 10.54 M | - | - | - | 92.37 k | 676.08 k |
| EXPONENTIATION | 50000 / 781 | 973.65 k | 2.75 M | - | 56.63 M | - | - | - | 278.16 k | n/a |
| DIVIDE | 50000 / 6250 | 117.28 M | 176.23 M | - | 317.20 M | - | - | - | 21.84 M | 1.74 G |
| ISQRT | 50000 / 1562 | 9.90 M | 12.12 M | - | 38.02 M | - | - | - | 9.71 M | n/a |
| MODMUL_R2 | 50000 / 50000 | 325.14 M | 1.12 G | - | 981.22 M | - | - | - | 9.47 M | 1.65 G |

#### rsa256(composite) (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | 1.13 G | 1.91 G | 2.85 G | 2.27 G | - | - | - | 54.49 M | 1.74 G |
| SUBTRACT | 50000 / 50000 | 1.14 G | 1.90 G | 2.74 G | 2.26 G | - | - | - | 69.90 M | 1.74 G |
| ADDMOD | 50000 / 50000 | 849.01 M | 1.49 G | 2.50 G | 2.47 G | - | - | - | 20.22 M | 1.66 G |
| SUBTRACTMOD | 50000 / 50000 | 786.03 M | 1.32 G | 2.28 G | 2.41 G | - | - | - | 19.58 M | 1.74 G |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | 31.01 M | 114.52 M | 418.07 M | 1.25 G | - | - | - | 43.12 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | 102.64 M | 335.92 M | 1.13 G | 821.61 M | - | - | - | 43.15 M | 1.74 G |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | 393.98 M | 1.33 G | 2.96 G | 1.89 G | - | - | - | 5.41 M | 1.81 G |
| COMPARE | 50000 / 50000 | 1.19 G | 1.93 G | - | 2.56 G | - | - | - | 132.42 M | 1.74 G |
| REDUCE | 50000 / 6250 | 226.34 M | 352.73 M | - | 703.78 M | - | - | - | 26.93 M | 1.74 G |
| MODMUL | 50000 / 3125 | 85.14 M | 147.91 M | - | 246.01 M | - | - | - | 9.58 M | 1.44 G |
| MODEXP | 50000 / 781 | 1.49 M | 9.74 M | - | 10.58 M | - | - | - | 98.01 k | 694.87 k |
| EXPONENTIATION | 50000 / 781 | 965.82 k | 2.75 M | - | 56.72 M | - | - | - | 279.32 k | n/a |
| DIVIDE | 50000 / 6250 | 117.00 M | 174.07 M | - | 307.52 M | - | - | - | 20.38 M | 1.71 G |
| ISQRT | 50000 / 1562 | 9.35 M | 12.11 M | - | 37.81 M | - | - | - | 9.79 M | n/a |
| MODMUL_R2 | 50000 / 50000 | 311.71 M | 1.07 G | - | 980.83 M | - | - | - | 9.46 M | 1.63 G |

#### brainpoolP512r1 (512-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 25000 | 428.48 M | 777.34 M | 235.13 M | 462.00 M | - | - | - | 50.29 M | 1.51 G |
| SUBTRACT | 50000 / 25000 | 427.86 M | 771.69 M | 241.46 M | 465.53 M | - | - | - | 59.15 M | 1.50 G |
| ADDMOD | 50000 / 25000 | 325.01 M | 629.32 M | 228.97 M | 800.65 M | - | - | - | 18.74 M | 1.50 G |
| SUBTRACTMOD | 50000 / 25000 | 285.09 M | 582.25 M | 223.44 M | 795.68 M | - | - | - | 18.74 M | 1.45 G |
| MULTIPLYOPERANDSCANNING | 50000 / 25000 | 3.79 M | 16.73 M | 14.00 M | 328.65 M | - | - | - | 17.82 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 25000 | 14.64 M | 55.63 M | 71.48 M | 81.32 M | - | - | - | 17.82 M | 1.51 G |
| MONTGOMERYMULTIPLICATION | 50000 / 25000 | 103.58 M | 357.23 M | 852.20 M | 625.77 M | - | - | - | 2.32 M | 1.52 G |
| COMPARE | 50000 / 25000 | 519.63 M | 895.96 M | - | 1.17 G | - | - | - | 125.64 M | 1.51 G |
| REDUCE | 50000 / 3125 | 73.71 M | 96.00 M | - | 265.75 M | - | - | - | 25.35 M | 1.51 G |
| MODMUL | 50000 / 1562 | 26.55 M | 27.37 M | - | 72.30 M | - | - | - | 4.38 M | 519.45 M |
| MODEXP | 50000 / 390 | 111.15 k | 898.79 k | - | 1.05 M | - | - | - | 18.00 k | 82.78 k |
| EXPONENTIATION | 50000 / 390 | 113.33 k | 405.11 k | - | 1.09 M | - | - | - | 80.07 k | n/a |
| DIVIDE | 50000 / 3125 | 38.83 M | 42.23 M | - | 130.44 M | - | - | - | 19.34 M | 1.02 G |
| ISQRT | 50000 / 781 | 1.36 M | 1.51 M | - | 7.84 M | - | - | - | 5.51 M | n/a |
| MODMUL_R2 | 50000 / 25000 | 74.15 M | 333.32 M | - | 443.16 M | - | - | - | 4.71 M | 1.27 G |

#### p1024 (1024-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 12500 | 73.44 M | 147.07 M | 93.25 M | 233.33 M | - | - | - | 42.62 M | 942.97 M |
| SUBTRACT | 50000 / 12500 | 75.66 M | 142.76 M | 93.09 M | 238.67 M | - | - | - | 53.07 M | 940.70 M |
| ADDMOD | 50000 / 12500 | 70.84 M | 130.88 M | 94.27 M | 364.80 M | - | - | - | 14.22 M | 939.00 M |
| SUBTRACTMOD | 50000 / 12500 | 55.33 M | 129.97 M | 93.83 M | 360.77 M | - | - | - | 15.81 M | 943.54 M |
| MULTIPLYOPERANDSCANNING | 50000 / 12500 | 234.79 k | 828.16 k | 2.82 M | 157.70 M | - | - | - | 5.21 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 12500 | 1.81 M | 7.16 M | 19.81 M | 26.72 M | - | - | - | 5.22 M | 946.97 M |
| MONTGOMERYMULTIPLICATION | 50000 / 12500 | 19.43 M | 95.66 M | 349.43 M | 201.37 M | - | - | - | 769.64 k | 828.47 M |
| COMPARE | 50000 / 12500 | 212.61 M | 394.01 M | - | 543.09 M | - | - | - | 131.27 M | 948.69 M |
| REDUCE | 50000 / 1562 | 14.41 M | 26.23 M | - | 76.92 M | - | - | - | 37.23 M | 939.00 M |
| MODMUL | 50000 / 781 | 2.82 M | 9.59 M | - | 16.80 M | - | - | - | 1.66 M | 176.51 M |
| MODEXP | 50000 / 195 | 11.80 k | 78.48 k | - | 139.87 k | - | - | - | 2.71 k | 32.97 k |
| EXPONENTIATION | 50000 / 195 | 11.65 k | 51.29 k | - | 156.89 k | - | - | - | 18.64 k | n/a |
| DIVIDE | 50000 / 1562 | 333.33 k | 3.52 M | - | 24.73 M | - | - | - | 17.59 M | 718.72 M |
| ISQRT | 50000 / 390 | 28.85 k | 185.75 k | - | 1.12 M | - | - | - | 2.81 M | n/a |
| MODMUL_R2 | 50000 / 12500 | 13.62 M | 83.51 M | - | 139.57 M | - | - | - | 1.70 M | 454.74 M |

#### p2048 (2048-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 6250 | 31.79 M | 65.16 M | 118.17 M | over budget | - | - | - | 30.97 M | 498.25 M |
| SUBTRACT | 50000 / 6250 | 31.36 M | 66.19 M | 118.74 M | - | - | - | - | 37.14 M | 497.93 M |
| ADDMOD | 50000 / 6250 | 26.83 M | 52.56 M | 102.75 M | - | - | - | - | 11.09 M | 501.93 M |
| SUBTRACTMOD | 50000 / 6250 | 26.23 M | 53.24 M | 97.03 M | - | - | - | - | 12.95 M | 499.04 M |
| MULTIPLYOPERANDSCANNING | 50000 / 6250 | 48.64 k | 179.17 k | 928.12 k | - | - | - | - | 1.54 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 6250 | 420.88 k | 1.58 M | 6.12 M | - | - | - | - | 1.54 M | 342.50 M |
| MONTGOMERYMULTIPLICATION | 50000 / 6250 | 357.53 k | 23.25 M | 127.88 M | - | - | - | - | 235.94 k | 191.13 M |
| COMPARE | 50000 / 6250 | 104.62 M | 155.18 M | - | - | - | - | - | 130.02 M | 499.84 M |
| REDUCE | 50000 / 781 | 45.31 k | 4.94 M | - | - | - | - | - | 24.07 M | 504.52 M |
| MODMUL | 50000 / 390 | 26.22 k | 1.59 M | - | - | - | - | - | 551.48 k | 38.00 M |
| MODEXP | 50000 / 97 | 198.9 | over budget | - | - | - | - | - | 377.1 | 20.88 k |
| EXPONENTIATION | 50000 / 97 | 906.3 | 4.69 k | - | - | - | - | - | 3.15 k | n/a |
| DIVIDE | 50000 / 781 | 19.04 k | 91.70 k | - | - | - | - | - | 13.81 M | 508.63 M |
| ISQRT | 50000 / 195 | 1.47 k | over budget | - | - | - | - | - | 1.79 M | n/a |
| MODMUL_R2 | 50000 / 6250 | 463.16 k | 22.19 M | - | - | - | - | - | 548.54 k | 125.06 M |

### Device 1 - cpu-haswell-AMD EPYC 7282 16-Core Processor (CPU)

#### secp256k1 (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | - | - | - | - | - | - | - | 55.32 M | 1.64 G |
| SUBTRACT | 50000 / 50000 | - | - | - | - | - | - | - | 70.65 M | 1.64 G |
| ADDMOD | 50000 / 50000 | - | - | - | - | - | - | - | 18.14 M | 1.75 G |
| SUBTRACTMOD | 50000 / 50000 | - | - | - | - | - | - | - | 19.62 M | 1.74 G |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | - | - | - | - | - | - | - | 43.26 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | - | - | - | - | - | - | - | 43.19 M | 1.70 G |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | - | - | - | - | - | - | - | 5.41 M | 1.69 G |
| COMPARE | 50000 / 50000 | - | - | - | - | - | - | - | 130.63 M | 1.69 G |
| REDUCE | 50000 / 6250 | - | - | - | - | - | - | - | 44.23 M | 1.70 G |
| MODMUL | 50000 / 3125 | - | - | - | - | - | - | - | 9.52 M | 1.11 G |
| MODEXP | 50000 / 781 | - | - | - | - | - | - | - | 92.37 k | 676.08 k |
| EXPONENTIATION | 50000 / 781 | - | - | - | - | - | - | - | 278.16 k | n/a |
| DIVIDE | 50000 / 6250 | - | - | - | - | - | - | - | 21.84 M | 1.74 G |
| ISQRT | 50000 / 1562 | - | - | - | - | - | - | - | 9.71 M | n/a |
| MODMUL_R2 | 50000 / 50000 | - | - | - | - | - | - | - | 9.47 M | 1.65 G |

#### rsa256(composite) (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | - | - | - | - | - | - | - | 54.49 M | 1.74 G |
| SUBTRACT | 50000 / 50000 | - | - | - | - | - | - | - | 69.90 M | 1.74 G |
| ADDMOD | 50000 / 50000 | - | - | - | - | - | - | - | 20.22 M | 1.66 G |
| SUBTRACTMOD | 50000 / 50000 | - | - | - | - | - | - | - | 19.58 M | 1.74 G |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | - | - | - | - | - | - | - | 43.12 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | - | - | - | - | - | - | - | 43.15 M | 1.74 G |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | - | - | - | - | - | - | - | 5.41 M | 1.81 G |
| COMPARE | 50000 / 50000 | - | - | - | - | - | - | - | 132.42 M | 1.74 G |
| REDUCE | 50000 / 6250 | - | - | - | - | - | - | - | 26.93 M | 1.74 G |
| MODMUL | 50000 / 3125 | - | - | - | - | - | - | - | 9.58 M | 1.44 G |
| MODEXP | 50000 / 781 | - | - | - | - | - | - | - | 98.01 k | 694.87 k |
| EXPONENTIATION | 50000 / 781 | - | - | - | - | - | - | - | 279.32 k | n/a |
| DIVIDE | 50000 / 6250 | - | - | - | - | - | - | - | 20.38 M | 1.71 G |
| ISQRT | 50000 / 1562 | - | - | - | - | - | - | - | 9.79 M | n/a |
| MODMUL_R2 | 50000 / 50000 | - | - | - | - | - | - | - | 9.46 M | 1.63 G |

#### brainpoolP512r1 (512-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 25000 | - | - | - | - | - | - | - | 50.29 M | 1.51 G |
| SUBTRACT | 50000 / 25000 | - | - | - | - | - | - | - | 59.15 M | 1.50 G |
| ADDMOD | 50000 / 25000 | - | - | - | - | - | - | - | 18.74 M | 1.50 G |
| SUBTRACTMOD | 50000 / 25000 | - | - | - | - | - | - | - | 18.74 M | 1.45 G |
| MULTIPLYOPERANDSCANNING | 50000 / 25000 | - | - | - | - | - | - | - | 17.82 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 25000 | - | - | - | - | - | - | - | 17.82 M | 1.51 G |
| MONTGOMERYMULTIPLICATION | 50000 / 25000 | - | - | - | - | - | - | - | 2.32 M | 1.52 G |
| COMPARE | 50000 / 25000 | - | - | - | - | - | - | - | 125.64 M | 1.51 G |
| REDUCE | 50000 / 3125 | - | - | - | - | - | - | - | 25.35 M | 1.51 G |
| MODMUL | 50000 / 1562 | - | - | - | - | - | - | - | 4.38 M | 519.45 M |
| MODEXP | 50000 / 390 | - | - | - | - | - | - | - | 18.00 k | 82.78 k |
| EXPONENTIATION | 50000 / 390 | - | - | - | - | - | - | - | 80.07 k | n/a |
| DIVIDE | 50000 / 3125 | - | - | - | - | - | - | - | 19.34 M | 1.02 G |
| ISQRT | 50000 / 781 | - | - | - | - | - | - | - | 5.51 M | n/a |
| MODMUL_R2 | 50000 / 25000 | - | - | - | - | - | - | - | 4.71 M | 1.27 G |

#### p1024 (1024-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 12500 | - | - | - | - | - | - | - | 42.62 M | 942.97 M |
| SUBTRACT | 50000 / 12500 | - | - | - | - | - | - | - | 53.07 M | 940.70 M |
| ADDMOD | 50000 / 12500 | - | - | - | - | - | - | - | 14.22 M | 939.00 M |
| SUBTRACTMOD | 50000 / 12500 | - | - | - | - | - | - | - | 15.81 M | 943.54 M |
| MULTIPLYOPERANDSCANNING | 50000 / 12500 | - | - | - | - | - | - | - | 5.21 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 12500 | - | - | - | - | - | - | - | 5.22 M | 946.97 M |
| MONTGOMERYMULTIPLICATION | 50000 / 12500 | - | - | - | - | - | - | - | 769.64 k | 828.47 M |
| COMPARE | 50000 / 12500 | - | - | - | - | - | - | - | 131.27 M | 948.69 M |
| REDUCE | 50000 / 1562 | - | - | - | - | - | - | - | 37.23 M | 939.00 M |
| MODMUL | 50000 / 781 | - | - | - | - | - | - | - | 1.66 M | 176.51 M |
| MODEXP | 50000 / 195 | - | - | - | - | - | - | - | 2.71 k | 32.97 k |
| EXPONENTIATION | 50000 / 195 | - | - | - | - | - | - | - | 18.64 k | n/a |
| DIVIDE | 50000 / 1562 | - | - | - | - | - | - | - | 17.59 M | 718.72 M |
| ISQRT | 50000 / 390 | - | - | - | - | - | - | - | 2.81 M | n/a |
| MODMUL_R2 | 50000 / 12500 | - | - | - | - | - | - | - | 1.70 M | 454.74 M |

#### p2048 (2048-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 6250 | - | - | - | - | - | - | - | 30.97 M | 498.25 M |
| SUBTRACT | 50000 / 6250 | - | - | - | - | - | - | - | 37.14 M | 497.93 M |
| ADDMOD | 50000 / 6250 | - | - | - | - | - | - | - | 11.09 M | 501.93 M |
| SUBTRACTMOD | 50000 / 6250 | - | - | - | - | - | - | - | 12.95 M | 499.04 M |
| MULTIPLYOPERANDSCANNING | 50000 / 6250 | - | - | - | - | - | - | - | 1.54 M | n/a |
| MULTIPLYPRODUCTSCANNING | 50000 / 6250 | - | - | - | - | - | - | - | 1.54 M | 342.50 M |
| MONTGOMERYMULTIPLICATION | 50000 / 6250 | - | - | - | - | - | - | - | 235.94 k | 191.13 M |
| COMPARE | 50000 / 6250 | - | - | - | - | - | - | - | 130.02 M | 499.84 M |
| REDUCE | 50000 / 781 | - | - | - | - | - | - | - | 24.07 M | 504.52 M |
| MODMUL | 50000 / 390 | - | - | - | - | - | - | - | 551.48 k | 38.00 M |
| MODEXP | 50000 / 97 | - | - | - | - | - | - | - | 377.1 | 20.88 k |
| EXPONENTIATION | 50000 / 97 | - | - | - | - | - | - | - | 3.15 k | n/a |
| DIVIDE | 50000 / 781 | - | - | - | - | - | - | - | 13.81 M | 508.63 M |
| ISQRT | 50000 / 195 | - | - | - | - | - | - | - | 1.79 M | n/a |
| MODMUL_R2 | 50000 / 6250 | - | - | - | - | - | - | - | 548.54 k | 125.06 M |

## 5. Head to head

Best OpenCL GPU result against best OpenCL CPU result and the CPU libraries.
Ratios above 1.00x mean the GPU is faster than that baseline.

### secp256k1 (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w32-opt | 2.74 G | none | n/a | 55.32 M | 1.64 G | n/a | 1.67x |
| SUBTRACT | w32-opt | 2.83 G | none | n/a | 70.65 M | 1.64 G | n/a | 1.73x |
| ADDMOD | w32-opt | 3.07 G | none | n/a | 18.14 M | 1.75 G | n/a | 1.76x |
| SUBTRACTMOD | w32-opt | 3.08 G | none | n/a | 19.62 M | 1.74 G | n/a | 1.77x |
| MULTIPLYOPERANDSCANNING | w32-opt | 1.67 G | none | n/a | 43.26 M | n/a | n/a | n/a |
| MULTIPLYPRODUCTSCANNING | w32 | 1.13 G | none | n/a | 43.19 M | 1.70 G | n/a | 0.67x |
| MONTGOMERYMULTIPLICATION | w32 | 3.12 G | none | n/a | 5.41 M | 1.69 G | n/a | 1.84x |
| COMPARE | w32-opt | 3.26 G | none | n/a | 130.63 M | 1.69 G | n/a | 1.93x |
| REDUCE | w32-opt | 113.15 M | none | n/a | 44.23 M | 1.70 G | n/a | 0.07x |
| MODMUL | w32-opt | 14.75 M | none | n/a | 9.52 M | 1.11 G | n/a | 0.01x |
| MODEXP | w32-opt | 164.59 k | none | n/a | 92.37 k | 676.08 k | n/a | 0.24x |
| EXPONENTIATION | w32-opt | 884.53 k | none | n/a | 278.16 k | n/a | n/a | n/a |
| DIVIDE | w32-opt | 39.65 M | none | n/a | 21.84 M | 1.74 G | n/a | 0.02x |
| ISQRT | w32-opt | 1.19 M | none | n/a | 9.71 M | n/a | n/a | n/a |
| MODMUL_R2 | w16 | 1.12 G | none | n/a | 9.47 M | 1.65 G | n/a | 0.68x |

### rsa256(composite) (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w32 | 2.85 G | none | n/a | 54.49 M | 1.74 G | n/a | 1.64x |
| SUBTRACT | w32 | 2.74 G | none | n/a | 69.90 M | 1.74 G | n/a | 1.57x |
| ADDMOD | w32 | 2.50 G | none | n/a | 20.22 M | 1.66 G | n/a | 1.50x |
| SUBTRACTMOD | w32-opt | 2.41 G | none | n/a | 19.58 M | 1.74 G | n/a | 1.38x |
| MULTIPLYOPERANDSCANNING | w32-opt | 1.25 G | none | n/a | 43.12 M | n/a | n/a | n/a |
| MULTIPLYPRODUCTSCANNING | w32 | 1.13 G | none | n/a | 43.15 M | 1.74 G | n/a | 0.65x |
| MONTGOMERYMULTIPLICATION | w32 | 2.96 G | none | n/a | 5.41 M | 1.81 G | n/a | 1.63x |
| COMPARE | w32-opt | 2.56 G | none | n/a | 132.42 M | 1.74 G | n/a | 1.47x |
| REDUCE | w32-opt | 87.97 M | none | n/a | 26.93 M | 1.74 G | n/a | 0.05x |
| MODMUL | w32-opt | 15.38 M | none | n/a | 9.58 M | 1.44 G | n/a | 0.01x |
| MODEXP | w32-opt | 165.23 k | none | n/a | 98.01 k | 694.87 k | n/a | 0.24x |
| EXPONENTIATION | w32-opt | 885.96 k | none | n/a | 279.32 k | n/a | n/a | n/a |
| DIVIDE | w32-opt | 38.44 M | none | n/a | 20.38 M | 1.71 G | n/a | 0.02x |
| ISQRT | w32-opt | 1.18 M | none | n/a | 9.79 M | n/a | n/a | n/a |
| MODMUL_R2 | w16 | 1.07 G | none | n/a | 9.46 M | 1.63 G | n/a | 0.66x |

### brainpoolP512r1 (512-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w16 | 388.67 M | none | n/a | 50.29 M | 1.51 G | n/a | 0.26x |
| SUBTRACT | w16 | 385.84 M | none | n/a | 59.15 M | 1.50 G | n/a | 0.26x |
| ADDMOD | w32-opt | 400.33 M | none | n/a | 18.74 M | 1.50 G | n/a | 0.27x |
| SUBTRACTMOD | w32-opt | 397.84 M | none | n/a | 18.74 M | 1.45 G | n/a | 0.27x |
| MULTIPLYOPERANDSCANNING | w32-opt | 164.32 M | none | n/a | 17.82 M | n/a | n/a | n/a |
| MULTIPLYPRODUCTSCANNING | w32-opt | 40.66 M | none | n/a | 17.82 M | 1.51 G | n/a | 0.03x |
| MONTGOMERYMULTIPLICATION | w32 | 426.10 M | none | n/a | 2.32 M | 1.52 G | n/a | 0.28x |
| COMPARE | w32-opt | 582.87 M | none | n/a | 125.64 M | 1.51 G | n/a | 0.38x |
| REDUCE | w32-opt | 16.61 M | none | n/a | 25.35 M | 1.51 G | n/a | 0.01x |
| MODMUL | w32-opt | 2.26 M | none | n/a | 4.38 M | 519.45 M | n/a | 0.00x |
| MODEXP | w32-opt | 8.18 k | none | n/a | 18.00 k | 82.78 k | n/a | 0.10x |
| EXPONENTIATION | w32-opt | 8.52 k | none | n/a | 80.07 k | n/a | n/a | n/a |
| DIVIDE | w32-opt | 8.15 M | none | n/a | 19.34 M | 1.02 G | n/a | 0.01x |
| ISQRT | w32-opt | 122.52 k | none | n/a | 5.51 M | n/a | n/a | n/a |
| MODMUL_R2 | w32-opt | 221.58 M | none | n/a | 4.71 M | 1.27 G | n/a | 0.18x |

### p1024 (1024-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w32-opt | 58.33 M | none | n/a | 42.62 M | 942.97 M | n/a | 0.06x |
| SUBTRACT | w32-opt | 59.67 M | none | n/a | 53.07 M | 940.70 M | n/a | 0.06x |
| ADDMOD | w32-opt | 91.20 M | none | n/a | 14.22 M | 939.00 M | n/a | 0.10x |
| SUBTRACTMOD | w32-opt | 90.19 M | none | n/a | 15.81 M | 943.54 M | n/a | 0.10x |
| MULTIPLYOPERANDSCANNING | w32-opt | 39.42 M | none | n/a | 5.21 M | n/a | n/a | n/a |
| MULTIPLYPRODUCTSCANNING | w32-opt | 6.68 M | none | n/a | 5.22 M | 946.97 M | n/a | 0.01x |
| MONTGOMERYMULTIPLICATION | w32 | 87.36 M | none | n/a | 769.64 k | 828.47 M | n/a | 0.11x |
| COMPARE | w32-opt | 135.77 M | none | n/a | 131.27 M | 948.69 M | n/a | 0.14x |
| REDUCE | w32-opt | 2.40 M | none | n/a | 37.23 M | 939.00 M | n/a | 0.00x |
| MODMUL | w32-opt | 262.35 k | none | n/a | 1.66 M | 176.51 M | n/a | 0.00x |
| MODEXP | w32-opt | 545.5 | none | n/a | 2.71 k | 32.97 k | n/a | 0.02x |
| EXPONENTIATION | w32-opt | 611.9 | none | n/a | 18.64 k | n/a | n/a | n/a |
| DIVIDE | w32-opt | 772.47 k | none | n/a | 17.59 M | 718.72 M | n/a | 0.00x |
| ISQRT | w32-opt | 8.71 k | none | n/a | 2.81 M | n/a | n/a | n/a |
| MODMUL_R2 | w32-opt | 34.89 M | none | n/a | 1.70 M | 454.74 M | n/a | 0.08x |

### p2048 (2048-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w32 | 14.77 M | none | n/a | 30.97 M | 498.25 M | n/a | 0.03x |
| SUBTRACT | w32 | 14.84 M | none | n/a | 37.14 M | 497.93 M | n/a | 0.03x |
| ADDMOD | w32 | 12.84 M | none | n/a | 11.09 M | 501.93 M | n/a | 0.03x |
| SUBTRACTMOD | w32 | 12.13 M | none | n/a | 12.95 M | 499.04 M | n/a | 0.02x |
| MULTIPLYOPERANDSCANNING | w32 | 116.01 k | none | n/a | 1.54 M | n/a | n/a | n/a |
| MULTIPLYPRODUCTSCANNING | w32 | 765.22 k | none | n/a | 1.54 M | 342.50 M | n/a | 0.00x |
| MONTGOMERYMULTIPLICATION | w32 | 15.98 M | none | n/a | 235.94 k | 191.13 M | n/a | 0.08x |
| COMPARE | w16 | 19.40 M | none | n/a | 130.02 M | 499.84 M | n/a | 0.04x |
| REDUCE | w16 | 77.18 k | none | n/a | 24.07 M | 504.52 M | n/a | 0.00x |
| MODMUL | w16 | 12.41 k | none | n/a | 551.48 k | 38.00 M | n/a | 0.00x |
| MODEXP | w8 | 0.4 | none | n/a | 377.1 | 20.88 k | n/a | 0.00x |
| EXPONENTIATION | w16 | 9.1 | none | n/a | 3.15 k | n/a | n/a | n/a |
| DIVIDE | w16 | 1.43 k | none | n/a | 13.81 M | 508.63 M | n/a | 0.00x |
| ISQRT | w8 | 5.7 | none | n/a | 1.79 M | n/a | n/a | n/a |
| MODMUL_R2 | w16 | 2.77 M | none | n/a | 548.54 k | 125.06 M | n/a | 0.02x |

## 6. Raw data

Also written to `Quadro_RTX_8000_Report.csv` for analysis.

```csv
kind,device,device_type,kernel,modulus,bits,operation,items,seconds,ops_per_sec,mismatches
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,ADD,50000,0.000903907,55315425.200,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,ADD,50000,0.000102164,489409228.539,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,ADD,50000,0.000105170,475421495.193,0
library,Quadro RTX 8000,gpu,cgbn,secp256k1,256,ADD,50000,0.000030464,1641281512.605,0
opencl-kernel,Quadro RTX 8000,GPU,w8,secp256k1,256,ADD,50000,0.000058882,849153468.804,0
opencl-e2e,Quadro RTX 8000,GPU,w8,secp256k1,256,ADD,50000,0.000752169,66474427.905,0
opencl-kernel,Quadro RTX 8000,GPU,w16,secp256k1,256,ADD,50000,0.000027442,1822016788.984,0
opencl-e2e,Quadro RTX 8000,GPU,w16,secp256k1,256,ADD,50000,0.000675503,74018923.758,0
opencl-kernel,Quadro RTX 8000,GPU,w32,secp256k1,256,ADD,50000,0.000019016,2629367903.714,0
opencl-e2e,Quadro RTX 8000,GPU,w32,secp256k1,256,ADD,50000,0.000660906,75653714.652,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,secp256k1,256,ADD,50000,0.000018215,2744984403.001,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,secp256k1,256,ADD,50000,0.000699129,71517564.248,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,SUBTRACT,50000,0.000707674,70653994.408,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,SUBTRACT,50000,0.000131239,380984455.006,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,SUBTRACT,50000,0.000119777,417442655.480,0
library,Quadro RTX 8000,gpu,cgbn,secp256k1,256,SUBTRACT,50000,0.000030528,1637840670.860,0
opencl-kernel,Quadro RTX 8000,GPU,w8,secp256k1,256,SUBTRACT,50000,0.000058761,850903072.784,0
opencl-e2e,Quadro RTX 8000,GPU,w8,secp256k1,256,SUBTRACT,50000,0.000723675,69091798.620,0
opencl-kernel,Quadro RTX 8000,GPU,w16,secp256k1,256,SUBTRACT,50000,0.000025919,1929073900.936,0
opencl-e2e,Quadro RTX 8000,GPU,w16,secp256k1,256,SUBTRACT,50000,0.000676525,73907117.852,0
opencl-kernel,Quadro RTX 8000,GPU,w32,secp256k1,256,SUBTRACT,50000,0.000018295,2732967214.325,0
opencl-e2e,Quadro RTX 8000,GPU,w32,secp256k1,256,SUBTRACT,50000,0.000666507,75017960.712,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,secp256k1,256,SUBTRACT,50000,0.000017663,2830776473.070,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,secp256k1,256,SUBTRACT,50000,0.000684772,73017012.771,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,ADDMOD,50000,0.002756126,18141405.778,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,ADDMOD,50000,0.000319286,156599404.225,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,ADDMOD,50000,0.000848301,58941343.598,0
library,Quadro RTX 8000,gpu,cgbn,secp256k1,256,ADDMOD,50000,0.000028608,1747762863.535,0
opencl-kernel,Quadro RTX 8000,GPU,w8,secp256k1,256,ADDMOD,50000,0.000086724,576541275.837,0
opencl-e2e,Quadro RTX 8000,GPU,w8,secp256k1,256,ADDMOD,50000,0.000748091,66836796.465,0
opencl-kernel,Quadro RTX 8000,GPU,w16,secp256k1,256,ADDMOD,50000,0.000036629,1365041728.960,0
opencl-e2e,Quadro RTX 8000,GPU,w16,secp256k1,256,ADDMOD,50000,0.000683689,73132664.447,0
opencl-kernel,Quadro RTX 8000,GPU,w32,secp256k1,256,ADDMOD,50000,0.000021982,2274587603.271,0
opencl-e2e,Quadro RTX 8000,GPU,w32,secp256k1,256,ADDMOD,50000,0.000660375,75714556.873,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,secp256k1,256,ADDMOD,50000,0.000016291,3069192996.899,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,secp256k1,256,ADDMOD,50000,0.000676626,73896055.005,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,SUBTRACTMOD,50000,0.002548301,19620916.870,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,SUBTRACTMOD,50000,0.000293137,170568686.925,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,SUBTRACTMOD,50000,0.000902915,55376203.858,0
library,Quadro RTX 8000,gpu,cgbn,secp256k1,256,SUBTRACTMOD,50000,0.000028672,1743861607.143,0
opencl-kernel,Quadro RTX 8000,GPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.000086584,577473041.551,0
opencl-e2e,Quadro RTX 8000,GPU,w8,secp256k1,256,SUBTRACTMOD,50000,0.000749093,66747385.192,0
opencl-kernel,Quadro RTX 8000,GPU,w16,secp256k1,256,SUBTRACTMOD,50000,0.000036399,1373659846.609,0
opencl-e2e,Quadro RTX 8000,GPU,w16,secp256k1,256,SUBTRACTMOD,50000,0.000685061,72986191.077,0
opencl-kernel,Quadro RTX 8000,GPU,w32,secp256k1,256,SUBTRACTMOD,50000,0.000021891,2284046806.564,0
opencl-e2e,Quadro RTX 8000,GPU,w32,secp256k1,256,SUBTRACTMOD,50000,0.000668050,74844694.699,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,secp256k1,256,SUBTRACTMOD,50000,0.000016241,3078653049.287,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,secp256k1,256,SUBTRACTMOD,50000,0.000676666,73891681.643,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001155915,43255768.077,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000146849,340485617.524,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000281325,177730482.231,0
opencl-kernel,Quadro RTX 8000,GPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.002256537,22157847.704,0
opencl-e2e,Quadro RTX 8000,GPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.003129074,15979167.680,0
opencl-kernel,Quadro RTX 8000,GPU,w16,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000436980,114421701.020,0
opencl-e2e,Quadro RTX 8000,GPU,w16,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.001211572,41268699.773,0
opencl-kernel,Quadro RTX 8000,GPU,w32,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000119207,419438592.160,0
opencl-e2e,Quadro RTX 8000,GPU,w32,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000891955,56056633.110,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000029927,1670738435.446,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,secp256k1,256,MULTIPLYOPERANDSCANNING,50000,0.000806122,62025362.227,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001157579,43193595.834,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000150436,332367098.270,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000280522,178239112.641,0
library,Quadro RTX 8000,gpu,cgbn,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000029472,1696525515.744,0
opencl-kernel,Quadro RTX 8000,GPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000682567,73252880.694,0
opencl-e2e,Quadro RTX 8000,GPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.001504418,33235444.436,0
opencl-kernel,Quadro RTX 8000,GPU,w16,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000146519,341252297.889,0
opencl-e2e,Quadro RTX 8000,GPU,w16,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000927251,53922841.865,0
opencl-kernel,Quadro RTX 8000,GPU,w32,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000044124,1133170975.827,0
opencl-e2e,Quadro RTX 8000,GPU,w32,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000820520,60936976.584,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000044755,1117195128.524,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,secp256k1,256,MULTIPLYPRODUCTSCANNING,50000,0.000828194,60372330.679,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.009234942,5414219.151,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.005993595,8342238.514,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000239775,208528785.251,0
library,Quadro RTX 8000,gpu,cgbn,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000029536,1692849404.117,0
opencl-kernel,Quadro RTX 8000,GPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000175353,285138892.312,0
opencl-e2e,Quadro RTX 8000,GPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000841999,59382496.091,0
opencl-kernel,Quadro RTX 8000,GPU,w16,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000036990,1351715321.235,0
opencl-e2e,Quadro RTX 8000,GPU,w16,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000656658,76143154.049,0
opencl-kernel,Quadro RTX 8000,GPU,w32,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000016051,3115049025.950,0
opencl-e2e,Quadro RTX 8000,GPU,w32,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000669813,74647698.348,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000022463,2225879110.263,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,secp256k1,256,MONTGOMERYMULTIPLICATION,50000,0.000670324,74590786.031,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,COMPARE,50000,0.000382747,130634584.637,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,COMPARE,50000,0.000061717,810148015.075,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,COMPARE,50000,0.000101823,491048702.341,0
library,Quadro RTX 8000,gpu,cgbn,secp256k1,256,COMPARE,50000,0.000029600,1689189189.189,0
opencl-kernel,Quadro RTX 8000,GPU,w8,secp256k1,256,COMPARE,50000,0.000056558,884047361.423,0
opencl-e2e,Quadro RTX 8000,GPU,w8,secp256k1,256,COMPARE,50000,0.000738192,67733045.178,0
opencl-kernel,Quadro RTX 8000,GPU,w16,secp256k1,256,COMPARE,50000,0.000025368,1970982192.648,0
opencl-e2e,Quadro RTX 8000,GPU,w16,secp256k1,256,COMPARE,50000,0.000693057,72144141.834,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,secp256k1,256,COMPARE,50000,0.000015349,3257514179.965,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,secp256k1,256,COMPARE,50000,0.000673320,74258905.967,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,REDUCE,6250,0.000141309,44229279.079,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,REDUCE,6250,0.000027612,226350169.066,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,REDUCE,6250,0.000112003,55801986.488,0
library,Quadro RTX 8000,gpu,cgbn,secp256k1,256,REDUCE,50000,0.000029344,1703925845.147,0
opencl-kernel,Quadro RTX 8000,GPU,w8,secp256k1,256,REDUCE,50000,0.000306883,162928474.652,0
opencl-e2e,Quadro RTX 8000,GPU,w8,secp256k1,256,REDUCE,50000,0.000978749,51085619.865,0
opencl-kernel,Quadro RTX 8000,GPU,w16,secp256k1,256,REDUCE,50000,0.000139375,358744434.700,0
opencl-e2e,Quadro RTX 8000,GPU,w16,secp256k1,256,REDUCE,50000,0.000782206,63921783.787,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,secp256k1,256,REDUCE,50000,0.000055235,905225116.342,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,secp256k1,256,REDUCE,50000,0.000700581,71369346.342,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,MODMUL,3125,0.000328273,9519514.413,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,MODMUL,3125,0.000049073,63680617.934,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,MODMUL,3125,0.000126691,24666300.580,0
library,Quadro RTX 8000,gpu,cgbn,secp256k1,256,MODMUL,50000,0.000045056,1109730113.636,0
opencl-kernel,Quadro RTX 8000,GPU,w8,secp256k1,256,MODMUL,50000,0.000818545,61083996.696,0
opencl-e2e,Quadro RTX 8000,GPU,w8,secp256k1,256,MODMUL,50000,0.001490131,33554096.459,0
opencl-kernel,Quadro RTX 8000,GPU,w16,secp256k1,256,MODMUL,50000,0.000337220,148271145.145,0
opencl-e2e,Quadro RTX 8000,GPU,w16,secp256k1,256,MODMUL,50000,0.000977938,51127982.302,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,secp256k1,256,MODMUL,50000,0.000211893,235968284.742,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,secp256k1,256,MODMUL,50000,0.000861988,58005452.095,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,MODEXP,781,0.008455371,92367.325,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,MODEXP,781,0.000723975,1078766.514,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,MODEXP,781,0.001398877,558304.980,0
library,Quadro RTX 8000,gpu,cgbn,secp256k1,256,MODEXP,50000,0.073955581,676081.498,0
opencl-kernel,Quadro RTX 8000,GPU,w8,secp256k1,256,MODEXP,50000,0.033577923,1489073.641,0
opencl-e2e,Quadro RTX 8000,GPU,w8,secp256k1,256,MODEXP,50000,0.034199123,1462025.794,0
opencl-kernel,Quadro RTX 8000,GPU,w16,secp256k1,256,MODEXP,50000,0.005142410,9723067.489,0
opencl-e2e,Quadro RTX 8000,GPU,w16,secp256k1,256,MODEXP,50000,0.005908196,8462820.323,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,secp256k1,256,MODEXP,50000,0.004745126,10537127.984,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,secp256k1,256,MODEXP,50000,0.005506794,9079693.246,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,EXPONENTIATION,781,0.002807774,278156.279,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,EXPONENTIATION,781,0.000263831,2960227.065,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,EXPONENTIATION,781,0.004177235,186965.777,0
opencl-kernel,Quadro RTX 8000,GPU,w8,secp256k1,256,EXPONENTIATION,50000,0.051353359,973646.147,0
opencl-e2e,Quadro RTX 8000,GPU,w8,secp256k1,256,EXPONENTIATION,50000,0.051900598,963380.037,0
opencl-kernel,Quadro RTX 8000,GPU,w16,secp256k1,256,EXPONENTIATION,50000,0.018182035,2749967.231,0
opencl-e2e,Quadro RTX 8000,GPU,w16,secp256k1,256,EXPONENTIATION,50000,0.018724036,2670364.454,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,secp256k1,256,EXPONENTIATION,50000,0.000882958,56627831.097,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,secp256k1,256,EXPONENTIATION,50000,0.001536019,32551682.507,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,DIVIDE,6250,0.000286193,21838402.200,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,DIVIDE,6250,0.000032471,192479281.812,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,DIVIDE,6250,0.000089109,70138862.876,0
library,Quadro RTX 8000,gpu,cgbn,secp256k1,256,DIVIDE,50000,0.000028672,1743861607.143,0
opencl-kernel,Quadro RTX 8000,GPU,w8,secp256k1,256,DIVIDE,50000,0.000426319,117283077.794,0
opencl-e2e,Quadro RTX 8000,GPU,w8,secp256k1,256,DIVIDE,50000,0.001221290,40940314.153,0
opencl-kernel,Quadro RTX 8000,GPU,w16,secp256k1,256,DIVIDE,50000,0.000283720,176230103.531,0
opencl-e2e,Quadro RTX 8000,GPU,w16,secp256k1,256,DIVIDE,50000,0.001059152,47207569.096,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,secp256k1,256,DIVIDE,50000,0.000157630,317198828.386,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,secp256k1,256,DIVIDE,50000,0.000940227,53178650.126,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,ISQRT,1562,0.000160836,9711761.136,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,ISQRT,1562,0.000015730,99300457.613,0
opencl-kernel,Quadro RTX 8000,GPU,w8,secp256k1,256,ISQRT,50000,0.005048189,9904541.849,0
opencl-e2e,Quadro RTX 8000,GPU,w8,secp256k1,256,ISQRT,50000,0.005843140,8557043.030,0
opencl-kernel,Quadro RTX 8000,GPU,w16,secp256k1,256,ISQRT,50000,0.004126010,12118244.996,0
opencl-e2e,Quadro RTX 8000,GPU,w16,secp256k1,256,ISQRT,50000,0.004859854,10288374.915,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,secp256k1,256,ISQRT,50000,0.001314969,38023709.056,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,secp256k1,256,ISQRT,50000,0.001968691,25397584.899,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,MODMUL_R2,50000,0.005281893,9466302.958,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,MODMUL_R2,50000,0.000286034,174804490.007,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,MODMUL_R2,50000,0.000716421,69791360.365,0
library,Quadro RTX 8000,gpu,cgbn,secp256k1,256,MODMUL_R2,50000,0.000030336,1648206751.055,0
opencl-kernel,Quadro RTX 8000,GPU,w8,secp256k1,256,MODMUL_R2,50000,0.000153782,325135414.519,0
opencl-e2e,Quadro RTX 8000,GPU,w8,secp256k1,256,MODMUL_R2,50000,0.000811732,61596673.898,0
opencl-kernel,Quadro RTX 8000,GPU,w16,secp256k1,256,MODMUL_R2,50000,0.000044584,1121477514.400,0
opencl-e2e,Quadro RTX 8000,GPU,w16,secp256k1,256,MODMUL_R2,50000,0.000683849,73115558.462,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,secp256k1,256,MODMUL_R2,50000,0.000050957,981222367.014,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,secp256k1,256,MODMUL_R2,50000,0.000696533,71784118.792,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,ADD,50000,0.000917623,54488604.957,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,ADD,50000,0.000128313,389672227.908,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,ADD,50000,0.000105280,474924176.976,0
library,Quadro RTX 8000,gpu,cgbn,rsa256(composite),256,ADD,50000,0.000028672,1743861607.143,0
opencl-kernel,Quadro RTX 8000,GPU,w8,rsa256(composite),256,ADD,50000,0.000044194,1131374016.395,0
opencl-e2e,Quadro RTX 8000,GPU,w8,rsa256(composite),256,ADD,50000,0.000745156,67100036.151,0
opencl-kernel,Quadro RTX 8000,GPU,w16,rsa256(composite),256,ADD,50000,0.000026170,1910589638.698,0
opencl-e2e,Quadro RTX 8000,GPU,w16,rsa256(composite),256,ADD,50000,0.000676806,73876404.368,0
opencl-kernel,Quadro RTX 8000,GPU,w32,rsa256(composite),256,ADD,50000,0.000017514,2854861141.687,0
opencl-e2e,Quadro RTX 8000,GPU,w32,rsa256(composite),256,ADD,50000,0.000664974,75190898.733,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,rsa256(composite),256,ADD,50000,0.000022012,2271483957.225,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,rsa256(composite),256,ADD,50000,0.000667678,74886375.772,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACT,50000,0.000715318,69898991.461,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACT,50000,0.000152670,327504094.481,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACT,50000,0.000151979,328992738.052,0
library,Quadro RTX 8000,gpu,cgbn,rsa256(composite),256,SUBTRACT,50000,0.000028672,1743861607.143,0
opencl-kernel,Quadro RTX 8000,GPU,w8,rsa256(composite),256,SUBTRACT,50000,0.000044044,1135225644.929,0
opencl-e2e,Quadro RTX 8000,GPU,w8,rsa256(composite),256,SUBTRACT,50000,0.000746789,66953323.061,0
opencl-kernel,Quadro RTX 8000,GPU,w16,rsa256(composite),256,SUBTRACT,50000,0.000026320,1899688305.439,0
opencl-e2e,Quadro RTX 8000,GPU,w16,rsa256(composite),256,SUBTRACT,50000,0.000678359,73707302.999,0
opencl-kernel,Quadro RTX 8000,GPU,w32,rsa256(composite),256,SUBTRACT,50000,0.000018265,2737461309.402,0
opencl-e2e,Quadro RTX 8000,GPU,w32,rsa256(composite),256,SUBTRACT,50000,0.000662609,75459282.529,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,rsa256(composite),256,SUBTRACT,50000,0.000022102,2262247461.734,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,rsa256(composite),256,SUBTRACT,50000,0.000658611,75917366.208,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,ADDMOD,50000,0.002472267,20224352.623,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,ADDMOD,50000,0.000312223,160141898.958,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,ADDMOD,50000,0.000793337,63024908.984,0
library,Quadro RTX 8000,gpu,cgbn,rsa256(composite),256,ADDMOD,50000,0.000030048,1664004259.851,0
opencl-kernel,Quadro RTX 8000,GPU,w8,rsa256(composite),256,ADDMOD,50000,0.000058892,849012468.619,0
opencl-e2e,Quadro RTX 8000,GPU,w8,rsa256(composite),256,ADDMOD,50000,0.000764492,65402910.825,0
opencl-kernel,Quadro RTX 8000,GPU,w16,rsa256(composite),256,ADDMOD,50000,0.000033454,1494587878.957,0
opencl-e2e,Quadro RTX 8000,GPU,w16,rsa256(composite),256,ADDMOD,50000,0.000683859,73114488.048,0
opencl-kernel,Quadro RTX 8000,GPU,w32,rsa256(composite),256,ADDMOD,50000,0.000020038,2495246096.464,0
opencl-e2e,Quadro RTX 8000,GPU,w32,rsa256(composite),256,ADDMOD,50000,0.000667438,74913335.268,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,rsa256(composite),256,ADDMOD,50000,0.000020269,2466812530.010,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,rsa256(composite),256,ADDMOD,50000,0.000667418,74915556.635,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACTMOD,50000,0.002554132,19576122.874,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACTMOD,50000,0.000299589,166895307.700,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACTMOD,50000,0.000899939,55559315.111,0
library,Quadro RTX 8000,gpu,cgbn,rsa256(composite),256,SUBTRACTMOD,50000,0.000028672,1743861607.143,0
opencl-kernel,Quadro RTX 8000,GPU,w8,rsa256(composite),256,SUBTRACTMOD,50000,0.000063611,786028047.598,0
opencl-e2e,Quadro RTX 8000,GPU,w8,rsa256(composite),256,SUBTRACTMOD,50000,0.000748492,66800974.135,0
opencl-kernel,Quadro RTX 8000,GPU,w16,rsa256(composite),256,SUBTRACTMOD,50000,0.000037822,1321983974.785,0
opencl-e2e,Quadro RTX 8000,GPU,w16,rsa256(composite),256,SUBTRACTMOD,50000,0.000682025,73311097.388,0
opencl-kernel,Quadro RTX 8000,GPU,w32,rsa256(composite),256,SUBTRACTMOD,50000,0.000021882,2284970311.652,0
opencl-e2e,Quadro RTX 8000,GPU,w32,rsa256(composite),256,SUBTRACTMOD,50000,0.000665805,75097081.350,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,rsa256(composite),256,SUBTRACTMOD,50000,0.000020749,2409762161.677,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,rsa256(composite),256,SUBTRACTMOD,50000,0.000664893,75200087.965,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001159573,43119321.619,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000146589,341089692.262,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000284891,175505795.454,0
opencl-kernel,Quadro RTX 8000,GPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001612533,31007115.129,0
opencl-e2e,Quadro RTX 8000,GPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.002476405,20190559.415,0
opencl-kernel,Quadro RTX 8000,GPU,w16,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000436600,114521283.949,0
opencl-e2e,Quadro RTX 8000,GPU,w16,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.001211682,41264956.821,0
opencl-kernel,Quadro RTX 8000,GPU,w32,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000119597,418070853.183,0
opencl-e2e,Quadro RTX 8000,GPU,w32,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000899159,55607510.476,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000040026,1249190650.922,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,50000,0.000808026,61879201.313,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001158701,43151769.959,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000147380,339259208.331,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000281375,177698862.719,0
library,Quadro RTX 8000,gpu,cgbn,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000028800,1736111111.111,0
opencl-kernel,Quadro RTX 8000,GPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000487145,102638837.015,0
opencl-e2e,Quadro RTX 8000,GPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.001366045,36602012.811,0
opencl-kernel,Quadro RTX 8000,GPU,w16,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000148843,335924334.403,0
opencl-e2e,Quadro RTX 8000,GPU,w16,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000926530,53964794.030,0
opencl-kernel,Quadro RTX 8000,GPU,w32,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000044094,1133936862.337,0
opencl-e2e,Quadro RTX 8000,GPU,w32,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000817394,61170002.287,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000060856,821610195.313,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,50000,0.000834616,59907797.028,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.009238168,5412328.564,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.001019126,49061645.153,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000303566,164708833.063,0
library,Quadro RTX 8000,gpu,cgbn,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000027584,1812645011.601,0
opencl-kernel,Quadro RTX 8000,GPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000126911,393977311.137,0
opencl-e2e,Quadro RTX 8000,GPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000843082,59306205.948,0
opencl-kernel,Quadro RTX 8000,GPU,w16,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000037632,1328658182.988,0
opencl-e2e,Quadro RTX 8000,GPU,w16,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000684090,73089777.666,0
opencl-kernel,Quadro RTX 8000,GPU,w32,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000016902,2958210937.543,0
opencl-e2e,Quadro RTX 8000,GPU,w32,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000670103,74615407.130,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000026410,1893223704.487,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,50000,0.000665514,75129922.333,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,COMPARE,50000,0.000377598,132415960.576,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,COMPARE,50000,0.000061828,808692768.970,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,COMPARE,50000,0.000094670,528151136.733,0
library,Quadro RTX 8000,gpu,cgbn,rsa256(composite),256,COMPARE,50000,0.000028672,1743861607.143,0
opencl-kernel,Quadro RTX 8000,GPU,w8,rsa256(composite),256,COMPARE,50000,0.000042070,1188490590.514,0
opencl-e2e,Quadro RTX 8000,GPU,w8,rsa256(composite),256,COMPARE,50000,0.000754594,66260788.033,0
opencl-kernel,Quadro RTX 8000,GPU,w16,rsa256(composite),256,COMPARE,50000,0.000025960,1926028850.742,0
opencl-e2e,Quadro RTX 8000,GPU,w16,rsa256(composite),256,COMPARE,50000,0.000674362,74144147.263,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,rsa256(composite),256,COMPARE,50000,0.000019507,2563210809.133,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,rsa256(composite),256,COMPARE,50000,0.000659613,75802030.345,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,REDUCE,6250,0.000232110,26926911.317,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,REDUCE,6250,0.000034956,178797245.128,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,REDUCE,6250,0.000113306,55160426.183,0
library,Quadro RTX 8000,gpu,cgbn,rsa256(composite),256,REDUCE,50000,0.000028672,1743861607.143,0
opencl-kernel,Quadro RTX 8000,GPU,w8,rsa256(composite),256,REDUCE,50000,0.000220909,226337525.098,0
opencl-e2e,Quadro RTX 8000,GPU,w8,rsa256(composite),256,REDUCE,50000,0.000928895,53827392.298,0
opencl-kernel,Quadro RTX 8000,GPU,w16,rsa256(composite),256,REDUCE,50000,0.000141750,352733462.547,0
opencl-e2e,Quadro RTX 8000,GPU,w16,rsa256(composite),256,REDUCE,50000,0.000784681,63720146.936,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,rsa256(composite),256,REDUCE,50000,0.000071045,703779183.053,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,rsa256(composite),256,REDUCE,50000,0.000706583,70763092.194,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MODMUL,3125,0.000326080,9583537.343,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MODMUL,3125,0.000046949,66561726.607,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MODMUL,3125,0.000128174,24380925.808,0
library,Quadro RTX 8000,gpu,cgbn,rsa256(composite),256,MODMUL,50000,0.000034752,1438766114.180,0
opencl-kernel,Quadro RTX 8000,GPU,w8,rsa256(composite),256,MODMUL,50000,0.000587256,85141754.588,0
opencl-e2e,Quadro RTX 8000,GPU,w8,rsa256(composite),256,MODMUL,50000,0.001295280,38601693.958,0
opencl-kernel,Quadro RTX 8000,GPU,w16,rsa256(composite),256,MODMUL,50000,0.000338052,147906166.374,0
opencl-e2e,Quadro RTX 8000,GPU,w16,rsa256(composite),256,MODMUL,50000,0.000973780,51346290.821,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,rsa256(composite),256,MODMUL,50000,0.000203246,246007279.809,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,rsa256(composite),256,MODMUL,50000,0.000843233,59295594.670,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MODEXP,781,0.007968897,98006.035,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MODEXP,781,0.000913826,854648.431,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MODEXP,781,0.001937340,403130.062,0
library,Quadro RTX 8000,gpu,cgbn,rsa256(composite),256,MODEXP,50000,0.071955517,694873.751,0
opencl-kernel,Quadro RTX 8000,GPU,w8,rsa256(composite),256,MODEXP,50000,0.033662615,1485327.264,0
opencl-e2e,Quadro RTX 8000,GPU,w8,rsa256(composite),256,MODEXP,50000,0.034263004,1459299.956,0
opencl-kernel,Quadro RTX 8000,GPU,w16,rsa256(composite),256,MODEXP,50000,0.005134546,9737959.372,0
opencl-e2e,Quadro RTX 8000,GPU,w16,rsa256(composite),256,MODEXP,50000,0.005908726,8462061.003,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,rsa256(composite),256,MODEXP,50000,0.004726801,10577978.486,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,rsa256(composite),256,MODEXP,50000,0.005376085,9300447.889,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,EXPONENTIATION,781,0.002796102,279317.413,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,EXPONENTIATION,781,0.000389880,2003180.286,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,EXPONENTIATION,781,0.006254852,124863.070,0
opencl-kernel,Quadro RTX 8000,GPU,w8,rsa256(composite),256,EXPONENTIATION,50000,0.051769409,965821.339,0
opencl-e2e,Quadro RTX 8000,GPU,w8,rsa256(composite),256,EXPONENTIATION,50000,0.052417540,953879.179,0
opencl-kernel,Quadro RTX 8000,GPU,w16,rsa256(composite),256,EXPONENTIATION,50000,0.018192456,2748391.955,0
opencl-e2e,Quadro RTX 8000,GPU,w16,rsa256(composite),256,EXPONENTIATION,50000,0.018728785,2669687.331,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,rsa256(composite),256,EXPONENTIATION,50000,0.000881525,56719889.512,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,rsa256(composite),256,EXPONENTIATION,50000,0.001519958,32895645.303,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,DIVIDE,6250,0.000306743,20375365.460,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,DIVIDE,6250,0.000043623,143273318.069,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,DIVIDE,6250,0.000089520,69816704.900,0
library,Quadro RTX 8000,gpu,cgbn,rsa256(composite),256,DIVIDE,50000,0.000029216,1711391018.620,0
opencl-kernel,Quadro RTX 8000,GPU,w8,rsa256(composite),256,DIVIDE,50000,0.000427342,117002301.816,0
opencl-e2e,Quadro RTX 8000,GPU,w8,rsa256(composite),256,DIVIDE,50000,0.001281244,39024573.548,0
opencl-kernel,Quadro RTX 8000,GPU,w16,rsa256(composite),256,DIVIDE,50000,0.000287235,174073484.131,0
opencl-e2e,Quadro RTX 8000,GPU,w16,rsa256(composite),256,DIVIDE,50000,0.001063892,46997255.164,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,rsa256(composite),256,DIVIDE,50000,0.000162589,307523631.599,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,rsa256(composite),256,DIVIDE,50000,0.000927412,53913473.877,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,ISQRT,1562,0.000159513,9792307.020,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,ISQRT,1562,0.000020389,76610013.890,0
opencl-kernel,Quadro RTX 8000,GPU,w8,rsa256(composite),256,ISQRT,50000,0.005346896,9351220.033,0
opencl-e2e,Quadro RTX 8000,GPU,w8,rsa256(composite),256,ISQRT,50000,0.006070632,8236374.700,0
opencl-kernel,Quadro RTX 8000,GPU,w16,rsa256(composite),256,ISQRT,50000,0.004127192,12114774.858,0
opencl-e2e,Quadro RTX 8000,GPU,w16,rsa256(composite),256,ISQRT,50000,0.004773580,10474319.057,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,rsa256(composite),256,ISQRT,50000,0.001322423,37809388.314,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,rsa256(composite),256,ISQRT,50000,0.001969242,25390480.159,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MODMUL_R2,50000,0.005286472,9458103.774,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MODMUL_R2,50000,0.000283969,176075639.555,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MODMUL_R2,50000,0.000720128,69432105.650,0
library,Quadro RTX 8000,gpu,cgbn,rsa256(composite),256,MODMUL_R2,50000,0.000030720,1627604166.667,0
opencl-kernel,Quadro RTX 8000,GPU,w8,rsa256(composite),256,MODMUL_R2,50000,0.000160405,311711085.242,0
opencl-e2e,Quadro RTX 8000,GPU,w8,rsa256(composite),256,MODMUL_R2,50000,0.000838934,59599445.603,0
opencl-kernel,Quadro RTX 8000,GPU,w16,rsa256(composite),256,MODMUL_R2,50000,0.000046879,1066574443.738,0
opencl-e2e,Quadro RTX 8000,GPU,w16,rsa256(composite),256,MODMUL_R2,50000,0.000686624,72820049.454,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,rsa256(composite),256,MODMUL_R2,50000,0.000050977,980832468.428,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,rsa256(composite),256,MODMUL_R2,50000,0.000685342,72956287.657,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,ADD,25000,0.000497114,50290284.241,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,ADD,25000,0.000071646,348936956.564,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,ADD,25000,0.000070154,356359028.107,0
library,Quadro RTX 8000,gpu,cgbn,brainpoolP512r1,512,ADD,50000,0.000033024,1514050387.597,0
opencl-kernel,Quadro RTX 8000,GPU,w8,brainpoolP512r1,512,ADD,50000,0.000116692,428478664.294,0
opencl-e2e,Quadro RTX 8000,GPU,w8,brainpoolP512r1,512,ADD,50000,0.001231609,40597295.169,0
opencl-kernel,Quadro RTX 8000,GPU,w16,brainpoolP512r1,512,ADD,50000,0.000064322,777338693.482,0
opencl-e2e,Quadro RTX 8000,GPU,w16,brainpoolP512r1,512,ADD,50000,0.001126801,44373402.452,0
opencl-kernel,Quadro RTX 8000,GPU,w32,brainpoolP512r1,512,ADD,50000,0.000212644,235134785.870,0
opencl-e2e,Quadro RTX 8000,GPU,w32,brainpoolP512r1,512,ADD,50000,0.001248753,40039946.056,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,brainpoolP512r1,512,ADD,50000,0.000108226,461995329.025,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,brainpoolP512r1,512,ADD,50000,0.001193798,41883127.495,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACT,25000,0.000422653,59150171.461,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACT,25000,0.000083899,297976884.311,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACT,25000,0.005961084,4193868.067,0
library,Quadro RTX 8000,gpu,cgbn,brainpoolP512r1,512,SUBTRACT,50000,0.000033376,1498082454.458,0
opencl-kernel,Quadro RTX 8000,GPU,w8,brainpoolP512r1,512,SUBTRACT,50000,0.000116862,427855475.994,0
opencl-e2e,Quadro RTX 8000,GPU,w8,brainpoolP512r1,512,SUBTRACT,50000,0.001236339,40441979.075,0
opencl-kernel,Quadro RTX 8000,GPU,w16,brainpoolP512r1,512,SUBTRACT,50000,0.000064793,771687789.453,0
opencl-e2e,Quadro RTX 8000,GPU,w16,brainpoolP512r1,512,SUBTRACT,50000,0.001117814,44730166.170,0
opencl-kernel,Quadro RTX 8000,GPU,w32,brainpoolP512r1,512,SUBTRACT,50000,0.000207074,241459590.544,0
opencl-e2e,Quadro RTX 8000,GPU,w32,brainpoolP512r1,512,SUBTRACT,50000,0.001347020,37118972.753,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,brainpoolP512r1,512,SUBTRACT,50000,0.000107404,465532683.570,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,brainpoolP512r1,512,SUBTRACT,50000,0.001136290,44002845.874,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,ADDMOD,25000,0.001333723,18744524.096,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,ADDMOD,25000,0.000133433,187359742.868,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,ADDMOD,25000,0.000375493,66579102.405,0
library,Quadro RTX 8000,gpu,cgbn,brainpoolP512r1,512,ADDMOD,50000,0.000033280,1502403846.154,0
opencl-kernel,Quadro RTX 8000,GPU,w8,brainpoolP512r1,512,ADDMOD,50000,0.000153842,325008459.819,0
opencl-e2e,Quadro RTX 8000,GPU,w8,brainpoolP512r1,512,ADDMOD,50000,0.001306602,38267198.996,0
opencl-kernel,Quadro RTX 8000,GPU,w16,brainpoolP512r1,512,ADDMOD,50000,0.000079451,629317679.053,0
opencl-e2e,Quadro RTX 8000,GPU,w16,brainpoolP512r1,512,ADDMOD,50000,0.001126451,44387196.713,0
opencl-kernel,Quadro RTX 8000,GPU,w32,brainpoolP512r1,512,ADDMOD,50000,0.000218365,228974309.687,0
opencl-e2e,Quadro RTX 8000,GPU,w32,brainpoolP512r1,512,ADDMOD,50000,0.001359393,36781125.170,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,brainpoolP512r1,512,ADDMOD,50000,0.000062449,800651579.877,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,brainpoolP512r1,512,ADDMOD,50000,0.001119517,44662116.774,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACTMOD,25000,0.001333964,18741137.909,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000128995,193806068.635,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACTMOD,25000,0.000573259,43610316.929,0
library,Quadro RTX 8000,gpu,cgbn,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000034432,1452137546.468,0
opencl-kernel,Quadro RTX 8000,GPU,w8,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000175383,285090439.474,0
opencl-e2e,Quadro RTX 8000,GPU,w8,brainpoolP512r1,512,SUBTRACTMOD,50000,0.001248872,40036131.565,0
opencl-kernel,Quadro RTX 8000,GPU,w16,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000085873,582254759.206,0
opencl-e2e,Quadro RTX 8000,GPU,w16,brainpoolP512r1,512,SUBTRACTMOD,50000,0.001230217,40643242.127,0
opencl-kernel,Quadro RTX 8000,GPU,w32,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000223775,223438793.224,0
opencl-e2e,Quadro RTX 8000,GPU,w32,brainpoolP512r1,512,SUBTRACTMOD,50000,0.001367398,36565796.766,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,brainpoolP512r1,512,SUBTRACTMOD,50000,0.000062839,795682587.109,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,brainpoolP512r1,512,SUBTRACTMOD,50000,0.001092566,45763828.038,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.001403015,17818767.508,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.000186224,134246926.707,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,25000,0.000199579,125263720.667,0
opencl-kernel,Quadro RTX 8000,GPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.013188502,3791181.121,0
opencl-e2e,Quadro RTX 8000,GPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.015006536,3331881.488,0
opencl-kernel,Quadro RTX 8000,GPU,w16,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.002989379,16725881.916,0
opencl-e2e,Quadro RTX 8000,GPU,w16,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.004553621,10980272.215,0
opencl-kernel,Quadro RTX 8000,GPU,w32,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.003571366,14000245.441,0
opencl-e2e,Quadro RTX 8000,GPU,w32,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.004960145,10080350.265,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.000152139,328646844.354,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,50000,0.001432181,34911791.545,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.001402975,17819276.131,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000172688,144769623.967,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,25000,0.000195291,128014009.101,0
library,Quadro RTX 8000,gpu,cgbn,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000033184,1506750241.080,0
opencl-kernel,Quadro RTX 8000,GPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.003416330,14635588.106,0
opencl-e2e,Quadro RTX 8000,GPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.004727790,10575766.069,0
opencl-kernel,Quadro RTX 8000,GPU,w16,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000898808,55629232.886,0
opencl-e2e,Quadro RTX 8000,GPU,w16,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.002194870,22780392.144,0
opencl-kernel,Quadro RTX 8000,GPU,w32,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000699529,71476669.238,0
opencl-e2e,Quadro RTX 8000,GPU,w32,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.002100581,23802937.532,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.000614838,81322237.626,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,50000,0.002022323,24724043.017,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.010771160,2321012.776,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.001143212,21868206.146,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,25000,0.000342751,72939264.074,0
library,Quadro RTX 8000,gpu,cgbn,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000032928,1518464528.669,0
opencl-kernel,Quadro RTX 8000,GPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000482717,103580345.100,0
opencl-e2e,Quadro RTX 8000,GPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.001595881,31330660.075,0
opencl-kernel,Quadro RTX 8000,GPU,w16,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000139966,357229846.178,0
opencl-e2e,Quadro RTX 8000,GPU,w16,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.001269412,39388315.402,0
opencl-kernel,Quadro RTX 8000,GPU,w32,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000058672,852196341.183,0
opencl-e2e,Quadro RTX 8000,GPU,w32,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.001100491,45434263.852,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.000079902,625767433.620,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,50000,0.001221121,40945981.362,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,COMPARE,25000,0.000198988,125635709.714,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,COMPARE,25000,0.000034495,724742211.873,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,COMPARE,25000,0.000064051,390313898.007,0
library,Quadro RTX 8000,gpu,cgbn,brainpoolP512r1,512,COMPARE,50000,0.000033024,1514050387.597,0
opencl-kernel,Quadro RTX 8000,GPU,w8,brainpoolP512r1,512,COMPARE,50000,0.000096222,519630859.170,0
opencl-e2e,Quadro RTX 8000,GPU,w8,brainpoolP512r1,512,COMPARE,50000,0.001156777,43223545.992,0
opencl-kernel,Quadro RTX 8000,GPU,w16,brainpoolP512r1,512,COMPARE,50000,0.000055806,895960801.886,0
opencl-e2e,Quadro RTX 8000,GPU,w16,brainpoolP512r1,512,COMPARE,50000,0.001104338,45275999.068,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,brainpoolP512r1,512,COMPARE,50000,0.000042891,1165748526.450,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,brainpoolP512r1,512,COMPARE,50000,0.001166106,42877749.822,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,REDUCE,3125,0.000123254,25354141.889,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,REDUCE,3125,0.000024136,129475056.674,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,REDUCE,3125,0.000079662,39228203.165,0
library,Quadro RTX 8000,gpu,cgbn,brainpoolP512r1,512,REDUCE,50000,0.000033088,1511121856.867,0
opencl-kernel,Quadro RTX 8000,GPU,w8,brainpoolP512r1,512,REDUCE,50000,0.000678349,73708390.842,0
opencl-e2e,Quadro RTX 8000,GPU,w8,brainpoolP512r1,512,REDUCE,50000,0.001778338,28116141.931,0
opencl-kernel,Quadro RTX 8000,GPU,w16,brainpoolP512r1,512,REDUCE,50000,0.000520809,96004512.045,0
opencl-e2e,Quadro RTX 8000,GPU,w16,brainpoolP512r1,512,REDUCE,50000,0.001560064,32049966.510,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,brainpoolP512r1,512,REDUCE,50000,0.000188148,265748078.556,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,brainpoolP512r1,512,REDUCE,50000,0.001255596,39821731.686,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL,1562,0.000356327,4383613.388,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL,1562,0.000044355,35215737.813,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL,1562,0.000118245,13209871.492,0
library,Quadro RTX 8000,gpu,cgbn,brainpoolP512r1,512,MODMUL,50000,0.000096256,519448138.298,0
opencl-kernel,Quadro RTX 8000,GPU,w8,brainpoolP512r1,512,MODMUL,50000,0.001883438,26547199.215,0
opencl-e2e,Quadro RTX 8000,GPU,w8,brainpoolP512r1,512,MODMUL,50000,0.002970934,16829724.741,0
opencl-kernel,Quadro RTX 8000,GPU,w16,brainpoolP512r1,512,MODMUL,50000,0.001826861,27369351.904,0
opencl-e2e,Quadro RTX 8000,GPU,w16,brainpoolP512r1,512,MODMUL,50000,0.002879912,17361641.403,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,brainpoolP512r1,512,MODMUL,50000,0.000691524,72304069.563,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,brainpoolP512r1,512,MODMUL,50000,0.001817403,27511786.254,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MODEXP,390,0.021668981,17998.077,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MODEXP,390,0.002985862,130615.544,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MODEXP,390,0.002980852,130835.075,0
library,Quadro RTX 8000,gpu,cgbn,brainpoolP512r1,512,MODEXP,50000,0.604014993,82779.402,0
opencl-kernel,Quadro RTX 8000,GPU,w8,brainpoolP512r1,512,MODEXP,50000,0.449828984,111153.353,0
opencl-e2e,Quadro RTX 8000,GPU,w8,brainpoolP512r1,512,MODEXP,50000,0.452650514,110460.495,0
opencl-kernel,Quadro RTX 8000,GPU,w16,brainpoolP512r1,512,MODEXP,50000,0.055630518,898787.247,0
opencl-e2e,Quadro RTX 8000,GPU,w16,brainpoolP512r1,512,MODEXP,50000,0.056598547,883414.902,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,brainpoolP512r1,512,MODEXP,50000,0.047663960,1049010.612,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,brainpoolP512r1,512,MODEXP,50000,0.048205689,1037221.980,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,EXPONENTIATION,390,0.004870812,80068.786,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,EXPONENTIATION,390,0.000735597,530181.617,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,EXPONENTIATION,390,0.008966543,43495.024,0
opencl-kernel,Quadro RTX 8000,GPU,w8,brainpoolP512r1,512,EXPONENTIATION,50000,0.441202231,113326.716,0
opencl-e2e,Quadro RTX 8000,GPU,w8,brainpoolP512r1,512,EXPONENTIATION,50000,0.442941035,112881.842,0
opencl-kernel,Quadro RTX 8000,GPU,w16,brainpoolP512r1,512,EXPONENTIATION,50000,0.123422743,405111.722,0
opencl-e2e,Quadro RTX 8000,GPU,w16,brainpoolP512r1,512,EXPONENTIATION,50000,0.124419227,401867.149,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,brainpoolP512r1,512,EXPONENTIATION,50000,0.045765483,1092526.440,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,brainpoolP512r1,512,EXPONENTIATION,50000,0.047250344,1058193.355,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,DIVIDE,3125,0.000161587,19339423.727,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,DIVIDE,3125,0.000027392,114085125.843,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,DIVIDE,3125,0.000086013,36331719.465,0
library,Quadro RTX 8000,gpu,cgbn,brainpoolP512r1,512,DIVIDE,50000,0.000048896,1022578534.031,0
opencl-kernel,Quadro RTX 8000,GPU,w8,brainpoolP512r1,512,DIVIDE,50000,0.001287816,38825419.480,0
opencl-e2e,Quadro RTX 8000,GPU,w8,brainpoolP512r1,512,DIVIDE,50000,0.002603726,19203249.528,0
opencl-kernel,Quadro RTX 8000,GPU,w16,brainpoolP512r1,512,DIVIDE,50000,0.001183889,42233693.599,0
opencl-e2e,Quadro RTX 8000,GPU,w16,brainpoolP512r1,512,DIVIDE,50000,0.002492185,20062716.183,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,brainpoolP512r1,512,DIVIDE,50000,0.000383328,130436615.097,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,brainpoolP512r1,512,DIVIDE,50000,0.001775684,28158162.199,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,ISQRT,781,0.000141729,5510511.298,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,ISQRT,781,0.000023234,33614621.433,0
opencl-kernel,Quadro RTX 8000,GPU,w8,brainpoolP512r1,512,ISQRT,50000,0.036823489,1357828.966,0
opencl-e2e,Quadro RTX 8000,GPU,w8,brainpoolP512r1,512,ISQRT,50000,0.037808751,1322445.167,0
opencl-kernel,Quadro RTX 8000,GPU,w16,brainpoolP512r1,512,ISQRT,50000,0.033031319,1513714.906,0
opencl-e2e,Quadro RTX 8000,GPU,w16,brainpoolP512r1,512,ISQRT,50000,0.034039465,1468883.250,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,brainpoolP512r1,512,ISQRT,50000,0.006374693,7843514.962,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,brainpoolP512r1,512,ISQRT,50000,0.007866638,6355955.354,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL_R2,25000,0.005305578,4712021.950,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL_R2,25000,0.000255816,97726522.384,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL_R2,25000,0.000874311,28593945.414,0
library,Quadro RTX 8000,gpu,cgbn,brainpoolP512r1,512,MODMUL_R2,50000,0.000039520,1265182186.235,0
opencl-kernel,Quadro RTX 8000,GPU,w8,brainpoolP512r1,512,MODMUL_R2,50000,0.000674351,74145350.438,0
opencl-e2e,Quadro RTX 8000,GPU,w8,brainpoolP512r1,512,MODMUL_R2,50000,0.001826710,27371612.433,0
opencl-kernel,Quadro RTX 8000,GPU,w16,brainpoolP512r1,512,MODMUL_R2,50000,0.000150005,333322517.097,0
opencl-e2e,Quadro RTX 8000,GPU,w16,brainpoolP512r1,512,MODMUL_R2,50000,0.001201724,41606893.680,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,brainpoolP512r1,512,MODMUL_R2,50000,0.000112825,443163904.412,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,brainpoolP512r1,512,MODMUL_R2,50000,0.001160665,43078754.139,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p1024,1024,ADD,12500,0.000293318,42615871.363,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p1024,1024,ADD,12500,0.000040036,312217751.259,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p1024,1024,ADD,12500,0.000053893,231941742.271,0
library,Quadro RTX 8000,gpu,cgbn,p1024,1024,ADD,50000,0.000053024,942969221.485,0
opencl-kernel,Quadro RTX 8000,GPU,w8,p1024,1024,ADD,50000,0.000680783,73444835.173,0
opencl-e2e,Quadro RTX 8000,GPU,w8,p1024,1024,ADD,50000,0.002556486,19558096.121,0
opencl-kernel,Quadro RTX 8000,GPU,w16,p1024,1024,ADD,50000,0.000339976,147069285.604,0
opencl-e2e,Quadro RTX 8000,GPU,w16,p1024,1024,ADD,50000,0.002207895,22646006.406,0
opencl-kernel,Quadro RTX 8000,GPU,w32,p1024,1024,ADD,50000,0.000536189,93250709.666,0
opencl-e2e,Quadro RTX 8000,GPU,w32,p1024,1024,ADD,50000,0.002371506,21083648.360,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,p1024,1024,ADD,50000,0.000214287,233332099.896,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,p1024,1024,ADD,50000,0.002051108,24377067.787,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p1024,1024,SUBTRACT,12500,0.000235518,53074505.260,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p1024,1024,SUBTRACT,12500,0.000045566,274326620.168,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p1024,1024,SUBTRACT,12500,0.000066647,187555122.674,0
library,Quadro RTX 8000,gpu,cgbn,p1024,1024,SUBTRACT,50000,0.000053152,940698374.473,0
opencl-kernel,Quadro RTX 8000,GPU,w8,p1024,1024,SUBTRACT,50000,0.000660845,75660698.134,0
opencl-e2e,Quadro RTX 8000,GPU,w8,p1024,1024,SUBTRACT,50000,0.002509027,19928045.502,0
opencl-kernel,Quadro RTX 8000,GPU,w16,p1024,1024,SUBTRACT,50000,0.000350245,142757195.782,0
opencl-e2e,Quadro RTX 8000,GPU,w16,p1024,1024,SUBTRACT,50000,0.002162528,23121088.003,0
opencl-kernel,Quadro RTX 8000,GPU,w32,p1024,1024,SUBTRACT,50000,0.000537090,93094266.650,0
opencl-e2e,Quadro RTX 8000,GPU,w32,p1024,1024,SUBTRACT,50000,0.002399109,20841069.963,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,p1024,1024,SUBTRACT,50000,0.000209498,238665778.456,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,p1024,1024,SUBTRACT,50000,0.002114408,23647281.489,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p1024,1024,ADDMOD,12500,0.000878849,14223147.606,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p1024,1024,ADDMOD,12500,0.000109528,114126111.403,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p1024,1024,ADDMOD,12500,0.000365774,34174115.508,0
library,Quadro RTX 8000,gpu,cgbn,p1024,1024,ADDMOD,50000,0.000053248,939002403.846,0
opencl-kernel,Quadro RTX 8000,GPU,w8,p1024,1024,ADDMOD,50000,0.000705831,70838488.270,0
opencl-e2e,Quadro RTX 8000,GPU,w8,p1024,1024,ADDMOD,50000,0.002672196,18711203.021,0
opencl-kernel,Quadro RTX 8000,GPU,w16,p1024,1024,ADDMOD,50000,0.000382016,130884588.352,0
opencl-e2e,Quadro RTX 8000,GPU,w16,p1024,1024,ADDMOD,50000,0.002273439,21993112.985,0
opencl-kernel,Quadro RTX 8000,GPU,w32,p1024,1024,ADDMOD,50000,0.000530418,94265276.051,0
opencl-e2e,Quadro RTX 8000,GPU,w32,p1024,1024,ADDMOD,50000,0.002394540,20880837.216,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,p1024,1024,ADDMOD,50000,0.000137060,364803870.433,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,p1024,1024,ADDMOD,50000,0.002010351,24871280.997,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p1024,1024,SUBTRACTMOD,12500,0.000790452,15813739.229,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p1024,1024,SUBTRACTMOD,12500,0.006087263,2053468.076,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p1024,1024,SUBTRACTMOD,12500,0.000370644,33725080.909,0
library,Quadro RTX 8000,gpu,cgbn,p1024,1024,SUBTRACTMOD,50000,0.000052992,943538647.343,0
opencl-kernel,Quadro RTX 8000,GPU,w8,p1024,1024,SUBTRACTMOD,50000,0.000903677,55329506.091,0
opencl-e2e,Quadro RTX 8000,GPU,w8,p1024,1024,SUBTRACTMOD,50000,0.002760754,18110994.533,0
opencl-kernel,Quadro RTX 8000,GPU,w16,p1024,1024,SUBTRACTMOD,50000,0.000384701,129971085.253,0
opencl-e2e,Quadro RTX 8000,GPU,w16,p1024,1024,SUBTRACTMOD,50000,0.002212053,22603437.412,0
opencl-kernel,Quadro RTX 8000,GPU,w32,p1024,1024,SUBTRACTMOD,50000,0.000532883,93829188.418,0
opencl-e2e,Quadro RTX 8000,GPU,w32,p1024,1024,SUBTRACTMOD,50000,0.002377017,21034768.097,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,p1024,1024,SUBTRACTMOD,50000,0.000138594,360765794.942,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,p1024,1024,SUBTRACTMOD,50000,0.001998798,25015034.156,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.002397314,5214168.973,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.000226881,55094988.537,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p1024,1024,MULTIPLYOPERANDSCANNING,12500,0.000381785,32740942.361,0
opencl-kernel,Quadro RTX 8000,GPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.212960670,234785.137,0
opencl-e2e,Quadro RTX 8000,GPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.216165498,231304.258,0
opencl-kernel,Quadro RTX 8000,GPU,w16,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.060374532,828163.771,0
opencl-e2e,Quadro RTX 8000,GPU,w16,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.062737331,796973.655,0
opencl-kernel,Quadro RTX 8000,GPU,w32,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.017716584,2822214.481,0
opencl-e2e,Quadro RTX 8000,GPU,w32,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.020085757,2489326.155,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.000317063,157697391.126,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.002651288,18858758.878,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.002394178,5220998.205,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.000309338,40408860.776,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,12500,0.000382366,32691200.421,0
library,Quadro RTX 8000,gpu,cgbn,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000052800,946969696.970,0
opencl-kernel,Quadro RTX 8000,GPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.027564461,1813929.893,0
opencl-e2e,Quadro RTX 8000,GPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.029955965,1669116.652,0
opencl-kernel,Quadro RTX 8000,GPU,w16,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.006983087,7160157.175,0
opencl-e2e,Quadro RTX 8000,GPU,w16,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.009528303,5247524.110,0
opencl-kernel,Quadro RTX 8000,GPU,w32,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.002524256,19807816.579,0
opencl-e2e,Quadro RTX 8000,GPU,w32,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.005022764,9954678.394,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.001871387,26718150.640,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.004215781,11860198.962,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.016241422,769637.047,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.001728884,7230097.798,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p1024,1024,MONTGOMERYMULTIPLICATION,12500,0.005979208,2090577.949,0
library,Quadro RTX 8000,gpu,cgbn,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000060352,828472958.643,0
opencl-kernel,Quadro RTX 8000,GPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.002573869,19426008.398,0
opencl-e2e,Quadro RTX 8000,GPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.004506110,11096045.195,0
opencl-kernel,Quadro RTX 8000,GPU,w16,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000522683,95660292.591,0
opencl-e2e,Quadro RTX 8000,GPU,w16,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.002341750,21351554.738,0
opencl-kernel,Quadro RTX 8000,GPU,w32,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000143092,349425806.126,0
opencl-e2e,Quadro RTX 8000,GPU,w32,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.002148853,23268225.666,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000248301,201368438.879,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.002121221,23571329.450,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p1024,1024,COMPARE,12500,0.000095221,131273589.570,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p1024,1024,COMPARE,12500,0.000019076,655272011.815,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p1024,1024,COMPARE,12500,0.000038974,320726744.408,0
library,Quadro RTX 8000,gpu,cgbn,p1024,1024,COMPARE,50000,0.000052704,948694596.236,0
opencl-kernel,Quadro RTX 8000,GPU,w8,p1024,1024,COMPARE,50000,0.000235177,212605724.310,0
opencl-e2e,Quadro RTX 8000,GPU,w8,p1024,1024,COMPARE,50000,0.002094548,23871499.192,0
opencl-kernel,Quadro RTX 8000,GPU,w16,p1024,1024,COMPARE,50000,0.000126901,394007670.686,0
opencl-e2e,Quadro RTX 8000,GPU,w16,p1024,1024,COMPARE,50000,0.001989250,25135100.267,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,p1024,1024,COMPARE,50000,0.000092065,543094778.158,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,p1024,1024,COMPARE,50000,0.001943483,25727006.386,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p1024,1024,REDUCE,1562,0.000041950,37234778.139,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p1024,1024,REDUCE,1562,0.000014277,109407180.749,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p1024,1024,REDUCE,1562,0.000065795,23740353.151,0
library,Quadro RTX 8000,gpu,cgbn,p1024,1024,REDUCE,50000,0.000053248,939002403.846,0
opencl-kernel,Quadro RTX 8000,GPU,w8,p1024,1024,REDUCE,50000,0.003470663,14406469.389,0
opencl-e2e,Quadro RTX 8000,GPU,w8,p1024,1024,REDUCE,50000,0.005292653,9447057.836,0
opencl-kernel,Quadro RTX 8000,GPU,w16,p1024,1024,REDUCE,50000,0.001905991,26233072.491,0
opencl-e2e,Quadro RTX 8000,GPU,w16,p1024,1024,REDUCE,50000,0.003925358,12737691.118,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,p1024,1024,REDUCE,50000,0.000650055,76916567.650,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,p1024,1024,REDUCE,50000,0.002502435,19980537.998,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p1024,1024,MODMUL,781,0.000470193,1661019.686,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p1024,1024,MODMUL,781,0.000060766,12852581.184,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p1024,1024,MODMUL,781,0.000161417,4838398.078,0
library,Quadro RTX 8000,gpu,cgbn,p1024,1024,MODMUL,50000,0.000283264,176513782.196,0
opencl-kernel,Quadro RTX 8000,GPU,w8,p1024,1024,MODMUL,50000,0.017717636,2822046.921,0
opencl-e2e,Quadro RTX 8000,GPU,w8,p1024,1024,MODMUL,50000,0.019751751,2531421.126,0
opencl-kernel,Quadro RTX 8000,GPU,w16,p1024,1024,MODMUL,50000,0.005212353,9592596.523,0
opencl-e2e,Quadro RTX 8000,GPU,w16,p1024,1024,MODMUL,50000,0.007453692,6708085.146,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,p1024,1024,MODMUL,50000,0.002976967,16795617.699,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,p1024,1024,MODMUL,50000,0.004816923,10380070.307,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p1024,1024,MODEXP,195,0.071865662,2713.396,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p1024,1024,MODEXP,195,0.007026938,27750.351,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p1024,1024,MODEXP,195,0.005906990,33011.737,0
library,Quadro RTX 8000,gpu,cgbn,p1024,1024,MODEXP,50000,1.516476393,32971.169,0
opencl-kernel,Quadro RTX 8000,GPU,w8,p1024,1024,MODEXP,50000,4.237869749,11798.381,0
opencl-e2e,Quadro RTX 8000,GPU,w8,p1024,1024,MODEXP,50000,4.276200058,11692.624,0
opencl-kernel,Quadro RTX 8000,GPU,w16,p1024,1024,MODEXP,50000,0.637082738,78482.742,0
opencl-e2e,Quadro RTX 8000,GPU,w16,p1024,1024,MODEXP,50000,0.640883328,78017.321,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,p1024,1024,MODEXP,50000,0.357484475,139866.214,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,p1024,1024,MODEXP,50000,0.359954639,138906.391,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p1024,1024,EXPONENTIATION,195,0.010461612,18639.575,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p1024,1024,EXPONENTIATION,195,0.000860544,226600.868,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p1024,1024,EXPONENTIATION,195,0.014110404,13819.590,0
opencl-kernel,Quadro RTX 8000,GPU,w8,p1024,1024,EXPONENTIATION,50000,4.293255999,11646.173,0
opencl-e2e,Quadro RTX 8000,GPU,w8,p1024,1024,EXPONENTIATION,50000,4.301126532,11624.861,0
opencl-kernel,Quadro RTX 8000,GPU,w16,p1024,1024,EXPONENTIATION,50000,0.974911446,51286.709,0
opencl-e2e,Quadro RTX 8000,GPU,w16,p1024,1024,EXPONENTIATION,50000,0.975195035,51271.795,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,p1024,1024,EXPONENTIATION,50000,0.318696902,156888.880,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,p1024,1024,EXPONENTIATION,50000,0.319693487,156399.808,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p1024,1024,DIVIDE,1562,0.000088779,17594246.351,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p1024,1024,DIVIDE,1562,0.000013185,118468256.836,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p1024,1024,DIVIDE,1562,0.000035187,44391398.733,0
library,Quadro RTX 8000,gpu,cgbn,p1024,1024,DIVIDE,50000,0.000069568,718721251.150,0
opencl-kernel,Quadro RTX 8000,GPU,w8,p1024,1024,DIVIDE,50000,0.150001700,333329.556,0
opencl-e2e,Quadro RTX 8000,GPU,w8,p1024,1024,DIVIDE,50000,0.152640222,327567.658,0
opencl-kernel,Quadro RTX 8000,GPU,w16,p1024,1024,DIVIDE,50000,0.014190702,3523433.865,0
opencl-e2e,Quadro RTX 8000,GPU,w16,p1024,1024,DIVIDE,50000,0.016461588,3037374.041,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,p1024,1024,DIVIDE,50000,0.002022092,24726867.054,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,p1024,1024,DIVIDE,50000,0.004391355,11386006.743,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p1024,1024,ISQRT,390,0.000138914,2807492.815,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p1024,1024,ISQRT,390,0.000015009,25984475.520,0
opencl-kernel,Quadro RTX 8000,GPU,w8,p1024,1024,ISQRT,50000,1.733258285,28847.403,0
opencl-e2e,Quadro RTX 8000,GPU,w8,p1024,1024,ISQRT,50000,1.735218852,28814.809,0
opencl-kernel,Quadro RTX 8000,GPU,w16,p1024,1024,ISQRT,50000,0.269173608,185753.724,0
opencl-e2e,Quadro RTX 8000,GPU,w16,p1024,1024,ISQRT,50000,0.272754772,183314.850,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,p1024,1024,ISQRT,50000,0.044776123,1116666.577,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,p1024,1024,ISQRT,50000,0.046523864,1074717.264,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p1024,1024,MODMUL_R2,12500,0.007345173,1701797.946,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p1024,1024,MODMUL_R2,12500,0.000349774,35737359.246,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p1024,1024,MODMUL_R2,12500,0.000884040,14139634.999,0
library,Quadro RTX 8000,gpu,cgbn,p1024,1024,MODMUL_R2,50000,0.000109952,454743888.242,0
opencl-kernel,Quadro RTX 8000,GPU,w8,p1024,1024,MODMUL_R2,50000,0.003671635,13617911.580,0
opencl-e2e,Quadro RTX 8000,GPU,w8,p1024,1024,MODMUL_R2,50000,0.005796983,8625176.026,0
opencl-kernel,Quadro RTX 8000,GPU,w16,p1024,1024,MODMUL_R2,50000,0.000598697,83514700.566,0
opencl-e2e,Quadro RTX 8000,GPU,w16,p1024,1024,MODMUL_R2,50000,0.002429506,20580315.140,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,p1024,1024,MODMUL_R2,50000,0.000358241,139570881.945,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,p1024,1024,MODMUL_R2,50000,0.002195242,22776531.181,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p2048,2048,ADD,6250,0.000201793,30972345.026,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p2048,2048,ADD,6250,0.000030327,206086151.673,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p2048,2048,ADD,6250,0.000059162,105642074.939,0
library,Quadro RTX 8000,gpu,cgbn,p2048,2048,ADD,50000,0.000100352,498246173.469,0
opencl-kernel,Quadro RTX 8000,GPU,w8,p2048,2048,ADD,50000,0.001572788,31790680.101,0
opencl-e2e,Quadro RTX 8000,GPU,w8,p2048,2048,ADD,50000,0.005266164,9494577.143,0
opencl-kernel,Quadro RTX 8000,GPU,w16,p2048,2048,ADD,50000,0.000767388,65156095.107,0
opencl-e2e,Quadro RTX 8000,GPU,w16,p2048,2048,ADD,50000,0.004470975,11183242.683,0
opencl-kernel,Quadro RTX 8000,GPU,w32,p2048,2048,ADD,50000,0.000423104,118174309.838,0
opencl-e2e,Quadro RTX 8000,GPU,w32,p2048,2048,ADD,50000,0.004665125,10717826.589,0
opencl-kernel,Quadro RTX 8000,GPU,w32-opt,p2048,2048,ADD,50000,0.000000000,inf,0
opencl-e2e,Quadro RTX 8000,GPU,w32-opt,p2048,2048,ADD,50000,0.000000000,inf,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p2048,2048,SUBTRACT,6250,0.000168300,37136066.338,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p2048,2048,SUBTRACT,6250,0.000031950,195618445.753,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p2048,2048,SUBTRACT,6250,0.000068320,91481316.285,0
library,Quadro RTX 8000,gpu,cgbn,p2048,2048,SUBTRACT,50000,0.000100416,497928616.953,0
opencl-kernel,Quadro RTX 8000,GPU,w8,p2048,2048,SUBTRACT,50000,0.001594179,31364105.024,0
opencl-e2e,Quadro RTX 8000,GPU,w8,p2048,2048,SUBTRACT,50000,0.005277326,9474495.052,0
opencl-kernel,Quadro RTX 8000,GPU,w16,p2048,2048,SUBTRACT,50000,0.000755405,66189675.565,0
opencl-e2e,Quadro RTX 8000,GPU,w16,p2048,2048,SUBTRACT,50000,0.004453742,11226514.830,0
opencl-kernel,Quadro RTX 8000,GPU,w32,p2048,2048,SUBTRACT,50000,0.000421090,118739512.996,0
opencl-e2e,Quadro RTX 8000,GPU,w32,p2048,2048,SUBTRACT,50000,0.004632613,10793044.512,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p2048,2048,ADDMOD,6250,0.000563441,11092557.563,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p2048,2048,ADDMOD,6250,0.000074521,83869044.100,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p2048,2048,ADDMOD,6250,0.000283268,22063907.508,0
library,Quadro RTX 8000,gpu,cgbn,p2048,2048,ADDMOD,50000,0.000099616,501927401.221,0
opencl-kernel,Quadro RTX 8000,GPU,w8,p2048,2048,ADDMOD,50000,0.001863310,26833969.227,0
opencl-e2e,Quadro RTX 8000,GPU,w8,p2048,2048,ADDMOD,50000,0.005468920,9142573.151,0
opencl-kernel,Quadro RTX 8000,GPU,w16,p2048,2048,ADDMOD,50000,0.000951247,52562593.359,0
opencl-e2e,Quadro RTX 8000,GPU,w16,p2048,2048,ADDMOD,50000,0.004646639,10760466.082,0
opencl-kernel,Quadro RTX 8000,GPU,w32,p2048,2048,ADDMOD,50000,0.000486615,102750659.595,0
opencl-e2e,Quadro RTX 8000,GPU,w32,p2048,2048,ADDMOD,50000,0.004717003,10599950.660,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p2048,2048,SUBTRACTMOD,6250,0.000482747,12946737.578,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p2048,2048,SUBTRACTMOD,6250,0.000070604,88522151.028,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p2048,2048,SUBTRACTMOD,6250,0.000283999,22007109.227,0
library,Quadro RTX 8000,gpu,cgbn,p2048,2048,SUBTRACTMOD,50000,0.000100192,499041839.668,0
opencl-kernel,Quadro RTX 8000,GPU,w8,p2048,2048,SUBTRACTMOD,50000,0.001906051,26232245.740,0
opencl-e2e,Quadro RTX 8000,GPU,w8,p2048,2048,SUBTRACTMOD,50000,0.005588537,8946885.252,0
opencl-kernel,Quadro RTX 8000,GPU,w16,p2048,2048,SUBTRACTMOD,50000,0.000939104,53242238.983,0
opencl-e2e,Quadro RTX 8000,GPU,w16,p2048,2048,SUBTRACTMOD,50000,0.004646459,10760882.344,0
opencl-kernel,Quadro RTX 8000,GPU,w32,p2048,2048,SUBTRACTMOD,50000,0.000515319,97027289.297,0
opencl-e2e,Quadro RTX 8000,GPU,w32,p2048,2048,SUBTRACTMOD,50000,0.004742391,10543205.023,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.004070532,1535425.820,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.000467689,13363581.086,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p2048,2048,MULTIPLYOPERANDSCANNING,6250,0.006003594,1041043.077,0
opencl-kernel,Quadro RTX 8000,GPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,50000,1.027948603,48640.564,0
opencl-e2e,Quadro RTX 8000,GPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,50000,1.032911182,48406.873,0
opencl-kernel,Quadro RTX 8000,GPU,w16,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.279068280,179167.622,0
opencl-e2e,Quadro RTX 8000,GPU,w16,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.283363189,176451.995,0
opencl-kernel,Quadro RTX 8000,GPU,w32,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.053872517,928117.020,0
opencl-e2e,Quadro RTX 8000,GPU,w32,p2048,2048,MULTIPLYOPERANDSCANNING,50000,0.058911372,848732.568,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.004068919,1536034.512,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.000366005,17076274.938,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,6250,0.000444074,14074233.137,0
library,Quadro RTX 8000,gpu,cgbn,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.000145984,342503288.032,0
opencl-kernel,Quadro RTX 8000,GPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.118799448,420877.376,0
opencl-e2e,Quadro RTX 8000,GPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.124082295,402958.375,0
opencl-kernel,Quadro RTX 8000,GPU,w16,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.031734918,1575551.573,0
opencl-e2e,Quadro RTX 8000,GPU,w16,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.036956078,1352957.416,0
opencl-kernel,Quadro RTX 8000,GPU,w32,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.008167550,6121786.851,0
opencl-e2e,Quadro RTX 8000,GPU,w32,p2048,2048,MULTIPLYPRODUCTSCANNING,50000,0.013481878,3708682.135,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.026490290,235935.507,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.002204006,2835745.363,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p2048,2048,MONTGOMERYMULTIPLICATION,6250,0.001102374,5669582.525,0
library,Quadro RTX 8000,gpu,cgbn,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000261600,191131498.471,0
opencl-kernel,Quadro RTX 8000,GPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.139848278,357530.323,0
opencl-e2e,Quadro RTX 8000,GPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.143140132,349308.048,0
opencl-kernel,Quadro RTX 8000,GPU,w16,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.002150666,23248612.629,0
opencl-e2e,Quadro RTX 8000,GPU,w16,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.005956497,8394195.386,0
opencl-kernel,Quadro RTX 8000,GPU,w32,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.000391003,127876264.139,0
opencl-e2e,Quadro RTX 8000,GPU,w32,p2048,2048,MONTGOMERYMULTIPLICATION,50000,0.004610521,10844761.314,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p2048,2048,COMPARE,6250,0.000048071,130016252.791,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p2048,2048,COMPARE,6250,0.000017243,362466520.835,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p2048,2048,COMPARE,6250,0.000053221,117434566.875,0
library,Quadro RTX 8000,gpu,cgbn,p2048,2048,COMPARE,50000,0.000100032,499840051.184,0
opencl-kernel,Quadro RTX 8000,GPU,w8,p2048,2048,COMPARE,50000,0.000477928,104618275.478,0
opencl-e2e,Quadro RTX 8000,GPU,w8,p2048,2048,COMPARE,50000,0.004242110,11786587.151,0
opencl-kernel,Quadro RTX 8000,GPU,w16,p2048,2048,COMPARE,50000,0.000322203,155181598.032,0
opencl-e2e,Quadro RTX 8000,GPU,w16,p2048,2048,COMPARE,50000,0.004110440,12164148.134,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p2048,2048,REDUCE,781,0.000032442,24073788.428,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p2048,2048,REDUCE,781,0.000014587,53540557.344,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p2048,2048,REDUCE,781,0.000067468,11575852.333,0
library,Quadro RTX 8000,gpu,cgbn,p2048,2048,REDUCE,50000,0.000099104,504520503.713,0
opencl-kernel,Quadro RTX 8000,GPU,w8,p2048,2048,REDUCE,50000,1.103548834,45308.371,0
opencl-e2e,Quadro RTX 8000,GPU,w8,p2048,2048,REDUCE,50000,1.108026324,45125.282,0
opencl-kernel,Quadro RTX 8000,GPU,w16,p2048,2048,REDUCE,50000,0.010119417,4940996.115,0
opencl-e2e,Quadro RTX 8000,GPU,w16,p2048,2048,REDUCE,50000,0.013730819,3641443.356,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p2048,2048,MODMUL,390,0.000707193,551476.012,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p2048,2048,MODMUL,390,0.000083809,4653449.992,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p2048,2048,MODMUL,390,0.000223414,1745637.821,0
library,Quadro RTX 8000,gpu,cgbn,p2048,2048,MODMUL,50000,0.001315840,37998540.856,0
opencl-kernel,Quadro RTX 8000,GPU,w8,p2048,2048,MODMUL,50000,1.906933753,26220.103,0
opencl-e2e,Quadro RTX 8000,GPU,w8,p2048,2048,MODMUL,50000,1.910398239,26172.553,0
opencl-kernel,Quadro RTX 8000,GPU,w16,p2048,2048,MODMUL,50000,0.031422193,1591232.035,0
opencl-e2e,Quadro RTX 8000,GPU,w16,p2048,2048,MODMUL,50000,0.035258442,1418100.099,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p2048,2048,MODEXP,97,0.257221635,377.107,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p2048,2048,MODEXP,97,0.025677636,3777.606,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p2048,2048,MODEXP,97,0.019657340,4934.544,0
library,Quadro RTX 8000,gpu,cgbn,p2048,2048,MODEXP,50000,2.394568205,20880.591,0
opencl-kernel,Quadro RTX 8000,GPU,w8,p2048,2048,MODEXP,50000,251.407893080,198.880,0
opencl-e2e,Quadro RTX 8000,GPU,w8,p2048,2048,MODEXP,50000,251.428107239,198.864,0
opencl-kernel,Quadro RTX 8000,GPU,w16,p2048,2048,MODEXP,50000,0.000000000,inf,0
opencl-e2e,Quadro RTX 8000,GPU,w16,p2048,2048,MODEXP,50000,0.000000000,inf,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p2048,2048,EXPONENTIATION,97,0.030840443,3145.221,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p2048,2048,EXPONENTIATION,97,0.002685611,36118.408,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p2048,2048,EXPONENTIATION,97,0.035726394,2715.079,0
opencl-kernel,Quadro RTX 8000,GPU,w8,p2048,2048,EXPONENTIATION,50000,55.167197083,906.336,0
opencl-e2e,Quadro RTX 8000,GPU,w8,p2048,2048,EXPONENTIATION,50000,55.154648437,906.542,0
opencl-kernel,Quadro RTX 8000,GPU,w16,p2048,2048,EXPONENTIATION,50000,10.651931911,4693.984,0
opencl-e2e,Quadro RTX 8000,GPU,w16,p2048,2048,EXPONENTIATION,50000,10.816278889,4622.662,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p2048,2048,DIVIDE,781,0.000056548,13811264.609,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p2048,2048,DIVIDE,781,0.007037898,110970.633,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p2048,2048,DIVIDE,781,0.003754272,208029.685,0
library,Quadro RTX 8000,gpu,cgbn,p2048,2048,DIVIDE,50000,0.000098304,508626302.083,0
opencl-kernel,Quadro RTX 8000,GPU,w8,p2048,2048,DIVIDE,50000,2.626582508,19036.143,0
opencl-e2e,Quadro RTX 8000,GPU,w8,p2048,2048,DIVIDE,50000,2.632506732,18993.304,0
opencl-kernel,Quadro RTX 8000,GPU,w16,p2048,2048,DIVIDE,50000,0.545283400,91695.438,0
opencl-e2e,Quadro RTX 8000,GPU,w16,p2048,2048,DIVIDE,50000,0.549332343,91019.582,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p2048,2048,ISQRT,195,0.000108746,1793167.536,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p2048,2048,ISQRT,195,0.000014247,13687181.283,0
opencl-kernel,Quadro RTX 8000,GPU,w8,p2048,2048,ISQRT,50000,34.069101874,1467.605,0
opencl-e2e,Quadro RTX 8000,GPU,w8,p2048,2048,ISQRT,50000,34.077358067,1467.250,0
opencl-kernel,Quadro RTX 8000,GPU,w16,p2048,2048,ISQRT,50000,0.000000000,inf,0
opencl-e2e,Quadro RTX 8000,GPU,w16,p2048,2048,ISQRT,50000,0.000000000,inf,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-1t,p2048,2048,MODMUL_R2,6250,0.011393853,548541.393,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,gmp-nt,p2048,2048,MODMUL_R2,6250,0.000542280,11525410.560,0
library,AMD EPYC 7282 16-Core Processor,host-cpu,openssl-nt,p2048,2048,MODMUL_R2,6250,0.001326831,4710471.589,0
library,Quadro RTX 8000,gpu,cgbn,p2048,2048,MODMUL_R2,50000,0.000399808,125060028.814,0
opencl-kernel,Quadro RTX 8000,GPU,w8,p2048,2048,MODMUL_R2,50000,0.107955074,463155.626,0
opencl-e2e,Quadro RTX 8000,GPU,w8,p2048,2048,MODMUL_R2,50000,0.111962508,446578.063,0
opencl-kernel,Quadro RTX 8000,GPU,w16,p2048,2048,MODMUL_R2,50000,0.002252881,22193802.492,0
opencl-e2e,Quadro RTX 8000,GPU,w16,p2048,2048,MODMUL_R2,50000,0.006039446,8278904.751,0
```
