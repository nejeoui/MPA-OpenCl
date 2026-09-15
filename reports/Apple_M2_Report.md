# MPA-OpenCL benchmark report - Apple M2

> **Note.** The multi-threaded GMP and OpenSSL baseline columns have been
> removed from this report: they predate the 2026-09-12 timing fix and were
> understated (see `reports/README.md`). The single-threaded GMP column, the
> OpenCL-on-CPU rows and all MPA measurements are unaffected and were verified
> against GMP before timing.


## 1. System under test

1 OpenCL device(s) exercised with the identical kernels and operands.

### Device 0 - Apple M2 (GPU)

| Property | Value |
|---|---|
| Model | Apple M2 |
| Type | GPU |
| Vendor | Apple |
| Device memory | 11.84 GiB |
| Max single allocation | 2.22 GiB |
| Local memory | 32 KiB |
| Global cache | 0 KiB |
| Compute units | 10 |
| Max clock | 1000 MHz |
| Max work-group size | 256 |
| OpenCL version | OpenCL 1.2  |
| Driver | 1.2 1.0 |

### Host

| Property | Value |
|---|---|
| CPU | Apple M2 |
| Logical cores | 8 |
| OpenMP threads used | 8 |
| RAM | 16.0 GB |
| OS | Darwin 25.3.0 |
| Kernel | 25.3.0 |
| Arch | arm64 |
| GMP | 6.3.0 |
| OpenSSL | OpenSSL 3.6.4 25 Aug 2026 |
| CGBN | not measured |

## 2. Method

- Base workload 20000 items, scaled down per operator by its cost weight and by modulus size. Device rows honour --min-items (64) so the GPU is not left idle; the CPU libraries keep the smaller count because a full-width MODEXP there costs minutes. Both counts appear in every row as dev/cpu, and throughput is per-second so they remain comparable.
- 5 timed repetitions, **minimum** reported. Two untimed warm-up launches precede them.
- `kernel` times `clEnqueueNDRangeKernel` + `clFinish` only. `e2e` adds the host->device operand writes and the device->host result read.
- Every OpenCL device runs the same kernels on the same operands, so GPU and CPU-OpenCL columns are directly comparable.
- CPU library baselines (GMP, OpenSSL) run those same operands, with temporaries preallocated outside the timed region, so the figure is the arithmetic and not marshalling. The generator is reseeded per modulus and operation so every backend sees identical inputs.
- OpenSSL rows time the nearest BN primitive, which is not always semantically identical (its Montgomery routine expects Montgomery-domain inputs); they measure comparable work, not identical results. Correctness is judged against GMP only.
- Every device cell is checked word-for-word against GMP before it is timed. A cell that mismatches is reported and excluded from the speedup tables.
- Total wall time 1259.4 s.

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

### Device 0 - Apple M2 (GPU)

#### secp256k1 (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 20000 / 20000 | 37.11 M | 80.97 M | 98.52 M | 80.97 M | 83.68 M | 95.69 M | 90.50 M | 68.49 M | n/a |
| SUBTRACT | 20000 / 20000 | 36.23 M | 75.76 M | 88.50 M | 82.99 M | 82.64 M | 93.02 M | 88.50 M | 137.93 M | n/a |
| ADDMOD | 20000 / 20000 | 36.10 M | 76.63 M | 86.21 M | 96.15 M | 92.17 M | 94.79 M | 90.50 M | 33.28 M | n/a |
| SUBTRACTMOD | 20000 / 20000 | 36.10 M | 74.35 M | 95.24 M | 87.34 M | 86.58 M | 92.59 M | 96.15 M | 38.24 M | n/a |
| MULTIPLYOPERANDSCANNING | 20000 / 20000 | 3.95 M | 14.97 M | 40.49 M | 86.96 M | 82.99 M | 90.50 M | 91.74 M | 67.11 M | n/a |
| MULTIPLYPRODUCTSCANNING | 20000 / 20000 | 7.88 M | 24.54 M | 54.50 M | 45.25 M | 53.62 M | 78.12 M | 71.43 M | 66.89 M | n/a |
| MONTGOMERYMULTIPLICATION | 20000 / 20000 | 43.76 M | 75.76 M | 89.29 M | 88.11 M | 86.58 M | 91.32 M | 88.50 M | 7.15 M | n/a |
| COMPARE | 20000 / 20000 | 99.01 M | 83.33 M | - | 85.47 M | 87.34 M | 103.09 M | 91.32 M | 263.16 M | n/a |
| REDUCE | 2500 / 2500 | 8.53 M | 9.77 M | - | 8.33 M | 10.42 M | 10.64 M | 10.33 M | 86.21 M | n/a |
| MODMUL | 1250 / 1250 | 3.44 M | 4.06 M | - | 3.80 M | 4.45 M | 4.73 M | 4.34 M | 13.74 M | n/a |
| MODEXP | 312 / 312 | 38.65 k | 186.05 k | - | 157.42 k | 253.45 k | 157.34 k | 255.53 k | 100.13 k | n/a |
| EXPONENTIATION | 312 / 312 | 5.30 k | 32.99 k | - | 146.27 k | 173.04 k | 172.19 k | 177.68 k | 381.42 k | n/a |
| DIVIDE | 2500 / 2500 | 3.31 M | 4.04 M | - | 6.72 M | 6.81 M | 6.93 M | 6.96 M | 40.32 M | n/a |
| ISQRT | 625 / 625 | 178.47 k | 255.21 k | - | 463.99 k | 479.29 k | 464.34 k | 478.56 k | 17.36 M | n/a |
| MODMUL_R2 | 20000 / 20000 | 33.06 M | 74.07 M | - | 73.80 M | 80.00 M | 80.32 M | 82.30 M | 14.00 M | n/a |

#### rsa256(composite) (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 20000 / 20000 | 66.01 M | 76.34 M | 89.29 M | 87.72 M | 82.64 M | 100.00 M | 94.34 M | 70.42 M | n/a |
| SUBTRACT | 20000 / 20000 | 69.20 M | 75.76 M | 95.69 M | 86.58 M | 88.11 M | 92.59 M | 89.29 M | 137.93 M | n/a |
| ADDMOD | 20000 / 20000 | 65.79 M | 72.20 M | 96.62 M | 69.93 M | 86.58 M | 92.17 M | 101.01 M | 38.39 M | n/a |
| SUBTRACTMOD | 20000 / 20000 | 66.01 M | 75.47 M | 90.50 M | 70.67 M | 86.58 M | 100.50 M | 89.69 M | 36.17 M | n/a |
| MULTIPLYOPERANDSCANNING | 20000 / 20000 | 3.94 M | 14.98 M | 40.65 M | 67.57 M | 84.03 M | 94.34 M | 88.89 M | 67.11 M | n/a |
| MULTIPLYPRODUCTSCANNING | 20000 / 20000 | 7.91 M | 25.13 M | 54.95 M | 54.79 M | 52.36 M | 74.07 M | 74.35 M | 67.11 M | n/a |
| MONTGOMERYMULTIPLICATION | 20000 / 20000 | 43.67 M | 74.63 M | 90.91 M | 62.70 M | 83.68 M | 85.11 M | 87.34 M | 7.08 M | n/a |
| COMPARE | 20000 / 20000 | 95.69 M | 80.00 M | - | 74.63 M | 85.84 M | 93.90 M | 88.50 M | 263.16 M | n/a |
| REDUCE | 2500 / 2500 | 8.42 M | 9.36 M | - | 6.96 M | 9.96 M | 9.77 M | 8.90 M | 53.19 M | n/a |
| MODMUL | 1250 / 1250 | 3.33 M | 4.14 M | - | 2.66 M | 4.19 M | 3.91 M | 3.97 M | 14.20 M | n/a |
| MODEXP | 312 / 312 | 46.52 k | 258.28 k | - | 196.47 k | 311.69 k | 193.55 k | 315.15 k | 107.03 k | n/a |
| EXPONENTIATION | 312 / 312 | 5.30 k | 32.98 k | - | 145.59 k | 174.40 k | 172.85 k | 178.29 k | 381.42 k | n/a |
| DIVIDE | 2500 / 2500 | 3.34 M | 4.19 M | - | 6.85 M | 6.83 M | 7.02 M | 7.10 M | 39.06 M | n/a |
| ISQRT | 625 / 625 | 154.89 k | 222.26 k | - | 411.45 k | 441.07 k | 415.28 k | 446.11 k | 17.36 M | n/a |
| MODMUL_R2 | 20000 / 20000 | 33.33 M | 72.73 M | - | 79.05 M | 76.92 M | 77.82 M | 86.21 M | 14.09 M | n/a |

#### brainpoolP512r1 (512-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 10000 / 10000 | 24.69 M | 31.95 M | 36.50 M | 38.46 M | 40.32 M | 45.66 M | 48.54 M | 63.29 M | n/a |
| SUBTRACT | 10000 / 10000 | 24.94 M | 35.34 M | 35.34 M | 38.91 M | 40.49 M | 46.51 M | 43.29 M | 120.48 M | n/a |
| ADDMOD | 10000 / 10000 | 24.15 M | 31.15 M | 30.21 M | 44.64 M | 44.64 M | 47.85 M | 45.87 M | 37.04 M | n/a |
| SUBTRACTMOD | 10000 / 10000 | 24.45 M | 31.55 M | 27.93 M | 45.05 M | 41.84 M | 46.08 M | 44.64 M | 36.63 M | n/a |
| MULTIPLYOPERANDSCANNING | 10000 / 10000 | 690.37 k | 2.69 M | 8.97 M | 28.74 M | 22.42 M | 24.94 M | 22.73 M | 30.40 M | n/a |
| MULTIPLYPRODUCTSCANNING | 10000 / 10000 | 928.68 k | 3.96 M | 12.22 M | 9.12 M | 8.32 M | 29.24 M | 26.60 M | 30.40 M | n/a |
| MONTGOMERYMULTIPLICATION | 10000 / 10000 | 1.41 M | 33.00 M | 40.00 M | 36.63 M | 36.76 M | 25.77 M | 21.83 M | 3.16 M | n/a |
| COMPARE | 10000 / 10000 | 44.05 M | 43.48 M | - | 41.32 M | 40.98 M | 40.82 M | 39.68 M | 270.27 M | n/a |
| REDUCE | 1250 / 1250 | 2.91 M | 3.63 M | - | 2.84 M | 2.27 M | 2.45 M | 4.01 M | 54.35 M | n/a |
| MODMUL | 625 / 625 | 747.61 k | 1.26 M | - | 1.10 M | 1.42 M | 699.11 k | 1.42 M | 6.79 M | n/a |
| MODEXP | 156 / 156 | 1.14 k | 21.08 k | - | 22.03 k | 11.16 k | 21.85 k | 11.65 k | 26.32 k | n/a |
| EXPONENTIATION | 156 / 156 | 330.0 | 1.71 k | - | 6.42 k | 7.29 k | 6.38 k | 7.24 k | 126.42 k | n/a |
| DIVIDE | 1250 / 1250 | 310.71 k | 1.10 M | - | 1.25 M | 1.22 M | 1.26 M | 1.23 M | 37.88 M | n/a |
| ISQRT | 312 / 312 | 11.55 k | 47.37 k | - | 58.70 k | 57.47 k | 58.99 k | 57.45 k | 9.45 M | n/a |
| MODMUL_R2 | 10000 / 10000 | 6.35 M | 31.85 M | - | 33.11 M | 34.97 M | 34.36 M | 29.24 M | 6.85 M | n/a |

#### p1024 (1024-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 5000 / 5000 | 12.17 M | 15.63 M | 17.48 M | 19.08 M | 19.16 M | 24.04 M | 22.52 M | 54.95 M | n/a |
| SUBTRACT | 5000 / 5000 | 12.41 M | 15.87 M | 15.11 M | 19.31 M | 19.46 M | 23.92 M | 24.27 M | 100.00 M | n/a |
| ADDMOD | 5000 / 5000 | 10.78 M | 15.97 M | 12.41 M | 24.04 M | 21.65 M | 23.70 M | 21.93 M | 24.27 M | n/a |
| SUBTRACTMOD | 5000 / 5000 | 10.78 M | 16.34 M | 12.89 M | 23.92 M | 22.83 M | 25.77 M | 22.03 M | 29.07 M | n/a |
| MULTIPLYOPERANDSCANNING | 5000 / 5000 | 174.28 k | 691.75 k | 2.56 M | 13.19 M | 12.76 M | 11.93 M | 11.96 M | 8.67 M | n/a |
| MULTIPLYPRODUCTSCANNING | 5000 / 5000 | 289.91 k | 1.10 M | 3.82 M | 3.83 M | 8.04 M | 8.18 M | 9.29 M | 8.68 M | n/a |
| MONTGOMERYMULTIPLICATION | 5000 / 5000 | 581.40 k | 1.51 M | 13.66 M | 16.78 M | 18.18 M | 16.29 M | 17.30 M | 1.04 M | n/a |
| COMPARE | 5000 / 5000 | 22.83 M | 22.94 M | - | 25.38 M | 23.26 M | 26.04 M | 27.17 M | 277.78 M | n/a |
| REDUCE | 625 / 625 | 331.39 k | 1.07 M | - | 1.50 M | 1.41 M | 1.43 M | 1.41 M | 62.50 M | n/a |
| MODMUL | 312 / 312 | 34.97 k | 91.58 k | - | 385.66 k | 359.45 k | 401.54 k | 369.67 k | 2.52 M | n/a |
| MODEXP | 78 / 78 | 15.7 | 189.9 | - | 752.9 | 1.73 k | 645.1 | 1.80 k | 4.98 k | n/a |
| EXPONENTIATION | 78 / 78 | 24.4 | 82.7 | - | 403.6 | 345.9 | 396.9 | 359.3 | 39.84 k | n/a |
| DIVIDE | 625 / 625 | 28.38 k | 66.64 k | - | 48.32 k | 84.85 k | 125.63 k | 137.30 k | 32.89 M | n/a |
| ISQRT | 156 / 156 | 1.11 k | 4.74 k | - | 5.54 k | 4.33 k | 5.42 k | 4.34 k | 4.73 M | n/a |
| MODMUL_R2 | 5000 / 5000 | 456.75 k | 2.93 M | - | 10.99 M | 14.84 M | 11.29 M | 14.08 M | 2.46 M | n/a |

#### p2048 (2048-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 2500 / 2500 | 6.22 M | 8.20 M | 6.51 M | 9.58 M | 9.29 M | 11.57 M | 11.21 M | 41.67 M | n/a |
| SUBTRACT | 2500 / 2500 | 6.11 M | 8.25 M | 6.76 M | 10.78 M | 10.42 M | 11.26 M | 10.73 M | 69.44 M | n/a |
| ADDMOD | 2500 / 2500 | 5.36 M | 6.85 M | 8.74 M | 11.01 M | 10.68 M | 8.36 M | 7.79 M | 21.19 M | n/a |
| SUBTRACTMOD | 2500 / 2500 | 6.04 M | 7.72 M | 5.81 M | 11.06 M | 10.20 M | 9.26 M | 8.20 M | 25.77 M | n/a |
| MULTIPLYOPERANDSCANNING | 2500 / 2500 | 45.07 k | 178.62 k | 677.32 k | 4.01 M | 3.65 M | 3.40 M | 3.63 M | 2.69 M | n/a |
| MULTIPLYPRODUCTSCANNING | 2500 / 2500 | 73.32 k | 286.76 k | 1.06 M | 1.06 M | 3.14 M | 2.54 M | 2.52 M | 2.71 M | n/a |
| MONTGOMERYMULTIPLICATION | 2500 / 2500 | 55.34 k | 329.08 k | 1.01 M | 1.54 M | 1.39 M | 1.60 M | 1.63 M | 368.84 k | n/a |
| COMPARE | 2500 / 2500 | 10.64 M | 11.42 M | - | 11.57 M | 11.11 M | 11.90 M | 11.63 M | 277.78 M | n/a |
| REDUCE | 312 / 312 | 22.46 k | 45.84 k | - | 70.96 k | 70.68 k | 70.84 k | 56.04 k | 44.57 M | n/a |
| MODMUL | 156 / 156 | 2.52 k | 3.87 k | - | 6.18 k | 6.47 k | 5.67 k | 6.40 k | 896.55 k | n/a |
| MODEXP | 64 / 64 | 3.1 | 5.3 | - | 31.8 | 32.9 | 32.0 | 32.7 | 718.7 | n/a |
| EXPONENTIATION | 64 / 64 | 2.5 | 9.2 | - | 34.4 | 34.8 | 34.4 | 34.6 | 6.40 k | n/a |
| DIVIDE | 312 / 312 | 2.94 k | 5.10 k | - | 9.20 k | 9.07 k | 9.35 k | 9.15 k | 28.36 M | n/a |
| ISQRT | 78 / 78 | 96.2 | 172.3 | - | 269.8 | 325.6 | 287.5 | 319.2 | 3.00 M | n/a |
| MODMUL_R2 | 2500 / 2500 | 154.88 k | 309.37 k | - | 1.01 M | 1.19 M | 1.31 M | 1.31 M | 821.29 k | n/a |

## 5. Head to head

Best OpenCL GPU result against best OpenCL CPU result and the CPU libraries.
Ratios above 1.00x mean the GPU is faster than that baseline.

### secp256k1 (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w32 | 98.52 M | none | n/a | 68.49 M | n/a | n/a | n/a |
| SUBTRACT | w32-il | 93.02 M | none | n/a | 137.93 M | n/a | n/a | n/a |
| ADDMOD | w32-opt | 96.15 M | none | n/a | 33.28 M | n/a | n/a | n/a |
| SUBTRACTMOD | w32-il64 | 96.15 M | none | n/a | 38.24 M | n/a | n/a | n/a |
| MULTIPLYOPERANDSCANNING | w32-il64 | 91.74 M | none | n/a | 67.11 M | n/a | n/a | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il | 78.12 M | none | n/a | 66.89 M | n/a | n/a | n/a |
| MONTGOMERYMULTIPLICATION | w32-il | 91.32 M | none | n/a | 7.15 M | n/a | n/a | n/a |
| COMPARE | w32-il | 103.09 M | none | n/a | 263.16 M | n/a | n/a | n/a |
| REDUCE | w32-il | 10.64 M | none | n/a | 86.21 M | n/a | n/a | n/a |
| MODMUL | w32-il | 4.73 M | none | n/a | 13.74 M | n/a | n/a | n/a |
| MODEXP | w32-il64 | 255.53 k | none | n/a | 100.13 k | n/a | n/a | n/a |
| EXPONENTIATION | w32-il64 | 177.68 k | none | n/a | 381.42 k | n/a | n/a | n/a |
| DIVIDE | w32-il64 | 6.96 M | none | n/a | 40.32 M | n/a | n/a | n/a |
| ISQRT | w32-o64 | 479.29 k | none | n/a | 17.36 M | n/a | n/a | n/a |
| MODMUL_R2 | w32-il64 | 82.30 M | none | n/a | 14.00 M | n/a | n/a | n/a |

### rsa256(composite) (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w32-il | 100.00 M | none | n/a | 70.42 M | n/a | n/a | n/a |
| SUBTRACT | w32 | 95.69 M | none | n/a | 137.93 M | n/a | n/a | n/a |
| ADDMOD | w32-il64 | 101.01 M | none | n/a | 38.39 M | n/a | n/a | n/a |
| SUBTRACTMOD | w32-il | 100.50 M | none | n/a | 36.17 M | n/a | n/a | n/a |
| MULTIPLYOPERANDSCANNING | w32-il | 94.34 M | none | n/a | 67.11 M | n/a | n/a | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 74.35 M | none | n/a | 67.11 M | n/a | n/a | n/a |
| MONTGOMERYMULTIPLICATION | w32 | 90.91 M | none | n/a | 7.08 M | n/a | n/a | n/a |
| COMPARE | w8 | 95.69 M | none | n/a | 263.16 M | n/a | n/a | n/a |
| REDUCE | w32-o64 | 9.96 M | none | n/a | 53.19 M | n/a | n/a | n/a |
| MODMUL | w32-o64 | 4.19 M | none | n/a | 14.20 M | n/a | n/a | n/a |
| MODEXP | w32-il64 | 315.15 k | none | n/a | 107.03 k | n/a | n/a | n/a |
| EXPONENTIATION | w32-il64 | 178.29 k | none | n/a | 381.42 k | n/a | n/a | n/a |
| DIVIDE | w32-il64 | 7.10 M | none | n/a | 39.06 M | n/a | n/a | n/a |
| ISQRT | w32-il64 | 446.11 k | none | n/a | 17.36 M | n/a | n/a | n/a |
| MODMUL_R2 | w32-il64 | 86.21 M | none | n/a | 14.09 M | n/a | n/a | n/a |

### brainpoolP512r1 (512-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w32-il64 | 48.54 M | none | n/a | 63.29 M | n/a | n/a | n/a |
| SUBTRACT | w32-il | 46.51 M | none | n/a | 120.48 M | n/a | n/a | n/a |
| ADDMOD | w32-il | 47.85 M | none | n/a | 37.04 M | n/a | n/a | n/a |
| SUBTRACTMOD | w32-il | 46.08 M | none | n/a | 36.63 M | n/a | n/a | n/a |
| MULTIPLYOPERANDSCANNING | w32-opt | 28.74 M | none | n/a | 30.40 M | n/a | n/a | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il | 29.24 M | none | n/a | 30.40 M | n/a | n/a | n/a |
| MONTGOMERYMULTIPLICATION | w32 | 40.00 M | none | n/a | 3.16 M | n/a | n/a | n/a |
| COMPARE | w8 | 44.05 M | none | n/a | 270.27 M | n/a | n/a | n/a |
| REDUCE | w32-il64 | 4.01 M | none | n/a | 54.35 M | n/a | n/a | n/a |
| MODMUL | w32-il64 | 1.42 M | none | n/a | 6.79 M | n/a | n/a | n/a |
| MODEXP | w32-opt | 22.03 k | none | n/a | 26.32 k | n/a | n/a | n/a |
| EXPONENTIATION | w32-o64 | 7.29 k | none | n/a | 126.42 k | n/a | n/a | n/a |
| DIVIDE | w32-il | 1.26 M | none | n/a | 37.88 M | n/a | n/a | n/a |
| ISQRT | w32-il | 58.99 k | none | n/a | 9.45 M | n/a | n/a | n/a |
| MODMUL_R2 | w32-o64 | 34.97 M | none | n/a | 6.85 M | n/a | n/a | n/a |

### p1024 (1024-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w32-il | 24.04 M | none | n/a | 54.95 M | n/a | n/a | n/a |
| SUBTRACT | w32-il64 | 24.27 M | none | n/a | 100.00 M | n/a | n/a | n/a |
| ADDMOD | w32-opt | 24.04 M | none | n/a | 24.27 M | n/a | n/a | n/a |
| SUBTRACTMOD | w32-il | 25.77 M | none | n/a | 29.07 M | n/a | n/a | n/a |
| MULTIPLYOPERANDSCANNING | w32-opt | 13.19 M | none | n/a | 8.67 M | n/a | n/a | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 9.29 M | none | n/a | 8.68 M | n/a | n/a | n/a |
| MONTGOMERYMULTIPLICATION | w32-o64 | 18.18 M | none | n/a | 1.04 M | n/a | n/a | n/a |
| COMPARE | w32-il64 | 27.17 M | none | n/a | 277.78 M | n/a | n/a | n/a |
| REDUCE | w32-opt | 1.50 M | none | n/a | 62.50 M | n/a | n/a | n/a |
| MODMUL | w32-il | 401.54 k | none | n/a | 2.52 M | n/a | n/a | n/a |
| MODEXP | w32-il64 | 1.80 k | none | n/a | 4.98 k | n/a | n/a | n/a |
| EXPONENTIATION | w32-opt | 403.6 | none | n/a | 39.84 k | n/a | n/a | n/a |
| DIVIDE | w32-il64 | 137.30 k | none | n/a | 32.89 M | n/a | n/a | n/a |
| ISQRT | w32-opt | 5.54 k | none | n/a | 4.73 M | n/a | n/a | n/a |
| MODMUL_R2 | w32-o64 | 14.84 M | none | n/a | 2.46 M | n/a | n/a | n/a |

### p2048 (2048-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w32-il | 11.57 M | none | n/a | 41.67 M | n/a | n/a | n/a |
| SUBTRACT | w32-il | 11.26 M | none | n/a | 69.44 M | n/a | n/a | n/a |
| ADDMOD | w32-opt | 11.01 M | none | n/a | 21.19 M | n/a | n/a | n/a |
| SUBTRACTMOD | w32-opt | 11.06 M | none | n/a | 25.77 M | n/a | n/a | n/a |
| MULTIPLYOPERANDSCANNING | w32-opt | 4.01 M | none | n/a | 2.69 M | n/a | n/a | n/a |
| MULTIPLYPRODUCTSCANNING | w32-o64 | 3.14 M | none | n/a | 2.71 M | n/a | n/a | n/a |
| MONTGOMERYMULTIPLICATION | w32-il64 | 1.63 M | none | n/a | 368.84 k | n/a | n/a | n/a |
| COMPARE | w32-il | 11.90 M | none | n/a | 277.78 M | n/a | n/a | n/a |
| REDUCE | w32-opt | 70.96 k | none | n/a | 44.57 M | n/a | n/a | n/a |
| MODMUL | w32-o64 | 6.47 k | none | n/a | 896.55 k | n/a | n/a | n/a |
| MODEXP | w32-o64 | 32.9 | none | n/a | 718.7 | n/a | n/a | n/a |
| EXPONENTIATION | w32-o64 | 34.8 | none | n/a | 6.40 k | n/a | n/a | n/a |
| DIVIDE | w32-il | 9.35 k | none | n/a | 28.36 M | n/a | n/a | n/a |
| ISQRT | w32-o64 | 325.6 | none | n/a | 3.00 M | n/a | n/a | n/a |
| MODMUL_R2 | w32-il64 | 1.31 M | none | n/a | 821.29 k | n/a | n/a | n/a |

## 6. CGBN

Not measured on this run. CGBN is CUDA-only and is not built into this host.
To populate the CGBN columns, produce `cgbn_results.tsv` next to the binary with one
whitespace-separated row per measurement and re-run:

```
# modulus_name  operation_name  items  seconds
secp256k1  MODMUL  20000  0.00123
```

`modulus_name` and `operation_name` must match the spellings used in the tables above.

## 7. Raw data

Also written to `Apple_M2_Report.csv` for analysis.

```csv
kind,device,device_type,kernel,modulus,bits,operation,items,seconds,ops_per_sec,mismatches
library,Apple M2,host-cpu,gmp-1t,secp256k1,256,ADD,20000,0.000292000,68493161.631,0
library,Apple M2,host-cpu,gmp-nt,secp256k1,256,ADD,20000,0.000160000,125000003.158,0
library,Apple M2,host-cpu,openssl-nt,secp256k1,256,ADD,20000,0.000150000,133333352.870,0
opencl-kernel,Apple M2,GPU,w8,secp256k1,256,ADD,20000,0.000539000,37105753.356,0
opencl-e2e,Apple M2,GPU,w8,secp256k1,256,ADD,20000,0.000634000,31545741.506,0
opencl-kernel,Apple M2,GPU,w16,secp256k1,256,ADD,20000,0.000247000,80971653.974,0
opencl-e2e,Apple M2,GPU,w16,secp256k1,256,ADD,20000,0.000363000,55096423.272,0
opencl-kernel,Apple M2,GPU,w32,secp256k1,256,ADD,20000,0.000203000,98522180.041,0
opencl-e2e,Apple M2,GPU,w32,secp256k1,256,ADD,20000,0.000288000,69444449.006,0
opencl-kernel,Apple M2,GPU,w32-opt,secp256k1,256,ADD,20000,0.000247000,80971653.974,0
opencl-e2e,Apple M2,GPU,w32-opt,secp256k1,256,ADD,20000,0.000326000,61349697.129,0
opencl-kernel,Apple M2,GPU,w32-o64,secp256k1,256,ADD,20000,0.000239000,83682023.347,0
opencl-e2e,Apple M2,GPU,w32-o64,secp256k1,256,ADD,20000,0.000335000,59701494.370,0
opencl-kernel,Apple M2,GPU,w32-il,secp256k1,256,ADD,20000,0.000209000,95693797.480,0
opencl-e2e,Apple M2,GPU,w32-il,secp256k1,256,ADD,20000,0.000281000,71174380.773,0
opencl-kernel,Apple M2,GPU,w32-il64,secp256k1,256,ADD,20000,0.000221000,90497739.694,0
opencl-e2e,Apple M2,GPU,w32-il64,secp256k1,256,ADD,20000,0.000309000,64724913.947,0
library,Apple M2,host-cpu,gmp-1t,secp256k1,256,SUBTRACT,20000,0.000145000,137931036.237,0
library,Apple M2,host-cpu,gmp-nt,secp256k1,256,SUBTRACT,20000,0.000093000,215053741.112,0
library,Apple M2,host-cpu,openssl-nt,secp256k1,256,SUBTRACT,20000,0.000099000,202020207.866,0
opencl-kernel,Apple M2,GPU,w8,secp256k1,256,SUBTRACT,20000,0.000552000,36231883.254,0
opencl-e2e,Apple M2,GPU,w8,secp256k1,256,SUBTRACT,20000,0.000654000,30581040.970,0
opencl-kernel,Apple M2,GPU,w16,secp256k1,256,SUBTRACT,20000,0.000264000,75757566.814,0
opencl-e2e,Apple M2,GPU,w16,secp256k1,256,SUBTRACT,20000,0.000362000,55248615.627,0
opencl-kernel,Apple M2,GPU,w32,secp256k1,256,SUBTRACT,20000,0.000226000,88495585.149,0
opencl-e2e,Apple M2,GPU,w32,secp256k1,256,SUBTRACT,20000,0.000312000,64102566.320,0
opencl-kernel,Apple M2,GPU,w32-opt,secp256k1,256,SUBTRACT,20000,0.000241000,82987561.355,0
opencl-e2e,Apple M2,GPU,w32-opt,secp256k1,256,SUBTRACT,20000,0.000317000,63091488.805,0
opencl-kernel,Apple M2,GPU,w32-o64,secp256k1,256,SUBTRACT,20000,0.000242000,82644615.030,0
opencl-e2e,Apple M2,GPU,w32-o64,secp256k1,256,SUBTRACT,20000,0.000326000,61349686.175,0
opencl-kernel,Apple M2,GPU,w32-il,secp256k1,256,SUBTRACT,20000,0.000215000,93023252.655,0
opencl-e2e,Apple M2,GPU,w32-il,secp256k1,256,SUBTRACT,20000,0.000290000,68965518.118,0
opencl-kernel,Apple M2,GPU,w32-il64,secp256k1,256,SUBTRACT,20000,0.000226000,88495585.149,0
opencl-e2e,Apple M2,GPU,w32-il64,secp256k1,256,SUBTRACT,20000,0.000313000,63897764.227,0
library,Apple M2,host-cpu,gmp-1t,secp256k1,256,ADDMOD,20000,0.000601000,33277871.440,0
library,Apple M2,host-cpu,gmp-nt,secp256k1,256,ADDMOD,20000,0.000203000,98522180.041,0
library,Apple M2,host-cpu,openssl-nt,secp256k1,256,ADDMOD,20000,0.000650000,30769230.685,0
opencl-kernel,Apple M2,GPU,w8,secp256k1,256,ADDMOD,20000,0.000554000,36101081.242,0
opencl-e2e,Apple M2,GPU,w8,secp256k1,256,ADDMOD,20000,0.000667000,29985006.831,0
opencl-kernel,Apple M2,GPU,w16,secp256k1,256,ADDMOD,20000,0.000261000,76628367.137,0
opencl-e2e,Apple M2,GPU,w16,secp256k1,256,ADDMOD,20000,0.000344000,58139534.877,0
opencl-kernel,Apple M2,GPU,w32,secp256k1,256,ADDMOD,20000,0.000232000,86206888.996,0
opencl-e2e,Apple M2,GPU,w32,secp256k1,256,ADDMOD,20000,0.000318000,62893086.372,0
opencl-kernel,Apple M2,GPU,w32-opt,secp256k1,256,ADDMOD,20000,0.000208000,96153840.510,0
opencl-e2e,Apple M2,GPU,w32-opt,secp256k1,256,ADDMOD,20000,0.000289000,69204154.956,0
opencl-kernel,Apple M2,GPU,w32-o64,secp256k1,256,ADDMOD,20000,0.000217000,92165889.048,0
opencl-e2e,Apple M2,GPU,w32-o64,secp256k1,256,ADDMOD,20000,0.000302000,66225171.863,0
opencl-kernel,Apple M2,GPU,w32-il,secp256k1,256,ADDMOD,20000,0.000211000,94786740.260,0
opencl-e2e,Apple M2,GPU,w32-il,secp256k1,256,ADDMOD,20000,0.000288000,69444449.006,0
opencl-kernel,Apple M2,GPU,w32-il64,secp256k1,256,ADDMOD,20000,0.000221000,90497739.694,0
opencl-e2e,Apple M2,GPU,w32-il64,secp256k1,256,ADDMOD,20000,0.000298000,67114094.180,0
library,Apple M2,host-cpu,gmp-1t,secp256k1,256,SUBTRACTMOD,20000,0.000523000,38240916.008,0
library,Apple M2,host-cpu,gmp-nt,secp256k1,256,SUBTRACTMOD,20000,0.000189000,105820114.808,0
library,Apple M2,host-cpu,openssl-nt,secp256k1,256,SUBTRACTMOD,20000,0.000627000,31897926.571,0
opencl-kernel,Apple M2,GPU,w8,secp256k1,256,SUBTRACTMOD,20000,0.000554000,36101085.035,0
opencl-e2e,Apple M2,GPU,w8,secp256k1,256,SUBTRACTMOD,20000,0.000659000,30349013.100,0
opencl-kernel,Apple M2,GPU,w16,secp256k1,256,SUBTRACTMOD,20000,0.000269000,74349439.330,0
opencl-e2e,Apple M2,GPU,w16,secp256k1,256,SUBTRACTMOD,20000,0.000359000,55710306.741,0
opencl-kernel,Apple M2,GPU,w32,secp256k1,256,SUBTRACTMOD,20000,0.000210000,95238082.795,0
opencl-e2e,Apple M2,GPU,w32,secp256k1,256,SUBTRACTMOD,20000,0.000310000,64516134.448,0
opencl-kernel,Apple M2,GPU,w32-opt,secp256k1,256,SUBTRACTMOD,20000,0.000229000,87336245.499,0
opencl-e2e,Apple M2,GPU,w32-opt,secp256k1,256,SUBTRACTMOD,20000,0.000300000,66666663.500,0
opencl-kernel,Apple M2,GPU,w32-o64,secp256k1,256,SUBTRACTMOD,20000,0.000231000,86580081.813,0
opencl-e2e,Apple M2,GPU,w32-o64,secp256k1,256,SUBTRACTMOD,20000,0.000306000,65359476.755,0
opencl-kernel,Apple M2,GPU,w32-il,secp256k1,256,SUBTRACTMOD,20000,0.000216000,92592586.199,0
opencl-e2e,Apple M2,GPU,w32-il,secp256k1,256,SUBTRACTMOD,20000,0.000290000,68965518.118,0
opencl-kernel,Apple M2,GPU,w32-il64,secp256k1,256,SUBTRACTMOD,20000,0.000208000,96153840.510,0
opencl-e2e,Apple M2,GPU,w32-il64,secp256k1,256,SUBTRACTMOD,20000,0.000309000,64724913.947,0
library,Apple M2,host-cpu,gmp-1t,secp256k1,256,MULTIPLYOPERANDSCANNING,20000,0.000298000,67114094.180,0
library,Apple M2,host-cpu,gmp-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,20000,0.000141000,141843945.570,0
library,Apple M2,host-cpu,openssl-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,20000,0.000179000,111731828.458,0
opencl-kernel,Apple M2,GPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,20000,0.005061000,3951788.181,0
opencl-e2e,Apple M2,GPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,20000,0.005180000,3861003.877,0
opencl-kernel,Apple M2,GPU,w16,secp256k1,256,MULTIPLYOPERANDSCANNING,20000,0.001336000,14970060.030,0
opencl-e2e,Apple M2,GPU,w16,secp256k1,256,MULTIPLYOPERANDSCANNING,20000,0.001453000,13764624.872,0
opencl-kernel,Apple M2,GPU,w32,secp256k1,256,MULTIPLYOPERANDSCANNING,20000,0.000494000,40485831.757,0
opencl-e2e,Apple M2,GPU,w32,secp256k1,256,MULTIPLYOPERANDSCANNING,20000,0.000589000,33955857.216,0
opencl-kernel,Apple M2,GPU,w32-opt,secp256k1,256,MULTIPLYOPERANDSCANNING,20000,0.000230000,86956519.810,0
opencl-e2e,Apple M2,GPU,w32-opt,secp256k1,256,MULTIPLYOPERANDSCANNING,20000,0.000322000,62111799.864,0
opencl-kernel,Apple M2,GPU,w32-o64,secp256k1,256,MULTIPLYOPERANDSCANNING,20000,0.000241000,82987561.355,0
opencl-e2e,Apple M2,GPU,w32-o64,secp256k1,256,MULTIPLYOPERANDSCANNING,20000,0.000329000,60790273.143,0
opencl-kernel,Apple M2,GPU,w32-il,secp256k1,256,MULTIPLYOPERANDSCANNING,20000,0.000221000,90497739.694,0
opencl-e2e,Apple M2,GPU,w32-il,secp256k1,256,MULTIPLYOPERANDSCANNING,20000,0.000316000,63291135.302,0
opencl-kernel,Apple M2,GPU,w32-il64,secp256k1,256,MULTIPLYOPERANDSCANNING,20000,0.000218000,91743131.076,0
opencl-e2e,Apple M2,GPU,w32-il64,secp256k1,256,MULTIPLYOPERANDSCANNING,20000,0.000315000,63492072.795,0
library,Apple M2,host-cpu,gmp-1t,secp256k1,256,MULTIPLYPRODUCTSCANNING,20000,0.000299000,66889630.623,0
library,Apple M2,host-cpu,gmp-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,20000,0.000144000,138888898.012,0
library,Apple M2,host-cpu,openssl-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,20000,0.000185000,108108101.273,0
opencl-kernel,Apple M2,GPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,20000,0.002537000,7883326.786,0
opencl-e2e,Apple M2,GPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,20000,0.002610000,7662835.176,0
opencl-kernel,Apple M2,GPU,w16,secp256k1,256,MULTIPLYPRODUCTSCANNING,20000,0.000815000,24539877.975,0
opencl-e2e,Apple M2,GPU,w16,secp256k1,256,MULTIPLYPRODUCTSCANNING,20000,0.000924000,21645021.817,0
opencl-kernel,Apple M2,GPU,w32,secp256k1,256,MULTIPLYPRODUCTSCANNING,20000,0.000367000,54495912.725,0
opencl-e2e,Apple M2,GPU,w32,secp256k1,256,MULTIPLYPRODUCTSCANNING,20000,0.000455000,43956046.648,0
opencl-kernel,Apple M2,GPU,w32-opt,secp256k1,256,MULTIPLYPRODUCTSCANNING,20000,0.000442000,45248863.888,0
opencl-e2e,Apple M2,GPU,w32-opt,secp256k1,256,MULTIPLYPRODUCTSCANNING,20000,0.000556000,35971220.259,0
opencl-kernel,Apple M2,GPU,w32-o64,secp256k1,256,MULTIPLYPRODUCTSCANNING,20000,0.000373000,53619304.670,0
opencl-e2e,Apple M2,GPU,w32-o64,secp256k1,256,MULTIPLYPRODUCTSCANNING,20000,0.000477000,41928724.248,0
opencl-kernel,Apple M2,GPU,w32-il,secp256k1,256,MULTIPLYPRODUCTSCANNING,20000,0.000256000,78124991.315,0
opencl-e2e,Apple M2,GPU,w32-il,secp256k1,256,MULTIPLYPRODUCTSCANNING,20000,0.000368000,54347833.477,0
opencl-kernel,Apple M2,GPU,w32-il64,secp256k1,256,MULTIPLYPRODUCTSCANNING,20000,0.000280000,71428576.945,0
opencl-e2e,Apple M2,GPU,w32-il64,secp256k1,256,MULTIPLYPRODUCTSCANNING,20000,0.000376000,53191487.823,0
library,Apple M2,host-cpu,gmp-1t,secp256k1,256,MONTGOMERYMULTIPLICATION,20000,0.002796000,7153075.858,0
library,Apple M2,host-cpu,gmp-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,20000,0.000763000,26212320.165,0
library,Apple M2,host-cpu,openssl-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,20000,0.000175000,114285718.361,0
opencl-kernel,Apple M2,GPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,20000,0.000457000,43763677.359,0
opencl-e2e,Apple M2,GPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,20000,0.000547000,36563073.024,0
opencl-kernel,Apple M2,GPU,w16,secp256k1,256,MONTGOMERYMULTIPLICATION,20000,0.000264000,75757583.518,0
opencl-e2e,Apple M2,GPU,w16,secp256k1,256,MONTGOMERYMULTIPLICATION,20000,0.000352000,56818181.374,0
opencl-kernel,Apple M2,GPU,w32,secp256k1,256,MONTGOMERYMULTIPLICATION,20000,0.000224000,89285707.261,0
opencl-e2e,Apple M2,GPU,w32,secp256k1,256,MONTGOMERYMULTIPLICATION,20000,0.000302000,66225171.863,0
opencl-kernel,Apple M2,GPU,w32-opt,secp256k1,256,MONTGOMERYMULTIPLICATION,20000,0.000227000,88105733.758,0
opencl-e2e,Apple M2,GPU,w32-opt,secp256k1,256,MONTGOMERYMULTIPLICATION,20000,0.000300000,66666663.500,0
opencl-kernel,Apple M2,GPU,w32-o64,secp256k1,256,MONTGOMERYMULTIPLICATION,20000,0.000231000,86580103.630,0
opencl-e2e,Apple M2,GPU,w32-o64,secp256k1,256,MONTGOMERYMULTIPLICATION,20000,0.000310000,64516134.448,0
opencl-kernel,Apple M2,GPU,w32-il,secp256k1,256,MONTGOMERYMULTIPLICATION,20000,0.000219000,91324209.440,0
opencl-e2e,Apple M2,GPU,w32-il,secp256k1,256,MONTGOMERYMULTIPLICATION,20000,0.000302000,66225171.863,0
opencl-kernel,Apple M2,GPU,w32-il64,secp256k1,256,MONTGOMERYMULTIPLICATION,20000,0.000226000,88495585.149,0
opencl-e2e,Apple M2,GPU,w32-il64,secp256k1,256,MONTGOMERYMULTIPLICATION,20000,0.000322000,62111799.864,0
library,Apple M2,host-cpu,gmp-1t,secp256k1,256,COMPARE,20000,0.000076000,263157906.424,0
library,Apple M2,host-cpu,gmp-nt,secp256k1,256,COMPARE,20000,0.000063000,317460246.653,0
library,Apple M2,host-cpu,openssl-nt,secp256k1,256,COMPARE,20000,0.000121000,165289269.816,0
opencl-kernel,Apple M2,GPU,w8,secp256k1,256,COMPARE,20000,0.000202000,99009888.869,0
opencl-e2e,Apple M2,GPU,w8,secp256k1,256,COMPARE,20000,0.000282000,70921987.424,0
opencl-kernel,Apple M2,GPU,w16,secp256k1,256,COMPARE,20000,0.000240000,83333325.333,0
opencl-e2e,Apple M2,GPU,w16,secp256k1,256,COMPARE,20000,0.000329000,60790273.143,0
opencl-kernel,Apple M2,GPU,w32-opt,secp256k1,256,COMPARE,20000,0.000234000,85470093.742,0
opencl-e2e,Apple M2,GPU,w32-opt,secp256k1,256,COMPARE,20000,0.000329000,60790273.143,0
opencl-kernel,Apple M2,GPU,w32-o64,secp256k1,256,COMPARE,20000,0.000229000,87336245.499,0
opencl-e2e,Apple M2,GPU,w32-o64,secp256k1,256,COMPARE,20000,0.000319000,62695916.410,0
opencl-kernel,Apple M2,GPU,w32-il,secp256k1,256,COMPARE,20000,0.000194000,103092771.803,0
opencl-e2e,Apple M2,GPU,w32-il,secp256k1,256,COMPARE,20000,0.000274000,72992703.155,0
opencl-kernel,Apple M2,GPU,w32-il64,secp256k1,256,COMPARE,20000,0.000219000,91324209.440,0
opencl-e2e,Apple M2,GPU,w32-il64,secp256k1,256,COMPARE,20000,0.000297000,67340069.289,0
library,Apple M2,host-cpu,gmp-1t,secp256k1,256,REDUCE,2500,0.000029000,86206932.254,0
library,Apple M2,host-cpu,gmp-nt,secp256k1,256,REDUCE,2500,0.000063000,39682567.496,0
library,Apple M2,host-cpu,openssl-nt,secp256k1,256,REDUCE,2500,0.000131000,19083968.305,0
opencl-kernel,Apple M2,GPU,w8,secp256k1,256,REDUCE,2500,0.000293000,8532422.650,0
opencl-e2e,Apple M2,GPU,w8,secp256k1,256,REDUCE,2500,0.000320000,7812500.197,0
opencl-kernel,Apple M2,GPU,w16,secp256k1,256,REDUCE,2500,0.000256000,9765626.135,0
opencl-e2e,Apple M2,GPU,w16,secp256k1,256,REDUCE,2500,0.000280000,8928572.118,0
opencl-kernel,Apple M2,GPU,w32-opt,secp256k1,256,REDUCE,2500,0.000300000,8333332.938,0
opencl-e2e,Apple M2,GPU,w32-opt,secp256k1,256,REDUCE,2500,0.000336000,7440476.250,0
opencl-kernel,Apple M2,GPU,w32-o64,secp256k1,256,REDUCE,2500,0.000240000,10416668.193,0
opencl-e2e,Apple M2,GPU,w32-o64,secp256k1,256,REDUCE,2500,0.000276000,9057970.813,0
opencl-kernel,Apple M2,GPU,w32-il,secp256k1,256,REDUCE,2500,0.000235000,10638298.553,0
opencl-e2e,Apple M2,GPU,w32-il,secp256k1,256,REDUCE,2500,0.000265000,9433960.884,0
opencl-kernel,Apple M2,GPU,w32-il64,secp256k1,256,REDUCE,2500,0.000242000,10330579.364,0
opencl-e2e,Apple M2,GPU,w32-il64,secp256k1,256,REDUCE,2500,0.000270000,9259260.616,0
library,Apple M2,host-cpu,gmp-1t,secp256k1,256,MODMUL,1250,0.000091000,13736255.791,0
library,Apple M2,host-cpu,gmp-nt,secp256k1,256,MODMUL,1250,0.000082000,15243905.394,0
library,Apple M2,host-cpu,openssl-nt,secp256k1,256,MODMUL,1250,0.000147000,8503400.586,0
opencl-kernel,Apple M2,GPU,w8,secp256k1,256,MODMUL,1250,0.000363000,3443525.902,0
opencl-e2e,Apple M2,GPU,w8,secp256k1,256,MODMUL,1250,0.000389000,3213367.663,0
opencl-kernel,Apple M2,GPU,w16,secp256k1,256,MODMUL,1250,0.000308000,4058442.102,0
opencl-e2e,Apple M2,GPU,w16,secp256k1,256,MODMUL,1250,0.000333000,3753753.385,0
opencl-kernel,Apple M2,GPU,w32-opt,secp256k1,256,MODMUL,1250,0.000329000,3799392.071,0
opencl-e2e,Apple M2,GPU,w32-opt,secp256k1,256,MODMUL,1250,0.000360000,3472222.731,0
opencl-kernel,Apple M2,GPU,w32-o64,secp256k1,256,MODMUL,1250,0.000281000,4448397.877,0
opencl-e2e,Apple M2,GPU,w32-o64,secp256k1,256,MODMUL,1250,0.000310000,4032258.403,0
opencl-kernel,Apple M2,GPU,w32-il,secp256k1,256,MODMUL,1250,0.000264000,4734847.926,0
opencl-e2e,Apple M2,GPU,w32-il,secp256k1,256,MODMUL,1250,0.000302000,4139072.444,0
opencl-kernel,Apple M2,GPU,w32-il64,secp256k1,256,MODMUL,1250,0.000288000,4340278.063,0
opencl-e2e,Apple M2,GPU,w32-il64,secp256k1,256,MODMUL,1250,0.000312000,4006410.395,0
library,Apple M2,host-cpu,gmp-1t,secp256k1,256,MODEXP,312,0.003116000,100128.370,0
library,Apple M2,host-cpu,gmp-nt,secp256k1,256,MODEXP,312,0.000771000,404669.265,0
library,Apple M2,host-cpu,openssl-nt,secp256k1,256,MODEXP,312,0.000751000,415446.064,0
opencl-kernel,Apple M2,GPU,w8,secp256k1,256,MODEXP,312,0.008073000,38647.343,0
opencl-e2e,Apple M2,GPU,w8,secp256k1,256,MODEXP,312,0.008111000,38466.280,0
opencl-kernel,Apple M2,GPU,w16,secp256k1,256,MODEXP,312,0.001677000,186046.509,0
opencl-e2e,Apple M2,GPU,w16,secp256k1,256,MODEXP,312,0.001678000,185935.634,0
opencl-kernel,Apple M2,GPU,w32-opt,secp256k1,256,MODEXP,312,0.001982000,157416.749,0
opencl-e2e,Apple M2,GPU,w32-opt,secp256k1,256,MODEXP,312,0.001987000,157020.634,0
opencl-kernel,Apple M2,GPU,w32-o64,secp256k1,256,MODEXP,312,0.001231000,253452.477,0
opencl-e2e,Apple M2,GPU,w32-o64,secp256k1,256,MODEXP,312,0.001239000,251815.980,0
opencl-kernel,Apple M2,GPU,w32-il,secp256k1,256,MODEXP,312,0.001983000,157337.366,0
opencl-e2e,Apple M2,GPU,w32-il,secp256k1,256,MODEXP,312,0.002010000,155223.881,0
opencl-kernel,Apple M2,GPU,w32-il64,secp256k1,256,MODEXP,312,0.001221000,255528.247,0
opencl-e2e,Apple M2,GPU,w32-il64,secp256k1,256,MODEXP,312,0.001256000,248407.651,0
library,Apple M2,host-cpu,gmp-1t,secp256k1,256,EXPONENTIATION,312,0.000818000,381418.093,0
library,Apple M2,host-cpu,gmp-nt,secp256k1,256,EXPONENTIATION,312,0.000265000,1177358.577,0
library,Apple M2,host-cpu,openssl-nt,secp256k1,256,EXPONENTIATION,312,0.001925000,162077.918,0
opencl-kernel,Apple M2,GPU,w8,secp256k1,256,EXPONENTIATION,312,0.058876000,5299.273,0
opencl-e2e,Apple M2,GPU,w8,secp256k1,256,EXPONENTIATION,312,0.058903000,5296.844,0
opencl-kernel,Apple M2,GPU,w16,secp256k1,256,EXPONENTIATION,312,0.009458000,32987.947,0
opencl-e2e,Apple M2,GPU,w16,secp256k1,256,EXPONENTIATION,312,0.009481000,32907.921,0
opencl-kernel,Apple M2,GPU,w32-opt,secp256k1,256,EXPONENTIATION,312,0.002133000,146272.855,0
opencl-e2e,Apple M2,GPU,w32-opt,secp256k1,256,EXPONENTIATION,312,0.002156000,144712.430,0
opencl-kernel,Apple M2,GPU,w32-o64,secp256k1,256,EXPONENTIATION,312,0.001803000,173044.926,0
opencl-e2e,Apple M2,GPU,w32-o64,secp256k1,256,EXPONENTIATION,312,0.001830000,170491.801,0
opencl-kernel,Apple M2,GPU,w32-il,secp256k1,256,EXPONENTIATION,312,0.001812000,172185.430,0
opencl-e2e,Apple M2,GPU,w32-il,secp256k1,256,EXPONENTIATION,312,0.001838000,169749.726,0
opencl-kernel,Apple M2,GPU,w32-il64,secp256k1,256,EXPONENTIATION,312,0.001756000,177676.533,0
opencl-e2e,Apple M2,GPU,w32-il64,secp256k1,256,EXPONENTIATION,312,0.001783000,174985.977,0
library,Apple M2,host-cpu,gmp-1t,secp256k1,256,DIVIDE,2500,0.000062000,40322576.459,0
library,Apple M2,host-cpu,gmp-nt,secp256k1,256,DIVIDE,2500,0.000067000,37313440.464,0
library,Apple M2,host-cpu,openssl-nt,secp256k1,256,DIVIDE,2500,0.000122000,20491803.674,0
opencl-kernel,Apple M2,GPU,w8,secp256k1,256,DIVIDE,2500,0.000755000,3311258.338,0
opencl-e2e,Apple M2,GPU,w8,secp256k1,256,DIVIDE,2500,0.000791000,3160556.147,0
opencl-kernel,Apple M2,GPU,w16,secp256k1,256,DIVIDE,2500,0.000619000,4038772.223,0
opencl-e2e,Apple M2,GPU,w16,secp256k1,256,DIVIDE,2500,0.000628000,3980891.663,0
opencl-kernel,Apple M2,GPU,w32-opt,secp256k1,256,DIVIDE,2500,0.000372000,6720430.461,0
opencl-e2e,Apple M2,GPU,w32-opt,secp256k1,256,DIVIDE,2500,0.000400000,6249999.703,0
opencl-kernel,Apple M2,GPU,w32-o64,secp256k1,256,DIVIDE,2500,0.000367000,6811989.091,0
opencl-e2e,Apple M2,GPU,w32-o64,secp256k1,256,DIVIDE,2500,0.000415000,6024096.987,0
opencl-kernel,Apple M2,GPU,w32-il,secp256k1,256,DIVIDE,2500,0.000361000,6925207.505,0
opencl-e2e,Apple M2,GPU,w32-il,secp256k1,256,DIVIDE,2500,0.000396000,6313131.496,0
opencl-kernel,Apple M2,GPU,w32-il64,secp256k1,256,DIVIDE,2500,0.000359000,6963788.343,0
opencl-e2e,Apple M2,GPU,w32-il64,secp256k1,256,DIVIDE,2500,0.000399000,6265663.982,0
library,Apple M2,host-cpu,gmp-1t,secp256k1,256,ISQRT,625,0.000036000,17361119.269,0
library,Apple M2,host-cpu,gmp-nt,secp256k1,256,ISQRT,625,0.000058000,10775855.717,0
opencl-kernel,Apple M2,GPU,w8,secp256k1,256,ISQRT,625,0.003502000,178469.446,0
opencl-e2e,Apple M2,GPU,w8,secp256k1,256,ISQRT,625,0.003507000,178214.996,0
opencl-kernel,Apple M2,GPU,w16,secp256k1,256,ISQRT,625,0.002449000,255206.204,0
opencl-e2e,Apple M2,GPU,w16,secp256k1,256,ISQRT,625,0.002479000,252117.794,0
opencl-kernel,Apple M2,GPU,w32-opt,secp256k1,256,ISQRT,625,0.001347000,463994.057,0
opencl-e2e,Apple M2,GPU,w32-opt,secp256k1,256,ISQRT,625,0.001398000,447067.241,0
opencl-kernel,Apple M2,GPU,w32-o64,secp256k1,256,ISQRT,625,0.001304000,479294.466,0
opencl-e2e,Apple M2,GPU,w32-o64,secp256k1,256,ISQRT,625,0.001329000,470278.408,0
opencl-kernel,Apple M2,GPU,w32-il,secp256k1,256,ISQRT,625,0.001346000,464338.780,0
opencl-e2e,Apple M2,GPU,w32-il,secp256k1,256,ISQRT,625,0.001378000,453555.873,0
opencl-kernel,Apple M2,GPU,w32-il64,secp256k1,256,ISQRT,625,0.001306000,478560.493,0
opencl-e2e,Apple M2,GPU,w32-il64,secp256k1,256,ISQRT,625,0.001328000,470632.536,0
library,Apple M2,host-cpu,gmp-1t,secp256k1,256,MODMUL_R2,20000,0.001429000,13995801.296,0
library,Apple M2,host-cpu,gmp-nt,secp256k1,256,MODMUL_R2,20000,0.000440000,45454542.093,0
library,Apple M2,host-cpu,openssl-nt,secp256k1,256,MODMUL_R2,20000,0.001215000,16460905.396,0
opencl-kernel,Apple M2,GPU,w8,secp256k1,256,MODMUL_R2,20000,0.000605000,33057850.783,0
opencl-e2e,Apple M2,GPU,w8,secp256k1,256,MODMUL_R2,20000,0.000694000,28818444.321,0
opencl-kernel,Apple M2,GPU,w16,secp256k1,256,MODMUL_R2,20000,0.000270000,74074084.928,0
opencl-e2e,Apple M2,GPU,w16,secp256k1,256,MODMUL_R2,20000,0.000363000,55096423.272,0
opencl-kernel,Apple M2,GPU,w32-opt,secp256k1,256,MODMUL_R2,20000,0.000271000,73800746.708,0
opencl-e2e,Apple M2,GPU,w32-opt,secp256k1,256,MODMUL_R2,20000,0.000336000,59523809.996,0
opencl-kernel,Apple M2,GPU,w32-o64,secp256k1,256,MODMUL_R2,20000,0.000250000,80000005.513,0
opencl-e2e,Apple M2,GPU,w32-o64,secp256k1,256,MODMUL_R2,20000,0.000333000,60060054.163,0
opencl-kernel,Apple M2,GPU,w32-il,secp256k1,256,MODMUL_R2,20000,0.000249000,80321293.155,0
opencl-e2e,Apple M2,GPU,w32-il,secp256k1,256,MODMUL_R2,20000,0.000330000,60606058.796,0
opencl-kernel,Apple M2,GPU,w32-il64,secp256k1,256,MODMUL_R2,20000,0.000243000,82304530.923,0
opencl-e2e,Apple M2,GPU,w32-il64,secp256k1,256,MODMUL_R2,20000,0.000330000,60606058.796,0
library,Apple M2,host-cpu,gmp-1t,rsa256(composite),256,ADD,20000,0.000284000,70422533.021,0
library,Apple M2,host-cpu,gmp-nt,rsa256(composite),256,ADD,20000,0.000129000,155038768.416,0
library,Apple M2,host-cpu,openssl-nt,rsa256(composite),256,ADD,20000,0.000141000,141843945.570,0
opencl-kernel,Apple M2,GPU,w8,rsa256(composite),256,ADD,20000,0.000303000,66006592.580,0
opencl-e2e,Apple M2,GPU,w8,rsa256(composite),256,ADD,20000,0.000391000,51150894.006,0
opencl-kernel,Apple M2,GPU,w16,rsa256(composite),256,ADD,20000,0.000262000,76335873.219,0
opencl-e2e,Apple M2,GPU,w16,rsa256(composite),256,ADD,20000,0.000360000,55555554.713,0
opencl-kernel,Apple M2,GPU,w32,rsa256(composite),256,ADD,20000,0.000224000,89285730.462,0
opencl-e2e,Apple M2,GPU,w32,rsa256(composite),256,ADD,20000,0.000299000,66889630.623,0
opencl-kernel,Apple M2,GPU,w32-opt,rsa256(composite),256,ADD,20000,0.000228000,87719302.141,0
opencl-e2e,Apple M2,GPU,w32-opt,rsa256(composite),256,ADD,20000,0.000331000,60422957.536,0
opencl-kernel,Apple M2,GPU,w32-o64,rsa256(composite),256,ADD,20000,0.000242000,82644615.030,0
opencl-e2e,Apple M2,GPU,w32-o64,rsa256(composite),256,ADD,20000,0.000338000,59171605.624,0
opencl-kernel,Apple M2,GPU,w32-il,rsa256(composite),256,ADD,20000,0.000200000,99999995.250,0
opencl-e2e,Apple M2,GPU,w32-il,rsa256(composite),256,ADD,20000,0.000291000,68728521.409,0
opencl-kernel,Apple M2,GPU,w32-il64,rsa256(composite),256,ADD,20000,0.000212000,94339629.558,0
opencl-e2e,Apple M2,GPU,w32-il64,rsa256(composite),256,ADD,20000,0.000306000,65359476.755,0
library,Apple M2,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACT,20000,0.000145000,137931036.237,0
library,Apple M2,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACT,20000,0.000090000,222222254.784,0
library,Apple M2,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACT,20000,0.000130000,153846153.427,0
opencl-kernel,Apple M2,GPU,w8,rsa256(composite),256,SUBTRACT,20000,0.000289000,69204154.956,0
opencl-e2e,Apple M2,GPU,w8,rsa256(composite),256,SUBTRACT,20000,0.000387000,51679589.472,0
opencl-kernel,Apple M2,GPU,w16,rsa256(composite),256,SUBTRACT,20000,0.000264000,75757566.814,0
opencl-e2e,Apple M2,GPU,w16,rsa256(composite),256,SUBTRACT,20000,0.000352000,56818181.374,0
opencl-kernel,Apple M2,GPU,w32,rsa256(composite),256,SUBTRACT,20000,0.000209000,95693797.480,0
opencl-e2e,Apple M2,GPU,w32,rsa256(composite),256,SUBTRACT,20000,0.000307000,65146577.822,0
opencl-kernel,Apple M2,GPU,w32-opt,rsa256(composite),256,SUBTRACT,20000,0.000231000,86580081.813,0
opencl-e2e,Apple M2,GPU,w32-opt,rsa256(composite),256,SUBTRACT,20000,0.000314000,63694266.616,0
opencl-kernel,Apple M2,GPU,w32-o64,rsa256(composite),256,SUBTRACT,20000,0.000227000,88105733.758,0
opencl-e2e,Apple M2,GPU,w32-o64,rsa256(composite),256,SUBTRACT,20000,0.000315000,63492061.063,0
opencl-kernel,Apple M2,GPU,w32-il,rsa256(composite),256,SUBTRACT,20000,0.000216000,92592586.199,0
opencl-e2e,Apple M2,GPU,w32-il,rsa256(composite),256,SUBTRACT,20000,0.000290000,68965518.118,0
opencl-kernel,Apple M2,GPU,w32-il64,rsa256(composite),256,SUBTRACT,20000,0.000224000,89285730.462,0
opencl-e2e,Apple M2,GPU,w32-il64,rsa256(composite),256,SUBTRACT,20000,0.000314000,63694266.616,0
library,Apple M2,host-cpu,gmp-1t,rsa256(composite),256,ADDMOD,20000,0.000521000,38387715.266,0
library,Apple M2,host-cpu,gmp-nt,rsa256(composite),256,ADDMOD,20000,0.000196000,102040827.238,0
library,Apple M2,host-cpu,openssl-nt,rsa256(composite),256,ADDMOD,20000,0.000496000,40322581.191,0
opencl-kernel,Apple M2,GPU,w8,rsa256(composite),256,ADDMOD,20000,0.000304000,65789476.606,0
opencl-e2e,Apple M2,GPU,w8,rsa256(composite),256,ADDMOD,20000,0.000381000,52493439.596,0
opencl-kernel,Apple M2,GPU,w16,rsa256(composite),256,ADDMOD,20000,0.000277000,72202162.484,0
opencl-e2e,Apple M2,GPU,w16,rsa256(composite),256,ADDMOD,20000,0.000355000,56338024.107,0
opencl-kernel,Apple M2,GPU,w32,rsa256(composite),256,ADDMOD,20000,0.000207000,96618355.344,0
opencl-e2e,Apple M2,GPU,w32,rsa256(composite),256,ADDMOD,20000,0.000312000,64102566.320,0
opencl-kernel,Apple M2,GPU,w32-opt,rsa256(composite),256,ADDMOD,20000,0.000286000,69930078.279,0
opencl-e2e,Apple M2,GPU,w32-opt,rsa256(composite),256,ADDMOD,20000,0.000377000,53050395.276,0
opencl-kernel,Apple M2,GPU,w32-o64,rsa256(composite),256,ADDMOD,20000,0.000231000,86580081.813,0
opencl-e2e,Apple M2,GPU,w32-o64,rsa256(composite),256,ADDMOD,20000,0.000305000,65573771.757,0
opencl-kernel,Apple M2,GPU,w32-il,rsa256(composite),256,ADDMOD,20000,0.000217000,92165913.771,0
opencl-e2e,Apple M2,GPU,w32-il,rsa256(composite),256,ADDMOD,20000,0.000296000,67567571.268,0
opencl-kernel,Apple M2,GPU,w32-il64,rsa256(composite),256,ADDMOD,20000,0.000198000,101010103.933,0
opencl-e2e,Apple M2,GPU,w32-il64,rsa256(composite),256,ADDMOD,20000,0.000305000,65573771.757,0
library,Apple M2,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACTMOD,20000,0.000553000,36166363.981,0
library,Apple M2,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACTMOD,20000,0.000193000,103626935.270,0
library,Apple M2,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACTMOD,20000,0.000583000,34305316.414,0
opencl-kernel,Apple M2,GPU,w8,rsa256(composite),256,SUBTRACTMOD,20000,0.000303000,66006605.260,0
opencl-e2e,Apple M2,GPU,w8,rsa256(composite),256,SUBTRACTMOD,20000,0.000384000,52083331.491,0
opencl-kernel,Apple M2,GPU,w16,rsa256(composite),256,SUBTRACTMOD,20000,0.000265000,75471687.069,0
opencl-e2e,Apple M2,GPU,w16,rsa256(composite),256,SUBTRACTMOD,20000,0.000360000,55555554.713,0
opencl-kernel,Apple M2,GPU,w32,rsa256(composite),256,SUBTRACTMOD,20000,0.000221000,90497739.694,0
opencl-e2e,Apple M2,GPU,w32,rsa256(composite),256,SUBTRACTMOD,20000,0.000298000,67114094.180,0
opencl-kernel,Apple M2,GPU,w32-opt,rsa256(composite),256,SUBTRACTMOD,20000,0.000283000,70671377.788,0
opencl-e2e,Apple M2,GPU,w32-opt,rsa256(composite),256,SUBTRACTMOD,20000,0.000382000,52356021.168,0
opencl-kernel,Apple M2,GPU,w32-o64,rsa256(composite),256,SUBTRACTMOD,20000,0.000231000,86580081.813,0
opencl-e2e,Apple M2,GPU,w32-o64,rsa256(composite),256,SUBTRACTMOD,20000,0.000306000,65359489.188,0
opencl-kernel,Apple M2,GPU,w32-il,rsa256(composite),256,SUBTRACTMOD,20000,0.000199000,100502511.611,0
opencl-e2e,Apple M2,GPU,w32-il,rsa256(composite),256,SUBTRACTMOD,20000,0.000301000,66445177.898,0
opencl-kernel,Apple M2,GPU,w32-il64,rsa256(composite),256,SUBTRACTMOD,20000,0.000223000,89686118.039,0
opencl-e2e,Apple M2,GPU,w32-il64,rsa256(composite),256,SUBTRACTMOD,20000,0.000301000,66445190.747,0
library,Apple M2,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYOPERANDSCANNING,20000,0.000298000,67114094.180,0
library,Apple M2,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,20000,0.000137000,145985344.285,0
library,Apple M2,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,20000,0.000138000,144927533.016,0
opencl-kernel,Apple M2,GPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,20000,0.005071000,3943995.251,0
opencl-e2e,Apple M2,GPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,20000,0.005165000,3872216.859,0
opencl-kernel,Apple M2,GPU,w16,rsa256(composite),256,MULTIPLYOPERANDSCANNING,20000,0.001335000,14981272.991,0
opencl-e2e,Apple M2,GPU,w16,rsa256(composite),256,MULTIPLYOPERANDSCANNING,20000,0.001437000,13917883.928,0
opencl-kernel,Apple M2,GPU,w32,rsa256(composite),256,MULTIPLYOPERANDSCANNING,20000,0.000492000,40650404.766,0
opencl-e2e,Apple M2,GPU,w32,rsa256(composite),256,MULTIPLYOPERANDSCANNING,20000,0.000585000,34188032.394,0
opencl-kernel,Apple M2,GPU,w32-opt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,20000,0.000296000,67567571.268,0
opencl-e2e,Apple M2,GPU,w32-opt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,20000,0.000400000,49999997.625,0
opencl-kernel,Apple M2,GPU,w32-o64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,20000,0.000238000,84033610.687,0
opencl-e2e,Apple M2,GPU,w32-o64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,20000,0.000340000,58823524.604,0
opencl-kernel,Apple M2,GPU,w32-il,rsa256(composite),256,MULTIPLYOPERANDSCANNING,20000,0.000212000,94339629.558,0
opencl-e2e,Apple M2,GPU,w32-il,rsa256(composite),256,MULTIPLYOPERANDSCANNING,20000,0.000313000,63897764.227,0
opencl-kernel,Apple M2,GPU,w32-il64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,20000,0.000225000,88888901.914,0
opencl-e2e,Apple M2,GPU,w32-il64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,20000,0.000321000,62305296.041,0
library,Apple M2,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,20000,0.000298000,67114094.180,0
library,Apple M2,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,20000,0.000116000,172413777.993,0
library,Apple M2,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,20000,0.000168000,119047619.993,0
opencl-kernel,Apple M2,GPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,20000,0.002529000,7908264.166,0
opencl-e2e,Apple M2,GPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,20000,0.002654000,7535794.996,0
opencl-kernel,Apple M2,GPU,w16,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,20000,0.000796000,25125627.903,0
opencl-e2e,Apple M2,GPU,w16,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,20000,0.000925000,21621621.615,0
opencl-kernel,Apple M2,GPU,w32,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,20000,0.000364000,54945058.310,0
opencl-e2e,Apple M2,GPU,w32,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,20000,0.000466000,42918455.149,0
opencl-kernel,Apple M2,GPU,w32-opt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,20000,0.000365000,54794514.013,0
opencl-e2e,Apple M2,GPU,w32-opt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,20000,0.000468000,42735046.871,0
opencl-kernel,Apple M2,GPU,w32-o64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,20000,0.000382000,52356021.168,0
opencl-e2e,Apple M2,GPU,w32-o64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,20000,0.000475000,42105262.448,0
opencl-kernel,Apple M2,GPU,w32-il,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,20000,0.000270000,74074084.928,0
opencl-e2e,Apple M2,GPU,w32-il,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,20000,0.000344000,58139534.877,0
opencl-kernel,Apple M2,GPU,w32-il64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,20000,0.000269000,74349439.330,0
opencl-e2e,Apple M2,GPU,w32-il64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,20000,0.000357000,56022411.692,0
library,Apple M2,host-cpu,gmp-1t,rsa256(composite),256,MONTGOMERYMULTIPLICATION,20000,0.002823000,7084661.665,0
library,Apple M2,host-cpu,gmp-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,20000,0.000769000,26007803.133,0
library,Apple M2,host-cpu,openssl-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,20000,0.000165000,121212138.973,0
opencl-kernel,Apple M2,GPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,20000,0.000458000,43668122.750,0
opencl-e2e,Apple M2,GPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,20000,0.000522000,38314175.024,0
opencl-kernel,Apple M2,GPU,w16,rsa256(composite),256,MONTGOMERYMULTIPLICATION,20000,0.000268000,74626864.720,0
opencl-e2e,Apple M2,GPU,w16,rsa256(composite),256,MONTGOMERYMULTIPLICATION,20000,0.000355000,56338033.345,0
opencl-kernel,Apple M2,GPU,w32,rsa256(composite),256,MONTGOMERYMULTIPLICATION,20000,0.000220000,90909096.212,0
opencl-e2e,Apple M2,GPU,w32,rsa256(composite),256,MONTGOMERYMULTIPLICATION,20000,0.000289000,69204154.956,0
opencl-kernel,Apple M2,GPU,w32-opt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,20000,0.000319000,62695916.410,0
opencl-e2e,Apple M2,GPU,w32-opt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,20000,0.000420000,47619047.997,0
opencl-kernel,Apple M2,GPU,w32-o64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,20000,0.000239000,83682002.967,0
opencl-e2e,Apple M2,GPU,w32-o64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,20000,0.000318000,62893086.372,0
opencl-kernel,Apple M2,GPU,w32-il,rsa256(composite),256,MONTGOMERYMULTIPLICATION,20000,0.000235000,85106388.422,0
opencl-e2e,Apple M2,GPU,w32-il,rsa256(composite),256,MONTGOMERYMULTIPLICATION,20000,0.000304000,65789476.606,0
opencl-kernel,Apple M2,GPU,w32-il64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,20000,0.000229000,87336245.499,0
opencl-e2e,Apple M2,GPU,w32-il64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,20000,0.000324000,61728390.799,0
library,Apple M2,host-cpu,gmp-1t,rsa256(composite),256,COMPARE,20000,0.000076000,263157906.424,0
library,Apple M2,host-cpu,gmp-nt,rsa256(composite),256,COMPARE,20000,0.000084000,238095239.985,0
library,Apple M2,host-cpu,openssl-nt,rsa256(composite),256,COMPARE,20000,0.000104000,192307734.837,0
opencl-kernel,Apple M2,GPU,w8,rsa256(composite),256,COMPARE,20000,0.000209000,95693770.828,0
opencl-e2e,Apple M2,GPU,w8,rsa256(composite),256,COMPARE,20000,0.000291000,68728521.409,0
opencl-kernel,Apple M2,GPU,w16,rsa256(composite),256,COMPARE,20000,0.000250000,79999986.887,0
opencl-e2e,Apple M2,GPU,w16,rsa256(composite),256,COMPARE,20000,0.000337000,59347180.138,0
opencl-kernel,Apple M2,GPU,w32-opt,rsa256(composite),256,COMPARE,20000,0.000268000,74626864.720,0
opencl-e2e,Apple M2,GPU,w32-opt,rsa256(composite),256,COMPARE,20000,0.000362000,55248624.511,0
opencl-kernel,Apple M2,GPU,w32-o64,rsa256(composite),256,COMPARE,20000,0.000233000,85836899.576,0
opencl-e2e,Apple M2,GPU,w32-o64,rsa256(composite),256,COMPARE,20000,0.000327000,61162081.940,0
opencl-kernel,Apple M2,GPU,w32-il,rsa256(composite),256,COMPARE,20000,0.000213000,93896717.110,0
opencl-e2e,Apple M2,GPU,w32-il,rsa256(composite),256,COMPARE,20000,0.000316000,63291135.302,0
opencl-kernel,Apple M2,GPU,w32-il64,rsa256(composite),256,COMPARE,20000,0.000226000,88495585.149,0
opencl-e2e,Apple M2,GPU,w32-il64,rsa256(composite),256,COMPARE,20000,0.000301000,66445190.747,0
library,Apple M2,host-cpu,gmp-1t,rsa256(composite),256,REDUCE,2500,0.000047000,53191545.465,0
library,Apple M2,host-cpu,gmp-nt,rsa256(composite),256,REDUCE,2500,0.000073000,34246587.642,0
library,Apple M2,host-cpu,openssl-nt,rsa256(composite),256,REDUCE,2500,0.000099000,25252525.983,0
opencl-kernel,Apple M2,GPU,w8,rsa256(composite),256,REDUCE,2500,0.000297000,8417508.661,0
opencl-e2e,Apple M2,GPU,w8,rsa256(composite),256,REDUCE,2500,0.000327000,7645258.882,0
opencl-kernel,Apple M2,GPU,w16,rsa256(composite),256,REDUCE,2500,0.000267000,9363296.027,0
opencl-e2e,Apple M2,GPU,w16,rsa256(composite),256,REDUCE,2500,0.000297000,8417508.661,0
opencl-kernel,Apple M2,GPU,w32-opt,rsa256(composite),256,REDUCE,2500,0.000359000,6963788.343,0
opencl-e2e,Apple M2,GPU,w32-opt,rsa256(composite),256,REDUCE,2500,0.000405000,6172839.523,0
opencl-kernel,Apple M2,GPU,w32-o64,rsa256(composite),256,REDUCE,2500,0.000251000,9960159.744,0
opencl-e2e,Apple M2,GPU,w32-o64,rsa256(composite),256,REDUCE,2500,0.000278000,8992806.948,0
opencl-kernel,Apple M2,GPU,w32-il,rsa256(composite),256,REDUCE,2500,0.000256000,9765623.914,0
opencl-e2e,Apple M2,GPU,w32-il,rsa256(composite),256,REDUCE,2500,0.000287000,8710800.432,0
opencl-kernel,Apple M2,GPU,w32-il64,rsa256(composite),256,REDUCE,2500,0.000281000,8896797.597,0
opencl-e2e,Apple M2,GPU,w32-il64,rsa256(composite),256,REDUCE,2500,0.000305000,8196721.470,0
library,Apple M2,host-cpu,gmp-1t,rsa256(composite),256,MODMUL,1250,0.000088000,14204550.041,0
library,Apple M2,host-cpu,gmp-nt,rsa256(composite),256,MODMUL,1250,0.000080000,15624994.710,0
library,Apple M2,host-cpu,openssl-nt,rsa256(composite),256,MODMUL,1250,0.000135000,9259256.624,0
opencl-kernel,Apple M2,GPU,w8,rsa256(composite),256,MODMUL,1250,0.000375000,3333333.304,0
opencl-e2e,Apple M2,GPU,w8,rsa256(composite),256,MODMUL,1250,0.000401000,3117207.228,0
opencl-kernel,Apple M2,GPU,w16,rsa256(composite),256,MODMUL,1250,0.000302000,4139073.241,0
opencl-e2e,Apple M2,GPU,w16,rsa256(composite),256,MODMUL,1250,0.000334000,3742515.171,0
opencl-kernel,Apple M2,GPU,w32-opt,rsa256(composite),256,MODMUL,1250,0.000470000,2659574.638,0
opencl-e2e,Apple M2,GPU,w32-opt,rsa256(composite),256,MODMUL,1250,0.000499000,2505009.939,0
opencl-kernel,Apple M2,GPU,w32-o64,rsa256(composite),256,MODMUL,1250,0.000298000,4194630.886,0
opencl-e2e,Apple M2,GPU,w32-o64,rsa256(composite),256,MODMUL,1250,0.000334000,3742515.171,0
opencl-kernel,Apple M2,GPU,w32-il,rsa256(composite),256,MODMUL,1250,0.000320000,3906250.099,0
opencl-e2e,Apple M2,GPU,w32-il,rsa256(composite),256,MODMUL,1250,0.000343000,3644314.949,0
opencl-kernel,Apple M2,GPU,w32-il64,rsa256(composite),256,MODMUL,1250,0.000315000,3968254.550,0
opencl-e2e,Apple M2,GPU,w32-il64,rsa256(composite),256,MODMUL,1250,0.000329000,3799392.071,0
library,Apple M2,host-cpu,gmp-1t,rsa256(composite),256,MODEXP,312,0.002915000,107032.589,0
library,Apple M2,host-cpu,gmp-nt,rsa256(composite),256,MODEXP,312,0.000738000,422764.243,0
library,Apple M2,host-cpu,openssl-nt,rsa256(composite),256,MODEXP,312,0.000689000,452830.184,0
opencl-kernel,Apple M2,GPU,w8,rsa256(composite),256,MODEXP,312,0.006707000,46518.563,0
opencl-e2e,Apple M2,GPU,w8,rsa256(composite),256,MODEXP,312,0.006725000,46394.052,0
opencl-kernel,Apple M2,GPU,w16,rsa256(composite),256,MODEXP,312,0.001208000,258278.145,0
opencl-e2e,Apple M2,GPU,w16,rsa256(composite),256,MODEXP,312,0.001226000,254486.141,0
opencl-kernel,Apple M2,GPU,w32-opt,rsa256(composite),256,MODEXP,312,0.001588000,196473.554,0
opencl-e2e,Apple M2,GPU,w32-opt,rsa256(composite),256,MODEXP,312,0.001607000,194150.590,0
opencl-kernel,Apple M2,GPU,w32-o64,rsa256(composite),256,MODEXP,312,0.001001000,311688.313,0
opencl-e2e,Apple M2,GPU,w32-o64,rsa256(composite),256,MODEXP,312,0.001010000,308910.889,0
opencl-kernel,Apple M2,GPU,w32-il,rsa256(composite),256,MODEXP,312,0.001612000,193548.388,0
opencl-e2e,Apple M2,GPU,w32-il,rsa256(composite),256,MODEXP,312,0.001636000,190709.046,0
opencl-kernel,Apple M2,GPU,w32-il64,rsa256(composite),256,MODEXP,312,0.000990000,315151.524,0
opencl-e2e,Apple M2,GPU,w32-il64,rsa256(composite),256,MODEXP,312,0.001007000,309831.169,0
library,Apple M2,host-cpu,gmp-1t,rsa256(composite),256,EXPONENTIATION,312,0.000818000,381418.093,0
library,Apple M2,host-cpu,gmp-nt,rsa256(composite),256,EXPONENTIATION,312,0.000276000,1130434.758,0
library,Apple M2,host-cpu,openssl-nt,rsa256(composite),256,EXPONENTIATION,312,0.001933000,161407.140,0
opencl-kernel,Apple M2,GPU,w8,rsa256(composite),256,EXPONENTIATION,312,0.058835000,5302.966,0
opencl-e2e,Apple M2,GPU,w8,rsa256(composite),256,EXPONENTIATION,312,0.058973000,5290.557,0
opencl-kernel,Apple M2,GPU,w16,rsa256(composite),256,EXPONENTIATION,312,0.009459000,32984.459,0
opencl-e2e,Apple M2,GPU,w16,rsa256(composite),256,EXPONENTIATION,312,0.009503000,32831.737,0
opencl-kernel,Apple M2,GPU,w32-opt,rsa256(composite),256,EXPONENTIATION,312,0.002143000,145590.297,0
opencl-e2e,Apple M2,GPU,w32-opt,rsa256(composite),256,EXPONENTIATION,312,0.002169000,143845.091,0
opencl-kernel,Apple M2,GPU,w32-o64,rsa256(composite),256,EXPONENTIATION,312,0.001789000,174399.105,0
opencl-e2e,Apple M2,GPU,w32-o64,rsa256(composite),256,EXPONENTIATION,312,0.001808000,172566.369,0
opencl-kernel,Apple M2,GPU,w32-il,rsa256(composite),256,EXPONENTIATION,312,0.001805000,172853.185,0
opencl-e2e,Apple M2,GPU,w32-il,rsa256(composite),256,EXPONENTIATION,312,0.001863000,167471.821,0
opencl-kernel,Apple M2,GPU,w32-il64,rsa256(composite),256,EXPONENTIATION,312,0.001750000,178285.715,0
opencl-e2e,Apple M2,GPU,w32-il64,rsa256(composite),256,EXPONENTIATION,312,0.001768000,176470.587,0
library,Apple M2,host-cpu,gmp-1t,rsa256(composite),256,DIVIDE,2500,0.000064000,39062522.303,0
library,Apple M2,host-cpu,gmp-nt,rsa256(composite),256,DIVIDE,2500,0.000058000,43103466.127,0
library,Apple M2,host-cpu,openssl-nt,rsa256(composite),256,DIVIDE,2500,0.000093000,26881717.639,0
opencl-kernel,Apple M2,GPU,w8,rsa256(composite),256,DIVIDE,2500,0.000749000,3337783.716,0
opencl-e2e,Apple M2,GPU,w8,rsa256(composite),256,DIVIDE,2500,0.000786000,3180661.620,0
opencl-kernel,Apple M2,GPU,w16,rsa256(composite),256,DIVIDE,2500,0.000596000,4194630.886,0
opencl-e2e,Apple M2,GPU,w16,rsa256(composite),256,DIVIDE,2500,0.000634000,3943217.688,0
opencl-kernel,Apple M2,GPU,w32-opt,rsa256(composite),256,DIVIDE,2500,0.000365000,6849315.344,0
opencl-e2e,Apple M2,GPU,w32-opt,rsa256(composite),256,DIVIDE,2500,0.000398000,6281406.976,0
opencl-kernel,Apple M2,GPU,w32-o64,rsa256(composite),256,DIVIDE,2500,0.000366000,6830601.225,0
opencl-e2e,Apple M2,GPU,w32-o64,rsa256(composite),256,DIVIDE,2500,0.000406000,6157635.370,0
opencl-kernel,Apple M2,GPU,w32-il,rsa256(composite),256,DIVIDE,2500,0.000356000,7022471.255,0
opencl-e2e,Apple M2,GPU,w32-il,rsa256(composite),256,DIVIDE,2500,0.000377000,6631300.433,0
opencl-kernel,Apple M2,GPU,w32-il64,rsa256(composite),256,DIVIDE,2500,0.000352000,7102272.672,0
opencl-e2e,Apple M2,GPU,w32-il64,rsa256(composite),256,DIVIDE,2500,0.000389000,6426735.327,0
library,Apple M2,host-cpu,gmp-1t,rsa256(composite),256,ISQRT,625,0.000036000,17361091.198,0
library,Apple M2,host-cpu,gmp-nt,rsa256(composite),256,ISQRT,625,0.000059000,10593223.285,0
opencl-kernel,Apple M2,GPU,w8,rsa256(composite),256,ISQRT,625,0.004035000,154894.672,0
opencl-e2e,Apple M2,GPU,w8,rsa256(composite),256,ISQRT,625,0.004052000,154244.817,0
opencl-kernel,Apple M2,GPU,w16,rsa256(composite),256,ISQRT,625,0.002812000,222261.736,0
opencl-e2e,Apple M2,GPU,w16,rsa256(composite),256,ISQRT,625,0.002883000,216788.067,0
opencl-kernel,Apple M2,GPU,w32-opt,rsa256(composite),256,ISQRT,625,0.001519000,411454.909,0
opencl-e2e,Apple M2,GPU,w32-opt,rsa256(composite),256,ISQRT,625,0.001550000,403225.810,0
opencl-kernel,Apple M2,GPU,w32-o64,rsa256(composite),256,ISQRT,625,0.001417000,441072.682,0
opencl-e2e,Apple M2,GPU,w32-o64,rsa256(composite),256,ISQRT,625,0.001454000,429848.690,0
opencl-kernel,Apple M2,GPU,w32-il,rsa256(composite),256,ISQRT,625,0.001505000,415282.394,0
opencl-e2e,Apple M2,GPU,w32-il,rsa256(composite),256,ISQRT,625,0.001526000,409567.503,0
opencl-kernel,Apple M2,GPU,w32-il64,rsa256(composite),256,ISQRT,625,0.001401000,446109.916,0
opencl-e2e,Apple M2,GPU,w32-il64,rsa256(composite),256,ISQRT,625,0.001439000,434329.391,0
library,Apple M2,host-cpu,gmp-1t,rsa256(composite),256,MODMUL_R2,20000,0.001419000,14094432.336,0
library,Apple M2,host-cpu,gmp-nt,rsa256(composite),256,MODMUL_R2,20000,0.000422000,47393370.130,0
library,Apple M2,host-cpu,openssl-nt,rsa256(composite),256,MODMUL_R2,20000,0.001241000,16116035.196,0
opencl-kernel,Apple M2,GPU,w8,rsa256(composite),256,MODMUL_R2,20000,0.000600000,33333331.750,0
opencl-e2e,Apple M2,GPU,w8,rsa256(composite),256,MODMUL_R2,20000,0.000726000,27548209.427,0
opencl-kernel,Apple M2,GPU,w16,rsa256(composite),256,MODMUL_R2,20000,0.000275000,72727273.121,0
opencl-e2e,Apple M2,GPU,w16,rsa256(composite),256,MODMUL_R2,20000,0.000369000,54200539.688,0
opencl-kernel,Apple M2,GPU,w32-opt,rsa256(composite),256,MODMUL_R2,20000,0.000253000,79051381.645,0
opencl-e2e,Apple M2,GPU,w32-opt,rsa256(composite),256,MODMUL_R2,20000,0.000344000,58139534.877,0
opencl-kernel,Apple M2,GPU,w32-o64,rsa256(composite),256,MODMUL_R2,20000,0.000260000,76923076.714,0
opencl-e2e,Apple M2,GPU,w32-o64,rsa256(composite),256,MODMUL_R2,20000,0.000325000,61538466.882,0
opencl-kernel,Apple M2,GPU,w32-il,rsa256(composite),256,MODMUL_R2,20000,0.000257000,77821000.750,0
opencl-e2e,Apple M2,GPU,w32-il,rsa256(composite),256,MODMUL_R2,20000,0.000345000,57971013.206,0
opencl-kernel,Apple M2,GPU,w32-il64,rsa256(composite),256,MODMUL_R2,20000,0.000232000,86206888.996,0
opencl-e2e,Apple M2,GPU,w32-il64,rsa256(composite),256,MODMUL_R2,20000,0.000305000,65573771.757,0
library,Apple M2,host-cpu,gmp-1t,brainpoolP512r1,512,ADD,10000,0.000158000,63291146.960,0
library,Apple M2,host-cpu,gmp-nt,brainpoolP512r1,512,ADD,10000,0.000103000,97087352.632,0
library,Apple M2,host-cpu,openssl-nt,brainpoolP512r1,512,ADD,10000,0.000121000,82644595.151,0
opencl-kernel,Apple M2,GPU,w8,brainpoolP512r1,512,ADD,10000,0.000405000,24691358.094,0
opencl-e2e,Apple M2,GPU,w8,brainpoolP512r1,512,ADD,10000,0.000503000,19880716.165,0
opencl-kernel,Apple M2,GPU,w16,brainpoolP512r1,512,ADD,10000,0.000313000,31948882.113,0
opencl-e2e,Apple M2,GPU,w16,brainpoolP512r1,512,ADD,10000,0.000420000,23809523.999,0
opencl-kernel,Apple M2,GPU,w32,brainpoolP512r1,512,ADD,10000,0.000274000,36496351.578,0
opencl-e2e,Apple M2,GPU,w32,brainpoolP512r1,512,ADD,10000,0.000365000,27397257.007,0
opencl-kernel,Apple M2,GPU,w32-opt,brainpoolP512r1,512,ADD,10000,0.000260000,38461538.357,0
opencl-e2e,Apple M2,GPU,w32-opt,brainpoolP512r1,512,ADD,10000,0.000350000,28571429.590,0
opencl-kernel,Apple M2,GPU,w32-o64,brainpoolP512r1,512,ADD,10000,0.000248000,40322585.923,0
opencl-e2e,Apple M2,GPU,w32-o64,brainpoolP512r1,512,ADD,10000,0.000335000,29850747.185,0
opencl-kernel,Apple M2,GPU,w32-il,brainpoolP512r1,512,ADD,10000,0.000219000,45662104.720,0
opencl-e2e,Apple M2,GPU,w32-il,brainpoolP512r1,512,ADD,10000,0.000294000,34013602.345,0
opencl-kernel,Apple M2,GPU,w32-il64,brainpoolP512r1,512,ADD,10000,0.000206000,48543690.032,0
opencl-e2e,Apple M2,GPU,w32-il64,brainpoolP512r1,512,ADD,10000,0.000277000,36101081.242,0
library,Apple M2,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACT,10000,0.000083000,120481939.732,0
library,Apple M2,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACT,10000,0.000080000,124999957.683,0
library,Apple M2,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACT,10000,0.000103000,97087407.498,0
opencl-kernel,Apple M2,GPU,w8,brainpoolP512r1,512,SUBTRACT,10000,0.000401000,24937654.205,0
opencl-e2e,Apple M2,GPU,w8,brainpoolP512r1,512,SUBTRACT,10000,0.000489000,20449897.826,0
opencl-kernel,Apple M2,GPU,w16,brainpoolP512r1,512,SUBTRACT,10000,0.000283000,35335688.894,0
opencl-e2e,Apple M2,GPU,w16,brainpoolP512r1,512,SUBTRACT,10000,0.000415000,24096384.567,0
opencl-kernel,Apple M2,GPU,w32,brainpoolP512r1,512,SUBTRACT,10000,0.000283000,35335688.894,0
opencl-e2e,Apple M2,GPU,w32,brainpoolP512r1,512,SUBTRACT,10000,0.000377000,26525201.734,0
opencl-kernel,Apple M2,GPU,w32-opt,brainpoolP512r1,512,SUBTRACT,10000,0.000257000,38910509.188,0
opencl-e2e,Apple M2,GPU,w32-opt,brainpoolP512r1,512,SUBTRACT,10000,0.000347000,28818441.903,0
opencl-kernel,Apple M2,GPU,w32-o64,brainpoolP512r1,512,SUBTRACT,10000,0.000247000,40485826.987,0
opencl-e2e,Apple M2,GPU,w32-o64,brainpoolP512r1,512,SUBTRACT,10000,0.000339000,29498528.383,0
opencl-kernel,Apple M2,GPU,w32-il,brainpoolP512r1,512,SUBTRACT,10000,0.000215000,46511626.327,0
opencl-e2e,Apple M2,GPU,w32-il,brainpoolP512r1,512,SUBTRACT,10000,0.000298000,33557047.090,0
opencl-kernel,Apple M2,GPU,w32-il64,brainpoolP512r1,512,SUBTRACT,10000,0.000231000,43290040.907,0
opencl-e2e,Apple M2,GPU,w32-il64,brainpoolP512r1,512,SUBTRACT,10000,0.000296000,33783785.634,0
library,Apple M2,host-cpu,gmp-1t,brainpoolP512r1,512,ADDMOD,10000,0.000270000,37037042.464,0
library,Apple M2,host-cpu,gmp-nt,brainpoolP512r1,512,ADDMOD,10000,0.000124000,80645152.917,0
library,Apple M2,host-cpu,openssl-nt,brainpoolP512r1,512,ADDMOD,10000,0.000314000,31847139.212,0
opencl-kernel,Apple M2,GPU,w8,brainpoolP512r1,512,ADDMOD,10000,0.000414000,24154592.232,0
opencl-e2e,Apple M2,GPU,w8,brainpoolP512r1,512,ADDMOD,10000,0.000501000,19960078.591,0
opencl-kernel,Apple M2,GPU,w16,brainpoolP512r1,512,ADDMOD,10000,0.000321000,31152648.021,0
opencl-e2e,Apple M2,GPU,w16,brainpoolP512r1,512,ADDMOD,10000,0.000406000,24630541.479,0
opencl-kernel,Apple M2,GPU,w32,brainpoolP512r1,512,ADDMOD,10000,0.000331000,30211484.081,0
opencl-e2e,Apple M2,GPU,w32,brainpoolP512r1,512,ADDMOD,10000,0.000407000,24570023.719,0
opencl-kernel,Apple M2,GPU,w32-opt,brainpoolP512r1,512,ADDMOD,10000,0.000224000,44642853.630,0
opencl-e2e,Apple M2,GPU,w32-opt,brainpoolP512r1,512,ADDMOD,10000,0.000308000,32467530.680,0
opencl-kernel,Apple M2,GPU,w32-o64,brainpoolP512r1,512,ADDMOD,10000,0.000224000,44642853.630,0
opencl-e2e,Apple M2,GPU,w32-o64,brainpoolP512r1,512,ADDMOD,10000,0.000311000,32154342.739,0
opencl-kernel,Apple M2,GPU,w32-il,brainpoolP512r1,512,ADDMOD,10000,0.000209000,47846898.740,0
opencl-e2e,Apple M2,GPU,w32-il,brainpoolP512r1,512,ADDMOD,10000,0.000315000,31746030.531,0
opencl-kernel,Apple M2,GPU,w32-il64,brainpoolP512r1,512,ADDMOD,10000,0.000218000,45871553.290,0
opencl-e2e,Apple M2,GPU,w32-il64,brainpoolP512r1,512,ADDMOD,10000,0.000302000,33112579.549,0
library,Apple M2,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACTMOD,10000,0.000273000,36630038.873,0
library,Apple M2,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACTMOD,10000,0.000130000,76923076.714,0
library,Apple M2,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACTMOD,10000,0.000357000,28011205.846,0
opencl-kernel,Apple M2,GPU,w8,brainpoolP512r1,512,SUBTRACTMOD,10000,0.000409000,24449879.477,0
opencl-e2e,Apple M2,GPU,w8,brainpoolP512r1,512,SUBTRACTMOD,10000,0.000493000,20283974.480,0
opencl-kernel,Apple M2,GPU,w16,brainpoolP512r1,512,SUBTRACTMOD,10000,0.000317000,31545738.610,0
opencl-e2e,Apple M2,GPU,w16,brainpoolP512r1,512,SUBTRACTMOD,10000,0.000405000,24691358.094,0
opencl-kernel,Apple M2,GPU,w32,brainpoolP512r1,512,SUBTRACTMOD,10000,0.000358000,27932961.656,0
opencl-e2e,Apple M2,GPU,w32,brainpoolP512r1,512,SUBTRACTMOD,10000,0.000475000,21052631.224,0
opencl-kernel,Apple M2,GPU,w32-opt,brainpoolP512r1,512,SUBTRACTMOD,10000,0.000222000,45045044.559,0
opencl-e2e,Apple M2,GPU,w32-opt,brainpoolP512r1,512,SUBTRACTMOD,10000,0.000299000,33444815.311,0
opencl-kernel,Apple M2,GPU,w32-o64,brainpoolP512r1,512,SUBTRACTMOD,10000,0.000239000,41841011.674,0
opencl-e2e,Apple M2,GPU,w32-o64,brainpoolP512r1,512,SUBTRACTMOD,10000,0.000328000,30487799.967,0
opencl-kernel,Apple M2,GPU,w32-il,brainpoolP512r1,512,SUBTRACTMOD,10000,0.000217000,46082956.885,0
opencl-e2e,Apple M2,GPU,w32-il,brainpoolP512r1,512,SUBTRACTMOD,10000,0.000320000,31250000.789,0
opencl-kernel,Apple M2,GPU,w32-il64,brainpoolP512r1,512,SUBTRACTMOD,10000,0.000224000,44642853.630,0
opencl-e2e,Apple M2,GPU,w32-il64,brainpoolP512r1,512,SUBTRACTMOD,10000,0.000309000,32362463.070,0
library,Apple M2,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,10000,0.000329000,30395136.571,0
library,Apple M2,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,10000,0.000148000,67567557.981,0
library,Apple M2,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,10000,0.000191000,52356021.168,0
opencl-kernel,Apple M2,GPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,10000,0.014485000,690369.348,0
opencl-e2e,Apple M2,GPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,10000,0.014904000,670960.814,0
opencl-kernel,Apple M2,GPU,w16,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,10000,0.003721000,2687449.620,0
opencl-e2e,Apple M2,GPU,w16,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,10000,0.003837000,2606202.725,0
opencl-kernel,Apple M2,GPU,w32,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,10000,0.001115000,8968609.463,0
opencl-e2e,Apple M2,GPU,w32,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,10000,0.001233000,8110300.351,0
opencl-kernel,Apple M2,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,10000,0.000348000,28735629.665,0
opencl-e2e,Apple M2,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,10000,0.000444000,22522522.280,0
opencl-kernel,Apple M2,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,10000,0.000446000,22421523.657,0
opencl-e2e,Apple M2,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,10000,0.000551000,18148820.174,0
opencl-kernel,Apple M2,GPU,w32-il,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,10000,0.000401000,24937657.825,0
opencl-e2e,Apple M2,GPU,w32-il,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,10000,0.000598000,16722407.656,0
opencl-kernel,Apple M2,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,10000,0.000440000,22727274.053,0
opencl-e2e,Apple M2,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,10000,0.000552000,18115941.627,0
library,Apple M2,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,10000,0.000329000,30395136.571,0
library,Apple M2,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,10000,0.000157000,63694254.809,0
library,Apple M2,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,10000,0.000186000,53763435.278,0
opencl-kernel,Apple M2,GPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,10000,0.010768000,928677.560,0
opencl-e2e,Apple M2,GPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,10000,0.010893000,918020.747,0
opencl-kernel,Apple M2,GPU,w16,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,10000,0.002526000,3958828.147,0
opencl-e2e,Apple M2,GPU,w16,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,10000,0.002660000,3759398.499,0
opencl-kernel,Apple M2,GPU,w32,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,10000,0.000818000,12224938.869,0
opencl-e2e,Apple M2,GPU,w32,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,10000,0.000934000,10706637.994,0
opencl-kernel,Apple M2,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,10000,0.001097000,9115770.038,0
opencl-e2e,Apple M2,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,10000,0.001288000,7763974.983,0
opencl-kernel,Apple M2,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,10000,0.001202000,8319467.457,0
opencl-e2e,Apple M2,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,10000,0.001316000,7598784.143,0
opencl-kernel,Apple M2,GPU,w32-il,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,10000,0.000342000,29239762.404,0
opencl-e2e,Apple M2,GPU,w32-il,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,10000,0.000434000,23041475.352,0
opencl-kernel,Apple M2,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,10000,0.000376000,26595748.029,0
opencl-e2e,Apple M2,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,10000,0.000476000,21008402.672,0
library,Apple M2,host-cpu,gmp-1t,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,10000,0.003168000,3156565.690,0
library,Apple M2,host-cpu,gmp-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,10000,0.000866000,11547344.665,0
library,Apple M2,host-cpu,openssl-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,10000,0.000237000,42194092.792,0
opencl-kernel,Apple M2,GPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,10000,0.007076000,1413227.816,0
opencl-e2e,Apple M2,GPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,10000,0.007245000,1380262.241,0
opencl-kernel,Apple M2,GPU,w16,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,10000,0.000303000,33003296.290,0
opencl-e2e,Apple M2,GPU,w16,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,10000,0.000392000,25510203.021,0
opencl-kernel,Apple M2,GPU,w32,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,10000,0.000250000,39999993.443,0
opencl-e2e,Apple M2,GPU,w32,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,10000,0.000335000,29850741.998,0
opencl-kernel,Apple M2,GPU,w32-opt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,10000,0.000273000,36630038.873,0
opencl-e2e,Apple M2,GPU,w32-opt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,10000,0.000344000,29069767.438,0
opencl-kernel,Apple M2,GPU,w32-o64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,10000,0.000272000,36764701.304,0
opencl-e2e,Apple M2,GPU,w32-o64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,10000,0.000350000,28571429.590,0
opencl-kernel,Apple M2,GPU,w32-il,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,10000,0.000388000,25773196.817,0
opencl-e2e,Apple M2,GPU,w32-il,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,10000,0.000487000,20533881.622,0
opencl-kernel,Apple M2,GPU,w32-il64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,10000,0.000458000,21834061.375,0
opencl-e2e,Apple M2,GPU,w32-il64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,10000,0.000544000,18382352.619,0
library,Apple M2,host-cpu,gmp-1t,brainpoolP512r1,512,COMPARE,10000,0.000037000,270270338.218,0
library,Apple M2,host-cpu,gmp-nt,brainpoolP512r1,512,COMPARE,10000,0.000060000,166666529.400,0
library,Apple M2,host-cpu,openssl-nt,brainpoolP512r1,512,COMPARE,10000,0.000085000,117647049.207,0
opencl-kernel,Apple M2,GPU,w8,brainpoolP512r1,512,COMPARE,10000,0.000227000,44052866.879,0
opencl-e2e,Apple M2,GPU,w8,brainpoolP512r1,512,COMPARE,10000,0.000299000,33444815.311,0
opencl-kernel,Apple M2,GPU,w16,brainpoolP512r1,512,COMPARE,10000,0.000230000,43478259.905,0
opencl-e2e,Apple M2,GPU,w16,brainpoolP512r1,512,COMPARE,10000,0.000312000,32051283.160,0
opencl-kernel,Apple M2,GPU,w32-opt,brainpoolP512r1,512,COMPARE,10000,0.000242000,41322317.454,0
opencl-e2e,Apple M2,GPU,w32-opt,brainpoolP512r1,512,COMPARE,10000,0.000329000,30395136.571,0
opencl-kernel,Apple M2,GPU,w32-o64,brainpoolP512r1,512,COMPARE,10000,0.000244000,40983607.348,0
opencl-e2e,Apple M2,GPU,w32-o64,brainpoolP512r1,512,COMPARE,10000,0.000350000,28571429.590,0
opencl-kernel,Apple M2,GPU,w32-il,brainpoolP512r1,512,COMPARE,10000,0.000245000,40816326.047,0
opencl-e2e,Apple M2,GPU,w32-il,brainpoolP512r1,512,COMPARE,10000,0.000347000,28818446.738,0
opencl-kernel,Apple M2,GPU,w32-il64,brainpoolP512r1,512,COMPARE,10000,0.000252000,39682539.998,0
opencl-e2e,Apple M2,GPU,w32-il64,brainpoolP512r1,512,COMPARE,10000,0.000341000,29325510.152,0
library,Apple M2,host-cpu,gmp-1t,brainpoolP512r1,512,REDUCE,1250,0.000023000,54347824.881,0
library,Apple M2,host-cpu,gmp-nt,brainpoolP512r1,512,REDUCE,1250,0.000054000,23148146.550,0
library,Apple M2,host-cpu,openssl-nt,brainpoolP512r1,512,REDUCE,1250,0.000106000,11792453.695,0
opencl-kernel,Apple M2,GPU,w8,brainpoolP512r1,512,REDUCE,1250,0.000430000,2906976.645,0
opencl-e2e,Apple M2,GPU,w8,brainpoolP512r1,512,REDUCE,1250,0.000440000,2840908.881,0
opencl-kernel,Apple M2,GPU,w16,brainpoolP512r1,512,REDUCE,1250,0.000344000,3633720.930,0
opencl-e2e,Apple M2,GPU,w16,brainpoolP512r1,512,REDUCE,1250,0.000375000,3333333.304,0
opencl-kernel,Apple M2,GPU,w32-opt,brainpoolP512r1,512,REDUCE,1250,0.000440000,2840909.257,0
opencl-e2e,Apple M2,GPU,w32-opt,brainpoolP512r1,512,REDUCE,1250,0.000477000,2620545.265,0
opencl-kernel,Apple M2,GPU,w32-o64,brainpoolP512r1,512,REDUCE,1250,0.000550000,2272727.045,0
opencl-e2e,Apple M2,GPU,w32-o64,brainpoolP512r1,512,REDUCE,1250,0.000578000,2162629.842,0
opencl-kernel,Apple M2,GPU,w32-il,brainpoolP512r1,512,REDUCE,1250,0.000510000,2450980.472,0
opencl-e2e,Apple M2,GPU,w32-il,brainpoolP512r1,512,REDUCE,1250,0.000565000,2212389.401,0
opencl-kernel,Apple M2,GPU,w32-il64,brainpoolP512r1,512,REDUCE,1250,0.000312000,4006410.395,0
opencl-e2e,Apple M2,GPU,w32-il64,brainpoolP512r1,512,REDUCE,1250,0.000340000,3676470.288,0
library,Apple M2,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL,625,0.000092000,6793478.110,0
library,Apple M2,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL,625,0.000083000,7530121.233,0
library,Apple M2,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL,625,0.000129000,4844961.513,0
opencl-kernel,Apple M2,GPU,w8,brainpoolP512r1,512,MODMUL,625,0.000836000,747607.637,0
opencl-e2e,Apple M2,GPU,w8,brainpoolP512r1,512,MODMUL,625,0.000866000,721709.042,0
opencl-kernel,Apple M2,GPU,w16,brainpoolP512r1,512,MODMUL,625,0.000498000,1255020.059,0
opencl-e2e,Apple M2,GPU,w16,brainpoolP512r1,512,MODMUL,625,0.000518000,1206563.739,0
opencl-kernel,Apple M2,GPU,w32-opt,brainpoolP512r1,512,MODMUL,625,0.000569000,1098418.229,0
opencl-e2e,Apple M2,GPU,w32-opt,brainpoolP512r1,512,MODMUL,625,0.000594000,1052188.583,0
opencl-kernel,Apple M2,GPU,w32-o64,brainpoolP512r1,512,MODMUL,625,0.000441000,1417233.618,0
opencl-e2e,Apple M2,GPU,w32-o64,brainpoolP512r1,512,MODMUL,625,0.000466000,1341201.723,0
opencl-kernel,Apple M2,GPU,w32-il,brainpoolP512r1,512,MODMUL,625,0.000894000,699105.148,0
opencl-e2e,Apple M2,GPU,w32-il,brainpoolP512r1,512,MODMUL,625,0.000927000,674217.896,0
opencl-kernel,Apple M2,GPU,w32-il64,brainpoolP512r1,512,MODMUL,625,0.000440000,1420454.440,0
opencl-e2e,Apple M2,GPU,w32-il64,brainpoolP512r1,512,MODMUL,625,0.000468000,1335470.049,0
library,Apple M2,host-cpu,gmp-1t,brainpoolP512r1,512,MODEXP,156,0.005927000,26320.229,0
library,Apple M2,host-cpu,gmp-nt,brainpoolP512r1,512,MODEXP,156,0.001480000,105405.403,0
library,Apple M2,host-cpu,openssl-nt,brainpoolP512r1,512,MODEXP,156,0.001171000,133219.469,0
opencl-kernel,Apple M2,GPU,w8,brainpoolP512r1,512,MODEXP,156,0.136530000,1142.606,0
opencl-e2e,Apple M2,GPU,w8,brainpoolP512r1,512,MODEXP,156,0.136519000,1142.698,0
opencl-kernel,Apple M2,GPU,w16,brainpoolP512r1,512,MODEXP,156,0.007399000,21083.930,0
opencl-e2e,Apple M2,GPU,w16,brainpoolP512r1,512,MODEXP,156,0.007433000,20987.488,0
opencl-kernel,Apple M2,GPU,w32-opt,brainpoolP512r1,512,MODEXP,156,0.007081000,22030.787,0
opencl-e2e,Apple M2,GPU,w32-opt,brainpoolP512r1,512,MODEXP,156,0.007105000,21956.369,0
opencl-kernel,Apple M2,GPU,w32-o64,brainpoolP512r1,512,MODEXP,156,0.013979000,11159.597,0
opencl-e2e,Apple M2,GPU,w32-o64,brainpoolP512r1,512,MODEXP,156,0.014015000,11130.931,0
opencl-kernel,Apple M2,GPU,w32-il,brainpoolP512r1,512,MODEXP,156,0.007140000,21848.739,0
opencl-e2e,Apple M2,GPU,w32-il,brainpoolP512r1,512,MODEXP,156,0.007153000,21809.031,0
opencl-kernel,Apple M2,GPU,w32-il64,brainpoolP512r1,512,MODEXP,156,0.013386000,11653.967,0
opencl-e2e,Apple M2,GPU,w32-il64,brainpoolP512r1,512,MODEXP,156,0.013430000,11615.786,0
library,Apple M2,host-cpu,gmp-1t,brainpoolP512r1,512,EXPONENTIATION,156,0.001234000,126418.150,0
library,Apple M2,host-cpu,gmp-nt,brainpoolP512r1,512,EXPONENTIATION,156,0.000386000,404145.048,0
library,Apple M2,host-cpu,openssl-nt,brainpoolP512r1,512,EXPONENTIATION,156,0.003521000,44305.595,0
opencl-kernel,Apple M2,GPU,w8,brainpoolP512r1,512,EXPONENTIATION,156,0.472733000,329.996,0
opencl-e2e,Apple M2,GPU,w8,brainpoolP512r1,512,EXPONENTIATION,156,0.473253000,329.633,0
opencl-kernel,Apple M2,GPU,w16,brainpoolP512r1,512,EXPONENTIATION,156,0.091327000,1708.148,0
opencl-e2e,Apple M2,GPU,w16,brainpoolP512r1,512,EXPONENTIATION,156,0.091445000,1705.943,0
opencl-kernel,Apple M2,GPU,w32-opt,brainpoolP512r1,512,EXPONENTIATION,156,0.024288000,6422.925,0
opencl-e2e,Apple M2,GPU,w32-opt,brainpoolP512r1,512,EXPONENTIATION,156,0.024366000,6402.364,0
opencl-kernel,Apple M2,GPU,w32-o64,brainpoolP512r1,512,EXPONENTIATION,156,0.021385000,7294.833,0
opencl-e2e,Apple M2,GPU,w32-o64,brainpoolP512r1,512,EXPONENTIATION,156,0.021456000,7270.694,0
opencl-kernel,Apple M2,GPU,w32-il,brainpoolP512r1,512,EXPONENTIATION,156,0.024456000,6378.803,0
opencl-e2e,Apple M2,GPU,w32-il,brainpoolP512r1,512,EXPONENTIATION,156,0.024514000,6363.711,0
opencl-kernel,Apple M2,GPU,w32-il64,brainpoolP512r1,512,EXPONENTIATION,156,0.021533000,7244.694,0
opencl-e2e,Apple M2,GPU,w32-il64,brainpoolP512r1,512,EXPONENTIATION,156,0.021588000,7226.237,0
library,Apple M2,host-cpu,gmp-1t,brainpoolP512r1,512,DIVIDE,1250,0.000033000,37878766.704,0
library,Apple M2,host-cpu,gmp-nt,brainpoolP512r1,512,DIVIDE,1250,0.000059000,21186446.569,0
library,Apple M2,host-cpu,openssl-nt,brainpoolP512r1,512,DIVIDE,1250,0.000088000,14204540.646,0
opencl-kernel,Apple M2,GPU,w8,brainpoolP512r1,512,DIVIDE,1250,0.004023000,310713.397,0
opencl-e2e,Apple M2,GPU,w8,brainpoolP512r1,512,DIVIDE,1250,0.004072000,306974.461,0
opencl-kernel,Apple M2,GPU,w16,brainpoolP512r1,512,DIVIDE,1250,0.001136000,1100352.135,0
opencl-e2e,Apple M2,GPU,w16,brainpoolP512r1,512,DIVIDE,1250,0.001158000,1079447.297,0
opencl-kernel,Apple M2,GPU,w32-opt,brainpoolP512r1,512,DIVIDE,1250,0.001000000,1250000.013,0
opencl-e2e,Apple M2,GPU,w32-opt,brainpoolP512r1,512,DIVIDE,1250,0.001027000,1217137.338,0
opencl-kernel,Apple M2,GPU,w32-o64,brainpoolP512r1,512,DIVIDE,1250,0.001028000,1215953.274,0
opencl-e2e,Apple M2,GPU,w32-o64,brainpoolP512r1,512,DIVIDE,1250,0.001054000,1185958.262,0
opencl-kernel,Apple M2,GPU,w32-il,brainpoolP512r1,512,DIVIDE,1250,0.000996000,1255020.059,0
opencl-e2e,Apple M2,GPU,w32-il,brainpoolP512r1,512,DIVIDE,1250,0.001023000,1221896.395,0
opencl-kernel,Apple M2,GPU,w32-il64,brainpoolP512r1,512,DIVIDE,1250,0.001015000,1231527.109,0
opencl-e2e,Apple M2,GPU,w32-il64,brainpoolP512r1,512,DIVIDE,1250,0.001043000,1198465.934,0
library,Apple M2,host-cpu,gmp-1t,brainpoolP512r1,512,ISQRT,312,0.000033000,9454556.846,0
library,Apple M2,host-cpu,gmp-nt,brainpoolP512r1,512,ISQRT,312,0.000053000,5886792.884,0
opencl-kernel,Apple M2,GPU,w8,brainpoolP512r1,512,ISQRT,312,0.027005000,11553.416,0
opencl-e2e,Apple M2,GPU,w8,brainpoolP512r1,512,ISQRT,312,0.027085000,11519.291,0
opencl-kernel,Apple M2,GPU,w16,brainpoolP512r1,512,ISQRT,312,0.006586000,47373.216,0
opencl-e2e,Apple M2,GPU,w16,brainpoolP512r1,512,ISQRT,312,0.006629000,47065.922,0
opencl-kernel,Apple M2,GPU,w32-opt,brainpoolP512r1,512,ISQRT,312,0.005315000,58701.787,0
opencl-e2e,Apple M2,GPU,w32-opt,brainpoolP512r1,512,ISQRT,312,0.005333000,58503.657,0
opencl-kernel,Apple M2,GPU,w32-o64,brainpoolP512r1,512,ISQRT,312,0.005429000,57469.147,0
opencl-e2e,Apple M2,GPU,w32-o64,brainpoolP512r1,512,ISQRT,312,0.005435000,57405.704,0
opencl-kernel,Apple M2,GPU,w32-il,brainpoolP512r1,512,ISQRT,312,0.005289000,58990.357,0
opencl-e2e,Apple M2,GPU,w32-il,brainpoolP512r1,512,ISQRT,312,0.005381000,57981.788,0
opencl-kernel,Apple M2,GPU,w32-il64,brainpoolP512r1,512,ISQRT,312,0.005431000,57447.984,0
opencl-e2e,Apple M2,GPU,w32-il64,brainpoolP512r1,512,ISQRT,312,0.005503000,56696.348,0
library,Apple M2,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL_R2,10000,0.001459000,6854009.634,0
library,Apple M2,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL_R2,10000,0.000422000,23696681.797,0
library,Apple M2,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL_R2,10000,0.001458000,6858710.636,0
opencl-kernel,Apple M2,GPU,w8,brainpoolP512r1,512,MODMUL_R2,10000,0.001574000,6353239.940,0
opencl-e2e,Apple M2,GPU,w8,brainpoolP512r1,512,MODMUL_R2,10000,0.001661000,6020469.747,0
opencl-kernel,Apple M2,GPU,w16,brainpoolP512r1,512,MODMUL_R2,10000,0.000314000,31847133.308,0
opencl-e2e,Apple M2,GPU,w16,brainpoolP512r1,512,MODMUL_R2,10000,0.000406000,24630541.479,0
opencl-kernel,Apple M2,GPU,w32-opt,brainpoolP512r1,512,MODMUL_R2,10000,0.000302000,33112585.931,0
opencl-e2e,Apple M2,GPU,w32-opt,brainpoolP512r1,512,MODMUL_R2,10000,0.000382000,26178010.584,0
opencl-kernel,Apple M2,GPU,w32-o64,brainpoolP512r1,512,MODMUL_R2,10000,0.000286000,34965032.023,0
opencl-e2e,Apple M2,GPU,w32-o64,brainpoolP512r1,512,MODMUL_R2,10000,0.000373000,26809652.335,0
opencl-kernel,Apple M2,GPU,w32-il,brainpoolP512r1,512,MODMUL_R2,10000,0.000291000,34364260.705,0
opencl-e2e,Apple M2,GPU,w32-il,brainpoolP512r1,512,MODMUL_R2,10000,0.000378000,26455024.628,0
opencl-kernel,Apple M2,GPU,w32-il64,brainpoolP512r1,512,MODMUL_R2,10000,0.000342000,29239767.380,0
opencl-e2e,Apple M2,GPU,w32-il64,brainpoolP512r1,512,MODMUL_R2,10000,0.000420000,23809523.999,0
library,Apple M2,host-cpu,gmp-1t,p1024,1024,ADD,5000,0.000091000,54945058.310,0
library,Apple M2,host-cpu,gmp-nt,p1024,1024,ADD,5000,0.000068000,73529418.343,0
library,Apple M2,host-cpu,openssl-nt,p1024,1024,ADD,5000,0.000097000,51546401.368,0
opencl-kernel,Apple M2,GPU,w8,p1024,1024,ADD,5000,0.000411000,12165450.526,0
opencl-e2e,Apple M2,GPU,w8,p1024,1024,ADD,5000,0.000495000,10101010.393,0
opencl-kernel,Apple M2,GPU,w16,p1024,1024,ADD,5000,0.000320000,15625000.395,0
opencl-e2e,Apple M2,GPU,w16,p1024,1024,ADD,5000,0.000414000,12077294.418,0
opencl-kernel,Apple M2,GPU,w32,p1024,1024,ADD,5000,0.000286000,17482519.570,0
opencl-e2e,Apple M2,GPU,w32,p1024,1024,ADD,5000,0.000381000,13123359.899,0
opencl-kernel,Apple M2,GPU,w32-opt,p1024,1024,ADD,5000,0.000262000,19083968.305,0
opencl-e2e,Apple M2,GPU,w32-opt,p1024,1024,ADD,5000,0.000336000,14880952.499,0
opencl-kernel,Apple M2,GPU,w32-o64,p1024,1024,ADD,5000,0.000261000,19157087.512,0
opencl-e2e,Apple M2,GPU,w32-o64,p1024,1024,ADD,5000,0.000345000,14492753.302,0
opencl-kernel,Apple M2,GPU,w32-il,p1024,1024,ADD,5000,0.000208000,24038460.128,0
opencl-e2e,Apple M2,GPU,w32-il,p1024,1024,ADD,5000,0.000282000,17730496.856,0
opencl-kernel,Apple M2,GPU,w32-il64,p1024,1024,ADD,5000,0.000222000,22522522.280,0
opencl-e2e,Apple M2,GPU,w32-il64,p1024,1024,ADD,5000,0.000308000,16233765.340,0
library,Apple M2,host-cpu,gmp-1t,p1024,1024,SUBTRACT,5000,0.000050000,100000053.458,0
library,Apple M2,host-cpu,gmp-nt,p1024,1024,SUBTRACT,5000,0.000080000,62500024.316,0
library,Apple M2,host-cpu,openssl-nt,p1024,1024,SUBTRACT,5000,0.000093000,53763468.928,0
opencl-kernel,Apple M2,GPU,w8,p1024,1024,SUBTRACT,5000,0.000403000,12406948.395,0
opencl-e2e,Apple M2,GPU,w8,p1024,1024,SUBTRACT,5000,0.000497000,10060362.155,0
opencl-kernel,Apple M2,GPU,w16,p1024,1024,SUBTRACT,5000,0.000315000,15873015.266,0
opencl-e2e,Apple M2,GPU,w16,p1024,1024,SUBTRACT,5000,0.000412000,12135922.508,0
opencl-kernel,Apple M2,GPU,w32,p1024,1024,SUBTRACT,5000,0.000331000,15105739.384,0
opencl-e2e,Apple M2,GPU,w32,p1024,1024,SUBTRACT,5000,0.000417000,11990408.427,0
opencl-kernel,Apple M2,GPU,w32-opt,p1024,1024,SUBTRACT,5000,0.000259000,19305019.820,0
opencl-e2e,Apple M2,GPU,w32-opt,p1024,1024,SUBTRACT,5000,0.000339000,14749264.192,0
opencl-kernel,Apple M2,GPU,w32-o64,p1024,1024,SUBTRACT,5000,0.000257000,19455250.188,0
opencl-e2e,Apple M2,GPU,w32-o64,p1024,1024,SUBTRACT,5000,0.000348000,14367817.236,0
opencl-kernel,Apple M2,GPU,w32-il,p1024,1024,SUBTRACT,5000,0.000209000,23923442.707,0
opencl-e2e,Apple M2,GPU,w32-il,p1024,1024,SUBTRACT,5000,0.000310000,16129033.612,0
opencl-kernel,Apple M2,GPU,w32-il64,p1024,1024,SUBTRACT,5000,0.000206000,24271845.016,0
opencl-e2e,Apple M2,GPU,w32-il64,p1024,1024,SUBTRACT,5000,0.000302000,16556292.966,0
library,Apple M2,host-cpu,gmp-1t,p1024,1024,ADDMOD,5000,0.000206000,24271838.158,0
library,Apple M2,host-cpu,gmp-nt,p1024,1024,ADDMOD,5000,0.000121000,41322317.454,0
library,Apple M2,host-cpu,openssl-nt,p1024,1024,ADDMOD,5000,0.000260000,19230769.178,0
opencl-kernel,Apple M2,GPU,w8,p1024,1024,ADDMOD,5000,0.000464000,10775861.125,0
opencl-e2e,Apple M2,GPU,w8,p1024,1024,ADDMOD,5000,0.000555000,9009009.384,0
opencl-kernel,Apple M2,GPU,w16,p1024,1024,ADDMOD,5000,0.000313000,15974441.057,0
opencl-e2e,Apple M2,GPU,w16,p1024,1024,ADDMOD,5000,0.000407000,12285011.860,0
opencl-kernel,Apple M2,GPU,w32,p1024,1024,ADDMOD,5000,0.000403000,12406946.603,0
opencl-e2e,Apple M2,GPU,w32,p1024,1024,ADDMOD,5000,0.000499000,10020040.925,0
opencl-kernel,Apple M2,GPU,w32-opt,p1024,1024,ADDMOD,5000,0.000208000,24038460.128,0
opencl-e2e,Apple M2,GPU,w32-opt,p1024,1024,ADDMOD,5000,0.000281000,17793595.193,0
opencl-kernel,Apple M2,GPU,w32-o64,p1024,1024,ADDMOD,5000,0.000231000,21645020.453,0
opencl-e2e,Apple M2,GPU,w32-o64,p1024,1024,ADDMOD,5000,0.000308000,16233765.340,0
opencl-kernel,Apple M2,GPU,w32-il,p1024,1024,ADDMOD,5000,0.000211000,23696685.065,0
opencl-e2e,Apple M2,GPU,w32-il,p1024,1024,ADDMOD,5000,0.000293000,17064848.691,0
opencl-kernel,Apple M2,GPU,w32-il64,p1024,1024,ADDMOD,5000,0.000228000,21929825.535,0
opencl-e2e,Apple M2,GPU,w32-il64,p1024,1024,ADDMOD,5000,0.000303000,16501648.145,0
library,Apple M2,host-cpu,gmp-1t,p1024,1024,SUBTRACTMOD,5000,0.000172000,29069762.520,0
library,Apple M2,host-cpu,gmp-nt,p1024,1024,SUBTRACTMOD,5000,0.000094000,53191512.527,0
library,Apple M2,host-cpu,openssl-nt,p1024,1024,SUBTRACTMOD,5000,0.000235000,21276597.106,0
opencl-kernel,Apple M2,GPU,w8,p1024,1024,SUBTRACTMOD,5000,0.000464000,10775862.476,0
opencl-e2e,Apple M2,GPU,w8,p1024,1024,SUBTRACTMOD,5000,0.000561000,8912655.613,0
opencl-kernel,Apple M2,GPU,w16,p1024,1024,SUBTRACTMOD,5000,0.000306000,16339869.189,0
opencl-e2e,Apple M2,GPU,w16,p1024,1024,SUBTRACTMOD,5000,0.000395000,12658228.459,0
opencl-kernel,Apple M2,GPU,w32,p1024,1024,SUBTRACTMOD,5000,0.000388000,12886598.409,0
opencl-e2e,Apple M2,GPU,w32,p1024,1024,SUBTRACTMOD,5000,0.000487000,10266940.811,0
opencl-kernel,Apple M2,GPU,w32-opt,p1024,1024,SUBTRACTMOD,5000,0.000209000,23923442.707,0
opencl-e2e,Apple M2,GPU,w32-opt,p1024,1024,SUBTRACTMOD,5000,0.000290000,17241379.530,0
opencl-kernel,Apple M2,GPU,w32-o64,p1024,1024,SUBTRACTMOD,5000,0.000219000,22831052.360,0
opencl-e2e,Apple M2,GPU,w32-o64,p1024,1024,SUBTRACTMOD,5000,0.000307000,16286644.455,0
opencl-kernel,Apple M2,GPU,w32-il,p1024,1024,SUBTRACTMOD,5000,0.000194000,25773200.684,0
opencl-e2e,Apple M2,GPU,w32-il,p1024,1024,SUBTRACTMOD,5000,0.000274000,18248175.789,0
opencl-kernel,Apple M2,GPU,w32-il64,p1024,1024,SUBTRACTMOD,5000,0.000227000,22026427.791,0
opencl-e2e,Apple M2,GPU,w32-il64,p1024,1024,SUBTRACTMOD,5000,0.000312000,16025641.580,0
library,Apple M2,host-cpu,gmp-1t,p1024,1024,MULTIPLYOPERANDSCANNING,5000,0.000577000,8665511.719,0
library,Apple M2,host-cpu,gmp-nt,p1024,1024,MULTIPLYOPERANDSCANNING,5000,0.000190000,26315790.642,0
library,Apple M2,host-cpu,openssl-nt,p1024,1024,MULTIPLYOPERANDSCANNING,5000,0.000310000,16129030.583,0
opencl-kernel,Apple M2,GPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,5000,0.028689000,174282.826,0
opencl-e2e,Apple M2,GPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,5000,0.028888000,173082.249,0
opencl-kernel,Apple M2,GPU,w16,p1024,1024,MULTIPLYOPERANDSCANNING,5000,0.007228000,691754.291,0
opencl-e2e,Apple M2,GPU,w16,p1024,1024,MULTIPLYOPERANDSCANNING,5000,0.007325000,682593.855,0
opencl-kernel,Apple M2,GPU,w32,p1024,1024,MULTIPLYOPERANDSCANNING,5000,0.001955000,2557544.776,0
opencl-e2e,Apple M2,GPU,w32,p1024,1024,MULTIPLYOPERANDSCANNING,5000,0.002099000,2382086.669,0
opencl-kernel,Apple M2,GPU,w32-opt,p1024,1024,MULTIPLYOPERANDSCANNING,5000,0.000379000,13192610.964,0
opencl-e2e,Apple M2,GPU,w32-opt,p1024,1024,MULTIPLYOPERANDSCANNING,5000,0.000475000,10526315.612,0
opencl-kernel,Apple M2,GPU,w32-o64,p1024,1024,MULTIPLYOPERANDSCANNING,5000,0.000392000,12755103.405,0
opencl-e2e,Apple M2,GPU,w32-o64,p1024,1024,MULTIPLYOPERANDSCANNING,5000,0.000495000,10101010.393,0
opencl-kernel,Apple M2,GPU,w32-il,p1024,1024,MULTIPLYOPERANDSCANNING,5000,0.000419000,11933174.536,0
opencl-e2e,Apple M2,GPU,w32-il,p1024,1024,MULTIPLYOPERANDSCANNING,5000,0.000509000,9823183.177,0
opencl-kernel,Apple M2,GPU,w32-il64,p1024,1024,MULTIPLYOPERANDSCANNING,5000,0.000418000,11961723.019,0
opencl-e2e,Apple M2,GPU,w32-il64,p1024,1024,MULTIPLYOPERANDSCANNING,5000,0.000515000,9708737.458,0
library,Apple M2,host-cpu,gmp-1t,p1024,1024,MULTIPLYPRODUCTSCANNING,5000,0.000576000,8680555.249,0
library,Apple M2,host-cpu,gmp-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,5000,0.000252000,19841269.999,0
library,Apple M2,host-cpu,openssl-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,5000,0.000315000,15873015.266,0
opencl-kernel,Apple M2,GPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,5000,0.017247000,289905.491,0
opencl-e2e,Apple M2,GPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,5000,0.017425000,286944.045,0
opencl-kernel,Apple M2,GPU,w16,p1024,1024,MULTIPLYPRODUCTSCANNING,5000,0.004529000,1103996.461,0
opencl-e2e,Apple M2,GPU,w16,p1024,1024,MULTIPLYPRODUCTSCANNING,5000,0.004655000,1074113.857,0
opencl-kernel,Apple M2,GPU,w32,p1024,1024,MULTIPLYPRODUCTSCANNING,5000,0.001309000,3819709.662,0
opencl-e2e,Apple M2,GPU,w32,p1024,1024,MULTIPLYPRODUCTSCANNING,5000,0.001413000,3538570.440,0
opencl-kernel,Apple M2,GPU,w32-opt,p1024,1024,MULTIPLYPRODUCTSCANNING,5000,0.001304000,3834355.728,0
opencl-e2e,Apple M2,GPU,w32-opt,p1024,1024,MULTIPLYPRODUCTSCANNING,5000,0.001414000,3536067.896,0
opencl-kernel,Apple M2,GPU,w32-o64,p1024,1024,MULTIPLYPRODUCTSCANNING,5000,0.000622000,8038584.933,0
opencl-e2e,Apple M2,GPU,w32-o64,p1024,1024,MULTIPLYPRODUCTSCANNING,5000,0.000733000,6821281.920,0
opencl-kernel,Apple M2,GPU,w32-il,p1024,1024,MULTIPLYPRODUCTSCANNING,5000,0.000611000,8183306.111,0
opencl-e2e,Apple M2,GPU,w32-il,p1024,1024,MULTIPLYPRODUCTSCANNING,5000,0.000712000,7022472.403,0
opencl-kernel,Apple M2,GPU,w32-il64,p1024,1024,MULTIPLYPRODUCTSCANNING,5000,0.000538000,9293679.916,0
opencl-e2e,Apple M2,GPU,w32-il64,p1024,1024,MULTIPLYPRODUCTSCANNING,5000,0.000642000,7788162.005,0
library,Apple M2,host-cpu,gmp-1t,p1024,1024,MONTGOMERYMULTIPLICATION,5000,0.004802000,1041232.817,0
library,Apple M2,host-cpu,gmp-nt,p1024,1024,MONTGOMERYMULTIPLICATION,5000,0.001206000,4145937.029,0
library,Apple M2,host-cpu,openssl-nt,p1024,1024,MONTGOMERYMULTIPLICATION,5000,0.000379000,13192610.964,0
opencl-kernel,Apple M2,GPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,5000,0.008600000,581395.349,0
opencl-e2e,Apple M2,GPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,5000,0.008705000,574382.540,0
opencl-kernel,Apple M2,GPU,w16,p1024,1024,MONTGOMERYMULTIPLICATION,5000,0.003301000,1514692.522,0
opencl-e2e,Apple M2,GPU,w16,p1024,1024,MONTGOMERYMULTIPLICATION,5000,0.003393000,1473622.142,0
opencl-kernel,Apple M2,GPU,w32,p1024,1024,MONTGOMERYMULTIPLICATION,5000,0.000366000,13661202.449,0
opencl-e2e,Apple M2,GPU,w32,p1024,1024,MONTGOMERYMULTIPLICATION,5000,0.000452000,11061946.719,0
opencl-kernel,Apple M2,GPU,w32-opt,p1024,1024,MONTGOMERYMULTIPLICATION,5000,0.000298000,16778523.545,0
opencl-e2e,Apple M2,GPU,w32-opt,p1024,1024,MONTGOMERYMULTIPLICATION,5000,0.000390000,12820512.786,0
opencl-kernel,Apple M2,GPU,w32-o64,p1024,1024,MONTGOMERYMULTIPLICATION,5000,0.000275000,18181818.280,0
opencl-e2e,Apple M2,GPU,w32-o64,p1024,1024,MONTGOMERYMULTIPLICATION,5000,0.000372000,13440858.820,0
opencl-kernel,Apple M2,GPU,w32-il,p1024,1024,MONTGOMERYMULTIPLICATION,5000,0.000307000,16286644.455,0
opencl-e2e,Apple M2,GPU,w32-il,p1024,1024,MONTGOMERYMULTIPLICATION,5000,0.000388000,12886598.409,0
opencl-kernel,Apple M2,GPU,w32-il64,p1024,1024,MONTGOMERYMULTIPLICATION,5000,0.000289000,17301035.254,0
opencl-e2e,Apple M2,GPU,w32-il64,p1024,1024,MONTGOMERYMULTIPLICATION,5000,0.000369000,13550134.922,0
library,Apple M2,host-cpu,gmp-1t,p1024,1024,COMPARE,5000,0.000018000,277777459.174,0
library,Apple M2,host-cpu,gmp-nt,p1024,1024,COMPARE,5000,0.000048000,104166757.721,0
library,Apple M2,host-cpu,openssl-nt,p1024,1024,COMPARE,5000,0.000070000,71428562.096,0
opencl-kernel,Apple M2,GPU,w8,p1024,1024,COMPARE,5000,0.000219000,22831052.360,0
opencl-e2e,Apple M2,GPU,w8,p1024,1024,COMPARE,5000,0.000314000,15923569.606,0
opencl-kernel,Apple M2,GPU,w16,p1024,1024,COMPARE,5000,0.000218000,22935782.769,0
opencl-e2e,Apple M2,GPU,w16,p1024,1024,COMPARE,5000,0.000312000,16025638.590,0
opencl-kernel,Apple M2,GPU,w32-opt,p1024,1024,COMPARE,5000,0.000197000,25380712.379,0
opencl-e2e,Apple M2,GPU,w32-opt,p1024,1024,COMPARE,5000,0.000292000,17123286.994,0
opencl-kernel,Apple M2,GPU,w32-o64,p1024,1024,COMPARE,5000,0.000215000,23255813.164,0
opencl-e2e,Apple M2,GPU,w32-o64,p1024,1024,COMPARE,5000,0.000302000,16556292.966,0
opencl-kernel,Apple M2,GPU,w32-il,p1024,1024,COMPARE,5000,0.000192000,26041665.746,0
opencl-e2e,Apple M2,GPU,w32-il,p1024,1024,COMPARE,5000,0.000273000,18315019.437,0
opencl-kernel,Apple M2,GPU,w32-il64,p1024,1024,COMPARE,5000,0.000184000,27173912.440,0
opencl-e2e,Apple M2,GPU,w32-il64,p1024,1024,COMPARE,5000,0.000292000,17123290.408,0
library,Apple M2,host-cpu,gmp-1t,p1024,1024,REDUCE,625,0.000010000,62500251.691,0
library,Apple M2,host-cpu,gmp-nt,p1024,1024,REDUCE,625,0.000052000,12019219.973,0
library,Apple M2,host-cpu,openssl-nt,p1024,1024,REDUCE,625,0.000087000,7183906.215,0
opencl-kernel,Apple M2,GPU,w8,p1024,1024,REDUCE,625,0.001886000,331389.186,0
opencl-e2e,Apple M2,GPU,w8,p1024,1024,REDUCE,625,0.001934000,323164.426,0
opencl-kernel,Apple M2,GPU,w16,p1024,1024,REDUCE,625,0.000584000,1070205.437,0
opencl-e2e,Apple M2,GPU,w16,p1024,1024,REDUCE,625,0.000614000,1017915.278,0
opencl-kernel,Apple M2,GPU,w32-opt,p1024,1024,REDUCE,625,0.000416000,1502403.758,0
opencl-e2e,Apple M2,GPU,w32-opt,p1024,1024,REDUCE,625,0.000445000,1404494.343,0
opencl-kernel,Apple M2,GPU,w32-o64,p1024,1024,REDUCE,625,0.000442000,1414027.183,0
opencl-e2e,Apple M2,GPU,w32-o64,p1024,1024,REDUCE,625,0.000468000,1335470.049,0
opencl-kernel,Apple M2,GPU,w32-il,p1024,1024,REDUCE,625,0.000436000,1433486.232,0
opencl-e2e,Apple M2,GPU,w32-il,p1024,1024,REDUCE,625,0.000464000,1346982.810,0
opencl-kernel,Apple M2,GPU,w32-il64,p1024,1024,REDUCE,625,0.000444000,1407657.642,0
opencl-e2e,Apple M2,GPU,w32-il64,p1024,1024,REDUCE,625,0.000469000,1332622.708,0
library,Apple M2,host-cpu,gmp-1t,p1024,1024,MODMUL,312,0.000124000,2516128.771,0
library,Apple M2,host-cpu,gmp-nt,p1024,1024,MODMUL,312,0.000085000,3670587.935,0
library,Apple M2,host-cpu,openssl-nt,p1024,1024,MODMUL,312,0.000176000,1772727.259,0
opencl-kernel,Apple M2,GPU,w8,p1024,1024,MODMUL,312,0.008922000,34969.738,0
opencl-e2e,Apple M2,GPU,w8,p1024,1024,MODMUL,312,0.008947000,34872.024,0
opencl-kernel,Apple M2,GPU,w16,p1024,1024,MODMUL,312,0.003407000,91576.167,0
opencl-e2e,Apple M2,GPU,w16,p1024,1024,MODMUL,312,0.003441000,90671.316,0
opencl-kernel,Apple M2,GPU,w32-opt,p1024,1024,MODMUL,312,0.000809000,385661.315,0
opencl-e2e,Apple M2,GPU,w32-opt,p1024,1024,MODMUL,312,0.000804000,388059.697,0
opencl-kernel,Apple M2,GPU,w32-o64,p1024,1024,MODMUL,312,0.000868000,359447.015,0
opencl-e2e,Apple M2,GPU,w32-o64,p1024,1024,MODMUL,312,0.000896000,348214.281,0
opencl-kernel,Apple M2,GPU,w32-il,p1024,1024,MODMUL,312,0.000777000,401544.412,0
opencl-e2e,Apple M2,GPU,w32-il,p1024,1024,MODMUL,312,0.000805000,387577.659,0
opencl-kernel,Apple M2,GPU,w32-il64,p1024,1024,MODMUL,312,0.000844000,369668.236,0
opencl-e2e,Apple M2,GPU,w32-il64,p1024,1024,MODMUL,312,0.000873000,357388.311,0
library,Apple M2,host-cpu,gmp-1t,p1024,1024,MODEXP,78,0.015647000,4984.981,0
library,Apple M2,host-cpu,gmp-nt,p1024,1024,MODEXP,78,0.004107000,18991.965,0
library,Apple M2,host-cpu,openssl-nt,p1024,1024,MODEXP,78,0.003504000,22260.274,0
opencl-kernel,Apple M2,GPU,w8,p1024,1024,MODEXP,78,4.966830000,15.704,0
opencl-e2e,Apple M2,GPU,w8,p1024,1024,MODEXP,78,4.967594000,15.702,0
opencl-kernel,Apple M2,GPU,w16,p1024,1024,MODEXP,78,0.410812000,189.868,0
opencl-e2e,Apple M2,GPU,w16,p1024,1024,MODEXP,78,0.412499000,189.091,0
opencl-kernel,Apple M2,GPU,w32-opt,p1024,1024,MODEXP,78,0.103597000,752.918,0
opencl-e2e,Apple M2,GPU,w32-opt,p1024,1024,MODEXP,78,0.103689000,752.250,0
opencl-kernel,Apple M2,GPU,w32-o64,p1024,1024,MODEXP,78,0.044986000,1733.873,0
opencl-e2e,Apple M2,GPU,w32-o64,p1024,1024,MODEXP,78,0.045021000,1732.525,0
opencl-kernel,Apple M2,GPU,w32-il,p1024,1024,MODEXP,78,0.120913000,645.092,0
opencl-e2e,Apple M2,GPU,w32-il,p1024,1024,MODEXP,78,0.120954000,644.873,0
opencl-kernel,Apple M2,GPU,w32-il64,p1024,1024,MODEXP,78,0.043339000,1799.765,0
opencl-e2e,Apple M2,GPU,w32-il64,p1024,1024,MODEXP,78,0.043373000,1798.354,0
library,Apple M2,host-cpu,gmp-1t,p1024,1024,EXPONENTIATION,78,0.001958000,39836.568,0
library,Apple M2,host-cpu,gmp-nt,p1024,1024,EXPONENTIATION,78,0.000514000,151750.969,0
library,Apple M2,host-cpu,openssl-nt,p1024,1024,EXPONENTIATION,78,0.010206000,7642.563,0
opencl-kernel,Apple M2,GPU,w8,p1024,1024,EXPONENTIATION,78,3.196950000,24.398,0
opencl-e2e,Apple M2,GPU,w8,p1024,1024,EXPONENTIATION,78,3.196277000,24.403,0
opencl-kernel,Apple M2,GPU,w16,p1024,1024,EXPONENTIATION,78,0.942950000,82.719,0
opencl-e2e,Apple M2,GPU,w16,p1024,1024,EXPONENTIATION,78,0.942911000,82.723,0
opencl-kernel,Apple M2,GPU,w32-opt,p1024,1024,EXPONENTIATION,78,0.193262000,403.597,0
opencl-e2e,Apple M2,GPU,w32-opt,p1024,1024,EXPONENTIATION,78,0.193300000,403.518,0
opencl-kernel,Apple M2,GPU,w32-o64,p1024,1024,EXPONENTIATION,78,0.225479000,345.930,0
opencl-e2e,Apple M2,GPU,w32-o64,p1024,1024,EXPONENTIATION,78,0.225523000,345.863,0
opencl-kernel,Apple M2,GPU,w32-il,p1024,1024,EXPONENTIATION,78,0.196513000,396.920,0
opencl-e2e,Apple M2,GPU,w32-il,p1024,1024,EXPONENTIATION,78,0.196604000,396.737,0
opencl-kernel,Apple M2,GPU,w32-il64,p1024,1024,EXPONENTIATION,78,0.217087000,359.303,0
opencl-e2e,Apple M2,GPU,w32-il64,p1024,1024,EXPONENTIATION,78,0.216933000,359.558,0
library,Apple M2,host-cpu,gmp-1t,p1024,1024,DIVIDE,625,0.000019000,32894788.691,0
library,Apple M2,host-cpu,gmp-nt,p1024,1024,DIVIDE,625,0.000075000,8333334.554,0
library,Apple M2,host-cpu,openssl-nt,p1024,1024,DIVIDE,625,0.000080000,7812503.040,0
opencl-kernel,Apple M2,GPU,w8,p1024,1024,DIVIDE,625,0.022023000,28379.422,0
opencl-e2e,Apple M2,GPU,w8,p1024,1024,DIVIDE,625,0.022093000,28289.503,0
opencl-kernel,Apple M2,GPU,w16,p1024,1024,DIVIDE,625,0.009379000,66638.234,0
opencl-e2e,Apple M2,GPU,w16,p1024,1024,DIVIDE,625,0.009425000,66312.997,0
opencl-kernel,Apple M2,GPU,w32-opt,p1024,1024,DIVIDE,625,0.012935000,48318.516,0
opencl-e2e,Apple M2,GPU,w32-opt,p1024,1024,DIVIDE,625,0.012952000,48255.096,0
opencl-kernel,Apple M2,GPU,w32-o64,p1024,1024,DIVIDE,625,0.007366000,84849.308,0
opencl-e2e,Apple M2,GPU,w32-o64,p1024,1024,DIVIDE,625,0.007519000,83122.756,0
opencl-kernel,Apple M2,GPU,w32-il,p1024,1024,DIVIDE,625,0.004975000,125628.140,0
opencl-e2e,Apple M2,GPU,w32-il,p1024,1024,DIVIDE,625,0.005014000,124650.977,0
opencl-kernel,Apple M2,GPU,w32-il64,p1024,1024,DIVIDE,625,0.004552000,137302.286,0
opencl-e2e,Apple M2,GPU,w32-il64,p1024,1024,DIVIDE,625,0.004635000,134843.581,0
library,Apple M2,host-cpu,gmp-1t,p1024,1024,ISQRT,156,0.000033000,4727278.423,0
library,Apple M2,host-cpu,gmp-nt,p1024,1024,ISQRT,156,0.000067000,2328358.685,0
opencl-kernel,Apple M2,GPU,w8,p1024,1024,ISQRT,156,0.140045000,1113.928,0
opencl-e2e,Apple M2,GPU,w8,p1024,1024,ISQRT,156,0.140109000,1113.419,0
opencl-kernel,Apple M2,GPU,w16,p1024,1024,ISQRT,156,0.032897000,4742.074,0
opencl-e2e,Apple M2,GPU,w16,p1024,1024,ISQRT,156,0.032899000,4741.785,0
opencl-kernel,Apple M2,GPU,w32-opt,p1024,1024,ISQRT,156,0.028179000,5536.037,0
opencl-e2e,Apple M2,GPU,w32-opt,p1024,1024,ISQRT,156,0.028113000,5549.034,0
opencl-kernel,Apple M2,GPU,w32-o64,p1024,1024,ISQRT,156,0.036023000,4330.567,0
opencl-e2e,Apple M2,GPU,w32-o64,p1024,1024,ISQRT,156,0.036131000,4317.622,0
opencl-kernel,Apple M2,GPU,w32-il,p1024,1024,ISQRT,156,0.028763000,5423.635,0
opencl-e2e,Apple M2,GPU,w32-il,p1024,1024,ISQRT,156,0.028792000,5418.172,0
opencl-kernel,Apple M2,GPU,w32-il64,p1024,1024,ISQRT,156,0.035904000,4344.920,0
opencl-e2e,Apple M2,GPU,w32-il64,p1024,1024,ISQRT,156,0.035823000,4354.744,0
library,Apple M2,host-cpu,gmp-1t,p1024,1024,MODMUL_R2,5000,0.002034000,2458210.417,0
library,Apple M2,host-cpu,gmp-nt,p1024,1024,MODMUL_R2,5000,0.000594000,8417508.661,0
library,Apple M2,host-cpu,openssl-nt,p1024,1024,MODMUL_R2,5000,0.001593000,3138731.908,0
opencl-kernel,Apple M2,GPU,w8,p1024,1024,MODMUL_R2,5000,0.010947000,456746.140,0
opencl-e2e,Apple M2,GPU,w8,p1024,1024,MODMUL_R2,5000,0.011082000,451182.098,0
opencl-kernel,Apple M2,GPU,w16,p1024,1024,MODMUL_R2,5000,0.001706000,2930832.339,0
opencl-e2e,Apple M2,GPU,w16,p1024,1024,MODMUL_R2,5000,0.001782000,2805836.129,0
opencl-kernel,Apple M2,GPU,w32-opt,p1024,1024,MODMUL_R2,5000,0.000455000,10989011.662,0
opencl-e2e,Apple M2,GPU,w32-opt,p1024,1024,MODMUL_R2,5000,0.000535000,9345794.406,0
opencl-kernel,Apple M2,GPU,w32-o64,p1024,1024,MODMUL_R2,5000,0.000337000,14836795.034,0
opencl-e2e,Apple M2,GPU,w32-o64,p1024,1024,MODMUL_R2,5000,0.000425000,11764706.532,0
opencl-kernel,Apple M2,GPU,w32-il,p1024,1024,MODMUL_R2,5000,0.000443000,11286681.788,0
opencl-e2e,Apple M2,GPU,w32-il,p1024,1024,MODMUL_R2,5000,0.000531000,9416196.411,0
opencl-kernel,Apple M2,GPU,w32-il64,p1024,1024,MODMUL_R2,5000,0.000355000,14084506.027,0
opencl-e2e,Apple M2,GPU,w32-il64,p1024,1024,MODMUL_R2,5000,0.000433000,11547344.665,0
library,Apple M2,host-cpu,gmp-1t,p2048,2048,ADD,2500,0.000060000,41666672.772,0
library,Apple M2,host-cpu,gmp-nt,p2048,2048,ADD,2500,0.000075000,33333338.218,0
library,Apple M2,host-cpu,openssl-nt,p2048,2048,ADD,2500,0.000105000,23809527.298,0
opencl-kernel,Apple M2,GPU,w8,p2048,2048,ADD,2500,0.000402000,6218905.844,0
opencl-e2e,Apple M2,GPU,w8,p2048,2048,ADD,2500,0.000490000,5102040.756,0
opencl-kernel,Apple M2,GPU,w16,p2048,2048,ADD,2500,0.000305000,8196721.470,0
opencl-e2e,Apple M2,GPU,w16,p2048,2048,ADD,2500,0.000404000,6188118.946,0
opencl-kernel,Apple M2,GPU,w32,p2048,2048,ADD,2500,0.000384000,6510417.423,0
opencl-e2e,Apple M2,GPU,w32,p2048,2048,ADD,2500,0.000488000,5122950.919,0
opencl-kernel,Apple M2,GPU,w32-opt,p2048,2048,ADD,2500,0.000261000,9578543.756,0
opencl-e2e,Apple M2,GPU,w32-opt,p2048,2048,ADD,2500,0.000347000,7204611.684,0
opencl-kernel,Apple M2,GPU,w32-o64,p2048,2048,ADD,2500,0.000269000,9293679.916,0
opencl-e2e,Apple M2,GPU,w32-o64,p2048,2048,ADD,2500,0.000367000,6811989.091,0
opencl-kernel,Apple M2,GPU,w32-il,p2048,2048,ADD,2500,0.000216000,11574073.275,0
opencl-e2e,Apple M2,GPU,w32-il,p2048,2048,ADD,2500,0.000295000,8474575.283,0
opencl-kernel,Apple M2,GPU,w32-il64,p2048,2048,ADD,2500,0.000223000,11210761.829,0
opencl-e2e,Apple M2,GPU,w32-il64,p2048,2048,ADD,2500,0.000297000,8417508.661,0
library,Apple M2,host-cpu,gmp-1t,p2048,2048,SUBTRACT,2500,0.000036000,69444477.077,0
library,Apple M2,host-cpu,gmp-nt,p2048,2048,SUBTRACT,2500,0.000072000,34722210.468,0
library,Apple M2,host-cpu,openssl-nt,p2048,2048,SUBTRACT,2500,0.000098000,25510206.809,0
opencl-kernel,Apple M2,GPU,w8,p2048,2048,SUBTRACT,2500,0.000409000,6112469.869,0
opencl-e2e,Apple M2,GPU,w8,p2048,2048,SUBTRACT,2500,0.000509000,4911591.589,0
opencl-kernel,Apple M2,GPU,w16,p2048,2048,SUBTRACT,2500,0.000303000,8250825.657,0
opencl-e2e,Apple M2,GPU,w16,p2048,2048,SUBTRACT,2500,0.000405000,6172839.523,0
opencl-kernel,Apple M2,GPU,w32,p2048,2048,SUBTRACT,2500,0.000370000,6756757.392,0
opencl-e2e,Apple M2,GPU,w32,p2048,2048,SUBTRACT,2500,0.000467000,5353318.997,0
opencl-kernel,Apple M2,GPU,w32-opt,p2048,2048,SUBTRACT,2500,0.000232000,10775861.125,0
opencl-e2e,Apple M2,GPU,w32-opt,p2048,2048,SUBTRACT,2500,0.000332000,7530119.913,0
opencl-kernel,Apple M2,GPU,w32-o64,p2048,2048,SUBTRACT,2500,0.000240000,10416665.667,0
opencl-e2e,Apple M2,GPU,w32-o64,p2048,2048,SUBTRACT,2500,0.000332000,7530121.233,0
opencl-kernel,Apple M2,GPU,w32-il,p2048,2048,SUBTRACT,2500,0.000222000,11261261.140,0
opencl-e2e,Apple M2,GPU,w32-il,p2048,2048,SUBTRACT,2500,0.000296000,8445946.408,0
opencl-kernel,Apple M2,GPU,w32-il64,p2048,2048,SUBTRACT,2500,0.000233000,10729612.447,0
opencl-e2e,Apple M2,GPU,w32-il64,p2048,2048,SUBTRACT,2500,0.000315000,7936509.099,0
library,Apple M2,host-cpu,gmp-1t,p2048,2048,ADDMOD,2500,0.000118000,21186446.569,0
library,Apple M2,host-cpu,gmp-nt,p2048,2048,ADDMOD,2500,0.000089000,28089875.835,0
library,Apple M2,host-cpu,openssl-nt,p2048,2048,ADDMOD,2500,0.000267000,9363296.027,0
opencl-kernel,Apple M2,GPU,w8,p2048,2048,ADDMOD,2500,0.000466000,5364806.894,0
opencl-e2e,Apple M2,GPU,w8,p2048,2048,ADDMOD,2500,0.000550000,4545454.570,0
opencl-kernel,Apple M2,GPU,w16,p2048,2048,ADDMOD,2500,0.000365000,6849315.344,0
opencl-e2e,Apple M2,GPU,w16,p2048,2048,ADDMOD,2500,0.000446000,5605381.646,0
opencl-kernel,Apple M2,GPU,w32,p2048,2048,ADDMOD,2500,0.000286000,8741259.785,0
opencl-e2e,Apple M2,GPU,w32,p2048,2048,ADDMOD,2500,0.000370000,6756757.392,0
opencl-kernel,Apple M2,GPU,w32-opt,p2048,2048,ADDMOD,2500,0.000227000,11013216.720,0
opencl-e2e,Apple M2,GPU,w32-opt,p2048,2048,ADDMOD,2500,0.000312000,8012820.790,0
opencl-kernel,Apple M2,GPU,w32-o64,p2048,2048,ADDMOD,2500,0.000234000,10683761.718,0
opencl-e2e,Apple M2,GPU,w32-o64,p2048,2048,ADDMOD,2500,0.000318000,7861635.796,0
opencl-kernel,Apple M2,GPU,w32-il,p2048,2048,ADDMOD,2500,0.000299000,8361203.828,0
opencl-e2e,Apple M2,GPU,w32-il,p2048,2048,ADDMOD,2500,0.000387000,6459947.712,0
opencl-kernel,Apple M2,GPU,w32-il64,p2048,2048,ADDMOD,2500,0.000321000,7788162.005,0
opencl-e2e,Apple M2,GPU,w32-il64,p2048,2048,ADDMOD,2500,0.000418000,5980861.510,0
library,Apple M2,host-cpu,gmp-1t,p2048,2048,SUBTRACTMOD,2500,0.000097000,25773200.684,0
library,Apple M2,host-cpu,gmp-nt,p2048,2048,SUBTRACTMOD,2500,0.000075000,33333338.218,0
library,Apple M2,host-cpu,openssl-nt,p2048,2048,SUBTRACTMOD,2500,0.000232000,10775861.125,0
opencl-kernel,Apple M2,GPU,w8,p2048,2048,SUBTRACTMOD,2500,0.000414000,6038647.209,0
opencl-e2e,Apple M2,GPU,w8,p2048,2048,SUBTRACTMOD,2500,0.000516000,4844960.966,0
opencl-kernel,Apple M2,GPU,w16,p2048,2048,SUBTRACTMOD,2500,0.000324000,7716050.236,0
opencl-e2e,Apple M2,GPU,w16,p2048,2048,SUBTRACTMOD,2500,0.000409000,6112469.869,0
opencl-kernel,Apple M2,GPU,w32,p2048,2048,SUBTRACTMOD,2500,0.000430000,5813953.291,0
opencl-e2e,Apple M2,GPU,w32,p2048,2048,SUBTRACTMOD,2500,0.000533000,4690431.661,0
opencl-kernel,Apple M2,GPU,w32-opt,p2048,2048,SUBTRACTMOD,2500,0.000226000,11061945.295,0
opencl-e2e,Apple M2,GPU,w32-opt,p2048,2048,SUBTRACTMOD,2500,0.000303000,8250825.657,0
opencl-kernel,Apple M2,GPU,w32-o64,p2048,2048,SUBTRACTMOD,2500,0.000245000,10204081.512,0
opencl-e2e,Apple M2,GPU,w32-o64,p2048,2048,SUBTRACTMOD,2500,0.000329000,7598784.143,0
opencl-kernel,Apple M2,GPU,w32-il,p2048,2048,SUBTRACTMOD,2500,0.000270000,9259258.620,0
opencl-e2e,Apple M2,GPU,w32-il,p2048,2048,SUBTRACTMOD,2500,0.000358000,6983240.414,0
opencl-kernel,Apple M2,GPU,w32-il64,p2048,2048,SUBTRACTMOD,2500,0.000305000,8196721.470,0
opencl-e2e,Apple M2,GPU,w32-il64,p2048,2048,SUBTRACTMOD,2500,0.000387000,6459948.684,0
library,Apple M2,host-cpu,gmp-1t,p2048,2048,MULTIPLYOPERANDSCANNING,2500,0.000930000,2688171.932,0
library,Apple M2,host-cpu,gmp-nt,p2048,2048,MULTIPLYOPERANDSCANNING,2500,0.000319000,7836989.551,0
library,Apple M2,host-cpu,openssl-nt,p2048,2048,MULTIPLYOPERANDSCANNING,2500,0.000505000,4950495.014,0
opencl-kernel,Apple M2,GPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,2500,0.055467000,45071.845,0
opencl-e2e,Apple M2,GPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,2500,0.056053000,44600.646,0
opencl-kernel,Apple M2,GPU,w16,p2048,2048,MULTIPLYOPERANDSCANNING,2500,0.013996000,178622.464,0
opencl-e2e,Apple M2,GPU,w16,p2048,2048,MULTIPLYOPERANDSCANNING,2500,0.014348000,174240.312,0
opencl-kernel,Apple M2,GPU,w32,p2048,2048,MULTIPLYOPERANDSCANNING,2500,0.003691000,677323.220,0
opencl-e2e,Apple M2,GPU,w32,p2048,2048,MULTIPLYOPERANDSCANNING,2500,0.003874000,645327.829,0
opencl-kernel,Apple M2,GPU,w32-opt,p2048,2048,MULTIPLYOPERANDSCANNING,2500,0.000624000,4006410.021,0
opencl-e2e,Apple M2,GPU,w32-opt,p2048,2048,MULTIPLYOPERANDSCANNING,2500,0.000735000,3401360.504,0
opencl-kernel,Apple M2,GPU,w32-o64,p2048,2048,MULTIPLYOPERANDSCANNING,2500,0.000685000,3649634.848,0
opencl-e2e,Apple M2,GPU,w32-o64,p2048,2048,MULTIPLYOPERANDSCANNING,2500,0.000772000,3238341.971,0
opencl-kernel,Apple M2,GPU,w32-il,p2048,2048,MULTIPLYOPERANDSCANNING,2500,0.000735000,3401360.504,0
opencl-e2e,Apple M2,GPU,w32-il,p2048,2048,MULTIPLYOPERANDSCANNING,2500,0.000821000,3045066.905,0
opencl-kernel,Apple M2,GPU,w32-il64,p2048,2048,MULTIPLYOPERANDSCANNING,2500,0.000689000,3628446.984,0
opencl-e2e,Apple M2,GPU,w32-il64,p2048,2048,MULTIPLYOPERANDSCANNING,2500,0.000797000,3136763.030,0
library,Apple M2,host-cpu,gmp-1t,p2048,2048,MULTIPLYPRODUCTSCANNING,2500,0.000923000,2708559.090,0
library,Apple M2,host-cpu,gmp-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,2500,0.000321000,7788162.005,0
library,Apple M2,host-cpu,openssl-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,2500,0.000499000,5010019.878,0
opencl-kernel,Apple M2,GPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,2500,0.034099000,73315.933,0
opencl-e2e,Apple M2,GPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,2500,0.034588000,72279.403,0
opencl-kernel,Apple M2,GPU,w16,p2048,2048,MULTIPLYPRODUCTSCANNING,2500,0.008718000,286763.020,0
opencl-e2e,Apple M2,GPU,w16,p2048,2048,MULTIPLYPRODUCTSCANNING,2500,0.008904000,280772.685,0
opencl-kernel,Apple M2,GPU,w32,p2048,2048,MULTIPLYPRODUCTSCANNING,2500,0.002355000,1061571.123,0
opencl-e2e,Apple M2,GPU,w32,p2048,2048,MULTIPLYPRODUCTSCANNING,2500,0.002503000,998801.428,0
opencl-kernel,Apple M2,GPU,w32-opt,p2048,2048,MULTIPLYPRODUCTSCANNING,2500,0.002357000,1060670.335,0
opencl-e2e,Apple M2,GPU,w32-opt,p2048,2048,MULTIPLYPRODUCTSCANNING,2500,0.002498000,1000800.646,0
opencl-kernel,Apple M2,GPU,w32-o64,p2048,2048,MULTIPLYPRODUCTSCANNING,2500,0.000796000,3140703.488,0
opencl-e2e,Apple M2,GPU,w32-o64,p2048,2048,MULTIPLYPRODUCTSCANNING,2500,0.000909000,2750274.867,0
opencl-kernel,Apple M2,GPU,w32-il,p2048,2048,MULTIPLYPRODUCTSCANNING,2500,0.000985000,2538070.938,0
opencl-e2e,Apple M2,GPU,w32-il,p2048,2048,MULTIPLYPRODUCTSCANNING,2500,0.001085000,2304147.473,0
opencl-kernel,Apple M2,GPU,w32-il64,p2048,2048,MULTIPLYPRODUCTSCANNING,2500,0.000991000,2522704.393,0
opencl-e2e,Apple M2,GPU,w32-il64,p2048,2048,MULTIPLYPRODUCTSCANNING,2500,0.001088000,2297794.077,0
library,Apple M2,host-cpu,gmp-1t,p2048,2048,MONTGOMERYMULTIPLICATION,2500,0.006778000,368840.367,0
library,Apple M2,host-cpu,gmp-nt,p2048,2048,MONTGOMERYMULTIPLICATION,2500,0.001771000,1411631.861,0
library,Apple M2,host-cpu,openssl-nt,p2048,2048,MONTGOMERYMULTIPLICATION,2500,0.000694000,3602305.540,0
opencl-kernel,Apple M2,GPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,2500,0.045179000,55335.443,0
opencl-e2e,Apple M2,GPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,2500,0.045158000,55361.176,0
opencl-kernel,Apple M2,GPU,w16,p2048,2048,MONTGOMERYMULTIPLICATION,2500,0.007597000,329077.265,0
opencl-e2e,Apple M2,GPU,w16,p2048,2048,MONTGOMERYMULTIPLICATION,2500,0.007766000,321916.044,0
opencl-kernel,Apple M2,GPU,w32,p2048,2048,MONTGOMERYMULTIPLICATION,2500,0.002481000,1007658.189,0
opencl-e2e,Apple M2,GPU,w32,p2048,2048,MONTGOMERYMULTIPLICATION,2500,0.002534000,986582.490,0
opencl-kernel,Apple M2,GPU,w32-opt,p2048,2048,MONTGOMERYMULTIPLICATION,2500,0.001625000,1538461.507,0
opencl-e2e,Apple M2,GPU,w32-opt,p2048,2048,MONTGOMERYMULTIPLICATION,2500,0.001703000,1467997.662,0
opencl-kernel,Apple M2,GPU,w32-o64,p2048,2048,MONTGOMERYMULTIPLICATION,2500,0.001793000,1394311.185,0
opencl-e2e,Apple M2,GPU,w32-o64,p2048,2048,MONTGOMERYMULTIPLICATION,2500,0.001870000,1336898.411,0
opencl-kernel,Apple M2,GPU,w32-il,p2048,2048,MONTGOMERYMULTIPLICATION,2500,0.001566000,1596424.019,0
opencl-e2e,Apple M2,GPU,w32-il,p2048,2048,MONTGOMERYMULTIPLICATION,2500,0.001671000,1496110.083,0
opencl-kernel,Apple M2,GPU,w32-il64,p2048,2048,MONTGOMERYMULTIPLICATION,2500,0.001537000,1626545.214,0
opencl-e2e,Apple M2,GPU,w32-il64,p2048,2048,MONTGOMERYMULTIPLICATION,2500,0.001620000,1543209.881,0
library,Apple M2,host-cpu,gmp-1t,p2048,2048,COMPARE,2500,0.000009000,277779255.714,0
library,Apple M2,host-cpu,gmp-nt,p2048,2048,COMPARE,2500,0.000051000,49019626.215,0
library,Apple M2,host-cpu,openssl-nt,p2048,2048,COMPARE,2500,0.000068000,36764709.171,0
opencl-kernel,Apple M2,GPU,w8,p2048,2048,COMPARE,2500,0.000235000,10638298.553,0
opencl-e2e,Apple M2,GPU,w8,p2048,2048,COMPARE,2500,0.000334000,7485029.037,0
opencl-kernel,Apple M2,GPU,w16,p2048,2048,COMPARE,2500,0.000219000,11415526.180,0
opencl-e2e,Apple M2,GPU,w16,p2048,2048,COMPARE,2500,0.000304000,8223683.001,0
opencl-kernel,Apple M2,GPU,w32-opt,p2048,2048,COMPARE,2500,0.000216000,11574073.275,0
opencl-e2e,Apple M2,GPU,w32-opt,p2048,2048,COMPARE,2500,0.000298000,8389261.773,0
opencl-kernel,Apple M2,GPU,w32-o64,p2048,2048,COMPARE,2500,0.000225000,11111112.739,0
opencl-e2e,Apple M2,GPU,w32-o64,p2048,2048,COMPARE,2500,0.000319000,7836989.551,0
opencl-kernel,Apple M2,GPU,w32-il,p2048,2048,COMPARE,2500,0.000210000,11904763.649,0
opencl-e2e,Apple M2,GPU,w32-il,p2048,2048,COMPARE,2500,0.000289000,8650519.369,0
opencl-kernel,Apple M2,GPU,w32-il64,p2048,2048,COMPARE,2500,0.000215000,11627906.582,0
opencl-e2e,Apple M2,GPU,w32-il64,p2048,2048,COMPARE,2500,0.000300000,8333332.938,0
library,Apple M2,host-cpu,gmp-1t,p2048,2048,REDUCE,312,0.000007000,44571459.811,0
library,Apple M2,host-cpu,gmp-nt,p2048,2048,REDUCE,312,0.000062000,5032257.542,0
library,Apple M2,host-cpu,openssl-nt,p2048,2048,REDUCE,312,0.000083000,3759036.520,0
opencl-kernel,Apple M2,GPU,w8,p2048,2048,REDUCE,312,0.013889000,22463.820,0
opencl-e2e,Apple M2,GPU,w8,p2048,2048,REDUCE,312,0.013970000,22333.572,0
opencl-kernel,Apple M2,GPU,w16,p2048,2048,REDUCE,312,0.006807000,45835.169,0
opencl-e2e,Apple M2,GPU,w16,p2048,2048,REDUCE,312,0.006864000,45454.545,0
opencl-kernel,Apple M2,GPU,w32-opt,p2048,2048,REDUCE,312,0.004397000,70957.471,0
opencl-e2e,Apple M2,GPU,w32-opt,p2048,2048,REDUCE,312,0.004423000,70540.357,0
opencl-kernel,Apple M2,GPU,w32-o64,p2048,2048,REDUCE,312,0.004414000,70684.187,0
opencl-e2e,Apple M2,GPU,w32-o64,p2048,2048,REDUCE,312,0.004449000,70128.119,0
opencl-kernel,Apple M2,GPU,w32-il,p2048,2048,REDUCE,312,0.004404000,70844.687,0
opencl-e2e,Apple M2,GPU,w32-il,p2048,2048,REDUCE,312,0.004440000,70270.270,0
opencl-kernel,Apple M2,GPU,w32-il64,p2048,2048,REDUCE,312,0.005567000,56044.548,0
opencl-e2e,Apple M2,GPU,w32-il64,p2048,2048,REDUCE,312,0.005575000,55964.125,0
library,Apple M2,host-cpu,gmp-1t,p2048,2048,MODMUL,156,0.000174000,896551.796,0
library,Apple M2,host-cpu,gmp-nt,p2048,2048,MODMUL,156,0.000095000,1642105.839,0
library,Apple M2,host-cpu,openssl-nt,p2048,2048,MODMUL,156,0.000228000,684210.382,0
opencl-kernel,Apple M2,GPU,w8,p2048,2048,MODMUL,156,0.062022000,2515.237,0
opencl-e2e,Apple M2,GPU,w8,p2048,2048,MODMUL,156,0.061698000,2528.445,0
opencl-kernel,Apple M2,GPU,w16,p2048,2048,MODMUL,156,0.040263000,3874.525,0
opencl-e2e,Apple M2,GPU,w16,p2048,2048,MODMUL,156,0.040337000,3867.417,0
opencl-kernel,Apple M2,GPU,w32-opt,p2048,2048,MODMUL,156,0.025263000,6175.039,0
opencl-e2e,Apple M2,GPU,w32-opt,p2048,2048,MODMUL,156,0.025315000,6162.354,0
opencl-kernel,Apple M2,GPU,w32-o64,p2048,2048,MODMUL,156,0.024122000,6467.125,0
opencl-e2e,Apple M2,GPU,w32-o64,p2048,2048,MODMUL,156,0.024177000,6452.413,0
opencl-kernel,Apple M2,GPU,w32-il,p2048,2048,MODMUL,156,0.027489000,5674.997,0
opencl-e2e,Apple M2,GPU,w32-il,p2048,2048,MODMUL,156,0.027508000,5671.078,0
opencl-kernel,Apple M2,GPU,w32-il64,p2048,2048,MODMUL,156,0.024391000,6395.802,0
opencl-e2e,Apple M2,GPU,w32-il64,p2048,2048,MODMUL,156,0.024416000,6389.253,0
library,Apple M2,host-cpu,gmp-1t,p2048,2048,MODEXP,64,0.089052000,718.681,0
library,Apple M2,host-cpu,gmp-nt,p2048,2048,MODEXP,64,0.020247000,3160.962,0
library,Apple M2,host-cpu,openssl-nt,p2048,2048,MODEXP,64,0.018702000,3422.094,0
opencl-kernel,Apple M2,GPU,w8,p2048,2048,MODEXP,64,20.476306000,3.126,0
opencl-e2e,Apple M2,GPU,w8,p2048,2048,MODEXP,64,20.476034000,3.126,0
opencl-kernel,Apple M2,GPU,w16,p2048,2048,MODEXP,64,12.009537000,5.329,0
opencl-e2e,Apple M2,GPU,w16,p2048,2048,MODEXP,64,12.010037000,5.329,0
opencl-kernel,Apple M2,GPU,w32-opt,p2048,2048,MODEXP,64,2.012495000,31.801,0
opencl-e2e,Apple M2,GPU,w32-opt,p2048,2048,MODEXP,64,2.012901000,31.795,0
opencl-kernel,Apple M2,GPU,w32-o64,p2048,2048,MODEXP,64,1.947243000,32.867,0
opencl-e2e,Apple M2,GPU,w32-o64,p2048,2048,MODEXP,64,1.947201000,32.868,0
opencl-kernel,Apple M2,GPU,w32-il,p2048,2048,MODEXP,64,1.999607000,32.006,0
opencl-e2e,Apple M2,GPU,w32-il,p2048,2048,MODEXP,64,2.000955000,31.985,0
opencl-kernel,Apple M2,GPU,w32-il64,p2048,2048,MODEXP,64,1.955654000,32.726,0
opencl-e2e,Apple M2,GPU,w32-il64,p2048,2048,MODEXP,64,1.958319000,32.681,0
library,Apple M2,host-cpu,gmp-1t,p2048,2048,EXPONENTIATION,64,0.010001000,6399.360,0
library,Apple M2,host-cpu,gmp-nt,p2048,2048,EXPONENTIATION,64,0.002750000,23272.727,0
library,Apple M2,host-cpu,openssl-nt,p2048,2048,EXPONENTIATION,64,0.052329000,1223.031,0
opencl-kernel,Apple M2,GPU,w8,p2048,2048,EXPONENTIATION,64,25.736167000,2.487,0
opencl-e2e,Apple M2,GPU,w8,p2048,2048,EXPONENTIATION,64,25.693380000,2.491,0
opencl-kernel,Apple M2,GPU,w16,p2048,2048,EXPONENTIATION,64,6.955807000,9.201,0
opencl-e2e,Apple M2,GPU,w16,p2048,2048,EXPONENTIATION,64,6.955732000,9.201,0
opencl-kernel,Apple M2,GPU,w32-opt,p2048,2048,EXPONENTIATION,64,1.859638000,34.415,0
opencl-e2e,Apple M2,GPU,w32-opt,p2048,2048,EXPONENTIATION,64,1.859856000,34.411,0
opencl-kernel,Apple M2,GPU,w32-o64,p2048,2048,EXPONENTIATION,64,1.840396000,34.775,0
opencl-e2e,Apple M2,GPU,w32-o64,p2048,2048,EXPONENTIATION,64,1.841140000,34.761,0
opencl-kernel,Apple M2,GPU,w32-il,p2048,2048,EXPONENTIATION,64,1.862383000,34.365,0
opencl-e2e,Apple M2,GPU,w32-il,p2048,2048,EXPONENTIATION,64,1.863766000,34.339,0
opencl-kernel,Apple M2,GPU,w32-il64,p2048,2048,EXPONENTIATION,64,1.851776000,34.561,0
opencl-e2e,Apple M2,GPU,w32-il64,p2048,2048,EXPONENTIATION,64,1.851749000,34.562,0
library,Apple M2,host-cpu,gmp-1t,p2048,2048,DIVIDE,312,0.000011000,28363570.478,0
library,Apple M2,host-cpu,gmp-nt,p2048,2048,DIVIDE,312,0.000053000,5886792.884,0
library,Apple M2,host-cpu,openssl-nt,p2048,2048,DIVIDE,312,0.000074000,4216217.276,0
opencl-kernel,Apple M2,GPU,w8,p2048,2048,DIVIDE,312,0.106022000,2942.785,0
opencl-e2e,Apple M2,GPU,w8,p2048,2048,DIVIDE,312,0.106024000,2942.730,0
opencl-kernel,Apple M2,GPU,w16,p2048,2048,DIVIDE,312,0.061149000,5102.291,0
opencl-e2e,Apple M2,GPU,w16,p2048,2048,DIVIDE,312,0.061219000,5096.457,0
opencl-kernel,Apple M2,GPU,w32-opt,p2048,2048,DIVIDE,312,0.033923000,9197.300,0
opencl-e2e,Apple M2,GPU,w32-opt,p2048,2048,DIVIDE,312,0.033895000,9204.897,0
opencl-kernel,Apple M2,GPU,w32-o64,p2048,2048,DIVIDE,312,0.034409000,9067.395,0
opencl-e2e,Apple M2,GPU,w32-o64,p2048,2048,DIVIDE,312,0.034436000,9060.286,0
opencl-kernel,Apple M2,GPU,w32-il,p2048,2048,DIVIDE,312,0.033364000,9351.397,0
opencl-e2e,Apple M2,GPU,w32-il,p2048,2048,DIVIDE,312,0.033247000,9384.305,0
opencl-kernel,Apple M2,GPU,w32-il64,p2048,2048,DIVIDE,312,0.034103000,9148.755,0
opencl-e2e,Apple M2,GPU,w32-il64,p2048,2048,DIVIDE,312,0.034149000,9136.432,0
library,Apple M2,host-cpu,gmp-1t,p2048,2048,ISQRT,78,0.000026000,2999997.305,0
library,Apple M2,host-cpu,gmp-nt,p2048,2048,ISQRT,78,0.000055000,1418181.526,0
opencl-kernel,Apple M2,GPU,w8,p2048,2048,ISQRT,78,0.810953000,96.183,0
opencl-e2e,Apple M2,GPU,w8,p2048,2048,ISQRT,78,0.811366000,96.134,0
opencl-kernel,Apple M2,GPU,w16,p2048,2048,ISQRT,78,0.452774000,172.271,0
opencl-e2e,Apple M2,GPU,w16,p2048,2048,ISQRT,78,0.452858000,172.239,0
opencl-kernel,Apple M2,GPU,w32-opt,p2048,2048,ISQRT,78,0.289073000,269.828,0
opencl-e2e,Apple M2,GPU,w32-opt,p2048,2048,ISQRT,78,0.289094000,269.808,0
opencl-kernel,Apple M2,GPU,w32-o64,p2048,2048,ISQRT,78,0.239578000,325.572,0
opencl-e2e,Apple M2,GPU,w32-o64,p2048,2048,ISQRT,78,0.239603000,325.538,0
opencl-kernel,Apple M2,GPU,w32-il,p2048,2048,ISQRT,78,0.271338000,287.464,0
opencl-e2e,Apple M2,GPU,w32-il,p2048,2048,ISQRT,78,0.271485000,287.309,0
opencl-kernel,Apple M2,GPU,w32-il64,p2048,2048,ISQRT,78,0.244328000,319.243,0
opencl-e2e,Apple M2,GPU,w32-il64,p2048,2048,ISQRT,78,0.244791000,318.639,0
library,Apple M2,host-cpu,gmp-1t,p2048,2048,MODMUL_R2,2500,0.003044000,821287.776,0
library,Apple M2,host-cpu,gmp-nt,p2048,2048,MODMUL_R2,2500,0.000855000,2923976.639,0
library,Apple M2,host-cpu,openssl-nt,p2048,2048,MODMUL_R2,2500,0.002651000,943040.366,0
opencl-kernel,Apple M2,GPU,w8,p2048,2048,MODMUL_R2,2500,0.016142000,154875.480,0
opencl-e2e,Apple M2,GPU,w8,p2048,2048,MODMUL_R2,2500,0.016143000,154865.886,0
opencl-kernel,Apple M2,GPU,w16,p2048,2048,MODMUL_R2,2500,0.008081000,309367.652,0
opencl-e2e,Apple M2,GPU,w16,p2048,2048,MODMUL_R2,2500,0.008171000,305960.103,0
opencl-kernel,Apple M2,GPU,w32-opt,p2048,2048,MODMUL_R2,2500,0.002468000,1012965.968,0
opencl-e2e,Apple M2,GPU,w32-opt,p2048,2048,MODMUL_R2,2500,0.002563000,975419.424,0
opencl-kernel,Apple M2,GPU,w32-o64,p2048,2048,MODMUL_R2,2500,0.002109000,1185395.926,0
opencl-e2e,Apple M2,GPU,w32-o64,p2048,2048,MODMUL_R2,2500,0.002179000,1147315.281,0
opencl-kernel,Apple M2,GPU,w32-il,p2048,2048,MODMUL_R2,2500,0.001906000,1311647.456,0
opencl-e2e,Apple M2,GPU,w32-il,p2048,2048,MODMUL_R2,2500,0.002008000,1245019.896,0
opencl-kernel,Apple M2,GPU,w32-il64,p2048,2048,MODMUL_R2,2500,0.001903000,1313715.189,0
opencl-e2e,Apple M2,GPU,w32-il64,p2048,2048,MODMUL_R2,2500,0.001974000,1266464.061,0
```
