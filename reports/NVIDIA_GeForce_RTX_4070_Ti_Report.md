# MPA-OpenCL benchmark report - NVIDIA GeForce RTX 4070 Ti


> **Note.** Two column groups have been removed from this report: the CGBN
> reference column, which was invalid, and the multi-threaded GMP and OpenSSL
> baselines, which predate the 2026-09-12 timing fix and were understated.
> See `reports/README.md`. The single-threaded GMP column, the OpenCL-on-CPU
> rows and every MPA measurement are unaffected, and every configuration was
> verified word-for-word against GMP before it was timed.


## 1. System under test

1 OpenCL device(s) exercised with the identical kernels and operands.

### Device 0 - NVIDIA GeForce RTX 4070 Ti (GPU)

| Property | Value |
|---|---|
| Model | NVIDIA GeForce RTX 4070 Ti |
| Type | GPU |
| Vendor | NVIDIA Corporation |
| Device memory | 11.61 GiB |
| Max single allocation | 2.90 GiB |
| Local memory | 48 KiB |
| Global cache | 1680 KiB |
| Compute units | 60 |
| Max clock | 2715 MHz |
| Max work-group size | 1024 |
| OpenCL version | OpenCL 3.0 CUDA |
| Driver | 570.195.03 |

### Host

| Property | Value |
|---|---|
| CPU | Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz |
| Logical cores | 72 |
| OpenMP threads used | 72 |
| RAM | 125.7 GB |
| OS | Ubuntu 24.04.4 LTS |
| Kernel | 6.8.0-87-generic |
| Arch | x86_64 |
| GMP | 6.3.0 |
| OpenSSL | OpenSSL 3.0.13 30 Jan 2024 |

## 2. Method

- Base workload 200000 items, scaled down per operator by its cost weight and by modulus size; the exact count is in every row.
- 5 timed repetitions, **minimum** reported. Two untimed warm-up launches precede them.
- `kernel` times `clEnqueueNDRangeKernel` + `clFinish` only. `e2e` adds the host->device operand writes and the device->host result read.
- Every OpenCL device runs the same kernels on the same operands, so GPU and CPU-OpenCL columns are directly comparable.
- CPU library baselines (GMP, OpenSSL) run those same operands, with temporaries preallocated outside the timed region, so the figure is the arithmetic and not marshalling. The generator is reseeded per modulus and operation so every backend sees identical inputs.
- OpenSSL rows time the nearest BN primitive, which is not always semantically identical (its Montgomery routine expects Montgomery-domain inputs); they measure comparable work, not identical results. Correctness is judged against GMP only.
- Every device cell is checked word-for-word against GMP before it is timed. A cell that mismatches is reported and excluded from the speedup tables.
- Total wall time 2228.2 s.

## 3. Correctness

| Device | Kernel | Configs run | Passed | Mismatched | Build/launch failed |
|---|---|---|---|---|---|
| [0] GPU | `mpaKernels_8bits.cl` (w8) | 75 | 75 | 0 | 0 |
| [0] GPU | `mpaKernel_16bits.cl` (w16) | 75 | 75 | 0 | 0 |
| [0] GPU | `mpaKernel_32bits.cl` (w32) | 35 | 35 | 0 | 0 |
| [0] GPU | `mpaKernel_32bits_opt.cl` (w32-opt) | 75 | 75 | 0 | 0 |
| [0] GPU | `mpaKernel_32bits_opt.cl` (w32-o64) | 75 | 75 | 0 | 0 |

**All configurations correct** - 335 configurations, 0 problems.

## 4. Throughput per device

Operations per second, higher is better. Kernel-only timings.

### Device 0 - NVIDIA GeForce RTX 4070 Ti (GPU)

#### secp256k1 (256-bit)

| Operation | items | w8 | w16 | w32 | w32-opt | w32-o64 | GMP 1T |
|---|---|---|---|---|---|---|---|
| ADD | 200000 | 1.22 G | 2.19 G | 3.87 G | 3.95 G | 4.05 G | 36.98 M |
| SUBTRACT | 200000 | 1.23 G | 2.27 G | 4.09 G | 3.82 G | 4.09 G | 39.72 M |
| ADDMOD | 200000 | 777.67 M | 1.54 G | 3.25 G | 5.99 G | 6.50 G | 12.26 M |
| SUBTRACTMOD | 200000 | 778.68 M | 1.54 G | 3.19 G | 5.86 G | 7.09 G | 16.81 M |
| MULTIPLYOPERANDSCANNING | 200000 | 28.89 M | 106.22 M | 379.25 M | 1.94 G | 1.96 G | 31.21 M |
| MULTIPLYPRODUCTSCANNING | 200000 | 255.91 M | 771.85 M | 1.84 G | 1.86 G | 1.99 G | 29.89 M |
| MONTGOMERYMULTIPLICATION | 200000 | 556.43 M | 2.30 G | 8.42 G | 5.70 G | 5.88 G | 3.96 M |
| COMPARE | 200000 | 1.30 G | 2.20 G | - | 4.22 G | 4.17 G | 94.06 M |
| REDUCE | 25000 | 223.96 M | 341.67 M | - | 867.36 M | 890.12 M | 33.96 M |
| MODMUL | 20000 | 73.80 M | 121.78 M | - | 254.12 M | 317.14 M | 7.13 M |
| MODEXP | 20000 | 1.65 M | 8.44 M | - | 12.96 M | 22.34 M | 52.49 k |
| EXPONENTIATION | 20000 | 886.14 k | 3.23 M | - | 51.53 M | 58.11 M | 159.16 k |
| DIVIDE | 25000 | 106.45 M | 121.47 M | - | 334.66 M | 348.56 M | 16.78 M |
| ISQRT | 20000 | 7.94 M | 10.21 M | - | 42.83 M | 47.17 M | 7.75 M |
| MODMUL_R2 | 200000 | 481.46 M | 2.37 G | - | 3.77 G | 3.89 G | 7.17 M |

#### rsa256(composite) (256-bit)

| Operation | items | w8 | w16 | w32 | w32-opt | w32-o64 | GMP 1T |
|---|---|---|---|---|---|---|---|
| ADD | 200000 | 1.27 G | 2.35 G | 4.01 G | 4.06 G | 3.91 G | 36.48 M |
| SUBTRACT | 200000 | 1.26 G | 2.38 G | 3.92 G | 4.00 G | 4.08 G | 46.10 M |
| ADDMOD | 200000 | 866.19 M | 1.70 G | 3.46 G | 6.44 G | 6.38 G | 15.76 M |
| SUBTRACTMOD | 200000 | 776.17 M | 1.48 G | 3.10 G | 6.72 G | 6.51 G | 16.68 M |
| MULTIPLYOPERANDSCANNING | 200000 | 28.70 M | 105.64 M | 381.05 M | 1.99 G | 1.96 G | 32.18 M |
| MULTIPLYPRODUCTSCANNING | 200000 | 255.29 M | 772.09 M | 1.83 G | 1.90 G | 1.97 G | 31.31 M |
| MONTGOMERYMULTIPLICATION | 200000 | 554.33 M | 2.30 G | 8.08 G | 6.70 G | 5.85 G | 4.16 M |
| COMPARE | 200000 | 1.26 G | 2.20 G | - | 4.20 G | 3.97 G | 90.85 M |
| REDUCE | 25000 | 219.36 M | 323.57 M | - | 864.96 M | 863.14 M | 19.75 M |
| MODMUL | 20000 | 71.73 M | 120.52 M | - | 271.77 M | 313.16 M | 7.42 M |
| MODEXP | 20000 | 1.66 M | 8.45 M | - | 13.01 M | 22.38 M | 72.78 k |
| EXPONENTIATION | 20000 | 888.08 k | 3.22 M | - | 51.78 M | 58.47 M | 164.40 k |
| DIVIDE | 25000 | 106.41 M | 120.34 M | - | 322.54 M | 333.08 M | 15.03 M |
| ISQRT | 20000 | 7.95 M | 9.96 M | - | 41.23 M | 43.48 M | 7.89 M |
| MODMUL_R2 | 200000 | 482.51 M | 2.34 G | - | 3.92 G | 3.79 G | 6.88 M |

#### brainpoolP512r1 (512-bit)

| Operation | items | w8 | w16 | w32 | w32-opt | w32-o64 | GMP 1T |
|---|---|---|---|---|---|---|---|
| ADD | 100000 | 531.36 M | 1.02 G | 1.49 G | 1.59 G | 1.59 G | 34.82 M |
| SUBTRACT | 100000 | 534.50 M | 1.02 G | 1.48 G | 1.65 G | 1.59 G | 43.42 M |
| ADDMOD | 100000 | 380.27 M | 729.83 M | 1.25 G | 1.53 G | 1.53 G | 14.11 M |
| SUBTRACTMOD | 100000 | 320.50 M | 645.92 M | 1.09 G | 1.50 G | 1.51 G | 15.91 M |
| MULTIPLYOPERANDSCANNING | 100000 | 6.38 M | 21.99 M | 67.11 M | 448.26 M | 445.75 M | 14.19 M |
| MULTIPLYPRODUCTSCANNING | 100000 | 35.98 M | 132.03 M | 442.91 M | 438.85 M | 399.81 M | 14.07 M |
| MONTGOMERYMULTIPLICATION | 100000 | 141.52 M | 443.01 M | 2.11 G | 1.56 G | 1.70 G | 1.58 M |
| COMPARE | 100000 | 511.19 M | 966.97 M | - | 1.75 G | 1.79 G | 90.01 M |
| REDUCE | 20000 | 62.73 M | 82.62 M | - | 294.49 M | 293.48 M | 17.08 M |
| MODMUL | 20000 | 22.55 M | 33.21 M | - | 84.08 M | 104.30 M | 3.52 M |
| MODEXP | 20000 | 113.93 k | 1.24 M | - | 1.61 M | 3.36 M | 13.17 k |
| EXPONENTIATION | 20000 | 109.12 k | 425.08 k | - | 1.35 M | 1.38 M | 56.97 k |
| DIVIDE | 20000 | 25.65 M | 25.63 M | - | 77.50 M | 79.12 M | 7.11 M |
| ISQRT | 20000 | 1.71 M | 1.67 M | - | 8.26 M | 8.82 M | 4.38 M |
| MODMUL_R2 | 100000 | 95.33 M | 481.53 M | - | 936.35 M | 1.37 G | 3.05 M |

#### p1024 (1024-bit)

| Operation | items | w8 | w16 | w32 | w32-opt | w32-o64 | GMP 1T |
|---|---|---|---|---|---|---|---|
| ADD | 50000 | 182.76 M | 348.24 M | 448.29 M | 631.63 M | 634.05 M | 15.22 M |
| SUBTRACT | 50000 | 185.63 M | 348.55 M | 454.71 M | 639.30 M | 642.07 M | 35.86 M |
| ADDMOD | 50000 | 117.85 M | 235.63 M | 346.22 M | 486.24 M | 484.51 M | 7.40 M |
| SUBTRACTMOD | 50000 | 118.47 M | 227.66 M | 346.21 M | 483.75 M | 482.55 M | 10.02 M |
| MULTIPLYOPERANDSCANNING | 50000 | 1.39 M | 5.15 M | 13.32 M | 160.85 M | 163.02 M | 4.29 M |
| MULTIPLYPRODUCTSCANNING | 50000 | 4.17 M | 16.33 M | 60.94 M | 61.39 M | 61.33 M | 4.22 M |
| MONTGOMERYMULTIPLICATION | 50000 | 22.99 M | 115.20 M | 380.52 M | 270.39 M | 330.68 M | 631.43 k |
| COMPARE | 50000 | 192.15 M | 363.42 M | - | 654.07 M | 648.60 M | 53.25 M |
| REDUCE | 20000 | 12.48 M | 22.98 M | - | 86.24 M | 86.68 M | 26.89 M |
| MODMUL | 20000 | 3.13 M | 8.71 M | - | 19.19 M | 25.51 M | 1.38 M |
| MODEXP | 20000 | 14.16 k | 92.69 k | - | 210.08 k | 330.36 k | 1.95 k |
| EXPONENTIATION | 20000 | 14.64 k | 54.87 k | - | 213.00 k | 216.04 k | 14.36 k |
| DIVIDE | 20000 | 1.53 M | 3.47 M | - | 17.93 M | 18.31 M | 6.42 M |
| ISQRT | 20000 | 111.91 k | 235.24 k | - | 1.18 M | 1.27 M | 2.40 M |
| MODMUL_R2 | 50000 | 14.64 M | 100.68 M | - | 183.43 M | 225.82 M | 1.33 M |

#### p2048 (2048-bit)

| Operation | items | w8 | w16 | w32 | w32-opt | w32-o64 | GMP 1T |
|---|---|---|---|---|---|---|---|
| ADD | 25000 | 85.68 M | 167.99 M | 275.25 M | 286.59 M | 280.93 M | 21.52 M |
| SUBTRACT | 25000 | 83.14 M | 174.03 M | 286.62 M | 294.16 M | 287.10 M | 24.64 M |
| ADDMOD | 25000 | 59.22 M | 122.49 M | 240.38 M | 233.17 M | 249.85 M | 7.02 M |
| SUBTRACTMOD | 25000 | 55.14 M | 114.60 M | 199.87 M | 236.63 M | 245.56 M | 10.11 M |
| MULTIPLYOPERANDSCANNING | 25000 | 333.79 k | 1.27 M | 4.76 M | 55.06 M | 55.52 M | 1.32 M |
| MULTIPLYPRODUCTSCANNING | 25000 | 1.05 M | 4.11 M | 16.28 M | 15.51 M | 15.56 M | 1.26 M |
| MONTGOMERYMULTIPLICATION | 25000 | 953.72 k | 24.90 M | 125.10 M | 84.63 M | 100.37 M | 156.13 k |
| COMPARE | 25000 | 98.28 M | 179.64 M | - | 347.62 M | 349.31 M | 39.40 M |
| REDUCE | 20000 | 136.53 k | 5.11 M | - | 27.16 M | 25.84 M | 18.74 M |
| MODMUL | 20000 | 74.80 k | 1.94 M | - | 4.79 M | 6.08 M | 415.85 k |
| MODEXP | 20000 | 558.8 | 3.74 k | - | 24.15 k | 13.38 k | 276.2 |
| EXPONENTIATION | 20000 | 1.49 k | 7.28 k | - | 27.44 k | 27.95 k | 2.56 k |
| DIVIDE | 20000 | 31.77 k | 191.22 k | - | 1.86 M | 2.14 M | 12.14 M |
| ISQRT | 20000 | 2.40 k | 12.93 k | - | 476.55 k | 553.95 k | 1.35 M |
| MODMUL_R2 | 25000 | 1.22 M | 22.12 M | - | 51.18 M | 65.36 M | 385.01 k |

## 5. Head to head

Best OpenCL GPU result against best OpenCL CPU result and the CPU libraries.
Ratios above 1.00x mean the GPU is faster than that baseline.

### secp256k1 (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w32-o64 | 4.05 G | none | n/a | 36.98 M | n/a |
| SUBTRACT | w32-o64 | 4.09 G | none | n/a | 39.72 M | n/a |
| ADDMOD | w32-o64 | 6.50 G | none | n/a | 12.26 M | n/a |
| SUBTRACTMOD | w32-o64 | 7.09 G | none | n/a | 16.81 M | n/a |
| MULTIPLYOPERANDSCANNING | w32-o64 | 1.96 G | none | n/a | 31.21 M | n/a |
| MULTIPLYPRODUCTSCANNING | w32-o64 | 1.99 G | none | n/a | 29.89 M | n/a |
| MONTGOMERYMULTIPLICATION | w32 | 8.42 G | none | n/a | 3.96 M | n/a |
| COMPARE | w32-opt | 4.22 G | none | n/a | 94.06 M | n/a |
| REDUCE | w32-o64 | 890.12 M | none | n/a | 33.96 M | n/a |
| MODMUL | w32-o64 | 317.14 M | none | n/a | 7.13 M | n/a |
| MODEXP | w32-o64 | 22.34 M | none | n/a | 52.49 k | n/a |
| EXPONENTIATION | w32-o64 | 58.11 M | none | n/a | 159.16 k | n/a |
| DIVIDE | w32-o64 | 348.56 M | none | n/a | 16.78 M | n/a |
| ISQRT | w32-o64 | 47.17 M | none | n/a | 7.75 M | n/a |
| MODMUL_R2 | w32-o64 | 3.89 G | none | n/a | 7.17 M | n/a |

### rsa256(composite) (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w32-opt | 4.06 G | none | n/a | 36.48 M | n/a |
| SUBTRACT | w32-o64 | 4.08 G | none | n/a | 46.10 M | n/a |
| ADDMOD | w32-opt | 6.44 G | none | n/a | 15.76 M | n/a |
| SUBTRACTMOD | w32-opt | 6.72 G | none | n/a | 16.68 M | n/a |
| MULTIPLYOPERANDSCANNING | w32-opt | 1.99 G | none | n/a | 32.18 M | n/a |
| MULTIPLYPRODUCTSCANNING | w32-o64 | 1.97 G | none | n/a | 31.31 M | n/a |
| MONTGOMERYMULTIPLICATION | w32 | 8.08 G | none | n/a | 4.16 M | n/a |
| COMPARE | w32-opt | 4.20 G | none | n/a | 90.85 M | n/a |
| REDUCE | w32-opt | 864.96 M | none | n/a | 19.75 M | n/a |
| MODMUL | w32-o64 | 313.16 M | none | n/a | 7.42 M | n/a |
| MODEXP | w32-o64 | 22.38 M | none | n/a | 72.78 k | n/a |
| EXPONENTIATION | w32-o64 | 58.47 M | none | n/a | 164.40 k | n/a |
| DIVIDE | w32-o64 | 333.08 M | none | n/a | 15.03 M | n/a |
| ISQRT | w32-o64 | 43.48 M | none | n/a | 7.89 M | n/a |
| MODMUL_R2 | w32-opt | 3.92 G | none | n/a | 6.88 M | n/a |

### brainpoolP512r1 (512-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w32-o64 | 1.59 G | none | n/a | 34.82 M | n/a |
| SUBTRACT | w32-opt | 1.65 G | none | n/a | 43.42 M | n/a |
| ADDMOD | w32-opt | 1.53 G | none | n/a | 14.11 M | n/a |
| SUBTRACTMOD | w32-o64 | 1.51 G | none | n/a | 15.91 M | n/a |
| MULTIPLYOPERANDSCANNING | w32-opt | 448.26 M | none | n/a | 14.19 M | n/a |
| MULTIPLYPRODUCTSCANNING | w32 | 442.91 M | none | n/a | 14.07 M | n/a |
| MONTGOMERYMULTIPLICATION | w32 | 2.11 G | none | n/a | 1.58 M | n/a |
| COMPARE | w32-o64 | 1.79 G | none | n/a | 90.01 M | n/a |
| REDUCE | w32-opt | 294.49 M | none | n/a | 17.08 M | n/a |
| MODMUL | w32-o64 | 104.30 M | none | n/a | 3.52 M | n/a |
| MODEXP | w32-o64 | 3.36 M | none | n/a | 13.17 k | n/a |
| EXPONENTIATION | w32-o64 | 1.38 M | none | n/a | 56.97 k | n/a |
| DIVIDE | w32-o64 | 79.12 M | none | n/a | 7.11 M | n/a |
| ISQRT | w32-o64 | 8.82 M | none | n/a | 4.38 M | n/a |
| MODMUL_R2 | w32-o64 | 1.37 G | none | n/a | 3.05 M | n/a |

### p1024 (1024-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w32-o64 | 634.05 M | none | n/a | 15.22 M | n/a |
| SUBTRACT | w32-o64 | 642.07 M | none | n/a | 35.86 M | n/a |
| ADDMOD | w32-opt | 486.24 M | none | n/a | 7.40 M | n/a |
| SUBTRACTMOD | w32-opt | 483.75 M | none | n/a | 10.02 M | n/a |
| MULTIPLYOPERANDSCANNING | w32-o64 | 163.02 M | none | n/a | 4.29 M | n/a |
| MULTIPLYPRODUCTSCANNING | w32-opt | 61.39 M | none | n/a | 4.22 M | n/a |
| MONTGOMERYMULTIPLICATION | w32 | 380.52 M | none | n/a | 631.43 k | n/a |
| COMPARE | w32-opt | 654.07 M | none | n/a | 53.25 M | n/a |
| REDUCE | w32-o64 | 86.68 M | none | n/a | 26.89 M | n/a |
| MODMUL | w32-o64 | 25.51 M | none | n/a | 1.38 M | n/a |
| MODEXP | w32-o64 | 330.36 k | none | n/a | 1.95 k | n/a |
| EXPONENTIATION | w32-o64 | 216.04 k | none | n/a | 14.36 k | n/a |
| DIVIDE | w32-o64 | 18.31 M | none | n/a | 6.42 M | n/a |
| ISQRT | w32-o64 | 1.27 M | none | n/a | 2.40 M | n/a |
| MODMUL_R2 | w32-o64 | 225.82 M | none | n/a | 1.33 M | n/a |

### p2048 (2048-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w32-opt | 286.59 M | none | n/a | 21.52 M | n/a |
| SUBTRACT | w32-opt | 294.16 M | none | n/a | 24.64 M | n/a |
| ADDMOD | w32-o64 | 249.85 M | none | n/a | 7.02 M | n/a |
| SUBTRACTMOD | w32-o64 | 245.56 M | none | n/a | 10.11 M | n/a |
| MULTIPLYOPERANDSCANNING | w32-o64 | 55.52 M | none | n/a | 1.32 M | n/a |
| MULTIPLYPRODUCTSCANNING | w32 | 16.28 M | none | n/a | 1.26 M | n/a |
| MONTGOMERYMULTIPLICATION | w32 | 125.10 M | none | n/a | 156.13 k | n/a |
| COMPARE | w32-o64 | 349.31 M | none | n/a | 39.40 M | n/a |
| REDUCE | w32-opt | 27.16 M | none | n/a | 18.74 M | n/a |
| MODMUL | w32-o64 | 6.08 M | none | n/a | 415.85 k | n/a |
| MODEXP | w32-opt | 24.15 k | none | n/a | 276.2 | n/a |
| EXPONENTIATION | w32-o64 | 27.95 k | none | n/a | 2.56 k | n/a |
| DIVIDE | w32-o64 | 2.14 M | none | n/a | 12.14 M | n/a |
| ISQRT | w32-o64 | 553.95 k | none | n/a | 1.35 M | n/a |
| MODMUL_R2 | w32-o64 | 65.36 M | none | n/a | 385.01 k | n/a |

## 6. Raw data

Also written to `NVIDIA_GeForce_RTX_4070_Ti_Report.csv` for analysis.

```csv
kind,device,device_type,kernel,modulus,bits,operation,items,seconds,ops_per_sec,mismatches
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,secp256k1,256,ADD,200000,0.005409051,36975062.162,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,secp256k1,256,ADD,200000,0.006175538,32385842.634,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,secp256k1,256,ADD,200000,0.005328600,37533310.376,0
library,NVIDIA GeForce RTX 4070 Ti,gpu,cgbn,secp256k1,256,ADD,2000,0.000007040,284090909.091,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,secp256k1,256,ADD,200000,0.000163328,1224529638.333,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,secp256k1,256,ADD,200000,0.004103136,48743205.006,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,secp256k1,256,ADD,200000,0.000091416,2187802106.062,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,secp256k1,256,ADD,200000,0.003824367,52296235.236,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32,secp256k1,256,ADD,200000,0.000051723,3866753361.828,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32,secp256k1,256,ADD,200000,0.003703160,54007928.676,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,secp256k1,256,ADD,200000,0.000050670,3947108613.860,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,secp256k1,256,ADD,200000,0.003961497,50485965.453,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,secp256k1,256,ADD,200000,0.000049425,4046530552.716,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,secp256k1,256,ADD,200000,0.003870401,51674232.881,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,secp256k1,256,SUBTRACT,200000,0.005035514,39717891.664,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,secp256k1,256,SUBTRACT,200000,0.004499210,44452247.842,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,secp256k1,256,SUBTRACT,200000,0.004957832,40340212.842,0
library,NVIDIA GeForce RTX 4070 Ti,gpu,cgbn,secp256k1,256,SUBTRACT,2000,0.000005376,372023809.524,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,secp256k1,256,SUBTRACT,200000,0.000162009,1234499063.627,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,secp256k1,256,SUBTRACT,200000,0.003896964,51322003.285,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,secp256k1,256,SUBTRACT,200000,0.000088251,2266264221.510,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,secp256k1,256,SUBTRACT,200000,0.003982747,50216597.068,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32,secp256k1,256,SUBTRACT,200000,0.000048956,4085299511.805,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32,secp256k1,256,SUBTRACT,200000,0.003991560,50105722.082,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,secp256k1,256,SUBTRACT,200000,0.000052420,3815340803.582,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,secp256k1,256,SUBTRACT,200000,0.003971611,50357398.630,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,secp256k1,256,SUBTRACT,200000,0.000048883,4091409447.465,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,secp256k1,256,SUBTRACT,200000,0.003749138,53345596.033,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,secp256k1,256,ADDMOD,200000,0.016317094,12257084.432,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,secp256k1,256,ADDMOD,200000,0.005466289,36587893.519,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,secp256k1,256,ADDMOD,200000,0.008312519,24060095.552,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,secp256k1,256,ADDMOD,200000,0.000257179,777668397.679,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,secp256k1,256,ADDMOD,200000,0.004069400,49147294.201,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,secp256k1,256,ADDMOD,200000,0.000129663,1542459638.875,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,secp256k1,256,ADDMOD,200000,0.003923073,50980442.335,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32,secp256k1,256,ADDMOD,200000,0.000061567,3248496600.965,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32,secp256k1,256,ADDMOD,200000,0.003856911,51854969.252,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,secp256k1,256,ADDMOD,200000,0.000033378,5991967348.761,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,secp256k1,256,ADDMOD,200000,0.003942843,50724820.343,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,secp256k1,256,ADDMOD,200000,0.000030783,6497091483.375,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,secp256k1,256,ADDMOD,200000,0.003739274,53486318.603,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,secp256k1,256,SUBTRACTMOD,200000,0.011896028,16812334.405,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,secp256k1,256,SUBTRACTMOD,200000,0.005082609,39349869.209,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,secp256k1,256,SUBTRACTMOD,200000,0.007103478,28155221.712,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,secp256k1,256,SUBTRACTMOD,200000,0.000256846,778676483.927,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,secp256k1,256,SUBTRACTMOD,200000,0.004056131,49308072.132,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,secp256k1,256,SUBTRACTMOD,200000,0.000130100,1537279690.752,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,secp256k1,256,SUBTRACTMOD,200000,0.004255676,46996059.534,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32,secp256k1,256,SUBTRACTMOD,200000,0.000062728,3188371298.501,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32,secp256k1,256,SUBTRACTMOD,200000,0.003969512,50384026.230,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,secp256k1,256,SUBTRACTMOD,200000,0.000034136,5858917416.600,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,secp256k1,256,SUBTRACTMOD,200000,0.003887186,51451101.881,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,secp256k1,256,SUBTRACTMOD,200000,0.000028218,7087697175.626,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,secp256k1,256,SUBTRACTMOD,200000,0.003716110,53819719.808,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,secp256k1,256,MULTIPLYOPERANDSCANNING,200000,0.006407681,31212539.891,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,200000,0.005663278,35315235.865,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,200000,0.006503451,30752903.185,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,200000,0.006921809,28894180.488,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,200000,0.011709681,17079884.565,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,secp256k1,256,MULTIPLYOPERANDSCANNING,200000,0.001882973,106215011.510,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,secp256k1,256,MULTIPLYOPERANDSCANNING,200000,0.006755039,29607527.152,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32,secp256k1,256,MULTIPLYOPERANDSCANNING,200000,0.000527356,379250411.461,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32,secp256k1,256,MULTIPLYOPERANDSCANNING,200000,0.005556983,35990752.912,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,secp256k1,256,MULTIPLYOPERANDSCANNING,200000,0.000103081,1940221871.820,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,secp256k1,256,MULTIPLYOPERANDSCANNING,200000,0.004844409,41284706.173,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,secp256k1,256,MULTIPLYOPERANDSCANNING,200000,0.000102291,1955206223.931,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,secp256k1,256,MULTIPLYOPERANDSCANNING,200000,0.004773840,41894994.746,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,secp256k1,256,MULTIPLYPRODUCTSCANNING,200000,0.006691410,29889066.848,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,200000,0.005132823,38964912.541,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,200000,0.006735300,29694296.932,0
library,NVIDIA GeForce RTX 4070 Ti,gpu,cgbn,secp256k1,256,MULTIPLYPRODUCTSCANNING,2000,0.000004992,400641025.641,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,200000,0.000781518,255912192.909,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,200000,0.005612452,35635048.328,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,secp256k1,256,MULTIPLYPRODUCTSCANNING,200000,0.000259118,771849351.830,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,secp256k1,256,MULTIPLYPRODUCTSCANNING,200000,0.005202717,38441453.431,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32,secp256k1,256,MULTIPLYPRODUCTSCANNING,200000,0.000108801,1838219330.126,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32,secp256k1,256,MULTIPLYPRODUCTSCANNING,200000,0.004841712,41307703.121,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,secp256k1,256,MULTIPLYPRODUCTSCANNING,200000,0.000107251,1864784629.195,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,secp256k1,256,MULTIPLYPRODUCTSCANNING,200000,0.004942962,40461569.657,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,secp256k1,256,MULTIPLYPRODUCTSCANNING,200000,0.000100529,1989474601.465,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,secp256k1,256,MULTIPLYPRODUCTSCANNING,200000,0.004834225,41371677.538,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,secp256k1,256,MONTGOMERYMULTIPLICATION,200000,0.050475594,3962310.973,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,200000,0.006547188,30547465.395,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,200000,0.005127565,39004868.753,0
library,NVIDIA GeForce RTX 4070 Ti,gpu,cgbn,secp256k1,256,MONTGOMERYMULTIPLICATION,2000,0.000012256,163185378.590,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,200000,0.000359435,556428864.951,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,200000,0.004220440,47388424.674,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,secp256k1,256,MONTGOMERYMULTIPLICATION,200000,0.000087003,2298771550.679,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,secp256k1,256,MONTGOMERYMULTIPLICATION,200000,0.003877472,51579999.850,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32,secp256k1,256,MONTGOMERYMULTIPLICATION,200000,0.000023752,8420348767.816,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32,secp256k1,256,MONTGOMERYMULTIPLICATION,200000,0.003929423,50898057.389,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,secp256k1,256,MONTGOMERYMULTIPLICATION,200000,0.000035095,5698831758.458,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,secp256k1,256,MONTGOMERYMULTIPLICATION,200000,0.003787465,52805767.885,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,secp256k1,256,MONTGOMERYMULTIPLICATION,200000,0.000033999,5882509564.801,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,secp256k1,256,MONTGOMERYMULTIPLICATION,200000,0.004000748,49990651.223,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,secp256k1,256,COMPARE,200000,0.002126259,94061916.777,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,secp256k1,256,COMPARE,200000,0.004720540,42368033.997,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,secp256k1,256,COMPARE,200000,0.006214632,32182114.360,0
library,NVIDIA GeForce RTX 4070 Ti,gpu,cgbn,secp256k1,256,COMPARE,2000,0.000004832,413907284.768,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,secp256k1,256,COMPARE,200000,0.000153682,1301388904.931,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,secp256k1,256,COMPARE,200000,0.003994569,50067980.432,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,secp256k1,256,COMPARE,200000,0.000091002,2197751724.314,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,secp256k1,256,COMPARE,200000,0.003988592,50143008.435,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,secp256k1,256,COMPARE,200000,0.000047432,4216560356.177,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,secp256k1,256,COMPARE,200000,0.004074700,49083368.777,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,secp256k1,256,COMPARE,200000,0.000047921,4173538202.012,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,secp256k1,256,COMPARE,200000,0.003759539,53198011.952,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,secp256k1,256,REDUCE,25000,0.000736127,33961529.920,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,secp256k1,256,REDUCE,25000,0.004244640,5889780.989,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,secp256k1,256,REDUCE,25000,0.004431988,5640809.502,0
library,NVIDIA GeForce RTX 4070 Ti,gpu,cgbn,secp256k1,256,REDUCE,2000,0.000005120,390625000.000,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,secp256k1,256,REDUCE,25000,0.000111628,223958197.767,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,secp256k1,256,REDUCE,25000,0.000808836,30908615.848,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,secp256k1,256,REDUCE,25000,0.000073170,341670362.834,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,secp256k1,256,REDUCE,25000,0.000772983,32342234.928,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,secp256k1,256,REDUCE,25000,0.000028823,867361765.514,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,secp256k1,256,REDUCE,25000,0.000778017,32132974.638,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,secp256k1,256,REDUCE,25000,0.000028086,890122834.985,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,secp256k1,256,REDUCE,25000,0.000835039,29938718.886,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,secp256k1,256,MODMUL,20000,0.002803400,7134194.381,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,secp256k1,256,MODMUL,20000,0.004400756,4544673.711,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,secp256k1,256,MODMUL,20000,0.006439646,3105760.779,0
library,NVIDIA GeForce RTX 4070 Ti,gpu,cgbn,secp256k1,256,MODMUL,2000,0.000009216,217013888.889,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,secp256k1,256,MODMUL,20000,0.000271013,73797196.132,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,secp256k1,256,MODMUL,20000,0.000929931,21506973.275,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,secp256k1,256,MODMUL,20000,0.000164234,121777404.973,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,secp256k1,256,MODMUL,20000,0.000820611,24372083.090,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,secp256k1,256,MODMUL,20000,0.000078703,254119777.178,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,secp256k1,256,MODMUL,20000,0.000709910,28172587.411,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,secp256k1,256,MODMUL,20000,0.000063064,317137931.758,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,secp256k1,256,MODMUL,20000,0.000642445,31131072.570,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,secp256k1,256,MODEXP,20000,0.381001183,52493.275,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,secp256k1,256,MODEXP,20000,0.021526369,929093.056,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,secp256k1,256,MODEXP,20000,0.036302281,550929.569,0
library,NVIDIA GeForce RTX 4070 Ti,gpu,cgbn,secp256k1,256,MODEXP,2000,0.001370112,1459734.679,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,secp256k1,256,MODEXP,20000,0.012089494,1654328.952,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,secp256k1,256,MODEXP,20000,0.012715905,1572833.381,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,secp256k1,256,MODEXP,20000,0.002369381,8441023.541,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,secp256k1,256,MODEXP,20000,0.002940224,6802202.784,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,secp256k1,256,MODEXP,20000,0.001543187,12960192.506,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,secp256k1,256,MODEXP,20000,0.002247913,8897141.315,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,secp256k1,256,MODEXP,20000,0.000895100,22343870.604,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,secp256k1,256,MODEXP,20000,0.001518029,13174978.688,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,secp256k1,256,EXPONENTIATION,20000,0.125658148,159161.983,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,secp256k1,256,EXPONENTIATION,20000,0.009554424,2093271.144,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,secp256k1,256,EXPONENTIATION,20000,0.182293111,109713.416,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,secp256k1,256,EXPONENTIATION,20000,0.022569814,886139.336,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,secp256k1,256,EXPONENTIATION,20000,0.023149803,863938.240,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,secp256k1,256,EXPONENTIATION,20000,0.006200240,3225681.590,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,secp256k1,256,EXPONENTIATION,20000,0.006802160,2940242.485,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,secp256k1,256,EXPONENTIATION,20000,0.000388134,51528590.739,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,secp256k1,256,EXPONENTIATION,20000,0.000992002,20161248.007,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,secp256k1,256,EXPONENTIATION,20000,0.000344169,58110970.420,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,secp256k1,256,EXPONENTIATION,20000,0.001006691,19867069.081,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,secp256k1,256,DIVIDE,25000,0.001489775,16781057.289,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,secp256k1,256,DIVIDE,25000,0.004176745,5985522.112,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,secp256k1,256,DIVIDE,25000,0.006294370,3971803.344,0
library,NVIDIA GeForce RTX 4070 Ti,gpu,cgbn,secp256k1,256,DIVIDE,2000,0.000005952,336021505.376,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,secp256k1,256,DIVIDE,25000,0.000234846,106452727.590,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,secp256k1,256,DIVIDE,25000,0.001223946,20425737.900,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,secp256k1,256,DIVIDE,25000,0.000205813,121469503.597,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,secp256k1,256,DIVIDE,25000,0.001084251,23057392.117,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,secp256k1,256,DIVIDE,25000,0.000074702,334663214.506,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,secp256k1,256,DIVIDE,25000,0.000905846,27598509.733,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,secp256k1,256,DIVIDE,25000,0.000071724,348558627.764,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,secp256k1,256,DIVIDE,25000,0.000914107,27349093.595,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,secp256k1,256,ISQRT,20000,0.002582119,7745576.167,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,secp256k1,256,ISQRT,20000,0.004406822,4538417.903,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,secp256k1,256,ISQRT,20000,0.002518092,7942521.405,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,secp256k1,256,ISQRT,20000,0.003152526,6344119.003,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,secp256k1,256,ISQRT,20000,0.001958847,10210087.634,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,secp256k1,256,ISQRT,20000,0.002539662,7875063.858,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,secp256k1,256,ISQRT,20000,0.000467008,42825820.687,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,secp256k1,256,ISQRT,20000,0.001164068,17181126.842,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,secp256k1,256,ISQRT,20000,0.000423987,47171265.353,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,secp256k1,256,ISQRT,20000,0.001135400,17614936.990,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,secp256k1,256,MODMUL_R2,200000,0.027910421,7165782.296,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,secp256k1,256,MODMUL_R2,200000,0.005167110,38706356.103,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,secp256k1,256,MODMUL_R2,200000,0.011630077,17196790.654,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,secp256k1,256,MODMUL_R2,200000,0.000415404,481459038.930,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,secp256k1,256,MODMUL_R2,200000,0.004641402,43090428.923,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,secp256k1,256,MODMUL_R2,200000,0.000084484,2367312264.921,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,secp256k1,256,MODMUL_R2,200000,0.003773403,53002554.922,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,secp256k1,256,MODMUL_R2,200000,0.000053023,3771948866.211,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,secp256k1,256,MODMUL_R2,200000,0.003976805,50291628.339,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,secp256k1,256,MODMUL_R2,200000,0.000051446,3887569454.967,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,secp256k1,256,MODMUL_R2,200000,0.003872468,51646649.882,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,rsa256(composite),256,ADD,200000,0.005482509,36479648.828,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,rsa256(composite),256,ADD,200000,0.004663243,42888608.160,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,rsa256(composite),256,ADD,200000,0.005136693,38935556.300,0
library,NVIDIA GeForce RTX 4070 Ti,gpu,cgbn,rsa256(composite),256,ADD,2000,0.000006304,317258883.249,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,rsa256(composite),256,ADD,200000,0.000157843,1267081889.958,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,rsa256(composite),256,ADD,200000,0.004480272,44640147.352,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,rsa256(composite),256,ADD,200000,0.000084944,2354490448.206,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,rsa256(composite),256,ADD,200000,0.003920547,51013290.174,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32,rsa256(composite),256,ADD,200000,0.000049875,4010024948.240,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32,rsa256(composite),256,ADD,200000,0.003936912,50801237.975,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,rsa256(composite),256,ADD,200000,0.000049280,4058441334.996,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,rsa256(composite),256,ADD,200000,0.004110047,48661243.882,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,rsa256(composite),256,ADD,200000,0.000051138,3910986628.786,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,rsa256(composite),256,ADD,200000,0.003751915,53306112.415,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACT,200000,0.004338571,46098127.841,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACT,200000,0.004598936,43488319.283,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACT,200000,0.004836341,41353576.867,0
library,NVIDIA GeForce RTX 4070 Ti,gpu,cgbn,rsa256(composite),256,SUBTRACT,2000,0.000005248,381097560.976,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,rsa256(composite),256,SUBTRACT,200000,0.000159014,1257750760.220,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,rsa256(composite),256,SUBTRACT,200000,0.003951046,50619506.402,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,rsa256(composite),256,SUBTRACT,200000,0.000084017,2380469084.011,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,rsa256(composite),256,SUBTRACT,200000,0.004119745,48546694.148,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32,rsa256(composite),256,SUBTRACT,200000,0.000051055,3917345028.605,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32,rsa256(composite),256,SUBTRACT,200000,0.003762680,53153603.432,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,rsa256(composite),256,SUBTRACT,200000,0.000050042,3996638234.952,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,rsa256(composite),256,SUBTRACT,200000,0.004293137,46585981.331,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,rsa256(composite),256,SUBTRACT,200000,0.000049034,4078801037.989,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,rsa256(composite),256,SUBTRACT,200000,0.003792906,52730018.298,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,rsa256(composite),256,ADDMOD,200000,0.012689658,15760866.084,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,rsa256(composite),256,ADDMOD,200000,0.006790240,29454040.150,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,rsa256(composite),256,ADDMOD,200000,0.007971097,25090649.327,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,rsa256(composite),256,ADDMOD,200000,0.000230895,866194567.739,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,rsa256(composite),256,ADDMOD,200000,0.004145430,48245899.100,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,rsa256(composite),256,ADDMOD,200000,0.000117412,1703402387.161,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,rsa256(composite),256,ADDMOD,200000,0.004277457,46756752.712,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32,rsa256(composite),256,ADDMOD,200000,0.000057860,3456616752.043,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32,rsa256(composite),256,ADDMOD,200000,0.003779017,52923815.661,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,rsa256(composite),256,ADDMOD,200000,0.000031050,6441235160.863,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,rsa256(composite),256,ADDMOD,200000,0.003997579,50030281.026,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,rsa256(composite),256,ADDMOD,200000,0.000031335,6382629746.699,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,rsa256(composite),256,ADDMOD,200000,0.003741631,53452624.603,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACTMOD,200000,0.011993287,16675995.476,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACTMOD,200000,0.004856436,41182464.431,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACTMOD,200000,0.006473859,30893475.024,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,rsa256(composite),256,SUBTRACTMOD,200000,0.000257676,776168516.775,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,rsa256(composite),256,SUBTRACTMOD,200000,0.004083706,48975121.316,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,rsa256(composite),256,SUBTRACTMOD,200000,0.000135086,1480538201.624,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,rsa256(composite),256,SUBTRACTMOD,200000,0.003875082,51611810.927,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32,rsa256(composite),256,SUBTRACTMOD,200000,0.000064508,3100387855.338,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32,rsa256(composite),256,SUBTRACTMOD,200000,0.003738139,53502559.201,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,rsa256(composite),256,SUBTRACTMOD,200000,0.000029779,6716133379.203,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,rsa256(composite),256,SUBTRACTMOD,200000,0.004002806,49964950.641,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,rsa256(composite),256,SUBTRACTMOD,200000,0.000030706,6513373439.034,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,rsa256(composite),256,SUBTRACTMOD,200000,0.003743879,53420529.436,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYOPERANDSCANNING,200000,0.006214952,32180457.813,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,200000,0.005721653,34954933.498,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,200000,0.004592899,43545481.838,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,200000,0.006968366,28701133.165,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,200000,0.012003581,16661694.646,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,rsa256(composite),256,MULTIPLYOPERANDSCANNING,200000,0.001893136,105644816.373,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,rsa256(composite),256,MULTIPLYOPERANDSCANNING,200000,0.006633202,30151350.902,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32,rsa256(composite),256,MULTIPLYOPERANDSCANNING,200000,0.000524866,381049606.338,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32,rsa256(composite),256,MULTIPLYOPERANDSCANNING,200000,0.005462938,36610336.860,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,200000,0.000100347,1993084383.520,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,200000,0.004988414,40092903.255,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,200000,0.000102209,1956774007.283,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,200000,0.004815463,41532869.453,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,200000,0.006387676,31310291.901,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,200000,0.004971544,40228950.800,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,200000,0.006640838,30116681.087,0
library,NVIDIA GeForce RTX 4070 Ti,gpu,cgbn,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,2000,0.000005120,390625000.000,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,200000,0.000783424,255289597.405,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,200000,0.005604382,35686360.995,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,200000,0.000259037,772090434.320,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,200000,0.004982561,40140000.378,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,200000,0.000109026,1834425225.141,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,200000,0.005005241,39958116.653,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,200000,0.000105197,1901193531.909,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,200000,0.005251655,38083232.767,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,200000,0.000101768,1965255346.039,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,200000,0.004886051,40932851.015,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,rsa256(composite),256,MONTGOMERYMULTIPLICATION,200000,0.048117425,4156498.394,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,200000,0.007635101,26194807.454,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,200000,0.006179729,32363878.683,0
library,NVIDIA GeForce RTX 4070 Ti,gpu,cgbn,rsa256(composite),256,MONTGOMERYMULTIPLICATION,2000,0.000011264,177556818.182,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,200000,0.000360797,554328268.214,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,200000,0.004385151,45608463.976,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,rsa256(composite),256,MONTGOMERYMULTIPLICATION,200000,0.000087073,2296921351.103,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,rsa256(composite),256,MONTGOMERYMULTIPLICATION,200000,0.004174697,47907668.738,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32,rsa256(composite),256,MONTGOMERYMULTIPLICATION,200000,0.000024743,8083084761.999,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32,rsa256(composite),256,MONTGOMERYMULTIPLICATION,200000,0.003738040,53503975.516,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,200000,0.000029856,6698823284.632,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,200000,0.004226203,47323803.870,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,200000,0.000034178,5851713183.486,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,200000,0.003883022,51506275.266,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,rsa256(composite),256,COMPARE,200000,0.002201382,90852016.401,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,rsa256(composite),256,COMPARE,200000,0.004503412,44410770.661,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,rsa256(composite),256,COMPARE,200000,0.005433020,36811939.193,0
library,NVIDIA GeForce RTX 4070 Ti,gpu,cgbn,rsa256(composite),256,COMPARE,2000,0.000004992,400641025.641,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,rsa256(composite),256,COMPARE,200000,0.000158889,1258740486.766,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,rsa256(composite),256,COMPARE,200000,0.003837086,52122886.630,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,rsa256(composite),256,COMPARE,200000,0.000090788,2202934256.914,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,rsa256(composite),256,COMPARE,200000,0.003922164,50992257.173,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,rsa256(composite),256,COMPARE,200000,0.000047652,4197091130.829,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,rsa256(composite),256,COMPARE,200000,0.004224190,47346354.904,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,rsa256(composite),256,COMPARE,200000,0.000050422,3966528795.417,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,rsa256(composite),256,COMPARE,200000,0.003826986,52260447.930,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,rsa256(composite),256,REDUCE,25000,0.001265593,19753586.032,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,rsa256(composite),256,REDUCE,25000,0.004309018,5801786.051,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,rsa256(composite),256,REDUCE,25000,0.004447814,5620738.687,0
library,NVIDIA GeForce RTX 4070 Ti,gpu,cgbn,rsa256(composite),256,REDUCE,2000,0.000006176,323834196.891,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,rsa256(composite),256,REDUCE,25000,0.000113966,219363554.626,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,rsa256(composite),256,REDUCE,25000,0.000812790,30758252.162,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,rsa256(composite),256,REDUCE,25000,0.000077263,323570200.276,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,rsa256(composite),256,REDUCE,25000,0.001000111,24997225.525,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,rsa256(composite),256,REDUCE,25000,0.000028903,864961694.895,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,rsa256(composite),256,REDUCE,25000,0.000870514,28718662.458,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,rsa256(composite),256,REDUCE,25000,0.000028964,863139983.682,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,rsa256(composite),256,REDUCE,25000,0.000735016,34012866.987,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,rsa256(composite),256,MODMUL,20000,0.002696308,7417550.336,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,rsa256(composite),256,MODMUL,20000,0.005386507,3712981.329,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,rsa256(composite),256,MODMUL,20000,0.006568301,3044927.446,0
library,NVIDIA GeForce RTX 4070 Ti,gpu,cgbn,rsa256(composite),256,MODMUL,2000,0.000008192,244140625.000,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,rsa256(composite),256,MODMUL,20000,0.000278814,71732410.119,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,rsa256(composite),256,MODMUL,20000,0.000938550,21309467.025,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,rsa256(composite),256,MODMUL,20000,0.000165946,120521157.545,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,rsa256(composite),256,MODMUL,20000,0.000905113,22096686.969,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,rsa256(composite),256,MODMUL,20000,0.000073591,271772310.399,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,rsa256(composite),256,MODMUL,20000,0.000787817,25386604.406,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,rsa256(composite),256,MODMUL,20000,0.000063865,313160671.678,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,rsa256(composite),256,MODMUL,20000,0.000643289,31090227.829,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,rsa256(composite),256,MODEXP,20000,0.274784831,72784.221,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,rsa256(composite),256,MODEXP,20000,0.016929333,1181381.453,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,rsa256(composite),256,MODEXP,20000,0.021710928,921195.080,0
library,NVIDIA GeForce RTX 4070 Ti,gpu,cgbn,rsa256(composite),256,MODEXP,2000,0.001345536,1486396.499,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,rsa256(composite),256,MODEXP,20000,0.012082205,1655326.982,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,rsa256(composite),256,MODEXP,20000,0.012662262,1579496.608,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,rsa256(composite),256,MODEXP,20000,0.002365593,8454539.787,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,rsa256(composite),256,MODEXP,20000,0.003122604,6404910.801,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,rsa256(composite),256,MODEXP,20000,0.001537811,13005499.282,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,rsa256(composite),256,MODEXP,20000,0.002216259,9024215.961,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,rsa256(composite),256,MODEXP,20000,0.000893485,22384258.540,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,rsa256(composite),256,MODEXP,20000,0.001483137,13484930.737,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,rsa256(composite),256,EXPONENTIATION,20000,0.121653766,164400.994,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,rsa256(composite),256,EXPONENTIATION,20000,0.010069839,1986129.059,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,rsa256(composite),256,EXPONENTIATION,20000,0.176256950,113470.703,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,rsa256(composite),256,EXPONENTIATION,20000,0.022520417,888083.023,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,rsa256(composite),256,EXPONENTIATION,20000,0.023252822,860110.657,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,rsa256(composite),256,EXPONENTIATION,20000,0.006203652,3223907.461,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,rsa256(composite),256,EXPONENTIATION,20000,0.006781036,2949401.817,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,rsa256(composite),256,EXPONENTIATION,20000,0.000386268,51777520.079,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,rsa256(composite),256,EXPONENTIATION,20000,0.001090485,18340462.843,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,rsa256(composite),256,EXPONENTIATION,20000,0.000342027,58474926.843,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,rsa256(composite),256,EXPONENTIATION,20000,0.000981517,20376623.096,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,rsa256(composite),256,DIVIDE,25000,0.001663869,15025222.261,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,rsa256(composite),256,DIVIDE,25000,0.004314149,5794885.577,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,rsa256(composite),256,DIVIDE,25000,0.004535557,5512002.160,0
library,NVIDIA GeForce RTX 4070 Ti,gpu,cgbn,rsa256(composite),256,DIVIDE,2000,0.000005120,390625000.000,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,rsa256(composite),256,DIVIDE,25000,0.000234932,106413745.143,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,rsa256(composite),256,DIVIDE,25000,0.001165429,21451327.970,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,rsa256(composite),256,DIVIDE,25000,0.000207753,120335184.259,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,rsa256(composite),256,DIVIDE,25000,0.001056956,23652829.203,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,rsa256(composite),256,DIVIDE,25000,0.000077511,322535069.554,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,rsa256(composite),256,DIVIDE,25000,0.001101337,22699681.800,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,rsa256(composite),256,DIVIDE,25000,0.000075056,333084697.686,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,rsa256(composite),256,DIVIDE,25000,0.000912329,27402396.934,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,rsa256(composite),256,ISQRT,20000,0.002536314,7885459.109,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,rsa256(composite),256,ISQRT,20000,0.004198823,4763239.594,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,rsa256(composite),256,ISQRT,20000,0.002515448,7950869.907,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,rsa256(composite),256,ISQRT,20000,0.003082224,6488821.196,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,rsa256(composite),256,ISQRT,20000,0.002008836,9956014.240,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,rsa256(composite),256,ISQRT,20000,0.002673529,7480749.305,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,rsa256(composite),256,ISQRT,20000,0.000485051,41232781.184,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,rsa256(composite),256,ISQRT,20000,0.001205364,16592498.071,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,rsa256(composite),256,ISQRT,20000,0.000459944,43483553.155,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,rsa256(composite),256,ISQRT,20000,0.001136802,17593212.925,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,rsa256(composite),256,MODMUL_R2,200000,0.029072860,6879268.178,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,rsa256(composite),256,MODMUL_R2,200000,0.007096223,28184006.914,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,rsa256(composite),256,MODMUL_R2,200000,0.008528828,23449881.018,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,rsa256(composite),256,MODMUL_R2,200000,0.000414501,482507812.277,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,rsa256(composite),256,MODMUL_R2,200000,0.004340832,46074116.611,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,rsa256(composite),256,MODMUL_R2,200000,0.000085576,2337104632.345,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,rsa256(composite),256,MODMUL_R2,200000,0.003756537,53240524.653,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,rsa256(composite),256,MODMUL_R2,200000,0.000051057,3917193185.218,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,rsa256(composite),256,MODMUL_R2,200000,0.003864981,51746696.798,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,rsa256(composite),256,MODMUL_R2,200000,0.000052806,3787449114.638,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,rsa256(composite),256,MODMUL_R2,200000,0.003836733,52127681.824,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,ADD,100000,0.002871797,34821403.319,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,ADD,100000,0.004517655,22135377.734,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,ADD,100000,0.006208042,16108138.360,0
library,NVIDIA GeForce RTX 4070 Ti,gpu,cgbn,brainpoolP512r1,512,ADD,2000,0.000006208,322164948.454,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,brainpoolP512r1,512,ADD,100000,0.000188197,531358072.003,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,brainpoolP512r1,512,ADD,100000,0.004042573,24736720.648,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,brainpoolP512r1,512,ADD,100000,0.000098278,1017521232.122,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,brainpoolP512r1,512,ADD,100000,0.004071118,24563277.449,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32,brainpoolP512r1,512,ADD,100000,0.000067177,1488604823.879,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32,brainpoolP512r1,512,ADD,100000,0.003741180,26729534.945,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,brainpoolP512r1,512,ADD,100000,0.000063001,1587277768.703,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,brainpoolP512r1,512,ADD,100000,0.004138451,24163629.899,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,brainpoolP512r1,512,ADD,100000,0.000062998,1587354031.060,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,brainpoolP512r1,512,ADD,100000,0.004022323,24861255.892,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACT,100000,0.002302989,43421831.096,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACT,100000,0.006142137,16280978.146,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACT,100000,0.006281789,15919032.169,0
library,NVIDIA GeForce RTX 4070 Ti,gpu,cgbn,brainpoolP512r1,512,SUBTRACT,2000,0.000005120,390625000.000,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,brainpoolP512r1,512,SUBTRACT,100000,0.000187091,534499404.953,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,brainpoolP512r1,512,SUBTRACT,100000,0.004206582,23772268.765,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,brainpoolP512r1,512,SUBTRACT,100000,0.000098285,1017448918.937,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,brainpoolP512r1,512,SUBTRACT,100000,0.003781715,26443029.558,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32,brainpoolP512r1,512,SUBTRACT,100000,0.000067570,1479946416.917,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32,brainpoolP512r1,512,SUBTRACT,100000,0.003880222,25771721.450,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,brainpoolP512r1,512,SUBTRACT,100000,0.000060745,1646227846.568,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,brainpoolP512r1,512,SUBTRACT,100000,0.004349771,22989715.629,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,brainpoolP512r1,512,SUBTRACT,100000,0.000062964,1588211019.569,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,brainpoolP512r1,512,SUBTRACT,100000,0.003763586,26570403.513,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,ADDMOD,100000,0.007084740,14114844.064,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,ADDMOD,100000,0.006344536,15761593.987,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,ADDMOD,100000,0.005478772,18252265.337,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,brainpoolP512r1,512,ADDMOD,100000,0.000262968,380274392.508,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,brainpoolP512r1,512,ADDMOD,100000,0.004072239,24556515.226,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,brainpoolP512r1,512,ADDMOD,100000,0.000137019,729825944.278,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,brainpoolP512r1,512,ADDMOD,100000,0.004181025,23917580.204,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32,brainpoolP512r1,512,ADDMOD,100000,0.000079988,1250186960.698,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32,brainpoolP512r1,512,ADDMOD,100000,0.003989922,25063146.603,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,brainpoolP512r1,512,ADDMOD,100000,0.000065289,1531652234.197,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,brainpoolP512r1,512,ADDMOD,100000,0.004119068,24277337.004,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,brainpoolP512r1,512,ADDMOD,100000,0.000065374,1529661155.571,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,brainpoolP512r1,512,ADDMOD,100000,0.003752915,26645953.366,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACTMOD,100000,0.006283480,15914747.941,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACTMOD,100000,0.004590819,21782605.638,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACTMOD,100000,0.007352769,13600318.418,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,brainpoolP512r1,512,SUBTRACTMOD,100000,0.000312011,320501530.584,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,brainpoolP512r1,512,SUBTRACTMOD,100000,0.004019662,24877714.076,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,brainpoolP512r1,512,SUBTRACTMOD,100000,0.000154819,645915593.604,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,brainpoolP512r1,512,SUBTRACTMOD,100000,0.003845147,26006808.241,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32,brainpoolP512r1,512,SUBTRACTMOD,100000,0.000092131,1085411352.807,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32,brainpoolP512r1,512,SUBTRACTMOD,100000,0.003896135,25666461.618,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,brainpoolP512r1,512,SUBTRACTMOD,100000,0.000066476,1504300966.157,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,brainpoolP512r1,512,SUBTRACTMOD,100000,0.003908945,25582350.129,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,brainpoolP512r1,512,SUBTRACTMOD,100000,0.000066411,1505775031.641,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,brainpoolP512r1,512,SUBTRACTMOD,100000,0.003921294,25501785.872,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,100000,0.007044863,14194740.050,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,100000,0.006406878,15608226.023,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,100000,0.005641227,17726639.948,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,100000,0.015683930,6375952.959,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,100000,0.020717017,4826949.748,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,100000,0.004547741,21988938.973,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,100000,0.009223305,10842100.482,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,100000,0.001490193,67105403.878,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,100000,0.006196305,16138650.574,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,100000,0.000223087,448255785.988,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,100000,0.004991287,20034913.241,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,100000,0.000224339,445754048.375,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,100000,0.004982649,20069645.687,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,100000,0.007107704,14069241.031,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,100000,0.006516458,15345759.918,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,100000,0.006353727,15738793.751,0
library,NVIDIA GeForce RTX 4070 Ti,gpu,cgbn,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,2000,0.000004096,488281250.000,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,100000,0.002779296,35980334.365,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,100000,0.007744239,12912824.709,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,100000,0.000757399,132030802.790,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,100000,0.005439085,18385445.500,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,100000,0.000225780,442909125.183,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,100000,0.004921978,20317035.183,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,100000,0.000227870,438846575.310,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,100000,0.005178624,19310148.986,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,100000,0.000250119,399809662.714,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,100000,0.005004994,19980043.976,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,100000,0.063167629,1583089.340,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,100000,0.006490989,15405972.789,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,100000,0.005001314,19994745.443,0
library,NVIDIA GeForce RTX 4070 Ti,gpu,cgbn,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,2000,0.000025600,78125000.000,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,100000,0.000706616,141519585.817,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,100000,0.004610369,21690237.894,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,100000,0.000225727,443013058.000,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,100000,0.003904285,25612884.136,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,100000,0.000047300,2114164416.791,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,100000,0.003761170,26587470.867,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,100000,0.000064089,1560328233.074,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,100000,0.003811443,26236782.294,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,100000,0.000058817,1700189336.426,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,100000,0.003905742,25603329.138,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,COMPARE,100000,0.001111004,90008683.358,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,COMPARE,100000,0.004886257,20465562.919,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,COMPARE,100000,0.004415661,22646666.139,0
library,NVIDIA GeForce RTX 4070 Ti,gpu,cgbn,brainpoolP512r1,512,COMPARE,2000,0.000004864,411184210.526,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,brainpoolP512r1,512,COMPARE,100000,0.000195621,511192742.079,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,brainpoolP512r1,512,COMPARE,100000,0.004003287,24979472.742,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,brainpoolP512r1,512,COMPARE,100000,0.000103416,966969096.378,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,brainpoolP512r1,512,COMPARE,100000,0.004005073,24968333.988,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,brainpoolP512r1,512,COMPARE,100000,0.000057093,1751528692.868,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,brainpoolP512r1,512,COMPARE,100000,0.003787226,26404550.962,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,brainpoolP512r1,512,COMPARE,100000,0.000055799,1792145916.422,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,brainpoolP512r1,512,COMPARE,100000,0.003869317,25844354.401,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,REDUCE,20000,0.001170761,17082904.907,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,REDUCE,20000,0.004232871,4724925.474,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,REDUCE,20000,0.006306886,3171137.084,0
library,NVIDIA GeForce RTX 4070 Ti,gpu,cgbn,brainpoolP512r1,512,REDUCE,2000,0.000005824,343406593.407,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,brainpoolP512r1,512,REDUCE,20000,0.000318851,62725228.126,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,brainpoolP512r1,512,REDUCE,20000,0.001521171,13147766.163,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,brainpoolP512r1,512,REDUCE,20000,0.000242067,82621741.437,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,brainpoolP512r1,512,REDUCE,20000,0.001279640,15629395.627,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,brainpoolP512r1,512,REDUCE,20000,0.000067913,294494541.764,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,brainpoolP512r1,512,REDUCE,20000,0.001094324,18276122.165,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,brainpoolP512r1,512,REDUCE,20000,0.000068148,293478829.971,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,brainpoolP512r1,512,REDUCE,20000,0.001320353,15147465.565,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL,20000,0.005682876,3519344.792,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL,20000,0.006333541,3157791.190,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL,20000,0.006846312,2921280.794,0
library,NVIDIA GeForce RTX 4070 Ti,gpu,cgbn,brainpoolP512r1,512,MODMUL,2000,0.000017408,114889705.882,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,brainpoolP512r1,512,MODMUL,20000,0.000886837,22552058.087,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,brainpoolP512r1,512,MODMUL,20000,0.002109565,9480627.192,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,brainpoolP512r1,512,MODMUL,20000,0.000602177,33212826.404,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,brainpoolP512r1,512,MODMUL,20000,0.001638210,12208446.806,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,brainpoolP512r1,512,MODMUL,20000,0.000237866,84080968.847,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,brainpoolP512r1,512,MODMUL,20000,0.001427106,14014375.028,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,brainpoolP512r1,512,MODMUL,20000,0.000191746,104304627.591,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,brainpoolP512r1,512,MODMUL,20000,0.001378426,14509301.634,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,MODEXP,20000,1.518640262,13169.676,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,MODEXP,20000,0.170375162,117388.003,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,MODEXP,20000,0.101436233,197168.205,0
library,NVIDIA GeForce RTX 4070 Ti,gpu,cgbn,brainpoolP512r1,512,MODEXP,2000,0.003442688,580941.404,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,brainpoolP512r1,512,MODEXP,20000,0.175540639,113933.731,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,brainpoolP512r1,512,MODEXP,20000,0.176796283,113124.550,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,brainpoolP512r1,512,MODEXP,20000,0.016090394,1242977.643,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,brainpoolP512r1,512,MODEXP,20000,0.017280775,1157355.500,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,brainpoolP512r1,512,MODEXP,20000,0.012424958,1609663.382,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,brainpoolP512r1,512,MODEXP,20000,0.013448796,1487121.970,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,brainpoolP512r1,512,MODEXP,20000,0.005944385,3364519.583,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,brainpoolP512r1,512,MODEXP,20000,0.007040187,2840833.624,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,EXPONENTIATION,20000,0.351067291,56969.135,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,EXPONENTIATION,20000,0.022579843,885745.749,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,EXPONENTIATION,20000,0.472826309,42298.831,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,brainpoolP512r1,512,EXPONENTIATION,20000,0.183281499,109121.761,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,brainpoolP512r1,512,EXPONENTIATION,20000,0.184841206,108200.982,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,brainpoolP512r1,512,EXPONENTIATION,20000,0.047050193,425077.959,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,brainpoolP512r1,512,EXPONENTIATION,20000,0.048306583,414022.247,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,brainpoolP512r1,512,EXPONENTIATION,20000,0.014811759,1350278.519,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,brainpoolP512r1,512,EXPONENTIATION,20000,0.016086034,1243314.541,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,brainpoolP512r1,512,EXPONENTIATION,20000,0.014539443,1375568.518,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,brainpoolP512r1,512,EXPONENTIATION,20000,0.015761984,1268875.802,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,DIVIDE,20000,0.002814307,7106545.250,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,DIVIDE,20000,0.000075627,264455809.283,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,DIVIDE,20000,0.005377314,3719328.981,0
library,NVIDIA GeForce RTX 4070 Ti,gpu,cgbn,brainpoolP512r1,512,DIVIDE,2000,0.000007264,275330396.476,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,brainpoolP512r1,512,DIVIDE,20000,0.000779806,25647402.828,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,brainpoolP512r1,512,DIVIDE,20000,0.002093275,9554406.550,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,brainpoolP512r1,512,DIVIDE,20000,0.000780199,25634483.184,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,brainpoolP512r1,512,DIVIDE,20000,0.002077454,9627168.628,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,brainpoolP512r1,512,DIVIDE,20000,0.000258081,77495046.838,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,brainpoolP512r1,512,DIVIDE,20000,0.001487747,13443145.301,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,brainpoolP512r1,512,DIVIDE,20000,0.000252795,79115507.605,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,brainpoolP512r1,512,DIVIDE,20000,0.001675539,11936456.978,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,ISQRT,20000,0.004566587,4379638.519,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,ISQRT,20000,0.003356499,5958589.769,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,brainpoolP512r1,512,ISQRT,20000,0.011669534,1713864.492,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,brainpoolP512r1,512,ISQRT,20000,0.012789121,1563829.131,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,brainpoolP512r1,512,ISQRT,20000,0.011976404,1669950.346,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,brainpoolP512r1,512,ISQRT,20000,0.013135238,1522621.824,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,brainpoolP512r1,512,ISQRT,20000,0.002422196,8256970.204,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,brainpoolP512r1,512,ISQRT,20000,0.003480558,5746204.979,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,brainpoolP512r1,512,ISQRT,20000,0.002267282,8821134.317,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,brainpoolP512r1,512,ISQRT,20000,0.003521898,5678756.257,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL_R2,100000,0.032816075,3047287.036,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL_R2,100000,0.005454601,18333146.569,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL_R2,100000,0.009338421,10708448.540,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,brainpoolP512r1,512,MODMUL_R2,100000,0.001048978,95330888.286,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,brainpoolP512r1,512,MODMUL_R2,100000,0.004796980,20846449.446,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,brainpoolP512r1,512,MODMUL_R2,100000,0.000207672,481528536.361,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,brainpoolP512r1,512,MODMUL_R2,100000,0.004003317,24979286.059,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,brainpoolP512r1,512,MODMUL_R2,100000,0.000106798,936346953.743,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,brainpoolP512r1,512,MODMUL_R2,100000,0.003849866,25974930.137,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,brainpoolP512r1,512,MODMUL_R2,100000,0.000073241,1365354166.985,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,brainpoolP512r1,512,MODMUL_R2,100000,0.003901846,25628894.576,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p1024,1024,ADD,50000,0.003284623,15222446.996,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p1024,1024,ADD,50000,0.005876524,8508431.139,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p1024,1024,ADD,50000,0.000121916,410118624.588,0
library,NVIDIA GeForce RTX 4070 Ti,gpu,cgbn,p1024,1024,ADD,2000,0.000009216,217013888.889,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p1024,1024,ADD,50000,0.000273590,182755155.841,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p1024,1024,ADD,50000,0.004190805,11930882.292,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p1024,1024,ADD,50000,0.000143581,348235603.328,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p1024,1024,ADD,50000,0.003966619,12605193.403,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32,p1024,1024,ADD,50000,0.000111534,448293683.760,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32,p1024,1024,ADD,50000,0.003989106,12534136.583,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p1024,1024,ADD,50000,0.000079160,631632343.940,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p1024,1024,ADD,50000,0.004176253,11972454.437,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p1024,1024,ADD,50000,0.000078858,634051137.241,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p1024,1024,ADD,50000,0.003856742,12964310.494,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p1024,1024,SUBTRACT,50000,0.001394279,35860826.671,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p1024,1024,SUBTRACT,50000,0.004414179,11327134.639,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p1024,1024,SUBTRACT,50000,0.006266891,7978437.761,0
library,NVIDIA GeForce RTX 4070 Ti,gpu,cgbn,p1024,1024,SUBTRACT,2000,0.000005120,390625000.000,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p1024,1024,SUBTRACT,50000,0.000269358,185626612.187,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p1024,1024,SUBTRACT,50000,0.003965925,12607399.041,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p1024,1024,SUBTRACT,50000,0.000143451,348551273.216,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p1024,1024,SUBTRACT,50000,0.003931943,12716359.388,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32,p1024,1024,SUBTRACT,50000,0.000109961,454706571.328,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32,p1024,1024,SUBTRACT,50000,0.004085406,12238685.814,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p1024,1024,SUBTRACT,50000,0.000078211,639296618.921,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p1024,1024,SUBTRACT,50000,0.004243962,11781443.846,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p1024,1024,SUBTRACT,50000,0.000077873,642071041.852,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p1024,1024,SUBTRACT,50000,0.003874741,12904088.203,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p1024,1024,ADDMOD,50000,0.006758713,7397858.091,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p1024,1024,ADDMOD,50000,0.004820149,10373123.479,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p1024,1024,ADDMOD,50000,0.005091951,9819418.930,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p1024,1024,ADDMOD,50000,0.000424252,117854480.466,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p1024,1024,ADDMOD,50000,0.004073702,12273848.114,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p1024,1024,ADDMOD,50000,0.000212193,235634539.074,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p1024,1024,ADDMOD,50000,0.003936482,12701696.554,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32,p1024,1024,ADDMOD,50000,0.000144416,346222023.596,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32,p1024,1024,ADDMOD,50000,0.003973899,12582101.536,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p1024,1024,ADDMOD,50000,0.000102830,486239410.033,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p1024,1024,ADDMOD,50000,0.004142591,12069741.142,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p1024,1024,ADDMOD,50000,0.000103197,484509924.463,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p1024,1024,ADDMOD,50000,0.003866184,12932648.710,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p1024,1024,SUBTRACTMOD,50000,0.004990413,10019210.661,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p1024,1024,SUBTRACTMOD,50000,0.006246360,8004661.891,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p1024,1024,SUBTRACTMOD,50000,0.006158913,8118315.657,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p1024,1024,SUBTRACTMOD,50000,0.000422048,118469927.745,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p1024,1024,SUBTRACTMOD,50000,0.004307329,11608121.996,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p1024,1024,SUBTRACTMOD,50000,0.000219630,227655642.961,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p1024,1024,SUBTRACTMOD,50000,0.003940247,12689559.819,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32,p1024,1024,SUBTRACTMOD,50000,0.000144423,346205278.804,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32,p1024,1024,SUBTRACTMOD,50000,0.004161298,12015481.586,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p1024,1024,SUBTRACTMOD,50000,0.000103360,483745933.575,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p1024,1024,SUBTRACTMOD,50000,0.004203877,11893782.963,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p1024,1024,SUBTRACTMOD,50000,0.000103616,482550774.503,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p1024,1024,SUBTRACTMOD,50000,0.003954624,12643427.221,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.011657760,4288988.617,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.004746547,10533973.503,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.006710196,7451347.140,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.035850597,1394676.911,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.040963871,1220587.773,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.009701378,5153907.031,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.014438397,3462988.332,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.003753393,13321279.978,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.008495130,5885725.120,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.000310839,160854957.011,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.004981425,10037288.509,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.000306707,163022041.583,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.005150195,9708370.306,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.011851872,4218742.839,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.002935522,17032745.776,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.005065319,9871046.598,0
library,NVIDIA GeForce RTX 4070 Ti,gpu,cgbn,p1024,1024,MULTIPLYPRODUCTSCANNING,2000,0.000006016,332446808.511,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.011977689,4174427.934,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.017149365,2915559.840,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.003061257,16333160.123,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.007900057,6329068.250,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000820455,60941792.639,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.005489180,9108828.673,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000814400,61394889.591,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.005812026,8602852.119,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000815315,61325986.381,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.005679605,8803429.014,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.079185849,631425.950,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.008268393,6047124.276,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.006345749,7879290.508,0
library,NVIDIA GeForce RTX 4070 Ti,gpu,cgbn,p1024,1024,MONTGOMERYMULTIPLICATION,2000,0.000048992,40822991.509,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.002174597,22992766.009,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.006000375,8332812.497,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000434026,115200470.354,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.004348316,11498704.334,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000131399,380520337.483,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.003862013,12946616.207,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000184918,270390208.211,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.003955432,12640844.348,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000151202,330683530.449,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.003974108,12581439.578,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p1024,1024,COMPARE,50000,0.000939005,53247856.289,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p1024,1024,COMPARE,50000,0.004768374,10485754.624,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p1024,1024,COMPARE,50000,0.005041087,9918495.891,0
library,NVIDIA GeForce RTX 4070 Ti,gpu,cgbn,p1024,1024,COMPARE,2000,0.000005120,390625000.000,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p1024,1024,COMPARE,50000,0.000260209,192153255.624,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p1024,1024,COMPARE,50000,0.004125270,12120418.605,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p1024,1024,COMPARE,50000,0.000137582,363419750.554,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p1024,1024,COMPARE,50000,0.003910987,12784496.416,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p1024,1024,COMPARE,50000,0.000076444,654073530.303,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p1024,1024,COMPARE,50000,0.003981650,12557608.052,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p1024,1024,COMPARE,50000,0.000077089,648601425.877,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p1024,1024,COMPARE,50000,0.003859691,12954404.991,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p1024,1024,REDUCE,20000,0.000743806,26888731.273,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p1024,1024,REDUCE,20000,0.004091320,4888397.764,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p1024,1024,REDUCE,20000,0.002898717,6899604.046,0
library,NVIDIA GeForce RTX 4070 Ti,gpu,cgbn,p1024,1024,REDUCE,2000,0.000006272,318877551.020,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p1024,1024,REDUCE,20000,0.001603182,12475190.436,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p1024,1024,REDUCE,20000,0.003483276,5741721.300,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p1024,1024,REDUCE,20000,0.000870507,22975117.389,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p1024,1024,REDUCE,20000,0.002611536,7658328.424,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p1024,1024,REDUCE,20000,0.000231912,86239604.279,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p1024,1024,REDUCE,20000,0.002093058,9555396.573,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p1024,1024,REDUCE,20000,0.000230726,86682872.404,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p1024,1024,REDUCE,20000,0.001975843,10122261.222,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p1024,1024,MODMUL,20000,0.014496329,1379659.636,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p1024,1024,MODMUL,20000,0.004940852,4047884.907,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p1024,1024,MODMUL,20000,0.006468269,3092017.390,0
library,NVIDIA GeForce RTX 4070 Ti,gpu,cgbn,p1024,1024,MODMUL,2000,0.000026720,74850299.401,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p1024,1024,MODMUL,20000,0.006390161,3129811.598,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p1024,1024,MODMUL,20000,0.008276400,2416509.616,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p1024,1024,MODMUL,20000,0.002295875,8711275.693,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p1024,1024,MODMUL,20000,0.004080358,4901530.745,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p1024,1024,MODMUL,20000,0.001042145,19191186.089,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p1024,1024,MODMUL,20000,0.002794847,7156026.765,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p1024,1024,MODMUL,20000,0.000783870,25514434.900,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p1024,1024,MODMUL,20000,0.002623693,7622843.157,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p1024,1024,MODEXP,20000,10.248072672,1951.586,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p1024,1024,MODEXP,20000,1.096651143,18237.340,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p1024,1024,MODEXP,20000,0.782086604,25572.615,0
library,NVIDIA GeForce RTX 4070 Ti,gpu,cgbn,p1024,1024,MODEXP,2000,0.011388832,175610.633,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p1024,1024,MODEXP,20000,1.412693940,14157.348,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p1024,1024,MODEXP,20000,1.414425981,14140.012,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p1024,1024,MODEXP,20000,0.215782906,92685.748,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p1024,1024,MODEXP,20000,0.217769911,91840.052,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p1024,1024,MODEXP,20000,0.095201832,210079.991,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p1024,1024,MODEXP,20000,0.097140355,205887.656,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p1024,1024,MODEXP,20000,0.060539140,330364.785,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p1024,1024,MODEXP,20000,0.062585254,319564.094,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p1024,1024,EXPONENTIATION,20000,1.392895232,14358.582,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p1024,1024,EXPONENTIATION,20000,0.127180196,157257.188,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p1024,1024,EXPONENTIATION,20000,2.558784243,7816.212,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p1024,1024,EXPONENTIATION,20000,1.366174006,14639.424,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p1024,1024,EXPONENTIATION,20000,1.367073118,14629.795,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p1024,1024,EXPONENTIATION,20000,0.364526334,54865.721,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p1024,1024,EXPONENTIATION,20000,0.366423364,54581.672,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p1024,1024,EXPONENTIATION,20000,0.093897393,212998.459,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p1024,1024,EXPONENTIATION,20000,0.095686911,209015.003,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p1024,1024,EXPONENTIATION,20000,0.092574623,216041.927,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p1024,1024,EXPONENTIATION,20000,0.095085666,210336.645,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p1024,1024,DIVIDE,20000,0.003114918,6420714.738,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p1024,1024,DIVIDE,20000,0.005829135,3431040.755,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p1024,1024,DIVIDE,20000,0.005702001,3507540.644,0
library,NVIDIA GeForce RTX 4070 Ti,gpu,cgbn,p1024,1024,DIVIDE,2000,0.000007136,280269058.296,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p1024,1024,DIVIDE,20000,0.013062349,1531118.185,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p1024,1024,DIVIDE,20000,0.015357681,1302279.940,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p1024,1024,DIVIDE,20000,0.005765676,3468803.985,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p1024,1024,DIVIDE,20000,0.007939791,2518958.006,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p1024,1024,DIVIDE,20000,0.001115485,17929421.187,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p1024,1024,DIVIDE,20000,0.003433144,5825563.855,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p1024,1024,DIVIDE,20000,0.001092085,18313594.210,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p1024,1024,DIVIDE,20000,0.003353232,5964394.855,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p1024,1024,ISQRT,20000,0.008334475,2399671.249,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p1024,1024,ISQRT,20000,0.004647372,4303507.391,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p1024,1024,ISQRT,20000,0.178707416,111914.774,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p1024,1024,ISQRT,20000,0.179555279,111386.310,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p1024,1024,ISQRT,20000,0.085021177,235235.511,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p1024,1024,ISQRT,20000,0.086994196,229900.395,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p1024,1024,ISQRT,20000,0.016908552,1182833.395,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p1024,1024,ISQRT,20000,0.018768049,1065640.869,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p1024,1024,ISQRT,20000,0.015718178,1272412.117,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p1024,1024,ISQRT,20000,0.017596053,1136618.539,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p1024,1024,MODMUL_R2,50000,0.037503844,1333196.676,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p1024,1024,MODMUL_R2,50000,0.008170693,6119431.926,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p1024,1024,MODMUL_R2,50000,0.007739331,6460506.736,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p1024,1024,MODMUL_R2,50000,0.003416194,14636171.138,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p1024,1024,MODMUL_R2,50000,0.007349581,6803108.919,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p1024,1024,MODMUL_R2,50000,0.000496642,100676147.940,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p1024,1024,MODMUL_R2,50000,0.004444245,11250504.825,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p1024,1024,MODMUL_R2,50000,0.000272587,183427637.182,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p1024,1024,MODMUL_R2,50000,0.004510649,11084879.611,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p1024,1024,MODMUL_R2,50000,0.000221420,225815230.525,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p1024,1024,MODMUL_R2,50000,0.003995880,12512888.216,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p2048,2048,ADD,25000,0.001161527,21523391.098,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p2048,2048,ADD,25000,0.004487497,5571034.468,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p2048,2048,ADD,25000,0.004288620,5829381.059,0
library,NVIDIA GeForce RTX 4070 Ti,gpu,cgbn,p2048,2048,ADD,2000,0.000006240,320512820.513,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p2048,2048,ADD,25000,0.000291783,85680107.661,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p2048,2048,ADD,25000,0.004254934,5875531.755,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p2048,2048,ADD,25000,0.000148816,167992650.354,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p2048,2048,ADD,25000,0.003947485,6333146.286,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32,p2048,2048,ADD,25000,0.000090826,275251368.000,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32,p2048,2048,ADD,25000,0.003808974,6563447.174,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p2048,2048,ADD,25000,0.000087232,286591760.508,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p2048,2048,ADD,25000,0.004245097,5889146.868,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p2048,2048,ADD,25000,0.000088991,280927236.741,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p2048,2048,ADD,25000,0.003796587,6584861.614,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p2048,2048,SUBTRACT,25000,0.001014599,24640276.735,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p2048,2048,SUBTRACT,25000,0.006389112,3912906.774,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p2048,2048,SUBTRACT,25000,0.002973248,8408313.013,0
library,NVIDIA GeForce RTX 4070 Ti,gpu,cgbn,p2048,2048,SUBTRACT,2000,0.000005120,390625000.000,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p2048,2048,SUBTRACT,25000,0.000300700,83139352.217,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p2048,2048,SUBTRACT,25000,0.004390029,5694723.196,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p2048,2048,SUBTRACT,25000,0.000143652,174031727.745,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p2048,2048,SUBTRACT,25000,0.003988159,6268556.486,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32,p2048,2048,SUBTRACT,25000,0.000087225,286615093.086,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32,p2048,2048,SUBTRACT,25000,0.003953505,6323502.883,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p2048,2048,SUBTRACT,25000,0.000084988,294159320.915,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p2048,2048,SUBTRACT,25000,0.004283605,5836205.682,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p2048,2048,SUBTRACT,25000,0.000087077,287102117.954,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p2048,2048,SUBTRACT,25000,0.003934885,6353425.789,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p2048,2048,ADDMOD,25000,0.003563147,7016269.761,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p2048,2048,ADDMOD,25000,0.004443967,5625604.457,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p2048,2048,ADDMOD,25000,0.004895248,5106993.574,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p2048,2048,ADDMOD,25000,0.000422128,59223741.083,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p2048,2048,ADDMOD,25000,0.004332955,5769734.518,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p2048,2048,ADDMOD,25000,0.000204098,122490139.283,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p2048,2048,ADDMOD,25000,0.004060091,6157497.478,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32,p2048,2048,ADDMOD,25000,0.000104004,240375251.067,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32,p2048,2048,ADDMOD,25000,0.003911074,6392106.077,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p2048,2048,ADDMOD,25000,0.000107220,233165653.433,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p2048,2048,ADDMOD,25000,0.004141944,6035813.143,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p2048,2048,ADDMOD,25000,0.000100059,249852664.107,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p2048,2048,ADDMOD,25000,0.003782055,6610163.008,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p2048,2048,SUBTRACTMOD,25000,0.002472386,10111689.692,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p2048,2048,SUBTRACTMOD,25000,0.004487930,5570496.744,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p2048,2048,SUBTRACTMOD,25000,0.004830298,5175664.095,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p2048,2048,SUBTRACTMOD,25000,0.000453373,55142239.765,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p2048,2048,SUBTRACTMOD,25000,0.004144499,6032092.233,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p2048,2048,SUBTRACTMOD,25000,0.000218153,114598469.945,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p2048,2048,SUBTRACTMOD,25000,0.004135112,6045785.392,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32,p2048,2048,SUBTRACTMOD,25000,0.000125080,199872085.359,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32,p2048,2048,SUBTRACTMOD,25000,0.006675265,3745169.650,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p2048,2048,SUBTRACTMOD,25000,0.000105648,236634709.039,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p2048,2048,SUBTRACTMOD,25000,0.004140298,6038212.710,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p2048,2048,SUBTRACTMOD,25000,0.000101807,245562533.790,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p2048,2048,SUBTRACTMOD,25000,0.003850208,6493155.672,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p2048,2048,MULTIPLYOPERANDSCANNING,25000,0.018972312,1317709.717,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p2048,2048,MULTIPLYOPERANDSCANNING,25000,0.006331547,3948482.072,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p2048,2048,MULTIPLYOPERANDSCANNING,25000,0.005307807,4710043.145,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,25000,0.074897980,333787.373,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,25000,0.080008752,312465.816,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p2048,2048,MULTIPLYOPERANDSCANNING,25000,0.019613007,1274664.310,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p2048,2048,MULTIPLYOPERANDSCANNING,25000,0.024636789,1014742.628,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32,p2048,2048,MULTIPLYOPERANDSCANNING,25000,0.005250459,4761488.513,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32,p2048,2048,MULTIPLYOPERANDSCANNING,25000,0.009936098,2516078.241,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p2048,2048,MULTIPLYOPERANDSCANNING,25000,0.000454061,55058685.648,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p2048,2048,MULTIPLYOPERANDSCANNING,25000,0.005191227,4815817.080,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p2048,2048,MULTIPLYOPERANDSCANNING,25000,0.000450270,55522240.039,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p2048,2048,MULTIPLYOPERANDSCANNING,25000,0.005148397,4855880.437,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p2048,2048,MULTIPLYPRODUCTSCANNING,25000,0.019919947,1255023.418,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,25000,0.005146456,4857711.729,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,25000,0.007112852,3514764.510,0
library,NVIDIA GeForce RTX 4070 Ti,gpu,cgbn,p2048,2048,MULTIPLYPRODUCTSCANNING,2000,0.000009216,217013888.889,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,25000,0.023802573,1050306.621,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,25000,0.028596647,874228.368,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p2048,2048,MULTIPLYPRODUCTSCANNING,25000,0.006086521,4107436.688,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p2048,2048,MULTIPLYPRODUCTSCANNING,25000,0.010938901,2285421.545,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32,p2048,2048,MULTIPLYPRODUCTSCANNING,25000,0.001535600,16280282.014,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32,p2048,2048,MULTIPLYPRODUCTSCANNING,25000,0.006224760,4016219.118,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p2048,2048,MULTIPLYPRODUCTSCANNING,25000,0.001612159,15507155.717,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p2048,2048,MULTIPLYPRODUCTSCANNING,25000,0.006750751,3703291.649,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p2048,2048,MULTIPLYPRODUCTSCANNING,25000,0.001606592,15560889.358,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p2048,2048,MULTIPLYPRODUCTSCANNING,25000,0.006342972,3941370.042,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p2048,2048,MONTGOMERYMULTIPLICATION,25000,0.160118537,156134.327,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p2048,2048,MONTGOMERYMULTIPLICATION,25000,0.010332441,2419563.762,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p2048,2048,MONTGOMERYMULTIPLICATION,25000,0.005989216,4174169.015,0
library,NVIDIA GeForce RTX 4070 Ti,gpu,cgbn,p2048,2048,MONTGOMERYMULTIPLICATION,2000,0.000080896,24723101.266,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,25000,0.026213168,953719.141,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,25000,0.030016158,832884.741,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p2048,2048,MONTGOMERYMULTIPLICATION,25000,0.001004127,24897249.781,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p2048,2048,MONTGOMERYMULTIPLICATION,25000,0.004732906,5282167.020,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32,p2048,2048,MONTGOMERYMULTIPLICATION,25000,0.000199847,125095746.188,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32,p2048,2048,MONTGOMERYMULTIPLICATION,25000,0.003912472,6389822.196,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p2048,2048,MONTGOMERYMULTIPLICATION,25000,0.000295393,84633032.633,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p2048,2048,MONTGOMERYMULTIPLICATION,25000,0.003995914,6256390.884,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p2048,2048,MONTGOMERYMULTIPLICATION,25000,0.000249084,100367713.085,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p2048,2048,MONTGOMERYMULTIPLICATION,25000,0.004068004,6145520.063,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p2048,2048,COMPARE,25000,0.000634570,39396755.370,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p2048,2048,COMPARE,25000,0.004054305,6166284.932,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p2048,2048,COMPARE,25000,0.006135310,4074773.697,0
library,NVIDIA GeForce RTX 4070 Ti,gpu,cgbn,p2048,2048,COMPARE,2000,0.000005120,390625000.000,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p2048,2048,COMPARE,25000,0.000254363,98284758.767,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p2048,2048,COMPARE,25000,0.003903538,6404446.446,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p2048,2048,COMPARE,25000,0.000139165,179642907.167,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p2048,2048,COMPARE,25000,0.003835621,6517849.294,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p2048,2048,COMPARE,25000,0.000071917,347623137.713,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p2048,2048,COMPARE,25000,0.004150773,6022974.444,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p2048,2048,COMPARE,25000,0.000071570,349308152.809,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p2048,2048,COMPARE,25000,0.003792071,6592703.536,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p2048,2048,REDUCE,20000,0.001067454,18736170.446,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p2048,2048,REDUCE,20000,0.003407518,5869374.666,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p2048,2048,REDUCE,20000,0.006854222,2917909.558,0
library,NVIDIA GeForce RTX 4070 Ti,gpu,cgbn,p2048,2048,REDUCE,2000,0.000007072,282805429.864,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p2048,2048,REDUCE,20000,0.146485686,136532.111,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p2048,2048,REDUCE,20000,0.149666810,133630.162,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p2048,2048,REDUCE,20000,0.003914104,5109726.314,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p2048,2048,REDUCE,20000,0.007106113,2814478.177,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p2048,2048,REDUCE,20000,0.000736329,27161771.214,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p2048,2048,REDUCE,20000,0.004019725,4975464.715,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p2048,2048,REDUCE,20000,0.000774004,25839662.596,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p2048,2048,REDUCE,20000,0.003826332,5226938.062,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p2048,2048,MODMUL,20000,0.048094718,415846.081,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p2048,2048,MODMUL,20000,0.006055987,3302517.011,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p2048,2048,MODMUL,20000,0.008196761,2439988.181,0
library,NVIDIA GeForce RTX 4070 Ti,gpu,cgbn,p2048,2048,MODMUL,2000,0.000043744,45720555.962,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p2048,2048,MODMUL,20000,0.267393849,74796.036,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p2048,2048,MODMUL,20000,0.270536099,73927.288,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p2048,2048,MODMUL,20000,0.010293162,1943037.530,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p2048,2048,MODMUL,20000,0.013398301,1492726.581,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p2048,2048,MODMUL,20000,0.004171192,4794792.532,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p2048,2048,MODMUL,20000,0.007397982,2703439.917,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p2048,2048,MODMUL,20000,0.003288350,6082077.727,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p2048,2048,MODMUL,20000,0.006448269,3101607.575,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p2048,2048,MODEXP,20000,72.413732084,276.191,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p2048,2048,MODEXP,20000,7.943209608,2517.874,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p2048,2048,MODEXP,20000,5.138911485,3891.875,0
library,NVIDIA GeForce RTX 4070 Ti,gpu,cgbn,p2048,2048,MODEXP,2000,0.055068672,36318.290,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p2048,2048,MODEXP,20000,35.791975987,558.784,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p2048,2048,MODEXP,20000,35.791458696,558.793,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p2048,2048,MODEXP,20000,5.344981153,3741.828,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p2048,2048,MODEXP,20000,5.347548628,3740.031,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p2048,2048,MODEXP,20000,0.828317936,24145.318,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p2048,2048,MODEXP,20000,0.832313903,24029.396,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p2048,2048,MODEXP,20000,1.494747906,13380.183,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p2048,2048,MODEXP,20000,1.497645794,13354.293,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p2048,2048,EXPONENTIATION,20000,7.806737083,2561.890,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p2048,2048,EXPONENTIATION,20000,0.880183947,22722.523,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p2048,2048,EXPONENTIATION,20000,11.801997208,1694.628,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p2048,2048,EXPONENTIATION,20000,13.460439140,1485.836,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p2048,2048,EXPONENTIATION,20000,13.459027972,1485.991,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p2048,2048,EXPONENTIATION,20000,2.747510969,7279.316,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p2048,2048,EXPONENTIATION,20000,2.753692979,7262.974,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p2048,2048,EXPONENTIATION,20000,0.728779224,27443.153,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p2048,2048,EXPONENTIATION,20000,0.732388024,27307.929,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p2048,2048,EXPONENTIATION,20000,0.715441467,27954.768,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p2048,2048,EXPONENTIATION,20000,0.718459985,27837.319,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p2048,2048,DIVIDE,20000,0.001647948,12136304.870,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p2048,2048,DIVIDE,20000,0.004236105,4721318.270,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p2048,2048,DIVIDE,20000,0.004728346,4229808.888,0
library,NVIDIA GeForce RTX 4070 Ti,gpu,cgbn,p2048,2048,DIVIDE,2000,0.000008032,249003984.064,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p2048,2048,DIVIDE,20000,0.629487340,31771.886,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p2048,2048,DIVIDE,20000,0.634924051,31499.831,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p2048,2048,DIVIDE,20000,0.104593754,191216.007,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p2048,2048,DIVIDE,20000,0.108730918,183940.321,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p2048,2048,DIVIDE,20000,0.010779107,1855441.274,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p2048,2048,DIVIDE,20000,0.014639620,1366155.682,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p2048,2048,DIVIDE,20000,0.009335711,2142311.385,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p2048,2048,DIVIDE,20000,0.013253862,1508994.142,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p2048,2048,ISQRT,20000,0.014815251,1349960.254,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p2048,2048,ISQRT,20000,0.004824107,4145844.983,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p2048,2048,ISQRT,20000,8.333632041,2399.914,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p2048,2048,ISQRT,20000,8.305543016,2408.030,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p2048,2048,ISQRT,20000,1.547090956,12927.488,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p2048,2048,ISQRT,20000,1.566107197,12770.518,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p2048,2048,ISQRT,20000,0.041968024,476553.292,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p2048,2048,ISQRT,20000,0.045499359,439566.633,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p2048,2048,ISQRT,20000,0.036104160,553952.785,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p2048,2048,ISQRT,20000,0.039222648,509909.479,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-1t,p2048,2048,MODMUL_R2,25000,0.064933560,385008.923,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,gmp-nt,p2048,2048,MODMUL_R2,25000,0.006553631,3814679.196,0
library,Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz,host-cpu,openssl-nt,p2048,2048,MODMUL_R2,25000,0.010982504,2276347.908,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p2048,2048,MODMUL_R2,25000,0.020468648,1221380.129,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w8,p2048,2048,MODMUL_R2,25000,0.024479600,1021258.517,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p2048,2048,MODMUL_R2,25000,0.001129967,22124539.057,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w16,p2048,2048,MODMUL_R2,25000,0.004798240,5210243.725,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p2048,2048,MODMUL_R2,25000,0.000488431,51184306.570,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-opt,p2048,2048,MODMUL_R2,25000,0.004329273,5774641.593,0
opencl-kernel,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p2048,2048,MODMUL_R2,25000,0.000382507,65358283.233,0
opencl-e2e,NVIDIA GeForce RTX 4070 Ti,GPU,w32-o64,p2048,2048,MODMUL_R2,25000,0.004300492,5813288.339,0
```

## CGBN comparison

Merged from `cgbn_results.tsv` after the sweep, by `cgbn_merge.py`.
MPA columns are the fastest correct kernel on the GPU, kernel time only.

### secp256k1 (256-bit)

| Operation | best MPA kernel | MPA items | MPA ops/s | CGBN ops/s | CGBN / MPA |
|---|---|---|---|---|---|
| ADD | w32-o64 | 200000 | 4.05 G | 6.89 G | 1.70x |
| SUBTRACT | w32-o64 | 200000 | 4.09 G | 6.98 G | 1.70x |
| MULTIPLYPRODUCTSCANNING | w32-o64 | 200000 | 1.99 G | 6.98 G | 3.51x |
| MONTGOMERYMULTIPLICATION | w32 | 200000 | 8.42 G | 540.00 M | 0.06x |
| COMPARE | w32-opt | 200000 | 4.22 G | 6.98 G | 1.65x |
| REDUCE | w32-o64 | 25000 | 890.12 M | 3.37 G | 3.78x |
| MODMUL | w32-o64 | 20000 | 317.14 M | 991.43 M | 3.13x |
| MODEXP | w32-o64 | 20000 | 22.34 M | 3.06 M | 0.14x |
| DIVIDE | w32-o64 | 25000 | 348.56 M | 2.79 G | 8.00x |

### rsa256(composite) (256-bit)

| Operation | best MPA kernel | MPA items | MPA ops/s | CGBN ops/s | CGBN / MPA |
|---|---|---|---|---|---|
| ADD | w32-opt | 200000 | 4.06 G | 8.07 G | 1.99x |
| SUBTRACT | w32-o64 | 200000 | 4.08 G | 7.23 G | 1.77x |
| MULTIPLYPRODUCTSCANNING | w32-o64 | 200000 | 1.97 G | 7.21 G | 3.67x |
| MONTGOMERYMULTIPLICATION | w32 | 200000 | 8.08 G | 556.45 M | 0.07x |
| COMPARE | w32-opt | 200000 | 4.20 G | 7.21 G | 1.72x |
| REDUCE | w32-opt | 25000 | 864.96 M | 3.49 G | 4.04x |
| MODMUL | w32-o64 | 20000 | 313.16 M | 1.02 G | 3.25x |
| MODEXP | w32-o64 | 20000 | 22.38 M | 3.13 M | 0.14x |
| DIVIDE | w32-o64 | 25000 | 333.08 M | 2.71 G | 8.14x |

### brainpoolP512r1 (512-bit)

| Operation | best MPA kernel | MPA items | MPA ops/s | CGBN ops/s | CGBN / MPA |
|---|---|---|---|---|---|
| ADD | w32-o64 | 100000 | 1.59 G | 3.08 G | 1.94x |
| SUBTRACT | w32-opt | 100000 | 1.65 G | 3.07 G | 1.86x |
| MULTIPLYPRODUCTSCANNING | w32 | 100000 | 442.91 M | 3.05 G | 6.89x |
| MONTGOMERYMULTIPLICATION | w32 | 100000 | 2.11 G | 170.62 M | 0.08x |
| COMPARE | w32-o64 | 100000 | 1.79 G | 3.10 G | 1.73x |
| REDUCE | w32-opt | 20000 | 294.49 M | 1.98 G | 6.72x |
| MODMUL | w32-o64 | 20000 | 104.30 M | 316.66 M | 3.04x |
| MODEXP | w32-o64 | 20000 | 3.36 M | 1.26 M | 0.38x |
| DIVIDE | w32-o64 | 20000 | 79.12 M | 1.46 G | 18.42x |

### p1024 (1024-bit)

| Operation | best MPA kernel | MPA items | MPA ops/s | CGBN ops/s | CGBN / MPA |
|---|---|---|---|---|---|
| ADD | w32-o64 | 50000 | 634.05 M | 880.78 M | 1.39x |
| SUBTRACT | w32-o64 | 50000 | 642.07 M | 883.77 M | 1.38x |
| MULTIPLYPRODUCTSCANNING | w32-opt | 50000 | 61.39 M | 880.03 M | 14.33x |
| MONTGOMERYMULTIPLICATION | w32 | 50000 | 380.52 M | 66.66 M | 0.18x |
| COMPARE | w32-opt | 50000 | 654.07 M | 885.02 M | 1.35x |
| REDUCE | w32-o64 | 20000 | 86.68 M | 883.77 M | 10.20x |
| MODMUL | w32-o64 | 20000 | 25.51 M | 126.03 M | 4.94x |
| MODEXP | w32-o64 | 20000 | 330.36 k | 247.80 k | 0.75x |
| DIVIDE | w32-o64 | 20000 | 18.31 M | 875.84 M | 47.82x |

### p2048 (2048-bit)

| Operation | best MPA kernel | MPA items | MPA ops/s | CGBN ops/s | CGBN / MPA |
|---|---|---|---|---|---|
| ADD | w32-opt | 25000 | 286.59 M | 452.47 M | 1.58x |
| SUBTRACT | w32-opt | 25000 | 294.16 M | 453.16 M | 1.54x |
| MULTIPLYPRODUCTSCANNING | w32 | 25000 | 16.28 M | 452.11 M | 27.77x |
| MONTGOMERYMULTIPLICATION | w32 | 25000 | 125.10 M | 35.77 M | 0.29x |
| COMPARE | w32-o64 | 25000 | 349.31 M | 453.19 M | 1.30x |
| REDUCE | w32-opt | 20000 | 27.16 M | 453.16 M | 16.68x |
| MODMUL | w32-o64 | 20000 | 6.08 M | 68.53 M | 11.27x |
| MODEXP | w32-opt | 20000 | 24.15 k | 41.83 k | 1.73x |
| DIVIDE | w32-o64 | 20000 | 2.14 M | 452.11 M | 211.04x |
