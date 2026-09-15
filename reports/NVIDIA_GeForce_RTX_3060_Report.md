# MPA-OpenCL benchmark report - NVIDIA GeForce RTX 3060


> **Note.** Two column groups have been removed from this report: the CGBN
> reference column, which was invalid, and the multi-threaded GMP and OpenSSL
> baselines, which predate the 2026-09-12 timing fix and were understated.
> See `reports/README.md`. The single-threaded GMP column, the OpenCL-on-CPU
> rows and every MPA measurement are unaffected, and every configuration was
> verified word-for-word against GMP before it was timed.


## 1. System under test

1 OpenCL device(s) exercised with the identical kernels and operands.

### Device 0 - NVIDIA GeForce RTX 3060 (GPU)

| Property | Value |
|---|---|
| Model | NVIDIA GeForce RTX 3060 |
| Type | GPU |
| Vendor | NVIDIA Corporation |
| Device memory | 11.63 GiB |
| Max single allocation | 2.91 GiB |
| Local memory | 48 KiB |
| Global cache | 784 KiB |
| Compute units | 28 |
| Max clock | 1807 MHz |
| Max work-group size | 1024 |
| OpenCL version | OpenCL 3.0 CUDA |
| Driver | 580.142 |

### Host

| Property | Value |
|---|---|
| CPU | Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz |
| Logical cores | 8 |
| OpenMP threads used | 8 |
| RAM | 31.3 GB |
| OS | Ubuntu 24.04.4 LTS |
| Kernel | 5.15.0-179-generic |
| Arch | x86_64 |
| GMP | 6.3.0 |
| OpenSSL | OpenSSL 3.0.13 30 Jan 2024 |

## 2. Method

- Base workload 200000 items, scaled down per operator by its cost weight and by modulus size. Device rows honour --min-items (20000) so the GPU is not left idle; the CPU libraries keep the smaller count because a full-width MODEXP there costs minutes. Both counts appear in every row as dev/cpu, and throughput is per-second so they remain comparable.
- 5 timed repetitions, **minimum** reported. Two untimed warm-up launches precede them.
- `kernel` times `clEnqueueNDRangeKernel` + `clFinish` only. `e2e` adds the host->device operand writes and the device->host result read.
- Every OpenCL device runs the same kernels on the same operands, so GPU and CPU-OpenCL columns are directly comparable.
- CPU library baselines (GMP, OpenSSL) run those same operands, with temporaries preallocated outside the timed region, so the figure is the arithmetic and not marshalling. The generator is reseeded per modulus and operation so every backend sees identical inputs.
- OpenSSL rows time the nearest BN primitive, which is not always semantically identical (its Montgomery routine expects Montgomery-domain inputs); they measure comparable work, not identical results. Correctness is judged against GMP only.
- Every device cell is checked word-for-word against GMP before it is timed. A cell that mismatches is reported and excluded from the speedup tables.
- Total wall time 3594.9 s.

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

### Device 0 - NVIDIA GeForce RTX 3060 (GPU)

#### secp256k1 (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 200000 / 200000 | 576.19 M | 1.25 G | 2.13 G | 2.17 G | 2.20 G | 2.96 G | 2.99 G | 45.80 M |
| SUBTRACT | 200000 / 200000 | 575.85 M | 1.25 G | 2.14 G | 2.22 G | 2.19 G | 2.96 G | 2.99 G | 63.13 M |
| ADDMOD | 200000 / 200000 | 368.78 M | 808.52 M | 1.70 G | 2.82 G | 2.83 G | 3.08 G | 3.10 G | 17.48 M |
| SUBTRACTMOD | 200000 / 200000 | 368.70 M | 804.66 M | 1.70 G | 2.85 G | 2.85 G | 3.08 G | 3.10 G | 24.52 M |
| MULTIPLYOPERANDSCANNING | 200000 / 200000 | 14.22 M | 54.07 M | 194.81 M | 1.10 G | 1.08 G | 2.35 G | 2.35 G | 42.53 M |
| MULTIPLYPRODUCTSCANNING | 200000 / 200000 | 88.34 M | 279.57 M | 471.39 M | 843.51 M | 975.90 M | 1.44 G | 1.53 G | 42.82 M |
| MONTGOMERYMULTIPLICATION | 200000 / 200000 | 195.84 M | 918.00 M | 2.86 G | 1.88 G | 1.89 G | 1.87 G | 1.91 G | 6.69 M |
| COMPARE | 200000 / 200000 | 670.93 M | 1.15 G | - | 2.31 G | 2.30 G | 6.27 G | 6.34 G | 96.00 M |
| REDUCE | 25000 / 25000 | 88.90 M | 130.85 M | - | 372.78 M | 386.67 M | 371.87 M | 393.36 M | 52.56 M |
| MODMUL | 20000 / 12500 | 35.59 M | 56.25 M | - | 127.10 M | 154.75 M | 127.95 M | 155.97 M | 11.68 M |
| MODEXP | 20000 / 3125 | 784.69 k | 3.81 M | - | 5.91 M | 10.23 M | 5.91 M | 10.53 M | 108.00 k |
| EXPONENTIATION | 20000 / 3125 | 515.07 k | 1.60 M | - | 29.44 M | 35.57 M | 29.53 M | 37.24 M | 308.53 k |
| DIVIDE | 25000 / 25000 | 46.90 M | 73.66 M | - | 195.76 M | 192.29 M | 216.70 M | 222.84 M | 22.67 M |
| ISQRT | 20000 / 6250 | 3.85 M | 5.01 M | - | 20.22 M | 23.35 M | 22.98 M | 25.15 M | 10.67 M |
| MODMUL_R2 | 200000 / 200000 | 172.49 M | 867.27 M | - | 1.26 G | 1.23 G | 1.29 G | 1.31 G | 11.12 M |

#### rsa256(composite) (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 200000 / 200000 | 633.29 M | 1.23 G | 1.97 G | 2.03 G | 2.03 G | 2.92 G | 2.96 G | 48.87 M |
| SUBTRACT | 200000 / 200000 | 634.93 M | 1.24 G | 1.97 G | 2.03 G | 2.04 G | 2.93 G | 2.97 G | 60.99 M |
| ADDMOD | 200000 / 200000 | 447.57 M | 886.72 M | 1.69 G | 2.75 G | 2.84 G | 3.07 G | 3.11 G | 23.05 M |
| SUBTRACTMOD | 200000 / 200000 | 407.06 M | 802.00 M | 1.56 G | 2.76 G | 2.75 G | 3.09 G | 3.09 G | 23.88 M |
| MULTIPLYOPERANDSCANNING | 200000 / 200000 | 14.21 M | 54.13 M | 177.85 M | 1.01 G | 986.10 M | 2.33 G | 2.36 G | 42.94 M |
| MULTIPLYPRODUCTSCANNING | 200000 / 200000 | 88.35 M | 255.93 M | 466.92 M | 848.06 M | 885.29 M | 1.34 G | 1.44 G | 42.36 M |
| MONTGOMERYMULTIPLICATION | 200000 / 200000 | 196.59 M | 844.59 M | 2.79 G | 1.85 G | 1.87 G | 1.87 G | 1.92 G | 6.67 M |
| COMPARE | 200000 / 200000 | 670.25 M | 1.15 G | - | 2.29 G | 2.28 G | 6.30 G | 6.32 G | 96.35 M |
| REDUCE | 25000 / 25000 | 88.76 M | 129.08 M | - | 370.68 M | 383.74 M | 370.38 M | 392.10 M | 32.47 M |
| MODMUL | 20000 / 12500 | 35.53 M | 56.08 M | - | 127.11 M | 154.77 M | 127.61 M | 155.74 M | 11.97 M |
| MODEXP | 20000 / 3125 | 790.86 k | 3.81 M | - | 5.92 M | 10.24 M | 5.91 M | 10.56 M | 113.62 k |
| EXPONENTIATION | 20000 / 3125 | 513.82 k | 1.45 M | - | 29.45 M | 35.55 M | 29.53 M | 37.24 M | 313.66 k |
| DIVIDE | 25000 / 25000 | 46.77 M | 73.63 M | - | 191.91 M | 187.79 M | 207.84 M | 212.95 M | 22.59 M |
| ISQRT | 20000 / 6250 | 3.83 M | 5.01 M | - | 19.97 M | 23.15 M | 22.86 M | 25.19 M | 12.12 M |
| MODMUL_R2 | 200000 / 200000 | 170.53 M | 855.06 M | - | 1.27 G | 1.23 G | 1.29 G | 1.31 G | 11.44 M |

#### brainpoolP512r1 (512-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 100000 / 100000 | 260.41 M | 496.09 M | 257.53 M | 827.25 M | 813.52 M | 1.32 G | 955.48 M | 42.82 M |
| SUBTRACT | 100000 / 100000 | 260.37 M | 493.52 M | 259.01 M | 832.29 M | 824.03 M | 1.31 G | 951.22 M | 51.81 M |
| ADDMOD | 100000 / 100000 | 189.33 M | 362.62 M | 346.11 M | 771.93 M | 767.87 M | 1.15 G | 1.36 G | 19.84 M |
| SUBTRACTMOD | 100000 / 100000 | 166.00 M | 323.25 M | 321.25 M | 743.22 M | 749.10 M | 1.16 G | 1.37 G | 22.11 M |
| MULTIPLYOPERANDSCANNING | 100000 / 100000 | 2.68 M | 10.23 M | 10.14 M | 219.79 M | 219.38 M | 394.05 M | 427.15 M | 20.46 M |
| MULTIPLYPRODUCTSCANNING | 100000 / 100000 | 12.03 M | 44.12 M | 76.79 M | 147.11 M | 147.16 M | 317.07 M | 220.99 M | 20.44 M |
| MONTGOMERYMULTIPLICATION | 100000 / 100000 | 49.84 M | 172.77 M | 971.95 M | 564.09 M | 726.86 M | 613.90 M | 670.44 M | 3.02 M |
| COMPARE | 100000 / 100000 | 263.75 M | 502.33 M | - | 926.98 M | 926.78 M | 3.44 G | 3.30 G | 78.32 M |
| REDUCE | 20000 / 12500 | 32.53 M | 40.61 M | - | 140.27 M | 146.84 M | 150.90 M | 127.57 M | 31.01 M |
| MODMUL | 20000 / 6250 | 11.36 M | 15.91 M | - | 37.77 M | 46.80 M | 37.55 M | 38.96 M | 6.12 M |
| MODEXP | 20000 / 1562 | 54.21 k | 588.68 k | - | 774.66 k | 1.51 M | 792.94 k | 1.42 M | 20.36 k |
| EXPONENTIATION | 20000 / 1562 | 64.21 k | 242.64 k | - | 909.85 k | 925.06 k | 887.22 k | 902.81 k | 92.96 k |
| DIVIDE | 20000 / 12500 | 16.56 M | 18.33 M | - | 54.47 M | 56.27 M | 56.81 M | 57.29 M | 21.25 M |
| ISQRT | 20000 / 3125 | 783.32 k | 788.61 k | - | 3.85 M | 4.08 M | 4.05 M | 3.35 M | 5.92 M |
| MODMUL_R2 | 100000 / 100000 | 33.55 M | 192.69 M | - | 333.76 M | 523.67 M | 352.23 M | 518.05 M | 5.95 M |

#### p1024 (1024-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 50000 / 50000 | 81.33 M | 157.31 M | 86.55 M | 316.11 M | 318.26 M | 476.19 M | 472.37 M | 30.96 M |
| SUBTRACT | 50000 / 50000 | 81.59 M | 156.56 M | 86.11 M | 319.23 M | 315.04 M | 476.02 M | 471.03 M | 34.46 M |
| ADDMOD | 50000 / 50000 | 55.57 M | 107.55 M | 97.56 M | 226.78 M | 223.25 M | 568.34 M | 566.26 M | 13.29 M |
| SUBTRACTMOD | 50000 / 50000 | 50.82 M | 108.32 M | 97.55 M | 223.25 M | 225.22 M | 571.96 M | 570.26 M | 17.86 M |
| MULTIPLYOPERANDSCANNING | 50000 / 50000 | 546.34 k | 2.05 M | 2.00 M | 88.08 M | 84.89 M | 129.80 M | 135.66 M | 6.50 M |
| MULTIPLYPRODUCTSCANNING | 50000 / 50000 | 1.49 M | 5.78 M | 22.12 M | 22.82 M | 22.76 M | 56.72 M | 60.76 M | 6.48 M |
| MONTGOMERYMULTIPLICATION | 50000 / 50000 | 8.66 M | 44.21 M | 193.49 M | 118.16 M | 147.84 M | 160.41 M | 235.55 M | 988.67 k |
| COMPARE | 50000 / 50000 | 94.38 M | 181.88 M | - | 320.38 M | 317.02 M | 1.53 G | 1.57 G | 47.38 M |
| REDUCE | 20000 / 6250 | 6.03 M | 11.05 M | - | 38.18 M | 37.32 M | 44.69 M | 41.54 M | 43.62 M |
| MODMUL | 20000 / 3125 | 1.49 M | 4.16 M | - | 8.59 M | 11.49 M | 8.84 M | 11.84 M | 2.25 M |
| MODEXP | 20000 / 781 | 6.66 k | 44.14 k | - | 100.86 k | 162.64 k | 99.35 k | 158.34 k | 3.11 k |
| EXPONENTIATION | 20000 / 781 | 8.30 k | 31.48 k | - | 114.47 k | 119.54 k | 113.13 k | 118.07 k | 23.63 k |
| DIVIDE | 20000 / 6250 | 368.06 k | 2.90 M | - | 12.86 M | 13.11 M | 12.98 M | 13.04 M | 20.03 M |
| ISQRT | 20000 / 1562 | 30.11 k | 111.46 k | - | 576.41 k | 607.59 k | 570.12 k | 615.17 k | 3.11 M |
| MODMUL_R2 | 50000 / 50000 | 5.35 M | 39.14 M | - | 71.81 M | 88.16 M | 82.43 M | 104.68 M | 2.20 M |

#### p2048 (2048-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 25000 / 25000 | 39.34 M | 77.33 M | 153.21 M | 155.26 M | 156.16 M | 233.44 M | 239.43 M | 20.04 M |
| SUBTRACT | 25000 / 25000 | 39.09 M | 77.38 M | 153.14 M | 156.03 M | 156.02 M | 233.63 M | 238.79 M | 20.32 M |
| ADDMOD | 25000 / 25000 | 29.00 M | 57.27 M | 122.39 M | 109.56 M | 114.06 M | 251.33 M | 253.65 M | 10.89 M |
| SUBTRACTMOD | 25000 / 25000 | 24.72 M | 52.87 M | 109.45 M | 111.64 M | 113.37 M | 234.09 M | 254.42 M | 12.41 M |
| MULTIPLYOPERANDSCANNING | 25000 / 25000 | 123.29 k | 306.00 k | 1.39 M | 21.52 M | 21.42 M | 20.40 M | 20.45 M | 1.99 M |
| MULTIPLYPRODUCTSCANNING | 25000 / 25000 | 376.21 k | 1.47 M | 5.87 M | 5.16 M | 5.69 M | 9.75 M | 11.27 M | 1.99 M |
| MONTGOMERYMULTIPLICATION | 25000 / 25000 | 261.56 k | 9.48 M | 52.48 M | 31.29 M | 37.00 M | 35.41 M | 41.38 M | 315.04 k |
| COMPARE | 25000 / 25000 | 45.31 M | 91.03 M | - | 155.98 M | 159.67 M | 580.63 M | 586.40 M | 70.95 M |
| REDUCE | 20000 / 3125 | 39.34 k | 2.55 M | - | 12.55 M | 11.14 M | 13.94 M | 11.64 M | 30.47 M |
| MODMUL | 20000 / 1562 | 23.59 k | 933.70 k | - | 2.16 M | 2.76 M | 2.40 M | 2.85 M | 753.00 k |
| MODEXP | 20000 / 390 | 151.7 | 1.33 k | - | 11.42 k | 4.11 k | 12.44 k | 3.73 k | 434.9 |
| EXPONENTIATION | 20000 / 390 | 692.6 | 3.84 k | - | 13.46 k | 13.88 k | 13.44 k | 14.04 k | 3.80 k |
| DIVIDE | 20000 / 3125 | 14.94 k | 65.16 k | - | 762.50 k | 1.05 M | 900.92 k | 870.06 k | 16.55 M |
| ISQRT | 20000 / 781 | 1.16 k | 3.50 k | - | 93.74 k | 262.31 k | 92.79 k | 264.09 k | 1.92 M |
| MODMUL_R2 | 25000 / 25000 | 334.58 k | 8.01 M | - | 19.20 M | 24.39 M | 21.12 M | 25.75 M | 736.53 k |

## 5. Head to head

Best OpenCL GPU result against best OpenCL CPU result and the CPU libraries.
Ratios above 1.00x mean the GPU is faster than that baseline.

### secp256k1 (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w32-il64 | 2.99 G | none | n/a | 45.80 M | n/a |
| SUBTRACT | w32-il64 | 2.99 G | none | n/a | 63.13 M | n/a |
| ADDMOD | w32-il64 | 3.10 G | none | n/a | 17.48 M | n/a |
| SUBTRACTMOD | w32-il64 | 3.10 G | none | n/a | 24.52 M | n/a |
| MULTIPLYOPERANDSCANNING | w32-il64 | 2.35 G | none | n/a | 42.53 M | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 1.53 G | none | n/a | 42.82 M | n/a |
| MONTGOMERYMULTIPLICATION | w32 | 2.86 G | none | n/a | 6.69 M | n/a |
| COMPARE | w32-il64 | 6.34 G | none | n/a | 96.00 M | n/a |
| REDUCE | w32-il64 | 393.36 M | none | n/a | 52.56 M | n/a |
| MODMUL | w32-il64 | 97.48 M | none | n/a | 11.68 M | n/a |
| MODEXP | w32-il64 | 1.64 M | none | n/a | 108.00 k | n/a |
| EXPONENTIATION | w32-il64 | 5.82 M | none | n/a | 308.53 k | n/a |
| DIVIDE | w32-il64 | 222.84 M | none | n/a | 22.67 M | n/a |
| ISQRT | w32-il64 | 7.86 M | none | n/a | 10.67 M | n/a |
| MODMUL_R2 | w32-il64 | 1.31 G | none | n/a | 11.12 M | n/a |

### rsa256(composite) (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w32-il64 | 2.96 G | none | n/a | 48.87 M | n/a |
| SUBTRACT | w32-il64 | 2.97 G | none | n/a | 60.99 M | n/a |
| ADDMOD | w32-il64 | 3.11 G | none | n/a | 23.05 M | n/a |
| SUBTRACTMOD | w32-il64 | 3.09 G | none | n/a | 23.88 M | n/a |
| MULTIPLYOPERANDSCANNING | w32-il64 | 2.36 G | none | n/a | 42.94 M | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 1.44 G | none | n/a | 42.36 M | n/a |
| MONTGOMERYMULTIPLICATION | w32 | 2.79 G | none | n/a | 6.67 M | n/a |
| COMPARE | w32-il64 | 6.32 G | none | n/a | 96.35 M | n/a |
| REDUCE | w32-il64 | 392.10 M | none | n/a | 32.47 M | n/a |
| MODMUL | w32-il64 | 97.34 M | none | n/a | 11.97 M | n/a |
| MODEXP | w32-il64 | 1.65 M | none | n/a | 113.62 k | n/a |
| EXPONENTIATION | w32-il64 | 5.82 M | none | n/a | 313.66 k | n/a |
| DIVIDE | w32-il64 | 212.95 M | none | n/a | 22.59 M | n/a |
| ISQRT | w32-il64 | 7.87 M | none | n/a | 12.12 M | n/a |
| MODMUL_R2 | w32-il64 | 1.31 G | none | n/a | 11.44 M | n/a |

### brainpoolP512r1 (512-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w32-il | 1.32 G | none | n/a | 42.82 M | n/a |
| SUBTRACT | w32-il | 1.31 G | none | n/a | 51.81 M | n/a |
| ADDMOD | w32-il64 | 1.36 G | none | n/a | 19.84 M | n/a |
| SUBTRACTMOD | w32-il64 | 1.37 G | none | n/a | 22.11 M | n/a |
| MULTIPLYOPERANDSCANNING | w32-il64 | 427.15 M | none | n/a | 20.46 M | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il | 317.07 M | none | n/a | 20.44 M | n/a |
| MONTGOMERYMULTIPLICATION | w32 | 971.95 M | none | n/a | 3.02 M | n/a |
| COMPARE | w32-il | 3.44 G | none | n/a | 78.32 M | n/a |
| REDUCE | w32-il | 94.31 M | none | n/a | 31.01 M | n/a |
| MODMUL | w32-o64 | 14.62 M | none | n/a | 6.12 M | n/a |
| MODEXP | w32-o64 | 118.23 k | none | n/a | 20.36 k | n/a |
| EXPONENTIATION | w32-o64 | 72.25 k | none | n/a | 92.96 k | n/a |
| DIVIDE | w32-il64 | 35.81 M | none | n/a | 21.25 M | n/a |
| ISQRT | w32-o64 | 638.09 k | none | n/a | 5.92 M | n/a |
| MODMUL_R2 | w32-o64 | 523.67 M | none | n/a | 5.95 M | n/a |

### p1024 (1024-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w32-il | 476.19 M | none | n/a | 30.96 M | n/a |
| SUBTRACT | w32-il | 476.02 M | none | n/a | 34.46 M | n/a |
| ADDMOD | w32-il | 568.34 M | none | n/a | 13.29 M | n/a |
| SUBTRACTMOD | w32-il | 571.96 M | none | n/a | 17.86 M | n/a |
| MULTIPLYOPERANDSCANNING | w32-il64 | 135.66 M | none | n/a | 6.50 M | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 60.76 M | none | n/a | 6.48 M | n/a |
| MONTGOMERYMULTIPLICATION | w32-il64 | 235.55 M | none | n/a | 988.67 k | n/a |
| COMPARE | w32-il64 | 1.57 G | none | n/a | 47.38 M | n/a |
| REDUCE | w32-il | 13.96 M | none | n/a | 43.62 M | n/a |
| MODMUL | w32-il64 | 1.85 M | none | n/a | 2.25 M | n/a |
| MODEXP | w32-o64 | 6.35 k | none | n/a | 3.11 k | n/a |
| EXPONENTIATION | w32-o64 | 4.67 k | none | n/a | 23.63 k | n/a |
| DIVIDE | w32-o64 | 4.10 M | none | n/a | 20.03 M | n/a |
| ISQRT | w32-il64 | 48.04 k | none | n/a | 3.11 M | n/a |
| MODMUL_R2 | w32-il64 | 104.68 M | none | n/a | 2.20 M | n/a |

### p2048 (2048-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w32-il64 | 239.43 M | none | n/a | 20.04 M | n/a |
| SUBTRACT | w32-il64 | 238.79 M | none | n/a | 20.32 M | n/a |
| ADDMOD | w32-il64 | 253.65 M | none | n/a | 10.89 M | n/a |
| SUBTRACTMOD | w32-il64 | 254.42 M | none | n/a | 12.41 M | n/a |
| MULTIPLYOPERANDSCANNING | w32-opt | 21.52 M | none | n/a | 1.99 M | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 11.27 M | none | n/a | 1.99 M | n/a |
| MONTGOMERYMULTIPLICATION | w32 | 52.48 M | none | n/a | 315.04 k | n/a |
| COMPARE | w32-il64 | 586.40 M | none | n/a | 70.95 M | n/a |
| REDUCE | w32-il | 2.18 M | none | n/a | 30.47 M | n/a |
| MODMUL | w32-il64 | 222.65 k | none | n/a | 753.00 k | n/a |
| MODEXP | w32-il | 242.7 | none | n/a | 434.9 | n/a |
| EXPONENTIATION | w32-il64 | 273.8 | none | n/a | 3.80 k | n/a |
| DIVIDE | w32-o64 | 164.72 k | none | n/a | 16.55 M | n/a |
| ISQRT | w32-il64 | 10.31 k | none | n/a | 1.92 M | n/a |
| MODMUL_R2 | w32-il64 | 25.75 M | none | n/a | 736.53 k | n/a |

## 6. Raw data

Also written to `NVIDIA_GeForce_RTX_3060_Report.csv` for analysis.

```csv
kind,device,device_type,kernel,modulus,bits,operation,items,seconds,ops_per_sec,mismatches
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,secp256k1,256,ADD,200000,0.004366861,45799488.297,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,secp256k1,256,ADD,200000,0.000907999,220264560.316,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,secp256k1,256,ADD,200000,0.001075022,186042707.566,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,secp256k1,256,ADD,200000,0.000122816,1628452318.916,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,secp256k1,256,ADD,200000,0.000347108,576189546.474,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,secp256k1,256,ADD,200000,0.002409690,82998227.936,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,secp256k1,256,ADD,200000,0.000160237,1248151039.341,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,secp256k1,256,ADD,200000,0.002459190,81327592.565,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32,secp256k1,256,ADD,200000,0.000093831,2131492128.307,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32,secp256k1,256,ADD,200000,0.002255584,88668830.385,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,secp256k1,256,ADD,200000,0.000091960,2174858317.789,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,secp256k1,256,ADD,200000,0.002252545,88788460.265,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,secp256k1,256,ADD,200000,0.000090805,2202521917.568,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,secp256k1,256,ADD,200000,0.002344084,85321175.506,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,secp256k1,256,ADD,200000,0.000067497,2963093711.398,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,secp256k1,256,ADD,200000,0.002255231,88682710.431,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,secp256k1,256,ADD,200000,0.000066810,2993563994.713,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,secp256k1,256,ADD,200000,0.002332720,85736823.198,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,secp256k1,256,SUBTRACT,200000,0.003168166,63128004.630,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,secp256k1,256,SUBTRACT,200000,0.000810506,246759432.699,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,secp256k1,256,SUBTRACT,200000,0.001071842,186594669.510,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,secp256k1,256,SUBTRACT,200000,0.000120832,1655190677.966,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,secp256k1,256,SUBTRACT,200000,0.000347310,575854459.353,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,secp256k1,256,SUBTRACT,200000,0.002361027,84708899.983,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,secp256k1,256,SUBTRACT,200000,0.000160092,1249281489.260,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,secp256k1,256,SUBTRACT,200000,0.002444023,81832289.324,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32,secp256k1,256,SUBTRACT,200000,0.000093612,2136478115.110,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32,secp256k1,256,SUBTRACT,200000,0.002273603,87966105.162,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,secp256k1,256,SUBTRACT,200000,0.000090090,2220002814.954,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,secp256k1,256,SUBTRACT,200000,0.002203359,90770501.782,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,secp256k1,256,SUBTRACT,200000,0.000091120,2194907365.261,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,secp256k1,256,SUBTRACT,200000,0.002386834,83793007.571,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,secp256k1,256,SUBTRACT,200000,0.000067626,2957441979.814,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,secp256k1,256,SUBTRACT,200000,0.002242194,89198346.775,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,secp256k1,256,SUBTRACT,200000,0.000066821,2993071141.557,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,secp256k1,256,SUBTRACT,200000,0.002271236,88057780.819,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,secp256k1,256,ADDMOD,200000,0.011438684,17484528.882,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,secp256k1,256,ADDMOD,200000,0.002809057,71198271.618,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,secp256k1,256,ADDMOD,200000,0.007301690,27390919.093,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,secp256k1,256,ADDMOD,200000,0.000122400,1633986928.105,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,secp256k1,256,ADDMOD,200000,0.000542325,368782538.980,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,secp256k1,256,ADDMOD,200000,0.002579973,77520192.219,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,secp256k1,256,ADDMOD,200000,0.000247366,808518460.245,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,secp256k1,256,ADDMOD,200000,0.002478187,80704160.834,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32,secp256k1,256,ADDMOD,200000,0.000117400,1703577210.564,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32,secp256k1,256,ADDMOD,200000,0.002260358,88481560.822,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,secp256k1,256,ADDMOD,200000,0.000070803,2824737387.690,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,secp256k1,256,ADDMOD,200000,0.002250167,88882291.163,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,secp256k1,256,ADDMOD,200000,0.000070554,2834706975.033,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,secp256k1,256,ADDMOD,200000,0.002318076,86278449.111,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,secp256k1,256,ADDMOD,200000,0.000064839,3084562023.689,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,secp256k1,256,ADDMOD,200000,0.002201528,90845994.345,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,secp256k1,256,ADDMOD,200000,0.000064487,3101400909.486,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,secp256k1,256,ADDMOD,200000,0.002277766,87805331.869,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,secp256k1,256,SUBTRACTMOD,200000,0.008157344,24517784.187,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,secp256k1,256,SUBTRACTMOD,200000,0.002156957,92723220.962,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,secp256k1,256,SUBTRACTMOD,200000,0.007274933,27491662.180,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,secp256k1,256,SUBTRACTMOD,200000,0.000122592,1631427825.633,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,secp256k1,256,SUBTRACTMOD,200000,0.000542450,368697577.351,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,secp256k1,256,SUBTRACTMOD,200000,0.002594872,77075091.873,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,secp256k1,256,SUBTRACTMOD,200000,0.000248553,804657361.547,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,secp256k1,256,SUBTRACTMOD,200000,0.002298041,87030648.618,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32,secp256k1,256,SUBTRACTMOD,200000,0.000117746,1698572237.723,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32,secp256k1,256,SUBTRACTMOD,200000,0.002312732,86477811.132,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,secp256k1,256,SUBTRACTMOD,200000,0.000070082,2853798867.774,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,secp256k1,256,SUBTRACTMOD,200000,0.002221887,90013578.552,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,secp256k1,256,SUBTRACTMOD,200000,0.000070089,2853514464.339,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,secp256k1,256,SUBTRACTMOD,200000,0.002255218,88683223.149,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,secp256k1,256,SUBTRACTMOD,200000,0.000064974,3078156600.887,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,secp256k1,256,SUBTRACTMOD,200000,0.002208726,90549936.319,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,secp256k1,256,SUBTRACTMOD,200000,0.000064421,3104578689.200,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,secp256k1,256,SUBTRACTMOD,200000,0.002130825,93860358.620,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,secp256k1,256,MULTIPLYOPERANDSCANNING,200000,0.004702380,42531653.791,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,200000,0.001204632,166025808.454,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,200000,0.002265832,88267798.236,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,200000,0.014064780,14219916.729,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,200000,0.016756898,11935383.249,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,secp256k1,256,MULTIPLYOPERANDSCANNING,200000,0.003699218,54065480.291,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,secp256k1,256,MULTIPLYOPERANDSCANNING,200000,0.006445546,31029178.937,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32,secp256k1,256,MULTIPLYOPERANDSCANNING,200000,0.001026622,194813666.725,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32,secp256k1,256,MULTIPLYOPERANDSCANNING,200000,0.003893701,51365012.225,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,secp256k1,256,MULTIPLYOPERANDSCANNING,200000,0.000181924,1099360004.044,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,secp256k1,256,MULTIPLYOPERANDSCANNING,200000,0.003031405,65976008.260,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,secp256k1,256,MULTIPLYOPERANDSCANNING,200000,0.000185576,1077725704.424,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,secp256k1,256,MULTIPLYOPERANDSCANNING,200000,0.003009357,66459379.331,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,secp256k1,256,MULTIPLYOPERANDSCANNING,200000,0.000085254,2345931891.535,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,secp256k1,256,MULTIPLYOPERANDSCANNING,200000,0.002859669,69938163.218,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,secp256k1,256,MULTIPLYOPERANDSCANNING,200000,0.000085208,2347197916.744,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,secp256k1,256,MULTIPLYOPERANDSCANNING,200000,0.002772566,72135343.495,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,secp256k1,256,MULTIPLYPRODUCTSCANNING,200000,0.004670507,42821903.640,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,200000,0.001178313,169734189.434,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,200000,0.002260206,88487510.477,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,secp256k1,256,MULTIPLYPRODUCTSCANNING,200000,0.000121856,1641281512.605,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,200000,0.002264008,88338911.080,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,200000,0.004957909,40339586.731,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,secp256k1,256,MULTIPLYPRODUCTSCANNING,200000,0.000715378,279572488.165,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,secp256k1,256,MULTIPLYPRODUCTSCANNING,200000,0.003454562,57894459.800,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32,secp256k1,256,MULTIPLYPRODUCTSCANNING,200000,0.000424279,471387912.787,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32,secp256k1,256,MULTIPLYPRODUCTSCANNING,200000,0.003237449,61777035.287,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,secp256k1,256,MULTIPLYPRODUCTSCANNING,200000,0.000237105,843508126.143,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,secp256k1,256,MULTIPLYPRODUCTSCANNING,200000,0.002988245,66928915.415,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,secp256k1,256,MULTIPLYPRODUCTSCANNING,200000,0.000204938,975904916.371,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,secp256k1,256,MULTIPLYPRODUCTSCANNING,200000,0.002965001,67453601.969,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,secp256k1,256,MULTIPLYPRODUCTSCANNING,200000,0.000139276,1435997863.874,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,secp256k1,256,MULTIPLYPRODUCTSCANNING,200000,0.002925696,68359803.000,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,secp256k1,256,MULTIPLYPRODUCTSCANNING,200000,0.000130854,1528421397.707,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,secp256k1,256,MULTIPLYPRODUCTSCANNING,200000,0.002806165,71271647.435,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,secp256k1,256,MONTGOMERYMULTIPLICATION,200000,0.029909195,6686906.813,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,200000,0.008418457,23757322.748,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,200000,0.001948040,102667294.412,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,secp256k1,256,MONTGOMERYMULTIPLICATION,200000,0.000124928,1600922131.148,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,200000,0.001021243,195839777.265,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,200000,0.003044759,65686644.234,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,secp256k1,256,MONTGOMERYMULTIPLICATION,200000,0.000217865,917999727.697,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,secp256k1,256,MONTGOMERYMULTIPLICATION,200000,0.002234981,89486219.832,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32,secp256k1,256,MONTGOMERYMULTIPLICATION,200000,0.000069832,2864015601.169,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32,secp256k1,256,MONTGOMERYMULTIPLICATION,200000,0.002181748,91669615.941,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,secp256k1,256,MONTGOMERYMULTIPLICATION,200000,0.000106481,1878269003.843,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,secp256k1,256,MONTGOMERYMULTIPLICATION,200000,0.002257802,88581726.838,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,secp256k1,256,MONTGOMERYMULTIPLICATION,200000,0.000105811,1890161752.429,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,secp256k1,256,MONTGOMERYMULTIPLICATION,200000,0.002280146,87713681.493,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,secp256k1,256,MONTGOMERYMULTIPLICATION,200000,0.000106926,1870453167.437,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,secp256k1,256,MONTGOMERYMULTIPLICATION,200000,0.002280329,87706642.144,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,secp256k1,256,MONTGOMERYMULTIPLICATION,200000,0.000104527,1913381430.035,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,secp256k1,256,MONTGOMERYMULTIPLICATION,200000,0.002116481,94496479.726,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,secp256k1,256,COMPARE,200000,0.002083362,95998678.047,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,secp256k1,256,COMPARE,200000,0.000757287,264100656.383,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,secp256k1,256,COMPARE,200000,0.001060844,188529134.903,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,secp256k1,256,COMPARE,200000,0.000121760,1642575558.476,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,secp256k1,256,COMPARE,200000,0.000298092,670933879.666,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,secp256k1,256,COMPARE,200000,0.002383493,83910462.381,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,secp256k1,256,COMPARE,200000,0.000173668,1151622837.565,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,secp256k1,256,COMPARE,200000,0.002426730,82415431.356,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,secp256k1,256,COMPARE,200000,0.000086743,2305662159.407,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,secp256k1,256,COMPARE,200000,0.002277928,87799087.705,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,secp256k1,256,COMPARE,200000,0.000087077,2296818478.968,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,secp256k1,256,COMPARE,200000,0.002425455,82458756.372,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,secp256k1,256,COMPARE,200000,0.000031886,6272348765.692,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,secp256k1,256,COMPARE,200000,0.002223512,89947794.981,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,secp256k1,256,COMPARE,200000,0.000031560,6337131774.613,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,secp256k1,256,COMPARE,200000,0.002058033,97180170.773,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,secp256k1,256,REDUCE,25000,0.000475638,52560985.232,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,secp256k1,256,REDUCE,25000,0.000130562,191479916.791,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,secp256k1,256,REDUCE,25000,0.000738905,33833846.708,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,secp256k1,256,REDUCE,200000,0.000153472,1303169307.756,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,secp256k1,256,REDUCE,25000,0.000281202,88904054.014,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,secp256k1,256,REDUCE,25000,0.000677511,36899771.005,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,secp256k1,256,REDUCE,25000,0.000191059,130849617.944,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,secp256k1,256,REDUCE,25000,0.000618439,40424355.394,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,secp256k1,256,REDUCE,25000,0.000067064,372778483.357,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,secp256k1,256,REDUCE,25000,0.000473524,52795636.871,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,secp256k1,256,REDUCE,25000,0.000064654,386673769.634,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,secp256k1,256,REDUCE,25000,0.000497749,50226116.259,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,secp256k1,256,REDUCE,25000,0.000067227,371874416.837,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,secp256k1,256,REDUCE,25000,0.000470674,53115324.648,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,secp256k1,256,REDUCE,25000,0.000063555,393359932.666,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,secp256k1,256,REDUCE,25000,0.000497000,50301810.774,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,secp256k1,256,MODMUL,12500,0.001069981,11682450.188,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,secp256k1,256,MODMUL,12500,0.000315646,39601325.451,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,secp256k1,256,MODMUL,12500,0.000970383,12881511.732,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,secp256k1,256,MODMUL,200000,0.000567296,352549638.989,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,secp256k1,256,MODMUL,20000,0.000561899,35593586.484,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,secp256k1,256,MODMUL,20000,0.000895281,22339354.126,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,secp256k1,256,MODMUL,20000,0.000355536,56253093.646,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,secp256k1,256,MODMUL,20000,0.000721138,27733943.868,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,secp256k1,256,MODMUL,20000,0.000157353,127102754.141,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,secp256k1,256,MODMUL,20000,0.000503378,39731574.363,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,secp256k1,256,MODMUL,20000,0.000129239,154752056.439,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,secp256k1,256,MODMUL,20000,0.000496013,40321525.977,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,secp256k1,256,MODMUL,20000,0.000156309,127951693.652,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,secp256k1,256,MODMUL,20000,0.000466003,42918176.383,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,secp256k1,256,MODMUL,20000,0.000128229,155970918.955,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,secp256k1,256,MODMUL,20000,0.000493715,40509201.198,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,secp256k1,256,MODEXP,3125,0.028935284,107999.631,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,secp256k1,256,MODEXP,3125,0.007131183,438216.212,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,secp256k1,256,MODEXP,3125,0.012566984,248667.460,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,secp256k1,256,MODEXP,200000,0.193794057,1032023.392,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,secp256k1,256,MODEXP,20000,0.025487894,784686.251,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,secp256k1,256,MODEXP,20000,0.025856310,773505.576,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,secp256k1,256,MODEXP,20000,0.005250946,3808837.509,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,secp256k1,256,MODEXP,20000,0.005632261,3550971.810,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,secp256k1,256,MODEXP,20000,0.003384430,5909414.550,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,secp256k1,256,MODEXP,20000,0.003758437,5321360.980,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,secp256k1,256,MODEXP,20000,0.001955291,10228656.683,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,secp256k1,256,MODEXP,20000,0.002332162,8575733.606,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,secp256k1,256,MODEXP,20000,0.003386745,5905375.226,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,secp256k1,256,MODEXP,20000,0.003731695,5359494.809,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,secp256k1,256,MODEXP,20000,0.001900006,10526282.719,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,secp256k1,256,MODEXP,20000,0.002286761,8745995.067,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,secp256k1,256,EXPONENTIATION,3125,0.010128817,308525.665,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,secp256k1,256,EXPONENTIATION,3125,0.002200888,1419881.427,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,secp256k1,256,EXPONENTIATION,3125,0.040668028,76841.690,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,secp256k1,256,EXPONENTIATION,20000,0.038829325,515074.625,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,secp256k1,256,EXPONENTIATION,20000,0.039345256,508320.495,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,secp256k1,256,EXPONENTIATION,20000,0.012526571,1596606.127,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,secp256k1,256,EXPONENTIATION,20000,0.014187550,1409686.662,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,secp256k1,256,EXPONENTIATION,20000,0.000679266,29443546.470,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,secp256k1,256,EXPONENTIATION,20000,0.001004146,19917422.276,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,secp256k1,256,EXPONENTIATION,20000,0.000562283,35569278.808,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,secp256k1,256,EXPONENTIATION,20000,0.000935925,21369233.734,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,secp256k1,256,EXPONENTIATION,20000,0.000677385,29525306.085,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,secp256k1,256,EXPONENTIATION,20000,0.001024351,19524557.698,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,secp256k1,256,EXPONENTIATION,20000,0.000537127,37235144.184,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,secp256k1,256,EXPONENTIATION,20000,0.000905107,22096834.757,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,secp256k1,256,DIVIDE,25000,0.001102810,22669363.449,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,secp256k1,256,DIVIDE,25000,0.000275132,90865471.027,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,secp256k1,256,DIVIDE,25000,0.000725770,34446173.462,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,secp256k1,256,DIVIDE,200000,0.000203776,981469849.246,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,secp256k1,256,DIVIDE,25000,0.000533043,46900531.525,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,secp256k1,256,DIVIDE,25000,0.001022504,24449782.048,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,secp256k1,256,DIVIDE,25000,0.000339381,73663525.026,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,secp256k1,256,DIVIDE,25000,0.000872069,28667455.096,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,secp256k1,256,DIVIDE,25000,0.000127705,195763691.416,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,secp256k1,256,DIVIDE,25000,0.000636246,39292978.049,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,secp256k1,256,DIVIDE,25000,0.000130012,192289955.583,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,secp256k1,256,DIVIDE,25000,0.000649559,38487653.085,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,secp256k1,256,DIVIDE,25000,0.000115369,216695919.375,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,secp256k1,256,DIVIDE,25000,0.000615687,40605048.149,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,secp256k1,256,DIVIDE,25000,0.000112188,222840251.430,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,secp256k1,256,DIVIDE,25000,0.000627514,39839746.493,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,secp256k1,256,ISQRT,6250,0.000585835,10668532.592,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,secp256k1,256,ISQRT,6250,0.000115097,54302034.591,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,secp256k1,256,ISQRT,20000,0.005201506,3845040.291,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,secp256k1,256,ISQRT,20000,0.005557511,3598733.317,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,secp256k1,256,ISQRT,20000,0.003990594,5011785.205,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,secp256k1,256,ISQRT,20000,0.004365573,4581300.098,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,secp256k1,256,ISQRT,20000,0.000988967,20223121.330,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,secp256k1,256,ISQRT,20000,0.001337328,14955194.301,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,secp256k1,256,ISQRT,20000,0.000856697,23345477.491,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,secp256k1,256,ISQRT,20000,0.001236195,16178677.606,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,secp256k1,256,ISQRT,20000,0.000870436,22976991.783,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,secp256k1,256,ISQRT,20000,0.001227856,16288554.590,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,secp256k1,256,ISQRT,20000,0.000795196,25151030.576,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,secp256k1,256,ISQRT,20000,0.001171486,17072333.119,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,secp256k1,256,MODMUL_R2,200000,0.017983020,11121602.494,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,secp256k1,256,MODMUL_R2,200000,0.005059747,39527667.885,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,secp256k1,256,MODMUL_R2,200000,0.015383292,13001118.316,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,secp256k1,256,MODMUL_R2,200000,0.000168960,1183712121.212,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,secp256k1,256,MODMUL_R2,200000,0.001159490,172489624.906,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,secp256k1,256,MODMUL_R2,200000,0.003275058,61067620.798,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,secp256k1,256,MODMUL_R2,200000,0.000230609,867268939.420,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,secp256k1,256,MODMUL_R2,200000,0.002488379,80373608.119,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,secp256k1,256,MODMUL_R2,200000,0.000158480,1261989204.254,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,secp256k1,256,MODMUL_R2,200000,0.002322002,86132569.231,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,secp256k1,256,MODMUL_R2,200000,0.000162925,1227558774.575,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,secp256k1,256,MODMUL_R2,200000,0.002484826,80488532.494,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,secp256k1,256,MODMUL_R2,200000,0.000154730,1292574238.679,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,secp256k1,256,MODMUL_R2,200000,0.002369417,84408947.682,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,secp256k1,256,MODMUL_R2,200000,0.000152569,1310882020.025,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,secp256k1,256,MODMUL_R2,200000,0.002192421,91223355.311,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,rsa256(composite),256,ADD,200000,0.004092756,48866826.924,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,rsa256(composite),256,ADD,200000,0.000911024,219533186.649,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,rsa256(composite),256,ADD,200000,0.001102188,181457243.563,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,rsa256(composite),256,ADD,200000,0.000121856,1641281512.605,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,rsa256(composite),256,ADD,200000,0.000315812,633288140.889,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,rsa256(composite),256,ADD,200000,0.002394305,83531548.127,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,rsa256(composite),256,ADD,200000,0.000162691,1229324344.743,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,rsa256(composite),256,ADD,200000,0.002508109,79741350.064,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32,rsa256(composite),256,ADD,200000,0.000101439,1971629117.902,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32,rsa256(composite),256,ADD,200000,0.002358133,84812859.268,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,rsa256(composite),256,ADD,200000,0.000098661,2027143538.785,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,rsa256(composite),256,ADD,200000,0.002293544,87201293.939,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,rsa256(composite),256,ADD,200000,0.000098291,2034773716.635,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,rsa256(composite),256,ADD,200000,0.002166264,92324850.868,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,rsa256(composite),256,ADD,200000,0.000068470,2920988890.466,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,rsa256(composite),256,ADD,200000,0.002206483,90641986.168,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,rsa256(composite),256,ADD,200000,0.000067531,2961602170.709,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,rsa256(composite),256,ADD,200000,0.002073516,96454526.294,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACT,200000,0.003279189,60990690.130,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACT,200000,0.000790504,253003153.124,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACT,200000,0.001105923,180844421.286,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,rsa256(composite),256,SUBTRACT,200000,0.000121856,1641281512.605,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,rsa256(composite),256,SUBTRACT,200000,0.000314997,634926711.566,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,rsa256(composite),256,SUBTRACT,200000,0.002394990,83507655.355,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,rsa256(composite),256,SUBTRACT,200000,0.000161808,1236032950.386,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,rsa256(composite),256,SUBTRACT,200000,0.002334835,85659157.642,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32,rsa256(composite),256,SUBTRACT,200000,0.000101338,1973593982.675,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32,rsa256(composite),256,SUBTRACT,200000,0.002372065,84314720.311,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,rsa256(composite),256,SUBTRACT,200000,0.000098299,2034608646.797,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,rsa256(composite),256,SUBTRACT,200000,0.002300316,86944576.226,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,rsa256(composite),256,SUBTRACT,200000,0.000097980,2041232398.573,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,rsa256(composite),256,SUBTRACT,200000,0.002104873,95017608.584,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,rsa256(composite),256,SUBTRACT,200000,0.000068196,2932721722.023,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,rsa256(composite),256,SUBTRACT,200000,0.002225391,89871847.858,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,rsa256(composite),256,SUBTRACT,200000,0.000067250,2973978216.742,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,rsa256(composite),256,SUBTRACT,200000,0.002124578,94136341.716,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,rsa256(composite),256,ADDMOD,200000,0.008675139,23054385.668,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,rsa256(composite),256,ADDMOD,200000,0.002360039,84744362.880,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,rsa256(composite),256,ADDMOD,200000,0.006568315,30449209.431,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,rsa256(composite),256,ADDMOD,200000,0.000122880,1627604166.667,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,rsa256(composite),256,ADDMOD,200000,0.000446853,447574486.914,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,rsa256(composite),256,ADDMOD,200000,0.002525824,79182080.599,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,rsa256(composite),256,ADDMOD,200000,0.000225550,886721232.129,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,rsa256(composite),256,ADDMOD,200000,0.002441765,81907964.659,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32,rsa256(composite),256,ADDMOD,200000,0.000118187,1692233660.439,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32,rsa256(composite),256,ADDMOD,200000,0.002383575,83907575.171,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,rsa256(composite),256,ADDMOD,200000,0.000072600,2754819838.382,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,rsa256(composite),256,ADDMOD,200000,0.002247542,88986102.518,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,rsa256(composite),256,ADDMOD,200000,0.000070355,2842727812.365,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,rsa256(composite),256,ADDMOD,200000,0.002113325,94637597.264,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,rsa256(composite),256,ADDMOD,200000,0.000065191,3067907741.347,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,rsa256(composite),256,ADDMOD,200000,0.002208873,90543908.907,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,rsa256(composite),256,ADDMOD,200000,0.000064358,3107616827.566,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,rsa256(composite),256,ADDMOD,200000,0.002090066,95690759.333,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACTMOD,200000,0.008376205,23877161.581,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACTMOD,200000,0.002166869,92299072.892,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACTMOD,200000,0.007282790,27462002.995,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,rsa256(composite),256,SUBTRACTMOD,200000,0.000122496,1632706374.086,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,rsa256(composite),256,SUBTRACTMOD,200000,0.000491333,407055890.976,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,rsa256(composite),256,SUBTRACTMOD,200000,0.002592817,77136180.244,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,rsa256(composite),256,SUBTRACTMOD,200000,0.000249377,801998624.921,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,rsa256(composite),256,SUBTRACTMOD,200000,0.002341724,85407161.530,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32,rsa256(composite),256,SUBTRACTMOD,200000,0.000128421,1557377673.670,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32,rsa256(composite),256,SUBTRACTMOD,200000,0.002248305,88955901.842,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,rsa256(composite),256,SUBTRACTMOD,200000,0.000072555,2756530433.018,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,rsa256(composite),256,SUBTRACTMOD,200000,0.002255221,88683104.124,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,rsa256(composite),256,SUBTRACTMOD,200000,0.000072659,2752584241.896,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,rsa256(composite),256,SUBTRACTMOD,200000,0.002141921,93374124.317,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,rsa256(composite),256,SUBTRACTMOD,200000,0.000064826,3085182424.742,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,rsa256(composite),256,SUBTRACTMOD,200000,0.002211211,90448174.560,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,rsa256(composite),256,SUBTRACTMOD,200000,0.000064711,3090663298.454,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,rsa256(composite),256,SUBTRACTMOD,200000,0.002142390,93353684.287,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYOPERANDSCANNING,200000,0.004658045,42936467.849,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,200000,0.001226880,163015134.799,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,200000,0.002241423,89229029.780,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,200000,0.014076040,14208541.597,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,200000,0.016790050,11911816.845,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,rsa256(composite),256,MULTIPLYOPERANDSCANNING,200000,0.003695105,54125660.432,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,rsa256(composite),256,MULTIPLYOPERANDSCANNING,200000,0.006472638,30899302.604,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32,rsa256(composite),256,MULTIPLYOPERANDSCANNING,200000,0.001124517,177854136.074,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32,rsa256(composite),256,MULTIPLYOPERANDSCANNING,200000,0.003929705,50894405.433,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,200000,0.000197473,1012796479.356,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,200000,0.003021261,66197525.766,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,200000,0.000202819,986100876.018,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,200000,0.002906642,68807923.817,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,rsa256(composite),256,MULTIPLYOPERANDSCANNING,200000,0.000085770,2331818253.064,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,rsa256(composite),256,MULTIPLYOPERANDSCANNING,200000,0.002886619,69285208.477,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,200000,0.000084924,2355047205.348,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,200000,0.002682399,74560123.102,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,200000,0.004721635,42358208.353,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,200000,0.001208737,165461961.989,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,200000,0.002259442,88517430.384,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,200000,0.000121856,1641281512.605,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,200000,0.002263635,88353467.289,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,200000,0.004959711,40324930.303,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,200000,0.000781464,255929901.249,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,200000,0.003436020,58206879.001,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,200000,0.000428340,466918813.370,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,200000,0.003198426,62530757.629,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,200000,0.000235833,848057754.456,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,200000,0.003062678,65302327.405,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,200000,0.000225914,885292625.915,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,200000,0.002870333,69678325.751,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,200000,0.000148972,1342534415.132,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,200000,0.002947342,67857751.154,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,200000,0.000139049,1438341656.801,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,200000,0.002826112,70768602.862,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,rsa256(composite),256,MONTGOMERYMULTIPLICATION,200000,0.029975037,6672218.626,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,200000,0.008416236,23763592.099,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,200000,0.001954520,102326916.821,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,rsa256(composite),256,MONTGOMERYMULTIPLICATION,200000,0.000121696,1643439389.955,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,200000,0.001017337,196591696.352,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,200000,0.003078748,64961471.217,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,rsa256(composite),256,MONTGOMERYMULTIPLICATION,200000,0.000236801,844591069.979,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,rsa256(composite),256,MONTGOMERYMULTIPLICATION,200000,0.002316865,86323544.132,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32,rsa256(composite),256,MONTGOMERYMULTIPLICATION,200000,0.000071612,2792827545.487,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32,rsa256(composite),256,MONTGOMERYMULTIPLICATION,200000,0.002207292,90608763.590,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,200000,0.000108300,1846721902.203,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,200000,0.002251812,88817362.655,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,200000,0.000107106,1867309090.810,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,200000,0.002105996,94966941.435,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,rsa256(composite),256,MONTGOMERYMULTIPLICATION,200000,0.000106919,1870575362.440,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,rsa256(composite),256,MONTGOMERYMULTIPLICATION,200000,0.002226064,89844677.264,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,200000,0.000104326,1917067688.069,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,200000,0.002124063,94159164.367,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,rsa256(composite),256,COMPARE,200000,0.002075850,96346073.604,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,rsa256(composite),256,COMPARE,200000,0.000723168,276560901.718,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,rsa256(composite),256,COMPARE,200000,0.001072645,186454986.016,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,rsa256(composite),256,COMPARE,200000,0.000121856,1641281512.605,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,rsa256(composite),256,COMPARE,200000,0.000298398,670245727.453,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,rsa256(composite),256,COMPARE,200000,0.002334222,85681652.404,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,rsa256(composite),256,COMPARE,200000,0.000173171,1154927837.233,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,rsa256(composite),256,COMPARE,200000,0.002236809,89413088.042,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,rsa256(composite),256,COMPARE,200000,0.000087433,2287465073.152,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,rsa256(composite),256,COMPARE,200000,0.002259101,88530793.164,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,rsa256(composite),256,COMPARE,200000,0.000087843,2276789647.022,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,rsa256(composite),256,COMPARE,200000,0.002099787,95247756.014,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,rsa256(composite),256,COMPARE,200000,0.000031764,6296428887.561,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,rsa256(composite),256,COMPARE,200000,0.002168652,92223186.664,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,rsa256(composite),256,COMPARE,200000,0.000031648,6319520472.018,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,rsa256(composite),256,COMPARE,200000,0.002064042,96897254.563,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,rsa256(composite),256,REDUCE,25000,0.000769849,32473901.000,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,rsa256(composite),256,REDUCE,25000,0.000198377,126022670.639,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,rsa256(composite),256,REDUCE,25000,0.000737984,33876072.071,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,rsa256(composite),256,REDUCE,200000,0.000153472,1303169307.756,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,rsa256(composite),256,REDUCE,25000,0.000281670,88756359.891,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,rsa256(composite),256,REDUCE,25000,0.000669373,37348382.855,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,rsa256(composite),256,REDUCE,25000,0.000193684,129076218.840,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,rsa256(composite),256,REDUCE,25000,0.000626332,39914932.648,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,rsa256(composite),256,REDUCE,25000,0.000067444,370677926.846,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,rsa256(composite),256,REDUCE,25000,0.000476092,52510861.040,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,rsa256(composite),256,REDUCE,25000,0.000065149,383735503.354,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,rsa256(composite),256,REDUCE,25000,0.000494055,50601655.270,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,rsa256(composite),256,REDUCE,25000,0.000067498,370381284.004,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,rsa256(composite),256,REDUCE,25000,0.000474963,52635680.327,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,rsa256(composite),256,REDUCE,25000,0.000063759,392101607.218,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,rsa256(composite),256,REDUCE,25000,0.000491398,50875260.091,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,rsa256(composite),256,MODMUL,12500,0.001044354,11969121.649,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,rsa256(composite),256,MODMUL,12500,0.000315061,39674854.892,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,rsa256(composite),256,MODMUL,12500,0.000971640,12864847.135,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,rsa256(composite),256,MODMUL,200000,0.000567296,352549638.989,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,rsa256(composite),256,MODMUL,20000,0.000562962,35526374.779,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,rsa256(composite),256,MODMUL,20000,0.000896252,22315152.603,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,rsa256(composite),256,MODMUL,20000,0.000356641,56078808.574,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,rsa256(composite),256,MODMUL,20000,0.000721476,27720950.500,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,rsa256(composite),256,MODMUL,20000,0.000157345,127109195.873,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,rsa256(composite),256,MODMUL,20000,0.000434827,45995303.742,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,rsa256(composite),256,MODMUL,20000,0.000129222,154772411.058,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,rsa256(composite),256,MODMUL,20000,0.000493586,40519787.413,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,rsa256(composite),256,MODMUL,20000,0.000156731,127607177.728,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,rsa256(composite),256,MODMUL,20000,0.000502972,39763641.044,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,rsa256(composite),256,MODMUL,20000,0.000128416,155743838.272,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,rsa256(composite),256,MODMUL,20000,0.000494022,40484028.620,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,rsa256(composite),256,MODEXP,3125,0.027502944,113624.200,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,rsa256(composite),256,MODEXP,3125,0.006886452,453789.559,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,rsa256(composite),256,MODEXP,3125,0.012555383,248897.225,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,rsa256(composite),256,MODEXP,200000,0.189743370,1054055.275,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,rsa256(composite),256,MODEXP,20000,0.025289075,790855.340,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,rsa256(composite),256,MODEXP,20000,0.025650606,779708.674,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,rsa256(composite),256,MODEXP,20000,0.005246038,3812400.883,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,rsa256(composite),256,MODEXP,20000,0.005627467,3553996.848,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,rsa256(composite),256,MODEXP,20000,0.003379910,5917317.290,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,rsa256(composite),256,MODEXP,20000,0.003710195,5390552.227,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,rsa256(composite),256,MODEXP,20000,0.001952976,10240781.310,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,rsa256(composite),256,MODEXP,20000,0.002329108,8586978.310,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,rsa256(composite),256,MODEXP,20000,0.003381296,5914891.800,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,rsa256(composite),256,MODEXP,20000,0.003748016,5336156.481,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,rsa256(composite),256,MODEXP,20000,0.001894200,10558547.212,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,rsa256(composite),256,MODEXP,20000,0.002267429,8820562.759,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,rsa256(composite),256,EXPONENTIATION,3125,0.009963038,313659.349,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,rsa256(composite),256,EXPONENTIATION,3125,0.002147799,1454977.858,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,rsa256(composite),256,EXPONENTIATION,3125,0.039670234,78774.428,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,rsa256(composite),256,EXPONENTIATION,20000,0.038924184,513819.378,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,rsa256(composite),256,EXPONENTIATION,20000,0.039352408,508228.111,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,rsa256(composite),256,EXPONENTIATION,20000,0.013818404,1447345.150,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,rsa256(composite),256,EXPONENTIATION,20000,0.014151520,1413275.749,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,rsa256(composite),256,EXPONENTIATION,20000,0.000679204,29446233.791,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,rsa256(composite),256,EXPONENTIATION,20000,0.001037829,19270996.849,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,rsa256(composite),256,EXPONENTIATION,20000,0.000562628,35547464.470,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,rsa256(composite),256,EXPONENTIATION,20000,0.000928173,21547706.784,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,rsa256(composite),256,EXPONENTIATION,20000,0.000677177,29534376.495,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,rsa256(composite),256,EXPONENTIATION,20000,0.001021059,19587506.751,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,rsa256(composite),256,EXPONENTIATION,20000,0.000537123,37235418.574,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,rsa256(composite),256,EXPONENTIATION,20000,0.000902732,22154969.107,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,rsa256(composite),256,DIVIDE,25000,0.001106465,22594478.762,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,rsa256(composite),256,DIVIDE,25000,0.000275606,90709201.320,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,rsa256(composite),256,DIVIDE,25000,0.000727904,34345186.999,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,rsa256(composite),256,DIVIDE,200000,0.000204384,978550180.053,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,rsa256(composite),256,DIVIDE,25000,0.000534554,46767958.772,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,rsa256(composite),256,DIVIDE,25000,0.001031407,24238734.538,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,rsa256(composite),256,DIVIDE,25000,0.000339535,73630110.299,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,rsa256(composite),256,DIVIDE,25000,0.000830839,30090064.178,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,rsa256(composite),256,DIVIDE,25000,0.000130271,191907616.266,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,rsa256(composite),256,DIVIDE,25000,0.000631318,39599696.993,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,rsa256(composite),256,DIVIDE,25000,0.000133129,187787768.736,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,rsa256(composite),256,DIVIDE,25000,0.000655210,38155706.103,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,rsa256(composite),256,DIVIDE,25000,0.000120286,207837960.453,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,rsa256(composite),256,DIVIDE,25000,0.000621387,40232577.882,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,rsa256(composite),256,DIVIDE,25000,0.000117398,212950741.119,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,rsa256(composite),256,DIVIDE,25000,0.000639380,39100377.633,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,rsa256(composite),256,ISQRT,6250,0.000515783,12117499.665,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,rsa256(composite),256,ISQRT,6250,0.000114912,54389421.804,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,rsa256(composite),256,ISQRT,20000,0.005217592,3833185.884,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,rsa256(composite),256,ISQRT,20000,0.005568183,3591835.955,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,rsa256(composite),256,ISQRT,20000,0.003991695,5010402.843,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,rsa256(composite),256,ISQRT,20000,0.004369098,4577603.895,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,rsa256(composite),256,ISQRT,20000,0.001001362,19972797.243,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,rsa256(composite),256,ISQRT,20000,0.001361443,14690295.292,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,rsa256(composite),256,ISQRT,20000,0.000863751,23154821.532,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,rsa256(composite),256,ISQRT,20000,0.001208258,16552756.420,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,rsa256(composite),256,ISQRT,20000,0.000874912,22859442.940,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,rsa256(composite),256,ISQRT,20000,0.001238284,16151383.790,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,rsa256(composite),256,ISQRT,20000,0.000793965,25190026.967,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,rsa256(composite),256,ISQRT,20000,0.001172690,17054805.443,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,rsa256(composite),256,MODMUL_R2,200000,0.017479251,11442137.904,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,rsa256(composite),256,MODMUL_R2,200000,0.005052850,39581622.256,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,rsa256(composite),256,MODMUL_R2,200000,0.015406299,12981703.154,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,rsa256(composite),256,MODMUL_R2,200000,0.000168960,1183712121.212,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,rsa256(composite),256,MODMUL_R2,200000,0.001172819,170529297.330,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,rsa256(composite),256,MODMUL_R2,200000,0.003233313,61856058.394,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,rsa256(composite),256,MODMUL_R2,200000,0.000233902,855059118.431,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,rsa256(composite),256,MODMUL_R2,200000,0.002243236,89156914.453,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,rsa256(composite),256,MODMUL_R2,200000,0.000157416,1270519014.963,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,rsa256(composite),256,MODMUL_R2,200000,0.002350054,85104425.902,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,rsa256(composite),256,MODMUL_R2,200000,0.000162866,1228003202.554,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,rsa256(composite),256,MODMUL_R2,200000,0.002198158,90985270.119,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,rsa256(composite),256,MODMUL_R2,200000,0.000155072,1289723325.515,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,rsa256(composite),256,MODMUL_R2,200000,0.002298782,87002595.381,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,rsa256(composite),256,MODMUL_R2,200000,0.000153192,1305551461.689,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,rsa256(composite),256,MODMUL_R2,200000,0.002224058,89925713.540,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,brainpoolP512r1,512,ADD,100000,0.002335167,42823489.267,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,brainpoolP512r1,512,ADD,100000,0.000540535,185001888.089,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,brainpoolP512r1,512,ADD,100000,0.000721992,138505691.063,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,brainpoolP512r1,512,ADD,200000,0.000184320,1085069444.444,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,brainpoolP512r1,512,ADD,100000,0.000384010,260409868.002,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,brainpoolP512r1,512,ADD,100000,0.002451652,40788823.834,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,brainpoolP512r1,512,ADD,100000,0.000201577,496088413.677,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,brainpoolP512r1,512,ADD,100000,0.002289907,43669896.545,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32,brainpoolP512r1,512,ADD,100000,0.000388306,257528867.047,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32,brainpoolP512r1,512,ADD,100000,0.002588461,38632994.983,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,brainpoolP512r1,512,ADD,100000,0.000120883,827246230.608,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,brainpoolP512r1,512,ADD,100000,0.002288528,43696209.592,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,brainpoolP512r1,512,ADD,100000,0.000122923,813517428.205,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,brainpoolP512r1,512,ADD,100000,0.002173365,46011598.771,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,brainpoolP512r1,512,ADD,100000,0.000075952,1316621451.897,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,brainpoolP512r1,512,ADD,100000,0.002232199,44798873.000,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,brainpoolP512r1,512,ADD,100000,0.000104659,955483962.652,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,brainpoolP512r1,512,ADD,100000,0.002131701,46910894.573,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACT,100000,0.001930134,51809873.758,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACT,100000,0.000506189,197554663.826,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACT,100000,0.000725057,137920195.652,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,brainpoolP512r1,512,SUBTRACT,200000,0.000184864,1081876406.439,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,brainpoolP512r1,512,SUBTRACT,100000,0.000384070,260369217.638,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,brainpoolP512r1,512,SUBTRACT,100000,0.002479691,40327605.332,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,brainpoolP512r1,512,SUBTRACT,100000,0.000202626,493520105.898,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,brainpoolP512r1,512,SUBTRACT,100000,0.002291120,43646775.249,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32,brainpoolP512r1,512,SUBTRACT,100000,0.000386082,259012344.098,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32,brainpoolP512r1,512,SUBTRACT,100000,0.002578031,38789293.360,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,brainpoolP512r1,512,SUBTRACT,100000,0.000120150,832293081.015,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,brainpoolP512r1,512,SUBTRACT,100000,0.002328401,42947929.274,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,brainpoolP512r1,512,SUBTRACT,100000,0.000121355,824028686.010,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,brainpoolP512r1,512,SUBTRACT,100000,0.002190030,45661474.847,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,brainpoolP512r1,512,SUBTRACT,100000,0.000076249,1311492454.968,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,brainpoolP512r1,512,SUBTRACT,100000,0.002257097,44304697.806,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,brainpoolP512r1,512,SUBTRACT,100000,0.000105128,951221519.696,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,brainpoolP512r1,512,SUBTRACT,100000,0.002159675,46303264.445,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,brainpoolP512r1,512,ADDMOD,100000,0.005041256,19836326.412,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,brainpoolP512r1,512,ADDMOD,100000,0.001307590,76476571.900,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,brainpoolP512r1,512,ADDMOD,100000,0.003800196,26314432.265,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,brainpoolP512r1,512,ADDMOD,200000,0.000184928,1081501989.964,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,brainpoolP512r1,512,ADDMOD,100000,0.000528169,189333324.340,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,brainpoolP512r1,512,ADDMOD,100000,0.002577313,38800099.336,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,brainpoolP512r1,512,ADDMOD,100000,0.000275773,362617084.925,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,brainpoolP512r1,512,ADDMOD,100000,0.002328290,42949976.831,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32,brainpoolP512r1,512,ADDMOD,100000,0.000288929,346105820.223,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32,brainpoolP512r1,512,ADDMOD,100000,0.002510499,39832717.967,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,brainpoolP512r1,512,ADDMOD,100000,0.000129546,771926689.950,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,brainpoolP512r1,512,ADDMOD,100000,0.002309792,43293940.967,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,brainpoolP512r1,512,ADDMOD,100000,0.000130230,767872007.723,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,brainpoolP512r1,512,ADDMOD,100000,0.002227036,44902731.876,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,brainpoolP512r1,512,ADDMOD,100000,0.000086712,1153242777.175,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,brainpoolP512r1,512,ADDMOD,100000,0.002223904,44965969.516,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,brainpoolP512r1,512,ADDMOD,100000,0.000073331,1363679739.518,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,brainpoolP512r1,512,ADDMOD,100000,0.002074202,48211312.138,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACTMOD,100000,0.004522715,22110612.779,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACTMOD,100000,0.001180714,84694515.619,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACTMOD,100000,0.004189954,23866610.338,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,brainpoolP512r1,512,SUBTRACTMOD,200000,0.000185280,1079447322.971,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,brainpoolP512r1,512,SUBTRACTMOD,100000,0.000602412,165999338.935,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,brainpoolP512r1,512,SUBTRACTMOD,100000,0.002651577,37713406.015,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,brainpoolP512r1,512,SUBTRACTMOD,100000,0.000309360,323248015.980,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,brainpoolP512r1,512,SUBTRACTMOD,100000,0.002395836,41739083.031,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32,brainpoolP512r1,512,SUBTRACTMOD,100000,0.000311285,321248990.882,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32,brainpoolP512r1,512,SUBTRACTMOD,100000,0.002492681,40117447.370,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,brainpoolP512r1,512,SUBTRACTMOD,100000,0.000134550,743218262.716,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,brainpoolP512r1,512,SUBTRACTMOD,100000,0.002326844,42976667.679,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,brainpoolP512r1,512,SUBTRACTMOD,100000,0.000133493,749103039.690,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,brainpoolP512r1,512,SUBTRACTMOD,100000,0.002184637,45774194.828,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,brainpoolP512r1,512,SUBTRACTMOD,100000,0.000085950,1163466940.086,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,brainpoolP512r1,512,SUBTRACTMOD,100000,0.002256372,44318932.986,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,brainpoolP512r1,512,SUBTRACTMOD,100000,0.000073208,1365970783.587,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,brainpoolP512r1,512,SUBTRACTMOD,100000,0.002109760,47398756.008,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,100000,0.004887303,20461182.589,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,100000,0.001177518,84924392.661,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,100000,0.001277622,78270411.131,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,100000,0.037312097,2680095.950,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,100000,0.040147491,2490815.679,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,100000,0.009771414,10233933.406,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,100000,0.012543657,7972156.740,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,100000,0.009861952,10139980.417,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,100000,0.012877755,7765328.672,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,100000,0.000454983,219788415.374,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,100000,0.003232353,30937214.919,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,100000,0.000455830,219380025.590,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,100000,0.003172563,31520256.834,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,100000,0.000253778,394045265.182,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,100000,0.003054704,32736396.136,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,100000,0.000234110,427149649.066,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,100000,0.002937085,34047363.343,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,100000,0.004891634,20443066.482,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,100000,0.001179566,84776941.602,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,100000,0.001282328,77983165.603,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,200000,0.000184320,1085069444.444,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,100000,0.008311707,12031222.945,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,100000,0.011098910,9009893.759,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,100000,0.002266748,44116064.199,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,100000,0.004982139,20071699.958,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,100000,0.001302335,76785158.737,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,100000,0.004104510,24363444.413,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,100000,0.000679777,147107065.166,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,100000,0.003453654,28954840.258,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,100000,0.000679511,147164653.351,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,100000,0.003338528,29953320.905,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,100000,0.000315390,317067753.982,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,100000,0.003103685,32219764.140,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,100000,0.000452506,220991541.594,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,100000,0.003145086,31795632.964,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,100000,0.033158376,3015829.243,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,100000,0.008681287,11519029.380,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,100000,0.001935011,51679292.545,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,200000,0.000181504,1101904090.268,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,100000,0.002006382,49840957.928,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,100000,0.004099673,24392189.272,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,100000,0.000578788,172774835.103,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,100000,0.002625173,38092728.216,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,100000,0.000102886,971949547.853,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,100000,0.002228840,44866388.637,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,100000,0.000177278,564085758.237,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,100000,0.002310501,43280656.370,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,100000,0.000137579,726855184.383,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,100000,0.002145256,46614483.871,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,100000,0.000162894,613896287.329,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,100000,0.002401877,41634104.980,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,100000,0.000149155,670443599.927,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,100000,0.002219491,45055375.350,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,brainpoolP512r1,512,COMPARE,100000,0.001276740,78324483.317,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,brainpoolP512r1,512,COMPARE,100000,0.000487628,205074366.094,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,brainpoolP512r1,512,COMPARE,100000,0.000670680,149102412.127,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,brainpoolP512r1,512,COMPARE,200000,0.000184320,1085069444.444,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,brainpoolP512r1,512,COMPARE,100000,0.000379143,263752732.295,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,brainpoolP512r1,512,COMPARE,100000,0.002458511,40675026.823,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,brainpoolP512r1,512,COMPARE,100000,0.000199071,502333299.922,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,brainpoolP512r1,512,COMPARE,100000,0.002283633,43789873.645,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,brainpoolP512r1,512,COMPARE,100000,0.000107877,926981528.970,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,brainpoolP512r1,512,COMPARE,100000,0.002286450,43735922.384,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,brainpoolP512r1,512,COMPARE,100000,0.000107901,926775502.515,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,brainpoolP512r1,512,COMPARE,100000,0.002160082,46294537.804,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,brainpoolP512r1,512,COMPARE,100000,0.000029075,3439378821.834,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,brainpoolP512r1,512,COMPARE,100000,0.002123973,47081578.154,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,brainpoolP512r1,512,COMPARE,100000,0.000030272,3303382663.454,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,brainpoolP512r1,512,COMPARE,100000,0.002051916,48734937.716,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,brainpoolP512r1,512,REDUCE,12500,0.000403158,31005212.496,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,brainpoolP512r1,512,REDUCE,12500,0.000106315,117575117.138,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,brainpoolP512r1,512,REDUCE,12500,0.000434425,28773668.925,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,brainpoolP512r1,512,REDUCE,200000,0.000281568,710307989.544,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,brainpoolP512r1,512,REDUCE,20000,0.000614839,32528839.619,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,brainpoolP512r1,512,REDUCE,20000,0.001156994,17286174.838,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,brainpoolP512r1,512,REDUCE,20000,0.000492521,40607407.254,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,brainpoolP512r1,512,REDUCE,20000,0.001067964,18727222.824,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,brainpoolP512r1,512,REDUCE,20000,0.000142585,140267246.116,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,brainpoolP512r1,512,REDUCE,20000,0.000706046,28326766.053,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,brainpoolP512r1,512,REDUCE,20000,0.000136204,146838622.305,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,brainpoolP512r1,512,REDUCE,20000,0.000695632,28750833.197,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,brainpoolP512r1,512,REDUCE,20000,0.000132539,150899025.859,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,brainpoolP512r1,512,REDUCE,20000,0.000711110,28125043.428,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,brainpoolP512r1,512,REDUCE,20000,0.000156781,127566481.422,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,brainpoolP512r1,512,REDUCE,20000,0.000728958,27436423.126,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL,6250,0.001020818,6122541.187,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL,6250,0.000284791,21945918.794,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL,6250,0.000795581,7855894.285,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,brainpoolP512r1,512,MODMUL,200000,0.001859584,107550936.123,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,brainpoolP512r1,512,MODMUL,20000,0.001761260,11355506.628,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,brainpoolP512r1,512,MODMUL,20000,0.002339405,8549182.350,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,brainpoolP512r1,512,MODMUL,20000,0.001256691,15914811.187,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,brainpoolP512r1,512,MODMUL,20000,0.001867737,10708145.850,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,brainpoolP512r1,512,MODMUL,20000,0.000529524,37769766.641,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,brainpoolP512r1,512,MODMUL,20000,0.001089765,18352580.741,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,brainpoolP512r1,512,MODMUL,20000,0.000427371,46797746.347,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,brainpoolP512r1,512,MODMUL,20000,0.001025693,19499012.511,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,brainpoolP512r1,512,MODMUL,20000,0.000532563,37554245.168,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,brainpoolP512r1,512,MODMUL,20000,0.001108793,18037632.567,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,brainpoolP512r1,512,MODMUL,20000,0.000513361,38958938.449,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,brainpoolP512r1,512,MODMUL,20000,0.001094433,18274302.544,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,brainpoolP512r1,512,MODEXP,1562,0.076734161,20355.992,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,brainpoolP512r1,512,MODEXP,1562,0.020775952,75183.077,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,brainpoolP512r1,512,MODEXP,1562,0.019055833,81969.652,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,brainpoolP512r1,512,MODEXP,200000,0.469232857,426227.612,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,brainpoolP512r1,512,MODEXP,20000,0.368964043,54205.824,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,brainpoolP512r1,512,MODEXP,20000,0.369620588,54109.540,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,brainpoolP512r1,512,MODEXP,20000,0.033974047,588684.651,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,brainpoolP512r1,512,MODEXP,20000,0.034591365,578178.977,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,brainpoolP512r1,512,MODEXP,20000,0.025817723,774661.654,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,brainpoolP512r1,512,MODEXP,20000,0.026427074,756799.635,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,brainpoolP512r1,512,MODEXP,20000,0.013211329,1513852.240,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,brainpoolP512r1,512,MODEXP,20000,0.013817884,1447399.619,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,brainpoolP512r1,512,MODEXP,20000,0.025222615,792939.194,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,brainpoolP512r1,512,MODEXP,20000,0.025939530,771023.993,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,brainpoolP512r1,512,MODEXP,20000,0.014081500,1420303.240,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,brainpoolP512r1,512,MODEXP,20000,0.014715318,1359127.952,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,brainpoolP512r1,512,EXPONENTIATION,1562,0.016803770,92955.331,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,brainpoolP512r1,512,EXPONENTIATION,1562,0.004230350,369236.586,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,brainpoolP512r1,512,EXPONENTIATION,1562,0.054472463,28675.039,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,brainpoolP512r1,512,EXPONENTIATION,20000,0.311475591,64210.489,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,brainpoolP512r1,512,EXPONENTIATION,20000,0.312455144,64009.188,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,brainpoolP512r1,512,EXPONENTIATION,20000,0.082425036,242644.723,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,brainpoolP512r1,512,EXPONENTIATION,20000,0.083061201,240786.309,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,brainpoolP512r1,512,EXPONENTIATION,20000,0.021981676,909848.731,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,brainpoolP512r1,512,EXPONENTIATION,20000,0.022505609,888667.353,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,brainpoolP512r1,512,EXPONENTIATION,20000,0.021620171,925062.065,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,brainpoolP512r1,512,EXPONENTIATION,20000,0.022291349,897209.049,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,brainpoolP512r1,512,EXPONENTIATION,20000,0.022542369,887218.199,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,brainpoolP512r1,512,EXPONENTIATION,20000,0.023160702,863531.683,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,brainpoolP512r1,512,EXPONENTIATION,20000,0.022153008,902811.935,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,brainpoolP512r1,512,EXPONENTIATION,20000,0.022996371,869702.442,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,brainpoolP512r1,512,DIVIDE,12500,0.000588230,21250192.198,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,brainpoolP512r1,512,DIVIDE,12500,0.000150060,83300018.813,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,brainpoolP512r1,512,DIVIDE,12500,0.000430032,29067602.323,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,brainpoolP512r1,512,DIVIDE,200000,0.000385696,518543101.303,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,brainpoolP512r1,512,DIVIDE,20000,0.001207552,16562434.064,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,brainpoolP512r1,512,DIVIDE,20000,0.001922871,10401113.826,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,brainpoolP512r1,512,DIVIDE,20000,0.001091256,18327505.506,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,brainpoolP512r1,512,DIVIDE,20000,0.001825572,10955470.446,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,brainpoolP512r1,512,DIVIDE,20000,0.000367148,54473941.754,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,brainpoolP512r1,512,DIVIDE,20000,0.001076691,18575432.111,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,brainpoolP512r1,512,DIVIDE,20000,0.000355457,56265593.914,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,brainpoolP512r1,512,DIVIDE,20000,0.001066900,18745899.777,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,brainpoolP512r1,512,DIVIDE,20000,0.000352044,56811088.577,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,brainpoolP512r1,512,DIVIDE,20000,0.001086280,18411458.599,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,brainpoolP512r1,512,DIVIDE,20000,0.000349094,57291152.265,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,brainpoolP512r1,512,DIVIDE,20000,0.001037849,19270626.131,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,brainpoolP512r1,512,ISQRT,3125,0.000528230,5915983.763,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,brainpoolP512r1,512,ISQRT,3125,0.000101258,30861762.510,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,brainpoolP512r1,512,ISQRT,20000,0.025532333,783320.507,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,brainpoolP512r1,512,ISQRT,20000,0.026104398,766154.423,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,brainpoolP512r1,512,ISQRT,20000,0.025361038,788611.253,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,brainpoolP512r1,512,ISQRT,20000,0.025941867,770954.534,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,brainpoolP512r1,512,ISQRT,20000,0.005195779,3849278.447,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,brainpoolP512r1,512,ISQRT,20000,0.005783269,3458251.721,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,brainpoolP512r1,512,ISQRT,20000,0.004897434,4083771.238,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,brainpoolP512r1,512,ISQRT,20000,0.005522713,3621408.570,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,brainpoolP512r1,512,ISQRT,20000,0.004937711,4050459.821,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,brainpoolP512r1,512,ISQRT,20000,0.005552814,3601777.372,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,brainpoolP512r1,512,ISQRT,20000,0.005977315,3345983.934,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,brainpoolP512r1,512,ISQRT,20000,0.006611423,3025067.379,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL_R2,100000,0.016813146,5947726.853,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL_R2,100000,0.004566001,21901002.615,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL_R2,100000,0.012608471,7931175.772,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,brainpoolP512r1,512,MODMUL_R2,200000,0.000271360,737028301.887,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,brainpoolP512r1,512,MODMUL_R2,100000,0.002980632,33549931.591,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,brainpoolP512r1,512,MODMUL_R2,100000,0.005081581,19678914.850,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,brainpoolP512r1,512,MODMUL_R2,100000,0.000518974,192687889.609,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,brainpoolP512r1,512,MODMUL_R2,100000,0.002603898,38403962.067,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,brainpoolP512r1,512,MODMUL_R2,100000,0.000299616,333760526.436,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,brainpoolP512r1,512,MODMUL_R2,100000,0.002499589,40006577.113,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,brainpoolP512r1,512,MODMUL_R2,100000,0.000190960,523670009.300,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,brainpoolP512r1,512,MODMUL_R2,100000,0.002244291,44557501.437,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,brainpoolP512r1,512,MODMUL_R2,100000,0.000283907,352228027.026,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,brainpoolP512r1,512,MODMUL_R2,100000,0.002362817,42322363.482,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,brainpoolP512r1,512,MODMUL_R2,100000,0.000193030,518054211.680,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,brainpoolP512r1,512,MODMUL_R2,100000,0.002243342,44576349.773,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,p1024,1024,ADD,50000,0.001614828,30963049.706,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,p1024,1024,ADD,50000,0.000401044,124674589.471,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,p1024,1024,ADD,50000,0.000523141,95576526.582,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,p1024,1024,ADD,200000,0.000338432,590960665.658,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,p1024,1024,ADD,50000,0.000614748,81334141.871,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,p1024,1024,ADD,50000,0.002813388,17772166.440,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,p1024,1024,ADD,50000,0.000317843,157310374.854,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,p1024,1024,ADD,50000,0.002408488,20759912.328,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32,p1024,1024,ADD,50000,0.000577725,86546366.556,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32,p1024,1024,ADD,50000,0.002788965,17927797.535,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,p1024,1024,ADD,50000,0.000158175,316105556.950,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,p1024,1024,ADD,50000,0.002333597,21426150.374,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,p1024,1024,ADD,50000,0.000157103,318262539.769,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,p1024,1024,ADD,50000,0.002207083,22654335.595,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,p1024,1024,ADD,50000,0.000105000,476190545.966,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,p1024,1024,ADD,50000,0.002149995,23255867.940,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,p1024,1024,ADD,50000,0.000105849,472371012.146,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,p1024,1024,ADD,50000,0.002166699,23076578.940,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,p1024,1024,SUBTRACT,50000,0.001451028,34458329.161,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,p1024,1024,SUBTRACT,50000,0.000393285,127134263.715,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,p1024,1024,SUBTRACT,50000,0.000514758,97133022.845,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,p1024,1024,SUBTRACT,200000,0.000337920,591856060.606,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,p1024,1024,SUBTRACT,50000,0.000612824,81589495.068,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,p1024,1024,SUBTRACT,50000,0.002700608,18514349.424,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,p1024,1024,SUBTRACT,50000,0.000319375,156555765.023,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,p1024,1024,SUBTRACT,50000,0.002396698,20862036.075,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32,p1024,1024,SUBTRACT,50000,0.000580625,86114105.889,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32,p1024,1024,SUBTRACT,50000,0.002729203,18320366.600,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,p1024,1024,SUBTRACT,50000,0.000156626,319231823.233,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,p1024,1024,SUBTRACT,50000,0.002360571,21181315.866,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,p1024,1024,SUBTRACT,50000,0.000158711,315038050.285,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,p1024,1024,SUBTRACT,50000,0.002205348,22672158.857,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,p1024,1024,SUBTRACT,50000,0.000105038,476018229.143,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,p1024,1024,SUBTRACT,50000,0.002173415,23005270.127,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,p1024,1024,SUBTRACT,50000,0.000106150,471031588.529,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,p1024,1024,SUBTRACT,50000,0.002181154,22923644.633,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,p1024,1024,ADDMOD,50000,0.003762366,13289509.990,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,p1024,1024,ADDMOD,50000,0.000884469,56531092.746,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,p1024,1024,ADDMOD,50000,0.002729500,18318373.304,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,p1024,1024,ADDMOD,200000,0.000338432,590960665.658,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,p1024,1024,ADDMOD,50000,0.000899794,55568275.280,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,p1024,1024,ADDMOD,50000,0.002929578,17067304.508,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,p1024,1024,ADDMOD,50000,0.000464914,107546781.541,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,p1024,1024,ADDMOD,50000,0.002554072,19576581.954,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32,p1024,1024,ADDMOD,50000,0.000512529,97555457.315,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32,p1024,1024,ADDMOD,50000,0.002703607,18493811.990,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,p1024,1024,ADDMOD,50000,0.000220481,226776954.640,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,p1024,1024,ADDMOD,50000,0.002400443,20829488.417,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,p1024,1024,ADDMOD,50000,0.000223962,223252150.852,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,p1024,1024,ADDMOD,50000,0.002263008,22094486.859,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,p1024,1024,ADDMOD,50000,0.000087975,568343276.130,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,p1024,1024,ADDMOD,50000,0.002310908,21636516.976,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,p1024,1024,ADDMOD,50000,0.000088298,566264274.847,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,p1024,1024,ADDMOD,50000,0.002135698,23411549.286,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,p1024,1024,SUBTRACTMOD,50000,0.002800305,17855197.758,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,p1024,1024,SUBTRACTMOD,50000,0.000679408,73593484.931,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,p1024,1024,SUBTRACTMOD,50000,0.002643290,18915820.987,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,p1024,1024,SUBTRACTMOD,200000,0.000337760,592136428.233,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,p1024,1024,SUBTRACTMOD,50000,0.000983870,50819721.799,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,p1024,1024,SUBTRACTMOD,50000,0.003087820,16192653.919,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,p1024,1024,SUBTRACTMOD,50000,0.000461588,108321713.045,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,p1024,1024,SUBTRACTMOD,50000,0.002645901,18897154.527,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32,p1024,1024,SUBTRACTMOD,50000,0.000512571,97547469.777,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32,p1024,1024,SUBTRACTMOD,50000,0.002677734,18672504.369,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,p1024,1024,SUBTRACTMOD,50000,0.000223967,223247160.971,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,p1024,1024,SUBTRACTMOD,50000,0.002414610,20707277.745,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,p1024,1024,SUBTRACTMOD,50000,0.000222007,225218136.622,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,p1024,1024,SUBTRACTMOD,50000,0.002304688,21694910.739,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,p1024,1024,SUBTRACTMOD,50000,0.000087418,571964509.432,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,p1024,1024,SUBTRACTMOD,50000,0.002139343,23371661.223,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,p1024,1024,SUBTRACTMOD,50000,0.000087679,570261883.777,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,p1024,1024,SUBTRACTMOD,50000,0.002098534,23826157.131,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.007691661,6500546.540,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.001828313,27347615.133,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.002357328,21210455.575,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.091517685,546342.491,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.094288249,530288.775,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.024419400,2047552.354,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.027426083,1823082.061,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.024970916,2002329.431,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.028081427,1780536.297,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.000567684,88077177.440,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.003423489,14604983.439,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.000588976,84893099.522,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.003249691,15386078.044,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.000385218,129796631.638,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.003127580,15986801.375,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.000368556,135664590.862,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,p1024,1024,MULTIPLYOPERANDSCANNING,50000,0.003073982,16265547.619,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.007710783,6484425.783,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.001825454,27390446.363,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.002360462,21182294.170,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,p1024,1024,MULTIPLYPRODUCTSCANNING,200000,0.000352256,567768895.349,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.033522652,1491528.773,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.036258675,1378980.341,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.008645620,5783275.255,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.011593512,4312756.987,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.002260092,22122993.168,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.005066547,9868654.098,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.002190690,22823859.527,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.005018257,9963618.818,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.002196716,22761248.772,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.004928963,10144121.659,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000881577,56716541.456,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.003601350,13883682.707,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.000822932,60758362.027,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,p1024,1024,MULTIPLYPRODUCTSCANNING,50000,0.003507556,14254939.963,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.050572994,988669.961,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.011628689,4299710.823,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.002958522,16900330.617,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,p1024,1024,MONTGOMERYMULTIPLICATION,200000,0.000563136,355153994.772,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.005774641,8658546.909,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.007856210,6364391.994,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.001131042,44207024.627,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.003202882,15610940.332,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000258409,193491704.729,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.002380237,21006311.704,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000423138,118164767.400,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.002559431,19535592.431,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000338204,147839752.054,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.002419128,20668604.843,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000311702,160409639.516,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.002413421,20717479.495,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.000212272,235546793.331,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,p1024,1024,MONTGOMERYMULTIPLICATION,50000,0.002188181,22850029.518,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,p1024,1024,COMPARE,50000,0.001055319,47379037.291,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,p1024,1024,COMPARE,50000,0.000306140,163323943.092,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,p1024,1024,COMPARE,50000,0.000387897,128900211.599,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,p1024,1024,COMPARE,200000,0.000338944,590067975.831,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,p1024,1024,COMPARE,50000,0.000529759,94382543.294,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,p1024,1024,COMPARE,50000,0.002638273,18951791.575,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,p1024,1024,COMPARE,50000,0.000274910,181877699.168,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,p1024,1024,COMPARE,50000,0.002348301,21291989.392,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,p1024,1024,COMPARE,50000,0.000156063,320383425.869,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,p1024,1024,COMPARE,50000,0.002364326,21147675.594,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,p1024,1024,COMPARE,50000,0.000157719,317019484.557,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,p1024,1024,COMPARE,50000,0.002466993,20267588.932,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,p1024,1024,COMPARE,50000,0.000032748,1526809050.012,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,p1024,1024,COMPARE,50000,0.002064741,24216112.526,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,p1024,1024,COMPARE,50000,0.000031822,1571242313.750,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,p1024,1024,COMPARE,50000,0.002056983,24307444.346,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,p1024,1024,REDUCE,6250,0.000143278,43621496.350,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,p1024,1024,REDUCE,6250,0.000038818,161007577.573,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,p1024,1024,REDUCE,6250,0.000285864,21863544.057,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,p1024,1024,REDUCE,200000,0.000454176,440357922.920,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,p1024,1024,REDUCE,20000,0.003315048,6033095.144,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,p1024,1024,REDUCE,20000,0.004351384,4596238.799,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,p1024,1024,REDUCE,20000,0.001810332,11047697.307,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,p1024,1024,REDUCE,20000,0.002834668,7055500.035,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,p1024,1024,REDUCE,20000,0.000523901,38175154.149,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,p1024,1024,REDUCE,20000,0.001467537,13628276.862,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,p1024,1024,REDUCE,20000,0.000535945,37317261.991,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,p1024,1024,REDUCE,20000,0.001614389,12388587.806,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,p1024,1024,REDUCE,20000,0.000447574,44685346.786,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,p1024,1024,REDUCE,20000,0.001461425,13685272.566,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,p1024,1024,REDUCE,20000,0.000481476,41538937.935,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,p1024,1024,REDUCE,20000,0.001468100,13623050.714,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,p1024,1024,MODMUL,3125,0.001390040,2248136.684,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,p1024,1024,MODMUL,3125,0.000338442,9233486.436,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,p1024,1024,MODMUL,3125,0.000952766,3279923.975,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,p1024,1024,MODMUL,200000,0.005202752,38441194.199,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,p1024,1024,MODMUL,20000,0.013434560,1488697.808,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,p1024,1024,MODMUL,20000,0.014460096,1383116.680,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,p1024,1024,MODMUL,20000,0.004812936,4155467.675,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,p1024,1024,MODMUL,20000,0.005842870,3422975.328,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,p1024,1024,MODMUL,20000,0.002329530,8585422.738,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,p1024,1024,MODMUL,20000,0.003339880,5988239.045,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,p1024,1024,MODMUL,20000,0.001740534,11490726.190,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,p1024,1024,MODMUL,20000,0.002807189,7124564.792,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,p1024,1024,MODMUL,20000,0.002262143,8841174.013,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,p1024,1024,MODMUL,20000,0.003254838,6144699.066,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,p1024,1024,MODMUL,20000,0.001689713,11836329.723,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,p1024,1024,MODMUL,20000,0.002696384,7417341.213,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,p1024,1024,MODEXP,781,0.250933913,3112.373,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,p1024,1024,MODEXP,781,0.067345260,11596.956,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,p1024,1024,MODEXP,781,0.046415128,16826.411,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,p1024,1024,MODEXP,200000,2.433285475,82193.397,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,p1024,1024,MODEXP,20000,3.004660793,6656.325,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,p1024,1024,MODEXP,20000,3.003860872,6658.098,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,p1024,1024,MODEXP,20000,0.453110787,44139.316,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,p1024,1024,MODEXP,20000,0.454167591,44036.608,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,p1024,1024,MODEXP,20000,0.198296667,100858.982,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,p1024,1024,MODEXP,20000,0.199936594,100031.713,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,p1024,1024,MODEXP,20000,0.122971631,162639.137,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,p1024,1024,MODEXP,20000,0.123921753,161392.165,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,p1024,1024,MODEXP,20000,0.201299561,99354.414,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,p1024,1024,MODEXP,20000,0.202273021,98876.261,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,p1024,1024,MODEXP,20000,0.126309147,158341.660,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,p1024,1024,MODEXP,20000,0.127309548,157097.408,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,p1024,1024,EXPONENTIATION,781,0.033045597,23634.011,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,p1024,1024,EXPONENTIATION,781,0.009383795,83228.587,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,p1024,1024,EXPONENTIATION,781,0.151203585,5165.221,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,p1024,1024,EXPONENTIATION,20000,2.409761306,8299.577,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,p1024,1024,EXPONENTIATION,20000,2.415450108,8280.030,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,p1024,1024,EXPONENTIATION,20000,0.635271281,31482.613,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,p1024,1024,EXPONENTIATION,20000,0.636619298,31415.950,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,p1024,1024,EXPONENTIATION,20000,0.174721370,114467.967,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,p1024,1024,EXPONENTIATION,20000,0.175586133,113904.211,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,p1024,1024,EXPONENTIATION,20000,0.167312639,119536.696,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,p1024,1024,EXPONENTIATION,20000,0.168365722,118789.025,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,p1024,1024,EXPONENTIATION,20000,0.176794227,113125.866,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,p1024,1024,EXPONENTIATION,20000,0.177764802,112508.212,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,p1024,1024,EXPONENTIATION,20000,0.169384578,118074.504,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,p1024,1024,EXPONENTIATION,20000,0.170319036,117426.686,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,p1024,1024,DIVIDE,6250,0.000312017,20030960.760,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,p1024,1024,DIVIDE,6250,0.000081855,76354539.779,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,p1024,1024,DIVIDE,6250,0.000281990,22163907.945,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,p1024,1024,DIVIDE,200000,0.000591872,337910899.654,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,p1024,1024,DIVIDE,20000,0.054339670,368055.235,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,p1024,1024,DIVIDE,20000,0.056040265,356886.249,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,p1024,1024,DIVIDE,20000,0.006906825,2895686.512,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,p1024,1024,DIVIDE,20000,0.008188486,2442453.941,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,p1024,1024,DIVIDE,20000,0.001555338,12858941.453,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,p1024,1024,DIVIDE,20000,0.002692482,7428090.600,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,p1024,1024,DIVIDE,20000,0.001525473,13110687.428,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,p1024,1024,DIVIDE,20000,0.002830494,7065904.345,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,p1024,1024,DIVIDE,20000,0.001541079,12977919.831,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,p1024,1024,DIVIDE,20000,0.002812015,7112337.571,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,p1024,1024,DIVIDE,20000,0.001533629,13040963.918,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,p1024,1024,DIVIDE,20000,0.002744605,7287023.066,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,p1024,1024,ISQRT,1562,0.000502624,3107690.683,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,p1024,1024,ISQRT,1562,0.000103787,15050055.139,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,p1024,1024,ISQRT,20000,0.664321015,30105.927,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,p1024,1024,ISQRT,20000,0.663178122,30157.810,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,p1024,1024,ISQRT,20000,0.179441243,111457.097,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,p1024,1024,ISQRT,20000,0.180491024,110808.834,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,p1024,1024,ISQRT,20000,0.034697487,576410.620,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,p1024,1024,ISQRT,20000,0.035711800,560038.979,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,p1024,1024,ISQRT,20000,0.032917082,607587.270,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,p1024,1024,ISQRT,20000,0.033962703,588881.279,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,p1024,1024,ISQRT,20000,0.035080114,570123.575,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,p1024,1024,ISQRT,20000,0.036081076,554307.194,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,p1024,1024,ISQRT,20000,0.032511318,615170.384,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,p1024,1024,ISQRT,20000,0.033481921,597337.291,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,p1024,1024,MODMUL_R2,50000,0.022707706,2201895.694,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,p1024,1024,MODMUL_R2,50000,0.005428430,9210766.308,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,p1024,1024,MODMUL_R2,50000,0.015079689,3315718.253,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,p1024,1024,MODMUL_R2,200000,0.001016480,196757437.431,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,p1024,1024,MODMUL_R2,50000,0.009352564,5346127.512,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,p1024,1024,MODMUL_R2,50000,0.011444578,4368881.052,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,p1024,1024,MODMUL_R2,50000,0.001277612,39135512.240,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,p1024,1024,MODMUL_R2,50000,0.003386443,14764754.575,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,p1024,1024,MODMUL_R2,50000,0.000696320,71806069.239,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,p1024,1024,MODMUL_R2,50000,0.002766307,18074639.005,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,p1024,1024,MODMUL_R2,50000,0.000567164,88157916.477,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,p1024,1024,MODMUL_R2,50000,0.002741345,18239221.790,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,p1024,1024,MODMUL_R2,50000,0.000606559,82432214.025,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,p1024,1024,MODMUL_R2,50000,0.002689657,18589731.123,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,p1024,1024,MODMUL_R2,50000,0.000477641,104681129.492,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,p1024,1024,MODMUL_R2,50000,0.002534379,19728699.135,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,p2048,2048,ADD,25000,0.001247258,20043967.774,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,p2048,2048,ADD,25000,0.000321188,77836045.341,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,p2048,2048,ADD,25000,0.000406409,61514380.954,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,p2048,2048,ADD,200000,0.000666528,300062412.982,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,p2048,2048,ADD,25000,0.000635433,39343253.918,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,p2048,2048,ADD,25000,0.002701414,9254412.715,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,p2048,2048,ADD,25000,0.000323308,77325648.168,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,p2048,2048,ADD,25000,0.002356328,10609728.303,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32,p2048,2048,ADD,25000,0.000163180,153205052.559,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32,p2048,2048,ADD,25000,0.002312545,10810600.601,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,p2048,2048,ADD,25000,0.000161021,155259292.077,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,p2048,2048,ADD,25000,0.002211733,11303353.408,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,p2048,2048,ADD,25000,0.000160094,156158255.726,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,p2048,2048,ADD,25000,0.002333406,10713952.007,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,p2048,2048,ADD,25000,0.000107096,233435456.697,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,p2048,2048,ADD,25000,0.002161532,11565870.872,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,p2048,2048,ADD,25000,0.000104416,239426803.742,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,p2048,2048,ADD,25000,0.002163307,11556381.123,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,p2048,2048,SUBTRACT,25000,0.001230317,20319966.428,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,p2048,2048,SUBTRACT,25000,0.000312292,80053294.018,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,p2048,2048,SUBTRACT,25000,0.000417018,59949448.007,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,p2048,2048,SUBTRACT,200000,0.000666624,300019201.229,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,p2048,2048,SUBTRACT,25000,0.000639492,39093530.149,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,p2048,2048,SUBTRACT,25000,0.002704920,9242417.670,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,p2048,2048,SUBTRACT,25000,0.000323094,77376857.640,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,p2048,2048,SUBTRACT,25000,0.002419333,10333426.722,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32,p2048,2048,SUBTRACT,25000,0.000163244,153144961.671,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32,p2048,2048,SUBTRACT,25000,0.002352252,10628112.882,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,p2048,2048,SUBTRACT,25000,0.000160226,156029592.300,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,p2048,2048,SUBTRACT,25000,0.002244890,11136402.985,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,p2048,2048,SUBTRACT,25000,0.000160233,156022847.295,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,p2048,2048,SUBTRACT,25000,0.002400599,10414067.687,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,p2048,2048,SUBTRACT,25000,0.000107008,233627448.470,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,p2048,2048,SUBTRACT,25000,0.002274926,10989368.469,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,p2048,2048,SUBTRACT,25000,0.000104693,238793430.490,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,p2048,2048,SUBTRACT,25000,0.002163796,11553769.463,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,p2048,2048,ADDMOD,25000,0.002295213,10892235.022,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,p2048,2048,ADDMOD,25000,0.000525165,47604086.335,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,p2048,2048,ADDMOD,25000,0.001874697,13335488.025,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,p2048,2048,ADDMOD,200000,0.000666144,300235384.541,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,p2048,2048,ADDMOD,25000,0.000862002,29002252.104,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,p2048,2048,ADDMOD,25000,0.002921273,8557912.936,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,p2048,2048,ADDMOD,25000,0.000436555,57266551.306,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,p2048,2048,ADDMOD,25000,0.002435972,10262843.816,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32,p2048,2048,ADDMOD,25000,0.000204273,122385219.058,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32,p2048,2048,ADDMOD,25000,0.002332588,10717709.444,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,p2048,2048,ADDMOD,25000,0.000228176,109564541.940,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,p2048,2048,ADDMOD,25000,0.002294307,10896536.549,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,p2048,2048,ADDMOD,25000,0.000219184,114059433.942,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,p2048,2048,ADDMOD,25000,0.002418371,10337536.983,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,p2048,2048,ADDMOD,25000,0.000099471,251329499.052,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,p2048,2048,ADDMOD,25000,0.002352941,10625000.686,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,p2048,2048,ADDMOD,25000,0.000098561,253650036.822,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,p2048,2048,ADDMOD,25000,0.002143491,11663216.680,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,p2048,2048,SUBTRACTMOD,25000,0.002013934,12413515.279,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,p2048,2048,SUBTRACTMOD,25000,0.000448082,55793360.293,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,p2048,2048,SUBTRACTMOD,25000,0.001850541,13509562.582,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,p2048,2048,SUBTRACTMOD,200000,0.000666048,300278658.595,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,p2048,2048,SUBTRACTMOD,25000,0.001011232,24722319.589,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,p2048,2048,SUBTRACTMOD,25000,0.003104114,8053827.837,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,p2048,2048,SUBTRACTMOD,25000,0.000472898,52865521.349,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,p2048,2048,SUBTRACTMOD,25000,0.002455414,10181582.489,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32,p2048,2048,SUBTRACTMOD,25000,0.000228412,109451322.230,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32,p2048,2048,SUBTRACTMOD,25000,0.002361393,10586971.290,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,p2048,2048,SUBTRACTMOD,25000,0.000223938,111638029.479,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,p2048,2048,SUBTRACTMOD,25000,0.002267335,11026160.866,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,p2048,2048,SUBTRACTMOD,25000,0.000220508,113374589.222,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,p2048,2048,SUBTRACTMOD,25000,0.002421801,10322896.151,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,p2048,2048,SUBTRACTMOD,25000,0.000106796,234091203.936,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,p2048,2048,SUBTRACTMOD,25000,0.002285711,10937515.715,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,p2048,2048,SUBTRACTMOD,25000,0.000098261,254424460.003,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,p2048,2048,SUBTRACTMOD,25000,0.002127141,11752864.563,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,p2048,2048,MULTIPLYOPERANDSCANNING,25000,0.012583896,1986666.134,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,p2048,2048,MULTIPLYOPERANDSCANNING,25000,0.002945331,8488010.226,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,p2048,2048,MULTIPLYOPERANDSCANNING,25000,0.003549233,7043775.332,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,25000,0.202767155,123294.130,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,25000,0.205736187,121514.841,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,p2048,2048,MULTIPLYOPERANDSCANNING,25000,0.081698787,306002.095,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,p2048,2048,MULTIPLYOPERANDSCANNING,25000,0.084451031,296029.542,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32,p2048,2048,MULTIPLYOPERANDSCANNING,25000,0.018046715,1385293.667,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32,p2048,2048,MULTIPLYOPERANDSCANNING,25000,0.020774404,1203403.955,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,p2048,2048,MULTIPLYOPERANDSCANNING,25000,0.001161582,21522371.867,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,p2048,2048,MULTIPLYOPERANDSCANNING,25000,0.003862604,6472317.583,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,p2048,2048,MULTIPLYOPERANDSCANNING,25000,0.001167324,21416504.347,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,p2048,2048,MULTIPLYOPERANDSCANNING,25000,0.004074631,6135524.871,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,p2048,2048,MULTIPLYOPERANDSCANNING,25000,0.001225775,20395260.286,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,p2048,2048,MULTIPLYOPERANDSCANNING,25000,0.004161290,6007752.474,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,p2048,2048,MULTIPLYOPERANDSCANNING,25000,0.001222660,20447221.492,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,p2048,2048,MULTIPLYOPERANDSCANNING,25000,0.003947010,6333908.403,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,p2048,2048,MULTIPLYPRODUCTSCANNING,25000,0.012585669,1986386.262,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,25000,0.002937595,8510363.038,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,25000,0.003552859,7036586.574,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,p2048,2048,MULTIPLYPRODUCTSCANNING,200000,0.001261440,158548959.919,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,25000,0.066451464,376214.435,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,25000,0.069474591,359843.788,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,p2048,2048,MULTIPLYPRODUCTSCANNING,25000,0.017005721,1470093.506,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,p2048,2048,MULTIPLYPRODUCTSCANNING,25000,0.019739155,1266518.247,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32,p2048,2048,MULTIPLYPRODUCTSCANNING,25000,0.004260188,5868285.588,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32,p2048,2048,MULTIPLYPRODUCTSCANNING,25000,0.007072181,3534977.412,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,p2048,2048,MULTIPLYPRODUCTSCANNING,25000,0.004849499,5155171.726,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,p2048,2048,MULTIPLYPRODUCTSCANNING,25000,0.007587440,3294918.962,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,p2048,2048,MULTIPLYPRODUCTSCANNING,25000,0.004390115,5694611.675,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,p2048,2048,MULTIPLYPRODUCTSCANNING,25000,0.007297843,3425669.751,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,p2048,2048,MULTIPLYPRODUCTSCANNING,25000,0.002564876,9747059.931,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,p2048,2048,MULTIPLYPRODUCTSCANNING,25000,0.005385161,4642386.712,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,p2048,2048,MULTIPLYPRODUCTSCANNING,25000,0.002218031,11271257.987,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,p2048,2048,MULTIPLYPRODUCTSCANNING,25000,0.004934318,5066556.355,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,p2048,2048,MONTGOMERYMULTIPLICATION,25000,0.079354739,315041.046,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,p2048,2048,MONTGOMERYMULTIPLICATION,25000,0.019021130,1314327.799,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,p2048,2048,MONTGOMERYMULTIPLICATION,25000,0.005343789,4678328.428,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,p2048,2048,MONTGOMERYMULTIPLICATION,200000,0.001767424,113159038.239,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,25000,0.095581873,261555.870,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,25000,0.098026645,255032.701,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,p2048,2048,MONTGOMERYMULTIPLICATION,25000,0.002637764,9477724.202,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,p2048,2048,MONTGOMERYMULTIPLICATION,25000,0.004672907,5349988.774,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32,p2048,2048,MONTGOMERYMULTIPLICATION,25000,0.000476396,52477350.518,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32,p2048,2048,MONTGOMERYMULTIPLICATION,25000,0.002582018,9682349.335,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,p2048,2048,MONTGOMERYMULTIPLICATION,25000,0.000799070,31286370.563,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,p2048,2048,MONTGOMERYMULTIPLICATION,25000,0.002807465,8904830.595,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,p2048,2048,MONTGOMERYMULTIPLICATION,25000,0.000675706,36998338.516,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,p2048,2048,MONTGOMERYMULTIPLICATION,25000,0.002902381,8613617.561,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,p2048,2048,MONTGOMERYMULTIPLICATION,25000,0.000706062,35407654.822,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,p2048,2048,MONTGOMERYMULTIPLICATION,25000,0.002909261,8593247.461,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,p2048,2048,MONTGOMERYMULTIPLICATION,25000,0.000604100,41383876.738,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,p2048,2048,MONTGOMERYMULTIPLICATION,25000,0.002655851,9413178.625,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,p2048,2048,COMPARE,25000,0.000352338,70954603.510,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,p2048,2048,COMPARE,25000,0.000063752,392144567.410,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,p2048,2048,COMPARE,25000,0.000153644,162713794.901,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,p2048,2048,COMPARE,200000,0.000666624,300019201.229,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,p2048,2048,COMPARE,25000,0.000551770,45308735.717,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,p2048,2048,COMPARE,25000,0.002907203,8599330.565,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,p2048,2048,COMPARE,25000,0.000274633,91030575.612,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,p2048,2048,COMPARE,25000,0.002256797,11077646.838,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,p2048,2048,COMPARE,25000,0.000160277,155979953.645,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,p2048,2048,COMPARE,25000,0.002248083,11120585.808,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,p2048,2048,COMPARE,25000,0.000156570,159672996.730,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,p2048,2048,COMPARE,25000,0.002358325,10600744.171,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,p2048,2048,COMPARE,25000,0.000043057,580626174.746,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,p2048,2048,COMPARE,25000,0.002197323,11377480.597,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,p2048,2048,COMPARE,25000,0.000042633,586399696.353,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,p2048,2048,COMPARE,25000,0.002067330,12092892.751,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,p2048,2048,REDUCE,3125,0.000102565,30468488.292,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,p2048,2048,REDUCE,3125,0.000027490,113677846.335,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,p2048,2048,REDUCE,3125,0.000213454,14640155.850,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,p2048,2048,REDUCE,200000,0.000659584,303221424.413,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,p2048,2048,REDUCE,20000,0.508415892,39337.873,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,p2048,2048,REDUCE,20000,0.510577258,39171.349,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,p2048,2048,REDUCE,20000,0.007853402,2546667.037,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,p2048,2048,REDUCE,20000,0.009546414,2095027.511,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,p2048,2048,REDUCE,20000,0.001593718,12549271.924,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,p2048,2048,REDUCE,20000,0.003400086,5882204.125,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,p2048,2048,REDUCE,20000,0.001794886,11142769.202,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,p2048,2048,REDUCE,20000,0.003490776,5729385.102,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,p2048,2048,REDUCE,20000,0.001435084,13936466.091,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,p2048,2048,REDUCE,20000,0.003151652,6345878.297,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,p2048,2048,REDUCE,20000,0.001718934,11635118.145,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,p2048,2048,REDUCE,20000,0.003402450,5878117.230,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,p2048,2048,MODMUL,1562,0.002074366,753001.148,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,p2048,2048,MODMUL,1562,0.000509442,3066099.743,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,p2048,2048,MODMUL,1562,0.001307487,1194658.129,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,p2048,2048,MODMUL,200000,0.009629696,20769087.622,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,p2048,2048,MODMUL,20000,0.847757633,23591.648,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,p2048,2048,MODMUL,20000,0.850148124,23525.312,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,p2048,2048,MODMUL,20000,0.021420110,933702.020,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,p2048,2048,MODMUL,20000,0.023148880,863972.685,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,p2048,2048,MODMUL,20000,0.009270714,2157331.137,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,p2048,2048,MODMUL,20000,0.010953001,1825983.585,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,p2048,2048,MODMUL,20000,0.007241829,2761733.255,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,p2048,2048,MODMUL,20000,0.008933458,2238774.725,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,p2048,2048,MODMUL,20000,0.008318424,2404301.586,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,p2048,2048,MODMUL,20000,0.010082813,1983573.437,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,p2048,2048,MODMUL,20000,0.007015392,2850874.188,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,p2048,2048,MODMUL,20000,0.008640668,2314635.858,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,p2048,2048,MODEXP,390,0.896729789,434.914,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,p2048,2048,MODEXP,390,0.239537489,1628.138,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,p2048,2048,MODEXP,390,0.148344923,2629.008,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,p2048,2048,MODEXP,200000,13.813170433,14478.935,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,p2048,2048,MODEXP,20000,131.869267422,151.665,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,p2048,2048,MODEXP,20000,131.844318098,151.694,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,p2048,2048,MODEXP,20000,15.083813892,1325.925,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,p2048,2048,MODEXP,20000,15.043730703,1329.457,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,p2048,2048,MODEXP,20000,1.751203719,11420.716,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,p2048,2048,MODEXP,20000,1.753364211,11406.643,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,p2048,2048,MODEXP,20000,4.870471476,4106.379,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,p2048,2048,MODEXP,20000,4.864875965,4111.102,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,p2048,2048,MODEXP,20000,1.607088357,12444.866,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,p2048,2048,MODEXP,20000,1.609037366,12429.792,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,p2048,2048,MODEXP,20000,5.355332584,3734.595,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,p2048,2048,MODEXP,20000,5.360876081,3730.734,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,p2048,2048,EXPONENTIATION,390,0.102649677,3799.330,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,p2048,2048,EXPONENTIATION,390,0.027374997,14246.577,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,p2048,2048,EXPONENTIATION,390,0.371798696,1048.955,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,p2048,2048,EXPONENTIATION,20000,28.874721708,692.647,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,p2048,2048,EXPONENTIATION,20000,28.847739546,693.295,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,p2048,2048,EXPONENTIATION,20000,5.206168970,3841.596,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,p2048,2048,EXPONENTIATION,20000,5.216705413,3833.837,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,p2048,2048,EXPONENTIATION,20000,1.485470173,13463.751,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,p2048,2048,EXPONENTIATION,20000,1.486140267,13457.680,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,p2048,2048,EXPONENTIATION,20000,1.440578103,13883.315,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,p2048,2048,EXPONENTIATION,20000,1.438746566,13900.989,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,p2048,2048,EXPONENTIATION,20000,1.487892463,13441.832,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,p2048,2048,EXPONENTIATION,20000,1.490678534,13416.709,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,p2048,2048,EXPONENTIATION,20000,1.424652141,14038.515,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,p2048,2048,EXPONENTIATION,20000,1.426505712,14020.273,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,p2048,2048,DIVIDE,3125,0.000188855,16547087.382,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,p2048,2048,DIVIDE,3125,0.000049107,63636536.392,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,p2048,2048,DIVIDE,3125,0.000209862,14890736.511,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,p2048,2048,DIVIDE,200000,0.000676832,295494302.870,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,p2048,2048,DIVIDE,20000,1.338371004,14943.540,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,p2048,2048,DIVIDE,20000,1.340777037,14916.723,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,p2048,2048,DIVIDE,20000,0.306928024,65161.857,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,p2048,2048,DIVIDE,20000,0.309121111,64699.560,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,p2048,2048,DIVIDE,20000,0.026229514,762499.832,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,p2048,2048,DIVIDE,20000,0.028707631,696678.873,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,p2048,2048,DIVIDE,20000,0.018971418,1054217.456,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,p2048,2048,DIVIDE,20000,0.021161919,945093.874,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,p2048,2048,DIVIDE,20000,0.022199600,900917.133,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,p2048,2048,DIVIDE,20000,0.024578379,813723.313,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,p2048,2048,DIVIDE,20000,0.022986804,870064.409,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,p2048,2048,DIVIDE,20000,0.025334290,789443.872,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,p2048,2048,ISQRT,781,0.000407643,1915892.205,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,p2048,2048,ISQRT,781,0.000086393,9040080.846,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,p2048,2048,ISQRT,20000,17.205783399,1162.400,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,p2048,2048,ISQRT,20000,17.203944969,1162.524,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,p2048,2048,ISQRT,20000,5.718122349,3497.652,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,p2048,2048,ISQRT,20000,5.703773955,3506.450,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,p2048,2048,ISQRT,20000,0.213349783,93742.772,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,p2048,2048,ISQRT,20000,0.215201533,92936.141,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,p2048,2048,ISQRT,20000,0.076244395,262314.364,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,p2048,2048,ISQRT,20000,0.077906060,256719.439,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,p2048,2048,ISQRT,20000,0.215538885,92790.681,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,p2048,2048,ISQRT,20000,0.217496811,91955.371,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,p2048,2048,ISQRT,20000,0.075731505,264090.883,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,p2048,2048,ISQRT,20000,0.077622531,257657.149,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-1t,p2048,2048,MODMUL_R2,25000,0.033943011,736528.648,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,gmp-nt,p2048,2048,MODMUL_R2,25000,0.008131855,3074329.292,0
library,Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz,host-cpu,openssl-nt,p2048,2048,MODMUL_R2,25000,0.020626474,1212034.589,0
library,NVIDIA GeForce RTX 3060,gpu,cgbn,p2048,2048,MODMUL_R2,200000,0.003708928,53923937.051,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w8,p2048,2048,MODMUL_R2,25000,0.074720552,334579.969,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w8,p2048,2048,MODMUL_R2,25000,0.077265769,323558.548,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w16,p2048,2048,MODMUL_R2,25000,0.003119439,8014261.531,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w16,p2048,2048,MODMUL_R2,25000,0.005348977,4673790.901,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-opt,p2048,2048,MODMUL_R2,25000,0.001302184,19198515.530,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-opt,p2048,2048,MODMUL_R2,25000,0.003624167,6898136.864,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-o64,p2048,2048,MODMUL_R2,25000,0.001025061,24388792.313,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-o64,p2048,2048,MODMUL_R2,25000,0.003129499,7988499.182,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il,p2048,2048,MODMUL_R2,25000,0.001183641,21121269.595,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il,p2048,2048,MODMUL_R2,25000,0.003502356,7138052.204,0
opencl-kernel,NVIDIA GeForce RTX 3060,GPU,w32-il64,p2048,2048,MODMUL_R2,25000,0.000971029,25745884.305,0
opencl-e2e,NVIDIA GeForce RTX 3060,GPU,w32-il64,p2048,2048,MODMUL_R2,25000,0.003229975,7739997.983,0
```

## CGBN comparison

Merged from `cgbn_results.tsv` after the sweep, by `cgbn_merge.py`.
MPA columns are the fastest correct kernel on the GPU, kernel time only.

### secp256k1 (256-bit)

| Operation | best MPA kernel | MPA items | MPA ops/s | CGBN ops/s | CGBN / MPA |
|---|---|---|---|---|---|
| ADD | w32-il64 | 200000 | 2.99 G | 1.63 G | 0.54x |
| SUBTRACT | w32-il64 | 200000 | 2.99 G | 1.66 G | 0.55x |
| ADDMOD | w32-il64 | 200000 | 3.10 G | 1.63 G | 0.53x |
| SUBTRACTMOD | w32-il64 | 200000 | 3.10 G | 1.63 G | 0.53x |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 200000 | 1.53 G | 1.64 G | 1.07x |
| MONTGOMERYMULTIPLICATION | w32 | 200000 | 2.86 G | 1.60 G | 0.56x |
| COMPARE | w32-il64 | 200000 | 6.34 G | 1.64 G | 0.26x |
| REDUCE | w32-il64 | 25000 | 393.36 M | 1.30 G | 3.31x |
| MODMUL | w32-il64 | 20000 | 155.97 M | 352.55 M | 2.26x |
| MODEXP | w32-il64 | 20000 | 10.53 M | 1.03 M | 0.10x |
| DIVIDE | w32-il64 | 25000 | 222.84 M | 981.47 M | 4.40x |
| MODMUL_R2 | w32-il64 | 200000 | 1.31 G | 1.18 G | 0.90x |

### rsa256(composite) (256-bit)

| Operation | best MPA kernel | MPA items | MPA ops/s | CGBN ops/s | CGBN / MPA |
|---|---|---|---|---|---|
| ADD | w32-il64 | 200000 | 2.96 G | 1.64 G | 0.55x |
| SUBTRACT | w32-il64 | 200000 | 2.97 G | 1.64 G | 0.55x |
| ADDMOD | w32-il64 | 200000 | 3.11 G | 1.63 G | 0.52x |
| SUBTRACTMOD | w32-il64 | 200000 | 3.09 G | 1.63 G | 0.53x |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 200000 | 1.44 G | 1.64 G | 1.14x |
| MONTGOMERYMULTIPLICATION | w32 | 200000 | 2.79 G | 1.64 G | 0.59x |
| COMPARE | w32-il64 | 200000 | 6.32 G | 1.64 G | 0.26x |
| REDUCE | w32-il64 | 25000 | 392.10 M | 1.30 G | 3.32x |
| MODMUL | w32-il64 | 20000 | 155.74 M | 352.55 M | 2.26x |
| MODEXP | w32-il64 | 20000 | 10.56 M | 1.05 M | 0.10x |
| DIVIDE | w32-il64 | 25000 | 212.95 M | 978.55 M | 4.60x |
| MODMUL_R2 | w32-il64 | 200000 | 1.31 G | 1.18 G | 0.91x |

### brainpoolP512r1 (512-bit)

| Operation | best MPA kernel | MPA items | MPA ops/s | CGBN ops/s | CGBN / MPA |
|---|---|---|---|---|---|
| ADD | w32-il | 100000 | 1.32 G | 1.09 G | 0.82x |
| SUBTRACT | w32-il | 100000 | 1.31 G | 1.08 G | 0.82x |
| ADDMOD | w32-il64 | 100000 | 1.36 G | 1.08 G | 0.79x |
| SUBTRACTMOD | w32-il64 | 100000 | 1.37 G | 1.08 G | 0.79x |
| MULTIPLYPRODUCTSCANNING | w32-il | 100000 | 317.07 M | 1.09 G | 3.42x |
| MONTGOMERYMULTIPLICATION | w32 | 100000 | 971.95 M | 1.10 G | 1.13x |
| COMPARE | w32-il | 100000 | 3.44 G | 1.09 G | 0.32x |
| REDUCE | w32-il | 20000 | 150.90 M | 710.31 M | 4.71x |
| MODMUL | w32-o64 | 20000 | 46.80 M | 107.55 M | 2.30x |
| MODEXP | w32-o64 | 20000 | 1.51 M | 426.23 k | 0.28x |
| DIVIDE | w32-il64 | 20000 | 57.29 M | 518.54 M | 9.05x |
| MODMUL_R2 | w32-o64 | 100000 | 523.67 M | 737.03 M | 1.41x |

### p1024 (1024-bit)

| Operation | best MPA kernel | MPA items | MPA ops/s | CGBN ops/s | CGBN / MPA |
|---|---|---|---|---|---|
| ADD | w32-il | 50000 | 476.19 M | 590.96 M | 1.24x |
| SUBTRACT | w32-il | 50000 | 476.02 M | 591.86 M | 1.24x |
| ADDMOD | w32-il | 50000 | 568.34 M | 590.96 M | 1.04x |
| SUBTRACTMOD | w32-il | 50000 | 571.96 M | 592.14 M | 1.04x |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 50000 | 60.76 M | 567.77 M | 9.34x |
| MONTGOMERYMULTIPLICATION | w32-il64 | 50000 | 235.55 M | 355.15 M | 1.51x |
| COMPARE | w32-il64 | 50000 | 1.57 G | 590.07 M | 0.38x |
| REDUCE | w32-il | 20000 | 44.69 M | 440.36 M | 9.85x |
| MODMUL | w32-il64 | 20000 | 11.84 M | 38.44 M | 3.25x |
| MODEXP | w32-o64 | 20000 | 162.64 k | 82.19 k | 0.51x |
| DIVIDE | w32-o64 | 20000 | 13.11 M | 337.91 M | 25.77x |
| MODMUL_R2 | w32-il64 | 50000 | 104.68 M | 196.76 M | 1.88x |

### p2048 (2048-bit)

| Operation | best MPA kernel | MPA items | MPA ops/s | CGBN ops/s | CGBN / MPA |
|---|---|---|---|---|---|
| ADD | w32-il64 | 25000 | 239.43 M | 300.06 M | 1.25x |
| SUBTRACT | w32-il64 | 25000 | 238.79 M | 300.02 M | 1.26x |
| ADDMOD | w32-il64 | 25000 | 253.65 M | 300.24 M | 1.18x |
| SUBTRACTMOD | w32-il64 | 25000 | 254.42 M | 300.28 M | 1.18x |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 25000 | 11.27 M | 158.55 M | 14.07x |
| MONTGOMERYMULTIPLICATION | w32 | 25000 | 52.48 M | 113.16 M | 2.16x |
| COMPARE | w32-il64 | 25000 | 586.40 M | 300.02 M | 0.51x |
| REDUCE | w32-il | 20000 | 13.94 M | 303.22 M | 21.76x |
| MODMUL | w32-il64 | 20000 | 2.85 M | 20.77 M | 7.29x |
| MODEXP | w32-il | 20000 | 12.44 k | 14.48 k | 1.16x |
| DIVIDE | w32-o64 | 20000 | 1.05 M | 295.49 M | 280.30x |
| MODMUL_R2 | w32-il64 | 25000 | 25.75 M | 53.92 M | 2.09x |
