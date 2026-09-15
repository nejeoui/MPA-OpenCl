# MPA-OpenCL benchmark report - gfx1010:xnack-

> **Note.** The multi-threaded GMP and OpenSSL baseline columns have been
> removed from this report: they predate the 2026-09-12 timing fix and were
> understated (see `reports/README.md`). The single-threaded GMP column, the
> OpenCL-on-CPU rows and all MPA measurements are unaffected and were verified
> against GMP before timing.


## 1. System under test

1 OpenCL device(s) exercised with the identical kernels and operands.

### Device 0 - gfx1010:xnack- (GPU)

| Property | Value |
|---|---|
| Model | gfx1010:xnack- |
| Type | GPU |
| Vendor | Advanced Micro Devices, Inc. |
| Device memory | 7.98 GiB |
| Max single allocation | 6.79 GiB |
| Local memory | 64 KiB |
| Global cache | 16 KiB |
| Compute units | 20 |
| Max clock | 2100 MHz |
| Max work-group size | 256 |
| OpenCL version | OpenCL 2.0  |
| Driver | 3649.0 (HSA1.1,LC) |

### Host

| Property | Value |
|---|---|
| CPU | AMD Ryzen 7 5700G with Radeon Graphics |
| Logical cores | 16 |
| OpenMP threads used | 16 |
| RAM | 11.7 GB |
| OS | Linux Mint 22.2 |
| Kernel | 6.17.10-061710-generic |
| Arch | x86_64 |
| GMP | 6.3.0 |
| OpenSSL | OpenSSL 3.0.13 30 Jan 2024 |
| CGBN | not measured |

## 2. Method

- Workload auto-sized from the device and host: --min-items from 700 x compute units, --items from ten times that capped by host RAM. Either flag, given explicitly, overrides its half.
- Base workload 20000 items, scaled down per operator by its cost weight and by modulus size. Device rows honour --min-items (2000) so the GPU is not left idle; the CPU libraries keep the smaller count because a full-width MODEXP there costs minutes. Both counts appear in every row as dev/cpu, and throughput is per-second so they remain comparable.
- 5 timed repetitions, **minimum** reported. Two untimed warm-up launches precede them.
- `kernel` times `clEnqueueNDRangeKernel` + `clFinish` only. `e2e` adds the host->device operand writes and the device->host result read.
- Every OpenCL device runs the same kernels on the same operands, so GPU and CPU-OpenCL columns are directly comparable.
- CPU library baselines (GMP, OpenSSL) run those same operands, with temporaries preallocated outside the timed region, so the figure is the arithmetic and not marshalling. The generator is reseeded per modulus and operation so every backend sees identical inputs.
- OpenSSL rows time the nearest BN primitive, which is not always semantically identical (its Montgomery routine expects Montgomery-domain inputs); they measure comparable work, not identical results. Correctness is judged against GMP only.
- Every device cell is checked word-for-word against GMP before it is timed. A cell that mismatches is reported and excluded from the speedup tables.
- Total wall time 1362.3 s.

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

### Device 0 - gfx1010:xnack- (GPU)

#### secp256k1 (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 20000 / 20000 | 225.79 M | 243.45 M | 292.10 M | 279.59 M | 284.73 M | 404.19 M | 317.62 M | 81.63 M | n/a |
| SUBTRACT | 20000 / 20000 | 218.50 M | 250.28 M | 292.96 M | 262.46 M | 325.23 M | 314.57 M | 315.27 M | 102.58 M | n/a |
| ADDMOD | 20000 / 20000 | 223.29 M | 226.15 M | 262.80 M | 284.16 M | 334.27 M | 335.62 M | 305.42 M | 30.32 M | n/a |
| SUBTRACTMOD | 20000 / 20000 | 224.35 M | 209.36 M | 268.86 M | 301.37 M | 321.00 M | 467.62 M | 318.13 M | 37.44 M | n/a |
| MULTIPLYOPERANDSCANNING | 20000 / 20000 | 10.58 M | 33.02 M | 94.79 M | 323.44 M | 301.37 M | 306.93 M | 320.12 M | 72.52 M | n/a |
| MULTIPLYPRODUCTSCANNING | 20000 / 20000 | 58.27 M | 81.52 M | 175.28 M | 230.19 M | 177.43 M | 258.35 M | 262.49 M | 72.28 M | n/a |
| MONTGOMERYMULTIPLICATION | 20000 / 20000 | 133.67 M | 250.79 M | 290.49 M | 280.61 M | 340.49 M | 295.70 M | 297.81 M | 8.51 M | n/a |
| COMPARE | 20000 / 20000 | 298.44 M | 371.68 M | - | 410.42 M | 279.51 M | 280.89 M | 319.56 M | 167.53 M | n/a |
| REDUCE | 2500 / 2500 | 19.74 M | 19.99 M | - | 28.57 M | 24.50 M | 29.15 M | 24.69 M | 48.48 M | n/a |
| MODMUL | 2000 / 1250 | 7.73 M | 12.54 M | - | 14.73 M | 15.33 M | 14.82 M | 15.15 M | 10.72 M | n/a |
| MODEXP | 2000 / 312 | 310.84 k | 869.83 k | - | 1.88 M | 2.49 M | 1.97 M | 2.42 M | 138.65 k | n/a |
| EXPONENTIATION | 2000 / 312 | 74.07 k | 5.82 M | - | 8.10 M | 10.09 M | 8.62 M | 10.91 M | 310.33 k | n/a |
| DIVIDE | 2500 / 2500 | 4.69 M | 10.33 M | - | 16.21 M | 16.80 M | 18.69 M | 15.90 M | 26.20 M | n/a |
| ISQRT | 2000 / 625 | 808.75 k | 1.64 M | - | 2.97 M | 2.90 M | 2.94 M | 3.09 M | 14.91 M | n/a |
| MODMUL_R2 | 20000 / 20000 | 134.62 M | 208.40 M | - | 272.71 M | 242.24 M | 230.81 M | 265.42 M | 15.08 M | n/a |

#### rsa256(composite) (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 20000 / 20000 | 224.17 M | 323.54 M | 301.09 M | 268.97 M | 283.97 M | 320.22 M | 311.28 M | 81.56 M | n/a |
| SUBTRACT | 20000 / 20000 | 227.10 M | 299.33 M | 282.19 M | 261.53 M | 257.38 M | 306.22 M | 309.65 M | 100.12 M | n/a |
| ADDMOD | 20000 / 20000 | 213.09 M | 258.18 M | 253.81 M | 312.99 M | 306.84 M | 332.93 M | 252.28 M | 35.65 M | n/a |
| SUBTRACTMOD | 20000 / 20000 | 209.80 M | 239.27 M | 259.39 M | 273.65 M | 401.50 M | 335.73 M | 251.23 M | 36.70 M | n/a |
| MULTIPLYOPERANDSCANNING | 20000 / 20000 | 10.58 M | 29.85 M | 108.29 M | 304.68 M | 299.83 M | 376.44 M | 291.55 M | 71.17 M | n/a |
| MULTIPLYPRODUCTSCANNING | 20000 / 20000 | 49.24 M | 82.34 M | 192.26 M | 226.41 M | 169.84 M | 255.67 M | 230.41 M | 71.58 M | n/a |
| MONTGOMERYMULTIPLICATION | 20000 / 20000 | 153.42 M | 303.10 M | 312.60 M | 291.47 M | 259.46 M | 299.78 M | 327.10 M | 8.51 M | n/a |
| COMPARE | 20000 / 20000 | 302.46 M | 305.19 M | - | 408.81 M | 361.97 M | 319.81 M | 281.88 M | 159.55 M | n/a |
| REDUCE | 2500 / 2500 | 18.94 M | 20.96 M | - | 24.31 M | 24.01 M | 24.71 M | 25.82 M | 32.32 M | n/a |
| MODMUL | 2000 / 1250 | 8.21 M | 11.35 M | - | 15.19 M | 16.12 M | 15.00 M | 14.44 M | 10.66 M | n/a |
| MODEXP | 2000 / 312 | 311.33 k | 866.78 k | - | 1.94 M | 2.50 M | 1.93 M | 2.40 M | 145.76 k | n/a |
| EXPONENTIATION | 2000 / 312 | 74.49 k | 6.46 M | - | 8.10 M | 10.12 M | 8.66 M | 9.49 M | 310.10 k | n/a |
| DIVIDE | 2500 / 2500 | 4.61 M | 9.85 M | - | 15.93 M | 16.06 M | 16.08 M | 15.64 M | 25.66 M | n/a |
| ISQRT | 2000 / 625 | 727.91 k | 1.52 M | - | 2.83 M | 2.85 M | 2.86 M | 2.85 M | 15.08 M | n/a |
| MODMUL_R2 | 20000 / 20000 | 159.34 M | 269.47 M | - | 249.09 M | 282.56 M | 306.74 M | 285.30 M | 15.09 M | n/a |

#### brainpoolP512r1 (512-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 10000 / 10000 | 80.38 M | 96.74 M | 128.76 M | 122.77 M | 128.28 M | 156.69 M | 158.11 M | 76.16 M | n/a |
| SUBTRACT | 10000 / 10000 | 89.47 M | 124.44 M | 112.45 M | 137.81 M | 125.71 M | 152.01 M | 178.53 M | 90.10 M | n/a |
| ADDMOD | 10000 / 10000 | 84.47 M | 102.16 M | 123.23 M | 158.08 M | 153.89 M | 171.44 M | 155.62 M | 31.23 M | n/a |
| SUBTRACTMOD | 10000 / 10000 | 86.53 M | 107.46 M | 124.33 M | 150.73 M | 198.32 M | 218.70 M | 157.33 M | 33.73 M | n/a |
| MULTIPLYOPERANDSCANNING | 10000 / 10000 | 1.35 M | 6.45 M | 26.17 M | 116.01 M | 122.80 M | 155.96 M | 159.88 M | 27.00 M | n/a |
| MULTIPLYPRODUCTSCANNING | 10000 / 10000 | 3.64 M | 12.88 M | 43.03 M | 42.05 M | 42.45 M | 92.93 M | 104.67 M | 26.83 M | n/a |
| MONTGOMERYMULTIPLICATION | 10000 / 10000 | 45.96 M | 97.29 M | 125.05 M | 130.65 M | 130.99 M | 117.44 M | 128.10 M | 3.47 M | n/a |
| COMPARE | 10000 / 10000 | 148.20 M | 149.35 M | - | 169.23 M | 153.44 M | 152.62 M | 152.88 M | 162.69 M | n/a |
| REDUCE | 2000 / 1250 | 7.69 M | 12.66 M | - | 13.57 M | 13.94 M | 13.77 M | 14.43 M | 30.18 M | n/a |
| MODMUL | 2000 / 625 | 2.79 M | 6.07 M | - | 7.71 M | 7.95 M | 8.00 M | 8.46 M | 6.07 M | n/a |
| MODEXP | 2000 / 156 | 47.92 k | 155.54 k | - | 252.15 k | 389.63 k | 244.77 k | 391.57 k | 25.92 k | n/a |
| EXPONENTIATION | 2000 / 156 | 6.67 k | 26.81 k | - | 424.87 k | 381.31 k | 423.93 k | 381.69 k | 114.58 k | n/a |
| DIVIDE | 2000 / 1250 | 532.01 k | 2.08 M | - | 4.90 M | 4.92 M | 4.45 M | 4.75 M | 23.22 M | n/a |
| ISQRT | 2000 / 312 | 67.12 k | 381.82 k | - | 892.93 k | 901.30 k | 894.36 k | 889.80 k | 8.73 M | n/a |
| MODMUL_R2 | 10000 / 10000 | 43.25 M | 106.53 M | - | 98.87 M | 123.27 M | 99.10 M | 141.82 M | 7.12 M | n/a |

#### p1024 (1024-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 5000 / 5000 | 29.23 M | 43.79 M | 55.08 M | 56.74 M | 58.20 M | 68.45 M | 88.03 M | 62.45 M | n/a |
| SUBTRACT | 5000 / 5000 | 28.88 M | 42.05 M | 55.64 M | 57.04 M | 66.44 M | 76.96 M | 84.19 M | 81.81 M | n/a |
| ADDMOD | 5000 / 5000 | 28.93 M | 44.21 M | 48.79 M | 77.28 M | 84.17 M | 79.19 M | 76.66 M | 23.85 M | n/a |
| SUBTRACTMOD | 5000 / 5000 | 28.44 M | 42.70 M | 55.66 M | 83.83 M | 83.78 M | 107.60 M | 77.78 M | 28.71 M | n/a |
| MULTIPLYOPERANDSCANNING | 5000 / 5000 | 260.21 k | 1.35 M | 4.29 M | 58.67 M | 26.49 M | 66.73 M | 54.87 M | 7.55 M | n/a |
| MULTIPLYPRODUCTSCANNING | 5000 / 5000 | 470.68 k | 1.82 M | 6.52 M | 6.55 M | 6.44 M | 26.09 M | 23.69 M | 7.66 M | n/a |
| MONTGOMERYMULTIPLICATION | 5000 / 5000 | 2.23 M | 37.99 M | 56.72 M | 46.89 M | 57.72 M | 50.47 M | 56.19 M | 1.11 M | n/a |
| COMPARE | 5000 / 5000 | 62.85 M | 66.03 M | - | 73.33 M | 75.12 M | 64.33 M | 86.37 M | 157.08 M | n/a |
| REDUCE | 2000 / 625 | 1.52 M | 6.47 M | - | 8.43 M | 8.26 M | 7.26 M | 8.73 M | 42.04 M | n/a |
| MODMUL | 2000 / 312 | 222.21 k | 1.94 M | - | 2.98 M | 3.31 M | 2.62 M | 3.27 M | 2.29 M | n/a |
| MODEXP | 2000 / 78 | 1.24 k | 20.14 k | - | 30.07 k | 52.03 k | 30.03 k | 49.06 k | 4.43 k | n/a |
| EXPONENTIATION | 2000 / 78 | 840.9 | 3.44 k | - | 13.39 k | 13.63 k | 13.33 k | 13.68 k | 32.09 k | n/a |
| DIVIDE | 2000 / 625 | 116.83 k | 360.76 k | - | 1.41 M | 1.34 M | 1.39 M | 1.31 M | 22.50 M | n/a |
| ISQRT | 2000 / 156 | 16.57 k | 29.92 k | - | 218.97 k | 216.82 k | 216.30 k | 213.24 k | 4.12 M | n/a |
| MODMUL_R2 | 5000 / 5000 | 4.59 M | 29.82 M | - | 44.57 M | 46.20 M | 38.00 M | 44.11 M | 2.52 M | n/a |

#### p2048 (2048-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T | CGBN |
|---|---|---|---|---|---|---|---|---|---|---|
| ADD | 2500 / 2500 | 3.96 M | 4.16 M | 21.43 M | 21.74 M | 19.41 M | 30.43 M | 33.83 M | 44.69 M | n/a |
| SUBTRACT | 2500 / 2500 | 4.25 M | 5.81 M | 18.55 M | 22.00 M | 19.62 M | 34.07 M | 37.67 M | 42.74 M | n/a |
| ADDMOD | 2500 / 2500 | 3.99 M | 4.34 M | 21.30 M | 32.14 M | 38.56 M | 35.49 M | 37.81 M | 18.64 M | n/a |
| SUBTRACTMOD | 2500 / 2500 | 3.38 M | 5.73 M | 19.15 M | 33.51 M | 42.76 M | 39.56 M | 38.46 M | 20.87 M | n/a |
| MULTIPLYOPERANDSCANNING | 2500 / 2500 | 20.96 k | 96.33 k | 629.47 k | 5.04 M | 3.79 M | 5.11 M | 12.45 M | 2.28 M | n/a |
| MULTIPLYPRODUCTSCANNING | 2500 / 2500 | 88.18 k | 289.78 k | 1.14 M | 921.54 k | 925.95 k | 3.68 M | 3.73 M | 2.28 M | n/a |
| MONTGOMERYMULTIPLICATION | 2500 / 2500 | 112.67 k | 452.25 k | 12.52 M | 12.94 M | 15.82 M | 12.87 M | 17.22 M | 334.67 k | n/a |
| COMPARE | 2500 / 2500 | 8.68 M | 8.60 M | - | 33.79 M | 33.13 M | 30.26 M | 29.36 M | 167.70 M | n/a |
| REDUCE | 2000 / 312 | 393.87 k | 1.29 M | - | 3.82 M | 3.06 M | 3.81 M | 3.21 M | 30.71 M | n/a |
| MODMUL | 2000 / 156 | 70.44 k | 287.14 k | - | 1.02 M | 969.46 k | 1.01 M | 997.69 k | 774.12 k | n/a |
| MODEXP | 2000 / 64 | 104.4 | 606.9 | - | 3.92 k | 6.92 k | 3.83 k | 6.91 k | 600.9 | n/a |
| EXPONENTIATION | 2000 / 64 | 100.2 | 420.6 | - | 1.72 k | 1.75 k | 1.71 k | 1.75 k | 4.79 k | n/a |
| DIVIDE | 2000 / 312 | 29.68 k | 52.39 k | - | 110.05 k | 139.49 k | 162.62 k | 115.85 k | 17.00 M | n/a |
| ISQRT | 2000 / 78 | 3.88 k | 6.88 k | - | 19.95 k | 21.95 k | 20.19 k | 21.54 k | 2.34 M | n/a |
| MODMUL_R2 | 2500 / 2500 | 139.09 k | 807.51 k | - | 7.99 M | 10.73 M | 8.14 M | 11.32 M | 774.34 k | n/a |

## 5. Head to head

Best OpenCL GPU result against best OpenCL CPU result and the CPU libraries.
Ratios above 1.00x mean the GPU is faster than that baseline.

### secp256k1 (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w32-il | 404.19 M | none | n/a | 81.63 M | n/a | n/a | n/a |
| SUBTRACT | w32-o64 | 325.23 M | none | n/a | 102.58 M | n/a | n/a | n/a |
| ADDMOD | w32-il | 335.62 M | none | n/a | 30.32 M | n/a | n/a | n/a |
| SUBTRACTMOD | w32-il | 467.62 M | none | n/a | 37.44 M | n/a | n/a | n/a |
| MULTIPLYOPERANDSCANNING | w32-opt | 323.44 M | none | n/a | 72.52 M | n/a | n/a | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 262.49 M | none | n/a | 72.28 M | n/a | n/a | n/a |
| MONTGOMERYMULTIPLICATION | w32-o64 | 340.49 M | none | n/a | 8.51 M | n/a | n/a | n/a |
| COMPARE | w32-opt | 410.42 M | none | n/a | 167.53 M | n/a | n/a | n/a |
| REDUCE | w32-il | 29.15 M | none | n/a | 48.48 M | n/a | n/a | n/a |
| MODMUL | w32-o64 | 9.58 M | none | n/a | 10.72 M | n/a | n/a | n/a |
| MODEXP | w32-o64 | 387.76 k | none | n/a | 138.65 k | n/a | n/a | n/a |
| EXPONENTIATION | w32-il64 | 1.70 M | none | n/a | 310.33 k | n/a | n/a | n/a |
| DIVIDE | w32-il | 18.69 M | none | n/a | 26.20 M | n/a | n/a | n/a |
| ISQRT | w32-il64 | 965.61 k | none | n/a | 14.91 M | n/a | n/a | n/a |
| MODMUL_R2 | w32-opt | 272.71 M | none | n/a | 15.08 M | n/a | n/a | n/a |

### rsa256(composite) (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w16 | 323.54 M | none | n/a | 81.56 M | n/a | n/a | n/a |
| SUBTRACT | w32-il64 | 309.65 M | none | n/a | 100.12 M | n/a | n/a | n/a |
| ADDMOD | w32-il | 332.93 M | none | n/a | 35.65 M | n/a | n/a | n/a |
| SUBTRACTMOD | w32-o64 | 401.50 M | none | n/a | 36.70 M | n/a | n/a | n/a |
| MULTIPLYOPERANDSCANNING | w32-il | 376.44 M | none | n/a | 71.17 M | n/a | n/a | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il | 255.67 M | none | n/a | 71.58 M | n/a | n/a | n/a |
| MONTGOMERYMULTIPLICATION | w32-il64 | 327.10 M | none | n/a | 8.51 M | n/a | n/a | n/a |
| COMPARE | w32-opt | 408.81 M | none | n/a | 159.55 M | n/a | n/a | n/a |
| REDUCE | w32-il64 | 25.82 M | none | n/a | 32.32 M | n/a | n/a | n/a |
| MODMUL | w32-o64 | 10.08 M | none | n/a | 10.66 M | n/a | n/a | n/a |
| MODEXP | w32-o64 | 390.45 k | none | n/a | 145.76 k | n/a | n/a | n/a |
| EXPONENTIATION | w32-o64 | 1.58 M | none | n/a | 310.10 k | n/a | n/a | n/a |
| DIVIDE | w32-il | 16.08 M | none | n/a | 25.66 M | n/a | n/a | n/a |
| ISQRT | w32-il | 893.98 k | none | n/a | 15.08 M | n/a | n/a | n/a |
| MODMUL_R2 | w32-il | 306.74 M | none | n/a | 15.09 M | n/a | n/a | n/a |

### brainpoolP512r1 (512-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w32-il64 | 158.11 M | none | n/a | 76.16 M | n/a | n/a | n/a |
| SUBTRACT | w32-il64 | 178.53 M | none | n/a | 90.10 M | n/a | n/a | n/a |
| ADDMOD | w32-il | 171.44 M | none | n/a | 31.23 M | n/a | n/a | n/a |
| SUBTRACTMOD | w32-il | 218.70 M | none | n/a | 33.73 M | n/a | n/a | n/a |
| MULTIPLYOPERANDSCANNING | w32-il64 | 159.88 M | none | n/a | 27.00 M | n/a | n/a | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 104.67 M | none | n/a | 26.83 M | n/a | n/a | n/a |
| MONTGOMERYMULTIPLICATION | w32-o64 | 130.99 M | none | n/a | 3.47 M | n/a | n/a | n/a |
| COMPARE | w32-opt | 169.23 M | none | n/a | 162.69 M | n/a | n/a | n/a |
| REDUCE | w32-il64 | 9.02 M | none | n/a | 30.18 M | n/a | n/a | n/a |
| MODMUL | w32-il64 | 2.64 M | none | n/a | 6.07 M | n/a | n/a | n/a |
| MODEXP | w32-il64 | 30.54 k | none | n/a | 25.92 k | n/a | n/a | n/a |
| EXPONENTIATION | w32-opt | 33.14 k | none | n/a | 114.58 k | n/a | n/a | n/a |
| DIVIDE | w32-o64 | 3.07 M | none | n/a | 23.22 M | n/a | n/a | n/a |
| ISQRT | w32-o64 | 140.60 k | none | n/a | 8.73 M | n/a | n/a | n/a |
| MODMUL_R2 | w32-il64 | 141.82 M | none | n/a | 7.12 M | n/a | n/a | n/a |

### p1024 (1024-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w32-il64 | 88.03 M | none | n/a | 62.45 M | n/a | n/a | n/a |
| SUBTRACT | w32-il64 | 84.19 M | none | n/a | 81.81 M | n/a | n/a | n/a |
| ADDMOD | w32-o64 | 84.17 M | none | n/a | 23.85 M | n/a | n/a | n/a |
| SUBTRACTMOD | w32-il | 107.60 M | none | n/a | 28.71 M | n/a | n/a | n/a |
| MULTIPLYOPERANDSCANNING | w32-il | 66.73 M | none | n/a | 7.55 M | n/a | n/a | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il | 26.09 M | none | n/a | 7.66 M | n/a | n/a | n/a |
| MONTGOMERYMULTIPLICATION | w32-o64 | 57.72 M | none | n/a | 1.11 M | n/a | n/a | n/a |
| COMPARE | w32-il64 | 86.37 M | none | n/a | 157.08 M | n/a | n/a | n/a |
| REDUCE | w32-il64 | 2.73 M | none | n/a | 42.04 M | n/a | n/a | n/a |
| MODMUL | w32-o64 | 516.79 k | none | n/a | 2.29 M | n/a | n/a | n/a |
| MODEXP | w32-o64 | 2.03 k | none | n/a | 4.43 k | n/a | n/a | n/a |
| EXPONENTIATION | w32-il64 | 533.4 | none | n/a | 32.09 k | n/a | n/a | n/a |
| DIVIDE | w32-opt | 441.65 k | none | n/a | 22.50 M | n/a | n/a | n/a |
| ISQRT | w32-opt | 17.08 k | none | n/a | 4.12 M | n/a | n/a | n/a |
| MODMUL_R2 | w32-o64 | 46.20 M | none | n/a | 2.52 M | n/a | n/a | n/a |

### p2048 (2048-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | CGBN | GPU vs CPU-CL | GPU vs CGBN |
|---|---|---|---|---|---|---|---|---|
| ADD | w32-il64 | 33.83 M | none | n/a | 44.69 M | n/a | n/a | n/a |
| SUBTRACT | w32-il64 | 37.67 M | none | n/a | 42.74 M | n/a | n/a | n/a |
| ADDMOD | w32-o64 | 38.56 M | none | n/a | 18.64 M | n/a | n/a | n/a |
| SUBTRACTMOD | w32-o64 | 42.76 M | none | n/a | 20.87 M | n/a | n/a | n/a |
| MULTIPLYOPERANDSCANNING | w32-il64 | 12.45 M | none | n/a | 2.28 M | n/a | n/a | n/a |
| MULTIPLYPRODUCTSCANNING | w32-il64 | 3.73 M | none | n/a | 2.28 M | n/a | n/a | n/a |
| MONTGOMERYMULTIPLICATION | w32-il64 | 17.22 M | none | n/a | 334.67 k | n/a | n/a | n/a |
| COMPARE | w32-opt | 33.79 M | none | n/a | 167.70 M | n/a | n/a | n/a |
| REDUCE | w32-opt | 595.71 k | none | n/a | 30.71 M | n/a | n/a | n/a |
| MODMUL | w32-opt | 79.37 k | none | n/a | 774.12 k | n/a | n/a | n/a |
| MODEXP | w32-o64 | 221.6 | none | n/a | 600.9 | n/a | n/a | n/a |
| EXPONENTIATION | w32-il64 | 56.2 | none | n/a | 4.79 k | n/a | n/a | n/a |
| DIVIDE | w32-il | 25.37 k | none | n/a | 17.00 M | n/a | n/a | n/a |
| ISQRT | w32-o64 | 855.9 | none | n/a | 2.34 M | n/a | n/a | n/a |
| MODMUL_R2 | w32-il64 | 11.32 M | none | n/a | 774.34 k | n/a | n/a | n/a |

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

Also written to `gfx1010_xnack-_Report.csv` for analysis.

```csv
kind,device,device_type,kernel,modulus,bits,operation,items,seconds,ops_per_sec,mismatches
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,secp256k1,256,ADD,20000,0.000245000,81632653.078,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,secp256k1,256,ADD,20000,0.000039784,502714658.866,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,secp256k1,256,ADD,20000,0.000049032,407896883.729,0
opencl-kernel,gfx1010:xnack-,GPU,w8,secp256k1,256,ADD,20000,0.000088576,225794797.553,0
opencl-e2e,gfx1010:xnack-,GPU,w8,secp256k1,256,ADD,20000,0.000354445,56426243.851,0
opencl-kernel,gfx1010:xnack-,GPU,w16,secp256k1,256,ADD,20000,0.000082154,243445237.043,0
opencl-e2e,gfx1010:xnack-,GPU,w16,secp256k1,256,ADD,20000,0.000378319,52865438.956,0
opencl-kernel,gfx1010:xnack-,GPU,w32,secp256k1,256,ADD,20000,0.000068469,292102995.098,0
opencl-e2e,gfx1010:xnack-,GPU,w32,secp256k1,256,ADD,20000,0.000378239,52876620.316,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,secp256k1,256,ADD,20000,0.000071534,279587329.214,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,secp256k1,256,ADD,20000,0.000371566,53826238.136,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,secp256k1,256,ADD,20000,0.000070241,284733986.338,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,secp256k1,256,ADD,20000,0.000386240,51781275.904,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,secp256k1,256,ADD,20000,0.000049482,404187383.986,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,secp256k1,256,ADD,20000,0.000346477,57723889.283,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,secp256k1,256,ADD,20000,0.000062968,317621647.613,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,secp256k1,256,ADD,20000,0.000367426,54432729.291,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,secp256k1,256,SUBTRACT,20000,0.000194966,102581988.664,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,secp256k1,256,SUBTRACT,20000,0.000045034,444108895.683,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,secp256k1,256,SUBTRACT,20000,0.000052639,379946428.043,0
opencl-kernel,gfx1010:xnack-,GPU,w8,secp256k1,256,SUBTRACT,20000,0.000091532,218502818.553,0
opencl-e2e,gfx1010:xnack-,GPU,w8,secp256k1,256,SUBTRACT,20000,0.000381906,52368907.528,0
opencl-kernel,gfx1010:xnack-,GPU,w16,secp256k1,256,SUBTRACT,20000,0.000079910,250281566.551,0
opencl-e2e,gfx1010:xnack-,GPU,w16,secp256k1,256,SUBTRACT,20000,0.000371136,53888601.481,0
opencl-kernel,gfx1010:xnack-,GPU,w32,secp256k1,256,SUBTRACT,20000,0.000068268,292963027.305,0
opencl-e2e,gfx1010:xnack-,GPU,w32,secp256k1,256,SUBTRACT,20000,0.000396553,50434620.315,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,secp256k1,256,SUBTRACT,20000,0.000076203,262456858.878,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,secp256k1,256,SUBTRACT,20000,0.000361217,55368379.658,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,secp256k1,256,SUBTRACT,20000,0.000061495,325229691.307,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,secp256k1,256,SUBTRACT,20000,0.000374278,53436215.799,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,secp256k1,256,SUBTRACT,20000,0.000063579,314569275.039,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,secp256k1,256,SUBTRACT,20000,0.000345715,57851120.127,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,secp256k1,256,SUBTRACT,20000,0.000063438,315268451.732,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,secp256k1,256,SUBTRACT,20000,0.000357668,55917778.482,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,secp256k1,256,ADDMOD,20000,0.000659617,30320625.454,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,secp256k1,256,ADDMOD,20000,0.000482114,41483964.377,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,secp256k1,256,ADDMOD,20000,0.000312816,63935348.572,0
opencl-kernel,gfx1010:xnack-,GPU,w8,secp256k1,256,ADDMOD,20000,0.000089568,223294033.530,0
opencl-e2e,gfx1010:xnack-,GPU,w8,secp256k1,256,ADDMOD,20000,0.000373350,53569037.094,0
opencl-kernel,gfx1010:xnack-,GPU,w16,secp256k1,256,ADDMOD,20000,0.000088436,226152245.841,0
opencl-e2e,gfx1010:xnack-,GPU,w16,secp256k1,256,ADDMOD,20000,0.000338805,59031006.004,0
opencl-kernel,gfx1010:xnack-,GPU,w32,secp256k1,256,ADDMOD,20000,0.000076103,262801728.948,0
opencl-e2e,gfx1010:xnack-,GPU,w32,secp256k1,256,ADDMOD,20000,0.000375774,53223480.092,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,secp256k1,256,ADDMOD,20000,0.000070382,284163564.130,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,secp256k1,256,ADDMOD,20000,0.000344676,58025508.023,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,secp256k1,256,ADDMOD,20000,0.000059832,334269285.880,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,secp256k1,256,ADDMOD,20000,0.000358078,55853752.556,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,secp256k1,256,ADDMOD,20000,0.000059591,335621150.388,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,secp256k1,256,ADDMOD,20000,0.000339644,58885185.634,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,secp256k1,256,ADDMOD,20000,0.000065483,305422781.153,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,secp256k1,256,ADDMOD,20000,0.000366023,54641374.932,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,secp256k1,256,SUBTRACTMOD,20000,0.000534241,37436288.119,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,secp256k1,256,SUBTRACTMOD,20000,0.000112751,177382018.668,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,secp256k1,256,SUBTRACTMOD,20000,0.000331632,60307811.071,0
opencl-kernel,gfx1010:xnack-,GPU,w8,secp256k1,256,SUBTRACTMOD,20000,0.000089147,224348547.804,0
opencl-e2e,gfx1010:xnack-,GPU,w8,secp256k1,256,SUBTRACTMOD,20000,0.000359745,55594935.306,0
opencl-kernel,gfx1010:xnack-,GPU,w16,secp256k1,256,SUBTRACTMOD,20000,0.000095529,209360507.908,0
opencl-e2e,gfx1010:xnack-,GPU,w16,secp256k1,256,SUBTRACTMOD,20000,0.000362650,55149593.280,0
opencl-kernel,gfx1010:xnack-,GPU,w32,secp256k1,256,SUBTRACTMOD,20000,0.000074389,268856954.587,0
opencl-e2e,gfx1010:xnack-,GPU,w32,secp256k1,256,SUBTRACTMOD,20000,0.000362790,55128311.136,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,secp256k1,256,SUBTRACTMOD,20000,0.000066364,301368211.333,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,secp256k1,256,SUBTRACTMOD,20000,0.000347882,57490758.375,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,secp256k1,256,SUBTRACTMOD,20000,0.000062306,320996372.648,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,secp256k1,256,SUBTRACTMOD,20000,0.000352327,56765447.975,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,secp256k1,256,SUBTRACTMOD,20000,0.000042770,467617487.104,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,secp256k1,256,SUBTRACTMOD,20000,0.000339633,58887092.782,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,secp256k1,256,SUBTRACTMOD,20000,0.000062868,318126869.064,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,secp256k1,256,SUBTRACTMOD,20000,0.000358439,55797499.763,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,secp256k1,256,MULTIPLYOPERANDSCANNING,20000,0.000275787,72519734.442,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,20000,0.000055074,363147764.904,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,20000,0.000090419,221192448.483,0
opencl-kernel,gfx1010:xnack-,GPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,20000,0.001889564,10584452.287,0
opencl-e2e,gfx1010:xnack-,GPU,w8,secp256k1,256,MULTIPLYOPERANDSCANNING,20000,0.002226295,8983535.425,0
opencl-kernel,gfx1010:xnack-,GPU,w16,secp256k1,256,MULTIPLYOPERANDSCANNING,20000,0.000605756,33016594.145,0
opencl-e2e,gfx1010:xnack-,GPU,w16,secp256k1,256,MULTIPLYOPERANDSCANNING,20000,0.000945953,21142699.480,0
opencl-kernel,gfx1010:xnack-,GPU,w32,secp256k1,256,MULTIPLYOPERANDSCANNING,20000,0.000210996,94788526.864,0
opencl-e2e,gfx1010:xnack-,GPU,w32,secp256k1,256,MULTIPLYOPERANDSCANNING,20000,0.000532298,37572938.463,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,secp256k1,256,MULTIPLYOPERANDSCANNING,20000,0.000061836,323436186.113,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,secp256k1,256,MULTIPLYOPERANDSCANNING,20000,0.000372578,53680034.785,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,secp256k1,256,MULTIPLYOPERANDSCANNING,20000,0.000066364,301368212.366,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,secp256k1,256,MULTIPLYOPERANDSCANNING,20000,0.000396750,50409577.799,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,secp256k1,256,MULTIPLYOPERANDSCANNING,20000,0.000065161,306932059.207,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,secp256k1,256,MULTIPLYOPERANDSCANNING,20000,0.000364290,54901314.870,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,secp256k1,256,MULTIPLYOPERANDSCANNING,20000,0.000062476,320122927.003,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,secp256k1,256,MULTIPLYOPERANDSCANNING,20000,0.000400227,49971641.122,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,secp256k1,256,MULTIPLYPRODUCTSCANNING,20000,0.000276719,72275485.240,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,20000,0.000051356,389438429.998,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,20000,0.000091692,218121537.156,0
opencl-kernel,gfx1010:xnack-,GPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,20000,0.000343214,58272681.175,0
opencl-e2e,gfx1010:xnack-,GPU,w8,secp256k1,256,MULTIPLYPRODUCTSCANNING,20000,0.000672892,29722451.745,0
opencl-kernel,gfx1010:xnack-,GPU,w16,secp256k1,256,MULTIPLYPRODUCTSCANNING,20000,0.000245340,81519523.929,0
opencl-e2e,gfx1010:xnack-,GPU,w16,secp256k1,256,MULTIPLYPRODUCTSCANNING,20000,0.000564819,35409573.692,0
opencl-kernel,gfx1010:xnack-,GPU,w32,secp256k1,256,MULTIPLYPRODUCTSCANNING,20000,0.000114103,175280229.286,0
opencl-e2e,gfx1010:xnack-,GPU,w32,secp256k1,256,MULTIPLYPRODUCTSCANNING,20000,0.000438552,45604626.146,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,secp256k1,256,MULTIPLYPRODUCTSCANNING,20000,0.000086883,230194629.574,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,secp256k1,256,MULTIPLYPRODUCTSCANNING,20000,0.000384691,51989778.791,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,secp256k1,256,MULTIPLYPRODUCTSCANNING,20000,0.000112721,177429228.195,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,secp256k1,256,MULTIPLYPRODUCTSCANNING,20000,0.000437946,45667730.709,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,secp256k1,256,MULTIPLYPRODUCTSCANNING,20000,0.000077415,258347864.382,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,secp256k1,256,MULTIPLYPRODUCTSCANNING,20000,0.000398363,50205465.878,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,secp256k1,256,MULTIPLYPRODUCTSCANNING,20000,0.000076192,262494749.928,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,secp256k1,256,MULTIPLYPRODUCTSCANNING,20000,0.000390549,51209963.436,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,secp256k1,256,MONTGOMERYMULTIPLICATION,20000,0.002351039,8506877.172,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,20000,0.000333215,60021307.572,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,20000,0.000103554,193135948.356,0
opencl-kernel,gfx1010:xnack-,GPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,20000,0.000149621,133671075.592,0
opencl-e2e,gfx1010:xnack-,GPU,w8,secp256k1,256,MONTGOMERYMULTIPLICATION,20000,0.000421581,47440468.135,0
opencl-kernel,gfx1010:xnack-,GPU,w16,secp256k1,256,MONTGOMERYMULTIPLICATION,20000,0.000079749,250786843.160,0
opencl-e2e,gfx1010:xnack-,GPU,w16,secp256k1,256,MONTGOMERYMULTIPLICATION,20000,0.000358772,55745710.363,0
opencl-kernel,gfx1010:xnack-,GPU,w32,secp256k1,256,MONTGOMERYMULTIPLICATION,20000,0.000068849,290490784.413,0
opencl-e2e,gfx1010:xnack-,GPU,w32,secp256k1,256,MONTGOMERYMULTIPLICATION,20000,0.000384771,51978969.313,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,secp256k1,256,MONTGOMERYMULTIPLICATION,20000,0.000071274,280607234.084,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,secp256k1,256,MONTGOMERYMULTIPLICATION,20000,0.000356228,56143818.024,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,secp256k1,256,MONTGOMERYMULTIPLICATION,20000,0.000058739,340489283.121,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,secp256k1,256,MONTGOMERYMULTIPLICATION,20000,0.000362326,55198909.276,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,secp256k1,256,MONTGOMERYMULTIPLICATION,20000,0.000067636,295700514.345,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,secp256k1,256,MONTGOMERYMULTIPLICATION,20000,0.000343110,58290344.144,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,secp256k1,256,MONTGOMERYMULTIPLICATION,20000,0.000067156,297814044.406,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,secp256k1,256,MONTGOMERYMULTIPLICATION,20000,0.000372044,53757082.501,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,secp256k1,256,COMPARE,20000,0.000119384,167526636.793,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,secp256k1,256,COMPARE,20000,0.000021550,928074244.121,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,secp256k1,256,COMPARE,20000,0.000028754,695555403.079,0
opencl-kernel,gfx1010:xnack-,GPU,w8,secp256k1,256,COMPARE,20000,0.000067015,298440647.557,0
opencl-e2e,gfx1010:xnack-,GPU,w8,secp256k1,256,COMPARE,20000,0.000335359,59637582.408,0
opencl-kernel,gfx1010:xnack-,GPU,w16,secp256k1,256,COMPARE,20000,0.000053810,371678127.940,0
opencl-e2e,gfx1010:xnack-,GPU,w16,secp256k1,256,COMPARE,20000,0.000328135,60950523.409,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,secp256k1,256,COMPARE,20000,0.000048731,410416368.425,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,secp256k1,256,COMPARE,20000,0.000350277,57097668.408,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,secp256k1,256,COMPARE,20000,0.000071553,279513088.208,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,secp256k1,256,COMPARE,20000,0.000354872,56358348.980,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,secp256k1,256,COMPARE,20000,0.000071203,280887041.034,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,secp256k1,256,COMPARE,20000,0.000331008,60421500.368,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,secp256k1,256,COMPARE,20000,0.000062587,319555180.340,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,secp256k1,256,COMPARE,20000,0.000364120,54926947.147,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,secp256k1,256,REDUCE,2500,0.000051566,48481557.699,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,secp256k1,256,REDUCE,2500,0.000007103,351963957.408,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,secp256k1,256,REDUCE,2500,0.000042129,59341546.193,0
opencl-kernel,gfx1010:xnack-,GPU,w8,secp256k1,256,REDUCE,2500,0.000126628,19742868.880,0
opencl-e2e,gfx1010:xnack-,GPU,w8,secp256k1,256,REDUCE,2500,0.000227978,10965970.401,0
opencl-kernel,gfx1010:xnack-,GPU,w16,secp256k1,256,REDUCE,2500,0.000125075,19988007.209,0
opencl-e2e,gfx1010:xnack-,GPU,w16,secp256k1,256,REDUCE,2500,0.000192642,12977440.018,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,secp256k1,256,REDUCE,2500,0.000087504,28570122.502,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,secp256k1,256,REDUCE,2500,0.000171121,14609545.293,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,secp256k1,256,REDUCE,2500,0.000102030,24502597.250,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,secp256k1,256,REDUCE,2500,0.000185096,13506504.750,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,secp256k1,256,REDUCE,2500,0.000085770,29147720.591,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,secp256k1,256,REDUCE,2500,0.000163926,15250783.898,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,secp256k1,256,REDUCE,2500,0.000101269,24686725.461,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,secp256k1,256,REDUCE,2500,0.000179595,13920209.345,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,secp256k1,256,MODMUL,1250,0.000116638,10716919.014,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,secp256k1,256,MODMUL,1250,0.000013415,93179276.844,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,secp256k1,256,MODMUL,1250,0.000040155,31129373.704,0
opencl-kernel,gfx1010:xnack-,GPU,w8,secp256k1,256,MODMUL,2000,0.000258745,7729617.964,0
opencl-e2e,gfx1010:xnack-,GPU,w8,secp256k1,256,MODMUL,2000,0.000328736,6083909.278,0
opencl-kernel,gfx1010:xnack-,GPU,w16,secp256k1,256,MODMUL,2000,0.000159529,12536905.513,0
opencl-e2e,gfx1010:xnack-,GPU,w16,secp256k1,256,MODMUL,2000,0.000238266,8393979.838,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,secp256k1,256,MODMUL,2000,0.000135815,14725913.925,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,secp256k1,256,MODMUL,2000,0.000218209,9165524.796,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,secp256k1,256,MODMUL,2000,0.000130433,15333542.907,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,secp256k1,256,MODMUL,2000,0.000204091,9799550.206,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,secp256k1,256,MODMUL,2000,0.000134982,14816790.391,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,secp256k1,256,MODMUL,2000,0.000208219,9605271.357,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,secp256k1,256,MODMUL,2000,0.000131976,15154270.483,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,secp256k1,256,MODMUL,2000,0.000205744,9720818.097,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,secp256k1,256,MODEXP,312,0.002250190,138654.958,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,secp256k1,256,MODEXP,312,0.000292168,1067878.755,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,secp256k1,256,MODEXP,312,0.000528431,590427.132,0
opencl-kernel,gfx1010:xnack-,GPU,w8,secp256k1,256,MODEXP,2000,0.006434226,310837.698,0
opencl-e2e,gfx1010:xnack-,GPU,w8,secp256k1,256,MODEXP,2000,0.006554713,305123.962,0
opencl-kernel,gfx1010:xnack-,GPU,w16,secp256k1,256,MODEXP,2000,0.002299299,869830.327,0
opencl-e2e,gfx1010:xnack-,GPU,w16,secp256k1,256,MODEXP,2000,0.002363730,846120.327,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,secp256k1,256,MODEXP,2000,0.001065457,1877128.781,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,secp256k1,256,MODEXP,2000,0.001175293,1701703.320,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,secp256k1,256,MODEXP,2000,0.000804631,2485611.416,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,secp256k1,256,MODEXP,2000,0.000908334,2201833.246,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,secp256k1,256,MODEXP,2000,0.001017799,1965024.528,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,secp256k1,256,MODEXP,2000,0.001129226,1771124.647,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,secp256k1,256,MODEXP,2000,0.000827986,2415499.779,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,secp256k1,256,MODEXP,2000,0.000946668,2112673.080,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,secp256k1,256,EXPONENTIATION,312,0.001005386,310328.570,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,secp256k1,256,EXPONENTIATION,312,0.000097884,3187446.363,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,secp256k1,256,EXPONENTIATION,312,0.001428569,218400.371,0
opencl-kernel,gfx1010:xnack-,GPU,w8,secp256k1,256,EXPONENTIATION,2000,0.026999865,74074.444,0
opencl-e2e,gfx1010:xnack-,GPU,w8,secp256k1,256,EXPONENTIATION,2000,0.027162861,73629.946,0
opencl-kernel,gfx1010:xnack-,GPU,w16,secp256k1,256,EXPONENTIATION,2000,0.000343513,5822195.958,0
opencl-e2e,gfx1010:xnack-,GPU,w16,secp256k1,256,EXPONENTIATION,2000,0.000415189,4817083.304,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,secp256k1,256,EXPONENTIATION,2000,0.000247053,8095428.914,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,secp256k1,256,EXPONENTIATION,2000,0.000321652,6217900.092,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,secp256k1,256,EXPONENTIATION,2000,0.000198280,10086745.995,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,secp256k1,256,EXPONENTIATION,2000,0.000278759,7174656.252,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,secp256k1,256,EXPONENTIATION,2000,0.000232093,8617235.340,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,secp256k1,256,EXPONENTIATION,2000,0.000280163,7138701.391,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,secp256k1,256,EXPONENTIATION,2000,0.000183302,10910955.712,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,secp256k1,256,EXPONENTIATION,2000,0.000265275,7539345.963,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,secp256k1,256,DIVIDE,2500,0.000095429,26197487.138,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,secp256k1,256,DIVIDE,2500,0.000010530,237416905.426,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,secp256k1,256,DIVIDE,2500,0.000032871,76054881.132,0
opencl-kernel,gfx1010:xnack-,GPU,w8,secp256k1,256,DIVIDE,2500,0.000532819,4692024.871,0
opencl-e2e,gfx1010:xnack-,GPU,w8,secp256k1,256,DIVIDE,2500,0.000633879,3943970.379,0
opencl-kernel,gfx1010:xnack-,GPU,w16,secp256k1,256,DIVIDE,2500,0.000242044,10328700.564,0
opencl-e2e,gfx1010:xnack-,GPU,w16,secp256k1,256,DIVIDE,2500,0.000333344,7499760.009,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,secp256k1,256,DIVIDE,2500,0.000154189,16213867.385,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,secp256k1,256,DIVIDE,2500,0.000249468,10021325.381,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,secp256k1,256,DIVIDE,2500,0.000148798,16801301.061,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,secp256k1,256,DIVIDE,2500,0.000237763,10514672.180,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,secp256k1,256,DIVIDE,2500,0.000133779,18687536.915,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,secp256k1,256,DIVIDE,2500,0.000212627,11757678.959,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,secp256k1,256,DIVIDE,2500,0.000157194,15903914.891,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,secp256k1,256,DIVIDE,2500,0.000240308,10403315.735,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,secp256k1,256,ISQRT,625,0.000041919,14909706.791,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,secp256k1,256,ISQRT,625,0.000005050,123762378.361,0
opencl-kernel,gfx1010:xnack-,GPU,w8,secp256k1,256,ISQRT,2000,0.002472938,808754.607,0
opencl-e2e,gfx1010:xnack-,GPU,w8,secp256k1,256,ISQRT,2000,0.002577574,775923.407,0
opencl-kernel,gfx1010:xnack-,GPU,w16,secp256k1,256,ISQRT,2000,0.001220708,1638393.457,0
opencl-e2e,gfx1010:xnack-,GPU,w16,secp256k1,256,ISQRT,2000,0.001329282,1504571.641,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,secp256k1,256,ISQRT,2000,0.000673102,2971317.868,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,secp256k1,256,ISQRT,2000,0.000762179,2624055.504,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,secp256k1,256,ISQRT,2000,0.000690767,2895332.290,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,secp256k1,256,ISQRT,2000,0.000804810,2485058.586,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,secp256k1,256,ISQRT,2000,0.000680179,2940402.453,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,secp256k1,256,ISQRT,2000,0.000763685,2618880.821,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,secp256k1,256,ISQRT,2000,0.000647258,3089957.944,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,secp256k1,256,ISQRT,2000,0.000725685,2756016.729,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,secp256k1,256,MODMUL_R2,20000,0.001326017,15082762.891,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,secp256k1,256,MODMUL_R2,20000,0.000198883,100561636.715,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,secp256k1,256,MODMUL_R2,20000,0.000500419,39966508.067,0
opencl-kernel,gfx1010:xnack-,GPU,w8,secp256k1,256,MODMUL_R2,20000,0.000148568,134618491.209,0
opencl-e2e,gfx1010:xnack-,GPU,w8,secp256k1,256,MODMUL_R2,20000,0.000428985,46621676.750,0
opencl-kernel,gfx1010:xnack-,GPU,w16,secp256k1,256,MODMUL_R2,20000,0.000095970,208398457.689,0
opencl-e2e,gfx1010:xnack-,GPU,w16,secp256k1,256,MODMUL_R2,20000,0.000393237,50859913.993,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,secp256k1,256,MODMUL_R2,20000,0.000073338,272709918.486,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,secp256k1,256,MODMUL_R2,20000,0.000383318,52176000.100,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,secp256k1,256,MODMUL_R2,20000,0.000082564,242236325.855,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,secp256k1,256,MODMUL_R2,20000,0.000394366,50714311.083,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,secp256k1,256,MODMUL_R2,20000,0.000086652,230808290.803,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,secp256k1,256,MODMUL_R2,20000,0.000375410,53275085.930,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,secp256k1,256,MODMUL_R2,20000,0.000075351,265424479.025,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,secp256k1,256,MODMUL_R2,20000,0.000387463,51617831.937,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,rsa256(composite),256,ADD,20000,0.000245220,81559416.031,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,rsa256(composite),256,ADD,20000,0.000038021,526025091.078,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,rsa256(composite),256,ADD,20000,0.000046598,429202970.764,0
opencl-kernel,gfx1010:xnack-,GPU,w8,rsa256(composite),256,ADD,20000,0.000089218,224170010.609,0
opencl-e2e,gfx1010:xnack-,GPU,w8,rsa256(composite),256,ADD,20000,0.000408937,48907288.901,0
opencl-kernel,gfx1010:xnack-,GPU,w16,rsa256(composite),256,ADD,20000,0.000061816,323540829.818,0
opencl-e2e,gfx1010:xnack-,GPU,w16,rsa256(composite),256,ADD,20000,0.000369703,54097478.223,0
opencl-kernel,gfx1010:xnack-,GPU,w32,rsa256(composite),256,ADD,20000,0.000066425,301091457.045,0
opencl-e2e,gfx1010:xnack-,GPU,w32,rsa256(composite),256,ADD,20000,0.000384731,51984373.505,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,rsa256(composite),256,ADD,20000,0.000074359,268965424.922,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,rsa256(composite),256,ADD,20000,0.000360095,55540898.909,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,rsa256(composite),256,ADD,20000,0.000070431,283965868.974,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,rsa256(composite),256,ADD,20000,0.000397060,50370221.117,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,rsa256(composite),256,ADD,20000,0.000062457,320220310.343,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,rsa256(composite),256,ADD,20000,0.000368798,54230228.988,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,rsa256(composite),256,ADD,20000,0.000064250,311284047.065,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,rsa256(composite),256,ADD,20000,0.000349983,57145632.791,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACT,20000,0.000199764,100118139.375,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACT,20000,0.000039504,506277844.916,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACT,20000,0.000053480,373971578.023,0
opencl-kernel,gfx1010:xnack-,GPU,w8,rsa256(composite),256,SUBTRACT,20000,0.000088065,227104979.277,0
opencl-e2e,gfx1010:xnack-,GPU,w8,rsa256(composite),256,SUBTRACT,20000,0.000391083,51140039.333,0
opencl-kernel,gfx1010:xnack-,GPU,w16,rsa256(composite),256,SUBTRACT,20000,0.000066815,299333981.523,0
opencl-e2e,gfx1010:xnack-,GPU,w16,rsa256(composite),256,SUBTRACT,20000,0.000369613,54110650.874,0
opencl-kernel,gfx1010:xnack-,GPU,w32,rsa256(composite),256,SUBTRACT,20000,0.000070873,282194911.466,0
opencl-e2e,gfx1010:xnack-,GPU,w32,rsa256(composite),256,SUBTRACT,20000,0.000368120,54330109.736,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,rsa256(composite),256,SUBTRACT,20000,0.000076473,261530213.156,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,rsa256(composite),256,SUBTRACT,20000,0.000363140,55075177.619,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,rsa256(composite),256,SUBTRACT,20000,0.000077705,257383694.852,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,rsa256(composite),256,SUBTRACT,20000,0.000380559,52554268.859,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,rsa256(composite),256,SUBTRACT,20000,0.000065312,306222439.696,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,rsa256(composite),256,SUBTRACT,20000,0.000386671,51723558.286,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,rsa256(composite),256,SUBTRACT,20000,0.000064590,309645456.119,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,rsa256(composite),256,SUBTRACT,20000,0.000363278,55054255.916,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,rsa256(composite),256,ADDMOD,20000,0.000560972,35652403.327,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,rsa256(composite),256,ADDMOD,20000,0.000126507,158094018.507,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,rsa256(composite),256,ADDMOD,20000,0.000257824,77572297.359,0
opencl-kernel,gfx1010:xnack-,GPU,w8,rsa256(composite),256,ADDMOD,20000,0.000093856,213092397.042,0
opencl-e2e,gfx1010:xnack-,GPU,w8,rsa256(composite),256,ADDMOD,20000,0.000423535,47221599.161,0
opencl-kernel,gfx1010:xnack-,GPU,w16,rsa256(composite),256,ADDMOD,20000,0.000077465,258181114.038,0
opencl-e2e,gfx1010:xnack-,GPU,w16,rsa256(composite),256,ADDMOD,20000,0.000372909,53632387.507,0
opencl-kernel,gfx1010:xnack-,GPU,w32,rsa256(composite),256,ADDMOD,20000,0.000078798,253813548.554,0
opencl-e2e,gfx1010:xnack-,GPU,w32,rsa256(composite),256,ADDMOD,20000,0.000365665,54694871.005,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,rsa256(composite),256,ADDMOD,20000,0.000063900,312989045.450,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,rsa256(composite),256,ADDMOD,20000,0.000347582,57540378.939,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,rsa256(composite),256,ADDMOD,20000,0.000065181,306837881.820,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,rsa256(composite),256,ADDMOD,20000,0.000378455,52846441.464,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,rsa256(composite),256,ADDMOD,20000,0.000060072,332933812.777,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,rsa256(composite),256,ADDMOD,20000,0.000368788,54231699.471,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,rsa256(composite),256,ADDMOD,20000,0.000079278,252276799.049,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,rsa256(composite),256,ADDMOD,20000,0.000343932,58151029.880,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACTMOD,20000,0.000544892,36704521.263,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACTMOD,20000,0.000104776,190883408.310,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACTMOD,20000,0.000290215,68914425.524,0
opencl-kernel,gfx1010:xnack-,GPU,w8,rsa256(composite),256,SUBTRACTMOD,20000,0.000095329,209799746.071,0
opencl-e2e,gfx1010:xnack-,GPU,w8,rsa256(composite),256,SUBTRACTMOD,20000,0.000411432,48610706.013,0
opencl-kernel,gfx1010:xnack-,GPU,w16,rsa256(composite),256,SUBTRACTMOD,20000,0.000083586,239274519.446,0
opencl-e2e,gfx1010:xnack-,GPU,w16,rsa256(composite),256,SUBTRACTMOD,20000,0.000389049,51407406.268,0
opencl-kernel,gfx1010:xnack-,GPU,w32,rsa256(composite),256,SUBTRACTMOD,20000,0.000077104,259389915.067,0
opencl-e2e,gfx1010:xnack-,GPU,w32,rsa256(composite),256,SUBTRACTMOD,20000,0.000367179,54469346.010,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,rsa256(composite),256,SUBTRACTMOD,20000,0.000073087,273646475.570,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,rsa256(composite),256,SUBTRACTMOD,20000,0.000351459,56905641.901,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,rsa256(composite),256,SUBTRACTMOD,20000,0.000049813,401501617.616,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,rsa256(composite),256,SUBTRACTMOD,20000,0.000377734,52947312.092,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,rsa256(composite),256,SUBTRACTMOD,20000,0.000059571,335733828.570,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,rsa256(composite),256,SUBTRACTMOD,20000,0.000380119,52615102.168,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,rsa256(composite),256,SUBTRACTMOD,20000,0.000079609,251227875.314,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,rsa256(composite),256,SUBTRACTMOD,20000,0.000344904,57987150.042,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYOPERANDSCANNING,20000,0.000281027,71167539.055,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,20000,0.000044313,451334822.846,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,20000,0.000089348,223843846.632,0
opencl-kernel,gfx1010:xnack-,GPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,20000,0.001891137,10575648.407,0
opencl-e2e,gfx1010:xnack-,GPU,w8,rsa256(composite),256,MULTIPLYOPERANDSCANNING,20000,0.002282611,8761895.917,0
opencl-kernel,gfx1010:xnack-,GPU,w16,rsa256(composite),256,MULTIPLYOPERANDSCANNING,20000,0.000669976,29851815.581,0
opencl-e2e,gfx1010:xnack-,GPU,w16,rsa256(composite),256,MULTIPLYOPERANDSCANNING,20000,0.001009272,19816263.605,0
opencl-kernel,gfx1010:xnack-,GPU,w32,rsa256(composite),256,MULTIPLYOPERANDSCANNING,20000,0.000184696,108286048.481,0
opencl-e2e,gfx1010:xnack-,GPU,w32,rsa256(composite),256,MULTIPLYOPERANDSCANNING,20000,0.000500117,39990642.190,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,20000,0.000065643,304678335.982,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,20000,0.000376876,53067852.554,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,20000,0.000066705,299827599.665,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,20000,0.000393003,50890196.796,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,rsa256(composite),256,MULTIPLYOPERANDSCANNING,20000,0.000053129,376442242.880,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,rsa256(composite),256,MULTIPLYOPERANDSCANNING,20000,0.000401629,49797200.905,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,20000,0.000068598,291553690.162,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,rsa256(composite),256,MULTIPLYOPERANDSCANNING,20000,0.000381001,52493300.504,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,20000,0.000279414,71578374.733,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,20000,0.000051998,384630178.058,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,20000,0.000092544,216113416.493,0
opencl-kernel,gfx1010:xnack-,GPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,20000,0.000406202,49236586.720,0
opencl-e2e,gfx1010:xnack-,GPU,w8,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,20000,0.000772478,25890704.977,0
opencl-kernel,gfx1010:xnack-,GPU,w16,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,20000,0.000242905,82336716.060,0
opencl-e2e,gfx1010:xnack-,GPU,w16,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,20000,0.000559759,35729662.230,0
opencl-kernel,gfx1010:xnack-,GPU,w32,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,20000,0.000104025,192261475.740,0
opencl-e2e,gfx1010:xnack-,GPU,w32,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,20000,0.000419015,47730988.147,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,20000,0.000088336,226408259.726,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,20000,0.000392696,50929981.453,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,20000,0.000117760,169836956.242,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,20000,0.000448316,44611390.173,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,20000,0.000078226,255669469.303,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,20000,0.000398313,50211768.131,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,20000,0.000086802,230409438.329,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,20000,0.000387294,51640355.950,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,rsa256(composite),256,MONTGOMERYMULTIPLICATION,20000,0.002350429,8509084.937,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,20000,0.000337042,59339785.554,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,20000,0.000380744,52528733.212,0
opencl-kernel,gfx1010:xnack-,GPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,20000,0.000130364,153416587.368,0
opencl-e2e,gfx1010:xnack-,GPU,w8,rsa256(composite),256,MONTGOMERYMULTIPLICATION,20000,0.000465133,42998454.202,0
opencl-kernel,gfx1010:xnack-,GPU,w16,rsa256(composite),256,MONTGOMERYMULTIPLICATION,20000,0.000065984,303103782.723,0
opencl-e2e,gfx1010:xnack-,GPU,w16,rsa256(composite),256,MONTGOMERYMULTIPLICATION,20000,0.000349014,57304291.519,0
opencl-kernel,gfx1010:xnack-,GPU,w32,rsa256(composite),256,MONTGOMERYMULTIPLICATION,20000,0.000063980,312597687.644,0
opencl-e2e,gfx1010:xnack-,GPU,w32,rsa256(composite),256,MONTGOMERYMULTIPLICATION,20000,0.000355787,56213408.586,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,20000,0.000068618,291468710.350,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,20000,0.000360055,55547069.182,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,20000,0.000077083,259460580.851,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,20000,0.000357657,55919498.327,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,rsa256(composite),256,MONTGOMERYMULTIPLICATION,20000,0.000066715,299782658.583,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,rsa256(composite),256,MONTGOMERYMULTIPLICATION,20000,0.000398073,50242041.056,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,20000,0.000061144,327096689.629,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,rsa256(composite),256,MONTGOMERYMULTIPLICATION,20000,0.000322252,62063230.032,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,rsa256(composite),256,COMPARE,20000,0.000125355,159546886.957,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,rsa256(composite),256,COMPARE,20000,0.000021761,919075411.368,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,rsa256(composite),256,COMPARE,20000,0.000027662,723013520.889,0
opencl-kernel,gfx1010:xnack-,GPU,w8,rsa256(composite),256,COMPARE,20000,0.000066124,302462041.014,0
opencl-e2e,gfx1010:xnack-,GPU,w8,rsa256(composite),256,COMPARE,20000,0.000344085,58125172.570,0
opencl-kernel,gfx1010:xnack-,GPU,w16,rsa256(composite),256,COMPARE,20000,0.000065533,305189751.741,0
opencl-e2e,gfx1010:xnack-,GPU,w16,rsa256(composite),256,COMPARE,20000,0.000338204,59135906.149,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,rsa256(composite),256,COMPARE,20000,0.000048922,408814031.288,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,rsa256(composite),256,COMPARE,20000,0.000346329,57748557.013,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,rsa256(composite),256,COMPARE,20000,0.000055253,361971295.495,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,rsa256(composite),256,COMPARE,20000,0.000350493,57062480.556,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,rsa256(composite),256,COMPARE,20000,0.000062537,319810670.530,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,rsa256(composite),256,COMPARE,20000,0.000359822,55583038.270,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,rsa256(composite),256,COMPARE,20000,0.000070952,281880708.135,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,rsa256(composite),256,COMPARE,20000,0.000348671,57360663.743,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,rsa256(composite),256,REDUCE,2500,0.000077356,32318113.648,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,rsa256(composite),256,REDUCE,2500,0.000009638,259389913.537,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,rsa256(composite),256,REDUCE,2500,0.000033183,75339782.356,0
opencl-kernel,gfx1010:xnack-,GPU,w8,rsa256(composite),256,REDUCE,2500,0.000132027,18935520.756,0
opencl-e2e,gfx1010:xnack-,GPU,w8,rsa256(composite),256,REDUCE,2500,0.000214352,11663058.897,0
opencl-kernel,gfx1010:xnack-,GPU,w16,rsa256(composite),256,REDUCE,2500,0.000119293,20956803.814,0
opencl-e2e,gfx1010:xnack-,GPU,w16,rsa256(composite),256,REDUCE,2500,0.000209002,11961608.022,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,rsa256(composite),256,REDUCE,2500,0.000102833,24311261.942,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,rsa256(composite),256,REDUCE,2500,0.000197851,12635771.367,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,rsa256(composite),256,REDUCE,2500,0.000104124,24009834.477,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,rsa256(composite),256,REDUCE,2500,0.000184504,13549841.729,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,rsa256(composite),256,REDUCE,2500,0.000101159,24713569.783,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,rsa256(composite),256,REDUCE,2500,0.000182480,13700131.510,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,rsa256(composite),256,REDUCE,2500,0.000096831,25818178.046,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,rsa256(composite),256,REDUCE,2500,0.000181329,13787094.174,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,rsa256(composite),256,MODMUL,1250,0.000117210,10664619.061,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,rsa256(composite),256,MODMUL,1250,0.000012573,99419390.370,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,rsa256(composite),256,MODMUL,1250,0.000036819,33949862.841,0
opencl-kernel,gfx1010:xnack-,GPU,w8,rsa256(composite),256,MODMUL,2000,0.000243687,8207249.465,0
opencl-e2e,gfx1010:xnack-,GPU,w8,rsa256(composite),256,MODMUL,2000,0.000326142,6132298.200,0
opencl-kernel,gfx1010:xnack-,GPU,w16,rsa256(composite),256,MODMUL,2000,0.000176270,11346230.223,0
opencl-e2e,gfx1010:xnack-,GPU,w16,rsa256(composite),256,MODMUL,2000,0.000241503,8281470.620,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,rsa256(composite),256,MODMUL,2000,0.000131657,15190988.693,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,rsa256(composite),256,MODMUL,2000,0.000218199,9165944.855,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,rsa256(composite),256,MODMUL,2000,0.000124051,16122401.263,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,rsa256(composite),256,MODMUL,2000,0.000189934,10529973.570,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,rsa256(composite),256,MODMUL,2000,0.000133349,14998237.712,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,rsa256(composite),256,MODMUL,2000,0.000204562,9776986.933,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,rsa256(composite),256,MODMUL,2000,0.000138518,14438556.737,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,rsa256(composite),256,MODMUL,2000,0.000213639,9361586.588,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,rsa256(composite),256,MODEXP,312,0.002140544,145757.340,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,rsa256(composite),256,MODEXP,312,0.000246532,1265555.790,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,rsa256(composite),256,MODEXP,312,0.000603643,516861.788,0
opencl-kernel,gfx1010:xnack-,GPU,w8,rsa256(composite),256,MODEXP,2000,0.006423968,311334.054,0
opencl-e2e,gfx1010:xnack-,GPU,w8,rsa256(composite),256,MODEXP,2000,0.006541288,305750.183,0
opencl-kernel,gfx1010:xnack-,GPU,w16,rsa256(composite),256,MODEXP,2000,0.002307395,866778.337,0
opencl-e2e,gfx1010:xnack-,GPU,w16,rsa256(composite),256,MODEXP,2000,0.002410378,829745.376,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,rsa256(composite),256,MODEXP,2000,0.001030622,1940575.691,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,rsa256(composite),256,MODEXP,2000,0.001114629,1794319.007,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,rsa256(composite),256,MODEXP,2000,0.000799070,2502909.632,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,rsa256(composite),256,MODEXP,2000,0.000857329,2332826.721,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,rsa256(composite),256,MODEXP,2000,0.001036394,1929768.023,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,rsa256(composite),256,MODEXP,2000,0.001125741,1776607.586,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,rsa256(composite),256,MODEXP,2000,0.000833746,2398812.107,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,rsa256(composite),256,MODEXP,2000,0.000943282,2120256.720,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,rsa256(composite),256,EXPONENTIATION,312,0.001006137,310096.935,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,rsa256(composite),256,EXPONENTIATION,312,0.000102402,3046815.494,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,rsa256(composite),256,EXPONENTIATION,312,0.001499983,208002.357,0
opencl-kernel,gfx1010:xnack-,GPU,w8,rsa256(composite),256,EXPONENTIATION,2000,0.026849522,74489.222,0
opencl-e2e,gfx1010:xnack-,GPU,w8,rsa256(composite),256,EXPONENTIATION,2000,0.027025112,74005.244,0
opencl-kernel,gfx1010:xnack-,GPU,w16,rsa256(composite),256,EXPONENTIATION,2000,0.000309800,6455777.922,0
opencl-e2e,gfx1010:xnack-,GPU,w16,rsa256(composite),256,EXPONENTIATION,2000,0.000386485,5174845.078,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,rsa256(composite),256,EXPONENTIATION,2000,0.000246963,8098379.106,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,rsa256(composite),256,EXPONENTIATION,2000,0.000300854,6647742.759,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,rsa256(composite),256,EXPONENTIATION,2000,0.000197669,10117924.407,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,rsa256(composite),256,EXPONENTIATION,2000,0.000275974,7247059.505,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,rsa256(composite),256,EXPONENTIATION,2000,0.000231031,8656846.912,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,rsa256(composite),256,EXPONENTIATION,2000,0.000301803,6626839.358,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,rsa256(composite),256,EXPONENTIATION,2000,0.000210844,9485686.106,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,rsa256(composite),256,EXPONENTIATION,2000,0.000286014,6992664.699,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,rsa256(composite),256,DIVIDE,2500,0.000097433,25658657.757,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,rsa256(composite),256,DIVIDE,2500,0.000011281,221611560.115,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,rsa256(composite),256,DIVIDE,2500,0.000032481,76968073.765,0
opencl-kernel,gfx1010:xnack-,GPU,w8,rsa256(composite),256,DIVIDE,2500,0.000542537,4607980.654,0
opencl-e2e,gfx1010:xnack-,GPU,w8,rsa256(composite),256,DIVIDE,2500,0.000652453,3831693.624,0
opencl-kernel,gfx1010:xnack-,GPU,w16,rsa256(composite),256,DIVIDE,2500,0.000253776,9851207.368,0
opencl-e2e,gfx1010:xnack-,GPU,w16,rsa256(composite),256,DIVIDE,2500,0.000350246,7137840.262,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,rsa256(composite),256,DIVIDE,2500,0.000156965,15927117.524,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,rsa256(composite),256,DIVIDE,2500,0.000256531,9745410.882,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,rsa256(composite),256,DIVIDE,2500,0.000155670,16059613.294,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,rsa256(composite),256,DIVIDE,2500,0.000253803,9850159.381,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,rsa256(composite),256,DIVIDE,2500,0.000155501,16077067.041,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,rsa256(composite),256,DIVIDE,2500,0.000240468,10396393.703,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,rsa256(composite),256,DIVIDE,2500,0.000159849,15639760.040,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,rsa256(composite),256,DIVIDE,2500,0.000253414,9865279.745,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,rsa256(composite),256,ISQRT,625,0.000041458,15075498.116,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,rsa256(composite),256,ISQRT,625,0.000005060,123517788.483,0
opencl-kernel,gfx1010:xnack-,GPU,w8,rsa256(composite),256,ISQRT,2000,0.002747593,727909.847,0
opencl-e2e,gfx1010:xnack-,GPU,w8,rsa256(composite),256,ISQRT,2000,0.002854152,700733.528,0
opencl-kernel,gfx1010:xnack-,GPU,w16,rsa256(composite),256,ISQRT,2000,0.001313732,1522380.516,0
opencl-e2e,gfx1010:xnack-,GPU,w16,rsa256(composite),256,ISQRT,2000,0.001426644,1401891.432,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,rsa256(composite),256,ISQRT,2000,0.000706634,2830319.515,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,rsa256(composite),256,ISQRT,2000,0.000812263,2462256.683,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,rsa256(composite),256,ISQRT,2000,0.000700656,2854467.813,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,rsa256(composite),256,ISQRT,2000,0.000768092,2603854.747,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,rsa256(composite),256,ISQRT,2000,0.000699124,2860722.847,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,rsa256(composite),256,ISQRT,2000,0.000801796,2494400.072,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,rsa256(composite),256,ISQRT,2000,0.000702652,2846359.221,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,rsa256(composite),256,ISQRT,2000,0.000792240,2524487.529,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,rsa256(composite),256,MODMUL_R2,20000,0.001325756,15085732.217,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,rsa256(composite),256,MODMUL_R2,20000,0.000703609,28424878.023,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,rsa256(composite),256,MODMUL_R2,20000,0.000478146,41828228.200,0
opencl-kernel,gfx1010:xnack-,GPU,w8,rsa256(composite),256,MODMUL_R2,20000,0.000125516,159342235.220,0
opencl-e2e,gfx1010:xnack-,GPU,w8,rsa256(composite),256,MODMUL_R2,20000,0.000452559,44193132.830,0
opencl-kernel,gfx1010:xnack-,GPU,w16,rsa256(composite),256,MODMUL_R2,20000,0.000074219,269472776.519,0
opencl-e2e,gfx1010:xnack-,GPU,w16,rsa256(composite),256,MODMUL_R2,20000,0.000366947,54503783.929,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,rsa256(composite),256,MODMUL_R2,20000,0.000080291,249093920.456,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,rsa256(composite),256,MODMUL_R2,20000,0.000383990,52084689.699,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,rsa256(composite),256,MODMUL_R2,20000,0.000070782,282557712.548,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,rsa256(composite),256,MODMUL_R2,20000,0.000389647,51328510.132,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,rsa256(composite),256,MODMUL_R2,20000,0.000065202,306739057.744,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,rsa256(composite),256,MODMUL_R2,20000,0.000389596,51335229.331,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,rsa256(composite),256,MODMUL_R2,20000,0.000070101,285302635.209,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,rsa256(composite),256,MODMUL_R2,20000,0.000359681,55604827.582,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,brainpoolP512r1,512,ADD,10000,0.000131306,76157982.077,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,brainpoolP512r1,512,ADD,10000,0.000021270,470145746.644,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,brainpoolP512r1,512,ADD,10000,0.000030267,330392836.467,0
opencl-kernel,gfx1010:xnack-,GPU,w8,brainpoolP512r1,512,ADD,10000,0.000124413,80377452.456,0
opencl-e2e,gfx1010:xnack-,GPU,w8,brainpoolP512r1,512,ADD,10000,0.000428223,23352318.770,0
opencl-kernel,gfx1010:xnack-,GPU,w16,brainpoolP512r1,512,ADD,10000,0.000103374,96736123.211,0
opencl-e2e,gfx1010:xnack-,GPU,w16,brainpoolP512r1,512,ADD,10000,0.000391283,25556949.829,0
opencl-kernel,gfx1010:xnack-,GPU,w32,brainpoolP512r1,512,ADD,10000,0.000077666,128756469.761,0
opencl-e2e,gfx1010:xnack-,GPU,w32,brainpoolP512r1,512,ADD,10000,0.000383028,26107751.915,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,brainpoolP512r1,512,ADD,10000,0.000081453,122770186.306,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,brainpoolP512r1,512,ADD,10000,0.000343724,29093109.588,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,brainpoolP512r1,512,ADD,10000,0.000077956,128277489.799,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,brainpoolP512r1,512,ADD,10000,0.000379618,26342270.363,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,brainpoolP512r1,512,ADD,10000,0.000063819,156693147.972,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,brainpoolP512r1,512,ADD,10000,0.000361855,27635378.817,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,brainpoolP512r1,512,ADD,10000,0.000063248,158107767.239,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,brainpoolP512r1,512,ADD,10000,0.000362226,27607074.029,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACT,10000,0.000110988,90099830.645,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACT,10000,0.000431309,23185233.790,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACT,10000,0.000344456,29031284.112,0
opencl-kernel,gfx1010:xnack-,GPU,w8,brainpoolP512r1,512,SUBTRACT,10000,0.000111769,89470246.669,0
opencl-e2e,gfx1010:xnack-,GPU,w8,brainpoolP512r1,512,SUBTRACT,10000,0.000394390,25355612.465,0
opencl-kernel,gfx1010:xnack-,GPU,w16,brainpoolP512r1,512,SUBTRACT,10000,0.000080360,124440019.781,0
opencl-e2e,gfx1010:xnack-,GPU,w16,brainpoolP512r1,512,SUBTRACT,10000,0.000380814,26259538.779,0
opencl-kernel,gfx1010:xnack-,GPU,w32,brainpoolP512r1,512,SUBTRACT,10000,0.000088927,112451786.323,0
opencl-e2e,gfx1010:xnack-,GPU,w32,brainpoolP512r1,512,SUBTRACT,10000,0.000378840,26396367.864,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,brainpoolP512r1,512,SUBTRACT,10000,0.000072566,137805583.936,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,brainpoolP512r1,512,SUBTRACT,10000,0.000342562,29191795.936,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,brainpoolP512r1,512,SUBTRACT,10000,0.000079548,125710262.937,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,brainpoolP512r1,512,SUBTRACT,10000,0.000402380,24852129.815,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,brainpoolP512r1,512,SUBTRACT,10000,0.000065783,152014958.083,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,brainpoolP512r1,512,SUBTRACT,10000,0.000362516,27584989.335,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,brainpoolP512r1,512,SUBTRACT,10000,0.000056014,178526796.345,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,brainpoolP512r1,512,SUBTRACT,10000,0.000355124,28159178.224,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,brainpoolP512r1,512,ADDMOD,10000,0.000320170,31233407.252,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,brainpoolP512r1,512,ADDMOD,10000,0.000318017,31444859.873,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,brainpoolP512r1,512,ADDMOD,10000,0.000188604,53021144.844,0
opencl-kernel,gfx1010:xnack-,GPU,w8,brainpoolP512r1,512,ADDMOD,10000,0.000118392,84465166.572,0
opencl-e2e,gfx1010:xnack-,GPU,w8,brainpoolP512r1,512,ADDMOD,10000,0.000422713,23656712.709,0
opencl-kernel,gfx1010:xnack-,GPU,w16,brainpoolP512r1,512,ADDMOD,10000,0.000097884,102161742.398,0
opencl-e2e,gfx1010:xnack-,GPU,w16,brainpoolP512r1,512,ADDMOD,10000,0.000386144,25897074.666,0
opencl-kernel,gfx1010:xnack-,GPU,w32,brainpoolP512r1,512,ADDMOD,10000,0.000081152,123225552.212,0
opencl-e2e,gfx1010:xnack-,GPU,w32,brainpoolP512r1,512,ADDMOD,10000,0.000374863,26676412.454,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,brainpoolP512r1,512,ADDMOD,10000,0.000063258,158082772.253,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,brainpoolP512r1,512,ADDMOD,10000,0.000324368,30829181.686,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,brainpoolP512r1,512,ADDMOD,10000,0.000064981,153891136.514,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,brainpoolP512r1,512,ADDMOD,10000,0.000362085,27617824.561,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,brainpoolP512r1,512,ADDMOD,10000,0.000058329,171441307.701,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,brainpoolP512r1,512,ADDMOD,10000,0.000368588,27130563.113,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,brainpoolP512r1,512,ADDMOD,10000,0.000064260,155617803.119,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,brainpoolP512r1,512,ADDMOD,10000,0.000345505,28943141.196,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACTMOD,10000,0.000296446,33732956.425,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACTMOD,10000,0.000059522,168005107.313,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACTMOD,10000,0.000172603,57936420.569,0
opencl-kernel,gfx1010:xnack-,GPU,w8,brainpoolP512r1,512,SUBTRACTMOD,10000,0.000115567,86529891.776,0
opencl-e2e,gfx1010:xnack-,GPU,w8,brainpoolP512r1,512,SUBTRACTMOD,10000,0.000418445,23898003.320,0
opencl-kernel,gfx1010:xnack-,GPU,w16,brainpoolP512r1,512,SUBTRACTMOD,10000,0.000093055,107463328.163,0
opencl-e2e,gfx1010:xnack-,GPU,w16,brainpoolP512r1,512,SUBTRACTMOD,10000,0.000387606,25799394.230,0
opencl-kernel,gfx1010:xnack-,GPU,w32,brainpoolP512r1,512,SUBTRACTMOD,10000,0.000080431,124330171.180,0
opencl-e2e,gfx1010:xnack-,GPU,w32,brainpoolP512r1,512,SUBTRACTMOD,10000,0.000382417,26149465.107,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,brainpoolP512r1,512,SUBTRACTMOD,10000,0.000066344,150729530.792,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,brainpoolP512r1,512,SUBTRACTMOD,10000,0.000336942,29678698.402,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,brainpoolP512r1,512,SUBTRACTMOD,10000,0.000050424,198318262.286,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,brainpoolP512r1,512,SUBTRACTMOD,10000,0.000356004,28089572.030,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,brainpoolP512r1,512,SUBTRACTMOD,10000,0.000045725,218698743.265,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,brainpoolP512r1,512,SUBTRACTMOD,10000,0.000327110,30570756.020,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,brainpoolP512r1,512,SUBTRACTMOD,10000,0.000063559,157334129.479,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,brainpoolP512r1,512,SUBTRACTMOD,10000,0.000334294,29913788.469,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,10000,0.000370405,26997475.743,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,10000,0.000054462,183614262.840,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,10000,0.000059342,168514711.374,0
opencl-kernel,gfx1010:xnack-,GPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,10000,0.007392343,1352751.083,0
opencl-e2e,gfx1010:xnack-,GPU,w8,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,10000,0.007650427,1307116.583,0
opencl-kernel,gfx1010:xnack-,GPU,w16,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,10000,0.001550847,6448089.334,0
opencl-e2e,gfx1010:xnack-,GPU,w16,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,10000,0.001888539,5295098.486,0
opencl-kernel,gfx1010:xnack-,GPU,w32,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,10000,0.000382096,26171433.352,0
opencl-e2e,gfx1010:xnack-,GPU,w32,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,10000,0.000732253,13656482.117,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,10000,0.000086201,116007935.036,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,10000,0.000393397,25419614.288,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,10000,0.000081432,122801846.766,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,10000,0.000409323,24430584.161,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,10000,0.000064120,155957579.700,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,10000,0.000418180,23913147.449,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,10000,0.000062546,159882326.075,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,10000,0.000376152,26584997.542,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,10000,0.000372729,26829143.963,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,10000,0.000053690,186254423.291,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,10000,0.000060223,166049515.990,0
opencl-kernel,gfx1010:xnack-,GPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,10000,0.002747202,3640067.239,0
opencl-e2e,gfx1010:xnack-,GPU,w8,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,10000,0.003123708,3201323.555,0
opencl-kernel,gfx1010:xnack-,GPU,w16,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,10000,0.000776496,12878366.408,0
opencl-e2e,gfx1010:xnack-,GPU,w16,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,10000,0.001133955,8818692.099,0
opencl-kernel,gfx1010:xnack-,GPU,w32,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,10000,0.000232406,43028148.992,0
opencl-e2e,gfx1010:xnack-,GPU,w32,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,10000,0.000551975,18116762.534,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,10000,0.000237795,42053028.874,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,10000,0.000553839,18055788.774,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,10000,0.000235559,42452209.430,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,10000,0.000560946,17827027.915,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,10000,0.000107611,92927302.795,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,10000,0.000412740,24228327.774,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,10000,0.000095538,104670393.326,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,10000,0.000378217,26439848.044,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,10000,0.002878228,3474359.919,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,10000,0.000345117,28975680.716,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,10000,0.000113994,87723915.291,0
opencl-kernel,gfx1010:xnack-,GPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,10000,0.000217558,45964754.220,0
opencl-e2e,gfx1010:xnack-,GPU,w8,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,10000,0.000534693,18702320.772,0
opencl-kernel,gfx1010:xnack-,GPU,w16,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,10000,0.000102783,97292353.912,0
opencl-e2e,gfx1010:xnack-,GPU,w16,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,10000,0.000411872,24279387.766,0
opencl-kernel,gfx1010:xnack-,GPU,w32,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,10000,0.000079970,125046892.760,0
opencl-e2e,gfx1010:xnack-,GPU,w32,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,10000,0.000394660,25338265.850,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,10000,0.000076543,130645519.550,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,10000,0.000358352,27905523.067,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,10000,0.000076342,130989494.935,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,10000,0.000389686,25661686.600,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,10000,0.000085149,117441190.814,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,10000,0.000387022,25838324.449,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,10000,0.000078066,128096739.006,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,10000,0.000387073,25834920.046,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,brainpoolP512r1,512,COMPARE,10000,0.000061465,162694216.404,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,brainpoolP512r1,512,COMPARE,10000,0.000012113,825559314.717,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,brainpoolP512r1,512,COMPARE,10000,0.000018044,554200844.106,0
opencl-kernel,gfx1010:xnack-,GPU,w8,brainpoolP512r1,512,COMPARE,10000,0.000067476,148200841.722,0
opencl-e2e,gfx1010:xnack-,GPU,w8,brainpoolP512r1,512,COMPARE,10000,0.000380513,26280311.055,0
opencl-kernel,gfx1010:xnack-,GPU,w16,brainpoolP512r1,512,COMPARE,10000,0.000066955,149354043.837,0
opencl-e2e,gfx1010:xnack-,GPU,w16,brainpoolP512r1,512,COMPARE,10000,0.000378680,26407520.859,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,brainpoolP512r1,512,COMPARE,10000,0.000059091,169230509.129,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,brainpoolP512r1,512,COMPARE,10000,0.000341039,29322159.638,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,brainpoolP512r1,512,COMPARE,10000,0.000065171,153442481.644,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,brainpoolP512r1,512,COMPARE,10000,0.000361604,27654561.354,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,brainpoolP512r1,512,COMPARE,10000,0.000065523,152618165.433,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,brainpoolP512r1,512,COMPARE,10000,0.000342419,29203986.934,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,brainpoolP512r1,512,COMPARE,10000,0.000065412,152877148.128,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,brainpoolP512r1,512,COMPARE,10000,0.000343291,29129805.342,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,brainpoolP512r1,512,REDUCE,1250,0.000041418,30180114.909,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,brainpoolP512r1,512,REDUCE,1250,0.000005410,231053605.651,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,brainpoolP512r1,512,REDUCE,1250,0.000026219,47675349.927,0
opencl-kernel,gfx1010:xnack-,GPU,w8,brainpoolP512r1,512,REDUCE,2000,0.000259948,7693846.461,0
opencl-e2e,gfx1010:xnack-,GPU,w8,brainpoolP512r1,512,REDUCE,2000,0.000361137,5538064.501,0
opencl-kernel,gfx1010:xnack-,GPU,w16,brainpoolP512r1,512,REDUCE,2000,0.000158006,12657747.175,0
opencl-e2e,gfx1010:xnack-,GPU,w16,brainpoolP512r1,512,REDUCE,2000,0.000248566,8046152.733,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,brainpoolP512r1,512,REDUCE,2000,0.000147396,13568889.262,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,brainpoolP512r1,512,REDUCE,2000,0.000238688,8379139.292,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,brainpoolP512r1,512,REDUCE,2000,0.000143458,13941362.593,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,brainpoolP512r1,512,REDUCE,2000,0.000224910,8892445.875,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,brainpoolP512r1,512,REDUCE,2000,0.000145211,13773061.268,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,brainpoolP512r1,512,REDUCE,2000,0.000214400,9328358.209,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,brainpoolP512r1,512,REDUCE,2000,0.000138639,14425955.185,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,brainpoolP512r1,512,REDUCE,2000,0.000240549,8314314.341,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL,625,0.000103023,6066606.489,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL,625,0.000011231,55649541.297,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL,625,0.000034946,17884736.436,0
opencl-kernel,gfx1010:xnack-,GPU,w8,brainpoolP512r1,512,MODMUL,2000,0.000717766,2786423.431,0
opencl-e2e,gfx1010:xnack-,GPU,w8,brainpoolP512r1,512,MODMUL,2000,0.000846818,2361782.579,0
opencl-kernel,gfx1010:xnack-,GPU,w16,brainpoolP512r1,512,MODMUL,2000,0.000329687,6066359.912,0
opencl-e2e,gfx1010:xnack-,GPU,w16,brainpoolP512r1,512,MODMUL,2000,0.000429044,4661526.557,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,brainpoolP512r1,512,MODMUL,2000,0.000259566,7705169.396,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,brainpoolP512r1,512,MODMUL,2000,0.000353652,5655276.938,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,brainpoolP512r1,512,MODMUL,2000,0.000251489,7952634.119,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,brainpoolP512r1,512,MODMUL,2000,0.000341156,5862420.709,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,brainpoolP512r1,512,MODMUL,2000,0.000250156,7995011.121,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,brainpoolP512r1,512,MODMUL,2000,0.000344252,5809697.553,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,brainpoolP512r1,512,MODMUL,2000,0.000236331,8462706.966,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,brainpoolP512r1,512,MODMUL,2000,0.000326920,6117704.635,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,brainpoolP512r1,512,MODEXP,156,0.006019178,25917.160,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,brainpoolP512r1,512,MODEXP,156,0.000747311,208748.433,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,brainpoolP512r1,512,MODEXP,156,0.000714821,218236.454,0
opencl-kernel,gfx1010:xnack-,GPU,w8,brainpoolP512r1,512,MODEXP,2000,0.041737350,47918.711,0
opencl-e2e,gfx1010:xnack-,GPU,w8,brainpoolP512r1,512,MODEXP,2000,0.041856314,47782.516,0
opencl-kernel,gfx1010:xnack-,GPU,w16,brainpoolP512r1,512,MODEXP,2000,0.012858199,155542.779,0
opencl-e2e,gfx1010:xnack-,GPU,w16,brainpoolP512r1,512,MODEXP,2000,0.013017147,153643.498,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,brainpoolP512r1,512,MODEXP,2000,0.007931907,252146.174,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,brainpoolP512r1,512,MODEXP,2000,0.007988563,250357.918,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,brainpoolP512r1,512,MODEXP,2000,0.005133052,389631.743,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,brainpoolP512r1,512,MODEXP,2000,0.005260699,380177.615,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,brainpoolP512r1,512,MODEXP,2000,0.008170835,244773.025,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,brainpoolP512r1,512,MODEXP,2000,0.008288775,241290.179,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,brainpoolP512r1,512,MODEXP,2000,0.005107688,391566.595,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,brainpoolP512r1,512,MODEXP,2000,0.005278876,378868.532,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,brainpoolP512r1,512,EXPONENTIATION,156,0.001361463,114582.622,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,brainpoolP512r1,512,EXPONENTIATION,156,0.000266860,584576.182,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,brainpoolP512r1,512,EXPONENTIATION,156,0.002292179,68057.512,0
opencl-kernel,gfx1010:xnack-,GPU,w8,brainpoolP512r1,512,EXPONENTIATION,2000,0.299823038,6670.601,0
opencl-e2e,gfx1010:xnack-,GPU,w8,brainpoolP512r1,512,EXPONENTIATION,2000,0.300015930,6666.313,0
opencl-kernel,gfx1010:xnack-,GPU,w16,brainpoolP512r1,512,EXPONENTIATION,2000,0.074594716,26811.551,0
opencl-e2e,gfx1010:xnack-,GPU,w16,brainpoolP512r1,512,EXPONENTIATION,2000,0.074723448,26765.360,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,brainpoolP512r1,512,EXPONENTIATION,2000,0.004707294,424872.549,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,brainpoolP512r1,512,EXPONENTIATION,2000,0.004834342,413706.767,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,brainpoolP512r1,512,EXPONENTIATION,2000,0.005245020,381314.085,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,brainpoolP512r1,512,EXPONENTIATION,2000,0.005369552,372470.552,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,brainpoolP512r1,512,EXPONENTIATION,2000,0.004717713,423934.224,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,brainpoolP512r1,512,EXPONENTIATION,2000,0.004836654,413509.009,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,brainpoolP512r1,512,EXPONENTIATION,2000,0.005239794,381694.395,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,brainpoolP512r1,512,EXPONENTIATION,2000,0.005326356,375491.236,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,brainpoolP512r1,512,DIVIDE,1250,0.000053841,23216507.837,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,brainpoolP512r1,512,DIVIDE,1250,0.000006272,199298467.348,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,brainpoolP512r1,512,DIVIDE,1250,0.000023454,53295813.273,0
opencl-kernel,gfx1010:xnack-,GPU,w8,brainpoolP512r1,512,DIVIDE,2000,0.003759310,532012.524,0
opencl-e2e,gfx1010:xnack-,GPU,w8,brainpoolP512r1,512,DIVIDE,2000,0.003905314,512122.713,0
opencl-kernel,gfx1010:xnack-,GPU,w16,brainpoolP512r1,512,DIVIDE,2000,0.000959688,2084010.637,0
opencl-e2e,gfx1010:xnack-,GPU,w16,brainpoolP512r1,512,DIVIDE,2000,0.001107946,1805142.128,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,brainpoolP512r1,512,DIVIDE,2000,0.000407954,4902513.517,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,brainpoolP512r1,512,DIVIDE,2000,0.000535324,3736055.173,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,brainpoolP512r1,512,DIVIDE,2000,0.000406799,4916432.933,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,brainpoolP512r1,512,DIVIDE,2000,0.000523637,3819439.801,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,brainpoolP512r1,512,DIVIDE,2000,0.000449599,4448408.469,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,brainpoolP512r1,512,DIVIDE,2000,0.000572538,3493217.917,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,brainpoolP512r1,512,DIVIDE,2000,0.000421237,4747921.005,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,brainpoolP512r1,512,DIVIDE,2000,0.000544206,3675078.922,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,brainpoolP512r1,512,ISQRT,312,0.000035757,8725564.231,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,brainpoolP512r1,512,ISQRT,312,0.000004479,69658405.269,0
opencl-kernel,gfx1010:xnack-,GPU,w8,brainpoolP512r1,512,ISQRT,2000,0.029797670,67119.342,0
opencl-e2e,gfx1010:xnack-,GPU,w8,brainpoolP512r1,512,ISQRT,2000,0.029984169,66701.865,0
opencl-kernel,gfx1010:xnack-,GPU,w16,brainpoolP512r1,512,ISQRT,2000,0.005238087,381818.782,0
opencl-e2e,gfx1010:xnack-,GPU,w16,brainpoolP512r1,512,ISQRT,2000,0.005379311,371794.827,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,brainpoolP512r1,512,ISQRT,2000,0.002239828,892925.707,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,brainpoolP512r1,512,ISQRT,2000,0.002358521,847989.058,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,brainpoolP512r1,512,ISQRT,2000,0.002219018,901299.584,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,brainpoolP512r1,512,ISQRT,2000,0.002350031,851052.603,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,brainpoolP512r1,512,ISQRT,2000,0.002236243,894357.187,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,brainpoolP512r1,512,ISQRT,2000,0.002365573,845461.121,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,brainpoolP512r1,512,ISQRT,2000,0.002247707,889795.690,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,brainpoolP512r1,512,ISQRT,2000,0.002380976,839991.667,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL_R2,10000,0.001404794,7118481.429,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL_R2,10000,0.000176451,56672957.360,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL_R2,10000,0.000468348,21351644.502,0
opencl-kernel,gfx1010:xnack-,GPU,w8,brainpoolP512r1,512,MODMUL_R2,10000,0.000231193,43253904.748,0
opencl-e2e,gfx1010:xnack-,GPU,w8,brainpoolP512r1,512,MODMUL_R2,10000,0.000526517,18992739.075,0
opencl-kernel,gfx1010:xnack-,GPU,w16,brainpoolP512r1,512,MODMUL_R2,10000,0.000093866,106534847.500,0
opencl-e2e,gfx1010:xnack-,GPU,w16,brainpoolP512r1,512,MODMUL_R2,10000,0.000411621,24294192.958,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,brainpoolP512r1,512,MODMUL_R2,10000,0.000101140,98872849.545,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,brainpoolP512r1,512,MODMUL_R2,10000,0.000398196,25113260.800,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,brainpoolP512r1,512,MODMUL_R2,10000,0.000081121,123272642.457,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,brainpoolP512r1,512,MODMUL_R2,10000,0.000402190,24863870.314,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,brainpoolP512r1,512,MODMUL_R2,10000,0.000100908,99100170.305,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,brainpoolP512r1,512,MODMUL_R2,10000,0.000396029,25250676.083,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,brainpoolP512r1,512,MODMUL_R2,10000,0.000070512,141819831.898,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,brainpoolP512r1,512,MODMUL_R2,10000,0.000403463,24785420.218,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,p1024,1024,ADD,5000,0.000080070,62445360.335,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,p1024,1024,ADD,5000,0.000312627,15993500.243,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,p1024,1024,ADD,5000,0.000352561,14181942.983,0
opencl-kernel,gfx1010:xnack-,GPU,w8,p1024,1024,ADD,5000,0.000171030,29234637.190,0
opencl-e2e,gfx1010:xnack-,GPU,w8,p1024,1024,ADD,5000,0.000465522,10740630.947,0
opencl-kernel,gfx1010:xnack-,GPU,w16,p1024,1024,ADD,5000,0.000114194,43785137.598,0
opencl-e2e,gfx1010:xnack-,GPU,w16,p1024,1024,ADD,5000,0.000418305,11953000.806,0
opencl-kernel,gfx1010:xnack-,GPU,w32,p1024,1024,ADD,5000,0.000090780,55078210.960,0
opencl-e2e,gfx1010:xnack-,GPU,w32,p1024,1024,ADD,5000,0.000394400,12677484.787,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,p1024,1024,ADD,5000,0.000088125,56737588.622,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,p1024,1024,ADD,5000,0.000401392,12456650.856,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,p1024,1024,ADD,5000,0.000085910,58200442.288,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,p1024,1024,ADD,5000,0.000382243,13080684.265,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,p1024,1024,ADD,5000,0.000073047,68449081.048,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,p1024,1024,ADD,5000,0.000365913,13664450.301,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,p1024,1024,ADD,5000,0.000056796,88034368.492,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,p1024,1024,ADD,5000,0.000362627,13788272.800,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,p1024,1024,SUBTRACT,5000,0.000061114,81814314.226,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,p1024,1024,SUBTRACT,5000,0.000014216,351716377.050,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,p1024,1024,SUBTRACT,5000,0.000025879,193206847.524,0
opencl-kernel,gfx1010:xnack-,GPU,w8,p1024,1024,SUBTRACT,5000,0.000173144,28877697.181,0
opencl-e2e,gfx1010:xnack-,GPU,w8,p1024,1024,SUBTRACT,5000,0.000465100,10750376.265,0
opencl-kernel,gfx1010:xnack-,GPU,w16,p1024,1024,SUBTRACT,5000,0.000118913,42047547.429,0
opencl-e2e,gfx1010:xnack-,GPU,w16,p1024,1024,SUBTRACT,5000,0.000404387,12364393.511,0
opencl-kernel,gfx1010:xnack-,GPU,w32,p1024,1024,SUBTRACT,5000,0.000089868,55637156.699,0
opencl-e2e,gfx1010:xnack-,GPU,w32,p1024,1024,SUBTRACT,5000,0.000386605,12933097.092,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,p1024,1024,SUBTRACT,5000,0.000087654,57042462.392,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,p1024,1024,SUBTRACT,5000,0.000399980,12500625.035,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,p1024,1024,SUBTRACT,5000,0.000075261,66435471.322,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,p1024,1024,SUBTRACT,5000,0.000383154,13049583.203,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,p1024,1024,SUBTRACT,5000,0.000064972,76956227.276,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,p1024,1024,SUBTRACT,5000,0.000360913,13853754.233,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,p1024,1024,SUBTRACT,5000,0.000059391,84187840.018,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,p1024,1024,SUBTRACT,5000,0.000370021,13512746.575,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,p1024,1024,ADDMOD,5000,0.000209603,23854620.411,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,p1024,1024,ADDMOD,5000,0.000042730,117013807.579,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,p1024,1024,ADDMOD,5000,0.000124353,40208117.200,0
opencl-kernel,gfx1010:xnack-,GPU,w8,p1024,1024,ADDMOD,5000,0.000172824,28931166.976,0
opencl-e2e,gfx1010:xnack-,GPU,w8,p1024,1024,ADDMOD,5000,0.000442128,11308942.212,0
opencl-kernel,gfx1010:xnack-,GPU,w16,p1024,1024,ADDMOD,5000,0.000113102,44207883.105,0
opencl-e2e,gfx1010:xnack-,GPU,w16,p1024,1024,ADDMOD,5000,0.000414697,12056995.831,0
opencl-kernel,gfx1010:xnack-,GPU,w32,p1024,1024,ADDMOD,5000,0.000102482,48789055.590,0
opencl-e2e,gfx1010:xnack-,GPU,w32,p1024,1024,ADDMOD,5000,0.000381956,13090513.045,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,p1024,1024,ADDMOD,5000,0.000064702,77277363.741,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,p1024,1024,ADDMOD,5000,0.000346299,14438389.944,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,p1024,1024,ADDMOD,5000,0.000059401,84173666.725,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,p1024,1024,ADDMOD,5000,0.000347638,14382777.478,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,p1024,1024,ADDMOD,5000,0.000063138,79191612.201,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,p1024,1024,ADDMOD,5000,0.000353960,14125889.923,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,p1024,1024,ADDMOD,5000,0.000065221,76662424.546,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,p1024,1024,ADDMOD,5000,0.000329996,15151698.813,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,p1024,1024,SUBTRACTMOD,5000,0.000174126,28714838.672,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,p1024,1024,SUBTRACTMOD,5000,0.000038603,129523612.130,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,p1024,1024,SUBTRACTMOD,5000,0.000119063,41994574.281,0
opencl-kernel,gfx1010:xnack-,GPU,w8,p1024,1024,SUBTRACTMOD,5000,0.000175809,28439954.726,0
opencl-e2e,gfx1010:xnack-,GPU,w8,p1024,1024,SUBTRACTMOD,5000,0.000469359,10652826.517,0
opencl-kernel,gfx1010:xnack-,GPU,w16,p1024,1024,SUBTRACTMOD,5000,0.000117090,42702194.892,0
opencl-e2e,gfx1010:xnack-,GPU,w16,p1024,1024,SUBTRACTMOD,5000,0.000394189,12684270.744,0
opencl-kernel,gfx1010:xnack-,GPU,w32,p1024,1024,SUBTRACTMOD,5000,0.000089828,55661931.692,0
opencl-e2e,gfx1010:xnack-,GPU,w32,p1024,1024,SUBTRACTMOD,5000,0.000381776,13096684.966,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,p1024,1024,SUBTRACTMOD,5000,0.000059642,83833540.128,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,p1024,1024,SUBTRACTMOD,5000,0.000347742,14378476.005,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,p1024,1024,SUBTRACTMOD,5000,0.000059681,83778757.642,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,p1024,1024,SUBTRACTMOD,5000,0.000343841,14541604.990,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,p1024,1024,SUBTRACTMOD,5000,0.000046467,107603244.909,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,p1024,1024,SUBTRACTMOD,5000,0.000342209,14610954.130,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,p1024,1024,SUBTRACTMOD,5000,0.000064280,77784692.450,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,p1024,1024,SUBTRACTMOD,5000,0.000333974,14971225.320,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,p1024,1024,MULTIPLYOPERANDSCANNING,5000,0.000662322,7549198.124,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,p1024,1024,MULTIPLYOPERANDSCANNING,5000,0.000092643,53970618.386,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,p1024,1024,MULTIPLYOPERANDSCANNING,5000,0.000487765,10250838.006,0
opencl-kernel,gfx1010:xnack-,GPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,5000,0.019215611,260205.101,0
opencl-e2e,gfx1010:xnack-,GPU,w8,p1024,1024,MULTIPLYOPERANDSCANNING,5000,0.019658199,254346.800,0
opencl-kernel,gfx1010:xnack-,GPU,w16,p1024,1024,MULTIPLYOPERANDSCANNING,5000,0.003709923,1347736.867,0
opencl-e2e,gfx1010:xnack-,GPU,w16,p1024,1024,MULTIPLYOPERANDSCANNING,5000,0.004069297,1228713.461,0
opencl-kernel,gfx1010:xnack-,GPU,w32,p1024,1024,MULTIPLYOPERANDSCANNING,5000,0.001166447,4286521.376,0
opencl-e2e,gfx1010:xnack-,GPU,w32,p1024,1024,MULTIPLYOPERANDSCANNING,5000,0.001494611,3345352.068,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,p1024,1024,MULTIPLYOPERANDSCANNING,5000,0.000085229,58665477.728,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,p1024,1024,MULTIPLYOPERANDSCANNING,5000,0.000417382,11979433.708,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,p1024,1024,MULTIPLYOPERANDSCANNING,5000,0.000188752,26489785.509,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,p1024,1024,MULTIPLYOPERANDSCANNING,5000,0.000493260,10136641.934,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,p1024,1024,MULTIPLYOPERANDSCANNING,5000,0.000074930,66728946.830,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,p1024,1024,MULTIPLYOPERANDSCANNING,5000,0.000389116,12849638.670,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,p1024,1024,MULTIPLYOPERANDSCANNING,5000,0.000091130,54866674.058,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,p1024,1024,MULTIPLYOPERANDSCANNING,5000,0.000393575,12704058.944,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,p1024,1024,MULTIPLYPRODUCTSCANNING,5000,0.000652414,7663845.350,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,5000,0.000461055,10844693.150,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,5000,0.000519804,9619010.241,0
opencl-kernel,gfx1010:xnack-,GPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,5000,0.010623031,470675.460,0
opencl-e2e,gfx1010:xnack-,GPU,w8,p1024,1024,MULTIPLYPRODUCTSCANNING,5000,0.011045020,452692.707,0
opencl-kernel,gfx1010:xnack-,GPU,w16,p1024,1024,MULTIPLYPRODUCTSCANNING,5000,0.002750465,1817874.432,0
opencl-e2e,gfx1010:xnack-,GPU,w16,p1024,1024,MULTIPLYPRODUCTSCANNING,5000,0.003109989,1607722.728,0
opencl-kernel,gfx1010:xnack-,GPU,w32,p1024,1024,MULTIPLYPRODUCTSCANNING,5000,0.000766597,6522331.810,0
opencl-e2e,gfx1010:xnack-,GPU,w32,p1024,1024,MULTIPLYPRODUCTSCANNING,5000,0.001105041,4524718.993,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,p1024,1024,MULTIPLYPRODUCTSCANNING,5000,0.000763671,6547322.080,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,p1024,1024,MULTIPLYPRODUCTSCANNING,5000,0.001089823,4587900.971,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,p1024,1024,MULTIPLYPRODUCTSCANNING,5000,0.000775827,6444735.748,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,p1024,1024,MULTIPLYPRODUCTSCANNING,5000,0.001115450,4482495.853,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,p1024,1024,MULTIPLYPRODUCTSCANNING,5000,0.000191628,26092220.307,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,p1024,1024,MULTIPLYPRODUCTSCANNING,5000,0.000505223,9896619.911,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,p1024,1024,MULTIPLYPRODUCTSCANNING,5000,0.000211084,23687252.462,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,p1024,1024,MULTIPLYPRODUCTSCANNING,5000,0.000523638,9548581.276,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,p1024,1024,MONTGOMERYMULTIPLICATION,5000,0.004521910,1105727.447,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,p1024,1024,MONTGOMERYMULTIPLICATION,5000,0.001240035,4032144.254,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,p1024,1024,MONTGOMERYMULTIPLICATION,5000,0.000609393,8204885.846,0
opencl-kernel,gfx1010:xnack-,GPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,5000,0.002241866,2230284.950,0
opencl-e2e,gfx1010:xnack-,GPU,w8,p1024,1024,MONTGOMERYMULTIPLICATION,5000,0.002589987,1930511.620,0
opencl-kernel,gfx1010:xnack-,GPU,w16,p1024,1024,MONTGOMERYMULTIPLICATION,5000,0.000131597,37994787.097,0
opencl-e2e,gfx1010:xnack-,GPU,w16,p1024,1024,MONTGOMERYMULTIPLICATION,5000,0.000421841,11852807.096,0
opencl-kernel,gfx1010:xnack-,GPU,w32,p1024,1024,MONTGOMERYMULTIPLICATION,5000,0.000088145,56724714.918,0
opencl-e2e,gfx1010:xnack-,GPU,w32,p1024,1024,MONTGOMERYMULTIPLICATION,5000,0.000389300,12843565.372,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,p1024,1024,MONTGOMERYMULTIPLICATION,5000,0.000106640,46886721.675,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,p1024,1024,MONTGOMERYMULTIPLICATION,5000,0.000373390,13390824.605,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,p1024,1024,MONTGOMERYMULTIPLICATION,5000,0.000086622,57722056.583,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,p1024,1024,MONTGOMERYMULTIPLICATION,5000,0.000350133,14280287.786,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,p1024,1024,MONTGOMERYMULTIPLICATION,5000,0.000099074,50467327.457,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,p1024,1024,MONTGOMERYMULTIPLICATION,5000,0.000378876,13196929.866,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,p1024,1024,MONTGOMERYMULTIPLICATION,5000,0.000088986,56188613.772,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,p1024,1024,MONTGOMERYMULTIPLICATION,5000,0.000358399,13950931.779,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,p1024,1024,COMPARE,5000,0.000031830,157084511.409,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,p1024,1024,COMPARE,5000,0.000204554,24443423.253,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,p1024,1024,COMPARE,5000,0.000277169,18039535.445,0
opencl-kernel,gfx1010:xnack-,GPU,w8,p1024,1024,COMPARE,5000,0.000079549,62854341.334,0
opencl-e2e,gfx1010:xnack-,GPU,w8,p1024,1024,COMPARE,5000,0.000383177,13048799.902,0
opencl-kernel,gfx1010:xnack-,GPU,w16,p1024,1024,COMPARE,5000,0.000075722,66031008.144,0
opencl-e2e,gfx1010:xnack-,GPU,w16,p1024,1024,COMPARE,5000,0.000362921,13777103.005,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,p1024,1024,COMPARE,5000,0.000068188,73326684.925,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,p1024,1024,COMPARE,5000,0.000345718,14462654.531,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,p1024,1024,COMPARE,5000,0.000066564,75115678.135,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,p1024,1024,COMPARE,5000,0.000335265,14913575.823,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,p1024,1024,COMPARE,5000,0.000077725,64329366.558,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,p1024,1024,COMPARE,5000,0.000353880,14129083.299,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,p1024,1024,COMPARE,5000,0.000057888,86373687.267,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,p1024,1024,COMPARE,5000,0.000335146,14918871.179,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,p1024,1024,REDUCE,625,0.000014868,42036588.664,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,p1024,1024,REDUCE,625,0.000241744,2585379.575,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,p1024,1024,REDUCE,625,0.000317556,1968156.797,0
opencl-kernel,gfx1010:xnack-,GPU,w8,p1024,1024,REDUCE,2000,0.001313799,1522302.879,0
opencl-e2e,gfx1010:xnack-,GPU,w8,p1024,1024,REDUCE,2000,0.001453219,1376255.059,0
opencl-kernel,gfx1010:xnack-,GPU,w16,p1024,1024,REDUCE,2000,0.000309229,6467698.694,0
opencl-e2e,gfx1010:xnack-,GPU,w16,p1024,1024,REDUCE,2000,0.000442820,4516507.836,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,p1024,1024,REDUCE,2000,0.000237275,8429038.035,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,p1024,1024,REDUCE,2000,0.000369883,5407115.222,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,p1024,1024,REDUCE,2000,0.000242232,8256547.446,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,p1024,1024,REDUCE,2000,0.000382063,5234738.774,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,p1024,1024,REDUCE,2000,0.000275424,7261531.313,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,p1024,1024,REDUCE,2000,0.000404074,4949588.442,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,p1024,1024,REDUCE,2000,0.000229178,8726841.139,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,p1024,1024,REDUCE,2000,0.000382094,5234314.074,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,p1024,1024,MODMUL,312,0.000136156,2291489.173,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,p1024,1024,MODMUL,312,0.000329127,947962.337,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,p1024,1024,MODMUL,312,0.000406993,766597.951,0
opencl-kernel,gfx1010:xnack-,GPU,w8,p1024,1024,MODMUL,2000,0.009000614,222207.063,0
opencl-e2e,gfx1010:xnack-,GPU,w8,p1024,1024,MODMUL,2000,0.009190279,217621.250,0
opencl-kernel,gfx1010:xnack-,GPU,w16,p1024,1024,MODMUL,2000,0.001032074,1937845.542,0
opencl-e2e,gfx1010:xnack-,GPU,w16,p1024,1024,MODMUL,2000,0.001193417,1675860.156,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,p1024,1024,MODMUL,2000,0.000670798,2981523.498,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,p1024,1024,MODMUL,2000,0.000848520,2357045.208,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,p1024,1024,MODMUL,2000,0.000603726,3312761.087,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,p1024,1024,MODMUL,2000,0.000748897,2670594.221,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,p1024,1024,MODMUL,2000,0.000762162,2624114.034,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,p1024,1024,MODMUL,2000,0.000912233,2192422.330,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,p1024,1024,MODMUL,2000,0.000611191,3272299.495,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,p1024,1024,MODMUL,2000,0.000785938,2544729.991,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,p1024,1024,MODEXP,78,0.017615015,4428.041,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,p1024,1024,MODEXP,78,0.002640972,29534.580,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,p1024,1024,MODEXP,78,0.002142668,36403.213,0
opencl-kernel,gfx1010:xnack-,GPU,w8,p1024,1024,MODEXP,2000,1.606523958,1244.924,0
opencl-e2e,gfx1010:xnack-,GPU,w8,p1024,1024,MODEXP,2000,1.607040274,1244.524,0
opencl-kernel,gfx1010:xnack-,GPU,w16,p1024,1024,MODEXP,2000,0.099312113,20138.530,0
opencl-e2e,gfx1010:xnack-,GPU,w16,p1024,1024,MODEXP,2000,0.099565207,20087.338,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,p1024,1024,MODEXP,2000,0.066519813,30066.230,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,p1024,1024,MODEXP,2000,0.066722393,29974.944,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,p1024,1024,MODEXP,2000,0.038437606,52032.377,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,p1024,1024,MODEXP,2000,0.038610189,51799.798,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,p1024,1024,MODEXP,2000,0.066609080,30025.936,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,p1024,1024,MODEXP,2000,0.066873122,29907.382,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,p1024,1024,MODEXP,2000,0.040763426,49063.590,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,p1024,1024,MODEXP,2000,0.040807749,49010.299,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,p1024,1024,EXPONENTIATION,78,0.002430799,32088.215,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,p1024,1024,EXPONENTIATION,78,0.000981020,79509.082,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,p1024,1024,EXPONENTIATION,78,0.007859089,9924.814,0
opencl-kernel,gfx1010:xnack-,GPU,w8,p1024,1024,EXPONENTIATION,2000,2.378530077,840.855,0
opencl-e2e,gfx1010:xnack-,GPU,w8,p1024,1024,EXPONENTIATION,2000,2.378803374,840.759,0
opencl-kernel,gfx1010:xnack-,GPU,w16,p1024,1024,EXPONENTIATION,2000,0.581247359,3440.876,0
opencl-e2e,gfx1010:xnack-,GPU,w16,p1024,1024,EXPONENTIATION,2000,0.581332991,3440.369,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,p1024,1024,EXPONENTIATION,2000,0.149341703,13392.107,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,p1024,1024,EXPONENTIATION,2000,0.149550193,13373.436,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,p1024,1024,EXPONENTIATION,2000,0.146763818,13627.337,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,p1024,1024,EXPONENTIATION,2000,0.146918459,13612.993,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,p1024,1024,EXPONENTIATION,2000,0.150082121,13326.038,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,p1024,1024,EXPONENTIATION,2000,0.150562266,13283.541,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,p1024,1024,EXPONENTIATION,2000,0.146238219,13676.315,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,p1024,1024,EXPONENTIATION,2000,0.146528751,13649.198,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,p1024,1024,DIVIDE,625,0.000027782,22496580.459,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,p1024,1024,DIVIDE,625,0.000269896,2315706.791,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,p1024,1024,DIVIDE,625,0.000239899,2605263.048,0
opencl-kernel,gfx1010:xnack-,GPU,w8,p1024,1024,DIVIDE,2000,0.017118507,116832.619,0
opencl-e2e,gfx1010:xnack-,GPU,w8,p1024,1024,DIVIDE,2000,0.017255563,115904.651,0
opencl-kernel,gfx1010:xnack-,GPU,w16,p1024,1024,DIVIDE,2000,0.005543840,360760.772,0
opencl-e2e,gfx1010:xnack-,GPU,w16,p1024,1024,DIVIDE,2000,0.005800863,344776.286,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,p1024,1024,DIVIDE,2000,0.001415143,1413284.735,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,p1024,1024,DIVIDE,2000,0.001625626,1230295.283,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,p1024,1024,DIVIDE,2000,0.001490610,1341732.579,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,p1024,1024,DIVIDE,2000,0.001686495,1185891.450,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,p1024,1024,DIVIDE,2000,0.001437302,1391496.011,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,p1024,1024,DIVIDE,2000,0.001640992,1218774.985,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,p1024,1024,DIVIDE,2000,0.001532131,1305371.408,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,p1024,1024,DIVIDE,2000,0.001717948,1164179.591,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,p1024,1024,ISQRT,156,0.000037861,4120334.908,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,p1024,1024,ISQRT,156,0.000330230,472398.026,0
opencl-kernel,gfx1010:xnack-,GPU,w8,p1024,1024,ISQRT,2000,0.120736264,16565.031,0
opencl-e2e,gfx1010:xnack-,GPU,w8,p1024,1024,ISQRT,2000,0.120972165,16532.729,0
opencl-kernel,gfx1010:xnack-,GPU,w16,p1024,1024,ISQRT,2000,0.066854371,29915.770,0
opencl-e2e,gfx1010:xnack-,GPU,w16,p1024,1024,ISQRT,2000,0.066931375,29881.352,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,p1024,1024,ISQRT,2000,0.009133540,218973.147,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,p1024,1024,ISQRT,2000,0.009322694,214530.263,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,p1024,1024,ISQRT,2000,0.009224312,216818.338,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,p1024,1024,ISQRT,2000,0.009405389,212644.049,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,p1024,1024,ISQRT,2000,0.009246324,216302.176,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,p1024,1024,ISQRT,2000,0.009442560,211806.968,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,p1024,1024,ISQRT,2000,0.009379043,213241.372,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,p1024,1024,ISQRT,2000,0.009555883,209295.154,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,p1024,1024,MODMUL_R2,5000,0.001982547,2522008.305,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,p1024,1024,MODMUL_R2,5000,0.000917240,5451136.017,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,p1024,1024,MODMUL_R2,5000,0.001367415,3656534.410,0
opencl-kernel,gfx1010:xnack-,GPU,w8,p1024,1024,MODMUL_R2,5000,0.001089901,4587572.633,0
opencl-e2e,gfx1010:xnack-,GPU,w8,p1024,1024,MODMUL_R2,5000,0.001447239,3454854.381,0
opencl-kernel,gfx1010:xnack-,GPU,w16,p1024,1024,MODMUL_R2,5000,0.000167664,29821547.867,0
opencl-e2e,gfx1010:xnack-,GPU,w16,p1024,1024,MODMUL_R2,5000,0.000478156,10456838.357,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,p1024,1024,MODMUL_R2,5000,0.000112180,44571224.873,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,p1024,1024,MODMUL_R2,5000,0.000441427,11326901.165,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,p1024,1024,MODMUL_R2,5000,0.000108231,46197485.003,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,p1024,1024,MODMUL_R2,5000,0.000412008,12135686.692,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,p1024,1024,MODMUL_R2,5000,0.000131585,37998252.091,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,p1024,1024,MODMUL_R2,5000,0.000447044,11184581.375,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,p1024,1024,MODMUL_R2,5000,0.000113362,44106490.587,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,p1024,1024,MODMUL_R2,5000,0.000425495,11751019.402,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,p2048,2048,ADD,2500,0.000055945,44686745.923,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,p2048,2048,ADD,2500,0.000301666,8287311.131,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,p2048,2048,ADD,2500,0.000472346,5292730.330,0
opencl-kernel,gfx1010:xnack-,GPU,w8,p2048,2048,ADD,2500,0.000631251,3960389.766,0
opencl-e2e,gfx1010:xnack-,GPU,w8,p2048,2048,ADD,2500,0.000957683,2610467.138,0
opencl-kernel,gfx1010:xnack-,GPU,w16,p2048,2048,ADD,2500,0.000600756,4161423.275,0
opencl-e2e,gfx1010:xnack-,GPU,w16,p2048,2048,ADD,2500,0.000930233,2687498.724,0
opencl-kernel,gfx1010:xnack-,GPU,w32,p2048,2048,ADD,2500,0.000116638,21433838.027,0
opencl-e2e,gfx1010:xnack-,GPU,w32,p2048,2048,ADD,2500,0.000448560,5573390.405,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,p2048,2048,ADD,2500,0.000115006,21737996.238,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,p2048,2048,ADD,2500,0.000430547,5806566.993,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,p2048,2048,ADD,2500,0.000128830,19405417.994,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,p2048,2048,ADD,2500,0.000433318,5769434.919,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,p2048,2048,ADD,2500,0.000082143,30434729.558,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,p2048,2048,ADD,2500,0.000378757,6600538.070,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,p2048,2048,ADD,2500,0.000073898,33830415.008,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,p2048,2048,ADD,2500,0.000387113,6458062.629,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,p2048,2048,SUBTRACT,2500,0.000058500,42735042.801,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,p2048,2048,SUBTRACT,2500,0.000214532,11653273.173,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,p2048,2048,SUBTRACT,2500,0.000268704,9303918.064,0
opencl-kernel,gfx1010:xnack-,GPU,w8,p2048,2048,SUBTRACT,2500,0.000587921,4252271.989,0
opencl-e2e,gfx1010:xnack-,GPU,w8,p2048,2048,SUBTRACT,2500,0.000966660,2586224.733,0
opencl-kernel,gfx1010:xnack-,GPU,w16,p2048,2048,SUBTRACT,2500,0.000430527,5806836.738,0
opencl-e2e,gfx1010:xnack-,GPU,w16,p2048,2048,SUBTRACT,2500,0.000798397,3131274.291,0
opencl-kernel,gfx1010:xnack-,GPU,w32,p2048,2048,SUBTRACT,2500,0.000134743,18553839.533,0
opencl-e2e,gfx1010:xnack-,GPU,w32,p2048,2048,SUBTRACT,2500,0.000406342,6152452.861,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,p2048,2048,SUBTRACT,2500,0.000113643,21998715.260,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,p2048,2048,SUBTRACT,2500,0.000412884,6054969.433,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,p2048,2048,SUBTRACT,2500,0.000127408,19622001.754,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,p2048,2048,SUBTRACT,2500,0.000404705,6177339.054,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,p2048,2048,SUBTRACT,2500,0.000073377,34070621.674,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,p2048,2048,SUBTRACT,2500,0.000376683,6636880.347,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,p2048,2048,SUBTRACT,2500,0.000066374,37665350.778,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,p2048,2048,SUBTRACT,2500,0.000353019,7081771.799,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,p2048,2048,ADDMOD,2500,0.000134102,18642525.846,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,p2048,2048,ADDMOD,2500,0.000246492,10142316.993,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,p2048,2048,ADDMOD,2500,0.000378860,6598743.598,0
opencl-kernel,gfx1010:xnack-,GPU,w8,p2048,2048,ADDMOD,2500,0.000626974,3987406.175,0
opencl-e2e,gfx1010:xnack-,GPU,w8,p2048,2048,ADDMOD,2500,0.001072419,2331178.392,0
opencl-kernel,gfx1010:xnack-,GPU,w16,p2048,2048,ADDMOD,2500,0.000575639,4342999.694,0
opencl-e2e,gfx1010:xnack-,GPU,w16,p2048,2048,ADDMOD,2500,0.000934752,2674506.179,0
opencl-kernel,gfx1010:xnack-,GPU,w32,p2048,2048,ADDMOD,2500,0.000117390,21296532.922,0
opencl-e2e,gfx1010:xnack-,GPU,w32,p2048,2048,ADDMOD,2500,0.000413616,6044253.606,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,p2048,2048,ADDMOD,2500,0.000077796,32135328.330,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,p2048,2048,ADDMOD,2500,0.000383830,6513300.160,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,p2048,2048,ADDMOD,2500,0.000064831,38561799.323,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,p2048,2048,ADDMOD,2500,0.000362436,6897769.543,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,p2048,2048,ADDMOD,2500,0.000070452,35485153.104,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,p2048,2048,ADDMOD,2500,0.000363148,6884245.544,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,p2048,2048,ADDMOD,2500,0.000066114,37813473.571,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,p2048,2048,ADDMOD,2500,0.000371043,6737763.546,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,p2048,2048,SUBTRACTMOD,2500,0.000119795,20868984.527,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,p2048,2048,SUBTRACTMOD,2500,0.000278332,8982078.956,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,p2048,2048,SUBTRACTMOD,2500,0.000452007,5530887.796,0
opencl-kernel,gfx1010:xnack-,GPU,w8,p2048,2048,SUBTRACTMOD,2500,0.000739665,3379908.472,0
opencl-e2e,gfx1010:xnack-,GPU,w8,p2048,2048,SUBTRACTMOD,2500,0.001053833,2372292.384,0
opencl-kernel,gfx1010:xnack-,GPU,w16,p2048,2048,SUBTRACTMOD,2500,0.000436588,5726222.435,0
opencl-e2e,gfx1010:xnack-,GPU,w16,p2048,2048,SUBTRACTMOD,2500,0.000819957,3048940.371,0
opencl-kernel,gfx1010:xnack-,GPU,w32,p2048,2048,SUBTRACTMOD,2500,0.000130525,19153418.881,0
opencl-e2e,gfx1010:xnack-,GPU,w32,p2048,2048,SUBTRACTMOD,2500,0.000418925,5967655.310,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,p2048,2048,SUBTRACTMOD,2500,0.000074610,33507572.731,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,p2048,2048,SUBTRACTMOD,2500,0.000367048,6811098.275,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,p2048,2048,SUBTRACTMOD,2500,0.000058469,42757700.572,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,p2048,2048,SUBTRACTMOD,2500,0.000349833,7146266.929,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,p2048,2048,SUBTRACTMOD,2500,0.000063188,39564474.402,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,p2048,2048,SUBTRACTMOD,2500,0.000341378,7323260.432,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,p2048,2048,SUBTRACTMOD,2500,0.000065011,38455030.651,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,p2048,2048,SUBTRACTMOD,2500,0.000350485,7132972.878,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,p2048,2048,MULTIPLYOPERANDSCANNING,2500,0.001095044,2283013.285,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,p2048,2048,MULTIPLYOPERANDSCANNING,2500,0.000725841,3444280.497,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,p2048,2048,MULTIPLYOPERANDSCANNING,2500,0.000667923,3742946.417,0
opencl-kernel,gfx1010:xnack-,GPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,2500,0.119266041,20961.541,0
opencl-e2e,gfx1010:xnack-,GPU,w8,p2048,2048,MULTIPLYOPERANDSCANNING,2500,0.131448610,19018.839,0
opencl-kernel,gfx1010:xnack-,GPU,w16,p2048,2048,MULTIPLYOPERANDSCANNING,2500,0.025952684,96329.150,0
opencl-e2e,gfx1010:xnack-,GPU,w16,p2048,2048,MULTIPLYOPERANDSCANNING,2500,0.023477985,106482.733,0
opencl-kernel,gfx1010:xnack-,GPU,w32,p2048,2048,MULTIPLYOPERANDSCANNING,2500,0.003971624,629465.428,0
opencl-e2e,gfx1010:xnack-,GPU,w32,p2048,2048,MULTIPLYOPERANDSCANNING,2500,0.004555890,548740.202,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,p2048,2048,MULTIPLYOPERANDSCANNING,2500,0.000496240,5037884.896,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,p2048,2048,MULTIPLYOPERANDSCANNING,2500,0.000840735,2973588.586,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,p2048,2048,MULTIPLYOPERANDSCANNING,2500,0.000659821,3788906.387,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,p2048,2048,MULTIPLYOPERANDSCANNING,2500,0.001008972,2477769.453,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,p2048,2048,MULTIPLYOPERANDSCANNING,2500,0.000488762,5114963.931,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,p2048,2048,MULTIPLYOPERANDSCANNING,2500,0.000850678,2938832.320,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,p2048,2048,MULTIPLYOPERANDSCANNING,2500,0.000200855,12446789.973,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,p2048,2048,MULTIPLYOPERANDSCANNING,2500,0.000518358,4822921.609,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,p2048,2048,MULTIPLYPRODUCTSCANNING,2500,0.001098030,2276804.823,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,2500,0.000535634,4667366.149,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,2500,0.000607249,4116927.323,0
opencl-kernel,gfx1010:xnack-,GPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,2500,0.028351543,88178.622,0
opencl-e2e,gfx1010:xnack-,GPU,w8,p2048,2048,MULTIPLYPRODUCTSCANNING,2500,0.033462241,74711.075,0
opencl-kernel,gfx1010:xnack-,GPU,w16,p2048,2048,MULTIPLYPRODUCTSCANNING,2500,0.008627341,289776.421,0
opencl-e2e,gfx1010:xnack-,GPU,w16,p2048,2048,MULTIPLYPRODUCTSCANNING,2500,0.009473777,263886.304,0
opencl-kernel,gfx1010:xnack-,GPU,w32,p2048,2048,MULTIPLYPRODUCTSCANNING,2500,0.002192900,1140042.866,0
opencl-e2e,gfx1010:xnack-,GPU,w32,p2048,2048,MULTIPLYPRODUCTSCANNING,2500,0.002128770,1174387.087,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,p2048,2048,MULTIPLYPRODUCTSCANNING,2500,0.002712855,921538.379,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,p2048,2048,MULTIPLYPRODUCTSCANNING,2500,0.003087637,809680.672,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,p2048,2048,MULTIPLYPRODUCTSCANNING,2500,0.002699936,925947.874,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,p2048,2048,MULTIPLYPRODUCTSCANNING,2500,0.002536041,985788.479,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,p2048,2048,MULTIPLYPRODUCTSCANNING,2500,0.000679378,3679836.555,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,p2048,2048,MULTIPLYPRODUCTSCANNING,2500,0.001002061,2494858.097,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,p2048,2048,MULTIPLYPRODUCTSCANNING,2500,0.000671133,3725044.066,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,p2048,2048,MULTIPLYPRODUCTSCANNING,2500,0.000996962,2507618.144,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,p2048,2048,MONTGOMERYMULTIPLICATION,2500,0.007470080,334668.437,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,p2048,2048,MONTGOMERYMULTIPLICATION,2500,0.001898861,1316578.728,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,p2048,2048,MONTGOMERYMULTIPLICATION,2500,0.000965401,2589597.483,0
opencl-kernel,gfx1010:xnack-,GPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,2500,0.022188002,112673.507,0
opencl-e2e,gfx1010:xnack-,GPU,w8,p2048,2048,MONTGOMERYMULTIPLICATION,2500,0.023669797,105619.833,0
opencl-kernel,gfx1010:xnack-,GPU,w16,p2048,2048,MONTGOMERYMULTIPLICATION,2500,0.005527951,452247.135,0
opencl-e2e,gfx1010:xnack-,GPU,w16,p2048,2048,MONTGOMERYMULTIPLICATION,2500,0.005850505,427313.540,0
opencl-kernel,gfx1010:xnack-,GPU,w32,p2048,2048,MONTGOMERYMULTIPLICATION,2500,0.000199644,12522289.671,0
opencl-e2e,gfx1010:xnack-,GPU,w32,p2048,2048,MONTGOMERYMULTIPLICATION,2500,0.000489418,5108107.999,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,p2048,2048,MONTGOMERYMULTIPLICATION,2500,0.000193192,12940494.424,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,p2048,2048,MONTGOMERYMULTIPLICATION,2500,0.000503313,4967088.074,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,p2048,2048,MONTGOMERYMULTIPLICATION,2500,0.000158005,15822284.127,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,p2048,2048,MONTGOMERYMULTIPLICATION,2500,0.000442225,5653230.822,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,p2048,2048,MONTGOMERYMULTIPLICATION,2500,0.000194183,12874453.474,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,p2048,2048,MONTGOMERYMULTIPLICATION,2500,0.000502960,4970574.203,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,p2048,2048,MONTGOMERYMULTIPLICATION,2500,0.000145201,17217512.261,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,p2048,2048,MONTGOMERYMULTIPLICATION,2500,0.000442266,5652706.743,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,p2048,2048,COMPARE,2500,0.000014908,167695198.201,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,p2048,2048,COMPARE,2500,0.000233448,10709022.994,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,p2048,2048,COMPARE,2500,0.000308289,8109274.089,0
opencl-kernel,gfx1010:xnack-,GPU,w8,p2048,2048,COMPARE,2500,0.000287979,8681188.560,0
opencl-e2e,gfx1010:xnack-,GPU,w8,p2048,2048,COMPARE,2500,0.000601096,4159069.433,0
opencl-kernel,gfx1010:xnack-,GPU,w16,p2048,2048,COMPARE,2500,0.000290685,8600374.977,0
opencl-e2e,gfx1010:xnack-,GPU,w16,p2048,2048,COMPARE,2500,0.000600025,4166493.063,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,p2048,2048,COMPARE,2500,0.000073979,33793373.675,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,p2048,2048,COMPARE,2500,0.000368190,6789972.564,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,p2048,2048,COMPARE,2500,0.000075451,33134086.954,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,p2048,2048,COMPARE,2500,0.000371613,6727428.804,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,p2048,2048,COMPARE,2500,0.000082604,30264878.214,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,p2048,2048,COMPARE,2500,0.000384328,6504860.430,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,p2048,2048,COMPARE,2500,0.000085159,29356850.193,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,p2048,2048,COMPARE,2500,0.000377204,6627713.387,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,p2048,2048,REDUCE,312,0.000010159,30711684.225,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,p2048,2048,REDUCE,312,0.000351399,887879.590,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,p2048,2048,REDUCE,312,0.000399650,780683.098,0
opencl-kernel,gfx1010:xnack-,GPU,w8,p2048,2048,REDUCE,2000,0.005077757,393874.697,0
opencl-e2e,gfx1010:xnack-,GPU,w8,p2048,2048,REDUCE,2000,0.005313349,376410.433,0
opencl-kernel,gfx1010:xnack-,GPU,w16,p2048,2048,REDUCE,2000,0.001545496,1294082.935,0
opencl-e2e,gfx1010:xnack-,GPU,w16,p2048,2048,REDUCE,2000,0.001875356,1066464.181,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,p2048,2048,REDUCE,2000,0.000523742,3818674.080,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,p2048,2048,REDUCE,2000,0.000775284,2579699.827,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,p2048,2048,REDUCE,2000,0.000653138,3062140.007,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,p2048,2048,REDUCE,2000,0.000898797,2225196.568,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,p2048,2048,REDUCE,2000,0.000524629,3812217.779,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,p2048,2048,REDUCE,2000,0.000769757,2598222.556,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,p2048,2048,REDUCE,2000,0.000622693,3211855.602,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,p2048,2048,REDUCE,2000,0.000886145,2256966.975,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,p2048,2048,MODMUL,156,0.000201518,774124.396,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,p2048,2048,MODMUL,156,0.000420198,371253.552,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,p2048,2048,MODMUL,156,0.000351479,443838.750,0
opencl-kernel,gfx1010:xnack-,GPU,w8,p2048,2048,MODMUL,2000,0.028394002,70437.411,0
opencl-e2e,gfx1010:xnack-,GPU,w8,p2048,2048,MODMUL,2000,0.028887406,69234.323,0
opencl-kernel,gfx1010:xnack-,GPU,w16,p2048,2048,MODMUL,2000,0.006965315,287137.050,0
opencl-e2e,gfx1010:xnack-,GPU,w16,p2048,2048,MODMUL,2000,0.007165120,279130.008,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,p2048,2048,MODMUL,2000,0.001965554,1017524.830,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,p2048,2048,MODMUL,2000,0.002217527,901905.591,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,p2048,2048,MODMUL,2000,0.002062999,969462.419,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,p2048,2048,MODMUL,2000,0.002296033,871067.620,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,p2048,2048,MODMUL,2000,0.001975278,1012515.707,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,p2048,2048,MODMUL,2000,0.002232508,895853.453,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,p2048,2048,MODMUL,2000,0.002004624,997693.333,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,p2048,2048,MODMUL,2000,0.002246746,890176.282,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,p2048,2048,MODEXP,64,0.106510039,600.882,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,p2048,2048,MODEXP,64,0.015027833,4258.764,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,p2048,2048,MODEXP,64,0.011761217,5441.614,0
opencl-kernel,gfx1010:xnack-,GPU,w8,p2048,2048,MODEXP,2000,19.151992975,104.428,0
opencl-e2e,gfx1010:xnack-,GPU,w8,p2048,2048,MODEXP,2000,19.142296621,104.481,0
opencl-kernel,gfx1010:xnack-,GPU,w16,p2048,2048,MODEXP,2000,3.295250215,606.934,0
opencl-e2e,gfx1010:xnack-,GPU,w16,p2048,2048,MODEXP,2000,3.293566359,607.244,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,p2048,2048,MODEXP,2000,0.510311720,3919.173,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,p2048,2048,MODEXP,2000,0.510754019,3915.779,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,p2048,2048,MODEXP,2000,0.288849167,6924.029,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,p2048,2048,MODEXP,2000,0.289184433,6916.002,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,p2048,2048,MODEXP,2000,0.521837170,3832.613,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,p2048,2048,MODEXP,2000,0.522020531,3831.267,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,p2048,2048,MODEXP,2000,0.289350390,6912.035,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,p2048,2048,MODEXP,2000,0.289707796,6903.508,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,p2048,2048,EXPONENTIATION,64,0.013356429,4791.700,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,p2048,2048,EXPONENTIATION,64,0.002689494,23796.298,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,p2048,2048,EXPONENTIATION,64,0.027939828,2290.637,0
opencl-kernel,gfx1010:xnack-,GPU,w8,p2048,2048,EXPONENTIATION,2000,19.952752551,100.237,0
opencl-e2e,gfx1010:xnack-,GPU,w8,p2048,2048,EXPONENTIATION,2000,19.953358319,100.234,0
opencl-kernel,gfx1010:xnack-,GPU,w16,p2048,2048,EXPONENTIATION,2000,4.754889572,420.620,0
opencl-e2e,gfx1010:xnack-,GPU,w16,p2048,2048,EXPONENTIATION,2000,4.754282369,420.673,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,p2048,2048,EXPONENTIATION,2000,1.164559971,1717.387,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,p2048,2048,EXPONENTIATION,2000,1.164885867,1716.906,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,p2048,2048,EXPONENTIATION,2000,1.142684715,1750.264,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,p2048,2048,EXPONENTIATION,2000,1.143151016,1749.550,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,p2048,2048,EXPONENTIATION,2000,1.168038656,1712.272,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,p2048,2048,EXPONENTIATION,2000,1.168221400,1712.004,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,p2048,2048,EXPONENTIATION,2000,1.139765726,1754.747,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,p2048,2048,EXPONENTIATION,2000,1.139942856,1754.474,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,p2048,2048,DIVIDE,312,0.000018354,16999019.217,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,p2048,2048,DIVIDE,312,0.000236634,1318491.848,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,p2048,2048,DIVIDE,312,0.000341882,912595.574,0
opencl-kernel,gfx1010:xnack-,GPU,w8,p2048,2048,DIVIDE,2000,0.067384882,29680.248,0
opencl-e2e,gfx1010:xnack-,GPU,w8,p2048,2048,DIVIDE,2000,0.067670598,29554.933,0
opencl-kernel,gfx1010:xnack-,GPU,w16,p2048,2048,DIVIDE,2000,0.038172169,52394.193,0
opencl-e2e,gfx1010:xnack-,GPU,w16,p2048,2048,DIVIDE,2000,0.038380629,52109.620,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,p2048,2048,DIVIDE,2000,0.018172804,110054.563,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,p2048,2048,DIVIDE,2000,0.018540560,107871.607,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,p2048,2048,DIVIDE,2000,0.014338071,139488.778,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,p2048,2048,DIVIDE,2000,0.014855436,134630.852,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,p2048,2048,DIVIDE,2000,0.012298272,162624.473,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,p2048,2048,DIVIDE,2000,0.012668183,157875.837,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,p2048,2048,DIVIDE,2000,0.017264290,115846.061,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,p2048,2048,DIVIDE,2000,0.017626636,113464.645,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,p2048,2048,ISQRT,78,0.000033353,2338620.219,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,p2048,2048,ISQRT,78,0.000431288,180853.629,0
opencl-kernel,gfx1010:xnack-,GPU,w8,p2048,2048,ISQRT,2000,0.515017663,3883.362,0
opencl-e2e,gfx1010:xnack-,GPU,w8,p2048,2048,ISQRT,2000,0.515803697,3877.444,0
opencl-kernel,gfx1010:xnack-,GPU,w16,p2048,2048,ISQRT,2000,0.290598568,6882.346,0
opencl-e2e,gfx1010:xnack-,GPU,w16,p2048,2048,ISQRT,2000,0.291126207,6869.873,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,p2048,2048,ISQRT,2000,0.100241457,19951.825,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,p2048,2048,ISQRT,2000,0.099123072,20176.937,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,p2048,2048,ISQRT,2000,0.091127519,21947.267,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,p2048,2048,ISQRT,2000,0.090760054,22036.126,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,p2048,2048,ISQRT,2000,0.099038188,20194.231,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,p2048,2048,ISQRT,2000,0.099662104,20067.808,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,p2048,2048,ISQRT,2000,0.092835208,21543.551,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,p2048,2048,ISQRT,2000,0.092520481,21616.835,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-1t,p2048,2048,MODMUL_R2,2500,0.003228575,774335.427,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,gmp-nt,p2048,2048,MODMUL_R2,2500,0.000833373,2999857.207,0
library,AMD Ryzen 7 5700G with Radeon Graphics,host-cpu,openssl-nt,p2048,2048,MODMUL_R2,2500,0.000910808,2744815.592,0
opencl-kernel,gfx1010:xnack-,GPU,w8,p2048,2048,MODMUL_R2,2500,0.017974247,139087.885,0
opencl-e2e,gfx1010:xnack-,GPU,w8,p2048,2048,MODMUL_R2,2500,0.018180264,137511.755,0
opencl-kernel,gfx1010:xnack-,GPU,w16,p2048,2048,MODMUL_R2,2500,0.003095952,807506.060,0
opencl-e2e,gfx1010:xnack-,GPU,w16,p2048,2048,MODMUL_R2,2500,0.003448333,724987.987,0
opencl-kernel,gfx1010:xnack-,GPU,w32-opt,p2048,2048,MODMUL_R2,2500,0.000313074,7985332.538,0
opencl-e2e,gfx1010:xnack-,GPU,w32-opt,p2048,2048,MODMUL_R2,2500,0.000625006,3999961.600,0
opencl-kernel,gfx1010:xnack-,GPU,w32-o64,p2048,2048,MODMUL_R2,2500,0.000232995,10729843.976,0
opencl-e2e,gfx1010:xnack-,GPU,w32-o64,p2048,2048,MODMUL_R2,2500,0.000538835,4639639.221,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il,p2048,2048,MODMUL_R2,2500,0.000307073,8141386.578,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il,p2048,2048,MODMUL_R2,2500,0.000639845,3907196.274,0
opencl-kernel,gfx1010:xnack-,GPU,w32-il64,p2048,2048,MODMUL_R2,2500,0.000220933,11315647.752,0
opencl-e2e,gfx1010:xnack-,GPU,w32-il64,p2048,2048,MODMUL_R2,2500,0.000531594,4702837.126,0
```
