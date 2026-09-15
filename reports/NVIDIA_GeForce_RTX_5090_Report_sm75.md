# MPA-OpenCL benchmark report - NVIDIA GeForce RTX 5090


> **Note.** Two column groups have been removed from this report: the CGBN
> reference column, which was invalid, and the multi-threaded GMP and OpenSSL
> baselines, which predate the 2026-09-12 timing fix and were understated.
> See `reports/README.md`. The single-threaded GMP column, the OpenCL-on-CPU
> rows and every MPA measurement are unaffected, and every configuration was
> verified word-for-word against GMP before it was timed.


> **Not a vendor-runtime result.** The OpenCL rows were produced by a third-party runtime (`OpenCL 3.0 PoCL HSTR: CUDA-sm_75`), not the GPU vendor's own OpenCL implementation, so the kernels went through a different compiler than on any vendor-ICD host. PoCL caps its PTX target at `sm_75` unless `POCL_CUDA_GPU_ARCH` says otherwise (here: unset), so newer hardware is addressed through forward JIT rather than native codegen. The work-group size was forced to 64 rather than derived from kernel register usage. 
>
> The CGBN column is native CUDA and the GMP/OpenSSL columns are native CPU code, so **only the OpenCL columns carry this handicap**: treat the gap to CGBN as an upper bound and the margin over GMP and OpenSSL as a lower bound. Correctness results are unaffected - a kernel that matches GMP on every word is correct whichever compiler built it.

## 1. System under test

1 OpenCL device(s) exercised with the identical kernels and operands.

### Device 0 - NVIDIA GeForce RTX 5090 (GPU)

| Property | Value |
|---|---|
| Model | NVIDIA GeForce RTX 5090 |
| Type | GPU |
| Vendor | NVIDIA Corporation |
| Device memory | 31.84 GiB |
| Max single allocation | 30.20 GiB |
| Local memory | 48 KiB |
| Global cache | 0 KiB |
| Compute units | 170 |
| Max clock | 2422 MHz |
| Max work-group size | 1024 |
| OpenCL version | OpenCL 3.0 PoCL HSTR: CUDA-sm_75 |
| Driver | 7.0 |

### Host

| Property | Value |
|---|---|
| CPU | AMD Ryzen 9 7950X 16-Core Processor |
| Logical cores | 32 |
| OpenMP threads used | 32 |
| RAM | 120.9 GB |
| OS | Ubuntu 24.04.2 LTS |
| Kernel | 6.18.33.2-microsoft-standard-WSL2 |
| Arch | x86_64 |
| GMP | 6.3.0 |
| OpenSSL | OpenSSL 3.0.13 30 Jan 2024 |

## 2. Method

- Work-group size forced to 64 via MPA_LOCAL_SIZE; the device item count is trimmed to a multiple of it. Runtimes that derive a launchable size from kernel register usage do not need this, and their timings are not directly comparable with these.
- Workload auto-sized from the device and host: --min-items from 700 x compute units, --items from ten times that capped by host RAM. Either flag, given explicitly, overrides its half.
- Base workload 604343 items, scaled down per operator by its cost weight and by modulus size. Device rows honour --min-items (119000) so the GPU is not left idle; the CPU libraries keep the smaller count because a full-width MODEXP there costs minutes. Both counts appear in every row as dev/cpu, and throughput is per-second so they remain comparable.
- 5 timed repetitions, **minimum** reported. Two untimed warm-up launches precede them.
- `kernel` times `clEnqueueNDRangeKernel` + `clFinish` only. `e2e` adds the host->device operand writes and the device->host result read.
- Every OpenCL device runs the same kernels on the same operands, so GPU and CPU-OpenCL columns are directly comparable.
- CPU library baselines (GMP, OpenSSL) run those same operands, with temporaries preallocated outside the timed region, so the figure is the arithmetic and not marshalling. The generator is reseeded per modulus and operation so every backend sees identical inputs.
- OpenSSL rows time the nearest BN primitive, which is not always semantically identical (its Montgomery routine expects Montgomery-domain inputs); they measure comparable work, not identical results. Correctness is judged against GMP only.
- Every device cell is checked word-for-word against GMP before it is timed. A cell that mismatches is reported and excluded from the speedup tables.
- Total wall time 515.1 s.

## 3. Correctness

| Device | Kernel | Configs run | Passed | Mismatched | Build/launch failed |
|---|---|---|---|---|---|
| [0] GPU | `mpaKernel_32bits_opt.cl` (w32-opt) | 75 | 75 | 0 | 0 |

**All configurations correct** - 75 configurations, 0 problems.

## 4. Throughput per device

Operations per second, higher is better. Kernel-only timings.

### Device 0 - NVIDIA GeForce RTX 5090 (GPU)

#### secp256k1 (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 604288 / 604343 | - | - | - | 2.72 G | - | - | - | 85.01 M |
| SUBTRACT | 604288 / 604343 | - | - | - | 2.61 G | - | - | - | 95.14 M |
| ADDMOD | 604288 / 604343 | - | - | - | 2.83 G | - | - | - | 34.96 M |
| SUBTRACTMOD | 604288 / 604343 | - | - | - | 2.84 G | - | - | - | 44.14 M |
| MULTIPLYOPERANDSCANNING | 604288 / 604343 | - | - | - | 2.05 G | - | - | - | 66.48 M |
| MULTIPLYPRODUCTSCANNING | 604288 / 604343 | - | - | - | 2.20 G | - | - | - | 67.26 M |
| MONTGOMERYMULTIPLICATION | 604288 / 604343 | - | - | - | 2.75 G | - | - | - | 9.65 M |
| COMPARE | 604288 / 604343 | - | - | - | 2.85 G | - | - | - | 115.47 M |
| REDUCE | 118976 / 75542 | - | - | - | 553.43 M | - | - | - | 88.54 M |
| MODMUL | 118976 / 37771 | - | - | - | 412.90 M | - | - | - | 16.84 M |
| MODEXP | 118976 / 9442 | - | - | - | 11.36 M | - | - | - | 156.43 k |
| EXPONENTIATION | 118976 / 9442 | - | - | - | 13.01 M | - | - | - | 499.45 k |
| DIVIDE | 118976 / 75542 | - | - | - | 469.03 M | - | - | - | 48.71 M |
| ISQRT | 118976 / 18885 | - | - | - | 68.10 M | - | - | - | 21.45 M |
| MODMUL_R2 | 604288 / 604343 | - | - | - | 2.39 G | - | - | - | 16.78 M |

#### rsa256(composite) (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 604288 / 604343 | - | - | - | 2.61 G | - | - | - | 83.93 M |
| SUBTRACT | 604288 / 604343 | - | - | - | 2.58 G | - | - | - | 95.61 M |
| ADDMOD | 604288 / 604343 | - | - | - | 2.79 G | - | - | - | 40.47 M |
| SUBTRACTMOD | 604288 / 604343 | - | - | - | 2.84 G | - | - | - | 44.27 M |
| MULTIPLYOPERANDSCANNING | 604288 / 604343 | - | - | - | 2.07 G | - | - | - | 67.05 M |
| MULTIPLYPRODUCTSCANNING | 604288 / 604343 | - | - | - | 2.19 G | - | - | - | 65.64 M |
| MONTGOMERYMULTIPLICATION | 604288 / 604343 | - | - | - | 2.73 G | - | - | - | 9.69 M |
| COMPARE | 604288 / 604343 | - | - | - | 3.02 G | - | - | - | 120.30 M |
| REDUCE | 118976 / 75542 | - | - | - | 554.60 M | - | - | - | 55.83 M |
| MODMUL | 118976 / 37771 | - | - | - | 417.22 M | - | - | - | 16.80 M |
| MODEXP | 118976 / 9442 | - | - | - | 12.19 M | - | - | - | 166.94 k |
| EXPONENTIATION | 118976 / 9442 | - | - | - | 13.69 M | - | - | - | 494.94 k |
| DIVIDE | 118976 / 75542 | - | - | - | 452.33 M | - | - | - | 46.33 M |
| ISQRT | 118976 / 18885 | - | - | - | 67.95 M | - | - | - | 21.06 M |
| MODMUL_R2 | 604288 / 604343 | - | - | - | 2.35 G | - | - | - | 16.55 M |

#### brainpoolP512r1 (512-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 302144 / 302171 | - | - | - | 1.23 G | - | - | - | 65.22 M |
| SUBTRACT | 302144 / 302171 | - | - | - | 1.23 G | - | - | - | 73.34 M |
| ADDMOD | 302144 / 302171 | - | - | - | 1.24 G | - | - | - | 36.14 M |
| SUBTRACTMOD | 302144 / 302171 | - | - | - | 1.27 G | - | - | - | 38.27 M |
| MULTIPLYOPERANDSCANNING | 302144 / 302171 | - | - | - | 715.59 M | - | - | - | 31.77 M |
| MULTIPLYPRODUCTSCANNING | 302144 / 302171 | - | - | - | 768.94 M | - | - | - | 31.68 M |
| MONTGOMERYMULTIPLICATION | 302144 / 302171 | - | - | - | 1.12 G | - | - | - | 4.10 M |
| COMPARE | 302144 / 302171 | - | - | - | 1.34 G | - | - | - | 85.63 M |
| REDUCE | 118976 / 37771 | - | - | - | 45.23 M | - | - | - | 51.12 M |
| MODMUL | 118976 / 18885 | - | - | - | 16.11 M | - | - | - | 8.40 M |
| MODEXP | 118976 / 4721 | - | - | - | 1.48 M | - | - | - | 27.22 k |
| EXPONENTIATION | 118976 / 4721 | - | - | - | 1.81 M | - | - | - | 157.91 k |
| DIVIDE | 118976 / 37771 | - | - | - | 78.75 M | - | - | - | 44.52 M |
| ISQRT | 118976 / 9442 | - | - | - | 3.03 M | - | - | - | 11.64 M |
| MODMUL_R2 | 302144 / 302171 | - | - | - | 844.36 M | - | - | - | 8.50 M |

#### p1024 (1024-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 151040 / 151085 | - | - | - | 574.34 M | - | - | - | 46.17 M |
| SUBTRACT | 151040 / 151085 | - | - | - | 577.75 M | - | - | - | 53.23 M |
| ADDMOD | 151040 / 151085 | - | - | - | 226.45 M | - | - | - | 18.02 M |
| SUBTRACTMOD | 151040 / 151085 | - | - | - | 583.36 M | - | - | - | 31.89 M |
| MULTIPLYOPERANDSCANNING | 151040 / 151085 | - | - | - | 279.12 M | - | - | - | 8.66 M |
| MULTIPLYPRODUCTSCANNING | 151040 / 151085 | - | - | - | 158.35 M | - | - | - | 8.66 M |
| MONTGOMERYMULTIPLICATION | 151040 / 151085 | - | - | - | 375.99 M | - | - | - | 1.30 M |
| COMPARE | 151040 / 151085 | - | - | - | 616.23 M | - | - | - | 124.84 M |
| REDUCE | 118976 / 18885 | - | - | - | 12.05 M | - | - | - | 70.37 M |
| MODMUL | 118976 / 9442 | - | - | - | 3.71 M | - | - | - | 2.96 M |
| MODEXP | 118976 / 2360 | - | - | - | 195.55 k | - | - | - | 4.30 k |
| EXPONENTIATION | 118976 / 2360 | - | - | - | 189.40 k | - | - | - | 32.35 k |
| DIVIDE | 118976 / 18885 | - | - | - | 10.94 M | - | - | - | 39.85 M |
| ISQRT | 118976 / 4721 | - | - | - | 745.72 k | - | - | - | 5.96 M |
| MODMUL_R2 | 151040 / 151085 | - | - | - | 267.58 M | - | - | - | 2.96 M |

#### p2048 (2048-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 118976 / 75542 | - | - | - | 372.18 M | - | - | - | 32.84 M |
| SUBTRACT | 118976 / 75542 | - | - | - | 380.95 M | - | - | - | 28.65 M |
| ADDMOD | 118976 / 75542 | - | - | - | 385.82 M | - | - | - | 19.21 M |
| SUBTRACTMOD | 118976 / 75542 | - | - | - | 397.24 M | - | - | - | 21.19 M |
| MULTIPLYOPERANDSCANNING | 118976 / 75542 | - | - | - | 104.64 M | - | - | - | 2.68 M |
| MULTIPLYPRODUCTSCANNING | 118976 / 75542 | - | - | - | 19.09 M | - | - | - | 2.64 M |
| MONTGOMERYMULTIPLICATION | 118976 / 75542 | - | - | - | 137.73 M | - | - | - | 394.65 k |
| COMPARE | 118976 / 75542 | - | - | - | 434.50 M | - | - | - | 200.63 M |
| REDUCE | 118976 / 9442 | - | - | - | 2.73 M | - | - | - | 49.51 M |
| MODMUL | 118976 / 4721 | - | - | - | 812.88 k | - | - | - | 903.62 k |
| MODEXP | 118976 / 1180 | - | - | - | 16.59 k | - | - | - | 572.7 |
| EXPONENTIATION | 118976 / 1180 | - | - | - | 39.49 k | - | - | - | 5.46 k |
| DIVIDE | 118976 / 9442 | - | - | - | 5.68 M | - | - | - | 24.14 M |
| ISQRT | 118976 / 2360 | - | - | - | 332.08 k | - | - | - | 3.52 M |
| MODMUL_R2 | 118976 / 75542 | - | - | - | 34.24 M | - | - | - | 903.98 k |

## 5. Head to head

Best OpenCL GPU result against best OpenCL CPU result and the CPU libraries.
Ratios above 1.00x mean the GPU is faster than that baseline.

### secp256k1 (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w32-opt | 2.72 G | none | n/a | 85.01 M | n/a |
| SUBTRACT | w32-opt | 2.61 G | none | n/a | 95.14 M | n/a |
| ADDMOD | w32-opt | 2.83 G | none | n/a | 34.96 M | n/a |
| SUBTRACTMOD | w32-opt | 2.84 G | none | n/a | 44.14 M | n/a |
| MULTIPLYOPERANDSCANNING | w32-opt | 2.05 G | none | n/a | 66.48 M | n/a |
| MULTIPLYPRODUCTSCANNING | w32-opt | 2.20 G | none | n/a | 67.26 M | n/a |
| MONTGOMERYMULTIPLICATION | w32-opt | 2.75 G | none | n/a | 9.65 M | n/a |
| COMPARE | w32-opt | 2.85 G | none | n/a | 115.47 M | n/a |
| REDUCE | w32-opt | 351.39 M | none | n/a | 88.54 M | n/a |
| MODMUL | w32-opt | 131.08 M | none | n/a | 16.84 M | n/a |
| MODEXP | w32-opt | 901.41 k | none | n/a | 156.43 k | n/a |
| EXPONENTIATION | w32-opt | 1.03 M | none | n/a | 499.45 k | n/a |
| DIVIDE | w32-opt | 297.81 M | none | n/a | 48.71 M | n/a |
| ISQRT | w32-opt | 10.81 M | none | n/a | 21.45 M | n/a |
| MODMUL_R2 | w32-opt | 2.39 G | none | n/a | 16.78 M | n/a |

### rsa256(composite) (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w32-opt | 2.61 G | none | n/a | 83.93 M | n/a |
| SUBTRACT | w32-opt | 2.58 G | none | n/a | 95.61 M | n/a |
| ADDMOD | w32-opt | 2.79 G | none | n/a | 40.47 M | n/a |
| SUBTRACTMOD | w32-opt | 2.84 G | none | n/a | 44.27 M | n/a |
| MULTIPLYOPERANDSCANNING | w32-opt | 2.07 G | none | n/a | 67.05 M | n/a |
| MULTIPLYPRODUCTSCANNING | w32-opt | 2.19 G | none | n/a | 65.64 M | n/a |
| MONTGOMERYMULTIPLICATION | w32-opt | 2.73 G | none | n/a | 9.69 M | n/a |
| COMPARE | w32-opt | 3.02 G | none | n/a | 120.30 M | n/a |
| REDUCE | w32-opt | 352.13 M | none | n/a | 55.83 M | n/a |
| MODMUL | w32-opt | 132.45 M | none | n/a | 16.80 M | n/a |
| MODEXP | w32-opt | 967.75 k | none | n/a | 166.94 k | n/a |
| EXPONENTIATION | w32-opt | 1.09 M | none | n/a | 494.94 k | n/a |
| DIVIDE | w32-opt | 287.20 M | none | n/a | 46.33 M | n/a |
| ISQRT | w32-opt | 10.78 M | none | n/a | 21.06 M | n/a |
| MODMUL_R2 | w32-opt | 2.35 G | none | n/a | 16.55 M | n/a |

### brainpoolP512r1 (512-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w32-opt | 1.23 G | none | n/a | 65.22 M | n/a |
| SUBTRACT | w32-opt | 1.23 G | none | n/a | 73.34 M | n/a |
| ADDMOD | w32-opt | 1.24 G | none | n/a | 36.14 M | n/a |
| SUBTRACTMOD | w32-opt | 1.27 G | none | n/a | 38.27 M | n/a |
| MULTIPLYOPERANDSCANNING | w32-opt | 715.65 M | none | n/a | 31.77 M | n/a |
| MULTIPLYPRODUCTSCANNING | w32-opt | 769.01 M | none | n/a | 31.68 M | n/a |
| MONTGOMERYMULTIPLICATION | w32-opt | 1.12 G | none | n/a | 4.10 M | n/a |
| COMPARE | w32-opt | 1.34 G | none | n/a | 85.63 M | n/a |
| REDUCE | w32-opt | 14.36 M | none | n/a | 51.12 M | n/a |
| MODMUL | w32-opt | 2.56 M | none | n/a | 8.40 M | n/a |
| MODEXP | w32-opt | 58.76 k | none | n/a | 27.22 k | n/a |
| EXPONENTIATION | w32-opt | 71.75 k | none | n/a | 157.91 k | n/a |
| DIVIDE | w32-opt | 25.00 M | none | n/a | 44.52 M | n/a |
| ISQRT | w32-opt | 240.25 k | none | n/a | 11.64 M | n/a |
| MODMUL_R2 | w32-opt | 844.43 M | none | n/a | 8.50 M | n/a |

### p1024 (1024-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w32-opt | 574.51 M | none | n/a | 46.17 M | n/a |
| SUBTRACT | w32-opt | 577.93 M | none | n/a | 53.23 M | n/a |
| ADDMOD | w32-opt | 226.52 M | none | n/a | 18.02 M | n/a |
| SUBTRACTMOD | w32-opt | 583.54 M | none | n/a | 31.89 M | n/a |
| MULTIPLYOPERANDSCANNING | w32-opt | 279.20 M | none | n/a | 8.66 M | n/a |
| MULTIPLYPRODUCTSCANNING | w32-opt | 158.40 M | none | n/a | 8.66 M | n/a |
| MONTGOMERYMULTIPLICATION | w32-opt | 376.10 M | none | n/a | 1.30 M | n/a |
| COMPARE | w32-opt | 616.41 M | none | n/a | 124.84 M | n/a |
| REDUCE | w32-opt | 1.91 M | none | n/a | 70.37 M | n/a |
| MODMUL | w32-opt | 294.61 k | none | n/a | 2.96 M | n/a |
| MODEXP | w32-opt | 3.88 k | none | n/a | 4.30 k | n/a |
| EXPONENTIATION | w32-opt | 3.76 k | none | n/a | 32.35 k | n/a |
| DIVIDE | w32-opt | 1.74 M | none | n/a | 39.85 M | n/a |
| ISQRT | w32-opt | 29.59 k | none | n/a | 5.96 M | n/a |
| MODMUL_R2 | w32-opt | 267.66 M | none | n/a | 2.96 M | n/a |

### p2048 (2048-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w32-opt | 236.31 M | none | n/a | 32.84 M | n/a |
| SUBTRACT | w32-opt | 241.88 M | none | n/a | 28.65 M | n/a |
| ADDMOD | w32-opt | 244.97 M | none | n/a | 19.21 M | n/a |
| SUBTRACTMOD | w32-opt | 252.22 M | none | n/a | 21.19 M | n/a |
| MULTIPLYOPERANDSCANNING | w32-opt | 66.44 M | none | n/a | 2.68 M | n/a |
| MULTIPLYPRODUCTSCANNING | w32-opt | 12.12 M | none | n/a | 2.64 M | n/a |
| MONTGOMERYMULTIPLICATION | w32-opt | 87.45 M | none | n/a | 394.65 k | n/a |
| COMPARE | w32-opt | 275.88 M | none | n/a | 200.63 M | n/a |
| REDUCE | w32-opt | 216.82 k | none | n/a | 49.51 M | n/a |
| MODMUL | w32-opt | 32.26 k | none | n/a | 903.62 k | n/a |
| MODEXP | w32-opt | 164.6 | none | n/a | 572.7 | n/a |
| EXPONENTIATION | w32-opt | 391.7 | none | n/a | 5.46 k | n/a |
| DIVIDE | w32-opt | 450.63 k | none | n/a | 24.14 M | n/a |
| ISQRT | w32-opt | 6.59 k | none | n/a | 3.52 M | n/a |
| MODMUL_R2 | w32-opt | 21.74 M | none | n/a | 903.98 k | n/a |

## 6. Raw data

Also written to `NVIDIA_GeForce_RTX_5090_Report.csv` for analysis.

```csv
kind,device,device_type,kernel,modulus,bits,operation,items,seconds,ops_per_sec,mismatches
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,ADD,604343,0.007108913,85012011.954,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,ADD,604343,0.010467887,57733046.253,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,ADD,604343,0.009356954,64587577.462,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,secp256k1,256,ADD,700000,0.000095840,7303839732.888,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,ADD,604288,0.000222282,2718564479.682,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,ADD,604288,0.004270685,141496734.235,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,SUBTRACT,604343,0.006352136,95140120.781,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,SUBTRACT,604343,0.008633075,70003213.983,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,SUBTRACT,604343,0.008646667,69893175.720,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,secp256k1,256,SUBTRACT,700000,0.000096128,7281957390.146,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,SUBTRACT,604288,0.000231529,2609989649.496,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,SUBTRACT,604288,0.004473741,135074428.535,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,ADDMOD,604343,0.017285767,34961885.093,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,ADDMOD,604343,0.008957983,67464182.547,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,ADDMOD,604343,0.007745310,78026958.832,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,secp256k1,256,ADDMOD,700000,0.000090880,7702464788.732,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,ADDMOD,604288,0.000213415,2831514347.270,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,ADDMOD,604288,0.004379722,137974052.268,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,SUBTRACTMOD,604343,0.013691766,44139156.483,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,SUBTRACTMOD,604343,0.009133254,66169516.761,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,SUBTRACTMOD,604343,0.010555022,57256441.908,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,secp256k1,256,SUBTRACTMOD,700000,0.000096768,7233796296.296,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,SUBTRACTMOD,604288,0.000212533,2843264463.500,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,SUBTRACTMOD,604288,0.004351959,138854246.185,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,MULTIPLYOPERANDSCANNING,604343,0.009090373,66481649.986,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,604343,0.009187833,65776446.298,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,604343,0.007959948,75922984.167,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,MULTIPLYOPERANDSCANNING,604288,0.000295070,2047947944.605,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,MULTIPLYOPERANDSCANNING,604288,0.005822831,103779071.231,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,MULTIPLYPRODUCTSCANNING,604343,0.008984923,67261899.881,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,604343,0.008416254,71806650.732,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,604343,0.010268249,58855506.435,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,secp256k1,256,MULTIPLYPRODUCTSCANNING,700000,0.000094528,7405213270.142,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,MULTIPLYPRODUCTSCANNING,604288,0.000275263,2195310618.404,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,MULTIPLYPRODUCTSCANNING,604288,0.005861875,103087836.236,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,MONTGOMERYMULTIPLICATION,604343,0.062626562,9649946.922,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,604343,0.011996906,50374904.921,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,604343,0.007649728,79001895.326,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,secp256k1,256,MONTGOMERYMULTIPLICATION,700000,0.000096992,7217090069.284,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,MONTGOMERYMULTIPLICATION,604288,0.000219717,2750302484.182,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,MONTGOMERYMULTIPLICATION,604288,0.004387116,137741519.452,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,COMPARE,604343,0.005233663,115472282.296,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,COMPARE,604343,0.008780690,68826368.565,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,COMPARE,604343,0.009047838,66794187.657,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,secp256k1,256,COMPARE,700000,0.000091456,7653953813.856,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,COMPARE,604288,0.000212093,2849166784.530,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,COMPARE,604288,0.004288098,140922149.923,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,REDUCE,75542,0.000853219,88537636.695,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,REDUCE,75542,0.006645608,11367206.741,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,REDUCE,75542,0.007861622,9608958.383,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,secp256k1,256,REDUCE,700000,0.000092704,7550914739.386,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,REDUCE,118976,0.000214978,553433661.903,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,REDUCE,118976,0.001074619,110714575.209,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,MODMUL,37771,0.002243558,16835312.895,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,MODMUL,37771,0.008002588,4719848.094,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,MODMUL,37771,0.005698611,6628106.360,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,secp256k1,256,MODMUL,700000,0.000221216,3164328077.535,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,MODMUL,118976,0.000288147,412900340.348,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,MODMUL,118976,0.001192944,99733106.349,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,MODEXP,9442,0.060359780,156428.668,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,MODEXP,9442,0.011788545,800947.025,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,MODEXP,9442,0.014472778,652397.213,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,secp256k1,256,MODEXP,700000,0.063936286,10948399.474,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,MODEXP,118976,0.010474660,11358459.402,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,MODEXP,118976,0.011464729,10377567.605,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,EXPONENTIATION,9442,0.018904690,499452.777,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,EXPONENTIATION,9442,0.008693090,1086150.012,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,EXPONENTIATION,9442,0.024938842,378606.195,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,EXPONENTIATION,118976,0.009142381,13013677.530,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,EXPONENTIATION,118976,0.009494740,12530727.473,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,DIVIDE,75542,0.001550924,48707737.063,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,DIVIDE,75542,0.007026467,10751064.501,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,DIVIDE,75542,0.007668083,9851484.306,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,secp256k1,256,DIVIDE,700000,0.000090560,7729681978.799,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,DIVIDE,118976,0.000253662,469033593.407,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,DIVIDE,118976,0.002057365,57829313.352,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,ISQRT,18885,0.000880371,21451181.682,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,ISQRT,18885,0.007722036,2445598.550,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,ISQRT,118976,0.001746996,68103186.856,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,ISQRT,118976,0.001859970,63966625.692,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,MODMUL_R2,604343,0.036021306,16777376.170,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,MODMUL_R2,604343,0.007896367,76534310.069,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,MODMUL_R2,604343,0.011926111,50673936.929,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,secp256k1,256,MODMUL_R2,700000,0.000106624,6565126050.420,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,MODMUL_R2,604288,0.000252399,2394178468.020,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,MODMUL_R2,604288,0.004286145,140986361.147,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,ADD,604343,0.007200978,83925127.973,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,ADD,604343,0.008847126,68309527.563,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,ADD,604343,0.008877586,68075149.333,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,rsa256(composite),256,ADD,700000,0.000085824,8156226696.495,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,ADD,604288,0.000231319,2612356616.520,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,ADD,604288,0.004267731,141594679.330,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACT,604343,0.006320856,95610941.874,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACT,604343,0.008856590,68236534.199,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACT,604343,0.009543734,63323538.362,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,rsa256(composite),256,SUBTRACT,700000,0.000096096,7284382284.382,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,SUBTRACT,604288,0.000234234,2579845886.362,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,SUBTRACT,604288,0.004435769,136230717.439,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,ADDMOD,604343,0.014933202,40469753.374,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,ADDMOD,604343,0.009242998,65383872.412,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,ADDMOD,604343,0.010709437,56430883.875,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,rsa256(composite),256,ADDMOD,700000,0.000090208,7759843916.282,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,ADDMOD,604288,0.000216932,2785609775.259,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,ADDMOD,604288,0.004405521,137166074.047,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACTMOD,604343,0.013649966,44274323.159,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACTMOD,604343,0.008193502,73758815.489,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACTMOD,604343,0.010784126,56040054.274,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,rsa256(composite),256,SUBTRACTMOD,700000,0.000092896,7535308301.757,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,SUBTRACTMOD,604288,0.000212764,2840177932.476,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,SUBTRACTMOD,604288,0.004392316,137578440.701,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYOPERANDSCANNING,604343,0.009013076,67051804.539,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,604343,0.009635459,62720726.598,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,604343,0.008024509,75312146.216,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,604288,0.000292364,2066902603.154,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,604288,0.005887193,102644504.701,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,604343,0.009206925,65640047.737,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,604343,0.009972001,60603985.339,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,604343,0.009504169,63587147.961,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,700000,0.000093920,7453151618.399,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,604288,0.000275584,2192754791.526,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,604288,0.005836898,103528963.584,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MONTGOMERYMULTIPLICATION,604343,0.062369513,9689718.138,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,604343,0.011677053,51754753.273,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,604343,0.010002640,60418349.097,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,rsa256(composite),256,MONTGOMERYMULTIPLICATION,700000,0.000089792,7795794725.588,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,604288,0.000221721,2725441645.479,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,604288,0.004465765,135315671.164,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,COMPARE,604343,0.005023706,120298241.834,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,COMPARE,604343,0.008130791,74327701.639,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,COMPARE,604343,0.009442502,64002422.685,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,rsa256(composite),256,COMPARE,700000,0.000090432,7740622788.393,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,COMPARE,604288,0.000200110,3019778789.813,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,COMPARE,604288,0.004325470,139704588.443,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,REDUCE,75542,0.001353169,55825996.466,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,REDUCE,75542,0.008160094,9257491.360,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,REDUCE,75542,0.007034421,10738908.057,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,rsa256(composite),256,REDUCE,700000,0.000084288,8304859529.233,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,REDUCE,118976,0.000214527,554597127.598,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,REDUCE,118976,0.001096080,108546815.723,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MODMUL,37771,0.002248858,16795636.809,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MODMUL,37771,0.007295923,5177000.856,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MODMUL,37771,0.007752814,4871908.384,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,rsa256(composite),256,MODMUL,700000,0.000218752,3199970743.125,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,MODMUL,118976,0.000285162,417222649.161,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,MODMUL,118976,0.002001689,59437810.248,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MODEXP,9442,0.056559887,166938.099,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MODEXP,9442,0.010666270,885220.435,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MODEXP,9442,0.014379401,656633.750,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,rsa256(composite),256,MODEXP,700000,0.062912092,11126636.832,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,MODEXP,118976,0.009756616,12194392.176,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,MODEXP,118976,0.010551726,11275501.300,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,EXPONENTIATION,9442,0.019077096,494939.059,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,EXPONENTIATION,9442,0.007949974,1187676.852,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,EXPONENTIATION,9442,0.024465985,385923.557,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,EXPONENTIATION,118976,0.008688791,13693043.993,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,EXPONENTIATION,118976,0.009472148,12560614.515,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,DIVIDE,75542,0.001630364,46334438.600,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,DIVIDE,75542,0.006714534,11250520.116,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,DIVIDE,75542,0.006702141,11271323.756,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,rsa256(composite),256,DIVIDE,700000,0.000089856,7790242165.242,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,DIVIDE,118976,0.000263030,452328330.565,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,DIVIDE,118976,0.002037587,58390632.439,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,ISQRT,18885,0.000896892,21056045.833,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,ISQRT,18885,0.006851428,2756359.620,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,ISQRT,118976,0.001751063,67945007.812,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,ISQRT,118976,0.001855722,64113045.789,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MODMUL_R2,604343,0.036509601,16552988.380,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MODMUL_R2,604343,0.008960716,67443604.612,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MODMUL_R2,604343,0.013150950,45954322.823,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,rsa256(composite),256,MODMUL_R2,700000,0.000101280,6911532385.466,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,MODMUL_R2,604288,0.000257518,2346586304.873,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,MODMUL_R2,604288,0.004388589,137695283.699,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,ADD,302171,0.004632932,65222413.463,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,ADD,302171,0.008209811,36806084.558,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,ADD,302171,0.008418583,35893332.309,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,brainpoolP512r1,512,ADD,700000,0.000146560,4776200873.362,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,ADD,302144,0.000245215,1232158837.934,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,ADD,302144,0.004451818,67869799.605,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACT,302171,0.004120070,73341229.255,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACT,302171,0.007889148,38302107.770,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACT,302171,0.008022972,37663224.346,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,brainpoolP512r1,512,SUBTRACT,700000,0.000148640,4709364908.504,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,SUBTRACT,302144,0.000244855,1233971380.698,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,SUBTRACT,302144,0.004460796,67733202.541,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,ADDMOD,302171,0.008361960,36136384.021,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,ADDMOD,302171,0.007927117,38118649.377,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,ADDMOD,302171,0.010191522,29649251.391,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,brainpoolP512r1,512,ADDMOD,700000,0.000147872,4733823847.652,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,ADDMOD,302144,0.000244424,1236147124.893,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,ADDMOD,302144,0.004427965,68235407.873,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACTMOD,302171,0.007894844,38274474.710,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACTMOD,302171,0.007547434,40036256.276,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACTMOD,302171,0.009927467,30437875.108,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,brainpoolP512r1,512,SUBTRACTMOD,700000,0.000142688,4905808477.237,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,SUBTRACTMOD,302144,0.000237791,1270629115.995,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,SUBTRACTMOD,302144,0.004347211,69502949.135,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,302171,0.009511903,31767670.704,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,302171,0.008576678,35231706.758,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,302171,0.008783210,34403254.122,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,302144,0.000422231,715589071.929,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,302144,0.005935124,50907783.167,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,302171,0.009537953,31680906.180,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,302171,0.008635840,34990342.925,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,302171,0.008948363,33768298.803,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,700000,0.000153184,4569667850.428,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,302144,0.000392935,768941416.332,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,302144,0.005975450,50564223.966,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,302171,0.073623691,4104263.124,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,302171,0.011014063,27435016.709,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,302171,0.008654421,34915217.688,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,700000,0.002374688,294775566.306,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,302144,0.000268610,1124842870.006,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,302144,0.004517515,66882785.623,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,COMPARE,302171,0.003528776,85630544.491,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,COMPARE,302171,0.007261331,41613720.908,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,COMPARE,302171,0.008443254,35788452.929,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,brainpoolP512r1,512,COMPARE,700000,0.000247200,2831715210.356,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,COMPARE,302144,0.000225748,1338411745.571,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,COMPARE,302144,0.004414328,68446206.597,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,REDUCE,37771,0.000738842,51121894.398,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,REDUCE,37771,0.006786380,5565706.637,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,REDUCE,37771,0.006935488,5446047.946,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,brainpoolP512r1,512,REDUCE,700000,0.000243744,2871865563.870,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,REDUCE,118976,0.002630312,45232656.883,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,REDUCE,118976,0.003289303,36170580.337,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL,18885,0.002248287,8399728.446,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL,18885,0.007160014,2637564.672,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL,18885,0.006770094,2789473.759,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,brainpoolP512r1,512,MODMUL,700000,0.002380640,294038577.861,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,MODMUL,118976,0.007386959,16106221.952,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,MODMUL,118976,0.006618592,17976028.856,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MODEXP,4721,0.173411150,27224.316,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MODEXP,4721,0.018221722,259086.381,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MODEXP,4721,0.019557403,241391.969,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,brainpoolP512r1,512,MODEXP,700000,0.211561188,3308735.438,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,MODEXP,118976,0.080342464,1480860.729,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,MODEXP,118976,0.082368989,1444427.101,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,EXPONENTIATION,4721,0.029895960,157914.314,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,EXPONENTIATION,4721,0.008171459,577742.607,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,EXPONENTIATION,4721,0.031465514,150037.275,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,EXPONENTIATION,118976,0.065801566,1808102.865,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,EXPONENTIATION,118976,0.069387262,1714666.303,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,DIVIDE,37771,0.000848491,44515499.863,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,DIVIDE,37771,0.006510096,5801911.456,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,DIVIDE,37771,0.006898984,5474864.061,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,brainpoolP512r1,512,DIVIDE,700000,0.002377216,294462093.474,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,DIVIDE,118976,0.001510787,78751005.622,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,DIVIDE,118976,0.003818999,31153712.952,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,ISQRT,9442,0.000811390,11636822.482,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,ISQRT,9442,0.006834266,1381567.529,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,ISQRT,118976,0.039300620,3027331.386,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,ISQRT,118976,0.040279568,2953755.616,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL_R2,302171,0.035549359,8500040.763,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL_R2,302171,0.010484665,28820281.994,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL_R2,302171,0.013554185,22293557.418,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,brainpoolP512r1,512,MODMUL_R2,700000,0.000170240,4111842105.263,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,MODMUL_R2,302144,0.000357839,844357269.947,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,MODMUL_R2,302144,0.004470445,67587007.816,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p1024,1024,ADD,151085,0.003272570,46167076.818,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p1024,1024,ADD,151085,0.007626264,19811142.411,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p1024,1024,ADD,151085,0.008074570,18711213.148,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p1024,1024,ADD,700000,0.000460992,1518464528.669,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,ADD,151040,0.000262979,574342033.496,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,ADD,151040,0.004347662,34740509.440,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p1024,1024,SUBTRACT,151085,0.002838317,53230491.365,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p1024,1024,SUBTRACT,151085,0.007183043,21033564.588,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p1024,1024,SUBTRACT,151085,0.007480433,20197360.293,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p1024,1024,SUBTRACT,700000,0.004714080,148491328.106,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,SUBTRACT,151040,0.000261426,577754378.623,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,SUBTRACT,151040,0.004433655,34066701.581,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p1024,1024,ADDMOD,151085,0.008382809,18023195.102,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p1024,1024,ADDMOD,151085,0.006513954,23194054.202,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p1024,1024,ADDMOD,151085,0.009228486,16371591.412,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p1024,1024,ADDMOD,700000,0.000246464,2840171384.056,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,ADDMOD,151040,0.000666976,226454933.391,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,ADDMOD,151040,0.004383950,34452947.460,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p1024,1024,SUBTRACTMOD,151085,0.004737481,31891420.716,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p1024,1024,SUBTRACTMOD,151085,0.007828548,19299236.380,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p1024,1024,SUBTRACTMOD,151085,0.007424540,20349409.164,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p1024,1024,SUBTRACTMOD,700000,0.000246464,2840171384.056,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,SUBTRACTMOD,151040,0.000258912,583364007.618,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,SUBTRACTMOD,151040,0.004324618,34925627.883,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p1024,1024,MULTIPLYOPERANDSCANNING,151085,0.017454758,8655806.029,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p1024,1024,MULTIPLYOPERANDSCANNING,151085,0.009346745,16164450.873,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p1024,1024,MULTIPLYOPERANDSCANNING,151085,0.008485034,17806056.885,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,MULTIPLYOPERANDSCANNING,151040,0.000541137,279115975.532,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,MULTIPLYOPERANDSCANNING,151040,0.006138590,24604998.209,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p1024,1024,MULTIPLYPRODUCTSCANNING,151085,0.017449558,8658385.495,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,151085,0.008717231,17331765.125,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,151085,0.008123738,18597965.633,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p1024,1024,MULTIPLYPRODUCTSCANNING,700000,0.000449600,1556939501.779,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,MULTIPLYPRODUCTSCANNING,151040,0.000953830,158351055.550,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,MULTIPLYPRODUCTSCANNING,151040,0.006454880,23399350.991,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p1024,1024,MONTGOMERYMULTIPLICATION,151085,0.116587743,1295890.939,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p1024,1024,MONTGOMERYMULTIPLICATION,151085,0.014911641,10132016.902,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p1024,1024,MONTGOMERYMULTIPLICATION,151085,0.008146972,18544927.036,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p1024,1024,MONTGOMERYMULTIPLICATION,700000,0.008761887,79891466.302,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,MONTGOMERYMULTIPLICATION,151040,0.000401712,375990738.282,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,MONTGOMERYMULTIPLICATION,151040,0.004494961,33602070.413,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p1024,1024,COMPARE,151085,0.001210197,124843299.507,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p1024,1024,COMPARE,151085,0.007684411,19661233.593,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p1024,1024,COMPARE,151085,0.007230142,20896546.327,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p1024,1024,COMPARE,700000,0.000458752,1525878906.250,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,COMPARE,151040,0.000245105,616225675.241,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,COMPARE,151040,0.004333475,34854244.032,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p1024,1024,REDUCE,18885,0.000268369,70369533.300,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p1024,1024,REDUCE,18885,0.007480458,2524577.999,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p1024,1024,REDUCE,18885,0.007451097,2534526.119,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p1024,1024,REDUCE,700000,0.000454048,1541687222.496,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,REDUCE,118976,0.009870704,12053446.201,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,REDUCE,118976,0.011281161,10546432.354,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p1024,1024,MODMUL,9442,0.003190316,2959581.624,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p1024,1024,MODMUL,9442,0.006439150,1466342.590,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p1024,1024,MODMUL,9442,0.007825398,1206583.979,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p1024,1024,MODMUL,700000,0.004740128,147675337.037,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,MODMUL,118976,0.032049057,3712308.906,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,MODMUL,118976,0.037574704,3166385.572,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p1024,1024,MODEXP,2360,0.548246237,4304.635,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p1024,1024,MODEXP,2360,0.045490842,51878.574,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p1024,1024,MODEXP,2360,0.039022155,60478.464,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p1024,1024,MODEXP,700000,1.212423444,577356.041,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,MODEXP,118976,0.608427604,195546.683,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,MODEXP,118976,0.610097935,195011.314,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p1024,1024,EXPONENTIATION,2360,0.072960551,32346.247,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p1024,1024,EXPONENTIATION,2360,0.012956760,182144.303,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p1024,1024,EXPONENTIATION,2360,0.087323602,27025.912,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,EXPONENTIATION,118976,0.628186918,189395.858,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,EXPONENTIATION,118976,0.650650762,182856.929,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p1024,1024,DIVIDE,18885,0.000473959,39845226.600,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p1024,1024,DIVIDE,18885,0.006528020,2892913.870,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p1024,1024,DIVIDE,18885,0.006222460,3034973.329,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p1024,1024,DIVIDE,700000,0.004723648,148190551.032,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,DIVIDE,118976,0.010871322,10944023.290,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,DIVIDE,118976,0.014599209,8149482.646,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p1024,1024,ISQRT,4721,0.000791512,5964533.973,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p1024,1024,ISQRT,4721,0.003599581,1311541.510,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,ISQRT,118976,0.159545113,745720.114,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,ISQRT,118976,0.164124604,724912.641,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p1024,1024,MODMUL_R2,151085,0.050957772,2964905.925,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p1024,1024,MODMUL_R2,151085,0.010304142,14662550.232,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p1024,1024,MODMUL_R2,151085,0.014260435,10594697.839,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p1024,1024,MODMUL_R2,700000,0.000556064,1258847902.400,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,MODMUL_R2,151040,0.000564472,267577464.008,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,MODMUL_R2,151040,0.004668190,32355153.255,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p2048,2048,ADD,75542,0.002300086,32843120.271,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p2048,2048,ADD,75542,0.007151363,10563300.737,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p2048,2048,ADD,75542,0.007403446,10203626.620,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p2048,2048,ADD,700000,0.000459456,1523540883.131,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,ADD,118976,0.000319677,372175735.366,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,ADD,118976,0.006888875,17270744.399,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p2048,2048,SUBTRACT,75542,0.002636534,28652011.327,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p2048,2048,SUBTRACT,75542,0.007078815,10671560.113,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p2048,2048,SUBTRACT,75542,0.006679983,11308711.620,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p2048,2048,SUBTRACT,700000,0.000459776,1522480512.249,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,SUBTRACT,118976,0.000312312,380952331.578,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,SUBTRACT,118976,0.006880128,17292701.386,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p2048,2048,ADDMOD,75542,0.003932313,19210576.270,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p2048,2048,ADDMOD,75542,0.007872546,9595624.988,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p2048,2048,ADDMOD,75542,0.008382169,9012225.764,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p2048,2048,ADDMOD,700000,0.000464160,1508100654.947,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,ADDMOD,118976,0.000308375,385815836.789,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,ADDMOD,118976,0.006915985,17203045.059,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p2048,2048,SUBTRACTMOD,75542,0.003565777,21185285.271,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p2048,2048,SUBTRACTMOD,75542,0.008196064,9216863.076,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p2048,2048,SUBTRACTMOD,75542,0.006895683,10954969.931,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p2048,2048,SUBTRACTMOD,700000,0.009610528,72836788.988,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,SUBTRACTMOD,118976,0.000299508,397238157.250,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,SUBTRACTMOD,118976,0.007011226,16969357.639,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p2048,2048,MULTIPLYOPERANDSCANNING,75542,0.028222352,2676672.728,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p2048,2048,MULTIPLYOPERANDSCANNING,75542,0.009204408,8207154.522,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p2048,2048,MULTIPLYOPERANDSCANNING,75542,0.011024579,6852143.779,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,MULTIPLYOPERANDSCANNING,118976,0.001136958,104644134.645,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,MULTIPLYOPERANDSCANNING,118976,0.011128141,10691453.240,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p2048,2048,MULTIPLYPRODUCTSCANNING,75542,0.028587306,2642501.543,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,75542,0.008222330,9187420.227,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,75542,0.009255847,8161543.758,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p2048,2048,MULTIPLYPRODUCTSCANNING,700000,0.009617024,72787590.007,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,MULTIPLYPRODUCTSCANNING,118976,0.006233571,19086330.796,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,MULTIPLYPRODUCTSCANNING,118976,0.014007426,8493780.310,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p2048,2048,MONTGOMERYMULTIPLICATION,75542,0.191416944,394646.359,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p2048,2048,MONTGOMERYMULTIPLICATION,75542,0.021301047,3546398.456,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p2048,2048,MONTGOMERYMULTIPLICATION,75542,0.008204411,9207486.313,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p2048,2048,MONTGOMERYMULTIPLICATION,700000,0.000977280,716273739.358,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,MONTGOMERYMULTIPLICATION,118976,0.000863849,137727781.744,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,MONTGOMERYMULTIPLICATION,118976,0.007501067,15861209.874,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p2048,2048,COMPARE,75542,0.000376515,200634844.730,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p2048,2048,COMPARE,75542,0.006816237,11082654.549,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p2048,2048,COMPARE,75542,0.006597271,11450492.257,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p2048,2048,COMPARE,700000,0.000457248,1530897893.484,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,COMPARE,118976,0.000273820,434504398.638,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,COMPARE,118976,0.006913491,17209250.577,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p2048,2048,REDUCE,9442,0.000190702,49511794.888,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p2048,2048,REDUCE,9442,0.006657746,1418197.707,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p2048,2048,REDUCE,9442,0.007089034,1331916.274,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p2048,2048,REDUCE,700000,0.000461600,1516464471.404,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,REDUCE,118976,0.043547289,2732110.374,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,REDUCE,118976,0.049952365,2381789.131,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p2048,2048,MODMUL,4721,0.005224515,903624.549,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p2048,2048,MODMUL,4721,0.006923720,681858.884,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p2048,2048,MODMUL,4721,0.007828839,603026.825,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p2048,2048,MODMUL,700000,0.009618752,72774513.783,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,MODMUL,118976,0.146363295,812881.399,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,MODMUL,118976,0.154938621,767891.177,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p2048,2048,MODEXP,1180,2.060571857,572.657,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p2048,2048,MODEXP,1180,0.154817005,7621.902,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p2048,2048,MODEXP,1180,0.125006579,9439.503,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p2048,2048,MODEXP,700000,6.987298489,100181.780,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,MODEXP,118976,7.169845908,16593.941,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,MODEXP,118976,7.182643102,16564.376,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p2048,2048,EXPONENTIATION,1180,0.215960840,5463.954,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p2048,2048,EXPONENTIATION,1180,0.021940471,53781.890,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p2048,2048,EXPONENTIATION,1180,0.195477184,6036.510,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,EXPONENTIATION,118976,3.012633087,39492.363,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,EXPONENTIATION,118976,3.068739812,38770.312,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p2048,2048,DIVIDE,9442,0.000391202,24135865.497,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p2048,2048,DIVIDE,9442,0.007017849,1345426.495,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p2048,2048,DIVIDE,9442,0.007317903,1290260.341,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p2048,2048,DIVIDE,700000,0.000884864,791082019.384,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,DIVIDE,118976,0.020952796,5678287.555,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,DIVIDE,118976,0.030086872,3954415.711,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p2048,2048,ISQRT,2360,0.000670953,3517385.208,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p2048,2048,ISQRT,2360,0.005883947,401091.307,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,ISQRT,118976,0.358277460,332077.826,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,ISQRT,118976,0.365727409,325313.326,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p2048,2048,MODMUL_R2,75542,0.083565991,903980.186,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p2048,2048,MODMUL_R2,75542,0.012996842,5812334.959,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p2048,2048,MODMUL_R2,75542,0.019142070,3946386.136,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p2048,2048,MODMUL_R2,700000,0.001894976,369397818.231,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,MODMUL_R2,118976,0.003475115,34236565.464,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,MODMUL_R2,118976,0.010010749,11884824.881,0
```
