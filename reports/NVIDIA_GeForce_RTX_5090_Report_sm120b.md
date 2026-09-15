# MPA-OpenCL benchmark report - NVIDIA GeForce RTX 5090


> **Note.** Two column groups have been removed from this report: the CGBN
> reference column, which was invalid, and the multi-threaded GMP and OpenSSL
> baselines, which predate the 2026-09-12 timing fix and were understated.
> See `reports/README.md`. The single-threaded GMP column, the OpenCL-on-CPU
> rows and every MPA measurement are unaffected, and every configuration was
> verified word-for-word against GMP before it was timed.


> **Not a vendor-runtime result.** The OpenCL rows were produced by a third-party runtime (`OpenCL 3.0 PoCL HSTR: CUDA-sm_120`), not the GPU vendor's own OpenCL implementation, so the kernels went through a different compiler than on any vendor-ICD host. PoCL caps its PTX target at `sm_75` unless `POCL_CUDA_GPU_ARCH` says otherwise (here: sm_120), so newer hardware is addressed through forward JIT rather than native codegen. The work-group size was forced to 64 rather than derived from kernel register usage. 
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
| OpenCL version | OpenCL 3.0 PoCL HSTR: CUDA-sm_120 |
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
- Total wall time 485.5 s.

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
| ADD | 604288 / 604343 | - | - | - | 6.80 G | - | - | - | 84.11 M |
| SUBTRACT | 604288 / 604343 | - | - | - | 6.33 G | - | - | - | 95.66 M |
| ADDMOD | 604288 / 604343 | - | - | - | 9.87 G | - | - | - | 37.03 M |
| SUBTRACTMOD | 604288 / 604343 | - | - | - | 10.05 G | - | - | - | 45.29 M |
| MULTIPLYOPERANDSCANNING | 604288 / 604343 | - | - | - | 2.60 G | - | - | - | 67.39 M |
| MULTIPLYPRODUCTSCANNING | 604288 / 604343 | - | - | - | 4.41 G | - | - | - | 66.75 M |
| MONTGOMERYMULTIPLICATION | 604288 / 604343 | - | - | - | 8.79 G | - | - | - | 9.93 M |
| COMPARE | 604288 / 604343 | - | - | - | 9.81 G | - | - | - | 118.86 M |
| REDUCE | 118976 / 75542 | - | - | - | 1.80 G | - | - | - | 93.38 M |
| MODMUL | 118976 / 37771 | - | - | - | 819.76 M | - | - | - | 17.30 M |
| MODEXP | 118976 / 9442 | - | - | - | 30.49 M | - | - | - | 160.97 k |
| EXPONENTIATION | 118976 / 9442 | - | - | - | 26.86 M | - | - | - | 513.57 k |
| DIVIDE | 118976 / 75542 | - | - | - | 1.10 G | - | - | - | 49.49 M |
| ISQRT | 118976 / 18885 | - | - | - | 146.17 M | - | - | - | 21.93 M |
| MODMUL_R2 | 604288 / 604343 | - | - | - | 6.37 G | - | - | - | 17.15 M |

#### rsa256(composite) (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 604288 / 604343 | - | - | - | 7.84 G | - | - | - | 82.93 M |
| SUBTRACT | 604288 / 604343 | - | - | - | 6.93 G | - | - | - | 97.23 M |
| ADDMOD | 604288 / 604343 | - | - | - | 9.45 G | - | - | - | 41.98 M |
| SUBTRACTMOD | 604288 / 604343 | - | - | - | 10.80 G | - | - | - | 45.23 M |
| MULTIPLYOPERANDSCANNING | 604288 / 604343 | - | - | - | 2.70 G | - | - | - | 66.78 M |
| MULTIPLYPRODUCTSCANNING | 604288 / 604343 | - | - | - | 4.10 G | - | - | - | 66.81 M |
| MONTGOMERYMULTIPLICATION | 604288 / 604343 | - | - | - | 9.18 G | - | - | - | 9.87 M |
| COMPARE | 604288 / 604343 | - | - | - | 10.61 G | - | - | - | 114.66 M |
| REDUCE | 118976 / 75542 | - | - | - | 1.66 G | - | - | - | 56.92 M |
| MODMUL | 118976 / 37771 | - | - | - | 842.38 M | - | - | - | 17.37 M |
| MODEXP | 118976 / 9442 | - | - | - | 30.52 M | - | - | - | 171.47 k |
| EXPONENTIATION | 118976 / 9442 | - | - | - | 26.39 M | - | - | - | 513.78 k |
| DIVIDE | 118976 / 75542 | - | - | - | 1.03 G | - | - | - | 49.73 M |
| ISQRT | 118976 / 18885 | - | - | - | 144.70 M | - | - | - | 22.05 M |
| MODMUL_R2 | 604288 / 604343 | - | - | - | 6.62 G | - | - | - | 17.16 M |

#### brainpoolP512r1 (512-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 302144 / 302171 | - | - | - | 2.79 G | - | - | - | 67.61 M |
| SUBTRACT | 302144 / 302171 | - | - | - | 2.91 G | - | - | - | 74.73 M |
| ADDMOD | 302144 / 302171 | - | - | - | 3.41 G | - | - | - | 37.16 M |
| SUBTRACTMOD | 302144 / 302171 | - | - | - | 2.75 G | - | - | - | 39.20 M |
| MULTIPLYOPERANDSCANNING | 302144 / 302171 | - | - | - | 548.56 M | - | - | - | 32.45 M |
| MULTIPLYPRODUCTSCANNING | 302144 / 302171 | - | - | - | 1.22 G | - | - | - | 32.38 M |
| MONTGOMERYMULTIPLICATION | 302144 / 302171 | - | - | - | 1.88 G | - | - | - | 4.23 M |
| COMPARE | 302144 / 302171 | - | - | - | 3.37 G | - | - | - | 88.96 M |
| REDUCE | 118976 / 37771 | - | - | - | 87.37 M | - | - | - | 52.64 M |
| MODMUL | 118976 / 18885 | - | - | - | 30.31 M | - | - | - | 8.65 M |
| MODEXP | 118976 / 4721 | - | - | - | 3.68 M | - | - | - | 28.23 k |
| EXPONENTIATION | 118976 / 4721 | - | - | - | 4.15 M | - | - | - | 161.47 k |
| DIVIDE | 118976 / 37771 | - | - | - | 84.72 M | - | - | - | 46.86 M |
| ISQRT | 118976 / 9442 | - | - | - | 7.05 M | - | - | - | 11.96 M |
| MODMUL_R2 | 302144 / 302171 | - | - | - | 1.71 G | - | - | - | 8.73 M |

#### p1024 (1024-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 151040 / 151085 | - | - | - | 1.35 G | - | - | - | 58.33 M |
| SUBTRACT | 151040 / 151085 | - | - | - | 1.36 G | - | - | - | 62.72 M |
| ADDMOD | 151040 / 151085 | - | - | - | 1.44 G | - | - | - | 27.02 M |
| SUBTRACTMOD | 151040 / 151085 | - | - | - | 1.36 G | - | - | - | 32.65 M |
| MULTIPLYOPERANDSCANNING | 151040 / 151085 | - | - | - | 356.68 M | - | - | - | 8.95 M |
| MULTIPLYPRODUCTSCANNING | 151040 / 151085 | - | - | - | 188.90 M | - | - | - | 8.94 M |
| MONTGOMERYMULTIPLICATION | 151040 / 151085 | - | - | - | 586.25 M | - | - | - | 1.32 M |
| COMPARE | 151040 / 151085 | - | - | - | 1.48 G | - | - | - | 128.79 M |
| REDUCE | 118976 / 18885 | - | - | - | 21.66 M | - | - | - | 73.70 M |
| MODMUL | 118976 / 9442 | - | - | - | 7.35 M | - | - | - | 3.04 M |
| MODEXP | 118976 / 2360 | - | - | - | 454.95 k | - | - | - | 4.47 k |
| EXPONENTIATION | 118976 / 2360 | - | - | - | 421.43 k | - | - | - | 33.42 k |
| DIVIDE | 118976 / 18885 | - | - | - | 21.80 M | - | - | - | 40.92 M |
| ISQRT | 118976 / 4721 | - | - | - | 1.56 M | - | - | - | 6.22 M |
| MODMUL_R2 | 151040 / 151085 | - | - | - | 445.90 M | - | - | - | 3.04 M |

#### p2048 (2048-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 118976 / 75542 | - | - | - | 209.69 M | - | - | - | 34.69 M |
| SUBTRACT | 118976 / 75542 | - | - | - | 210.09 M | - | - | - | 28.17 M |
| ADDMOD | 118976 / 75542 | - | - | - | 229.04 M | - | - | - | 20.25 M |
| SUBTRACTMOD | 118976 / 75542 | - | - | - | 226.26 M | - | - | - | 21.84 M |
| MULTIPLYOPERANDSCANNING | 118976 / 75542 | - | - | - | 54.78 M | - | - | - | 2.77 M |
| MULTIPLYPRODUCTSCANNING | 118976 / 75542 | - | - | - | 50.55 M | - | - | - | 2.75 M |
| MONTGOMERYMULTIPLICATION | 118976 / 75542 | - | - | - | 54.26 M | - | - | - | 404.97 k |
| COMPARE | 118976 / 75542 | - | - | - | 785.12 M | - | - | - | 198.72 M |
| REDUCE | 118976 / 9442 | - | - | - | 5.82 M | - | - | - | 51.93 M |
| MODMUL | 118976 / 4721 | - | - | - | 1.79 M | - | - | - | 936.52 k |
| MODEXP | 118976 / 1180 | - | - | - | 16.65 k | - | - | - | 590.7 |
| EXPONENTIATION | 118976 / 1180 | - | - | - | 39.62 k | - | - | - | 5.61 k |
| DIVIDE | 118976 / 9442 | - | - | - | 5.62 M | - | - | - | 26.80 M |
| ISQRT | 118976 / 2360 | - | - | - | 335.05 k | - | - | - | 3.75 M |
| MODMUL_R2 | 118976 / 75542 | - | - | - | 35.03 M | - | - | - | 931.27 k |

## 5. Head to head

Best OpenCL GPU result against best OpenCL CPU result and the CPU libraries.
Ratios above 1.00x mean the GPU is faster than that baseline.

### secp256k1 (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w32-opt | 6.80 G | none | n/a | 84.11 M | n/a |
| SUBTRACT | w32-opt | 6.33 G | none | n/a | 95.66 M | n/a |
| ADDMOD | w32-opt | 9.88 G | none | n/a | 37.03 M | n/a |
| SUBTRACTMOD | w32-opt | 10.05 G | none | n/a | 45.29 M | n/a |
| MULTIPLYOPERANDSCANNING | w32-opt | 2.60 G | none | n/a | 67.39 M | n/a |
| MULTIPLYPRODUCTSCANNING | w32-opt | 4.41 G | none | n/a | 66.75 M | n/a |
| MONTGOMERYMULTIPLICATION | w32-opt | 8.79 G | none | n/a | 9.93 M | n/a |
| COMPARE | w32-opt | 9.81 G | none | n/a | 118.86 M | n/a |
| REDUCE | w32-opt | 1.14 G | none | n/a | 93.38 M | n/a |
| MODMUL | w32-opt | 260.25 M | none | n/a | 17.30 M | n/a |
| MODEXP | w32-opt | 2.42 M | none | n/a | 160.97 k | n/a |
| EXPONENTIATION | w32-opt | 2.13 M | none | n/a | 513.57 k | n/a |
| DIVIDE | w32-opt | 695.88 M | none | n/a | 49.49 M | n/a |
| ISQRT | w32-opt | 23.20 M | none | n/a | 21.93 M | n/a |
| MODMUL_R2 | w32-opt | 6.37 G | none | n/a | 17.15 M | n/a |

### rsa256(composite) (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w32-opt | 7.84 G | none | n/a | 82.93 M | n/a |
| SUBTRACT | w32-opt | 6.93 G | none | n/a | 97.23 M | n/a |
| ADDMOD | w32-opt | 9.46 G | none | n/a | 41.98 M | n/a |
| SUBTRACTMOD | w32-opt | 10.80 G | none | n/a | 45.23 M | n/a |
| MULTIPLYOPERANDSCANNING | w32-opt | 2.70 G | none | n/a | 66.78 M | n/a |
| MULTIPLYPRODUCTSCANNING | w32-opt | 4.10 G | none | n/a | 66.81 M | n/a |
| MONTGOMERYMULTIPLICATION | w32-opt | 9.18 G | none | n/a | 9.87 M | n/a |
| COMPARE | w32-opt | 10.61 G | none | n/a | 114.66 M | n/a |
| REDUCE | w32-opt | 1.05 G | none | n/a | 56.92 M | n/a |
| MODMUL | w32-opt | 267.43 M | none | n/a | 17.37 M | n/a |
| MODEXP | w32-opt | 2.42 M | none | n/a | 171.47 k | n/a |
| EXPONENTIATION | w32-opt | 2.09 M | none | n/a | 513.78 k | n/a |
| DIVIDE | w32-opt | 653.93 M | none | n/a | 49.73 M | n/a |
| ISQRT | w32-opt | 22.97 M | none | n/a | 22.05 M | n/a |
| MODMUL_R2 | w32-opt | 6.62 G | none | n/a | 17.16 M | n/a |

### brainpoolP512r1 (512-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w32-opt | 2.79 G | none | n/a | 67.61 M | n/a |
| SUBTRACT | w32-opt | 2.91 G | none | n/a | 74.73 M | n/a |
| ADDMOD | w32-opt | 3.41 G | none | n/a | 37.16 M | n/a |
| SUBTRACTMOD | w32-opt | 2.75 G | none | n/a | 39.20 M | n/a |
| MULTIPLYOPERANDSCANNING | w32-opt | 548.61 M | none | n/a | 32.45 M | n/a |
| MULTIPLYPRODUCTSCANNING | w32-opt | 1.22 G | none | n/a | 32.38 M | n/a |
| MONTGOMERYMULTIPLICATION | w32-opt | 1.88 G | none | n/a | 4.23 M | n/a |
| COMPARE | w32-opt | 3.37 G | none | n/a | 88.96 M | n/a |
| REDUCE | w32-opt | 27.74 M | none | n/a | 52.64 M | n/a |
| MODMUL | w32-opt | 4.81 M | none | n/a | 8.65 M | n/a |
| MODEXP | w32-opt | 145.98 k | none | n/a | 28.23 k | n/a |
| EXPONENTIATION | w32-opt | 164.87 k | none | n/a | 161.47 k | n/a |
| DIVIDE | w32-opt | 26.90 M | none | n/a | 46.86 M | n/a |
| ISQRT | w32-opt | 559.62 k | none | n/a | 11.96 M | n/a |
| MODMUL_R2 | w32-opt | 1.71 G | none | n/a | 8.73 M | n/a |

### p1024 (1024-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w32-opt | 1.35 G | none | n/a | 58.33 M | n/a |
| SUBTRACT | w32-opt | 1.36 G | none | n/a | 62.72 M | n/a |
| ADDMOD | w32-opt | 1.44 G | none | n/a | 27.02 M | n/a |
| SUBTRACTMOD | w32-opt | 1.36 G | none | n/a | 32.65 M | n/a |
| MULTIPLYOPERANDSCANNING | w32-opt | 356.78 M | none | n/a | 8.95 M | n/a |
| MULTIPLYPRODUCTSCANNING | w32-opt | 188.95 M | none | n/a | 8.94 M | n/a |
| MONTGOMERYMULTIPLICATION | w32-opt | 586.42 M | none | n/a | 1.32 M | n/a |
| COMPARE | w32-opt | 1.48 G | none | n/a | 128.79 M | n/a |
| REDUCE | w32-opt | 3.44 M | none | n/a | 73.70 M | n/a |
| MODMUL | w32-opt | 582.92 k | none | n/a | 3.04 M | n/a |
| MODEXP | w32-opt | 9.02 k | none | n/a | 4.47 k | n/a |
| EXPONENTIATION | w32-opt | 8.36 k | none | n/a | 33.42 k | n/a |
| DIVIDE | w32-opt | 3.46 M | none | n/a | 40.92 M | n/a |
| ISQRT | w32-opt | 62.05 k | none | n/a | 6.22 M | n/a |
| MODMUL_R2 | w32-opt | 446.03 M | none | n/a | 3.04 M | n/a |

### p2048 (2048-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w32-opt | 133.14 M | none | n/a | 34.69 M | n/a |
| SUBTRACT | w32-opt | 133.39 M | none | n/a | 28.17 M | n/a |
| ADDMOD | w32-opt | 145.43 M | none | n/a | 20.25 M | n/a |
| SUBTRACTMOD | w32-opt | 143.66 M | none | n/a | 21.84 M | n/a |
| MULTIPLYOPERANDSCANNING | w32-opt | 34.78 M | none | n/a | 2.77 M | n/a |
| MULTIPLYPRODUCTSCANNING | w32-opt | 32.10 M | none | n/a | 2.75 M | n/a |
| MONTGOMERYMULTIPLICATION | w32-opt | 34.45 M | none | n/a | 404.97 k | n/a |
| COMPARE | w32-opt | 498.50 M | none | n/a | 198.72 M | n/a |
| REDUCE | w32-opt | 461.86 k | none | n/a | 51.93 M | n/a |
| MODMUL | w32-opt | 70.83 k | none | n/a | 936.52 k | n/a |
| MODEXP | w32-opt | 165.1 | none | n/a | 590.7 | n/a |
| EXPONENTIATION | w32-opt | 392.9 | none | n/a | 5.61 k | n/a |
| DIVIDE | w32-opt | 445.97 k | none | n/a | 26.80 M | n/a |
| ISQRT | w32-opt | 6.65 k | none | n/a | 3.75 M | n/a |
| MODMUL_R2 | w32-opt | 22.24 M | none | n/a | 931.27 k | n/a |

## 6. Raw data

Also written to `NVIDIA_GeForce_RTX_5090_Report.csv` for analysis.

```csv
kind,device,device_type,kernel,modulus,bits,operation,items,seconds,ops_per_sec,mismatches
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,ADD,604343,0.007185287,84108400.079,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,ADD,604343,0.000598236,1010208367.608,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,ADD,604343,0.001783385,338874118.498,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,secp256k1,256,ADD,700000,0.000095840,7303839732.888,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,ADD,604288,0.000088899,6797454284.269,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,ADD,604288,0.016630186,36336815.473,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,SUBTRACT,604343,0.006317330,95664310.185,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,SUBTRACT,604343,0.000657688,918890258.730,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,SUBTRACT,604343,0.001848838,326877195.643,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,secp256k1,256,SUBTRACT,700000,0.000096128,7281957390.146,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,SUBTRACT,604288,0.000095531,6325562321.907,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,SUBTRACT,604288,0.016622570,36353464.263,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,ADDMOD,604343,0.016318892,37033335.588,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,ADDMOD,604343,0.001750762,345188576.656,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,ADDMOD,604343,0.002325583,259867308.848,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,secp256k1,256,ADDMOD,700000,0.000090880,7702464788.732,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,ADDMOD,604288,0.000061196,9874663085.313,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,ADDMOD,604288,0.016502854,36617181.086,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,SUBTRACTMOD,604343,0.013343285,45291920.793,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,SUBTRACTMOD,604343,0.000741628,814886979.861,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,SUBTRACTMOD,604343,0.002304804,262210143.689,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,secp256k1,256,SUBTRACTMOD,700000,0.000096768,7233796296.296,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,SUBTRACTMOD,604288,0.000060134,10049045766.941,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,SUBTRACTMOD,604288,0.016552377,36507627.124,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,MULTIPLYOPERANDSCANNING,604343,0.008967720,67390931.238,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,604343,0.000622502,970829275.073,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,604343,0.001854370,325902038.113,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,MULTIPLYOPERANDSCANNING,604288,0.000232401,2600194556.322,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,MULTIPLYOPERANDSCANNING,604288,0.022224453,27190230.454,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,MULTIPLYPRODUCTSCANNING,604343,0.009053155,66754961.128,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,604343,0.000573358,1054041219.225,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,604343,0.001945974,310560663.967,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,secp256k1,256,MULTIPLYPRODUCTSCANNING,700000,0.000094528,7405213270.142,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,MULTIPLYPRODUCTSCANNING,604288,0.000137061,4408902066.627,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,MULTIPLYPRODUCTSCANNING,604288,0.021997232,27471092.739,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,MONTGOMERYMULTIPLICATION,604343,0.060843073,9932815.184,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,604343,0.002986409,202364444.860,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,604343,0.001920816,314628258.404,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,secp256k1,256,MONTGOMERYMULTIPLICATION,700000,0.000096992,7217090069.284,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,MONTGOMERYMULTIPLICATION,604288,0.000068721,8793366166.019,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,MONTGOMERYMULTIPLICATION,604288,0.016652839,36287386.137,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,COMPARE,604343,0.005084381,118862650.772,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,COMPARE,604343,0.000772075,782751646.187,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,COMPARE,604343,0.001783555,338841825.013,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,secp256k1,256,COMPARE,700000,0.000091456,7653953813.856,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,COMPARE,604288,0.000061587,9811909333.555,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,COMPARE,604288,0.016552708,36506896.390,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,REDUCE,75542,0.000808995,93377587.587,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,REDUCE,75542,0.000141819,532664846.750,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,REDUCE,75542,0.000234786,321748376.608,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,secp256k1,256,REDUCE,700000,0.000092704,7550914739.386,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,REDUCE,118976,0.000066266,1795432447.942,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,REDUCE,118976,0.003584865,33188420.548,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,MODMUL,37771,0.002182873,17303343.158,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,MODMUL,37771,0.000113195,333680694.856,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,MODMUL,37771,0.000280692,134563859.157,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,secp256k1,256,MODMUL,700000,0.000221216,3164328077.535,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,MODMUL,118976,0.000145135,819761015.495,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,MODMUL,118976,0.003615663,32905722.296,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,MODEXP,9442,0.058656272,160971.703,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,MODEXP,9442,0.002888253,3269103.979,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,MODEXP,9442,0.005462609,1728478.096,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,secp256k1,256,MODEXP,700000,0.063936286,10948399.474,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,MODEXP,118976,0.003902307,30488631.367,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,MODEXP,118976,0.007405065,16066840.737,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,EXPONENTIATION,9442,0.018384947,513572.327,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,EXPONENTIATION,9442,0.000924725,10210602.208,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,EXPONENTIATION,9442,0.013743007,687040.321,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,EXPONENTIATION,118976,0.004429047,26862663.708,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,EXPONENTIATION,118976,0.007915936,15029934.657,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,DIVIDE,75542,0.001526257,49494947.689,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,DIVIDE,75542,0.000082777,912595230.924,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,DIVIDE,75542,0.000231850,325822457.844,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,secp256k1,256,DIVIDE,700000,0.000090560,7729681978.799,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,DIVIDE,118976,0.000108556,1095988428.825,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,DIVIDE,118976,0.002514142,47322707.848,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,ISQRT,18885,0.000861035,21932901.380,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,ISQRT,18885,0.000052029,362970412.932,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,ISQRT,118976,0.000813965,146168449.350,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,ISQRT,118976,0.002676751,44447914.672,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,MODMUL_R2,604343,0.035245544,17146649.779,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,MODMUL_R2,604343,0.001693805,356796093.241,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,MODMUL_R2,604343,0.004551760,132771281.985,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,secp256k1,256,MODMUL_R2,700000,0.000106624,6565126050.420,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,MODMUL_R2,604288,0.000094800,6374360076.248,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,MODMUL_R2,604288,0.008507257,71032062.370,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,ADD,604343,0.007287573,82927884.710,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,ADD,604343,0.000599288,1008434759.146,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,ADD,604343,0.001738240,347675232.513,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,rsa256(composite),256,ADD,700000,0.000085824,8156226696.495,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,ADD,604288,0.000077037,7844136226.400,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,ADD,604288,0.008592399,70328206.529,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACT,604343,0.006215598,97230065.886,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACT,604343,0.000604177,1000274929.696,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACT,604343,0.001841245,328225197.898,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,rsa256(composite),256,SUBTRACT,700000,0.000096096,7284382284.382,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,SUBTRACT,604288,0.000087165,6932687621.010,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,SUBTRACT,604288,0.008565127,70552136.094,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,ADDMOD,604343,0.014395556,41981219.914,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,ADDMOD,604343,0.000764361,790651154.973,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,ADDMOD,604343,0.002119834,285089766.596,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,rsa256(composite),256,ADDMOD,700000,0.000090208,7759843916.282,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,ADDMOD,604288,0.000063912,9454998897.505,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,ADDMOD,604288,0.008530782,70836178.329,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACTMOD,604343,0.013362755,45225928.107,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACTMOD,604343,0.000752569,803040047.399,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACTMOD,604343,0.002351223,257033464.479,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,rsa256(composite),256,SUBTRACTMOD,700000,0.000092896,7535308301.757,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,SUBTRACTMOD,604288,0.000055966,10797418989.588,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,SUBTRACTMOD,604288,0.016577835,36451563.270,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYOPERANDSCANNING,604343,0.009049708,66780386.552,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,604343,0.000559993,1079197248.476,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,604343,0.001957035,308805424.104,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,604288,0.000223604,2702494345.795,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,604288,0.022165661,27262349.657,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,604343,0.009045140,66814113.169,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,604343,0.000597595,1011291921.093,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,604343,0.001913102,315896924.750,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,700000,0.000093920,7453151618.399,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,604288,0.000147390,4099926856.122,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,604288,0.022010085,27455050.793,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MONTGOMERYMULTIPLICATION,604343,0.061219920,9871672.498,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,604343,0.002986738,202342154.382,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,604343,0.002127637,284044234.053,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,rsa256(composite),256,MONTGOMERYMULTIPLICATION,700000,0.000089792,7795794725.588,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,604288,0.000065845,9177403342.133,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,604288,0.016552696,36506923.092,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,COMPARE,604343,0.005270714,114660557.725,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,COMPARE,604343,0.000623604,969113334.160,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,COMPARE,604343,0.001878415,321730288.301,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,rsa256(composite),256,COMPARE,700000,0.000090432,7740622788.393,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,COMPARE,604288,0.000056938,10613124718.316,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,COMPARE,604288,0.016551565,36509417.367,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,REDUCE,75542,0.001327219,56917511.462,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,REDUCE,75542,0.000103917,726945732.825,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,REDUCE,75542,0.000231239,326683628.761,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,rsa256(composite),256,REDUCE,700000,0.000084288,8304859529.233,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,REDUCE,118976,0.000071837,1656196919.685,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,REDUCE,118976,0.003560167,33418656.682,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MODMUL,37771,0.002174267,17371830.882,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MODMUL,37771,0.000106272,355417988.842,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MODMUL,37771,0.000276986,136364271.167,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,rsa256(composite),256,MODMUL,700000,0.000218752,3199970743.125,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,MODMUL,118976,0.000141238,842380350.222,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,MODMUL,118976,0.003643145,32657497.076,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MODEXP,9442,0.055065827,171467.506,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MODEXP,9442,0.010324068,914561.961,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MODEXP,9442,0.005451157,1732109.352,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,rsa256(composite),256,MODEXP,700000,0.062912092,11126636.832,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,MODEXP,118976,0.003897869,30523344.715,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,MODEXP,118976,0.007414782,16045785.397,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,EXPONENTIATION,9442,0.018377532,513779.543,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,EXPONENTIATION,9442,0.000934975,10098665.859,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,EXPONENTIATION,9442,0.013953936,676654.957,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,EXPONENTIATION,118976,0.004508527,26389106.711,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,EXPONENTIATION,118976,0.007857503,15141705.696,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,DIVIDE,75542,0.001519164,49726038.843,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,DIVIDE,75542,0.000114728,658444331.084,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,DIVIDE,75542,0.000241859,312339156.348,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,rsa256(composite),256,DIVIDE,700000,0.000089856,7790242165.242,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,DIVIDE,118976,0.000115519,1029924416.173,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,DIVIDE,118976,0.004636941,25658294.917,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,ISQRT,18885,0.000856385,22051995.050,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,ISQRT,18885,0.000052781,357798685.380,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,ISQRT,118976,0.000822221,144700768.137,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,ISQRT,118976,0.004321422,27531678.471,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MODMUL_R2,604343,0.035217034,17160530.929,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MODMUL_R2,604343,0.001704685,354518837.370,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MODMUL_R2,604343,0.004453392,135703975.730,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,rsa256(composite),256,MODMUL_R2,700000,0.000101280,6911532385.466,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,MODMUL_R2,604288,0.000091304,6618412429.402,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,MODMUL_R2,604288,0.008609190,70191040.980,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,ADD,302171,0.004469614,67605614.110,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,ADD,302171,0.000183980,1642412391.481,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,ADD,302171,0.000474711,636536824.748,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,brainpoolP512r1,512,ADD,700000,0.000146560,4776200873.362,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,ADD,302144,0.000108146,2793855961.752,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,ADD,302144,0.016627389,18171463.955,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACT,302171,0.004043685,74726640.935,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACT,302171,0.000168991,1788088858.945,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACT,302171,0.000476554,634075308.814,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,brainpoolP512r1,512,SUBTRACT,700000,0.000148640,4709364908.504,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,SUBTRACT,302144,0.000103707,2913439956.811,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,SUBTRACT,302144,0.016596630,18205141.665,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,ADDMOD,302171,0.008132084,37157879.643,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,ADDMOD,302171,0.000407283,741918655.575,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,ADDMOD,302171,0.002406668,125555751.894,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,brainpoolP512r1,512,ADDMOD,700000,0.000147872,4733823847.652,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,ADDMOD,302144,0.000088709,3406006248.462,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,ADDMOD,302144,0.016628070,18170719.714,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACTMOD,302171,0.007709212,39196094.726,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACTMOD,302171,0.000390932,772950354.250,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACTMOD,302171,0.002578384,117193947.377,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,brainpoolP512r1,512,SUBTRACTMOD,700000,0.000142688,4905808477.237,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,SUBTRACTMOD,302144,0.000109738,2753322255.284,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,SUBTRACTMOD,302144,0.016567075,18237618.973,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,302171,0.009312445,32448084.031,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,302171,0.000564011,535753690.974,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,302171,0.000552970,546450997.984,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,302144,0.000550795,548559927.683,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,302144,0.022513641,13420485.802,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,302171,0.009330951,32383730.544,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,302171,0.000560113,539482215.452,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,302171,0.000555715,543751686.499,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,700000,0.000153184,4569667850.428,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,302144,0.000247640,1220093341.597,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,302144,0.022222177,13596507.592,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,302171,0.071410019,4231493.052,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,302171,0.003532235,85546685.357,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,302171,0.001006090,300341893.174,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,700000,0.002374688,294775566.306,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,302144,0.000160655,1880700801.560,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,302144,0.016787882,17997743.889,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,COMPARE,302171,0.003396547,88964176.232,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,COMPARE,302171,0.000143513,2105532700.934,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,COMPARE,302171,0.000592935,509618984.881,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,brainpoolP512r1,512,COMPARE,700000,0.000247200,2831715210.356,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,COMPARE,302144,0.000089711,3367968831.741,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,COMPARE,302144,0.016604113,18196937.188,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,REDUCE,37771,0.000717481,52643893.036,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,REDUCE,37771,0.000112754,334985720.380,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,REDUCE,37771,0.000285202,132435877.963,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,brainpoolP512r1,512,REDUCE,700000,0.000243744,2871865563.870,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,REDUCE,118976,0.001361824,87365173.989,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,REDUCE,118976,0.008033838,14809360.311,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL,18885,0.002182072,8654618.066,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL,18885,0.000107384,175864099.913,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL,18885,0.000301342,62669659.422,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,brainpoolP512r1,512,MODMUL,700000,0.002380640,294038577.861,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,MODMUL,118976,0.003924879,30313291.089,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,MODMUL,118976,0.010560283,11266364.767,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MODEXP,4721,0.167210842,28233.815,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MODEXP,4721,0.009341021,505405.146,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MODEXP,4721,0.008363084,564504.677,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,brainpoolP512r1,512,MODEXP,700000,0.211561188,3308735.438,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,MODEXP,118976,0.032340952,3678803.275,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,MODEXP,118976,0.034505041,3448075.883,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,EXPONENTIATION,4721,0.029237555,161470.410,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,EXPONENTIATION,4721,0.001693905,2787051.443,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,EXPONENTIATION,4721,0.020009760,235934.866,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,EXPONENTIATION,118976,0.028634438,4154996.886,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,EXPONENTIATION,118976,0.030668478,3879423.036,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,DIVIDE,37771,0.000805959,46864660.679,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,DIVIDE,37771,0.000062699,602418247.010,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,DIVIDE,37771,0.000142631,264816411.719,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,brainpoolP512r1,512,DIVIDE,700000,0.002377216,294462093.474,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,DIVIDE,118976,0.001404386,84717448.647,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,DIVIDE,118976,0.003513518,33862356.545,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,ISQRT,9442,0.000789298,11962529.070,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,ISQRT,9442,0.000041479,227633194.362,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,ISQRT,118976,0.016872075,7051651.949,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,ISQRT,118976,0.018726053,6353501.127,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL_R2,302171,0.034616544,8729092.042,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL_R2,302171,0.001701821,177557464.779,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL_R2,302171,0.003884713,77784640.780,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,brainpoolP512r1,512,MODMUL_R2,700000,0.000170240,4111842105.263,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,MODMUL_R2,302144,0.000176415,1712690312.055,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,MODMUL_R2,302144,0.004176949,72336053.404,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p1024,1024,ADD,151085,0.002590347,58326165.889,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p1024,1024,ADD,151085,0.000130668,1156250517.032,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p1024,1024,ADD,151085,0.000139205,1085341621.353,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p1024,1024,ADD,700000,0.000460992,1518464528.669,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,ADD,151040,0.000112003,1348535929.579,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,ADD,151040,0.004158463,36321111.070,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p1024,1024,SUBTRACT,151085,0.002408912,62719181.597,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p1024,1024,SUBTRACT,151085,0.000166155,909301926.651,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p1024,1024,SUBTRACT,151085,0.000145206,1040487343.028,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p1024,1024,SUBTRACT,700000,0.004714080,148491328.106,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,SUBTRACT,151040,0.000111161,1358749521.688,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,SUBTRACT,151040,0.004086376,36961846.488,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p1024,1024,ADDMOD,151085,0.005590792,27023899.071,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p1024,1024,ADDMOD,151085,0.000292745,516097710.161,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p1024,1024,ADDMOD,151085,0.001831095,82510739.588,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p1024,1024,ADDMOD,700000,0.000246464,2840171384.056,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,ADDMOD,151040,0.000105099,1437120451.372,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,ADDMOD,151040,0.008511976,17744410.712,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p1024,1024,SUBTRACTMOD,151085,0.004627954,32646175.953,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p1024,1024,SUBTRACTMOD,151085,0.000243032,621667398.518,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p1024,1024,SUBTRACTMOD,151085,0.000874730,172721856.639,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p1024,1024,SUBTRACTMOD,700000,0.000246464,2840171384.056,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,SUBTRACTMOD,151040,0.000111061,1359971531.450,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,SUBTRACTMOD,151040,0.008414261,17950477.080,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p1024,1024,MULTIPLYOPERANDSCANNING,151085,0.016873427,8954019.849,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p1024,1024,MULTIPLYOPERANDSCANNING,151085,0.000986333,153178525.791,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p1024,1024,MULTIPLYOPERANDSCANNING,151085,0.001047819,144189990.178,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,MULTIPLYOPERANDSCANNING,151040,0.000423464,356677315.137,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,MULTIPLYOPERANDSCANNING,151040,0.011590849,13030969.619,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p1024,1024,MULTIPLYPRODUCTSCANNING,151085,0.016902401,8938670.966,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,151085,0.000982284,153809920.013,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,151085,0.001057076,142927277.687,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p1024,1024,MULTIPLYPRODUCTSCANNING,700000,0.000449600,1556939501.779,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,MULTIPLYPRODUCTSCANNING,151040,0.000799588,188897299.167,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,MULTIPLYPRODUCTSCANNING,151040,0.011929754,12660780.836,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p1024,1024,MONTGOMERYMULTIPLICATION,151085,0.114112251,1324003.328,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p1024,1024,MONTGOMERYMULTIPLICATION,151085,0.005702214,26495848.421,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p1024,1024,MONTGOMERYMULTIPLICATION,151085,0.001745192,86572131.454,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p1024,1024,MONTGOMERYMULTIPLICATION,700000,0.008761887,79891466.302,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,MONTGOMERYMULTIPLICATION,151040,0.000257639,586246689.387,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,MONTGOMERYMULTIPLICATION,151040,0.008578072,17607686.539,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p1024,1024,COMPARE,151085,0.001173067,128794870.221,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p1024,1024,COMPARE,151085,0.000066065,2286914518.625,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p1024,1024,COMPARE,151085,0.000089550,1687159421.956,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p1024,1024,COMPARE,700000,0.000458752,1525878906.250,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,COMPARE,151040,0.000101944,1481598234.055,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,COMPARE,151040,0.016607882,9094476.954,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p1024,1024,REDUCE,18885,0.000256226,73704483.476,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p1024,1024,REDUCE,18885,0.000017292,1092131973.191,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p1024,1024,REDUCE,18885,0.000096373,195957318.666,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p1024,1024,REDUCE,700000,0.000454048,1541687222.496,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,REDUCE,118976,0.005492326,21662225.215,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,REDUCE,118976,0.008751081,13595577.398,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p1024,1024,MODMUL,9442,0.003110463,3035560.912,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p1024,1024,MODMUL,9442,0.000156367,60383508.453,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p1024,1024,MODMUL,9442,0.000337771,27953858.069,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p1024,1024,MODMUL,700000,0.004740128,147675337.037,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,MODMUL,118976,0.016197654,7345261.149,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,MODMUL,118976,0.019592819,6072428.944,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p1024,1024,MODEXP,2360,0.528405607,4466.266,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p1024,1024,MODEXP,2360,0.032534633,72538.086,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p1024,1024,MODEXP,2360,0.025675804,91915.330,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p1024,1024,MODEXP,700000,1.212423444,577356.041,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,MODEXP,118976,0.261515756,454947.732,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,MODEXP,118976,0.265100690,448795.512,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p1024,1024,EXPONENTIATION,2360,0.070613158,33421.533,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p1024,1024,EXPONENTIATION,2360,0.004689450,503257.296,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p1024,1024,EXPONENTIATION,2360,0.068343089,34531.655,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,EXPONENTIATION,118976,0.282315452,421429.288,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,EXPONENTIATION,118976,0.285922972,416112.071,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p1024,1024,DIVIDE,18885,0.000461556,40915943.652,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p1024,1024,DIVIDE,18885,0.000028173,670323278.830,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p1024,1024,DIVIDE,18885,0.000114758,164563642.789,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p1024,1024,DIVIDE,700000,0.004723648,148190551.032,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,DIVIDE,118976,0.005457929,21798744.132,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,DIVIDE,118976,0.009878388,12044070.461,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p1024,1024,ISQRT,4721,0.000758850,6221255.719,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p1024,1024,ISQRT,4721,0.000038213,123544031.369,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,ISQRT,118976,0.076082379,1563778.649,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,ISQRT,118976,0.080014853,1486923.933,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p1024,1024,MODMUL_R2,151085,0.049660319,3042368.697,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p1024,1024,MODMUL_R2,151085,0.002509624,60202245.465,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p1024,1024,MODMUL_R2,151085,0.006852947,22046719.647,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p1024,1024,MODMUL_R2,700000,0.000556064,1258847902.400,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,MODMUL_R2,151040,0.000338733,445896964.002,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,MODMUL_R2,151040,0.004335880,34834913.901,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p2048,2048,ADD,75542,0.002177944,34685006.611,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p2048,2048,ADD,75542,0.000092356,817943608.664,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p2048,2048,ADD,75542,0.000113926,663080144.765,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p2048,2048,ADD,700000,0.000459456,1523540883.131,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,ADD,118976,0.000567397,209687359.049,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,ADD,118976,0.013876997,8573612.761,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p2048,2048,SUBTRACT,75542,0.002682051,28165756.585,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p2048,2048,SUBTRACT,75542,0.000079932,945079213.744,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p2048,2048,SUBTRACT,75542,0.000112473,671645440.134,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p2048,2048,SUBTRACT,700000,0.000459776,1522480512.249,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,SUBTRACT,118976,0.000566315,210088060.131,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,SUBTRACT,118976,0.013871167,8577216.268,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p2048,2048,ADDMOD,75542,0.003730260,20251134.937,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p2048,2048,ADDMOD,75542,0.000191544,394384683.471,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p2048,2048,ADDMOD,75542,0.000685361,110222207.096,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p2048,2048,ADDMOD,700000,0.000464160,1508100654.947,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,ADDMOD,118976,0.000519456,229039665.255,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,ADDMOD,118976,0.013859665,8584334.457,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p2048,2048,SUBTRACTMOD,75542,0.003458254,21843970.727,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p2048,2048,SUBTRACTMOD,75542,0.000185923,406307582.989,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p2048,2048,SUBTRACTMOD,75542,0.000676203,111714964.337,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p2048,2048,SUBTRACTMOD,700000,0.009610528,72836788.988,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,SUBTRACTMOD,118976,0.000525838,226259799.947,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,SUBTRACTMOD,118976,0.013862310,8582696.549,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p2048,2048,MULTIPLYOPERANDSCANNING,75542,0.027287623,2768361.304,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p2048,2048,MULTIPLYOPERANDSCANNING,75542,0.001567214,48201464.401,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p2048,2048,MULTIPLYOPERANDSCANNING,75542,0.001689858,44703163.452,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,MULTIPLYOPERANDSCANNING,118976,0.002171692,54784933.357,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,MULTIPLYOPERANDSCANNING,118976,0.019784902,6013474.318,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p2048,2048,MULTIPLYPRODUCTSCANNING,75542,0.027460011,2750982.151,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,75542,0.001564199,48294363.399,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,75542,0.001709214,44196924.660,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p2048,2048,MULTIPLYPRODUCTSCANNING,700000,0.009617024,72787590.007,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,MULTIPLYPRODUCTSCANNING,118976,0.002353567,50551356.073,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,MULTIPLYPRODUCTSCANNING,118976,0.011152017,10668563.208,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p2048,2048,MONTGOMERYMULTIPLICATION,75542,0.186538416,404967.522,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p2048,2048,MONTGOMERYMULTIPLICATION,75542,0.009457350,7987649.809,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p2048,2048,MONTGOMERYMULTIPLICATION,75542,0.003332545,22667960.636,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p2048,2048,MONTGOMERYMULTIPLICATION,700000,0.000977280,716273739.358,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,MONTGOMERYMULTIPLICATION,118976,0.002192762,54258510.948,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,MONTGOMERYMULTIPLICATION,118976,0.015377908,7736813.089,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p2048,2048,COMPARE,75542,0.000380142,198720412.149,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p2048,2048,COMPARE,75542,0.000145866,517886857.510,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p2048,2048,COMPARE,75542,0.000059403,1271686327.477,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p2048,2048,COMPARE,700000,0.000457248,1530897893.484,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,COMPARE,118976,0.000151538,785124113.097,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,COMPARE,118976,0.013392048,8884078.103,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p2048,2048,REDUCE,9442,0.000181835,51926221.977,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p2048,2048,REDUCE,9442,0.000014157,666947589.120,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p2048,2048,REDUCE,9442,0.000080733,116953277.641,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p2048,2048,REDUCE,700000,0.000461600,1516464471.404,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,REDUCE,118976,0.020443472,5819755.056,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,REDUCE,118976,0.027392380,4343397.699,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p2048,2048,MODMUL,4721,0.005040977,936524.820,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p2048,2048,MODMUL,4721,0.000245807,19206125.606,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p2048,2048,MODMUL,4721,0.000544213,8674912.350,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p2048,2048,MODMUL,700000,0.009618752,72774513.783,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,MODMUL,118976,0.066651865,1785036.319,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,MODMUL,118976,0.073620182,1616078.591,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p2048,2048,MODEXP,1180,1.997553209,590.723,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p2048,2048,MODEXP,1180,0.119384366,9884.041,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p2048,2048,MODEXP,1180,0.093119777,12671.852,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p2048,2048,MODEXP,700000,6.987298489,100181.780,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,MODEXP,118976,7.145882232,16649.589,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,MODEXP,118976,7.158294992,16620.718,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p2048,2048,EXPONENTIATION,1180,0.210283312,5611.477,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p2048,2048,EXPONENTIATION,1180,0.015258181,77335.562,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p2048,2048,EXPONENTIATION,1180,0.169478082,6962.552,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,EXPONENTIATION,118976,3.003139334,39617.209,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,EXPONENTIATION,118976,3.008653507,39544.600,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p2048,2048,DIVIDE,9442,0.000352289,26801853.985,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p2048,2048,DIVIDE,9442,0.000027111,348271495.511,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p2048,2048,DIVIDE,9442,0.000102415,92193542.188,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p2048,2048,DIVIDE,700000,0.000884864,791082019.384,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,DIVIDE,118976,0.021171733,5619568.308,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,DIVIDE,118976,0.030015590,3963806.821,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p2048,2048,ISQRT,2360,0.000629475,3749155.403,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p2048,2048,ISQRT,2360,0.000032572,72454700.767,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,ISQRT,118976,0.355095594,335053.439,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,ISQRT,118976,0.360934437,329633.274,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p2048,2048,MODMUL_R2,75542,0.081116996,931272.159,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p2048,2048,MODMUL_R2,75542,0.004041951,18689489.693,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p2048,2048,MODMUL_R2,75542,0.008868613,8517904.614,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p2048,2048,MODMUL_R2,700000,0.001894976,369397818.231,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,MODMUL_R2,118976,0.003395936,35034818.231,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,MODMUL_R2,118976,0.009943261,11965490.946,0
```
