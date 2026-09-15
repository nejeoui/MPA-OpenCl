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
- Total wall time 503.5 s.

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
| ADD | 604288 / 604343 | - | - | - | 6.43 G | - | - | - | 83.43 M |
| SUBTRACT | 604288 / 604343 | - | - | - | 6.14 G | - | - | - | 97.04 M |
| ADDMOD | 604288 / 604343 | - | - | - | 9.23 G | - | - | - | 37.00 M |
| SUBTRACTMOD | 604288 / 604343 | - | - | - | 8.29 G | - | - | - | 45.23 M |
| MULTIPLYOPERANDSCANNING | 604288 / 604343 | - | - | - | 2.95 G | - | - | - | 67.55 M |
| MULTIPLYPRODUCTSCANNING | 604288 / 604343 | - | - | - | 3.84 G | - | - | - | 66.87 M |
| MONTGOMERYMULTIPLICATION | 604288 / 604343 | - | - | - | 8.61 G | - | - | - | 9.94 M |
| COMPARE | 604288 / 604343 | - | - | - | 11.48 G | - | - | - | 116.64 M |
| REDUCE | 118976 / 75542 | - | - | - | 1.28 G | - | - | - | 93.03 M |
| MODMUL | 118976 / 37771 | - | - | - | 837.68 M | - | - | - | 17.27 M |
| MODEXP | 118976 / 9442 | - | - | - | 30.37 M | - | - | - | 160.96 k |
| EXPONENTIATION | 118976 / 9442 | - | - | - | 26.89 M | - | - | - | 512.00 k |
| DIVIDE | 118976 / 75542 | - | - | - | 1.09 G | - | - | - | 48.27 M |
| ISQRT | 118976 / 18885 | - | - | - | 145.79 M | - | - | - | 22.12 M |
| MODMUL_R2 | 604288 / 604343 | - | - | - | 5.19 G | - | - | - | 17.29 M |

#### rsa256(composite) (256-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 604288 / 604343 | - | - | - | 7.80 G | - | - | - | 81.04 M |
| SUBTRACT | 604288 / 604343 | - | - | - | 7.44 G | - | - | - | 97.33 M |
| ADDMOD | 604288 / 604343 | - | - | - | 10.71 G | - | - | - | 42.20 M |
| SUBTRACTMOD | 604288 / 604343 | - | - | - | 10.07 G | - | - | - | 45.21 M |
| MULTIPLYOPERANDSCANNING | 604288 / 604343 | - | - | - | 3.10 G | - | - | - | 66.69 M |
| MULTIPLYPRODUCTSCANNING | 604288 / 604343 | - | - | - | 4.37 G | - | - | - | 66.98 M |
| MONTGOMERYMULTIPLICATION | 604288 / 604343 | - | - | - | 8.36 G | - | - | - | 9.86 M |
| COMPARE | 604288 / 604343 | - | - | - | 9.84 G | - | - | - | 117.46 M |
| REDUCE | 118976 / 75542 | - | - | - | 1.26 G | - | - | - | 57.37 M |
| MODMUL | 118976 / 37771 | - | - | - | 861.76 M | - | - | - | 17.22 M |
| MODEXP | 118976 / 9442 | - | - | - | 30.40 M | - | - | - | 172.72 k |
| EXPONENTIATION | 118976 / 9442 | - | - | - | 26.97 M | - | - | - | 515.95 k |
| DIVIDE | 118976 / 75542 | - | - | - | 1.03 G | - | - | - | 49.65 M |
| ISQRT | 118976 / 18885 | - | - | - | 144.21 M | - | - | - | 22.03 M |
| MODMUL_R2 | 604288 / 604343 | - | - | - | 5.82 G | - | - | - | 17.31 M |

#### brainpoolP512r1 (512-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 302144 / 302171 | - | - | - | 2.77 G | - | - | - | 66.91 M |
| SUBTRACT | 302144 / 302171 | - | - | - | 2.84 G | - | - | - | 72.12 M |
| ADDMOD | 302144 / 302171 | - | - | - | 3.30 G | - | - | - | 37.44 M |
| SUBTRACTMOD | 302144 / 302171 | - | - | - | 3.19 G | - | - | - | 40.36 M |
| MULTIPLYOPERANDSCANNING | 302144 / 302171 | - | - | - | 559.71 M | - | - | - | 32.30 M |
| MULTIPLYPRODUCTSCANNING | 302144 / 302171 | - | - | - | 1.20 G | - | - | - | 31.95 M |
| MONTGOMERYMULTIPLICATION | 302144 / 302171 | - | - | - | 1.99 G | - | - | - | 4.23 M |
| COMPARE | 302144 / 302171 | - | - | - | 3.24 G | - | - | - | 85.31 M |
| REDUCE | 118976 / 37771 | - | - | - | 80.12 M | - | - | - | 53.55 M |
| MODMUL | 118976 / 18885 | - | - | - | 30.45 M | - | - | - | 8.77 M |
| MODEXP | 118976 / 4721 | - | - | - | 3.67 M | - | - | - | 28.30 k |
| EXPONENTIATION | 118976 / 4721 | - | - | - | 4.20 M | - | - | - | 160.40 k |
| DIVIDE | 118976 / 37771 | - | - | - | 84.69 M | - | - | - | 46.97 M |
| ISQRT | 118976 / 9442 | - | - | - | 7.04 M | - | - | - | 12.11 M |
| MODMUL_R2 | 302144 / 302171 | - | - | - | 1.70 G | - | - | - | 8.74 M |

#### p1024 (1024-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 151040 / 151085 | - | - | - | 1.20 G | - | - | - | 47.56 M |
| SUBTRACT | 151040 / 151085 | - | - | - | 1.24 G | - | - | - | 52.12 M |
| ADDMOD | 151040 / 151085 | - | - | - | 1.35 G | - | - | - | 27.27 M |
| SUBTRACTMOD | 151040 / 151085 | - | - | - | 735.24 M | - | - | - | 32.50 M |
| MULTIPLYOPERANDSCANNING | 151040 / 151085 | - | - | - | 170.78 M | - | - | - | 8.99 M |
| MULTIPLYPRODUCTSCANNING | 151040 / 151085 | - | - | - | 182.51 M | - | - | - | 9.05 M |
| MONTGOMERYMULTIPLICATION | 151040 / 151085 | - | - | - | 344.04 M | - | - | - | 1.34 M |
| COMPARE | 151040 / 151085 | - | - | - | 1.45 G | - | - | - | 117.45 M |
| REDUCE | 118976 / 18885 | - | - | - | 21.59 M | - | - | - | 75.64 M |
| MODMUL | 118976 / 9442 | - | - | - | 7.34 M | - | - | - | 3.08 M |
| MODEXP | 118976 / 2360 | - | - | - | 452.78 k | - | - | - | 4.48 k |
| EXPONENTIATION | 118976 / 2360 | - | - | - | 421.46 k | - | - | - | 33.14 k |
| DIVIDE | 118976 / 18885 | - | - | - | 21.79 M | - | - | - | 40.56 M |
| ISQRT | 118976 / 4721 | - | - | - | 1.56 M | - | - | - | 6.21 M |
| MODMUL_R2 | 151040 / 151085 | - | - | - | 444.14 M | - | - | - | 3.07 M |

#### p2048 (2048-bit)

| Operation | items dev/cpu | w8 | w16 | w32 | w32-opt | w32-o64 | w32-il | w32-il64 | GMP 1T |
|---|---|---|---|---|---|---|---|---|---|
| ADD | 118976 / 75542 | - | - | - | 92.11 M | - | - | - | 34.32 M |
| SUBTRACT | 118976 / 75542 | - | - | - | 731.54 M | - | - | - | 28.79 M |
| ADDMOD | 118976 / 75542 | - | - | - | 818.29 M | - | - | - | 19.79 M |
| SUBTRACTMOD | 118976 / 75542 | - | - | - | 839.64 M | - | - | - | 21.85 M |
| MULTIPLYOPERANDSCANNING | 118976 / 75542 | - | - | - | 127.05 M | - | - | - | 2.77 M |
| MULTIPLYPRODUCTSCANNING | 118976 / 75542 | - | - | - | 50.54 M | - | - | - | 2.76 M |
| MONTGOMERYMULTIPLICATION | 118976 / 75542 | - | - | - | 55.04 M | - | - | - | 407.23 k |
| COMPARE | 118976 / 75542 | - | - | - | 832.10 M | - | - | - | 188.92 M |
| REDUCE | 118976 / 9442 | - | - | - | 5.80 M | - | - | - | 51.74 M |
| MODMUL | 118976 / 4721 | - | - | - | 1.79 M | - | - | - | 938.23 k |
| MODEXP | 118976 / 1180 | - | - | - | 16.63 k | - | - | - | 594.1 |
| EXPONENTIATION | 118976 / 1180 | - | - | - | 39.64 k | - | - | - | 5.63 k |
| DIVIDE | 118976 / 9442 | - | - | - | 5.62 M | - | - | - | 32.38 M |
| ISQRT | 118976 / 2360 | - | - | - | 336.72 k | - | - | - | 3.74 M |
| MODMUL_R2 | 118976 / 75542 | - | - | - | 34.97 M | - | - | - | 932.59 k |

## 5. Head to head

Best OpenCL GPU result against best OpenCL CPU result and the CPU libraries.
Ratios above 1.00x mean the GPU is faster than that baseline.

### secp256k1 (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w32-opt | 6.43 G | none | n/a | 83.43 M | n/a |
| SUBTRACT | w32-opt | 6.14 G | none | n/a | 97.04 M | n/a |
| ADDMOD | w32-opt | 9.23 G | none | n/a | 37.00 M | n/a |
| SUBTRACTMOD | w32-opt | 8.29 G | none | n/a | 45.23 M | n/a |
| MULTIPLYOPERANDSCANNING | w32-opt | 2.95 G | none | n/a | 67.55 M | n/a |
| MULTIPLYPRODUCTSCANNING | w32-opt | 3.84 G | none | n/a | 66.87 M | n/a |
| MONTGOMERYMULTIPLICATION | w32-opt | 8.61 G | none | n/a | 9.94 M | n/a |
| COMPARE | w32-opt | 11.48 G | none | n/a | 116.64 M | n/a |
| REDUCE | w32-opt | 815.39 M | none | n/a | 93.03 M | n/a |
| MODMUL | w32-opt | 265.94 M | none | n/a | 17.27 M | n/a |
| MODEXP | w32-opt | 2.41 M | none | n/a | 160.96 k | n/a |
| EXPONENTIATION | w32-opt | 2.13 M | none | n/a | 512.00 k | n/a |
| DIVIDE | w32-opt | 692.56 M | none | n/a | 48.27 M | n/a |
| ISQRT | w32-opt | 23.14 M | none | n/a | 22.12 M | n/a |
| MODMUL_R2 | w32-opt | 5.19 G | none | n/a | 17.29 M | n/a |

### rsa256(composite) (256-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w32-opt | 7.80 G | none | n/a | 81.04 M | n/a |
| SUBTRACT | w32-opt | 7.44 G | none | n/a | 97.33 M | n/a |
| ADDMOD | w32-opt | 10.71 G | none | n/a | 42.20 M | n/a |
| SUBTRACTMOD | w32-opt | 10.07 G | none | n/a | 45.21 M | n/a |
| MULTIPLYOPERANDSCANNING | w32-opt | 3.10 G | none | n/a | 66.69 M | n/a |
| MULTIPLYPRODUCTSCANNING | w32-opt | 4.37 G | none | n/a | 66.98 M | n/a |
| MONTGOMERYMULTIPLICATION | w32-opt | 8.36 G | none | n/a | 9.86 M | n/a |
| COMPARE | w32-opt | 9.84 G | none | n/a | 117.46 M | n/a |
| REDUCE | w32-opt | 800.16 M | none | n/a | 57.37 M | n/a |
| MODMUL | w32-opt | 273.58 M | none | n/a | 17.22 M | n/a |
| MODEXP | w32-opt | 2.41 M | none | n/a | 172.72 k | n/a |
| EXPONENTIATION | w32-opt | 2.14 M | none | n/a | 515.95 k | n/a |
| DIVIDE | w32-opt | 654.84 M | none | n/a | 49.65 M | n/a |
| ISQRT | w32-opt | 22.89 M | none | n/a | 22.03 M | n/a |
| MODMUL_R2 | w32-opt | 5.82 G | none | n/a | 17.31 M | n/a |

### brainpoolP512r1 (512-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w32-opt | 2.77 G | none | n/a | 66.91 M | n/a |
| SUBTRACT | w32-opt | 2.84 G | none | n/a | 72.12 M | n/a |
| ADDMOD | w32-opt | 3.30 G | none | n/a | 37.44 M | n/a |
| SUBTRACTMOD | w32-opt | 3.20 G | none | n/a | 40.36 M | n/a |
| MULTIPLYOPERANDSCANNING | w32-opt | 559.76 M | none | n/a | 32.30 M | n/a |
| MULTIPLYPRODUCTSCANNING | w32-opt | 1.20 G | none | n/a | 31.95 M | n/a |
| MONTGOMERYMULTIPLICATION | w32-opt | 1.99 G | none | n/a | 4.23 M | n/a |
| COMPARE | w32-opt | 3.24 G | none | n/a | 85.31 M | n/a |
| REDUCE | w32-opt | 25.44 M | none | n/a | 53.55 M | n/a |
| MODMUL | w32-opt | 4.83 M | none | n/a | 8.77 M | n/a |
| MODEXP | w32-opt | 145.45 k | none | n/a | 28.30 k | n/a |
| EXPONENTIATION | w32-opt | 166.48 k | none | n/a | 160.40 k | n/a |
| DIVIDE | w32-opt | 26.89 M | none | n/a | 46.97 M | n/a |
| ISQRT | w32-opt | 558.49 k | none | n/a | 12.11 M | n/a |
| MODMUL_R2 | w32-opt | 1.70 G | none | n/a | 8.74 M | n/a |

### p1024 (1024-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w32-opt | 1.20 G | none | n/a | 47.56 M | n/a |
| SUBTRACT | w32-opt | 1.24 G | none | n/a | 52.12 M | n/a |
| ADDMOD | w32-opt | 1.35 G | none | n/a | 27.27 M | n/a |
| SUBTRACTMOD | w32-opt | 735.46 M | none | n/a | 32.50 M | n/a |
| MULTIPLYOPERANDSCANNING | w32-opt | 170.84 M | none | n/a | 8.99 M | n/a |
| MULTIPLYPRODUCTSCANNING | w32-opt | 182.57 M | none | n/a | 9.05 M | n/a |
| MONTGOMERYMULTIPLICATION | w32-opt | 344.14 M | none | n/a | 1.34 M | n/a |
| COMPARE | w32-opt | 1.45 G | none | n/a | 117.45 M | n/a |
| REDUCE | w32-opt | 3.43 M | none | n/a | 75.64 M | n/a |
| MODMUL | w32-opt | 582.50 k | none | n/a | 3.08 M | n/a |
| MODEXP | w32-opt | 8.98 k | none | n/a | 4.48 k | n/a |
| EXPONENTIATION | w32-opt | 8.36 k | none | n/a | 33.14 k | n/a |
| DIVIDE | w32-opt | 3.46 M | none | n/a | 40.56 M | n/a |
| ISQRT | w32-opt | 61.75 k | none | n/a | 6.21 M | n/a |
| MODMUL_R2 | w32-opt | 444.27 M | none | n/a | 3.07 M | n/a |

### p2048 (2048-bit)

| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GPU vs CPU-CL |
|---|---|---|---|---|---|---|
| ADD | w32-opt | 58.48 M | none | n/a | 34.32 M | n/a |
| SUBTRACT | w32-opt | 464.48 M | none | n/a | 28.79 M | n/a |
| ADDMOD | w32-opt | 519.56 M | none | n/a | 19.79 M | n/a |
| SUBTRACTMOD | w32-opt | 533.12 M | none | n/a | 21.85 M | n/a |
| MULTIPLYOPERANDSCANNING | w32-opt | 80.67 M | none | n/a | 2.77 M | n/a |
| MULTIPLYPRODUCTSCANNING | w32-opt | 32.09 M | none | n/a | 2.76 M | n/a |
| MONTGOMERYMULTIPLICATION | w32-opt | 34.95 M | none | n/a | 407.23 k | n/a |
| COMPARE | w32-opt | 528.33 M | none | n/a | 188.92 M | n/a |
| REDUCE | w32-opt | 460.35 k | none | n/a | 51.74 M | n/a |
| MODMUL | w32-opt | 70.86 k | none | n/a | 938.23 k | n/a |
| MODEXP | w32-opt | 164.9 | none | n/a | 594.1 | n/a |
| EXPONENTIATION | w32-opt | 393.1 | none | n/a | 5.63 k | n/a |
| DIVIDE | w32-opt | 446.23 k | none | n/a | 32.38 M | n/a |
| ISQRT | w32-opt | 6.68 k | none | n/a | 3.74 M | n/a |
| MODMUL_R2 | w32-opt | 22.20 M | none | n/a | 932.59 k | n/a |

## 6. Raw data

Also written to `NVIDIA_GeForce_RTX_5090_Report.csv` for analysis.

```csv
kind,device,device_type,kernel,modulus,bits,operation,items,seconds,ops_per_sec,mismatches
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,ADD,604343,0.007243788,83429138.379,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,ADD,604343,0.000512213,1179866087.454,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,ADD,604343,0.001753458,344657802.545,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,secp256k1,256,ADD,700000,0.000095840,7303839732.888,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,ADD,604288,0.000093928,6433518577.171,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,ADD,604288,0.016768235,36037662.847,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,SUBTRACT,604343,0.006227570,97043150.206,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,SUBTRACT,604343,0.000524696,1151796550.764,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,SUBTRACT,604343,0.001888835,319955444.050,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,secp256k1,256,SUBTRACT,700000,0.000096128,7281957390.146,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,SUBTRACT,604288,0.000098477,6136343576.421,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,SUBTRACT,604288,0.016748519,36080085.328,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,ADDMOD,604343,0.016335193,36996379.381,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,ADDMOD,604343,0.000886813,681477350.298,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,ADDMOD,604343,0.005235127,115439989.682,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,secp256k1,256,ADDMOD,700000,0.000090880,7702464788.732,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,ADDMOD,604288,0.000065454,9232247796.748,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,ADDMOD,604288,0.016701388,36181903.095,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,SUBTRACTMOD,604343,0.013362142,45228003.032,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,SUBTRACTMOD,604343,0.000739905,816784529.932,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,SUBTRACTMOD,604343,0.005200520,116208188.550,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,secp256k1,256,SUBTRACTMOD,700000,0.000096768,7233796296.296,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,SUBTRACTMOD,604288,0.000072899,8289381944.258,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,SUBTRACTMOD,604288,0.016672403,36244805.582,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,MULTIPLYOPERANDSCANNING,604343,0.008947222,67545324.553,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,604343,0.000619456,975602882.782,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,MULTIPLYOPERANDSCANNING,604343,0.001746626,346005935.795,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,MULTIPLYOPERANDSCANNING,604288,0.000204558,2954115434.587,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,MULTIPLYOPERANDSCANNING,604288,0.022350701,27036646.366,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,MULTIPLYPRODUCTSCANNING,604343,0.009037983,66867021.805,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,604343,0.000670192,901745979.721,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,MULTIPLYPRODUCTSCANNING,604343,0.002267874,266479962.729,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,secp256k1,256,MULTIPLYPRODUCTSCANNING,700000,0.000094528,7405213270.142,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,MULTIPLYPRODUCTSCANNING,604288,0.000157559,3835316749.221,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,MULTIPLYPRODUCTSCANNING,604288,0.022163776,27264668.280,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,MONTGOMERYMULTIPLICATION,604343,0.060822580,9936161.865,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,604343,0.003670777,164636262.666,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,MONTGOMERYMULTIPLICATION,604343,0.001952766,309480504.931,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,secp256k1,256,MONTGOMERYMULTIPLICATION,700000,0.000096992,7217090069.284,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,MONTGOMERYMULTIPLICATION,604288,0.000070224,8605171587.603,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,MONTGOMERYMULTIPLICATION,604288,0.016761503,36052136.460,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,COMPARE,604343,0.005181324,116638720.823,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,COMPARE,604343,0.000566255,1067262853.295,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,COMPARE,604343,0.001982693,304809187.345,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,secp256k1,256,COMPARE,700000,0.000091456,7653953813.856,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,COMPARE,604288,0.000052650,11477456318.778,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,COMPARE,604288,0.016700418,36184004.567,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,REDUCE,75542,0.000812042,93027203.497,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,REDUCE,75542,0.000131359,575080106.055,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,REDUCE,75542,0.000611802,123474568.524,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,secp256k1,256,REDUCE,700000,0.000092704,7550914739.386,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,REDUCE,118976,0.000092645,1284214726.076,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,REDUCE,118976,0.003545078,33560897.831,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,MODMUL,37771,0.002186920,17271322.002,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,MODMUL,37771,0.000107864,350172705.748,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,MODMUL,37771,0.000451898,83582983.530,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,secp256k1,256,MODMUL,700000,0.000221216,3164328077.535,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,MODMUL,118976,0.000142030,837682461.401,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,MODMUL,118976,0.003647914,32614804.651,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,MODEXP,9442,0.058660877,160959.066,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,MODEXP,9442,0.002893572,3263094.916,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,MODEXP,9442,0.005406512,1746412.491,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,secp256k1,256,MODEXP,700000,0.063936286,10948399.474,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,MODEXP,118976,0.003917715,30368722.224,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,MODEXP,118976,0.007462514,15943152.459,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,EXPONENTIATION,9442,0.018441250,512004.336,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,EXPONENTIATION,9442,0.001718602,5493999.898,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,EXPONENTIATION,9442,0.013892346,679654.834,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,EXPONENTIATION,118976,0.004424568,26889856.618,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,EXPONENTIATION,118976,0.007926064,15010729.012,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,DIVIDE,75542,0.001565111,48266229.225,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,DIVIDE,75542,0.000091764,823221226.557,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,DIVIDE,75542,0.000237941,317482317.568,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,secp256k1,256,DIVIDE,700000,0.000090560,7729681978.799,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,DIVIDE,118976,0.000109077,1090752748.257,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,DIVIDE,118976,0.004597115,25880580.212,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,ISQRT,18885,0.000853881,22116664.786,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,ISQRT,18885,0.000051978,363326483.063,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,ISQRT,118976,0.000816059,145793409.944,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,ISQRT,118976,0.002668976,44577394.551,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,secp256k1,256,MODMUL_R2,604343,0.034959221,17287084.273,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,secp256k1,256,MODMUL_R2,604343,0.001717590,351855209.916,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,secp256k1,256,MODMUL_R2,604343,0.004965403,121710762.597,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,secp256k1,256,MODMUL_R2,700000,0.000106624,6565126050.420,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,MODMUL_R2,604288,0.000116350,5193713261.584,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,secp256k1,256,MODMUL_R2,604288,0.008469777,71346388.889,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,ADD,604343,0.007457244,81041065.822,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,ADD,604343,0.000682906,884957793.405,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,ADD,604343,0.001804324,334941488.654,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,rsa256(composite),256,ADD,700000,0.000085824,8156226696.495,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,ADD,604288,0.000077457,7801599751.605,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,ADD,604288,0.008422466,71747158.179,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACT,604343,0.006209074,97332229.800,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACT,604343,0.000578738,1044242726.811,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACT,604343,0.001757606,343844424.266,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,rsa256(composite),256,SUBTRACT,700000,0.000096096,7284382284.382,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,SUBTRACT,604288,0.000081264,7436114207.438,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,SUBTRACT,604288,0.008448365,71527213.410,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,ADDMOD,604343,0.014322273,42196025.824,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,ADDMOD,604343,0.000828623,729334228.524,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,ADDMOD,604343,0.002104905,287111774.215,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,rsa256(composite),256,ADDMOD,700000,0.000090208,7759843916.282,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,ADDMOD,604288,0.000056447,10705405907.346,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,ADDMOD,604288,0.008417136,71792590.428,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,SUBTRACTMOD,604343,0.013366509,45213226.593,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,SUBTRACTMOD,604343,0.000782986,771843771.913,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,SUBTRACTMOD,604343,0.002289976,263908006.791,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,rsa256(composite),256,SUBTRACTMOD,700000,0.000092896,7535308301.757,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,SUBTRACTMOD,604288,0.000059994,10072484388.545,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,SUBTRACTMOD,604288,0.016691020,36204378.181,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYOPERANDSCANNING,604343,0.009061949,66690179.917,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,604343,0.000555734,1087467811.504,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,604343,0.001824062,331317164.138,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,604288,0.000194700,3103688216.597,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,MULTIPLYOPERANDSCANNING,604288,0.022284656,27116774.800,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,604343,0.009022124,66984560.190,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,604343,0.000610700,989590832.779,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,604343,0.001967875,307104355.401,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,700000,0.000093920,7453151618.399,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,604288,0.000138343,4368046607.860,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,MULTIPLYPRODUCTSCANNING,604288,0.022171492,27255179.772,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MONTGOMERYMULTIPLICATION,604343,0.061276711,9862523.454,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,604343,0.003670797,164635364.609,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,604343,0.001920546,314672504.112,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,rsa256(composite),256,MONTGOMERYMULTIPLICATION,700000,0.000089792,7795794725.588,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,604288,0.000072257,8363022721.273,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,MONTGOMERYMULTIPLICATION,604288,0.016733239,36113032.402,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,COMPARE,604343,0.005145135,117459117.258,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,COMPARE,604343,0.000536919,1125576004.919,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,COMPARE,604343,0.001775510,340377148.205,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,rsa256(composite),256,COMPARE,700000,0.000090432,7740622788.393,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,COMPARE,604288,0.000061436,9836042526.757,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,COMPARE,604288,0.016701580,36181486.966,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,REDUCE,75542,0.001316809,57367465.232,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,REDUCE,75542,0.000206182,366385016.543,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,REDUCE,75542,0.000232160,325387861.905,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,rsa256(composite),256,REDUCE,700000,0.000084288,8304859529.233,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,REDUCE,118976,0.000094409,1260217490.724,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,REDUCE,118976,0.003561189,33409067.099,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MODMUL,37771,0.002193373,17220509.788,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MODMUL,37771,0.000107123,352594513.546,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MODMUL,37771,0.000326219,115784176.531,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,rsa256(composite),256,MODMUL,700000,0.000218752,3199970743.125,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,MODMUL,118976,0.000138062,861757433.756,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,MODMUL,118976,0.003664085,32470860.881,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MODEXP,9442,0.054666853,172718.924,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MODEXP,9442,0.002784626,3390760.725,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MODEXP,9442,0.005422090,1741394.865,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,rsa256(composite),256,MODEXP,700000,0.062912092,11126636.832,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,MODEXP,118976,0.003913478,30401601.892,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,MODEXP,118976,0.007425834,16021904.160,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,EXPONENTIATION,9442,0.018300052,515954.817,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,EXPONENTIATION,9442,0.000954121,9896018.249,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,EXPONENTIATION,9442,0.014049425,672055.984,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,EXPONENTIATION,118976,0.004412034,26966246.942,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,EXPONENTIATION,118976,0.007749008,15353707.181,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,DIVIDE,75542,0.001521568,49647467.397,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,DIVIDE,75542,0.000092025,820884310.932,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,DIVIDE,75542,0.000247600,305096902.558,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,rsa256(composite),256,DIVIDE,700000,0.000089856,7790242165.242,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,DIVIDE,118976,0.000115359,1031352487.787,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,DIVIDE,118976,0.004605692,25832382.427,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,ISQRT,18885,0.000857397,22025969.747,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,ISQRT,18885,0.000060084,314309740.737,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,ISQRT,118976,0.000825016,144210525.888,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,ISQRT,118976,0.004327594,27492413.679,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,rsa256(composite),256,MODMUL_R2,604343,0.034914755,17309100.393,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,rsa256(composite),256,MODMUL_R2,604343,0.001702361,355002885.909,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,rsa256(composite),256,MODMUL_R2,604343,0.004957659,121900882.585,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,rsa256(composite),256,MODMUL_R2,700000,0.000101280,6911532385.466,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,MODMUL_R2,604288,0.000103827,5820150554.826,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,rsa256(composite),256,MODMUL_R2,604288,0.008733547,69191589.152,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,ADD,302171,0.004516341,66906152.199,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,ADD,302171,0.000176886,1708282355.449,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,ADD,302171,0.000470934,641641969.381,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,brainpoolP512r1,512,ADD,700000,0.000146560,4776200873.362,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,ADD,302144,0.000109107,2769244850.596,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,ADD,302144,0.016721096,18069629.048,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACT,302171,0.004189952,72118013.197,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACT,302171,0.000162348,1861255771.075,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACT,302171,0.000446267,677108159.119,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,brainpoolP512r1,512,SUBTRACT,700000,0.000148640,4709364908.504,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,SUBTRACT,302144,0.000106242,2843921839.035,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,SUBTRACT,302144,0.016726977,18063275.661,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,ADDMOD,302171,0.008069897,37444219.841,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,ADDMOD,302171,0.000419466,720370962.349,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,ADDMOD,302171,0.002429060,124398330.025,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,brainpoolP512r1,512,ADDMOD,700000,0.000147872,4733823847.652,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,ADDMOD,302144,0.000091664,3296211267.336,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,ADDMOD,302144,0.016753427,18034757.919,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,SUBTRACTMOD,302171,0.007487250,40358075.253,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,SUBTRACTMOD,302171,0.000402264,751175724.875,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,SUBTRACTMOD,302171,0.001391321,217182813.866,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,brainpoolP512r1,512,SUBTRACTMOD,700000,0.000142688,4905808477.237,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,SUBTRACTMOD,302144,0.000094569,3194956308.261,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,SUBTRACTMOD,302144,0.016718962,18071935.330,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,302171,0.009355115,32300083.945,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,302171,0.000569811,530300313.322,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,302171,0.000964240,313377395.017,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,302144,0.000539825,559707486.503,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYOPERANDSCANNING,302144,0.022647864,13340949.091,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,302171,0.009458612,31946653.645,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,302171,0.000568358,531656108.958,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,302171,0.000621099,486510182.486,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,700000,0.000153184,4569667850.428,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,302144,0.000251527,1201239099.030,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,MULTIPLYPRODUCTSCANNING,302144,0.022372011,13505446.547,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,302171,0.071442822,4229550.171,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,302171,0.003531553,85563204.859,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,302171,0.001009266,299396752.258,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,700000,0.002374688,294775566.306,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,302144,0.000152028,1987424246.208,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,MONTGOMERYMULTIPLICATION,302144,0.016906397,17871578.238,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,COMPARE,302171,0.003542133,85307634.678,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,COMPARE,302171,0.000124246,2432040362.646,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,COMPARE,302171,0.000427722,706465882.729,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,brainpoolP512r1,512,COMPARE,700000,0.000247200,2831715210.356,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,COMPARE,302144,0.000093187,3242341621.004,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,COMPARE,302144,0.016751614,18036709.848,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,REDUCE,37771,0.000705288,53554007.273,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,REDUCE,37771,0.000113054,334097140.304,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,REDUCE,37771,0.000311762,121153369.428,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,brainpoolP512r1,512,REDUCE,700000,0.000243744,2871865563.870,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,REDUCE,118976,0.001484989,80119114.379,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,REDUCE,118976,0.008081207,14722553.108,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL,18885,0.002154329,8766070.170,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL,18885,0.000109518,172437491.252,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL,18885,0.000268039,70456148.991,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,brainpoolP512r1,512,MODMUL,700000,0.002380640,294038577.861,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,MODMUL,118976,0.003907386,30449001.094,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,MODMUL,118976,0.010526939,11302050.795,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MODEXP,4721,0.166833075,28297.746,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MODEXP,4721,0.009314479,506845.306,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MODEXP,4721,0.008417596,560848.973,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,brainpoolP512r1,512,MODEXP,700000,0.211561188,3308735.438,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,MODEXP,118976,0.032458494,3665481.220,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,MODEXP,118976,0.034562827,3442311.004,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,EXPONENTIATION,4721,0.029433283,160396.650,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,EXPONENTIATION,4721,0.001703844,2770793.668,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,EXPONENTIATION,4721,0.019961876,236500.819,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,EXPONENTIATION,118976,0.028357591,4195560.897,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,EXPONENTIATION,118976,0.030788164,3864342.139,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,DIVIDE,37771,0.000804136,46970921.566,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,DIVIDE,37771,0.000062088,608346057.380,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,DIVIDE,37771,0.000162969,231768178.884,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,brainpoolP512r1,512,DIVIDE,700000,0.002377216,294462093.474,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,DIVIDE,118976,0.001404776,84693929.545,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,DIVIDE,118976,0.003504732,33947246.905,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,ISQRT,9442,0.000779369,12114926.579,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,ISQRT,9442,0.000063119,149590478.612,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,ISQRT,118976,0.016906168,7037431.610,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,ISQRT,118976,0.018563442,6409156.256,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,brainpoolP512r1,512,MODMUL_R2,302171,0.034582265,8737744.574,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,brainpoolP512r1,512,MODMUL_R2,302171,0.001680940,179763116.867,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,brainpoolP512r1,512,MODMUL_R2,302171,0.003911905,77243951.660,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,brainpoolP512r1,512,MODMUL_R2,700000,0.000170240,4111842105.263,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,MODMUL_R2,302144,0.000177667,1700619597.421,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,brainpoolP512r1,512,MODMUL_R2,302144,0.004131201,73137086.093,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p1024,1024,ADD,151085,0.003176740,47559763.946,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p1024,1024,ADD,151085,0.000128504,1175723173.496,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p1024,1024,ADD,151085,0.000148922,1014522985.535,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p1024,1024,ADD,700000,0.000460992,1518464528.669,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,ADD,151040,0.000125649,1202079577.043,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,ADD,151040,0.016752756,9015830.010,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p1024,1024,SUBTRACT,151085,0.002898842,52119085.821,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p1024,1024,SUBTRACT,151085,0.000125678,1202159263.845,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p1024,1024,SUBTRACT,151085,0.000148041,1020562278.209,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p1024,1024,SUBTRACT,700000,0.004714080,148491328.106,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,SUBTRACT,151040,0.000122042,1237607667.014,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,SUBTRACT,151040,0.016821947,8978746.661,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p1024,1024,ADDMOD,151085,0.005540195,27270700.002,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p1024,1024,ADDMOD,151085,0.000303036,498571015.586,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p1024,1024,ADDMOD,151085,0.000965642,156460668.721,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p1024,1024,ADDMOD,700000,0.000246464,2840171384.056,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,ADDMOD,151040,0.000112253,1345531867.154,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,ADDMOD,151040,0.016750541,9017022.170,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p1024,1024,SUBTRACTMOD,151085,0.004649204,32496960.986,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p1024,1024,SUBTRACTMOD,151085,0.000249693,605083193.774,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p1024,1024,SUBTRACTMOD,151085,0.000882926,171118544.187,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p1024,1024,SUBTRACTMOD,700000,0.000246464,2840171384.056,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,SUBTRACTMOD,151040,0.000205430,735238390.357,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,SUBTRACTMOD,151040,0.016750932,9016811.736,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p1024,1024,MULTIPLYOPERANDSCANNING,151085,0.016797009,8994756.190,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p1024,1024,MULTIPLYOPERANDSCANNING,151085,0.001194017,126535048.491,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p1024,1024,MULTIPLYOPERANDSCANNING,151085,0.001043250,144821491.776,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,MULTIPLYOPERANDSCANNING,151040,0.000884389,170784568.010,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,MULTIPLYOPERANDSCANNING,151040,0.022952362,6580586.357,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p1024,1024,MULTIPLYPRODUCTSCANNING,151085,0.016693833,9050348.086,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,151085,0.000985931,153240961.397,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p1024,1024,MULTIPLYPRODUCTSCANNING,151085,0.001049362,143977943.063,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p1024,1024,MULTIPLYPRODUCTSCANNING,700000,0.000449600,1556939501.779,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,MULTIPLYPRODUCTSCANNING,151040,0.000827561,182512247.788,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,MULTIPLYPRODUCTSCANNING,151040,0.022983310,6571725.295,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p1024,1024,MONTGOMERYMULTIPLICATION,151085,0.113068937,1336220.223,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p1024,1024,MONTGOMERYMULTIPLICATION,151085,0.005766735,26199400.770,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p1024,1024,MONTGOMERYMULTIPLICATION,151085,0.001737878,86936484.339,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p1024,1024,MONTGOMERYMULTIPLICATION,700000,0.008761887,79891466.302,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,MONTGOMERYMULTIPLICATION,151040,0.000439023,344036722.949,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,MONTGOMERYMULTIPLICATION,151040,0.017274847,8743348.049,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p1024,1024,COMPARE,151085,0.001286421,117446002.038,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p1024,1024,COMPARE,151085,0.000076445,1976386898.252,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p1024,1024,COMPARE,151085,0.000174251,867054248.808,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p1024,1024,COMPARE,700000,0.000458752,1525878906.250,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,COMPARE,151040,0.000103887,1453885010.013,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,COMPARE,151040,0.016676732,9056930.281,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p1024,1024,REDUCE,18885,0.000249654,75644675.698,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p1024,1024,REDUCE,18885,0.000023245,812428832.847,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p1024,1024,REDUCE,18885,0.000133183,141797284.670,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p1024,1024,REDUCE,700000,0.000454048,1541687222.496,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,REDUCE,118976,0.005510459,21590942.054,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,REDUCE,118976,0.018767789,6339372.206,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p1024,1024,MODMUL,9442,0.003065378,3080207.227,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p1024,1024,MODMUL,9442,0.000154082,61278980.578,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p1024,1024,MODMUL,9442,0.000339324,27825902.821,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p1024,1024,MODMUL,700000,0.004740128,147675337.037,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,MODMUL,118976,0.016209344,7339963.951,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,MODMUL,118976,0.019710871,6036060.022,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p1024,1024,MODEXP,2360,0.526814813,4479.753,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p1024,1024,MODEXP,2360,0.032843103,71856.792,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p1024,1024,MODEXP,2360,0.025942618,90970.001,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p1024,1024,MODEXP,700000,1.212423444,577356.041,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,MODEXP,118976,0.262765911,452783.238,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,MODEXP,118976,0.266055824,447184.347,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p1024,1024,EXPONENTIATION,2360,0.071207857,33142.410,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p1024,1024,EXPONENTIATION,2360,0.004487437,525912.680,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p1024,1024,EXPONENTIATION,2360,0.066702937,35380.751,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,EXPONENTIATION,118976,0.282291729,421464.704,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,EXPONENTIATION,118976,0.286294348,415572.298,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p1024,1024,DIVIDE,18885,0.000465584,40561961.978,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p1024,1024,DIVIDE,18885,0.000028525,662050519.818,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p1024,1024,DIVIDE,18885,0.000100521,187870933.292,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p1024,1024,DIVIDE,700000,0.004723648,148190551.032,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,DIVIDE,118976,0.005461115,21786026.585,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,DIVIDE,118976,0.010013595,11881447.184,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p1024,1024,ISQRT,4721,0.000760083,6211163.072,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p1024,1024,ISQRT,4721,0.000052881,89275986.494,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,ISQRT,118976,0.076455969,1556137.496,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,ISQRT,118976,0.080119221,1484986.979,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p1024,1024,MODMUL_R2,151085,0.049265544,3066747.831,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p1024,1024,MODMUL_R2,151085,0.002529682,59724900.753,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p1024,1024,MODMUL_R2,151085,0.007007972,21559018.465,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p1024,1024,MODMUL_R2,700000,0.000556064,1258847902.400,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,MODMUL_R2,151040,0.000340075,444137323.619,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p1024,1024,MODMUL_R2,151040,0.004405221,34286588.802,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p2048,2048,ADD,75542,0.002201347,34316263.167,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p2048,2048,ADD,75542,0.000084511,893871757.255,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p2048,2048,ADD,75542,0.000293046,257782051.873,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p2048,2048,ADD,700000,0.000459456,1523540883.131,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,ADD,118976,0.001291732,92105787.782,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,ADD,118976,0.011058698,10758590.036,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p2048,2048,SUBTRACT,75542,0.002624180,28786897.952,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p2048,2048,SUBTRACT,75542,0.000112734,670090438.805,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p2048,2048,SUBTRACT,75542,0.000132462,570291815.149,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p2048,2048,SUBTRACT,700000,0.000459776,1522480512.249,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,SUBTRACT,118976,0.000162638,731538640.720,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,SUBTRACT,118976,0.006760020,17599947.820,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p2048,2048,ADDMOD,75542,0.003816644,19792781.496,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p2048,2048,ADDMOD,75542,0.000193337,390727115.758,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p2048,2048,ADDMOD,75542,0.000678488,111338770.884,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p2048,2048,ADDMOD,700000,0.000464160,1508100654.947,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,ADDMOD,118976,0.000145396,818289446.602,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,ADDMOD,118976,0.006714884,17718251.086,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p2048,2048,SUBTRACTMOD,75542,0.003457012,21851819.563,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p2048,2048,SUBTRACTMOD,75542,0.000196934,383590482.043,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p2048,2048,SUBTRACTMOD,75542,0.000665754,113468364.234,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p2048,2048,SUBTRACTMOD,700000,0.009610528,72836788.988,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,SUBTRACTMOD,118976,0.000141699,839639740.596,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,SUBTRACTMOD,118976,0.006676331,17820566.604,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p2048,2048,MULTIPLYOPERANDSCANNING,75542,0.027292249,2767892.083,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p2048,2048,MULTIPLYOPERANDSCANNING,75542,0.001554520,48595068.798,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p2048,2048,MULTIPLYOPERANDSCANNING,75542,0.001692433,44635151.884,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,MULTIPLYOPERANDSCANNING,118976,0.000936448,127050317.592,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,MULTIPLYOPERANDSCANNING,118976,0.009719747,12240647.971,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p2048,2048,MULTIPLYPRODUCTSCANNING,75542,0.027321796,2764898.748,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,75542,0.001550463,48722219.507,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p2048,2048,MULTIPLYPRODUCTSCANNING,75542,0.001715295,44040237.546,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p2048,2048,MULTIPLYPRODUCTSCANNING,700000,0.009617024,72787590.007,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,MULTIPLYPRODUCTSCANNING,118976,0.002354118,50539521.763,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,MULTIPLYPRODUCTSCANNING,118976,0.011336365,10495074.905,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p2048,2048,MONTGOMERYMULTIPLICATION,75542,0.185501233,407231.795,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p2048,2048,MONTGOMERYMULTIPLICATION,75542,0.009479762,7968765.414,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p2048,2048,MONTGOMERYMULTIPLICATION,75542,0.003324270,22724387.503,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p2048,2048,MONTGOMERYMULTIPLICATION,700000,0.000977280,716273739.358,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,MONTGOMERYMULTIPLICATION,118976,0.002161613,55040380.408,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,MONTGOMERYMULTIPLICATION,118976,0.015446616,7702399.088,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p2048,2048,COMPARE,75542,0.000399869,188916790.820,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p2048,2048,COMPARE,75542,0.000031089,2429848790.690,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p2048,2048,COMPARE,75542,0.000072318,1044582375.168,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p2048,2048,COMPARE,700000,0.000457248,1530897893.484,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,COMPARE,118976,0.000142982,832104759.314,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,COMPARE,118976,0.013619188,8735909.962,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p2048,2048,REDUCE,9442,0.000182476,51743820.827,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p2048,2048,REDUCE,9442,0.000018956,498097194.763,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p2048,2048,REDUCE,9442,0.000116651,80942335.131,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p2048,2048,REDUCE,700000,0.000461600,1516464471.404,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,REDUCE,118976,0.020510637,5800697.406,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,REDUCE,118976,0.027297538,4358488.299,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p2048,2048,MODMUL,4721,0.005031831,938227.041,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p2048,2048,MODMUL,4721,0.000245666,19217138.209,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p2048,2048,MODMUL,4721,0.000564020,8370272.127,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p2048,2048,MODMUL,700000,0.009618752,72774513.783,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,MODMUL,118976,0.066627913,1785678.019,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,MODMUL,118976,0.073721037,1613867.694,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p2048,2048,MODEXP,1180,1.986313961,594.065,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p2048,2048,MODEXP,1180,0.119374885,9884.826,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p2048,2048,MODEXP,1180,0.092965990,12692.814,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p2048,2048,MODEXP,700000,6.987298489,100181.780,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,MODEXP,118976,7.153930234,16630.858,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,MODEXP,118976,7.165998883,16602.849,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p2048,2048,EXPONENTIATION,1180,0.209555117,5630.977,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p2048,2048,EXPONENTIATION,1180,0.014762435,79932.613,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p2048,2048,EXPONENTIATION,1180,0.172788786,6829.147,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,EXPONENTIATION,118976,3.001722644,39635.907,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,EXPONENTIATION,118976,3.008412902,39547.763,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p2048,2048,DIVIDE,9442,0.000291614,32378433.377,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p2048,2048,DIVIDE,9442,0.000027252,346468352.105,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p2048,2048,DIVIDE,9442,0.005703445,1655490.723,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p2048,2048,DIVIDE,700000,0.000884864,791082019.384,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,DIVIDE,118976,0.021159369,5622851.977,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,DIVIDE,118976,0.030520655,3898212.541,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p2048,2048,ISQRT,2360,0.000631409,3737672.419,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p2048,2048,ISQRT,2360,0.006681521,353212.991,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,ISQRT,118976,0.353343094,336715.227,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,ISQRT,118976,0.358985864,331422.521,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-1t,p2048,2048,MODMUL_R2,75542,0.081002801,932585.035,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,gmp-nt,p2048,2048,MODMUL_R2,75542,0.004120490,18333256.268,0
library,AMD Ryzen 9 7950X 16-Core Processor,host-cpu,openssl-nt,p2048,2048,MODMUL_R2,75542,0.008865216,8521168.530,0
library,NVIDIA GeForce RTX 5090,gpu,cgbn,p2048,2048,MODMUL_R2,700000,0.001894976,369397818.231,0
opencl-kernel,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,MODMUL_R2,118976,0.003402397,34968288.157,0
opencl-e2e,NVIDIA GeForce RTX 5090,GPU,w32-opt,p2048,2048,MODMUL_R2,118976,0.010042179,11847627.722,0
```
