# Benchmark reports

One `<Device>_Report.md` (human-readable) and `<Device>_Report.csv`
(machine-readable) per device, written by `GPU_Host`. Every configuration in
every report was verified word-for-word against GMP before it was timed; a
configuration that fails verification is never timed and never appears.

Reproduce a report with:

```sh
./GPU_Host --devices gpu --budget 10800
```

`GPU_Host` writes into the working directory, so move the resulting pair into
this folder afterwards.

## Device index

`GPU_Host` names each report from the device name the OpenCL runtime reports,
which is not always the marketing name. ROCm, for instance, reports the LLVM
target triple, so the AMD card appears as `gfx1010:xnack-` rather than by model.
This index identifies every report in the directory.

| Device                        | Vendor | Architecture               | CU  | Memory     | Runtime  | Report file                              |
|-------------------------------|--------|----------------------------|-----|------------|----------|------------------------------------------|
| Apple M2                      | Apple  | Apple silicon (M2)         | 10  | 11.84 GiB  | Apple    | Apple_M2_Report.md                       |
| Intel(R) Graphics             | Intel  | Intel Xe (Core Ultra iGPU) | 64  | 28.01 GiB  | NEO      | Intel_R_Graphics_Report.md               |
| NVIDIA A100-SXM4-40GB         | NVIDIA | Ampere                     | 108 | 39.49 GiB  | CUDA ICD | NVIDIA_A100-SXM4-40GB_Report.md          |
| NVIDIA A40                    | NVIDIA | Ampere                     | 84  | 44.43 GiB  | CUDA ICD | NVIDIA_A40_Report.md                     |
| NVIDIA B200                   | NVIDIA | Blackwell                  | 148 | 178.34 GiB | CUDA ICD | NVIDIA_B200_Report.md                    |
| NVIDIA B300 SXM6 AC           | NVIDIA | Blackwell                  | 148 | 267.68 GiB | CUDA ICD | NVIDIA_B300_SXM6_AC_Report.md            |
| NVIDIA GB10                   | NVIDIA | Grace-Blackwell            | 48  | 121.63 GiB | CUDA ICD | NVIDIA_GB10_Report.md                    |
| NVIDIA GeForce GTX 1060 3GB   | NVIDIA | Pascal                     | 9   | 2.94 GiB   | CUDA ICD | NVIDIA_GeForce_GTX_1060_3GB_Report.md    |
| NVIDIA GeForce GTX 1660 SUPER | NVIDIA | Turing                     | 22  | 5.61 GiB   | CUDA ICD | NVIDIA_GeForce_GTX_1660_SUPER_Report.md  |
| NVIDIA GeForce RTX 2060 SUPER | NVIDIA | Turing                     | 34  | 7.60 GiB   | CUDA ICD | NVIDIA_GeForce_RTX_2060_SUPER_Report.md  |
| NVIDIA GeForce RTX 2070       | NVIDIA | Turing                     | 36  | 7.60 GiB   | CUDA ICD | NVIDIA_GeForce_RTX_2070_Report.md        |
| NVIDIA GeForce RTX 3060       | NVIDIA | Ampere                     | 28  | 11.63 GiB  | CUDA ICD | NVIDIA_GeForce_RTX_3060_Report.md        |
| NVIDIA GeForce RTX 3080 Ti    | NVIDIA | Ampere                     | 80  | 11.63 GiB  | CUDA ICD | NVIDIA_GeForce_RTX_3080_Ti_Report.md     |
| NVIDIA GeForce RTX 4070 Ti    | NVIDIA | Ada Lovelace               | 60  | 11.61 GiB  | CUDA ICD | NVIDIA_GeForce_RTX_4070_Ti_Report.md     |
| NVIDIA GeForce RTX 5080       | NVIDIA | Blackwell                  | 84  | 15.45 GiB  | CUDA ICD | NVIDIA_GeForce_RTX_5080_Report.md        |
| NVIDIA GeForce RTX 5090       | NVIDIA | Blackwell                  | 170 | 31.84 GiB  | PoCL     | NVIDIA_GeForce_RTX_5090_Report.md        |
| NVIDIA GeForce RTX 5090       | NVIDIA | Blackwell                  | 170 | 31.84 GiB  | PoCL     | NVIDIA_GeForce_RTX_5090_Report_sm120.md  |
| NVIDIA GeForce RTX 5090       | NVIDIA | Blackwell                  | 170 | 31.84 GiB  | PoCL     | NVIDIA_GeForce_RTX_5090_Report_sm120b.md |
| NVIDIA GeForce RTX 5090       | NVIDIA | Blackwell                  | 170 | 31.84 GiB  | PoCL     | NVIDIA_GeForce_RTX_5090_Report_sm75.md   |
| NVIDIA H100 80GB HBM3         | NVIDIA | Hopper                     | 132 | 79.19 GiB  | CUDA ICD | NVIDIA_H100_80GB_HBM3_Report.md          |
| NVIDIA H100 80GB HBM3         | NVIDIA | Hopper                     | 132 | 79.19 GiB  | CUDA ICD | NVIDIA_H100_80GB_HBM3_Report_rerun.md    |
| NVIDIA H100 NVL               | NVIDIA | Hopper                     | 132 | 93.09 GiB  | CUDA ICD | NVIDIA_H100_NVL_Report.md                |
| NVIDIA H200 NVL               | NVIDIA | Hopper                     | 132 | 139.80 GiB | CUDA ICD | NVIDIA_H200_NVL_Report.md                |
| NVIDIA RTX A2000              | NVIDIA | Ampere                     | 26  | 5.66 GiB   | CUDA ICD | NVIDIA_RTX_A2000_Report.md               |
| NVIDIA TITAN Xp               | NVIDIA | Pascal                     | 30  | 11.89 GiB  | CUDA ICD | NVIDIA_TITAN_Xp_Report.md                |
| Quadro RTX 6000               | NVIDIA | Turing                     | 72  | 21.97 GiB  | CUDA ICD | Quadro_RTX_6000_Report.md                |
| Quadro RTX 8000               | NVIDIA | Turing                     | 72  | 47.27 GiB  | CUDA ICD | Quadro_RTX_8000_Report.md                |
| Tesla P40                     | NVIDIA | Pascal                     | 30  | 23.87 GiB  | CUDA ICD | Tesla_P40_Report.md                      |
| Tesla V100-PCIE-32GB          | NVIDIA | Volta                      | 80  | 31.73 GiB  | CUDA ICD | Tesla_V100-PCIE-32GB_Report.md           |
| gfx1010:xnack-                | AMD    | RDNA1 (Navi 10)            | 20  | 7.98 GiB   | ROCm     | gfx1010_xnack-_Report.md                 |

`CU` is the compute-unit count as reported by the runtime; RDNA runtimes report
work-group processors, so the AMD figure corresponds to roughly twice as many
physical compute units. Reports with a `_rerun`, `_new`, `_sm75` or `_sm120`
suffix are repeat runs of the same card, kept because they are the evidence for
the reproducibility figure and for the PoCL-versus-vendor-ICD comparison.


## Available results table

`ok` means the column is present and valid. `-` means the column is not present.

`cells` is the number of configurations that were verified against GMP and then
timed on that device. `vars` counts how many of the seven build-time variants
completed the *entire* cross-product of operations and moduli, so it is a
coverage figure, not a pass rate: a variant that ran 40 configurations
correctly out of 75 contributes 0. Two different things reduce it. On most
devices the time budget ran out, which is why the widest and slowest cells are
missing. 

```
report                           vars  cells  MPA  CPU  CGBN
Apple_M2                           7/7    485  ok    -     -
NVIDIA_A100-SXM4-40GB              7/7    485  ok    -    ok
NVIDIA_A40                         7/7    485  ok    -    ok
NVIDIA_B200                        7/7    485  ok   ok    ok
NVIDIA_B300_SXM6_AC                7/7    485  ok    -    ok
NVIDIA_GB10                        7/7    485  ok    -     -
NVIDIA_GeForce_GTX_1060_3GB        7/7    485  ok    -     -
NVIDIA_GeForce_GTX_1660_SUPER      7/7    485  ok    -     -
NVIDIA_GeForce_RTX_2060_SUPER      7/7    485  ok    -     -
NVIDIA_GeForce_RTX_2070            7/7    485  ok   ok     -
NVIDIA_GeForce_RTX_3060            7/7    485  ok    -     -
NVIDIA_GeForce_RTX_5080            7/7    485  ok    -    ok
NVIDIA_H100_80GB_HBM3              7/7    485  ok   ok    ok
NVIDIA_H200_NVL                    7/7    485  ok    -    ok
NVIDIA_RTX_A2000                   7/7    485  ok    -     -
gfx1010_xnack-                     7/7    485  ok    -     -
NVIDIA_GeForce_RTX_3080_Ti         5/7    335  ok    -     -
NVIDIA_GeForce_RTX_4070_Ti         5/7    335  ok    -     -
NVIDIA_H100_NVL                    4/7    300  ok   ok    ok
NVIDIA_TITAN_Xp                    4/7    332  ok    -     -
Tesla_V100-PCIE-32GB               4/7    331  ok    -    ok
Quadro_RTX_8000                    2/7    246  ok    -    ok
NVIDIA_GeForce_RTX_5090            1/7     75  ok    -     -
Tesla_P40                          1/7    149  ok    -     -
Intel_R_Graphics                   0/7     40  ok    -     -
Quadro_RTX_6000                    0/7    147  ok    -    ok

every column valid (MPA + CPU + CGBN): NVIDIA_B200, NVIDIA_H100_80GB_HBM3, NVIDIA_H100_NVL
...of those, with all seven variants: NVIDIA_B200, NVIDIA_H100_80GB_HBM3
MPA data usable: 26/26 reports
```

## Which reports the paper's figures come from

- CGBN comparison (10.8x, 3.25x, 1.64x, 0.97x at 256, 512, 1024 and 2048 bits)
  and the CPU comparison (502x single-threaded GMP, 40x across the host's 40
  threads): `NVIDIA_H100_NVL_Report.md`. The sweep and CGBN both ran at 500000
  items, with `cgbn/cgbn_results_NVIDIA_H100_NVL.tsv` written by the same run.
  All columns valid; covers the four `w32-*` variants.
- Seven-variant coverage and the interleaved-layout result: the reports listed
  as 7/7 in the table above.
- Portability and correctness across 26 devices: every report, since the MPA
  measurements are valid throughout.

## `cgbn/`

CGBN reference results, one file per device, as produced by `cgbn_bench.cu`. The
device and item count are recorded in each file's header, and `GPU_Host` refuses
a file whose device does not match the report it is building.

## `ablation/`

Controls behind specific claims, each row carrying its own `device` column:

- `order_mulhi_first.csv`, `order_mulhi_second.csv` -- the same cells with the
  `mul_hi` variants run in both orders, to show the effect is not an artifact of
  run order, clock ramp or thermal drift.
- `reduce_old_GeForce_RTX_5060.csv`, `reduce_new_GeForce_RTX_5060.csv` -- the
  REDUCE operation before and after the Montgomery-based rewrite.
