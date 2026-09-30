# ChainHash-512 measurements

ChainHash-512 (`chainhash512.h`) against ChainHash-128
(`include/chainhash128.h`) and XXH3-128, one-shot and streaming, 64 B to
16 MiB, all in the same binary.
Development name: PH-512 v1.1 (GF(2^512) field chain, 16 KiB blocks,
two-level outer stage over regions of 8 blocks, 576-byte key). v1.1
changed only the schedule for inputs of at most one block; the bulk
figures below are v1's and the digests are identical.

## Method

Three hosts: an Intel Xeon Platinum 8375C (Ice Lake-SP), an AMD EPYC 9R14
(Zen 4) and an Apple M2 Pro. Compilers: clang 21 and gcc 11.5 on the Xeon,
gcc 11 on Zen 4, Apple clang on the M2. The benchmark translation units are
built with `-O3 -march=native` on x86. ChainHash-512 picks its backend at
run time (AVX-512 VPCLMULQDQ on both x86 hosts, NEON PMULL on the M2).

Timing is `clock_gettime(CLOCK_MONOTONIC)` on every host; there is no TSC
conversion. Each length is timed in 15 round-robin rounds across all
hashes. Cells are the median GB/s over the rounds; the ×ChainHash-128
column is the median of the per-round ratios, so it need not equal the
ratio of the two medians. On x86 the benchmark is pinned to one core under
`nice -n 10`.

## Bulk summary (GB/s, 1 MiB, aligned)

| Host | ChainHash-512 | ChainHash-128 | XXH3-128 |
| --- | ---: | ---: | ---: |
| Xeon, clang 21 | 23.96 | 41.34 | 61.78 |
| Xeon, gcc 11.5 | 22.67 | 38.37 | 51.55 |
| Zen 4, gcc 11 | 23.24 | 40.88 | 44.40 |
| M2 Pro | 26.70 | 42.12 | 40.54 |

With the input one byte off a 64-byte boundary, ChainHash-512 at 1 MiB
measured 23.85 (Xeon clang), 22.91 (Xeon gcc), 23.34 (Zen 4) and 26.14 (M2)
GB/s, within 3% of aligned.

Built without `-march` (the production configuration, runtime dispatch
only), the 1 MiB figure is 23.43 GB/s on the Xeon under clang and 23.74 on
Zen 4. The corresponding Xeon gcc run coincided with a contention spike
that lowered every column by about 25% and is not reported.

## Full sweeps, aligned input (GB/s)

### Xeon 8375C, clang 21

| Length | ChainHash-512 | ×ChainHash-128 | ChainHash-512 streaming | ChainHash-128 | XXH3-128 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 64 B | 0.58 | 0.427× | 0.55 | 1.35 | 4.99 |
| 256 B | 2.18 | 0.552× | 2.07 | 3.95 | 9.84 |
| 1 KiB | 6.36 | 0.455× | 5.52 | 13.96 | 25.59 |
| 4 KiB | 14.65 | 0.550× | 11.89 | 26.65 | 46.93 |
| 8 KiB | 18.72 | 0.566× | 14.34 | 33.07 | 54.52 |
| 16 KiB | 21.65 | 0.575× | 19.54 | 37.70 | 59.42 |
| 64 KiB | 23.04 | 0.569× | 22.42 | 40.51 | 62.53 |
| 128 KiB | 23.73 | 0.567× | 23.31 | 41.85 | 63.46 |
| 256 KiB | 23.95 | 0.564× | 23.60 | 42.46 | 63.86 |
| 1 MiB | 23.96 | 0.580× | 23.71 | 41.34 | 61.78 |
| 16 MiB | 21.75 | 0.742× | 21.74 | 29.37 | 32.94 |

### Xeon 8375C, gcc 11.5

| Length | ChainHash-512 | ×ChainHash-128 | ChainHash-512 streaming | ChainHash-128 | XXH3-128 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 64 B | 0.49 | 0.369× | 0.48 | 1.33 | 4.61 |
| 256 B | 1.88 | 0.534× | 1.76 | 3.53 | 9.80 |
| 1 KiB | 5.16 | 0.459× | 5.30 | 11.53 | 30.37 |
| 4 KiB | 13.80 | 0.567× | 11.66 | 24.44 | 48.41 |
| 8 KiB | 18.11 | 0.588× | 14.79 | 30.74 | 52.70 |
| 16 KiB | 21.27 | 0.599× | 18.36 | 35.55 | 54.23 |
| 64 KiB | 22.27 | 0.574× | 21.49 | 38.57 | 53.93 |
| 128 KiB | 22.87 | 0.584× | 22.18 | 39.12 | 54.06 |
| 256 KiB | 23.14 | 0.582× | 22.68 | 39.73 | 54.25 |
| 1 MiB | 22.67 | 0.597× | 22.67 | 38.37 | 51.55 |
| 16 MiB | 20.09 | 0.813× | 19.76 | 24.58 | 26.58 |

### Zen 4 EPYC 9R14, gcc 11

| Length | ChainHash-512 | ×ChainHash-128 | ChainHash-512 streaming | ChainHash-128 | XXH3-128 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 64 B | 0.46 | 0.353× | 0.43 | 1.29 | 5.10 |
| 256 B | 1.72 | 0.508× | 1.65 | 3.39 | 11.84 |
| 1 KiB | 5.79 | 0.516× | 5.26 | 11.23 | 31.81 |
| 4 KiB | 14.47 | 0.599× | 12.46 | 24.17 | 42.95 |
| 8 KiB | 18.90 | 0.601× | 15.66 | 31.52 | 44.69 |
| 16 KiB | 21.99 | 0.591× | 19.16 | 37.24 | 45.71 |
| 64 KiB | 23.40 | 0.549× | 22.70 | 42.63 | 45.37 |
| 128 KiB | 23.56 | 0.539× | 23.34 | 43.71 | 44.93 |
| 256 KiB | 23.54 | 0.532× | 23.43 | 44.27 | 45.14 |
| 1 MiB | 23.24 | 0.569× | 23.21 | 40.88 | 44.40 |
| 16 MiB | 23.10 | 0.563× | 23.01 | 41.04 | 44.39 |

### M2 Pro, Apple clang

| Length | ChainHash-512 | ×ChainHash-128 | ChainHash-512 streaming | ChainHash-128 | XXH3-128 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 64 B | 0.65 | 0.386× | 0.61 | 1.69 | 9.24 |
| 256 B | 2.57 | 0.480× | 2.39 | 5.33 | 20.13 |
| 1 KiB | 8.39 | 0.622× | 7.14 | 13.50 | 34.16 |
| 4 KiB | 17.83 | 0.597× | 13.99 | 29.87 | 38.87 |
| 8 KiB | 21.80 | 0.624× | 16.87 | 34.93 | 39.71 |
| 16 KiB | 24.46 | 0.641× | 22.53 | 38.13 | 39.90 |
| 64 KiB | 26.11 | 0.636× | 25.73 | 41.04 | 40.45 |
| 128 KiB | 26.51 | 0.637× | 26.35 | 41.63 | 40.45 |
| 256 KiB | 26.61 | 0.636× | 26.49 | 41.86 | 40.41 |
| 1 MiB | 26.70 | 0.634× | 26.70 | 42.12 | 40.54 |
| 16 MiB | 26.52 | 0.644× | 26.58 | 41.08 | 39.85 |

Streaming matches one-shot in bulk. At 1–8 KiB it is 15–25% lower, from
the buffer copy and a separate finalizer call.

## Validation

Built without `-march` and dispatching at run time, every backend matches
the reference implementation: 14,498 one-shot hashes and 2,718 streaming
runs on the Xeon (clang and gcc) and Zen 4 (gcc), each exercising all three
x86 backends (AVX-512, PCLMUL+SSE4.1, portable), and 20,249 and 5,062 on
the M2. The 72 known-answer vectors (0 to 384 KiB, every region shape in
regions 1–3) match on every backend, one-shot and streaming, and the
independent Python oracle reproduces them byte for byte.
