# ChainHash-128 v2 measurements

ChainHash-128 v2 (`include/chainhash128v2.h`) against ChainHash-128
(`include/chainhash128.h`, commit 2c61966) and XXH3-128, one-shot hashes
of 16 B to 1 MiB, all three in the same binary.
Development name: CH-128/P v2.2 (the PENCIL block with the two-level outer
stage, 2 KiB blocks, 8 blocks per 16 KiB region, 192-byte key). v2.2
changed only the schedule (short inputs; on Zen 4 a pre-broadcast mask
table; on NEON a region prefetch); the digests are v2's.

## Method

Three hosts: an Intel Xeon Platinum 8375C (Ice Lake-SP), an AMD EPYC 9R14
(Zen 4) and an Apple M2 Pro. On x86 each compiler builds its own binary
with `-O3 -march=native`: clang 21.1.8 and gcc 11.5 on the Xeon, clang 18
and gcc 11.5 on Zen 4. The M2 uses Apple clang with performance-core QoS.
All hashes are called through the same function-pointer harness.

x86 throughput is bytes per `rdtsc` tick times the TSC frequency (2.9 GHz
on the Xeon, 2.6 GHz on Zen 4). The TSC is not the core clock: the Xeon
core runs at about 3.36 GHz under a zmm CLMUL load. The M2 uses
`CLOCK_MONOTONIC`. Each cell is the median over invocations (5 on the Xeon,
3 on Zen 4 and the M2) of an 11-sample median. On the Xeon the benchmark
is pinned to an idle core and run under `nice`.

The sweep was run on the v2.0 header. v2.1 computes the same digests
(3,000 lengths cross-checked on the Xeon and Zen 4, 0 mismatches). It
changed only CPUID detection, gcc 9 support, the XMM product method, the
≤ 128 B path and key expansion. Its bulk speed is 1.00× v2.0 on Zen 4 and
within ±2% on the Xeon.

## Bulk summary (GB/s)

| Host | 256 KiB: v2 / ChainHash-128 / XXH3-128 | 1 MiB: v2 / ChainHash-128 / XXH3-128 | v2 / ChainHash-128 at 256 KiB | Worst v2 / ChainHash-128 (≤ 1 KiB) |
| --- | --- | --- | ---: | ---: |
| Xeon, clang | 72.4 / 44.8 / 65.4 | 58.6 / 40.8 / 57.6 | 1.62× | 0.84× (256 B) |
| Xeon, gcc | 71.5 / 42.2 / 67.0 | 68.3 / 41.4 / 64.2 | 1.69× | 0.88× (128 B) |
| Zen 4, clang | 68.6 / 29.4 / 71.5 | 64.7 / 28.9 / 68.2 | 2.33× | 0.87× (64 B) |
| Zen 4, gcc | 67.2 / 45.0 / 61.6 | 65.0 / 43.2 / 61.3 | 1.49× | 0.78× (16 B) |
| M2 Pro | 52.8 / 42.1 / 40.5 | 53.0 / 42.2 / 40.5 | 1.26× | 0.76× (128 B) |

The Zen 4 clang ratio is high because ChainHash-128's clang build of its
ZMM schoolbook loop is slow on Zen 4.

## Full sweeps, aligned input (GB/s)

### Xeon 8375C, clang 21

| Length | ChainHash-128 v2 | ChainHash-128 | XXH3-128 | ×ChainHash-128 | ×XXH3-128 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 16 B | 0.63 | 0.67 | 4.02 | 0.93× | 0.16× |
| 32 B | 0.98 | 1.04 | 6.22 | 0.94× | 0.16× |
| 64 B | 1.76 | 1.93 | 9.10 | 0.91× | 0.19× |
| 128 B | 2.90 | 3.38 | 13.17 | 0.86× | 0.22× |
| 256 B | 3.87 | 4.61 | 14.62 | 0.84× | 0.26× |
| 512 B | 8.01 | 8.48 | 23.77 | 0.94× | 0.34× |
| 1 KiB | 14.40 | 14.38 | 34.28 | 1.00× | 0.42× |
| 2 KiB | 23.60 | 21.67 | 38.79 | 1.09× | 0.61× |
| 4 KiB | 35.30 | 28.24 | 48.61 | 1.25× | 0.73× |
| 8 KiB | 45.96 | 34.83 | 56.00 | 1.32× | 0.82× |
| 16 KiB | 64.01 | 39.50 | 60.44 | 1.62× | 1.06× |
| 64 KiB | 70.46 | 43.47 | 64.17 | 1.62× | 1.10× |
| 256 KiB | 72.44 | 44.81 | 65.40 | 1.62× | 1.11× |
| 1 MiB | 58.59 | 40.80 | 57.58 | 1.44× | 1.02× |

### Xeon 8375C, gcc 11.5

| Length | ChainHash-128 v2 | ChainHash-128 | XXH3-128 | ×ChainHash-128 | ×XXH3-128 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 16 B | 0.55 | 0.58 | 4.49 | 0.94× | 0.12× |
| 32 B | 0.87 | 0.89 | 6.58 | 0.97× | 0.13× |
| 64 B | 1.53 | 1.66 | 9.56 | 0.92× | 0.16× |
| 128 B | 2.58 | 2.93 | 11.37 | 0.88× | 0.23× |
| 256 B | 3.54 | 3.69 | 14.27 | 0.96× | 0.25× |
| 512 B | 6.74 | 6.74 | 23.68 | 1.00× | 0.28× |
| 1 KiB | 12.35 | 11.68 | 35.60 | 1.06× | 0.35× |
| 2 KiB | 20.65 | 18.72 | 48.27 | 1.10× | 0.43× |
| 4 KiB | 30.67 | 25.34 | 56.36 | 1.21× | 0.54× |
| 8 KiB | 40.12 | 31.77 | 61.38 | 1.26× | 0.65× |
| 16 KiB | 60.00 | 36.43 | 64.23 | 1.65× | 0.93× |
| 64 KiB | 68.56 | 40.79 | 66.11 | 1.68× | 1.04× |
| 256 KiB | 71.48 | 42.20 | 66.96 | 1.69× | 1.07× |
| 1 MiB | 68.29 | 41.42 | 64.20 | 1.65× | 1.06× |

### Zen 4 EPYC 9R14, clang 18

| Length | ChainHash-128 v2 | ChainHash-128 | XXH3-128 | ×ChainHash-128 | ×XXH3-128 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 16 B | 0.68 | 0.77 | 5.38 | 0.88× | 0.13× |
| 32 B | 1.12 | 1.20 | 8.25 | 0.93× | 0.14× |
| 64 B | 1.89 | 2.17 | 11.93 | 0.87× | 0.16× |
| 128 B | 3.59 | 3.52 | 15.07 | 1.02× | 0.24× |
| 256 B | 5.45 | 5.97 | 21.79 | 0.91× | 0.25× |
| 512 B | 9.59 | 10.13 | 39.08 | 0.95× | 0.25× |
| 1 KiB | 16.28 | 16.86 | 54.21 | 0.97× | 0.30× |
| 2 KiB | 24.05 | 25.44 | 55.54 | 0.95× | 0.43× |
| 4 KiB | 34.78 | 19.73 | 62.69 | 1.76× | 0.55× |
| 8 KiB | 42.29 | 23.73 | 66.73 | 1.78× | 0.63× |
| 16 KiB | 62.59 | 26.55 | 68.99 | 2.36× | 0.91× |
| 64 KiB | 67.17 | 28.69 | 70.75 | 2.34× | 0.95× |
| 256 KiB | 68.55 | 29.38 | 71.46 | 2.33× | 0.96× |
| 1 MiB | 64.65 | 28.94 | 68.20 | 2.23× | 0.95× |

### Zen 4 EPYC 9R14, gcc 11.5

| Length | ChainHash-128 v2 | ChainHash-128 | XXH3-128 | ×ChainHash-128 | ×XXH3-128 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 16 B | 0.48 | 0.62 | 4.79 | 0.78× | 0.10× |
| 32 B | 0.80 | 0.99 | 6.83 | 0.81× | 0.12× |
| 64 B | 1.45 | 1.82 | 10.23 | 0.80× | 0.14× |
| 128 B | 3.14 | 2.96 | 12.79 | 1.06× | 0.25× |
| 256 B | 4.27 | 4.59 | 17.80 | 0.93× | 0.24× |
| 512 B | 6.54 | 8.13 | 29.07 | 0.80× | 0.23× |
| 1 KiB | 13.42 | 13.18 | 40.63 | 1.02× | 0.33× |
| 2 KiB | 20.75 | 20.40 | 54.18 | 1.02× | 0.38× |
| 4 KiB | 29.72 | 24.34 | 58.51 | 1.22× | 0.51× |
| 8 KiB | 38.78 | 31.99 | 60.65 | 1.21× | 0.64× |
| 16 KiB | 58.26 | 37.71 | 62.05 | 1.55× | 0.94× |
| 64 KiB | 65.19 | 43.14 | 61.52 | 1.51× | 1.06× |
| 256 KiB | 67.20 | 45.04 | 61.56 | 1.49× | 1.09× |
| 1 MiB | 65.02 | 43.21 | 61.26 | 1.50× | 1.06× |

### M2 Pro, Apple clang

| Length | ChainHash-128 v2 | ChainHash-128 | XXH3-128 | ×ChainHash-128 | ×XXH3-128 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 16 B | 1.07 | 0.95 | 5.41 | 1.13× | 0.20× |
| 32 B | 1.67 | 1.66 | 9.26 | 1.00× | 0.18× |
| 64 B | 2.71 | 3.13 | 12.71 | 0.87× | 0.21× |
| 128 B | 4.17 | 5.50 | 15.87 | 0.76× | 0.26× |
| 256 B | 5.97 | 5.46 | 20.87 | 1.09× | 0.29× |
| 512 B | 10.71 | 8.66 | 28.07 | 1.24× | 0.38× |
| 1 KiB | 14.12 | 13.98 | 34.43 | 1.01× | 0.41× |
| 2 KiB | 21.54 | 21.07 | 37.77 | 1.02× | 0.57× |
| 4 KiB | 32.45 | 30.76 | 39.91 | 1.05× | 0.81× |
| 8 KiB | 43.73 | 35.83 | 40.74 | 1.22× | 1.07× |
| 16 KiB | 54.22 | 39.23 | 40.85 | 1.38× | 1.33× |
| 64 KiB | 57.83 | 41.54 | 40.40 | 1.39× | 1.43× |
| 256 KiB | 52.78 | 42.05 | 40.54 | 1.26× | 1.30× |
| 1 MiB | 52.96 | 42.17 | 40.54 | 1.26× | 1.31× |

## Misaligned input

With the input one byte off a 64-byte boundary, v2 measured 67.28 (Xeon
clang), 67.56 (Xeon gcc), 68.26 (Zen 4 clang) and 67.22 (Zen 4 gcc) GB/s at
256 KiB. The largest x86 drop is at 1 MiB on the Xeon under gcc, 60.86
against 68.29 aligned. On the M2 misalignment costs 14% at 256 KiB–1 MiB
(45.20 and 45.39 GB/s); ChainHash-128 drops 8% there (38.74 and 38.79).

## Key setup (v2.1)

Key expansion uses in-register carry-less products in v2.1.

| Host | v2.0 | v2.1 |
| --- | ---: | ---: |
| Xeon, clang (ticks) | 30,218 | 7,248 |
| Xeon, gcc (ticks) | 42,454 | 9,160 |
| Zen 4, gcc (ticks) | 31,122 | 10,816 |
| M2 Pro (µs) | 9.5 | 1.8 |

## Validation

The v2.1 header matches the independent Python oracle's 151 known-answer
vectors (lengths 0 to 49,153 B) on every backend, one-shot and streaming,
plus 20,000 random inputs per backend: 24 builds on the Xeon (clang 21,
gcc 11, gcc 9.5; dispatched, each XMM product method forced, portable
only), the Zen 4 builds under gcc 11 and clang 18, and the M2 (NEON and
portable). Zero mismatches everywhere.
