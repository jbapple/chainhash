# ChainHash-256 measurements

ChainHash-256 (`chainhash256.h`) against ChainHash-128 (commit 8806d02 of
`include/chainhash128.h`) and XXH3-128, one-shot hashes of 64 B to 16 MiB,
all three in the same binary.
Development name: PH-256 v2 with the v1.2 schedule (field point r = n = 4
over GF(2^256), two-level outer stage, 4 KiB blocks). The bulk figures
below were measured on v1.1 (independent masks); v2 changes only the key
(the 128 masks are the powers of one element, derived at key setup) and
v1.2 only the short-input schedule, so the region kernels timed here are
the shipped ones. The v2 known-answer vectors differ from v1's.

## Method

Three hosts: an Intel Xeon Platinum 8375C (Ice Lake-SP), an AMD EPYC 9R14
(Zen 4) and an Apple M2 Pro. Compilers: clang 21 on the Xeon, gcc 11 on
Zen 4 (the only compiler there), Apple clang with `-mcpu=apple-m2` on the
M2. On x86 ChainHash-256 is built for baseline x86-64-v2 and selects its
AVX-512 VPCLMULQDQ backend at run time; the M2 uses the NEON backend.

The harness times each length in round-robin rounds across all hashes and
takes the median: 15 rounds on x86, 13 on the M2. The sweep runs twice per
input offset, and each cell is the mean of the two medians. x86 figures are
bytes per `rdtsc` tick; multiply by the TSC frequency (2.9 GHz on the Xeon,
2.6 GHz on Zen 4) for GB/s. The M2 is timed with `mach_absolute_time` and
reported in GB/s. On the Xeon the benchmark is pinned to CPUs 16–23 under
`nice -n 10`; the machine is shared and its load average was 24–27 during
the sweep. Zen 4 ran on cores 4–7 at load about 2.5, the M2 at load
about 5.

"Offset 1" is the same length with the input one byte off a 64-byte
boundary.

## Bulk summary

| Host | 256 KiB | 1 MiB | ×ChainHash-128 at 256 KiB / 1 MiB | ×XXH3-128 at 256 KiB / 1 MiB |
| --- | ---: | ---: | --- | --- |
| Xeon 8375C (B/tick) | 13.51 | 13.19 | 0.87× / 0.89× | 0.59× / 0.62× |
| Zen 4 (B/tick) | 15.60 | 15.46 | 0.87× / 0.89× | 0.60× / 0.62× |
| M2 Pro (GB/s) | 39.04 | 38.82 | 0.91× / 0.90× | 0.94× / 0.94× |

On the Xeon, ChainHash-256 overtakes ChainHash-128 once the input leaves
the caches (1.13–1.14× at 4–16 MiB), because its prefetched 8 KiB segments
stream better from L3 and DRAM. On the M2 it peaks at 44.31 GB/s at 64 KiB
(1.05× ChainHash-128).

## Full sweeps, input offset 0

### Xeon 8375C, clang 21 (B/tick; × 2.9 = GB/s)

| Length | ChainHash-256 | ChainHash-128 | XXH3-128 | ×ChainHash-128 | ×XXH3-128 | ChainHash-256, offset 1 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 64 B | 0.41 | 0.53 | 3.27 | 0.77× | 0.12× | 0.40 |
| 128 B | 0.56 | 0.79 | 4.28 | 0.71× | 0.13× | 0.56 |
| 256 B | 1.11 | 1.41 | 5.28 | 0.79× | 0.21× | 1.10 |
| 512 B | 2.26 | 2.87 | 8.94 | 0.78× | 0.25× | 2.23 |
| 1 KiB | 3.88 | 5.17 | 12.37 | 0.75× | 0.31× | 3.85 |
| 2 KiB | 6.07 | 7.97 | 13.73 | 0.76× | 0.44× | 6.01 |
| 3 KiB | 7.48 | 9.74 | 15.78 | 0.77× | 0.47× | 7.40 |
| 4 KiB | 8.51 | 9.73 | 17.14 | 0.87× | 0.50× | 8.38 |
| 8 KiB | 9.57 | 12.02 | 19.65 | 0.80× | 0.49× | 9.43 |
| 16 KiB | 10.67 | 13.59 | 21.16 | 0.78× | 0.50× | 10.57 |
| 64 KiB | 13.18 | 14.94 | 22.41 | 0.88× | 0.59× | 11.53 |
| 256 KiB | 13.51 | 15.45 | 22.90 | 0.87× | 0.59× | 11.82 |
| 1 MiB | 13.19 | 14.85 | 21.38 | 0.89× | 0.62× | 11.83 |
| 4 MiB | 10.44 | 9.17 | 11.41 | 1.14× | 0.91× | 10.05 |
| 16 MiB | 10.29 | 9.10 | 11.19 | 1.13× | 0.92× | 9.87 |

### Zen 4 EPYC 9R14, gcc 11 (B/tick; × 2.6 = GB/s)

| Length | ChainHash-256 | ChainHash-128 | XXH3-128 | ×ChainHash-128 | ×XXH3-128 | ChainHash-256, offset 1 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 64 B | 0.55 | 0.62 | 3.95 | 0.90× | 0.14× | 0.55 |
| 128 B | 0.70 | 0.83 | 5.25 | 0.85× | 0.13× | 0.71 |
| 256 B | 1.41 | 1.48 | 7.64 | 0.95× | 0.19× | 1.42 |
| 512 B | 2.76 | 2.75 | 12.64 | 1.00× | 0.22× | 2.77 |
| 1 KiB | 4.79 | 4.88 | 15.71 | 0.98× | 0.30× | 4.79 |
| 2 KiB | 7.26 | 7.46 | 18.31 | 0.97× | 0.40× | 7.36 |
| 3 KiB | 9.11 | 9.24 | 20.26 | 0.99× | 0.45× | 8.99 |
| 4 KiB | 8.49 | 9.66 | 21.94 | 0.88× | 0.39× | 8.45 |
| 8 KiB | 9.29 | 12.72 | 23.85 | 0.73× | 0.39× | 9.42 |
| 16 KiB | 10.53 | 15.00 | 24.85 | 0.70× | 0.42× | 10.53 |
| 64 KiB | 15.20 | 17.25 | 25.65 | 0.88× | 0.59× | 13.81 |
| 256 KiB | 15.60 | 17.93 | 25.91 | 0.87× | 0.60× | 14.30 |
| 1 MiB | 15.46 | 17.44 | 24.95 | 0.89× | 0.62× | 13.97 |
| 4 MiB | 15.38 | 17.18 | 25.02 | 0.90× | 0.61× | 13.84 |
| 16 MiB | 15.31 | 17.14 | 24.98 | 0.89× | 0.61× | 13.78 |

### M2 Pro, Apple clang (GB/s)

| Length | ChainHash-256 | ChainHash-128 | XXH3-128 | ×ChainHash-128 | ×XXH3-128 | ChainHash-256, offset 1 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 64 B | 2.05 | 1.73 | 12.31 | 1.18× | 0.17× | 2.05 |
| 128 B | 3.65 | 1.37 | 15.33 | 2.67× | 0.24× | 3.79 |
| 256 B | 6.61 | 5.29 | 20.40 | 1.25× | 0.32× | 6.62 |
| 512 B | 11.54 | 8.56 | 27.98 | 1.35× | 0.41× | 11.60 |
| 1 KiB | 17.72 | 13.98 | 34.48 | 1.27× | 0.51× | 17.67 |
| 2 KiB | 24.12 | 20.50 | 36.90 | 1.18× | 0.65× | 24.60 |
| 3 KiB | 28.07 | 24.32 | 38.30 | 1.15× | 0.73× | 28.62 |
| 4 KiB | 30.81 | 30.32 | 38.92 | 1.02× | 0.79× | 31.44 |
| 8 KiB | 32.26 | 34.83 | 39.75 | 0.93× | 0.81× | 33.66 |
| 16 KiB | 36.39 | 39.35 | 38.70 | 0.92× | 0.94× | 35.48 |
| 64 KiB | 44.31 | 42.32 | 40.92 | 1.05× | 1.08× | 43.81 |
| 256 KiB | 39.04 | 43.13 | 41.53 | 0.91× | 0.94× | 38.69 |
| 1 MiB | 38.82 | 43.35 | 41.50 | 0.90× | 0.94× | 38.61 |
| 4 MiB | 36.96 | 41.40 | 39.88 | 0.89× | 0.93× | 38.33 |
| 16 MiB | 36.44 | 41.80 | 40.70 | 0.87× | 0.90× | 35.88 |

## Region kernel against its bound

Per 256 B, 4 KiB blocks, timed with calibrated probes of the CLMUL and
ternary-logic throughput.

| Host / compiler | Bound | Measured | Measured / bound |
| --- | ---: | ---: | ---: |
| Xeon, clang 21 (L2-resident) | 18.2 ticks | 18.84–18.86 ticks | 96.6% |
| Xeon, gcc 11 `-O3` (L2-resident) | 18.2 ticks | 19.08 ticks | 95.4% |
| Zen 4, gcc 11 `-O3` (L2-resident) | 15.7 ticks | 16.07 ticks | 98% (fitted bound, ±6%) |
| M2 Pro | 44.0 GB/s | 41.4 GB/s | 94% |

gcc at `-O2` reaches only 89% (Xeon) and 91% (Zen 4); build with `-O3`.

## Validation

The header's test script checks the 54 known-answer vectors from the
independent Python implementation and, on each fast backend, 61,031
one-shot inputs (every length up to one region plus two blocks plus 64 B,
misaligned, plus 20,000 random) and 4,006 streaming runs against the
reference. It passes on the Xeon (clang 21, gcc 15, gcc 11; AVX-512,
PCLMUL+SSE and portable, and both forced x86 latency paths), on Zen 4
(gcc 11) and on the M2 (Apple clang 17, clang 22; NEON and portable).
