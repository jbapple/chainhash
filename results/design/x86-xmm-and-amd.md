# x86 XMM paths, short inputs and AMD dispatch (2026-09-28)

The machines behind the public SMHasher tables do not take the AVX-512
path. rurban/smhasher's main table comes from an AMD Ryzen 5 3350G (Zen+),
its alternates from an Intel i7-6820HQ (Skylake), an AMD EPYC 9554P
(Zen 4) and an Intel i5-2300 (Sandy Bridge); SMHasher3's from a Ryzen 9
3950X (Zen 2) built with GCC 9.3. Zen+, Zen 2, Skylake and Sandy Bridge
have PCLMULQDQ but no VPCLMULQDQ, so `chainhash()` and `chainhash128()`
run the XMM backend there, and GCC compiles it. These changes make that
path fast and make the headers build with GCC 9. Every digest is
unchanged: the frozen vectors, the property tests and an exhaustive
comparison of every backend with the portable evaluator (all lengths
0..9000, inputs ending at a guard page, GCC 9, GCC 11 and Clang, on both
hosts below) pass.

## What changed

- **GCC 9 build.** Both headers called `__get_cpuid` from `<cpuid.h>`,
  which in GCC 9 has no include guard, so a translation unit with both
  headers (both SMHasher registrations) failed to compile. CPUID is now
  inline assembly. `chainhash_calibrate.h` includes `<cpuid.h>` itself.
- **GCC 9 -O3 miscompile.** GCC 9.5 at -O3 (not -O2, not with
  `-fno-tree-vrp`, not GCC 11) compiled ChainHash-128's partial-region tail
  so that `k->yp[count-1-j]` was read with the index wrapped to 2^32-1:
  a segfault for one full region plus a short tail (4145 bytes), which
  SMHasher3's thread-safety sanity test hits. `count` is now a `size_t` in
  the x86 and NEON tails; every length 0..9000 then matches the portable
  evaluator with GCC 9 -O3.
- **ChainHash partial regions on XMM and YMM.** Only ZMM and NEON had a
  vector tail; XMM and YMM built every word of a partial region through
  the portable byte loader with one hardware product per pair. The x86
  tail now follows the NEON one (raw lane products, Horner weights
  `(y^e, X^64 y^e)` on raw lanes, one reduction) and serves every x86
  backend; ZMM keeps its own tail past 64 bytes where PCLMULQDQ is no
  faster than VPCLMULQDQ on ZMM (AMD) and past 768 bytes where it issues
  every cycle (Intel). Inputs up to
  64 bytes take a register-only path; loads never leave the input (whole
  16-byte words, the last partial word as the 16 bytes ending at the input
  end shifted down, inputs under 16 bytes from overlapping 8/4/1-byte
  loads). In the finalizer the reductions are shifts instead of two
  dependent PCLMULQDQs, and the last reduction is folded into the final
  product: with `z = x^c2`, `R` the raw middle product and
  `w = z*X^64 mod p`, `z*(R^c3) = z*(R.lo^c3) ^ w*R.hi mod p`, and `w`
  depends on `x` only.
- **ChainHash-128 XMM bulk kernel.** The kernel was one ~2,000-instruction
  loop body per 4 KiB region with sixteen raw lane registers live; GCC 9
  and 11 ran it at half of Clang's speed on the Xeon. It is now
  chunk-outer over lane pairs: a chunk's two keys serve both lanes, each
  data word is keyed by one register XOR, at most 15 vector registers are
  live, and the lane states stay in memory, touched once per region. The
  Karatsuba middle operands are `unpacklo(a,b)^unpackhi(a,b)` with one
  `0x10` product. The loop body is 86-96 instructions.
- **Product method by PCLMULQDQ rate.** AMD cores issue PCLMULQDQ once per
  two cycles at every vector width, Intel Haswell once per two and Sandy
  and Ivy Bridge once per eight, and shuffles run on other pipes, so for
  XMM Karatsuba (three products) beats schoolbook (four) there.
  `chainhash128()` keeps schoolbook for XMM only where PCLMULQDQ issues
  every cycle and shares its port with the shuffles: Intel cores with ADX
  (Broadwell and later). The same predicate sets where ZMM switches to its
  own tail.
- **ChainHash-128 short inputs.** The integer twist `+tau` is done in a
  register (lane sums plus a carry from an unsigned comparison) instead of
  a store/scalar add/load round trip; inputs under 16 bytes are loaded
  without a stack word; two words use `V = n*y^2 + A_0*(kb*y) + A_1*kb`,
  where `kb*y` depends on the key only.
- **Key expansion** multiplies with the dispatched hardware product
  instead of the bit-serial one: 8,736 → 572 TSC ticks for ChainHash and
  21,762 → 1,066 for ChainHash-128 on the Zen 4 host.

## Instruction costs

From uops.info (measured reciprocal throughput, latency): PCLMULQDQ xmm
is 2 cycles, 4 on Zen+, Zen 2, Zen 3 and Zen 4 (microcoded, 4 uops, on
Zen+ and Zen 2); VPCLMULQDQ ymm and zmm are also 2 and 4 on Zen 3/Zen 4.
On Ice Lake PCLMULQDQ xmm is 1 and 6, VPCLMULQDQ ymm and zmm 2 and 8
(three uops, two on port 5); on Skylake xmm is 1 and 7 (port 5, which also
executes the shuffles), on Broadwell 1 and 5, on Haswell 2 and 7, on Sandy
and Ivy Bridge 8 and 13-14 (18 uops). So:

- On AMD one 512-bit product does the work of four 128-bit ones at the
  same issue rate. ChainHash's XMM kernel (72 products per KiB) is bound
  at 144 cycles per KiB, 7.1 B/cycle, which it reaches on Zen 4 under every
  compiler; XMM/ZMM is at most 0.25 by product count.
- On Ice Lake YMM issues no faster than XMM per lane, which is why the
  YMM and XMM rows coincide there, and XMM/ZMM is bounded by 0.5 by
  product count.
- ChainHash-128's XMM kernel is bound at 4.74 B/cycle with Karatsuba and
  3.56 with schoolbook on AMD.

llvm-mca's `znver1`/`znver2` models carry no PCLMULQDQ entry (latency
100, throughput 0.25), so they cannot model these loops. Its
`znver3`/`znver4` models (throughput 2, latency 4) predict the ChainHash
XMM loop at 151-153 cycles per KiB (measured 144 on Zen 4) and the new
ChainHash-128 XMM loop at 2.45 (Karatsuba) and 2.2 (schoolbook) cycles
per PCLMULQDQ (measured 2.2 and 2.0). The model rated the old ChainHash-128
XMM loop the same for GCC and Clang, so it does not explain the GCC loss
on the Xeon, which disappeared with the smaller loop.

## Measurement protocol

A standalone harness in the style of the SMHasher speed tests, one binary
per compiler with every backend forced explicitly: bulk is 256 KiB at
alignments 0..7, minimum over 25 batches of 20 calls, averaged over the
alignments; short inputs are 1..31 bytes hashed 2,000 times in a chain in
which each key depends on the previous hash (as SMHasher's own loop),
minimum over 25 batches, averaged over the lengths. `-O3 -march=native
-DNDEBUG`, one pinned core. Zen 4: AMD EPYC 9R14, GCC 9.5 (container),
GCC 11.5, Clang 18; core cycles from `perf_event_open` user-mode cycles.
Xeon: Platinum 8375C (Ice Lake), GCC 9.5 (container), GCC 11.5, GCC 14,
Clang 21; invariant-TSC ticks (no PMU on that host). "Before" is
0a03c63 with only the CPUID change (0a03c63 itself does not compile with
GCC 9); the ChainHash-128 XMM rows use the dispatched product method of
each version on each host. Xeon bulk rows are the best of three
interleaved runs (the host is shared).

### Zen 4 (EPYC 9R14), core cycles

| | GCC 9 before | GCC 9 after | GCC 11 before | GCC 11 after | Clang before | Clang after |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| ChainHash XMM, 1-31 B, cycles | 272.2 | 94.1 | 290.7 | 93.7 | 166.8 | 82.0 |
| ChainHash ZMM, 1-31 B, cycles | 122.5 | 94.2 | 121.2 | 93.9 | 117.0 | 82.0 |
| ChainHash XMM, 256 KiB, B/cycle | 7.10 | 7.10 | 7.08 | 7.08 | 7.10 | 7.10 |
| ChainHash ZMM, 256 KiB, B/cycle | 24.43 | 24.51 | 24.12 | 24.13 | 19.5-21.8 | 19.5-21.8 |
| ChainHash-128 XMM, 1-31 B, cycles | 192.7 | 158.9 | 176.3 | 146.2 | 147.6 | 123.9 |
| ChainHash-128 XMM, 256 KiB, B/cycle | 3.15 | 4.33 | 3.15 | 4.35 | 3.35 | 4.72 |
| ChainHash-128 ZMM, 256 KiB, B/cycle | 13.37 | 13.37 | 13.45 | 13.41 | 13.34 | 13.38 |

YMM equals XMM on short inputs and runs at 14.16 (ChainHash) and
8.5-9.2 (ChainHash-128) B/cycle in bulk, unchanged. Clang's ChainHash ZMM
bulk alternates between two levels from run to run in both versions.

### Xeon 8375C (Ice Lake), TSC ticks

| | GCC 9 before | GCC 9 after | GCC 11 before | GCC 11 after | Clang before | Clang after |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| ChainHash XMM, 1-31 B, ticks | 323.9 | 84.6 | 296.7 | 84.2 | 168.0 | 82.0 |
| ChainHash ZMM, 1-31 B, ticks | 124.9 | 83.6 | 125.9 | 83.7 | 122.0 | 82.1 |
| ChainHash XMM, 256 KiB, B/tick | 16.20 | 16.19 | 16.29 | 16.29 | 17.03 | 17.03 |
| ChainHash ZMM, 256 KiB, B/tick | 28.25 | 28.13 | 28.41 | 28.30 | 28.29 | 28.17 |
| ChainHash-128 XMM, 1-31 B, ticks | 178.6 | 131.9 | 131.1 | 110.6 | 125.9 | 108.4 |
| ChainHash-128 XMM, 256 KiB, B/tick | 4.11 | 8.16 | 4.20 | 8.24 | 8.35 | 8.43 |
| ChainHash-128 ZMM, 256 KiB, B/tick | 14.81 | 14.57 | 14.92 | 14.88 | 14.98 | 14.99 |

GCC 14: ChainHash 1-31 B 269.2 → 76.1 ticks, ChainHash-128 XMM bulk
7.99 → 8.51 B/tick, ChainHash-128 1-31 B 187.2 → 104.3 ticks. With the
Ice Lake port model of the ceiling study, the ChainHash-128 XMM loop was
bound by the legacy decoder under GCC 9 and 11 (3,400-3,500 instructions,
22 KB per iteration); it is now 82-99 instructions and port-bound, and
GCC 11 runs it at 97% (schoolbook) and 93% (Karatsuba) of that bound.

### Crossover of the ChainHash tails (GCC 11, ticks per chained call)

| bytes | 1 | 16 | 31 | 64 | 128 | 256 | 512 | 640 | 896 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Zen 4 XMM tail | 59 | 63 | 69 | 73 | 79 | 79 | 99 | 112 | 139 |
| Zen 4 ZMM tail | 82 | 82 | 82 | 78 | 77 | 81 | 85 | 88 | 91 |
| Xeon XMM tail | 78 | 85 | 91 | 93 | 98 | 100 | 107 | 116 | 139 |
| Xeon ZMM tail | 122 | 122 | 122 | 119 | 118 | 120 | 119 | 120 | 123 |

## In the SMHasher harnesses

rurban/smhasher 2688c165 and SMHasher3 3de870c7 with the draft
registrations, CMake Release (`-O3 -march=native`), GCC 11.5 and GCC 9.5,
ticks as the tools print them. SMHasher's small keys (1-31 bytes) and its
256 KiB bulk blocks (no partial region) run the code measured above.

| SMHasher3, 1-31 B ticks / 256 KiB bytes per tick | Zen 4 GCC 11 | Zen 4 GCC 9 | Xeon GCC 11 | Xeon GCC 9 |
| --- | ---: | ---: | ---: | ---: |
| chainhash, before | 110.87 / 35.08 | 110.11 / 35.26 | 154.84 / 28.31 | 156.97 / 28.16 |
| chainhash, after | 62.46 / 35.13 | 61.83 / 35.27 | 80.08 / 28.23 | 81.46 / 28.16 |
| chainhash forced XMM, before | 157.07 / 10.08 | 131.20 / 10.10 | 176.87 / 16.28 | 178.16 / 16.18 |
| chainhash forced XMM, after | 62.48 / 10.08 | 61.81 / 10.10 | 81.76 / 16.28 | 83.59 / 16.18 |
| chainhash-128, before | 134.80 / 17.47 | 140.03 / 17.58 | 146.21 / 14.07 | 176.34 / 14.31 |
| chainhash-128, after | 117.17 / 17.39 | 130.48 / 17.35 | 129.37 / 14.23 | 143.42 / 14.27 |
| chainhash-128 forced XMM, after | 127.79 / 5.03 | 110.10 / 5.03 | 132.92 / 8.20 | 134.40 / 8.18 |
| XXH3-64, same binary | 20.38 / 18.77 | 20.48 / 18.64 | 29.48 / 19.67 | 30.15 / 20.27 |

"Before" for GCC 9 is 30c0111 with the CPUID change: without it the
registrations do not compile with GCC 9, and at 30c0111 chainhash-128
then crashed in SMHasher3's thread-safety sanity test under GCC 9 -O3
(the lane-count miscompile above); after, Sanity passes on both hosts for
both byte orders. The forced-XMM ChainHash-128 rows are schoolbook; on
AMD the dispatch now takes Karatsuba there. SMHasher3 prints GiB/s at an
assumed 3.5 GHz; its reference table was measured on a Ryzen 9 3950X
(Zen 2) pinned at 3.5 GHz, where ticks are core cycles and the XMM rows
apply.

rurban's harness passes a new 32-bit seed to every small-key call and to
every bulk trial. The draft registration expanded the key on each seed
change, so its figures were key expansion: 8,787 ticks per 1-31-byte
hash and 16.0 bytes per tick in bulk for chainhash on Zen 4 (22,820 and
7.05 for chainhash-128), against 58.6 and 33.9 with a fixed key.

## Run-time dispatch without the hardware

qemu-user 7.2 with its CPU models exercises the CPUID paths (it cannot
time them). A binary built with `-march=znver2` (GCC 9 and 11) and one
with `-march=x86-64-v2 -mavx -mpclmul` select XMM for both functions
under the EPYC (Zen 1), EPYC-Rome (Zen 2), Skylake-Client, Haswell and
SandyBridge models; Karatsuba on Zen, Sandy Bridge and Haswell, schoolbook
on Broadwell and Skylake. Both self-tests pass and the dispatched digests
equal the portable ones at every length 0..1500 under each model. qemu's
TCG does not implement VPCLMULQDQ, so the Zen 3 (YMM) selection was not
exercised there.

## Limits

The Zen 4 host stands in for Zen 2 and Zen+, which were not available.
Their PCLMULQDQ has the same measured issue rate and latency, but it is
microcoded there (4 uops from the microcode sequencer) where Zen 4 runs it
on the FP pipes; Zen 2 has a 4K-op cache and 256-bit AVX2, Zen+ 128-bit
AVX2 datapaths. The short-input kernels are latency-bound chains of
products and shifts and the XMM bulk kernels are product-bound, so the
same shape is expected, but the absolute Zen 2 and Zen+ figures, and
whether microcode sequencing lowers the product issue rate inside these
loops, are unmeasured.
