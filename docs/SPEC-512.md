# ChainHash-512 specification

This defines the function in [`include/chainhash512.h`](../include/chainhash512.h) (with
[`include/chainhash512_body.inc`](../include/chainhash512_body.inc)), a keyed 512-bit hash: a
pseudo-dot-product over GF(2^512) with a power key, 16 KiB blocks and the version 2 family's
two-level outer stage ([FAMILY.md](FAMILY.md)). The collision bound is in
[THEOREM-512.md](THEOREM-512.md). Development name: PH-512 v1.1 (the v1 function; v1.1 changed only the
schedule for inputs of at most one block, so the digests are v1's).

The independent oracle [`test/512/pyref512.py`](../test/512/pyref512.py) is written from this
document with Python integers and shares no code with the header; it regenerates
[`test/512/vectors.txt`](../test/512/vectors.txt) byte for byte.

## 1. Definition

**Field.** L = GF(2)[x]/(f), f = x^512 + x^8 + x^5 + x^2 + 1. f is irreducible (Rabin,
[`test/512/cert512.py`](../test/512/cert512.py)). An element is a 512-bit string: bit i is the
coefficient of x^i, stored as 8 little-endian 64-bit limbs. Addition is XOR.

**Key.** 576 bytes, 9 independent uniform 512-bit words, each little endian, in this order:

| bytes | word |
|---|---|
| 0–63 | s |
| 64–127 | y |
| 128–191 | τ |
| 192–511 | c0, c1, c2, c3, c4 |
| 512–575 | z |

No key is rejected. `chainhash512_key_from_words` takes the same bytes as 72 little-endian 64-bit
words. `chainhash512_key_from_seed` fills the 576 bytes from SplitMix64; it is for benchmarks and
tests, and the bound is not asserted for it.

**Pairs.** A message of n bytes is cut into 1 KiB chunks. Pairs are indexed g = 0, 1, … in message
order.

- **Full chunk.** It holds pairs 8c..8c+7. Pair i of the chunk has
  - u = the 512-bit word whose limb j is the 8 bytes at offset 64j + 8i, and
  - v = the word whose limb j is the 8 bytes at 512 + 64j + 8i
  (limb-major layout).
- **Final partial chunk** (n mod 1024 ≠ 0). It holds ⌈(n mod 1024)/128⌉ pairs back to back, 128
  bytes each: u is the first 64 bytes and v the next 64. The last pair is zero-padded.
- **Count.** G = 8⌊n/1024⌋ + ⌈(n mod 1024)/128⌉ pairs.

**Blocks.** Block t (1-based) holds pairs 128(t−1)..128t−1, at positions q = 0..127. Let
m = max(1, ⌈G/128⌉). Every block except the last is exactly 16 KiB of message. The block value is

    b_t = Σ_q (u_q + s^(2q+1)) · (v_q + s^(2q+2))        (product in L)

The empty message has m = 1 and b_1 = 0.

**Regions.** Regions hold R = 8 block values each, so there are m′ = ⌈m/8⌉ regions. Region ρ
holds a = (b_{8ρ−7}, …) with q = 8 values, except that the last region holds
q = m − 8(m′−1) ∈ {1..8}. With h = ⌈q/2⌉ and f = ⌊q/2⌋:

    c(a) = Σ_{i=1}^{f} (a_i + y^(2i−1)) · (a_{h+i} + y^(2i))  +  [q odd] · a_h

**Outer, then output.** n is embedded in L as the integer n < 2^64 in limb 0.

    V = n · z^(m′) + Σ_{ρ=1}^{m′} c(a^(ρ)) · z^(m′−ρ)
    v = (V + τ) mod 2^512                          (integer addition of the 512-bit strings)
    H = (v + c2) · ((v² + c0)(v + v² + c1) + c3) + c4

The output is H as 64 little-endian bytes.

- **≤ 1 block** (n ≤ 16384): V = n·z + b_1.
- **Streaming** computes V as W + n·z^(m′), where W is the Horner sum from 0 (section 2).

## 2. Interface and implementation notes

These are not part of the definition.

**Backends.** All give identical digests.

| backend | ISA | selection |
|---|---|---|
| CH512_AVX512 | AVX-512 F, DQ, BW, VL + VPCLMULQDQ | CPUID leaf 7, XGETBV XCR0 & 0xE6 = 0xE6 |
| CH512_PCLMUL | PCLMULQDQ + SSE4.1 | CPUID leaf 1 |
| CH512_NEON | AArch64 + PMULL | compile time: `__ARM_FEATURE_AES` (Apple's default target has it) |
| CH512_PORTABLE | bit-serial | fallback |

The backend is chosen when the key is initialized (`*_with_backend` variants pick one; an
unavailable one falls back to the best). Every x86 kernel carries its own target attribute, so the
header builds and dispatches correctly without `-march`. `CHAINHASH512_PORTABLE` omits all
hardware code. The expanded key is about 110 KiB of tables (initialize it in place, static or on
the heap); the stream state buffers one 16 KiB block.

**Blocks.** Karatsuba³ with 27 point sums accumulated unreduced across the block.
- x86: a 3-pass ZMM kernel (9 accumulators per pass, no in-loop spills).
- NEON: 3 passes with the Apple-fused `pmull; eor`.
- Tail pairs: register-resident vector evaluation.

**Region step.** Recombine+reduce R is F2-linear and R(E(a)⊙E(b)) = a·b mod f. So

    V·z + c  =  R( E(V)⊙E(z) + Σ_i E(a_i + y^(2i−1))⊙E(a_{h+i} + y^(2i)) )  [+ a_h]

- **Cost.** One reduction per region: 5 point-products of 27 clmul each, plus 8 block reductions.
- **Odd q.** The bare block a_h is not reduced separately; its unreduced point sums seed the region
  accumulator.
- **≤ 1 block.** The path is V = R(P_1 + E(n)⊙E(z)), with one reduction.

**Streaming** (`chainhash512_init/update/final`). It buffers one 16 KiB block and keeps at most 8
block values.
- A full block is laid out identically whether or not it is the last one, and a full region has
  the same c. So both are folded as soon as they complete, into W by Horner in z from 0.
- `final`:
  - If nothing has been folded yet (n ≤ 16 KiB), it applies the one-shot function to the buffer.
  - Otherwise it hashes the partial block, folds the last region with its actual q, and adds
    n·z^(m′), computing z^(m′) by square-and-multiply.

**Finalizer.** v² is computed as 8 self-clmuls plus a reduction (characteristic 2).

## 3. Bound

For messages M ≠ M′ of at most 8L bytes, fixed before the key, Pr[H(M) = H(M′)] ≤ N(L)/2^512 with
N(L) ≤ 2L, so the score is 511 bits. The statement, the numerators and the proof status are in
[THEOREM-512.md](THEOREM-512.md).

## 4. Test vectors

- Key: key[i] = (73i + 11) & 255 for 576 bytes. Message: msg[i] = (137i + 29) & 255.
- [`test/512/vectors.txt`](../test/512/vectors.txt): 72 lengths from 0 to 393,216 B. They cover
  the ≤ 1-block path, chunk and block edges, and every q = 1..8 in regions 1, 2 and 3.
- The C reference, every backend (one-shot and streaming), the key given as words, and a
  `CHAINHASH512_PORTABLE` build all reproduce them (`make test-512`), and `pyref512.py`
  regenerates the file.
