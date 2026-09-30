# ChainHash-128 v2 specification

This defines the function in [`include/chainhash128v2.h`](../include/chainhash128v2.h), a keyed
128-bit hash. It is a different function from ChainHash-128 ([SPEC-128.md](SPEC-128.md)), whose
header and digests are unchanged. The shared construction of the version 2 family is in
[FAMILY.md](FAMILY.md), and the collision bound is in [THEOREM-128v2.md](THEOREM-128v2.md).

The block is PENCIL's nine-product block, with the masks taken as powers of one key element in an
8-word algebra. The outer stage is the **two-level formula**: the block formula is applied again to
the block values, and regions of 8 block values are chained by a Horner polynomial in an
independent key word. The twist and the degree-5 finalizer are ChainHash-128's.

Bytes are little endian, and bit i of a word is X^i. The test vectors
[`test/128v2/vectors.txt`](../test/128v2/vectors.txt) come from the independent Python oracle
[`test/128v2/pyref_2l.py`](../test/128v2/pyref_2l.py). Development name: CH-128/P v2.2 (v2.1 with a
short-input schedule, a pre-broadcast mask table on Zen 4 and a NEON region prefetch; the digests
are v2's).

## 1. Fields, ring and key

- F = GF(2)[X]/(X^64+X^4+X^3+X+1) holds block coefficients, 64 bits each.
- K = GF(2)[X]/(X^128+X^7+X^2+X+1) is the outer field (ChainHash-128's).
- The roles are α_i = i (the F-element with bit pattern i), for i = 0..7, and ξ = X.
- A = F[T]/g(T), with g(T) = T^8 + 0x17d T^4 + 0x60c T^2 + 0x770 T + 0xee0. Equivalently
  g = ∏_i (T + α_i) + ξ·c with c = 0x770. E: A → F^8, E(P) = (P(α_0), …, P(α_7)), is the
  evaluation map.

**Key: 192 bytes, all words independent and uniform; no key is rejected.**

| bytes | word | role |
|---|---|---|
| 0..63 | s = Σ_{k<8} s_k T^k ∈ A (s_k = 64-bit word k) | block key (level 1) |
| 64..79 | y ∈ K | level-2 key (region formula) |
| 80..159 | c0, c1, c2, c3, c4 ∈ K | finalizer |
| 160..175 | τ (128-bit integer) | twist |
| 176..191 | z ∈ K | region chain key |

`chainhash128v2_key_from_words` takes the same 192 bytes as twelve 128-bit words (word i is bytes
16i..16i+15). `chainhash128v2_key_from_seed` fills the 192 bytes from SplitMix64; it is for
benchmarks and tests, and the bound is not asserted for it.

## 2. Blocks

**Masks.** For m = 0..15 and i = 0..7, with powers taken in A: μ_{m,i} = E(s^(2m+1))_i and
ν_{m,i} = E(s^(2m+2))_i.

**Layout.** A block is B = 2048 bytes: 16 pair-vectors of 8 roles. A layout region is 8
interleaved blocks, 16384 bytes. Word a (0..127) of block j (0..7) in layout region Q is the 16
bytes at 16384·Q + 128·a + 16·j. Write a = 8m + i (pair-vector m, role i). The word's bytes 0..7
are x_{m,i} and bytes 8..15 are z_{m,i}. Then u_{m,i} = x_{m,i} ⊕ μ_{m,i} and
v_{m,i} = z_{m,i} ⊕ ν_{m,i}.

**Block count.** Write the message length as ℓ = 16384·Q + r, 0 ≤ r < 16384. Then
m = p(ℓ) = 1 if ℓ = 0, else 8Q + (r > 0 ? min(8, ⌈r/16⌉) : 0). Blocks are numbered
t = 8Q' + j + 1.

**Presence rule.** Every pair-vector of a full layout region is present. In the final partial
layout region, pair-vector (m, j) is present iff 1024·m + 16·j < r, i.e. iff its first word
(a = 8m) has a byte. Missing bytes are zero. A present pair-vector always contributes all 8
roles, masks included. Presence depends only on ℓ.

**Block value.** For block t, sum over its present pair-vectors, in F:
Q0 = Σ_i u_i v_i and Q1 = Σ_i α_i u_i v_i + ξ (Σ_i u_i)(Σ_i v_i).
The packed value is b_t = Q0 + X^64·Q1 ∈ K, with Q0 as the low 64 bits. The empty message has one
block, b_1 = 0.

**Equivalent ring form.** Let U_m = E^{-1}(x_m) + s^(2m+1) and V_m = E^{-1}(z_m) + s^(2m+2). Then
(Q0, Q1) = c·(W_7, W_6), where W = Σ_{present m} U_m V_m mod g. This identity is F-bilinear and
is checked on a basis in [`test/128v2/cert/cert64.py`](../test/128v2/cert/cert64.py). The Python
oracle computes blocks this way.

## 3. Two-level outer stage

**Regions.** R = 8, so a level-2 region is exactly one 16 KiB layout region.
- m' = ⌈m/8⌉.
- Region ρ (1 ≤ ρ ≤ m') holds a^(ρ) = (b_{8ρ−7}, …, b_{8ρ}).
- Every region has q = 8 values except possibly the last, which has q = m − 8(m' − 1) ∈ {1..8}.

**Region value.** For a_1..a_q ∈ K, with h = ⌈q/2⌉ and f = ⌊q/2⌋:

    c(a_1..a_q) = Σ_{i=1}^{f} (a_i + y^{2i−1}) · (a_{h+i} + y^{2i})   +   [q odd] · a_h

This is the block formula read over K with key y. The first half of the region is paired with the
second half (lanes 0–3 with lanes 4–7 on ZMM). An odd q happens only in the last region, and its
middle value enters bare.

**Chain, twist, finalizer.** ℓ < 2^64 is embedded in K as an integer in the low limb.

    V = ℓ·z^{m'} + Σ_{ρ=1}^{m'} c(a^(ρ))·z^{m'−ρ}        (length-leading Horner in z)
    v = (V + τ) mod 2^128;   qv = v·v;   rv = (qv ⊕ c0)·(v ⊕ qv ⊕ c1);   H = (v ⊕ c2)·(rv ⊕ c3) ⊕ c4

The output is the 16 little-endian bytes of H (`chainhash128_store`).

**The ≤ 1-block path** (0..16 bytes: m = 1, one region, q = 1). The formula gives V = ℓ·z + b_1,
and V = 0 for the empty message. This is not a special case of the definition. It is the odd-q
rule, and it costs one multiply.

## 4. Interface and implementation notes

These are not part of the definition.

- `chainhash128v2(&key, data, len)` dispatches at run time: CH128_ZMM (AVX-512 VPCLMULQDQ) >
  CH128_XMM (PCLMULQDQ, with AVX2 tails when present) > CH128_PORTABLE on x86; CH128_NEON when
  the AArch64 CPU has PMULL (detected at run time, as in `chainhash128.h`). The numbering is
  `chainhash128.h`'s; there is no YMM backend. `chainhash128v2_with_backend` requires
  `chainhash128v2_has_backend`. `CHAINHASH128_PORTABLE` omits all hardware code.
- The header includes `chainhash128.h` for the GF(2^128) helpers, the finalizer and CPU
  detection; both headers can be included in one translation unit with every other ChainHash
  header.
- Any evaluation schedule gives the same V. This includes the unreduced 256-bit state, the lane
  split on ZMM (lane c carries pair c's chain) and streaming. Streaming
  (`chainhash128v2_init/update/final`) buffers one 16 KiB region, folds each full region at once
  and adds ℓ·z^{m'} at the end.
- The expanded key caches the mask table (2 KiB) and tail constants, y^1..y^8 and the ZMM lane
  tables, z, 0x87·z and the limb swaps: about 7 KiB, initialized in place.
- On x86, the outer 128×128 products use schoolbook on Intel cores with ADX and Karatsuba elsewhere
  (`-DCH128P_XSCHOOL=0/1` forces one for tests); same digest.
- `chainhash128v2_reference` is the definition read literally (direct blocks, bit-serial field
  arithmetic); the tests compare every backend with it.

## 5. Bound

For messages M ≠ M′ of at most 8L bytes, fixed before the key, Pr[H(M) = H(M′)] ≤ N(L)/2^128 with
N(L) ≤ 2L, so the score is 127 bits. The statement, the numerators and the proof status are in
[THEOREM-128v2.md](THEOREM-128v2.md).

## 6. Test vectors

- Key: key[i] = (i·73 + 11) mod 256, 192 bytes. Message: msg[i] = (i·137 + 29) mod 256.
- [`test/128v2/vectors.txt`](../test/128v2/vectors.txt) lists `len digest` for 151 lengths from 0
  to 49153 B. They cover:
  - every last-region shape q = 1..8 in regions 1, 2 and 3;
  - block boundaries (B·k ± 1) and pair-vector presence boundaries (1024m, 1024m + 112 ± 1);
  - layout-region boundaries (16384k ± 1);
  - the short-path boundaries (16/32/64/128 ± 1).
- They are produced by `pyref_2l.py`. That oracle computes blocks by the ring form
  (`algebra.py`, `pyref.py`), and its outer stage and finalizer are its own GF(2^128) code, read
  literally from section 3. `make test-128v2` checks the vectors on the reference, every backend,
  one-shot and streaming, and recomputes the short ones in Python; `make vectors` regenerates the
  whole file in Python and compares it byte for byte.
