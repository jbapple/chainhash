# ChainHash-256 specification

This defines the function in [`include/chainhash256.h`](../include/chainhash256.h), a keyed
256-bit hash: a pseudo-dot-product over GF(2^256) with a power key, 4 KiB blocks and the version 2
family's two-level outer stage ([FAMILY.md](FAMILY.md)). The collision bound is in
[THEOREM-256.md](THEOREM-256.md). Development name: PH-256 v2 with the v1.2 schedule (v2 changed
the key from 128 independent masks to the powers of one element; v1.2 changed only the schedule).

The normative reference is `ph2l_hash_ref` in the header (bit-serial), exposed as
`chainhash256_reference`. [`test/256/pyref.py`](../test/256/pyref.py) is an independent Python
implementation written from this document, and it reproduces every vector in
[`test/256/vectors.txt`](../test/256/vectors.txt).

## 1. Field

L = GF(2)[x]/(f), with f = x^256 + x^10 + x^5 + x^2 + 1. f is irreducible: a Rabin test
([`test/256/ph_field.py`](../test/256/ph_field.py)) checks x^(2^256) ≡ x mod f and
gcd(x^(2^128) − x, f) = 1, and rejects two negative controls.

An element is 32 bytes: four little-endian 64-bit limbs. Limb 0 holds the coefficients of
x^0..x^63; bit i of limb l is the coefficient of x^(64l+i). Addition is XOR, and multiplication is
polynomial multiplication mod f.

## 2. Key (288 bytes)

The key is nine elements of L, 32 bytes each, all independent and uniform, in this order:

| element | role |
|---|---|
| s | block mask generator: for pair i = 1..64, k_i = s^(2i−1) masks x_i and l_i = s^(2i) masks y_i |
| Y | level-2 key y (region formula) |
| z | region chain key |
| tau | twist |
| c_0, c_1, c_2, c_3, c_4 | finalizer |

`chainhash256_key_from_bytes` reads them in this order, each as four little-endian limbs, limb 0
first; `chainhash256_key_from_words` takes the same 72 limbs. No key is rejected. The 128 masks
are computed once at key setup (the expanded key holds them; they are derived, not key material).

`chainhash256_key_from_seed` (splitmix64 of `seed ^ 0x5048323536763200`, four outputs per element,
limb 0 first, in the same order) is for tests and vectors only; the bound is not asserted for it.
It keeps the development rule that a zero Y is replaced by 1. For key bytes, Y is used as given:
the outer stage uses Y only through the masks y^1..y^8, so Y = 0 is a valid key, and the bound is
over a uniform Y including 0. For every seed whose Y is nonzero (all but a 2^-256 fraction), the
seed key equals the key given by the same splitmix bytes.

## 3. Message layout

- Block = m = 64 pairs (x_i, y_i) of elements (4096 bytes). Region = 8 blocks (32768 bytes). Let
  n_full = ⌊n / REGION⌋ and rem = n − n_full·REGION.
- **Full regions.** Region r is the bytes [r·REGION, (r+1)·REGION). Its blocks 0..7 are
  interleaved limb-major:
  - limb l of x_i of block b is the 8 bytes at 512i + 64l + 8b;
  - limb l of y_i of block b is the 8 bytes at 512i + 256 + 64l + 8b.
- **Tail.** The last rem bytes form ⌈rem/BLOCK⌉ blocks, each zero-padded to BLOCK. Within a tail
  block, pair i = 8g + t has:
  - x limb l at 512g + 64l + 8t;
  - y limb l at 512g + 256 + 64l + 8t.

  A tail block therefore has the same 512-byte row shape as a region, with lanes = pairs instead
  of blocks.

## 4. Block value

c = Σ_{i=1}^{m} (x_i + k_i)·(y_i + l_i) in L, with the power-key masks k_i = s^(2i−1), l_i = s^(2i)
(pairs numbered 1..64 here; the layout of section 3 numbers them 0..63).

## 5. Outer stage (two-level) and finalizer

1. **Block values.** b_1..b_m are the block values of section 4 in message order: region 0 blocks
   0..7, region 1 blocks 0..7, …, then the tail blocks. So m = 8·n_full + ⌈rem/BLOCK⌉.
   **Exception:** the empty message has m = 1 and b_1 = 0.
2. **Regions.** R = 8 values each; m′ = ⌈m/8⌉. Every region has q = 8 values except possibly the
   last, which has q = m − 8(m′ − 1) ∈ {1..8}.
3. **Region value** (the block formula applied to the block values, with key y = Y). For values
   a_1..a_q, with h = ⌈q/2⌉ and f = ⌊q/2⌋:

       c(a_1..a_q) = Σ_{i=1}^{f} (a_i + Y^{2i−1})·(a_{h+i} + Y^{2i})  +  [q odd]·a_h

   For q = 8, value i is paired with value 4 + i. For odd q, the middle value a_h is unpaired and
   enters bare.
4. **Chain.** V = n·z^{m′} + Σ_{ρ=1}^{m′} c_ρ·z^{m′−ρ}: a Horner polynomial in the independent key
   z with the length leading, where n is the byte length as an element (limb 0 = n). For a message
   of at most one block this is V = n·z + b_1, and V = 0 for n = 0.
5. **Finalizer.** X = V +_Z tau (256-bit integer addition mod 2^256 of the limb vectors), G = X²,
   t = (G + c_0)·(X + G + c_1), out = (X + c_2)·(t + c_3) + c_4.
6. **Digest.** out as 32 bytes, limb 0 first, little-endian.

## 6. Interface and implementation notes

These are not part of the definition.

- Backends, all with identical digests: CH256_AVX512 (AVX-512 F/VL/BW/DQ/VBMI2 + VPCLMULQDQ +
  GFNI), CH256_PCLMUL (PCLMULQDQ + SSE4.1), CH256_NEON (AArch64 PMULL, when the compiler targets
  +crypto), CH256_PORTABLE (the bit-serial reference). The backend is chosen when the key is
  initialized (`*_with_backend` variants pick one; an unavailable one falls back to the best).
  x86 kernels carry their own target attributes, so the header builds without `-march`; build
  with `-O3` under gcc. `CHAINHASH256_PORTABLE` omits all hardware code.
- The expanded key holds pointers into itself in some backends: initialize it in place and do not
  copy it.
- Field products in the kernels are Karatsuba² over 64-bit limbs (9 carry-less products),
  accumulated unreduced; the reference multiplies schoolbook with a bit-level reduction and never
  uses that evaluation. On x86 the finalizer and tail have a Zen 4 latency path;
  `-DPHX_ZEN=0/1` forces either for tests (same digest).
- Short inputs (at most one block) take a latency schedule: the tail pairs' products are formed
  directly and folded once, and the finalizer's independent products are started early. Same
  digests; the property and page-edge tests cover it.
- Streaming (`chainhash256_init/update/final`) folds every full region into V (V ← V·z + c) as
  soon as it is complete; `final` folds the tail region and adds n·z^{m′}. Before any region
  completes, `final` applies the one-shot function to the buffered bytes.
- The header needs a little-endian host and GCC or Clang (`unsigned __int128`, statement
  expressions in the AVX-512 kernel).

## 7. Conformance

- [`test/256/vectors.txt`](../test/256/vectors.txt): lines `64 key n digest`. The message is
  m[j] = (137·j + 29) mod 256; `key` is a decimal seed (`chainhash256_key_from_seed`) or `R`, the
  raw key b[i] = (73·i + 11) mod 256, i = 0..287, which the C test gives both as bytes and as words.
- `make test-256` runs every backend against the reference (one-shot and streaming), the vectors,
  and the Python cross-check; `make certs` runs the certificates of
  [THEOREM-256.md](THEOREM-256.md).
