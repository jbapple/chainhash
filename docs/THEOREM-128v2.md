# ChainHash-128 v2 collision bound

**Proof status: paper proof and machine-checked certificates; not yet in Lean.** The argument is
below; the exact computations it relies on are run by `make certs` (listed at the end). Unlike
ChainHash and ChainHash-128 ([THEOREM.md](THEOREM.md), [THEOREM-128.md](THEOREM-128.md)), nothing
here is checked by a proof assistant.

## Statement

Let `q = 2^128` and let L be an integer with `1 <= L <= 2^61-1`. Fix two distinct byte strings
m, m', each at most 8L bytes long, independently of the key. The key is 192 uniformly random
bytes ([SPEC-128v2.md](SPEC-128v2.md) section 1). Then

```text
Pr[ chainhash128v2(key, m) = chainhash128v2(key, m') ] <= min(1, N(L) / 2^128),
N(L) = max( d(L) + E(p(8L)) + 1,  m'(p(8L)) + 1 ).
```

L counts 8-byte words: the bound is for messages of at most 8L bytes, so messages of up to n bytes
use L = ⌈n/8⌉.

The probability is over the key alone; the messages are chosen without seeing the key or any hash
value (the guarantee is not adaptive, and the function is not a MAC). The event is equality of all
128 output bits; empty messages, partial final words and unequal lengths are included.

- `p(ℓ)` is the block count of an ℓ-byte message: with `ℓ = 16384Q + r`, `p(0) = 1` and
  `p(ℓ) = 8Q + (r > 0 ? min(8, ceil(r/16)) : 0)`. Write `m = p(8L)` and `m' = ceil(m/8)`.
- `E(m) = 8 + m' − 1` if `m' >= 2`, and `E(m) = 2 floor(m/2)` if `m' = 1`.
- `d(L)` bounds the root count of a level-1 block difference in s: `d(1) = 1`, and
  `d(L) = min(32, 2 ceil(L/128))` for `L >= 2`.

All three are nondecreasing in L, so the bound for "at most 8L bytes" is the bound at the length
limit.

| Message limit | L | ChainHash-128 `p + d` | ChainHash-128 v2 `N(L)` |
| --- | ---: | ---: | ---: |
| 16 bytes | 2 | 2 | 3 |
| 128 bytes | 16 | 9 | 11 |
| 1 KiB | 128 | 16 | 11 |
| 4 KiB | 512 | 40 | 17 |
| 16 KiB | 2048 | 64 | 41 |
| 64 KiB | 8192 | 160 | 44 |
| 1 MiB | 131072 | 2080 | 104 |
| 16 MiB | 2097152 | 32800 | 1064 |
| 1 GiB | 134217728 | 2097184 | 65576 |
| 2^40 bytes | 137438953472 | 2147483680 | 67108904 |

Divide by 2^128. The ChainHash-128 column is [THEOREM-128.md](THEOREM-128.md)'s bound for the
released `chainhash128.h`, for comparison. These are upper bounds; they are not claimed to be
attained. The SplitMix64 seed constructor is a different key distribution and the bound is not
asserted for it.

## Score

The strength score `min_L log2(L / epsilon(L))` with `epsilon(L) = N(L)/2^128`, L in 8-byte words,
is **127 bits**: `N(L) <= 2L` for every L, with equality only at L = 1 (`N(1) = 2`).
`test/128v2/cert128p.py` checks the integer inequality for every L < 400,000 and in ±3000 windows
around 2^19 .. 2^61.

## Proof outline

1. **Level 1** (per block). The block difference is Π applied to a nonzero polynomial in s of
   degree at most 2·(present pair-vectors): a changed x-word of pair-vector m survives as the
   monomial `dU·s^(2m+2)` and a changed z-word as `dV·s^(2m+1)`, and these exponents are
   distinct. The Chinese remainder theorem splits A into four fields K_c (g is a product of four
   distinct irreducible quadratics over F), Π restricted to each K_c is injective, and a root
   count in the component s_c gives at most `d(L)/2^128`. At L = 1 only the x-word of pair-vector
   0 is present, the difference is `dU·s^2`, and squaring is a bijection: `d = 1`.
2. **Lemma A** (region). The linear-term exponents of the region polynomial in y are exactly
   {1..2f}, one per value; data × data terms and the bare value sit at y^0. So two differing
   regions with the same q give a nonzero polynomial of degree at most 2f ≤ 8 in y.
3. **Lemma B** (chain). For equal lengths, the leading z-coefficient of the difference is the
   first differing region's difference; for unequal lengths it is ℓ − ℓ' or ℓ', a nonzero
   constant. y and z are independent of s and of each other, which gives the terms `E + 1` and
   `m' + 1`.
4. **Finalizer.** The integer twist is a bijection for fixed τ, and the degree-5 finalizer adds at
   most 1/2^128, as in [THEOREM-128.md](THEOREM-128.md) step 4 (the argument holds over any
   characteristic-two field).

## Machine-checked certificates

`make certs` runs [`test/128v2/cert128p.py`](../test/128v2/cert128p.py), which checks:

0. the block certificates of [`test/128v2/cert/cert64.py`](../test/128v2/cert/cert64.py), an
   independent derivation of the block algebra: F irreducible; L, c and g; (Q0, Q1) = c·(W7, W6)
   exactly on a basis; E invertible; g a product of four distinct irreducible quadratics; the CRT
   bijection; Π restricted to each K_c injective; the literal pencil equal to the ring form;
1. Rabin irreducibility of X^128 + X^7 + X^2 + X + 1 (the outer field K);
2. the level-1 exponent classes (2M distinct exponents ≥ 1) and the root counts d(L);
3. Lemma A's exponent bookkeeping for q = 1..8;
4. the numerators and the score (N(L) ≤ 2L, equality only at L = 1).

The vectors in [`test/128v2/vectors.txt`](../test/128v2/vectors.txt) come from an independent
Python oracle and are reproduced by the C reference and every backend (`make test-128v2`).

**Gaps.** The level-1 lemma for the power key in A, Lemmas A and B and the composition are proved
on paper and backed by the exact computations above; none of it is formalized in Lean yet. The
bound needs the 192 key bytes to be uniform; the seed constructor is not covered.
