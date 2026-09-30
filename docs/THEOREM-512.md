# ChainHash-512 collision bound

**Proof status: paper proof and machine-checked certificates; not yet in Lean.** The argument is
below; the exact computations it relies on are run by `make certs` (listed at the end). Nothing
here is checked by a proof assistant.

## Statement

Let `Q = 2^512` and let L be an integer with `1 <= L <= 2^61-1`. Fix two distinct byte strings
m, m', each at most 8L bytes long, independently of the key. The key is 576 uniformly random
bytes: (s, y, z) and (τ, c0..c4) independent and uniform ([SPEC-512.md](SPEC-512.md) section 1).
Then

```text
Pr[ chainhash512(key, m) = chainhash512(key, m') ] <= min(1, N(L) / 2^512),
N(L) = max( d(L) + E(m) + 1,  m' + 1 ),   m = p(8L),   m' = ceil(m/8).
```

L counts 8-byte words: the bound is for messages of at most 8L bytes, so messages of up to n bytes
use L = ⌈n/8⌉.

The probability is over the key alone; the messages are chosen without seeing the key or any hash
value (the guarantee is not adaptive, and the function is not a MAC). The event is equality of all
512 output bits; empty messages, partial pairs and unequal lengths are included.

- `G(n) = 8 floor(n/1024) + ceil((n mod 1024)/128)` is the pair count and `p(n) = max(1, ceil(G(n)/128))`
  the block count of an n-byte message.
- `E(m) = 8 + m' − 1` if `m' >= 2`, and `E(m) = 2 floor(m/2)` if `m' = 1`.
- `d(L)` is the level-1 root count: `d(L) = 1` if `8L <= 64` (only u of pair 0 is present and
  v = 0, so the difference is `(u − u')·s^2`, whose only root is 0); otherwise
  `d(L) = min(2·G(8L), 256)`.

| Message limit | L | `N(L)` |
| --- | ---: | ---: |
| 64 bytes | 8 | 2 |
| 128 bytes | 16 | 3 |
| 1 KiB | 128 | 17 |
| 16 KiB | 2048 | 257 |
| 64 KiB | 8192 | 261 |
| 1 MiB | 131072 | 272 |
| 16 MiB | 2097152 | 392 |
| 1 GiB | 134217728 | 8456 |
| 2^40 bytes | 137438953472 | 8388872 |

Divide by 2^512. These are upper bounds; they are not claimed to be attained. Up to 128 KiB (one
region) the two-level outer stage is at most 1 worse than a direct Horner chain in y; beyond that
the numerator grows as m/8 + 8 instead of m. The seed constructor is a different key distribution
and the bound is not asserted for it.

## Score

`min_L log2(L·Q/N(L))` = **511 bits**, attained at L = 1 (`N(1) = 2`). `test/512/cert512.py` checks
`N(L) <= 2L` for every L < 200,000 and in ±3000 windows around 2^18 .. 2^61.

## Proof outline

1. **Level 1** (block, root count in s). For a block where the two messages differ, the difference
   of b_t is a polynomial in s: `Σ_q [(u_q − u'_q) s^(2q+2) + (v_q − v'_q) s^(2q+1)]` plus a constant.
   The 256 (operand, position) classes map to 256 distinct exponents 1..256, so a differing pair
   gives a nonzero coefficient on its own monomial and the polynomial is nonzero, with at most
   `d(L)` roots.
2. **Lemma A** (region). The linear exponents of the region polynomial in y are {1..2f}, each used
   once; data × data terms and the bare value sit at `y^0`.
3. **Lemma B** (chain). A length-leading Horner polynomial in the independent z: for equal lengths
   the first differing region's difference leads; for unequal lengths the leading coefficient is a
   nonzero constant.
4. **Finalizer.** The twist and the quintic add `1/Q`, as in [THEOREM.md](THEOREM.md) step 4 (the
   argument holds over any characteristic-two field).

ChainHash-512 has no pencil block, so the level-1 step is a plain root count in s; the power-key lemma in
an algebra that ChainHash-128 v2 needs does not arise.

## Machine-checked certificates

`make certs` runs [`test/512/cert512.py`](../test/512/cert512.py), which checks:

- Rabin irreducibility of `x^512 + x^8 + x^5 + x^2 + 1`: `x^(2^512) = x mod f` and
  `gcd(x^(2^256) − x, f) = 1`;
- the level-1 exponent classes: 256 (operand, position) classes → 256 distinct exponents in
  1..256, and the root counts d(1..8) = 1, d(9) = 2, d(large) = 256;
- Lemma A's bookkeeping: for q = 1..8 the linear exponents are exactly {1..2f}, data × data and the
  bare value at 0;
- the score (511 bits at L = 1) and `N(L) <= 2L` on 464,000 values of L up to 2^61, with the
  numerator table above.

The vectors in [`test/512/vectors.txt`](../test/512/vectors.txt) come from an independent Python
oracle and are reproduced by the reference and every backend (`make test-512`).

**Gaps.** The composition, Lemma A and Lemma B are proved on paper (field-generic) and backed by
the exact computations above; none of it is formalized in Lean yet. The bound needs the 576 key
bytes to be uniform; the seed constructor is not covered.
