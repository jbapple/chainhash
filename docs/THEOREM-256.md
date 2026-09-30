# ChainHash-256 collision bound

**Proof status: paper proof and machine-checked certificates; not yet in Lean.** The argument is
below; the exact computations it relies on are run by `make certs` (listed at the end). Nothing
here is checked by a proof assistant.

## Statement

Let `q = 2^256` and let L be an integer with `1 <= L <= 2^61-1`. Fix two distinct byte strings
m, m', each at most 8L bytes long, independently of the key. The key is 288 uniformly random
bytes ([SPEC-256.md](SPEC-256.md) section 2): nine elements of GF(2^256), of which s generates the
128 block masks s^1..s^128. Then

```text
Pr[ chainhash256(key, m) = chainhash256(key, m') ] <= min(1, N(L) / 2^256),
N(L) = max( d(8L) + E(p(8L)) + 1,  m'(p(8L)) + 1 ).
```

L counts 8-byte words: the bound is for messages of at most 8L bytes, so messages of up to n bytes
use L = ⌈n/8⌉.

The probability is over the key alone; the messages are chosen without seeing the key or any hash
value (the guarantee is not adaptive, and the function is not a MAC). The event is equality of all
256 output bits; empty messages, partial blocks and unequal lengths are included.

- `p(ℓ)` is the block count of an ℓ-byte message: with 4096-byte blocks and 32768-byte regions,
  `p(0) = 1` and `p(ℓ) = 8 floor(ℓ/32768) + ceil((ℓ mod 32768)/4096)`. Write `m = p(8L)` and
  `m' = ceil(m/8)`.
- `E(m) = 8 + m' − 1` if `m' >= 2`, and `E(m) = 2 floor(m/2)` if `m' = 1`.
- `d(n)` is the level-1 root count for n-byte messages: with `G = floor((n−1)/512)` and
  `P(n) = 8G + min(8, ceil((n − 512G)/8))` the highest pair index that carries message bytes,
  `d(n) = P(n)` for `n <= 256`, `d(n) = 2 P(n)` for `256 < n < 4096`, and `d(n) = 128` from one full
  block on (`d(0) = 1`).

| Message limit | L | `d(8L)` | `N(L)` |
| --- | ---: | ---: | ---: |
| 8 bytes | 1 | 1 | 2 |
| 64 bytes | 8 | 8 | 9 |
| 256 bytes | 32 | 8 | 9 |
| 1 KiB | 128 | 32 | 33 |
| 4 KiB | 512 | 128 | 129 |
| 16 KiB | 2048 | 128 | 133 |
| 64 KiB | 8192 | 128 | 138 |
| 1 MiB | 131072 | 128 | 168 |
| 16 MiB | 2097152 | 128 | 648 |
| 1 GiB | 134217728 | 128 | 32904 |
| 2^40 bytes | 137438953472 | 128 | 33554568 |

Divide by 2^256. These are upper bounds; they are not claimed to be attained. The seed
constructor is a different key distribution and the bound is not asserted for it.

## Score

The strength score `min_L log2(L / epsilon(L))` with `epsilon(L) = N(L)/2^256`, L in 8-byte words,
is **255 bits**: `N(L) <= 2L` for every L, with equality only at L = 1 (`N(1) = 2`).
`test/256/ph256_bounds.py` checks this for every L < 2^20 and in windows around 2^20..2^61, and
prints the numerators next to those of the independent-mask variant (d = 1 at every length) that
the power key replaced: the power key costs at most a factor 128 in the numerator (7 bits at 1 MiB)
and saves 4 KiB of key.

## Proof outline

1. **Level 1** (power key). The block value is `c(s) = Σ_i (x_i + s^(2i−1))(y_i + s^(2i))`, pairs
   i = 1..64. Expanded, `c(s) = Σ x_i y_i + Σ (x_i s^(2i) + y_i s^(2i−1)) + Σ s^(4i−1)`: every data
   word owns one exponent (x_i → 2i, y_i → 2i−1, all distinct and ≥ 1) and every data × data term
   sits at `s^0`. For two differing blocks the difference `Δc(s) = Σ Δ(x_i y_i) + Σ (Δx_i s^(2i) +
   Δy_i s^(2i−1))` is therefore a nonzero polynomial of degree at most `2P`, P the highest pair index
   holding message bytes, so it has at most `2P` roots s, and at most 128 for a full block. When only
   x words carry bytes (n ≤ 256: the first tail group's x rows), `Δc = (Σ √Δx_i s^i)²` has at most
   `P` roots, and a 1–8-byte difference is `Δx s²` with the single root s = 0. Hence `d(n)` above.
   s is uniform and independent of y, z and the finalizer key.
2. **Lemma A** (region). For region values `a ≠ a'` with the same q, `Δc(y)` is a nonzero polynomial
   of degree at most 2f: each paired value has its own linear exponent in 1..2f, and every
   data × data term, together with the bare value, sits at `y^0`.
3. **Lemma B** (chain). For equal lengths, the first differing region's `Δc` is the leading
   z-coefficient. For unequal lengths, the `z^{m'}` coefficient is a nonzero constant (the length
   or the length difference). y and z are independent of the masks and of each other.
4. **Finalizer.** The integer twist is a bijection for fixed tau, and the degree-5 finalizer adds at
   most `1/2^256`, as in [THEOREM.md](THEOREM.md) step 4 (the argument holds over any
   characteristic-two field).

Using an independent z rather than a power of y matters: `w = y^8` is not injective across lengths,
and `w = y^9` costs score.

## Machine-checked certificates

`make certs` runs [`test/256/certs.sh`](../test/256/certs.sh), which checks, on the as-built
reference ([`test/256/ph_ref.h`](../test/256/ph_ref.h), the specification section of the header
compiled at block sizes of 1, 2 and 4 KiB, and at m = 8 pairs for the interpolation):

- **field**: Rabin irreducibility of `x^256 + x^10 + x^5 + x^2 + 1` (`ph_field.py`), which rejects two
  negative controls; the C field product agrees with an independent Python product on 300
  products (`ph_oracle.py`);
- **level 1** (`cert_pk.c`): the as-built block value (`ph2_derive` + `ph_block`), as a polynomial
  in s over GF(2^256), is interpolated from 4m nodes and checked at 8 more, and must equal the
  exponent-class expansion above coefficient by coefficient; adversarial differences (single x or
  y words at four positions, a full pair, several pairs, dense, x-only with y = 0) give nonzero
  difference polynomials of degree ≤ 2m with the predicted coefficients, and the 1-word difference
  is exactly `Δx s²`; at m = 8 and at the released m = 64. A negative control (a derivation with
  colliding exponents, `l_j = k_j`) is rejected;
- **level 2** (`cert_2l.c`): the region polynomial equals Lemma A's expansion for every q = 1..8,
  including adversarial differences (single values, the bare value, a zero partner, a full pair);
  the outer polynomial in z is length-leading with the region values as coefficients, for
  m' = 1..4 with odd last q; the ≤ 1-block path; a negative control (colliding exponents) is
  rejected;
- **bound** (`ph256_bounds.py`): `N(L) <= 2L`, score 255.

The vectors in [`test/256/vectors.txt`](../test/256/vectors.txt) are reproduced by an independent
Python implementation (`pyref.py`) and by the reference and every backend (`make test-256`).

**Gaps.** Lemmas A and B and the composition are proved on paper and backed by the exact
computations above; none of it is formalized in Lean yet. The bound needs the 288 key bytes to be
uniform; the seed constructor is not covered.
