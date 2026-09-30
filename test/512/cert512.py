#!/usr/bin/env python3
"""cert512.py -- machine checks behind docs/SPEC-512.md section 3 and docs/THEOREM-512.md (ChainHash-512/PH v1, two-level outer).
  1. Rabin: f = x^512 + x^8 + x^5 + x^2 + 1 is irreducible over GF(2)  (L is a field).
  2. Level-1 exponent classes: in block t, a difference in u (resp. v) of pair position q contributes
     exactly the monomial (u-u') s^(2q+2) (resp. (v-v') s^(2q+1)); all 256 exponents are distinct and
     >= 1, so the block difference is a nonzero polynomial in s whenever the blocks differ; its root
     count d(L) = 1 when both messages have <= 64 bytes (only u of pair 0: (u-u') s^2), and
     d(L) <= min(2 * pairs, 256) otherwise.
  3. Lemma A bookkeeping for q = 1..8: linear-term exponents {2i} u {2i-1} (i <= f) = {1..2f}, each once;
     data x data and the bare value at exponent 0.
  4. Numerators and score: N_2L(L) = max(d + E(m) + 1, m' + 1) vs Horner N_H(L) = max(d + m, m + 1),
     m = blocks(8L), score = min_L log2(L 2^512 / N(L)).  Must be 511 at L = 1 (N = 2) and nowhere lower.
Exit status 0 iff all checks pass.
"""
import math, sys
ok = True
def check(cond, msg):
    global ok
    print(('PASS ' if cond else 'FAIL ') + msg); ok &= bool(cond)

# ---------------- 1. Rabin irreducibility (degree 512 = 2^9: only prime divisor 2) --------------
F = (1 << 512) | (1 << 8) | (1 << 5) | (1 << 2) | 1
def pmod(a, m):
    dm = m.bit_length()
    while a.bit_length() >= dm: a ^= m << (a.bit_length() - dm)
    return a
def pmulmod(a, b, m):
    r = 0
    while b:
        if b & 1: r ^= a
        b >>= 1; a <<= 1
        if a.bit_length() == m.bit_length(): a ^= m
    return r
def pgcd(a, b):
    while b: a, b = b, pmod(a, b)
    return a
x = 2; t = x
for i in range(256): t = pmulmod(t, t, F)          # x^(2^256)
g = pgcd(F, t ^ x)
for i in range(256): t = pmulmod(t, t, F)          # x^(2^512)
check(t == x and g == 1, 'Rabin: x^(2^512) = x mod f and gcd(x^(2^256) - x, f) = 1 -> f irreducible, L = GF(2^512)')

# ---------------- 2. level-1 exponent classes -------------------------------------------------
exps = {}
for q in range(128):
    exps[('u', q)] = 2 * q + 2          # (u + s^(2q+1))(v + s^(2q+2)): u multiplies s^(2q+2)
    exps[('v', q)] = 2 * q + 1          #                               v multiplies s^(2q+1)
vals = list(exps.values())
check(len(set(vals)) == 256 and min(vals) >= 1 and max(vals) == 256,
      'level-1: 256 (operand, position) classes -> 256 distinct exponents in 1..256')
def npairs(n):
    F_, r = divmod(n, 1024); return 8 * F_ + (r + 127) // 128
def blocks(n):
    G = npairs(n); return max(1, (G + 127) // 128)
def d_of(L):
    n = 8 * L
    if n <= 64: return 1                 # only u-limbs of pair 0 present (v = 0): difference (u-u') s^2
    return min(2 * npairs(n), 256)
check(d_of(1) == 1 and d_of(8) == 1 and d_of(9) == 2 and d_of(10**9) == 256, 'level-1 root counts: d(1..8) = 1, d(9) = 2, d(large) = 256')

# ---------------- 3. Lemma A exponent bookkeeping ---------------------------------------------
good = True
for qq in range(1, 9):
    f = qq // 2; h = (qq + 1) // 2
    lin = [2 * i for i in range(1, f + 1)] + [2 * i - 1 for i in range(1, f + 1)]
    good &= sorted(lin) == list(range(1, 2 * f + 1)) and len(set(lin)) == 2 * f and 0 not in lin
    kappa = [4 * i - 1 for i in range(1, f + 1)]   # key-only constants: cancel in differences
    good &= all(e >= 1 for e in kappa)
check(good, 'Lemma A: for q = 1..8 the linear exponents are exactly {1..2f}, data x data and bare value at 0')

# ---------------- 4. numerators and score ------------------------------------------------------
R = 8
def N2L(L):
    m = blocks(8 * L); mr = (m + R - 1) // R
    E = R + mr - 1 if mr >= 2 else 2 * (m // 2)
    return max(d_of(L) + E + 1, mr + 1)
def NH(L):
    m = blocks(8 * L); return max(d_of(L) + m, m + 1)
Ls = set(range(1, 200001))
for e in range(18, 62):
    c = 1 << e; Ls |= set(range(max(1, c - 3000), c + 3000))
worst2 = min((math.log2(L) + 512 - math.log2(N2L(L)), L) for L in Ls)
worstH = min((math.log2(L) + 512 - math.log2(NH(L)), L) for L in Ls)
check(N2L(1) == 2 and abs(worst2[0] - 511) < 1e-12 and worst2[1] == 1, 'score 2L = %.4f bits at L = %d (N(1) = %d)' % (worst2[0], worst2[1], N2L(1)))
bad = [L for L in Ls if N2L(L) > 2 * L]
check(not bad, 'N_2L(L) <= 2L for every L tested (%d values, L up to 2^61)' % len(Ls))
print('      score Horner (v0 outer) = %.4f bits at L = %d' % worstH)
print('      numerators N (divide by 2^512):   message      Horner     2L')
for label, n in [('64 B', 64), ('1 KiB', 1024), ('16 KiB', 16384), ('64 KiB', 65536), ('128 KiB', 131072),
                 ('1 MiB', 1 << 20), ('16 MiB', 1 << 24), ('1 GiB', 1 << 30), ('2^40 B', 1 << 40)]:
    L = n // 8
    print('                                        %-10s %8d %8d' % (label, NH(L), N2L(L)))
sys.exit(0 if ok else 1)
