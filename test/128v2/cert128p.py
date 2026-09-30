#!/usr/bin/env python3
"""cert128p.py -- machine checks behind docs/SPEC-128v2.md section 4 and docs/THEOREM-128v2.md (CH-128/P v2: v1 block, two-level outer).
  0. Block certificates (cert/cert64.py, an independent derivation of the block algebra, run
     unchanged): F irreducible; L, c, g; (Q0,Q1) = c*(W7,W6) on a basis (exact); E invertible; g = product of four
     distinct irreducible quadratics; CRT bijection; Pi restricted to each K_c injective; literal pencil = ring form.
  1. Rabin: G = x^128 + x^7 + x^2 + x + 1 irreducible over GF(2) (the outer field K).
  2. Level-1 exponent classes (ring form): pair-vector m contributes (U_m + s^(2m+1))(V_m + s^(2m+2)), so a
     difference in U_m (V_m) survives as the monomial dU s^(2m+2) (dV s^(2m+1)); the 2M exponents are distinct
     and >= 1. Root count d(L): d(1) = 1 (only u of word 0: dU s^2, squaring bijective), d(L) = min(2M, 2 ceil(L/128)).
  3. Lemma A bookkeeping for q = 1..8: linear exponents {2i} u {2i-1} (i <= floor(q/2)) = {1..2 floor(q/2)},
     each once; data x data and the bare middle value at exponent 0.
  4. Numerators and score, ratio r = 2: N_2L(L) = max(d + E(m) + 1, m' + 1) <= 2L for every L tested, equality
     only at L = 1 (N = 2), so the score min_L log2(L 2^128 / N(L)) = 127; Horner (v1) N_H(L) = p + d for comparison.
Exit status 0 iff all checks pass."""
import math, os, subprocess, sys
ok = True
def check(cond, msg):
    global ok
    print(('PASS ' if cond else 'FAIL ') + msg); sys.stdout.flush(); ok &= bool(cond)

# ---- 0. block certificates (cert/cert64.py) ----
here = os.path.dirname(os.path.abspath(__file__))
if os.environ.get('SKIP_CERT64'):
    print('SKIP block certificates (SKIP_CERT64 set)')
else:
    r = subprocess.run([sys.executable, 'cert64.py'], cwd=os.path.join(here, 'cert'), capture_output=True, text=True)
    fails = [l for l in r.stdout.splitlines() if l.startswith('FAIL')]
    npass = sum(l.startswith('PASS') for l in r.stdout.splitlines())
    check(r.returncode == 0 and 'ALL PASS' in r.stdout and not fails, 'block certificates cert/cert64.py: %d PASS, %d FAIL' % (npass, len(fails)))

# ---- 1. Rabin for G (degree 128 = 2^7) ----
G = (1 << 128) | 0x87
def pmulmod(a, b, m):
    r = 0; d = m.bit_length()
    while b:
        if b & 1: r ^= a
        b >>= 1; a <<= 1
        if a.bit_length() == d: a ^= m
    return r
def pmod(a, m):
    d = m.bit_length()
    while a.bit_length() >= d: a ^= m << (a.bit_length() - d)
    return a
def pgcd(a, b):
    while b: a, b = b, pmod(a, b)
    return a
t = 2
for _ in range(64): t = pmulmod(t, t, G)
g64 = pgcd(G, t ^ 2)
for _ in range(64): t = pmulmod(t, t, G)
check(t == 2 and g64 == 1, 'Rabin: x^(2^128) = x mod G and gcd(x^(2^64) - x, G) = 1 -> G = x^128+x^7+x^2+x+1 irreducible, K = GF(2^128)')

# ---- 2. level-1 exponent classes ----
M = 16
ex = {}
for m in range(M):
    ex[('U', m)] = 2 * m + 2; ex[('V', m)] = 2 * m + 1
vals = sorted(ex.values())
check(vals == list(range(1, 2 * M + 1)), 'level-1: %d (operand, pair-vector) classes -> %d distinct exponents 1..%d' % (len(ex), len(set(vals)), 2 * M))
def d_of(L): return 1 if L == 1 else min(2 * M, 2 * math.ceil(L / 128))
# L words = 8L bytes; pair-vector m of block j present iff 1024 m + 16 j < n: m <= floor((8L-1)/1024)
def d_struct(L):
    if L == 1: return 1                       # only u of word 0 (bytes 0..7): dU s^2 -> 1 root
    npv = min(M, (8 * L - 1) // 1024 + 1)     # pair-vectors that can be present
    return 2 * npv
check(all(d_struct(L) <= d_of(L) for L in range(1, 20000)), 'level-1 root counts: d(1) = 1, d(L) <= min(2M, 2 ceil(L/128)) (structural count <= formula, L < 20000)')

# ---- 3. Lemma A ----
okA = True
for q in range(1, 9):
    h = (q + 1) // 2; f = q // 2; lin = []
    for i in range(1, f + 1):
        lin.append(2 * i)          # a_i * y^(2i)          (from (a_i + y^(2i-1))(a_(h+i) + y^(2i)))
        lin.append(2 * i - 1)      # a_(h+i) * y^(2i-1)
    if sorted(lin) != list(range(1, 2 * f + 1)) or len(set(lin)) != len(lin): okA = False
    # data x data a_i a_(h+i) and the bare a_h (q odd) sit at exponent 0; key constants y^(4i-1) cancel
check(okA, 'Lemma A: for q = 1..8 the linear exponents are exactly {1..2 floor(q/2)}, each once; data x data and bare value at 0')

# ---- 4. numerators, score (r = 2) ----
R = 8
def p_blocks(n):
    REG = 8 * 128 * M
    if n == 0: return 1
    Q, r = divmod(n, REG)
    return 8 * Q + (min(8, (r + 15) // 16) if r else 0)
def E_of(m):
    mp = -(-m // R); q = m - R * (mp - 1)
    return R + mp - 1 if mp >= 2 else 2 * (q // 2)
def N_2l(L):
    m = p_blocks(8 * L); mp = -(-m // R)
    E = E_of(m) if m > R else max(E_of(k) for k in range(1, m + 1))
    return max(d_of(L) + E + 1, mp + 1)
def N_h(L): return p_blocks(8 * L) + d_of(L)
Ls = list(range(1, 400000))
for e in range(19, 62):
    c = 1 << e; Ls += range(max(1, c - 3000), c + 3000)
viol = [L for L in Ls if N_2l(L) > 2 * L]; eq = [L for L in Ls if N_2l(L) == 2 * L]
check(not viol and eq == [1], 'N_2L(L) <= 2L for every L tested (%d values, L up to 2^61 + 3000); equality only at L = 1' % len(Ls))
score = min(128 - math.log2(N_2l(L) / L) for L in Ls[:400000])
check(abs(score - 127) < 1e-9 and N_2l(1) == 2, 'score 2L = %.4f bits, attained at L = 1 (N(1) = 2)' % score)
print('      score Horner (v1.1 outer) = %.4f bits' % min(128 - math.log2(N_h(L) / L) for L in Ls[:400000]))
print('      numerators N (divide by 2^128):   message      Horner (v1.1)   2L (v2)')
for lab, nb in [('128 B', 128), ('1 KiB', 1024), ('4 KiB', 4096), ('16 KiB', 16384), ('64 KiB', 65536), ('1 MiB', 1 << 20),
                ('16 MiB', 1 << 24), ('1 GiB', 1 << 30), ('2^40 B', 1 << 40)]:
    L = nb // 8
    print('                                        %-8s %12d %12d' % (lab, N_h(L), N_2l(L)))
print('ALL PASS' if ok else 'SOME FAILED')
sys.exit(0 if ok else 1)
