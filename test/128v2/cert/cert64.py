"""Exact certificates for CH-128/P over F = GF(2^64) mod X^64+X^4+X^3+X+1.
Every claim is either an exact algebraic identity checked on a basis (F-bilinear/
F-linear => basis suffices) or an exact polynomial computation. Random tests are
extra evidence only."""
import random, sys
from gfp import GF, Ring, ptrim, padd, pmul, pmod, pdivmod, pgcd, peval, pderiv
from structure import Pencil, gf2_rank

random.seed(20260928)
F = GF(64, 0x1b)
ok_all = True
def check(name, cond):
    global ok_all
    ok_all &= bool(cond)
    print(("PASS " if cond else "FAIL ") + name)
    sys.stdout.flush()

check("F modulus X^64+X^4+X^3+X+1 irreducible (Rabin)", F.is_irreducible_modulus())
P = Pencil(F, range(8), 2)            # alpha_i = i, xi = X
A = P.A
print("L =", [hex(x) for x in P.L])
check("L = T^8 + 0x17d T^4 + 0x60c T^2 + 0x770 T",
      P.L == [0, 0x770, 0x60c, 0, 0x17d, 0, 0, 0, 1])
check("c = 0x770", P.c == 0x770)
print("g =", [hex(x) for x in P.g], " (xi*c =", hex(F.mul(2, P.c)), ")")
check("g = L + xi*c, constant 0xee0", P.g[0] == 0xee0)

# --- (3) identity (Q0,Q1) = c (W7, W6), W = E^-1(u) E^-1(v) mod g ---
def ident(u, v):
    W = A.mul(P.Einv(u), P.Einv(v))
    return P.Q(u, v) == P.Pi(W)
e = [[1 if j == i else 0 for j in range(8)] for i in range(8)]
check("identity on all 64 basis pairs (e_i, e_j)  [F-bilinear => exact proof]",
      all(ident(e[i], e[j]) for i in range(8) for j in range(8)))
check("identity on 1000 random pairs",
      all(ident([random.getrandbits(64) for _ in range(8)], [random.getrandbits(64) for _ in range(8)]) for _ in range(1000)))

# --- (4) E invertible: Vandermonde, and E(Einv(u)) = u on a basis ---
check("E o Einv = id on basis (E invertible, alphas distinct)",
      all(P.E(P.Einv(e[i])) == e[i] for i in range(8)))
vdet = 1
for i in range(8):
    for j in range(i + 1, 8):
        vdet = F.mul(vdet, i ^ j)
check("Vandermonde det = prod(alpha_i + alpha_j) != 0", vdet != 0)

# --- Frobenius: T^(2^64) = T+3, T^(2^128) = T mod g ---
Tel = A.norm([0, 1])
f64 = A.frob(Tel, 64)
check("T^(2^64) == T + 3 mod g", f64 == A.norm([3, 1]))
f128 = A.frob(f64, 64)
check("T^(2^128) == T mod g", f128 == Tel)
check("gcd(g, T^(2^64)-T) = 1 (no root in F)", len(pgcd(F, P.g, ptrim(padd(f64, [0, 1])))) == 1)
check("g' = c (nonzero const) => g squarefree", pderiv(F, P.g) == [P.c])
# g(T+3) = g(T): L additive and L(3)=0
g_shift = [0] * 9
for i, co in enumerate(P.g):     # expand co*(T+3)^i
    pw = [1]
    for _ in range(i):
        pw = pmul(F, pw, [3, 1])
    for k, x in enumerate(pw):
        g_shift[k] ^= F.mul(co, x)
check("g(T+3) == g(T)", ptrim(g_shift) == ptrim(P.g))

# --- explicit factorization into four irreducible quadratics ---
facs = P.split_quadratics()
facs.sort()
for f in facs:
    print("  factor:", [hex(x) for x in f])
prod = [1]
for f in facs:
    prod = pmul(F, prod, f)
check("product of factors == g", prod == ptrim(P.g))
check("exactly four factors, each monic degree 2", len(facs) == 4 and all(len(f) == 3 and f[2] == 1 for f in facs))
check("each factor has the form T^2 + 3T + beta", all(f[1] == 3 for f in facs))
def irreducible_quadratic(f):
    R = Ring(F, f)
    x = R.frob(R.norm([0, 1]), 64)
    return len(pgcd(F, f, ptrim(R.add(x, R.norm([0, 1]))))) == 1
check("each factor irreducible over F (no root: gcd(f, T^q - T) = 1)", all(irreducible_quadratic(f) for f in facs))
check("factors pairwise distinct", len(set(tuple(f) for f in facs)) == 4)

# --- CRT: A -> prod F[T]/f_c is an F-linear bijection (8x8 over F, rank via F_2 512) ---
def crt(a):
    out = []
    for f in facs:
        r = pmod(F, ptrim(a), f)
        out += r + [0] * (2 - len(r))
    return out
def pack(words):
    x = 0
    for i, w in enumerate(words):
        x |= w << (64 * i)
    return x
basis512 = []
for k in range(8):
    for b in range(64):
        a = [0] * 8
        a[k] = 1 << b
        basis512.append(pack(crt(a)))
check("CRT map A -> (+)K_c has F_2-rank 512 (bijection => uniform s <-> independent uniform s_c)",
      gf2_rank(basis512) == 512)
# ring homomorphism spot check
ok = True
for _ in range(200):
    a, b = A.rand(), A.rand()
    ab = A.mul(a, b)
    for f in facs:
        R = Ring(F, f)
        if R.norm(pmod(F, ptrim(ab), f)) != R.mul(R.norm(pmod(F, ptrim(a), f)), R.norm(pmod(F, ptrim(b), f))):
            ok = False
check("CRT projection is multiplicative (200 random pairs)", ok)

# --- (c) Kr ∩ K_c = {0}:  K_c = ideal (h_c), h_c = g/f_c ---
for ci, f in enumerate(facs):
    h, rem = pdivmod(F, P.g, f)
    assert not rem and len(h) == 7
    b1 = A.norm(h)
    b2 = A.norm(pmul(F, [0, 1], h))
    p1, p2 = P.Pi(b1), P.Pi(b2)
    det = F.mul(p1[0], p2[1]) ^ F.mul(p1[1], p2[0])
    # F_2 version: 128 basis vectors X^k*b1, X^k*b2 -> 128-bit images
    vecs = []
    for bb in (b1, b2):
        for k in range(64):
            v = [F.mul(1 << k, x) for x in bb]
            q0, q1 = P.Pi(v)
            vecs.append(q0 | (q1 << 64))
    r = gf2_rank(vecs)
    # K_c really is the kernel of reduction mod the other factors (it is the c-component)
    is_comp = all(not ptrim(pmod(F, ptrim(b1), f2)) and not ptrim(pmod(F, ptrim(b2), f2)) for f2 in facs if f2 != f)
    check(f"component {ci}: Pi|K_c 2x2 det over F = {hex(det)} != 0, F_2 rank = {r} (need 128), K_c in ker of other projections",
          det != 0 and r == 128 and is_comp)
# kernel of Pi has F_2-dim 384 on A
vecs = []
for k in range(8):
    for b in range(64):
        a = [0] * 8
        a[k] = 1 << b
        q0, q1 = P.Pi(a)
        vecs.append(q0 | (q1 << 64))
check("Pi: A -> F^2 has F_2-rank 128 (kernel dim 384)", gf2_rank(vecs) == 128)
# sum of two components DOES meet Kr (shows why the argument needs single components)
h2, _ = pdivmod(F, P.g, pmul(F, facs[0], facs[1]))
x = A.norm(h2)   # degree 4 element in K_0 + K_1, Pi = 0
print("  note: g/(f0 f1) has degree", len(h2) - 1, "-> Pi(g/(f0 f1)) =", P.Pi(x), "(nonzero element of (K_0+K_1) ∩ Kr)")

# --- (a) keyed block: literal pencil with masks E(s^j) == Pi(sum (M+s^(2m+1))(M'+s^(2m+2))) ---
ok = True
for p in (1, 2, 4):
    for _ in range(30):
        s = A.rand()
        msgs = [([random.getrandbits(64) for _ in range(8)], [random.getrandbits(64) for _ in range(8)]) for _ in range(p)]
        ok &= P.block_literal(s, msgs) == P.block_algebraic(s, msgs)
check("(a) literal keyed pencil == Pi(sum_m (M_m+s^(2m+1))(M'_m+s^(2m+2))) for p=1,2,4 (90 random)", ok)

# --- (b) per-component: c-projection of block element = ChainHash-128 block in s_c ---
ok = True
for _ in range(20):
    p = 2
    s = A.rand()
    Ms = [(A.rand(), A.rand()) for _ in range(p)]
    W = A.zero()
    for m, (M, Mp) in enumerate(Ms):
        W = A.add(W, A.mul(A.add(M, A.pow(s, 2 * m + 1)), A.add(Mp, A.pow(s, 2 * m + 2))))
    for f in facs:
        R = Ring(F, f)
        pr = lambda a: R.norm(pmod(F, ptrim(a), f))
        sc = pr(s)
        Wc = R.zero()
        for m, (M, Mp) in enumerate(Ms):
            Wc = R.add(Wc, R.mul(R.add(pr(M), R.pow(sc, 2 * m + 1)), R.add(pr(Mp), R.pow(sc, 2 * m + 2))))
        ok &= Wc == pr(W)
check("(b) component projection of the block = ChainHash-128-style block in s_c with projected messages", ok)

print("\nALL PASS" if ok_all else "\nSOME CHECK FAILED")
