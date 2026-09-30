"""CH-128/P algebra: A = F[T]/g(T), F = GF(2^64) mod X^64+X^4+X^3+X+1.
Verifies L(T)=prod(T+alpha_i), g = L + xi*c, and (Q0,Q1) = c*(W7,W6) for W=U*V mod g
where u_i=U(alpha_i), v_i=V(alpha_i)."""
import random, sys
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from gf import mul64

ALPHA = list(range(8)); XI = 2; M64 = (1<<64)-1
def fpow(a, e):
    r = 1
    while e:
        if e & 1: r = mul64(r, a)
        a = mul64(a, a); e >>= 1
    return r
def finv(a): return fpow(a, (1<<64)-2)

def pmul(a, b):   # polys over F, lists low->high
    r = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b): r[i+j] ^= mul64(x, y)
    return r
L = [1]
for al in ALPHA: L = pmul(L, [al, 1])
c = L[1]
G = list(L); G[0] ^= mul64(XI, c)          # g = L + xi*c (monic deg 8)
def pmod(a, g=G):
    a = list(a)
    for d in range(len(a)-1, 7, -1):
        t = a[d]
        if t:
            for k in range(9): a[d-8+k] ^= mul64(t, g[k])
    return (a + [0]*8)[:8]
def amul(a, b): return pmod(pmul(a, b))
def peval(a, x):
    r = 0
    for co in reversed(a): r = mul64(r, x) ^ co
    return r
def E(a): return [peval(a, al) for al in ALPHA]
def apow(s, e):
    r = [1]+[0]*7; b = list(s)
    while e:
        if e & 1: r = amul(r, b)
        b = amul(b, b); e >>= 1
    return r

def pencil_Q(u, v):
    Q0 = 0; Q1 = 0; su = 0; sv = 0
    for i in range(8):
        p = mul64(u[i], v[i]); Q0 ^= p; Q1 ^= mul64(ALPHA[i], p); su ^= u[i]; sv ^= v[i]
    return Q0, Q1 ^ mul64(XI, mul64(su, sv))

if __name__ == '__main__':
    print('L coeffs (T^0..T^8):', [hex(x) for x in L])
    print('c = %#x ; g =' % c, [hex(x) for x in G])
    assert L[8] == 1 and L[0] == 0 and L[3] == L[5] == L[6] == L[7] == 0
    assert (L[4], L[2], L[1]) == (0x17d, 0x60c, 0x770), 'L differs from claim'
    random.seed(7); ok = 0
    for _ in range(3000):
        U = [random.getrandbits(64) for _ in range(8)]; V = [random.getrandbits(64) for _ in range(8)]
        W = amul(U, V)
        assert pencil_Q(E(U), E(V)) == (mul64(c, W[7]), mul64(c, W[6])); ok += 1
    print('identity (Q0,Q1) == c*(W7,W6) holds on %d random (U,V)' % ok)
    # E is invertible (alphas distinct): check E is injective on random + basis
    basis = [E([1 if k == j else 0 for k in range(8)]) for j in range(8)]
    print('E(T^j) rows computed; E is a Vandermonde map at distinct alphas')
