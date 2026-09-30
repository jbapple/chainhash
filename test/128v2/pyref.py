"""Independent CH-128/P oracle via the ring formula b_t = c*(W7,W6), W = sum_m U_m V_m mod g,
U_m = E^-1(u_data) + s^(2m+1), V_m = E^-1(v_data) + s^(2m+2). Own GF(2^128) outer + finalizer.
Usage: python3 pyref.py M  -> prints vectors for the test lengths (message byte i = (i*137+29)&255, seed 123)."""
import sys
from algebra import amul, apow, E, c, mul64, fpow, finv, ALPHA
M64 = (1 << 64) - 1; M128 = (1 << 128) - 1

def mul128(a, b):
    r = 0
    while b:
        if b & 1: r ^= a
        b >>= 1; a <<= 1
        if a >> 128: a = (a & M128) ^ 0x87
    return r

def splitmix_bytes(seed, nwords):
    out = b''
    for _ in range(nwords):
        seed = (seed + 0x9e3779b97f4a7c15) & M64; z = seed
        z = ((z ^ (z >> 30)) * 0xbf58476d1ce4e5b9) & M64; z = ((z ^ (z >> 27)) * 0x94d049bb133111eb) & M64; z ^= z >> 31
        out += z.to_bytes(8, 'little')
    return out

# inverse Vandermonde: E(U)_i = sum_k U_k alpha_i^k ; solve by Gaussian elimination over F
V = [[fpow(al, k) if not (al == 0 and k == 0) else 1 for k in range(8)] for al in ALPHA]
def solve(rhs):
    A = [row[:] + [r] for row, r in zip(V, rhs)]
    for col in range(8):
        piv = next(r for r in range(col, 8) if A[r][col])
        A[col], A[piv] = A[piv], A[col]
        inv = finv(A[col][col]); A[col] = [mul64(x, inv) for x in A[col]]
        for r in range(8):
            if r != col and A[r][col]:
                f = A[r][col]; A[r] = [x ^ mul64(f, y) for x, y in zip(A[r], A[col])]
    return [A[k][8] for k in range(8)]

def hash_py(Mpv, key, msg):
    s = [int.from_bytes(key[8*i:8*i+8], 'little') for i in range(8)]
    w = [int.from_bytes(key[64+16*i:80+16*i], 'little') for i in range(7)]
    y, c0, c1, c2, c3, c4, tau = w
    pw = [apow(s, e) for e in range(2 * Mpv + 1)]
    B = 128 * Mpv; REG = 8 * B; n = len(msg)
    q, r = divmod(n, REG)
    p = 1 if n == 0 else 8 * q + (min(8, (r + 15) // 16) if r else 0)
    def byte(o): return msg[o] if o < n else 0
    def word(o): return int.from_bytes(bytes(byte(o + t) for t in range(8)), 'little')
    Vh = n
    for t in range(p):
        base = (t // 8) * REG; j = t % 8; Wsum = [0] * 8
        for m in range(Mpv):
            if base + 1024 * m + 16 * j >= n: continue
            ud = [word(base + 1024*m + 128*i + 16*j) for i in range(8)]
            vd = [word(base + 1024*m + 128*i + 16*j + 8) for i in range(8)]
            U = [a ^ b for a, b in zip(solve(ud), pw[2*m+1])]
            Vv = [a ^ b for a, b in zip(solve(vd), pw[2*m+2])]
            Wsum = [a ^ b for a, b in zip(Wsum, amul(U, Vv))]
        bt = mul64(c, Wsum[7]) | (mul64(c, Wsum[6]) << 64)
        Vh = mul128(Vh, y) ^ bt
    v = (Vh + tau) & M128
    qv = mul128(v, v); rv = mul128(qv ^ c0, v ^ qv ^ c1)
    return mul128(v ^ c2, rv ^ c3) ^ c4

if __name__ == '__main__':
    Mpv = int(sys.argv[1]); REG = 1024 * Mpv
    key = splitmix_bytes(123, 22)
    lens = [0, 1, 16, 17, 100, 128, 129, 1000, 1024, 1025, 1500, REG - 1, REG, REG + 17, 2 * REG + 1100]
    buf = bytes((i * 137 + 29) & 255 for i in range(max(lens)))
    for L in lens:
        h = hash_py(Mpv, key, buf[:L])
        print('%d %016x %016x' % (L, h & M64, h >> 64), flush=True)
