"""CH-128/P v2 canonical oracle (independent of the C code), written from docs/SPEC-128v2.md (development name CH-128/P v2).
Block values: v1's ring formula (b_t = c*(W7, W6), W = sum_m U_m V_m mod g; algebra.py / pyref.solve,
the v1 oracle's code, which does not share code with the C kernels). Outer stage: this file's own
GF(2^128) arithmetic (carry-less product then reduction by X^128 = X^7+X^2+X+1), the two-level
formula, then the twist + degree-5 finalizer.
Usage: python3 pyref_2l.py [out]  -> vectors.txt (M = 16, the released instance); python3 pyref_2l.py --m M [out] for other block sizes."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from algebra import amul, apow, c as CC, mul64
from pyref import solve

M128 = (1 << 128) - 1
R = 8                       # blocks per level-2 region

def clmul(a, b):
    r = 0
    while b:
        if b & 1: r ^= a
        a <<= 1; b >>= 1
    return r
def gmul(a, b):             # GF(2^128), X^128 = X^7 + X^2 + X + 1
    r = clmul(a, b)
    for i in range(r.bit_length() - 1, 127, -1):
        if (r >> i) & 1: r ^= (1 << i) ^ (0x87 << (i - 128))
    return r
def gpow(a, e):
    r = 1
    for _ in range(e): r = gmul(r, a)
    return r

def block_values(Mpv, s, msg):
    """v1's block values b_1..b_m (docs/SPEC-128v2.md section 2), via the ring formula."""
    pw = [apow(s, e) for e in range(2 * Mpv + 1)]
    B = 128 * Mpv; REG = 8 * B; n = len(msg)
    q, r = divmod(n, REG)
    m = 1 if n == 0 else 8 * q + (min(8, (r + 15) // 16) if r else 0)
    def word(o): return int.from_bytes(bytes(msg[o + t] if o + t < n else 0 for t in range(8)), 'little')
    out = []
    for t in range(m):
        base = (t // 8) * REG; j = t % 8; Wsum = [0] * 8
        for mm in range(Mpv):
            if base + 1024 * mm + 16 * j >= n: continue
            ud = [word(base + 1024 * mm + 128 * i + 16 * j) for i in range(8)]
            vd = [word(base + 1024 * mm + 128 * i + 16 * j + 8) for i in range(8)]
            U = [a ^ b for a, b in zip(solve(ud), pw[2 * mm + 1])]
            V = [a ^ b for a, b in zip(solve(vd), pw[2 * mm + 2])]
            Wsum = [a ^ b for a, b in zip(Wsum, amul(U, V))]
        out.append(mul64(CC, Wsum[7]) | (mul64(CC, Wsum[6]) << 64))
    return out

def region_value(a, y):
    """c(a_1..a_q) = sum_{i=1}^{floor(q/2)} (a_i + y^(2i-1)) (a_(ceil(q/2)+i) + y^(2i)) + [q odd] a_ceil(q/2)."""
    q = len(a); h = (q + 1) // 2; c = 0
    for i in range(1, q // 2 + 1):
        c ^= gmul(a[i - 1] ^ gpow(y, 2 * i - 1), a[h + i - 1] ^ gpow(y, 2 * i))
    if q % 2: c ^= a[h - 1]
    return c

def hash_2l(Mpv, key, msg):
    s = [int.from_bytes(key[8 * i:8 * i + 8], 'little') for i in range(8)]
    y, c0, c1, c2, c3, c4, tau = [int.from_bytes(key[64 + 16 * i:80 + 16 * i], 'little') for i in range(7)]
    z = int.from_bytes(key[176:192], 'little')
    b = block_values(Mpv, s, msg)
    regions = [b[i:i + R] for i in range(0, len(b), R)]
    mp = len(regions)
    V = gmul(len(msg), gpow(z, mp))                        # length-leading: l * z^m'
    for rho, reg in enumerate(regions, 1):
        V ^= gmul(region_value(reg, y), gpow(z, mp - rho))
    v = (V + tau) & M128
    qv = gmul(v, v); rv = gmul(qv ^ c0, v ^ qv ^ c1)
    return gmul(v ^ c2, rv ^ c3) ^ c4

def lengths(Mpv):
    """Every region shape q = 1..8 in regions 1, 2 and 3; every block-size boundary (B*k -1/0/+1 inside the first
    region: block k's first byte is at 16k, a block's pair-vector m starts at 1024m + 16j); region boundaries;
    the short-path boundaries 16/32/64/128."""
    B = 128 * Mpv; REG = 8 * B
    s = {0, 1, 7, 8, 9, 15, 16, 17, 31, 32, 33, 47, 48, 49, 63, 64, 65, 80, 96, 100, 111, 112, 113, 127, 128, 129,
         255, 256, 257, 511, 512, 513, 767, 768, 1000, 1023, 1024, 1025, 1040, 1041, 1136, 1137, 2047, 2048, 2049,
         3000, 4095, 4096, 4097, 8191, 8192, 8193, 12288, 20000, 32768, 40000}
    for reg in range(3):
        for q in range(1, 9):
            s.add(REG * reg + 16 * (q - 1) + 1)           # last region has exactly q blocks (1 word each)
            s.add(REG * reg + 16 * q)                     # q blocks, last word full
        s.update({REG * reg + B - 1, REG * reg + B, REG * reg + B + 1})
    for k in range(1, 4):
        s.update({k * REG - 1, k * REG, k * REG + 1})
    for m in range(1, Mpv):                               # pair-vector presence boundaries in the first region
        s.update({1024 * m, 1024 * m + 1, 1024 * m + 16 * 7, 1024 * m + 16 * 7 + 1})
    return sorted(L for L in s if L <= 3 * REG + 200)

if __name__ == '__main__':
    args = sys.argv[1:]; Mpv = 16
    if args[:1] == ['--m']: Mpv = int(args[1]); args = args[2:]
    out = open(args[0], 'w') if args else sys.stdout
    key = bytes((i * 73 + 11) & 255 for i in range(192))
    lens = lengths(Mpv)
    msg = bytes((i * 137 + 29) & 255 for i in range(max(lens)))
    print('# CH-128/P v2 (M=%d: %d B blocks, 8 blocks per region, level-2 regions of R=8). key[i]=(i*73+11)&255 (192 bytes), msg[i]=(i*137+29)&255' % (Mpv, 128 * Mpv), file=out)
    print('# len digest(16 bytes, little-endian hex)', file=out)
    for L in lens:
        print('%d %s' % (L, hash_2l(Mpv, key, msg[:L]).to_bytes(16, 'little').hex()), file=out, flush=True)
