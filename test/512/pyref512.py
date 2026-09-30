#!/usr/bin/env python3
"""pyref512.py -- independent oracle for ChainHash-512 (docs/SPEC-512.md; development name PH-512 v1), written
from the definition with Python integers (no code shared with the C header).
usage: pyref512.py [out]  (default: standard output) -- the vectors of test/512/vectors.txt.

  field      L = GF(2)[x]/(x^512 + x^8 + x^5 + x^2 + 1), elements = ints < 2^512 (bit i = x^i)
  key        576 bytes = 9 little-endian 512-bit words: s, y, tau, c0, c1, c2, c3, c4, z
  pairs      a message of n bytes is cut into 1 KiB chunks.  In a full chunk, u-limb j of pair i
             (i < 8) is the 8 bytes at 64j + 8i and v-limb j is at 512 + 64j + 8i; the final partial
             chunk (if any) holds pairs back to back, 128 bytes each (u = first 64, v = next 64),
             zero padded.  Pair g of the message belongs to block g // 128 at position q = g % 128.
  block      b_t = sum over its pairs (u + s^(2q+1)) (v + s^(2q+2)); the empty message has one
             block, b_1 = 0.  m = max(1, ceil(pairs / 128)).
  region     regions of 8 block values (last one q = m - 8(m'-1) values); h = ceil(q/2), f = q // 2:
             c = sum_{i=1..f} (a_i + y^(2i-1)) (a_{h+i} + y^(2i)) + [q odd] a_h
  outer      V = n z^m' + sum_r c_r z^(m'-r)
  final      v = (V + tau) mod 2^512 (integer addition), H = (v+c2)((v^2+c0)(v+v^2+c1)+c3)+c4,
             output = H as 64 little-endian bytes.
"""
import sys
F = (1 << 512) | (1 << 8) | (1 << 5) | (1 << 2) | 1

def clmul(a, b):
    # carry-less product via an 8-bit window table of a
    T = [0] * 256
    for x in range(1, 256):
        hb = x.bit_length() - 1
        T[x] = T[x ^ (1 << hb)] ^ (a << hb)
    r, sh = 0, 0
    while b:
        r ^= T[b & 255] << sh
        b >>= 8; sh += 8
    return r

def mod(r):
    while r.bit_length() > 512:
        d = r.bit_length() - 513
        r ^= F << d
    return r

def mul(a, b):
    return mod(clmul(a, b))

def word(bs):
    return int.from_bytes(bs, 'little')

def pairs(msg):
    n = len(msg); full = n // 1024; out = []
    for c in range(full):
        ch = msg[1024 * c:1024 * c + 1024]
        for i in range(8):
            u = sum(word(ch[64 * j + 8 * i:64 * j + 8 * i + 8]) << (64 * j) for j in range(8))
            v = sum(word(ch[512 + 64 * j + 8 * i:512 + 64 * j + 8 * i + 8]) << (64 * j) for j in range(8))
            out.append((u, v))
    rest = msg[1024 * full:]
    for o in range(0, len(rest), 128):
        pr = rest[o:o + 128] + bytes(128 - len(rest[o:o + 128]))
        out.append((word(pr[:64]), word(pr[64:])))
    return out

def hash512(key, msg):
    w = [word(key[64 * e:64 * e + 64]) for e in range(9)]
    s, y, tau, c0, c1, c2, c3, c4, z = w
    ku, kv, p = [], [], s
    for q in range(128):
        ku.append(p); p = mul(p, s); kv.append(p); p = mul(p, s)
    P = pairs(msg)
    m = max(1, (len(P) + 127) // 128)
    b = []
    for t in range(m):
        acc = 0
        for q, (u, v) in enumerate(P[128 * t:128 * t + 128]):
            acc ^= clmul(u ^ ku[q], v ^ kv[q])
        b.append(mod(acc))
    ypow = [1]
    for e in range(8):
        ypow.append(mul(ypow[-1], y))
    mr = (m + 7) // 8
    V = len(msg)
    for r in range(mr):
        a = b[8 * r:8 * r + 8]; q = len(a); h = (q + 1) // 2; f = q // 2
        c = 0
        for i in range(1, f + 1):
            c ^= mul(a[i - 1] ^ ypow[2 * i - 1], a[h + i - 1] ^ ypow[2 * i])
        if q & 1:
            c ^= a[h - 1]
        V = mul(V, z) ^ c
    v = (V + tau) % (1 << 512)
    v2 = mul(v, v)
    H = mul(v ^ c2, mul(v2 ^ c0, v ^ v2 ^ c1) ^ c3) ^ c4
    return H.to_bytes(64, 'little')

if __name__ == '__main__':
    key = bytes((i * 73 + 11) & 255 for i in range(576))
    B = 16384
    lens = [0, 1, 7, 8, 9, 63, 64, 65, 127, 128, 129, 255, 1000, 1023, 1024, 1025, 1151, 2048, 4095, 4096,
            8191, 8192, 16383, B, B + 1, B + 1024, B + 1025]
    for r in range(3):
        for q in range(1, 9):
            base = r * 8 * B + (q - 1) * B
            lens += [base + 1, base + B]
    lens = sorted(set(lens))
    maxn = max(lens)
    msg = bytes((i * 137 + 29) & 255 for i in range(maxn))
    with (open(sys.argv[1], 'w') if len(sys.argv) > 1 else sys.stdout) as fo:
        fo.write('# ChainHash-512/PH v1 vectors from pyref512.py. key[i]=(i*73+11)&255 (576 B), msg[i]=(i*137+29)&255\n')
        for n in lens:
            fo.write('%d %s\n' % (n, hash512(key, msg[:n]).hex()))
    print('wrote %d vectors, max length %d' % (len(lens), maxn), file=sys.stderr)
