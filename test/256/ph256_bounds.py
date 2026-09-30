#!/usr/bin/env python3
"""ph256_bounds_v2.py -- collision numerators and score for PH-256 v2 (POWER KEY, two-level outer), against v1/v1.1
(independent masks).  q = 2^256, messages of at most 8L bytes, PH_M = 64 (4 KiB blocks).
Level 1 (v2): for a differing block the difference is a nonzero polynomial in s with exponent classes x_i -> s^2i,
y_i -> s^(2i-1) (cert_pk.c); its degree is <= 2P, P = the highest 1-based pair index holding message data.
If only x words carry data (n <= 256 bytes: the first tail group's x rows), dc = sum dx_i s^2i = (sum sqrt(dx_i) s^i)^2,
so the root count is <= P.  Hence
   d(n) = P(n)            if n <= 256
          2 P(n)           if 256 < n < BLOCK        (P(n) = 8 G + min(8, ceil((n - 512 G)/8)), G = (n-1) // 512)
          2 * 64 = 128     if n >= BLOCK            (a full block may differ)
v1/v1.1 (independent masks): d = 1 for every length.
Outer (both): N(L) = max(d(8L) + E(m) + 1, m' + 1), E(m) = 8 + m' - 1 if m' >= 2 else 2 floor(m/2) (SPEC-2L section 3).
Claim: N_v2(L) <= 2L for every L, equality only at L = 1 -> score 255.  Checked for every L < 2^20 and in +-3000
windows around 2^20..2^61."""
import math
M=64; BLOCK=64*M; REGION=8*BLOCK
def nblocks(n):
    if n==0: return 1
    return (n//REGION)*8 + -(-(n%REGION)//BLOCK)
def P(n):
    G=(n-1)//512; return 8*G+min(8,-(-(n-512*G)//8))
def d_v2(n):
    if n<=0: return 1
    if n>=BLOCK: return 2*M
    return P(n) if n<=256 else 2*P(n)
def N(L,d):
    n=8*L; m=nblocks(n); mp=-(-m//8); E= 8+mp-1 if mp>=2 else 2*(m//2)
    return max(d+E+1, mp+1)
def Ls():
    for L in range(1,1<<20): yield L
    for k in range(20,62):
        for L in range(max(1,(1<<k)-3000),(1<<k)+3000): yield L
worst=(1e9,None); bad=0; eq=[]
for L in Ls():
    n2=N(L,d_v2(8*L)); s=math.log2(L)+256-math.log2(n2)
    if s<worst[0]: worst=(s,L)
    if n2>2*L: bad+=1
    if n2==2*L: eq.append(L)
print("PH-256 v2 (power key), 4 KiB blocks: score %.2f (at L=%d); cases N > 2L: %d; N = 2L at L in %s"%(worst[0],worst[1],bad,eq[:5]))
print("\nNumerators N(L) (divide by 2^256): v1/v1.1 independent masks vs v2 power key")
print("| message limit | v1.1 (d = 1) | **v2** (d(L)) | d(L) |"); print("|---|---:|---:|---:|")
for name,nb in (("8 B",8),("16 B",16),("64 B",64),("256 B",256),("264 B",264),("512 B",512),("1 KiB",1024),("4 KiB",4096),("16 KiB",16384),("64 KiB",65536),("1 MiB",1<<20),("1 GiB",1<<30)):
    L=nb//8; print("| %s | %d | %d | %d |"%(name,N(L,1),N(L,d_v2(nb)),d_v2(nb)))
ok= bad==0 and eq==[1]
print("\nCLAIM (N_v2 <= 2L for all tested L, equality only at L = 1; score 255):", "PASS" if ok else "FAIL")
