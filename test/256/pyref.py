#!/usr/bin/env python3
"""ph256_pyref_v2.py VECTORS [maxlen] -- an independent Python implementation of PH-256 v2 (POWER KEY) written from SPEC.md
(field, key expansion from a seed, region/tail layout, two-level outer stage, finalizer); checks every vector line
"PH_M key n digest" (key = a decimal seed, or R = raw key bytes (73 i + 11) mod 256) with n <= maxlen.  Shares no code with the C reference."""
import sys
M64=(1<<64)-1
F=(1<<256)|(1<<10)|(1<<5)|(1<<2)|1
def mul(a,b):
    r=0
    while b:
        if b&1: r^=a
        a<<=1; b>>=1
    for i in range(r.bit_length()-1,255,-1):
        if (r>>i)&1: r^=F<<(i-256)
    return r
def splitmix(st):
    st[0]=(st[0]+0x9E3779B97F4A7C15)&M64; z=st[0]
    z=((z^(z>>30))*0xBF58476D1CE4E5B9)&M64; z=((z^(z>>27))*0x94D049BB133111EB)&M64
    return z^(z>>31)
def el(st): return sum(splitmix(st)<<(64*i) for i in range(4))
def key(seed,m):
    if seed=='R':                                              # raw key bytes: s | y | z | t | c0..c4, each 32 B LE
        b=bytes(((73*i+11)&255) for i in range(288)); el9=[int.from_bytes(b[32*i:32*i+32],'little') for i in range(9)]
    else:
        st=[int(seed)^0x5048323536763200]; el9=[el(st) for _ in range(9)]
        if el9[1]==0: el9[1]=1
    s,Y,z,tau=el9[0],el9[1],el9[2],el9[3]; c=el9[4:9]
    k=[];l=[]; pw=s
    for j in range(m): k.append(pw); pw=mul(pw,s); l.append(pw); pw=mul(pw,s)      # k_j = s^(2j+1), l_j = s^(2j+2)
    return k,l,Y,c,tau,z
def w64(b,o): return int.from_bytes(b[o:o+8],'little')
def elem(b,o,stride): return sum(w64(b,o+stride*l)<<(64*l) for l in range(4))
def hash_(seed,m,msg):
    k,l,Y,c,tau,z=key(seed,m); n=len(msg); BLOCK=64*m; REGION=8*BLOCK
    def block(xs,ys):
        acc=0
        for i in range(m): acc^=mul(xs[i]^k[i],ys[i]^l[i])
        return acc
    bv=[]; nfull=n//REGION                                   # level 1: block values in message order
    for r in range(nfull):
        R=msg[r*REGION:(r+1)*REGION]
        for b in range(8):
            xs=[elem(R,512*i+8*b,64) for i in range(m)]; ys=[elem(R,512*i+256+8*b,64) for i in range(m)]
            bv.append(block(xs,ys))
    T=msg[nfull*REGION:]; rem=len(T); b=0
    while b*BLOCK<rem:
        blk=T[b*BLOCK:(b+1)*BLOCK]; blk=blk+bytes(BLOCK-len(blk))
        xs=[elem(blk,512*(i//8)+8*(i%8),64) for i in range(m)]; ys=[elem(blk,512*(i//8)+256+8*(i%8),64) for i in range(m)]
        bv.append(block(xs,ys)); b+=1
    if n==0: bv=[0]
    yp=[1]
    for e in range(8): yp.append(mul(yp[-1],Y))
    def region(a):                                            # level 2: the block formula over K with key y
        q=len(a); h=(q+1)//2; f=q//2; cc=0
        for i in range(1,f+1): cc^=mul(a[i-1]^yp[2*i-1],a[h+i-1]^yp[2*i])
        if q%2: cc^=a[h-1]
        return cc
    v=n                                                       # Horner in z, length leading
    for r in range(0,len(bv),8): v=mul(v,z)^region(bv[r:r+8])
    X=(v+tau)&((1<<256)-1); G=mul(X,X); t=mul(G^c[0],X^G^c[1]); out=mul(X^c[2],t^c[3])^c[4]
    return out.to_bytes(32,'little').hex()
maxlen=int(sys.argv[2]) if len(sys.argv)>2 else 70000
n_ok=n_all=0
for line in open(sys.argv[1]):
    m,seed,n,d=line.split(); m=int(m); n=int(n)
    if n>maxlen: continue
    msg=bytes(((j*137+29)&255) for j in range(n)); n_all+=1; got=hash_(seed,m,msg)
    if got==d: n_ok+=1
    else: print("MISMATCH m=%d seed=%d n=%d"%(m,seed,n))
print("PH-256 v2 vectors vs independent Python implementation: %d/%d match (n <= %d) -> %s"%(n_ok,n_all,maxlen,"PASS" if n_all and n_ok==n_all else "FAIL"))
sys.exit(0 if n_all and n_ok==n_all else 1)
