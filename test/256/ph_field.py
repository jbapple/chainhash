#!/usr/bin/env python3
"""ph_field.py -- field certificate for PH-256: f = x^256 + x^10 + x^5 + x^2 + 1 is irreducible over GF(2).
Rabin: deg f = 256 = 2^8, the only prime divisor of 256 is 2, so f is irreducible iff
  (i)  x^(2^256) == x  (mod f)   and   (ii) gcd(x^(2^128) - x, f) == 1.
Polynomials are Python ints (bit i = coefficient of x^i)."""
N=256; F=(1<<256)|(1<<10)|(1<<5)|(1<<2)|1
def pmod(a,m):
    dm=m.bit_length()-1
    while a.bit_length()-1>=dm: a^=m<<(a.bit_length()-1-dm)
    return a
def pmul(a,b):
    r=0
    while b:
        if b&1: r^=a
        a<<=1; b>>=1
    return r
def sqmod(a): 
    # squaring in char 2: spread bits
    r=0; i=0
    while a:
        if a&1: r|=1<<(2*i)
        a>>=1; i+=1
    return pmod(r,F)
def pgcd(a,b):
    while b: a,b=b,pmod(a,b)
    return a
x=2; t=x; x2_128=None
for k in range(1,257):
    t=sqmod(t)
    if k==128: x2_128=t
ok1 = (t==x)
g = pgcd(F, x2_128^x)
ok2 = (g==1)
print("f = x^256 + x^10 + x^5 + x^2 + 1  (low word 0x425)")
print("(i)  x^(2^256) mod f == x :", ok1)
print("(ii) gcd(x^(2^128) - x, f) == 1 :", ok2, "(gcd degree %d)"%(g.bit_length()-1))
print("FIELD CERT", "PASS: f irreducible, GF(2)[x]/f = GF(2^256)" if ok1 and ok2 else "FAIL")
# negative control: a reducible pentanomial must fail
Fbad=(1<<256)|(1<<4)|(1<<3)|(1<<1)|1  # x^256+x^4+x^3+x+1 (check whatever result, report)
def rabin256(Fp):
    global F
    F0=F; F=Fp; t=2; h=None
    for k in range(1,257):
        t=sqmod(t)
        if k==128: h=t
    r=(t==2, pgcd(Fp,h^2)==1); F=F0; return r
G=(1<<128)|0x87
def compose_x_plus_1(p):   # p(x+1)
    r=0; xp1=1; d=0; q=p
    powv=1; res=0
    for i in range(p.bit_length()):
        if (p>>i)&1: res^=powv
        powv=pmul(powv,3)
    return res
G1=compose_x_plus_1(G)
for name,Fp in [("x^256+1",(1<<256)|1),("GCM(x)*GCM(x+1) (two degree-128 irreducibles)",pmul(G,G1))]:
    a,b=rabin256(Fp); print("negative control %-70s (i)=%s (ii)=%s -> %s"%(name,a,b,"rejected (OK)" if not(a and b) else "ACCEPTED (tool broken)"))
