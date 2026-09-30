#!/usr/bin/env python3
"""ph_oracle.py DUMP -- check C ph_mul products (cert_pk -DDUMP) against an independent Python GF(2^256) multiply."""
import sys
F=(1<<256)|(1<<10)|(1<<5)|(1<<2)|1
def mul(a,b):
    r=0
    while b:
        if b&1: r^=a
        a<<=1; b>>=1
    for i in range(511,255,-1):
        if (r>>i)&1: r^=F<<(i-256)
    return r
n=bad=0
for line in open(sys.argv[1]):
    a,b,c=(int(t,16) for t in line.split()); n+=1; bad+= mul(a,b)!=c
print("ph_mul vs Python oracle: %d products, %d mismatches -> %s"%(n,bad,"PASS" if n and not bad else "FAIL"))
