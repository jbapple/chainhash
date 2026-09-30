"""GF(2^64) with modulus X^64 + X^4 + X^3 + X + 1  (low reduction word 0x1b).
Rebuilt independently. add=XOR, mul=clmul then reduce. Also GF(2^128) helper."""

MOD64_LOW = 0x1b            # X^4+X^3+X+1 ; full modulus = (1<<64)|0x1b
POLY64 = (1 << 64) | MOD64_LOW

def clmul(a, b):
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
    return r

def reduce_mod(p, deg, low):
    """reduce polynomial p modulo X^deg + low (low has degree < deg)."""
    full = (1 << deg) | low
    for i in range(p.bit_length() - 1, deg - 1, -1):
        if (p >> i) & 1:
            p ^= full << (i - deg)
    return p

def reduce64(p):
    return reduce_mod(p, 64, MOD64_LOW)

def mul64(a, b):
    return reduce64(clmul(a, b))

def xtimes(a):
    a <<= 1
    if a & (1 << 64):
        a ^= POLY64            # equivalently (a ^ (1<<64)) ^ 0x1b, then it's <64 bits
    return a & ((1 << 64) - 1)

# sanity: xtimes via mul by X (X = 0b10 = 2)
def _selftest():
    import random
    X = 2
    for _ in range(2000):
        a = random.getrandbits(64)
        assert xtimes(a) == mul64(a, X), "xtimes != mul by X"
    # field axioms spot check: distributivity, associativity
    for _ in range(2000):
        a = random.getrandbits(64); b = random.getrandbits(64); c = random.getrandbits(64)
        assert mul64(a, b ^ c) == (mul64(a,b) ^ mul64(a,c))
        assert mul64(mul64(a,b),c) == mul64(a, mul64(b,c))
        assert mul64(a,b) == mul64(b,a)
    # X^128 == 0x145 claim: X^64 = 0x1b, so X^128 = 0x1b^2 (carryless) reduced? Actually
    # X^128 mod (X^64+0x1b) as a 64-bit element:
    x128 = reduce_mod(1 << 128, 64, MOD64_LOW)
    assert x128 == mul64(reduce_mod(1<<64,64,MOD64_LOW), reduce_mod(1<<64,64,MOD64_LOW))
    print("gf selftest OK; X^128 mod poly =", hex(x128),
          "; clmul(0x1b,0x1b)=", hex(clmul(0x1b,0x1b)))

if __name__ == "__main__":
    _selftest()
