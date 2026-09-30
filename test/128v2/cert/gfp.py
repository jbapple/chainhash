"""Generic GF(2^n) + polynomials over it + the quotient ring A = F[T]/g.
Written independently for the CH-128/P check. Field elements are ints (bit i = X^i).
Polynomials over F are lists of ints, index = degree (little-endian)."""
import random


def clmul(a, b):
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
    return r


class GF:
    def __init__(self, n, modlow):
        self.n = n
        self.q = 1 << n
        self.mod = (1 << n) | modlow
        self.mask = self.q - 1
        self.table = None
        if n <= 8:
            self.table = [[self._mul(a, b) for b in range(self.q)] for a in range(self.q)]

    def _mul(self, a, b):
        p = clmul(a, b)
        n = self.n
        for i in range(p.bit_length() - 1, n - 1, -1):
            if (p >> i) & 1:
                p ^= self.mod << (i - n)
        return p

    def mul(self, a, b):
        if self.table is not None:
            return self.table[a][b]
        return self._mul(a, b)

    def pow(self, a, e):
        r = 1
        while e:
            if e & 1:
                r = self.mul(r, a)
            a = self.mul(a, a)
            e >>= 1
        return r

    def inv(self, a):
        assert a != 0
        return self.pow(a, self.q - 2)

    def is_irreducible_modulus(self):
        # Rabin: X^(2^n) == X mod m, and gcd(X^(2^(n/r)) - X, m) = 1 for prime r | n
        n = self.n
        def xpow2k(k):
            x = 2
            for _ in range(k):
                x = self._mul(x, x)
            return x
        if xpow2k(n) != 2:
            return False
        # gcd over GF(2)[X]
        def pgcd(a, b):
            while b:
                while a and a.bit_length() >= b.bit_length():
                    a ^= b << (a.bit_length() - b.bit_length())
                a, b = b, a
            return a
        for r in {r for r in range(2, n + 1) if n % r == 0 and all(r % d for d in range(2, r))}:
            h = xpow2k(n // r) ^ 2
            if pgcd(self.mod, h) != 1:
                return False
        return True


# ---------------- polynomials over F ----------------
def ptrim(a):
    a = list(a)
    while a and a[-1] == 0:
        a.pop()
    return a


def padd(a, b):
    n = max(len(a), len(b))
    return ptrim([(a[i] if i < len(a) else 0) ^ (b[i] if i < len(b) else 0) for i in range(n)])


def pmul(F, a, b):
    if not a or not b:
        return []
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x == 0:
            continue
        for j, y in enumerate(b):
            if y:
                r[i + j] ^= F.mul(x, y)
    return ptrim(r)


def pdivmod(F, a, b):
    a = ptrim(a)
    b = ptrim(b)
    assert b
    inv_lead = F.inv(b[-1])
    q = [0] * max(0, len(a) - len(b) + 1)
    a = list(a)
    while len(a) >= len(b) and a:
        coef = F.mul(a[-1], inv_lead)
        shift = len(a) - len(b)
        q[shift] = coef
        for i, y in enumerate(b):
            a[i + shift] ^= F.mul(coef, y)
        a = ptrim(a)
    return ptrim(q), a


def pmod(F, a, b):
    return pdivmod(F, a, b)[1]


def pgcd(F, a, b):
    a, b = ptrim(a), ptrim(b)
    while b:
        a, b = b, pmod(F, a, b)
    if a:
        il = F.inv(a[-1])
        a = [F.mul(x, il) for x in a]
    return a


def peval(F, a, x):
    r = 0
    for c in reversed(a):
        r = F.mul(r, x) ^ c
    return r


def pderiv(F, a):
    # char 2: d/dT T^i = i T^(i-1): only odd i survive
    return ptrim([a[i] if i % 2 == 1 else 0 for i in range(1, len(a))])


class Ring:
    """A = F[T]/g, elements as length-deg(g) coefficient lists."""
    def __init__(self, F, g):
        self.F = F
        self.g = ptrim(g)
        self.d = len(self.g) - 1
        assert self.g[-1] == 1

    def norm(self, a):
        a = pmod(self.F, a, self.g)
        return a + [0] * (self.d - len(a))

    def mul(self, a, b):
        return self.norm(pmul(self.F, ptrim(a), ptrim(b)))

    def add(self, a, b):
        return [x ^ y for x, y in zip(a, b)]

    def pow(self, a, e):
        r = self.norm([1])
        while e:
            if e & 1:
                r = self.mul(r, a)
            a = self.mul(a, a)
            e >>= 1
        return r

    def sq(self, a):
        return self.mul(a, a)

    def frob(self, a, k):
        for _ in range(k):
            a = self.sq(a)
        return a

    def rand(self):
        return [random.getrandbits(self.F.n) for _ in range(self.d)]

    def zero(self):
        return [0] * self.d
