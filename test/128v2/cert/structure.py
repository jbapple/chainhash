"""Generic CH-128/P structure over F = GF(2^n): pencil, L, g, A, E, Pi, CRT."""
import random
from gfp import GF, Ring, ptrim, padd, pmul, pmod, pdivmod, pgcd, peval, pderiv


class Pencil:
    def __init__(self, F, alphas, xi):
        self.F = F
        self.al = list(alphas)
        self.N = len(alphas)
        assert len(set(self.al)) == self.N
        self.xi = xi
        L = [1]
        for a in self.al:
            L = pmul(F, L, [a, 1])
        self.L = L
        self.c = L[1]                               # coefficient of T
        cc = 1
        for a in self.al:
            if a:
                cc = F.mul(cc, a)
        assert cc == self.c, "c != prod of nonzero alphas"
        self.g = padd(L, [F.mul(xi, self.c)])
        self.A = Ring(F, self.g)
        # evaluation matrix E (rows = alpha_i), and its inverse (Lagrange basis)
        self.lag = []
        for i, ai in enumerate(self.al):
            num = [1]
            den = 1
            for j, aj in enumerate(self.al):
                if j != i:
                    num = pmul(F, num, [aj, 1])
                    den = F.mul(den, ai ^ aj)
            il = F.inv(den)
            ell = [F.mul(x, il) for x in num]
            self.lag.append(ell + [0] * (self.N - len(ell)))

    # --- the literal pencil block (evaluation coordinates) ---
    def Q(self, u, v):
        F = self.F
        q0 = q1 = 0
        su = sv = 0
        for i in range(self.N):
            p = F.mul(u[i], v[i])
            q0 ^= p
            q1 ^= F.mul(self.al[i], p)
            su ^= u[i]
            sv ^= v[i]
        q1 ^= F.mul(self.xi, F.mul(su, sv))
        return (q0, q1)

    def E(self, U):
        return [peval(self.F, ptrim(U), a) for a in self.al]

    def Einv(self, u):
        r = [0] * self.N
        for i in range(self.N):
            if u[i]:
                for k in range(self.N):
                    r[k] ^= self.F.mul(u[i], self.lag[i][k])
        return r

    def Pi(self, W):
        F = self.F
        return (F.mul(self.c, W[self.N - 1]), F.mul(self.c, W[self.N - 2]))

    # --- keyed block: p positions, masks E(s^(2m+1)), E(s^(2m+2)) ---
    def block_literal(self, s, msgs):
        """msgs = list of (m, m') word vectors (evaluation coords). Returns summed (Q0,Q1)."""
        A = self.A
        q0 = q1 = 0
        spow = A.norm([1])
        pw = []
        for j in range(1, 2 * len(msgs) + 1):
            spow = A.mul(spow, s)
            pw.append(spow)
        for m, (mu, mv) in enumerate(msgs):
            k = self.E(pw[2 * m])
            kp = self.E(pw[2 * m + 1])
            u = [a ^ b for a, b in zip(mu, k)]
            v = [a ^ b for a, b in zip(mv, kp)]
            a0, a1 = self.Q(u, v)
            q0 ^= a0
            q1 ^= a1
        return (q0, q1)

    def block_algebraic(self, s, msgs):
        A = self.A
        W = A.zero()
        for m, (mu, mv) in enumerate(msgs):
            U = A.add(self.Einv(mu), A.pow(s, 2 * m + 1))
            V = A.add(self.Einv(mv), A.pow(s, 2 * m + 2))
            W = A.add(W, A.mul(U, V))
        return self.Pi(W)

    # --- factorization of g into irreducibles (equal-degree 2, char 2 trace split) ---
    def split_quadratics(self, rng=random):
        F, A = self.F, self.A
        n = F.n
        todo = [ptrim(self.g)]
        done = []
        while todo:
            f = todo.pop()
            if len(f) - 1 == 2:
                done.append(f)
                continue
            R = Ring(F, f)
            while True:
                r = [rng.getrandbits(n) for _ in range(len(f) - 1)]
                t = R.zero()
                x = list(r)
                for _ in range(2 * n):
                    t = R.add(t, x)
                    x = R.sq(x)
                h = pgcd(F, f, ptrim(t))
                if 0 < len(h) - 1 < len(f) - 1:
                    q, rem = pdivmod(F, f, h)
                    assert not rem
                    todo += [h, q]
                    break
        return done


def gf2_rank(vecs):
    basis = []
    for v in vecs:
        for b in basis:
            v = min(v, v ^ b)
        if v:
            basis.append(v)
            basis.sort(reverse=True)
    return len(basis)
