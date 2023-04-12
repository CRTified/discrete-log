from sage.rings.finite_rings.finite_field_constructor import FiniteField, GF
from sage.rings.integer import Integer

from util import countCalls, timeall

class StrongOracle():
    def __init__(self, P, k):
        self.F = GF(P.order())
        self._base = P
        self.public = k * P
        self._secret = k
        self._cache = {}
        self.__addToCache(1, self._base)
        self.__addToCache(k, self.public)
        self.__addToCache(self._base.order(), self._base.curve()(0))
        self._calls = 0
        print(self.F)

    def __addToCache(self, x, P):
        if P.is_zero():
            self._cache[(0, 0)] = self._base.order()
            return
        self._cache[P.xy()] = self.F(x)

    def __getFromCache(self, Q):
        if Q.is_zero():
            return self._base.order()
        if res := self._cache.get(Q.xy(), False):
            return Integer(res)
        print("Cache miss, computing naively: discrete_log({}, {})".format(Q, self._base))
        print("This might take a *long* time (or is not feasible)")
        res = discrete_log(Q, self._base, operation="+")
        self.__addToCache(res, Q)
        return Integer(res)

    def dlogSolution(self, Q):
        """
          Oracle call!

          Returns dlog_P(Q) = q for Q = q*P.
          Uses `self._cache` to do this efficiently, and `discrete_log` for cache
          misses.
        """
        return self.__getFromCache(Q)

    def createScalarMultiple(self, v):
        """
          Not an oracle call!
          Can be computed with public knowledge.

          Returns a point vP.
          Side effect: Caches dlog(vP) = v
        """
        vP = Integer(v) * self._base
        self.__addToCache(v, vP)
        return vP

    def add(self, Q, R):
        """
          Not an oracle call!
          Can be computed with public knowledge.

          Returns Q + R = qP + rP = (q+r) P
          Side effect: Caches dlog( (q+r)P ) = (q+r)
        """
        q = self.dlogSolution(Q)
        r = self.dlogSolution(R)
        self.__addToCache(q + r, Q+R)
        return Q + R

    def mult(self, v, Q):
        """
          Not an oracle call!
          Can be computed with public (adaptive) knowledge.

          Returns vQ = (vq) P
          Side effect: Caches dlog( vqP ) = vq
        """
        v = Integer(v)
        q = self.dlogSolution(Q)
        self.__addToCache(v * q, v * Q)
        return v * Q

    @timeall
    def oracle(self, Q, R):
        """
          Oracle call!

          Computes DH(qP, rP) = qrP.
          Side effect: Caches dlog(qrP) = qr
        """
        self._calls = self._calls + 1
        q = self.dlogSolution(Q)
        r = self.dlogSolution(R)
        res = self.F(q) * R
        self.__addToCache(q * r, res)
        return res
