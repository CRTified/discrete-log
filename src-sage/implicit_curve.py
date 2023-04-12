from abc import ABC, abstractmethod

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.finite_rings.finite_field_constructor import FiniteField, GF
from sage.schemes.elliptic_curves.constructor import EllipticCurve
from sage.rings.integer import Integer

class ImplicitCurve(ABC):
    def __init__(self, P, oracle):
        self._P = P
        self._o = oracle

    def Fq_add(self, aP, bP):
        """
        Given aP and bP
        Returns (a+b)P
        """
        return self._o.add(aP, bP)

    def Fq_sub(self, aP, bP):
        """
        Given aP and bP
        Return (a-b)P
        """
        return self._o.add(aP, self.Fq_smult(-1, bP))

    def Fq_smult(self, s, aP):
        """
        Given s and aP
        Return saP
        """
        return self._o.mult(s, aP)

    def Fq_pmult(self, aP, bP):
        """
        Given aP and bP
        Return abP
        """
        return self._o.oracle(aP, bP)

    def Fq_psquare(self, aP):
        """
        Given aP
        Return a^2 P
        """
        return self._o.oracle(aP, aP)

    def Fq_pow(self, aP, e):
        """
        Generic scalar multiplication using
        double-and-add.
        """
        tmp = aP
        for ei in Integer(e).bits()[-2::-1]:
            tmp = self.Fq_psquare(tmp)
            if ei:
                tmp = self.Fq_pmult(tmp, aP)
        return tmp

    def Fq_inverse(self, aP):
        return self.Fq_pow(aP, self._P.order() - 2)

    @abstractmethod
    def compute_rhs(self, XP):
        pass

    @abstractmethod
    def is_valid_x(self, XP):
        pass
    @abstractmethod
    def is_valid_y(self, YP):
        pass

    @abstractmethod
    def lift(self, P):
        pass

    def is_square(self, XP):
        """
        Compute the implicit Legendre symbol and check if it is 1
        """
        b_r1half = self.Fq_pow(XP, (self._P.order() - 1) / 2)
        return b_r1half == self._P

    def scalarMult(self, aP_I, b):
        """
        Generic scalar multiplication using
        double-and-add.
        """
        tmp = aP_I
        for b in b.bits()[-2::-1]:
            tmp = self.pointDouble(tmp)
            if b:
                tmp = self.pointAdd(tmp, aP_I)
        return tmp

    #@abstractmethod
    #def pointAdd(self, aP_I, bP_I):
    #    pass

class ImplicitWeierstrass(ImplicitCurve):
    def __init__(self, P, oracle, a, b):
        super().__init__(P, oracle)
        self.a = a
        self.b = b
        self.E_e = EllipticCurve(GF(P.order()), [a, b])
        nullP = self._o.createScalarMultiple(0)
        self.impl_zero = nullP, nullP

    def hashP(self, aP_I):
        XP, YP = aP_I
        if XP.is_zero():
            return 0
        x, y = XP.xy()
        return int(x) << 1 | int(Integer(y).bits()[0])

    def compute_rhs(self, XP):
        xCubeP = self.Fq_pmult(XP, self.Fq_psquare(XP))
        AxP = self.Fq_smult(self.a, XP)
        BP = self._o.createScalarMultiple(self.b)
        return self.Fq_add(xCubeP, self.Fq_add(AxP, BP))

    def compute_lhs(self, YP):
        return self.Fq_psquare(YP)


    def is_valid_x(self, XP):
        """
        This can be checked by implicitly evaluating the Legendre symbol of the rhs
        """
        tmp = self.compute_rhs(XP)
        return self.is_square(XP)

    def is_valid_y(self, YP):
        """
        lhs is valid if it is a square
        """
        return self.is_square(YP)

    def lift(self, P):
        if P.is_zero():
            return self.impl_zero
        X, Y = P.xy()
        return (self._o.createScalarMultiple(X),
                self._o.createScalarMultiple(Y))

    def pointDouble(self, aP_I):
        # Doubling: 2M + 1.5log(q) M
        # Both: 1M + 1M
        X1, Y1 = aP_I

        # # λ = (3x^2 + A) / 2y
        y1_2_inv = self.Fq_pow(self.Fq_smult(2, Y1),
                               self._P.order() - 2)
        lambd = self.Fq_pmult(self.Fq_add(self.Fq_smult(3,
                                                        self.Fq_psquare(X1)),
                                          self._o.createScalarMultiple(self.a)
                                          ),
                               y1_2_inv)
        X3 = self.Fq_sub(self.Fq_sub(self.Fq_psquare(lambd),
                                     X1),
                         X1)
        Y3 = self.Fq_sub(self.Fq_pmult(lambd,
                                       self.Fq_sub(X1,
                                                   X3)),
                         Y1)
        return X3, Y3

    def pointAdd(self, aP_I, bP_I):
        # Addition: 1M + 1.5log(q) M
        # Both: 1M + 1M
        X1, Y1 = aP_I
        X2, Y2 = bP_I

        if X1 == X2:
            if Y1 != Y2:
                return self.impl_zero
            # Perform doubling instead
            return self.pointDouble(aP_I)

        # Addition, as -op1 != op2
        # λ = (y2 - y1) / (x2 - x1)
        lambd = self.Fq_pmult(self.Fq_sub(Y2, Y1),
                              self.Fq_pow(self.Fq_sub(X2, X1),
                                          self._P.order() - 2))
        X3 = self.Fq_sub(self.Fq_sub(self.Fq_psquare(lambd),
                                     X1),
                         X2)
        Y3 = self.Fq_sub(self.Fq_pmult(lambd,
                                       self.Fq_sub(X1,
                                                   X3)),
                         Y1)
        return X3, Y3

class ImplicitProjWeierstrass(ImplicitCurve):
    def __init__(self, P, oracle, a, b):
        super().__init__(P, oracle)
        self.a = a
        self.b = b
        self.E_e = EllipticCurve(GF(P.order()), [a, b])
        nullP = self._o.createScalarMultiple(0)
        self.impl_zero = self.lift(self.E_e(0))

    def hashP(self, aP_I):
        XP, YP, ZP = self.scale(aP_I)
        if XP.is_zero():
            return 0
        x, y = XP.xy()
        return int(x) << 1 | int(Integer(y).bits()[0])

    def compute_rhs(self, XP):
        xCubeP = self.Fq_pmult(XP, self.Fq_psquare(XP))
        AxP = self.Fq_smult(self.a, XP)
        BP = self._o.createScalarMultiple(self.b)
        return self.Fq_add(xCubeP, self.Fq_add(AxP, BP))

    def compute_lhs(self, YP):
        return self.Fq_psquare(YP)

    def is_valid_x(self, XP):
        """
        This can be checked by implicitly evaluating the Legendre symbol of the rhs
        """
        tmp = self.compute_rhs(XP)
        return self.is_square(XP)

    def is_valid_y(self, YP):
        """
        lhs is valid if it is a square
        """
        return self.is_square(YP)

    def lift(self, P_E):
        X, Y, Z = P_E
        return (self._o.createScalarMultiple(X),
                self._o.createScalarMultiple(Y),
                self._o.createScalarMultiple(Z))

    def scale(self, aP_I):
        X1, Y1, Z1 = aP_I
        if Z1 == self._P:
            return X1, Y1, Z1
        A = self.Fq_inverse(Z1)
        X3 = self.Fq_pmult(X1, A)
        Y3 = self.Fq_pmult(Y1, A)
        Z3 = self._P
        return X3, Y3, Z3

    def pointInv(self, aP_I):
        X1, Y1, Z1 = aP_I
        Y2 = self.Fq_smult(-1, Y1)
        return X1, Y2, Z1

    def pointAdd(self, aP_I, bP_I):
        """
        add-1998-cmo-2
        12M + 2S + 6add + 1*2
        """

        X1, Y1, Z1 = aP_I
        X2, Y2, Z2 = bP_I

        if Z1.is_zero():
            return bP_I
        if Z2.is_zero():
            return aP_I

        if X1 == X2:
            if Y1 != Y2:
                return self.impl_zero
            # Perform doubling instead
            return self.pointDouble(aP_I)

        Y1Z2 = self.Fq_pmult(Y1, Z2)    # O
        X1Z2 = self.Fq_pmult(X1, Z2)    # O
        Z1Z2 = self.Fq_pmult(Z1, Z2)    # O
        t0   = self.Fq_pmult(Y2, Z1)    # O
        u    = self.Fq_sub(t0, Y1Z2)
        uu   = self.Fq_psquare(u)       # O
        t1   = self.Fq_pmult(X2, Z1)    # O
        v    = self.Fq_sub(t1, X1Z2)
        vv   = self.Fq_psquare(v)       # O
        vvv  = self.Fq_pmult(v, vv)     # O
        R    = self.Fq_pmult(vv, X1Z2)  # O
        t2   = self.Fq_smult(2, R)
        t3   = self.Fq_pmult(uu, Z1Z2)  # O
        t4   = self.Fq_sub(t3, vvv)
        A    = self.Fq_sub(t4, t2)
        X3   = self.Fq_pmult(v, A)      # O
        t5   = self.Fq_sub(R, A)
        t6   = self.Fq_pmult(vvv, Y1Z2) # O
        t7   = self.Fq_pmult(u, t5)     # O
        Y3   = self.Fq_sub(t7, t6)
        Z3   = self.Fq_pmult(vvv, Z1Z2) # O

        return X3, Y3, Z3

    def pointDouble(self, aP_I):
        """
        dbl-2007-bl
        5M + 6S + 1*a + 7add + 3*2 + 1*3
        """
        X1, Y1, Z1 = aP_I

        XX  = self.Fq_psquare(X1)   # O
        ZZ  = self.Fq_psquare(Z1)   # O
        t0  = self.Fq_smult(3, XX)
        t1  = self.Fq_smult(self.a, ZZ)
        w   = self.Fq_add(t1, t0)
        t2  = self.Fq_pmult(Y1, Z1) # O
        s   = self.Fq_smult(2, t2)
        ss  = self.Fq_psquare(s)    # O
        sss = self.Fq_pmult(s, ss)  # O
        R   = self.Fq_pmult(Y1, s)  # O
        RR  = self.Fq_psquare(R)    # O
        t3  = self.Fq_add(X1, R)
        t4  = self.Fq_psquare(t3)   # O
        t5  = self.Fq_sub(t4, XX)
        B   = self.Fq_sub(t5, RR)
        t6  = self.Fq_psquare(w)    # O
        t7  = self.Fq_smult(2, B)
        h   = self.Fq_sub(t6, t7)
        X3  = self.Fq_pmult(h, s)   # O
        t8  = self.Fq_sub(B, h)
        t9  = self.Fq_smult(2, RR)
        t10 = self.Fq_pmult(w, t8)  # O
        Y3  = self.Fq_sub(t10, t9)
        Z3  = sss

        return X3, Y3, Z3
