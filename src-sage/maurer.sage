import itertools
from operator import itemgetter
from collections import defaultdict
from os.path import exists
import os
import errno
import binascii
import sys

from curve_list import CURVES, getChallenge, getCurveChallenge, pickImplCurve
from oracles import StrongOracle
from util import countCalls, repeatUntilNotNone, ppDict, timeit
from numeric_util import smoothness, genericSquareAndMultiply
from implicit_curve import ImplicitWeierstrass, ImplicitProjWeierstrass

@timeit
def maurer_dlog(E, P, Q, o, METHOD):
    scalarSquareMultiply = genericSquareAndMultiply(o.oracle)

    # Prime Order of generator point P
    q = P.order(); assert q.is_prime()
    F = GF(q)

    impl_param, factors, impl_gen, impl_type = pickImplCurve(NAME)
    N = prod(pi ^ ei for pi, ei in factors)
    print("Using auxiliary curve: ", impl_type, impl_param)
    print("With Generators ", impl_gen)
    print("Order: ", N, "=", factors)
    print("Smoothness: 2^", len(smoothness(factors).bits()))

    A, B = impl_param
    impl_curve = ImplicitProjWeierstrass(P, o, A, B)
    impl_doubleAndAdd = impl_curve.scalarMult


    while True:
        print(o._calls)
        # Build implicit representation of (k+α)P
        alpha = F.random_element()
        print(f"Trying α={alpha}")
        alphaP = o.createScalarMultiple(alpha)
        XP = o.add(Q, alphaP)

        # Split is_valid_x apart, to not recompute the rhs
        betaP = impl_curve.compute_rhs(XP)
        if impl_curve.is_square(betaP):
            break
        print(o._calls)


    print("α = ", alpha)
    print("βP=", betaP)

    def cipolla(betaP):
        print("B", o._calls)
        while impl_curve.is_square(asqn := (o.add(o.createScalarMultiple((a := F.random_element())^2), o.mult(-1, betaP)))):
            print(f"Trying a={a}")
            pass
        impl_a = o.createScalarMultiple(a)
        print(o._calls)

        def implicit_mult_ext(op1, op2):
            if op1 == op2:
                # Doubling is slightly cheaper, as we can reuse o.oracle(op1[0], op2[1])
                # This saves log_2(q) - 1 oracle calls
                new_left = o.add(o.oracle(op1[0], op2[0]), o.oracle(o.oracle(op1[1], op2[1]), asqn))
                v = o.oracle(op1[0], op2[1])
                new_right = o.add(v, v)
            else:
                new_left = o.add(o.oracle(op1[0], op2[0]), o.oracle(o.oracle(op1[1], op2[1]), asqn))
                new_right = o.add(o.oracle(op1[0], op2[1]), o.oracle(op1[1], op2[0]))
            return new_left, new_right

        tmp = (impl_a, P)
        for b in Integer((P.order() + 1) / 2).bits()[-2::-1]:
            # Square
            tmp = implicit_mult_ext(tmp, tmp)
            if b:
                # Multiply
                tmp = implicit_mult_ext(tmp, (impl_a, P))
        assert(tmp[1].is_zero())
        return tmp

    print("Calls before calculating sqrt(β)P:" ,o._calls)
    YP, _ = cipolla(betaP)
    print("Calls after calculating sqrt(β)P:" ,o._calls)
    print("Implicit representation of (X, Y):")
    print("(", XP, ",")
    print(" ", YP, ")")
    R = (XP, YP, P)

    Omega = impl_gen
    print("Generator for implicit curve: ", Omega)

    print("Calls before calculating list:" ,o._calls)

    @timeit
    def calc_subR():
        sub_R = {}
        min_f = min(pi ^ ei for pi, ei in factors)
        Rpow2i_lut = []
        tmp = R
        for i in range(len(Integer(N / min_f).bits())):
            Rpow2i_lut.append(tmp)
            tmp = impl_curve.pointDouble(tmp)

        print("Doubling", o._calls)
        for pi, ei in factors:
            print(f"Calculating N/(pi^ei) Q for pi^ei = {pi^ei}, so exp = {N / (pi^ei)}")
            tmp = reduce(impl_curve.pointAdd, itertools.compress(Rpow2i_lut, Integer(N / pi^ei).bits()))
            sub_R[pi^ei] = impl_curve.scale(tmp)
        return sub_R

    sub_R = calc_subR()

    print("Calls after calculating list:", o._calls)

    @timeit
    def solve_factorgroup(arg):
        """
        This solves u_i * P_aux_i = Q_aux_i by pure brute-force.
        Prime powers are slightly optimized as in the prime-power
        Pohlig-Hellman, i.e. we search p_i dlogs for every e_i.
        """
        pi, ei = arg

        u_i = 0
        print("p_i =", pi, ", e_i =", ei)
        for j in range(1, ei + 1):
            OmegaNull = Integer(N / (pi ^ j)) * Omega
            print("j = ", j)
            print("u_i = ", u_i)

            # This needs oracle calls!
            Rnull = impl_curve.pointAdd(sub_R[pi^j], impl_curve.lift(-u_i * OmegaNull))
            Rnull = impl_curve.scale(Rnull)

            if Rnull[2]:
                OmegaNull = Integer(N / (pi)) * Omega
                n = 1
                T = OmegaNull
                impl_T = impl_curve.lift(T)
                while Rnull[0] != impl_T[0]: # TODO: Inverse is also fine
                    n += 1
                    T += OmegaNull
                    impl_T = impl_curve.lift(T)
                else:
                    if Rnull[1] == -impl_T[1]:
                        print("Inverse!")
                        n = pi - n
                    u_i += n * pi ^ (j-1)
            print(pi, j, u_i)
        return (u_i, pi ^ ei)

    def solve_cheat(arg):
        """

        """
        explicit_curve = EllipticCurve(GF(P.order()), impl_param)
        pi, ei = arg
        p = pi^ei
        Omegap = Integer(N/p) * Omega
        R = sub_R[p]

        R_ex, R_ey = o.dlogSolution(R[0]), o.dlogSolution(R[1])
        R_E = explicit_curve(R_ex, R_ey)
        u = discrete_log(R_E, Omegap, ord=p, operation="+")
        return u, p

    @timeit
    def solve_externalcb(arg):
        """
        This solves u_i * P_aux_i = Q_aux_i with a pre-computed codebook.
        Looking up values is free in terms of oracle calls
        """
        pi, ei = arg
        p = pi^ei
        Omegap = Integer(N/p) * Omega
        print(f"Solving {p}")
        entries = ceil(p / 2)
        size_d = ceil((len(Integer(entries).bits()) + 1) / 8)
        size_hash = ceil((len(Omega[0].base_ring().cardinality().bits()) + 1) / 8)
        size_entry = size_d + size_hash


        target = impl_curve.hashP(sub_R[p])
        if target == b"\x00" * size_hash:
            return 0, p
        pos = 0
        if not exists(f"./codebook/{NAME}/{p}.bin"):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT),
                                    f"./codebook/{NAME}/{p}.bin")

        print(f"Opening ./codebook/{NAME}/{p}.bin")
        with open(f"./codebook/{NAME}/{p}.bin", "rb") as codebook:
            for i in srange(len(entries.bits()) - 1, -2, -1):
                codebook.seek(size_entry * floor(pos + 2**i), 0)
                hashVal_b = codebook.read(size_hash)
                hashVal = int.from_bytes(hashVal_b, byteorder="big")

                if target > hashVal:
                    if pos + 2**i < entries:
                        pos += 2**i

                if hashVal == target:
                    d = int.from_bytes(codebook.read(size_d), byteorder="big")
                    if impl_curve.lift(d * Omegap) == sub_R[p]:
                        return d, p
                    else:
                        return p-d, p
            print("Not found. This should never happen. Is the codebook correct?")

    if METHOD == "naive":
        u_list = list(map(solve_factorgroup, factors))
    if METHOD == "externalcb":
        u_list = list(map(solve_externalcb, factors))

    print("Solutions from subgroups", u_list)


    res = CRT([d_i for d_i, _ in u_list], [m_i for _, m_i in u_list])
    print("Explicit Point", (res * Omega).xy())

    # Removing randomization
    return (res * Omega).xy()[0] - alpha

if __name__ == "__main__":
    NAME = (sys.argv[1:2] or ["NIST P-256"])[0]
    METHOD = (sys.argv[2:3] or ["externalcb"])[0]
    METHOD = METHOD if METHOD in ["naive", "externalcb"] else "externalcb"
    E, P, k = getCurveChallenge(NAME)
    print(f"Trying to solve a dlog instance on '{NAME}' with method '{METHOD}'")
    print("Solution will be ", k)

    assert(P.order().is_prime())

    # Instantiate the oracle, which knows k (and all intermediate dlogs)
    o = StrongOracle(P, k)

    res = maurer_dlog(E, P, o.public, o, METHOD)

    print()
    print("Found solution ", res, " is correct" if res == o._secret else "is incorrect!")
    print("Expected: ", k)
    print("Oracle calls:", o._calls)
    print("Time spent in oracle:", round(o.oracle.t, 3), "ms")
    print()

    exit(0 if res == o._secret else -1)
