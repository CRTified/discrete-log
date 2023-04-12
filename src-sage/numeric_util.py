from collections import Counter

def ecmfToTupels(l):
    return [(Integer(a), Integer(b)) for a, b in Counter(l).items()]

def smoothness(factors):
    return max([pi for pi, _ in factors])

def genericSquareAndMultiply(operation):
    def squareAndMultiply(base, exponent):
        tmp = base
        for b in exponent.bits()[-2::-1]:
            # Square
            tmp = operation(tmp, tmp)
            if b:
                # Multiply
                tmp = operation(tmp, base)
        return tmp
    return squareAndMultiply
