import multiprocessing as mp
from collections import Counter
from curve_list import CURVES
import os

def smoothness(factors):
    return max([pi for pi, _ in factors])

def log_smoothness(s):
    return len(s.bits())

def random_curve(F):
  # Choose random values for a and b
  a, b = F.random_element(), F.random_element()
  # Keep choosing new values for a and b until the curve is non-singular
  while 4*a**3 + 27*b**2 == 0:
    a, b = F.random_element(), F.random_element()
  # Construct an elliptic curve in short Weierstrass form
  E = EllipticCurve(F, [a, b])
  # Compute the order of the curve
  order = E.order()
  factors = factor(order)
  return E, factors, smoothness(factors)

def many_random_curves(F):
    def timeout_rc():
        return random_curve(F)
    while True:
        try:
            res = timeout_rc()
            if res:
                yield res
        except:
            print(" T", end="")
            return

def readlastline(f):
    f.seek(-2, 2)              # Jump to the second last byte.
    while f.read(1) != b"\n":  # Until EOL is found ...
        f.seek(-2, 1)          # ... jump back, over the read byte plus one more.
    return f.read()            # Read all data from this point on.


def proc(F, q, cqu):
    set_random_seed()
    for c in many_random_curves(F):
        q.put(c)
        if not cqu.get():
            q.put("magic")
            break

for NAME in CURVES.keys():
    target = CURVES.get(NAME, None)
    if not target:
        exit(1)
    THRESHOLD = 40
    E, B, expected_order = target

    print("Targeting ", NAME)

    q = B.order()
    if expected_order:
        assert (q == expected_order)

    print("Targeting {}".format(E))
    print("Basepoint B={}".format(B))
    print("with order q={}".format(q))
    print("q is prime?",  "Yes" if q.is_prime() else "No")

    F = GF(q)

    os.makedirs("data", exist_ok=True)
    try:
        with open("data/Curve-{}.txt".format(NAME), "rb") as f:
            line = readlastline(f).decode().split(" ")
        A = line[0]
        B = line[1]
        print("Read {} {}".format(A, B))
        best_c = EllipticCurve(F, [F(A), F(B)])
        best_f = best_c.order().factor()
        best_s = smoothness(best_f)
    except:
        print("Choosing random starting point")
        best_c, best_f, best_s = random_curve(F)

    print("Starting with:")
    print(best_c)
    print("# Order 2^{}-smooth, Factors: {}".format(str(log_smoothness(best_s)), str(best_f)))

    qu = mp.Queue()
    controlqu = mp.Queue()

    procs = [mp.Process(target=proc, args=(F, qu,controlqu)) for _ in range(mp.cpu_count() - 2)]
    for p in procs:
        controlqu.put(True)
        p.start()

    i = 1
    with open("data/log-{}.txt".format(NAME), "a") as fulllog:
        while True:
            try:
                res = qu.get(True, 5)
                if res == "magic":
                    continue
            except:
                if log_smoothness(best_s) < THRESHOLD:
                    break
                else:
                    continue
            c, f, s = res

            # i A B log2(s) factors
            fulllog.write("{:6}\t{}\t{}\t{}\t\"{}\"\n".format(i, c.a4(), c.a6(), log_smoothness(s), f))
            if i & 0xFF == 0:
                fulllog.flush()
                print("\r{:10}          ".format(i), end='', flush=True)

            print("\r{:10}".format(i), end='', flush=True)
            i += 1
            if s < best_s:
                best_c, best_f, best_s = c, f, s
                print()
                print(c,f,s)
                print()
                print("New best:")
                print("# Order 2^{}-smooth, Factors: {}".format(str(log_smoothness(best_s)), str(best_f)))
                print("A = F({})".format(best_c.a4()))
                print("B = F({})".format(best_c.a6()))
                print()
                with open("data/Curve-{}.txt".format(NAME), "a") as logfile:
                    logfile.write("{} {}\n".format(
                        best_c.a4(),
                        best_c.a6()
                    ))
            controlqu.put(log_smoothness(best_s) >= THRESHOLD)
    print("Found a good curve")
    for p in procs:
        controlqu.put(False)
