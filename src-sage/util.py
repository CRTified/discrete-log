import time
from functools import wraps

def countCalls(f):
    """
    Decorator to count function calls
    """
    @wraps(f)
    def w(*args, **kwargs):
        w.calls += 1
        return f(*args, **kwargs)
    def resetCounter():
        w.calls = 0
        w.calls = 0
        w.resetCounter = resetCounter
    return w

def repeatUntilNotNone(f):
    def wrapped(*args, **kw):
        res = f(*args, **kw)
        while res == None:
            res = f(*args, **kw)
        return res
    return wrapped

def ppDict(d):
    """
    prettyprinter for dicts
    """
    print("{")
    for k, v in d.items():
        print("  {:30}: {},".format(str(k), str(v)))
        print("}")

def timeit(method):
    """
    Decorator to measure wall time of code execution
    """
    @wraps(method)
    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()
        print('%r  %2.2f ms' % \
              (method.__name__, (te - ts) * 1000))
        return result
    return timed

def timeall(method):
    """
    Decorator to measure total wall time of all calls
    """
    @wraps(method)
    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()
        timed.t = timed.t + (te - ts) * 1000
        return result
    timed.t = 0
    timed.calls = 0
    return timed
