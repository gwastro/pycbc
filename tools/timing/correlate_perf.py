from time import time

from numpy.random import uniform

from pycbc.filter import correlate
from pycbc.filter.matchedfilter import BatchCorrelator, Correlator
from pycbc.types import Array, complex64, complex128, zeros

niter = 2000


for N in [2**10, 2**15, 2**18]:
    a = zeros(N, dtype=complex64)
    a.data += uniform(-1, 1, size=len(a))
    b = a * 0.5
    c = a * 1.5

    xs = []
    zs = []
    for i in range(50):
        xs.append(a * i)
        zs.append(c * i)
    corr = BatchCorrelator(xs, zs, N)
    t1 = time()
    for i in range(niter):
        corr.execute(b)

    t2 = time()
    print(f"Batch Correlate Perf Size:{N} Time:{(t2 - t1) * 1000 / niter:3.3f}")


for dtp in [complex64, complex128]:
    for N in [2**10, 2**15, 2**20]:
        a = zeros(N, dtype=dtp)
        a += Array(uniform(-1, 1, size=N) * (1 + -0.5j), dtype=a.dtype)
        b = zeros(N, dtype=dtp)
        c = zeros(N, dtype=dtp)
        correlate(a, b, c)

        t1 = time()
        for i in range(niter):
            correlate(a, b, c)
        t2 = time()
        print(
            f"Correlate Perf Type:{dtp!r} Size:{N} Time:{(t2 - t1) * 1000 / niter:3.3f}"
        )

        if dtp is complex64:
            corr = Correlator(a, b, c)
            t1 = time()
            for i in range(niter):
                corr.correlate()
            t2 = time()
            print(
                f"Correlator Perf Type:{dtp!r} Size:{N} Time:{(t2 - t1) * 1000 / niter:3.3f}"
            )
