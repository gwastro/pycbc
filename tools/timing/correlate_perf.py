from pycbc.filter import correlate
from pycbc.filter.matchedfilter import BatchCorrelator, Correlator
from pycbc.types import zeros, complex64, complex128, Array
from time import time
from numpy.random import uniform
niter = 2000


for N in [2**10, 2**15, 2**18]:
    a = zeros(N, dtype=complex64)
    a.data += uniform(-1, 1, size=len(a))
    b = a * 0.5
    c = a * 1.5

    xs = [a*1, a*2, a*3]
    zs = [c*1, c*2, c*3]
    corr = BatchCorrelator(xs, zs, N)
    t1 = time()
    for i in range(niter):
        corr.execute(b)

    t2 = time()
    print("Batch Correlate Perf Size:{} Time:{:3.3f}".format(
          N, (t2-t1)*1000 / niter))


for dtp in [complex64, complex128]:
    for N in [2**10, 2**15, 2**20]:
        a = zeros(N, dtype=dtp)
        a += Array(uniform(-1, 1, size=N) * (1 + -.5j), dtype=a.dtype)
        b = zeros(N, dtype=dtp)
        c = zeros(N, dtype=dtp)
        correlate(a, b, c)

        t1 = time()
        for i in range(niter):
            correlate(a, b, c)
        t2 = time()
        print("Correlate Perf Type:{} Size:{} Time:{:3.3f}".format(repr(dtp),
              N, (t2-t1)*1000 / niter))

        if dtp is complex64:
            corr = Correlator(a, b, c)
            t1 = time()
            for i in range(niter):
                corr.correlate()
            t2 = time()
            print("Correlator Perf Type:{} Size:{} Time:{:3.3f}".format(repr(dtp),
                  N, (t2-t1)*1000 / niter))

