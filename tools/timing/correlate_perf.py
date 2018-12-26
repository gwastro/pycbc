from pycbc.filter import correlate
from pycbc.types import zeros, complex64, complex128
from time import time
from numpy.random import uniform
niter = 1000

for dtp in [complex64, complex128]:
    for N in [2**10, 2**15, 2**20]:
        a = zeros(N, dtype=dtp) 
        a.data += uniform(-1, 1, size=len(a))
        b = a * 0.5
        c = a * 1.5
        correlate(a, b, c)
        
        t1 = time()
        for i in range(niter):
            correlate(a, b, c)
        t2 = time()
        print("Correlate Perf Type:{} Size:{} Time:{:3.3f}".format(repr(dtp),
              N, (t2-t1)*1000 / niter))


