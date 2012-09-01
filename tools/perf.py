from pycbc.scheme import *
from pycbc.types import *
from pycbc.filter import *
from pycbc.psd import *
import pycbc
from math import log
import numpy
import sys
from numpy.random import normal
from math import sin, log
import gc
ctx = CUDAScheme()
print "SETUP"

for N in [2**18, 2**19, 2**20, 2**21, 2**22]:
    print "         SIZE    ",  int(log(N,2))
    n = N/2 +1
    noise = normal(0.0,2,N) +5
    nplus2 = TimeSeries(noise,delta_t=1.0/4096,dtype=float32)
    ntilde2 = make_frequency_series(nplus2)
    psd2 = ntilde2.squared_norm()    

    o,ind = match(ntilde2,ntilde2,psd=psd2)
    with ctx:
        o,ind = match(ntilde2,ntilde2,psd=psd2)
        o,ind = match(ntilde2,ntilde2,psd=None)

    #o,ind = match(ntilde2,ntilde2,psd=psd2)

    with ctx:
        o,ind = matched_filter(ntilde2,ntilde2)

    ngpu = 5000
    ncpu = 50

    def cpu_s():
        for i in range(0,ncpu):
            o,ind = match(ntilde2,ntilde2,psd=psd2)

    def gpu_s():
        with ctx:
            for i in range(0,ngpu):
                o,ind = match(ntilde2,ntilde2,psd=psd2)

    def cpu_st():
        for i in range(0,ncpu):
            o,ind = match(ntilde2,ntilde2,psd=None,h_norm=1,s_norm=1)

    def gpu_st():
        with ctx:
            for i in range(0,ngpu):
                o,ind = match(ntilde2,ntilde2,psd=None,h_norm=1,s_norm=1)

    def cpu_s2():
        for i in range(0,ncpu):
            o,ind = matched_filter(ntilde2,ntilde2,psd=None,h_norm=1)


    def gpu_s2():
        with ctx:
            for i in range(0,ngpu):
                o,ind = matched_filter(ntilde2,ntilde2,psd=None,h_norm=1)

    def sig():
        with ctx:
            for i in range(0,ngpu):
                a = sigmasq(ntilde2)

    def sq():
        with ctx:
            for i in range(0,ngpu):
                ntilde2.squared_norm()

    def sm():
        with ctx:
            for i in range(0,ngpu):
                ntilde2.sum()

    import timeit
    gt = timeit.Timer(sig)
    t = (1000 * gt.timeit(number=1)/ngpu)
    print "GPU sigma %.2f msec/op" % t, " %5.1f op/min " % (1000 *60 /t)

    #import timeit
    #gt = timeit.Timer(sq)
    # = (1000 * gt.timeit(number=1)/ngpu)
    #print "GPU square %.2f msec/op" % t, " %5.1f op/min " % (1000 *60 /t)

    #import timeit
    #gt = timeit.Timer(sm)
    #t = (1000 * gt.timeit(number=1)/ngpu)
    #print "GPU sum %.2f msec/op" % t, " %5.1f op/min " % (1000 *60 /t)

    gt = timeit.Timer(gpu_s)
    t = (1000 * gt.timeit(number=1)/ngpu)
    print "GPU single %.2f msec/match" % t, " %5.1f matches/min " % (1000 *60 /t)
    pt = t
    gt = timeit.Timer(cpu_s)
    t = (1000 * gt.timeit(number=1)/ncpu)
    print "CPU single %.2f msec/match" % t, " %5.1f matches/min " % (1000 *60 /t)
    print "          %3.0f X Speedup" % (t/pt)

    gt = timeit.Timer(gpu_st)
    t = (1000 * gt.timeit(number=1)/ngpu)
    print "GPU single(nopsd) %.2f msec/match" % t, " %5.1f matches/min " % (1000 *60 /t)
    pt = t
    gt = timeit.Timer(cpu_st)
    t = (1000 * gt.timeit(number=1)/ncpu)
    print "CPU single(nopsd) %.2f msec/match" % t, " %5.1f matches/min " % (1000 *60 /t)
    print "          %3.0f X Speedup" % (t/pt)

    gt = timeit.Timer(gpu_s2)
    t = (1000 * gt.timeit(number=1)/ngpu)
    print "GPU single snrfast %.2f msec/op" % t, " %5.1f op/min " % (1000 *60 /t)
    pt = t
    gt = timeit.Timer(cpu_s2)
    t = (1000 * gt.timeit(number=1)/ncpu)
    print "CPU single snrfast %.2f msec/op" % t, " %5.1f op/min " % (1000 *60 /t)
    print "          %3.0f X Speedup" % (t/pt)



