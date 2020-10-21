#!/usr/bin/env python
from pycbc.scheme import *
from pycbc.types import *
from pycbc.filter import *
from pycbc.psd import *
import pycbc
from math import log
import numpy
import numpy.random
import sys
from optparse import OptionParser
from math import sin
import gc
parser = OptionParser()

import logging
logging.basicConfig(format='%(asctime)s : %(message)s', level=logging.DEBUG)

parser.add_option('--scheme','-s',  type = 'choice',
                    choices = ('cpu','cuda','opencl'),
                    default = 'cpu', dest = 'scheme',
                    help = 'specifies processing scheme, can be cpu [default], cuda, or opencl')

parser.add_option('--device-num','-d', action='store', type = 'int',
                    dest = 'devicenum', default=0,
                    help = 'specifies a GPU device to use for CUDA or OpenCL, 0 by default')

parser.add_option('--size', type=int, help='fft size in log2')
parser.add_option('--iterations', type=int, help='Number of iterations to perform')

(options, args) = parser.parse_args()

#Changing the optvalues to a dict makes them easier to read
_options = vars(options)

if _options['scheme'] == 'cpu':
    ctx = CPUScheme()
if _options['scheme'] == 'cuda':
    ctx = CUDAScheme(device_num=_options['devicenum'])
if _options['scheme'] == 'opencl':
    ctx = OpenCLScheme(device_num=_options['devicenum'])




size = options.size
niter = options.iterations


if type(ctx) is CUDAScheme:
    print("RUNNING ON ", ctx.device.name())
else:
    print("RUNNING ON CPU")

N = 2**size
print("         SIZE    ",  int(log(N,2)))
n = N/2 +1
a = numpy.zeros(N) + 1000
noise = numpy.random.normal(a).astype(numpy.float32)

with ctx:
    nplus2 = TimeSeries(noise,delta_t=1.0/4096,dtype=float32)
    ntilde2 = make_frequency_series(nplus2)
    psd2 = ntilde2.squared_norm()
    o = match(ntilde2,ntilde2,psd=psd2)
    o = match(ntilde2,ntilde2,psd=None, v1_norm=1, v2_norm=1)
    o = matched_filter_core(ntilde2, ntilde2)
    out=zeros(N,dtype=complex64)
    o = overlap_cplx(ntilde2, ntilde2, normalized=False)
    ntilde3 = ntilde2 +10j

def matcht():
    with ctx:
        for i in range(0,niter):
            o,ind = match(ntilde2,ntilde2,psd=psd2)

def match_fast():
    with ctx:
        for i in range(0,niter):
            o,ind = match(ntilde2,ntilde2,psd=None,v1_norm=1,v2_norm=1)

def ovlp():
    with ctx:
        for i in range(0,niter):
            o = overlap_cplx(ntilde2,ntilde3, normalized=False)

def filter_fast():
    with ctx:
        for i in range(0,niter):
            snr, corr, norm = matched_filter_core(ntilde2, ntilde2, psd=None, h_norm=1, out=out)

import timeit
gt = timeit.Timer(ovlp)
t = (1000 * gt.timeit(number=1)/niter)
print("Foverlap %.2f msec" % t, " %5.1f op/min " % (1000 *60 /t))

gt = timeit.Timer(matcht)
t = (1000 * gt.timeit(number=1)/niter)
print("MATCH %.2f msec" % t, " %5.1f op/min " % (1000 *60 /t))


gt = timeit.Timer(match_fast)
t = (1000 * gt.timeit(number=1)/niter)
print("MATCH FAST %.2f msec" % t, " %5.1f op/min " % (1000 *60 /t))


gt = timeit.Timer(filter_fast)
t = (1000 * gt.timeit(number=1)/niter)
print("FILTER FAST %.2f msec" % t, " %5.1f op/min " % (1000 *60 /t))




