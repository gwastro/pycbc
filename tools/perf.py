#!/usr/bin/python
from pycbc.scheme import *
from pycbc.types import *
from pycbc.filter import *
from pycbc.psd import *
import pycbc
from math import log
import numpy
import sys
from optparse import OptionParser
from math import sin, log
import gc
parser = OptionParser()

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
    ctx = DefaultScheme()
if _options['scheme'] == 'cuda':
    ctx = CUDAScheme(device_num=_options['devicenum'])
if _options['scheme'] == 'opencl':
    ctx = OpenCLScheme(device_num=_options['devicenum'])




size = options.size
niter = options.iterations


if type(ctx) is CUDAScheme:
    print "RUNNING ON ", ctx.device.name()
else:
    print "RUNNING ON CPU"

N = 2**size
print "         SIZE    ",  int(log(N,2))
n = N/2 +1
noise = numpy.arange(0,N,1)
nplus2 = TimeSeries(noise,delta_t=1.0/4096,dtype=float32)
ntilde2 = make_frequency_series(nplus2)
psd2 = ntilde2.squared_norm()    

o,ind = match(ntilde2,ntilde2,psd=psd2)
with ctx:
    o,ind = match(ntilde2,ntilde2,psd=psd2)
    o,ind = match(ntilde2,ntilde2,psd=None, h_norm=1, s_norm=1)
    o,ind = matched_filter(ntilde2,ntilde2)

def matcht():
    with ctx:
        for i in range(0,niter):
            o,ind = match(ntilde2,ntilde2,psd=psd2)

def match_fast():
    with ctx:
        for i in range(0,niter):
            o,ind = match(ntilde2,ntilde2,psd=None,h_norm=1,s_norm=1)

def filter_fast():
    with ctx:
        for i in range(0,niter):
            o,ind = matched_filter(ntilde2,ntilde2,psd=None,h_norm=1)

import timeit
gt = timeit.Timer(matcht)
t = (1000 * gt.timeit(number=1)/niter)
print "MATCH %.2f msec" % t, " %5.1f op/min " % (1000 *60 /t)


gt = timeit.Timer(match_fast)
t = (1000 * gt.timeit(number=1)/niter)
print "MATCH FAST %.2f msec" % t, " %5.1f op/min " % (1000 *60 /t)


gt = timeit.Timer(filter_fast)
t = (1000 * gt.timeit(number=1)/niter)
print "FILTER FAST %.2f msec" % t, " %5.1f op/min " % (1000 *60 /t)




