#!/usr/bin/python
from pycbc.scheme import *
from pycbc.types import *
from pycbc.fft import *
from pycbc.events import *
import pycbc
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


parser.add_option('--size',type=float, help='FFT size')
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


niter = options.iterations


if type(ctx) is CUDAScheme:
    print("RUNNING ON ", ctx.device.name())
else:
    print("RUNNING ON CPU")
print(type(ctx))

N = 2**options.size

v = zeros(N, dtype=complex64) + 1

with ctx:
    v+v
    v*v
    v.squared_norm()

def addc():
    with ctx:
        for i in range(0, niter):
            v + 3
        v[0]

def add():
    with ctx:
	for i in range(0, niter):
	    v+v
	v[0]

def mul():
    with ctx:
	for i in range(0, niter):
	    v*v
	v[0]

def sqnm():
    with ctx:
	for i in range(0, niter):
	    v.squared_norm()
	v[0]

import timeit

gt = timeit.Timer(addc)
t = (1000 * gt.timeit(number=1)/niter)
print("ADDC  %.2f msec" % t, " %5.1f /min " % (1000 *60 /t))

gt = timeit.Timer(add)
t = (1000 * gt.timeit(number=1)/niter)
print("ADD  %.2f msec" % t, " %5.1f /min " % (1000 *60 /t))

gt = timeit.Timer(mul)
t = (1000 * gt.timeit(number=1)/niter)
print("MUL %.2f msec" % t, " %5.1f /min " % (1000 *60 /t))

gt = timeit.Timer(sqnm)
t = (1000 * gt.timeit(number=1)/niter)
print("SQNRM  %.2f msec" % t, " %5.1f /min " % (1000 *60 /t))


