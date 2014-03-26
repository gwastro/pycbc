#!/usr/bin/python
from pycbc.scheme import *
from pycbc.types import *
from pycbc.fft import *
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
parser.add_option('--measure-level', type=int, help='Set the measure level (only applies to FFTW- cpu scheme)', default=1)
parser.add_option('--backend', type=str, help='set the backend type for this scheme', default='Default')
         
          
(options, args) = parser.parse_args()   

#Changing the optvalues to a dict makes them easier to read
_options = vars(options)

if _options['scheme'] == 'cpu':
    ctx = CPUScheme()
    from pycbc.fft.lalfft import set_measure_level
    set_measure_level(options.measure_level)
if _options['scheme'] == 'cuda':
    ctx = CUDAScheme(device_num=_options['devicenum'])
if _options['scheme'] == 'opencl':
    ctx = OpenCLScheme(device_num=_options['devicenum'])


niter = options.iterations


if type(ctx) is CUDAScheme:
    print "RUNNING ON ", ctx.device.name()
else:
    print "RUNNING ON CPU"
print type(ctx)

N = 2**options.size

vecin = zeros(N, dtype=complex64) + 1
vecout = vecin * 1

vecdin = zeros(N, dtype=complex128) + 1
vecdout = vecdin * 1

with ctx:
    fft(vecin, vecout, backend=options.backend)
    ifft(vecin, vecout, backend=options.backend)

def tfft():
    with ctx:
	for i in range(0, niter):
	    fft(vecin, vecout, backend=options.backend)
        b = vecout[9]
def tifft():
    with ctx:
	for i in range(0, niter):
	    ifft(vecin, vecout, backend=options.backend)
        b = vecout[9]

def dtifft():
    with ctx:
        for i in range(0, niter):
            ifft(vecdin, vecdout, backend=options.backend)
        b = vecout[9]

import timeit
gt = timeit.Timer(tfft)
t = (1000 * gt.timeit(number=1)/niter)
print "C2C FFT %.2f msec" % t, " %5.1f /min " % (1000 *60 /t)

gt = timeit.Timer(tifft)
t = (1000 * gt.timeit(number=1)/niter)
print "C2C iFFT %.2f msec" % t, " %5.1f /min " % (1000 *60 /t)

gt = timeit.Timer(dtifft)
t = (1000 * gt.timeit(number=1)/niter)
print "C2C DOUBLE iFFT %.2f msec" % t, " %5.1f /min " % (1000 *60 /t)

