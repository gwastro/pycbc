#!/usr/bin/python
from pycbc.scheme import *
from pycbc.types import *
from pycbc.waveform import *
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

parser.add_option('--approximant', type=str, default="TaylorF2")
parser.add_option('--deltaf',type=float, help='frequency step')
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
    print "RUNNING ON ", ctx.device.name()
else:
    print "RUNNING ON CPU"


with ctx:
    wf_taylor = get_fd_waveform(mass1=1, mass2=1, f_lower=14, 
                                approximant=options.approximant, delta_f=options.deltaf)

def taylorf2():
    with ctx:
        for i in range(0,niter):
            wf_taylor = get_fd_waveform(mass1=1, mass2=1, f_lower=14, 
                                        approximant=options.approximant, delta_f=options.deltaf)


import timeit
gt = timeit.Timer(taylorf2)
t = (1000 * gt.timeit(number=1)/niter)
print "Waveform Generation %.2f msec" % t, " %5.1f gen/min " % (1000 *60 /t)

if type(ctx) is CUDAScheme:
    def SPAtmplt():
        with ctx:
            n = int(1.0 / options.deltaf * 4096)
            out = zeros(n, dtype=complex64)
            for i in range(0,niter):
                wf_taylor = get_fd_waveform(mass1=1, mass2=1, f_lower=14, 
                                            approximant="SPAtmplt", delta_f=options.deltaf, out=out, amplitude_order=0)

    gt = timeit.Timer(SPAtmplt)
    t = (1000 * gt.timeit(number=1)/niter)
    print "SPAtmplt Generation %.2f msec" % t, " %5.1f gen/min " % (1000 *60 /t)


