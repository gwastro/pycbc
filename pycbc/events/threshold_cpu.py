# Copyright (C) 2012  Alex Nitz
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#
import numpy
import events, pycbc
from scipy.weave import inline
from .simd_threshold import thresh_cluster_support, ThreshClusterObject, default_segsize
from .events import _BaseThresholdCluster

if pycbc.HAVE_OMP:
    omp_libs = ['gomp']
    omp_flags = ['-fopenmp']
else:
    omp_libs = []
    omp_flags = []

def threshold_numpy(series, value):
    arr = series.data
    locs = numpy.where(arr.real**2 + arr.imag**2 > value**2)[0]
    vals = arr[locs]
    return locs, vals

outl = None
outv = None
count = None
def threshold_inline(series, value):
    arr = numpy.array(series.data.view(dtype=numpy.float32), copy=False)
    global outl, outv, count
    if outl is None or len(outl) < len(series):
        outl = numpy.zeros(len(series), dtype=numpy.uint32)
        outv = numpy.zeros(len(series), dtype=numpy.complex64)
        count = numpy.zeros(1, dtype=numpy.uint32)
        
    N = len(series)
    threshold = value**2.0
    code = """  
        float v = threshold;
        unsigned int num_parallel_regions = 16;
        unsigned int t=0;
     
        #pragma omp parallel for ordered shared(t)
        for (unsigned int p=0; p<num_parallel_regions; p++){
            unsigned int start  = (N * p) / num_parallel_regions;
            unsigned int end    = (N * (p+1)) / num_parallel_regions;
            unsigned int c = 0;
            
            for (unsigned int i=start; i<end; i++){
                float r = arr[i*2];
                float im = arr[i*2+1];
                if ((r * r + im * im) > v){
                    outl[c+start] = i;
                    outv[c+start] = std::complex<float>(r, im);
                    c++;
                }
            } 
            
            #pragma omp ordered
            {
                t+=c;
            }
            memmove(outl+t-c, outl+start, sizeof(unsigned int)*c);
            memmove(outv+t-c, outv+start, sizeof(std::complex<float>)*c);

        }       
        
        count[0] = t;
    """
    inline(code, ['N', 'arr', 'outv', 'outl', 'count', 'threshold'],
                    extra_compile_args=['-march=native -O3 -w'] + omp_flags,
                    libraries=omp_libs
          )
    num = count[0]
    if num > 0:
        return outl[0:num], outv[0:num]
    else:
        return numpy.array([], numpy.uint32), numpy.array([], numpy.float32)

threshold=threshold_inline

# The CUDA function we are trying to emulate seems to hard-code
# choices based on a 4096 window, so for now we do the same. This means
# we can have up to 256 points in a 2^20 size SNR time series.
simd_outv = numpy.zeros(256, dtype = numpy.complex64)
simd_outl = numpy.zeros(256, dtype = numpy.uint32)
simd_count = 0
def threshold_and_cluster(series, threshold, window):
    code = """
    return_val = parallel_thresh_cluster(series, (uint32_t) slen, values, locs, 
                                         (float) threshold, (uint32_t) window, (uint32_t) segsize);
    """
    series = numpy.array(series.data, copy = False)
    slen = len(series)
    values = simd_outv
    locs = simd_outl
    segsize = default_segsize
    simd_count = inline(code, ['series', 'slen', 'values', 'locs', 'threshold', 'window', 'segsize'],
                        extra_compile_args = ['-march=native -O3 -w'] + omp_flags,
                        support_code = thresh_cluster_support, libraries = omp_libs,
                        auto_downcast = 1)
    if simd_count > 0:
        return values[0:simd_count], locs[0:simd_count]
    else:
        return numpy.array([], dtype = numpy.complex64), numpy.array([], dtype = numpy.uint32)
    
def CPUThresholdCluster(_BaseThresholdCluster):
    def __init__(self, series, threshold, window):
        self.series = numpy.array(series.data, copy = False)
        self.slen = len(series)
        self.threshold = threshold
        self.window = window
        self.outlen = len(series)/window
        self.outv = numpy.zeros(self.outlen, numpy.complex64)
        self.outl = numpy.zeros(self.outlen, numpy.uint32)
        self.segsize = default_segsize
        self.code = """
             return_val = parallel_thresh_cluster(series, (uint32_t) slen, values, locs, 
                                         (float) threshold, (uint32_t) window, (uint32_t) segsize);
              """
        self.support = thresh_cluster_support

    def threshold_and_cluster(self):
        series = self.series
        slen = self.slen
        values = self.outv
        locs = self.outl
        threshold = self.threshold
        window = self.window
        segsize = self.segsize
        self.count = inline(self.code, ['series', 'slen', 'values', 'locs', 'threshold', 'window', 'segsize'],
                            extra_compile_args = ['-march=native -O3 -w'] + omp_flags,
                            support_code = self.support, libraries = omp_libs,
                            auto_downcast = 1)
        if self.count > 0:
            return values[0:self.count], locs[0:self.count]
        else:
            return numpy.array([], dtype = numpy.complex64), numpy.array([], dtype = numpy.uint32)

def _threshold_cluster_factory(series, threshold, window):
    return CPUThresholdCluster
