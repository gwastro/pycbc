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
from pycbc import WEAVE_FLAGS
from scipy.weave import inline
from .simd_threshold import thresh_cluster_support, default_segsize
from .events import _BaseThresholdCluster
from pycbc.opt import omp_libs, omp_flags

def threshold_numpy(series, value):
    arr = series.data
    locs = numpy.where(arr.real**2 + arr.imag**2 > value**2)[0]
    vals = arr[locs]
    return locs, vals

threshold_only = threshold_numpy

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
                    extra_compile_args=[WEAVE_FLAGS] + omp_flags,
                    libraries=omp_libs
          )
    num = count[0]
    if num > 0:
        return outl[0:num], outv[0:num]
    else:
        return numpy.array([], numpy.uint32), numpy.array([], numpy.float32)

threshold=threshold_inline

class CPUThresholdCluster(_BaseThresholdCluster):
    def __init__(self, series):
        self.series = numpy.array(series.data, copy=False)

        self.slen = len(series)
        self.outv = numpy.zeros(self.slen, numpy.complex64)
        self.outl = numpy.zeros(self.slen, numpy.uint32)
        self.segsize = default_segsize
        self.code = """
             return_val = parallel_thresh_cluster(series, (uint32_t) slen, values, locs,
                                         (float) threshold, (uint32_t) window, (uint32_t) segsize);
              """
        self.support = thresh_cluster_support

    def threshold_and_cluster(self, threshold, window):
        series = self.series
        slen = self.slen
        values = self.outv
        locs = self.outl
        segsize = self.segsize
        self.count = inline(self.code, ['series', 'slen', 'values', 'locs', 'threshold', 'window', 'segsize'],
                            extra_compile_args = [WEAVE_FLAGS] + omp_flags,
                            #extra_compile_args = ['-mno-avx -mno-sse2 -mno-sse3 -mno-ssse3 -mno-sse4 -mno-sse4.1 -mno-sse4.2 -mno-sse4a -O2 -w'] + omp_flags,
                            #extra_compile_args = ['-msse3 -O3 -w'] + omp_flags,
                            support_code = self.support, libraries = omp_libs,
                            auto_downcast = 1)
        if self.count > 0:
            return values[0:self.count], locs[0:self.count]
        else:
            return numpy.array([], dtype = numpy.complex64), numpy.array([], dtype = numpy.uint32)

def _threshold_cluster_factory(series):
    return CPUThresholdCluster
