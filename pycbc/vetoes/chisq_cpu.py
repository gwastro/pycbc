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
import numpy, pycbc
from pycbc.types import Array, real_same_precision_as, complex_same_precision_as
from scipy.weave import inline

if pycbc.HAVE_OMP:
    omp_libs = ['gomp']
    omp_flags = ['-fopenmp']
else:
    omp_libs = []
    omp_flags = []

def chisq_accum_bin_numpy(chisq, q):
    chisq += q.squared_norm()
    
def chisq_accum_bin_inline(chisq, q):
    
    chisq = numpy.array(chisq.data, copy=False)
    q = numpy.array(q.data, copy=False)
    N = len(chisq)
    code = """
        #pragma omp parallel for
        for (int i=0; i<N; i++){
            chisq[i] += q[i].real()*q[i].real()+q[i].imag()*q[i].imag();
        }
    """
    inline(code, ['chisq', 'q', 'N'], 
                    extra_compile_args=['-march=native -O3 -w'] + omp_flags,
                    libraries=omp_libs
          )
          
chisq_accum_bin = chisq_accum_bin_inline

point_chisq_code = """
    int num_parallel_regions = 256 / num_bins;
    TYPE* __restrict__ shift = shifts;
    TYPE* __restrict__ v = v1;
    TYPE* __restrict__ pr = (TYPE*) malloc(sizeof(TYPE)*n);
    TYPE* __restrict__  pi = (TYPE*) malloc(sizeof(TYPE)*n);
    TYPE* __restrict__  vsr = (TYPE*) malloc(sizeof(TYPE)*n);
    TYPE* __restrict__  vsi = (TYPE*) malloc(sizeof(TYPE)*n);
    TYPE* __restrict__ outr = (TYPE* ) malloc(sizeof(TYPE)*n);
    TYPE* __restrict__ outi = (TYPE*) malloc(sizeof(TYPE)*n);
    
    TYPE TWO_PI = 6.28318530718;
    
    for (int i=0; i<n; i++){
        vsr[i] = cos(TWO_PI * shift[i] / slen);
        vsi[i] = sin(TWO_PI * shift[i] / slen);
    }
    
    for (unsigned int r=0; r<num_bins; r++){     
        int bstart = bins[r];
        int bend = bins[r+1];
        int blen = bend - bstart;

        for (int i=0; i<n; i++){
            outr[i] = 0;
            outi[i] = 0;
        }
    
        for (unsigned int k=0; k<num_parallel_regions; k++){
            unsigned int start = blen * k / num_parallel_regions + bstart;
            unsigned int end = blen * (k + 1) / num_parallel_regions + bstart;
        
            //start the cumulative rotations at the offset point
            for (int i=0; i<n; i++){
                pr[i] = cos(TWO_PI * shift[i] * (start) / slen);
                pi[i] = sin(TWO_PI * shift[i] * (start) / slen);
            }
            
            TYPE t1, t2, k1, k2, k3, vs, va;          
            for (unsigned int j=start; j<end; j++){
                TYPE vr = v[j*2];
                TYPE vi = v[j*2+1];  
                vs = vr + vi;
                va = vi - vr;
                    
                for (int i=0; i<n; i++){            
                    t1 = pr[i];
                    t2 = pi[i];
                                
                    // Complex multiply pr[i] * v
                    k1 = vr * (t1 + t2);
                    k2 = t1 * va;
                    k3 = t2 * vs;
                                
                    outr[i] += k1 - k3;
                    outi[i] += k1 + k2;
                    
                    // phase shift for the next time point
                    pr[i] = t1 * vsr[i] - t2 * vsi[i];
                    pi[i] = t1 * vsi[i] + t2 * vsr[i]; 
                }                                              
            } 
               
        }    
        
        for (unsigned int i=0; i<n; i++){
            chisq[i] += outr[i]*outr[i] + outi[i]*outi[i];
        }
        
    }  
    
    free(pr);
    free(pi);  
    free(vsr);
    free(vsi);      
    free(outr);
    free(outi);  
"""

point_chisq_code_single = point_chisq_code.replace('TYPE', 'float')
point_chisq_code_double = point_chisq_code.replace('TYPE', 'double')

def shift_sum(v1, shifts, bins):
    real_type = real_same_precision_as(v1)
    shifts = numpy.array(shifts, dtype=real_type)
    
    bins = numpy.array(bins, dtype=numpy.uint32)
    num_bins = len(bins) - 1
    v1 = numpy.array(v1.data, copy=False).view(dtype=real_type)
    slen = len(v1)
    
    if v1.dtype.name == 'float32':
        code = point_chisq_code_single
    else:
        code = point_chisq_code_double
    
    n = int(len(shifts))
    
    # Create some output memory
    chisq =  numpy.zeros(n, dtype=real_type)
    
    inline(code, ['v1', 'n', 'chisq', 'slen', 'shifts', 'bins', 'num_bins'],
                    extra_compile_args=['-march=native -O3 -w -ffast-math -ftree-vectorizer-verbose=6'],
                    libraries=omp_libs
          )
          
    return  Array(chisq, copy=False)
