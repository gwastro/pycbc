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
from pycbc.types import real_same_precision_as
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
    int num_parallel_regions = 16;
    
    for (unsigned int r=0; r<blen; r++){     
        int bstart = bins[r];
        int bend = bins[r+1];
        int blen = bend - bstart;
        
        TYPE* outr = (TYPE*) malloc(sizeof(TYPE)*n);
        TYPE* outi = (TYPE*) malloc(sizeof(TYPE)*n);
        for (int i=0; i<n; i++){
            outr[i] = 0;
            outi[i] = 0;
        }
    
        #pragma omp parallel for
        for (unsigned int k=0; k<num_parallel_regions; k++){
            unsigned int start = blen * k / num_parallel_regions + bstart;
            unsigned int end = blen * (k + 1) / num_parallel_regions + bstart;
        
            //start the cumulative rotations at the offset point
            TYPE* pr = (TYPE*) malloc(sizeof(TYPE)*n);
            TYPE* pi = (TYPE*) malloc(sizeof(TYPE)*n);
            TYPE* vsr = (TYPE*) malloc(sizeof(TYPE)*n);
            TYPE* vsi = (TYPE*) malloc(sizeof(TYPE)*n);
            TYPE* outr_tmp = (TYPE*) malloc(sizeof(TYPE)*n);
            TYPE* outi_tmp = (TYPE*) malloc(sizeof(TYPE)*n);

            
            for (int i=0; i<n; i++){
                pr[i] = cos(2 * 3.141592653 * shifts[i] * (start) / slen);
                pi[i] = sin(2 * 3.141592653 * shifts[i] * (start) / slen);
                vsr[i] = cos(2 * 3.141592653 * shifts[i] / slen);
                vsi[i] = sin(2 * 3.141592653 * shifts[i] / slen);
                outr_tmp[i] = 0;
                outi_tmp[i] = 0;
            }
            
            TYPE t1, t2, k1, k2, k3, vs, va;
            
            for (unsigned int j=start; j<end; j++){
                std::complex<TYPE> v = v1[j];
                TYPE vr = v.real();
                TYPE vi = v.imag();  
                vs = vr + vi;
                va = vi - vr;
                
                for (int i=0; i<n; i++){
                    t1 = pr[i];
                    t2 = pi[i];
                    
                    // Complex multiply pr[i] * v
                    k1 = vr * (t1 + t2);
                    k2 = t1 * va;
                    k3 = t2 * vs;
                                
                    outr_tmp[i] += k1 - k3;
                    outi_tmp[i] += k1 + k2;
                    
                    // phase shift for the next time point
                    pr[i] = t1 * vsr[i] - t2 * vsi[i];
                    pi[i] = t1 * vsi[i] + t2 * vsr[i]; 
                }                                              
            } 
            
            #pragma omp critical
            {
                for (unsigned int i=0; i<n; i++){
                    outr[i] += outr_tmp[i];
                    outi[i] += outi_tmp[i];
                }
            }
            
            free(pr);
            free(pi);  
            free(outr_tmp);
            free(outi_tmp); 
            free(vsr);
            free(vsi);
            
        }    
        
        for (unsigned int i=0; i<n; i++){
            chisq[i] += outr[i]*outr[i] + outi[i]*outi[i];
        }
        free(outr);
        free(outi);
    }    
"""

point_chisq_code_single = point_chisq_code.replace('TYPE', 'float')
point_chisq_code_double = point_chisq_code.replace('TYPE', 'double')

def shift_sum(v1, shifts, bins):
    real_type = real_same_precision_as(v1)
    shifts = numpy.array(shifts, dtype=real_type)
    
    bins = numpy.array(bins, dtype=numpy.uint32)
    blen = len(bins) - 1
    v1 = numpy.array(v1.data, copy=False)
    slen = len(v1)

    if v1.dtype.name == 'complex64':
        code = point_chisq_code_single
    else:
        code = point_chisq_code_double
    
    n = int(len(shifts))
    
    # Create some output memory
    chisq =  numpy.zeros(n, dtype=real_type)
    
    inline(code, ['v1', 'n', 'chisq', 'slen', 'shifts', 'bins', 'blen'],
                    extra_compile_args=['-march=native -O3 -w'] + omp_flags,
                    libraries=omp_libs
          )
          
    return  chisq
