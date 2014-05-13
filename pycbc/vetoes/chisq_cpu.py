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
from pycbc.types import Array
from scipy.weave import inline

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
           extra_compile_args=['-march=native  -O3  -fopenmp'],
           libraries=['gomp']
          )
          
chisq_accum_bin = chisq_accum_bin_inline

def shift_sum(v1, shifts, slen=None, offset=0):
    v1 = numpy.array(v1.data, copy=False)
    shifts = numpy.array(shifts, dtype=numpy.float32)
    vlen = len(v1)
    if slen is None:
        slen = vlen

    code = """
        int num_parallel_regions = 16;
        
        #pragma omp parallel for
        for (unsigned int k=0; k<num_parallel_regions; k++){
            unsigned int start = vlen * k / num_parallel_regions;
            unsigned int end = vlen * (k + 1) / num_parallel_regions;
        
            //start the cumulative rotations at the offset point
            float* pr = (float*) malloc(sizeof(float)*n);
            float* pi = (float*) malloc(sizeof(float)*n);
            float* vsr = (float*) malloc(sizeof(float)*n);
            float* vsi = (float*) malloc(sizeof(float)*n);
            float* outr_tmp = (float*) malloc(sizeof(float)*n);
            float* outi_tmp = (float*) malloc(sizeof(float)*n);
            for (int i=0; i<n; i++){
                pr[i] = cos(2 * 3.141592653 * shifts[i] * (start + offset) / slen);
                pi[i] = sin(2 * 3.141592653 * shifts[i] * (start + offset) / slen);
                vsr[i] = cos(2 * 3.141592653 * shifts[i] / slen);
                vsi[i] = sin(2 * 3.141592653 * shifts[i] / slen);
                outr_tmp[i] = 0;
                outi_tmp[i] = 0;
            }
            float t1, t2, k1, k2, k3, vs, va;
            
            for (unsigned int j=start; j<end; j++){
                std::complex<float> v = v1[j];
                float vr = v.real();
                float vi = v.imag();  
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
    """
    n = int(len(shifts))
    
    # Create some output memory
    outr =  numpy.zeros(n, dtype=numpy.float32)
    outi =  numpy.zeros(n, dtype=numpy.float32)
    
    inline(code, ['v1', 'n', 'vlen', 'outi', 'outr', 'slen', 'shifts', 'offset'],
                       extra_compile_args=['-march=native -O3 -fopenmp'],
                    libraries=['gomp'] )
    return  Array(outr + 1.0j * outi, dtype=numpy.complex64)
