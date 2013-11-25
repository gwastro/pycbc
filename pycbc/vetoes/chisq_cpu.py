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

def chisq_accum_bin(chisq, q):
    chisq += q.squared_norm()

def shift_sum(v1, shifts, slen=None, offset=0):
    from scipy.weave import inline
    v1 = v1.data
    shifts = numpy.array(shifts, dtype=numpy.float32)
    vlen = len(v1)
    if slen is None:
        slen = vlen
        
    code1 = """
        float t1, t2;
        for (int j=0; j<vlen; j++){
            std::complex<float> v = v1[j];
            float vr = v.real();
            float vi = v.imag();  
                       
            for (int i=0; i<n; i++){
                outr[i] += vr * pr[i] - vi * pi[i];
                outi[i] += vr * pi[i] + vi * pr[i];
                t1 = pr[i];
                t2 = pi[i];
                pr[i] = t1 * vsr[i] - t2 * vsi[i];
                pi[i] = t1 * vsi[i] + t2 * vsr[i]; 
            }                                              
        }            
    """
    code = """
        
        float t1, t2, k1, k2, k3, vs, va;
        for (int j=0; j<vlen; j++){
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
                            
                outr[i] += k1 - k3;
                outi[i] += k1 + k2;
                
                // phase shift for the next time point
                pr[i] = t1 * vsr[i] - t2 * vsi[i];
                pi[i] = t1 * vsi[i] + t2 * vsr[i]; 
            }                                              
        }            
    """
    n = int(len(shifts))
    
    #Calculate the incremental rotation for each time shift
    vs = numpy.exp(numpy.pi * 2j * shifts / slen )
    vsr = vs.real*1
    vsi = vs.imag*1
    
    # Create some output memory
    outr =  numpy.zeros(n, dtype=numpy.float32)
    outi =  numpy.zeros(n, dtype=numpy.float32)
    
    # Create memory for storing the cumulative rotation for each time shift
    p = numpy.exp(numpy.pi * 2j *  offset * shifts / slen)
    pi = numpy.zeros(n, dtype=numpy.float32) + p.imag
    pr = numpy.zeros(n, dtype=numpy.float32) + p.real
    inline(code, ['v1', 'n', 'vlen', 'pr', 'pi', 'outi', 'outr', 'vsr', 'vsi'] )
    return  Array(outr + 1.0j * outi, dtype=numpy.complex64)
