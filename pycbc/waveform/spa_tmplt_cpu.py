#  Copyright (C) 2013 Alex Nitz
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with with program; see the file COPYING. If not, write to the
#  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
#  MA  02111-1307  USA
import numpy
import lal
import pycbc
from pycbc.types import Array, float32, FrequencySeries
from pycbc.waveform.spa_tmplt import spa_tmplt_precondition
from scipy.weave import inline

support = """
    #include <stdio.h>
    #include <math.h>
"""

if pycbc.HAVE_OMP:
    omp_libs = ['gomp']
    omp_flags = ['-fopenmp']
else:
    omp_libs = []
    omp_flags = []

# Precompute cbrt(f) ###########################################################

def cbrt_lookup(vmax, delta):
    vec = numpy.arange(0, vmax*1.2, delta)
    return FrequencySeries(vec**(1.0/3.0), delta_f=delta).astype(float32)
    
_cbrt_vec = None
    
def get_cbrt(vmax, delta):
    global _cbrt_vec
    if _cbrt_vec is None or (_cbrt_vec.delta_f != delta) or (len(_cbrt_vec) < int(vmax/delta)):
        _cbrt_vec = cbrt_lookup(vmax, delta)
    return _cbrt_vec   
    
# Precompute log(v) ############################################################
    
def logv_lookup(vmax, delta):
    vec = numpy.arange(0, vmax*1.2, delta)
    vec[1:len(vec)] = numpy.log(vec[1:len(vec)])
    return FrequencySeries(vec, delta_f=delta).astype(float32)
    
_logv_vec = None
    
def get_log(vmax, delta):
    global _logv_vec
    if _logv_vec is None or (_logv_vec.delta_f != delta) or (len(_logv_vec) < int(vmax/delta)):
        _logv_vec = logv_lookup(vmax, delta)
    return _logv_vec   

# Precompute the sine function #################################################
def sin_cos_lookup():
    vec = numpy.arange(0, lal.TWOPI*3, lal.TWOPI/10000)
    return Array(numpy.sin(vec)).astype(float32)
sin_cos = Array([], dtype=float32)

def spa_tmplt_engine(htilde,  kmin,  phase_order, delta_f, piM,  pfaN, 
                    pfa2,  pfa3,  pfa4,  pfa5,  pfl5,
                    pfa6,  pfl6,  pfa7, amp_factor):
    """ Calculate the spa tmplt phase 
    """
    kfac = numpy.array(spa_tmplt_precondition(len(htilde), delta_f, kmin).data, copy=False)
    htilde = numpy.array(htilde.data, copy=False)
    cbrt_vec = numpy.array(get_cbrt(len(htilde)*delta_f + kmin, delta_f).data, copy=False)
    logv_vec = numpy.array(get_log(len(htilde)*delta_f + kmin, delta_f).data, copy=False)
    length = len(htilde)
    
    code = """ 
    float piM13 = cbrtf(piM);
    float logpiM13 = log(piM13);
    float log4 = log(4.);
    const float _pfaN=pfaN;
    const float _pfa2=pfa2;
    const float _pfa3=pfa3;
    const float _pfa4=pfa4;
    const float _pfa5=pfa5;
    const float _pfl5=pfl5;
    const float _pfa6=pfa6;
    const float _pfl6=pfl6;
    const float _pfa7=pfa7;
    const float ampc = amp_factor;
    const float two_pi = 2 * M_PI;
    const float inv_two_pi = 1 / (2 * M_PI);
    
    #pragma omp parallel for schedule(dynamic, 1024)
    for (unsigned int i=0; i<length; i++){
        int index = i + kmin;
        const float v =  piM13 * cbrt_vec[index];
        const float logv = logv_vec[index] * 1.0/3.0 + logpiM13;
        const float v5 = v * v * v * v * v;
        float phasing = 0;
        float sinp, cosp;

        switch (phase_order)
        {   
            case -1:
            case 7:
                phasing = _pfa7 * v;
            case 6:
                phasing = (phasing + _pfa6 + _pfl6 * (logv + log4) ) * v;
            case 5:
                phasing = (phasing + _pfa5 + _pfl5 * (logv) ) * v;
            case 4:
                phasing = (phasing + _pfa4) * v;
            case 3:
                phasing = (phasing + _pfa3) * v;
            case 2:
                phasing = (phasing + _pfa2) * v * v;
            case 0:
                phasing += 1.;
                break;
            default:
                break;
        }
        float amp = ampc * kfac[i];
        phasing *= _pfaN / v5;
        phasing -= M_PI_4;
        
        phasing -= int(phasing / two_pi) * two_pi;
         
        while (phasing < -M_PI){
            phasing += two_pi;
        }
        while (phasing > M_PI){
            phasing -= two_pi;
        }
        
        // compute sine
        if (phasing < 0)
        {
            sinp = 1.27323954 * phasing + .405284735 * phasing * phasing;
            
            if (sinp < 0)
                sinp = .225 * (sinp *-sinp - sinp) + sinp;
            else
                sinp = .225 * (sinp * sinp - sinp) + sinp;
        }
        else
        {
            sinp = 1.27323954 * phasing - 0.405284735 * phasing * phasing;          
            if (sinp < 0)
                sinp = .225 * (sinp *-sinp - sinp) + sinp;
            else
                sinp = .225 * (sinp * sinp - sinp) + sinp;
        }
        
        
        //compute cosine
        
        phasing += M_PI_2;
        if (phasing >  M_PI)
            phasing -= two_pi;
        if (phasing < 0)
        {
            cosp = 1.27323954 * phasing + .405284735 * phasing * phasing;
                     
            if (cosp < 0)
                cosp = .225 * (cosp *-cosp - cosp) + cosp;
            else
                cosp = .225 * (cosp * cosp - cosp) + cosp;
        }
        else
        {
            cosp = 1.27323954 * phasing - 0.405284735 * phasing * phasing;          
            if (cosp < 0)
                cosp = .225 * (cosp *-cosp - cosp) + cosp;
            else
                cosp = .225 * (cosp * cosp - cosp) + cosp;
        }
           
            
        
        //printf("%f %f %f \\n", sinp, sin(phasing), phasing);
        
        htilde[i] = std::complex<float>(cosp, - sinp) * amp;
    }
    """
    inline(code, ['htilde', 'cbrt_vec', 'logv_vec', 'kmin', 'phase_order', 
                   'piM',  'pfaN', 'amp_factor', 'kfac',
                   'pfa2',  'pfa3',  'pfa4',  'pfa5',  'pfl5',
                   'pfa6',  'pfl6',  'pfa7', 'length'],
                    extra_compile_args=['-march=native -O3 -w'] + omp_flags,
                    support_code = support,
                    libraries=omp_libs
                )
