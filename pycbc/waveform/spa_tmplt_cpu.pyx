#  Copyright (C) 2013 Alex Nitz
#
#  This program is free software you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with with program see the file COPYING. If not, write to the
#  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
#  MA  02111-1307  USA
cimport numpy, cython
import numpy
import lal
import pycbc
from pycbc.types import Array, float32, FrequencySeries
from pycbc.waveform.spa_tmplt import spa_tmplt_precondition
from libc.math cimport cbrt, log, M_PI, M_PI_2, M_PI_4, floor, fabs

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
    
@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
cdef spa_tmplt_inline(float piM, float pfaN,
                      float pfa2, float pfa3,
                      float pfa4, float pfa5,
                      float pfl5, float pfa6,
                      float pfl6, float pfa7,
                      float ampc, int kmin,
                      numpy.ndarray[numpy.float32_t, ndim=1] _logv_vec,
                      numpy.ndarray[numpy.float32_t, ndim=1] _cbrt_vec,
                      numpy.ndarray[numpy.float32_t, ndim=1] _kfac,
                      numpy.ndarray[numpy.complex64_t, ndim=1] _htilde,
                      ):
    cdef float piM13 = cbrt(piM)
    cdef float logpiM13 = log(piM13)
    cdef float log4 = log(4.)
    cdef float two_pi = 2 * M_PI
    cdef float v, logv, v5, phasing, amp
    
    cdef float complex* htilde = &_htilde[0]
    cdef float* kfac = &_kfac[0]
    cdef float* cbrt_vec = &_cbrt_vec[kmin]
    cdef float* logv_vec = &_logv_vec[kmin]
    cdef unsigned int i, xmax = _htilde.shape[0]

    for i in range(xmax):
        v = piM13 * cbrt_vec[i]
        logv = logv_vec[i] * 1.0/3.0 + logpiM13
        amp = ampc * kfac[i]
        v5 = v * v * v * v * v
          
        phasing = pfa7 * v
        phasing = (phasing + pfa6 + pfl6 * (logv + log4) ) * v
        phasing = (phasing + pfa5 + pfl5 * logv) * v
        phasing = (phasing + pfa4) * v
        phasing = (phasing + pfa3) * v
        phasing = (phasing + pfa2) * v * v + 1

        phasing = phasing * pfaN / v5 - M_PI_4
        phasing -= <int>(phasing / two_pi) * two_pi
        
        if (phasing < -M_PI):
            phasing += two_pi
        if (phasing > M_PI):
            phasing -= two_pi
       
        sinp = 1.273239545 * phasing - .405284735 * phasing * fabs(phasing)
        sinp = .225 * (sinp * fabs(sinp) - sinp) + sinp  
        
        phasing += M_PI_2
        if phasing > M_PI:
            phasing -= two_pi

        cosp = 1.273239545 * phasing - .405284735 * phasing * fabs(phasing)
        cosp = .225 * (cosp * fabs(cosp) - cosp) + cosp         

        htilde[i] = (cosp - sinp * 1j) * amp

def spa_tmplt_engine(htilde,  kmin,  phase_order, delta_f, piM,  pfaN, 
                    pfa2,  pfa3,  pfa4,  pfa5,  pfl5,
                    pfa6,  pfl6,  pfa7, amp_factor):
    """ Calculate the spa tmplt phase 
    """
    kfac = spa_tmplt_precondition(len(htilde), delta_f, kmin).data
    cbrt_vec = get_cbrt(len(htilde)*delta_f + kmin, delta_f).data
    logv_vec = get_log(len(htilde)*delta_f + kmin, delta_f).data
    spa_tmplt_inline(piM, pfaN, pfa2, pfa3, pfa4, pfa5, pfl5,
                      pfa6, pfl6, pfa7, amp_factor,
                      kmin, logv_vec, cbrt_vec, kfac, htilde.data,
                      )
