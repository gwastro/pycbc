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
from scipy.weave import inline

support = """
    #include <stdio.h>
    #include <omp.h>
    #include <math.h>
"""

def correlate_numpy(x, y, z):
    z.data[:] = numpy.conjugate(x.data)[:]
    z *= y

def correlate_inline(x, y, z):
    za = numpy.array(z.data, copy=False)
    xa = numpy.array(x.data, copy=False)
    ya = numpy.array(y.data, copy=False)
    N = len(x)
    code = """
        #pragma omp parallel for
        for (int i=0; i<N; i++){
            float re, im;
            re = xa[i].real() * ya[i].real() + xa[i].imag() * ya[i].imag();
            im = xa[i].real() * ya[i].imag() - xa[i].imag() * ya[i].real();
            za[i] = std::complex<float>(re, im);
        }
    """
    inline(code, ['xa', 'ya', 'za', 'N'], 
           extra_compile_args=['-march=native  -O3  -fopenmp' ],
           support_code = support,
           libraries=['gomp']
          )
    
correlate = correlate_inline
    
