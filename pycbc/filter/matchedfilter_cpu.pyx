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
cimport numpy
import numpy
from .matchedfilter import _BaseCorrelator

ctypedef fused COMPLEXTYPE:
    float complex
    double complex

def batch_corr(numpy.ndarray [numpy.int64_t, ndim = 1] x,
               numpy.ndarray [COMPLEXTYPE, ndim = 1] y,
               numpy.ndarray [numpy.int64_t, ndim = 1] z,
               unsigned int numvec,
               unsigned int xmax):
    cdef float xr, yr, xi, yi, re, im
    cdef COMPLEXTYPE* xp
    cdef COMPLEXTYPE* zp
    for i in range(numvec):
        xp = <COMPLEXTYPE*> x[i]
        zp = <COMPLEXTYPE*> z[i]
        for j in range(xmax):
            xr = xp[j].real
            xi = xp[j].imag
            yr = y[j].real
            yi = y[j].imag
            re = xr*yr + xi*yi
            im = xr*yi - xi*yr
            zp[j] = re + 1.0j * im

def batch_correlate_execute(self, y):
    return batch_corr(self.x.data, y.data, self.z.data,
                      self.size, self.num_vectors)

def cor(numpy.ndarray [COMPLEXTYPE, ndim = 1] x,
                     numpy.ndarray [COMPLEXTYPE, ndim = 1] y,
                     numpy.ndarray [COMPLEXTYPE, ndim = 1] z):
    cdef unsigned int xmax = x.shape[0]
    cdef float xr, yr, xi, yi, re, im
    for i in range(xmax):
        xr = x[i].real
        xi = x[i].imag
        yr = y[i].real
        yi = y[i].imag

        re = xr*yr + xi*yi
        im = xr*yi - xi*yr

        z[i] = re + im*1.0j

def correlate_inline(x, y, z):
    return cor(x.data, y.data, z.data)

def correlate_numpy(x, y, z):
    z.data[:] = numpy.conjugate(x.data)[:]
    z *= y

correlate = correlate_inline

class CPUCorrelator(_BaseCorrelator):
    def __init__(self, x, y, z):
        self.x = numpy.array(x.data, copy=False)
        self.y = numpy.array(y.data, copy=False)
        self.z = numpy.array(z.data, copy=False)

    def correlate(self):
        correlate_inline(self.x, self.y, self.z)
        
def _correlate_factory(x, y, z):
    return CPUCorrelator
