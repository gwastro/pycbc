# Copyright (C) 2018  Alex Nitz, Josh Willis
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
# cython: embedsignature=True
from __future__ import absolute_import
import numpy
from .matchedfilter import _BaseCorrelator
cimport numpy, cython

ctypedef fused COMPLEXTYPE:
    float complex
    double complex

def _batch_correlate(numpy.ndarray [long, ndim=1] x,
                     numpy.ndarray [float complex, ndim=1] y,
                     numpy.ndarray [long, ndim=1] z,
                     size, num_vectors):
    cdef unsigned int nvec = num_vectors
    cdef unsigned int vsize = size

    cdef float complex* xp
    cdef float complex* zp

    for i in range(nvec):
        xp = <float complex*> x[i]
        zp = <float complex*> z[i]
        for j in range(vsize):
            zp[j] = xp[j].conjugate() * y[j]

def batch_correlate_execute(self, y):
    num_vectors = self.num_vectors # pylint:disable=unused-variable
    size = self.size # pylint:disable=unused-variable
    _batch_correlate(self.x.data, y.data, self.z.data, size, num_vectors)

def correlate_numpy(x, y, z):
    z.data[:] = numpy.conjugate(x.data)[:]
    z *= y

@cython.boundscheck(False)
@cython.wraparound(False)
def _correlate(numpy.ndarray [COMPLEXTYPE, ndim=1] x,
               numpy.ndarray [COMPLEXTYPE, ndim=1] y,
               numpy.ndarray [COMPLEXTYPE, ndim=1] z):
    cdef unsigned int xmax = x.shape[0]
    for i in range(xmax):
        z[i] = x[i].conjugate() * y[i]

def correlate(x, y, z):
    _correlate(x.data, y.data, z.data)

class CPUCorrelator(_BaseCorrelator):
    def __init__(self, x, y, z):
        self.x = numpy.array(x.data, copy=False)
        self.y = numpy.array(y.data, copy=False)
        self.z = numpy.array(z.data, copy=False)

    def correlate(self):
        _correlate(self.x, self.y, self.z)

def _correlate_factory(x, y, z):
    return CPUCorrelator
