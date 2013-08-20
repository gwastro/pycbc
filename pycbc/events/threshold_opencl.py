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
import pyopencl
from pycbc.types import zeros, Array
from pyopencl.array import to_device
from pyopencl.array import zeros as pzeros
from pyopencl.tools import get_or_register_dtype, dtype_to_ctype
from pyopencl.elementwise import ElementwiseKernel
from pycbc.scheme import mgr

threshold_op = """
    if (i == 0)
        bn[0] = 0;

    cfloat_t val = in[i];
    if ( cfloat_abs(val) > threshold){
        int n_w = atomic_add(bn, 1);
        outv[n_w] = val;
        outl[n_w] = i;
    }

"""

threshold_kernel = ElementwiseKernel(mgr.state.context,
            " %(tp_in)s *in, %(tp_out1)s *outv, %(tp_out2)s *outl, %(tp_th)s threshold, %(tp_n)s *bn" % {
                "tp_in": dtype_to_ctype(numpy.complex64),
                "tp_out1": dtype_to_ctype(numpy.complex64),
                "tp_out2": dtype_to_ctype(numpy.uint32),
                "tp_th": dtype_to_ctype(numpy.float32),
                "tp_n": dtype_to_ctype(numpy.uint32),
                },
            threshold_op,
            "getstuff")

n = pzeros(mgr.state.queue, 1, numpy.uint32)
val = pzeros(mgr.state.queue, 4096*256, numpy.complex64)
loc = pzeros(mgr.state.queue, 4096*256, numpy.uint32)


def threshold(series, value):
    threshold_kernel(series.data, val, loc, value, n)
    n0 = n.get()[0]
    return loc[0:n0].get(), val[0:n0].get()

