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
from pycbc.types import zeros, Array
from pyopencl.array import to_device
from pyopencl.tools import get_or_register_dtype, dtype_to_ctype
from pyopencl.elementwise import ElementwiseKernel
from events import complex64_subset
from pycbc.scheme import mgr

complex64_subset = get_or_register_dtype("event", dtype=complex64_subset)

preamble = """
    typedef struct{
        cfloat_t val;
        long loc;
    } event;
    
    """

threshold_op = """
    if (i == 0)
        bn[0] = 0;

    cfloat_t val = in[i];
    event nv;
    if ( cfloat_abs(val) > threshold){
        nv.val = val;
        nv.loc = i;
        int n_w = atomic_add(bn, 1);
        out[n_w] = nv;
    }

"""

threshold_kernel = ElementwiseKernel(mgr.state.context,
            " %(tp_in)s *in, %(tp_out)s *out, %(tp_th)s threshold, %(tp_n)s *bn" % {
                "tp_in": dtype_to_ctype(numpy.complex64),
                "tp_out": dtype_to_ctype(complex64_subset),
                "tp_th": dtype_to_ctype(numpy.float32),
                "tp_n": dtype_to_ctype(numpy.int32),
                },
            threshold_op,
            "getstuff", preamble=preamble)
            
n_events = numpy.zeros(1, dtype=numpy.int64)
n_events = to_device(mgr.state.queue, n_events)
buffer_vec = numpy.zeros(4096*2048, dtype=complex64_subset)
buffer_vec = to_device(mgr.state.queue, buffer_vec)
            
def threshold(series, value):
    threshold_kernel(series.data, buffer_vec, value, n_events)
    n = n_events.get()[0]
    return numpy.sort(buffer_vec[0:n].get(), order='loc')

