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
import events
from scipy.weave import inline

def threshold_numpy(series, value):
    arr = series.data
    locs = numpy.where(arr.real**2 + arr.imag**2 > value**2)[0]
    vals = arr[locs]
    return locs, vals

outl = outv = count = None
def threshold_inline(series, value):
    arr = series.data.view(dtype=numpy.float32)
    global outl, outv, count
    if outl is None:
        outl = numpy.zeros(len(series), dtype=numpy.uint32)
        outv = numpy.zeros(len(series), dtype=series.dtype)
        count = numpy.zeros(1, dtype=numpy.uint32)
        
    N = len(series)
    threshold = value**2
    code = """
        unsigned int c = 0;
        float r, im;
        for (int i=0; i<N; i++){
            r = arr[i*2];
            im = arr[i*2+1];
            if ((r*r+im*im) > threshold){
                outl[c] = i;
                outv[c] = std::complex<float>(r, im);
                c++;
            }
        }
        count[0] = c;          
    """
    inline(code, ['N', 'arr', 'outv', 'outl', 'count', 'threshold'])
    return outl[0:count], outv[0:count]

threshold=threshold_numpy
