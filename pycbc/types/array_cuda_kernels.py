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

from pycuda.elementwise import ElementwiseKernel


squared_norm_double = ElementwiseKernel(
        "pycuda::complex<double> *x, double *z",
        "z[i] = x[i].real() * x[i].real() + x[i].imag() * x[i].imag()",
        "squared_norm_double")
        
squared_norm_single = ElementwiseKernel(
        "pycuda::complex<float> *x, float *z",
        "z[i] = x[i].real() * x[i].real() + x[i].imag() * x[i].imag()",
        "squared_norm_single")
        
squared_norm = {'single':squared_norm_single,'double':squared_norm_double}        
        
correlate_single = ElementwiseKernel(
        "pycuda::complex<float> *x, pycuda::complex<float> *y, pycuda::complex<float> *z",
        "z[i] = conj(x[i]) * y[i]",
        "correlate_single")
        
correlate_double = ElementwiseKernel(
        "pycuda::complex<double> *x, pycuda::complex<double> *y, pycuda::complex<double> *z",
        "z[i] = conj(x[i]) * y[i]",
        "correlate_double")

correlate = {'single':correlate_single,'double':correlate_double}


        
        
