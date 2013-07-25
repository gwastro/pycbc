# Copyright (C) 2012  Josh Willis, Andrew Miller
#
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
"""
These are the unittests for the pycbc.fft subpackage
"""

import pycbc
import pycbc.scheme
import pycbc.types
import numpy
from numpy import dtype
import pycbc.fft
import unittest
import base_test
import sys
from utils import parse_args_all_schemes

_scheme, _context = parse_args_all_schemes("FFT")

# We will need these imported in order to check that things are in the current scheme
if _scheme == 'cuda':
    import pycuda
    import pycuda.gpuarray
elif _scheme == 'opencl':
    import pyopencl
    import pyopencl.array


class _BaseTestFFTClass(base_test.function_base):
    """
    This is the base class from which unit tests for all FFT backends
    are derived.
    """
    def setUp(self):
        # Number of decimal places to compare for single precision
        self.splaces = 6
        self.smsg = 'FFT output differs by more than {0} digits from expected'.format(self.splaces)
        # Number of decimal places to compare for double precision
        self.dplaces = 14
        self.dmsg = 'FFT output differs by more than {0} digits from expected'.format(self.dplaces)
        # Message if overwrote input
        self.omsg = 'FFT overwrote input array'
        self.scheme = _scheme
        self.context = _context

    def test_forward_real_single(self):
        # First, test case when input array length is even
        self.in_even = pycbc.types.Array([1.0,-1.0,2.0,-2.0],
                                         dtype=dtype('float32'))
        self.in_pristine = pycbc.types.Array([1.0,-1.0,2.0,-2.0],
                                             dtype=dtype('float32'))
        self.out_even = pycbc.types.Array([0.0+0.0j,-1.0-1.0j,6.0+0.0j],
                                          dtype=dtype('complex64'))
        self.out_even_test = pycbc.types.zeros(3,dtype=dtype('complex64'))

        if self.scheme != 'cpu':
            self.scheme_test(pycbc.fft.fft,(self.in_even,self.out_even_test),
                                            (self.in_pristine,self.out_even),self.splaces,backend=self.backend)
            for back in pycbc.fft.cpu_backends:
                self.cpu_test(pycbc.fft.fft,(self.in_even,self.out_even_test),
                                            (self.in_pristine,self.out_even),self.splaces,backend=back)
        with self.context:
            pycbc.fft.fft(self.in_even,self.out_even_test,backend=self.backend)
            # First, check that we have not overwritten the input array
            self.assertEqual(self.in_even[0],self.in_pristine[0],msg=self.omsg)
            self.assertEqual(self.in_even[1],self.in_pristine[1],msg=self.omsg)
            self.assertEqual(self.in_even[2],self.in_pristine[2],msg=self.omsg)
            self.assertEqual(self.in_even[3],self.in_pristine[3],msg=self.omsg)
            # Check that output is correct. Note that we compare most
            # entries to be AlmostEqual, but the imaginary parts of DC and
            # Nyquist to be exactly equal to zero.
            self.assertAlmostEqual(self.out_even[0].real,
                              self.out_even_test[0].real,
                              places=self.splaces,msg=self.smsg)
            self.assertEqual(self.out_even[0].imag,self.out_even_test[0].imag,
                        msg="Imaginary part of DC was not exactly zero")
            self.assertAlmostEqual(self.out_even[1].real,
                              self.out_even_test[1].real,
                              places=self.splaces,msg=self.smsg)
            self.assertAlmostEqual(self.out_even[1].imag,
                              self.out_even_test[1].imag,
                              places=self.splaces,msg=self.smsg)
            self.assertAlmostEqual(self.out_even[2].real,
                              self.out_even_test[2].real,
                              places=self.splaces,msg=self.smsg)
            self.assertEqual(self.out_even[2].imag,self.out_even_test[2].imag,
                        msg="Imaginary part of Nyquist was not exactly zero")
        # Now, another test case when input array length is odd
        self.in_odd = pycbc.types.Array([1.0,2.0,2.0],
                                        dtype=dtype('float32'))
        self.out_odd = pycbc.types.Array([5.0+0.0j,-1.0+0.0j],
                                         dtype=dtype('complex64'))
        self.out_odd_test = pycbc.types.zeros(2,dtype=dtype('complex64'))
        if self.scheme != 'cpu':
            self.scheme_test(pycbc.fft.fft,(self.in_odd,self.out_odd_test), 
                                            (self.in_odd,self.out_odd),self.splaces,backend=self.backend)
            for back in pycbc.fft.cpu_backends:
                self.cpu_test(pycbc.fft.fft,(self.in_odd,self.out_odd_test), 
                                            (self.in_odd,self.out_odd),self.splaces,backend=back)
        with self.context:
            pycbc.fft.fft(self.in_odd,self.out_odd_test,backend=self.backend)
            # Compare again.  Now only imaginary part of DC is strictly compared
            # with zero.
            self.assertAlmostEqual(self.out_odd[0].real,
                              self.out_odd_test[0].real,
                              places=self.splaces,msg=self.smsg)
            self.assertEqual(self.out_odd[0].imag,self.out_odd_test[0].imag,
                        msg="Imaginary part of DC was not exactly zero")
            self.assertAlmostEqual(self.out_odd[1].real,
                              self.out_odd_test[1].real,
                              places=self.splaces,msg=self.smsg)
            self.assertAlmostEqual(self.out_odd[1].imag,
                              self.out_odd_test[1].imag,
                              places=self.splaces,msg=self.smsg)
            # Now test that the proper exceptions are raised when we give
            # erroneous arguments
            self.out_badlen = pycbc.types.zeros(3,dtype=dtype('complex64'))
            args = [self.in_odd,self.out_badlen,self.backend]
            self.assertRaises(ValueError,pycbc.fft.fft,*args)
            self.out_badprec = pycbc.types.zeros(2,dtype=dtype('complex128'))
            args = [self.in_odd,self.out_badprec,self.backend]
            self.assertRaises(ValueError,pycbc.fft.fft,*args)
            self.out_baddtype = pycbc.types.zeros(2,dtype=dtype('float32'))
            args = [self.in_odd,self.out_baddtype,self.backend]
            self.assertRaises(ValueError,pycbc.fft.fft,*args)
            self.out_badarray = numpy.zeros(2,dtype=dtype('complex64'))
            args = [self.in_odd,self.out_badarray,self.backend]
            self.assertRaises(TypeError,pycbc.fft.fft,*args)
            self.in_badarray = numpy.zeros(3,dtype=dtype('float32'))
            args = [self.in_badarray,self.out_odd_test,self.backend]
            self.assertRaises(TypeError,pycbc.fft.fft,*args)


    def test_forward_real_double(self):
        # First, test case when input array length is even
        self.in_even = pycbc.types.Array([1.0,-1.0,2.0,-2.0],
                                         dtype=dtype('float64'))
        self.in_pristine = pycbc.types.Array([1.0,-1.0,2.0,-2.0],
                                             dtype=dtype('float64'))
        self.out_even = pycbc.types.Array([0.0+0.0j,-1.0-1.0j,6.0+0.0j],
                                          dtype=dtype('complex128'))
        self.out_even_test = pycbc.types.zeros(3,dtype=dtype('complex128'))
        if self.scheme != 'cpu':
            self.scheme_test(pycbc.fft.fft,(self.in_even,self.out_even_test),
                                            (self.in_pristine,self.out_even),self.dplaces,backend=self.backend)
            for back in pycbc.fft.cpu_backends:
                self.cpu_test(pycbc.fft.fft,(self.in_even,self.out_even_test),
                                            (self.in_pristine,self.out_even),self.dplaces,backend=back)
        with self.context:
            pycbc.fft.fft(self.in_even,self.out_even_test,backend=self.backend)
            # First, check that we have not overwritten the input array
            self.assertEqual(self.in_even[0],self.in_pristine[0],msg=self.omsg)
            self.assertEqual(self.in_even[1],self.in_pristine[1],msg=self.omsg)
            self.assertEqual(self.in_even[2],self.in_pristine[2],msg=self.omsg)
            self.assertEqual(self.in_even[3],self.in_pristine[3],msg=self.omsg)
            # Check that output is correct. Note that we compare most
            # entries to be AlmostEqual, but the imaginary parts of DC and
            # Nyquist to be exactly equal to zero.
            self.assertAlmostEqual(self.out_even[0].real,
                              self.out_even_test[0].real,
                              places=self.dplaces,msg=self.dmsg)
            self.assertEqual(self.out_even[0].imag,self.out_even_test[0].imag,
                        msg="Imaginary part of DC was not exactly zero")
            self.assertAlmostEqual(self.out_even[1].real,
                              self.out_even_test[1].real,
                              places=self.dplaces,msg=self.dmsg)
            self.assertAlmostEqual(self.out_even[1].imag,
                              self.out_even_test[1].imag,
                              places=self.dplaces,msg=self.dmsg)
            self.assertAlmostEqual(self.out_even[2].real,
                              self.out_even_test[2].real,
                              places=self.dplaces,msg=self.dmsg)
            self.assertEqual(self.out_even[2].imag,self.out_even_test[2].imag,
                        msg="Imaginary part of Nyquist was not exactly zero")
        # Now, another test case when input array length is odd
        self.in_odd = pycbc.types.Array([1.0,2.0,2.0],
                                        dtype=dtype('float64'))
        self.out_odd = pycbc.types.Array([5.0+0.0j,-1.0+0.0j],
                                         dtype=dtype('complex128'))
        self.out_odd_test = pycbc.types.zeros(2,dtype=dtype('complex128'))
        if self.scheme != 'cpu':
            self.scheme_test(pycbc.fft.fft,(self.in_odd,self.out_odd_test),
                                            (self.in_odd,self.out_odd),self.dplaces,backend=self.backend)
            for back in pycbc.fft.cpu_backends:
                self.cpu_test(pycbc.fft.fft,(self.in_odd,self.out_odd_test),
                                            (self.in_odd,self.out_odd),self.dplaces,backend=back)
        with self.context:
            pycbc.fft.fft(self.in_odd,self.out_odd_test,backend=self.backend)
            # Compare again.  Now only imaginary part of DC is strictly compared
            # with zero.
            self.assertAlmostEqual(self.out_odd[0].real,
                              self.out_odd_test[0].real,
                              places=self.dplaces,msg=self.dmsg)
            self.assertEqual(self.out_odd[0].imag,self.out_odd_test[0].imag,
                        msg="Imaginary part of DC was not exactly zero")
            self.assertAlmostEqual(self.out_odd[1].real,
                              self.out_odd_test[1].real,
                              places=self.dplaces,msg=self.dmsg)
            self.assertAlmostEqual(self.out_odd[1].imag,
                              self.out_odd_test[1].imag,
                              places=self.dplaces,msg=self.dmsg)
            # Now test that the proper exceptions are raised when we give
            # erroneous arguments
            self.out_badlen = pycbc.types.zeros(3,dtype=dtype('complex128'))
            args = [self.in_odd,self.out_badlen,self.backend]
            self.assertRaises(ValueError,pycbc.fft.fft,*args)
            self.out_badprec = pycbc.types.zeros(2,dtype=dtype('complex64'))
            args = [self.in_odd,self.out_badprec,self.backend]
            self.assertRaises(ValueError,pycbc.fft.fft,*args)
            self.out_baddtype = pycbc.types.zeros(2,dtype=dtype('float64'))
            args = [self.in_odd,self.out_baddtype,self.backend]
            self.assertRaises(ValueError,pycbc.fft.fft,*args)
            self.out_badarray = numpy.zeros(2,dtype=dtype('complex128'))
            args = [self.in_odd,self.out_badarray,self.backend]
            self.assertRaises(TypeError,pycbc.fft.fft,*args)
            self.in_badarray = numpy.zeros(3,dtype=dtype('float64'))
            args = [self.in_badarray,self.out_odd_test,self.backend]
            self.assertRaises(TypeError,pycbc.fft.fft,*args)


    def test_inverse_real_single(self):
        # First, test case when output array length is even
        self.in_even = pycbc.types.Array([0.0+0.0j,-1.0-1.0j,6.0+0.0j],
                                          dtype=dtype('complex64'))
        self.in_pristine = pycbc.types.Array([0.0+0.0j,-1.0-1.0j,6.0+0.0j],
                                          dtype=dtype('complex64'))
        self.out_even = pycbc.types.Array([4.0,-4.0,8.0,-8.0],
                                         dtype=dtype('float32'))
        self.out_even_test = pycbc.types.zeros(4,dtype=dtype('float32'))
        if self.scheme != 'cpu':
            self.scheme_test(pycbc.fft.ifft,(self.in_even,self.out_even_test),
                                            (self.in_pristine,self.out_even),self.splaces,backend=self.backend)
            for back in pycbc.fft.cpu_backends:
                self.cpu_test(pycbc.fft.ifft,(self.in_even,self.out_even_test),
                                            (self.in_pristine,self.out_even),self.splaces,backend=back)
        with self.context:
            pycbc.fft.ifft(self.in_even,self.out_even_test,backend=self.backend)
            # First, check that we have not overwritten the input array
            self.assertEqual(self.in_even[0].real,self.in_pristine[0].real,msg=self.omsg)
            self.assertEqual(self.in_even[0].imag,self.in_pristine[0].imag,msg=self.omsg)
            self.assertEqual(self.in_even[1].real,self.in_pristine[1].real,msg=self.omsg)
            self.assertEqual(self.in_even[1].imag,self.in_pristine[1].imag,msg=self.omsg)
            self.assertEqual(self.in_even[2].real,self.in_pristine[2].real,msg=self.omsg)
            self.assertEqual(self.in_even[2].imag,self.in_pristine[2].imag,msg=self.omsg)
            # Check that output is correct.
            self.assertAlmostEqual(self.out_even[0],self.out_even_test[0],
                              places=self.splaces,msg=self.smsg)
            self.assertAlmostEqual(self.out_even[1],self.out_even_test[1],
                              places=self.splaces,msg=self.smsg)
            self.assertAlmostEqual(self.out_even[2],self.out_even_test[2],
                              places=self.splaces,msg=self.smsg)
            self.assertAlmostEqual(self.out_even[3],self.out_even_test[3],
                              places=self.splaces,msg=self.smsg)

        # Now, another test case when output array length is odd
        self.in_odd = pycbc.types.Array([5.0+0.0j,-1.0+0.0j],
                                         dtype=dtype('complex64'))
        self.out_odd = pycbc.types.Array([3.0,6.0,6.0],
                                        dtype=dtype('float32'))
        self.out_odd_test = pycbc.types.zeros(3,dtype=dtype('float32'))
        if self.scheme != 'cpu':
            self.scheme_test(pycbc.fft.ifft,(self.in_odd,self.out_odd_test),
                                            (self.in_odd,self.out_odd),self.splaces,backend=self.backend)
            for back in pycbc.fft.cpu_backends:
                self.cpu_test(pycbc.fft.ifft,(self.in_odd,self.out_odd_test),
                                            (self.in_odd,self.out_odd),self.splaces,backend=back)
        with self.context:
            pycbc.fft.ifft(self.in_odd,self.out_odd_test,backend=self.backend)
            # Compare again.
            self.assertAlmostEqual(self.out_odd[0],self.out_odd_test[0],
                              places=self.splaces,msg=self.smsg)
            self.assertAlmostEqual(self.out_odd[1],self.out_odd_test[1],
                              places=self.splaces,msg=self.smsg)
            self.assertAlmostEqual(self.out_odd[2],self.out_odd_test[2],
                              places=self.splaces,msg=self.smsg)
            # Now test that the proper exceptions are raised when we give
            # erroneous arguments
            self.out_badlen = pycbc.types.zeros(5,dtype=dtype('float32'))
            args = [self.in_odd,self.out_badlen,self.backend]
            self.assertRaises(ValueError,pycbc.fft.ifft,*args)
            self.out_badprec = pycbc.types.zeros(3,dtype=dtype('float64'))
            args = [self.in_odd,self.out_badprec,self.backend]
            self.assertRaises(ValueError,pycbc.fft.ifft,*args)
            self.in_baddtype = pycbc.types.zeros(3,dtype=dtype('float32'))
            args = [self.in_baddtype,self.out_odd,self.backend]
            self.assertRaises(ValueError,pycbc.fft.ifft,*args)
            self.out_badarray = numpy.zeros(3,dtype=dtype('float32'))
            args = [self.in_odd,self.out_badarray,self.backend]
            self.assertRaises(TypeError,pycbc.fft.ifft,*args)
            self.in_badarray = numpy.zeros(2,dtype=dtype('complex64'))
            args = [self.in_badarray,self.out_odd_test,self.backend]
            self.assertRaises(TypeError,pycbc.fft.ifft,*args)

    def test_inverse_real_double(self):
        # First, test case when output array length is even
        self.in_even = pycbc.types.Array([0.0+0.0j,-1.0-1.0j,6.0+0.0j],
                                          dtype=dtype('complex128'))
        self.in_pristine = pycbc.types.Array([0.0+0.0j,-1.0-1.0j,6.0+0.0j],
                                          dtype=dtype('complex128'))
        self.out_even = pycbc.types.Array([4.0,-4.0,8.0,-8.0],
                                         dtype=dtype('float64'))
        self.out_even_test = pycbc.types.zeros(4,dtype=dtype('float64'))
        if self.scheme != 'cpu':
            self.scheme_test(pycbc.fft.ifft,(self.in_even,self.out_even_test),
                                            (self.in_pristine,self.out_even),self.dplaces,backend=self.backend)
            for back in pycbc.fft.cpu_backends:
                self.cpu_test(pycbc.fft.ifft,(self.in_even,self.out_even_test),
                                            (self.in_pristine,self.out_even),self.dplaces,backend=back)
        with self.context:
            pycbc.fft.ifft(self.in_even,self.out_even_test,backend=self.backend)
            # First, check that we have not overwritten the input array
            self.assertEqual(self.in_even[0].real,self.in_pristine[0].real,msg=self.omsg)
            self.assertEqual(self.in_even[0].imag,self.in_pristine[0].imag,msg=self.omsg)
            self.assertEqual(self.in_even[1].real,self.in_pristine[1].real,msg=self.omsg)
            self.assertEqual(self.in_even[1].imag,self.in_pristine[1].imag,msg=self.omsg)
            self.assertEqual(self.in_even[2].real,self.in_pristine[2].real,msg=self.omsg)
            self.assertEqual(self.in_even[2].imag,self.in_pristine[2].imag,msg=self.omsg)
            # Check that output is correct.
            self.assertAlmostEqual(self.out_even[0],self.out_even_test[0],
                              places=self.dplaces,msg=self.dmsg)
            self.assertAlmostEqual(self.out_even[1],self.out_even_test[1],
                              places=self.dplaces,msg=self.dmsg)
            self.assertAlmostEqual(self.out_even[2],self.out_even_test[2],
                              places=self.dplaces,msg=self.dmsg)
            self.assertAlmostEqual(self.out_even[3],self.out_even_test[3],
                              places=self.dplaces,msg=self.dmsg)
        # Now, another test case when output array length is odd
        self.in_odd = pycbc.types.Array([5.0+0.0j,-1.0+0.0j],
                                         dtype=dtype('complex128'))
        self.out_odd = pycbc.types.Array([3.0,6.0,6.0],
                                        dtype=dtype('float64'))
        self.out_odd_test = pycbc.types.zeros(3,dtype=dtype('float64'))
        if self.scheme != 'cpu':
            self.scheme_test(pycbc.fft.ifft,(self.in_odd,self.out_odd_test),
                                            (self.in_odd,self.out_odd),self.dplaces,backend=self.backend)
            for back in pycbc.fft.cpu_backends:
                self.cpu_test(pycbc.fft.ifft,(self.in_odd,self.out_odd_test),
                                            (self.in_odd,self.out_odd),self.dplaces,backend=back)
        with self.context:
            pycbc.fft.ifft(self.in_odd,self.out_odd_test,backend=self.backend)
            # Compare again.
            self.assertAlmostEqual(self.out_odd[0],self.out_odd_test[0],
                              places=self.dplaces,msg=self.dmsg)
            self.assertAlmostEqual(self.out_odd[1],self.out_odd_test[1],
                              places=self.dplaces,msg=self.dmsg)
            self.assertAlmostEqual(self.out_odd[2],self.out_odd_test[2],
                              places=self.dplaces,msg=self.dmsg)
            # Now test that the proper exceptions are raised when we give
            # erroneous arguments
            self.out_badlen = pycbc.types.zeros(5,dtype=dtype('float64'))
            args = [self.in_odd,self.out_badlen,self.backend]
            self.assertRaises(ValueError,pycbc.fft.ifft,*args)
            self.out_badprec = pycbc.types.zeros(3,dtype=dtype('float32'))
            args = [self.in_odd,self.out_badprec,self.backend]
            self.assertRaises(ValueError,pycbc.fft.ifft,*args)
            self.in_baddtype = pycbc.types.zeros(3,dtype=dtype('float64'))
            args = [self.in_baddtype,self.out_odd,self.backend]
            self.assertRaises(ValueError,pycbc.fft.ifft,*args)
            self.out_badarray = numpy.zeros(3,dtype=dtype('float64'))
            args = [self.in_odd,self.out_badarray,self.backend]
            self.assertRaises(TypeError,pycbc.fft.ifft,*args)
            self.in_badarray = numpy.zeros(2,dtype=dtype('complex128'))
            args = [self.in_badarray,self.out_odd_test,self.backend]
            self.assertRaises(TypeError,pycbc.fft.ifft,*args)

    def test_forward_complex_single(self):
        # A forward complex test case
        self.in_cmplx = pycbc.types.Array([1.0+1.0j,2.0-2.0j],
                                          dtype=dtype('complex64'))
        self.in_pristine = pycbc.types.Array([1.0+1.0j,2.0-2.0j],
                                             dtype=dtype('complex64'))
        self.out_cmplx = pycbc.types.Array([3.0-1.0j,-1.0+3.0j],
                                           dtype=dtype('complex64'))
        self.out_cmplx_test = pycbc.types.zeros(2,dtype=dtype('complex64'))
        if self.scheme != 'cpu':
            self.scheme_test(pycbc.fft.fft,(self.in_cmplx,self.out_cmplx_test),
                                            (self.in_pristine,self.out_cmplx),self.splaces,backend=self.backend)
            for back in pycbc.fft.cpu_backends:
                self.cpu_test(pycbc.fft.fft,(self.in_cmplx,self.out_cmplx_test),
                                            (self.in_pristine,self.out_cmplx),self.splaces,backend=back)
        with self.context:
            pycbc.fft.fft(self.in_cmplx,self.out_cmplx_test,backend=self.backend)
            # First, check that we have not overwritten the input array
            self.assertEqual(self.in_cmplx[0].real,self.in_pristine[0].real,msg=self.omsg)
            self.assertEqual(self.in_cmplx[0].imag,self.in_pristine[0].imag,msg=self.omsg)
            self.assertEqual(self.in_cmplx[1].real,self.in_pristine[1].real,msg=self.omsg)
            self.assertEqual(self.in_cmplx[1].imag,self.in_pristine[1].imag,msg=self.omsg)
            # Check that output is correct.
            self.assertAlmostEqual(self.out_cmplx[0].real,self.out_cmplx_test[0].real,
                              places=self.splaces,msg=self.smsg)
            self.assertAlmostEqual(self.out_cmplx[0].imag,self.out_cmplx_test[0].imag,
                              places=self.splaces,msg=self.smsg)
            self.assertAlmostEqual(self.out_cmplx[1].real,self.out_cmplx_test[1].real,
                              places=self.splaces,msg=self.smsg)
            self.assertAlmostEqual(self.out_cmplx[1].imag,self.out_cmplx_test[1].imag,
                              places=self.splaces,msg=self.smsg)
            # Now test that the proper exceptions are raised when we give
            # erroneous arguments
            self.out_badlen = pycbc.types.zeros(3,dtype=dtype('complex64'))
            args = [self.in_cmplx,self.out_badlen,self.backend]
            self.assertRaises(ValueError,pycbc.fft.fft,*args)
            self.out_badprec = pycbc.types.zeros(2,dtype=dtype('complex128'))
            args = [self.in_cmplx,self.out_badprec,self.backend]
            self.assertRaises(ValueError,pycbc.fft.fft,*args)
            self.out_badarray = numpy.zeros(2,dtype=dtype('complex64'))
            args = [self.in_cmplx,self.out_badarray,self.backend]
            self.assertRaises(TypeError,pycbc.fft.fft,*args)
            self.in_badarray = numpy.zeros(2,dtype=dtype('complex64'))
            args = [self.in_badarray,self.out_cmplx_test,self.backend]
            self.assertRaises(TypeError,pycbc.fft.fft,*args)

    def test_forward_complex_double(self):
        # A forward complex test case
        self.in_cmplx = pycbc.types.Array([1.0+1.0j,2.0-2.0j],
                                          dtype=dtype('complex128'))
        self.in_pristine = pycbc.types.Array([1.0+1.0j,2.0-2.0j],
                                             dtype=dtype('complex128'))
        self.out_cmplx = pycbc.types.Array([3.0-1.0j,-1.0+3.0j],
                                           dtype=dtype('complex128'))
        self.out_cmplx_test = pycbc.types.zeros(2,dtype=dtype('complex128'))
        if self.scheme != 'cpu':
            self.scheme_test(pycbc.fft.fft,(self.in_cmplx,self.out_cmplx_test),
                                            (self.in_pristine,self.out_cmplx),self.dplaces,backend=self.backend)
            for back in pycbc.fft.cpu_backends:
                self.cpu_test(pycbc.fft.fft,(self.in_cmplx,self.out_cmplx_test),
                                            (self.in_pristine,self.out_cmplx),self.dplaces,backend=back)
        with self.context:
            pycbc.fft.fft(self.in_cmplx,self.out_cmplx_test,backend=self.backend)
            # First, check that we have not overwritten the input array
            self.assertEqual(self.in_cmplx[0].real,self.in_pristine[0].real,msg=self.omsg)
            self.assertEqual(self.in_cmplx[0].imag,self.in_pristine[0].imag,msg=self.omsg)
            self.assertEqual(self.in_cmplx[1].real,self.in_pristine[1].real,msg=self.omsg)
            self.assertEqual(self.in_cmplx[1].imag,self.in_pristine[1].imag,msg=self.omsg)
            # Check that output is correct.
            self.assertAlmostEqual(self.out_cmplx[0].real,self.out_cmplx_test[0].real,
                              places=self.dplaces,msg=self.dmsg)
            self.assertAlmostEqual(self.out_cmplx[0].imag,self.out_cmplx_test[0].imag,
                              places=self.dplaces,msg=self.dmsg)
            self.assertAlmostEqual(self.out_cmplx[1].real,self.out_cmplx_test[1].real,
                              places=self.dplaces,msg=self.dmsg)
            self.assertAlmostEqual(self.out_cmplx[1].imag,self.out_cmplx_test[1].imag,
                              places=self.dplaces,msg=self.dmsg)
            # Now test that the proper exceptions are raised when we give
            # erroneous arguments
            self.out_badlen = pycbc.types.zeros(3,dtype=dtype('complex128'))
            args = [self.in_cmplx,self.out_badlen,self.backend]
            self.assertRaises(ValueError,pycbc.fft.fft,*args)
            self.out_badprec = pycbc.types.zeros(2,dtype=dtype('complex64'))
            args = [self.in_cmplx,self.out_badprec,self.backend]
            self.assertRaises(ValueError,pycbc.fft.fft,*args)
            self.out_badarray = numpy.zeros(2,dtype=dtype('complex128'))
            args = [self.in_cmplx,self.out_badarray,self.backend]
            self.assertRaises(TypeError,pycbc.fft.fft,*args)
            self.in_badarray = numpy.zeros(2,dtype=dtype('complex128'))
            args = [self.in_badarray,self.out_cmplx_test,self.backend]
            self.assertRaises(TypeError,pycbc.fft.fft,*args)

    def test_inverse_complex_single(self):
        # A reverse complex test case
        self.in_cmplx = pycbc.types.Array([3.0-1.0j,-1.0+3.0j],
                                         dtype=dtype('complex64'))
        self.in_pristine = pycbc.types.Array([3.0-1.0j,-1.0+3.0j],
                                         dtype=dtype('complex64'))
        self.out_cmplx = pycbc.types.Array([2.0+2.0j,4.0-4.0j],
                                          dtype=dtype('complex64'))
        self.out_cmplx_test = pycbc.types.zeros(2,dtype=dtype('complex64'))
        if self.scheme != 'cpu':
            self.scheme_test(pycbc.fft.ifft,(self.in_cmplx,self.out_cmplx_test),
                                            (self.in_pristine,self.out_cmplx),self.splaces,backend=self.backend)
            for back in pycbc.fft.cpu_backends:
                self.cpu_test(pycbc.fft.ifft,(self.in_cmplx,self.out_cmplx_test),
                                            (self.in_pristine,self.out_cmplx),self.splaces,backend=back)
        with self.context:
            pycbc.fft.ifft(self.in_cmplx,self.out_cmplx_test,backend=self.backend)
            # First, check that we have not overwritten the input array
            self.assertEqual(self.in_cmplx[0].real,self.in_pristine[0].real,msg=self.omsg)
            self.assertEqual(self.in_cmplx[0].imag,self.in_pristine[0].imag,msg=self.omsg)
            self.assertEqual(self.in_cmplx[1].real,self.in_pristine[1].real,msg=self.omsg)
            self.assertEqual(self.in_cmplx[1].imag,self.in_pristine[1].imag,msg=self.omsg)
            # Check that output is correct.
            self.assertAlmostEqual(self.out_cmplx[0].real,self.out_cmplx_test[0].real,
                              places=self.splaces,msg=self.smsg)
            self.assertAlmostEqual(self.out_cmplx[0].imag,self.out_cmplx_test[0].imag,
                              places=self.splaces,msg=self.smsg)
            self.assertAlmostEqual(self.out_cmplx[1].real,self.out_cmplx_test[1].real,
                              places=self.splaces,msg=self.smsg)
            self.assertAlmostEqual(self.out_cmplx[1].imag,self.out_cmplx_test[1].imag,
                              places=self.splaces,msg=self.smsg)
            # Now test that the proper exceptions are raised when we give
            # erroneous arguments
            self.out_badlen = pycbc.types.zeros(3,dtype=dtype('complex64'))
            args = [self.in_cmplx,self.out_badlen,self.backend]
            self.assertRaises(ValueError,pycbc.fft.ifft,*args)
            self.out_badprec = pycbc.types.zeros(2,dtype=dtype('complex128'))
            args = [self.in_cmplx,self.out_badprec,self.backend]
            self.assertRaises(ValueError,pycbc.fft.ifft,*args)
            self.out_badarray = numpy.zeros(2,dtype=dtype('complex64'))
            args = [self.in_cmplx,self.out_badarray,self.backend]
            self.assertRaises(TypeError,pycbc.fft.ifft,*args)
            self.in_badarray = numpy.zeros(2,dtype=dtype('complex64'))
            args = [self.in_badarray,self.out_cmplx_test,self.backend]
            self.assertRaises(TypeError,pycbc.fft.ifft,*args)

    def test_inverse_complex_double(self):
        # A reverse complex test case
        self.in_cmplx = pycbc.types.Array([3.0-1.0j,-1.0+3.0j],
                                         dtype=dtype('complex128'))
        self.in_pristine = pycbc.types.Array([3.0-1.0j,-1.0+3.0j],
                                         dtype=dtype('complex128'))
        self.out_cmplx = pycbc.types.Array([2.0+2.0j,4.0-4.0j],
                                          dtype=dtype('complex128'))
        self.out_cmplx_test = pycbc.types.zeros(2,dtype=dtype('complex128'))
        if self.scheme != 'cpu':
            self.scheme_test(pycbc.fft.ifft,(self.in_cmplx,self.out_cmplx_test),
                                            (self.in_pristine,self.out_cmplx),self.dplaces,backend=self.backend)
            for back in pycbc.fft.cpu_backends:
                self.cpu_test(pycbc.fft.ifft,(self.in_cmplx,self.out_cmplx_test),
                                            (self.in_pristine,self.out_cmplx),self.dplaces,backend=back)
        with self.context:
            pycbc.fft.ifft(self.in_cmplx,self.out_cmplx_test,backend=self.backend)
            # First, check that we have not overwritten the input array
            self.assertEqual(self.in_cmplx[0].real,self.in_pristine[0].real,msg=self.omsg)
            self.assertEqual(self.in_cmplx[0].imag,self.in_pristine[0].imag,msg=self.omsg)
            self.assertEqual(self.in_cmplx[1].real,self.in_pristine[1].real,msg=self.omsg)
            self.assertEqual(self.in_cmplx[1].imag,self.in_pristine[1].imag,msg=self.omsg)
            # Check that output is correct.
            self.assertAlmostEqual(self.out_cmplx[0].real,self.out_cmplx_test[0].real,
                              places=self.dplaces,msg=self.dmsg)
            self.assertAlmostEqual(self.out_cmplx[0].imag,self.out_cmplx_test[0].imag,
                              places=self.dplaces,msg=self.dmsg)
            self.assertAlmostEqual(self.out_cmplx[1].real,self.out_cmplx_test[1].real,
                              places=self.dplaces,msg=self.dmsg)
            self.assertAlmostEqual(self.out_cmplx[1].imag,self.out_cmplx_test[1].imag,
                              places=self.dplaces,msg=self.dmsg)
            # Now test that the proper exceptions are raised when we give
            # erroneous arguments
            self.out_badlen = pycbc.types.zeros(3,dtype=dtype('complex128'))
            args = [self.in_cmplx,self.out_badlen,self.backend]
            self.assertRaises(ValueError,pycbc.fft.ifft,*args)
            self.out_badprec = pycbc.types.zeros(2,dtype=dtype('complex64'))
            args = [self.in_cmplx,self.out_badprec,self.backend]
            self.assertRaises(ValueError,pycbc.fft.ifft,*args)
            self.out_badarray = numpy.zeros(2,dtype=dtype('complex128'))
            args = [self.in_cmplx,self.out_badarray,self.backend]
            self.assertRaises(TypeError,pycbc.fft.ifft,*args)
            self.in_badarray = numpy.zeros(2,dtype=dtype('complex128'))
            args = [self.in_badarray,self.out_cmplx_test,self.backend]
            self.assertRaises(TypeError,pycbc.fft.ifft,*args)

    def test_time_frequency(self):
        with self.context:
            # In these tests we look only at exceptions, metadata, and scaling
            self.in_ts = pycbc.types.TimeSeries([1.0,2.0,3.0],
                                dtype=dtype('float32'),delta_t=1.0/4096.0)
            self.out_ts = pycbc.types.TimeSeries([1.0,2.0,3.0],
                                dtype=dtype('float32'),delta_t=1.0/4096.0)
            self.fs_deltaf = 4096.0/3.0
            self.out_ts_test = pycbc.types.TimeSeries([0.0,0.0,0.0],
                                dtype=dtype('float32'),delta_t=0.01)
            self.out_fs_test = pycbc.types.FrequencySeries([0.0,0.0],
                                dtype=dtype('complex64'),delta_f = 0.01)
            # First, check for appropriate exceptions to be raised by forward fft:
            self.out_ts_badtype = pycbc.types.TimeSeries([0.0,0.0],
                                   dtype=dtype('complex64'),delta_t=0.01)
            args = [self.in_ts,self.out_ts_badtype,self.backend]
            self.assertRaises(TypeError,pycbc.fft.fft,*args)
            self.out_fs_badtype = pycbc.types.FrequencySeries([0.0,0.0],
                                   dtype=dtype('complex64'),delta_f=0.01)
            args = [self.out_fs_test,self.out_fs_badtype,self.backend]
            self.assertRaises(TypeError,pycbc.fft.fft,*args)

            self.out_badarray = pycbc.types.zeros(2,dtype=dtype('complex64'))
            args = [self.in_ts,self.out_badarray,self.backend]
            self.assertRaises(TypeError,pycbc.fft.fft,*args)
            args = [self.out_fs_test,self.out_badarray,self.backend]
            self.assertRaises(TypeError,pycbc.fft.fft,*args)

            # Next, check for appropriate exceptions to be raised by inverse fft:
            self.in_ts_badtype = pycbc.types.TimeSeries([0.0,0.0],
                                   dtype=dtype('complex64'),delta_t=0.01)
            args = [self.in_ts_badtype,self.out_ts_test,self.backend]
            self.assertRaises(TypeError,pycbc.fft.ifft,*args)
            self.in_fs_badtype = pycbc.types.FrequencySeries([0.0,0.0],
                                   dtype=dtype('complex64'),delta_f=0.01)
            args = [self.in_fs_badtype,self.out_fs_test,self.backend]
            self.assertRaises(TypeError,pycbc.fft.ifft,*args)

            self.in_badarray = pycbc.types.zeros(2,dtype=dtype('complex64'))
            args = [self.in_badarray,self.out_ts_test,self.backend]
            self.assertRaises(TypeError,pycbc.fft.ifft,*args)
            args = [self.in_badarray,self.out_fs_test,self.backend]
            self.assertRaises(TypeError,pycbc.fft.ifft,*args)

            # Finally, check that we get the correct scaled values.  Since we're
            # testing both forward and backward, we check that ifft(fft(input))
            # is input, and that the intermediate delta_f is set correctly.
            pycbc.fft.fft(self.in_ts,self.out_fs_test,backend=self.backend)
            self.assertAlmostEqual(self.out_fs_test._delta_f,self.fs_deltaf,
                                   places=self.splaces,msg=self.smsg)
            pycbc.fft.ifft(self.out_fs_test,self.out_ts_test,backend=self.backend)
            self.assertAlmostEqual(self.out_ts_test[0],self.out_ts[0],
                                   places=self.splaces,msg=self.smsg)
            self.assertAlmostEqual(self.out_ts_test[1],self.out_ts[1],
                                   places=self.splaces,msg=self.smsg)
            self.assertAlmostEqual(self.out_ts_test[2],self.out_ts[2],
                                   places=self.splaces,msg=self.smsg)
            self.assertAlmostEqual(self.out_ts_test._delta_t,self.out_ts._delta_t,
                                   places=self.splaces,msg=self.smsg)

# Now, factories to create test cases for each available backend.
# The automation means that the default for each scheme will get created
# and run twice, once as 'Default' and once under its own name.

# Get our list of backends:
if _scheme == 'cpu':
    backends = pycbc.fft.cpu_backends
elif _scheme == 'cuda':
    backends == pycbc.fft.cuda_backends
elif _scheme == 'opencl':
    backends == pycbc.fft.opencl_backends

FFTTestClasses = []
for backend in backends:
    # This creates, for each backend, a new class derived from
    # both _BaseTestFFTClass and unittest.TestCase, and with
    # the additional property 'self.backend' set to the value
    # of backend.  One such class for each backend is appended
    # to the list
    klass = type('{0}_{1}_test'.format(_scheme,backend),
                 (_BaseTestFFTClass,unittest.TestCase),
                 {'backend': backend})
    FFTTestClasses.append(klass)


# Finally, we create suites and run them

if __name__ == '__main__':

    suite = unittest.TestSuite()
    for klass in FFTTestClasses:
        suite.addTest(unittest.TestLoader().loadTestsFromTestCase(klass))

    results = unittest.TextTestRunner(verbosity=2).run(suite)
