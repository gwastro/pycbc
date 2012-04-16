# Copyright (C) 2012  Josh Willis
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
import pycbc.array
import numpy
from numpy import dtype
import pycbc.fft
import unittest

class _BaseTestFFTClass(object):
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

    def test_forward_real_single(self):
        # First, test case when input array length is even
        self.in_even = pycbc.array.Array([1.0,-1.0,2.0,-2.0],
                                         dtype=dtype('float32'))
        self.in_pristine = pycbc.array.Array([1.0,-1.0,2.0,-2.0],
                                             dtype=dtype('float32'))
        self.out_even = pycbc.array.Array([0.0+0.0j,-1.0-1.0j,6.0+0.0j],
                                          dtype=dtype('complex64'))
        self.out_even_test = pycbc.array.zeros(3,dtype=dtype('complex64'))
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
        self.in_odd = pycbc.array.Array([1.0,2.0,2.0],
                                        dtype=dtype('float32'))
        self.out_odd = pycbc.array.Array([5.0+0.0j,-1.0+0.0j],
                                         dtype=dtype('complex64'))
        self.out_odd_test = pycbc.array.zeros(2,dtype=dtype('complex64'))
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
        self.out_badlen = pycbc.array.zeros(3,dtype=dtype('complex64'))
        args = [self.in_odd,self.out_badlen,self.backend]
        self.assertRaises(ValueError,pycbc.fft.fft,*args)
        self.out_badprec = pycbc.array.zeros(2,dtype=dtype('complex128'))
        args = [self.in_odd,self.out_badprec,self.backend]
        self.assertRaises(ValueError,pycbc.fft.fft,*args)
        self.out_baddtype = pycbc.array.zeros(2,dtype=dtype('float32'))
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
        self.in_even = pycbc.array.Array([1.0,-1.0,2.0,-2.0],
                                         dtype=dtype('float64'))
        self.in_pristine = pycbc.array.Array([1.0,-1.0,2.0,-2.0],
                                             dtype=dtype('float64'))
        self.out_even = pycbc.array.Array([0.0+0.0j,-1.0-1.0j,6.0+0.0j],
                                          dtype=dtype('complex128'))
        self.out_even_test = pycbc.array.zeros(3,dtype=dtype('complex128'))
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
        self.in_odd = pycbc.array.Array([1.0,2.0,2.0],
                                        dtype=dtype('float64'))
        self.out_odd = pycbc.array.Array([5.0+0.0j,-1.0+0.0j],
                                         dtype=dtype('complex128'))
        self.out_odd_test = pycbc.array.zeros(2,dtype=dtype('complex128'))
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
        self.out_badlen = pycbc.array.zeros(3,dtype=dtype('complex128'))
        args = [self.in_odd,self.out_badlen,self.backend]
        self.assertRaises(ValueError,pycbc.fft.fft,*args)
        self.out_badprec = pycbc.array.zeros(2,dtype=dtype('complex64'))
        args = [self.in_odd,self.out_badprec,self.backend]
        self.assertRaises(ValueError,pycbc.fft.fft,*args)
        self.out_baddtype = pycbc.array.zeros(2,dtype=dtype('float64'))
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
        self.in_even = pycbc.array.Array([0.0+0.0j,-1.0-1.0j,6.0+0.0j],
                                          dtype=dtype('complex64'))
        self.in_pristine = pycbc.array.Array([0.0+0.0j,-1.0-1.0j,6.0+0.0j],
                                          dtype=dtype('complex64'))
        self.out_even = pycbc.array.Array([4.0,-4.0,8.0,-8.0],
                                         dtype=dtype('float32'))
        self.out_even_test = pycbc.array.zeros(4,dtype=dtype('float32'))
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
        self.in_odd = pycbc.array.Array([5.0+0.0j,-1.0+0.0j],
                                         dtype=dtype('complex64'))
        self.out_odd = pycbc.array.Array([3.0,6.0,6.0],
                                        dtype=dtype('float32'))
        self.out_odd_test = pycbc.array.zeros(3,dtype=dtype('float32'))
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
        self.out_badlen = pycbc.array.zeros(2,dtype=dtype('float32'))
        args = [self.in_odd,self.out_badlen,self.backend]
        self.assertRaises(ValueError,pycbc.fft.ifft,*args)
        self.out_badprec = pycbc.array.zeros(3,dtype=dtype('float64'))
        args = [self.in_odd,self.out_badprec,self.backend]
        self.assertRaises(ValueError,pycbc.fft.ifft,*args)
        self.in_baddtype = pycbc.array.zeros(3,dtype=dtype('float32'))
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
        self.in_even = pycbc.array.Array([0.0+0.0j,-1.0-1.0j,6.0+0.0j],
                                          dtype=dtype('complex128'))
        self.in_pristine = pycbc.array.Array([0.0+0.0j,-1.0-1.0j,6.0+0.0j],
                                          dtype=dtype('complex128'))
        self.out_even = pycbc.array.Array([4.0,-4.0,8.0,-8.0],
                                         dtype=dtype('float64'))
        self.out_even_test = pycbc.array.zeros(4,dtype=dtype('float64'))
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
        self.in_odd = pycbc.array.Array([5.0+0.0j,-1.0+0.0j],
                                         dtype=dtype('complex128'))
        self.out_odd = pycbc.array.Array([3.0,6.0,6.0],
                                        dtype=dtype('float64'))
        self.out_odd_test = pycbc.array.zeros(3,dtype=dtype('float64'))
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
        self.out_badlen = pycbc.array.zeros(2,dtype=dtype('float64'))
        args = [self.in_odd,self.out_badlen,self.backend]
        self.assertRaises(ValueError,pycbc.fft.ifft,*args)
        self.out_badprec = pycbc.array.zeros(3,dtype=dtype('float32'))
        args = [self.in_odd,self.out_badprec,self.backend]
        self.assertRaises(ValueError,pycbc.fft.ifft,*args)
        self.in_baddtype = pycbc.array.zeros(3,dtype=dtype('float64'))
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
        self.in_cmplx = pycbc.array.Array([1.0+1.0j,2.0-2.0j],
                                          dtype=dtype('complex64'))
        self.in_pristine = pycbc.array.Array([1.0+1.0j,2.0-2.0j],
                                             dtype=dtype('complex64'))
        self.out_cmplx = pycbc.array.Array([3.0-1.0j,-1.0+3.0j],
                                           dtype=dtype('complex64'))
        self.out_cmplx_test = pycbc.array.zeros(2,dtype=dtype('complex64'))
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
        self.out_badlen = pycbc.array.zeros(3,dtype=dtype('complex64'))
        args = [self.in_cmplx,self.out_badlen,self.backend]
        self.assertRaises(ValueError,pycbc.fft.fft,*args)
        self.out_badprec = pycbc.array.zeros(2,dtype=dtype('complex128'))
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
        self.in_cmplx = pycbc.array.Array([1.0+1.0j,2.0-2.0j],
                                          dtype=dtype('complex128'))
        self.in_pristine = pycbc.array.Array([1.0+1.0j,2.0-2.0j],
                                             dtype=dtype('complex128'))
        self.out_cmplx = pycbc.array.Array([3.0-1.0j,-1.0+3.0j],
                                           dtype=dtype('complex128'))
        self.out_cmplx_test = pycbc.array.zeros(2,dtype=dtype('complex128'))
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
        self.out_badlen = pycbc.array.zeros(3,dtype=dtype('complex128'))
        args = [self.in_cmplx,self.out_badlen,self.backend]
        self.assertRaises(ValueError,pycbc.fft.fft,*args)
        self.out_badprec = pycbc.array.zeros(2,dtype=dtype('complex64'))
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
        self.in_cmplx = pycbc.array.Array([3.0-1.0j,-1.0+3.0j],
                                         dtype=dtype('complex64'))
        self.in_pristine = pycbc.array.Array([3.0-1.0j,-1.0+3.0j],
                                         dtype=dtype('complex64'))
        self.out_cmplx = pycbc.array.Array([2.0+2.0j,4.0-4.0j],
                                          dtype=dtype('complex64'))
        self.out_cmplx_test = pycbc.array.zeros(2,dtype=dtype('complex64'))
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
        self.out_badlen = pycbc.array.zeros(3,dtype=dtype('complex64'))
        args = [self.in_cmplx,self.out_badlen,self.backend]
        self.assertRaises(ValueError,pycbc.fft.ifft,*args)
        self.out_badprec = pycbc.array.zeros(2,dtype=dtype('complex128'))
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
        self.in_cmplx = pycbc.array.Array([3.0-1.0j,-1.0+3.0j],
                                         dtype=dtype('complex128'))
        self.in_pristine = pycbc.array.Array([3.0-1.0j,-1.0+3.0j],
                                         dtype=dtype('complex128'))
        self.out_cmplx = pycbc.array.Array([2.0+2.0j,4.0-4.0j],
                                          dtype=dtype('complex128'))
        self.out_cmplx_test = pycbc.array.zeros(2,dtype=dtype('complex128'))
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
        self.out_badlen = pycbc.array.zeros(3,dtype=dtype('complex128'))
        args = [self.in_cmplx,self.out_badlen,self.backend]
        self.assertRaises(ValueError,pycbc.fft.ifft,*args)
        self.out_badprec = pycbc.array.zeros(2,dtype=dtype('complex64'))
        args = [self.in_cmplx,self.out_badprec,self.backend]
        self.assertRaises(ValueError,pycbc.fft.ifft,*args)
        self.out_badarray = numpy.zeros(2,dtype=dtype('complex128'))
        args = [self.in_cmplx,self.out_badarray,self.backend]
        self.assertRaises(TypeError,pycbc.fft.ifft,*args)
        self.in_badarray = numpy.zeros(2,dtype=dtype('complex128'))
        args = [self.in_badarray,self.out_cmplx_test,self.backend]
        self.assertRaises(TypeError,pycbc.fft.ifft,*args)


# Now, factories to create test cases for each available backend.
# The automation means that the default for each scheme will get created
# and run twice, once as 'Default' and once under its own name.

CPUTestClasses = []

for backend in pycbc.fft.cpu_backends:
    klass = type('CPU_{0}Test'.format(backend),
                 (_BaseTestFFTClass,unittest.TestCase),
                 {'backend': backend})
    CPUTestClasses.append(klass)

if pycbc.HAVE_CUDA:
    CUDATestClasses = []
    for backend in pycbc.fft.cuda_backends:
        CUDATestClasses.append(type('CUDA_{0}Test'.format(backend),
                                    (_BaseTestFFTClass,unittest.TestCase),
                                    {'backend': backend}))

if pycbc.HAVE_OPENCL:
    OpenCLTestClasses = []
    for backend in pycbc.fft.opencl_backends:
        OpenCLTestClasses.append(type('OpenCL_{0}Test'.format(backend),
                                      (_BaseTestFFTClass,unittest.TestCase),
                                      {'backend': backend}))

# Finally, we create suites and run them, for every available backend

if __name__ == '__main__':

    suiteCPU = unittest.TestSuite()
    for klass in CPUTestClasses:
        suiteCPU.addTest(unittest.makeSuite(klass))

    unittest.TextTestRunner().run(suiteCPU)

    if pycbc.HAVE_CUDA:
        suiteCUDA = unittest.TestSuite()
        for klass in CUDATestClasses:
            suiteCUDA.addTest(unittest.makeSuite(klass))

        with CUDAScheme():
            unittest.TextTestRunner().run(suiteCUDA)

    if pycbc.HAVE_OPENCL:
        suiteOpenCL = unittest.TestSuite()
        for klass in OpenCLTestClasses:
            suiteOpenCL.addTest(unittest.makeSuite(klass))

        with OpenCLScheme():
            unittest.TextTestRunner().run(suiteOpenCL)
