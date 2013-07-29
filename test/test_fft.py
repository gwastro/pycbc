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
These are the unit-tests for the pycbc.fft subpackage

For both forward and reverse complex FFT, and forward R2C and reverse C2R FFT,
these tests validate that the FFT produces the correct output and does not
overwrite its input on small test cases where the answer may be hand-computed.
They also verify that an FFT of (larger sized) random data followed by the
reverse FFT gives back the original input (modulo the scaling factor of the
length of the array, when testing just the basic Array class; that factor is
not present for Time/FrequencySeries transformations).  Finally, it also checks
that the correct exceptions are raised when several different kinds of erroneous
inputs are given.

For the R2C (resp. C2R), the tests are run with input arrays (resp. output arrays)
that are both even and odd in length, since the length of those output arrays is
not the same as that of the input.  However these unit tests no longer check that
the imaginary parts of DC (and for even-length also Nyquist) outputs are exactly
zero, both because that complicates the testing framework considerably, and because
it will not in general hold for GPU algorithms (and those cannot be made to enforce
that, without unacceptable computational overhead).

All of these tests are performed for inputs from each of the three basic types
(Array, TimeSeries, and FrequencySeries) and for each precision (single and double).
They are also checked for each possible backend of the current scheme, which means
in particular that whichever backend is the default will be tested twice, once as
'Default' and once under its own name.
"""

import pycbc
import pycbc.scheme
import pycbc.types
from pycbc.types import Array as ar, TimeSeries as ts, FrequencySeries as fs
import numpy
from numpy import dtype, float32, float64, complex64, complex128, zeros
from numpy.random import randn
import pycbc.fft
import unittest
import sys
from utils import parse_args_all_schemes
from lal import LIGOTimeGPS as LTG

_scheme, _context = parse_args_all_schemes("FFT")

# Because we run many similar tests where we only vary dtypes, precisions,
# or Array/TimeSeries/FrequencySeries, it is helpful to define the following
# dictionarys and functions and then call those from within the actual tests.

# Map the same kind (real or complex) to the other precision
_other_prec = {float32:float64, float64:float32,
               complex64:complex128, complex128:complex64}
# Map the same precision to the other kind
_other_kind = {float32:complex64, float64:complex128,
                complex64:float32, complex128:float64}

# Our actual helper functions.  Note that these perform the necessary operations
# within the appropriate context, so they should not themselves be called inside
# of a context block.

def _test_fft(test_case,inarr,expec,tol):
    # Basic test to see that the forward FFT doesn't
    # overwrite its input and computes its output to
    # within the required accuracy.
    tc = test_case
    inty = type(inarr)
    in_pristine = inty(inarr)
    outty = type(expec)
    # Make a copy...
    outarr = outty(expec)
    # Clear it and change values (so we test something meaningful):
    outarr.clear()
    if hasattr(outarr,'_epoch'):
        outarr._epoch *= 5*tol
    if hasattr(outarr,'delta_t'):
        outarr._delta_t *= 5*tol
    if hasattr(outarr,'delta_f'):
        outarr._delta_f *= 5*tol
    with tc.context:
        pycbc.fft.fft(inarr,outarr,tc.backend)
        # First, verify that the input hasn't been overwritten
        emsg = 'FFT overwrote input array'
        tc.assertEqual(inarr,in_pristine,emsg)
        # Next, check that the output is correct to within tolerance.
        # That will require exact equality of all other meta-data
        emsg = 'FFT output differs by more than a factor of {0} from expected'.format(tol)
        if isinstance(outarr,ts) or isinstance(outarr,fs):
            tc.assertTrue(outarr.almost_equal_elem(expec,tol=tol,dtol=tol),msg=emsg)
        else:
            tc.assertTrue(outarr.almost_equal_elem(expec,tol=tol),msg=emsg)

def _test_ifft(test_case,inarr,expec,tol):
    # Basic test to see that the reverse FFT doesn't
    # overwrite its input and computes its output to
    # within the required accuracy.
    tc = test_case
    inty = type(inarr)
    in_pristine = inty(inarr)
    outty = type(expec)
    # Make a copy...
    outarr = outty(expec)
    # Clear it and change values (so we test something meaningful):
    outarr.clear()
    if hasattr(outarr,'_epoch'):
        outarr._epoch *= 5*tol
    if hasattr(outarr,'delta_t'):
        outarr._delta_t *= 5*tol
    if hasattr(outarr,'delta_f'):
        outarr._delta_f *= 5*tol
    with tc.context:
        pycbc.fft.ifft(inarr,outarr,tc.backend)
        # First, verify that the input hasn't been overwritten
        emsg = 'Inverse FFT overwrote input array'
        tc.assertEqual(inarr,in_pristine,emsg)
        # Next, check that the output is correct to within tolerance.
        # That will require exact equality of all other meta-data
        emsg = 'Inverse FFT output differs by more than a factor of {0} from expected'.format(tol)
        if isinstance(outarr,ts) or isinstance(outarr,fs):
            tc.assertTrue(outarr.almost_equal_elem(expec,tol=tol,dtol=tol),msg=emsg)
        else:
            tc.assertTrue(outarr.almost_equal_elem(expec,tol=tol),msg=emsg)

def _test_random(test_case,inarr,outarr,tol):
    tc = test_case
    # Test that applying a transform and its inverse to reasonably long, random
    # input gives back the (appropriately scaled) input. We must allow for numerical
    # error, and we use 100*tol for the relative error on an element-by-element
    # comparison, but just tol for the norm comparison.
    #
    # First test IFFT(FFT(random))
    # The numpy randn(n) provides an array of n numbers drawn from standard normal
    if dtype(inarr).kind is 'c':
        inarr._data[:] = randn(len(inarr)) +1j*randn(len(inarr))
    else:
        inarr._data[:] = randn(len(inarr))
    incopy = type(inarr)(inarr)
    outarr.clear()
    # An FFT followed by IFFT gives Array scaled by len(Array), but for
    # Time/FrequencySeries there should be no scaling.
    if type(inarr) == pycbc.types.Array:
        incopy *= len(inarr)
    with tc.context:
        pycbc.fft.fft(inarr,outarr,tc.backend)
        pycbc.fft.ifft(outarr,inarr,tc.backend)
        emsg="IFFT(FFT(random)) did not reproduce original array to within tolerance {0}".format(tol)
        if isinstance(incopy,ts) or isinstance(incopy,fs):
            tc.assertTrue(incopy.almost_equal_elem(inarr,tol=100*tol,dtol=tol),
                          msg=emsg)
            tc.assertTrue(incopy.almost_equal_norm(inarr,tol=tol,dtol=tol),
                          msg=emsg)
        else:
            tc.assertTrue(incopy.almost_equal_elem(inarr,tol=100*tol),
                          msg=emsg)
            tc.assertTrue(incopy.almost_equal_norm(inarr,tol=tol),
                          msg=emsg)
    # Now the same for FFT(IFFT(random))
    if dtype(outarr).kind is 'c':
        outarr._data[:] = randn(len(outarr))+1j*randn(len(outarr))
    else:
        outarr._data[:] = randn(len(outarr))
    inarr.clear()
    outcopy = type(outarr)(outarr)
    if type(outarr) == pycbc.types.Array:
        outcopy *= len(outarr)
    with tc.context:
        pycbc.fft.ifft(outarr,inarr,tc.backend)
        pycbc.fft.fft(inarr,outarr,tc.backend)
        emsg="FFT(IFFT(random)) did not reproduce original array to within tolerance {0}".format(tol)
        if isinstance(outcopy,ts) or isinstance(outcopy,fs):
            tc.assertTrue(outcopy.almost_equal_elem(outarr,tol=100*tol,dtol=tol),
                          msg=emsg)
            tc.assertTrue(outcopy.almost_equal_norm(outarr,tol=tol,dtol=tol),
                          msg=emsg)
        else:
            tc.assertTrue(outcopy.almost_equal_elem(outarr,tol=100*tol),
                          msg=emsg)
            tc.assertTrue(outcopy.almost_equal_norm(outarr,tol=tol),
                          msg=emsg)

def _test_raise_excep_fft(test_case,inarr,outarr,other_args={}):
    # As far as can be told from the unittest module documentation, the
    # 'assertRaises' tests do not permit a custom message.  So more
    # comments than usual here, to help diagnose and test failures.
    #
    # The 'other_args' argument is needed to pass additional keywords to
    # the constructors of some types (T/F series); we cannot simply copy since
    # the whole point is to vary the input/output in some way that should cause
    # an exception.
    tc = test_case
    with tc.context:
        outty = type(outarr)
        outzer = pycbc.types.zeros(len(outarr))
        # If we give an output array that is wrong only in length, raise ValueError:
        out_badlen = outty(pycbc.types.zeros(len(outarr)+1),dtype=outarr.dtype,**other_args)
        args = [inarr,out_badlen,tc.backend]
        tc.assertRaises(ValueError,pycbc.fft.fft,*args)
        # If we give an output array that has the wrong precision, raise ValueError:
        out_badprec = outty(outzer,dtype=_other_prec[dtype(outarr).type],**other_args)
        args = [inarr,out_badprec,tc.backend]
        tc.assertRaises(ValueError,pycbc.fft.fft,*args)
        # If we give an output array that has the wrong kind (real or complex) but
        # correct precision, then raise a ValueError:
        out_badkind = outty(outzer,dtype=_other_kind[dtype(outarr).type],**other_args)
        args = [inarr,out_badkind,tc.backend]
        tc.assertRaises(ValueError,pycbc.fft.fft,*args)
        # If we give an output array that isn't a PyCBC type, raise TypeError:
        out_badtype = numpy.zeros(len(outarr),dtype=outarr.dtype)
        args = [inarr,out_badtype,tc.backend]
        tc.assertRaises(TypeError,pycbc.fft.fft,*args)
        # If we give an input array that isn't a PyCBC type, raise TypeError:
        in_badtype = numpy.zeros(len(inarr),dtype=inarr.dtype)
        args = [in_badtype,outarr,tc.backend]
        tc.assertRaises(TypeError,pycbc.fft.fft,*args)

def _test_raise_excep_ifft(test_case,inarr,outarr,other_args={}):
    # As far as can be told from the unittest module documentation, the
    # 'assertRaises' tests do not permit a custom message.  So more
    # comments than usual here, to help diagnose and test failures.
    #
    # The 'other_args' argument is needed to pass additional keywords to
    # the constructors of some types (T/F series); we cannot simply copy since
    # the whole point is to vary the input/output in some way that should cause
    # an exception.
    tc = test_case
    with tc.context:
        outty = type(outarr)
        outzer = pycbc.types.zeros(len(outarr))
        # If we give an output array that is wrong only in length, raise ValueError:
        out_badlen = outty(pycbc.types.zeros(len(outarr)+1),dtype=outarr.dtype,**other_args)
        args = [inarr,out_badlen,tc.backend]
        tc.assertRaises(ValueError,pycbc.fft.ifft,*args)
        # If we give an output array that has the wrong precision, raise ValueError:
        out_badprec = outty(outzer,dtype=_other_prec[dtype(outarr).type],**other_args)
        args = [inarr,out_badprec,tc.backend]
        tc.assertRaises(ValueError,pycbc.fft.ifft,*args)
        # If we give an output array that has the wrong kind (real or complex) but
        # correct precision, then raise a ValueError:
        out_badkind = outty(outzer,dtype=_other_kind[dtype(outarr).type],**other_args)
        args = [inarr,out_badkind,tc.backend]
        tc.assertRaises(ValueError,pycbc.fft.ifft,*args)
        # If we give an output array that isn't a PyCBC type, raise TypeError:
        out_badtype = numpy.zeros(len(outarr),dtype=outarr.dtype)
        args = [inarr,out_badtype,tc.backend]
        tc.assertRaises(TypeError,pycbc.fft.ifft,*args)
        # If we give an input array that isn't a PyCBC type, raise TypeError:
        in_badtype = numpy.zeros(len(inarr),dtype=inarr.dtype)
        args = [in_badtype,outarr,tc.backend]
        tc.assertRaises(TypeError,pycbc.fft.ifft,*args)

class _BaseTestFFTClass(unittest.TestCase):
    """
    This is the base class from which unit tests for all FFT backends
    are derived.
    """
    def setUp(self):
        # Dictionary to convert a dtype to a relative precision to test
        self.tdict = { float32: 1e-6, float64: 1e-14,
                       complex64: 1e-6, complex128: 1e-14}
        # Save our scheme and context
        self.scheme = _scheme
        self.context = _context
        # Next we set up various lists that are used to build our 'known'
        # test, which are repeated for a variety of different precisions
        # and basic types. All of the lists should be consistent with a
        # direct calculation in the formulas of the "What FFTW actually
        # computes" section of the FFTW manual.

        # First, R2C transforms.  We have both even and odd length inputs,
        # since the appropriate values for the output lengths vary.
        self.in_r2c_e = [1.0,-1.0,2.0,-2.0]
        self.out_r2c_e = [0.0+0.0j,-1.0-1.0j,6.0+0.0j]
        self.in_r2c_o = [1.0,2.0,2.0]
        self.out_r2c_o = [5.0+0.0j,-1.0+0.0j]
        # Next, C2R transforms, again for both even and odd lengths
        self.in_c2r_e = [0.0+0.0j,-1.0-1.0j,6.0+0.0j]
        self.out_c2r_e = [4.0, -4.0, 8.0, -8.0]
        self.in_c2r_o = [5.0+0.0j,-1.0+0.0j]
        self.out_c2r_o = [3.0,6.0,6.0]
        # Finally, C2C transforms, where we don't do both even and odd,
        # but do have different lists for fwd and rev (i.e., fft and ifft)
        self.in_c2c_fwd = [1.0+1.0j,2.0-2.0j]
        self.out_c2c_fwd = [3.0-1.0j,-1.0+3.0j]
        self.in_c2c_rev = [3.0-1.0j,-1.0+3.0j]
        self.out_c2c_rev = [2.0+2.0j,4.0-4.0j]
        # For Time/FrequencySeries, we want to test with a non-trivial epoch
        self.epoch = LTG(3,4)
        # When we need a delta_t or delta_f for input, use this.
        # Output-appropriate variable is computed.
        self.delta = 1.0/4096.0
        # Length of our random arrays, for both real and complex arrays.
        self.rand_len_r = 2046
        self.rand_len_c = 1024

    def test_fwd_real_arr(self):
        for fwd_dtype in [float32,float64]:
            # Even input
            inarr = ar(self.in_r2c_e,dtype=fwd_dtype)
            outexp = ar(self.out_r2c_e,dtype=_other_kind[fwd_dtype])
            _test_fft(self,inarr,outexp,self.tdict[fwd_dtype])
            # Odd input
            inarr = ar(self.in_r2c_o,dtype=fwd_dtype)
            outexp = ar(self.out_r2c_o,dtype=_other_kind[fwd_dtype])
            _test_fft(self,inarr,outexp,self.tdict[fwd_dtype])
            # Random
            rand_inarr = ar(zeros(self.rand_len_r,dtype=fwd_dtype))
            rand_outarr = ar(zeros(self.rand_len_c,dtype=_other_kind[fwd_dtype]))
            _test_random(self,rand_inarr,rand_outarr,self.tdict[fwd_dtype])
            # Clean these up since they could be big:
            del rand_inarr
            del rand_outarr
            # Check that exceptions are raised.  Need input and
            # output arrays; just reuse inarr and outexp (values won't
            # matter, we're just checking exceptions).
            _test_raise_excep_fft(self,inarr,outexp)

    def test_fwd_real_ts(self):
        for fwd_dtype in [float32,float64]:
            delta_t = self.delta
            # Even input
            inarr = ts(self.in_r2c_e,dtype=fwd_dtype,delta_t=delta_t,epoch=self.epoch)
            delta_f = 1.0/(inarr.delta_t * len(inarr))
            outexp = fs(self.out_r2c_e,dtype=_other_kind[fwd_dtype],delta_f=delta_f,epoch=self.epoch)
            _test_fft(self,inarr,outexp,self.tdict[fwd_dtype])
            # Odd input
            inarr = ts(self.in_r2c_o,dtype=fwd_dtype,delta_t=delta_t,epoch=self.epoch)
            delta_f = 1.0/(inarr.delta_t * len(inarr))
            outexp = fs(self.out_r2c_o,dtype=_other_kind[fwd_dtype],delta_f=delta_f,epoch=self.epoch)
            _test_fft(self,inarr,outexp,self.tdict[fwd_dtype])
            # Random
            rand_inarr = ts(zeros(self.rand_len_r,dtype=fwd_dtype),epoch=self.epoch,delta_t=delta_t)
            delta_f = 1.0/(rand_inarr.delta_t * len(rand_inarr))
            rand_outarr = fs(zeros(self.rand_len_c,dtype=_other_kind[fwd_dtype]),epoch=self.epoch,
                             delta_f=delta_f)
            _test_random(self,rand_inarr,rand_outarr,self.tdict[fwd_dtype])
            # Clean these up since they could be big:
            del rand_inarr
            del rand_outarr
            # Check that exceptions are raised.  Need input and
            # output arrays; just reuse inarr and outexp (values won't
            # matter, we're just checking exceptions).
            output_args = {"delta_f": self.delta, "epoch": self.epoch}
            _test_raise_excep_fft(self,inarr,outexp,output_args)

    def test_fwd_real_fs(self):
        for fwd_dtype in [float32,float64]:
            delta_f = self.delta
            # Even input
            inarr = fs(self.in_r2c_e,dtype=fwd_dtype,delta_f=delta_f,epoch=self.epoch)
            delta_t = 1.0/(inarr.delta_f * len(inarr))
            outexp = ts(self.out_r2c_e,dtype=_other_kind[fwd_dtype],delta_t=delta_t,epoch=self.epoch)
            _test_fft(self,inarr,outexp,self.tdict[fwd_dtype])
            # Odd input
            inarr = fs(self.in_r2c_o,dtype=fwd_dtype,delta_f=delta_f,epoch=self.epoch)
            delta_t = 1.0/(inarr.delta_f * len(inarr))
            outexp = ts(self.out_r2c_o,dtype=_other_kind[fwd_dtype],delta_t=delta_t,epoch=self.epoch)
            _test_fft(self,inarr,outexp,self.tdict[fwd_dtype])
            # Random
            rand_inarr = fs(zeros(self.rand_len_r,dtype=fwd_dtype),epoch=self.epoch,delta_f=delta_f)
            delta_t = 1.0/(rand_inarr.delta_f * len(rand_inarr))
            rand_outarr = fs(zeros(self.rand_len_c,dtype=_other_kind[fwd_dtype]),epoch=self.epoch,
                             delta_t=delta_t)
            _test_random(self,rand_inarr,rand_outarr,self.tdict[fwd_dtype])
            # Clean these up since they could be big:
            del rand_inarr
            del rand_outarr
            # Check that exceptions are raised.  Need input and
            # output arrays; just reuse inarr and outexp (values won't
            # matter, we're just checking exceptions).
            output_args = {"delta_t": self.delta, "epoch": self.epoch}
            _test_raise_excep_fft(self,inarr,outexp,output_args)

    def test_rev_real_arr(self):
        for rev_dtype in [float32,float64]:
            # Even input
            inarr = ar(self.in_c2r_e,dtype=_other_kind[rev_dtype])
            outexp = ar(self.out_c2r_e,dtype=rev_dtype)
            _test_ifft(self,inarr,outexp,self.tdict[rev_dtype])
            # Odd input
            inarr = ar(self.in_c2r_o,dtype=_other_kind[rev_dtype])
            outexp = ar(self.out_c2r_o,dtype=rev_dtype)
            _test_ifft(self,inarr,outexp,self.tdict[rev_dtype])
            # Random---we don't do that in 'reverse' tests, since both
            # directions are already tested in forward, and if we just passed
            # in arrays in the other order we'd only get exceptions
            #
            # Check that exceptions are raised.  Need input and
            # output arrays; just reuse inarr and outexp (values won't
            # matter, we're just checking exceptions).
            _test_raise_excep_ifft(self,inarr,outexp)

    def test_rev_real_ts(self):
        for rev_dtype in [float32,float64]:
            delta_t = self.delta
            # Even input
            inarr = ts(self.in_c2r_e,dtype=_other_kind[rev_dtype],delta_t=delta_t,epoch=self.epoch)
            delta_f = 1.0/(delta_t*len(self.out_c2r_e))
            outexp = fs(self.out_c2r_e,dtype=rev_dtype,delta_f=delta_f,epoch=self.epoch)
            _test_ifft(self,inarr,outexp,self.tdict[rev_dtype])
            # Odd input
            inarr = ts(self.in_c2r_o,dtype=_other_kind[rev_dtype])
            delta_f = 1.0/(delta_t*len(self.out_c2r_o))
            outexp = fs(self.out_c2r_o,dtype=rev_dtype,delta_f=delta_f,epoch=self.epoch)
            _test_ifft(self,inarr,outexp,self.tdict[rev_dtype])
            # Random---we don't do that in 'reverse' tests, since both
            # directions are already tested in forward, and if we just passed
            # in arrays in the other order we'd only get exceptions
            #
            # Check that exceptions are raised.  Need input and
            # output arrays; just reuse inarr and outexp (values won't
            # matter, we're just checking exceptions).
            output_args = {"delta_f": self.delta, "epoch": self.epoch}
            _test_raise_excep_ifft(self,inarr,outexp)

    def test_rev_real_fs(self):
        for rev_dtype in [float32,float64]:
            delta_f = self.delta
            # Even input
            inarr = fs(self.in_c2r_e,dtype=_other_kind[rev_dtype],delta_f=delta_f,epoch=self.epoch)
            delta_t = 1.0/(delta_f*len(self.out_c2r_e))
            outexp = ts(self.out_c2r_e,dtype=rev_dtype,delta_t=delta_t,epoch=self.epoch)
            _test_ifft(self,inarr,outexp,self.tdict[rev_dtype])
            # Odd input
            inarr = fs(self.in_c2r_o,dtype=_other_kind[rev_dtype])
            delta_t = 1.0/(delta_f*len(self.out_c2r_o))
            outexp = ts(self.out_c2r_o,dtype=rev_dtype,delta_t=delta_t,epoch=self.epoch)
            _test_ifft(self,inarr,outexp,self.tdict[rev_dtype])
            # Random---we don't do that in 'reverse' tests, since both
            # directions are already tested in forward, and if we just passed
            # in arrays in the other order we'd only get exceptions
            #
            # Check that exceptions are raised.  Need input and
            # output arrays; just reuse inarr and outexp (values won't
            # matter, we're just checking exceptions).
            output_args = {"delta_t": self.delta, "epoch": self.epoch}
            _test_raise_excep_ifft(self,inarr,outexp)

    def test_fwd_complex_arr(self):
        for fwd_dtype in [complex64,complex128]:
            # Don't do separate even/odd tests for complex
            inarr = ar(self.in_c2c_fwd,dtype=fwd_dtype)
            outexp = ar(self.out_c2c_fwd,dtype=fwd_dtype)
            _test_fft(self,inarr,outexp,self.tdict[fwd_dtype])
            # Random
            rand_inarr = ar(zeros(self.rand_len_c,dtype=fwd_dtype))
            rand_outarr = ar(zeros(self.rand_len_c,dtype=fwd_dtype))
            _test_random(self,rand_inarr,rand_outarr,self.tdict[fwd_dtype])
            # Clean these up since they could be big:
            del rand_inarr
            del rand_outarr
            # Check that exceptions are raised.  Need input and
            # output arrays; just reuse inarr and outexp (values won't
            # matter, we're just checking exceptions).
            _test_raise_excep_fft(self,inarr,outexp)

    def test_fwd_complex_ts(self):
        for fwd_dtype in [complex64,complex128]:
            delta_t = self.delta
            # Don't do separate even/odd tests for complex
            inarr = ts(self.in_c2c_fwd,dtype=fwd_dtype,delta_t=delta_t,epoch=self.epoch)
            delta_f = 1.0/(delta_t * len(inarr))
            outexp = fs(self.out_c2c_fwd,dtype=fwd_dtype,delta_f=delta_f,epoch=self.epoch)
            _test_fft(self,inarr,outexp,self.tdict[fwd_dtype])
            # Random
            rand_inarr = ts(zeros(self.rand_len_c,dtype=fwd_dtype),delta_t=delta_t,epoch=self.epoch)
            delta_f = 1.0/(delta_t*len(rand_inarr))
            rand_outarr = fs(zeros(self.rand_len_c,dtype=fwd_dtype),delta_f=delta_f,epoch=self.epoch)
            _test_random(self,rand_inarr,rand_outarr,self.tdict[fwd_dtype])
            # Clean these up since they could be big:
            del rand_inarr
            del rand_outarr
            # Check that exceptions are raised.  Need input and
            # output arrays; just reuse inarr and outexp (values won't
            # matter, we're just checking exceptions).
            output_args = {"delta_f": self.delta, "epoch": self.epoch}
            _test_raise_excep_fft(self,inarr,outexp,output_args)

    def test_fwd_complex_fs(self):
        for fwd_dtype in [complex64,complex128]:
            delta_f = self.delta
            # Don't do separate even/odd tests for complex
            inarr = fs(self.in_c2c_fwd,dtype=fwd_dtype,delta_f=delta_f,epoch=self.epoch)
            delta_t = 1.0/(delta_f * len(inarr))
            outexp = ts(self.out_c2c_fwd,dtype=fwd_dtype,delta_t=delta_t,epoch=self.epoch)
            _test_fft(self,inarr,outexp,self.tdict[fwd_dtype])
            # Random
            rand_inarr = fs(zeros(self.rand_len_c,dtype=fwd_dtype),delta_f=delta_f,epoch=self.epoch)
            delta_t = 1.0/(delta_t*len(rand_inarr))
            rand_outarr = ts(zeros(self.rand_len_c,dtype=fwd_dtype),delta_t=delta_t,epoch=self.epoch)
            _test_random(self,rand_inarr,rand_outarr,self.tdict[fwd_dtype])
            # Clean these up since they could be big:
            del rand_inarr
            del rand_outarr
            # Check that exceptions are raised.  Need input and
            # output arrays; just reuse inarr and outexp (values won't
            # matter, we're just checking exceptions).
            output_args = {"delta_t": self.delta, "epoch": self.epoch}
            _test_raise_excep_fft(self,inarr,outexp,output_args)

    def test_rev_complex_arr(self):
        for rev_dtype in [complex64,complex128]:
            # Don't do separate even/odd tests for complex
            inarr = ar(self.in_c2c_rev,dtype=rev_dtype)
            outexp = ar(self.out_c2c_rev,dtype=rev_dtype)
            _test_ifft(self,inarr,outexp,self.tdict[rev_dtype])
            # Random---we don't do that in 'reverse' tests, since both
            # directions are already tested in forward, and if we just passed
            # in arrays in the other order we'd only get exceptions
            #
            # Check that exceptions are raised.  Need input and
            # output arrays; just reuse inarr and outexp (values won't
            # matter, we're just checking exceptions).
            _test_raise_excep_ifft(self,inarr,outexp)

    def test_rev_complex_ts(self):
        for rev_dtype in [complex64,complex128]:
            delta_t = self.delta
            # Don't do separate even/odd tests for complex
            inarr = ts(self.in_c2c_rev,dtype=rev_dtype,delta_t=delta_t,epoch=self.epoch)
            delta_f = 1.0/(delta_t*len(self.out_c2c_rev))
            outexp = fs(self.out_c2c_rev,dtype=rev_dtype,delta_f=delta_f,epoch=self.epoch)
            _test_ifft(self,inarr,outexp,self.tdict[rev_dtype])
            # Random---we don't do that in 'reverse' tests, since both
            # directions are already tested in forward, and if we just passed
            # in arrays in the other order we'd only get exceptions
            #
            # Check that exceptions are raised.  Need input and
            # output arrays; just reuse inarr and outexp (values won't
            # matter, we're just checking exceptions).
            output_args = {"delta_f": self.delta, "epoch": self.epoch}
            _test_raise_excep_ifft(self,inarr,outexp,output_args)

    def test_rev_complex_fs(self):
        for rev_dtype in [complex64,complex128]:
            delta_f = self.delta
            # Don't do separate even/odd tests for complex
            inarr = fs(self.in_c2c_rev,dtype=rev_dtype,delta_f=delta_f,epoch=self.epoch)
            delta_t = 1.0/(delta_f*len(self.out_c2c_rev))
            outexp = ts(self.out_c2c_rev,dtype=rev_dtype,delta_t=delta_t,epoch=self.epoch)
            _test_ifft(self,inarr,outexp,self.tdict[rev_dtype])
            # Random---we don't do that in 'reverse' tests, since both
            # directions are already tested in forward, and if we just passed
            # in arrays in the other order we'd only get exceptions
            #
            # Check that exceptions are raised.  Need input and
            # output arrays; just reuse inarr and outexp (values won't
            # matter, we're just checking exceptions).
            output_args = {"delta_t": self.delta, "epoch": self.epoch}
            _test_raise_excep_ifft(self,inarr,outexp,output_args)

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
                 (_BaseTestFFTClass,),
                 {'backend': backend})
    FFTTestClasses.append(klass)


# Finally, we create suites and run them

if __name__ == '__main__':

    suite = unittest.TestSuite()
    for klass in FFTTestClasses:
        suite.addTest(unittest.TestLoader().loadTestsFromTestCase(klass))

    results = unittest.TextTestRunner(verbosity=2).run(suite)
