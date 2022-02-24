# Copyright (C) 2012--2013  Alex Nitz, Josh Willis, Andrew Miller
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
This module contains a few helper functions designed to make writing PyCBC
unit tests easier, while still allowing the tests to be run on CPU and CUDA

All tests starting with 'test_' in the test subdirectory of pycbc are run
whenever the command 'python setup.py test' is given.  That command will
attempt to call each test, passing it the argument '-s <scheme>' where
scheme is each of 'cpu', 'cuda' in turn. Unit tests
designed to validate code that should run under multiple schemes should accept
each of these options, rerunning the same tests under each successive scheme.

This will usually be done by putting something like:

   _scheme, _context = parse_args_all_schemes('MyFeature')

before the definition of the unit test class.  In the definition of 'setUp'
(which is a mandatory function that must be defined in the class, per the
Python unittest module design) one would then usually have:

   self.scheme = _scheme
   self.context = _context

and those properties of any instance can then be used by later tests defined
in the class; for example, beginning a block with 'with self.context:' to
ensure the appropriate context manager is used for that scheme.

Some other unit tests may be for features or sub-packages that are not GPU
capable, and cannot be meaningfully tested on the GPU.  Those tests, after
importing pycbc.test.utils, should instead just call:

   parse_args_cpu_only('MyFeature')

This call is needed because the tests must still be able to accept the arguments
specifying a GPU environment (since setup.py does not know which tests are GPU
capable and which are not) but when called with a GPU scheme will exit immediately.

Both functions take a single string as an argument.  That string is used to customize
the heading of all of the tests (according to feature and scheme) to make the output
of running all of the unit tests somewhat easier to parse when they are all run at
once.
"""

import pycbc
import optparse
from sys import exit as _exit
from optparse import OptionParser
from pycbc.scheme import CPUScheme, CUDAScheme
from numpy import float32, float64, complex64, complex128
from pycbc.types import Array


def _check_scheme_all(option, opt_str, scheme, parser):
    if scheme=='cuda' and not pycbc.HAVE_CUDA:
        raise optparse.OptionValueError("CUDA not found")

    setattr (parser.values, option.dest, scheme)


def parse_args_all_schemes(feature_str):
    _parser = OptionParser()
    _parser.add_option('--scheme','-s', action='callback', type = 'choice',
                       choices = ('cpu','cuda'),
                       default = 'cpu', dest = 'scheme', callback = _check_scheme_all,
                       help = 'specifies processing scheme, can be cpu [default], cuda')
    _parser.add_option('--device-num','-d', action='store', type = 'int',
                       dest = 'devicenum', default=0,
                       help = 'specifies a GPU device to use for CUDA, 0 by default')
    (_opt_list, _args) = _parser.parse_args()

    # Changing the optvalues to a dict makes them easier to read
    _options = vars(_opt_list)

    _scheme = _options['scheme']

    if _scheme == 'cpu':
        _context = CPUScheme()
    if _scheme == 'cuda':
        _context = CUDAScheme(device_num=_options['devicenum'])

    _scheme_dict = { 'cpu': 'CPU', 'cuda': 'CUDA'}

    print(72*'=')
    print("Running {0} unit tests for {1}:".format(_scheme_dict[_scheme],feature_str))

    return [_scheme,_context]

def _check_scheme_cpu(option, opt_str, scheme, parser):
    if scheme=='cuda':
        exit(0)

    setattr (parser.values, option.dest, scheme)


def parse_args_cpu_only(feature_str):
    _parser = OptionParser()
    _parser.add_option('--scheme','-s', action='callback', type = 'choice',
                       choices = ('cpu','cuda'),
                       default = 'cpu', dest = 'scheme', callback = _check_scheme_cpu,
                       help = 'specifies processing scheme, can be cpu [default], cuda')
    _parser.add_option('--device-num','-d', action='store', type = 'int',
                       dest = 'devicenum', default=0,
                       help = 'specifies a GPU device to use for CUDA, 0 by default')
    (_opt_list, _args) = _parser.parse_args()

    # In this case, the only reason we parsed the arguments was to exit if we were given
    # a GPU scheme.  So if we get here we're on the CPU, and should print out our message
    # and return.

    print(72*'=')
    print("Running {0} unit tests for {1}:".format('CPU', feature_str))

    return

def simple_exit(results):
    """
    A simpler version of exit_based_on_results(); this function causes the script
    to exit normally with return value of zero if and only if all tests within the
    script passed and had no errors. Otherwise it returns the number of failures
    plus the number of errors

    Parameters
    ----------
    results: an instance of unittest.TestResult, returned (for instance) from a call such as
        results = unittest.TextTestRunner(verbosity=2).run(suite)
    """
    if results.wasSuccessful():
        _exit(0)
    else:
        nfail = len(results.errors)+len(results.failures)
        _exit(nfail)

def exit_based_on_results(results):
    """
    A probably-obsolete function to exit from a unit test-script with a status that depends
    on whether or not the only errors or failures were NotImplemented errors.  Specifically,
    if the unit-test suite execution encoded in results was:
       All tests successful:                               Exit 0
       All tests successful or only NotImplemented errors: Exit 1
       Some tests either failed or had other errors:       Exit 2
    The intent was that failures due to missing features be treated differently (especially
    when that happens on one of the GPU schemes) and that these exit statuses could then be
    interpreted by NMI or some other automatic build/test system accordingly.

    Parameters
    ----------
    results: an instance of unittest.TestResult, returned (for instance) from a call such as
        results = unittest.TextTestRunner(verbosity=2).run(suite)
    """
    NotImpErrors = 0
    for error in results.errors:
        for errormsg in error:
            if type(errormsg) is str:
                if 'NotImplemented' in errormsg:
                    NotImpErrors +=1
                    break
    if results.wasSuccessful():
        _exit(0)
    elif len(results.failures)==0 and len(results.errors)==NotImpErrors:
        _exit(1)
    else:
        _exit(2)

# Copied over from base_array.py so we can refactor it to remove
# the scheme shuffling and take advantage of the new equality/almost
# equal methods of the types.

# The following dictionary converts dtypes into the corresponding real
# dtype; it is needed for functions that return arrays of the same
# precision but which are always real.

_real_dtype_dict = { float32: float32, complex64: float32,
                     float64: float64, complex128: float64 }

class array_base(object):
    def setNumbers(self):
        # We create instances of our types, and need to know the generic answer
        # type so that we can convert the many basic lists into the appropriate
        # precision and kind.  This logic essentially implements (for our limited
        # use cases) what is in the function numpy.result_type, but that function
        # doesn't become available until Numpy 1.6.0
        if self.kind == 'real':
            if self.okind == 'real':
                self.result_dtype = self.dtype
            else:
                self.result_dtype = self.odtype
        else:
            self.result_dtype = self.dtype

        self.rdtype = _real_dtype_dict[self.dtype]

        # The child class (testing one of Array, TimeSeries, or FrequencySeries)
        # should set the following in its setUp method before this method (setNumbers)
        # is called:
        #    self.type = one of [Array,TimeSeries,FrequencySeries]
        #    self.kwds = dict for other kwd args beyond 'dtype'; normally an
        #                epoch and one of delta_t or delta_f

        # These are the values that should be used to initialize the test arrays.
        if self.kind == 'real':
            self.a = self.type([5,3,1],dtype=self.dtype,**self.kwds)
            self.alist = [5,3,1]
        else:
            self.a = self.type([5+1j,3+3j,1+5j],dtype=self.dtype,**self.kwds)
            self.alist = [5+1j,3+3j,1+5j]
        if self.okind == 'real':
            self.b = self.type([10,8,6],dtype=self.odtype,**self.kwds)
            self.blist = [10,8,6]
        else:
            self.b = self.type([10+6j,8+4j,6+2j],dtype=self.odtype,**self.kwds)
            self.blist = [10+6j,8+4j,6+2j]
        # And the scalar to test on
        if self.okind == 'real':
            self.scalar = 5
        else:
            self.scalar = 5+2j

        # The weights used in the weighted inner product test are always an Array,
        # regardless of the types whose inner product is being tested.
        self.w = Array([1, 2, 1],dtype=self.dtype)

        # All the answers are stored here to make it easier to read in the actual tests.
        # Again, it makes a difference whether they are complex or real valued, so there
        # are four sets of possible answers, depending on the dtypes.
        if self.kind == 'real' and self.okind == 'real':
            self.cumsum=self.type([5,8,9],dtype=self.dtype,**self.kwds)

            self.mul = self.type([50, 24, 6],dtype=self.result_dtype,**self.kwds)
            self.mul_s = self.type([25, 15, 5],dtype=self.result_dtype,**self.kwds)

            self.add = self.type([15, 11, 7],dtype=self.result_dtype,**self.kwds)
            self.add_s = self.type([10, 8, 6],dtype=self.result_dtype,**self.kwds)

            #self.div = [.5, 3./8., 1./6.]
            self.div = self.type([.5, 0.375, .16666666666666666667],dtype=self.result_dtype,**self.kwds)
            #self.div_s = [1., 3./5., 1./5.]
            self.div_s = self.type([1., 0.6, 0.2],dtype=self.result_dtype,**self.kwds)

            #self.rdiv = [2., 8./3., 6.]
            self.rdiv = self.type([2., 2.66666666666666666667, 6.],dtype=self.result_dtype,**self.kwds)
            #self.rdiv_s = [1., 5./3., 5.]
            self.rdiv_s = self.type([1., 1.66666666666666666667, 5.],dtype=self.result_dtype,**self.kwds)

            self.sub = self.type([-5, -5, -5],dtype=self.result_dtype,**self.kwds)
            self.sub_s = self.type([0, -2, -4],dtype=self.result_dtype,**self.kwds)

            self.rsub = self.type([5, 5, 5],dtype=self.result_dtype,**self.kwds)
            self.rsub_s = self.type([0, 2, 4],dtype=self.result_dtype,**self.kwds)

            self.pow1 = self.type([25., 9., 1.],dtype=self.dtype,**self.kwds)
            #self.pow2 = [pow(5,-1.5), pow(3,-1.5), pow(1,-1.5)]
            self.pow2 = self.type([0.08944271909999158786,
                                   0.19245008972987525484, 1.],dtype=self.dtype,**self.kwds)

            self.abs = self.type([5, 3, 1],dtype=self.rdtype,**self.kwds)
            self.real = self.type([5,3,1],dtype=self.rdtype,**self.kwds)
            self.imag = self.type([0, 0, 0],dtype=self.rdtype,**self.kwds)
            self.conj = self.type([5, 3, 1],dtype=self.dtype,**self.kwds)

            self.sum = 9

            self.dot = 80
            self.inner = self.dot
            self.weighted_inner = 68

        if self.kind =='real' and self.okind == 'complex':
            self.cumsum= self.type([5,8,9],dtype=self.dtype,**self.kwds)
            self.mul = self.type([50+30j, 24+12j, 6+2j],dtype=self.result_dtype,**self.kwds)
            self.mul_s = self.type([25+10j, 15+6j, 5+2j],dtype=self.result_dtype,**self.kwds)

            self.add = self.type([15+6j, 11+4j, 7+2j],dtype=self.result_dtype,**self.kwds)
            self.add_s = self.type([10+2j, 8+2j, 6+2j],dtype=self.result_dtype,**self.kwds)

            #self.div = [25./68.-15.j/68., 3./10.-3.j/20., 3./20.-1.j/20.]
            self.div = self.type([0.36764705882352941176-0.22058823529411764706j,
                                  0.3-0.15j, 0.15-0.05j],dtype=self.result_dtype,**self.kwds)
            #self.div_s = [25./29.-10.j/29., 15./29.-6.j/29., 5./29.-2.j/29.]
            self.div_s = self.type([0.86206896551724137931-0.34482758620689655172j,
                          0.51724137931034482759-0.20689655172413793103j,
                          0.17241379310344827586-0.06896551724137931034j],
                          dtype=self.result_dtype,**self.kwds)

            #self.rdiv = [2.+6.j/5., 8./3.+4.j/3, 6.+2.j]
            self.rdiv = self.type([2.+1.2j, 2.66666666666666666667+1.33333333333333333333j,
                                   6.+2.j],dtype=self.result_dtype,**self.kwds)
            #self.rdiv_s = [1.+2.j/5., 5./3.+2.j/3., 5.+2.j]
            self.rdiv_s = self.type([1.+0.4j, 1.66666666666666666667+0.666666666666666666667j,
                                     5.+2.j],dtype=self.result_dtype,**self.kwds)

            self.sub = self.type([-5-6j, -5-4j, -5-2j],dtype=self.result_dtype,**self.kwds)
            self.sub_s = self.type([0-2j, -2-2j, -4-2j],dtype=self.result_dtype,**self.kwds)

            self.rsub = self.type([5+6j, 5+4j, 5+2j],dtype=self.result_dtype,**self.kwds)
            self.rsub_s = self.type([0+2j, 2+2j, 4+2j],dtype=self.result_dtype,**self.kwds)

            self.pow1 = self.type([25., 9., 1.],dtype=self.dtype,**self.kwds)
            #self.pow2 = [pow(5,-1.5), pow(3,-1.5), pow(1,-1.5)]
            self.pow2 = self.type([0.08944271909999158786, 0.19245008972987525484,
                                   1.],dtype=self.dtype,**self.kwds)

            self.abs = self.type([5, 3, 1],dtype=self.rdtype,**self.kwds)
            self.real = self.type([5,3,1],dtype=self.rdtype,**self.kwds)
            self.imag = self.type([0, 0, 0],dtype=self.rdtype,**self.kwds)
            self.conj = self.type([5, 3, 1],dtype=self.dtype,**self.kwds)

            self.sum = 9

            self.dot = 80+44j
            self.inner = self.dot
            self.weighted_inner = 68 + 38j

        if self.kind == 'complex' and self.okind == 'real':
            self.cumsum = self.type([5+1j,8+4j,9+9j],dtype=self.dtype,**self.kwds)
            self.mul = self.type([50+10j, 24+24j, 6+30j],dtype=self.result_dtype,**self.kwds)
            self.mul_s = self.type([25+5j, 15+15j, 5+25j],dtype=self.result_dtype,**self.kwds)

            self.add = self.type([15+1j, 11+3j, 7+5j],dtype=self.result_dtype,**self.kwds)
            self.add_s = self.type([10+1j, 8+3j, 6+5j],dtype=self.result_dtype,**self.kwds)

            #self.div = [1./2.+1.j/10., 3./8.+3.j/8., 1./6.+5.j/6.]
            self.div = self.type([0.5+0.1j, 0.375+0.375j,0.16666666666666666667+0.83333333333333333333j],
                                 dtype=self.result_dtype,**self.kwds)
            #self.div_s = [1.+1.j/5., 3./5.+3.j/5., 1./5.+1.j]
            self.div_s = self.type([1.+0.2j, 0.6+0.6j, 0.2+1.j],dtype=self.result_dtype,**self.kwds)

            #self.rdiv = [25./13.-5.j/13., 4./3.-4.j/3., 3./13.-15.j/13.]
            self.rdiv = self.type([1.92307692307692307692-0.38461538461538461538j,
                         1.33333333333333333333-1.33333333333333333333j,
                         0.23076923076923076923-1.15384615384615384615j],
                         dtype=self.result_dtype,**self.kwds)
            #self.rdiv_s = [25./26.-5.j/26., 5./6.-5.j/6., 5./26.-25.j/26.]
            self.rdiv_s = self.type([0.96153846153846153846-0.19230769230769230769j,
                           0.83333333333333333333-0.83333333333333333333j,
                           0.19230769230769230769-0.96153846153846153846j],
                           dtype=self.result_dtype,**self.kwds)

            self.sub = self.type([-5+1j, -5+3j, -5+5j],dtype=self.result_dtype,**self.kwds)
            self.sub_s = self.type([0+1j, -2+3j, -4+5j],dtype=self.result_dtype,**self.kwds)

            self.rsub = self.type([5-1j, 5-3j, 5-5j],dtype=self.result_dtype,**self.kwds)
            self.rsub_s = self.type([0-1j, 2-3j, 4-5j],dtype=self.result_dtype,**self.kwds)

            self.pow1 = self.type([24.+10.j, 0.+18.j, -24.+10.j],dtype=self.dtype,**self.kwds)
            #self.pow2 = [pow(5+1j,-1.5), pow(3+3j,-1.5), pow(1+5j,-1.5)]
            self.pow2 = self.type([0.08307064054041229214-0.0253416052125975132j,
                         0.04379104225017853491-0.1057209281108342370j,
                        -0.04082059235165559671-0.0766590341356157206j],
                         dtype=self.dtype,**self.kwds)

            #self.abs = [pow(26,.5), 3*pow(2,.5), pow(26,.5)]
            self.abs = self.type([5.09901951359278483003,
                        4.24264068711928514641,
                        5.09901951359278483003],dtype=self.rdtype,**self.kwds)
            self.real = self.type([5,3,1],dtype=self.rdtype,**self.kwds)
            self.imag = self.type([1, 3, 5],dtype=self.rdtype,**self.kwds)
            self.conj = self.type([5-1j, 3-3j, 1-5j],dtype=self.dtype,**self.kwds)

            self.sum = 9+9j

            self.dot = 80+64j
            self.inner = 80-64j
            self.weighted_inner = 68- 52j

        if self.kind =='complex' and self.okind =='complex':
            self.cumsum = self.type([5+1j,8+4j,9+9j],dtype=self.dtype,**self.kwds)
            self.mul = self.type([44+40j, 12+36j, -4+32j],dtype=self.result_dtype,**self.kwds)
            self.mul_s = self.type([23+15j, 9+21j, -5+27j],dtype=self.result_dtype,**self.kwds)

            self.add = self.type([15+7j, 11+7j, 7+7j],dtype=self.result_dtype,**self.kwds)
            self.add_s = self.type([10+3j, 8+5j, 6+7j],dtype=self.result_dtype,**self.kwds)

            #self.div = [7./17.-5.j/34., 9./20.+3.j/20., 2./5.+7.j/10.]
            self.div = self.type([0.41176470588235294118-0.14705882352941176471j,
                                  0.45+0.15j, 0.4+0.7j],dtype=self.result_dtype,**self.kwds)
            #self.div_s = [27./29.-5.j/29., 21./29.+9.j/29., 15./29.+23.j/29.]
            self.div_s = self.type([0.93103448275862068966-0.17241379310344827586j,
                          0.72413793103448275862+0.31034482758620689655j,
                          0.51724137931034482759+0.79310344827586206897j],
                          dtype=self.result_dtype,**self.kwds)

            #self.rdiv = [28./13.+10.j/13., 2.-2.j/3., 8./13.-14.j/13.]
            self.rdiv = self.type([2.15384615384615384615+0.76923076923076923077j,
                         2.                    -0.66666666666666666667j,
                         0.61538461538461538462-1.07692307692307692308j],
                         dtype=self.result_dtype,**self.kwds)
            #self.rdiv_s = [27./26.+5.j/26., 7./6.-1.j/2., 15./26.-23.j/26]
            self.rdiv_s = self.type([1.03846153846153846154+0.19230769230769230769j,
                           1.16666666666666666667-0.5j,
                           0.57692307692307692308-0.88461538461538461538j],
                           dtype=self.result_dtype,**self.kwds)

            self.sub = self.type([-5-5j, -5-1j, -5+3j],dtype=self.result_dtype,**self.kwds)
            self.sub_s = self.type([0-1j, -2+1j, -4+3j],dtype=self.result_dtype,**self.kwds)

            self.rsub = self.type([5+5j, 5+1j, 5-3j],dtype=self.result_dtype,**self.kwds)
            self.rsub_s = self.type([0+1j, 2-1j, 4-3j],dtype=self.result_dtype,**self.kwds)

            self.pow1 = self.type([24.+10.j, 0.+18.j, -24.+10.j],dtype=self.dtype,**self.kwds)
            #self.pow2 = [pow(5+1j,-1.5), pow(3+3j,-1.5), pow(1+5j,-1.5)]
            self.pow2 = self.type([0.08307064054041229214-0.0253416052125975132j,
                         0.04379104225017853491-0.1057209281108342370j,
                        -0.04082059235165559671-0.0766590341356157206j],
                         dtype=self.dtype,**self.kwds)

            #self.abs = [pow(26,.5), 3*pow(2,.5), pow(26,.5)]
            self.abs = self.type([5.09901951359278483003,
                        4.24264068711928514641,
                        5.09901951359278483003],dtype=self.rdtype,**self.kwds)
            self.real = self.type([5,3,1],dtype=self.rdtype,**self.kwds)
            self.imag = self.type([1, 3, 5],dtype=self.rdtype,**self.kwds)
            self.conj = self.type([5-1j, 3-3j, 1-5j],dtype=self.dtype,**self.kwds)

            self.sum = 9+9j

            self.dot = 52+108j
            self.inner = 108-20j
            self.weighted_inner= 90 -14j
        self.min = 1
        self.max = 5

    def test_mul(self):
        # Make copies to see we don't overwrite
        acopy = type(self.a)(self.a)
        bcopy = type(self.b)(self.b)
        with self.context:
            # Two of whichever type
            c = self.a * self.b
            self.assertEqual(self.a,acopy)
            self.assertEqual(self.b,bcopy)
            self.assertTrue(self.mul.almost_equal_elem(c,tol=self.tol))

            # Type with scalar
            c = self.a * self.s
            self.assertEqual(self.a,acopy)
            self.assertEqual(self.scalar,self.s)
            self.assertTrue(self.mul_s.almost_equal_elem(c,tol=self.tol))

            # Input that should raise an error
            self.assertRaises(TypeError, self.a.__mul__, self.bad)
            self.assertRaises(ValueError, self.a.__mul__, self.bad2)

    def test_rmul(self):
        # Make copies to see we don't overwrite
        acopy = type(self.a)(self.a)
        bcopy = type(self.b)(self.b)
        with self.context:
            # Two of whichever type
            c = self.a.__rmul__(self.b)
            self.assertEqual(self.a,acopy)
            self.assertEqual(self.b,bcopy)
            self.assertTrue(self.mul.almost_equal_elem(c,tol=self.tol))

            # Type with scalar
            c = self.s * self.a
            self.assertEqual(self.a,acopy)
            self.assertEqual(self.b,bcopy)
            self.assertTrue(self.mul_s.almost_equal_elem(c,tol=self.tol))

            # Input that should raise an error
            self.assertRaises(TypeError, self.a.__rmul__, self.bad)
            self.assertRaises(ValueError, self.a.__rmul__, self.bad2)

    def test_imul(self):
        if not (self.kind == 'real' and self.okind == 'complex'):
            # Make copy to see we don't overwrite
            acopy = type(self.a)(self.a)
            bcopy = type(self.b)(self.b)
            with self.context:
                # Type with itself
                self.a *= self.b
                self.assertEqual(bcopy,self.b)
                self.assertTrue(self.mul.almost_equal_elem(self.a,tol=self.tol))

                # Reset for next test
                self.a = type(self.a)(acopy)
                # Type with scalar
                self.a *= self.s
                self.assertEqual(self.scalar,self.s)
                self.assertTrue(self.mul_s.almost_equal_elem(self.a,tol=self.tol))

                # Input that should raise an error
                self.assertRaises(TypeError, self.a.__imul__, self.bad)
                self.assertRaises(ValueError, self.a.__imul__, self.bad2)

        else:
            with self.context:
                self.assertRaises(TypeError, self.a.__imul__,self.s)
                self.assertRaises(TypeError, self.a.__imul__,self.b)

    def test_add(self):
        # Make copies to see we don't overwrite
        acopy = type(self.a)(self.a)
        bcopy = type(self.b)(self.b)
        with self.context:
            # Type with itself
            c = self.a + self.b
            self.assertEqual(self.a,acopy)
            self.assertEqual(self.b,bcopy)
            self.assertTrue(self.add.almost_equal_elem(c,tol=self.tol))

            # Type with scalar
            c = self.a + self.s
            self.assertEqual(self.a,acopy)
            self.assertEqual(self.scalar,self.s)
            self.assertTrue(self.add_s.almost_equal_elem(c,tol=self.tol))

            # Input that should raise an error
            self.assertRaises(TypeError, self.a.__add__, self.bad)
            self.assertRaises(ValueError, self.a.__add__, self.bad2)


    def test_radd(self):
        # Make copies to see we don't overwrite
        acopy = type(self.a)(self.a)
        bcopy = type(self.b)(self.b)
        with self.context:
            # Type with itself
            c = self.a.__radd__(self.b)
            self.assertEqual(self.a,acopy)
            self.assertEqual(self.b,bcopy)
            self.assertTrue(self.add.almost_equal_elem(c,tol=self.tol))

            # Type with scalar
            c = self.s + self.a
            self.assertEqual(self.a,acopy)
            self.assertEqual(self.scalar,self.s)
            self.assertTrue(self.add_s.almost_equal_elem(c,tol=self.tol))

            # Input that should raise an error
            self.assertRaises(TypeError, self.a.__radd__, self.bad)
            self.assertRaises(ValueError, self.a.__radd__, self.bad2)

    def test_iadd(self):
        if not (self.kind == 'real' and self.okind == 'complex'):
            # Make copy to see we don't overwrite
            acopy = type(self.a)(self.a)
            bcopy = type(self.b)(self.b)
            with self.context:
                # Type with itself
                self.a += self.b
                self.assertEqual(bcopy,self.b)
                self.assertTrue(self.add.almost_equal_elem(self.a,tol=self.tol))

                # Reset for next test
                self.a = type(self.a)(acopy)
                # Type with scalar
                self.a += self.s
                self.assertEqual(self.scalar,self.s)
                self.assertTrue(self.add_s.almost_equal_elem(self.a,tol=self.tol))

                # Input that should raise an error
                self.assertRaises(TypeError, self.a.__iadd__, self.bad)
                self.assertRaises(ValueError, self.a.__iadd__, self.bad2)

        else:
            with self.context:
                self.assertRaises(TypeError, self.a.__iadd__,self.s)
                self.assertRaises(TypeError, self.a.__iadd__,self.b)

    def test_div(self):
        # Make copies to see we don't overwrite
        acopy = type(self.a)(self.a)
        bcopy = type(self.b)(self.b)
        with self.context:
            # Type with itself
            c = self.a / self.b
            self.assertEqual(self.a,acopy)
            self.assertEqual(self.b,bcopy)
            self.assertTrue(self.div.almost_equal_elem(c,tol=self.tol))

            # Type with scalar
            c = self.a / self.s
            self.assertEqual(self.a,acopy)
            self.assertEqual(self.scalar,self.s)
            self.assertTrue(self.div_s.almost_equal_elem(c,tol=self.tol))

            # Input that should raise an error
            self.assertRaises(TypeError, self.a.__div__, self.bad)
            self.assertRaises(ValueError, self.a.__div__, self.bad2)

    def test_rdiv(self):
        # Make copies to see we don't overwrite
        acopy = type(self.a)(self.a)
        bcopy = type(self.b)(self.b)
        with self.context:
            # Type with scalar
            c = self.s / self.a
            self.assertEqual(self.a,acopy)
            self.assertEqual(self.scalar,self.s)
            self.assertTrue(self.rdiv_s.almost_equal_elem(c,tol=self.tol))

            # Input that should raise an error
            self.assertRaises(TypeError, self.a.__rdiv__, self.bad)
            self.assertRaises(ValueError, self.a.__rdiv__, self.bad2)

    def test_idiv(self):
        if not (self.kind == 'real' and self.okind == 'complex'):
            # Make copy to see we don't overwrite
            acopy = type(self.a)(self.a)
            bcopy = type(self.b)(self.b)
            with self.context:
                # Type with itself
                self.a /= self.b
                self.assertEqual(bcopy,self.b)
                self.assertTrue(self.div.almost_equal_elem(self.a,tol=self.tol))

                # Reset for next test
                self.a = type(self.a)(acopy)
                # Type with scalar
                self.a /= self.s
                self.assertEqual(self.scalar,self.s)
                self.assertTrue(self.div_s.almost_equal_elem(self.a,tol=self.tol))

                # Input that should raise an error
                self.assertRaises(TypeError, self.a.__idiv__, self.bad)
                self.assertRaises(ValueError, self.a.__idiv__, self.bad2)

        else:
            with self.context:
                self.assertRaises(TypeError, self.a.__idiv__,self.s)
                self.assertRaises(TypeError, self.a.__idiv__,self.b)

    def test_sub(self):
        # Make copies to see we don't overwrite
        acopy = type(self.a)(self.a)
        bcopy = type(self.b)(self.b)
        with self.context:
            # Type with itself
            c = self.a - self.b
            self.assertEqual(self.a,acopy)
            self.assertEqual(self.b,bcopy)
            self.assertTrue(self.sub.almost_equal_elem(c,tol=self.tol))

            # Type with scalar
            c = self.a - self.s
            self.assertEqual(self.a,acopy)
            self.assertEqual(self.scalar,self.s)
            self.assertTrue(self.sub_s.almost_equal_elem(c,tol=self.tol))

            # Input that should raise an error
            self.assertRaises(TypeError, self.a.__sub__, self.bad)
            self.assertRaises(ValueError, self.a.__sub__, self.bad2)

    def test_rsub(self):
        # Make copies to see we don't overwrite
        acopy = type(self.a)(self.a)
        bcopy = type(self.b)(self.b)
        with self.context:
            # Type with scalar
            c = self.s - self.a
            self.assertEqual(self.a,acopy)
            self.assertEqual(self.scalar,self.s)
            self.assertTrue(self.rsub_s.almost_equal_elem(c,tol=self.tol))

            # Input that should raise an error
            self.assertRaises(TypeError, self.a.__rsub__, self.bad)
            self.assertRaises(ValueError, self.a.__rsub__, self.bad2)

    def test_isub(self):
        if not (self.kind == 'real' and self.okind == 'complex'):
            # Make copy to see we don't overwrite
            acopy = type(self.a)(self.a)
            bcopy = type(self.b)(self.b)
            with self.context:
                # Type with itself
                self.a -= self.b
                self.assertEqual(bcopy,self.b)
                self.assertTrue(self.sub.almost_equal_elem(self.a,tol=self.tol))

                # Reset for next test
                self.a = type(self.a)(acopy)
                # Type with scalar
                self.a -= self.s
                self.assertEqual(self.scalar,self.s)
                self.assertTrue(self.sub_s.almost_equal_elem(self.a,tol=self.tol))

                # Input that should raise an error
                self.assertRaises(TypeError, self.a.__isub__, self.bad)
                self.assertRaises(ValueError, self.a.__isub__, self.bad2)

        else:
            with self.context:
                self.assertRaises(TypeError, self.a.__isub__,self.s)
                self.assertRaises(TypeError, self.a.__isub__,self.b)

    def test_pow(self):
        # Make copy to see we don't overwrite
        acopy = type(self.a)(self.a)
        with self.context:
            # From CPU
            c1 = self.a ** 2
            c2 = self.a ** -1.5
            self.assertEqual(acopy,self.a)
            self.assertTrue(self.pow1.almost_equal_elem(c1,tol=self.tol))
            self.assertTrue(self.pow2.almost_equal_elem(c2,tol=self.tol))

    def test_abs(self):
        # Make copy to see we don't overwrite
        acopy = type(self.a)(self.a)
        # We want to check that absolute value behaves correctly no matter
        # what quadrant it's in.
        t1 = self.a * 1
        t2 = self.a * -1
        t3 = self.a * 1j
        t4 = self.a * -1j
        with self.context:
            c1 = abs(t1)
            c2 = abs(t2)
            c3 = abs(t3)
            c4 = abs(t4)
            self.assertEqual(self.a,acopy)
            # Because complex arrays can involve floating-point math, we
            # must use almost-equal comparisons, esp. on the GPU
            self.assertTrue(self.abs.almost_equal_norm(c1,tol=self.tol))
            self.assertTrue(self.abs.almost_equal_norm(c2,tol=self.tol))
            self.assertTrue(self.abs.almost_equal_norm(c3,tol=self.tol))
            self.assertTrue(self.abs.almost_equal_norm(c4,tol=self.tol))

    def test_real(self):
        # Make copy to see we don't overwrite
        acopy = type(self.a)(self.a)
        with self.context:
            c = self.a.real()
            self.assertEqual(self.a,acopy)
            self.assertEqual(self.real,c)

    def test_imag(self):
        # Make copy to see we don't overwrite
        acopy = type(self.a)(self.a)
        with self.context:
            c = self.a.imag()
            self.assertEqual(self.a,acopy)
            self.assertEqual(self.imag,c)

    def test_conj(self):
        # Make copy to see we don't overwrite
        acopy = type(self.a)(self.a)
        with self.context:
            c = self.a.conj()
            self.assertEqual(self.a,acopy)
            self.assertEqual(self.conj,c)

    def test_cumsum(self):
        # Make copy to see we don't overwrite
        acopy = type(self.a)(self.a)
        with self.context:
            c = self.a.cumsum()
            self.assertEqual(self.a,acopy)
            self.assertTrue(self.cumsum.almost_equal_elem(c,tol=self.tol))

    def test_sum(self):
        # Make copy to see we don't overwrite
        acopy = type(self.a)(self.a)
        with self.context:
            # From CPU
            c = self.a.sum()
            self.assertEqual(self.a,acopy)
            # Hand calculate the relative tolerance for a scalar answer
            self.assertTrue(abs(c-self.sum)<=self.tol*abs(self.sum))

    def test_dot(self):
        # Make copies to see we don't overwrite
        acopy = type(self.a)(self.a)
        bcopy = type(self.b)(self.b)
        with self.context:
            c = self.a.dot(self.b)
            self.assertEqual(self.a,acopy)
            self.assertEqual(self.b,bcopy)
            # Hand calculate the relative tolerance for a scalar answer
            self.assertTrue(abs(c-self.dot)<=self.tol*abs(self.dot))

    def test_inner(self):
        # Make copies to see we don't overwrite
        acopy = type(self.a)(self.a)
        bcopy = type(self.b)(self.b)
        with self.context:
            # CPU with CPU
            c = self.a.inner(self.b)
            self.assertEqual(self.a,acopy)
            self.assertEqual(self.b,bcopy)
            # Hand calculate the relative tolerance for a scalar answer
            self.assertTrue(abs(c-self.inner)<=self.tol*abs(self.inner))

            # Input that should raise an error
            self.assertRaises(TypeError, self.a.inner, self.bad)
            self.assertRaises(ValueError, self.a.inner, self.bad2)

    def test_weighted_inner(self):
        # Make copies to see we don't overwrite
        acopy = type(self.a)(self.a)
        bcopy = type(self.b)(self.b)
        wcopy = type(self.w)(self.w)
        with self.context:
            # CPU with CPU
            c = self.a.weighted_inner(self.b, self.w)
            self.assertEqual(self.a,acopy)
            self.assertEqual(self.b,bcopy)
            self.assertEqual(self.w,wcopy)
            # Hand calculate the relative tolerance for a scalar answer
            self.assertTrue(abs(c-self.weighted_inner)<=self.tol*abs(self.weighted_inner))

            # Input that should raise an error
            self.assertRaises(TypeError, self.a.weighted_inner, self.bad, self.w)
            self.assertRaises(ValueError, self.a.weighted_inner, self.bad2, self.w)

    def test_max(self):
        if self.kind == 'real':
            # Make a copy to see we don't overwrite
            acopy = type(self.a)(self.a)
            with self.context:
                c = self.a.max()
                self.assertEqual(self.a,acopy)
                self.assertEqual(self.max,c)

    def test_min(self):
        if self.kind == 'real':
            # Make a copy to see we don't overwrite
            acopy = type(self.a)(self.a)
            with self.context:
                c = self.a.min()
                self.assertEqual(self.a,acopy)
                self.assertEqual(self.min,c)

    def test_view(self):
        rtypes = { complex64: float32, complex128: float64}
        if self.kind == 'complex':
            rtype = rtypes[self.dtype]
            # Create an array that is the complex array
            # reinterpreted as real
            c_cmp = self.type([5,1,3,3,1,5],dtype=rtype,**self.kwds)
            d_cmp = self.type([5+2j,3+3j,1+5j],dtype=self.dtype,**self.kwds)
            with self.context:
                c = self.a.view(rtype)
                # Check that we correctly created the view
                self.assertEqual(c,c_cmp)
                # That the memory locations are the same
                self.assertEqual(self.a.ptr,c.ptr)
                # And that changing the view changes the original
                c[1] = 2.0
                self.assertEqual(self.a,d_cmp)
