# Copyright (C) 2012  Alex Nitz, Andrew Miller, Josh Willis
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
'''
The tests in this file are designed to ensure that the schemes behave as designed,
and that data are moved to and from the GPU as they should, or apprpriate exceptions
are raised.  We only attempt this for two representative functions: one basic
arithemtic operation, and one that should *not* move its data (regardless of scheme).

We do not specifically test that the lalwrapped functions raise exceptions from the
GPU, because that test is done in the test_lalwrap unit tests.
'''
import pycbc
import unittest
from pycbc.types import *
from pycbc.scheme import *
import numpy
from numpy import dtype, float32, float64, complex64, complex128
import lal
from utils import parse_args_all_schemes, simple_exit

_scheme, _context = parse_args_all_schemes("Scheme")

# By importing the current schemes array type, it will make it
# easier to check the  array types later
if isinstance(_context,CUDAScheme):
    import pycuda
    import pycuda.gpuarray
    from pycuda.gpuarray import GPUArray as SchemeArray
elif isinstance(_context,CPUScheme):
    from pycbc.types.aligned import ArrayWithAligned as SchemeArray

from pycbc.types.aligned import ArrayWithAligned as CPUArray


class SchemeTestBase(unittest.TestCase):
    def setUp(self):
        self.context = _context
        self.scheme = _scheme
        # Determine kind (real or complex) from dtype:
        if self.dtype == float32 or self.dtype == float64:
            self.kind = 'real'
        else:
            self.kind = 'complex'
        if self.odtype == float32 or self.odtype == float64:
            self.okind = 'real'
        else:
            self.okind = 'complex'
        # Now set up the arrays we'll need.  We run this from a factory
        # constructor that creates many different instances for the
        # various kind/precision combinations.
        if self.kind == 'real':
            self.a = Array([5,3,1],dtype=self.dtype)
            if self.okind == 'real':
                self.b = Array([10,8,6],dtype=self.odtype)
                self.answer = Array([50,24,6],dtype=self.dtype)
            else:
                self.b = Array([10+6j,8+4j,6+2j],dtype=self.odtype)
                self.answer = Array([50+30j,24+12j,6+2j],dtype=self.odtype)
        else:
            self.a = Array([5+1j,3+3j,1+5j],dtype=self.dtype)
            if self.okind == 'real':
                self.b = Array([10,8,6],dtype=self.odtype)
                self.answer = Array([50+10j,24+24j,6+30j],dtype=self.dtype)
            else:
                self.b = Array([10+6j,8+4j,6+2j],dtype=self.odtype)
                self.answer = Array([44+40j,12+36j,-4+32j],dtype=self.dtype)

    def test_move(self):
        '''
        This test uses the __mul__ special method to verify that arrays are moved
        on and off of the GPU automatically when they should be, and that the _scheme
        property and array types are correct for the executing architecture.
        '''
        # Make some copies
        a1 = type(self.a)(self.a)
        a2 = type(self.a)(self.a)
        b1 = type(self.b)(self.b)
        with self.context:
            # The following should move both of a1 and b1 onto the GPU (if self.context
            # isn't CPU)
            c = a1 * b1
            # Check that the data types are correct
            self.assertTrue(isinstance(a1._data, SchemeArray))
            self.assertTrue(isinstance(b1._data, SchemeArray))
            self.assertTrue(isinstance(c._data, SchemeArray))
            # Check that schemes are correct
            self.assertTrue(isinstance(a1._scheme, type(self.context)))
            self.assertTrue(isinstance(b1._scheme, type(self.context)))
            self.assertTrue(isinstance(c._scheme, type(self.context)))
            # And finally check that the values are correct
            self.assertEqual(a1,self.a)
            self.assertEqual(b1,self.b)
            self.assertEqual(c,self.answer)
            # Now check that nothing about a2 has changed, since it wasn't involved
            # in the computation
            self.assertTrue(isinstance(a2._data, CPUArray))
            self.assertTrue(isinstance(a2._scheme, DefaultScheme))
            self.assertEqual(a2,self.a)

        # Now move back to the CPU, and check that everything is correctly
        # transferred:
        c = a1 * b1
        # Check that schemes are correct
        self.assertTrue(isinstance(a1._scheme, DefaultScheme))
        self.assertTrue(isinstance(b1._scheme, DefaultScheme))
        self.assertTrue(isinstance(c._scheme, DefaultScheme))
        # Check that the data types are correct
        self.assertTrue(isinstance(a1._data, CPUArray))
        self.assertTrue(isinstance(b1._data, CPUArray))
        self.assertTrue(isinstance(c.data, CPUArray))
        # And finally check that the values are correct
        self.assertEqual(a1,self.a)
        self.assertEqual(b1,self.b)
        self.assertEqual(c,self.answer)

    def test_do_not_move(self):
        '''
        This test checks that the __eq__ special method (invoked via the
        '==' operator) does *not* change the scheme or type, since it
        does its comparisons by copying from the CPU to GPU, but should
        leave the original arrays in place, with their data properties and
        schemes unchanged.
        '''
        acopy = type(self.a)(self.a)
        with self.context:
            # Force a move to the GPU by trivially multiplying by one:
            a1 = acopy*1
            a2 = acopy*1
            truth = (a1 == a2)
            # Now verify that nothing moved
            self.assertTrue(isinstance(a1._scheme, type(self.context)))
            self.assertTrue(isinstance(a2._scheme, type(self.context)))
            self.assertTrue(isinstance(a1.data, SchemeArray))
            self.assertTrue(isinstance(a2.data, SchemeArray))

# Now the function that creates our various classes

def scheme_test_maker(dtype,odtype):
    class tests(SchemeTestBase):
        def __init__(self,*args):
            self.dtype = dtype
            self.odtype = odtype
            unittest.TestCase.__init__(self,*args)
    tests.__name__ = _scheme + " " + dtype.__name__ + " with " + odtype.__name__
    return tests

types = [ (float32,[float32,complex64]), (float64,[float64,complex128]),
        (complex64,[complex64,float32]), (complex128,[float64,complex128]) ]

suite = unittest.TestSuite()

ind = 0
for ty,oktype in types:
    for ot in oktype:
        na = 'test' + str(ind)
        vars()[na] = scheme_test_maker(ty,ot)
        suite.addTest(unittest.TestLoader().loadTestsFromTestCase(vars()[na]))
        ind += 1

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)


