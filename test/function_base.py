# Copyright (C) 2012 Josh Willis, Andrew Miller
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
This runs a function, with all permutations of CPU and GPU arrays for inputs, and checks that they end up in the correct scheme
"""
import pycbc
import numpy
import pycbc.scheme
import pycbc.types

class function_base(object):
    def checkScheme(self, inputs, results, places):
        # First importing the correct dtype the context
        if type(self.context) is pycbc.scheme.CUDAScheme:
            import pycuda
            import pycuda.gpuarray
            from pycuda.gpuarray import GPUArray as SchemeArray
        elif type(self.context) is pycbc.scheme.OpenCLScheme:
            import pyopencl
            import pyopencl.array
            from pyopencl.array import Array as SchemeArray
        else:
            from numpy import ndarray as SchemeArray
        # Check will specify all the arguments that should be tested, and what their values should be.
        for k in range(len(inputs)):
            val = inputs[k]
            trueval = results[k]
            if type(val) is pycbc.types.Array:
                # We will need to check the scheme, and data type
                self.assertTrue(type(val.data)==SchemeArray)
                if type(self.context) != pycbc.scheme.DefaultScheme:
                    self.assertTrue(type(val._scheme)==type(self.context))
                else:
                    self.assertTrue(val._scheme==None)
                self.assertTrue(val.dtype == trueval.dtype)
                
                # These can be almostequal tests, because the child class should perform
                # the numerical tests, this is just to be sure that it isn't garbled
                self.assertEqual(len(val), len(trueval))
                for i in range(len(trueval)):
                    self.assertAlmostEqual(val[i], trueval[i], places=places)
            else:
                self.assertTrue(val==trueval)
        
    def checkCPU(self,inputs, results, places):
        # Check will specify all the arguments that should be tested, and what their values should be.
        for k in range(len(inputs)):
            val = inputs[k]
            trueval = results[k]
            if type(val) is pycbc.types.Array:
                # We will need to check the scheme, and data type
                self.assertTrue(type(val.data)==numpy.ndarray)
                self.assertTrue(val._scheme==None)
                self.assertTrue(val.dtype == trueval.dtype)
                
                # These can be almostequal tests, because the child class should perform
                # the numerical tests, this is just to be sure that it isn't garbled
                self.assertEqual(len(val), len(trueval))
                for i in range(len(trueval)):
                    self.assertAlmostEqual(val[i], trueval[i], places=places)
            else:
                self.assertTrue(val==trueval)

    def scheme_swap(self, function, inputs, results, places, *args, **kwargs):                    
        # We will need to make lists of inputs, involving the CPU and current scheme, covering all possible permutations
        finalargs = []
        for i in range(pow(2,len(inputs))):
            finalargs.append([])
            for j in range(len(inputs)):
                finalargs[i].append(inputs[j] * 1)
                # Now that it has been added to the list (on the CPU) we will move it to the GPU if necessary
                with self.context:
                    if (i/(pow(2,j))%2 == 1):
                        finalargs[i][j]*=1
            # After these inputs, we will tack on all the other arguments
            for arg in args:
                finalargs[i].append(arg)
        # Now we'll enter the context, run and check all these
        with self.context:
            for i in range(pow(2,len(inputs))):
                function(*finalargs[i],**kwargs)
                self.checkScheme(finalargs[i][0:len(inputs)],results,places)
                
    # This function is here so that when different kwargs need to be specified for running on the cpu,
    # the tests can call them correctly.
    def cpu_swap(self, function, inputs, results, places, *args, **kwargs):        
        finalargs = []
        for i in range(pow(2,len(inputs))):
            finalargs.append([])
            for j in range(len(inputs)):
                finalargs[i].append(inputs[j] * 1)
                # Now that it has been added to the list (on the CPU) we will move it to the GPU if necessary
                with self.context:
                    if (i/(pow(2,j))%2 == 1):
                        finalargs[i][j]*=1
            # After these inputs, we will tack on all the other arguments
            for arg in args:
                finalargs[i].append(arg)
        
        for i in range(pow(2,len(inputs))):
            function(*finalargs[i],**kwargs)
            self.checkCPU(finalargs[i][0:len(inputs)],results,places)
        
    
