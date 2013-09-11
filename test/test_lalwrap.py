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
These are the unit tests for the custom SWIG wrapping in pycbc.lalwrap.
Note that these unit tests do not attempt to check any of the functions
actually wrapped in lalwrap (those should be tested through unit tests of
the subpackages that call them) but rather test that the SWIG typemaps used
in lalwrap behave as expected and documented.

The only way to do this is to write functions which are specifically designed
for this testing; these are all in the module testlalwrap.  Look in the source
for that module for the definitions of those functions.

Because of the wide variety of possible inputs, outputs, and ways each can be
malformed, there is extensive use in this unittest of automatic creation of
classes and functions using Python's introspection capabilities.  As a result, the
unittests in this script are not necessarily all that easy to read.
"""

import pycbc
import pycbc.scheme
import pycbc.types
import numpy
from numpy import dtype, float32, float64, complex64, complex128
import unittest
from utils import parse_args_all_schemes, simple_exit
import sys
from lal import LIGOTimeGPS as LTG

import pycbc.testlalwrap as tlw

_scheme, _context = parse_args_all_schemes("LALWrap")

# First we define several utility functions that the base typemap tests
# will share.


def DoublePyCBCType(self):
    """
    A convenience function to double all relevant attributes of
    a PyCBC type.  This is used in Input and Noneout typemap tests,
    which by design do the same.  So we do it also in Python, and
    compare.
    """

    returnobj = type(self)(self)
    returnobj *= 2
    if isinstance(returnobj,pycbc.types.TimeSeries):
        returnobj._delta_t *= 2
        returnobj._epoch *= 2
    if isinstance(self,pycbc.types.FrequencySeries):
        returnobj._delta_f *= 2
        returnobj._epoch *= 2

    return returnobj

def baddata(self,pinst,fn,argtuple,kwdict,key):
    """
    A function to ensure that the typemapped function raises
    the appropriate exceptions based on various ways in which
    the '_data' attribute can be erroneous.

    self is an instance of the base test class, and provides
    the 'assertRaises' functionality.

    fn is the function to be tested.

    argtuple is a tuple of arguments to fn

    kwdict is a dict of keyword arguments to fn.  It should
    *not* already contain pinst.

    pinst is the particular instance whose _data property will
    be (repeatedly) modified.

    key is the key with which kwdict will be updated, i.e. a
    call of kwdict.update({key:pinst})
    """

    if not isinstance(pinst,pycbc.types.Array):
        raise TypeError("pinst must be an instance of a PyCBC type")

    otypedict = {dtype('float32'): float64,
                 dtype('float64'): float32,
                 dtype('complex64'): complex128,
                 dtype('complex128'): complex64}

    saved_data = pinst._data.copy()
    pinst._data = None
    kwdict.update({key:pinst})
    self.assertRaises(TypeError,fn,*argtuple,**kwdict)

    pinst._data = numpy.array([saved_data,saved_data],dtype=saved_data.dtype,order='F')
    kwdict.update({key:pinst})
    self.assertRaises(TypeError,fn,*argtuple,**kwdict)

    pinst._data = numpy.array([saved_data,saved_data],dtype=saved_data.dtype)
    kwdict.update({key:pinst})
    self.assertRaises(ValueError,fn,*argtuple,**kwdict)

    pinst._data = numpy.array(saved_data,dtype=otypedict[saved_data.dtype])
    kwdict.update({key:pinst})
    self.assertRaises(ValueError,fn,*argtuple,**kwdict)

    pinst._data = numpy.array(saved_data)

def badepoch(self,pinst,fn,argtuple,kwdict,key):
    """
    A function to ensure that the typemapped function raises
    the appropriate exceptions if the '_epoch' attribute is
    erroneous.

    self is an instance of the base test class, and provides
    the 'assertRaises' functionality.

    fn is the function to be tested.

    argtuple is a tuple of arguments to fn

    kwdict is a dict of keyword arguments to fn.  It should
    *not* already contain pinst.

    pinst is the particular instance whose _epoch property will
    be modified.

    key is the key with which kwdict will be updated, i.e. a
    call of kwdict.update({key:pinst})
    """

    if not isinstance(pinst,pycbc.types.Array):
        raise TypeError("pinst must be an instance of a PyCBC type")

    saved_epoch = pinst._epoch
    pinst._epoch = "A string"
    kwdict.update({key:pinst})
    self.assertRaises(TypeError,fn,*argtuple,**kwdict)
    pinst._epoch = saved_epoch

def baddelta(self,pinst,fn,argtuple,kwdict,key):
    """
    A function to ensure that the typemapped function raises
    the appropriate exceptions if the '_delta_{t,f}' attribute is
    erroneous.

    self is an instance of the base test class, and provides
    the 'assertRaises' functionality.

    fn is the function to be tested.

    argtuple is a tuple of arguments to fn

    kwdict is a dict of keyword arguments to fn.  It should
    *not* already contain pinst.

    pinst is the particular instance whose _delta_{t,f} property will
    be modified.

    key is the key with which kwdict will be updated, i.e. a
    call of kwdict.update({key:pinst})
    """

    if isinstance(pinst,pycbc.types.TimeSeries):
        saved_dt = pinst._delta_t
        pinst._delta_t = "A string"
        kwdict.update({key:pinst})
        self.assertRaises(TypeError,fn,*argtuple,**kwdict)
        pinst._delta_t = saved_dt
    elif isinstance(pinst,pycbc.types.FrequencySeries):
        saved_df = pinst._delta_f
        pinst._delta_f = "A string"
        kwdict.update({key:pinst})
        self.assertRaises(TypeError,fn,*argtuple,**kwdict)
        pinst._delta_f = saved_df
    else:
        raise TypeError("pinst must be an instance of pycbc.type.{Time,Frequency}Series")


# Now we define several 'Base Test' cases, one for each of our four kinds
# of typemaps.  They organize and call the appropriate functions, and perform
# the necessary tests.  The take as input 'self', which should be an instance
# of the TestCase class, and will be passed on to other test functions, as
# well as a 'LALType' string, which should be chosen from the 12 possibilities:
# {REAL4,REAL8,COMPLEX8,COMPLEX16}{Vector,TimeSeries,FrequencySeries}

# First, some more helper functions, to parse these things:

def GetTestFunc(TMType,LALType):
    fnstr = "Test{0}{1}".format(TMType,LALType)
    return getattr(tlw,fnstr)

def GetDtype(LALType):
    if "REAL4" in LALType:
        return float32
    elif "REAL8" in LALType:
        return float64
    elif "COMPLEX8" in LALType:
        return complex64
    elif "COMPLEX16" in LALType:
        return complex128
    else:
        raise ValueError("LALType did not contain a valid string")

def isV(LALType):
    return ("Vector" in LALType)

def isTS(LALType):
    return ("TimeSeries" in LALType)

def isFS(LALType):
    return ("FrequencySeries" in LALType)

def isC(LALType):
    return ("COMPLEX" in LALType)


possible_laltypes = ["REAL4Vector","REAL8Vector",
                     "COMPLEX8Vector","COMPLEX16Vector",
                     "REAL4TimeSeries","REAL8TimeSeries",
                     "COMPLEX8TimeSeries","COMPLEX16TimeSeries",
                     "REAL4FrequencySeries","REAL8FrequencySeries",
                     "COMPLEX8FrequencySeries","COMPLEX16FrequencySeries"]

def ValidLALType(LALType):
    return (LALType in possible_laltypes)

class _BaseTestTMClass(unittest.TestCase):
    """
    This is the base class from which unit tests for all FFT backends
    are derived.
    """
    def setUp(self):
        # Various error messages
        self.tmfail = "SWIG-wrapped function indicated failure"
        self.outfail = "Wrapped function did not modify or generate PyCBC type correctly"
        self.infail = "Wrapped function modified read-only input"
        self.scheme = _scheme
        self.context = _context

    # Now the four typemap tests themselves

    def test_Input(self):
        if not ValidLALType(self.laltype):
            raise ValueError

        if isC(self.laltype):
            self.inbase1 = numpy.array([1+1j,2+2j],dtype=GetDtype(self.laltype))
            self.inbase2 = numpy.array([3+3j,4+4j],dtype=GetDtype(self.laltype))
        else:
            self.inbase1 = numpy.array([1,2],dtype=GetDtype(self.laltype))
            self.inbase2 = numpy.array([3,4],dtype=GetDtype(self.laltype))

        if isV(self.laltype):
            self.input1 = pycbc.types.Array(self.inbase1)
            self.key1 = "invec1"
            self.input2 = pycbc.types.Array(self.inbase2)
            self.key2 = "invec2"
        if isTS(self.laltype):
            self.input1 = pycbc.types.TimeSeries(self.inbase1,delta_t=1.0,epoch=LTG(1,2))
            self.key1 = "ints1"
            self.input2 = pycbc.types.TimeSeries(self.inbase2,delta_t=2.0,epoch=LTG(3,4))
            self.key2 = "ints2"
        if isFS(self.laltype):
            self.input1 = pycbc.types.FrequencySeries(self.inbase1,delta_f=1.0,epoch=LTG(1,2))
            self.key1 = "infs1"
            self.input2 = pycbc.types.FrequencySeries(self.inbase2,delta_f=2.0,epoch=LTG(3,4))
            self.key2 = "infs2"

        self.expectedout = DoublePyCBCType(self.input1)
        self.copy2 = type(self.input2)(self.input2)
        self.fn = GetTestFunc("Input",self.laltype)

        # Now begin our testing.  If we're on the GPU, we only check that we
        # raise TypeError if we call these functions, and that they can
        # be moved back to the CPU and the correct values obtained.
        if self.scheme != 'cpu':
            with self.context:
                intuple=[self.input1,self.input2]
                self.assertRaises(TypeError,self.fn,*intuple)
        else:
            self.retval=self.fn(self.input1,self.input2)
            self.assertEqual(self.retval,0,msg=self.tmfail)
            self.assertEqual(self.expectedout,self.input1,msg=self.outfail)
            self.assertEqual(self.copy2,self.input2,msg=self.infail)

            # Now test that the correct errors are raised when we
            # give erroneous input.  First save a copy so that we
            # can reset things, ensuring only one property of the
            # input at a time is broken.
            baddata(self,self.input1,self.fn,[],{self.key2:self.input2},self.key1)
            baddata(self,self.input2,self.fn,[],{self.key1:self.input1},self.key2)
            if isTS(self.laltype) or isFS(self.laltype):
                badepoch(self,self.input1,self.fn,[],{self.key2:self.input2},self.key1)
                badepoch(self,self.input2,self.fn,[],{self.key1:self.input1},self.key2)
                baddelta(self,self.input1,self.fn,[],{self.key2:self.input2},self.key1)
                baddelta(self,self.input2,self.fn,[],{self.key1:self.input1},self.key2)

    def test_Noneout(self):
        if not ValidLALType(self.laltype):
            raise ValueError

        if isC(self.laltype):
            self.ibase = numpy.array([1+1j,2+2j],dtype=GetDtype(self.laltype))
        else:
            self.ibase = numpy.array([1,2],dtype=GetDtype(self.laltype))
        self.obase = numpy.zeros([2],dtype=GetDtype(self.laltype))

        if isV(self.laltype):
            self.input = pycbc.types.Array(self.ibase)
            self.ikey = "invec"
            self.output = pycbc.types.Array(self.obase)
            self.okey = "outvec"
        if isTS(self.laltype):
            self.input = pycbc.types.TimeSeries(self.ibase,delta_t=1.0,epoch=LTG(1,2))
            self.ikey = "ints"
            self.output = pycbc.types.TimeSeries(self.obase,delta_t=1.0,epoch=LTG(0,0))
            self.okey = "outts"
        if isFS(self.laltype):
            self.input = pycbc.types.FrequencySeries(self.ibase,delta_f=1.0,epoch=LTG(1,2))
            self.ikey = "infs"
            self.output = pycbc.types.FrequencySeries(self.obase,delta_f=1.0,epoch=LTG(0,0))
            self.okey = "outfs"

        self.expectedout = DoublePyCBCType(self.input)
        self.savedinput = type(self.input)(self.input)
        self.fn = GetTestFunc("Noneout",self.laltype)

        # Now begin our testing.  If we're on the GPU, we only check that we
        # raise TypeError if we call these functions, and that they can
        # be moved back to the CPU and the correct values obtained.
        if self.scheme != 'cpu':
            with self.context:
                intuple=[self.input,self.output]
                self.assertRaises(TypeError,self.fn,*intuple)
        else:
            self.retval=self.fn(self.input,self.output)
            self.assertEqual(self.retval,None,msg="Noneout typemap did not return 'None'")
            self.assertEqual(self.expectedout,self.output,msg=self.outfail)
            self.assertEqual(self.savedinput,self.input,msg=self.infail)

            # Now test that the correct errors are raised when we
            # give erroneous input.  First save a copy so that we
            # can reset things, ensuring only one property of the
            # input at a time is broken.
            baddata(self,self.input,self.fn,[],{self.okey:self.output},self.ikey)
            baddata(self,self.output,self.fn,[],{self.ikey:self.input},self.okey)
            if isTS(self.laltype) or isFS(self.laltype):
                badepoch(self,self.input,self.fn,[],{self.okey:self.output},self.ikey)
                badepoch(self,self.output,self.fn,[],{self.ikey:self.input},self.okey)
                baddelta(self,self.input,self.fn,[],{self.okey:self.output},self.ikey)
                baddelta(self,self.output,self.fn,[],{self.ikey:self.input},self.okey)

    def test_Newout(self):
        if not ValidLALType(self.laltype):
            raise ValueError

        if isC(self.laltype):
            self.value = 1.2+3.4j
            self.ibase = numpy.array([self.value,self.value],dtype=GetDtype(self.laltype))
        else:
            self.value = 5.6
            self.ibase = numpy.array([self.value,self.value],dtype=GetDtype(self.laltype))

        self.kwdict = {'length':2,'value':self.value}
        if isV(self.laltype):
            self.expectedout = pycbc.types.Array(self.ibase)
        if isTS(self.laltype):
            self.expectedout = pycbc.types.TimeSeries(self.ibase,delta_t=1.0,epoch=LTG(1,2))
            self.kwdict.update({'epoch':self.expectedout._epoch,
                                'deltaT':self.expectedout._delta_t})
        if isFS(self.laltype):
            self.expectedout = pycbc.types.FrequencySeries(self.ibase,delta_f=1.0,epoch=LTG(1,2))
            self.kwdict.update({'epoch':self.expectedout._epoch,
                                'deltaF':self.expectedout._delta_f})

        self.fn = GetTestFunc("Newout",self.laltype)

        # Now begin our testing.  On the GPU, there's nothing to check
        if self.scheme != 'cpu':
            pass
        else:
            empty_args = []
            self.output=self.fn(*empty_args,**self.kwdict)
            self.assertEqual(self.expectedout,self.output,msg=self.outfail)

    def test_Argout(self):
        if not ValidLALType(self.laltype):
            raise ValueError

        if isC(self.laltype):
            self.value = 1.2+3.4j
            self.ibase = numpy.array([self.value,self.value],dtype=GetDtype(self.laltype))
        else:
            self.value = 5.6
            self.ibase = numpy.array([self.value,self.value],dtype=GetDtype(self.laltype))

        self.kwdict = {'length':2,'value':self.value}
        if isV(self.laltype):
            self.expectedout = pycbc.types.Array(self.ibase)
        if isTS(self.laltype):
            self.expectedout = pycbc.types.TimeSeries(self.ibase,delta_t=1.0,epoch=LTG(1,2))
            self.kwdict.update({'epoch':self.expectedout._epoch,
                                'deltaT':self.expectedout._delta_t})
        if isFS(self.laltype):
            self.expectedout = pycbc.types.FrequencySeries(self.ibase,delta_f=1.0,epoch=LTG(1,2))
            self.kwdict.update({'epoch':self.expectedout._epoch,
                                'deltaF':self.expectedout._delta_f})

        self.fn = GetTestFunc("Argout",self.laltype)

        # Now begin our testing.  On the GPU, there's nothing to check
        if self.scheme != 'cpu':
            pass
        else:
            empty_args = []
            self.output=self.fn(*empty_args,**self.kwdict)
            self.assertEqual(len(self.output),2,msg="Argout typemap did not return a list of values")
            self.assertEqual(self.output[0],0,msg=self.tmfail)
            self.assertEqual(self.expectedout,self.output[1],msg=self.outfail)


# Now, factories to create test cases.

LALWrapTestClasses = []
for laltype in possible_laltypes:
    klass = type('{0}_{1}Test'.format(_scheme,laltype),
                 (_BaseTestTMClass,),
                 {'laltype': laltype})
    LALWrapTestClasses.append(klass)

# Finally, we create suites and run them:

if __name__ == '__main__':

    suite = unittest.TestSuite()
    for klass in LALWrapTestClasses:
            suite.addTest(unittest.TestLoader().loadTestsFromTestCase(klass))

    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
