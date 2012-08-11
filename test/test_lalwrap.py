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
"""

import pycbc
import pycbc.scheme
import pycbc.types
import numpy
from numpy import dtype, float32, float64, complex64, complex128
import unittest
import base_test
import sys
from sys import getrefcount as grc
from lal import LIGOTimeGPS as LTG

import pycbc.testlalwrap as tlw

import optparse
from optparse import OptionParser

_parser = OptionParser()

def _check_scheme(option, opt_str, scheme, parser):
    if scheme=='cuda' and not pycbc.HAVE_CUDA:
        raise optparse.OptionValueError("CUDA not found")

    if scheme=='opencl' and not pycbc.HAVE_OPENCL:
        raise optparse.OptionValueError("OpenCL not found")
    setattr (parser.values, option.dest, scheme)

_parser.add_option('--scheme','-s', action='callback', type = 'choice', choices = ('cpu','cuda','opencl'), 
                    default = 'cpu', dest = 'scheme', callback = _check_scheme,
                    help = 'specifies processing scheme, can be cpu [default], cuda, or opencl')

_parser.add_option('--device-num','-d', action='store', type = 'int', dest = 'devicenum', default=0,
                    help = 'specifies a GPU device to use for CUDA or OpenCL, 0 by default')

(_opt_list, _args) = _parser.parse_args()

#Changing the optvalues to a dict makes them easier to read
_options = vars(_opt_list)

# We will need these imported in order to check that things are in the current scheme
if _options['scheme'] == 'cuda':
    import pycuda
    import pycuda.gpuarray
elif _options['scheme'] == 'opencl':
    import pyopencl
    import pyopencl.array

# First we define several utility functions that the base typemap tests
# will share.

def ourcmp(inst1,inst2):
    """
    A utility function to compare two pycbc.types.  The usual numpy ==
    statement returns an array of booleans, which is not what we want,
    and also does not complain if dtypes do not match.  It is of course
    unaware of the Time/FrequencySeries metadata.  We could make something
    like this a method of pycbc.types, but the numpy behavior might seem
    more natural to users, so we keep this here.
    """
    if (type(inst1) != type(inst2)):
        return False
    if isinstance(inst1,pycbc.types.TimeSeries):
        return ( (inst1.dtype==inst2.dtype) and (inst1._delta_t == inst2._delta_t)
                 and (inst1._epoch == inst2._epoch) and
                 bool((inst1._data == inst2._data).all()))
    if isinstance(inst1,pycbc.types.FrequencySeries):
        return ( (inst1.dtype==inst2.dtype) and (inst1._delta_f == inst2._delta_f)
                 and (inst1._epoch == inst2._epoch) and
                 bool((inst1._data == inst2._data).all()))
    if isinstance(inst1,pycbc.types.Array):
        return ( (inst1.dtype==inst2.dtype) and
                 bool((inst1._data == inst2._data).all()))

def ourcopy(other):
    """
    A convenience function to return an exact copy of a pycbc.type instance.
    """
    if not isinstance(other,pycbc.types.Array):
        raise TypeError("ourcopy() can only be used to duplicate PyCBC types")
    if isinstance(other,pycbc.types.TimeSeries):
        return pycbc.types.TimeSeries(initial_array=other._data,dtype=other.dtype,
                                      delta_t=other._delta_t,epoch=other._epoch,copy=True)
    if isinstance(other,pycbc.types.FrequencySeries):
        return pycbc.types.FrequencySeries(initial_array=other._data,dtype=other.dtype,
                                      delta_f=other._delta_f,epoch=other._epoch,copy=True)
    if isinstance(other,pycbc.types.Array):
        return pycbc.types.Array(initial_array=other._data,dtype=other.dtype,copy=True)

def saverefcnt(other):
    """
    A convenience function to save a copy of the reference count on all pertinent
    properties of a PyCBC type.  Useful for debugging typemaps.
    """
    if not isinstance(other,pycbc.types.Array):
        raise TypeError("saverefcnt() can only be used to monitor PyCBC types")
    refdict = {}
    refdict.update({'_data':grc(other._data)})
    if isinstance(other,pycbc.types.TimeSeries):
        refdict.update({'_epoch':grc(other._epoch)})
        refdict.update({'_delta_t':grc(other._delta_t)})
    if isinstance(other,pycbc.types.FrequencySeries):
        refdict.update({'_epoch':grc(other._epoch)})
        refdict.update({'_delta_f':grc(other._delta_f)})
    return (grc(other),refdict)

def cmprefcnt(cmptuple,other):
    """
    A convenience function to compare a saved copy of the reference counts of
    pertinent properties of a PyCBC type with the current reference count of
    the PyCBC type instance 'other'.  Useful for debugging typemaps.
    """
    if not isinstance(other,pycbc.types.Array):
        raise TypeError("cmprefcnt() can only be used to monitor PyCBC types")
    trutharray = [cmptuple[0]==grc(other)]
    refdict = cmptuple[1]
    keys = refdict.keys()
    for key in keys:
        trutharray.append(refdict[key]==grc(getattr(other,key)))
    return bool(numpy.array(trutharray).all())

def DoublePyCBCType(self):
    """
    A convenience function to double all relevant attributes of 
    a PyCBC type.  This is used in Input and Noneout typemap tests,
    which by design do the same.  So we do it also in Python, and
    compare.
    """
    self *= 2
    if isinstance(self,pycbc.types.TimeSeries):
        self._delta_t *= 2
        self._epoch *= 2
    if isinstance(self,pycbc.types.FrequencySeries):
        self._delta_f *= 2
        self._epoch *= 2

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

    otypedict = {float32: float64,
                 float64: float32,
                 complex64: complex128,
                 complex128: complex64}

    saved_data = pinst._data.copy()
    pinst._data = None
    kwdict.update({key:pinst})
    self.assertRaises(TypeError,fn,*argtuple,**kwdict)

    pinst._data = numpy.array(saved_data,dtype=pinst.dtype,order='F')
    kwdict.update({key:pinst})
    self.assertRaises(TypeError,fn,*argtuple,**kwdict)

    pinst._data = numpy.array([saved_data,saved_data],dtype=pinst.dtype)
    kwdict.update({key:pinst})
    self.assertRaises(ValueError,fn,*argtuple,**kwdict)

    pinst._data = numpy.array(saved_data,dtype=otypedict[pinst.dtype])
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

def GetPlaces(LALType):
    if GetDtype(LALType) in [float32,complex64]:
        return 7
    else:
        return 15

possible_laltypes = ["REAL4Vector","REAL8Vector",
                     "COMPLEX8Vector","COMPLEX16Vector",
                     "REAL4TimeSeries","REAL8TimeSeries",
                     "COMPLEX8TimeSeries","COMPLEX16TimeSeries",
                     "REAL4FrequencySeries","REAL8FrequencySeries",
                     "COMPLEX8FrequencySeries","COMPLEX16FrequencySeries"]

def ValidLALType(LALType):
    return (LALType in possible_laltypes)
                       


class _BaseTestTMClass(base_test.function_base):
    """
    This is the base class from which unit tests for all FFT backends
    are derived.
    """
    def setUp(self):
        # Various error messages
        self.tmfail = "SWIG-wrapped function indicated failure"
        self.referror = "Reference count incorrectly changed after function call"
        self.outfail = "Wrapped function did not modify or generate PyCBC type correctly"
        self.infail = "Wrapped function modified read-only input"

    # Now the four typemap tests themselves

    def test_Input(self):
        if not ValidLALType(self.laltype):
            raise ValueError

        # Only used for GPU checks:
        self.places = GetPlaces(self.laltype)

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

        self.expectedout = DoublePyCBCType(ourcopy(self.input1))
        self.copy2 = ourcopy(self.input2)
        self.fn = GetTestFunc("Input",self.laltype)

        # Now begin our testing.  If we're on the GPU, we only check that we 
        # raise TypeError if we call these functions, and that they can
        # be moved back to the CPU and the correct values obtained.
        if _options['scheme'] != 'cpu':
            self.cpu_test(self.fn,(self.input1,self.input2),(self.expectedout,self.copy2),self.places)
            with self.context:
                intuple=[self.input1,self.input2]
                self.assertRaises(TypeError,self.fn,*intuple)
        else:
            self.ref1 = saverefcnt(self.input1)
            self.ref2 = saverefcnt(self.input2)
            self.retval=self.fn(self.input1,self.input2)
            self.assertEqual(self.retval,0,msg=self.tmfail)
            self.assertTrue(cmprefcnt(self.ref1,self.input1),msg=self.referror)
            self.assertTrue(cmprefcnt(self.ref2,self.input2),msg=self.referror)
            self.assertTrue(ourcmp(self.expectedout,self.input1),msg=self.outfail)
            self.assertTrue(ourcmp(self.copy2,self.input2),msg=self.infail)

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

        # Only used for GPU checks:
        self.places = GetPlaces(self.laltype)

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
            self.key = "infs"
            self.output = pycbc.types.FrequencySeries(self.obase,delta_f=1.0,epoch=LTG(0,0))
            self.okey = "outfs"

        self.expectedout = DoublePyCBCType(ourcopy(self.input))
        self.savedinput = ourcopy(self.input)
        self.fn = GetTestFunc("Noneout",self.laltype)


        # Now begin our testing.  If we're on the GPU, we only check that we 
        # raise TypeError if we call these functions, and that they can
        # be moved back to the CPU and the correct values obtained.
        if _options['scheme'] != 'cpu':
            self.cpu_test(self.fn,(self.input,self.output),
                          (self.savedinput,self.expectedout),self.places)
            with self.context:
                intuple=[self.input,self.output]
                self.assertRaises(TypeError,self.fn,*intuple)
        else:
            self.iref = saverefcnt(self.input)
            self.oref = saverefcnt(self.output)
            self.retval=self.fn(self.input,self.output)
            self.assertEqual(self.retval,None,msg="Noneout typemap did not return 'None'")
            self.assertTrue(cmprefcnt(self.iref,self.input),msg=self.referror)
            self.assertTrue(cmprefcnt(self.oref,self.output),msg=self.referror)
            self.assertTrue(ourcmp(self.expectedout,self.output),msg=self.outfail)
            self.assertTrue(ourcmp(self.savedinput,self.input),msg=self.infail)

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

        # Only used for GPU checks:
        self.places = GetPlaces(self.laltype)

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
        if _options['scheme'] != 'cpu':
            pass
        else:
            empty_args = []
            self.output=self.fn(*empty_args,**self.kwdict)
            self.assertTrue(ourcmp(self.expectedout,self.output),msg=self.outfail)

    def test_Argout(self):
        if not ValidLALType(self.laltype):
            raise ValueError

        # Only used for GPU checks:
        self.places = GetPlaces(self.laltype)

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
        if _options['scheme'] != 'cpu':
            pass
        else:
            empty_args = []
            self.output=self.fn(*empty_args,**self.kwdict)
            self.assertEqual(len(self.output),2,msg="Argout typemap did not return a list of values")
            self.assertEqual(self.output[0],0,msg=self.tmfail)
            self.assertTrue(ourcmp(self.expectedout,self.output[1]),msg=self.outfail)


# Now, factories to create test cases for each available backend.
# The automation means that the default for each scheme will get created
# and run twice, once as 'Default' and once under its own name.

if _options['scheme']=='cpu':
    CPUTestClasses = []
    context = pycbc.scheme.DefaultScheme()
    for laltype in possible_laltypes:
        klass = type('CPU_{0}Test'.format(laltype),
                     (_BaseTestTMClass,unittest.TestCase),
                     {'laltype': laltype,'context' : context})
        CPUTestClasses.append(klass)

if _options['scheme']=='cuda':
    CUDATestClasses = []
    context = pycbc.scheme.CUDAScheme(device_num=_options['devicenum'])
    for laltype in possible_laltypes:
        CUDATestClasses.append(type('CUDA_{0}Test'.format(laltype),
                                    (_BaseTestTMClass,unittest.TestCase),
                                    {'laltype': laltype,'context' : context}))

if _options['scheme']=='opencl':
    OpenCLTestClasses = []
    context = pycbc.scheme.OpenCLScheme(device_num=_options['devicenum'])
    for laltype in possible_laltypes:
        OpenCLTestClasses.append(type('OpenCL_{0}Test'.format(laltype),
                                      (_BaseTestTMClass,unittest.TestCase),
                                      {'laltype': laltype,'context' : context}))

# Finally, we create suites and run them:

if __name__ == '__main__':

    if _options['scheme']=='cpu':
        suiteCPU = unittest.TestSuite()
        for klass in CPUTestClasses:
            suiteCPU.addTest(unittest.TestLoader().loadTestsFromTestCase(klass))
        results = unittest.TextTestRunner(verbosity=2).run(suiteCPU)

    if _options['scheme']=='cuda':
        suiteCUDA = unittest.TestSuite()
        for klass in CUDATestClasses:
            suiteCUDA.addTest(unittest.TestLoader().loadTestsFromTestCase(klass))
        results = unittest.TextTestRunner(verbosity=2).run(suiteCUDA)

    if _options['scheme']=='opencl':
        suiteOpenCL = unittest.TestSuite()
        for klass in OpenCLTestClasses:
            suiteOpenCL.addTest(unittest.TestLoader().loadTestsFromTestCase(klass))
        results = unittest.TextTestRunner(verbosity=2).run(suiteOpenCL)
        
    NotImpErrors = 0
    for error in results.errors:
        for errormsg in error:
            if type(errormsg) is str:
                if 'NotImplemented' in errormsg:
                    NotImpErrors +=1
                    break
    if results.wasSuccessful():
        sys.exit(0)
    elif len(results.failures)==0 and len(results.errors)==NotImpErrors:
        sys.exit(1)
    else:
        sys.exit(2)
