# Copyright (C) 2011 Josh Willis
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
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
FFTW implementation class of the Fast Fourier Transform

This is a derived class of FastFourierTransformBase, which
implements the constructor and performs the transform of
a Fast Fourier Transform using the FFTW3 library.  It is assumed
that the input and output datavectors are aligned memory, as is
the case in PyCBC.

When an instance of this class is instantiated, you must give
the length of the vectors to be transformed, their type, the
direction of the transform, and the precision of the vectors, as
below:

  MyFFT = FastFourierTransformFFTW(vector_length=1024,
                                   data_type='complex',
                                   transform_direction='forward',
                                   data_precision='double',
                                   measure_level=0)

where:

  data_type            must be  'real' or 'complex'
  transform_direction  must be  'forward' or 'reverse'
  data_precision       must be  'single' or 'double'
  measure_level        must be  0, 1, 2, or 3, and defaults to
                                1 if not given

Other keyword arguments may be given, but they will be ignored.

When the FFT is instantiated, an FFTW plan is created, drawing
on system-wide wisdom if present.  Once this plan is created, the
perform_transform() method of the instance may be called repeatedly
on new input and output vectors, so long as they are the same length,
type, precision, and the same direction of transform and level of
planning is desired.  For example:

 mylen=1024
 output_list = []
 MyFFT = FastFourierTransformFFTW(vector_length=mylen,data_type='complex',
                        transform_direction='forward',data_precision='single')
 for invec in vector_list:
      new_output = complex_vector_single_cpu_t(mylen)
      MyFFT.perform_transform(invec,new_output)
      output_list.append(new_output)

to loop over an existing list of input vectors, creating a list
of their FFT's.
"""

from abc import ABCMeta, abstractmethod, abstractproperty
from pycbc.cpu import CpuProcessingObj
from pycbc.fft.base import FastFourierTransformBase

# Need to check that SWIG-wrapped Python code provides '__all__',
# or else provide it ourselves.
import pycbc.fft.clayer.fftw as cfftw

class FastFourierTransformFFTW(FastFourierTransformBase,CpuProcessingObj):
    """
    This is the FFTW implementation of the FastFourierTransform
    base class.  When the constructor is called, a plan is created,
    and that instance may have its perform_transform() method
    invoked repeatedly on different input and output vectors, so
    long as a new plan is not needed.
    """

    # Constructor:
    def __init__(self,measure_level=1,**kwargs):
        # First, call our parent init, since it does most
        # of the heavy-lifting, except for sanity-checking
        # measure_level and creating the plan
        super(FastFourierTransformFFTW, self).__init__(**kwargs)

        # Next, responsibilities specific to FFTW:

        if measure_level is not in [0,1,2,3]:
            raise ValueError("Invalid measure_level; must be 0, 1, 2, or 3")

        # The values of all params not specific to FFTW have
        # already been sanity-checked, so the logic below
        # is safe; exactly one plan-creation 'if' statement
        # will be invoked.

        ordered_args = (self._data_type,self._data_precision)

        if ordered_args == ('real','single'):
            self._fftw_plan = cfftw.fftw_real_single_plan(
                self._vector_length,self._fwdflag,measure_level)
        if ordered_args == ('real','double'):
            self._fftw_plan = cfftw.fftw_real_double_plan(
                self._vector_length,self._fwdflag,measure_level)
        if ordered_args == ('complex','single'):
            self._fftw_plan = cfftw.fftw_complex_single_plan(
                self._vector_length,self._fwdflag,measure_level)
        if ordered_args == ('complex','double'):
            self._fftw_plan = cfftw.fftw_complex_double_plan(
                self._vector_length,self._fwdflag,measure_level)

    def perform_transform(self,input_vector,output_vector,**kwargs):
        # Here's where we do the work; the plan was created in
        # by __init__(). Again, sanity checking was already
        # performed, so exactly one of the six exectute functions
        # gets called.

        ordered_args = (self._data_type,self._data_precision,self._fwdflag)

        if ordered_args == ('real','single',1):
            cfftw.execute_real_single_forward_fft(output_vector,
                                           input_vector,self._fftw_plan)
        if ordered_args == ('real','single',0):
            cfftw.execute_real_single_reverse_fft(output_vector,
                                           input_vector,self._fftw_plan)
        if ordered_args == ('real','double',1):
            cfftw.execute_real_double_forward_fft(output_vector,
                                           input_vector,self._fftw_plan)
        if ordered_args == ('real','double',0):
            cfftw.execute_real_double_reverse_fft(output_vector,
                                           input_vector,self._fftw_plan)
        if ordered_args[0:2] == ('complex','single'):
            cfftw.execute_complex_single_fft(output_vector,input_vector,
                                           self._fftw_plan)
        if ordered_args[0:2] == ('complex','double'):
            cfftw.execute_complex_double_fft(output_vector,input_vector,
                                           self._fftw_plan)

        # Maybe call the method of the parent class, for
        # consistency? Note that this shouldn't do anything

        #super(FastFourierTransformFFTW,self).perform_transform(
        #    input_vector,output_vector,**kwargs)
