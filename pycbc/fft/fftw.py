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
                                   measure_level=0,
                                   device_context=MyContext)

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
                        transform_direction='forward',data_precision='single',
                        device_context=Mycontext)
 for invec in vector_list:
      new_output = complex_vector_single_cpu_t(mylen,...)
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

# The following is the dictionary mapping transform parameters to
# their corresponding implementation function. The tuple that is
# the key is (type,precision,fwdflag).
FFTWTransformDict = \
    {('real','single',0) : getattr(cfftw,'execute_real_single_reverse_fft'),
     ('real','single',1) : getattr(cfftw,'execute_real_single_forward_fft'),
     ('real','double',0) : getattr(cfftw,'execute_real_double_reverse_fft'),
     ('real','double',1) : getattr(cfftw,'execute_real_double_forward_fft'),
     ('complex','single',0) : getattr(cfftw,'execute_complex_single_fft'),
     ('complex','single',1) : getattr(cfftw,'execute_complex_single_fft'),
     ('complex','double',0) : getattr(cfftw,'execute_complex_double_fft'),
     ('complex','double',1) : getattr(cfftw,'execute_complex_double_fft')}

# The following is the dictionary mapping transform parameters to
# their corresponding plan constructor. The tuple that is
# the key is (type,precision).
FFTWPlanConstructorDict = \
    {('real','single') : getattr(cfftw,'fftw_real_single_plan'),
     ('real','double') : getattr(cfftw,'fftw_real_double_plan'),
     ('complex','single') : getattr(cfftw,'fftw_complex_single_plan'),
     ('complex','double') : getattr(cfftw,'fftw_complex_double_plan')}

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

        if measure_level not in [0,1,2,3]:
            raise ValueError("Invalid measure_level; must be 0, 1, 2, or 3")

        # The values of all params not specific to FFTW have
        # already been sanity-checked, so create the plan

        ordered_args = (self._data_type,self._data_precision)
        plan_constructor = FFTWPlanConstructorDict[ordered_args]
        self._fftw_plan = plan_constructor(self._vector_length,
                                           self._fwdflag,
                                           measure_level)

    def perform_transform(self,input_vector,output_vector,**kwargs):
        # Create the key, lookup and execute the transform

        ordered_args = (self._data_type,self._data_precision,self._fwdflag)
        transform_function = FFTWTransformDict[ordered_args]
        transform_function(output_vector,input_vector,self._fftw_plan)
