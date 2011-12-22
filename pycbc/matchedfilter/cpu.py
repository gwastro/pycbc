# Copyright (C) 2011 Karsten Wiesner
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
MatchedFilter Cpu implementation class for the pycbc package
"""

from pycbc.matchedfilter.base import MatchedFilterBase

# import modul data structure and processing functions from the clayer
from pycbc.matchedfilter.clayer.matchedfiltercpu import matched_filter_cpu_t 
from pycbc.matchedfilter.clayer.matchedfiltercpu import gen_snr_cpu
from pycbc.matchedfilter.clayer.matchedfiltercpu import find_max_cpu

# get functions from swig's pointer library to allow
# call by reference from python to c 
from pycbc.matchedfilter.clayer.matchedfiltercpu import copy_ulongp
from pycbc.matchedfilter.clayer.matchedfiltercpu import copy_doublep
from pycbc.matchedfilter.clayer.matchedfiltercpu import ulongp_value
from pycbc.matchedfilter.clayer.matchedfiltercpu import doublep_value
from pycbc.matchedfilter.clayer.matchedfiltercpu import delete_ulongp
from pycbc.matchedfilter.clayer.matchedfiltercpu import delete_doublep

# import member datavectors
from pycbc.datavector.cpu import complex_vector_single_cpu_t
from pycbc.datavector.cpu import real_vector_single_cpu_t

from pycbc.cpu import CpuProcessingObj

from pycbc.fft.fftw import FastFourierTransformFFTW 

import logging

class MatchedFilterCpu(MatchedFilterBase, CpuProcessingObj):

    def __init__(self, context, length=0, delta_x=1):

        self.__logger= logging.getLogger('pycbc.MatchedFilterCpu')
        self.__logger.debug("instanciated MatchedFilterCpu")
                
        super(MatchedFilterCpu, self).__init__(context, 
                                length, 
                                delta_x,
                                snr_vector_t=    real_vector_single_cpu_t, 
                                qtilde_vector_t= complex_vector_single_cpu_t, 
                                q_vector_t =     complex_vector_single_cpu_t)

        # instantiate the matched filter data structure in the clayer
        self._mf_clayer_state = matched_filter_cpu_t(self._context)
        
        # instantiate the matched filter clayer functions from functors
        self._gen_snr_cpu=  GenerateSnrCpu()
        self._find_max_cpu= FindMaximumCpu()

        # instantiate the fft to transform q into qtilde in clayer/gen_snr
        self._fft_qtilde= FastFourierTransformFFTW(vector_length=self._length,
                                                   data_type='complex',
                                                   transform_direction='reverse',
                                                   data_precision='single',
                                                   measure_level=1,
                                                   device_context=self._context)

        self._fft_plan= self._fft_qtilde._fftw_plan

        self.__logger.debug("instanciated FFT {0} with plan {1}".format(self._fft_qtilde, self._fft_plan))


    # implementation of ABC's abstractmethod
    def perform_generate_snr(self, stilde, htilde):
        """
        calls the generate_snr methode of the derived implementation object
        @type  context: Device Context
        @param context: Input: Device Context
        @type  stilde: DataVectorBase
        @param stilde: Input: Straindata frequency domain
        @type  htilde: DataVectorBase
        @param htilde: Input: Template waveform frequency domain
        @rtype  snr:   DataVectorBase
        @return snr:   Output: Signal to noise ratio series
        """
        
        stilde= self.data_in(stilde)
        htilde= self.data_in(htilde)

        self._gen_snr_cpu(self._context, self._snr, stilde, htilde, 
                                         self._q, self._qtilde, self._fft_plan,
                                         1.0, 1.0) 

        return self._snr


    # implementation of ABC's abstractmethod
    def perform_find_max(self, snr):
        """
        calls the max methode of the derived implementation object
        @rtype:  snr:  DataVectorBase
        @return: snr:  Signal to noise ratio series
        @rtype:  float
        @return: Maximum of snr series
        """
        return self._find_max_cpu(self._context, snr)
        


## Functors of the matched filter cpu clayer extension
#
class GenerateSnrCpu:
    """
    functor definition for gen_snr_cpu()
    """                                           #FIXME use kwargs

    def __init__(self):
        # status variables etc for/in clayer regrading this function goes here !
        pass
    
    ## Functor GenSnr 
    #  @param self The object pointer.
    #  @param output_snr FIXME describe purpose
    #  @param stilde     FIXME describe purpose
    #  @param htilde     FIXME describe purpose
    #  @param q          FIXME describe purpose
    #  @param qtilde     FIXME describe purpose
    #  @param f_min      FIXME describe purpose
    #  @param sigma_sq   FIXME describe purpose
    def __call__(self, context, output_snr, stilde, htilde, q, qtilde, fft_plan, f_min, sigma_sq ):
        
        gen_snr_cpu(context, output_snr, stilde, htilde, q, qtilde, fft_plan, f_min, sigma_sq)
        return 

class FindMaximumCpu:
    """
    functor definition for find_max_cpu()
    """
                                     #FIXME use kwargs


    def __init__(self):
        # status variables etc for/in clayer regrading this function goes here !
        pass
    
    ## Functor GenSnr 
    #  @param self The object pointer.
    #  @param snr        FIXME describe purpose
    #  @param max        FIXME describe purpose
    #  @param index      FIXME describe purpose
    def __call__(self, context, snr):
        max= copy_doublep(0.0)
        index= copy_ulongp(0)
        
        find_max_cpu(context, max, index, snr)
        
        max_ret= doublep_value(max)           ## FIXME probably put all pointer alloc/free to constructor
        index_ret= ulongp_value(index)
        delete_doublep(max)
        delete_ulongp(index)
        return max_ret, index_ret
