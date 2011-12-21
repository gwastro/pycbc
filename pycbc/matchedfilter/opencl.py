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
MatchedFilter OpenCl implementation class for the pycbc package
"""
from pycbc.matchedfilter.base import MatchedFilterBase

# import modul data structure and processing functions from the clayer
from pycbc.matchedfilter.clayer.matchedfilteropencl import matched_filter_opencl_t
from pycbc.matchedfilter.clayer.matchedfilteropencl import gen_snr_opencl
##### FIXME add to clayer from pycbc.matchedfilter.clayer.opencl import find_max_opencl

# get functions from swig's pointer library to allow
# call by reference from python to c 

##### FIXME add to .i file --> for find max
#from pycbc.matchedfilter.clayer.matchedfilteropencl import copy_ulongp
#from pycbc.matchedfilter.clayer.matchedfilteropencl import copy_doublep
#from pycbc.matchedfilter.clayer.matchedfilteropencl import ulongp_value
#from pycbc.matchedfilter.clayer.matchedfilteropencl import doublep_value
#from pycbc.matchedfilter.clayer.matchedfilteropencl import delete_ulongp
#from pycbc.matchedfilter.clayer.matchedfilteropencl import delete_doublep


# import member datavectors
from pycbc.datavector.opencl import complex_vector_single_opencl_t
from pycbc.datavector.opencl import real_vector_single_opencl_t

from pycbc.opencl import OpenClProcessingObj

import logging

class MatchedFilterOpenCl(MatchedFilterBase, OpenClProcessingObj):

    def __init__(self, context, length=0, delta_x=1):

        self.__logger= logging.getLogger('pycbc.MatchedFilterOpenCl')
        self.__logger.debug("instanciated MatchedFilterOpenCl") 

        super(MatchedFilterOpenCl, self).__init__(context, 
                length, 
                delta_x, 
                snr_vector_t=    real_vector_single_opencl_t, 
                qtilde_vector_t= complex_vector_single_opencl_t, 
                q_vector_t =     complex_vector_single_opencl_t)

        self.__logger.debug("A")                
        # instantiate the matched filter data structure in the clayer
        self._mf_clayer_state = matched_filter_opencl_t(self._context)
        self.__logger.debug("B")        

        # instantiate the matched filter clayer functions from functors
        self._gen_snr_opencl=  GenerateSnrOpenCl()
        self.__logger.debug("C")
#### FIXME        self._find_max_opencl= FindMaximumOpenCl()
                
                

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

        self._gen_snr_opencl(self._context, self._mf_clayer_state, 
                             self._snr, stilde, htilde)
                             # FIXME at gen_snr_cpu : self._q, self._qtilde, 1.0, 1.0) 
        
        return self._snr



    def perform_find_max(self, snr):
        """
        calls the max methode of the derived implementation object
        @rtype:  snr:  DataVectorBase
        @return: snr:  Signal to noise ratio series
        @rtype:  float
        @return: Maximum of snr series
        """
        pass
# FIXME add find max in clayer        return self._find_max_opencl(self._context, snr)

## Functors of the matched filter cpu clayer extension
#
class GenerateSnrOpenCl:
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
    def __call__(self, context, mf_clayer_state, output_snr, stilde, htilde): # FIXME at gen_snr_cpu : , q, qtilde, f_min, sigma_sq ):
        
        gen_snr_opencl(context,  mf_clayer_state, output_snr, stilde, htilde) # FIXME at gen_snr_cpu : , q, qtilde, f_min, sigma_sq)
        return 

#class FindMaximumOpenCl:
#    """
#    functor definition for find_max_cpu()
#    """
#                                     #FIXME use kwargs#
#
#
#    def __init__(self):
#        # status variables etc for/in clayer regrading this function goes here !
#        pass
    
    ## Functor GenSnr 
    #  @param self The object pointer.
    #  @param snr        FIXME describe purpose
    #  @param max        FIXME describe purpose
    #  @param index      FIXME describe purpose
#    def __call__(self, context, snr):
#        max= copy_doublep(0.0)
#        index= copy_ulongp(0)
#        
#        find_max_opencl(context, max, index, snr)
#        
#        max_ret= doublep_value(max)           ## FIXME probably put all pointer alloc/free to constructor
#        index_ret= ulongp_value(index)
#        delete_doublep(max)
#        delete_ulongp(index)
#        return max_ret, index_ret
