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
Base class of all data vectors for the pycbc package
"""

from abc import ABCMeta, abstractmethod, abstractproperty

class DataVectorBase:
    
    __metaclass__ = ABCMeta
    
    __slots__ = ("data_vector", "element_size", "element_type", 
                 "meta_type", "length" )
    
    def __init__(self, element_type=0, length=0, element_size=0,
                 meta_type='non_specified_memory', data_vector=0):
        """
        Initialize the DataVectorBase type
        """
        print "DataVectorBase.__init__ called"
        
        self.element_size = element_size     #: element size in bytes
        self.element_type = element_type     #: element type (float, double, int)
        self.meta_type = meta_type                    
        """specifies origin and meta type of the memory used: 
        cbc_memory_meta_type, 
        gpu_cuda_global_memory,
        gpu_cuda_global_memory,
        gpu_cuda_constant_memory,
        gpu_cuda_texture_memory,
        gpu_opencl_global_memory,
        gpu_opencl_constant_memory,
        gpu_opencl_local_memory,
        gpu_opencl_private_memory,
        cpugpu_cuda_zero_latency_memory,
        cpu_generic_memory,
        cpu_pinned_memory,
        cpu_fftw_aligned_memory,
        non_specified_memory
        """
        self.length = length           #: length of data vector
        self.data_vector = data_vector #: true allocated memory via swig from C .. GPU etc
        

