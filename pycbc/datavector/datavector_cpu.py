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
package for cpu data vectors
"""


from datavector_base import *
from datavectorcpu import real_vector_t

class DataVectorCpuGeneric(DataVectorBase):
    
    def __init__(self, length=1024):
        print "instanciated DataVectorCpuGeneric"
        
        self.data_vector= real_vector_t(length, 0)#: references by swig to C Layer
       
        # currently overwrites self.data_vector which is obviously wrong!
        # Would call the real_vector_t destructor then
        
        #super(DataVectorCpuGeneric, self).__init__('float', 'real', length, 
        #      4, 'cpu_generic_memory')


class DataVectorCpuFfftwAligned(DataVectorBase):
    
    def __init__(self, length=0):
        print "instanciated DataVectorCpuFfftwAligned"
        
        self.data_vector= [length] #: references by swig to C Layer
        
        super(DataVectorCpuFfftwAligned, self).__init__('float', 'real', length, 
              4, 'cpu_fftw_aligned_memory')
    
    