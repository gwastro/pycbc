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
package for opencl data vectors
"""


from datavector_base import *

class DataVectorOpenClGlobal(DataVectorBase):
    
    def __init__(self, length=0):
        print "instanciated DataVectorOpenClGlobal"
        
        self.data_vector= [length] #: references by swig to C Layer 
        
        super(DataVectorOpenClGlobal, self).__init__('float', 'real', length, 
              4, 'gpu_opencl_global_memory')


class DataVectorOpenClPrivate(DataVectorBase):
    
    def __init__(self, length=0):
        print "instanciated DataVectorOpenClPrivate"
        
        self.data_vector= [length] #: references by swig to C Layer
        
        super(DataVectorOpenClPrivate, self).__init__('float', 'real', length, 
              4, 'gpu_opencl_private_memory')
    
    