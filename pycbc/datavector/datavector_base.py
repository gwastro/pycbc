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


class DataVectorBase:
    
    __slots__ = ("data_vector", "element_size", "element_type", "length" )
    
    def __init__(self, element_size=0, element_type=0, length=0):
        """
        Initialize the DataVectorBase type
        """
        self.element_size = element_size
        self.element_type = element_type
        self.length = length
        
        self.data_vector = 0  # allocate true memory via swig to C .. GPU etc
                              # in derived class
        

