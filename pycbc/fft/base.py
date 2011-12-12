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
Base class of Fast Fourier Transform

This base class is not designed to be called directly; rather,
each implementation should derive from it, and therefore must
implement the abstract method perform_transform.  Their __init__()
should call the __init__ of this base class.
"""

from abc import ABCMeta, abstractmethod, abstractproperty

class FastFourierTransformBase(object):
    """
    This is the base class for all FastFourierTransform
    classes in PyCBC.  All particular implementations should derive from
    this class.
    """

    __metaclass__ = ABCMeta

    # Constructor:
    def __init__(self,vector_length,data_type,transform_direction,
                 data_precision,**kwargs):
        self._vector_length = vector_length

        # Note: all of the error checks below need to be
        # converted to pycbc_error calls.

        if data_type not in ['complex','real']:
            raise ValueError("Invalid value for data_type")
        else:
            self._data_type = data_type

        if transform_direction not in ['forward','reverse','backward']:
            raise ValueError("Invalid value for transform_direction")
        else:
            if transform_direction is 'forward':
                self._fwdflag = 1
            else:
                self._fwdflag = 0

        if data_precision not in ['single','double']:
            raise ValueError("Invalid value for data_precision")
        else:
            self._data_precision = data_precision

        super(FastFourierTransformBase, self).__init__(**kwargs)

    @abstractmethod
    def perform_transform(self,input_vector,output_vector,**kwargs):
        # At present, a placeholder that will just perform the
        # following check if invoked by a derived class.  Probably
        # just needs to be 'pass'.
        assert not hasattr(super(FastFourierTransformBase, self), 'perform_transform')
