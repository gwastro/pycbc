# Copyright (C) 2011 Karsten Wiesner, Josh Willis
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
pycbc processing object base class
"""

from abc import ABCMeta, abstractmethod, abstractproperty

import logging

# All PyCBC classes should ultimately inherit from the following class,
# as its __init__ method intercepts any remaining keyword arguments,
# so that unneeded arguments may safely be passed to __init__
class PyCBCRootClass(object):
    def __init__(self,**kwargs):
        pass

# ------------------- pycbc processing base classes section --------------------

# All processing objects have to inherit from this base class
class PyCbcProcessingObj(PyCBCRootClass):

    __metaclass__ = ABCMeta

    def __init__(self, device_context, **kwargs):
        self._logger= logging.getLogger('pycbc.pycbc_base')
        self._devicecontext = device_context
        super(PyCbcProcessingObj,self).__init__(**kwargs)

    @abstractmethod
    def data_in(self, datavector):
        pass
