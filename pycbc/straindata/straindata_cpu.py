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
Cpu version of strain data class
"""

from straindata_base import StrainDataBase

from pycbc.datavector.datavectorcpu import real_vector_double_t as InitialTimeSeriesDoublePreci
from pycbc.datavector.datavectorcpu import real_vector_single_t as TimeSeriesSinglePreci
from pycbc.datavector.datavectorcpu import complex_vector_single_t as FrequencySeries

class StrainDataCpu(StrainDataBase):

    def __init__(self, segments, length, ifo):
        
        super(StrainDataCpu, self).__init__(segments, length, ifo,
                                            InitialTimeSeriesDoublePreci,
                                            TimeSeriesSinglePreci,
                                            FrequencySeries)

    def render(self):
        pass         # nothing to render in cpu version (data remaines in place)
                                            