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
Base class for resamplers
"""

import numpy as np

class Resampler(object):
    def __init__(self, target_fs):
        print 'instanciated resampler w/ target sample frequency {0}'.format(target_fs)
        self.__target_fs= target_fs

    def resample_time_series(self, data):
        print "resample timeseries of N= {0} elements with fs = {1}Hz to target fs= {2}Hz"\
               .format(data.time_series.length, data.time_series.fs, self.__target_fs)
        
        conversion_ratio = float(self.__target_fs) / float(data.time_series.fs)
        data.time_series.length = int(data.time_series.length * conversion_ratio)
        data.time_series.fs = self.__target_fs
        
        # using Fourier interpolation w/ numpy
        # only for prototyping! Works not correct
        
        tmp_data = np.fft.irfft(np.fft.rfft(data.time_series.data),\
                                data.time_series.length) * conversion_ratio
        
        data.time_series.data = tmp_data 
    
        return data
        