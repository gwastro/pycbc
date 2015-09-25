# Copyright (C) 2015  Christopher M. Biwer
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
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

import logging
import numpy
import sys
from pycbc import frame

# import dependencies that are not standard to pycbc
import ROOT
ROOT.gSystem.Load('/usr/lib64/libdmtsigp.so')
ROOT.gSystem.Load('/usr/lib64/libgdsplot.so')
from foton import FilterFile, Filter, iir2z

def get_swstat_bits(frame_filenames, swstat_channel_name, start_time, end_time):
    ''' This function just checks the first time in the SWSTAT channel
    to see if the filter was on, it doesn't check times beyond that.

    This is just for a first test on a small chunck of data.

    To read the SWSTAT bits, reference: https://dcc.ligo.org/DocDB/0107/T1300711/001/LIGO-T1300711-v1.pdf

    Bit 0-9 = Filter on/off switches for the 10 filters in an SFM.
    Bit 10 = Filter module input switch on/off
    Bit 11 = Filter module offset switch on/off
    Bit 12 = Filter module output switch on/off
    Bit 13 = Filter module limit switch on/off 
    Bit 14 = Filter module history reset momentary switch
    '''

    # read frames
    swstat = frame.read_frame(frame_filenames, swstat_channel_name,
                      start_time=start_time, end_time=end_time)

    # convert number in channel to binary
    bits = bin(int(swstat[0]))

    # check if filterbank input or output was off
    filterbank_off = False
    if len(bits) < 14 or int(bits[-13]) == 0 or int(bits[-11]) == 0:
        filterbank_off = True

    return bits[-10:], filterbank_off


def filter_data(data, filter_name, filter_file, bits, filterbank_off=False,
                    swstat_channel_name=None):
    '''
    A naive function to determine if the filter was on at the time
    and then filter the data.
    '''

    # if filterbank is off then return a time series of zeroes
    if filterbank_off:
        return numpy.zeros(len(data))

    # loop over the 10 filters in the filterbank
    for i in range(10):

        # read the filter
        filter = Filter(filter_file[filter_name][i])

        # if bit is on then filter the data
        bit = int(bits[-(i+1)])
        if bit:
            logging.info('filtering with filter module %d', i)

            # if there are second-order sections then filter with them
            if len(filter.sections):
                data = filter.apply(data)

            # else it is a filter with only gain so apply the gain
            else:
                coeffs = iir2z(filter_file[filter_name][i])
                if len(coeffs) > 1:
                    logging.info('Gain-only filter module return more than one number')
                    sys.exit()
                gain = coeffs[0]
                data = gain * data

    return  data

def read_gain_from_frames(frame_filenames, gain_channel_name, start_time, end_time):
    '''
    Returns the gain from the file.
    '''

    # get timeseries from frame
    gain = frame.read_frame(frame_filenames, gain_channel_name,
                      start_time=start_time, end_time=end_time)

    return gain[0]
