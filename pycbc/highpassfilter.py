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
Base class for bandpass filters implemented as ForwardsBackwardsButterworthFilter
"""


class ForwardsBackwardsButterworthFilter(object):
    def __init__(self, param1, param2, param3):
        print 'instanciated bandpassfilter'

    def high_pass(self, data, frequency, attenuation):
        print "highpass filtering data at {0} Hz with {1} dB".format(frequency, attenuation)
        return data