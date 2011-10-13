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
Base class of injector
"""

class Injector(object):
    def __init__(self, length):
        print 'instanciated injector'
        self.__length= length

    def read_injections(self, url):
        """
        Read the injections from a ligolw xml file containing a sim_inspiral table
        """
        print 'reading injections from {0}'.format(url)

    def inject_to_time_series(self, data):
        """
        Create the simulated signals and add them to the time data
        """
        assert self.__length == data.time_series.length, "length must fit to time_series length "
        print 'inject into {0}'.format(data.time_series)
        # a very cheesy injection:
        data.time_series.data[1] = 10
        return data

    def inject_to_ferquency_series(self, data):
        """
        Create the simulated signals and add them to the frequency data
        """
        pass
