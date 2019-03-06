# Copyright (C) 2017  Alex Nitz
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


#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#
""" This package provides information about LIGO/Virgo detections of
compact binary mergers
"""
import numpy
from . import catalog

class Merger(object):
    """Informaton about a specific compact binary merger"""
    def __init__(self, name, source='gwtc-1'):
        """ Return the information of a merger

        Parameters
        ----------
        name: str
            The name (GW prefixed date) of the merger event.
        """
        self.data = catalog.get_source(source)[name]

        # Set some basic params from the dataset
        for key in self.data:
            if 'best' in self.data[key]:
                setattr(self, key, self.data[key]['best'])

        self.time = self.data['tc']['best']

    def median1d(self, name, return_errors=False):
        """ Return median 1d marginalized parameters

        Parameters
        ----------
        name: str
            The name of the parameter requested
        return_errors: Optional, {bool, False}
            If true, return a second and third parameter that represents the
            lower and upper 90% error on the parameter.

        Returns
        -------
        param: float or tuple
            The requested parameter
        """
        if return_errors:
            mid = self.data[name]['best']
            low, high = self.data[name]['err']
            return (mid, low, high)
        else:
            return self.data[name]['best']

    def strain(self, ifo, duration=32, sample_rate=4096):
        """ Return strain around the event

        Currently this will return the strain around the event in the smallest
        format available. Selection of other data is not yet available.

        Parameters
        ----------
        ifo: str
            The name of the observatory you want strain for. Ex. H1, L1, V1

        Returns
        -------
        strain: pycbc.types.TimeSeries
            Strain around the event.
        """
        from astropy.utils.data import download_file
        from pycbc.frame import read_frame

        # Information is currently wrong on GWOSC!
        # channels = self.data['files']['FrameChannels']
        # for channel in channels:
        #    if ifo in channel:
        #        break

        length = "{}sec".format(duration)
        if sample_rate == 4096:
            sampling = "4KHz"
        elif sample_rate == 16384:
            sampling = "16KHz"

        channel = "{}:GWOSC-{}_R1_STRAIN".format(ifo, sampling.upper())
        url = self.data['files'][ifo][length][sampling]['GWF']
        filename = download_file(url, cache=True)
        return read_frame(str(filename), str(channel))

class Catalog(object):
    """Manage a set of binary mergers"""
    def __init__(self, source='gwtc-1'):
        """ Return the set of detected mergers

        The set of detected mergers. At some point this may have some selection
        abilities.
        """
        self.data = catalog.get_source(source=source)
        self.mergers = {name: Merger(name,
                                     source=source) for name in self.data}
        self.names = self.mergers.keys()

    def __len__(self):
        return len(self.mergers)

    def __getitem__(self, key):
        return self.mergers[key]

    def __setitem__(self, key, value):
        self.mergers[key] = value

    def __delitem__(self, key):
        del self.mergers[key]

    def __iter__(self):
        return iter(self.mergers)

    def median1d(self, param, return_errors=False):
        """ Return median 1d marginalized parameters

        Parameters
        ----------
        name: str
            The name of the parameter requested
        return_errors: Optional, {bool, False}
            If true, return a second and third parameter that represents the
            lower and upper 90% error on the parameter.

        Returns
        -------
        param: nump.ndarray or tuple
            The requested parameter
        """
        v = [self.mergers[m].median1d(param, return_errors=return_errors) for m in self.mergers]
        if return_errors:
            value, merror, perror = zip(*v)
            return numpy.array(value), numpy.array(merror), numpy.array(perror)
        else:
            return numpy.array(v)
