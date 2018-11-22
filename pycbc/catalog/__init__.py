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
    def __init__(self, name):
        """ Return the information of a merger

        Parameters
        ----------
        name: str
            The name (GW prefixed date) of the merger event.
        """
        self.data = catalog.data[name]

        # Set some basic params from the dataset
        for key in self.data['median1d']:
            setattr(self, key, self.data['median1d'][key][0])

        self.time = self.data['time']

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
            return self.data['median1d'][name]
        else:
            return self.data['median1d'][name][0]

    def strain(self, ifo):
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

        channel = '%s:LOSC-STRAIN' % ifo
        url = self.data['frames'][ifo]
        filename = download_file(url, cache=True)
        return read_frame(filename, channel)

class Catalog(object):
    """Manage a set of binary mergers"""
    def __init__(self):
        """ Return the set of detected mergers

        The set of detected mergers. At some point this may have some selection
        abilities.
        """
        self.mergers = {name: Merger(name) for name in catalog.data}
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
