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

_aliases = {}
_aliases['mchirp'] = 'chirp_mass_source'
_aliases['mass1'] = 'mass_1_source'
_aliases['mass2'] = 'mass_2_source'
_aliases['snr'] = 'network_matched_filter_snr'
_aliases['z'] = _aliases['redshift'] = 'redshift'
_aliases['distance'] = 'luminosity_distance'


class Merger(object):
    """Informaton about a specific compact binary merger"""
    def __init__(self, name, source='gwtc-1'):
        """ Return the information of a merger

        Parameters
        ----------
        name: str
            The name (GW prefixed date) of the merger event.
        """
        try:
            self.data = catalog.get_source(source)[name]
        except KeyError:
            # Try common name
            data = catalog.get_source(source)
            for mname in data:
                cname = data[mname]['commonName']
                if cname == name:
                    name = mname
                    self.data = data[name]
                    break
            else:
                raise ValueError('Did not find merger matching'
                                 ' name: {}'.format(name))

        # Set some basic params from the dataset
        for key in self.data:
            setattr(self, '_raw_' + key, self.data[key])

        for key in _aliases:
            setattr(self, key, self.data[_aliases[key]])

        self.common_name = self.data['commonName']
        self.time = self.data['GPS']
        self.frame = 'source'

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
        if name in _aliases:
            name = _aliases[name]

        try:
            if return_errors:
                mid = self.data[name]
                high = self.data[name + '_upper']
                low = self.data[name + '_lower']
                return (mid, low, high)
            else:
                return self.data[name]
        except KeyError as e:
            print(e)
            raise RuntimeError("Cannot get parameter {}".format(name))

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
        from pycbc.io import get_file
        from pycbc.frame import read_frame

        for fdict in self.data['strain']:
            if (fdict['detector'] == ifo and fdict['duration'] == duration and
                    fdict['sampling_rate'] == sample_rate and
                    fdict['format'] == 'gwf'):
                url = fdict['url']
                break
        else:
            raise ValueError('no strain data is available as requested '
                             'for ' + self.common_name)

        ver = url.split('/')[-1].split('-')[1].split('_')[-1]
        sampling_map = {4096: "4KHZ",
                        16384: "16KHZ"}
        channel = "{}:GWOSC-{}_{}_STRAIN".format(
                ifo, sampling_map[sample_rate], ver)

        filename = get_file(url, cache=True)
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
        try:
            return self.mergers[key]
        except KeyError:
            # Try common name
            for m in self.mergers:
                if key == self.mergers[m].common_name:
                    break
            else:
                raise ValueError('Did not find merger matching'
                                 ' name: {}'.format(key))
            return self.mergers[m]

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
        v = [self.mergers[m].median1d(param, return_errors=return_errors)
             for m in self.mergers]
        if return_errors:
            value, merror, perror = zip(*v)
            return numpy.array(value), numpy.array(merror), numpy.array(perror)
        else:
            return numpy.array(v)
