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
""" This modules contains information about the announced LIGO/Virgo
compact binary mergers
"""
import json
from astropy.utils.data import download_file

# For the time being all quantities are the 1-d median value
# FIXME with posteriors when available and we can just post-process that

# GWTC-1 catalog
gwtc1_url = "https://www.gw-openscience.org/catalog/GWTC-1-confident/filelist/"


def get_source(source):
    """Get the source data for a particular GW catalog
    """
    if source == 'gwtc-1':
        fname = download_file(gwtc1_url, cache=True)
        data = json.load(open(fname, 'r'))
    else:
        raise ValueError('Unkown catalog source {}'.format(source))
    return data['data']
