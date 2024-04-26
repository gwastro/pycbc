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
import logging
import json

from pycbc.io import get_file

logger = logging.getLogger('pycbc.catalog.catalog')

# For the time being all quantities are the 1-d median value
# FIXME with posteriors when available and we can just post-process that

# LVC catalogs
base_lvc_url = "https://www.gwosc.org/eventapi/jsonfull/{}/"

def lvk_catalogs():
    _catalog_source = "https://gwosc.org/eventapi/json/"
    catalog_list = json.load(open(get_file(_catalog_source), 'r'))
    return catalog_list
    
def populate_catalogs():
    """ Refresh set of known catalogs
    """
    global _catalogs
    if _catalogs is None:
        # update the LVK catalog information
         _catalogs = {cname: 'LVK' for cname in lvk_catalogs().keys()}

_catalogs = None

# add some aliases
_aliases = {}
_aliases['gwtc-1'] = 'GWTC-1-confident'
_aliases['gwtc-2'] = 'GWTC-2'
_aliases['gwtc-2.1'] = 'GWTC-2.1-confident'
_aliases['gwtc-3'] = 'GWTC-3-confident'

def list_catalogs():
    """Return a list of possible GW catalogs to query"""
    populate_catalogs()    
    return list(_catalogs.keys())

def get_source(source):
    """Get the source data for a particular GW catalog
    """
    populate_catalogs()

    if source in _aliases:
        source = _aliases[source]

    if source in _catalogs:
        catalog_type = _catalogs[source]
        if catalog_type == 'LVK':
            fname = get_file(base_lvc_url.format(source), cache=True)
            data = json.load(open(fname, 'r'))
    else:
        raise ValueError('Unkown catalog source {}'.format(source))
    return data['events']
