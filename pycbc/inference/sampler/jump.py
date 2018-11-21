# Copyright (C) 2018 Alex Nitz Sebastian Khan
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
""" This module contains jump proposals for samplers to pick new parameters.
The format will follow the emcee format.
"""
from __future__ import absolute_import
from emcee.moves import (GaussianMove, StretchMove, WalkMove,
                        KDEMove, DEMove, DESnookerMove)

_jump_proposals = {'stretch': StretchMove,
                   'walk': WalkMove,
                   'de': DEMove,
                   'desnooker':DESnookerMove,
                   'kde': KDEMove,
                   'gaussian':GaussianMove}

def get_jump_from_config(section, cp):
    if not cp.has_option(section, 'jump'):
        return None

    name = cp.get(section, 'jump')

    jump_section = 'jump-{}'.format(name)
    options = {}
    if cp.has_section(jump_section):
        options = cp.items(jump_section)

    kwds = {}
    for key, value in options:
        try:
            kwds[key] = float(value)
        except:
            kwds[key] = value
    return _jump_proposals[name](**kwds)

