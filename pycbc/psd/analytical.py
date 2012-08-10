#!/usr/bin/python
# Copyright (C) 2012 Alex Nitz
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
""" Provides reference PSDs from LALSimulation.
"""

from pycbc.types import FrequencySeries, zeros
import lalsimulation

_name_prefix = 'SimNoisePSD'

def get_list():
    " Return a list of available PSD functions. "
    l = []
    for name in lalsimulation.__dict__:
        if name.startswith(_name_prefix) and name != _name_prefix:
            l.append(name[len(_name_prefix):])
    return l

for _name in get_list():
    exec("""
def %s(length, delta_f, low_freq_cutoff):
    " Return a FrequencySeries containing the %s PSD from LALSimulation. "
    return from_lalsimulation(lalsimulation.%s, length, delta_f, low_freq_cutoff)
""" % (_name, _name, _name_prefix + _name))

def from_lalsimulation(func, length, delta_f, low_freq_cutoff):
    " Return a FrequencySeries containing the specified LALSimulation PSD. "
    psd = FrequencySeries(zeros(length), delta_f=delta_f)
    kmin = int(low_freq_cutoff / delta_f)
    for k in xrange(kmin, length - 1):
        psd[k] = func(k * delta_f)
    return psd

