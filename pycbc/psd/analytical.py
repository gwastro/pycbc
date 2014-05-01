#!/usr/bin/python
# Copyright (C) 2012 Alex Nitz, Tito Dal Canton
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
"""Provides reference PSDs from LALSimulation.
"""

from pycbc.types import FrequencySeries, zeros
import lalsimulation
import numpy

# build a list of usable PSD functions from lalsimulation
_name_prefix = 'SimNoisePSD'
_psd_list = []
for _name in lalsimulation.__dict__:
    if _name.startswith(_name_prefix) and _name != _name_prefix:
        try:
            eval('lalsimulation.' + _name + '(100.)')
        except TypeError:
            # ignore fancy PSDs taking extra args
            continue
        _psd_list.append(_name[len(_name_prefix):])
_psd_list = sorted(_psd_list)

# add functions wrapping lalsimulation PSDs
for _name in _psd_list:
    exec("""
def %s(length, delta_f, low_freq_cutoff):
    \"\"\"Return a FrequencySeries containing the %s PSD from LALSimulation.
    \"\"\"
    return from_lalsimulation(lalsimulation.%s, length, delta_f, low_freq_cutoff)
""" % (_name, _name, _name_prefix + _name))

def get_lalsim_psd_list():
    """Return a list of available reference PSD functions.
    """
    return _psd_list

def from_lalsimulation(func, length, delta_f, low_freq_cutoff):
    """Generate a frequency series containing the specified LALSimulation PSD.

    Parameters
    ----------
    func : function
        LALSimulation PSD function.
    length : int
        Length of the frequency series in samples.
    delta_f : float
        Frequency resolution of the frequency series.
    low_freq_cutoff : float
        Frequencies below this value are set to zero.

    Returns
    -------
    psd : FrequencySeries
        The generated frequency series.
    """
    psd = FrequencySeries(zeros(length), delta_f=delta_f)
    kmin = int(low_freq_cutoff / delta_f)
    psd.data[kmin:] = map(func, numpy.arange(length)[kmin:] * delta_f)
    return psd

def from_string(psd_name, length, delta_f, low_freq_cutoff):
    """Generate a frequency series containing a LALSimulation PSD specified by name.

    Parameters
    ----------
    psd_name : string
        PSD name as found in LALSimulation, minus the SimNoisePSD prefix.
    length : int
        Length of the frequency series in samples.
    delta_f : float
        Frequency resolution of the frequency series.
    low_freq_cutoff : float
        Frequencies below this value are set to zero.
    
    Returns
    -------
    psd : FrequencySeries
        The generated frequency series.
    """
    try:
        func = lalsimulation.__dict__[_name_prefix + psd_name]
    except KeyError:
        raise ValueError(psd_name + ' not found among LALSimulation PSD functions.')
    return from_lalsimulation(func, length, delta_f, low_freq_cutoff)

