#!/usr/bin/python
# Copyright (C) 2012-2016 Alex Nitz, Tito Dal Canton, Leo Singer
#               2022 Shichao Wu
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
"""Provides reference PSDs from LALSimulation and pycbc.psd.analytical_space.

More information about how to use these ground-based detectors' PSD can be
found in the guide about :ref:`Analytic PSDs from lalsimulation`. For
space-borne ones, see `pycbc.psd.analytical_space` module.
"""
import numbers
from pycbc.types import FrequencySeries
from pycbc.psd.analytical_space import (
    analytical_psd_lisa_tdi_1p5_XYZ, analytical_psd_lisa_tdi_2p0_XYZ,
    analytical_psd_lisa_tdi_1p5_AE, analytical_psd_lisa_tdi_1p5_T,
    sh_transformed_psd_lisa_tdi_XYZ)
import lal
import numpy

# build a list of usable PSD functions from lalsimulation
_name_prefix = 'SimNoisePSD'
_name_suffix = 'Ptr'
_name_blacklist = ('FromFile', 'MirrorTherm', 'Quantum', 'Seismic', 'Shot', 'SuspTherm')
_psd_list = []

try:
    import lalsimulation
    for _name in lalsimulation.__dict__:
        if _name != _name_prefix and _name.startswith(_name_prefix) and not _name.endswith(_name_suffix):
            _name = _name[len(_name_prefix):]
            if _name not in _name_blacklist:
                _psd_list.append(_name)
except ImportError:
    pass

_psd_list = sorted(_psd_list)

# add functions wrapping lalsimulation PSDs
for _name in _psd_list:
    exec("""
def %s(length, delta_f, low_freq_cutoff):
    \"\"\"Return a FrequencySeries containing the %s PSD from LALSimulation.
    \"\"\"
    return from_string("%s", length, delta_f, low_freq_cutoff)
""" % (_name, _name, _name))

def get_psd_model_list():
    """ Returns a list of available reference PSD functions.

    Returns
    -------
    list
        Returns a list of names of reference PSD functions.
    """
    return get_lalsim_psd_list() + get_pycbc_psd_list()

def get_lalsim_psd_list():
    """Return a list of available reference PSD functions from LALSimulation.
    """
    return _psd_list

def get_pycbc_psd_list():
    """ Return a list of available reference PSD functions coded in PyCBC.

    Returns
    -------
    list
        Returns a list of names of all reference PSD functions coded in PyCBC.
    """
    pycbc_analytical_psd_list = pycbc_analytical_psds.keys()
    pycbc_analytical_psd_list = sorted(pycbc_analytical_psd_list)
    return pycbc_analytical_psd_list

def from_string(psd_name, length, delta_f, low_freq_cutoff, **kwargs):
    """Generate a frequency series containing a LALSimulation or
    built-in space-borne detectors' PSD specified by name.

    Parameters
    ----------
    psd_name : string
        PSD name as found in LALSimulation (minus the SimNoisePSD prefix)
        or pycbc.psd.analytical_space.
    length : int
        Length of the frequency series in samples.
    delta_f : float
        Frequency resolution of the frequency series.
    low_freq_cutoff : float
        Frequencies below this value are set to zero.
    **kwargs :
        All other keyword arguments are passed to the PSD model.

    Returns
    -------
    psd : FrequencySeries
        The generated frequency series.
    """

    # check if valid PSD model
    if psd_name not in get_psd_model_list():
        raise ValueError(psd_name + ' not found among analytical '
                         'PSD functions.')

    # make sure length has the right type for CreateREAL8FrequencySeries
    if not isinstance(length, numbers.Integral) or length <= 0:
        raise TypeError('length must be a positive integer')
    length = int(length)

    # if PSD model is in LALSimulation
    if psd_name in get_lalsim_psd_list():
        lalseries = lal.CreateREAL8FrequencySeries(
            '', lal.LIGOTimeGPS(0), 0, delta_f, lal.DimensionlessUnit, length)
        try:
            func = lalsimulation.__dict__[
                                        _name_prefix + psd_name + _name_suffix]
        except KeyError:
            func = lalsimulation.__dict__[_name_prefix + psd_name]
            func(lalseries, low_freq_cutoff)
        else:
            lalsimulation.SimNoisePSD(lalseries, 0, func)
        psd = FrequencySeries(lalseries.data.data, delta_f=delta_f)

    # if PSD model is coded in PyCBC
    else:
        func = pycbc_analytical_psds[psd_name]
        psd = func(length, delta_f, low_freq_cutoff, **kwargs)

    # zero-out content below low-frequency cutoff
    kmin = int(low_freq_cutoff / delta_f)
    psd.data[:kmin] = 0

    return psd

def flat_unity(length, delta_f, low_freq_cutoff):
    """ Returns a FrequencySeries of ones above the low_frequency_cutoff.

    Parameters
    ----------
    length : int
        Length of output Frequencyseries.
    delta_f : float
        Frequency step for output FrequencySeries.
    low_freq_cutoff : int
        Low-frequency cutoff for output FrequencySeries.

    Returns
    -------
    FrequencySeries
        Returns a FrequencySeries containing the unity PSD model.
    """
    fseries = FrequencySeries(numpy.ones(length), delta_f=delta_f)
    kmin = int(low_freq_cutoff / fseries.delta_f)
    fseries.data[:kmin] = 0
    return fseries

# dict of analytical PSDs coded in PyCBC
pycbc_analytical_psds = {
    'flat_unity' : flat_unity,
    'analytical_psd_lisa_tdi_1p5_XYZ' : analytical_psd_lisa_tdi_1p5_XYZ,
    'analytical_psd_lisa_tdi_2p0_XYZ' : analytical_psd_lisa_tdi_2p0_XYZ,
    'analytical_psd_lisa_tdi_1p5_AE' : analytical_psd_lisa_tdi_1p5_AE,
    'analytical_psd_lisa_tdi_1p5_T' : analytical_psd_lisa_tdi_1p5_T,
    'sh_transformed_psd_lisa_tdi_XYZ' : sh_transformed_psd_lisa_tdi_XYZ,
}
