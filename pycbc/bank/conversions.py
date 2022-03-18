# Copyright (C) 2022 Gareth Cabourn Davies
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
"""
This module is supplied to make a convenience function for converting into 
specific values from PyCBC template banks.
"""

from pycbc import conversions as conv
from pycbc import pnutils

# Convert from parameter name to helper function
# some multiple names are used for the same function
_bank_conversion_functions = {
    'mass1': get_mass1,
    'mass2': get_mass2,
    'spin1z': get_spin1z,
    'spin2z': get_spin2z,
    'duration': get_template_duration,
    'template_duration': get_template_duration,
    'mtotal': get_mtotal,
    'total_mass': get_mtotal,
    'q': get_q,
    'invq': get_invq,
    'eta': get_eta,
    'chirp_mass': get_mchirp,
    'mchirp': get_mchirp,
    'chieff': get_chi_eff,
    'chi_eff': get_chi_eff,
    'effective_spin': get_chi_eff,
    'chi_a': get_chi_a,
}

_conversion_options = list(_bank_conversion_functions.keys())

def bank_conversion(parameter, bank, template_ids):
    """
    Get a specific value from a hdf file object in standard PyCBC
    template bank format

    Parameters:
    -----------
    parameter: str
        the parameter to convert to, must be in _conversion_options

    bank: h5py File object or dictionary of arrays
        Template bank containing the parameters for use in conversions
        must contain mass1, mass2, spin1z, spin2z as a minimum
    """
    if parameter not in _conversion_options:
        raise NotImplementedError("Bank conversion function " + parameter
                                  + " not recognised: choose from '" + 
                                  "', '".join(_conversion_options) + "'.")
    return _bank_conversion_functions[parameter](bank, template_ids)

# These are functions to allow a generalisable function to just return
# things already in the bank
def get_mass1(bank, template_ids):
    return bank['mass1'][:][template_ids]


def get_mass2(bank, template_ids):
    return bank['mass2'][:][template_ids]


def get_spin1z(bank, template_ids):
    return bank['spin1z'][:][template_ids]


def get_spin2z(bank, template_ids):
    return bank['spin2z'][:][template_ids]


# These things may be in the bank, but if not, we need to calculate
def get_template_duration(bank, template_ids):
    if 'template_duration' in bank:
        return bank['template_duration'][:][template_ids]
    else:
        duration = pnutils.get_imr_duration(bank['mass1'][:][template_ids],
                                            bank['mass2'][:][template_ids],
                                            bank['spin1z'][:][template_ids],
                                            bank['spin2z'][:][template_ids],
                                            bank['f_lower'][:][template_ids],
                                            approximant="SEOBNRv4")
 
        return duration


# Basic conversions
def get_mtotal(bank, template_ids):
    return conv.mtotal_from_mass1_mass2(bank['mass1'][:][template_ids],
                                        bank['mass2'][:][template_ids])


def get_q(bank, template_ids):
    return conv.q_from_mass1_mass2(bank['mass1'][:][template_ids],
                                   bank['mass2'][:][template_ids])


def get_invq(bank, template_ids):
    return 1. / get_q(bank, template_ids)


def get_eta(bank, template_ids):
    return conv.eta_from_mass1_mass2(bank['mass1'][:][template_ids],
                                     bank['mass2'][:][template_ids])


def get_mchirp(bank, template_ids):
    return conv.mchirp_from_mass1_mass2(bank['mass1'][:][template_ids],
                                        bank['mass2'][:][template_ids])


def get_chi_eff(bank, template_ids):
    return conv.chi_eff(bank['mass1'][:][template_ids],
                        bank['mass2'][:][template_ids],
                        bank['spin1z'][:][template_ids],
                        bank['spin2z'][:][template_ids])


def get_chi_a(bank, template_ids)
    return conv.chi_a(bank['mass1'][:][template_ids],
                      bank['mass2'][:][template_ids],
                      bank['spin1z'][:][template_ids],
                      bank['spin2z'][:][template_ids])
