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
_conversion_options = ['mass1', 'mass2', 'spin1z', 'spin2z', 'duration',
                       'template_duration', 'mtotal', 'total_mass',
                       'q', 'invq', 'eta', 'chirp_mass','mchirp',
                       'chieff', 'chi_eff', 'effective_spin', 'chi_a']

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
    # These just return things already in the bank
    if parameter == 'mass1':
        return bank['mass1'][:][template_ids]
    if parameter == 'mass2':
        return bank['mass2'][:][template_ids]
    if parameter == 'spin1z':
        return bank['spin1z'][:][template_ids]
    if parameter == 'spin2z':
        return bank['spin2z'][:][template_ids]
    # These things may be in the bank, but if not, we need to calculate
    if parameter in ['template_duration', 'duration']:
        if 'template_duration' in bank:
            return bank['template_duration'][:][template_ids]
        duration = pnutils.get_imr_duration(bank['mass1'][:][template_ids],
                                            bank['mass2'][:][template_ids],
                                            bank['spin1z'][:][template_ids],
                                            bank['spin2z'][:][template_ids],
                                            bank['f_lower'][:][template_ids],
                                            approximant="SEOBNRv4")
        return duration
    # Basic conversions
    if parameter in ['mtotal', 'total_mass']:
        return conv.mtotal_from_mass1_mass2(bank['mass1'][:][template_ids],
                                            bank['mass2'][:][template_ids])
    if parameter == 'q':
        return conv.q_from_mass1_mass2(bank['mass1'][:][template_ids],
                                       bank['mass2'][:][template_ids])

    if parameter == 'invq':
        return 1. / get_bank_param('q', bank, template_ids)

    if parameter == 'eta':
        return conv.eta_from_mass1_mass2(bank['mass1'][:][template_ids],
                                         bank['mass2'][:][template_ids])


    if parameter in ['mchirp', 'chirp_mass']:
        return conv.mchirp_from_mass1_mass2(bank['mass1'][:][template_ids],
                                            bank['mass2'][:][template_ids])

    if parameter in ['chieff','chi_eff','effective_spin']:
        return conv.chi_eff(bank['mass1'][:][template_ids],
                            bank['mass2'][:][template_ids],
                            bank['spin1z'][:][template_ids],
                            bank['spin2z'][:][template_ids])

    if parameter == 'chi_a':
        return conv.chi_a(bank['mass1'][:][template_ids],
                          bank['mass2'][:][template_ids],
                          bank['spin1z'][:][template_ids],
                          bank['spin2z'][:][template_ids])

    # parameter not in the current conversion parameter list
    raise NotImplementedError("Bank conversion function " + parameter
                              + " not recognised: choose from '" +
                              "', '".join(_conversion_options) + "'.")    
