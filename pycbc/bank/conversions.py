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
                       'q', 'invq', 'eta', 'chirp_mass', 'mchirp',
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
    # These just give things already in the bank
    if parameter == 'mass1':
         values = bank['mass1'][:][template_ids]
    elif parameter == 'mass2':
         values = bank['mass2'][:][template_ids]
    elif parameter == 'spin1z':
         values = bank['spin1z'][:][template_ids]
    elif parameter == 'spin2z':
         values = bank['spin2z'][:][template_ids]
    # These things may be in the bank, but if not, we need to calculate
    elif parameter in ['template_duration', 'duration']:
        if 'template_duration' in bank:
             values = bank['template_duration'][:][template_ids]
        duration = pnutils.get_imr_duration(bank['mass1'][:][template_ids],
                                            bank['mass2'][:][template_ids],
                                            bank['spin1z'][:][template_ids],
                                            bank['spin2z'][:][template_ids],
                                            bank['f_lower'][:][template_ids],
                                            approximant="SEOBNRv4")
         values = duration
    # Basic conversions
    elif parameter in ['mtotal', 'total_mass']:
         values = conv.mtotal_from_mass1_mass2(bank['mass1'][:][template_ids],
                                            bank['mass2'][:][template_ids])
    elif parameter == 'q':
         values = conv.q_from_mass1_mass2(bank['mass1'][:][template_ids],
                                       bank['mass2'][:][template_ids])

    elif parameter == 'invq':
         values = 1. / bank_conversion('q', bank, template_ids)

    elif parameter == 'eta':
         values = conv.eta_from_mass1_mass2(bank['mass1'][:][template_ids],
                                         bank['mass2'][:][template_ids])


    elif parameter in ['mchirp', 'chirp_mass']:
         values = conv.mchirp_from_mass1_mass2(bank['mass1'][:][template_ids],
                                            bank['mass2'][:][template_ids])

    elif parameter in ['chieff','chi_eff','effective_spin']:
         values = conv.chi_eff(bank['mass1'][:][template_ids],
                            bank['mass2'][:][template_ids],
                            bank['spin1z'][:][template_ids],
                            bank['spin2z'][:][template_ids])

    elif parameter == 'chi_a':
         values = conv.chi_a(bank['mass1'][:][template_ids],
                          bank['mass2'][:][template_ids],
                          bank['spin1z'][:][template_ids],
                          bank['spin2z'][:][template_ids])

    else:
        # parameter not in the current conversion parameter list
        raise NotImplementedError("Bank conversion function " + parameter
                                  + " not recognised: choose from '" +
                                  "', '".join(_conversion_options) + "'.")
    return values
