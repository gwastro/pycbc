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
conversion_options = ['mass1', 'mass2', 'spin1z', 'spin2z', 'duration',
                      'template_duration', 'mtotal', 'total_mass',
                      'q', 'invq', 'eta', 'chirp_mass', 'mchirp',
                      'chieff', 'chi_eff', 'effective_spin', 'chi_a']


mass_conversions = {
    'mtotal': conv.mtotal_from_mass1_mass2,
    'total_mass': conv.mtotal_from_mass1_mass2,
    'q': conv.q_from_mass1_mass2,
    'invq': conv.invq_from_mass1_mass2,
    'eta': conv.eta_from_mass1_mass2,
    'mchirp': conv.mchirp_from_mass1_mass2,
    'chirp_mass': conv.mchirp_from_mass1_mass2,
}

spin_conversions = {
    'chieff': conv.chi_eff,
    'chi_eff': conv.chi_eff,
    'effective_spin': conv.chi_eff,
    'chi_a': conv.chi_a
}


def get_bank_property(parameter, bank, template_ids,
                      duration_approximant="SEOBNRv4"):
    """ Get a specific value from a hdf file object in standard PyCBC
    template bank format

    Parameters
    ----------
    parameter: str
        the parameter to convert to, must be in conversion_options

    bank: h5py File object or dictionary of arrays
        Template bank containing the parameters for use in conversions
        must contain mass1, mass2, spin1z, spin2z as a minimum

    template_ids: numpy array
        Array of template IDs for reading a set of templates from the bank

    Optional Parameters
    -------------------
    duration_approximant: string
        The approximant used to calculate the duration of the template if
        not already given

    Returns
    -------
    values: numpy array, same size as template_ids
        Array of whatever the requested parameter is calculated for
        the specified templates in the bank

    """
    # These just give things already in the bank
    if parameter in bank:
        values = bank[parameter][:][template_ids]

    # These things may be in the bank, but if not, we need to calculate
    elif parameter in ['template_duration', 'duration']:
        if 'template_duration' in bank:
            # This statement should be the reached only if 'duration'
            # is given, but 'template_duration' is in the bank
            values = bank['template_duration'][:][template_ids]
        else:
            values = pnutils.get_imr_duration(bank['mass1'][:][template_ids],
                                              bank['mass2'][:][template_ids],
                                              bank['spin1z'][:][template_ids],
                                              bank['spin2z'][:][template_ids],
                                              bank['f_lower'][:][template_ids],
                                              approximant=duration_approximant)
    # Basic conversions
    elif parameter in mass_conversions.keys():
        values = mass_conversions[parameter](bank['mass1'][:][template_ids],
                                             bank['mass2'][:][template_ids])

    elif parameter in spin_conversions.keys():
        values = spin_conversions[parameter](bank['mass1'][:][template_ids],
                                             bank['mass2'][:][template_ids],
                                             bank['spin1z'][:][template_ids],
                                             bank['spin2z'][:][template_ids])

    else:
        # parameter not in the current conversion parameter list
        raise NotImplementedError("Bank conversion function " + parameter
                                  + " not recognised: choose from '" +
                                  "', '".join(conversion_options) + "'.")
    return values
