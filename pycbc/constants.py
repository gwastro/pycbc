"""
Constants module for PyCBC.

This module provides constants from various packages (LAL, astropy, numpy)
based on a choice set by environment variable.

The constant choice is set by setting the PYCBC_CONSTANT_SOURCE environment
variable. Allowed values are 'default' or 'lal'.
"""

import logging
import os

import numpy as np
from astropy import (
    constants as aconstants,
    units as aunits
)

from pycbc.libutils import import_optional

lal = import_optional('lal')

# We define a global logger for this module
logger = logging.getLogger('pycbc.constants')

# Get the environment variable which defines which constants to use
# Allowed values are 'default' or 'lal'
_CONSTANTS = os.environ.get('PYCBC_CONSTANT_SOURCE', 'default').lower()


# first, do mappings for constants where the value is directly in astropy/numpy
# All are in SI units
_DEFAULT_MAPPING = {
    'C_SI': aconstants.c.value,             # Speed of light
    'G_SI': aconstants.G.value,             # Gravitational constant
    'MSUN_SI': aconstants.M_sun.value,      # Mass of the Sun
    'PC_SI': aconstants.pc.value,           # Parsec
    'REARTH_SI': aconstants.R_earth.value,  # Earth equatorial radius
    'YRJUL_SI': aunits.year.to(aunits.s),   # years in seconds
    'PI': np.pi,
    'TWOPI': 2 * np.pi,
    'PI_4': np.pi / 4,
    'GAMMA': np.euler_gamma,
    'LN2': np.log(2.)
}

# We need to define some constants from astropy values
MSUN_SI = _DEFAULT_MAPPING['MSUN_SI']
C_SI = _DEFAULT_MAPPING['C_SI']
G_SI = _DEFAULT_MAPPING['G_SI']

MTSUN_SI = MSUN_SI * G_SI / (C_SI ** 3)
MRSUN_SI = MSUN_SI * G_SI / (C_SI * C_SI)

# Add in these hybrid constants
_DEFAULT_MAPPING.update({
    'MTSUN_SI': MTSUN_SI,
    'MRSUN_SI': MRSUN_SI
})

# Define the Constant Lookup Function
def get_constant(name):
    """
    Retrieves a constant value by name from the preferred order of packages.

    The order of preference is: Astropy + NumPy > LAL

    Parameters
    ----------
    name : str
        The name of the constant to retrieve (e.g., 'C_SI' for the speed of light).

    Returns
    -------
    float or object
        The value of the constant.

    Raises
    ------
    NotImplementedError
        If the constant is not found in any of the available packages.
    """
    if _CONSTANTS.lower() == 'lal': # Allow LAL
        if lal is None:
            raise ImportError(
                "PYCBC_CONSTANT_SOURCE is set to 'lal', but the 'lal' module is not installed. "
                "Please install the 'lal' package to use LAL constants."
            )
        return getattr(lal, name)

    elif name in _DEFAULT_MAPPING:
        return _DEFAULT_MAPPING[name]

    raise NotImplementedError(
        'Contact the PyCBC team, this should never happen. You are '
        'trying to use a constant which is not defined.'
    )

# Expose the constants as attributes of this module:
for const_name in _DEFAULT_MAPPING:
    globals()[const_name] = get_constant(const_name)
