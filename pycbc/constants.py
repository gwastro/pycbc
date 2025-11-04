"""
Constants module for PyCBC.

This module provides constants from various packages (LAL, astropy, numpy)
based on a preferred order of availability. This helps minimize dependency bloat
by relying on smaller, more general packages (astropy/numpy) when the large
LALSuite dependency is not available.

The priority order for constant lookup is: LAL > Astropy / SciPy / NumPy
"""

import logging

# We define a global logger for this module
logger = logging.getLogger('pycbc.constants')

# Attempt to import LALSuite
try:
    import lal
    _LAL_AVAILABLE = True
    logger.debug("LAL constants module found.")
except ImportError:
    lal = None
    _LAL_AVAILABLE = False
    logger.debug("LAL constants module not found.")


import numpy as np
from astropy import (
    constants as aconstants,
    units as aunits
)
# first, do mappings for constants which the fallback is directly in astropy/numpy
# All are in SI units
_FALLBACK_MAPPING = {
    'C_SI': aconstants.c.value,             # Speed of light
    'G_SI': aconstants.G.value,             # Gravitational constant
    'MSUN_SI': aconstants.M_sun.value,      # Mass of the Sun
    'PC_SI': aconstants.pc.value,           # Parsec
    'REARTH_SI': aconstants.R_earth.value,  # Earth equatorial radius
    'YRJUL_SI': aunits.year.to(aunits.s),   # years in seconds
    'PI': np.pi,
    'TWOPI': 2 * np.pi,
    'GAMMA': np.euler_gamma,
}

# We need to define some constants from astropy values
MSUN_SI = _FALLBACK_MAPPING['MSUN_SI']
C_SI = _FALLBACK_MAPPING['C_SI']
G_SI = _FALLBACK_MAPPING['G_SI']

MTSUN_SI = MSUN_SI * G_SI / (C_SI ** 3)
MRSUN_SI = MSUN_SI * G_SI / (C_SI * C_SI)

# Add in these hybrid constants
_FALLBACK_MAPPING.update({
    'MTSUN_SI': MTSUN_SI,
    'MRSUN_SI': MRSUN_SI
})

# Define the Constant Lookup Function
def get_constant(name):
    """
    Retrieves a constant value by name from the preferred order of packages.

    The order of preference is: LAL > Astropy > NumPy

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
        

    if _LAL_AVAILABLE:
        return getattr(lal, name)

    elif name in _FALLBACK_MAPPING:
        return _FALLBACK_MAPPING[name]

    raise NotImplementedError(
        'Contact the PyCBC team, this should never happen. You are '
        'trying to use a LAL constant which doesnt have an astropy/numpy '
        'fallback defined.'
    )

# Expose the constants as attributes of the module

C_SI = get_constant('C_SI')
G_SI = get_constant('G_SI')

MTSUN_SI = get_constant('MTSUN_SI')
MRSUN_SI = get_constant('MTSUN_SI')
MSUN_SI = get_constant('MSUN_SI')

PC_SI = get_constant('PC_SI')
REARTH_SI = get_constant('REARTH_SI')
YRJUL_SI = get_constant('YRJUL_SI')

PI = get_constant('PI')
GAMMA = get_constant('PI')
