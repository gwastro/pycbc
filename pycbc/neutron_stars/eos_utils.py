# Copyright (C) 2022 Francesco Pannarale, Andrew Williamson,
# Samuel Higginbotham
#
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
"""
Utility functions for handling NS equations of state
"""
import os.path
import numpy as np
from scipy.interpolate import interp1d
import lalsimulation as lalsim
from . import NS_SEQUENCES, NS_DATA_DIRECTORY
from .pg_isso_solver import PG_ISSO_solver


def load_ns_sequence(eos_name):
    """
    Load the data of an NS non-rotating equilibrium sequence generated
    using the equation of state (EOS) chosen by the user.
    File format is: grav mass (Msun), baryonic mass (Msun), compactness

    Parameters
    -----------
    eos_name: string
        NS equation of state label ('2H' is the only supported
        choice at the moment)

    Returns
    ----------
    ns_sequence: numpy.array
        contains the sequence data in the form NS gravitational
        mass (in solar masses), NS baryonic mass (in solar
        masses), NS compactness (dimensionless)
    max_ns_g_mass: float
        the maximum NS gravitational mass (in solar masses) in
        the sequence (this is the mass of the most massive stable
        NS)
    """
    ns_sequence_file = os.path.join(
        NS_DATA_DIRECTORY, 'equil_{}.dat'.format(eos_name))
    if eos_name not in NS_SEQUENCES:
        raise NotImplementedError(
            f'{eos_name} does not have an implemented NS sequence file! '
            f'To implement, the file {ns_sequence_file} must exist and '
            'contain: NS gravitational mass (in solar masses), NS baryonic '
            'mass (in solar masses), NS compactness (dimensionless)')
    ns_sequence = np.loadtxt(ns_sequence_file)
    max_ns_g_mass = max(ns_sequence[:, 0])
    return (ns_sequence, max_ns_g_mass)


def interp_grav_mass_to_baryon_mass(ns_g_mass, ns_sequence):
    """
    Determines the baryonic mass of an NS given its gravitational
    mass and an NS equilibrium sequence (in solar masses).

    Parameters
    -----------
    ns_g_mass: float
        NS gravitational mass (in solar masses)
    ns_sequence: numpy.array
        contains the sequence data in the form NS gravitational
        mass (in solar masses), NS baryonic mass (in solar
        masses), NS compactness (dimensionless)

    Returns
    ----------
    float
    """
    x = ns_sequence[:, 0]
    y = ns_sequence[:, 1]
    f = interp1d(x, y)

    return f(ns_g_mass)


def interp_grav_mass_to_compactness(ns_g_mass, ns_sequence):
    """
    Determines the dimensionless compactness parameter of an NS given
    its gravitational mass and an NS equilibrium sequence.

    Parameters
    -----------
    ns_g_mass: float
        NS gravitational mass (in solar masses)
    ns_sequence: numpy.array
        contains the sequence data in the form NS gravitational
        mass (in solar masses), NS baryonic mass (in solar
        masses), NS compactness (dimensionless)

    Returns
    ----------
    float
    """
    x = ns_sequence[:, 0]
    y = ns_sequence[:, 2]
    f = interp1d(x, y)

    return f(ns_g_mass)


def initialize_eos(ns_mass, eos):
    """Load an equation of state and return the compactness and baryonic
    mass for a given neutron star mass

    Parameters
    ----------
    ns_mass : {float, array}
        The mass of the neutron star, in solar masses.
    eos : str
        Name of the equation of state.

    Returns
    -------
    ns_compactness : float
        Compactness parameter of the neutron star.
    ns_b_mass : float
        Baryonic mass of the neutron star.
    """
    if isinstance(ns_mass, np.ndarray):
        input_is_array = True
    if eos in NS_SEQUENCES:
        ns_seq, ns_max = load_ns_sequence(eos)
        try:
            if any(ns_mass > ns_max) and input_is_array:
                raise ValueError(
                    f'Maximum NS mass for {eos} is {ns_max}, received masses '
                    f'up to {max(ns_mass[ns_mass > ns_max])}')
        except TypeError:
            if ns_mass > ns_max and not input_is_array:
                raise ValueError(
                    f'Maximum NS mass for {eos} is {ns_max}, received '
                    f'{ns_mass}')
        # Interpolate NS compactness and rest mass
        ns_compactness = interp_grav_mass_to_compactness(ns_mass, ns_seq)
        ns_b_mass = interp_grav_mass_to_baryon_mass(ns_mass, ns_seq)
    elif eos in lalsim.SimNeutronStarEOSNames:
        #eos_obj = lalsim.SimNeutronStarEOSByName(eos)
        #eos_fam = lalsim.CreateSimNeutronStarFamily(eos_obj)
        #r_ns = lalsim.SimNeutronStarRadius(ns_mass * lal.MSUN_SI, eos_obj)
        #ns_compactness = lal.G_SI * ns_mass * lal.MSUN_SI / (r_ns * lal.C_SI**2)
        raise NotImplementedError(
            'LALSimulation EOS interface not yet implemented!')
    else:
        raise NotImplementedError(
            f'{eos} is not implemented! Available are: '
            f'{NS_SEQUENCES + list(lalsim.SimNeutronStarEOSNames)}')
    return (ns_compactness, ns_b_mass)


def foucart18(
        eta, ns_compactness, ns_b_mass, bh_spin_mag, bh_spin_pol):
    """Function that determines the remnant disk mass of an NS-BH system
    using the fit to numerical-relativity results discussed in
    `Foucart, Hinderer & Nissanke, PRD 98, 081501(R) (2018)`_.

    .. _Foucart, Hinderer & Nissanke, PRD 98, 081501(R) (2018):
        https://doi.org/10.1103/PhysRevD.98.081501

    Parameters
    ----------
    eta : {float, array}
        The symmetric mass ratio of the system
        (note: primary is assumed to be the BH).
    ns_compactness : {float, array}
        NS compactness parameter.
    ns_b_mass : {float, array}
        Baryonic mass of the NS.
    bh_spin_mag: {float, array}
        Dimensionless spin magnitude of the BH.
    bh_spin_pol : {float, array}
        The tilt angle of the BH spin.
    """
    isso = PG_ISSO_solver(bh_spin_mag, bh_spin_pol)
    # Fit parameters and tidal correction
    alpha = 0.406
    beta  = 0.139
    gamma = 0.255
    delta = 1.761
    fit = (
        alpha / eta ** (1/3) * (1 - 2 * ns_compactness)
        - beta * ns_compactness / eta * isso
        + gamma
        )
    return ns_b_mass * np.where(fit > 0.0, fit, 0.0) ** delta
