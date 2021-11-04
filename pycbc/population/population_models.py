# Copyright (C) 2021  Shichao Wu, Alex Nitz, Collin Capano
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
This module provides functions for star formation rate, merger rate density,
and population models of BBH/BNS/NSBH.
"""

import numpy as np
import scipy.integrate as scipy_integrate
import scipy.interpolate as scipy_interpolate
from scipy.misc import derivative as scipy_derivative
from sympy import symbols, sqrt, exp, log, integrate, lambdify
from pycbc.cosmology import get_cosmology, ComovingVolInterpolator
from astropy import units


z = symbols('z')


def sfr_grb_2008(z):
    r""" The star formation rate (SFR) calibrated by high-z GRBs data.

    Parameters
    ----------
    z : float
        The redshift.

    Returns
    -------
    rho_z : float
            The SFR at redshift z, in unit of "Msolar/Mpc^3/yr".

    Note
    ----
        Please see Eq.(5) in <arXiv:0804.4008> for more details.
    """

    rho_local = 0.02  # Msolar/yr/Mpc^3
    a = 3.4
    b = -0.3
    c = -3.5
    eta = -10
    B = 5000
    C = 9

    rho_z = rho_local*((1+z)**(a*eta) + ((1+z)/B)**(b*eta) +
                ((1+z)/C)**(c*eta))**(1./eta)

    return rho_z


def sfr_madau_dickinson_2014(z):
    r""" The madau-dickinson 2014 star formation rate (SFR).

    Parameters
    ----------
    z : float
        The redshift.

    Returns
    -------
    rho_z : float
            The SFR at redshift z, in unit of "Msolar/Mpc^3/yr".

    Notes
    -----
         Pease see Eq.(15) in <arXiv:1403.0007> for more details.
    """

    rho_z = 0.015 * (1+z)**2.7 / (1 + ((1+z)/2.9)**5.6)

    return rho_z


def sfr_madau_fragos_2017(z, k_imf=0.66, mode='high'):
    r""" The madau-fragos 2017 star formation rate (SFR),
         which updates madau-dickinson 2014 SFR by better reproducing
         a number of recent 4 < z < 10 results.

    Parameters
    ----------
    z : float
        The redshift.
    k_imf : float
        The correction factor KIMF adjusts the SFR for the assumed IMF,
        for the Salpeter IMF, k_imf=1.0, for the  three component broken
        power-law Kroupa IMF, k_imf=0.66, here we choose Kroupa IMF as default.
    model : string
        The model of SFR, choose from 'high' and 'low'. Default to 'high'.

    Returns
    -------
    rho_z : float
            The SFR at redshift z, in unit of "Msolar/Mpc^3/yr".

    Notes
    -----
         Pease see <arXiv:1606.07887> and <arXiv:1706.07053> for more details.
    """

    if mode == 'low':
        rho_z = k_imf * 0.015 * (1+z)**2.6 / (1 + ((1+z)/3.2)**6.2)
    elif mode == 'high':
        rho_z = k_imf * 0.015 * (1+z)**2.7 / (1 + ((1+z)/3.0)**5.35)
    else:
        raise ValueError("'mode' must choose from 'high' or 'low'.")

    return rho_z


H_0 = 67.8 * (3.0856776E+19)**(-1) / (1/24/3600/365*1e-9)  # Gyr^-1
Omega_Lambda = 0.692
Omega_m = 0.308
c = 3e10  # cm/s
R_0 = 35
M_Sun = 1.9891e33  # g


def derivative_lookback_time(z):
    r""" The derivative of lookback time t(z)
            with respect to redshit z.

    Parameters
    ----------
    z : float
         The redshift.

    Returns
    -------
    dt_dz : float
            The value of dt/dz ate redshift z.

    Notes
    -----
         Pease see Eq.(A3) in <arXiv:2011.02717v3> for more details.
    """

    dt_dz = 1/H_0/(1+z)/sqrt((Omega_Lambda+Omega_m*(1+z)**3))
    return dt_dz


def p_tau(tau, td_model="inverse"):
    r""" The probability distribution of the time delay.

    Parameters
    ----------
    tau : float
         The merger delay time from the
         formation of the binary system and the orbital
         decay timescale through gravitational wave radiation.
    td_model : str
         The time delay model.

    Returns
    -------
    p_t : float
          The probability at time delay tau.

    Notes
    -----
         Pease see the Appendix in <arXiv:2011.02717v3> for more details.
    """

    if td_model == "log_normal":
        t_LN = 2.9  # Gyr
        sigma_LN = 0.2
        p_t = exp(-(log(tau)-log(t_LN))**2/(2*sigma_LN**2)) / \
             (sqrt(2*np.pi)*sigma_LN)
    elif td_model == "gaussian":
        t_G = 2  # Gyr
        sigma_G = 0.3
        p_t = exp(-(tau-t_G)**2/(2*sigma_G**2)) / (sqrt(2*np.pi)*sigma_G)
    elif td_model == "power_law":
        alpha_t = 0.81
        p_t = tau**(-alpha_t)
    elif td_model == "inverse":
        p_t = tau**(-0.999)  # Try to avoid dividing zero.
    else:
        raise ValueError("'model' must choose from \
                ['log_normal', 'gaussian', 'power_law', 'inverse'].")

    return p_t


def convolution_converter(z_min, sfr_func, derivative_lookback_time, td_model):
    r""" This function is used in a symbolic integral, which to calculate
        the merger rate density of CBC sources. This function converts the
        convolution of the star formation rate SFR(tau) and the time delay
        probability P(tau) on the time delay 'tau' into the convolution on
        the redshift 'z'.

    Parameters
    ----------
    z_min : float
            The redshift.
    sfr_func : function
            The star formation rate function used in the convolution.
    derivative_lookback_time : function
            The derivative of lookback time t(z)
            with respect to redshit z.
    td_model : str
            The name of time delay model.

    Returns
    -------
    f : sympy.core.symbol.Symbol
          The product of SFR(z), P(tau(z)) and dt(z)/dz.

    Notes
    -----
         Pease see Eq.(A2) in <arXiv:2011.02717v3> for more details.
    """

    tau = integrate(derivative_lookback_time(z), (z, z_min, z))
    f = sfr_func(z) * p_tau(tau, td_model) * derivative_lookback_time(z)

    return f


def merger_rate_density(z_array, sfr_func, td_model, rho_local):
    r""" This function uses the symbolic integral to calculate
        the merger rate density of CBC sources. This function converts the
        convolution of the star formation rate SFR(tau) and the time delay
        probability P(tau) on the time delay 'tau' into the convolution on
        the redshift 'z'. This function relies on `convolution_converter`.

    Parameters
    ----------
    z_array : numpy.array
            The array of redshift.
    sfr_func : function
            The star formation rate function used in the convolution.
    td_model : str
            The name of time delay model.
    rho_local : float
            The local merger rate of a certain type of CBC source.

    Returns
    -------
    rho_z : scipy.interpolate.interp1d
          The merger rate density.

    Notes
    -----
         Pease see Eq.(A2) in <arXiv:2011.02717v3> for more details.
    """

    f_z = np.zeros(len(z_array))

    # Note: this 'for' loop is very slow.
    for i in range(len(z_array)):
        f = lambdify(z, convolution_converter(z_min=z_array[i], sfr_func=sfr_func,
                    derivative_lookback_time=derivative_lookback_time, 
                    td_model=td_model), 'scipy')
        f_z[i] = scipy_integrate.quad(f, z_array[i], np.inf, epsabs=1.49e-3)[0]

    f_z = f_z/f_z[0]*rho_local  # Normalization & Rescale
    rho_z = scipy_interpolate.interp1d(z_array, f_z)

    return rho_z


def coalescence_rate_per_z_bin(z_array, sfr_func, td_model, rho_local):
    r""" This function calculates the coalescence rate per redshift bin.

    Parameters
    ----------
    z_array : numpy.array
            The array of redshift.
    sfr_func : function
            The star formation rate function used in the convolution.
    td_model : str
            The name of time delay model.
    rho_local : float
            The local merger rate of a certain type of CBC source.

    Returns
    -------
    dR_dz : numpy.array
          The coalescence rate per redshift bin.

    Notes
    -----
         Pease see Eq.(1) in <arXiv:2011.02717v3> for more details.
    """

    density_func = merger_rate_density(z_array, sfr_func, td_model, rho_local)

    R_list = []
    for redshift in z_array:
        R_z = rate_from_redshift(z=redshift, rate_density=density_func)
        R_list.append(R_z)
    cumulative_rate_func = scipy_interpolate.interp1d(
                    z_array, R_list, fill_value='extrapolate')

    dR_dz = []
    for redshift in z_min:
        dR_dz.append(scipy_derivative(cumulative_rate_func, redshift, dx=1e-6))

    return np.array(dR_dz)


def norm_redshift_distribution(z_array, sfr_func, td_model, rho_local):
    r""" This function calculates the normalized redshift distribution
        of a certain type of CBC sources.

    Parameters
    ----------
    z_array : numpy.array
            The array of redshift.
    sfr_func : function
            The star formation rate function used in the convolution.
    td_model : str
            The name of time delay model.
    rho_local : float
            The local merger rate of a certain type of CBC source.

    Returns
    -------
    morm_dR_dz : numpy.array
            The normalized redshift distribution.

    Notes
    -----
         The can be used as a population-informed prior for redshift
         and luminosity distance of CBC sources.
    """

    dR_dz = coalescence_rate_per_z_bin(z_array, sfr_func, td_model, rho_local)
    density_func = merger_rate_density(z_array, sfr_func, td_model, rho_local)
    lambd = rate_from_redshift(z=z_array[-1], rate_density=density_func)
    morm_dR_dz = np.array(dR_dz)/lambd

    return morm_dR_dz


def rate_from_redshift(z, rate_density, **kwargs):
    """Total rate of occurances out to some redshift

    Parameters
    ----------
    z : float
        The redshift.
    rate_density : function
        The rate density as a function of redshift.
    \**kwargs :
        All other keyword args are passed to :py:func:`get_cosmology` to
        select a cosmology. If none provided, will use
        :py:attr:`DEFAULT_COSMOLOGY`.

    Returns
    -------
    rate: float
        The total rate of occurances out to some redshift.
    """

    cosmology = get_cosmology(**kwargs)
    def diff_rate(z):
        dr = cosmology.differential_comoving_volume(z) / (1+z)
        return (dr * 4 * np.pi * rate_density(z)).value
    return scipy_integrate.quad(diff_rate, 0, z)[0]


def distance_from_rate(vc, rate_density, **kwargs):
    r"""Returns the luminosity distance from the given total rate value

    Parameters
    ----------
    vc : float
        The total rate
    \**kwargs :
        All other keyword args are passed to :py:func:`get_cosmology` to
        select a cosmology. If none provided, will use
        :py:attr:`DEFAULT_COSMOLOGY`.

    Returns
    -------
    float :
        The luminosity distance at the given comoving volume.
    """

    cosmology = get_cosmology(**kwargs)
    if not hasattr(rate_density, 'dist_interp'):
        rate_density.dist_interp = {}

    if cosmology.name not in rate_density.dist_interp:
        def rate_func(z):
            return rate_from_redshift(z, rate_density) * rate_density(0.5).unit

        rate_density.dist_interp[cosmology.name] = ComovingVolInterpolator(
                                        'luminosity_distance',
                                        vol_func=rate_func,
                                        cosmology=cosmology)
    return rate_density.dist_interp[cosmology.name](vc)


# TODO: Add population models for BNS and NSBH.

# M_TOV = 2.22 # AP4 EoS
# A1 = 0.045
# A2 = 0.023

# M_TOV = 2.42 # DD2 EoS
# A1 = 0.046
# A2 = 0.014

# M_TOV = 2.77 # Ms1 EoS
# A1 = 0.042
# A2 = 0.010


__all__ = ['sfr_grb_2008', 'sfr_madau_dickinson_2014',
           'sfr_madau_fragos_2017', 'derivative_lookback_time',
           'p_tau', 'merger_rate_density', 'coalescence_rate_per_z_bin',
           'norm_redshift_distribution', 'rate_from_redshift', 
           'distance_from_rate',
          ]
