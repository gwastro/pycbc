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
This module provides functions for star formation rate models, time delay
models, merger rate density, and population models of BBH/BNS/NSBH.
"""

from functools import partial
import numpy as np
import scipy.integrate as scipy_integrate
import scipy.interpolate as scipy_interpolate
from astropy import units
from pycbc.cosmology import get_cosmology
from pycbc.cosmology import cosmological_quantity_from_redshift


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
    eta = -10

    rho_z = rho_local*((1+z)**(3.4*eta) + ((1+z)/5000)**(-0.3*eta) +
                       ((1+z)/9)**(-3.5*eta))**(1./eta)
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
        factor_a = 2.6
        factor_b = 3.2
        factor_c = 6.2
    elif mode == 'high':
        factor_a = 2.7
        factor_b = 3.0
        factor_c = 5.35
    else:
        raise ValueError("'mode' must choose from 'high' or 'low'.")
    rho_z = k_imf * 0.015 * (1+z)**factor_a / (1 + ((1+z)/factor_b)**factor_c)

    return rho_z


def diff_lookback_time(z, **kwargs):
    r""" The derivative of lookback time t(z)
         with respect to redshit z.

    Parameters
    ----------
    z : float
         The redshift.

    Returns
    -------
    dt_dz : float
            The value of dt/dz at the redshift z.
    \**kwargs :
        All other keyword args are passed to :py:func:`get_cosmology` to
        select a cosmology. If none provided, will use
        :py:attr:`DEFAULT_COSMOLOGY`.

    Notes
    -----
         Pease see Eq.(A3) in <arXiv:2011.02717v3> for more details.
    """
    from sympy import sqrt

    cosmology = get_cosmology(**kwargs)
    H0 = cosmology.H0.value * \
        (3.0856776E+19)**(-1)/(1/24/3600/365*1e-9)  # Gyr^-1
    dt_dz = 1/H0/(1+z)/sqrt((cosmology.Ode0+cosmology.Om0*(1+z)**3))
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
    from sympy import sqrt, exp, log

    if td_model == "log_normal":
        t_ln = 2.9  # Gyr
        sigma_ln = 0.2
        p_t = exp(-(log(tau)-log(t_ln))**2/(2*sigma_ln**2)) / \
                 (sqrt(2*np.pi)*sigma_ln)
    elif td_model == "gaussian":
        t_g = 2  # Gyr
        sigma_g = 0.3
        p_t = exp(-(tau-t_g)**2/(2*sigma_g**2)) / (sqrt(2*np.pi)*sigma_g)
    elif td_model == "power_law":
        alpha_t = 0.81
        p_t = tau**(-alpha_t)
    elif td_model == "inverse":
        p_t = tau**(-0.999)  # Try to avoid dividing zero.
    else:
        raise ValueError("'model' must choose from \
        ['log_normal', 'gaussian', 'power_law', 'inverse'].")

    return p_t


def convolution_trans(sfr, diff_lookback_t, model_td, **kwargs):
    r""" This function is used in a symbolic integral, which to calculate
        the merger rate density of CBC sources. This function converts the
        convolution of the star formation rate SFR(tau) and the time delay
        probability P(tau) on the time delay 'tau' into the convolution on
        the redshift 'z'.

    Parameters
    ----------
    sfr : function
            The star formation rate function used in the convolution.
    diff_lookback_t : function
            The derivative of lookback time t(z)
            with respect to redshit z.
    model_td : str
            The name of time delay model.
    \**kwargs :
        All other keyword args are passed to :py:func:`get_cosmology` to
        select a cosmology. If none provided, will use
        :py:attr:`DEFAULT_COSMOLOGY`.

    Returns
    -------
    func : sympy.core.symbol.Symbol
          The product of SFR(z), P(tau(z)) and dt(z)/dz.

    Notes
    -----
         Pease see Eq.(A2) in <arXiv:2011.02717v3> for more details.
    """
    from sympy import integrate, symbols

    if model_td not in ['log_normal', 'gaussian', 'power_law', 'inverse']:
        raise ValueError("'model_td' must choose from \
        ['log_normal', 'gaussian', 'power_law', 'inverse'].")

    # Fix the cosmology, set 'z/z_0' to be the only
    # parameter in the symbolic integration.
    diff_lookback_time_z = partial(diff_lookback_t, **kwargs)
    z = symbols('z')
    z_0 = symbols('z_0')
    tau = integrate(diff_lookback_time_z(z), (z, z_0, z))
    func = sfr(z) * p_tau(tau, model_td) * diff_lookback_time_z(z)
    return func


def merger_rate_density(sfr_func, td_model, rho_local, maxz=10.0,
                        npoints=10000, z_array=None, **kwargs):
    r""" This function uses the symbolic integral to calculate
        the merger rate density of CBC sources. This function converts the
        convolution of the star formation rate SFR(tau) and the time delay
        probability P(tau) on the time delay 'tau' into the convolution on
        the redshift 'z'. This function relies on `convolution_trans`.

    Parameters
    ----------
    sfr_func : function
            The star formation rate function used in the convolution.
    td_model : str
            The name of time delay model.
    rho_local : float
            The local merger rate of a certain type of CBC source.
            In the unit of "Mpc^-3yr^-1".
    maxz : float
            The max redshift. The default value is 10.
    npoints : int
            The number of points used in the interpolation. The default
            value is 10000.
    z_array : numpy.array
            The array of redshift. The default value is None.
    \**kwargs :
        All other keyword args are passed to :py:func:`get_cosmology` to
        select a cosmology. If none provided, will use
        :py:attr:`DEFAULT_COSMOLOGY`.

    Returns
    -------
    rho_z : scipy.interpolate.interp1d
          The merger rate density.

    Notes
    -----
         Pease see Eq.(A1), Eq.(A2) in <arXiv:2011.02717v3> for more details.
    """
    from sympy import symbols, lambdify

    if z_array is None:
        z_array = np.linspace(0, maxz, npoints)

    if td_model not in ['log_normal', 'gaussian', 'power_law', 'inverse']:
        raise ValueError("'td_model' must choose from \
        ['log_normal', 'gaussian', 'power_law', 'inverse'].")

    z = symbols('z')
    z_0 = symbols('z_0')
    f_z = np.zeros(len(z_array))

    func_1 = convolution_trans(
                sfr=sfr_func, diff_lookback_t=diff_lookback_time,
                model_td=td_model, **kwargs)
    for i in range(len(z_array)):
        func_2 = lambdify(z, func_1.subs(z_0, z_array[i]), 'scipy')
        f_z[i] = scipy_integrate.quad(
                    func_2, z_array[i], np.inf, epsabs=1.49e-3)[0]

    f_z = f_z/f_z[0]*rho_local  # Normalize & Rescale
    rho_z = scipy_interpolate.interp1d(z_array, f_z)
    return rho_z


def coalescence_rate(rate_den, maxz=10.0, npoints=10000,
                     z_array=None, **kwargs):
    r""" This function calculates the coalescence(merger) rate at the redshift z.

    Parameters
    ----------
    rate_den : function or scipy.interpolate.interp1d
        The merger rate density as a function of redshift. In the unit of
        "Mpc^-3yr^-1". Use `merger_rate_density` function to calculate.
    maxz : float
            The max redshift. The default value is 10.
    npoints : int
            The number of points used in the interpolation. The default
            value is 10000.
    z_array : numpy.array or list
            The redshift range. The default value is None.
    \**kwargs :
        All other keyword args are passed to :py:func:`get_cosmology` to
        select a cosmology. If none provided, will use
        :py:attr:`DEFAULT_COSMOLOGY`.

    Returns
    -------
    coalescence_rate_interp : scipy.interpolate.interp1d
          The coalescence rate.

    Notes
    -----
         Pease see Eq.(1) in <arXiv:2011.02717v3> for more details.
    """

    if z_array is None:
        z_array = np.linspace(0, maxz, npoints)

    dr_dz = []
    cosmology = get_cosmology(**kwargs)

    for z in z_array:
        dr = cosmology.differential_comoving_volume(z) / (1+z)
        dr_dz.append((dr*4*np.pi*units.sr*rate_den(z)*(units.Mpc)**(-3)).value)

    coalescence_rate_interp = scipy_interpolate.interp1d(
                z_array, dr_dz, fill_value='extrapolate')

    return coalescence_rate_interp


def total_rate_upto_redshift(z, merger_rate):
    r"""Total rate of occurrences out to some redshift.

    Parameters
    ----------
    z : int, float, tuple, numpy.ndarray or list
            The redshift.
    merger_rate : scipy.interpolate.interp1d or function
        The merger or coalescence rate. Function should take the
        redshift as a single argument. Provided by users or
        calculated by the `coalescence_rate` function.

    Returns
    -------
    rate: float or list
        The total rate of occurrences out to some redshift. In the
        unit of "yr^-1".
    """

    if isinstance(z, (float, int)):
        total_rate = scipy_integrate.quad(
                        merger_rate, 0, z,
                        epsabs=2.00e-4, epsrel=2.00e-4, limit=1000)[0]
    elif isinstance(z, (tuple, np.ndarray, list)):
        total_rate = []
        for redshift in z:
            total_rate.append(
                scipy_integrate.quad(
                    merger_rate, 0, redshift,
                    epsabs=2.00e-4, epsrel=2.00e-4, limit=1000)[0]
            )
    else:
        raise ValueError("'z' must be 'int', 'float', 'tuple', \
                            'numpy.ndarray' or 'list'.")

    return total_rate


def average_time_between_signals(z_array, merger_rate):
    r""" This function calculates the average time interval
        of a certain type of CBC source.

    Parameters
    ----------
    z_array : numpy.array
            The array of redshift.
    merger_rate : scipy.interpolate.interp1d or function
            The coalescence rate. Provided by users or calculated by
            the `coalescence_rate` function.

    Returns
    -------
    average_time : float
            The average time interval (s).
    """

    total_rate = total_rate_upto_redshift(
            z_array[-1], merger_rate)  # yr^-1
    average_time = 1./total_rate * 365*24*3600  # s
    return average_time


def norm_redshift_distribution(z_array, merger_rate):
    r""" This function calculates the normalized redshift distribution
        of a certain type of CBC source.

    Parameters
    ----------
    z_array : numpy.array
            The array of redshift.
    merger_rate : scipy.interpolate.interp1d or function
            The coalescence rate. Provided by users or calculated by
            the `coalescence_rate` function.

    Returns
    -------
    norm_coalescence_rate : numpy.array
            The normalized redshift distribution.

    Notes
    -----
         The can be used as a population-informed prior for redshift
         and luminosity distance of CBC sources.
    """

    lambd = average_time_between_signals(z_array, merger_rate)
    norm_coalescence_rate = lambd/(365*24*3600) * merger_rate(z_array)
    return norm_coalescence_rate


def distance_from_rate(
        total_rate, merger_rate, maxz=10, npoints=10000, **kwargs):
    r"""Returns the luminosity distance from the given total rate value.

    Parameters
    ----------
    total_rate : float
        The total rate.
    merger_rate : scipy.interpolate.interp1d or function
            The coalescence rate. Provided by users or calculated by
            the `coalescence_rate` function.
    maxz : float
        The max redshift in the interpolation, the default value is 10.
        Please use the same `maxz` as in `merger_rate`.
    npoints : int
        The number of points used in the interpolation, the default value
        is 10000.
    \**kwargs :
        All other keyword args are passed to :py:func:`get_cosmology` to
        select a cosmology. If none provided, will use
        :py:attr:`DEFAULT_COSMOLOGY`.

    Returns
    -------
    dl : float
        The luminosity distance at the given total rate value.
        In the unit of "Mpc".

    Notes
    -----
         This can be used in a population-informed prior for redshift
         and luminosity distance of CBC sources. When this used in
         high redshift range, please first use the `total_rate_upto_redshift`
         function to plot the curve and find the point where the curve
         starts to stay almost horizontal, then set `maxz` to the
         corresponding value and change `npoints` to a reasonable value.
    """
    cosmology = get_cosmology(**kwargs)

    if not hasattr(merger_rate, 'dist_interp'):
        merger_rate.dist_interp = {}

    if ((cosmology.name not in merger_rate.dist_interp) or
       (len(merger_rate.dist_interp[cosmology.name].x) != npoints)):
        def rate_func(redshift):
            return total_rate_upto_redshift(redshift, merger_rate)

        z_array = np.linspace(0, maxz, npoints)
        dists = cosmological_quantity_from_redshift(
                z_array, 'luminosity_distance', **kwargs)
        total_rates = rate_func(z_array)
        interp = scipy_interpolate.interp1d(total_rates, dists)
        merger_rate.dist_interp[cosmology.name] = interp

    dl = merger_rate.dist_interp[cosmology.name](total_rate)
    if np.isscalar(dl):
        dl = float(dl)
    return dl


__all__ = ['sfr_grb_2008', 'sfr_madau_dickinson_2014',
           'sfr_madau_fragos_2017', 'diff_lookback_time',
           'p_tau', 'merger_rate_density', 'coalescence_rate',
           'norm_redshift_distribution', 'total_rate_upto_redshift',
           'distance_from_rate', 'average_time_between_signals']
