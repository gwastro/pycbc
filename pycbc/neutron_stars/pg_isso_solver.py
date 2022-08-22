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
Innermost Stable Spherical Orbit (ISSO) solver in the Perez-Giz (PG)
formalism. See `Stone, Loeb, Berger, PRD 87, 084053 (2013)`_.

.. _Stone, Loeb, Berger, PRD 87, 084053 (2013):
    http://dx.doi.org/10.1103/PhysRevD.87.084053
"""
import pickle
import os.path
import numpy as np
from scipy.optimize import root_scalar
from scipy.interpolate import RectBivariateSpline
from . import NS_DATA_DIRECTORY


def ISCO_eq(r, chi):
    """Polynomial that enables the calculation of the Kerr innermost
    stable circular orbit (ISCO) radius via its roots,

    .. math:: Z(r) = [r (r-6)]^2 - \chi^2 [2r (3r+14) - 9 \chi^2]\,.

    Parameters
    -----------
    r: float
        the radial coordinate in BH mass units
    chi: float
        the BH dimensionless spin parameter

    Returns
    ----------
    float
    """
    return (r * (r - 6))**2 - chi**2 * (2 * r * (3 * r + 14) - 9 * chi**2)


def ISCO_eq_dr(r, chi):
    """Partial derivative of :func:`ISCO_eq` with respect to r.

    Parameters
    ----------
    r: float
        the radial coordinate in BH mass units
    chi: float
        the BH dimensionless spin parameter

    Returns
    -------
    float
    """
    return 4 * r**3 - 36 * r**2 + 12 * (6 - chi**2) * r - 28 * chi**2


def ISCO_eq_dr2(r, chi):
    """Double partial derivative of :func:`ISCO_eq` with respect to r.

    Parameters
    ----------
    r: float
        the radial coordinate in BH mass units
    chi: float
        the BH dimensionless spin parameter

    Returns
    -------
    float
    """
    return 12 * r**2 - 72 * r + 12 * (6 - chi**2)


def ISSO_eq_at_pole(r, chi):
    r"""Polynomial that enables the calculation of the Kerr polar
    (:math:`\iota = \pm \pi / 2`) innermost stable spherical orbit
    (ISSO) radius via the roots of

    .. math::

        P(r) &= r^3 [r^2 (r - 6) + \chi^2 (3 r + 4)] \\
             &\quad + \chi^4 [3 r (r - 2) + \chi^2] \, ,

    where :math:`\chi` is the BH dimensionless spin parameter. Physical
    solutions are between 6 and
    :math:`1 + \sqrt{3} + \sqrt{3 + 2 \sqrt{3}}`.

    Parameters
    ----------
    r: float
        the radial coordinate in BH mass units
    chi: float
        the BH dimensionless spin parameter

    Returns
    -------
    float
    """
    chi2 = chi * chi
    return (
        r**3 * (r**2 * (r - 6) + chi2 * (3 * r + 4))
        + chi2 * chi2 * (3 * r * (r - 2) + chi2))


def ISSO_eq_at_pole_dr(r, chi):
    """Partial derivative of :func:`ISSO_eq_at_pole` with respect to r.

    Parameters
    ----------
    r: float
        the radial coordinate in BH mass units
    chi: float
        the BH dimensionless spin parameter

    Returns
    -------
    float
    """
    chi2 = chi * chi
    twlvchi2 = 12 * chi2
    sxchi4 = 6 * chi2 * chi2
    return (
        6 * r**5 - 30 * r**4 + twlvchi2 * r**3 + twlvchi2 * r**2 + sxchi4 * r
        + sxchi4)


def ISSO_eq_at_pole_dr2(r, chi):
    """Double partial derivative of :func:`ISSO_eq_at_pole` with
    respect to r.

    Parameters
    ----------
    r: float
        the radial coordinate in BH mass units
    chi: float
        the BH dimensionless spin parameter

    Returns
    -------
    float
    """
    chi2 = chi * chi
    return (
        30 * r**4 - 120 * r**3 + 36 * chi2 * r**2 + 24 * chi2 * r
        + 6 * chi2 * chi2)


def PG_ISSO_eq(r, chi, incl):
    r"""Polynomial that enables the calculation of a generic innermost
    stable spherical orbit (ISSO) radius via the roots in :math:`r` of

    .. math::

        S(r) &= r^8 Z(r) + \chi^2 (1 - \cos(\iota)^2) \\
             &\quad * [\chi^2 (1 - \cos(\iota)^2) Y(r) - 2 r^4 X(r)]\,,

    where

    .. math::

        X(r) &= \chi^2 (\chi^2 (3 \chi^2 + 4 r (2 r - 3)) \\
             &\quad + r^2 (15 r (r - 4) + 28)) - 6 r^4 (r^2 - 4) \, ,

    .. math::

        Y(r) &= \chi^4 (\chi^4 + r^2 [7 r (3 r - 4) + 36]) \\
             &\quad + 6 r (r - 2) \\
             &\qquad * (\chi^6 + 2 r^3
             [\chi^2 (3 r + 2) + 3 r^2 (r - 2)]) \, ,

    and :math:`Z(r) =` :func:`ISCO_eq`. Physical solutions are between
    the equatorial ISSO (i.e. the ISCO) radius (:func:`ISCO_eq`) and
    the polar ISSO radius (:func:`ISSO_eq_at_pole`).
    See `Stone, Loeb, Berger, PRD 87, 084053 (2013)`_.

    .. _Stone, Loeb, Berger, PRD 87, 084053 (2013):
        http://dx.doi.org/10.1103/PhysRevD.87.084053

    Parameters
    ----------
    r: float
        the radial coordinate in BH mass units
    chi: float
        the BH dimensionless spin parameter
    incl: float
        inclination angle between the BH spin and the orbital angular
        momentum in radians

    Returns
    -------
    float
    """
    chi2 = chi * chi
    chi4 = chi2 * chi2
    r2 = r * r
    r4 = r2 * r2
    three_r = 3 * r
    r_minus_2 = r - 2
    sin_incl2 = (np.sin(incl))**2

    X = (
        chi2 * (
            chi2 * (3 * chi2 + 4 * r * (2 * r - 3))
            + r2 * (15 * r * (r - 4) + 28))
        - 6 * r4 * (r2 - 4))
    Y = (
        chi4 * (chi4 + r2 * (7 * r * (three_r - 4) + 36))
        + 6 * r * r_minus_2 * (
            chi4 * chi2 + 2 * r2 * r * (
                chi2 * (three_r + 2) + 3 * r2 * r_minus_2)))
    Z = ISCO_eq(r, chi)

    return r4 * r4 * Z + chi2 * sin_incl2 * (chi2 * sin_incl2 * Y - 2 * r4 * X)


def PG_ISSO_eq_dr(r, chi, incl):
    """Partial derivative of :func:`PG_ISSO_eq` with respect to r.

    Parameters
    ----------
    r: float
        the radial coordinate in BH mass units
    chi: float
        the BH dimensionless spin parameter
    incl: float
        inclination angle between the BH spin and the orbital angular
        momentum in radians

    Returns
    -------
    float
    """
    sini = np.sin(incl)
    sin2i = sini * sini
    sin4i = sin2i * sin2i
    chi2 = chi * chi
    chi4 = chi2 * chi2
    chi6 = chi4 * chi2
    chi8 = chi4 * chi4
    chi10 = chi6 * chi4
    return (
        12 * r**11 - 132 * r**10
        + r**9 * (120 * chi2 * sin2i - 60 * chi2 + 360) - r**8 * 252 * chi2
        + 8 * r**7 * (
            36 * chi4 * sin4i - 30 * chi4 * sin2i + 9 * chi4
            - 48 * chi2 * sin2i)
        + 7 * r**6 * (120 * chi4 * sin2i - 144 * chi4 * sin4i)
        + 6 * r**5 * (
            36 * chi6 * sin4i - 16 * chi6 * sin2i + 144 * chi4 * sin4i
            - 56 * chi4 * sin2i)
        + r**4 * (120 * chi6 * sin2i - 240 * chi6 * sin4i)
        + r**3 * (84 * chi8 * sin4i - 24 * chi8 * sin2i - 192 * chi6 * sin4i)
        - 84 * r**2 * chi8 * sin4i
        + r * (12 * chi10 * sin4i + 72 * chi8 * sin4i) - 12 * chi10 * sin4i)


def PG_ISSO_eq_dr2(r, chi, incl):
    """Second partial derivative of :func:`PG_ISSO_eq` with respect to
    r.

    Parameters
    ----------
    r: float
        the radial coordinate in BH mass units
    chi: float
        the BH dimensionless spin parameter
    incl: float
        inclination angle between the BH spin and the orbital angular
        momentum in radians

    Returns
    -------
    float
    """
    sini = np.sin(incl)
    sin2i = sini * sini
    sin4i = sin2i * sin2i
    chi2 = chi * chi
    chi4 = chi2 * chi2
    chi6 = chi4 * chi2
    chi8 = chi4 * chi4
    return (
        132 * r**10 - 1320 * r**9
        + 90 * r**8 * (12 * chi2 * sin2i - 6 * chi2 + 36) - 2016 * chi2 * r**7
        + 56 * r**6 * (
            36 * chi4 * sin4i - 30 * chi4 * sin2i + 9 * chi4
            - 48 * chi2 * sin2i)
        + 42 * r**5 * (120 * chi4 * sin2i - 144 * chi4 * sin4i)
        + 30 * r**4 * (
            36 * chi6 * sin4i - 16 * chi6 * sin2i + 144 * chi4 * sin4i
            - 56 * chi4 * sin2i)
        + r**3 * (480 * chi6 * sin2i - 960 * chi6 * sin4i)
        + r**2 * (
            252 * chi8 * sin4i - 72 * chi8 * sin2i - 576 * chi6 * sin4i)
        - r * 168 * chi8 * sin4i
        + 12 * chi8 * chi2 * sin4i + 72 * chi8 * sin4i)


def PG_ISSO_solver(chi, incl):
    """Function that determines the radius of the innermost stable
    spherical orbit (ISSO) for a Kerr BH and a generic inclination
    angle between the BH spin and the orbital angular momentum.
    This function finds the appropriate root of :func:`PG_ISSO_eq`.

    Parameters
    ----------
    chi: float
        the BH dimensionless spin parameter
    incl: float
        the inclination angle between the BH spin and the orbital
        angular momentum in radians

    Returns
    -------
    solution: float
        the radius of the orbit in BH mass units
    """
    # Auxiliary variables
    cos_incl = np.cos(incl)
    sgnchi = np.sign(cos_incl)*chi
    if np.isscalar(sgnchi):
        sgnchi = np.array(sgnchi, copy=False, ndmin=1)
        chi = np.array(chi, copy=False, ndmin=1)
        incl = np.array(incl, copy=False, ndmin=1)

    # ISCO radius for the given spin magnitude
    initial_guess = [2 if s > 0.99 else (9 if s < 0 else 5) for s in sgnchi]
    rISCO_limit = np.array([
        root_scalar(
            ISCO_eq, x0=g0, fprime=ISCO_eq_dr, fprime2=ISCO_eq_dr2,
            args=(sc)).root
        for g0, sc in zip(initial_guess, sgnchi)])
    # If the inclination is 0 or pi, just output the ISCO radius
    equatorial = np.isclose(incl, 0) | np.isclose(incl, np.pi)
    if all(equatorial):
        return rISCO_limit

    # ISSO radius for an inclination of pi/2
    initial_guess = [9 if c < 0 else 6 for c in chi]
    rISSO_at_pole_limit = np.array([
        root_scalar(
            ISSO_eq_at_pole, x0=g0, fprime=ISSO_eq_at_pole_dr,
            fprime2=ISSO_eq_at_pole_dr2, args=(c)).root
        for g0, c in zip(initial_guess, chi)])
    # If the inclination is pi/2, just output the ISSO radius at the pole(s)
    polar = np.isclose(incl, np.pi / 2)
    if all(polar):
        return rISSO_at_pole_limit

    # Otherwise, find the ISSO radius for a generic inclination
    initial_hi = np.maximum(rISCO_limit, rISSO_at_pole_limit)
    initial_lo = np.minimum(rISCO_limit, rISSO_at_pole_limit)
    brackets = [
        (bl, bh) if c == 1 else None
        for bl, bh, c in zip(initial_lo, initial_hi, chi)]
    solution = np.array([
        root_scalar(
            PG_ISSO_eq, x0=g0, fprime=PG_ISSO_eq_dr, bracket=bracket,
            fprime2=PG_ISSO_eq_dr2, args=(c, inc), xtol=1e-12).root
        for g0, bracket, c, inc in zip(initial_hi, brackets, chi, incl)])
    oob = (solution < 1) | (solution > 9)
    if any(oob):
        solution = np.array([
            root_scalar(
                PG_ISSO_eq, x0=g0, fprime=PG_ISSO_eq_dr, bracket=bracket,
                fprime2=PG_ISSO_eq_dr2, args=(c, inc)).root
            if ob else sol for g0, bracket, c, inc, ob, sol
            in zip(initial_lo, brackets, chi, incl, oob, solution)
            ])
        oob = (solution < 1) | (solution > 9)
        if any(oob):
            raise RuntimeError('Unable to obtain some solutions!')
    return solution


def concat_grid(bounds):
    '''Constructs non-uniform grid given specified bounds'''
    out = np.concatenate([
        np.linspace(0, lim[0], lim[1], endpoint=False) if ii == 0
        else (
            np.linspace(bounds[ii-1][0], lim[0], lim[1])
            if ii == len(bounds)
            else np.linspace(bounds[ii-1][0], lim[0], lim[1], endpoint=False)
            )
        for ii, lim in enumerate(bounds)])
    return out


def generate_isso_bivariate_interp():
    """Constructs a grid in spin magnitude and spin tilt angle then
    solves for the ISSO radius (in mass units). Uses the resulting
    grid to create a bivariate spine which is saved to disk.
    Note: the grid density need not be uniform, and so has been tuned
    by hand to be more dense in regions of steeper gradient in order to
    keep the fractional interpolation error over the full span of
    parameter values to less than 1e-5.
    """
    chis = concat_grid((
        (0.5, 25), (0.75, 20), (0.95, 30), (0.99, 20), (0.995, 50),
        (0.9975, 50), (1, 101)))
    incs = concat_grid((
        (0.125 * np.pi, 50), (0.375 * np.pi, 300), (0.5 * np.pi, 50),
        (0.9 * np.pi, 50), (np.pi, 26)))
    chig, incg = np.meshgrid(chis, incs)
    roots = np.empty_like(chig)
    for i, (chi, inc) in enumerate(zip(chig, incg)):
        # more dense at larger chi
        roots[i] = PG_ISSO_solver(chi, inc)
    bivar = RectBivariateSpline(incs, chis, roots)
    with open(os.path.join(NS_DATA_DIRECTORY, 'isso_inc_chi.pkl'), 'wb') as f:
        pickle.dump(bivar, f)

        
def load_isso_bivariate_interp():
    with open(os.path.join(NS_DATA_DIRECTORY, 'isso_inc_chi.pkl'), 'rb') as f:
        func = pickle.load(f)
    return func


def pg_isso_interp(incl, chi):
    return load_isso_bivariate_interp()(incl, chi, grid=False)
