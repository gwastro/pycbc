# Module with utilities for estimating candidate events source probabilities
# Initial code by A. Curiel Barroso, August 2019
# Modified by V. Villa-Ortega, January 2020, March 2021

"""Functions to compute the area corresponding to different CBC on the m1 & m2
plane when given a central mchirp value and uncertainty.
It also includes a function that calculates the source frame when given the
detector frame mass and redshift.
"""

import math
import numpy as np
from pycbc.conversions import mass2_from_mchirp_mass1 as m2mcm1
from scipy.integrate import quad
from pycbc.cosmology import _redshift


def insert_args(parser):
    mchirp_group = parser.add_argument_group("Arguments for estimating the "
                                             "source probabilities of a "
                                             "candidate event using the snr, "
                                             "mchirp, and effective distance.")
    mchirp_group.add_argument('--src-class-mass-range', type=float, nargs=2,
                              metavar=('MIN_M2', 'MAX_M1'),
                              default=[1.0, 45.0],
                              help="Minimum and maximum values for the mass "
                                   "of the binary components, used as limits "
                                   "of the mass plane when computing the area "
                                   "corresponding to different CBC sources.")
    mchirp_group.add_argument('--src-class-mass-gap', type=float, nargs=2,
                              metavar=('MAX_NS', 'MIN_BH'), default=[3.0, 5.0],
                              help="Limits of the mass gap, that correspond "
                                   "to the maximum mass of a neutron star "
                                   "and the minimum mass of a black hole. "
                                   "Used as limits of integration of the "
                                   "different CBC regions.")
    mchirp_group.add_argument('--src-class-mchirp-to-delta', type=float,
                              metavar='m0', default=0.01,
                              help='Coefficient to estimate the value of the '
                                   'mchirp uncertainty by mchirp_delta = '
                                   'm0 * mchirp.')
    mchirp_group.add_argument('--src-class-eff-to-lum-distance', type=float,
                              metavar='a0', default=0.759,
                              help='Coefficient to estimate the value of the '
                                   'luminosity distance from the minimum '
                                   'eff distance by D_lum = a0 * min(D_eff).')
    mchirp_group.add_argument('--src-class-lum-distance-to-delta', type=float,
                              nargs=2, metavar=('b0', 'b1'),
                              default=[-0.449, -0.342],
                              help='Coefficients to estimate the value of the '
                                   'uncertainty on the luminosity distance '
                                   'from the estimated luminosity distance and'
                                   ' the coinc snr by delta_lum = D_lum * '
                                   'exp(b0) * coinc_snr ** b1.')
    mchirp_group.add_argument('--src-class-mass-gap-separate',
                              action='store_true',
                              help='Gives separate probabilities for each kind'
                                   ' of mass gap CBC sources: GNS, GG, BHG.')


def from_cli(args):
    return {'mass_limits': {'max_m1': args.src_class_mass_range[1],
                            'min_m2': args.src_class_mass_range[0]},
            'mass_bdary': {'ns_max': args.src_class_mass_gap[0],
                           'gap_max': args.src_class_mass_gap[1]},
            'estimation_coeff': {'a0': args.src_class_eff_to_lum_distance,
                                 'b0': args.src_class_lum_distance_to_delta[0],
                                 'b1': args.src_class_lum_distance_to_delta[1],
                                 'm0': args.src_class_mchirp_to_delta},
            'mass_gap': args.src_class_mass_gap_separate}


def src_mass_from_z_det_mass(z, del_z, mdet, del_mdet):
    """Takes values of redshift, redshift uncertainty, detector mass and its
    uncertainty and computes the source mass and its uncertainty.
    """
    msrc = mdet / (1. + z)
    del_msrc = msrc * ((del_mdet / mdet) ** 2.
                       + (del_z / (1. + z)) ** 2.) ** 0.5
    return (msrc, del_msrc)


def intmc(mc, x_min, x_max):
    """Returns the integral of m2 over m1 between x_min and x_max,
       assuming that mchirp is fixed.
    """
    integral = quad(lambda x, mc: m2mcm1(mc, x), x_min, x_max, args=mc)
    return integral[0]


def get_area(trig_mc, lim_h1, lim_h2, lim_v1, lim_v2):
    """Returns the area of the chirp mass contour in each region of the m1m2
       plane taking horizontal and vertical limits of the region as arguments.
    """
    mcb = trig_mc[0] + trig_mc[1]
    mcs = trig_mc[0] - trig_mc[1]
    # The points where the equal mass line and a chirp mass
    # curve intersect is m1 = m2 = 2**0.2 * mchirp
    mib = (2.**0.2) * mcb
    mis = (2.**0.2) * mcs

    if lim_h1 == 'diagonal':
        Pb_h1 = mib
        Ps_h1 = mis
        fun_sup = lambda x: x
    else:
        Pb_h1 = m2mcm1(mcb, lim_h1)
        Ps_h1 = m2mcm1(mcs, lim_h1)
        fun_sup = lambda x: lim_h1

    Pb_h2 = m2mcm1(mcb, lim_h2)
    Ps_h2 = m2mcm1(mcs, lim_h2)
    fun_inf = lambda x: lim_h2

    limb1 = np.clip(Pb_h1, lim_v1, lim_v2)
    limb2 = np.clip(Pb_h2, lim_v1, lim_v2)
    lims1 = np.clip(Ps_h1, lim_v1, lim_v2)
    lims2 = np.clip(Ps_h2, lim_v1, lim_v2)

    intb = intmc(mcb, limb1, limb2)
    ints = intmc(mcs, lims1, lims2)
    intline_sup = quad(fun_sup, lims1, limb1)[0]
    intline_inf = quad(fun_inf, lims2, limb2)[0]
    area = intb + intline_sup - ints - intline_inf
    return area


def calc_areas(trig_mc_det, mass_limits, mass_bdary, z, mass_gap):
    """Computes the area inside the lines of the second component mass as a
    function of the first component mass for the two extreme values
    of mchirp: mchirp +/- mchirp_uncertainty, for each region of the source
    classifying diagram.
    """
    trig_mc = src_mass_from_z_det_mass(z["central"], z["delta"],
                                       trig_mc_det["central"],
                                       trig_mc_det["delta"])
    m2_min = mass_limits["min_m2"]
    m1_max = mass_limits["max_m1"]
    ns_max = mass_bdary["ns_max"]
    gap_max = mass_bdary["gap_max"]

    abbh = get_area(trig_mc, 'diagonal', gap_max, gap_max, m1_max)
    abhg = get_area(trig_mc, gap_max, ns_max, gap_max, m1_max)
    ansbh = get_area(trig_mc, ns_max, m2_min, gap_max, m1_max)
    agg = get_area(trig_mc, 'diagonal', ns_max, ns_max, gap_max)
    agns = get_area(trig_mc, ns_max, m2_min, ns_max, gap_max)
    abns = get_area(trig_mc, 'diagonal', m2_min, m2_min, ns_max)

    if mass_gap:
        return {
            "BNS": abns,
            "GNS": agns,
            "NSBH": ansbh,
            "GG": agg,
            "BHG": abhg,
            "BBH": abbh
            }
    return {
        "BNS": abns,
        "NSBH": ansbh,
        "BBH": abbh,
        "Mass Gap": agns + agg + abhg
        }


def calc_probabilities(mchirp, snr, eff_distance, src_args):
    """Computes the different probabilities that a candidate event belongs to
       each CBC source category taking as arguments the chirp mass, the
       coincident SNR and the effective distance, and estimating the
       chirp mass uncertainty, the luminosity distance (and its uncertainty)
       and the redshift (and its uncertainty). Probability estimation is done
       assuming it is directly proportional to the area laying in the
       correspondent CBC region.
    """
    mass_limits = src_args['mass_limits']
    mass_bdary = src_args['mass_bdary']
    coeff = src_args['estimation_coeff']
    trig_mc_det = {'central': mchirp, 'delta': mchirp * coeff['m0']}
    dist_estimation = coeff['a0'] * eff_distance
    dist_std_estimation = (dist_estimation * math.exp(coeff['b0']) *
                           snr ** coeff['b1'])
    z_estimation = _redshift(dist_estimation)
    z_est_max = _redshift(dist_estimation + dist_std_estimation)
    z_est_min = _redshift(dist_estimation - dist_std_estimation)
    z_std_estimation = 0.5 * (z_est_max - z_est_min)
    z = {'central': z_estimation, 'delta': z_std_estimation}
    mass_gap = src_args['mass_gap']

    # If the mchirp is greater than the mchirp corresponding to two masses
    # equal to the maximum mass, the probability for BBH is 100%
    mc_max = mass_limits['max_m1'] / (2 ** 0.2)
    if trig_mc_det['central'] > mc_max * (1 + z['central']):
        if mass_gap:
            probabilities = {"BNS": 0.0, "GNS": 0.0, "NSBH": 0.0, "GG": 0.0,
                             "BHG": 0.0, "BBH": 1.0}
        else:
            probabilities = {"BNS": 0.0, "NSBH": 0.0, "BBH": 1.0,
                             "Mass Gap": 0.0}
    else:
        areas = calc_areas(trig_mc_det, mass_limits, mass_bdary, z, mass_gap)
        total_area = sum(areas.values())
        probabilities = {key: areas[key]/total_area for key in areas}
    return probabilities

