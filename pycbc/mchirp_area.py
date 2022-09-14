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
from astropy.cosmology import FlatLambdaCDM


def insert_args(parser):
    mchirp_group = parser.add_argument_group("Arguments for estimating the "
                                             "source probabilities of a "
                                             "candidate event using the snr, "
                                             "mchirp, and effective distance.")
    mchirp_group.add_argument('--src-class-mass-limits', type=float, nargs=3,
                              metavar=('MIN_M2', 'MAX_NS', 'MAX_M1'),
                              default=[1.0, 3.0, 45.0],
                              help="Minimum and maximum values for the mass "
                                   "of the binary components and maximum mass "
                                   "of a neutron star, used as limits "
                                   "when computing the area corresponding"
                                   "to different CBC sources.")
    mchirp_group.add_argument('--src-class-mass-gap-max', type=float,
                              metavar=('MAX_GAP'),
                              help="Upper limit of the mass gap, corresponding"
                                   " to the minimum mass of a black hole. "
                                   "Used as limit of integration of the "
                                   "different CBC regions when considering "
                                   "the MassGap category.")
    mchirp_group.add_argument('--src-class-mchirp-to-delta', type=float,
                              metavar='m0', required=True,
                              help='Coefficient to estimate the value of the '
                                   'mchirp uncertainty by mchirp_delta = '
                                   'm0 * mchirp.')
    mchirp_group.add_argument('--src-class-eff-to-lum-distance', type=float,
                              metavar='a0', required=True,
                              help='Coefficient to estimate the value of the '
                                   'luminosity distance from the minimum '
                                   'eff distance by D_lum = a0 * min(D_eff).')
    mchirp_group.add_argument('--src-class-lum-distance-to-delta', type=float,
                              nargs=2, metavar=('b0', 'b1'), required=True,
                              help='Coefficients to estimate the value of the '
                                   'uncertainty on the luminosity distance '
                                   'from the estimated luminosity distance and'
                                   ' the coinc snr by delta_lum = D_lum * '
                                   'exp(b0) * coinc_snr ** b1.')
    mchirp_group.add_argument('--src-class-mass-gap-separate',
                              action='store_true',
                              help='Gives separate probabilities for each kind'
                                   ' of mass gap CBC sources: GNS, GG, BHG.')
    mchirp_group.add_argument('--src-class-lal-cosmology',
                              action='store_true',
                              help='Uses the Planck15 cosmology defined in '
                                   'lalsuite instead of the astropy Planck15 '
                                   'default model.')


def from_cli(args, parser):
    mass_limits_sorted = sorted(args.src_class_mass_limits)
    if args.src_class_mass_gap_max:
        if args.src_class_mass_gap_max < mass_limits_sorted[1]:
            parser.error('MAX_GAP value cannot be lower than MAX_NS limit')
        return {'mass_limits':
                 {'max_m1': mass_limits_sorted[2],
                  'min_m2': mass_limits_sorted[0]},
                'mass_bdary':
                 {'ns_max': mass_limits_sorted[1],
                  'gap_max': args.src_class_mass_gap_max},
                'estimation_coeff':
                 {'a0': args.src_class_eff_to_lum_distance,
                  'b0': args.src_class_lum_distance_to_delta[0],
                  'b1': args.src_class_lum_distance_to_delta[1],
                  'm0': args.src_class_mchirp_to_delta},
                'mass_gap': True,
                'mass_gap_separate': args.src_class_mass_gap_separate,
                'lal_cosmology': args.src_class_lal_cosmology}
    return {'mass_limits':
             {'max_m1': mass_limits_sorted[2],
              'min_m2': mass_limits_sorted[0]},
            'mass_bdary':
             {'ns_max': mass_limits_sorted[1],
              'gap_max': mass_limits_sorted[1]},
            'estimation_coeff':
             {'a0': args.src_class_eff_to_lum_distance,
              'b0': args.src_class_lum_distance_to_delta[0],
              'b1': args.src_class_lum_distance_to_delta[1],
              'm0': args.src_class_mchirp_to_delta},
            'mass_gap': False,
            'mass_gap_separate': args.src_class_mass_gap_separate,
            'lal_cosmology': args.src_class_lal_cosmology}


def redshift_estimation(distance, distance_std, lal_cosmology):
    """Takes values of distance and its uncertainty and returns a
       dictionary with estimates of the redshift and its uncertainty.
       If the argument 'lal_cosmology' is True, it uses Planck15 cosmology
       model as defined in lalsuite instead of the astropy default.
       Constants for lal_cosmology taken from Planck15_lal_cosmology() in
       https://git.ligo.org/lscsoft/pesummary/-/blob/master/pesummary/gw/
       cosmology.py.
    """
    if lal_cosmology:
        cosmology = FlatLambdaCDM(H0=67.90, Om0=0.3065)
    else:
        cosmology = None
    z_estimation = _redshift(distance, cosmology=cosmology)
    z_est_max = _redshift((distance + distance_std),
                          cosmology=cosmology)
    z_est_min = _redshift((distance - distance_std),
                          cosmology=cosmology)
    z_std_estimation = 0.5 * (z_est_max - z_est_min)
    z = {'central': z_estimation, 'delta': z_std_estimation}
    return z


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
    """Returns the area under the chirp mass contour in a region of the m1-m2
    plane (m1 > m2).

    Parameters
    ----------
    trig_mc : sequence of two values
        first represents central estimate of mchirp in source frame,
        second its uncertainty
    lim_h1, lim_h2 : floats or the string 'diagonal'
        upper and lower horizontal limits of the region (limits on m2)
    lim_v1, lim_v2 : floats
        right and left vertical limits of the region (limits on m1)

    Returns
    -------
    area : float
    """
    mc_max = trig_mc[0] + trig_mc[1]
    mc_min = trig_mc[0] - trig_mc[1]
    # The points where the equal mass line and a chirp mass
    # curve intersect is m1 = m2 = 2**0.2 * mchirp
    mi_max = (2.**0.2) * mc_max
    mi_min = (2.**0.2) * mc_min

    if lim_h1 == 'diagonal':
        max_h1 = mi_max
        min_h1 = mi_min
        fun_sup = lambda x: x
    else:
        max_h1 = m2mcm1(mc_max, lim_h1)
        min_h1 = m2mcm1(mc_min, lim_h1)
        fun_sup = lambda x: lim_h1

    max_h2 = m2mcm1(mc_max, lim_h2)
    min_h2 = m2mcm1(mc_min, lim_h2)
    fun_inf = lambda x: lim_h2

    lim_max1 = np.clip(max_h1, lim_v1, lim_v2)
    lim_max2 = np.clip(max_h2, lim_v1, lim_v2)
    lim_min1 = np.clip(min_h1, lim_v1, lim_v2)
    lim_min2 = np.clip(min_h2, lim_v1, lim_v2)

    int_max = intmc(mc_max, lim_max1, lim_max2)
    int_min = intmc(mc_min, lim_min1, lim_min2)
    intline_sup = quad(fun_sup, lim_min1, lim_max1)[0]
    intline_inf = quad(fun_inf, lim_min2, lim_max2)[0]
    area = int_max + intline_sup - int_min - intline_inf
    return area


def calc_areas(
        trig_mc_det,
        mass_limits,
        mass_bdary,
        z,
        mass_gap,
        mass_gap_separate):
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
        if mass_gap_separate:
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
    return {
        "BNS": abns,
        "NSBH": ansbh,
        "BBH": abbh
        }


def calc_probabilities(mchirp, snr, eff_distance, src_args):
    """Computes the different probabilities that a candidate event belongs to
       each CBC source category taking as arguments the chirp mass, the
       coincident SNR and the effective distance, and estimating the
       chirp mass uncertainty, the luminosity distance (and its uncertainty)
       and the redshift (and its uncertainty). Probability is estimated to be
       directly proportional to the area of the corresponding CBC region.
    """
    mass_limits = src_args['mass_limits']
    mass_bdary = src_args['mass_bdary']
    coeff = src_args['estimation_coeff']
    trig_mc_det = {'central': mchirp, 'delta': mchirp * coeff['m0']}
    dist_estimation = coeff['a0'] * eff_distance
    dist_std_estimation = (dist_estimation * math.exp(coeff['b0']) *
                           snr ** coeff['b1'])
    z = redshift_estimation(dist_estimation, dist_std_estimation,
                            src_args['lal_cosmology'])
    mass_gap = src_args['mass_gap']
    mass_gap_separate = src_args['mass_gap_separate']

    # If the mchirp is greater than the mchirp corresponding to two masses
    # equal to the maximum mass, the probability for BBH is 100%.
    # If it is less than the mchirp corresponding to two masses equal to the
    # minimum mass, the probability for BNS is 100%.
    mc_max = mass_limits['max_m1'] / (2 ** 0.2)
    mc_min = mass_limits['min_m2'] / (2 ** 0.2)

    if trig_mc_det['central'] > mc_max * (1 + z['central']):
        if mass_gap:
            if mass_gap_separate:
                probabilities = {"BNS": 0.0, "GNS": 0.0, "NSBH": 0.0,
                                 "GG": 0.0, "BHG": 0.0, "BBH": 1.0}
            else:
                probabilities = {"BNS": 0.0, "NSBH": 0.0, "BBH": 1.0,
                                 "Mass Gap": 0.0}
        else:
            probabilities = {"BNS": 0.0, "NSBH": 0.0, "BBH": 1.0}

    elif trig_mc_det['central'] < mc_min * (1 + z['central']):
        if mass_gap:
            if mass_gap_separate:
                probabilities = {"BNS": 1.0, "GNS": 0.0, "NSBH": 0.0,
                                 "GG": 0.0, "BHG": 0.0, "BBH": 0.0}
            else:
                probabilities = {"BNS": 1.0, "NSBH": 0.0, "BBH": 0.0,
                                 "Mass Gap": 0.0}
        else:
            probabilities = {"BNS": 1.0, "NSBH": 0.0, "BBH": 0.0}

    else:
        areas = calc_areas(trig_mc_det, mass_limits, mass_bdary, z,
            mass_gap, mass_gap_separate)
        total_area = sum(areas.values())
        probabilities = {key: areas[key] / total_area for key in areas}

    return probabilities
