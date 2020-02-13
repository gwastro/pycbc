# Module with utilities for estimating candidate events source probabilities
# Initial code by A. Curiel Barroso, August 2019
# Modified by V. Villa-Ortega, January 2020

"""Functions to compute the area corresponding to different CBC on the m1 & m2
plane when given a central mchirp value and uncertainty.
It also includes a function that calculates the source frame when given the
detector frame mass and redshift.
"""

import math
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
    """Returns the integral of a component mass as a function of the mass of
       the other component, taking mchirp as an argument.
    """
    integral = quad(lambda x, mc: m2mcm1(mc, x), x_min, x_max, args=mc)
    return integral[0]


def calc_areas(trig_mc_det, mass_limits, mass_bdary, z, mass_gap):
    """Computes the area inside the lines of the second component mass as a
    function of the first component mass for the two extreme values
    of mchirp: mchirp +/- mchirp_uncertainty, for each region of the source
    classifying diagram.
    """
    trig_mc = src_mass_from_z_det_mass(z["central"], z["delta"],
                                       trig_mc_det["central"],
                                       trig_mc_det["delta"])
    mcb = trig_mc[0] + trig_mc[1]
    mcs = trig_mc[0] - trig_mc[1]
    m2_min = mass_limits["min_m2"]
    m1_max = mass_limits["max_m1"]
    ns_max = mass_bdary["ns_max"]
    gap_max = mass_bdary["gap_max"]
    # The points where the equal mass line and a chirp mass
    # curve intersect is m1 = m2 = 2**0.2 * mchirp
    mib = (2.**0.2) * mcb
    mis = (2.**0.2) * mcs

    # AREA FOR BBH
    if mib < gap_max:
        abbh = 0.0
    else:
        limb_bbh = min(m1_max, m2mcm1(mcb, gap_max))
        intb_bbh = intmc(mcb, mib, limb_bbh)

        if mis < gap_max:
            lims1_bbh = gap_max
            lims2_bbh = lims1_bbh
        else:
            lims1_bbh = mis
            lims2_bbh = min(m1_max, m2mcm1(mcs, gap_max))

        ints_bbh = intmc(mcs, lims1_bbh, lims2_bbh)

        limdiag_bbh = max(m2mcm1(mcs, lims1_bbh), gap_max)
        intline_sup_bbh = 0.5 * (limdiag_bbh + mib) * (mib - lims1_bbh)
        intline_inf_bbh = (limb_bbh - lims2_bbh) * gap_max
        int_sup_bbh = intb_bbh + intline_sup_bbh
        int_inf_bbh = ints_bbh + intline_inf_bbh

        abbh = int_sup_bbh - int_inf_bbh

    # AREA FOR BHG
    if m2mcm1(mcb, gap_max) < ns_max or m2mcm1(mcs, m1_max) > gap_max:
        abhg = 0.0
    else:
        if m2mcm1(mcb, m1_max) > gap_max:
            limb2_bhg = m1_max
            limb1_bhg = limb2_bhg
        else:
            limb2_bhg = min(m1_max, m2mcm1(mcb, ns_max))
            limb1_bhg = max(gap_max, m2mcm1(mcb, gap_max))

        intb_bhg = intmc(mcb, limb1_bhg, limb2_bhg)

        if m2mcm1(mcs, gap_max) < ns_max:
            lims2_bhg = gap_max
            lims1_bhg = lims2_bhg
        else:
            lims1_bhg = max(gap_max, m2mcm1(mcs, gap_max))
            lims2_bhg = min(m1_max, m2mcm1(mcs, ns_max))

        intline_inf_bhg = (limb2_bhg - lims2_bhg) * ns_max
        intline_sup_bhg = (limb1_bhg - lims1_bhg) * gap_max
        ints_bhg = intmc(mcs, lims1_bhg, lims2_bhg)
        int_sup_bhg = intb_bhg + intline_sup_bhg
        int_inf_bhg = ints_bhg + intline_inf_bhg

        abhg = int_sup_bhg - int_inf_bhg

    # AREA FOR GG
    if m2mcm1(mcs, gap_max) > gap_max or m2mcm1(mcb, ns_max) < ns_max:
        agg = 0.0
    else:
        if m2mcm1(mcb, gap_max) > gap_max:
            limb2_gg = gap_max
            limb1_gg = limb2_gg
        else:
            limb1_gg = mib
            limb2_gg = min(gap_max, m2mcm1(mcb, ns_max))

        intb_gg = intmc(mcb, limb1_gg, limb2_gg)

        if m2mcm1(mcs, ns_max) < ns_max:
            lims2_gg = ns_max
            lims1_gg = lims2_gg
        else:
            lims1_gg = mis
            lims2_gg = min(gap_max, m2mcm1(mcs, ns_max))

        ints_gg = intmc(mcs, lims1_gg, lims2_gg)
        limdiag1_gg = max(m2mcm1(mcs, lims1_gg), ns_max)
        limdiag2_gg = min(m2mcm1(mcb, limb1_gg), gap_max)
        intline_sup_gg = (0.5 * (limb1_gg - lims1_gg)
                          * (limdiag1_gg + limdiag2_gg))
        intline_inf_gg = (limb2_gg - lims2_gg) * ns_max
        int_sup_gg = intb_gg + intline_sup_gg
        int_inf_gg = ints_gg + intline_inf_gg

        agg = int_sup_gg - int_inf_gg

    # AREA FOR BNS
    if m2mcm1(mcs, ns_max) > ns_max:
        abns = 0.0
    else:
        if m2mcm1(mcb, ns_max) > ns_max:
            limb2_bns = ns_max
            limb1_bns = limb2_bns
        else:
            limb2_bns = min(ns_max, m2mcm1(mcb, m2_min))
            limb1_bns = mib

        intb_bns = intmc(mcb, limb1_bns, limb2_bns)

        if mis < m2_min:
            lims2_bns = m2_min
            lims1_bns = lims2_bns
        else:
            lims2_bns = min(ns_max, m2mcm1(mcs, m2_min))
            lims1_bns = mis

        ints_bns = intmc(mcs, lims1_bns, lims2_bns)
        intline_inf_bns = (limb2_bns - lims2_bns) * m2_min
        limdiag1_bns = max(m2mcm1(mcs, lims1_bns), m2_min)
        limdiag2_bns = min(m2mcm1(mcb, limb1_bns), ns_max)
        intline_sup_bns = (0.5 * (limdiag1_bns + limdiag2_bns)
                           * (limb1_bns - lims1_bns))
        int_sup_bns = intb_bns + intline_sup_bns
        int_inf_bns = ints_bns + intline_inf_bns

        abns = int_sup_bns - int_inf_bns

    # AREA FOR GNS
    if m2mcm1(mcs, gap_max) > ns_max or m2mcm1(mcb, ns_max) < m2_min:
        agns = 0.0
    else:
        if m2mcm1(mcb, gap_max) > ns_max:
            limb2_gns = gap_max
            limb1_gns = limb2_gns
        else:
            limb2_gns = min(gap_max, m2mcm1(mcb, m2_min))
            limb1_gns = max(ns_max, m2mcm1(mcb, ns_max))

        intb_gns = intmc(mcb, limb1_gns, limb2_gns)

        if m2mcm1(mcs, ns_max) < m2_min:
            lims2_gns = ns_max
            lims1_gns = lims2_gns
        else:
            lims1_gns = max(ns_max, m2mcm1(mcs, ns_max))
            lims2_gns = min(gap_max, m2mcm1(mcs, m2_min))

        intline_inf_gns = (limb2_gns - lims2_gns) * m2_min
        intline_sup_gns = (limb1_gns - lims1_gns) * ns_max
        ints_gns = intmc(mcs, lims1_gns, lims2_gns)
        int_sup_gns = intb_gns + intline_sup_gns
        int_inf_gns = ints_gns + intline_inf_gns

        agns = int_sup_gns - int_inf_gns

    # AREA FOR NSBH
    if m2mcm1(mcs, m1_max) > ns_max or m2mcm1(mcb, gap_max) < m2_min:
        ansbh = 0.0
    else:
        if m2mcm1(mcb, m1_max) > ns_max:
            limb2_nsbh = m1_max
            limb1_nsbh = limb2_nsbh
        else:
            limb1_nsbh = max(gap_max, m2mcm1(mcb, ns_max))
            limb2_nsbh = min(m1_max, m2mcm1(mcb, m2_min))

        intb_nsbh = intmc(mcb, limb1_nsbh, limb2_nsbh)

        if m2mcm1(mcs, gap_max) < m2_min:
            lims1_nsbh = gap_max
            lims2_nsbh = lims1_nsbh
        else:
            lims1_nsbh = max(gap_max, m2mcm1(mcs, ns_max))
            lims2_nsbh = min(m1_max, m2mcm1(mcs, m2_min))

        intline_inf_nsbh = (limb2_nsbh - lims2_nsbh) * m2_min
        intline_sup_nsbh = (limb1_nsbh - lims1_nsbh) * ns_max
        ints_nsbh = intmc(mcs, lims1_nsbh, lims2_nsbh)
        int_sup_nsbh = intb_nsbh + intline_sup_nsbh
        int_inf_nsbh = ints_nsbh + intline_inf_nsbh

        ansbh = int_sup_nsbh - int_inf_nsbh
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
