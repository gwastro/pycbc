# Integration of the area laying in the different cbc regions
# By A. Curiel Barroso
# August 2019

"""Functions to compute the area corresponding to different CBC on the m1 & m2
plane when given a central mchirp value and uncertainty.
It also includes a function that calculates the source frame when given the
detector frame mass and redshift.
"""

from pycbc.conversions import mass2_from_mchirp_mass1 as m2mcm1
from scipy.integrate import quad


def src_mass_from_z_det_mass(z, del_z, mdet, del_mdet):
    """Takes values of redshift, redshift uncertainty, detector mass and its
    uncertainty and computes the source mass and its uncertainty.
    """
    msrc = mdet / (1. + z)
    del_msrc = msrc * ((del_mdet / mdet) ** 2.
                       + (del_z / (1. + z)) ** 2.) ** 0.5
    return (msrc, del_msrc)


# Integration function
def mchange(x, mc):
    """Returns a component mass as a function of mchirp and the other
    component mass.
    """
    return m2mcm1(mc, x)


def intmc(mc, x_min, x_max):
    """Returns the integral of mchange between the minimum and maximum values
    of a component mass taking mchirp as an argument.
    """
    integral = quad(mchange, x_min, x_max, args=mc)
    return integral[0]


def calc_areas(trig_mc_det, mass_limits, mass_bdary, z):
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

    return {
        "bns": abns,
        "gns": agns,
        "nsbh": ansbh,
        "gg": agg,
        "bhg": abhg,
        "bbh": abbh
        }
