# Copyright (C) 2012  Alex Nitz
#
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


#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#
"""This module contains convenience pN functions. This includes calculating conversions
between quantities.
"""
import logging
import numpy

import lal
from scipy.optimize import bisect, brentq, minimize

from pycbc import conversions, libutils

logger = logging.getLogger('pycbc.pnutils')

lalsim = libutils.import_optional('lalsimulation')

def nearest_larger_binary_number(input_len):
    """ Return the nearest binary number larger than input_len.
    """
    return int(2**numpy.ceil(numpy.log2(input_len)))

def chirp_distance(dist, mchirp, ref_mass=1.4):
    return conversions.chirp_distance(dist, mchirp, ref_mass=ref_mass)

def mass1_mass2_to_mtotal_eta(mass1, mass2):
    m_total = conversions.mtotal_from_mass1_mass2(mass1, mass2)
    eta = conversions.eta_from_mass1_mass2(mass1, mass2)
    return m_total,eta

def mtotal_eta_to_mass1_mass2(m_total, eta):
    mass1 = conversions.mass1_from_mtotal_eta(m_total, eta)
    mass2 = conversions.mass2_from_mtotal_eta(m_total, eta)
    return mass1,mass2

def mass1_mass2_to_mchirp_eta(mass1, mass2):
    m_chirp = conversions.mchirp_from_mass1_mass2(mass1, mass2)
    eta = conversions.eta_from_mass1_mass2(mass1, mass2)
    return m_chirp,eta

def mchirp_eta_to_mass1_mass2(m_chirp, eta):
    mtotal = conversions.mtotal_from_mchirp_eta(m_chirp, eta)
    mass1 = conversions.mass1_from_mtotal_eta(mtotal, eta)
    mass2 = conversions.mass2_from_mtotal_eta(mtotal, eta)
    return mass1, mass2

def mchirp_mass1_to_mass2(mchirp, mass1):
    """
    This function takes a value of mchirp and one component mass and returns
    the second component mass. As this is a cubic equation this requires
    finding the roots and returning the one that is real.
    Basically it can be shown that:

    m2^3 - a(m2 + m1) = 0

    where

    a = Mc^5 / m1^3

    this has 3 solutions but only one will be real.
    """
    return conversions.mass2_from_mchirp_mass1(mchirp, mass1)

def eta_mass1_to_mass2(eta, mass1, return_mass_heavier=False, force_real=True):
    """
    This function takes values for eta and one component mass and returns the
    second component mass. Similar to mchirp_mass1_to_mass2 this requires
    finding the roots of a quadratic equation. Basically:

    eta m2^2 + (2 eta - 1)m1 m2 + eta m1^2 = 0

    This has two solutions which correspond to mass1 being the heavier mass
    or it being the lighter mass. By default the value corresponding to
    mass1 > mass2 is returned. Use the return_mass_heavier kwarg to invert this
    behaviour.
    """
    return conversions.mass_from_knownmass_eta(mass1, eta,
        known_is_secondary=return_mass_heavier, force_real=force_real)

def mchirp_q_to_mass1_mass2(mchirp, q):
    """ This function takes a value of mchirp and the mass ratio
    mass1/mass2 and returns the two component masses.

    The map from q to eta is

        eta = (mass1*mass2)/(mass1+mass2)**2 = (q)/(1+q)**2

    Then we can map from (mchirp,eta) to (mass1,mass2).
    """
    eta = conversions.eta_from_q(q)
    mass1 = conversions.mass1_from_mchirp_eta(mchirp, eta)
    mass2 = conversions.mass2_from_mchirp_eta(mchirp, eta)
    return mass1, mass2

def A0(f_lower):
    """used in calculating chirp times: see Cokelaer, arxiv.org:0706.4437
       appendix 1, also lalinspiral/python/sbank/tau0tau3.py
    """
    return conversions._a0(f_lower)

def A3(f_lower):
    """another parameter used for chirp times"""
    return conversions._a3(f_lower)

def mass1_mass2_to_tau0_tau3(mass1, mass2, f_lower):
    tau0 = conversions.tau0_from_mass1_mass2(mass1, mass2, f_lower)
    tau3 = conversions.tau3_from_mass1_mass2(mass1, mass2, f_lower)
    return tau0,tau3

def tau0_tau3_to_mtotal_eta(tau0, tau3, f_lower):
    mtotal = conversions.mtotal_from_tau0_tau3(tau0, tau3, f_lower)
    eta = conversions.eta_from_tau0_tau3(tau0, tau3, f_lower)
    return mtotal, eta

def tau0_tau3_to_mass1_mass2(tau0, tau3, f_lower):
    m_total,eta = tau0_tau3_to_mtotal_eta(tau0, tau3, f_lower)
    return mtotal_eta_to_mass1_mass2(m_total, eta)

def mass1_mass2_spin1z_spin2z_to_beta_sigma_gamma(mass1, mass2,
                                                  spin1z, spin2z):
    _, eta = mass1_mass2_to_mtotal_eta(mass1, mass2)
    # get_beta_sigma_from_aligned_spins() takes
    # the spin of the heaviest body first
    heavy_spin = numpy.where(mass2 <= mass1, spin1z, spin2z)
    light_spin = numpy.where(mass2 > mass1, spin1z, spin2z)
    beta, sigma, gamma = get_beta_sigma_from_aligned_spins(
        eta, heavy_spin, light_spin)
    return beta, sigma, gamma

def get_beta_sigma_from_aligned_spins(eta, spin1z, spin2z):
    """
    Calculate the various PN spin combinations from the masses and spins.
    See <http://arxiv.org/pdf/0810.5336v3.pdf>.

    Parameters
    -----------
    eta : float or numpy.array
        Symmetric mass ratio of the input system(s)
    spin1z : float or numpy.array
        Spin(s) parallel to the orbit of the heaviest body(ies)
    spin2z : float or numpy.array
        Spin(s) parallel to the orbit of the smallest body(ies)

    Returns
    --------
    beta : float or numpy.array
        The 1.5PN spin combination
    sigma : float or numpy.array
        The 2PN spin combination
    gamma : float or numpy.array
        The 2.5PN spin combination
    chis : float or numpy.array
        (spin1z + spin2z) / 2.
    """
    chiS = 0.5 * (spin1z + spin2z)
    chiA = 0.5 * (spin1z - spin2z)
    delta = (1 - 4 * eta) ** 0.5
    spinspin = spin1z * spin2z
    beta = (113. / 12. - 19. / 3. * eta) * chiS
    beta += 113. / 12. * delta * chiA
    sigma = eta / 48. * (474 * spinspin)
    sigma += (1 - 2 * eta) * (81. / 16. * (chiS * chiS + chiA * chiA))
    sigma += delta * (81. / 8. * (chiS * chiA))
    gamma = (732985. / 2268. - 24260. / 81. * eta - \
            340. / 9. * eta * eta) * chiS
    gamma += (732985. / 2268. + 140. / 9. * eta) * delta * chiA
    return beta, sigma, gamma

def solar_mass_to_kg(solar_masses):
    return solar_masses * lal.MSUN_SI

def parsecs_to_meters(distance):
    return distance * lal.PC_SI

def megaparsecs_to_meters(distance):
    return parsecs_to_meters(distance) * 1e6

def velocity_to_frequency(v, M):
    return conversions.velocity_to_frequency(v, M)

def frequency_to_velocity(f, M):
    return conversions.frequency_to_velocity(f, M)

def f_SchwarzISCO(M):
    """
    Innermost stable circular orbit (ISCO) for a test particle
    orbiting a Schwarzschild black hole

    Parameters
    ----------
    M : float or numpy.array
        Total mass in solar mass units

    Returns
    -------
    f : float or numpy.array
        Frequency in Hz
    """
    return conversions.f_schwarzchild_isco(M)

def f_BKLISCO(m1, m2):
    """
    Mass ratio dependent ISCO derived from estimates of the final spin
    of a merged black hole in a paper by Buonanno, Kidder, Lehner
    (arXiv:0709.3839).  See also arxiv:0801.4297v2 eq.(5)

    Parameters
    ----------
    m1 : float or numpy.array
        First component mass in solar mass units
    m2 : float or numpy.array
        Second component mass in solar mass units

    Returns
    -------
    f : float or numpy.array
        Frequency in Hz
    """
    # q is defined to be in [0,1] for this formula
    q = numpy.minimum(m1/m2, m2/m1)
    return f_SchwarzISCO(m1+m2) * ( 1 + 2.8*q - 2.6*q*q + 0.8*q*q*q )

def f_LightRing(M):
    """
    Gravitational wave frequency corresponding to the light-ring orbit,
    equal to 1/(3**(3/2) pi M) : see InspiralBankGeneration.c

    Parameters
    ----------
    M : float or numpy.array
        Total mass in solar mass units

    Returns
    -------
    f : float or numpy.array
        Frequency in Hz
    """
    return 1.0 / (3.0**(1.5) * lal.PI * M * lal.MTSUN_SI)

def f_ERD(M):
    """
    Effective RingDown frequency studied in Pan et al. (arXiv:0704.1964)
    found to give good fit between stationary-phase templates and
    numerical relativity waveforms [NB equal-mass & nonspinning!]
    Equal to 1.07*omega_220/2*pi

    Parameters
    ----------
    M : float or numpy.array
        Total mass in solar mass units

    Returns
    -------
    f : float or numpy.array
        Frequency in Hz
    """
    return 1.07 * 0.5326 / (2*lal.PI * 0.955 * M * lal.MTSUN_SI)

def f_FRD(m1, m2):
    """
    Fundamental RingDown frequency calculated from the Berti, Cardoso and
    Will (gr-qc/0512160) value for the omega_220 QNM frequency using
    mass-ratio dependent fits to the final BH mass and spin from Buonanno
    et al. (arXiv:0706.3732) : see also InspiralBankGeneration.c

    Parameters
    ----------
    m1 : float or numpy.array
        First component mass in solar mass units
    m2 : float or numpy.array
        Second component mass in solar mass units

    Returns
    -------
    f : float or numpy.array
        Frequency in Hz
    """
    m_total, eta = mass1_mass2_to_mtotal_eta(m1, m2)
    tmp = ( (1. - 0.63*(1. - 3.4641016*eta + 2.9*eta**2)**(0.3)) /
    (1. - 0.057191*eta - 0.498*eta**2) )
    return tmp / (2.*lal.PI * m_total*lal.MTSUN_SI)

def f_LRD(m1, m2):
    """
    Lorentzian RingDown frequency = 1.2*FRD which captures part of
    the Lorentzian tail from the decay of the QNMs

    Parameters
    ----------
    m1 : float or numpy.array
        First component mass in solar mass units
    m2 : float or numpy.array
        Second component mass in solar mass units

    Returns
    -------
    f : float or numpy.array
        Frequency in Hz
    """
    return 1.2 * f_FRD(m1, m2)

def _get_freq(freqfunc, m1, m2, s1z, s2z):
    """Wrapper of the LALSimulation function returning the frequency
    for a given frequency function and template parameters.

    Parameters
    ----------
    freqfunc : lalsimulation FrequencyFunction wrapped object e.g.
        lalsimulation.fEOBNRv2RD
    m1 : float-ish, i.e. castable to float
        First component mass in solar masses
    m2 : float-ish
        Second component mass in solar masses
    s1z : float-ish
        First component dimensionless spin S_1/m_1^2 projected onto L
    s2z : float-ish
        Second component dimensionless spin S_2/m_2^2 projected onto L

    Returns
    -------
    f : float
        Frequency in Hz
    """
    return lalsim.SimInspiralGetFrequency(
        solar_mass_to_kg(m1),
        solar_mass_to_kg(m2),
        0,
        0,
        float(s1z),
        0,
        0,
        float(s2z),
        int(freqfunc)
    )

# vectorize to enable calls with numpy arrays
_vec_get_freq = numpy.vectorize(_get_freq)

def get_freq(freqfunc, m1, m2, s1z, s2z):
    """
    Returns the LALSimulation function which evaluates the frequency
    for the given frequency function and template parameters.

    Parameters
    ----------
    freqfunc : string
        Name of the frequency function to use, e.g., 'fEOBNRv2RD'
    m1 : float or numpy.array
        First component mass in solar masses
    m2 : float or numpy.array
        Second component mass in solar masses
    s1z : float or numpy.array
        First component dimensionless spin S_1/m_1^2 projected onto L
    s2z : float or numpy.array
        Second component dimensionless spin S_2/m_2^2 projected onto L

    Returns
    -------
    f : float or numpy.array
        Frequency in Hz
    """
    lalsim_ffunc = getattr(lalsim, freqfunc)
    return _vec_get_freq(lalsim_ffunc, m1, m2, s1z, s2z)

def _get_final_freq(approx, m1, m2, s1z, s2z):
    """Wrapper of the LALSimulation function returning the final (highest)
    frequency for a given approximant an template parameters

    Parameters
    ----------
    approx : lalsimulation approximant wrapped object e.g.
        lalsimulation.EOBNRv2
    m1 : float-ish, i.e. castable to float
        First component mass in solar masses
    m2 : float-ish
        Second component mass in solar masses
    s1z : float-ish
        First component dimensionless spin S_1/m_1^2 projected onto L
    s2z : float-ish
        Second component dimensionless spin S_2/m_2^2 projected onto L

    Returns
    -------
    f : float
        Frequency in Hz
    """
    return lalsim.SimInspiralGetFinalFreq(
        solar_mass_to_kg(m1),
        solar_mass_to_kg(m2),
        0,
        0,
        float(s1z),
        0,
        0,
        float(s2z),
        int(approx)
    )

# vectorize to enable calls with numpy arrays
_vec_get_final_freq = numpy.vectorize(_get_final_freq)

def get_final_freq(approx, m1, m2, s1z, s2z):
    """Returns the final (highest) frequency for a given approximant using
    given template parameters.

    NOTE: TaylorTx and TaylorFx are currently all given an ISCO cutoff !!

    Parameters
    ----------
    approx : string
        Name of the approximant e.g. 'EOBNRv2'
    m1 : float or numpy.array
        First component mass in solar masses
    m2 : float or numpy.array
        Second component mass in solar masses
    s1z : float or numpy.array
        First component dimensionless spin S_1/m_1^2 projected onto L
    s2z : float or numpy.array
        Second component dimensionless spin S_2/m_2^2 projected onto L

    Returns
    -------
    f : float or numpy.array
        Frequency in Hz
    """
    # Unfortunately we need a few special cases (quite hacky in the case of
    # IMRPhenomXAS) because some useful approximants are not understood by
    # GetApproximantFromString().
    if approx in ['IMRPhenomD', 'IMRPhenomXAS']:
        return frequency_cutoff_from_name('IMRPhenomDPeak', m1, m2, s1z, s2z)
    if approx == 'SEOBNRv5':
        return frequency_cutoff_from_name('SEOBNRv5RD', m1, m2, s1z, s2z)
    lalsim_approx = lalsim.GetApproximantFromString(approx)
    return _vec_get_final_freq(lalsim_approx, m1, m2, s1z, s2z)

# Dictionary of functions with uniform API taking a
# parameter dict indexed on mass1, mass2, spin1z, spin2z
named_frequency_cutoffs = {
    # functions depending on the total mass alone
    "SchwarzISCO": lambda p: f_SchwarzISCO(p["mass1"]+p["mass2"]),
    "LightRing"  : lambda p: f_LightRing(p["mass1"]+p["mass2"]),
    "ERD"        : lambda p: f_ERD(p["mass1"]+p["mass2"]),
    # functions depending on the 2 component masses
    "BKLISCO"    : lambda p: f_BKLISCO(p["mass1"], p["mass2"]),
    "FRD"        : lambda p: f_FRD(p["mass1"], p["mass2"]),
    "LRD"        : lambda p: f_LRD(p["mass1"], p["mass2"]),
    # functions depending on 2 component masses and aligned spins
    "MECO"       : lambda p: meco_frequency(p["mass1"], p["mass2"],
                                              p["spin1z"], p["spin2z"]),
    "HybridMECO" : lambda p: hybrid_meco_frequency(
        p["mass1"], p["mass2"], p["spin1z"], p["spin2z"], qm1=None, qm2=None),
    "IMRPhenomBFinal": lambda p: get_freq("fIMRPhenomBFinal",
                                              p["mass1"], p["mass2"],
                                              p["spin1z"], p["spin2z"]),
    "IMRPhenomCFinal": lambda p: get_freq("fIMRPhenomCFinal",
                                              p["mass1"], p["mass2"],
                                              p["spin1z"], p["spin2z"]),
    "IMRPhenomDPeak": lambda p: get_freq("fIMRPhenomDPeak",
                                              p["mass1"], p["mass2"],
                                              p["spin1z"], p["spin2z"]),
    "EOBNRv2RD"   : lambda p: get_freq("fEOBNRv2RD", p["mass1"], p["mass2"],
                                              p["spin1z"], p["spin2z"]),
    "EOBNRv2HMRD" : lambda p: get_freq("fEOBNRv2HMRD", p["mass1"], p["mass2"],
                                              p["spin1z"], p["spin2z"]),
    "SEOBNRv1RD"  : lambda p: get_freq("fSEOBNRv1RD",  p["mass1"], p["mass2"],
                                              p["spin1z"], p["spin2z"]),
    "SEOBNRv1Peak": lambda p: get_freq("fSEOBNRv1Peak", p["mass1"], p["mass2"],
                                              p["spin1z"], p["spin2z"]),
    "SEOBNRv2RD": lambda p: get_freq("fSEOBNRv2RD", p["mass1"], p["mass2"],
                                     p["spin1z"], p["spin2z"]),
    "SEOBNRv2Peak": lambda p: get_freq("fSEOBNRv2Peak", p["mass1"], p["mass2"],
                                       p["spin1z"], p["spin2z"]),
    "SEOBNRv4RD": lambda p: get_freq("fSEOBNRv4RD", p["mass1"], p["mass2"],
                                     p["spin1z"], p["spin2z"]),
    "SEOBNRv4Peak": lambda p: get_freq("fSEOBNRv4Peak", p["mass1"], p["mass2"],
                                       p["spin1z"], p["spin2z"]),
    "SEOBNRv5RD": lambda p: get_freq("fSEOBNRv5RD", p["mass1"], p["mass2"],
                                     p["spin1z"], p["spin2z"]),
    "SEOBNRv5Peak": lambda p: get_freq("fSEOBNRv5Peak", p["mass1"], p["mass2"],
                                       p["spin1z"], p["spin2z"])
}

def frequency_cutoff_from_name(name, m1, m2, s1z, s2z):
    """
    Returns the result of evaluating the frequency cutoff function
    specified by 'name' on a template with given parameters.

    Parameters
    ----------
    name : string
        Name of the cutoff function
    m1 : float or numpy.array
        First component mass in solar masses
    m2 : float or numpy.array
        Second component mass in solar masses
    s1z : float or numpy.array
        First component dimensionless spin S_1/m_1^2 projected onto L
    s2z : float or numpy.array
        Second component dimensionless spin S_2/m_2^2 projected onto L

    Returns
    -------
    f : float or numpy.array
        Frequency in Hz
    """
    params = {"mass1": m1, "mass2": m2, "spin1z": s1z, "spin2z": s2z}
    return named_frequency_cutoffs[name](params)

def _get_imr_duration(m1, m2, s1z, s2z, f_low, approximant="SEOBNRv4"):
    """Wrapper of lalsimulation template duration approximate formula"""
    m1, m2, s1z, s2z, f_low = float(m1), float(m2), float(s1z), float(s2z),\
                              float(f_low)
    if approximant == "SEOBNRv2":
        chi = lalsim.SimIMRPhenomBComputeChi(m1, m2, s1z, s2z)
        time_length = lalsim.SimIMRSEOBNRv2ChirpTimeSingleSpin(
                                m1 * lal.MSUN_SI, m2 * lal.MSUN_SI, chi, f_low)
    elif approximant == "IMRPhenomXAS":
        time_length = lalsim.SimIMRPhenomXASDuration(
                           m1 * lal.MSUN_SI, m2 * lal.MSUN_SI, s1z, s2z, f_low)
    elif approximant == "IMRPhenomD":
        time_length = lalsim.SimIMRPhenomDChirpTime(
                           m1 * lal.MSUN_SI, m2 * lal.MSUN_SI, s1z, s2z, f_low)
    elif approximant in ["SEOBNRv4", "SEOBNRv4_ROM"]:
        # NB the LALSim function has f_low as first argument
        time_length = lalsim.SimIMRSEOBNRv4ROMTimeOfFrequency(
                           f_low, m1 * lal.MSUN_SI, m2 * lal.MSUN_SI, s1z, s2z)
    elif approximant in ["SEOBNRv5", "SEOBNRv5_ROM"]:
        time_length = lalsim.SimIMRSEOBNRv5ROMTimeOfFrequency(
                           f_low, m1 * lal.MSUN_SI, m2 * lal.MSUN_SI, s1z, s2z)
    elif approximant in ["SPAtmplt", "TaylorF2"]:
        chi = lalsim.SimInspiralTaylorF2ReducedSpinComputeChi(
            m1, m2, s1z, s2z
        )
        time_length = lalsim.SimInspiralTaylorF2ReducedSpinChirpTime(
            f_low, m1 * lal.MSUN_SI, m2 * lal.MSUN_SI, chi, -1
        )
    else:
        raise RuntimeError("I can't calculate a duration for %s" % approximant)
    # FIXME Add an extra factor of 1.1 for 'safety' since the duration
    # functions are approximate
    return time_length * 1.1

get_imr_duration = numpy.vectorize(_get_imr_duration)

def get_inspiral_tf(tc, mass1, mass2, spin1, spin2, f_low, n_points=100,
        pn_2order=7, eccentricity=None, approximant='TaylorF2'):
    """Compute the time-frequency evolution of an inspiral signal.

    Return a tuple of time and frequency vectors tracking the evolution of an
    inspiral signal in the time-frequency plane.
    """
    # handle param-dependent approximant specification
    class Params:
        pass
    params = Params()
    params.mass1 = mass1
    params.mass2 = mass2
    params.spin1z = spin1
    params.spin2z = spin2
    if eccentricity is not None:
        params.eccentricity = eccentricity
    try:
        approximant = eval(approximant, {'__builtins__': None},
                           dict(params=params))
    except (NameError, TypeError):
        pass

    if approximant in ['TaylorF2', 'SPAtmplt']:
        from pycbc.waveform.spa_tmplt import findchirp_chirptime

        # FIXME spins are not taken into account
        f_high = f_SchwarzISCO(mass1 + mass2)
        def tof_func(f):
            return findchirp_chirptime(
                float(mass1),
                float(mass2),
                float(f),
                pn_2order
            )
    elif approximant.startswith('SEOBNRv'):
        approximant_prefix = approximant[:len('SEOBNRv*')]
        f_high = get_final_freq(approximant_prefix, mass1, mass2, spin1, spin2)
        f_high *= 0.999  # avoid errors due to rounding
        tof_func_map = {
            # use HI function for v2 as it has wider freq range validity
            'SEOBNRv2': lalsim.SimIMRSEOBNRv2ROMDoubleSpinHITimeOfFrequency,
            'SEOBNRv4': lalsim.SimIMRSEOBNRv4ROMTimeOfFrequency,
            'SEOBNRv5': lalsim.SimIMRSEOBNRv5ROMTimeOfFrequency
        }
        def tof_func(f):
            return tof_func_map[approximant_prefix](
                f,
                solar_mass_to_kg(mass1),
                solar_mass_to_kg(mass2),
                float(spin1),
                float(spin2)
            )
    elif approximant in ['IMRPhenomD', 'IMRPhenomXAS']:
        f_high = get_final_freq(approximant, mass1, mass2, spin1, spin2)
        tof_func_map = {
            'IMRPhenomD': lalsim.SimIMRPhenomDChirpTime,
            'IMRPhenomXAS': lalsim.SimIMRPhenomXASDuration
        }
        def tof_func(f):
            return tof_func_map[approximant](
                solar_mass_to_kg(mass1),
                solar_mass_to_kg(mass2),
                float(spin1),
                float(spin2),
                f
            )
    elif approximant=='TaylorF2Ecc':
        f_high = f_SchwarzISCO(mass1 + mass2)
        from pycbc.waveform.spa_tmplt import eccentric_chirp_time
        def tof_func(f):
            return eccentric_chirp_time(
                    float(mass1),
                    float(mass2),
                    float(eccentricity),
                    f)
    else:
        raise ValueError(f'Approximant {approximant} not supported')
    track_f = numpy.logspace(numpy.log10(f_low), numpy.log10(f_high), n_points)
    tof_func_vec = numpy.vectorize(tof_func)
    track_t = tc - tof_func_vec(track_f)
    return (track_t, track_f)

# Functions adapted from https://github.com/gmorras/pyEFPE/blob/main/pyEFPE/waveform/functions.py
#      Original function --->  Functions here
#      tLO_func          --->  eccentric_newtonian_time
#      F_tLO_series      --->  _eccentric_enhancement_function
#      F_tLO_series_at_0 --->  _enhancement_function_low
#      F_tLO_series_at_1 --->  _enhancement_function_high

def eccentric_newtonian_time(m1, m2, eccentricity, f_low):
    """Returns the Newtonian eccentric-binary time to coalescence.

    This computes the leading-order inspiral time for an eccentric compact
    binary using Eq. (119) of Morras et al., arXiv:2502.03929. The
    eccentricity correction is evaluated with
    :func:`_eccentric_enhancement_function`.

    Parameters
    ----------
    m1, m2 : float or array_like
        Component masses in solar masses. Inputs must be broadcastable to a
        common shape.
    eccentricity : float or array_like
        Orbital eccentricity defined at ``f_low``.
    f_low : float or array_like
        Starting gravitational-wave frequency in Hz.

    Returns
    -------
    float or numpy.ndarray
        Newtonian time to coalescence in seconds. Scalar inputs return a
        scalar; array-like inputs return a NumPy array.
    """
    m1, m2, eccentricity, f_low = numpy.broadcast_arrays(
        numpy.asarray(m1, dtype=float),
        numpy.asarray(m2, dtype=float),
        numpy.asarray(eccentricity, dtype=float),
        numpy.asarray(f_low, dtype=float),
    )

    M = m1 + m2
    ecc_sq = eccentricity * eccentricity
    one_minus_ecc_sq = 1 - ecc_sq
    one_minus_ecc_sq_sqrt = numpy.sqrt(one_minus_ecc_sq)
    velocity = frequency_to_velocity(f_low, M)
    # weighted_velocity is the `y` variable in Morras et al., arXiv:2502.03929
    weighted_velocity = velocity / one_minus_ecc_sq_sqrt
    eta = conversions.eta_from_mass1_mass2(m1, m2)
    M_sec = M * lal.MTSUN_SI
    t_N = (5 / (256 * one_minus_ecc_sq_sqrt * eta) * M_sec *
            (weighted_velocity**-8) * _eccentric_enhancement_function(ecc_sq))
    return t_N.item() if t_N.ndim == 0 else t_N


def _eccentric_enhancement_function(eccentricity_sq):

    eccentricity_sq_arr = numpy.asarray(eccentricity_sq, dtype=float)

    # Value of eccentricity^2 to decide which enhancement function to be used,
    # based on Morras et al., arXiv:2502.03929
    eccentricity_sq_switch = 0.4

    #if eccentricity^2 <= 0.4  use series expansion from `_enhancement_function_low`,
    #                otherwise use series expansion from `_enhancement_function_high`

    if eccentricity_sq_arr.ndim == 0:
        eccentricity_sq_scalar = float(eccentricity_sq_arr)
        if eccentricity_sq_scalar <= eccentricity_sq_switch:
            return _enhancement_function_low(eccentricity_sq_scalar)
        else:
            return _enhancement_function_high(eccentricity_sq_scalar)

    else:
        F = numpy.empty_like(eccentricity_sq_arr, dtype=float)
        low_ecc_index = eccentricity_sq_arr <= eccentricity_sq_switch
        high_ecc_index = numpy.logical_not(low_ecc_index)
        if numpy.any(low_ecc_index):
            F[low_ecc_index] = _enhancement_function_low(eccentricity_sq_arr[low_ecc_index])
        if numpy.any(high_ecc_index):
            F[high_ecc_index] = _enhancement_function_high(eccentricity_sq_arr[high_ecc_index])
        return F

#coefficeints of polynomials for `_enhancement_function_low`
_LOW_ECC_COEFFS = numpy.array([
    1.000000000000000, -0.1511627906976744, 0.2656836084021005, 0.007463780007501875,
    0.08800790590714085, 0.03153077124184580, 0.04185392210761341, 0.02761371124737777,
    0.02642686119635599, 0.02145414169943131, 0.01926563742403364, 0.01677217450181226,
    0.01503898450331388, 0.01347003088516929, 0.01220348405026613, 0.01110346195121149,
    0.01016749330223870, 0.009352447459389354, 0.008642114232972108, 0.008016709107501086,
    0.007463457161231162, 0.006970902964332086, 0.006530253812462707, 0.006134111124121483,
    0.005776451398258840, 0.005452229768143720, 0.005157232928828976, 0.004887901117379935,
    0.004641215172072290, 0.004414596725230578, 0.004205833180162198, 0.004013015784562471,
    0.003834490347048246, 0.003668816871424952, 0.003514736646109113, 0.003371145103099870,
    0.003237069368178693, 0.003111649567641214, 0.002994123194217794, 0.002883811964918853,
    0.002780110725151045, 0.002682478039498533, 0.002590428180410570, 0.002503524280447905,
    0.002421372457429333, 0.002343616756405161, 0.002269934780185545, 0.002200033902499887,
    0.002133647975961232, 0.002070534461714599, 0.002010471919657125,
])

#coefficeints of polynomials for `_enhancement_function_high`
_HIGH_ECC_COEFFS = numpy.array([
    0.40941176470588236, 0.02286320645905421, 0.008142925951557094, 0.004155878512401501,
    0.0024847568765827364, 0.001633290511588348, 0.0011455395465032605, 0.000842577094447224,
    0.0006427195446513564, 0.0005045715814217113, 0.00040544477831148924, 0.0003321116545345162,
    0.0002764636962467965, 0.00023331895675436905, 0.0001992473575067662, 0.00017190919726595898,
    0.00014966654751355843, 0.000131346469484909, 0.00011609207361015848, 0.00010326617525262309,
    9.238740896587563e-05, 8.308691874613337e-05, 7.50784086790562e-05, 6.813705814151314e-05,
    6.20844345948999e-05, 5.677753690123865e-05, 5.210072977646673e-05, 4.7959732146031274e-05,
    4.427708468062032e-05, 4.0988696117937945e-05, 3.80411855904995e-05, 3.538981870077757e-05,
    3.2996890965966614e-05, 3.083045152872919e-05, 2.886328796018351e-05, 2.707211306399751e-05,
    2.543690918050699e-05, 2.3940396192787112e-05, 2.256759736012561e-05, 2.1305483020854193e-05,
    2.014267666038337e-05, 1.906921121897629e-05, 1.8076326095558655e-05, 1.7156297290322297e-05,
    1.6302294667351972e-05, 1.550826151746742e-05, 1.4768812541410176e-05, 1.407914711455732e-05,
])

def _enhancement_function_low(eccentricity_sq):
    # Maclaurin seies of Eq. D6. Morras et al., arXiv:2502.03929
    return numpy.polyval(numpy.flip(_LOW_ECC_COEFFS), eccentricity_sq)

def _enhancement_function_high(eccentricity_sq):
    # Series expansion of Eq. D7, Morras et al., arXiv:2502.03929
    # around sqrt(1 - eccentricity^2) = 0
    cns = _HIGH_ECC_COEFFS
    pre_factor = (48/19)*((425/304)**(1181/2299))
    u = 1 - eccentricity_sq
    eccentric_factor = eccentricity_sq**(-24 / 19) * (1 + (121 / 304) * eccentricity_sq)**(-3480 / 2299)
    eccentric_factor *= 1 - 1.4555165803216864*numpy.sqrt(u) + u*numpy.polyval(numpy.flip(cns),u)
    return pre_factor * eccentric_factor

##############################This code was taken from Andy ###########


def _energy_coeffs(m1, m2, chi1, chi2):
    """ Return the center-of-mass energy coefficients up to 3.0pN (2.5pN spin)
    """
    mtot = m1 + m2
    eta = m1*m2 / (mtot*mtot)
    chi = (m1*chi1 + m2*chi2) / mtot
    chisym = (chi1 + chi2) / 2.
    beta = (113.*chi - 76.*eta*chisym)/12.
    sigma12 = 79.*eta*chi1*chi2/8.
    sigmaqm = 81.*m1*m1*chi1*chi1/(16.*mtot*mtot) \
            + 81.*m2*m2*chi2*chi2/(16.*mtot*mtot)

    energy0 = -0.5*eta
    energy2 = -0.75 - eta/12.
    energy3 = 0.
    energy4 = -3.375 + (19*eta)/8. - pow(eta,2)/24.
    energy5 = 0.
    energy6 = -10.546875 - (155*pow(eta,2))/96. - (35*pow(eta,3))/5184. \
                + eta*(59.80034722222222 - (205*pow(lal.PI,2))/96.)

    energy3 += (32*beta)/113. + (52*chisym*eta)/113.

    energy4 += (-16*sigma12)/79. - (16*sigmaqm)/81.
    energy5 += (96*beta)/113. + ((-124*beta)/339. - (522*chisym)/113.)*eta \
                - (710*chisym*pow(eta,2))/339.

    return (energy0, energy2, energy3, energy4, energy5, energy6)

def meco_velocity(m1, m2, chi1, chi2):
    """
    Returns the velocity of the minimum energy cutoff for 3.5pN (2.5pN spin)

    Parameters
    ----------
    m1 : float
        First component mass in solar masses
    m2 : float
        Second component mass in solar masses
    chi1 : float
        First component dimensionless spin S_1/m_1^2 projected onto L
    chi2 : float
        Second component dimensionless spin S_2/m_2^2 projected onto L

    Returns
    -------
    v : float
        Velocity (dimensionless)
    """
    _, energy2, energy3, energy4, energy5, energy6 = \
        _energy_coeffs(m1, m2, chi1, chi2)
    def eprime(v):
        return 2. + v * v * (4.*energy2 + v * (5.*energy3 \
                + v * (6.*energy4
                + v * (7.*energy5 + 8.*energy6 * v))))
    return bisect(eprime, 0.05, 1.0)

def _meco_frequency(m1, m2, chi1, chi2):
    """Returns the frequency of the minimum energy cutoff for 3.5pN (2.5pN spin)
    """
    return velocity_to_frequency(meco_velocity(m1, m2, chi1, chi2), m1+m2)

meco_frequency = numpy.vectorize(_meco_frequency)

def _dtdv_coeffs(m1, m2, chi1, chi2):
    """ Returns the dt/dv coefficients up to 3.5pN (2.5pN spin)
    """
    mtot = m1 + m2
    eta = m1*m2 / (mtot*mtot)
    chi = (m1*chi1 + m2*chi2) / mtot
    chisym = (chi1 + chi2) / 2.
    beta = (113.*chi - 76.*eta*chisym)/12.
    sigma12 = 79.*eta*chi1*chi2/8.
    sigmaqm = 81.*m1*m1*chi1*chi1/(16.*mtot*mtot) \
            + 81.*m2*m2*chi2*chi2/(16.*mtot*mtot)

    dtdv0 = 1. # FIXME: Wrong but doesn't matter for now.
    dtdv2 = (1./336.) * (743. + 924.*eta)
    dtdv3 = -4. * lal.PI + beta
    dtdv4 = (3058673. + 5472432.*eta + 4353552.*eta*eta)/1016064. - sigma12 - sigmaqm
    dtdv5 = (1./672.) * lal.PI * (-7729. + 1092.*eta) + (146597.*beta/18984. + 42.*beta*eta/113. - 417307.*chisym*eta/18984. - 1389.*chisym*eta*eta/226.)
    dtdv6 = 22.065 + 165.416*eta - 2.20067*eta*eta + 4.93152*eta*eta*eta
    dtdv6log = 1712./315.
    dtdv7 = (lal.PI/1016064.) * (-15419335. - 12718104.*eta + 4975824.*eta*eta)

    return (dtdv0, dtdv2, dtdv3, dtdv4, dtdv5, dtdv6, dtdv6log, dtdv7)

def _dtdv_cutoff_velocity(m1, m2, chi1, chi2):
    _, dtdv2, dtdv3, dtdv4, dtdv5, dtdv6, dtdv6log, dtdv7 = _dtdv_coeffs(m1, m2, chi1, chi2)

    def dtdv_func(v):
        x = dtdv7
        x = v * x + dtdv6 + dtdv6log * 3. * numpy.log(v)
        x = v * x + dtdv5
        x = v * x + dtdv4
        x = v * x + dtdv3
        x = v * x + dtdv2
        return v * v * x + 1.

    if dtdv_func(1.0) < 0.:
        return bisect(dtdv_func, 0.05, 1.0)
    else:
        return 1.0

def energy_coefficients(m1, m2, s1z=0, s2z=0, phase_order=-1, spin_order=-1):
    """ Return the energy coefficients. This assumes that the system has aligned spins only.
    """
    implemented_phase_order = 7
    implemented_spin_order = 7
    if phase_order > implemented_phase_order:
        raise ValueError("pN coeffiecients of that order have not been implemented")
    elif phase_order == -1:
        phase_order = implemented_phase_order

    if spin_order > implemented_spin_order:
        raise ValueError("pN coeffiecients of that order have not been implemented")
    elif spin_order == -1:
        spin_order = implemented_spin_order

    qmdef1 = 1.0
    qmdef2 = 1.0

    M = m1 + m2
    dm = (m1-m2)/M
    m1M = m1 / M
    m2M = m2 / M

    s1z = s1z * m1M * m1M
    s2z = s2z * m2M * m2M

    _, eta = mass1_mass2_to_mchirp_eta(m1, m2)

    ecof = numpy.zeros(phase_order+1)
    # Orbital terms
    if phase_order >= 0:
        ecof[0] = 1.0
    if phase_order >= 1:
        ecof[1] = 0
    if phase_order >= 2:
        ecof[2] = -(1.0/12.0) * (9.0 + eta)
    if phase_order >= 3:
        ecof[3] = 0
    if phase_order >= 4:
        ecof[4] = (-81.0 + 57.0*eta - eta*eta) / 24.0
    if phase_order >= 5:
        ecof[5] = 0
    if phase_order >= 6:
        ecof[6] = - 675.0/64.0 + ( 34445.0/576.0    \
              - 205.0/96.0 * lal.PI * lal.PI ) * eta  \
              - (155.0/96.0) *eta * eta - 35.0/5184.0 * eta * eta
    # Spin terms

    ESO15s1 = 8.0/3.0 + 2.0*m2/m1
    ESO15s2 = 8.0/3.0 + 2.0*m1/m2

    ESS2 = 1.0 / eta
    EQM2s1 = qmdef1/2.0/m1M/m1M
    EQM2s1L = -qmdef1*3.0/2.0/m1M/m1M
    #EQM2s2 = qmdef2/2.0/m2M/m2M
    EQM2s2L = -qmdef2*3.0/2.0/m2M/m2M

    ESO25s1 = 11.0 - 61.0*eta/9.0 + (dm/m1M) * (-3.0 + 10.*eta/3.0)
    ESO25s2 = 11.0 - 61.0*eta/9.0 + (dm/m2M) * (3.0 - 10.*eta/3.0)

    ESO35s1 = 135.0/4.0 - 367.0*eta/4.0 + 29.0*eta*eta/12.0 + (dm/m1M) * (-27.0/4.0 + 39.0*eta - 5.0*eta*eta/4.0)
    ESO35s2 = 135.0/4.0 - 367.0*eta/4.0 + 29.0*eta*eta/12.0 - (dm/m2M) * (-27.0/4.0 + 39.0*eta - 5.0*eta*eta/4.0)

    if spin_order >=3:
        ecof[3] += ESO15s1 * s1z + ESO15s2 * s2z
    if spin_order >=4:
        ecof[4] += ESS2 * (s1z*s2z - 3.0*s1z*s2z)
        ecof[4] += EQM2s1*s1z*s1z + EQM2s1*s2z*s2z + EQM2s1L*s1z*s1z + EQM2s2L*s2z*s2z
    if spin_order >=5:
        ecof[5] = ESO25s1*s1z + ESO25s2*s2z
    if spin_order >=7:
        ecof[7] += ESO35s1*s1z + ESO35s2*s2z

    return ecof

def energy(v, mass1, mass2, s1z=0, s2z=0, phase_order=-1, spin_order=-1):
    ecof = energy_coefficients(mass1, mass2, s1z, s2z, phase_order, spin_order)
    _, eta = mass1_mass2_to_mchirp_eta(mass1, mass2)
    amp = - (1.0/2.0) * eta
    e = 0.0
    for i in numpy.arange(0, len(ecof), 1):
        e += v**(i+2.0) * ecof[i]

    return e * amp

def meco2(m1, m2, s1z=0, s2z=0, phase_order=-1, spin_order=-1):
    ecof = energy_coefficients(m1, m2, s1z, s2z, phase_order, spin_order)

    def test(v):
        de = 0
        for i in numpy.arange(0, len(ecof), 1):
            de += v**(i+1.0)* ecof[i] * (i + 2)

        return de

    return bisect(test, 0.001, 1.0)


def t2_cutoff_velocity(m1, m2, chi1, chi2):
    return min(meco_velocity(m1,m2,chi1,chi2), _dtdv_cutoff_velocity(m1,m2,chi1,chi2))

def t2_cutoff_frequency(m1, m2, chi1, chi2):
    return velocity_to_frequency(t2_cutoff_velocity(m1, m2, chi1, chi2), m1 + m2)

t4_cutoff_velocity = meco_velocity
t4_cutoff_frequency = meco_frequency

# Hybrid MECO in arXiv:1602.03134
# To obtain the MECO, find minimum in v of eq. (6)


def kerr_lightring(v, chi):
    """Return the function whose first root defines the Kerr light ring"""
    return 1 + chi * v**3 - 3 * v**2 * (1 - chi * v**3)**(1./3)


def kerr_lightring_velocity(chi):
    """Return the velocity at the Kerr light ring"""
    # If chi > 0.9996, the algorithm cannot solve the function
    if chi >= 0.9996:
        return brentq(kerr_lightring, 0, 0.8, args=(0.9996))
    else:
        return brentq(kerr_lightring, 0, 0.8, args=(chi))


def hybridEnergy(v, m1, m2, chi1, chi2, qm1, qm2):
    """Return hybrid MECO energy.

    Return the hybrid energy [eq. (6)] whose minimum defines the hybrid MECO
    up to 3.5PN (including the 3PN spin-spin)

    Parameters
    ----------
    m1 : float
        Mass of the primary object in solar masses.
    m2 : float
        Mass of the secondary object in solar masses.
    chi1: float
        Dimensionless spin of the primary object.
    chi2: float
        Dimensionless spin of the secondary object.
    qm1: float
        Quadrupole-monopole term of the primary object (1 for black holes).
    qm2: float
        Quadrupole-monopole term of the secondary object (1 for black holes).

    Returns
    -------
    h_E: float
        The hybrid energy as a function of v
    """
    pi_sq = numpy.pi**2
    v2, v3, v4, v5, v6, v7 = v**2, v**3, v**4, v**5, v**6, v**7
    chi1_sq, chi2_sq = chi1**2, chi2**2
    m1, m2 = float(m1), float(m2)
    M = float(m1 + m2)
    M_2, M_4 = M**2, M**4
    eta = m1 * m2 / M_2
    eta2, eta3 = eta**2, eta**3
    m1_2, m1_4 = m1**2, m1**4
    m2_2, m2_4 = m2**2, m2**4

    chi = (chi1 * m1 + chi2 * m2) / M
    Kerr = -1. + (1. - 2. * v2 * (1. - chi * v3)**(1./3.)) / \
        numpy.sqrt((1. - chi * v3) * (1. + chi * v3 - 3. * v2 * (1 - chi * v3)**(1./3.)))

    h_E = Kerr - \
        (v2 / 2.) * \
        (
        - eta * v2 / 12. - 2 * (chi1 + chi2) * eta * v3 / 3. +
        (19. * eta / 8. - eta2 / 24. + chi1_sq * m1_2 * (1 - qm1) / M_2 +
            chi2_sq * m2_2 * (1 - qm2) / M_2) * v4
        - 1. / 9. * (120. * (chi1 + chi2) * eta2 +
            (76. * chi1 + 45. * chi2) * m1_2 * eta / M_2 +
            (45. * chi1 + 76. * chi2) * m2_2 * eta / M_2) * v5
        + (34445. * eta / 576. - 205. * pi_sq * eta / 96. - 155. * eta2 / 96. -
            35. * eta3 / 5184. +
            5. / 18. * (21. * chi1_sq * (1. - qm1) * m1_4 / M_4 +
            21. * chi2_sq * (1. - qm2) * m2_4 / M_4 +
            (chi1_sq * (56. - 27. * qm1) + 20. * chi1 * chi2) * eta * m1_2 / M_2 +
            (chi2_sq * (56. - 27. * qm2) + 20. * chi1 * chi2) * eta * m2_2 / M_2 +
            (chi1_sq * (31. - 9. * qm1) + 38. * chi1 * chi2 +
            chi2_sq * (31. - 9. * qm2)) * eta2)) * v6
        - eta / 12. * (3. * (292. * chi1 + 81. * chi2) * m1_4 / M_4 +
            3. * (81. * chi1 + 292. * chi2) * m2_4 / M_4 +
            4. * (673. * chi1 + 360. * chi2) * eta * m1_2 / M_2 +
            4. * (360. * chi1 + 673. * chi2) * eta * m2_2 / M_2 +
            3012. * eta2 * (chi1 + chi2)) * v7
        )

    return h_E


def hybrid_meco_velocity(m1, m2, chi1, chi2, qm1=None, qm2=None):
    """Return the velocity of the hybrid MECO

    Parameters
    ----------
    m1 : float
        Mass of the primary object in solar masses.
    m2 : float
        Mass of the secondary object in solar masses.
    chi1: float
        Dimensionless spin of the primary object.
    chi2: float
        Dimensionless spin of the secondary object.
    qm1: {None, float}, optional
        Quadrupole-monopole term of the primary object (1 for black holes).
        If None, will be set to qm1 = 1.
    qm2: {None, float}, optional
        Quadrupole-monopole term of the secondary object (1 for black holes).
        If None, will be set to qm2 = 1.

    Returns
    -------
    v: float
        The velocity (dimensionless) of the hybrid MECO
    """

    if qm1 is None:
        qm1 = 1
    if qm2 is None:
        qm2 = 1

    # Set bounds at 0.1 to skip v=0 and at the lightring velocity
    chi = (chi1 * m1 + chi2 * m2) / (m1 + m2)
    vmax = kerr_lightring_velocity(chi) - 0.01

    return minimize(hybridEnergy, 0.2, args=(m1, m2, chi1, chi2, qm1, qm2),
                    bounds=[(0.1, vmax)]).x.item()


def hybrid_meco_frequency(m1, m2, chi1, chi2, qm1=None, qm2=None):
    """Return the frequency of the hybrid MECO

    Parameters
    ----------
    m1 : float
        Mass of the primary object in solar masses.
    m2 : float
        Mass of the secondary object in solar masses.
    chi1: float
        Dimensionless spin of the primary object.
    chi2: float
        Dimensionless spin of the secondary object.
    qm1: {None, float}, optional
        Quadrupole-monopole term of the primary object (1 for black holes).
        If None, will be set to qm1 = 1.
    qm2: {None, float}, optional
        Quadrupole-monopole term of the secondary object (1 for black holes).
        If None, will be set to qm2 = 1.

    Returns
    -------
    f: float
        The frequency (in Hz) of the hybrid MECO
    """
    if qm1 is None:
        qm1 = 1
    if qm2 is None:
        qm2 = 1

    return velocity_to_frequency(hybrid_meco_velocity(m1, m2, chi1, chi2, qm1, qm2), m1 + m2)


def jframe_to_l0frame(mass1, mass2, f_ref, phiref=0., thetajn=0., phijl=0.,
                      spin1_a=0., spin2_a=0.,
                      spin1_polar=0., spin2_polar=0.,
                      spin12_deltaphi=0.):
    """Converts J-frame parameters into L0 frame.

    Parameters
    ----------
    mass1 : float
        The mass of the first component object in the
        binary (in solar masses)
    mass2 : float
        The mass of the second component object in the
        binary (in solar masses)
    f_ref : float
        The reference frequency.
    phiref : float
        The reference phase.
    thetajn : float
        Angle between the line of sight and the total angular momentume J.
    phijl : float
        Azimuthal angle of L on its cone about J.
    spin1_a : float
        The dimensionless spin magnitude :math:`|\\vec{{s}}_1/m^2_1|`.
    spin2_a : float
        The dimensionless spin magnitude :math:`|\\vec{{s}}_2/m^2_2|`.
    spin1_polar : float
        Angle between L and the spin magnitude of the larger object.
    spin2_polar : float
        Angle betwen L and the spin magnitude of the smaller object.
    spin12_deltaphi : float
        Difference between the azimuthal angles of the spin of the larger
        object (S1) and the spin of the smaller object (S2).

    Returns
    -------
    dict :
        Dictionary of:

        * inclination : float
            Inclination (rad), defined as the angle between
            the orbital angular momentum L and the
            line-of-sight at the reference frequency.
        * spin1x : float
            The x component of the first binary component's
            dimensionless spin.
        * spin1y : float
            The y component of the first binary component's
            dimensionless spin.
        * spin1z : float
            The z component of the first binary component's
            dimensionless spin.
        * spin2x : float
            The x component of the second binary component's
            dimensionless spin.
        * spin2y : float
            The y component of the second binary component's
            dimensionless spin.
        * spin2z : float
            The z component of the second binary component's
            dimensionless spin.
    """
    inclination, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z = \
        lalsim.SimInspiralTransformPrecessingNewInitialConditions(
            thetajn, phijl, spin1_polar, spin2_polar, spin12_deltaphi,
            spin1_a, spin2_a, mass1*lal.MSUN_SI, mass2*lal.MSUN_SI, f_ref,
            phiref)
    out = {'inclination': inclination,
           'spin1x': spin1x,
           'spin1y': spin1y,
           'spin1z': spin1z,
           'spin2x': spin2x,
           'spin2y': spin2y,
           'spin2z': spin2z}
    return out

def l0frame_to_jframe(mass1, mass2, f_ref, phiref=0., inclination=0.,
                      spin1x=0., spin1y=0., spin1z=0.,
                      spin2x=0., spin2y=0., spin2z=0.):
    """Converts L0-frame parameters to J-frame.

    Parameters
    ----------
    mass1 : float
        The mass of the first component object in the
        binary (in solar masses)
    mass2 : float
        The mass of the second component object in the
        binary (in solar masses)
    f_ref : float
        The reference frequency.
    phiref : float
        The orbital phase at ``f_ref``.
    inclination : float
        Inclination (rad), defined as the angle between
        the orbital angular momentum L and the
        line-of-sight at the reference frequency.
    spin1x : float
        The x component of the first binary component's
        dimensionless spin.
    spin1y : float
        The y component of the first binary component's
        dimensionless spin.
    spin1z : float
        The z component of the first binary component's
        dimensionless spin.
    spin2x : float
        The x component of the second binary component's
        dimensionless spin.
    spin2y : float
        The y component of the second binary component's
        dimensionless spin.
    spin2z : float
        The z component of the second binary component's
        dimensionless spin.

    Returns
    -------
    dict :
        Dictionary of:

        * thetajn : float
            Angle between the line of sight and the total angular momentume J.
        * phijl : float
            Azimuthal angle of L on its cone about J.
        * spin1_a : float
            The dimensionless spin magnitude :math:`|\\vec{{s}}_1/m^2_1|`.
        * spin2_a : float
            The dimensionless spin magnitude :math:`|\\vec{{s}}_2/m^2_2|`.
        * spin1_polar : float
            Angle between L and the spin magnitude of the larger object.
        * spin2_polar : float
            Angle betwen L and the spin magnitude of the smaller object.
        * spin12_deltaphi : float
            Difference between the azimuthal angles of the spin of the larger
            object (S1) and the spin of the smaller object (S2).
    """
    # Note: unlike other LALSimulation functions, this one takes masses in
    # solar masses
    thetajn, phijl, s1pol, s2pol, s12_deltaphi, spin1_a, spin2_a = \
        lalsim.SimInspiralTransformPrecessingWvf2PE(
            inclination, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z,
            mass1, mass2, f_ref, phiref)
    out = {'thetajn': thetajn,
           'phijl': phijl,
           'spin1_polar': s1pol,
           'spin2_polar': s2pol,
           'spin12_deltaphi': s12_deltaphi,
           'spin1_a': spin1_a,
           'spin2_a': spin2_a}
    return out
