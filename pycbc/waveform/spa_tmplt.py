#  Adapted from code in LALSimInspiralTaylorF2.c
#
#  Copyright (C) 2007 Jolien Creighton, B.S. Sathyaprakash, Thomas Cokelaer
#  Copyright (C) 2012 Leo Singer, Alex Nitz
#
#  This program is free software you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with with program see the file COPYING. If not, write to the
#  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
#  MA  02111-1307  USA

"""This module contains functions for generating common SPA template precalculated
   vectors.
"""
from math import sqrt, log
import warnings
import numpy, lal, pycbc.pnutils
from pycbc.scheme import schemed
from pycbc.types import FrequencySeries, Array, complex64, float32, zeros
from pycbc.waveform.utils import ceilpow2

lalsimulation = pycbc.libutils.import_optional('lalsimulation')

def findchirp_chirptime(m1, m2, fLower, porder):
    # variables used to compute chirp time
    m1 = float(m1)
    m2 = float(m2)
    m = m1 + m2
    eta = m1 * m2 / m / m
    c0T = c2T = c3T = c4T = c5T = c6T = c6LogT = c7T = 0.

    # All implemented option
    if porder == -1:
        porder = 7

    if porder >= 7:
        c7T = lal.PI * (14809.0 * eta * eta / 378.0 - 75703.0 * eta / 756.0 - 15419335.0 / 127008.0)

    if porder >= 6:
        c6T = lal.GAMMA * 6848.0 / 105.0 - 10052469856691.0 / 23471078400.0 +\
            lal.PI * lal.PI * 128.0 / 3.0 + \
            eta * (3147553127.0 / 3048192.0 - lal.PI * lal.PI * 451.0 / 12.0) -\
            eta * eta * 15211.0 / 1728.0 + eta * eta * eta * 25565.0 / 1296.0 +\
            eta * eta * eta * 25565.0 / 1296.0 + numpy.log(4.0) * 6848.0 / 105.0
        c6LogT = 6848.0 / 105.0

    if porder >= 5:
        c5T = 13.0 * lal.PI * eta / 3.0 - 7729.0 * lal.PI / 252.0

    if porder >= 4:
        c4T = 3058673.0 / 508032.0 + eta * (5429.0 / 504.0 + eta * 617.0 / 72.0)
        c3T = -32.0 * lal.PI / 5.0
        c2T = 743.0 / 252.0 + eta * 11.0 / 3.0
        c0T = 5.0 * m * lal.MTSUN_SI / (256.0 * eta)

    # This is the PN parameter v evaluated at the lower freq. cutoff
    xT = pow (lal.PI * m * lal.MTSUN_SI * fLower, 1.0 / 3.0)
    x2T = xT * xT
    x3T = xT * x2T
    x4T = x2T * x2T
    x5T = x2T * x3T
    x6T = x3T * x3T
    x7T = x3T * x4T
    x8T = x4T * x4T

    # Computes the chirp time as tC = t(v_low);
    # tC = t(v_low) - t(v_upper) would be more
    # correct, but the difference is negligible.

    # This formula works for any PN order, because
    # higher order coeffs will be set to zero.
    return c0T * (1 + c2T * x2T + c3T * x3T + c4T * x4T + c5T * x5T +
                  (c6T + c6LogT * numpy.log(xT)) * x6T + c7T * x7T) / x8T


def spa_length_in_time(**kwds):
    """
    Returns the length in time of the template,
    based on the masses, PN order, and low-frequency
    cut-off.
    """
    m1 = kwds['mass1']
    m2 = kwds['mass2']
    flow = kwds['f_lower']
    porder = int(kwds['phase_order'])

    # For now, we call the swig-wrapped function below in
    # lalinspiral.  Eventually would be nice to replace this
    # with a function using PN coeffs from lalsimulation.
    return findchirp_chirptime(m1, m2, flow, porder)


def spa_amplitude_factor(**kwds):
    m1 = kwds['mass1']
    m2 = kwds['mass2']

    _, eta = pycbc.pnutils.mass1_mass2_to_mchirp_eta(m1, m2)

    FTaN = 32. * eta * eta / 5.
    dETaN = 2. * -eta / 2.

    M = m1 + m2

    m_sec = M * lal.MTSUN_SI
    piM = lal.PI * m_sec

    amp0 = 4. * m1 * m2 / (1e6 * lal.PC_SI) * lal.MRSUN_SI * lal.MTSUN_SI * sqrt(lal.PI / 12.)

    fac = numpy.sqrt(-dETaN / FTaN) * amp0 * (piM ** (-7./6.))
    return -fac


_prec = None
def spa_tmplt_precondition(length, delta_f, kmin=0):
    """Return the amplitude portion of the TaylorF2 approximant, used to precondition
    the strain data. The result is cached, and so should not be modified, only read.
    """
    global _prec
    if _prec is None or _prec.delta_f != delta_f or len(_prec) < length:
        v = numpy.arange(0, (kmin + length*2), 1.) * delta_f
        v = numpy.power(v[1:len(v)], -7./6.)
        _prec = FrequencySeries(v, delta_f=delta_f, dtype=float32)
    return _prec[kmin:kmin + length]


def spa_tmplt_norm(psd, length, delta_f, f_lower):
    amp = spa_tmplt_precondition(length, delta_f)
    k_min = int(f_lower / delta_f)
    sigma = (amp[k_min:length].numpy() ** 2. / psd[k_min:length].numpy())
    norm_vec = numpy.zeros(length)
    norm_vec[k_min:length] = sigma.cumsum() * 4. * delta_f
    return norm_vec


def spa_tmplt_end(**kwds):
    return pycbc.pnutils.f_SchwarzISCO(kwds['mass1'] + kwds['mass2'])


def spa_distance(psd, mass1, mass2, lower_frequency_cutoff, snr=8):
    """ Return the distance at a given snr (default=8) of the SPA TaylorF2
    template.
    """
    kend = int(spa_tmplt_end(mass1=mass1, mass2=mass2) / psd.delta_f)
    norm1 = spa_tmplt_norm(psd, len(psd), psd.delta_f, lower_frequency_cutoff)
    norm2 = spa_amplitude_factor(mass1=mass1, mass2=mass2) ** 2.0

    if kend >= len(psd):
        kend = len(psd) - 2
    return sqrt(norm1[kend] * norm2) / snr


@schemed("pycbc.waveform.spa_tmplt_")
def spa_tmplt_engine(htilde, kmin, phase_order, delta_f, piM, pfaN,
                     pfa2, pfa3, pfa4, pfa5, pfl5,
                     pfa6, pfl6, pfa7, amp_factor):
    """ Calculate the spa tmplt phase
    """
    err_msg = "This function is a stub that should be overridden using the "
    err_msg += "scheme. You shouldn't be seeing this error!"
    raise ValueError(err_msg)


def spa_tmplt(**kwds):
    """ Generate a minimal TaylorF2 approximant with optimizations for the sin/cos
    """
    distance = kwds['distance']
    mass1 = kwds['mass1']
    mass2 = kwds['mass2']
    s1z = kwds['spin1z']
    s2z = kwds['spin2z']
    phase_order = int(kwds['phase_order'])
    #amplitude_order = int(kwds['amplitude_order'])
    spin_order = int(kwds['spin_order'])

    if 'out' in kwds:
        out = kwds['out']
    else:
        out = None

    amp_factor = spa_amplitude_factor(mass1=mass1, mass2=mass2) / distance

    lal_pars = lal.CreateDict()
    if phase_order != -1:
        lalsimulation.SimInspiralWaveformParamsInsertPNPhaseOrder(
            lal_pars, phase_order)

    if spin_order != -1:
        lalsimulation.SimInspiralWaveformParamsInsertPNSpinOrder(
            lal_pars, spin_order)

    # Calculate the PN terms
    phasing = lalsimulation.SimInspiralTaylorF2AlignedPhasing(
                                    float(mass1), float(mass2),
                                    float(s1z), float(s2z),
                                    lal_pars)

    pfaN = phasing.v[0]
    pfa2 = phasing.v[2] / pfaN
    pfa3 = phasing.v[3] / pfaN
    pfa4 = phasing.v[4] / pfaN
    pfa5 = phasing.v[5] / pfaN
    pfa6 = (phasing.v[6] - phasing.vlogv[6] * log(4)) / pfaN
    pfa7 = phasing.v[7] / pfaN

    pfl5 = phasing.vlogv[5] / pfaN
    pfl6 = phasing.vlogv[6] / pfaN

    piM = lal.PI * (mass1 + mass2) * lal.MTSUN_SI

    if 'sample_points' not in kwds:
        f_lower = kwds['f_lower']
        delta_f = kwds['delta_f']
        kmin = int(f_lower / float(delta_f))

        # Get max frequency one way or another
        # f_final is assigned default value 0 in parameters.py
        if 'f_final' in kwds and kwds['f_final'] > 0.:
            fstop = kwds['f_final']
        elif 'f_upper' in kwds:
            fstop = kwds['f_upper']
            warnings.warn('f_upper is deprecated in favour of f_final!',
                          DeprecationWarning)
        else:
            # Schwarzschild ISCO frequency
            vISCO = 1. / sqrt(6.)
            fstop = vISCO * vISCO * vISCO / piM
        if fstop <= f_lower:
            raise ValueError("cannot generate waveform! f_lower >= f_final"
                             f" ({f_lower}, {fstop})")

        kmax = int(fstop / delta_f)
        f_max = ceilpow2(fstop)
        n = int(f_max / delta_f) + 1

        if not out:
            htilde = FrequencySeries(zeros(n, dtype=numpy.complex64), delta_f=delta_f, copy=False)
        else:
            if type(out) is not Array:
                raise TypeError("Output must be an instance of Array")
            if len(out) < kmax:
                kmax = len(out)
            if out.dtype != complex64:
                raise TypeError("Output array is the wrong dtype")
            htilde = FrequencySeries(out, delta_f=delta_f, copy=False)

        spa_tmplt_engine(htilde[kmin:kmax], kmin, phase_order,
                         delta_f, piM, pfaN,
                         pfa2, pfa3, pfa4, pfa5, pfl5,
                         pfa6, pfl6, pfa7, amp_factor)
    else:
        from .spa_tmplt_cpu import spa_tmplt_inline_sequence
        htilde = numpy.empty(len(kwds['sample_points']), dtype=numpy.complex64)
        spa_tmplt_inline_sequence(
            piM, pfaN, pfa2, pfa3, pfa4, pfa5, pfl5, pfa6, pfl6, pfa7,
            amp_factor, kwds['sample_points'], htilde)

    return htilde


def spa_tmplt_batch(templates_params, filter_length, out, **common_kwds):
    """
    Generate multiple TaylorF2 templates in a single batched operation.
    
    Parameters
    ----------
    templates_params : list of dict
        List of parameter dictionaries, one per template. Each should contain
        mass1, mass2, spin1z, spin2z, and optionally distance.
    filter_length : int
        Length of the output frequency series
    out : list of FrequencySeries
        Pre-allocated output arrays to write templates into
    **common_kwds : dict
        Common parameters for all templates (delta_f, f_lower, etc.)
        
    Returns
    -------
    list of FrequencySeries
        List of generated templates (same as out parameter)
    """
    import cupy as cp
    from pycbc.types import FrequencySeries, zeros
    from pycbc.types.array import Array
    from .spa_tmplt_cupy import spa_tmplt_engine_batch
    from .waveform import get_waveform_filter_length_in_time
    import lalsimulation
    
    num_templates = len(templates_params)
    if num_templates == 0:
        return []
    
    # Extract common parameters
    delta_f = common_kwds['delta_f']
    f_lower = common_kwds['f_lower']
    phase_order = int(common_kwds.get('phase_order', -1))
    spin_order = int(common_kwds.get('spin_order', -1))
    
    # Pre-allocate arrays for parameters
    kmin_list = []
    kmax_list = []
    piM_list = []
    pfaN_list = []
    pfa2_list = []
    pfa3_list = []
    pfa4_list = []
    pfa5_list = []
    pfl5_list = []
    pfa6_list = []
    pfl6_list = []
    pfa7_list = []
    amp_list = []
    phase_order_list = []
    
    lal_pars = lal.CreateDict()
    if phase_order != -1:
        lalsimulation.SimInspiralWaveformParamsInsertPNPhaseOrder(lal_pars, phase_order)
    if spin_order != -1:
        lalsimulation.SimInspiralWaveformParamsInsertPNSpinOrder(lal_pars, spin_order)
    
    # Compute PN coefficients for each template
    for params in templates_params:
        mass1 = params['mass1']
        mass2 = params['mass2']
        s1z = params['spin1z']
        s2z = params['spin2z']
        distance = params.get('distance', common_kwds.get('distance', 1.0))
        
        amp_factor = spa_amplitude_factor(mass1=mass1, mass2=mass2) / distance
        
        # Calculate PN terms
        phasing = lalsimulation.SimInspiralTaylorF2AlignedPhasing(
            float(mass1), float(mass2), float(s1z), float(s2z), lal_pars)
        
        pfaN = phasing.v[0]
        pfa2 = phasing.v[2] / pfaN
        pfa3 = phasing.v[3] / pfaN
        pfa4 = phasing.v[4] / pfaN
        pfa5 = phasing.v[5] / pfaN
        pfa6 = (phasing.v[6] - phasing.vlogv[6] * log(4)) / pfaN
        pfa7 = phasing.v[7] / pfaN
        pfl5 = phasing.vlogv[5] / pfaN
        pfl6 = phasing.vlogv[6] / pfaN
        
        piM = lal.PI * (mass1 + mass2) * lal.MTSUN_SI
        
        kmin = int(f_lower / float(delta_f))
        
        # Get max frequency
        if 'f_final' in common_kwds and common_kwds['f_final'] > 0.:
            fstop = common_kwds['f_final']
        else:
            vISCO = 1. / sqrt(6.)
            fstop = vISCO * vISCO * vISCO / piM
            
        kmax = int(fstop / delta_f)
        # Ensure kmax doesn't exceed filter_length - 1 (final point must be 0)
        if kmax >= filter_length:
            kmax = filter_length - 1
        
        kmin_list.append(kmin)
        kmax_list.append(kmax)
        piM_list.append(piM)
        pfaN_list.append(pfaN)
        pfa2_list.append(pfa2)
        pfa3_list.append(pfa3)
        pfa4_list.append(pfa4)
        pfa5_list.append(pfa5)
        pfl5_list.append(pfl5)
        pfa6_list.append(pfa6)
        pfl6_list.append(pfl6)
        pfa7_list.append(pfa7)
        amp_list.append(amp_factor)
        phase_order_list.append(phase_order)
    
    # Use the provided filter_length for output size
    n = filter_length
    
    # Find the common frequency length for the batch kernel
    # Each template needs different kmax-kmin, so we use the maximum span
    freq_spans = [kmax_list[i] - kmin_list[i] for i in range(num_templates)]
    freq_length = max(freq_spans)
    
    # Transfer parameters to GPU
    kmin_gpu = cp.asarray(kmin_list, dtype=cp.int64)
    freq_spans_gpu = cp.asarray(freq_spans, dtype=cp.int32)
    phase_order_gpu = cp.asarray(phase_order_list, dtype=cp.int64)
    piM_gpu = cp.asarray(piM_list, dtype=cp.float32)
    pfaN_gpu = cp.asarray(pfaN_list, dtype=cp.float32)
    pfa2_gpu = cp.asarray(pfa2_list, dtype=cp.float32)
    pfa3_gpu = cp.asarray(pfa3_list, dtype=cp.float32)
    pfa4_gpu = cp.asarray(pfa4_list, dtype=cp.float32)
    pfa5_gpu = cp.asarray(pfa5_list, dtype=cp.float32)
    pfl5_gpu = cp.asarray(pfl5_list, dtype=cp.float32)
    pfa6_gpu = cp.asarray(pfa6_list, dtype=cp.float32)
    pfl6_gpu = cp.asarray(pfl6_list, dtype=cp.float32)
    pfa7_gpu = cp.asarray(pfa7_list, dtype=cp.float32)
    amp_gpu = cp.asarray(amp_list, dtype=cp.float32)
    
    # Allocate output array on GPU
    htilde_batch_gpu = cp.zeros(num_templates * freq_length, dtype=cp.complex64)
    
    # Call batched kernel
    spa_tmplt_engine_batch(
        htilde_batch_gpu, kmin_gpu, freq_spans_gpu, phase_order_gpu,
        delta_f, piM_gpu, pfaN_gpu,
        pfa2_gpu, pfa3_gpu, pfa4_gpu, pfa5_gpu, pfl5_gpu,
        pfa6_gpu, pfl6_gpu, pfa7_gpu, amp_gpu,
        num_templates, freq_length
    )
    
    # Extract individual templates from batch result
    templates = []
    for i, params in enumerate(templates_params):
        # Use pre-allocated output memory
        htilde_data = out[i]
        htilde = FrequencySeries(htilde_data, delta_f=delta_f, copy=False)
        
        # Clear the output array
        htilde_data.clear()
        
        # Copy the relevant portion from GPU batch to this template
        # Each template occupies freq_length elements in the batch array
        start_idx = i * freq_length
        end_idx = start_idx + freq_spans[i]
        # Extract the slice from GPU and assign to the PyCBC array
        htilde_data._data[kmin_list[i]:kmax_list[i]] = htilde_batch_gpu[start_idx:end_idx]
        
        # Set metadata
        htilde.chirp_length = get_waveform_filter_length_in_time(
            mass1=params['mass1'], mass2=params['mass2'],
            spin1z=params['spin1z'], spin2z=params['spin2z'],
            f_lower=f_lower, phase_order=phase_order,
            approximant='SPAtmplt'
        )
        htilde.length_in_time = htilde.chirp_length
        
        templates.append(htilde)
    
    return templates

