# Copyright (C) 2015  Alex Nitz
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
""" Utilities for handling frequency compressed an unequally spaced frequency
domain waveforms.
"""
import lalsimulation, lal, numpy
from pycbc import pnutils
from pycbc.types import Array
from pycbc.opt import omp_libs, omp_flags
from pycbc import WEAVE_FLAGS
from scipy.weave import inline
from scipy import interpolate
from pycbc.types import FrequencySeries, zeros, complex_same_precision_as
from pycbc.waveform import utils

def rough_time_estimate(m1, m2, flow, fudge_length=1.1, fudge_min=0.02):
    """ A very rough estimate of the duration of the waveform.

    An estimate of the waveform duration starting from flow. This is intended
    to be fast but not necessarily accurate. It should be an overestimate of
    the length. It is derived from a simplification of the 0PN post-newtonian
    terms and includes a fudge factor for possible ringdown, etc.

    Parameters
    ----------
    m1: float
        mass of first component object in solar masses
    m2: float
        mass of second component object in solar masses
    flow: float
        starting frequency of the waveform
    fudge_length: optional, {1.1, float}
        Factor to multiply length estimate by to ensure it is a convservative
        value
    fudge_min: optional, {0.2, float}
        Minimum signal duration that can be returned. This should be long
        enough to encompass the ringdown and errors in the precise end time.

    Returns
    -------
    time: float
        Time from flow untill the end of the waveform
    """
    m = m1 + m2
    msun = m * lal.MTSUN_SI
    t =  5.0 / 256.0 * m * m * msun / (m1 * m2) / \
        (numpy.pi * msun * flow) **  (8.0 / 3.0)

    # fudge factoriness
    return .022 if t < 0 else (t + fudge_min) * fudge_length 

def mchirp_compression(m1, m2, fmin, fmax, min_seglen=0.02, df_multiple=None):
    """Return the frequencies needed to compress a waveform with the given
    chirp mass. This is based on the estimate in rough_time_estimate.

    Parameters
    ----------
    m1: float
        mass of first component object in solar masses
    m2: float
        mass of second component object in solar masses
    fmin : float
        The starting frequency of the compressed waveform.
    fmax : float
        The ending frequency of the compressed waveform.
    min_seglen : float
        The inverse of this gives the maximum frequency step that is used.
    df_multiple : {None, float}
        Make the compressed sampling frequencies a multiple of the given value.
        If None provided, the returned sample points can have any floating
        point value.

    Returns
    -------
    array
        The frequencies at which to evaluate the compressed waveform.
    """
    sample_points = []
    f = fmin
    while f < fmax:
        if df_multiple is not None:
            f = int(f/df_multiple)*df_multiple
        sample_points.append(f)
        f += 1.0 / rough_time_estimate(m1, m2, f, fudge_min=min_seglen)
    # add the last point
    if sample_points[-1] < fmax:
        sample_points.append(fmax)
    return numpy.array(sample_points)

def spa_compression(htilde, fmin, fmax, min_seglen=0.02,
        sample_frequencies=None):
    """Returns the frequencies needed to compress the given frequency domain
    waveform. This is done by estimating t(f) of the waveform using the
    stationary phase approximation.

    Parameters
    ----------
    htilde : FrequencySeries
        The waveform to compress.
    fmin : float
        The starting frequency of the compressed waveform.
    fmax : float
        The ending frequency of the compressed waveform.
    min_seglen : float
        The inverse of this gives the maximum frequency step that is used.
    sample_frequencies : {None, array}
        The frequencies that the waveform is evaluated at. If None, will
        retrieve the frequencies from the waveform's sample_frequencies
        attribute.

    Returns
    -------
    array
        The frequencies at which to evaluate the compressed waveform.
    """
    if sample_frequencies is None:
        sample_frequencies = htilde.sample_frequencies.numpy()
    kmin = int(fmin/htilde.delta_f)
    kmax = int(fmax/htilde.delta_f)
    tf = abs(utils.time_from_frequencyseries(htilde,
            sample_frequencies=sample_frequencies).data[kmin:kmax])
    sample_frequencies = sample_frequencies[kmin:kmax]
    sample_points = []
    f = fmin
    while f < fmax:
        f = int(f/htilde.delta_f)*htilde.delta_f
        sample_points.append(f)
        jj = numpy.searchsorted(sample_frequencies, f)
        f += 1./(tf[jj:].max()+min_seglen)
    # add the last point
    if sample_points[-1] < fmax:
        sample_points.append(fmax)
    return numpy.array(sample_points)

compression_algorithms = {
        'mchirp': mchirp_compression,
        'spa': spa_compression
        }


def fd_decompress(amp, phase, sample_frequencies, out=None, df=None,
        f_lower=None, interpolation='linear'):
    """Decompresses an FD waveform using the given amplitude, phase, and the
    frequencies at which they are sampled at.

    Parameters
    ----------
    amp : array
        The amplitude of the waveform at the sample frequencies.
    phase : array
        The phase of the waveform at the sample frequencies.
    sample_frequencies : array
        The frequency (in Hz) of the waveform at the sample frequencies.
    out : {None, FrequencySeries}
        The output array to save the decompressed waveform to. If this contains
        slots for frequencies > the maximum frequency in sample_frequencies,
        the rest of the values are zeroed. If not provided, must provide a df.
    df : {None, float}
        The frequency step to use for the decompressed waveform. Must be
        provided if out is None.
    f_lower : {None, float}
        The frequency to start the decompression at. If None, will use whatever
        the lowest frequency is in sample_frequencies. All values at
        frequencies less than this will be 0 in the decompressed waveform.
    interpolation : {'linear', str}
        The interpolation to use for the amplitude and phase. Default is
        'linear'. If 'linear' a custom interpolater is used. Otherwise,
        ``scipy.interpolate.interp1d`` is used; for other options, see
        possible values for that function's ``kind`` argument.

    Returns
    -------
    out : FrqeuencySeries
        If out was provided, writes to that array. Otherwise, a new
        FrequencySeries with the decompressed waveform.
    """
        
    if out is None:
        if df is None:
            raise ValueError("Either provide output memory or a df")
        flen = int(numpy.ceil(sample_frequencies.max()/df+1))
        out = FrequencySeries(numpy.zeros(flen,
            dtype=numpy.complex128), copy=False, delta_f=df)
    else:
        df = out.delta_f
        flen = len(out)
    if f_lower is None:
        jmin = 0
        f_lower = sample_frequencies[0]
    else:
        if f_lower >= sample_frequencies.max():
            raise ValueError("f_lower is > than the maximum sample frequency")
        jmin = int(numpy.searchsorted(sample_frequencies, f_lower))
    imin = int(numpy.floor(f_lower/df))
    # interpolate the amplitude and the phase
    if interpolation == "linear":
        # use custom interpolation
        sflen = len(sample_frequencies)
        h = numpy.array(out.data, copy=False)
        # make sure df is a float
        df = float(df)
        code = r"""
        # include <math.h>
        # include <stdio.h>
        int j = jmin-1;
        double sf = 0.;
        double A = 0.;
        double nextA = 0.;
        double phi = 0.;
        double nextPhi = 0.;
        double next_sf = sample_frequencies[jmin];
        double f = 0.;
        double invsdf = 0.;
        double mAmp = 0.;
        double bAmp = 0.;
        double mPhi = 0.;
        double bPhi = 0.;
        double interpAmp = 0.;
        double interpPhi = 0.;
        // zero-out beginning of array
        std::fill(h, h+imin, std::complex<double>(0., 0.));
        // cycle over desired samples
        for (int i=imin; i<flen; i++){
            f = i*df;
            if (f >= next_sf){
                // update linear interpolations
                j += 1;
                // if we have gone beyond the sampled frequencies, just break
                if ((j+1) == sflen) {
                    // zero-out rest the rest of the array & exit
                    std::fill(h+i, h+flen, std::complex<double>(0., 0.));
                    break;
                }
                sf = (double) sample_frequencies[j];
                next_sf = (double) sample_frequencies[j+1];
                A = (double) amp[j];
                nextA = (double) amp[j+1];
                phi = (double) phase[j];
                nextPhi = (double) phase[j+1];
                invsdf = 1./(next_sf - sf);
                mAmp = (nextA - A)*invsdf;
                bAmp = A - mAmp*sf;
                mPhi = (nextPhi - phi)*invsdf;
                bPhi = phi - mPhi*sf;
            }
            interpAmp = mAmp * f + bAmp;
            interpPhi = mPhi * f + bPhi;
            h[i] = std::complex<double> (interpAmp*cos(interpPhi),
                                         interpAmp*sin(interpPhi));
        }
        """
        inline(code, ['flen', 'sflen', 'df', 'sample_frequencies',
                      'amp', 'phase', 'h', 'imin', 'jmin'],
               extra_compile_args=[WEAVE_FLAGS + '-march=native -O3 -w'] +\
                                  omp_flags,
               libraries=omp_libs)
    else:
        # use scipy for fancier interpolation
        outfreq = out.sample_frequencies.numpy()
        amp_interp = interpolate.interp1d(sample_frequencies, amp,
            kind=interpolation, bounds_error=False, fill_value=0.,
            assume_sorted=True)
        phase_interp = interpolate.interp1d(sample_frequencies, phase,
            kind=interpolation, bounds_error=False, fill_value=0.,
            assume_sorted=True)
        A = amp_interp(outfreq)
        phi = phase_interp(outfreq)
        out.data[:] = A*numpy.cos(phi) + (1j)*A*numpy.sin(phi)
    return out
