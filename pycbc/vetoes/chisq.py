# Copyright (C) 2012  Alex Nitz
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
import numpy, logging, math, pycbc.fft

from pycbc.types import zeros, real_same_precision_as, TimeSeries, complex_same_precision_as
from pycbc.filter import sigmasq_series, make_frequency_series, matched_filter_core, get_cutoff_indices
from pycbc.scheme import schemed
import pycbc.pnutils

BACKEND_PREFIX="pycbc.vetoes.chisq_"

def power_chisq_bins_from_sigmasq_series(sigmasq_series, num_bins, kmin, kmax):
    """Returns bins of equal power for use with the chisq functions

    Parameters
    ----------

    sigmasq_series: FrequencySeries
        A frequency series containing the cumulative power of a filter template
        preweighted by a psd.
    num_bins: int
        The number of chisq bins to calculate.
    kmin: int
        DOCUMENTME
    kmax: int
        DOCUMENTME

    Returns
    -------

    bins: List of ints
        A list of the edges of the chisq bins is returned.

    """
    sigmasq = sigmasq_series[kmax - 1]
    edge_vec = numpy.arange(0, num_bins) * sigmasq / num_bins
    bins = numpy.searchsorted(sigmasq_series[kmin:kmax], edge_vec, side='right')
    bins += kmin
    return numpy.append(bins, kmax)

def power_chisq_bins(htilde, num_bins, psd, low_frequency_cutoff=None,
                     high_frequency_cutoff=None):
    """Returns bins of equal power for use with the chisq functions

    Parameters
    ----------

    htilde: FrequencySeries
        A frequency series containing the template waveform
    num_bins: int
        The number of chisq bins to calculate.
    psd: FrequencySeries
        A frequency series containing the psd. Its length must be commensurate
        with the template waveform.
    low_frequency_cutoff: {None, float}, optional
        The low frequency cutoff to apply
    high_frequency_cutoff: {None, float}, optional
        The high frequency cutoff to apply

    Returns
    -------

    bins: List of ints
        A list of the edges of the chisq bins is returned.
    """
    sigma_vec = sigmasq_series(htilde, psd, low_frequency_cutoff,
                               high_frequency_cutoff).numpy()
    kmin, kmax = get_cutoff_indices(low_frequency_cutoff,
                                    high_frequency_cutoff,
                                    htilde.delta_f,
                                    (len(htilde)-1)*2)
    return power_chisq_bins_from_sigmasq_series(sigma_vec, num_bins, kmin, kmax)

@schemed(BACKEND_PREFIX)
def chisq_accum_bin(chisq, q):
    pass

@schemed(BACKEND_PREFIX)
def shift_sum(v1, shifts, bins):
    """ Calculate the time shifted sum of the FrequencySeries
    """
    pass

def power_chisq_at_points_from_precomputed(corr, snr, snr_norm, bins, indices):
    """Calculate the chisq timeseries from precomputed values for only select points.

    This function calculates the chisq at each point by explicitly time shifting
    and summing each bin. No FFT is involved.

    Parameters
    ----------
    corr: FrequencySeries
        The product of the template and data in the frequency domain.
    snr: numpy.ndarray
        The unnormalized array of snr values at only the selected points in `indices`.
    snr_norm: float
        The normalization of the snr (EXPLAINME : refer to Findchirp paper?)
    bins: List of integers
        The edges of the equal power bins
    indices: Array
        The indices where we will calculate the chisq. These must be relative
        to the given `corr` series.

    Returns
    -------
    chisq: Array
        An array containing only the chisq at the selected points.
    """
    logging.info('doing fast point chisq')
    num_bins = len(bins) - 1
    chisq = shift_sum(corr, indices, bins) # pylint:disable=assignment-from-no-return
    return (chisq * num_bins - (snr.conj() * snr).real) * (snr_norm ** 2.0)

_q_l = None
_qtilde_l = None
_chisq_l = None
def power_chisq_from_precomputed(corr, snr, snr_norm, bins, indices=None, return_bins=False):
    """Calculate the chisq timeseries from precomputed values.

    This function calculates the chisq at all times by performing an
    inverse FFT of each bin.

    Parameters
    ----------

    corr: FrequencySeries
        The produce of the template and data in the frequency domain.
    snr: TimeSeries
        The unnormalized snr time series.
    snr_norm:
        The snr normalization factor (true snr = snr * snr_norm) EXPLAINME - define 'true snr'?
    bins: List of integers
        The edges of the chisq bins.
    indices: {Array, None}, optional
        Index values into snr that indicate where to calculate
        chisq values. If none, calculate chisq for all possible indices.
    return_bins: {boolean, False}, optional
        Return a list of the SNRs for each chisq bin.

    Returns
    -------
    chisq: TimeSeries
    """
    # Get workspace memory
    global _q_l, _qtilde_l, _chisq_l
    
    bin_snrs = []

    if _q_l is None or len(_q_l) != len(snr):
        q = zeros(len(snr), dtype=complex_same_precision_as(snr))
        qtilde = zeros(len(snr), dtype=complex_same_precision_as(snr))
        _q_l = q
        _qtilde_l = qtilde
    else:
        q = _q_l
        qtilde = _qtilde_l

    if indices is not None:
        snr = snr.take(indices)

    if _chisq_l is None or len(_chisq_l) < len(snr):
        chisq = zeros(len(snr), dtype=real_same_precision_as(snr))
        _chisq_l = chisq
    else:
        chisq = _chisq_l[0:len(snr)]
        chisq.clear()

    num_bins = len(bins) - 1

    for j in range(num_bins):
        k_min = int(bins[j])
        k_max = int(bins[j+1])

        qtilde[k_min:k_max] = corr[k_min:k_max]
        pycbc.fft.ifft(qtilde, q)
        qtilde[k_min:k_max].clear()

        if return_bins:
            bin_snrs.append(TimeSeries(q  * snr_norm *  num_bins ** 0.5,
                                      delta_t=snr.delta_t,
                                      epoch=snr.start_time))

        if indices is not None:
            chisq_accum_bin(chisq, q.take(indices))
        else:
            chisq_accum_bin(chisq, q)

    chisq = (chisq * num_bins - snr.squared_norm()) * (snr_norm ** 2.0)

    if indices is None:
        chisq = TimeSeries(chisq, delta_t=snr.delta_t, epoch=snr.start_time, copy=False)

    if return_bins:
        return chisq, bin_snrs
    else:
        return chisq

def fastest_power_chisq_at_points(corr, snr, snrv, snr_norm, bins, indices):
    """Calculate the chisq values for only selected points.

    This function looks at the number of points to be evaluated and selects
    the fastest method (FFT, or direct time shift and sum). In either case,
    only the selected points are returned.

    Parameters
    ----------
    corr: FrequencySeries
        The product of the template and data in the frequency domain.
    snr: Array
        The unnormalized snr
    snr_norm: float
        The snr normalization factor  --- EXPLAINME
    bins: List of integers
        The edges of the equal power bins
    indices: Array
        The indices where we will calculate the chisq. These must be relative
        to the given `snr` series.

    Returns
    -------
    chisq: Array
        An array containing only the chisq at the selected points.
    """
    import pycbc.scheme
    if isinstance(pycbc.scheme.mgr.state, pycbc.scheme.CPUScheme):
        # We don't have that many points so do the direct time shift.
        return power_chisq_at_points_from_precomputed(corr, snrv,
                                                      snr_norm, bins, indices)
    else:
        # We have a lot of points so it is faster to use the fourier transform
        return power_chisq_from_precomputed(corr, snr, snr_norm, bins,
                                            indices=indices)


def power_chisq(template, data, num_bins, psd,
                low_frequency_cutoff=None,
                high_frequency_cutoff=None,
                return_bins=False):
    """Calculate the chisq timeseries

    Parameters
    ----------
    template: FrequencySeries or TimeSeries
        A time or frequency series that contains the filter template.
    data: FrequencySeries or TimeSeries
        A time or frequency series that contains the data to filter. The length
        must be commensurate with the template.
        (EXPLAINME - does this mean 'the same as' or something else?)
    num_bins: int
        The number of bins in the chisq. Note that the dof goes as 2*num_bins-2.
    psd: FrequencySeries
        The psd of the data.
    low_frequency_cutoff: {None, float}, optional
        The low frequency cutoff for the filter
    high_frequency_cutoff: {None, float}, optional
        The high frequency cutoff for the filter
    return_bins: {boolean, False}, optional
        Return a list of the individual chisq bins

    Returns
    -------
    chisq: TimeSeries
        TimeSeries containing the chisq values for all times.
    """
    htilde = make_frequency_series(template)
    stilde = make_frequency_series(data)

    bins = power_chisq_bins(htilde, num_bins, psd, low_frequency_cutoff,
                            high_frequency_cutoff)
    corra = zeros((len(htilde)-1)*2, dtype=htilde.dtype)
    total_snr, corr, tnorm = matched_filter_core(htilde, stilde, psd,
                           low_frequency_cutoff, high_frequency_cutoff,
                           corr_out=corra)

    return power_chisq_from_precomputed(corr, total_snr, tnorm, bins, return_bins=return_bins)


class SingleDetPowerChisq(object):
    """Class that handles precomputation and memory management for efficiently
    running the power chisq in a single detector inspiral analysis.
    """
    def __init__(self, num_bins=0, snr_threshold=None):
        if not (num_bins == "0" or num_bins == 0):
            self.do = True
            self.column_name = "chisq"
            self.table_dof_name = "chisq_dof"
            self.num_bins = num_bins
        else:
            self.do = False
        self.snr_threshold = snr_threshold

    @staticmethod
    def parse_option(row, arg):
        safe_dict = {}
        safe_dict.update(row.__dict__)
        safe_dict.update(math.__dict__)
        safe_dict.update(pycbc.pnutils.__dict__)
        return eval(arg, {"__builtins__":None}, safe_dict)

    def cached_chisq_bins(self, template, psd):
        key = id(psd)
        if not hasattr(psd, '_chisq_cached_key'):
            psd._chisq_cached_key = {}

        if not hasattr(template, '_bin_cache'):
            template._bin_cache = {}

        if key not in template._bin_cache or id(template.params) not in psd._chisq_cached_key:
            psd._chisq_cached_key[id(template.params)] = True
            num_bins = int(self.parse_option(template, self.num_bins))

            if hasattr(psd, 'sigmasq_vec') and template.approximant in psd.sigmasq_vec:
                logging.info("...Calculating fast power chisq bins")
                kmin = int(template.f_lower / psd.delta_f)
                kmax = template.end_idx
                bins = power_chisq_bins_from_sigmasq_series(
                    psd.sigmasq_vec[template.approximant], num_bins, kmin, kmax)
            else:
                logging.info("...Calculating power chisq bins")
                bins = power_chisq_bins(template, num_bins, psd, template.f_lower)
            template._bin_cache[key] = bins

        return template._bin_cache[key]

    def values(self, corr, snrv, snr_norm, psd, indices, template):
        """ Calculate the chisq at points given by indices.

        Returns
        -------
        chisq: Array
            Chisq values, one for each sample index

        chisq_dof: Array
            Number of statistical degrees of freedom for the chisq test
            in the given template
        """

        if self.do:
            logging.info("...Doing power chisq")

            num_above = len(indices)
            if self.snr_threshold:
                above = abs(snrv * snr_norm) > self.snr_threshold
                num_above = above.sum()
                logging.info('%s above chisq activation threshold' % num_above)
                above_indices = indices[above]
                above_snrv = snrv[above]
                rchisq = numpy.zeros(len(indices), dtype=numpy.float32)
                dof = -100
            else:
                above_indices = indices
                above_snrv = snrv

            if num_above > 0:
                bins = self.cached_chisq_bins(template, psd)
                dof = (len(bins) - 1) * 2 - 2
                chisq = power_chisq_at_points_from_precomputed(corr,
                                     above_snrv, snr_norm, bins, above_indices)

            if self.snr_threshold:
                if num_above > 0:
                    rchisq[above] = chisq
            else:
                rchisq = chisq

            return rchisq, numpy.repeat(dof, len(indices))# dof * numpy.ones_like(indices)
        else:
            return None, None

class SingleDetSkyMaxPowerChisq(SingleDetPowerChisq):
    """Class that handles precomputation and memory management for efficiently
    running the power chisq in a single detector inspiral analysis when
    maximizing analytically over sky location.
    """
    def __init__(self, **kwds):
        super(SingleDetSkyMaxPowerChisq, self).__init__(**kwds)
        self.template_mem = None
        self.corr_mem = None

    def calculate_chisq_bins(self, template, psd):
        """ Obtain the chisq bins for this template and PSD.
        """
        num_bins = int(self.parse_option(template, self.num_bins))
        if hasattr(psd, 'sigmasq_vec') and \
                                       template.approximant in psd.sigmasq_vec:
            logging.info("...Calculating fast power chisq bins")
            kmin = int(template.f_lower / psd.delta_f)
            kmax = template.end_idx
            bins = power_chisq_bins_from_sigmasq_series(
                   psd.sigmasq_vec[template.approximant], num_bins, kmin, kmax)
        else:
            logging.info("...Calculating power chisq bins")
            bins = power_chisq_bins(template, num_bins, psd, template.f_lower)
        return bins


    def values(self, corr_plus, corr_cross, snrv, psd,
               indices, template_plus, template_cross, u_vals,
               hplus_cross_corr, hpnorm, hcnorm):
        """ Calculate the chisq at points given by indices.

        Returns
        -------
        chisq: Array
            Chisq values, one for each sample index

        chisq_dof: Array
            Number of statistical degrees of freedom for the chisq test
            in the given template
        """
        if self.do:
            logging.info("...Doing power chisq")

            num_above = len(indices)
            if self.snr_threshold:
                above = abs(snrv) > self.snr_threshold
                num_above = above.sum()
                logging.info('%s above chisq activation threshold' % num_above)
                above_indices = indices[above]
                above_snrv = snrv[above]
                rchisq = numpy.zeros(len(indices), dtype=numpy.float32)
                dof = -100
            else:
                above_indices = indices
                above_snrv = snrv

            if num_above > 0:
                chisq = []
                curr_tmplt_mult_fac = 0.
                curr_corr_mult_fac = 0.
                if self.template_mem is None or \
                        (not len(self.template_mem) == len(template_plus)):
                    self.template_mem = zeros(len(template_plus),
                                dtype=complex_same_precision_as(corr_plus))
                if self.corr_mem is None or \
                                (not len(self.corr_mem) == len(corr_plus)):
                    self.corr_mem = zeros(len(corr_plus),
                                dtype=complex_same_precision_as(corr_plus))

                tmplt_data = template_cross.data
                corr_data = corr_cross.data
                numpy.copyto(self.template_mem.data, template_cross.data)
                numpy.copyto(self.corr_mem.data, corr_cross.data)
                template_cross._data = self.template_mem.data
                corr_cross._data = self.corr_mem.data

                for lidx, index in enumerate(above_indices):
                    above_local_indices = numpy.array([index])
                    above_local_snr = numpy.array([above_snrv[lidx]])
                    local_u_val = u_vals[lidx]
                    # Construct template from _plus and _cross
                    # Note that this modifies in place, so we store that and
                    # revert on the next pass.
                    template = template_cross.multiply_and_add(template_plus,
                                               local_u_val-curr_tmplt_mult_fac)
                    curr_tmplt_mult_fac = local_u_val

                    template.f_lower = template_plus.f_lower
                    template.params = template_plus.params
                    # Construct the corr vector
                    norm_fac = local_u_val*local_u_val + 1
                    norm_fac += 2 * local_u_val * hplus_cross_corr
                    norm_fac = hcnorm / (norm_fac**0.5)
                    hp_fac = local_u_val * hpnorm / hcnorm
                    corr = corr_cross.multiply_and_add(corr_plus,
                                                   hp_fac - curr_corr_mult_fac)
                    curr_corr_mult_fac = hp_fac

                    bins = self.calculate_chisq_bins(template, psd)
                    dof = (len(bins) - 1) * 2 - 2
                    curr_chisq = power_chisq_at_points_from_precomputed(corr,
                                          above_local_snr/ norm_fac, norm_fac,
                                          bins, above_local_indices)
                    chisq.append(curr_chisq[0])
                chisq = numpy.array(chisq)
                # Must reset corr and template to original values!
                template_cross._data = tmplt_data
                corr_cross._data = corr_data

            if self.snr_threshold:
                if num_above > 0:
                    rchisq[above] = chisq
            else:
                rchisq = chisq

            return rchisq, numpy.repeat(dof, len(indices))# dof * numpy.ones_like(indices)
        else:
            return None, None
