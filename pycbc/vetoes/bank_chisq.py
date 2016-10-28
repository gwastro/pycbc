# Copyright (C) 2013  Ian Harry, Alex Nitz
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
import logging, numpy
from pycbc.types import Array, zeros, real_same_precision_as, TimeSeries
from pycbc.filter import overlap_cplx, matched_filter_core
from pycbc.waveform import FilterBank
from math import sqrt

def segment_snrs(filters, stilde, psd, low_frequency_cutoff):
    """ This functions calculates the snr of each bank veto template against
    the segment

    Parameters
    ----------
    filters: list of FrequencySeries
        The list of bank veto templates filters.
    stilde: FrequencySeries
        The current segment of data.
    psd: FrequencySeries
    low_frequency_cutoff: float

    Returns
    -------
    snr (list): List of snr time series.
    norm (list): List of normalizations factors for the snr time series.
    """
    snrs = []
    norms = []

    for i, bank_template in enumerate(filters):
        # For every template compute the snr against the stilde segment
        snr, corr, norm = matched_filter_core(
                bank_template, stilde, h_norm=bank_template.sigmasq(psd),
                psd=None, low_frequency_cutoff=low_frequency_cutoff)
        # SNR time series stored here
        snrs.append(snr)
        # Template normalization factor stored here
        norms.append(norm)

    return snrs, norms

def template_overlaps(bank_filters, template, psd, low_frequency_cutoff):
    """ This functions calculates the overlaps between the template and the
    bank veto templates.

    Parameters
    ----------
    bank_filters: List of FrequencySeries
    template: FrequencySeries
    psd: FrequencySeries
    low_frequency_cutoff: float

    Returns
    -------
    overlaps: List of complex overlap values.
    """
    overlaps = []
    template_ow = template / psd
    for bank_template in bank_filters:
        overlap = overlap_cplx(template_ow, bank_template,
                low_frequency_cutoff=low_frequency_cutoff, normalized=False)
        norm = sqrt(1 / template.sigmasq(psd) / bank_template.sigmasq(psd))
        overlaps.append(overlap * norm)
        if (abs(overlaps[-1]) > 0.99):
            errMsg = "Overlap > 0.99 between bank template and filter. "
            errMsg += "This bank template will not be used to calculate "
            errMsg += "bank chisq for this filter template. The expected "
            errMsg += "value will be added to the chisq to account for "
            errMsg += "the removal of this template.\n"
            errMsg += "Masses of filter template: %e %e\n" \
                      %(template.params.mass1, template.params.mass2)
            errMsg += "Masses of bank filter template: %e %e\n" \
                      %(bank_template.params.mass1, bank_template.params.mass2)
            errMsg += "Overlap: %e" %(abs(overlaps[-1]))
            logging.debug(errMsg)
    return overlaps

def bank_chisq_from_filters(tmplt_snr, tmplt_norm, bank_snrs, bank_norms,
        tmplt_bank_matches, indices=None):
    """ This function calculates and returns a TimeSeries object containing the
    bank veto calculated over a segment.

    Parameters
    ----------
    tmplt_snr: TimeSeries
        The SNR time series from filtering the segment against the current
        search template
    tmplt_norm: float
        The normalization factor for the search template
    bank_snrs: list of TimeSeries
        The precomputed list of SNR time series between each of the bank veto
        templates and the segment
    bank_norms: list of floats
        The normalization factors for the list of bank veto templates
        (usually this will be the same for all bank veto templates)
    tmplt_bank_matches: list of floats
        The complex overlap between the search template and each
        of the bank templates
    indices: {None, Array}, optional
        Array of indices into the snr time series. If given, the bank chisq
        will only be calculated at these values.

    Returns
    -------
    bank_chisq: TimeSeries of the bank vetos
    """
    if indices is not None:
        tmplt_snr = Array(tmplt_snr, copy=False)
        bank_snrs_tmp = []
        for bank_snr in bank_snrs:
            bank_snrs_tmp.append(bank_snr.take(indices))
        bank_snrs=bank_snrs_tmp

    # Initialise bank_chisq as 0s everywhere
    bank_chisq = zeros(len(tmplt_snr), dtype=real_same_precision_as(tmplt_snr))

    # Loop over all the bank templates
    for i in range(len(bank_snrs)):
        bank_match = tmplt_bank_matches[i]
        if (abs(bank_match) > 0.99):
            # Not much point calculating bank_chisquared if the bank template
            # is very close to the filter template. Can also hit numerical
            # error due to approximations made in this calculation.
            # The value of 2 is the expected addition to the chisq for this
            # template
            bank_chisq += 2.
            continue
        bank_norm = sqrt((1 - bank_match*bank_match.conj()).real)

        bank_SNR = bank_snrs[i] * (bank_norms[i] / bank_norm)
        tmplt_SNR = tmplt_snr * (bank_match.conj() * tmplt_norm / bank_norm)

        bank_SNR = Array(bank_SNR, copy=False)
        tmplt_SNR = Array(tmplt_SNR, copy=False)

        bank_chisq += (bank_SNR - tmplt_SNR).squared_norm()

    if indices is not None:
        return bank_chisq
    else:
        return TimeSeries(bank_chisq, delta_t=tmplt_snr.delta_t,
                          epoch=tmplt_snr.start_time, copy=False)

class SingleDetBankVeto(object):
    """This class reads in a template bank file for a bank veto, handles the
       memory management of its filters internally, and calculates the bank
       veto TimeSeries.
    """
    def __init__(self, bank_file, flen, delta_f, f_low, cdtype, approximant=None, **kwds):
        if bank_file is not None:
            self.do = True

            self.column_name = "bank_chisq"
            self.table_dof_name = "bank_chisq_dof"

            self.cdtype = cdtype
            self.delta_f = delta_f
            self.f_low = f_low
            self.seg_len_freq = flen
            self.seg_len_time = (self.seg_len_freq-1)*2

            logging.info("Read in bank veto template bank")
            bank_veto_bank = FilterBank(bank_file,
                    self.seg_len_freq,
                    self.delta_f, self.cdtype,
                    low_frequency_cutoff=f_low,
                    approximant=approximant, **kwds)

            self.filters = list(bank_veto_bank)
            self.dof = len(bank_veto_bank) * 2

            self._overlaps_cache = {}
            self._segment_snrs_cache = {}
        else:
            self.do = False

    def cache_segment_snrs(self, stilde, psd):
        key = (id(stilde), id(psd))
        if key not in self._segment_snrs_cache:
            logging.info("Precalculate the bank veto template snrs")
            data = segment_snrs(self.filters, stilde, psd, self.f_low)
            self._segment_snrs_cache[key] = data
        return self._segment_snrs_cache[key]

    def cache_overlaps(self, template, psd):
        key = (id(template.params), id(psd))
        if key not in self._overlaps_cache:
            logging.info("...Calculate bank veto overlaps")
            o = template_overlaps(self.filters, template, psd, self.f_low)
            self._overlaps_cache[key] = o
        return self._overlaps_cache[key]

    def values(self, template, psd, stilde, snrv, norm, indices):
        """
        Returns
        -------
        bank_chisq_from_filters: TimeSeries of bank veto values - if indices
        is None then evaluated at all time samples, if not then only at 
        requested sample indices

        bank_chisq_dof: int, approx number of statistical degrees of freedom
        """
        if self.do:
            logging.info("...Doing bank veto")
            overlaps = self.cache_overlaps(template, psd)
            bank_veto_snrs, bank_veto_norms = self.cache_segment_snrs(stilde, psd)
            chisq = bank_chisq_from_filters(snrv, norm, bank_veto_snrs,
                                            bank_veto_norms, overlaps, indices) 
            dof = numpy.repeat(self.dof, len(chisq))
            return chisq, dof
        else:
            return None, None      

class SingleDetSkyMaxBankVeto(SingleDetBankVeto):
    """Stub for precessing bank veto if anyone ever wants to code it up.
    """
    def __init__(self, *args, **kwds):
        super(SingleDetSkyMaxBankVeto, self).__init__(*args, **kwds)

    def values(self, *args, **kwargs):
        if self.do:
            err_msg = "Precessing single detector sky-max bank veto has not "
            err_msg += "been written. If you want to use it, why not help "
            err_msg += "write it?"
            raise NotImplementedError(err_msg)
        else:
            return None, None
