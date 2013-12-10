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
import logging
from pycbc.types import zeros, real_same_precision_as, TimeSeries, complex_same_precision_as
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
                bank_template, stilde, psd,
                low_frequency_cutoff=low_frequency_cutoff)
        # SNR time series stored here
        snrs.append(snr)
        # Template normalization factor stored here
        norms.append(norm)
        
    return snrs, norms

def template_overlaps(bank_filters, template, template_sigmasq, psd, low_frequency_cutoff):
    """ This functions calculates the overlaps between the template and the
    bank veto templates.
    
    Parameters
    ----------
    bank_filters: List of FrequencySeries
    template: FrequencySeries
    template_sigmasq: float
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
        norm = sqrt(1 / template_sigmasq / bank_template.sigmasq)
        overlaps.append(overlap * norm)
    return overlaps

def bank_chisq_from_filters(tmplt_snr, tmplt_norm, bank_snrs, bank_norms,
        tmplt_bank_matchs, indices=None):
    """ This function calculates and returns a TimeSeries object containing the
    bank veto calcuated over a segment.
    
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
    tmplt_bank_matchs: list of floats
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
        tmplt_snr = tmplt_snr.take(indices)
        bank_snrs_tmp = []
        for bank_snr in bank_snrs:
            bank_snrs_tmp.append(bank_snr.take(indices))
        bank_snrs=bank_snrs_tmp
    
    # Initialise bank_chisq as 0s everywhere
    bank_chisq = zeros(len(tmplt_snr), dtype=real_same_precision_as(tmplt_snr))

    # Loop over all the bank templates
    for i in range(len(bank_snrs)):
        bank_match = tmplt_bank_matchs[i]
        
        bank_norm = sqrt((1 - bank_match*bank_match.conj()).real)
        
        bank_SNR = bank_snrs[i] * (bank_norms[i] / bank_norm)     
        tmplt_SNR = tmplt_snr * (bank_match.conj() * tmplt_norm / bank_norm)
        
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
    def __init__(self, bank_file, approximant, psd, segments, f_low, **kwds):
        if bank_file is not None:
            print bank_file
            self.do = True
            
            self.column_name = "bank_chisq"
            self.table_dof_name = "bank_chisq_dof"

            self.psd = psd    
            self.cdtype = complex_same_precision_as(psd)
            self.delta_f = psd.delta_f
            self.f_low = f_low
            self.seg_len_freq = len(psd)
            self.seg_len_time = (self.seg_len_freq-1)*2

            logging.info("Read in bank veto template bank")
            bank_veto_bank = FilterBank(bank_file,
                    approximant, self.seg_len_freq, 
                    self.delta_f, f_low, dtype=self.cdtype, psd=self.psd, **kwds)

            self.filters = list(bank_veto_bank)

            logging.info("Precalculate the bank veto template snrs")
            self.snr_data = []
            for seg in segments:
                self.snr_data.append(segment_snrs(self.filters, seg, psd, f_low))
                  
            self.dof = len(bank_veto_bank) * 2 - 2

            self._overlaps = None
            self._template = None
        else:
            self.do = False
        
    def values(self, template, s_num, snr, norm, indices):
        if self.do:
            logging.info("...Doing Bank Chisq")
            
            #Only calculate the overlaps if I haven't already and the template hasn't changed
            if self._overlaps is None or self._template != template:
                logging.info("...Calculate Bank Chisq Overlaps")
                self._overlaps = template_overlaps(self.filters, template, template.sigmasq, self.psd, self.f_low)
                self._template = template
                         
            bank_veto_snrs, bank_veto_norms = self.snr_data[s_num]
            return bank_chisq_from_filters(snr, norm, bank_veto_snrs, bank_veto_norms, self._overlaps, indices)

        
