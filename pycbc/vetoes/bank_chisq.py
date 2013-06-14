# Copyright (C) 2012  Ian Harry
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

from pycbc.types import zeros, real_same_precision_as, TimeSeries
from pycbc.filter import overlap_cplx, matched_filter_core
from pycbc.waveform import TemplateBank

def bank_chisq_from_filters(tmplt_snr, tmplt_norm, bank_snrs, bank_norms,\
        tmplt_bank_matchs):
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

    Returns
    -------
    The function returns the bank_veto TimeSeries object.
    """
    
    # Initialise bank_chisq as 0s everywhere
    bank_chisq = TimeSeries(zeros(len(tmplt_snr),\
            dtype=real_same_precision_as(tmplt_snr)), delta_t=tmplt_snr.delta_t,\
            epoch=tmplt_snr.start_time, copy=False)

    # Get normalized real and imaginary components of the SNR as numpy arrays
    tmplt_SNR = tmplt_snr * tmplt_norm

    # Loop over all the bank templates
    for i in range(len(bank_snrs)):
        # Get normalized real and imaginary components of the bank template SNRs
        bank_SNR = bank_snrs[i] * bank_norms[i]       
        # Get the overlap between the search template and the bank veto tmplt
        bank_match = tmplt_bank_matchs[i]  
        # Calculate components of bank veto
        numerator = bank_SNR - tmplt_SNR * bank_match.conj()
        bank_chisq_tmp = numerator.squared_norm()
        bank_chisq += bank_chisq_tmp / \
                (1 - bank_match*bank_match.conjugate()).real

    return bank_chisq
    
class BankVeto(object):
    """This class reads in a template bank file for a bank veto, handles the
       memory management of its filters internally, and calculates the bank
       veto TimeSeries.
    """
    def __init__(self, bank_file, approximant, seg_len_freq, delta_f, f_low, dtype, **kwds):
        self.filters = []
        self.snrs_work_mem = []
        
        self.dtype = dtype
        self.delta_f = delta_f
        self.f_low = f_low
        self.seg_len_freq = seg_len_freq
        self.seg_len_time = (seg_len_freq-1)*2
    
        # Read in the bank veto bank
        bank_veto_bank = TemplateBank(bank_file,
                approximant, self.seg_len_freq, 
                delta_f, f_low, dtype=dtype, **kwds)
                
        # The following command actually generates all the filters
        self.filters = list(bank_veto_bank)
        
        # Create memory for storing the bank veto filter time series
        for i in range(len(self)):
            self.snrs_work_mem.append(zeros(self.seg_len_time, dtype=dtype))
        
    def __len__(self):
        return len(self.filters)
        
    def segment_snrs(self, stilde, psd=None):
        """ Return the snrs for each bank filter against stilde.
        """
        snrs = []
        norms = []
        
        for i, bank_template in enumerate(self.filters):
            # For every template compute the snr against the stilde segment
            snr, corr, norm = matched_filter_core(
                    bank_template, stilde, psd,
                    low_frequency_cutoff=self.f_low,
                    out=self.snrs_work_mem[i])
            # SNR time series stored here
            snrs.append(snr)
            # Template normalization factor stored here
            norms.append(norm)
            
        return snrs, norms
        
    def template_overlaps(self, template, psd=None):
        overlaps = []
        for bank_template in self.filters:
            overlap = overlap_cplx(template, bank_template, psd=psd,
                    low_frequency_cutoff=self.f_low)
            overlaps.append(overlap)
        return overlaps

