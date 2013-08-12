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
from pycbc.types import Array, zeros, real_same_precision_as, TimeSeries, complex_same_precision_as, FrequencySeries
from pycbc.filter import sigmasq_series, make_frequency_series, sigmasq, matched_filter_core
import numpy
from pycbc.scheme import schemed
import pycbc.fft

BACKEND_PREFIX="pycbc.vetoes.chisq_"

def power_chisq_bins_from_sigmasq_series(sigmasq_series, num_bins, kmin, kmax):
    """Returns bins of equal power for use with the chisq functions
    """
    sigmasq = sigmasq_series[kmax - 1]                        
    edge_vec = numpy.arange(0, num_bins) * sigmasq / num_bins
    bins = numpy.searchsorted(sigmasq_series[kmin:kmax], edge_vec, side='right')
    bins += kmin
    return numpy.append(bins, kmax)

def power_chisq_bins(htilde, num_bins, psd, low_frequency_cutoff=None, 
                     high_frequency_cutoff=None):
    """Returns bins of equal power for use with the chisq functions
    """
    sigma_vec = sigmasq_series(htilde, psd, low_frequency_cutoff, 
                               high_frequency_cutoff).numpy() 
    kmin = int(low_frequency_cutoff / htilde.delta_f)
    kmax = len(sigma_vec)
    return power_chisq_bins_from_sigmasq_series(sigma_vec, num_bins, kmin, kmax)
    
@schemed(BACKEND_PREFIX)
def chisq_accum_bin(chisq, q):
    pass

def shift_sum(v1, shifts):
    v1 = v1.data

    from scipy.weave import inline
    import numpy

    pre = "float out0i=0, out0r=0, p0i=0, p0r=1, vs0r = vs[0].real(), vs0i = vs[0].imag();"

    op = """
            out0r += vr * p0r - vi * p0i;
            out0i += vr * p0i + vi * p0r;
            t1 = p0r;
            t2 = p0i;
            p0r = t1 * vs0r - t2 * vs0i;
            p0i = t1 * vs0i + t2 * vs0r; 
    """

    outstr = "out[0] = std::complex<float> (out0r, out0i);"
        
    outline = """
        std::complex<float> v;
        float t1, t2;
        
        %s;
        
        for (int j=0; j<vlen; j++){
            v = v1[j];
            float vr = v.real();
            float vi = v.imag();
             
             %s;          
                             
        }
        
        %s;
    """

    vs = numpy.zeros(len(shifts), dtype=numpy.complex64)
    out = numpy.zeros(len(shifts), dtype=numpy.complex64)
    n = int(len(out))
    vlen = int(len(v1))
    
    pre_u = ""
    outstr_u =""
    op_u = ""
    
    for i in range(len(shifts)):
        vs[i] = numpy.exp(1.0 * numpy.pi * 2j * float(shifts[i])/len(v1) )
        pre_u += pre.replace('0', str(i))
        outstr_u += outstr.replace('0', str(i))
        op_u += op.replace('0', str(i))
        
     
    code = outline % (pre_u, op_u, outstr_u)   

    inline(code, ['v1', 'vs', 'out', 'n', 'vlen'], )
    return Array(out, dtype=numpy.complex64)
    
def power_chisq_at_points_from_precomputed(corr, snr, h_norm, bins, indices):
    """ Returns the chisq time series
    """
    snr = Array(snr, copy=False)
    
    chisq = zeros(len(indices), dtype=real_same_precision_as(corr))     
    num_bins = len(bins) - 1
    
    norm = 4 * corr.delta_f / numpy.sqrt(h_norm)
    
    for j in range(len(bins)-1):
        k_min = int(bins[j]) - bins[0]
        k_max = int(bins[j+1]) - bins[0]           
        qi = shift_sum(corr[k_min:k_max], indices) * norm
        chisq += qi.squared_norm()

    return (chisq * num_bins - snr.squared_norm())

    
def power_chisq_from_precomputed(corr, snr, bins, snr_norm):
    """ Returns the chisq time series
    """
    q = zeros(len(snr), dtype=complex_same_precision_as(snr))
    qtilde = zeros(len(snr), dtype=complex_same_precision_as(snr))
    chisq = TimeSeries(zeros(len(snr), dtype=real_same_precision_as(snr)), 
                       delta_t=snr.delta_t, epoch=snr.start_time, copy=False)
    
    chisq_norm = snr_norm ** 2.0
    num_bins = len(bins) - 1
    
    for j in range(len(bins)-1): 
        k_min = int(bins[j])
        k_max = int(bins[j+1])
        
        qtilde[k_min:k_max] = corr[k_min:k_max]
        pycbc.fft.ifft(qtilde, q) 
        qtilde[k_min:k_max].clear()
        chisq_accum_bin(chisq, q)
        
    return (chisq * num_bins - snr.squared_norm()) * chisq_norm

def power_chisq(template, data, num_bins, psd, low_frequency_cutoff=None, high_frequency_cutoff=None):
    """ Return the chisq time series.
    """  
    htilde = make_frequency_series(template)
    stilde = make_frequency_series(data)
    
    bins = power_chisq_bins(htilde, num_bins, psd, low_frequency_cutoff, high_frequency_cutoff)
    
    corra = zeros((len(htilde)-1)*2, dtype=htilde.dtype)
    
    total_snr, corr, tnorm = matched_filter_core(htilde, stilde, psd,
                           low_frequency_cutoff, high_frequency_cutoff, corr_out=corra)

    return power_chisq_from_precomputed(corr, total_snr, bins, tnorm)

    
                
        


