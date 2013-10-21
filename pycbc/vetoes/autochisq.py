# Copyright (C) 2013  Stas Babak
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


import sys, re
import pycbc
from pycbc.filter import resample_to_delta_t, highpass, make_frequency_series
from pycbc.filter import  matched_filter_core, sigmasq, sigmasq_series
from pycbc.types import Array, TimeSeries, FrequencySeries, float32, complex64, zeros
from pycbc.types import complex_same_precision_as,real_same_precision_as
import numpy as np
from math import cos, sin, sqrt, pi, atan2, exp 


BACKEND_PREFIX="pycbc.vetoes.autochisq_"


def autochisq_from_precomputed(sn, corr, hautocorr, stride=1, num_points=None, indices=None):

    """ Compute correlation (two sided) between template and data
        and compares with autocorrelation of the template: C(t) = IFFT(A*A/S(f))

	Parameters
	----------
	sn: Array[complex]
	    normalized (!) array of complex snr 
	corr: FrequencySeries[complex] 
	    correlation of the template and the data in freq. domain ### not used but could be
	hautocorr: Array[real] 
	     time domain autocorrelation for the template
	num_points: [int, optional] number of points used for autochisq
	             on each side, if None allpoints are used
	stride: [int, optional] stride for points selection for autochisq
	        total length <= 2*num_points*stride, default=1
	indices: Array[int], optional], compute autochisquare only at the points specified 
		 in this array, default: compute autochisquare for max snr in the snr array

	Returns
	-------
	autochisq: [tupple] returns autochisq values and snr corresponding to the instances 
	          of time defined by indices 

    """
 	
     
    snr_v = Array(np.abs(sn), copy=True)
    Nsnr = len(sn) 

    ## normalizing the hautocorr
    hautocorr = hautocorr/hautocorr[0]

    indx = np.array([])
    if (indices == None): ### find the maximum snr 
	maxSNR, max_ind = snr_v.max_loc()
	indx = np.append(indx, max_ind)
	indices = Array(indx, copy=True)
    
    
    achisq = np.zeros(len(indices))
    num_points_all = int(Nsnr/stride)
    if num_points is None:
       num_points = num_points_all
    if (num_points > num_points_all):   
       num_points = num_points_all
    
    ip = 0
    for ind in indices:
        end_point = ind + stride*num_points+1
        phi = atan2(sn[ind].imag, sn[ind].real)
        cphi = cos(phi)
        sphi = sin(phi)
        k = 1   
        snr_ind =  sn[ind].real*cphi + sn[ind].imag*sphi
	# going right
	for i in xrange(int(ind)+1, int(end_point), stride):
            if (i>Nsnr-1):
		  i = i - Nsnr  # folding it
	    z = sn[i].real*cphi + sn[i].imag*sphi
	    dz = z - hautocorr[k]*snr_ind
	    chisq_norm = 1.0 - hautocorr[k]*hautocorr[k]
	    achisq[ip] += dz*dz/chisq_norm
	    k += stride
        # going left
        end_point = ind - stride*num_points -1
        k=1
	for i in xrange(int(ind)-1, int(end_point), -stride):
            if (i<0):
		  i = i + Nsnr  # folding it
	    z = sn[i].real*cphi + sn[i].imag*sphi
	    dz = z - hautocorr[k]*snr_ind
	    chisq_norm = 1.0 - hautocorr[k]*hautocorr[k]
	    achisq[ip] += dz*dz/chisq_norm
	    k += stride
        # Stas last two cycles can be combined in a single cycle in "k" 
        ip += 1

    dof = 2*num_points	    
    achisq_list=np.zeros((len(indices), 3))

    achisq_list[:,0] = indices
    achisq_list[:,2] = achisq
    for i in xrange(len(indices)):
	achisq_list[i, 1] = snr_v[indices[i]]
 
    return(dof, achisq_list)


def autochisq(template, data, psd, stride=1, snr_thr=8.0, num_points=None, \
		low_frequency_cutoff=None, high_frequency_cutoff=None, max_snr=True):

    
    """ Compute correlation (two sided) between template and data
        and compares with autocorrelation of the template: C(t) = IFFT(A*A/S(f))

	Parameters
	----------
        template: FrequencySeries or TimeSeries
	    A time or frequency series that contains the filter template. The length
            must be commensurate with the data. 
        data: FrequencySeries or TimeSeries
            A time ore frequency series that contains the data to filter. The length
            must be commensurate with the template.
        psd: FrequencySeries
            The psd of the data. 
	num_points: [int, optional] number of points used for autochisq
	             on each side, if None allpoints are used
	stride: [int, optional] stride for points selection for autochisq
	        total length <= 2*num_points*stride, default=1
	low_frequency_cutoff: [float, None], 
	        lower frequency for filtering, if None, start at df
	high_frequency_cutoff: [float, None], higher frequency for filtering, 
	                      if None, stops at Nyquist frequency
	snr_thr: float, default=8.0
	    Compute autochisq for every point with snr >= snr_thr, if snr_thr<=0.0
	    then autochisq is not computed unless max_snr=True
	max_snr: bool, default=True
	     if True we compute autochisq for a point with maximum snr (independent of the 
	     snr_thr) otherwise autochisq is computed only if snr>=snr_thr>0

	Returns
	-------
	autochisq: [tupple] returns autochisq values and snr corresponding to the instances 
	          of time defined by indices 

    """
 
    ### {{{

    htilde = make_frequency_series(template)
    xtilde = make_frequency_series(data)

    if len(htilde) != len(xtilde):
        raise ValueError("Length of template and data must match")

    N = (len(xtilde)-1) * 2 

    _Ptilde = zeros(N, dtype=complex_same_precision_as(template))
    Pt = zeros(N, dtype=complex_same_precision_as(template))

    Pt, _Ptilde, P_norm = matched_filter_core(htilde, htilde, psd=psd, \
		   low_frequency_cutoff=low_frequency_cutoff,\
                  high_frequency_cutoff=high_frequency_cutoff)
    _P = Array(Pt.real(), copy=True)
    #_P = _P/_P[0]  # normalize it
 
    
    sn, cor, hnrm = matched_filter_core(htilde, xtilde, psd=psd, \
		   low_frequency_cutoff=low_frequency_cutoff,\
                  high_frequency_cutoff=high_frequency_cutoff)
 
    sn = sn*hnrm 
    Nsnr = len(sn) 
    index_list = np.array([])
    if (max_snr or snr_thr > 0.0):
        # compute snr data set
        snr_v = np.abs(sn)
        if (max_snr):
	    # looking for max snr
            ind_max = np.argmax(snr_v)
	    maxSNR = snr_v[ind_max]
            index_list = np.append(index_list, int(ind_max))
        if (maxSNR >= snr_thr and snr_thr>0.0):
	    for i in xrange(Nsnr):
		if (snr_v[i] >= snr_thr and snr_v[i] != maxSNR):
			index_list = np.append(index_list, int(i))
    #print "Stas, size of index list = ", len(index_list)

    dof= 0
    achi_list = np.array([])
    if (len(index_list) > 0):
	index_list = Array(index_list, copy=False)
        dof, achi_list = autochisq_from_precomputed(sn, cor, _P, stride=stride, \
			num_points=num_points, indices=index_list)
    return(dof, achi_list)
 

		    
	    
 




 

