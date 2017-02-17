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

from pycbc.filter import make_frequency_series
from pycbc.filter import  matched_filter_core
from pycbc.types import Array
import numpy as np
import logging

BACKEND_PREFIX="pycbc.vetoes.autochisq_"


def autochisq_from_precomputed(sn, corr_sn, hautocorr, indices,
                       stride=1, num_points=None, oneside=None,
                       twophase=True, maxvalued=False):
    """ 
    Compute correlation (two sided) between template and data
    and compares with autocorrelation of the template: C(t) = IFFT(A*A/S(f))

    Parameters
    ----------
    sn: Array[complex]
        normalized (!) array of complex snr for the template that produced the
        trigger(s) being tested
    corr_sn : Array[complex]
        normalized (!) array of complex snr for the template that you want to
        produce a correlation chisq test for. In the [common] case that sn and
        corr_sn are the same, you are computing auto-correlation chisq.
    hautocorr: Array[complex] 
        time domain autocorrelation for the template
    indices: Array[int]
        compute correlation chisquare at the points specified in this array,
    num_points: [int, optional; default=None]
        Number of points used for autochisq on each side, if None all points
        are used.
    stride: [int, optional; default = 1]
        stride for points selection for autochisq
        total length <= 2*num_points*stride
    oneside: [str, optional; default=None]
        whether to use one or two sided autochisquare. If None (or not
        provided) twosided chi-squared will be used. If given, options are
        'left' or 'right', to do one-sided chi-squared on the left or right.
    twophase: Boolean, optional; default=True
        If True calculate the auto-chisq using both phases of the filter.
        If False only use the phase of the obtained trigger(s).
    maxvalued: Boolean, optional; default=False
        Return the largest auto-chisq at any of the points tested if True.
        If False, return the sum of auto-chisq at all points tested.

    Returns
    -------
    autochisq: [tuple]
        returns autochisq values and snr corresponding to the instances 
        of time defined by indices
    """
    Nsnr = len(sn)

    indx = np.array([])

    achisq = np.zeros(len(indices))
    num_points_all = int(Nsnr/stride)
    if num_points is None:
        num_points = num_points_all
    if (num_points > num_points_all):
        num_points = num_points_all

    snrabs = np.abs(sn[indices])
    cphi_array = (sn[indices]).real / snrabs
    sphi_array = (sn[indices]).imag / snrabs

    start_point = - stride*num_points
    end_point = stride*num_points+1
    if oneside == 'left':
        achisq_idx_list = np.arange(start_point, 0, stride)
    elif oneside == 'right':
        achisq_idx_list = np.arange(stride, end_point, stride)
    else:
        achisq_idx_list_pt1 = np.arange(start_point, 0, stride)
        achisq_idx_list_pt2 = np.arange(stride, end_point, stride)
        achisq_idx_list = np.append(achisq_idx_list_pt1,
                                    achisq_idx_list_pt2)

    hauto_corr_vec = hautocorr[achisq_idx_list]
    hauto_norm = hauto_corr_vec.real*hauto_corr_vec.real
    # REMOVE THIS LINE TO REPRODUCE OLD RESULTS
    hauto_norm += hauto_corr_vec.imag*hauto_corr_vec.imag
    chisq_norm = 1.0 - hauto_norm

    for ip,ind in enumerate(indices):
        curr_achisq_idx_list = achisq_idx_list + ind
        
        cphi = cphi_array[ip]
        sphi = sphi_array[ip]
        # By construction, the other "phase" of the SNR is 0
        snr_ind =  sn[ind].real*cphi + sn[ind].imag*sphi

        # Wrap index if needed (maybe should fail in this case?)
        if curr_achisq_idx_list[0] < 0:
            curr_achisq_idx_list[curr_achisq_idx_list < 0] += Nsnr
        if curr_achisq_idx_list[-1] > (Nsnr - 1):
            curr_achisq_idx_list[curr_achisq_idx_list > (Nsnr-1)] -= Nsnr

        z = corr_sn[curr_achisq_idx_list].real*cphi + \
             corr_sn[curr_achisq_idx_list].imag*sphi
        dz = z - hauto_corr_vec.real*snr_ind
        curr_achisq_list = dz*dz/chisq_norm

        if twophase:
            chisq_norm = 1.0 - hauto_norm
            z = -corr_sn[curr_achisq_idx_list].real*sphi + \
                 corr_sn[curr_achisq_idx_list].imag*cphi
            dz = z - hauto_corr_vec.imag*snr_ind
            curr_achisq_list += dz*dz/chisq_norm
          
        if maxvalued:
            achisq[ip] = curr_achisq_list.max()
        else:
            achisq[ip] = curr_achisq_list.sum()

    dof = num_points
    if oneside is None:
        dof = dof * 2
    if twophase:
        dof = dof * 2

    return dof, achisq, indices

class SingleDetAutoChisq(object):
    """Class that handles precomputation and memory management for efficiently
    running the auto chisq in a single detector inspiral analysis.
    """	
    def __init__(self, stride, num_points, onesided=None, twophase=False,
                 reverse_template=False, take_maximum_value=False,
                 maximal_value_dof=None):
        """
        Initialize autochisq calculation instance

        Parameters
        -----------
        stride : int
            Number of sample points between points at which auto-chisq is
            calculated.
        num_points : int
            Number of sample points at which to calculate auto-chisq in each
            direction from the trigger
        onesided : optional, default=None, choices=['left','right']
            If None (default), calculate auto-chisq in both directions from the
            trigger. If left (backwards in time) or right (forwards in time)
            calculate auto-chisq only in that direction.
        twophase : optional, default=False
            If False calculate auto-chisq using only the phase of the trigger.
            If True, compare also against the orthogonal phase.
        reverse_template : optional, default=False
            If true, time-reverse the template before calculating auto-chisq.
            In this case this is more of a cross-correlation chisq than auto.
        take_maximum_value : optional, default=False
            If provided, instead of adding the auto-chisq value at each sample
            point tested, return only the maximum value.
        maximal_value_dof : int, required if using take_maximum_value
            If using take_maximum_value the expected value is not known. This
            value specifies what to store in the cont_chisq_dof output.
        """
        if stride > 0:
            self.do = True
            self.column_name = "cont_chisq"
            self.table_dof_name = "cont_chisq_dof"
            self.dof = num_points
            self.num_points = num_points
            self.stride = stride
            self.one_sided = onesided
            if (onesided is not None):
                self.dof = self.dof * 2
            self.two_phase = twophase
            if self.two_phase:
                self.dof = self.dof * 2
            self.reverse_template = reverse_template
            self.take_maximum_value=take_maximum_value
            if self.take_maximum_value:
                if maximal_value_dof is None:
                    err_msg = "Must provide the maximal_value_dof keyword "
                    err_msg += "argument if using the take_maximum_value "
                    err_msg += "option."
                    raise ValueError(err_msg)
                self.dof = maximal_value_dof
            
            self._autocor = None
            self._autocor_id = None
        else:
            self.do = False

    def values(self, sn, indices, template, psd, norm, stilde=None,
               low_frequency_cutoff=None, high_frequency_cutoff=None):
        """
        Calculate the auto-chisq at the specified indices.

        Parameters
        -----------
        sn : Array[complex]
            SNR time series of the template for which auto-chisq is being
            computed. Provided unnormalized. 
        indices : Array[int]
            List of points at which to calculate auto-chisq
        template : Pycbc template object 
            The template for which we are calculating auto-chisq
        psd : Pycbc psd object
            The PSD of the data being analysed
        norm : float
            The normalization factor to apply to sn
        stilde : Pycbc data object, needed if using reverse-template
            The data being analysed. Only needed if using reverse-template,
            otherwise ignored
        low_frequency_cutoff : float
            The lower frequency to consider in matched-filters
        high_frequency_cutoff : float
            The upper frequency to consider in matched-filters
        """
        if self.do and (len(indices) > 0):
            htilde = make_frequency_series(template)

            # Check if we need to recompute the autocorrelation
            key = (id(template), id(psd))
            if key != self._autocor_id:
                logging.info("Calculating autocorrelation")

                if not self.reverse_template:
                    Pt, _Ptilde, P_norm = matched_filter_core(htilde,
                              htilde, psd=psd,
                              low_frequency_cutoff=low_frequency_cutoff,
                              high_frequency_cutoff=high_frequency_cutoff)
                    Pt = Pt * (1./ Pt[0])
                    self._autocor = Array(Pt, copy=True)
                else:
                    Pt, _Ptilde, P_norm = matched_filter_core(htilde.conj(),
                              htilde, psd=psd,
                              low_frequency_cutoff=low_frequency_cutoff,
                              high_frequency_cutoff=high_frequency_cutoff)

                    # T-reversed template has same norm as forward template
                    # so we can normalize using that
                    # FIXME: Here sigmasq has to be cast to a float or the
                    #        code is really slow ... why??
                    norm_fac = P_norm / float(((template.sigmasq(psd))**0.5))
                    Pt *= norm_fac
                    self._autocor = Array(Pt, copy=True)
                self._autocor_id = key
            
            logging.info("...Calculating autochisquare")
            sn = sn*norm
            if self.reverse_template:
                assert(stilde is not None)
                asn, acor, ahnrm = matched_filter_core(htilde.conj(), stilde,
                                 low_frequency_cutoff=low_frequency_cutoff,
                                 high_frequency_cutoff=high_frequency_cutoff,
                                 h_norm=template.sigmasq(psd))
                correlation_snr = asn * ahnrm
            else:
                correlation_snr = sn

            achi_list = np.array([])
            index_list = np.array(indices)
            dof, achi_list, _ = autochisq_from_precomputed(sn, correlation_snr,
                               self._autocor, index_list, stride=self.stride,
                               num_points=self.num_points,
                               oneside=self.one_sided, twophase=self.two_phase,
                               maxvalued=self.take_maximum_value)
            self.dof = dof
            return achi_list

class SingleDetSkyMaxAutoChisq(SingleDetAutoChisq):
    """Stub for precessing auto chisq if anyone ever wants to code it up.
    """
    def __init__(self, *args, **kwds):
        super(SingleDetSkyMaxAutoChisq, self).__init__(*args, **kwds)

    def values(self, *args, **kwargs):
        if self.do:
            err_msg = "Precessing single detector sky-max auto chisq has not "
            err_msg += "been written. If you want to use it, why not help "
            err_msg += "write it?"
            raise NotImplementedError(err_msg)
        else:
            return None
