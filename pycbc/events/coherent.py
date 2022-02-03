# Copyright (C) 2022 Andrew Williamson
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
""" This module contains functions for calculating and manipulating coherent
triggers.
"""

import numpy as np


def get_coinc_indexes(idx_dict, time_delay_idx):
    """Return the indexes corresponding to coincident triggers

    Parameters
    ----------
    idx_dict: dict
        Dictionary of indexes of triggers above threshold in each
        detector
    time_delay_idx: dict
        Dictionary giving time delay index (time_delay*sample_rate) for
        each ifo

    Returns
    -------
    coinc_idx: list
        List of indexes for triggers in geocent time that appear in
        multiple detectors
    """
    coinc_list = np.array([], dtype=int)
    for ifo in idx_dict.keys():
        # Create list of indexes above threshold in single detector in geocent
        # time. Can then search for triggers that appear in multiple detectors
        # later.
        if len(idx_dict[ifo]) != 0:
            coinc_list = np.hstack([coinc_list,
                                    idx_dict[ifo] - time_delay_idx[ifo]])
    # Search through coinc_idx for repeated indexes. These must have been loud
    # in at least 2 detectors.
    coinc_idx = np.unique(coinc_list, return_counts=True)[0][
        np.unique(coinc_list, return_counts=True)[1] > 1]
    return coinc_idx


def coincident_snr(snr_dict, index, threshold, time_delay_idx):
    """Calculate the coincident SNR for all coinciden triggers above
    threshold

    Parameters
    ----------
    snr_dict: dict
        Dictionary of individual detector SNRs
    index: list
        List of indexes (geocentric) for which to calculate coincident
        SNR
    threshold: float
        Coincident SNR threshold. Triggers below this are cut
    time_delay_idx: dict
        Dictionary of time delay from geocenter in indexes for each
        detector

    Returns
    -------
    rho_coinc: numpy.ndarray
        Coincident SNR values for surviving triggers
    index: list
        The subset of input indexes corresponding to triggers that
        survive the cuts
    coinc_triggers: dict
        Dictionary of individual detector SNRs for triggers that
        survive cuts
    """
    # Restrict the snr timeseries to just the interesting points
    coinc_triggers = {ifo: snr_dict[ifo][index + time_delay_idx[ifo]]
                      for ifo in snr_dict.keys()}
    # Calculate the coincident snr
    snr_array = np.array([coinc_triggers[ifo]
                          for ifo in coinc_triggers.keys()])
    rho_coinc = np.sqrt(np.sum(snr_array * snr_array.conj(), axis=0))
    # Apply threshold
    thresh_indexes = rho_coinc > threshold
    index = index[thresh_indexes]
    coinc_triggers = {ifo: snr_dict[ifo][index + time_delay_idx[ifo]]
                      for ifo in snr_dict.keys()}
    rho_coinc = rho_coinc[thresh_indexes]
    return rho_coinc, index, coinc_triggers


def get_projection_matrix(fp, fc, sigma, projection='standard'):
    """ Calculate the matrix that prjects the signal onto the network.

    Parameters
    ----------
    fp: dict
        Dictionary containing the plus antenna response factors for
        each IFO
    fc: dict
        Dictionary containing the cross antenna response factors for
        each IFO
    sigma: dict
        Dictionary of the sensitivity weights for each IFO
    projection: optional, {string, 'standard'}
        The signal polarization to project. Choice of 'standard'
        (unrestricted; default), 'right' or 'left' (circular
        polarizations)

    Returns
    -------
    projection_matrix: np.ndarray
        The matrix that projects the signal onto the detector network
    """
    # Calculate the weighted antenna responses
    keys = sorted(sigma.keys())
    wp = np.array([sigma[ifo] * fp[ifo] for ifo in keys])
    wc = np.array([sigma[ifo] * fc[ifo] for ifo in keys])

    # Get the projection matrix associated with the requested projection
    if projection == 'standard':
        denominator = np.dot(wp, wp) * np.dot(wc, wc) - np.dot(wp, wc)**2
        projection_matrix = (np.dot(wc, wc)*np.outer(wp, wp) +
                             np.dot(wp, wp)*np.outer(wc, wc) -
                             np.dot(wp, wc)*(np.outer(wp, wc) +
                             np.outer(wc, wp))) / denominator
    elif projection == 'left':
        projection_matrix = ((np.outer(wp, wp) + np.outer(wc, wc) +
                             (np.outer(wp, wc) - np.outer(wc, wp)) * 1j)
                             / (np.dot(wp, wp) + np.dot(wc, wc)))
    elif projectioni == 'right':
        projection_matrix = ((np.outer(wp, wp) + np.outer(wc, wc) +
                             (np.outer(wc, wp) - np.outer(wp, wc)) * 1j)
                             / (np.dot(wp, wp) + np.dot(wc, wc)))
    else:
        raise ValueError('Unknown projection: {}'.format(projection))

    return projection_matrix


def coherent_snr(snr_triggers, index, threshold, projection_matrix,
                 coinc_snr=[]):
    """Calculate the coherent SNR for a given set of triggers

    Parameters
    ----------
    snr_triggers: dict
        Dictionary of the normalised complex snr time series for each
        ifo
    index: numpy.ndarray
        Array of the indexes corresponding to triggers
    threshold: float
        Coherent SNR threshold. Triggers below this are cut
    projection_matrix: numpy.ndarray
        Matrix that projects the signal onto the network
    coinc_snr: Optional- The coincident snr for each trigger.

    Returns
    -------
    rho_coh: numpy.ndarray
        Array of coherent SNR for the detector network
    index: numpy.ndarray
        Indexes that survive cuts
    snrv: dict
        Dictionary of individual deector triggers that survive cuts
    coinc_snr: list
        The coincident SNR values for triggers surviving the coherent
        cut
    """
    # Calculate rho_coh
    snr_array = np.array([snr_triggers[ifo]
                          for ifo in sorted(snr_triggers.keys())])
    x = np.inner(snr_array.conj().transpose(), projection_matrix)
    rho_coh2 = sum(x.transpose() * snr_array)
    rho_coh = np.sqrt(rho_coh2)
    # Apply thresholds
    index = index[rho_coh > threshold]
    if coinc_snr:
        coinc_snr = coinc_snr[rho_coh > threshold]
    snrv = {ifo: snr_triggers[ifo][rho_coh > threshold]
            for ifo in snr_triggers.keys()}
    rho_coh = rho_coh[rho_coh > threshold]
    return rho_coh, index, snrv, coinc_snr


def network_chisq(chisq, chisq_dof, snr_dict):
    """Calculate the network chi-squared statistic

    Parameters
    ----------
    chisq: dict
        Dictionary of individual detector chi-squared statistics
    chisq_dof: dict
        Dictionary of the number of degrees of freedom of the
        chi-squared statistic
    snr_dict: dict
        Dictionary of complex individual detector SNRs

    Returns
    -------
    net_chisq: list
        Network chi-squared values
    """
    ifos = sorted(snr_dict.keys())
    chisq_per_dof = {}
    for ifo in ifos:
        chisq_per_dof[ifo] = chisq[ifo] / chisq_dof[ifo]
        chisq_per_dof[ifo][chisq_per_dof[ifo] < 1] = 1
    snr2 = {ifo: np.real(np.array(snr_dict[ifo]) *
                 np.array(snr_dict[ifo]).conj())
            for ifo in ifos}
    coinc_snr2 = sum(snr2.values())
    snr2_ratio = {ifo: snr2[ifo] / coinc_snr2 for ifo in ifos}
    net_chisq = sum([chisq_per_dof[ifo] * snr2_ratio[ifo] for ifo in ifos])
    return net_chisq


def reweighted_snr(netwk_snr, netwk_chisq, a=3, b=1./6.):
    """
    Output: reweighted_snr: Reweighted SNR for each trigger
    Input:  netwk_snr:  Dictionary of coincident or coherent SNR for each
                        trigger
            netwk_chisq: A chisq value for each trigger
    """
    denom = ((1 + netwk_chisq)**a) / 2
    reweighted_snr = netwk_snr / denom**b
    return reweighted_snr


def null_snr(rho_coh, rho_coinc, null_min=5.25, null_grad=0.2, null_step=20.,
             index={}, snrv={}):
    """
    Output: null: null snr for surviving triggers
            rho_coh: Coherent snr for surviving triggers
            rho_coinc: Coincident snr for suviving triggers
            index: Indexes for surviving triggers
            snrv: Single detector snr for surviving triggers
    Input:  rho_coh: Numpy array of coherent snr triggers
            rho_coinc: Numpy array of coincident snr triggers
            null_min: Any trigger with null snr below this is cut
            null_grad: Any trigger with null snr<(null_grad*rho_coh+null_min)
                       is cut
            null_step: The value for required for coherent snr to start
                       increasing the null threshold
            index: Optional- Indexes of triggers. If given, will remove
                   triggers that fail cuts
            snrv: Optional- Individual ifo snr for triggers. If given will
                  remove triggers that fail cut
    """
    null2 = rho_coinc**2 - rho_coh**2
    # Numerical errors may make this negative and break the sqrt, so set
    # negative values to 0.
    null2[null2 < 0] = 0
    null = null2**0.5
    # Make cut on null.
    keep1 = np.logical_and(null < null_min, rho_coh <= null_step)
    keep2 = np.logical_and(null < (rho_coh * null_grad + null_min),
                           rho_coh > null_step)
    keep = np.logical_or(keep1, keep2)
    index = index[keep]
    rho_coh = rho_coh[keep]
    snrv = {ifo: snrv[ifo][keep] for ifo in snrv}
    rho_coinc = rho_coinc[keep]
    null = null[keep]
    return null, rho_coh, rho_coinc, index, snrv


def reweight_snr_by_null(network_snr, nullsnr):
    """
    Output: reweighted_snr: Reweighted SNR for each trigger
    Input:  network_snr:  Dictionary of coincident, coherent, or reweighted
                          SNR for each trigger
            null: Null snr for each trigger
    """
    nullsnr = np.array(nullsnr)
    nullsnr[nullsnr <= 4.25] = 4.25
    reweighted_snr = network_snr / (nullsnr - 3.25)
    return reweighted_snr
