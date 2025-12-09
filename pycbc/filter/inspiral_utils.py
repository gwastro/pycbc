"""GPU-optimized utilities for inspiral search

This module contains GPU-specific optimized implementations for inspiral
template matching, designed to work directly with CuPy arrays and batched
operations.
"""

import cupy as cp
import numpy as np
from pycbc.types import Array, FrequencySeries, float32, complex64
from pycbc import DYN_RANGE_FAC


def template_triggers_gpu_batched(t_nums, bank, segments, matched_filter,
                                   power_chisq, sg_chisq, inj_filter_rejector,
                                   cluster_window, out_vals_ref, opt, psd_var=None):
    """GPU-optimized batched template processing for CUPY scheme
    
    This function is specifically designed for GPU processing with batched
    templates. It works directly with CuPy arrays and calls GPU kernels
    without unnecessary PyCBC Array wrapper overhead.
    
    Assumptions:
    - Using CUPYScheme
    - Processing multiple templates in batch (len(t_nums) > 1)
    - All templates are SPAtmplt approximant
    - batch_size evenly divides number of templates
    
    Parameters
    ----------
    t_nums : list of int
        Template indices to process in this batch
    bank : FilterBank
        Template bank
    segments : list
        List of frequency-domain data segments
    matched_filter : MatchedFilterControl
        Matched filter engine
    power_chisq : SingleDetPowerChisq
        Power chi-squared veto
    sg_chisq : SingleDetSGChisq
        Sine-Gaussian chi-squared veto
    inj_filter_rejector : InjFilterRejector
        Injection filter rejector
    cluster_window : int
        Clustering window in samples
    out_vals_ref : dict
        Reference dictionary for output values
    opt : Namespace
        Command-line options
    psd_var : array-like, optional
        PSD variation data
        
    Returns
    -------
    out_vals_all : list of list of dict
        Trigger data for each template
    tparams : list of dict
        Template parameters for each template
    """
    import copy
    
    num_templates = len(t_nums)
    tparams = []
    out_vals_all = [[] for _ in range(num_templates)]
    
    # Get template data as 2D cupy array (batch_size x freq_length)
    # This calls directly into batched template generation
    htilde_batch = _generate_templates_gpu(t_nums, bank)
    tparams = [bank.table[i] for i in t_nums]
    
    # Get the generated templates to extract proper parameters
    templates = bank[t_nums]
    
    # Get filter parameters for each template
    template_flow = cp.array([t.f_lower for t in templates], dtype=cp.float32)
    
    # Process each segment
    for s_num, stilde in enumerate(segments):
        # Skip if any template in batch should be rejected
        if not all(inj_filter_rejector.template_segment_checker(bank, t_num, stilde) 
                  for t_num in t_nums):
            continue
        
        # Get segment data as cupy array
        stilde_data = stilde._data  # Underlying cupy array
        psd_data = stilde.psd._data  # Underlying cupy array
        
        # Compute sigmasq for all templates
        sigmasqs = _compute_sigmasqs_gpu(htilde_batch, psd_data, template_flow, 
                                         stilde.psd.delta_f, templates)
        
        # Debug: compare with template method
        ref_sigmasq = templates[0].sigmasq(stilde.psd)
        print(f"DEBUG: Our sigmasqs[0]={float(sigmasqs[0]):.2e}, template.sigmasq()={float(ref_sigmasq):.2e}, ratio={float(sigmasqs[0]/ref_sigmasq):.2f}")
        
        # Batched matched filtering - direct kernel call
        # Note: stilde is already overwhitened (divided by PSD)
        snr_batch, corr_batch = _matched_filter_gpu(htilde_batch, stilde_data, 
                                                     matched_filter.kmin, 
                                                     matched_filter.kmax)
        
        # The SNR normalization for overwhitened data
        # SNR = IFFT(conj(h) * s) / sqrt(sigmasq / (4 * delta_f))
        # Simplifying: SNR = IFFT_result * sqrt(4 * delta_f) / sqrt(sigmasq)
        # But let's match what PyCBC does: just divide by sqrt(sigmasq)
        snr_batch = snr_batch / cp.sqrt(sigmasqs[:, cp.newaxis])
        
        # Debug: check SNR values
        max_snr = float(cp.max(cp.abs(snr_batch)))
        print(f"DEBUG: Max SNR: {max_snr:.2f}")
        
        # Debug: check SNR values
        max_snr = float(cp.max(cp.abs(snr_batch)))
        print(f"DEBUG: Max SNR in batch: {max_snr:.2f}, sigmasqs: {sigmasqs[:3]}")
        
        # Threshold and cluster - operates on 2D arrays
        trigger_indices, trigger_snrs = _threshold_and_cluster_gpu(
            snr_batch, opt.snr_threshold, cluster_window, stilde.analyze)
        
        print(f"DEBUG: Found {len(trigger_indices)} total triggers above threshold {opt.snr_threshold}")
        
        # Compute chi-squared for triggers
        chisqs_batch, chisq_dofs_batch = _compute_chisq_gpu(
            corr_batch, trigger_indices, trigger_snrs, sigmasqs,
            psd_data, power_chisq, stilde.analyze.start)
        
        # Extract results for each template
        for tidx in range(num_templates):
            # Get triggers for this template
            tidx_mask = trigger_indices[:, 0] == tidx
            if not cp.any(tidx_mask):
                continue
                
            idx = trigger_indices[tidx_mask, 1]
            snrv = trigger_snrs[tidx_mask]
            chisq = chisqs_batch[tidx_mask]
            chisq_dof = chisq_dofs_batch[tidx_mask]
            
            # Build output dictionary
            out_vals = out_vals_ref.copy()
            out_vals['chisq'] = chisq
            out_vals['chisq_dof'] = chisq_dof
            
            # TODO: Implement GPU SG chisq
            # For now, create dummy values
            out_vals['sg_chisq'] = cp.zeros(len(idx), dtype=float32)
            
            # Add segment offset to indices
            idx = idx + stilde.cumulative_index
            
            out_vals['time_index'] = idx
            out_vals['snr'] = snrv
            out_vals['sigmasq'] = cp.zeros(len(snrv), dtype=float32) + sigmasqs[tidx]
            
            if opt.psdvar_short_segment is not None:
                import pycbc.psd
                out_vals['psd_var_val'] = \
                    pycbc.psd.find_trigger_value(psd_var, idx,
                                               opt.gps_start_time,
                                               opt.sample_rate)
            
            out_vals_all[tidx].append(copy.deepcopy(out_vals))
    
    return out_vals_all, tparams


def _generate_templates_gpu(t_nums, bank):
    """Generate batch of templates directly as 2D cupy array
    
    Parameters
    ----------
    t_nums : list of int
        Template indices
    bank : FilterBank
        Template bank
        
    Returns
    -------
    cupy.ndarray
        2D array of templates (num_templates x freq_length)
    """
    # Use the bank's batched generation
    templates = bank[t_nums]
    
    # Extract underlying cupy data into 2D array
    htilde_batch = cp.stack([t._data for t in templates])
    
    return htilde_batch


def _compute_sigmasqs_gpu(htilde_batch, psd_data, template_flow, delta_f, templates):
    """Compute sigmasq for batch of templates using GPU
    
    Parameters
    ----------
    htilde_batch : cupy.ndarray
        2D array of templates (num_templates x freq_length)
    psd_data : cupy.ndarray
        PSD array
    template_flow : cupy.ndarray
        Lower frequency cutoff for each template
    delta_f : float
        Frequency spacing
    templates : list
        List of template FrequencySeries objects (for metadata)
        
    Returns
    -------
    cupy.ndarray
        Array of sigmasq values for each template
    """
    num_templates = htilde_batch.shape[0]
    freq_length = min(htilde_batch.shape[1], len(psd_data))
    sigmasqs = cp.zeros(num_templates, dtype=cp.float32)
    
    # Compute sigmasq for each template
    # sigmasq = 4 * delta_f * sum(|h[k]|^2 / S[k]) for k in [kmin, kmax)
    for i in range(num_templates):
        # Get frequency range from template metadata
        kmin = int(templates[i].f_lower / delta_f)
        kmax = templates[i].end_idx if hasattr(templates[i], 'end_idx') else freq_length
        kmax = min(kmax, freq_length)
        
        # Compute integrand
        htilde = htilde_batch[i, :freq_length]
        integrand = cp.abs(htilde[kmin:kmax])**2 / psd_data[kmin:kmax]
        sigmasqs[i] = 4.0 * delta_f * cp.sum(integrand)
    
    return sigmasqs


def _matched_filter_gpu(htilde_batch, stilde_data, kmin, kmax):
    """Batched matched filtering on GPU
    
    Parameters
    ----------
    htilde_batch : cupy.ndarray
        2D array of templates (num_templates x freq_length)
    stilde_data : cupy.ndarray
        Strain data in frequency domain
    kmin : int
        Minimum frequency index
    kmax : int
        Maximum frequency index
        
    Returns
    -------
    snr_batch : cupy.ndarray
        2D array of SNR time series (num_templates x time_length)
    corr_batch : cupy.ndarray
        2D array of correlation in frequency domain (num_templates x freq_length)
    """
    num_templates = htilde_batch.shape[0]
    freq_length = htilde_batch.shape[1]
    
    # Allocate output arrays
    corr_batch = cp.zeros((num_templates, freq_length), dtype=cp.complex64)
    
    # Ensure stilde_data length matches
    stilde_len = min(len(stilde_data), freq_length)
    
    # Batched correlation: corr = conj(h) * s
    corr_batch[:, kmin:kmax] = cp.conj(htilde_batch[:, kmin:kmax]) * stilde_data[kmin:kmax]
    
    # Debug: check correlation
    max_corr = float(cp.max(cp.abs(corr_batch)))
    print(f"DEBUG MF: max corr={max_corr:.2e}, kmin={kmin}, kmax={kmax}, freq_length={freq_length}")
    print(f"DEBUG MF: htilde nonzero={cp.count_nonzero(htilde_batch)}, stilde nonzero={cp.count_nonzero(stilde_data)}")
    
    # Batched IFFT to get SNR time series
    snr_batch = cp.fft.ifft(corr_batch, axis=1)
    delta_f = 1./256
    norm = (4.0 * delta_f)
    snr_batch *= norm
    
    # Debug: check IFFT result
    max_snr_before_norm = float(cp.max(cp.abs(snr_batch)))
    print(f"DEBUG MF: max SNR before norm={max_snr_before_norm:.2e}")
    
    return snr_batch, corr_batch


def _normalize_snr_gpu(snr_batch, sigmasqs):
    """Normalize SNR by sigmasq
    
    Parameters
    ----------
    snr_batch : cupy.ndarray
        2D array of SNR time series
    sigmasqs : cupy.ndarray
        Array of sigmasq values
        
    Returns
    -------
    cupy.ndarray
        Normalized SNR time series
    """
    # Reshape sigmasqs for broadcasting
    norm = cp.sqrt(sigmasqs[:, cp.newaxis])
    
    # Debug
    print(f"DEBUG NORM: norm shape={norm.shape}, snr shape={snr_batch.shape}, norm[0]={float(norm[0]):.2e}")
    
    result = snr_batch / norm
    
    max_result = float(cp.max(cp.abs(result)))
    print(f"DEBUG NORM: max normalized SNR={max_result:.2e}")
    
    return result


def _threshold_and_cluster_gpu(snr_batch, threshold, window, analyze_segment):
    """Threshold and cluster SNR time series
    
    Parameters
    ----------
    snr_batch : cupy.ndarray
        2D array of SNR time series (num_templates x time_length)
    threshold : float
        SNR threshold
    window : int
        Clustering window in samples
    analyze_segment : slice
        Segment to analyze
        
    Returns
    -------
    trigger_indices : cupy.ndarray
        2D array of trigger indices (num_triggers x 2) [template_idx, time_idx]
    trigger_snrs : cupy.ndarray
        Array of trigger SNR values
    """
    # Get absolute SNR
    abs_snr = cp.abs(snr_batch)
    
    # Extract the analyze segment
    if isinstance(analyze_segment, slice):
        snr_analyze = abs_snr[:, analyze_segment]
        segment_start = analyze_segment.start if analyze_segment.start is not None else 0
    else:
        snr_analyze = abs_snr
        segment_start = 0
    
    # Apply threshold
    above_thresh = snr_analyze > threshold
    
    # Get indices of triggers
    template_idxs, time_idxs = cp.where(above_thresh)
    
    if len(template_idxs) == 0:
        return cp.array([], dtype=cp.int32).reshape(0, 2), cp.array([], dtype=cp.complex64)
    
    # Adjust time indices to absolute positions
    time_idxs = time_idxs + segment_start
    
    # Get SNR values at trigger locations
    trigger_snrs = snr_batch[template_idxs, time_idxs]
    
    # Cluster triggers (simple maximum in window for now)
    # TODO: Implement proper clustering
    trigger_indices = cp.stack([template_idxs, time_idxs], axis=1)
    
    return trigger_indices, trigger_snrs


def _compute_chisq_gpu(corr_batch, trigger_indices, trigger_snrs, sigmasqs,
                       psd_data, power_chisq, analyze_start):
    """Compute chi-squared for triggers
    
    Parameters
    ----------
    corr_batch : cupy.ndarray
        2D array of correlations
    trigger_indices : cupy.ndarray
        Trigger indices (num_triggers x 2)
    trigger_snrs : cupy.ndarray
        Trigger SNR values
    sigmasqs : cupy.ndarray
        Template sigmasq values
    psd_data : cupy.ndarray
        PSD data
    power_chisq : SingleDetPowerChisq
        Power chisq calculator
    analyze_start : int
        Start of analyze segment
        
    Returns
    -------
    chisqs : cupy.ndarray
        Chi-squared values
    chisq_dofs : cupy.ndarray
        Chi-squared degrees of freedom
    """
    num_triggers = len(trigger_indices)
    
    if num_triggers == 0:
        return cp.array([], dtype=cp.float32), cp.array([], dtype=cp.int32)
    
    # For now, return dummy values
    # TODO: Implement batched chisq computation
    chisqs = cp.ones(num_triggers, dtype=cp.float32)
    chisq_dofs = cp.ones(num_triggers, dtype=cp.int32) * 10
    
    return chisqs, chisq_dofs


def _compute_sigmasqs_batch(templates, psd):
    """Compute sigmasq for a batch of SPAtmplt templates efficiently
    
    Parameters
    ----------
    templates : list
        List of template waveforms
    psd : FrequencySeries
        Power spectral density
        
    Returns
    -------
    list
        List of sigmasq values for each template
    """
    import pycbc.waveform
    from pycbc import DYN_RANGE_FAC
    
    # Ensure the sigmasq_vec is computed for SPAtmplt
    if not hasattr(psd, 'sigmasq_vec'):
        psd.sigmasq_vec = {}
    
    if 'SPAtmplt' not in psd.sigmasq_vec:
        psd.sigmasq_vec['SPAtmplt'] = \
            pycbc.waveform.get_waveform_filter_norm(
                'SPAtmplt',
                psd,
                len(psd),
                psd.delta_f,
                templates[0].min_f_lower
            )
    
    curr_sigmasq = psd.sigmasq_vec['SPAtmplt']
    sigmasqs = []
    
    for template in templates:
        # Compute sigma_scale if not already done
        if not hasattr(template, 'sigma_scale'):
            amp_norm = pycbc.waveform.get_template_amplitude_norm(
                template.params, approximant='SPAtmplt')
            amp_norm = 1 if amp_norm is None else amp_norm
            template.sigma_scale = (DYN_RANGE_FAC * amp_norm) ** 2.0
        
        kmin = int(template.f_lower / psd.delta_f)
        sigmasq = template.sigma_scale * \
            (curr_sigmasq[template.end_idx-1] - curr_sigmasq[kmin])
        sigmasqs.append(sigmasq)
    
    return sigmasqs
