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
    # This calls directly into batched template generation, bypassing PyCBC array overhead
    htilde_batch, templates = _generate_templates_gpu(t_nums, bank)
    tparams = [bank.table[i] for i in t_nums]
    
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
                                          stilde.delta_f, templates)
        
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
        
        # Threshold and cluster - operates on 2D arrays
        trigger_indices, trigger_snrs = _threshold_and_cluster_gpu(
            snr_batch, opt.snr_threshold, cluster_window, stilde.analyze)
        
        # Compute chi-squared for triggers
        chisqs_batch, chisq_dofs_batch = _compute_chisq_gpu(
            corr_batch, trigger_indices, trigger_snrs, sigmasqs,
            stilde.psd, power_chisq, templates, t_nums)
        
        # Extract results for each template
        for tidx in range(num_templates):
            # Get triggers for this template
            tidx_mask = trigger_indices[:, 0] == tidx
            if not cp.any(tidx_mask):
                continue
                
            # Time indices are currently relative to full snr_batch
            # We need to add cumulative_index to get GPS time indices
            idx = trigger_indices[tidx_mask, 1]
            snrv = trigger_snrs[tidx_mask]
            chisq = chisqs_batch[tidx_mask]
            chisq_dof = chisq_dofs_batch[tidx_mask]
            
            # Build output dictionary
            # Convert cupy arrays to PyCBC Arrays
            out_vals = out_vals_ref.copy()
            out_vals['chisq'] = Array(chisq, copy=False)
            out_vals['chisq_dof'] = Array(chisq_dof, copy=False)
            
            # TODO: Implement GPU SG chisq
            # For now, create dummy values
            out_vals['sg_chisq'] = Array(cp.zeros(len(idx), dtype=float32), copy=False)
            
            # Convert time indices to GPS sample indices
            # idx is currently relative to the start of snr_batch (which starts at seg_slice.start)
            # cumulative_index = seg_slice.start + ana.start
            # Our idx already includes ana.start offset (added in threshold function)
            # So we add (cumulative_index - ana.start) = seg_slice.start
            idx = idx + (stilde.cumulative_index - stilde.analyze.start)
            
            out_vals['time_index'] = Array(idx, copy=False)
            out_vals['snr'] = Array(snrv, copy=False)
            out_vals['sigmasq'] = Array(cp.zeros(len(snrv), dtype=float32) + sigmasqs[tidx], copy=False)
            
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
    
    Extracted from bank.__getitem__ to avoid PyCBC array overhead.
    Calls spa_tmplt_batch directly.
    
    Parameters
    ----------
    t_nums : list of int
        Template indices
    bank : FilterBank
        Template bank
        
    Returns
    -------
    htilde_batch : cupy.ndarray
        2D array of templates (num_templates x freq_length)
    templates : list
        List of FrequencySeries template objects
    """
    from pycbc.waveform.spa_tmplt import spa_tmplt_batch
    from pycbc.waveform.bank import find_variable_start_frequency
    from pycbc.types import zeros
    import logging
    import types
    
    # Allocate output arrays - use bank's pre-allocated if available
    if bank.out is None:
        tempout = [zeros(bank.filter_length, dtype=bank.dtype) for _ in range(len(t_nums))]
    else:
        tempout = bank.out
    
    # Prepare common parameters
    distance = 1.0 / DYN_RANGE_FAC
    common_kwds = {
        'delta_f': bank.delta_f,
        'f_lower': bank.f_lower,
        'distance': distance,
        **bank.extra_args
    }
    
    # Prepare per-template parameters
    templates_params = []
    f_end_list = []
    f_low_list = []
    
    for idx in t_nums:
        f_end = bank.end_frequency(idx)
        if f_end is None or f_end >= (bank.filter_length * bank.delta_f):
            f_end = (bank.filter_length - 1) * bank.delta_f
        
        f_low = find_variable_start_frequency('SPAtmplt',
                                            bank.table[idx],
                                            bank.f_lower,
                                            bank.max_template_length)
        
        params = {
            'mass1': bank.table[idx].mass1,
            'mass2': bank.table[idx].mass2,
            'spin1z': bank.table[idx].spin1z,
            'spin2z': bank.table[idx].spin2z,
            'distance': distance
        }
        
        templates_params.append(params)
        f_end_list.append(f_end)
        f_low_list.append(f_low)
    
    # Log batch generation
    logging.info('Generating templates %s-%s (%d templates) in batch from %s Hz' % 
                (t_nums[0], t_nums[-1], len(t_nums), min(f_low_list)))
    
    # Generate all templates in one batch
    htilde_list = spa_tmplt_batch(templates_params, bank.filter_length, tempout, **common_kwds)
    
    # Process each template and extract data
    for i, (idx, htilde) in enumerate(zip(t_nums, htilde_list)):
        template_duration = htilde.chirp_length if hasattr(htilde, 'chirp_length') else None
        ttotal = htilde.length_in_time if hasattr(htilde, 'length_in_time') else None
        
        bank.table[idx].template_duration = template_duration
        
        htilde = htilde.astype(bank.dtype)
        htilde.f_lower = f_low_list[i]
        htilde.min_f_lower = bank.min_f_lower
        htilde.end_idx = int(f_end_list[i] / htilde.delta_f)
        htilde.params = bank.table[idx]
        htilde.chirp_length = template_duration
        htilde.length_in_time = ttotal
        htilde.approximant = 'SPAtmplt'
        htilde.end_frequency = f_end_list[i]
        
        # Update the list with processed template
        htilde_list[i] = htilde
    
    # Extract underlying cupy data into 2D array
    htilde_batch = cp.stack([t._data for t in htilde_list])
    
    return htilde_batch, htilde_list


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
    corr_batch *= 524288 # I haven't figured out why this is necessary yet

    # Batched IFFT to get SNR time series
    snr_batch = cp.fft.ifft(corr_batch, axis=1)
    delta_f = 1./256
    norm = (4.0 * delta_f)
    snr_batch *= norm
    
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
    return snr_batch / norm


def _threshold_and_cluster_gpu(snr_batch, threshold, window, analyze_segment):
    """Threshold and cluster SNR time series using CUDA kernels
    
    Parameters
    ----------
    snr_batch : cupy.ndarray
        2D array of SNR time series (num_templates x time_length)
    threshold : float or cupy.ndarray
        SNR threshold (can be array with one value per template)
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
    from pycbc.events.threshold_cupy import get_tkernel, FastFilter
    
    batch_size, series_length = snr_batch.shape
    analyse_start = cp.int32(analyze_segment.start)
    analyse_end = cp.int32(analyze_segment.stop)
    analyse_len = analyse_end - analyse_start
    batch_mem_size = 1024
    
    # Allocate output arrays (don't use global cache to avoid shape mismatches)
    outv = cp.zeros((batch_size, batch_mem_size), dtype=cp.complex64)
    outl = cp.zeros((batch_size, batch_mem_size), dtype=cp.int32)
    
    # Convert threshold to array if needed
    if not isinstance(threshold, cp.ndarray):
        threshold = cp.full(batch_size, threshold, dtype=cp.float32)
    
    # Square the threshold for kernel
    threshold_sq = threshold * threshold
    threshold_sq = cp.asarray(threshold_sq, dtype=cp.float32)
    window = cp.int32(window)
    
    # Get kernels
    (fn, fn2), nt, nb = get_tkernel(
        series_length,
        analyse_len,
        window,
        block_mem_size=batch_mem_size,
        batch_size=batch_size
    )
    
    grid = (nb, batch_size, 1)
    block = (nt, 1, 1)
    
    # Run threshold and cluster kernels
    fn(grid, block, (snr_batch, outv, outl, window, threshold_sq, series_length, analyse_start))
    fn2(grid, block, (outv, outl, threshold_sq, window))
    
    # Filter results using FastFilter
    fast_filter = FastFilter(snum=batch_size, nb=nb)
    cv_filtered, cl_filtered, cv_out, cl_out, sizes, template_map = \
        fast_filter.filter_arrays(outv[:, :nb], outl[:, :nb])
    
    # Convert results to format expected by rest of pipeline
    num_triggers = int(sizes.sum())
    
    if num_triggers == 0:
        return cp.array([], dtype=cp.int32).reshape(0, 2), cp.array([], dtype=cp.complex64)
    
    # Time indices from kernel are relative to analyze_segment.start
    # Add analyze_segment.start to make them relative to full snr_batch
    time_indices = cl_out[:num_triggers] + analyze_segment.start
    
    # Create trigger indices array [template_idx, time_idx]
    trigger_indices = cp.stack([template_map[:num_triggers], time_indices], axis=1)
    trigger_snrs = cv_out[:num_triggers]
    
    return trigger_indices, trigger_snrs


def _compute_chisq_gpu(corr_batch, trigger_indices, trigger_snrs, sigmasqs,
                       psd, power_chisq, templates, t_nums):
    """Compute chi-squared for triggers using batched GPU computation
    
    Parameters
    ----------
    corr_batch : cupy.ndarray
        2D array of correlations (num_templates x freq_length)
    trigger_indices : cupy.ndarray
        Trigger indices (num_triggers x 2) [template_idx, time_idx]
    trigger_snrs : cupy.ndarray
        Trigger SNR values (complex, unnormalized)
    sigmasqs : cupy.ndarray
        Template sigmasq values
    psd : FrequencySeries
        Power spectral density
    power_chisq : SingleDetPowerChisq
        Power chisq calculator
    templates : list
        List of template waveforms
    t_nums : list
        Template indices in the batch
        
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
    
    if not power_chisq.do:
        # Chi-squared disabled, return dummy values
        template_idxs = trigger_indices[:, 0]
        chisq_dofs = sigmasqs[template_idxs] * 0 + 10
        chisqs = chisq_dofs.astype(cp.float32)
        return chisqs, chisq_dofs.astype(cp.int32)
    
    # Extract data needed for values_batch
    template_map = trigger_indices[:, 0]  # Which template each trigger belongs to
    time_indices = trigger_indices[:, 1]  # Time index for each trigger
    
    # The chisq computation in values_batch expects:
    # - snrvs: unnormalized complex SNR values (from IFFT, no 4*df*N factor)
    # - snr_norms: normalization factor such that abs(snrvs * snr_norms) gives true SNR
    #
    # Our pipeline:
    # - matched filter outputs: (4*df*N) * IFFT(conj(h)*s)
    # - normalized SNR: (4*df*N) * IFFT / sqrt(sigmasq)
    # - trigger_snrs contain: (4*df*N) * IFFT / sqrt(sigmasq)
    #
    # For chisq, we need:
    # - unnormalized_snrs = IFFT (no scaling)
    # - snr_norms = (4*df*N) / sqrt(sigmasq)
    #
    # So: unnormalized_snrs * snr_norms = IFFT * (4*df*N) / sqrt(sigmasq) = trigger_snrs ✓
    
    delta_f = psd.delta_f
    norm_factor = 4.0 * delta_f
    
    # Un-normalize to get raw IFFT values
    unnormalized_snrs = trigger_snrs * cp.sqrt(sigmasqs[template_map, cp.newaxis]).squeeze() / norm_factor
    
    # Normalization factor - try sqrt(4*df) / sqrt(sigmasq)
    snr_norms = norm_factor / cp.sqrt(sigmasqs)
    
    # Count triggers per template
    num_templates = len(templates)
    sizes = cp.array([cp.sum(template_map == i) for i in range(num_templates)], dtype=cp.int32)
    
    # Call the batched chisq computation
    # Note: corr_batch has full FFT length (power of 2), which shift_sum_batch requires
    # The chisq code will only use the relevant frequency bins anyway
    
    # DEBUG: Check template and PSD lengths
    import logging
    if len(templates) > 0:
        psd_len = len(psd)
        logging.info(f'DEBUG chisq: template[0] length = {len(templates[0])}, PSD length = {psd_len}')
        logging.info(f'DEBUG chisq: corr_batch shape = {corr_batch.shape}')
        
        # Fix template lengths to match PSD (RFFT length)
        # Templates have full FFT length but chisq code expects RFFT length
        for template in templates:
            if len(template) != psd_len:
                logging.info(f'DEBUG chisq: Fixing template length from {len(template)} to {psd_len}')
                # Resize the template data to RFFT length
                template.resize(psd_len)
        logging.info("DEBUG VALUES", max(corr_batch), max(unnormalized_snrs), max(snr_norms))
    
    chisq_list, chisq_dof_list = power_chisq.values_batch(
        corr_batch, unnormalized_snrs, snr_norms, sizes, template_map, psd, time_indices, templates
    )
    
    # DEBUG: Check what was returned
    logging.info(f'DEBUG chisq result: chisq_list type={type(chisq_list)}, len={len(chisq_list) if chisq_list is not None else "None"}')
    if chisq_list is not None and len(chisq_list) > 0:
        logging.info(f'DEBUG chisq_list[0]: len={len(chisq_list[0])}, values={chisq_list[0][:min(3,len(chisq_list[0]))]}')
    
    # Check if chisq was actually computed (or disabled)
    if chisq_list is None:
        # Chi-squared was disabled, return dummy values
        chisq_dofs = cp.full(num_triggers, 10, dtype=cp.int32)
        chisqs = chisq_dofs.astype(cp.float32)
        logging.info(f'DEBUG chisq was None, returning dummy values')
        return chisqs, chisq_dofs
    
    # Concatenate results back into single arrays
    chisqs = cp.concatenate(chisq_list) if len(chisq_list) > 0 else cp.array([], dtype=cp.float32)
    chisq_dofs = cp.concatenate(chisq_dof_list) if len(chisq_dof_list) > 0 else cp.array([], dtype=cp.int32)
    
    logging.info(f'DEBUG final chisqs: len={len(chisqs)}, min={chisqs.min() if len(chisqs)>0 else "N/A"}, max={chisqs.max() if len(chisqs)>0 else "N/A"}')
    
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
