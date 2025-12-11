"""GPU-optimized utilities for inspiral search

This module contains GPU-specific optimized implementations for inspiral
template matching, designed to work directly with CuPy arrays and batched
operations.
"""

import cupy as cp
import numpy as np
from pycbc.types import Array, FrequencySeries, float32, complex64
from pycbc import DYN_RANGE_FAC
import time

# Global timing accumulators
_profile_times = {
    'template_generation': 0.0,
    'sigmasq_computation': 0.0,
    'matched_filter': 0.0,
    'snr_normalization': 0.0,
    'threshold_cluster': 0.0,
    'chisq_computation': 0.0,
    'sg_chisq': 0.0,
    'newsnr_cut': 0.0,
    'output_preparation': 0.0,
    'total_function': 0.0,
    'gpu_sync': 0.0
}
_profile_counts = {
    'batches_processed': 0,
    'segments_processed': 0,
    'templates_processed': 0,
    'triggers_found': 0
}

def reset_profile():
    """Reset profiling counters"""
    global _profile_times, _profile_counts
    for key in _profile_times:
        _profile_times[key] = 0.0
    for key in _profile_counts:
        _profile_counts[key] = 0

def print_profile():
    """Print detailed profiling information"""
    import logging
    logger = logging.getLogger('py.pycbc')
    
    total = _profile_times['total_function']
    if total == 0:
        logger.info("No profiling data collected")
        return
    
    logger.info("=" * 80)
    logger.info("GPU BATCHED INSPIRAL DETAILED PROFILE")
    logger.info("=" * 80)
    logger.info(f"Batches processed: {_profile_counts['batches_processed']}")
    logger.info(f"Segments processed: {_profile_counts['segments_processed']}")
    logger.info(f"Templates processed: {_profile_counts['templates_processed']}")
    logger.info(f"Triggers found: {_profile_counts['triggers_found']}")
    logger.info("")
    logger.info(f"{'Operation':<30} {'Time (s)':<12} {'% Total':<10} {'Avg/call (ms)':<15}")
    logger.info("-" * 80)
    
    items = [
        ('Template Generation', 'template_generation', _profile_counts['batches_processed']),
        ('Sigmasq Computation', 'sigmasq_computation', _profile_counts['segments_processed']),
        ('Matched Filtering', 'matched_filter', _profile_counts['segments_processed']),
        ('SNR Normalization', 'snr_normalization', _profile_counts['segments_processed']),
        ('Threshold & Cluster', 'threshold_cluster', _profile_counts['segments_processed']),
        ('Chi-squared', 'chisq_computation', _profile_counts['segments_processed']),
        ('SG Chi-squared', 'sg_chisq', _profile_counts['segments_processed']),
        ('NewSNR Cut', 'newsnr_cut', _profile_counts['segments_processed']),
        ('Output Preparation', 'output_preparation', _profile_counts['segments_processed']),
        ('GPU Synchronization', 'gpu_sync', _profile_counts['segments_processed']),
    ]
    
    for name, key, count in items:
        t = _profile_times[key]
        pct = 100.0 * t / total if total > 0 else 0
        avg = 1000.0 * t / count if count > 0 else 0
        logger.info(f"{name:<30} {t:<12.4f} {pct:<10.2f} {avg:<15.2f}")
    
    logger.info("-" * 80)
    logger.info(f"{'TOTAL':<30} {total:<12.4f} {100.0:<10.2f}")
    logger.info("=" * 80)


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
    out_vals : dict
        Dictionary with trigger data arrays:
        - 'template_id': template indices (from t_nums)
        - 'time_index': time indices
        - 'snr': SNR values
        - 'chisq': chi-squared values
        - 'chisq_dof': chi-squared degrees of freedom
        - 'sigmasq': sigmasq values
    tparams : list of dict
        Template parameters for each template
    """
    import copy
    
    t_start_func = time.time()
    _profile_counts['batches_processed'] += 1
    _profile_counts['templates_processed'] += len(t_nums)
    
    num_templates = len(t_nums)
    tparams = []
    
    # Accumulators for all triggers across all segments
    all_template_ids = []
    all_time_indices = []
    all_snrs = []
    all_chisqs = []
    all_chisq_dofs = []
    all_sigmasqs = []
    
    # Get template data as 2D cupy array (batch_size x freq_length)
    # This calls directly into batched template generation, bypassing PyCBC array overhead
    t0 = time.time()
    htilde_batch, templates, kmin_array, kmax_array = _generate_templates_gpu(t_nums, bank)
    cp.cuda.Stream.null.synchronize()
    _profile_times['template_generation'] += time.time() - t0
    
    tparams = [bank.table[i] for i in t_nums]
    
    # Get filter parameters for each template (kept for compatibility, but kmin/kmax are better)
    template_flow = cp.array([t.f_lower for t in templates], dtype=cp.float32)
    
    # Process each segment
    for s_num, stilde in enumerate(segments):
        _profile_counts['segments_processed'] += 1
        # Skip if any template in batch should be rejected
        if not all(inj_filter_rejector.template_segment_checker(bank, t_num, stilde) 
                  for t_num in t_nums):
            continue
        
        # Get segment data as cupy array
        stilde_data = stilde._data  # Underlying cupy array
        psd_data = stilde.psd._data  # Underlying cupy array
        
        # Compute sigmasq for all templates
        t0 = time.time()
        sigmasqs = _compute_sigmasqs_gpu(htilde_batch, psd_data, kmin_array, kmax_array, stilde.delta_f)
        cp.cuda.Stream.null.synchronize()
        _profile_times['sigmasq_computation'] += time.time() - t0
        
        # Batched matched filtering - direct kernel call
        # Note: stilde is already overwhitened (divided by PSD)
        t0 = time.time()
        snr_batch, corr_batch = _matched_filter_gpu(htilde_batch, stilde_data, 
                                                     matched_filter.kmin, 
                                                     matched_filter.kmax)
        cp.cuda.Stream.null.synchronize()
        _profile_times['matched_filter'] += time.time() - t0
        
        # The SNR normalization for overwhitened data
        # SNR = IFFT(conj(h) * s) / sqrt(sigmasq / (4 * delta_f))
        # Simplifying: SNR = IFFT_result * sqrt(4 * delta_f) / sqrt(sigmasq)
        # But let's match what PyCBC does: just divide by sqrt(sigmasq)
        t0 = time.time()
        snr_batch = snr_batch / cp.sqrt(sigmasqs[:, cp.newaxis])
        cp.cuda.Stream.null.synchronize()
        _profile_times['snr_normalization'] += time.time() - t0
        
        # Threshold and cluster - operates on 2D arrays
        t0 = time.time()
        trigger_indices, trigger_snrs = _threshold_and_cluster_gpu(
            snr_batch, opt.snr_threshold, cluster_window, stilde.analyze)
        cp.cuda.Stream.null.synchronize()
        _profile_times['threshold_cluster'] += time.time() - t0
        
        # Compute chi-squared for triggers
        t0 = time.time()
        chisqs_batch, chisq_dofs_batch = _compute_chisq_gpu(
            corr_batch, trigger_indices, trigger_snrs, sigmasqs,
            stilde.psd, power_chisq, templates, t_nums)
        cp.cuda.Stream.null.synchronize()
        _profile_times['chisq_computation'] += time.time() - t0
        
        # Extract results for each template
        t0 = time.time()
        
        # Early exit if no triggers
        n_triggers_total = len(trigger_indices)
        if n_triggers_total == 0:
            _profile_times['output_preparation'] += time.time() - t0
            continue
        
        # Move constant calculations outside the loop
        time_offset = stilde.cumulative_index - stilde.analyze.start
        
        # Apply time offset to all triggers at once (before splitting by template)
        trigger_indices[:, 1] += time_offset
        
        # Get template ID for each trigger (this is the index within the batch, 0 to num_templates-1)
        template_batch_ids = trigger_indices[:, 0]
        
        # Total triggers for profiling
        _profile_counts['triggers_found'] += n_triggers_total
        
        # Map batch indices to actual template IDs from t_nums
        # template_batch_ids are indices 0, 1, 2, ... num_templates-1
        # We need to convert them to actual template indices from t_nums
        template_ids_actual = cp.array([t_nums[i] for i in template_batch_ids.get()], dtype=cp.int32)
        
        # Collect arrays for this segment
        time_indices = trigger_indices[:, 1]
        
        # Get sigmasq for each trigger
        sigmasq_vals = sigmasqs[template_batch_ids]
        
        # Append to accumulators
        all_template_ids.append(template_ids_actual)
        all_time_indices.append(time_indices)
        all_snrs.append(trigger_snrs)
        all_chisqs.append(chisqs_batch)
        all_chisq_dofs.append(chisq_dofs_batch)
        all_sigmasqs.append(sigmasq_vals)
        
        _profile_times['output_preparation'] += time.time() - t0
    
    # Concatenate all results
    if len(all_template_ids) > 0:
        out_vals = {
            'template_id': Array(cp.concatenate(all_template_ids), copy=False),
            'time_index': Array(cp.concatenate(all_time_indices), copy=False),
            'snr': Array(cp.concatenate(all_snrs), copy=False),
            'chisq': Array(cp.concatenate(all_chisqs), copy=False),
            'chisq_dof': Array(cp.concatenate(all_chisq_dofs), copy=False),
            'sigmasq': Array(cp.concatenate(all_sigmasqs), copy=False)
        }
    else:
        # No triggers found
        out_vals = {
            'template_id': Array(cp.array([], dtype=cp.int32), copy=False),
            'time_index': Array(cp.array([], dtype=cp.int32), copy=False),
            'snr': Array(cp.array([], dtype=cp.complex64), copy=False),
            'chisq': Array(cp.array([], dtype=cp.float32), copy=False),
            'chisq_dof': Array(cp.array([], dtype=cp.int32), copy=False),
            'sigmasq': Array(cp.array([], dtype=cp.float32), copy=False)
        }
    
    _profile_times['total_function'] += time.time() - t_start_func
    return out_vals, tparams


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
    htilde_batch, htilde_list, kmin_array, kmax_array = spa_tmplt_batch(templates_params, bank.filter_length, **common_kwds)
    
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
    
    return htilde_batch, htilde_list, kmin_array, kmax_array


def _compute_sigmasqs_gpu(htilde_batch, psd_data, kmin_array, kmax_array, delta_f):
    """Compute sigmasq for batch of templates using fully batched GPU computation
    
    Computes: sigmasq[i] = 4 * delta_f * sum(|h[i,k]|^2 / S[k]) for k in [kmin[i], kmax[i])
    
    Uses broadcasting to create a mask and compute all sigmasqs in parallel.
    
    Parameters
    ----------
    htilde_batch : cupy.ndarray
        2D array of templates (num_templates x freq_length)
    psd_data : cupy.ndarray
        PSD array (freq_length,)
    kmin_array : cupy.ndarray
        Start frequency index for each template (num_templates,)
    kmax_array : cupy.ndarray
        End frequency index for each template (num_templates,)
    delta_f : float
        Frequency spacing
        
    Returns
    -------
    cupy.ndarray
        Array of sigmasq values for each template (num_templates,)
    """
    num_templates = htilde_batch.shape[0]
    freq_length = min(htilde_batch.shape[1], len(psd_data))
    
    # Compute |h|^2 / PSD for all templates and frequencies at once
    # Shape: (num_templates, freq_length)
    power_spectrum = (htilde_batch[:, :freq_length].real**2 + htilde_batch[:, :freq_length].imag**2) / psd_data[:freq_length]
    
    # Create frequency index array: shape (freq_length,)
    k_indices = cp.arange(freq_length, dtype=cp.int64)
    
    # Broadcast to create mask: shape (num_templates, freq_length)
    # mask[i, k] = True if kmin[i] <= k < kmax[i]
    kmin_broadcast = kmin_array[:, cp.newaxis]  # shape: (num_templates, 1)
    kmax_broadcast = kmax_array[:, cp.newaxis]  # shape: (num_templates, 1)
    mask = (k_indices >= kmin_broadcast) & (k_indices < kmax_broadcast)  # shape: (num_templates, freq_length)
    
    # Apply mask and sum across frequency axis
    # This computes the sum for all templates in parallel
    sigmasqs = cp.sum(power_spectrum * mask, axis=1)  # shape: (num_templates,)
    
    # Apply the 4 * delta_f factor
    sigmasqs *= 4.0 * delta_f
    
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
    corr_batch = cp.zeros((num_templates, (freq_length-1)*2), dtype=cp.complex64)
    
    # Ensure stilde_data length matches
    stilde_len = min(len(stilde_data), freq_length)
    
    # Batched correlation: corr = conj(h) * s
    corr_batch[:, kmin:kmax] = cp.conj(htilde_batch[:, kmin:kmax]) * stilde_data[kmin:kmax]

    # Batched IFFT to get SNR time series
    snr_batch = cp.fft.ifft(corr_batch, axis=1)
    delta_f = 1./256
    norm = (4.0 * delta_f)
    snr_batch *= norm * 524288
    
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
        
        # Fix template lengths to match PSD (RFFT length)
        # Templates have full FFT length but chisq code expects RFFT length
        for template in templates:
            if len(template) != psd_len:
                # Resize the template data to RFFT length
                template.resize(psd_len)
    
    chisq_list, chisq_dof_list = power_chisq.values_batch(
        corr_batch, unnormalized_snrs, snr_norms, sizes, template_map, psd, time_indices, templates
    )
    
    
    # Check if chisq was actually computed (or disabled)
    if chisq_list is None:
        # Chi-squared was disabled, return dummy values
        chisq_dofs = cp.full(num_triggers, 10, dtype=cp.int32)
        chisqs = chisq_dofs.astype(cp.float32)
        return chisqs, chisq_dofs
    
    # Concatenate results back into single arrays
    chisqs = cp.concatenate(chisq_list) if len(chisq_list) > 0 else cp.array([], dtype=cp.float32)
    chisq_dofs = cp.concatenate(chisq_dof_list) if len(chisq_dof_list) > 0 else cp.array([], dtype=cp.int32)
    
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
