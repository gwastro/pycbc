import logging
import numpy as np
from pycbc.types import zeros, complex64
from pycbc.filter.matchedfilter import matched_filter_core
from pycbc.filter.matchedfilter_cpu import fast_multiply_analytic_cython, find_peaks_in_block_cython
import pycbc.fft

class RatioMatchedFilterControl(object):
    """
    High-performance engine for hierarchical "Ratio/FIR" matched filtering.
    
    This class owns the memory buffers and execution logic to:
    1. Calculate the SNR of a Reference Template against stored data segments.
    2. Batch-process many 'Target' templates by convolving that Reference SNR
       with provided FIR filters using AVX-optimized kernels.
    """

    def __init__(self, snr_threshold, tlen, delta_f, 
                 high_frequency_cutoff=None, fir_fft_length=4096, batch_size=64):
        """
        Initialize the controller with PyCBC FFT plans and memory buffers.
        
        Parameters
        ----------
        snr_threshold : float
            Threshold for recording triggers.
        tlen : int
            Length of the analysis segment in samples (time domain).
        delta_f : float
            Frequency resolution.
        high_frequency_cutoff : float, optional
            Upper frequency limit for the reference template calculation.
            If None, the full valid bandwidth (up to Nyquist) is used.
        fir_fft_length : int, optional
            The block size for the cache-blocked FIR engine (default: 4096).
        batch_size : int, optional
            Number of fine templates to process simultaneously.
        """
        self.tlen = tlen
        self.delta_f = delta_f
        self.delta_t = 1.0 / (tlen * delta_f)
        self.snr_threshold = snr_threshold
        self.f_high = high_frequency_cutoff
        
        self.threshold_sq = float(snr_threshold**2)
        
        # --- Buffers ---
        self.ref_snr_mem = zeros(tlen, dtype=complex64)
        
        self.fir_fft_len = fir_fft_length
        self.batch_size = batch_size
        
        # 1. Input Block & Plan
        self.block_data_in = zeros(fir_fft_length, dtype=complex64)
        self.block_data_f = zeros(fir_fft_length, dtype=complex64)
        self.block_fft_plan = pycbc.fft.FFT(self.block_data_in, self.block_data_f)
        
        # 2. Intermediate Buffers
        total_batch_size = batch_size * fir_fft_length
        self.temp_freq_mult = zeros(total_batch_size, dtype=complex64)
        self.corr_output_buffer = zeros(total_batch_size, dtype=complex64)

        self.batch_ifft_plan = pycbc.fft.IFFT(
            self.temp_freq_mult, self.corr_output_buffer, nbatch=batch_size
        )

        # 3. Prep Buffers
        self.filters_padded = zeros(total_batch_size, dtype=complex64)
        self.filters_f_buffer = zeros(total_batch_size, dtype=complex64)
        self.filter_batch_fft_plan = pycbc.fft.FFT(
            self.filters_padded, self.filters_f_buffer, nbatch=batch_size
        )

        logging.info("RatioControl Initialized: BlockSize=%d, BatchSize=%d", 
                     fir_fft_length, batch_size)

    def prepare_filters(self, fir_taps):
        """
        Prepare frequency-domain filters for a batch of taps.
        
        Parameters
        ----------
        fir_taps : numpy.ndarray
            2D array of FIR taps (shape: [N_templates, N_taps]).

        Returns
        -------
        filters_f : numpy.ndarray (Complex64)
            Frequency domain filters ready for processing.
        n_taps : int
            The tap count used (passed through for process_segment).
        """
        n_filters, n_taps = fir_taps.shape
        if n_taps >= self.fir_fft_len:
             raise ValueError("FIR Taps (%d) exceed FFT block length (%d)" % 
                              (n_taps, self.fir_fft_len))
        
        filters_f = self._fft_all_filters(fir_taps)
        return filters_f, n_taps

    def process_segment(self, stilde, psd, ref_template, filters_f, n_taps, 
                        valid_slice=None):
        """
        Process a single data segment against a prepared batch of filters.

        Parameters
        ----------
        stilde : FrequencySeries
            The frequency-domain strain data segment (already overwhitened).
        psd : FrequencySeries
            The PSD used to normalize the reference SNR.
        ref_template : FrequencySeries
            The coarse reference waveform.
        filters_f : numpy.ndarray
            Pre-calculated frequency-domain filters (from `prepare_filters`).
        n_taps : int
            Number of taps in the FIR filters.
        valid_slice : slice, optional
            The slice of the time-series that contains valid analysis data.
            If None, `stilde.analyze` is used.

        Returns
        -------
        triggers : tuple
            (local_filter_indices, time_indices, snr_values)
            Note: local_filter_indices are 0..N-1 relative to the input filters_f.
        """
        if valid_slice is None:
            valid_slice = getattr(stilde, 'analyze', None)

        # 1. Calculate Reference Normalization
        h_norm = ref_template.sigmasq(psd)

        # 2. Calculate Reference SNR
        # We pass h_norm explicitly to avoid re-whitening.
        snr, _, norm = matched_filter_core(
            ref_template, stilde, 
            low_frequency_cutoff=ref_template.f_lower,
            high_frequency_cutoff=self.f_high,
            out=self.ref_snr_mem,
            h_norm=h_norm 
        )
        snr *= norm
        
        if hasattr(snr, 'numpy'):
            ref_data = snr.numpy()
        else:
            ref_data = snr
        
        # 3. Execute Blocked Kernel
        # This returns raw indices relative to the filter batch (0..N)
        local_idxs, t_idxs, snr_vals = self._execute_blocked_kernel(
            ref_data, filters_f, n_taps, valid_slice
        )
        
        return local_idxs, t_idxs, snr_vals

    def _fft_all_filters(self, taps):
        """Helper to FFT all filters (Includes Roll and Conjugate)."""
        n_filters, n_taps = taps.shape
        filters_f = np.zeros((n_filters, self.fir_fft_len), dtype=np.complex64)
        
        padded_view = self.filters_padded.data
        fft_view = self.filters_f_buffer.data
        
        padded_reshaped = padded_view.reshape(self.batch_size, self.fir_fft_len)
        fft_reshaped = fft_view.reshape(self.batch_size, self.fir_fft_len)

        roll_amt = -(n_taps // 2)

        for start in range(0, n_filters, self.batch_size):
            end = min(start + self.batch_size, n_filters)
            batch_len = end - start
            
            padded_reshaped[:batch_len, :] = 0
            
            # Copy and Roll
            tmp_taps = taps[start:end]
            padded_reshaped[:batch_len, :n_taps] = tmp_taps
            
            # Apply Roll to center impulse at index 0
            padded_rolled = np.roll(padded_reshaped[:batch_len], roll_amt, axis=1)
            padded_reshaped[:batch_len] = padded_rolled

            self.filter_batch_fft_plan.execute()
            
            # Conjugate (Correlation = A * conj(B))
            filters_f[start:end] = np.conj(fft_reshaped[:batch_len])
            
        return filters_f

    def _execute_blocked_kernel(self, data, filters_f, n_taps, valid_slice):
        """
        Inner loop: Time-Blocking + Filter-Batching.
        Aligned perfectly to valid_slice to avoid post-filtering.
        """
        n_samples = len(data)
        n_filters = len(filters_f)
        
        N_FFT = self.fir_fft_len
        # Valid output samples per block (Overlap-Save)
        # Note: We lose (n_taps - 1) samples due to circular convolution corruption
        N_VALID = N_FFT - n_taps + 1
        STEP = N_VALID
        
        # Calculate indices to drop from block output due to filter roll
        bad_start = n_taps // 2
        
        all_f_idxs = []
        all_t_idxs = []
        all_snrs = []
        
        freq_mult_view = self.temp_freq_mult.data.reshape(self.batch_size, N_FFT)
        corr_out_view = self.corr_output_buffer.data.reshape(self.batch_size, N_FFT)
        
        if valid_slice:
            v_start = valid_slice.start
            v_stop = valid_slice.stop
        else:
            v_start = 0
            v_stop = n_samples

        # Determine Loop Bounds
        # We align the loop such that 'bad_start' (first valid point of block) 
        # aligns with the STEP grid covering v_start.
        
        first_block_idx = (v_start - bad_start) // STEP
        loop_start = first_block_idx * STEP
        
        # --- OUTER LOOP: Time Blocks ---
        for t_start in range(loop_start, n_samples, STEP):
            
            # The time index corresponding to the first valid sample in this block
            block_valid_t0 = t_start + bad_start
            
            # If this block's valid region starts after our analysis window, we are done.
            if block_valid_t0 >= v_stop:
                break
                
            # If this block ends before our analysis window, skip it.
            if block_valid_t0 + N_VALID <= v_start:
                continue

            # Calculate intersection of this block's valid region and the requested valid slice
            roi_start = max(v_start, block_valid_t0)
            roi_stop = min(v_stop, block_valid_t0 + N_VALID)
            
            roi_len = roi_stop - roi_start
            
            if roi_len <= 0: 
                continue

            # Determine buffer slice indices 
            # Time T corresponds to buffer index (T - t_start)
            buf_slice_start = roi_start - t_start
            buf_slice_stop = roi_stop - t_start
            
            # Load Data
            t_end = min(t_start + N_FFT, n_samples)
            self.block_data_in.clear()
            self.block_data_in[0:t_end-t_start] = data[t_start:t_end]
            self.block_fft_plan.execute()

            # --- INNER LOOP: Filter Batches ---
            for f_start in range(0, n_filters, self.batch_size):
                f_end = min(f_start + self.batch_size, n_filters)
                actual_batch_size = f_end - f_start
                
                filter_batch_f = filters_f[f_start:f_end]
                current_mult_view = freq_mult_view[:actual_batch_size]
                current_corr_view = corr_out_view[:actual_batch_size]
                
                fast_multiply_analytic_cython(
                    self.block_data_f.data, filter_batch_f, current_mult_view
                )
                
                self.batch_ifft_plan.execute()
                
                # Sliced Peak Finding
                # We slice the buffer exactly to the Region of Interest
                valid_window_view = current_corr_view[:, buf_slice_start : buf_slice_stop]
                
                f_list, t_list, s_list = find_peaks_in_block_cython(
                    valid_window_view,
                    roi_start,          # Time index corresponding to view[0]
                    roi_len,            # Number of valid points to check
                    self.threshold_sq, 
                    0, 
                    f_start 
                )
                
                if f_list:
                    all_f_idxs.extend(f_list)
                    all_t_idxs.extend(t_list)
                    all_snrs.extend(s_list)

        return (np.array(all_f_idxs, dtype=np.int32), 
                np.array(all_t_idxs, dtype=np.int64), 
                np.array(all_snrs, dtype=np.float32))
