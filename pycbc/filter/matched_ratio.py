import logging
import numpy as np
import mkl_fft
from pycbc.types import zeros, complex64
from pycbc.filter.matchedfilter import matched_filter_core
from pycbc.filter.matchedfilter_cpu import fast_multiply_analytic_cython, find_peaks_in_block_cython

class RatioMatchedFilterControl(object):
    """
    High-performance engine for hierarchical "Ratio/FIR" matched filtering.
    Uses mkl_fft for ALL FFT operations to maximize throughput and consistency.
    """

    def __init__(self, snr_threshold, delta_f, 
                 high_frequency_cutoff=None, fir_fft_length=4096, batch_size=64):
        self.delta_f = delta_f
        self.snr_threshold = snr_threshold
        self.f_high = high_frequency_cutoff
        
        self.threshold_sq = float(snr_threshold**2)
        
        self.fir_fft_len = fir_fft_length
        self.batch_size = batch_size
        
        # 2. Intermediate Buffers (Batch x Block)
        total_batch_size = batch_size * fir_fft_length
        self.temp_freq_mult = zeros(total_batch_size, dtype=complex64)
        self.corr_output_buffer = zeros(total_batch_size, dtype=complex64)

        # 3. Filter Preparation Buffers
        self.filters_padded = zeros(total_batch_size, dtype=complex64)
        self.filters_f_buffer = zeros(total_batch_size, dtype=complex64)
        
        # Direct MKL handle
        self.fft_lib = mkl_fft

    def prepare_filters(self, fir_taps, tap_counts):
        """
        Prepare frequency-domain filters for a batch of taps.
        """
        n_filters, n_taps = fir_taps.shape
        if n_taps >= self.fir_fft_len:
             raise ValueError("FIR Taps (%d) exceed FFT block length (%d)" % 
                              (n_taps, self.fir_fft_len))
        
        # Calculate max tap count for validity logic
        n_taps_max = int(np.max(tap_counts))
        
        filters_f = self._fft_all_filters(fir_taps, tap_counts)
        return filters_f, n_taps_max

    def process_segment(self, stilde, psd, ref_template, filters_f, n_taps, indices, 
                        valid_slice=None):
        """
        Process a single data segment.
        """
        if valid_slice is None:
            valid_slice = getattr(stilde, 'analyze', None)

        # 1. Calculate Reference Normalization
        h_norm = ref_template.sigmasq(psd)

        # 2. Calculate Reference SNR
        snr, _, norm = matched_filter_core(
            ref_template, stilde, psd=psd,
            low_frequency_cutoff=ref_template.f_lower,
            high_frequency_cutoff=self.f_high,
            h_norm=h_norm 
        )
        
        # 3. Execute Blocked Kernel
        local_idxs, t_idxs, snr_vals = self._execute_blocked_kernel(
            snr.numpy() * (norm * stilde.delta_t), filters_f, n_taps, valid_slice
        )
        
        # 4. Map indices
        if len(local_idxs) > 0:
            global_ids = indices[local_idxs]
            return global_ids, t_idxs, snr_vals, h_norm
        else:
            return [], [], []

    def _fft_all_filters(self, taps, counts):
        """Helper to FFT all filters using mkl_fft."""
        n_filters, n_taps_alloc = taps.shape
        filters_f = np.zeros((n_filters, self.fir_fft_len), dtype=np.complex64)
        
        padded_view = self.filters_padded.data
        padded_reshaped = padded_view.reshape(self.batch_size, self.fir_fft_len)

        for start in range(0, n_filters, self.batch_size):
            end = min(start + self.batch_size, n_filters)
            batch_len = end - start
            
            # 1. Zero out and Fill
            padded_reshaped[:batch_len, :] = 0
            
            tmp_taps = taps[start:end]
            padded_reshaped[:batch_len, :n_taps_alloc] = tmp_taps
            
            # 2. Variable Roll Logic
            current_counts = counts[start:end]
            roll_offsets = -(current_counts // 2)
            
            cols = np.arange(self.fir_fft_len)
            rows = np.arange(batch_len)[:, None]
            shifted_cols = (cols[None, :] - roll_offsets[:, None]) % self.fir_fft_len
            
            current_data = padded_reshaped[:batch_len].copy()
            padded_reshaped[:batch_len] = current_data[rows, shifted_cols]

            # 3. Execute FFT (Direct MKL call)
            fft_out = self.fft_lib.fft(padded_reshaped[:batch_len], axis=-1)
            
            # 4. Conjugate & Store
            filters_f[start:end] = np.conj(fft_out)
            
        return filters_f

    def _execute_blocked_kernel(self, data, filters_f, n_taps, valid_slice):
        """
        Inner loop: Time-Blocking + Filter-Batching using mkl_fft.
        """
        n_samples = len(data)
        n_filters = len(filters_f)
        
        N_FFT = self.fir_fft_len
        
        all_f_idxs = []
        all_t_idxs = []
        all_snrs = []
        
        # Reshape views for Cython kernel
        freq_mult_view = self.temp_freq_mult.data.reshape(self.batch_size, N_FFT)
        corr_out_view = self.corr_output_buffer.data.reshape(self.batch_size, N_FFT)

        if valid_slice:
            v_start = valid_slice.start
            v_stop = valid_slice.stop
        else:
            v_start = 0
            v_stop = n_samples
 
        block_f_cache = {}
        
        # --- OUTER LOOP: Time Blocks ---
        for f_start in range(0, n_filters, self.batch_size):  
        
            f_end = min(f_start + self.batch_size, n_filters)
            actual_batch_size = f_end - f_start
            
            # Valid output samples per block (Overlap-Save)
            n_taps_max = n_taps[f_start:f_end].max()
            N_VALID = N_FFT - n_taps_max + 1
            STEP = N_VALID
            bad_start = n_taps_max // 2

            # Determine Loop Bounds
            first_block_idx = (v_start - bad_start) // STEP
            loop_start = first_block_idx * STEP
         
            for t_start in range(loop_start, n_samples, STEP):

                block_valid_t0 = t_start + bad_start
                
                if block_valid_t0 >= v_stop:
                    break
                
                if block_valid_t0 + N_VALID <= v_start:
                    continue

                roi_start = max(v_start, block_valid_t0)
                roi_stop = min(v_stop, block_valid_t0 + N_VALID)
                
                roi_len = roi_stop - roi_start
                
                if roi_len <= 0: 
                    continue

                buf_slice_start = roi_start - t_start

                if t_start not in block_f_cache:
                    t_end = min(t_start + N_FFT, n_samples)
                    block_in_view = np.zeros(self.fir_fft_len, dtype=complex64)
                    block_in_view[0:t_end-t_start] = data[t_start:t_end]
                    block_f_view = self.fft_lib.fft(block_in_view)
                    block_f_cache[t_start] = block_f_view
                
                block_f_view = block_f_cache[t_start]
                filter_batch_f = filters_f[f_start:f_end]
                
                current_mult_view = freq_mult_view[:actual_batch_size]
                current_corr_view = corr_out_view[:actual_batch_size]
                
                fast_multiply_analytic_cython(
                    block_f_view, filter_batch_f, current_mult_view
                )
                
                self.fft_lib.ifft(
                    current_mult_view, 
                    axis=-1, 
                    out=current_corr_view
                )
                
                f_list, t_list, s_list = find_peaks_in_block_cython(
                    current_corr_view, 
                    roi_start,          
                    roi_len,            
                    self.threshold_sq, 
                    f_start,
                    input_offset=buf_slice_start
                )

                if f_list:
                    all_f_idxs.extend(f_list)
                    all_t_idxs.extend(t_list)
                    all_snrs.extend(s_list)

        return (np.array(all_f_idxs, dtype=np.int32), 
                np.array(all_t_idxs, dtype=np.int64), 
                np.array(all_snrs, dtype=np.complex64))
