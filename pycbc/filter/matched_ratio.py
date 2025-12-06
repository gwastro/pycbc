import logging
import numpy as np
from pycbc.types import TimeSeries, FrequencySeries, zeros, complex64
from pycbc.filter.matchedfilter import matched_filter_core
# Import the kernels we just merged into matchedfilter.pyx
from pycbc.filter.matchedfilter import fast_multiply_analytic_cython, find_peaks_in_block_cython
import pycbc.fft

class RatioMatchedFilterControl(object):
    """
    Controls the hierarchical execution of matched filtering using Ratio/FIR
    decomposition.
    
    This class manages the memory and execution logic to:
    1. Compute the SNR of a 'Coarse' reference template using standard FFT convolution.
    2. Batch-process associated 'Fine' templates by convolving the Reference SNR
       with short FIR filters (Ratio Filters) using cache-optimized kernels.
    """

    def __init__(self, high_frequency_cutoff, snr_threshold, tlen, delta_f,
                 segments, ratio_bank, 
                 fir_fft_length=4096, batch_size=64):
        """
        Parameters
        ----------
        high_frequency_cutoff : float
            Upper frequency limit for the reference template generation.
        snr_threshold : float
            Threshold for recording triggers.
        tlen : int
            Length of the analysis segment in samples (time domain).
        delta_f : float
            Frequency resolution.
        segments : list of FrequencySeries
            The data segments (strain) to be analyzed.
        ratio_bank : RatioFilterBank
            The hierarchical bank object.
        fir_fft_length : int, optional
            The block size for the cache-blocked FIR engine (default: 4096).
            This defines the size of the L1/L2 cache-resident operations.
        batch_size : int, optional
            Number of fine templates to process simultaneously in the Cython kernel.
        """
        self.tlen = tlen
        self.delta_f = delta_f
        self.delta_t = 1.0 / (tlen * delta_f)
        self.snr_threshold = snr_threshold
        self.f_high = high_frequency_cutoff
        self.segments = segments
        self.bank = ratio_bank
        
        # Threshold squared for the Cython kernel (avoid sqrts in tight loops)
        self.threshold_sq = float(snr_threshold**2)
        
        # --- Memory Allocation for Stage 1 (Reference SNR) ---
        # We use standard PyCBC buffers for the reference matched filter
        self.ref_snr_mem = zeros(tlen, dtype=complex64)
        
        # --- Memory Allocation for Stage 2 (FIR Engine) ---
        self.fir_fft_len = fir_fft_length
        self.batch_size = batch_size
        
        # 1. Input Block (Time Domain -> Freq Domain)
        # Used to hold a chunk of the Reference SNR time series
        self.block_data_in = np.zeros(fir_fft_length, dtype=np.complex64)
        self.block_data_f = np.zeros(fir_fft_length, dtype=np.complex64)
        
        # 2. Intermediate Frequency Multiplication (Batch x Block)
        # Used to hold the result of RefBlock(f) * FilterBatch(f)
        self.temp_freq_mult = np.zeros((batch_size, fir_fft_length), dtype=np.complex64)
        
        # 3. Output Time Domain Correlation Blocks
        # Used to hold the IFFT result before peak finding
        self.corr_output_buffer = np.zeros((batch_size, fir_fft_length), dtype=np.complex64)

        # 4. Pre-allocated filter buffer (for FFTing the taps)
        self.filters_padded = np.zeros((batch_size, fir_fft_length), dtype=np.complex64)
        
        # We need a standard FFT engine for the small blocks. 
        # numpy.fft is sufficient here as the block size is small (4096) and 
        # usually fits in cache, plus it links to MKL in standard PyCBC installs.
        self.fft_lib = np.fft

        logging.info("RatioControl Initialized: BlockSize=%d, BatchSize=%d", 
                     fir_fft_length, batch_size)

    def process_coarse_group(self, coarse_index):
        """
        Process an entire group: One Coarse Reference + Many Fine Targets.
        
        Parameters
        ----------
        coarse_index : int
            The index of the reference template in the bank.

        Yields
        -------
        segment_index : int
            The index of the data segment currently being processed.
        triggers : tuple
            (template_indices, time_indices, snr_values)
        """
        
        # 1. --- STAGE 1: Reference Calculation ---
        # Generate the Reference Waveform
        h_ref = self.bank.get_coarse_template(coarse_index)
        
        # 2. --- PREP STAGE: Get FIR Taps & Prepare Filters ---
        # We fetch the taps once per group.
        taps, tap_counts, fine_indices = self.bank.get_firs(coarse_index)
        n_filters, n_taps_max = taps.shape
        
        if n_taps_max >= self.fir_fft_len:
             raise ValueError("FIR Taps (%d) exceed FFT block length (%d)" % 
                              (n_taps_max, self.fir_fft_len))

        # FFT the filters for this group
        # This converts time-domain FIR taps into frequency-domain filters 
        # ready for the block convolution engine.
        filters_f = self._fft_all_filters(taps)

        # 3. --- EXECUTION LOOP (Per Segment) ---
        for i, stilde in enumerate(self.segments):
            # Calculate Reference SNR (Standard Matched Filter)
            # This fills self.ref_snr_mem with the complex SNR time series
            snr, _, norm = matched_filter_core(
                h_ref, stilde, 
                low_frequency_cutoff=h_ref.f_lower,
                high_frequency_cutoff=self.f_high,
                out=self.ref_snr_mem # Reuse memory
            )
            
            # Normalize Reference SNR immediately.
            # The Ratio method relies on linearity: Ratio * (Data*Ref/Norm) = Target_SNR
            snr *= norm
            
            # Convert to numpy for Cython consumption (zero-copy if possible)
            if hasattr(snr, 'numpy'):
                ref_data = snr.numpy()
            else:
                ref_data = snr
            
            # Execute the "Dechirper" kernel
            local_idxs, t_idxs, snr_vals = self._execute_blocked_kernel(
                ref_data, filters_f, n_taps_max
            )
            
            if len(local_idxs) > 0:
                # Map local batch indices (0..N) back to global bank indices
                global_idxs = fine_indices[local_idxs]
                yield i, (global_idxs, t_idxs, snr_vals)
            else:
                yield i, ([], [], [])

    def _fft_all_filters(self, taps):
        """Helper to FFT all filters in the group at once (or in chunks)."""
        n_filters, n_taps = taps.shape
        
        # Allocate buffer for this group's freq-domain filters
        filters_f = np.zeros((n_filters, self.fir_fft_len), dtype=np.complex64)
        
        # Process in batches to keep cache usage clean
        for start in range(0, n_filters, self.batch_size):
            end = min(start + self.batch_size, n_filters)
            batch_len = end - start
            
            # Clear padding buffer
            self.filters_padded[:batch_len, :] = 0
            # Copy taps into buffer
            self.filters_padded[:batch_len, :n_taps] = taps[start:end]
            
            # FFT (Complex Conjugate for Correlation!)
            # Standard correlation in freq domain is A(f) * conj(B(f)).
            # Since fast_multiply_analytic_cython does direct multiplication, 
            # we pre-conjugate the filter here.
            filters_f[start:end] = np.conj(self.fft_lib.fft(self.filters_padded[:batch_len], axis=-1))
            
        return filters_f

    def _execute_blocked_kernel(self, data, filters_f, n_taps):
        """
        The inner loop logic: Time-Blocking + Filter-Batching.
        
        This mimics the structure of `find_triggers_stateless` from the
        optimization notebook.
        """
        n_samples = len(data)
        n_filters = len(filters_f)
        
        # Constants for overlap-save method
        N_FFT = self.fir_fft_len
        # Valid output samples per block = Length - Filter + 1
        N_VALID = N_FFT - n_taps + 1
        STEP = N_VALID
        
        all_f_idxs = []
        all_t_idxs = []
        all_snrs = []
        
        # --- OUTER LOOP: Time Blocks (Cache Blocking) ---
        # We process one chunk of Reference SNR against ALL filters before moving on.
        # This keeps the Data chunk hot in L2 cache.
        for t_start in range(0, n_samples, STEP):
            t_end = min(t_start + N_FFT, n_samples)
            
            # 1. Prepare Data Block
            self.block_data_in[:] = 0j
            self.block_data_in[:t_end-t_start] = data[t_start:t_end]
            
            # 2. FFT Data Block
            # Write directly into pre-allocated buffer
            np.copyto(self.block_data_f, self.fft_lib.fft(self.block_data_in))
            
            # Validity check: if near end of file, valid region might shrink
            current_n_valid = min(N_VALID, n_samples - t_start - n_taps + 1)
            if current_n_valid <= 0:
                break

            # --- INNER LOOP: Filter Batches ---
            # Process batches of filters against this single data block
            for f_start in range(0, n_filters, self.batch_size):
                f_end = min(f_start + self.batch_size, n_filters)
                actual_batch_size = f_end - f_start
                
                # Create views into pre-allocated buffers
                filter_batch_f = filters_f[f_start:f_end]
                temp_freq_batch = self.temp_freq_mult[:actual_batch_size]
                corr_block_batch = self.corr_output_buffer[:actual_batch_size]
                
                # 3. Optimized Complex Multiply (Cython AVX)
                # Multiplies Block_Data(f) * Filter_Batch(f)
                fast_multiply_analytic_cython(self.block_data_f, filter_batch_f, temp_freq_batch)
                
                # 4. Batch IFFT
                # Transforms back to time domain
                self.fft_lib.ifft(temp_freq_batch, axis=-1, out=corr_block_batch)
                
                # 5. Fused Peak Finding (Cython)
                # Scans the time-domain block for peaks > threshold_sq.
                # Returns lists of indices/values.
                # f_start passed as offset so returned indices are relative to full filter list.
                f_list, t_list, s_list = find_peaks_in_block_cython(
                    corr_block_batch, 
                    t_start, 
                    current_n_valid, 
                    self.threshold_sq, 
                    f_start 
                )
                
                if f_list:
                    all_f_idxs.extend(f_list)
                    all_t_idxs.extend(t_list)
                    all_snrs.extend(s_list)

        # Final conversion to numpy arrays
        return (np.array(all_f_idxs, dtype=np.int32), 
                np.array(all_t_idxs, dtype=np.int64), 
                np.array(all_snrs, dtype=np.float32))
