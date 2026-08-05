import logging
import numpy as np
from pycbc.fft import FFT, IFFT
from pycbc.types import zeros, complex64
from pycbc.filter.matchedfilter import matched_filter_core
from pycbc.filter.matchedfilter_cpu import fast_multiply_analytic_cython, find_peaks_in_block_cython

class RatioMatchedFilterControl(object):
    """
    High-performance engine for hierarchical "Ratio/FIR" matched filtering.
    Uses pycbc's class-based FFT interface (pycbc.fft.FFT/IFFT) for all FFT
    operations, so the actual backend (MKL, FFTW, numpy, ...) follows
    whatever --processing-scheme the caller selected, with FFT/IFFT plans
    built once per batch size and reused.

    Note: pycbc's class-based IFFT does not normalize by 1/N the way
    numpy.fft.ifft does; see _forward_block_fft for where that's corrected.
    """

    def __init__(self, snr_threshold, delta_f,
                 high_frequency_cutoff=None, fir_fft_length=4096, batch_size=64, tap_sample_rate=2048, engine_sample_rate=2048):
        self.delta_f = delta_f
        self.snr_threshold = snr_threshold
        self.f_high = high_frequency_cutoff

        self.tap_sr = int(tap_sample_rate)
        self.engine_sr = int(engine_sample_rate)

        self.threshold_sq = float(snr_threshold**2)

        self.fir_fft_len = fir_fft_length
        self.batch_size = batch_size

        # The bank's taps are generated at tap_sample_rate but filtered
        # against data at engine_sample_rate; the ratio must be an exact
        # integer so the high-resolution FFT below preserves delta_f.
        exact_ratio = self.tap_sr / self.engine_sr
        self.decimation_factor = int(np.round(exact_ratio))
        if abs(exact_ratio - self.decimation_factor) > 1e-5 or self.decimation_factor < 1:
            raise ValueError(
                f"Multi-rate Error: The bank sample rate ({self.tap_sr} Hz) must be "
                f"an exact integer multiple of the engine sample "
                f"rate ({self.engine_sr} Hz).\n"
                f"Calculated ratio was {exact_ratio:.4f}. Please use standard power-of-2 "
                f"downsampling scales (e.g., 2048/512)."
            )
        # 512/4096 = 0.125   2048/(4*4096) = 0.125 preserving delta_f
        # 4096/4096 = 1      2048/(1/2*4096) = 1 for 4096 engine 2048 bank
        self.high_res_fft_len = self.fir_fft_len * self.decimation_factor

        # FFT/IFFT plans are built lazily, keyed by the (nbatch, size) they
        # were requested with, and reused for the life of this engine: the
        # dominant batch size is self.batch_size, but the last filter-batch
        # of any coarse-template group is usually smaller (n_filters is not
        # generally a multiple of batch_size), so a handful of extra plan
        # sizes get cached too.
        self._fft_plans = {}
        self._ifft_plans = {}

    def _get_plan(self, plans, cls, nbatch, size):
        key = (nbatch, size)
        cached = plans.get(key)
        if cached is None:
            in_arr = zeros(nbatch * size, dtype=complex64)
            out_arr = zeros(nbatch * size, dtype=complex64)
            if nbatch > 1:
                plan = cls(in_arr, out_arr, nbatch=nbatch, size=size)
            else:
                plan = cls(in_arr, out_arr)
            in_view = in_arr.data.reshape(nbatch, size)
            out_view = out_arr.data.reshape(nbatch, size)
            cached = plans[key] = (plan, in_view, out_view)
        return cached

    def _get_fft_plan(self, nbatch, size):
        return self._get_plan(self._fft_plans, FFT, nbatch, size)

    def _get_ifft_plan(self, nbatch, size):
        return self._get_plan(self._ifft_plans, IFFT, nbatch, size)

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

        decimate = int(np.round(self.tap_sr / self.engine_sr))
        self.ref_snr = snr.numpy() * (norm * stilde.delta_t)  / decimate

        # 3. Execute Blocked Kernel
        local_idxs, t_idxs, snr_vals, tstarts = self._execute_blocked_kernel(
            self.ref_snr, filters_f, n_taps, valid_slice
        )

        # 4. Map indices
        if len(local_idxs) > 0:
            global_ids = indices[local_idxs]
            return global_ids, t_idxs, snr_vals, tstarts, h_norm
        else:
            return [], [], [], tstarts, h_norm

    def _fft_all_filters(self, taps, counts):
        """Helper to FFT all filters, batch_size rows at a time."""
        n_filters, n_taps_alloc = taps.shape
        filters_f = np.zeros((n_filters, self.fir_fft_len), dtype=np.complex64)

        high_res_fft_len = self.high_res_fft_len

        for start in range(0, n_filters, self.batch_size):
            end = min(start + self.batch_size, n_filters)
            batch_len = end - start

            plan, in_view, out_view = self._get_fft_plan(batch_len, high_res_fft_len)

            # Zero out processing buffer, then copy raw 2048 Hz taps into
            # the start of the buffer.
            in_view[:] = 0.0
            tmp_taps = taps[start:end]
            in_view[:, :n_taps_alloc] = tmp_taps

            # Handle Variable Time-Domain Roll Logic at the native 2048 Hz rate
            current_counts = counts[start:end]
            roll_offsets = -(current_counts // 2)

            cols_high = np.arange(high_res_fft_len)
            rows = np.arange(batch_len)[:, None]
            shifted_cols_high = (cols_high[None, :] - roll_offsets[:, None]) % high_res_fft_len

            current_data = in_view.copy()
            in_view[:] = current_data[rows, shifted_cols_high]

            # Transform to Frequency Domain at native resolution
            plan.execute()

            # Brick-Wall Frequency Slicing (Anti-Aliasing & Decimation Match)
            # Because the data engine goes up to 256Hz (the 512Hz Nyquist limit), only need the first 4096 bins of that spectrum
            fft_sliced = out_view[:, :self.fir_fft_len]

            # Conjugate & Store back into the 512 Hz buffer block
            filters_f[start:end] = np.conj(fft_sliced)

        return filters_f

    def _forward_block_fft(self, segment):
        """FFT one time-domain block (nbatch=1), pre-dividing by fir_fft_len.

        pycbc's class-based IFFT (unlike numpy.fft.ifft/mkl_fft.ifft) does
        not normalize its output by 1/N. Since this block spectrum only
        ever gets multiplied by a filter spectrum and then inverse
        transformed (in _execute_blocked_kernel), pre-dividing it here by N
        once -- rather than dividing the batched product after every IFFT --
        gives an already-correctly-normalized correlation for free:
        ifft_raw(block_f/N * filter_f) == ifft_normalized(block_f * filter_f).
        """
        plan, in_view, out_view = self._get_fft_plan(1, self.fir_fft_len)
        in_view[0, :] = 0.0
        in_view[0, :len(segment)] = segment
        plan.execute()
        return out_view[0].copy() / self.fir_fft_len

    def _execute_blocked_kernel(self, data, filters_f, n_taps, valid_slice):
        """
        Inner loop: Time-Blocking + Filter-Batching.
        """
        tap_groups = 3
        nsizes = np.quantile(n_taps, np.linspace(0, 1, tap_groups+1)[1:]).astype(int)
        n_samples = len(data)
        n_filters = len(filters_f)

        N_FFT = self.fir_fft_len

        all_f_idxs = []
        all_t_idxs = []
        all_snrs = []
        all_tstarts = []

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

            ifft_plan, current_mult_view, current_corr_view = self._get_ifft_plan(
                actual_batch_size, N_FFT
            )

            # Valid output samples per block (Overlap-Save)
            n_taps_max = n_taps[f_start:f_end].max()
            i = np.searchsorted(nsizes, n_taps_max)
            n_taps_max = nsizes[i]

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

                t_end = min(t_start + N_FFT, n_samples)
                if t_start not in block_f_cache:
                    block_f_cache[t_start] = self._forward_block_fft(data[t_start:t_end])

                block_f_view = block_f_cache[t_start]
                filter_batch_f = filters_f[f_start:f_end]

                fast_multiply_analytic_cython(
                    block_f_view, filter_batch_f, current_mult_view
                )

                ifft_plan.execute()
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
                    all_tstarts.extend([t_start] * len(s_list))

        return (np.array(all_f_idxs, dtype=np.int32),
                np.array(all_t_idxs, dtype=np.int64),
                np.array(all_snrs, dtype=np.complex64),
                np.array(all_tstarts, dtype=np.int32))
