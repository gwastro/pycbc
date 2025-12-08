# Copyright (C) 2018  Alex Nitz, Josh Willis
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
# cython: embedsignature=True
import numpy
from .matchedfilter import _BaseCorrelator
cimport numpy, cython
from cython.parallel import prange
from libc.math cimport sqrt

# --- Typedefs for Ratio Filter Kernels ---
ctypedef numpy.complex64_t complex64_t
ctypedef numpy.float32_t float32_t
ctypedef numpy.int32_t int32_t
ctypedef numpy.int64_t int64_t

ctypedef fused COMPLEXTYPE:
    float complex
    double complex

@cython.boundscheck(False)
@cython.wraparound(False)
def _batch_correlate(numpy.ndarray [long, ndim=1] x,
                     numpy.ndarray [float complex, ndim=1] y,
                     numpy.ndarray [long, ndim=1] z,
                     size, num_vectors):
    cdef unsigned int nvec = num_vectors
    cdef unsigned int vsize = size

    cdef float complex* xp
    cdef float complex* zp

    cdef unsigned int i, j

    for i in prange(nvec, nogil=True):
        xp = <float complex*> x[i]
        zp = <float complex*> z[i]
        for j in range(vsize):
            zp[j] = xp[j].conjugate() * y[j]

def batch_correlate_execute(self, y):
    num_vectors = self.num_vectors # pylint:disable=unused-variable
    size = self.size # pylint:disable=unused-variable
    _batch_correlate(self.x.data, y.data, self.z.data, size, num_vectors)

def correlate_numpy(x, y, z):
    z.data[:] = numpy.conjugate(x.data)[:]
    z *= y

@cython.boundscheck(False)
@cython.wraparound(False)
def _correlate(COMPLEXTYPE[:] x,
               COMPLEXTYPE[:] y,
               COMPLEXTYPE[:] z):
    cdef unsigned int xmax = x.shape[0]
    cdef unsigned int i
    for i in prange(xmax, nogil=True):
        z[i] = x[i].conjugate() * y[i]

def correlate(x, y, z):
    _correlate(x.data, y.data, z.data)

class CPUCorrelator(_BaseCorrelator):
    def __init__(self, x, y, z):
        self.x = numpy.array(x.data, copy=False)
        self.y = numpy.array(y.data, copy=False)
        self.z = numpy.array(z.data, copy=False)

    def correlate(self):
        _correlate(self.x, self.y, self.z)

def _correlate_factory(x, y, z):
    return CPUCorrelator

# -----------------------------------------------------------------------------
# Ratio Filter Optimization Kernels
# -----------------------------------------------------------------------------

@cython.boundscheck(False) 
@cython.wraparound(False) 
@cython.cdivision(True)   
def fast_multiply_analytic_cython(
    numpy.ndarray[complex64_t, ndim=1, mode="c"] data_f,
    numpy.ndarray[complex64_t, ndim=2, mode="c"] filter_batch_f,
    numpy.ndarray[complex64_t, ndim=2, mode="c"] out_batch
):
    """
    Cython version of the "half-only" analytic signal multiply.
    
    This kernel is single-threaded and relies on the C compiler's
    autovectorizer (enabled by -march=native) to use AVX.
    """
    
    # --- C-level variable declarations ---
    cdef long batch_size = filter_batch_f.shape[0]
    cdef long n_fft = filter_batch_f.shape[1]
    
    # We only compute up to the Nyquist bin
    cdef long n_half_plus_one = (n_fft // 2) + 1 
    
    cdef long i, j # Loop iterators

    # This is a pure C-loop, no Python overhead.
    for i in range(batch_size):
        for j in range(n_half_plus_one):
            # Direct C-level complex multiplication
            out_batch[i, j] = data_f[j] * filter_batch_f[i, j]

# ... (Previous imports and fast_multiply remain the same) ...

@cython.boundscheck(False) 
@cython.wraparound(False) 
@cython.cdivision(True)   
def find_peaks_in_block_cython(
    numpy.ndarray[complex64_t, ndim=2, mode="c"] corr_output,
    long t_start,
    long n_valid,
    float threshold_sq,
    long f_start_offset,
    long input_offset=0
):
    """
    Cython version of the manually vectorized "max-reduction" kernel.
    Returns three separate, flat lists (f_idx, t_idx, snr)
    """
    
    cdef long n_filters_in_batch = corr_output.shape[0]
    cdef list f_idx_list = []
    cdef list t_idx_list = []
    cdef list snr_list = []
    cdef int VEC_WIDTH = 8
    
    cdef float32_t current_max_snr_sq_vec[8]
    cdef int64_t current_max_idx_vec[8]
    cdef complex64_t current_max_z_vec[8]
    
    cdef long f_batch_idx, i, idx, read_idx
    cdef int v_lane
    cdef int32_t f_global_idx
    cdef complex64_t z
    cdef float32_t mag_sq
    cdef float32_t final_max_snr_sq
    cdef int64_t final_max_idx
    
    for f_batch_idx in range(n_filters_in_batch):
        f_global_idx = <int32_t>(f_start_offset + f_batch_idx)

        # --- Initialize C stack arrays ---
        for v_lane in range(VEC_WIDTH):
            current_max_snr_sq_vec[v_lane] = threshold_sq
            current_max_idx_vec[v_lane] = -1

        # --- Main vectorized loop ---
        for i in range(n_valid // VEC_WIDTH):
            for v_lane in range(VEC_WIDTH):
                idx = i * VEC_WIDTH + v_lane
                
                # Apply Offset Here
                read_idx = idx + input_offset
                
                z = corr_output[f_batch_idx, read_idx]
                mag_sq = z.real * z.real + z.imag * z.imag

                if mag_sq > current_max_snr_sq_vec[v_lane]:
                    current_max_snr_sq_vec[v_lane] = mag_sq
                    # Return Global Time Index (t_start corresponds to idx=0)
                    current_max_idx_vec[v_lane] = t_start + idx
                    current_max_z_vec[v_lane] = z

        # --- Epilogue ---
        for i in range((n_valid // VEC_WIDTH) * VEC_WIDTH, n_valid):
            read_idx = i + input_offset
            z = corr_output[f_batch_idx, read_idx]
            mag_sq = z.real * z.real + z.imag * z.imag
            
            if mag_sq > current_max_snr_sq_vec[0]:
                current_max_snr_sq_vec[0] = mag_sq
                current_max_idx_vec[0] = t_start + i
                current_max_z_vec[0] = z
        
        # --- Final Reduction ---
        final_max_snr_sq = threshold_sq
        final_max_idx = -1
        final_max_z = 0 + 0j
        for v_lane in range(VEC_WIDTH):
            if current_max_snr_sq_vec[v_lane] > final_max_snr_sq:
                final_max_snr_sq = current_max_snr_sq_vec[v_lane]
                final_max_idx = current_max_idx_vec[v_lane]
                final_max_z = current_max_z_vec[v_lane]
        
        if final_max_idx != -1:
            f_idx_list.append(f_global_idx)
            t_idx_list.append(final_max_idx)
            snr_list.append(final_max_z)
            
    return (f_idx_list, t_idx_list, snr_list)
