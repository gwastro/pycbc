# Copyright (C) 2014 Josh Willis
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

from pycbc.types import zeros, complex64, float32
from scipy.weave import inline
import numpy as _np
import pycbc.opt
from pycbc.opt import omp_support, omp_libs, omp_flags

"""
This module defines several C functions that are weave compiled and compute
a combined thresholding and time clustering of a complex array using a
multithreaded, SIMD vectoried code.

This module also defines several classes that call this function, of which the
last, ThreshClusterObject, is used in the matched filtering control object to
perform the time clustering and thresholding on the output SNR time series
of the matched filtering.

There are three C functions defined:

max_simd: A single-threaded function that uses SIMD vectorization to compute
          the maximum ov a complex time series over a given window.

windowed_max: A single threaded (but not vectorized) function that finds the
              locations, norms, and complex values of the maxima in each of
              a set of predefined windows in a given complex array. It does
              this by calling max_simd on each window with the appropriate
              parameters.

parallel_thresh_cluster: A multithreaded function that finds the maxima in
                         each of a set of contiguous, fixed-length windows
                         by calling windowed_max in parallel. It then sweeps
                         through the results of that (single-threaded) and
                         tests for above threshold, and time-clusters the
                         surviving triggers.

A user calls only the last function; the other two exist to conveniently
compartmentalize SIMD code from OpenMP code.
"""

tc_common_support = omp_support + pycbc.opt.simd_intel_intrin_support + """
#include <stdint.h> // For uint32_t, int64_t
#include <error.h>
#include <complex> // Must use C++ header with weave
#include <math.h> // For M_SQRT2

"""

# The following maximizes over an interval that can be no longer than
# the clustering window.  We do *not* assume alignment, because the window length
# itself may not be a multiple of the alignment.

thresh_cluster_support = tc_common_support + """
void max_simd(float * __restrict inarr, float * __restrict mval,
              float * __restrict norm, int64_t *mloc,
              int64_t nstart, int64_t howmany){

  /*

   This function, using SIMD vectorization, takes an input float
   array (which consists of alternating real and imaginary parts
   of a complex array) and writes the two-float value of the maximum
   into 'mval', the single-float norm of the maximum into 'norm', and
   the location of the maximum (as an index into a *complex* array)
   into 'mloc'. The number 'nstart' is added to whatever location is
   found for the maximum (returned in mloc), and 'howmany' is the length
   of inarr (as a *real* array).

  */

  int64_t i, curr_mloc;
  float re, im, curr_norm, curr, curr_mval[2], *arrptr;

  // Set everything up.  Doesn't depend on any SIMD.
  curr_norm = 0.0;
  curr_mval[0] = 0.0;
  curr_mval[1] = 0.0;
  curr_mloc = 0;
  arrptr = inarr;

  /*

   Note that at most one of _HAVE_AVX and _HAVE_SSE4_1 will be defined (in
   'pycbc.opt.simd_intel_intrin_support' prepended above); if neither is,
   then a non-vectorized code path will be executed.

   As of this writing, documentation on the SIMD instrinsic functions may
   be found at:

      https://software.intel.com/sites/landingpage/IntrinsicsGuide/

   though Intel websites change frequently.

  */

  /*

   The basic vectorized algorithm is to read in complex elements of the input
   array several at a time, using vectorized reads. We then compute the norm
   of the array on an SIMD vector chunk, and compare to the current maximum,
   obtaining a mask for where this vector is larger. Iterating give us the maximum
   over elements 0, V, 2V, 3V, etc; 1, V+1, 2V+1, 3V+1, etc; and so forth, where
   V is the vector length (or a multiple of it, when we are able to unroll the loop).
   After this vectorized maximum, a non-vectorized comparison takes the max over the
   V distinct elements that are the respective maxima for the different array elements
   modulo V.

  */

#if _HAVE_AVX

  __m256 norm_lo, cval_lo, arr_lo, reg0_lo, reg1_lo;
  __m256d mloc_lo, count_lo, incr_lo;
  __m256 norm_hi, cval_hi, arr_hi, reg0_hi, reg1_hi;
  __m256d mloc_hi, count_hi, incr_hi;
  float output_vals[2*ALGN_FLT] __attribute__ ((aligned(ALGN)));
  float output_norms[2*ALGN_FLT] __attribute__ ((aligned(ALGN)));
  double output_locs[2*ALGN_DBL] __attribute__ ((aligned(ALGN)));
  double curr_mloc_dbl;
  int64_t peel, misalgn, j;

  // We calculate size of our various pointers modulo our alignment size

  misalgn = (int64_t)  (((uintptr_t) inarr) % ALGN);

  // Some kinds of misalignment are impossible to handle

  if (misalgn % 2*sizeof(float)) {
    error(EXIT_FAILURE, 0, "Array given to max_simd must be aligned on a least a complex float boundary\\n");
  }

  // 'peel' is how many elements must be handled before we get to
  // something aligned on an SIMD boundary.
  peel = ( misalgn ? ((ALGN - misalgn) / (sizeof(float))) : 0 );
  peel = (peel > howmany ? howmany : peel);

  // Below is the only place i gets initialized! It should never be initialized
  // in the 'for' loops.
  i = 0;

  // Peel off any unaligned beginning for the array.
  for ( ; i < peel; i += 2){
    re = *arrptr;
    im = *(arrptr+1);
    curr = re*re + im*im;
    if (curr > curr_norm){
        curr_mval[0] = re;
        curr_mval[1] = im;
        curr_mloc = i;
        curr_norm = curr;
    }
    arrptr += 2;
  }

  // Now a loop, unrolled once (which is all the AVX registers we can use)

  // AVX---as opposed to AVX2---has very little support for packed
  // integer types.  For instance, one cannot add two packed int
  // SIMD vectors.  So we compromise by storing the array indices
  // as doubles, which should be more than enough to exactly represent
  // a 32 bit integer.

  _mm256_zeroupper();

  // Note that the "set_p{s,d}" functions take their arguments from
  // most-significant value to least.

  incr_lo = _mm256_set_pd(6.0, 4.0, 2.0, 0.0); // incr_lo = [0, 2, 4, 6]
  count_lo = _mm256_set1_pd( (double) i);      // count_lo = [i, i, i, i]
  count_lo = _mm256_add_pd(count_lo, incr_lo); // count_lo = [i, i+2, i+4, i+6]
  incr_lo = _mm256_set_pd(1.0*ALGN_FLT, 1.0*ALGN_FLT, 1.0*ALGN_FLT, 1.0*ALGN_FLT);
                                               // incr_lo = [8, 8, 8, 8]
  count_hi = _mm256_add_pd(count_lo, incr_lo); // count_hi = [i+8, i+10, i+12, i+14]
  incr_lo = _mm256_set_pd(2.0*ALGN_FLT, 2.0*ALGN_FLT, 2.0*ALGN_FLT, 2.0*ALGN_FLT);
                                               // incr_lo = [16, 16, 16, 16]
  incr_hi = _mm256_set_pd(2.0*ALGN_FLT, 2.0*ALGN_FLT, 2.0*ALGN_FLT, 2.0*ALGN_FLT);
                                               // incr_hi = [16, 16, 16, 16]
  // Now count_lo and count_hi have the current indices into the array

  // We don't need to initialize to what we found in the peel-off loop,
  // since we'll store the results of the high and low SIMD loops into
  // an array that we then loop over comparing with the peel-off values.
  mloc_lo = _mm256_setzero_pd();
  norm_lo = _mm256_setzero_ps();
  cval_lo = _mm256_setzero_ps();

  mloc_hi = _mm256_setzero_pd();
  norm_hi = _mm256_setzero_ps();
  cval_hi = _mm256_setzero_ps();

  for (; i <= howmany - 2*ALGN_FLT; i += 2*ALGN_FLT){
      // Load everything into a register
      arr_lo =  _mm256_load_ps(arrptr);
      arr_hi = _mm256_load_ps(arrptr + ALGN_FLT);

      reg0_lo = _mm256_mul_ps(arr_lo, arr_lo);               // 4 x [re*re, im*im]
      reg1_lo = _mm256_shuffle_ps(reg0_lo, reg0_lo, 0xB1);   // 4 x [im*im, re*re]
      reg0_lo = _mm256_add_ps(reg0_lo, reg1_lo);             // 4 x [re^2 +im^2, re^2 +im^2]
      reg1_lo = _mm256_cmp_ps(reg0_lo, norm_lo, _CMP_GT_OQ); // Now a mask for where > curr_norm

      // Now use the mask to selectively update complex value, norm, and location
      mloc_lo = _mm256_blendv_pd(mloc_lo, count_lo, _mm256_castps_pd(reg1_lo) );
      norm_lo = _mm256_blendv_ps(norm_lo, reg0_lo, reg1_lo);
      cval_lo = _mm256_blendv_ps(cval_lo, arr_lo, reg1_lo);

      reg0_hi = _mm256_mul_ps(arr_hi, arr_hi);               // 4 x [re*re, im*im]
      reg1_hi = _mm256_shuffle_ps(reg0_hi, reg0_hi, 0xB1);   // 4 x [im*im, re*re]
      reg0_hi = _mm256_add_ps(reg0_hi, reg1_hi);             // 4 x [re^2 +im^2, re^2 +im^2]
      reg1_hi = _mm256_cmp_ps(reg0_hi, norm_hi, _CMP_GT_OQ); // Now a mask for where > curr_norm

      // Now use the mask to selectively update complex value, norm, and location
      mloc_hi = _mm256_blendv_pd(mloc_hi, count_hi, _mm256_castps_pd(reg1_hi) );
      norm_hi = _mm256_blendv_ps(norm_hi, reg0_hi, reg1_hi);
      cval_hi = _mm256_blendv_ps(cval_hi, arr_hi, reg1_hi);

      count_lo = _mm256_add_pd(count_lo, incr_lo); // count_lo += [16, 16, 16, 16]
      count_hi = _mm256_add_pd(count_hi, incr_hi); // count_hi += [16, 16, 16, 16]
      arrptr += 2*ALGN_FLT;
  }

  // Finally, one last SIMD loop that is not unrolled, just in case we can.
  // We don't reset increments because we won't use them after this loop, and
  // this loop executes at most once.

  for (; i <= howmany - ALGN_FLT; i += ALGN_FLT){
      // Load everything into a register
      arr_lo = _mm256_load_ps(arrptr);

      reg0_lo = _mm256_mul_ps(arr_lo, arr_lo);               // 4 x [re*re, im*im]
      reg1_lo = _mm256_shuffle_ps(reg0_lo, reg0_lo, 0xB1);   // 4 x [im*im, re*re]
      reg0_lo = _mm256_add_ps(reg0_lo, reg1_lo);             // 4 x [re^2 +im^2, re^2 +im^2]
      reg1_lo = _mm256_cmp_ps(reg0_lo, norm_lo, _CMP_GT_OQ); // Now a mask for where > curr_norm

      // Now use the mask to selectively update complex value, norm, and location
      mloc_lo = _mm256_blendv_pd(mloc_lo, count_lo, _mm256_castps_pd(reg1_lo) );
      norm_lo = _mm256_blendv_ps(norm_lo, reg0_lo, reg1_lo);
      cval_lo = _mm256_blendv_ps(cval_lo, arr_lo, reg1_lo);

      arrptr += ALGN_FLT;
  }

  // Now write out the results to our temporary tables:
  _mm256_store_ps(output_vals, cval_lo);
  _mm256_store_ps(output_vals + ALGN_FLT, cval_hi);
  _mm256_store_ps(output_norms, norm_lo);
  _mm256_store_ps(output_norms + ALGN_FLT, norm_hi);
  _mm256_store_pd(output_locs, mloc_lo);
  _mm256_store_pd(output_locs + ALGN_DBL, mloc_hi);

  _mm256_zeroupper();

  // Now loop over our output arrays
  // When we start, curr_norm, curr_mloc, and
  // curr_mval all have the values they had at
  // the end of the *first* peeling loop

  curr_mloc_dbl = (double) curr_mloc;

  for (j = 0; j < 2*ALGN_FLT; j += 2){
    if (output_norms[j] > curr_norm) {
      curr_norm = output_norms[j];
      curr_mloc_dbl = output_locs[j/2];
      curr_mval[0] = output_vals[j];
      curr_mval[1] = output_vals[j+1];
    }
  }

  curr_mloc = (int64_t) curr_mloc_dbl;

#elif _HAVE_SSE4_1

  __m128 norm_lo, cval_lo, arr_lo, reg0_lo, reg1_lo;
  __m128d mloc_lo, count_lo, incr_lo;
  float output_vals[ALGN_FLT] __attribute__ ((aligned(ALGN)));
  float output_norms[ALGN_FLT] __attribute__ ((aligned(ALGN)));
  double output_locs[ALGN_DBL] __attribute__ ((aligned(ALGN)));
  double curr_mloc_dbl;
  int64_t peel, misalgn, j;

  // We calculate size of our various pointers modulo our alignment size

  misalgn = (int64_t)  (((uintptr_t) inarr) % ALGN);

  // Some kinds of misalignment are impossible to handle

  if (misalgn % 2*sizeof(float)) {
    error(EXIT_FAILURE, 0, "Array given to max_simd must be aligned on a least a complex float boundary");
  }

  // 'peel' is how many elements must be handled before we get to
  // something aligned on an SIMD boundary.
  peel = ( misalgn ? ((ALGN - misalgn) / (sizeof(float))) : 0 );
  peel = (peel > howmany ? howmany : peel);

  // Below is the only place i gets initialized! It should never be initialized
  // in the for loops.
  i = 0;

  // Peel off any unaligned beginning for the array.
  for ( ; i < peel; i += 2){
    re = *arrptr;
    im = *(arrptr+1);
    curr = re*re + im*im;
    if (curr > curr_norm){
        curr_mval[0] = re;
        curr_mval[1] = im;
        curr_mloc = i;
        curr_norm = curr;
    }
    arrptr += 2;
  }

  // Note that the "set_p{s,d}" functions take their arguments from
  // most-significant value to least.

  incr_lo = _mm_set_pd(2.0, 0.0);                   // incr_lo = [0, 2]
  count_lo = _mm_set1_pd( (double) i);              // count_lo = [i, i]
  count_lo = _mm_add_pd(count_lo, incr_lo);         // count_lo = [i, i+2]
  incr_lo = _mm_set_pd(1.0*ALGN_FLT, 1.0*ALGN_FLT); // incr_lo = [4, 4]

  // Now count_lo has the current indices into the array

  // We don't need to initialize to what we found in the peel-off loop,
  // since we'll store the results of the high and low SIMD loops into
  // an array that we then loop over comparing with the peel-off values.
  mloc_lo = _mm_setzero_pd();
  norm_lo = _mm_setzero_ps();
  cval_lo = _mm_setzero_ps();

  for (; i <= howmany - ALGN_FLT; i += ALGN_FLT){
      // Load everything into a register
      arr_lo =  _mm_load_ps(arrptr);

      reg0_lo = _mm_mul_ps(arr_lo, arr_lo);               // 2 x [re*re, im*im]
      reg1_lo = _mm_shuffle_ps(reg0_lo, reg0_lo, 0xB1);   // 2 x [im*im, re*re]
      reg0_lo = _mm_add_ps(reg0_lo, reg1_lo);             // 2 x [re^2 +im^2, re^2 +im^2]
      reg1_lo = _mm_cmpgt_ps(reg0_lo, norm_lo);           // Now a mask for where > curr_norm

      // Now use the mask to selectively update complex value, norm, and location
      mloc_lo = _mm_blendv_pd(mloc_lo, count_lo, _mm_castps_pd(reg1_lo) );
      norm_lo = _mm_blendv_ps(norm_lo, reg0_lo, reg1_lo);
      cval_lo = _mm_blendv_ps(cval_lo, arr_lo, reg1_lo);

      count_lo = _mm_add_pd(count_lo, incr_lo); // count_lo += [4, 4]
      arrptr += ALGN_FLT;
  }

  // Now write out the results to our temporary tables:
  _mm_store_ps(output_vals, cval_lo);
  _mm_store_ps(output_norms, norm_lo);
  _mm_store_pd(output_locs, mloc_lo);

  // Now loop over our output arrays
  // When we start, curr_norm, curr_mloc, and
  // curr_mval all have the values they had at
  // the end of the *first* peeling loop

  curr_mloc_dbl = (double) curr_mloc;

  for (j = 0; j < ALGN_FLT; j += 2){
    if (output_norms[j] > curr_norm) {
      curr_norm = output_norms[j];
      curr_mloc_dbl = output_locs[j/2];
      curr_mval[0] = output_vals[j];
      curr_mval[1] = output_vals[j+1];
    }
  }

  curr_mloc = (int64_t) curr_mloc_dbl;

#else
 // If we have no SSE, all we have to do is initialize
 // our loop counter, and the last "cleanup" loop
 // will in fact do all the work.

 i = 0;

#endif

  for ( ; i < howmany; i += 2){
    re = *arrptr;
    im = *(arrptr+1);
    curr = re*re + im*im;
    if (curr > curr_norm){
        curr_mval[0] = re;
        curr_mval[1] = im;
        curr_mloc = i;
        curr_norm = curr;
    }
    arrptr += 2;
  }

  // Store our answers and return
  *mval = curr_mval[0];
  *(mval+1) = curr_mval[1];
  *norm = curr_norm;

  // Note that curr_mloc is a real array index, but we
  // must return the index into the complex array.
  *mloc = (curr_mloc/2) + nstart;

  return;

}

void max_simple(float * __restrict inarr, float * __restrict mval,
                float * __restrict norm, int64_t *mloc,
                int64_t nstart, int64_t howmany){

  /*

   This function does *NOT* use explicit SIMD vectorization, and
   takes an input float array (which consists of alternating real and
   imaginary parts of a complex array) and writes the two-float value
   of the maximum into 'mval', the single-float norm of the maximum into
   'norm', and the location of the maximum (as an index into a *complex* array)
   into 'mloc'. The number 'nstart' is added to whatever location is
   found for the maximum (returned in mloc), and 'howmany' is the length
   of inarr (as a *real* array).

  */

  int64_t i, curr_mloc;
  float re, im, curr_norm, curr, curr_mval[2], *arrptr;

  // Set everything up.
  curr_norm = 0.0;
  curr_mval[0] = 0.0;
  curr_mval[1] = 0.0;
  curr_mloc = 0;
  arrptr = inarr;

  for (i = 0; i < howmany; i += 2){
    re = *arrptr;
    im = *(arrptr+1);
    curr = re*re + im*im;
    if (curr > curr_norm){
        curr_mval[0] = re;
        curr_mval[1] = im;
        curr_mloc = i;
        curr_norm = curr;
    }
    arrptr += 2;
  }

  // Store our answers and return
  *mval = curr_mval[0];
  *(mval+1) = curr_mval[1];
  *norm = curr_norm;

  // Note that curr_mloc is a real array index, but we
  // must return the index into the complex array.
  *mloc = (curr_mloc/2) + nstart;

  return;

}


void windowed_max(std::complex<float> * __restrict inarr, const int64_t arrlen,
                  std::complex<float> * __restrict cvals, float * __restrict norms,
                  int64_t * __restrict locs, const int64_t winsize,
                  const int64_t startoffset){


  /*

   This function fills in the arrays cvals, norms, and locs, with the max (as
   complex value, norm, and location, resp.) of the array inarr.  The length of
   inarr is arrlen, and the function assumes that it computes the max over successive
   windows of length winsize, starting at the beginning of the array and continuing
   to the end.  If winsize does not evenly divide arrlen, then the last partial window
   will be maximized over as well.  If winsize is greater than arrlen, then just one
   maximization will be performed over the (partial) array inarr.

   Thus, in all cases, the length of cvals, norms, and locs should be:
      nwindows = ( (arrlen % winsize) ? (arrlen/winsize) + 1 : (arrlen/winsize) )

   Note that all input sizes are for *complex* arrays; the function this function calls
   often requires *real* arrays, and lengths will be converted where appropriate.

  */

  int64_t i, nwindows;

  nwindows = ( (arrlen % winsize) ? (arrlen/winsize) + 1 : (arrlen/winsize) );

  // Everything but the last window, which may not be full length

  for (i = 0; i < nwindows-1; i++){
    // The factor of 2 multiplying lengths[i] is because max_simd needs its length as a real
    // length, not complex.  But startpts and startoffset are complex values.
    max_simple((float *) &inarr[i*winsize], (float *) &cvals[i],
               &norms[i], &locs[i], startoffset + i*winsize, 2*winsize);
  }
  // Now the last window (which will be the only window if arrlen <= winzise)
  max_simple((float *) &inarr[i*winsize], (float *) &cvals[i],
              &norms[i], &locs[i], startoffset + i*winsize, 2*(arrlen - i*winsize));

  return;
}

int parallel_thresh_cluster(std::complex<float> * __restrict inarr, const uint32_t arrlen,
                            std::complex<float> * __restrict values, uint32_t * __restrict locs,
                            const float thresh, const uint32_t winsize, const uint32_t segsize){

  /*

  This function takes a complex input array 'inarr', of length 'arrlen', and returns
  in the complex array 'values' the time-clustered values of all maxima within a window
  of size 'winsize'. The locations (as indices into the original array) are returned in
  the array 'locs'. Both 'values' and 'locs' must be pre-allocated. The time-clustered
  triggers are only returned when their norm is above the value 'thresh'. The last argument,
  'segsize', specifies in what size chunks the array should be processed, and should be no
  larger than what can fit in the cache local to a single processor core (the parallelization
  will call 'windowed_max' in parallel on chunks of this size).

  */

  int64_t i, nsegs, nwins_ps, last_arrlen, last_nwins_ps, outlen, curr_mloc;
  int64_t *seglens, *mlocs, cnt, s_segsize, s_arrlen, s_winsize, curr_mark, cluster_win;
  float *norms, thr_sqr, curr_norm;
  std::complex<float> *cvals, curr_cval;
  short int *marks;

  thr_sqr = (thresh * thresh);
  // Signed versions of all our array length inputs, so loop
  // logic and calls to other functions are safe.
  s_segsize = (int64_t) segsize; 
  s_arrlen = (int64_t) arrlen;
  // We divide segsize by two since our initial pass must be with a segment
  // half the length of that we will use for clustering.
  s_winsize = (int64_t) winsize/2;
  cluster_win = (int64_t) winsize;

  // If segsize divides arrlen evenly, then the number of segments 'nsegs' is that
  // ratio; if not, it is one more than the floor of that ratio and the last segment
  // will be shorter than all of the others.
  nsegs = ( (s_arrlen % s_segsize) ? (s_arrlen/s_segsize) + 1 : (s_arrlen/s_segsize) );

  // If winsize divides segsize evenly, then the number of windows per segment
  // 'nwins_ps' is that ratio; if not, it is one more than the floor of that ratio
  // and the last window in each segment will be shorter than all of the others.
  nwins_ps = ( (s_segsize % s_winsize) ? (s_segsize/s_winsize) + 1 : (s_segsize/s_winsize) );

  // Our code will always handle the last segment separately.  However if
  // segsize evenly divides arrlen, then the last segment will be no different.
  // The following ternary operator captures that logic:
  last_arrlen = ( (s_arrlen % s_segsize) ? (s_arrlen - (nsegs-1) * s_segsize) : (s_segsize) );
  last_nwins_ps = ( (last_arrlen % s_winsize) ? (last_arrlen/s_winsize) + 1 : (last_arrlen/s_winsize) );

  // Then the total length of the working arrays we must dynamically allocate is:
  outlen = (nsegs-1) * nwins_ps + last_nwins_ps;

  // Now dynamic allocation.  No reason really to align this memory; it will be parceled
  // out to different cores anyway.

  cvals = (std::complex<float> *) malloc(outlen * sizeof(std::complex<float>) );
  norms = (float *) malloc(outlen * sizeof(float) );
  mlocs = (int64_t *) malloc(outlen * sizeof(int64_t) );

  // The next array will be used in our forward/reverse algorithm for
  // clustering.  Note the calloc, rather than malloc!
  marks = (short int *) calloc((size_t) outlen, sizeof(short int));

  // The next array allows us to dynamically communicate possibly changed sizes to the
  // many parallel calls to windowed_max:

  seglens = (int64_t *) malloc(nsegs * sizeof(int64_t) );

  // check to see if anything failed
  if ( (cvals == NULL) || (norms == NULL) || (mlocs == NULL) || (seglens == NULL) || (marks == NULL) ){
    error(EXIT_FAILURE, ENOMEM, "Could not allocate temporary memory needed by parallel_thresh_cluster");
  }

  for (i = 0; i < (nsegs-1); i++){
    seglens[i] = segsize;
  }
  seglens[i] = last_arrlen;

  // Now the real work, in an OpenMP parallel for loop:
#pragma omp parallel for schedule(dynamic,1)
  for (i = 0; i < nsegs; i++){
    windowed_max(&inarr[i*segsize], seglens[i], &cvals[i*nwins_ps],
                 &norms[i*nwins_ps], &mlocs[i*nwins_ps],
                 s_winsize, i*s_segsize);
  }

  /*

  We have now the requisite maxima over our windows in cvals,
  norms, and mlocs. We want to apply the threshold and cluster over
  time.

  Our time clustering algorithm is that we should keep a candidate
  trigger if it is above threshold and louder than anything before
  or after it within 'window' in index. We test this by sliding through
  our candidate triggers and comparing to what comes before and after.

  */

  // First, go through the data *backwards*, and mark as valid anything
  // that survives a sliding window in this direction. These candidates
  // are guaranteed to be larger than anything that comes *earlier* than
  // them within a window size.

  curr_norm = norms[outlen-1];
  curr_mloc = mlocs[outlen-1];
  curr_mark = outlen - 1;

  // For this pass we just use 'cnt' as a logical, noting
  // whether we have found a point larger than something
  // after it.  Otherwise writing out 'marks' after the loop
  // could be mistaken.
  //
  // Note that in this reverse loop, 'curr_mark' records
  // where in the arrays (of lengths 'outlen') we found
  // our last potential local maximum.
  cnt = 0;
  for (i = outlen-2; i >= 0; i--){
    if ( (curr_mloc - mlocs[i]) > cluster_win){
      marks[curr_mark] = 1;
      curr_mark = i;
      curr_norm = norms[i];
      curr_mloc = mlocs[i];
      cnt = 1;
    } else if (norms[i] > curr_norm) {
      // Note that we required strictly greater than:
      // if there is a sequence of several equal values
      // all within a window, only the greatest in index
      // will be marked (rightmost).
      curr_mark = i;
      curr_norm = norms[i];
      curr_mloc = mlocs[i];
      cnt = 1;
    }
  }
  // We may not have marked the last point, so do so:
  if (cnt) marks[curr_mark] = 1;

  // Now we have a sliding forward window; if something
  // is marked and survives this, then it is a clustered
  // trigger. This is also the only place where we apply
  // the threshold condition.

  cnt = 0;
  curr_norm = norms[0];
  curr_mloc = mlocs[0];
  curr_cval = cvals[0];
  // Note that in this pass, we treat 'curr_mark' as a
  // boolean, to know whether our forward potential
  // maximum was also marked on the reverse loop earlier.
  curr_mark = marks[0];

  for (i = 1; i < outlen; i++){    
    if ( (mlocs[i] - curr_mloc) > cluster_win){
      // The last one is a maximum for all points following,
      // so if also for points preceding (curr_mark) and
      // if above threshold, then write it out.
      if ( (curr_norm > thr_sqr) && curr_mark){
        cnt += 1;
        values[cnt-1] = curr_cval;
        locs[cnt-1] = (uint32_t) curr_mloc;
      }
      // Even if we didn't write it out, we still update
      // our current max
      curr_cval = cvals[i];
      curr_norm = norms[i];
      curr_mloc = mlocs[i];
      curr_mark = marks[i];
    } else if (norms[i] >= curr_norm) {
      // Here we allow greater-than or equal-to, since
      // only the rightmost may have been marked in the
      // first pass, we want the rightmost in a sequence
      // of equal values to count as the maximum.
      curr_cval = cvals[i];
      curr_norm = norms[i];
      curr_mloc = mlocs[i];
      curr_mark = marks[i];
    }
  }

  // It's possible that the last point would survive as a trigger,
  // so we need a separate test for that.

  if ((cnt > 0) && (curr_mloc != locs[cnt-1])){
    if ((curr_norm > thr_sqr) && curr_mark){
      cnt += 1;
      values[cnt-1] = curr_cval;
      locs[cnt-1] = (uint32_t) curr_mloc;
    }
  }

  free(cvals);
  free(norms);
  free(mlocs);
  free(seglens);
  free(marks);

  return cnt;
}

"""

### Now some actual code that just implements the different
### correlations in a parallelized fashion.

# First, a simple thing, that just does the max...

max_only_code = """
max_simd(inarr, mval, norm, (int64_t *) mloc, (int64_t) nstart[0], (int64_t) howmany[0]);
"""

class MaxOnlyObject(object):
    def __init__(self, inarray, verbose=0):
        self.inarr = _np.array(inarray.data, copy=False).view(dtype = float32)
        self.howmany = _np.zeros(1, dtype = _np.int64)
        self.howmany[0] = len(self.inarr)
        self.nstart = _np.zeros(1, dtype = _np.int64)
        self.nstart[0] = 0
        self.cmplx_mval = zeros(1, dtype = complex64)
        self.mval = _np.array(self.cmplx_mval.data, copy = False).view(dtype = float32)
        self.norm = _np.zeros(1, dtype = float32)
        self.mloc = _np.zeros(1, dtype = _np.int64)
        self.code = max_only_code
        self.support = thresh_cluster_support
        self.verbose = verbose

    def execute(self):
        inarr = self.inarr
        mval = self.mval
        norm = self.norm
        mloc = self.mloc
        nstart = self.nstart
        howmany = self.howmany
        inline(self.code, ['inarr', 'mval', 'norm', 'mloc', 'nstart', 'howmany'],
               extra_compile_args = ['-march=native -O3 -w'],
               #extra_compile_args = ['-mno-avx -mno-sse2 -mno-sse3 -mno-ssse3 -mno-sse4 -mno-sse4.1 -mno-sse4.2 -mno-sse4a -O2 -w'],
               #extra_compile_args = ['-msse4.1 -O3 -w'],
               support_code = self.support, auto_downcast = 1, verbose = self.verbose)

windowed_max_code = """
windowed_max(inarr, (int64_t) arrlen[0], cvals, norms, (int64_t *) locs, (int64_t ) winsize[0],
             (int64_t) startoffset[0]);
"""

class WindowedMaxObject(object):
    def __init__(self, inarray, winsize, verbose=0):
        self.inarr = _np.array(inarray.data, copy=False)
        self.arrlen = _np.zeros(1, dtype = _np.int64)
        self.arrlen[0] = len(self.inarr)
        self.len_win = winsize
        nwindows = int( len(self.inarr) / winsize)
        if (nwindows * winsize < len(self.inarr)):
            nwindows = nwindows + 1
        self.nwindows = nwindows
        self.cvals = _np.zeros(self.nwindows, dtype = complex64)
        self.norms = _np.zeros(self.nwindows, dtype = float32)
        self.locs = _np.zeros(self.nwindows, dtype = _np.int64)
        self.winsize = _np.zeros(1, dtype = _np.int64)
        self.winsize[0] = self.len_win
        self.startoffset = _np.zeros(1, dtype = _np.int64)
        self.code = windowed_max_code
        self.support = thresh_cluster_support
        self.verbose = verbose

    def execute(self):
        inarr = self.inarr
        arrlen = self.arrlen
        cvals = self.cvals
        norms = self.norms
        locs = self.locs
        winsize = self.winsize
        startoffset = self.startoffset
        inline(self.code, ['inarr', 'arrlen', 'cvals', 'norms', 'locs', 'winsize', 'startoffset'],
               extra_compile_args = ['-march=native -O3 -w'],
               #extra_compile_args = ['-mno-avx -mno-sse2 -mno-sse3 -mno-ssse3 -mno-sse4 -mno-sse4.1 -mno-sse4.2 -mno-sse4a -O2 -w'],
               #extra_compile_args = ['-msse4.1 -O3 -w'],
               support_code = self.support, auto_downcast = 1, verbose = self.verbose)



thresh_cluster_code = """
return_val = parallel_thresh_cluster(series, (uint32_t) slen, values, locs,
                                     (float) thresh, (uint32_t) window, (uint32_t) segsize);
"""

if pycbc.opt.HAVE_GETCONF:
    default_segsize = pycbc.opt.LEVEL2_CACHE_SIZE / _np.dtype( _np.complex64).itemsize
else:
    # Seems to work for Sandy Bridge/Ivy Bridge/Haswell, for now?
    default_segsize = 32768

class ThreshClusterObject(object):
    """
    This class takes a complex SNR time series, a real threshold value, and a window size
    (expressed as a number of complex sample points), and optionally a segment size and
    verbosity indicator (for debugging weave).

    The execute method returns two numpy arrays: one giving the location of clustered
    maxima above the threshold, and the other the corresponding (complex) values of the
    SNR time series at those clustered maxima.  The expectation is that the memory of
    the SNR time series will be filled new inputs many times, and execute() called
    repeatedly.
    """
    def __init__(self, series, window, segsize = default_segsize, verbose=0):
        self.series = _np.array(series.data, copy=False)
        self.slen = len(self.series)
        nwindows = int( self.slen / window)
        if (nwindows * window < self.slen):
            nwindows = nwindows + 1
        self.nwindows = nwindows
        self.values = _np.zeros(self.nwindows, dtype = complex64)
        self.locs = _np.zeros(self.nwindows, dtype = _np.uint32)
        self.window = window
        self.segsize = segsize
        self.code = thresh_cluster_code
        self.support = thresh_cluster_support
        self.verbose = verbose

    def execute(self, thresh):
        series = self.series
        slen = self.slen
        values = self.values
        locs = self.locs
        window = self.window
        segsize = self.segsize
        nthr = inline(self.code, ['series', 'slen', 'values', 'locs', 'thresh', 'window', 'segsize'],
                      extra_compile_args = ['-march=native -O3 -w'] + omp_flags,
                      #extra_compile_args = ['-mno-avx -mno-sse2 -mno-sse3 -mno-ssse3 -mno-sse4 -mno-sse4.1 -mno-sse4.2 -mno-sse4a -O2 -w'],
                      #extra_compile_args = ['-msse4.1 -O3 -w'],
                      support_code = self.support, libraries = omp_libs,
                      auto_downcast = 1, verbose = self.verbose)
        if nthr > 0:
            return self.values[0:nthr], self.locs[0:nthr]
        else:
            return _np.array([], dtype = complex64), _np.array([], dtype = _np.uint32)

