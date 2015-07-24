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

from pycbc.types import float32
from scipy.weave import inline
import numpy as _np
import pycbc.opt
from pycbc.opt import omp_support, omp_libs, omp_flags

"""
This module contains C functions (to be compiled by weave) for multiplying
the complex conjugate of one vector by a second vector, writing the output
to a third vector. They do this multi-threaded and with SIMD vectorization.

The long string of code defined here, and the other calling that function,
are imported and used in the CPUCorrelator class defined in
matchedfilter_cpu.py.

Two functions are defined in the 'support' long code string:

ccorrf_simd: Runs on a single core, but vectorized
ccorrf_parallel: Runs multicore, but not explicitly vectorized.
                 Parallelized using OpenMP, and calls ccorrf_simd
"""

# Common support, so defined here


corr_common_support = omp_support + pycbc.opt.simd_intel_intrin_support + """
#include <stdint.h> // For uint32_t, int64_t
#include <error.h>
#include <complex> // Must use C++ header with weave
"""

corr_support = corr_common_support + """
void ccorrf_simd(float * __restrict inconj, float * __restrict innoconj,
                 float * __restrict out, const int64_t len){

  /*
   This function takes two input vectors, and puts in the third argument the output
   vector that is the elementwise product of the conjugate of the first input and
   the second input.  All vectors are of type float, but are interpreted as complex
   vectors; the length should be the length of the *real* vectors (which must be
   the same). The last argument is the common length of all three vectors.

   The vectors do not have to be aligned on any particular SIMD boundary, but the
   *relative* misalignment of all three vectors must be the same.  In practice this
   means that the vectors will often be allocated as the same index range of three
   vectors that *are* aligned.
  */

  int64_t i;
  float ar, ai, br, bi, cr, ci;
  float *aptr, *bptr, *cptr;

  aptr = inconj;
  bptr = innoconj;
  cptr = out;

  /*

   Note that at most one of _HAVE_AVX ad _HAVE_SSE3 will be defined (in
   'pycbc.opt.simd_intel_intrin_support' prepended above); if neither is,
   then a non-vectorized code path will be executed.

   As of this writing, documentation on the SIMD instrinsic functions may
   be found at:

      https://software.intel.com/sites/landingpage/IntrinsicsGuide/

   though Intel websites change frequently.

  */

#if _HAVE_AVX
  __m256 ymm0, ymm1, ymm2, ymm3, ymm4, ymm5, signed_zero;
  int64_t peel, misalgna, misalgnb, misalgnc;

  // We calculate size of our various pointers modulo our alignment size

  misalgna = (int64_t) (((uintptr_t) aptr) % ALGN);
  misalgnb = (int64_t) (((uintptr_t) bptr) % ALGN);
  misalgnc = (int64_t) (((uintptr_t) cptr) % ALGN);

  // Some kinds of misalignment are impossible to handle

  if ((misalgna != misalgnb) || (misalgnb != misalgnc)){
    error(EXIT_FAILURE, 0, "Arrays given to ccorrf_simd must all three have same alignment\\n");
  }

  if (misalgna % 2*sizeof(float)) {
    error(EXIT_FAILURE, 0, "Arrays given to ccorrf_simd must be aligned on a least a complex float boundary\\n");
  }

  // If we were not aligned, we must loop through a few iterations
  // before we may use SIMD functions.

  peel = ( misalgna ? ((ALGN - misalgna) / (sizeof(float))) : 0 );
  peel = (peel > len ? len : peel);

  // Below is the only place i gets initialized! It should never be initialized
  // in the 'for' loops.
  i = 0;

  // Peel off any unaligned beginning for the array.
  for ( ; i < peel; i += 2){
    ar = *aptr;
    ai = *(aptr+1);
    br = *bptr;
    bi = *(bptr+1);

    // Note that 'a' is the input we conjugate

    cr = ar*br + ai*bi;
    ci = ar*bi - ai*br;

    // Store output and increment pointers

    *cptr = cr;
    *(cptr+1) = ci;

    aptr += 2;
    bptr += 2;
    cptr += 2;
  }

  _mm256_zeroupper();

  // We use the following to flip sign bits of imaginary parts
  // of the first vector.
  signed_zero = _mm256_castsi256_ps( _mm256_set1_epi32(0x80000000) );

  // Main loop using AVX; unrolled once.
  for (; i <= len - 2*ALGN_FLT ; i += 2*ALGN_FLT){

      // Load everything into registers
      ymm1 = _mm256_load_ps(aptr);
      ymm4 = _mm256_load_ps(aptr+ALGN_FLT);
      ymm0 = _mm256_load_ps(bptr);
      ymm3 = _mm256_load_ps(bptr+ALGN_FLT);

      // First iteration
      ymm2 = _mm256_movehdup_ps(ymm1);  // only imag parts (duplicated) of inconj
      ymm1 = _mm256_moveldup_ps(ymm1);  // only real parts (duplicated) of inconj
      ymm1 = _mm256_mul_ps(ymm1, ymm0); // 4 x [inconj_re * innoconj_re, inconj_re * innoconj_im]
      ymm0 = _mm256_shuffle_ps(ymm0, ymm0, 0xB1); // Swap re and im in innoconj
      ymm2 = _mm256_xor_ps(ymm2, signed_zero);    // Change sign of inconj imag
      ymm2 = _mm256_mul_ps(ymm2, ymm0);    // 4 x [-inconj_im * innoconj_im, -inconj_im * innoconj_re]
      ymm0 = _mm256_addsub_ps(ymm1, ymm2); // 4 x [inconj_re * innoconj_re + inconj_im * innoconj_im,
                                           //      inconj_re * innoconj_im - inconj_im * innoconj_re]

      // Second iteration
      ymm5 = _mm256_movehdup_ps(ymm4);  // only imag parts (duplicated) of inconj
      ymm4 = _mm256_moveldup_ps(ymm4);  // only real parts (duplicated) of inconj
      ymm4 = _mm256_mul_ps(ymm4, ymm3); // 4 x [inconj_re * innoconj_re, inconj_re * innoconj_im]
      ymm3 = _mm256_shuffle_ps(ymm3, ymm3, 0xB1); // Swap re and im in innoconj
      ymm5 = _mm256_xor_ps(ymm5, signed_zero);    // Change sign of inconj imag
      ymm5 = _mm256_mul_ps(ymm5, ymm3);    // 4 x [-inconj_im * innoconj_im, -inconj_im * innoconj_re]
      ymm3 = _mm256_addsub_ps(ymm4, ymm5); // 4 x [inconj_re * innoconj_re + inconj_im * innoconj_im,
                                           //      inconj_re * innoconj_im - inconj_im * innoconj_re]

      // Write out results
      _mm256_store_ps(cptr, ymm0);
      _mm256_store_ps(cptr+ALGN_FLT, ymm3);

      aptr += 2*ALGN_FLT;
      bptr += 2*ALGN_FLT;
      cptr += 2*ALGN_FLT;
  }

  // One last SIMD loop, not unrolled (which we may not always be able to do)
  for (; i <= len - ALGN_FLT; i += ALGN_FLT){
      // Load everything into registers
      ymm1 = _mm256_load_ps(aptr);
      ymm0 = _mm256_load_ps(bptr);

      // Do the work
      ymm2 = _mm256_movehdup_ps(ymm1);  // only imag parts (duplicated) of inconj
      ymm1 = _mm256_moveldup_ps(ymm1);  // only real parts (duplicated) of inconj
      ymm1 = _mm256_mul_ps(ymm1, ymm0); // 4 x [inconj_re * innoconj_re, inconj_re * innoconj_im]
      ymm0 = _mm256_shuffle_ps(ymm0, ymm0, 0xB1); // Swap re and im in innoconj
      ymm2 = _mm256_xor_ps(ymm2, signed_zero);    // Change sign of inconj imag
      ymm2 = _mm256_mul_ps(ymm2, ymm0);    // 4 x [-inconj_im * innoconj_im, -inconj_im * innoconj_re]
      ymm0 = _mm256_addsub_ps(ymm1, ymm2); // 4 x [inconj_re * innoconj_re + inconj_im * innoconj_im,
                                           //      inconj_re * innoconj_im - inconj_im * innoconj_re]

      // Write out results
      _mm256_store_ps(cptr, ymm0);

      aptr += ALGN_FLT;
      bptr += ALGN_FLT;
      cptr += ALGN_FLT;
  }
  _mm256_zeroupper();

#elif _HAVE_SSE3  // Need at least SSE3 for move{h,l}dup and addsub instructions

  __m128 ymm0, ymm1, ymm2, ymm3, ymm4, ymm5, signed_zero;
  int64_t peel, misalgna, misalgnb, misalgnc;

  // We calculate size of our various pointers modulo our alignment size

  misalgna = (int64_t) (((uintptr_t) aptr) % ALGN);
  misalgnb = (int64_t) (((uintptr_t) bptr) % ALGN);
  misalgnc = (int64_t) (((uintptr_t) cptr) % ALGN);

  // Some kinds of misalignment are impossible to handle

  if ((misalgna != misalgnb) || (misalgnb != misalgnc)){
    error(EXIT_FAILURE, 0, "Arrays given to ccorrf_simd must all three have same alignment\\n");
  }

  if (misalgna % 2*sizeof(float)) {
    error(EXIT_FAILURE, 0, "Arrays given to ccorrf_simd must be aligned on a least a complex float boundary\\n");
  }

  // If we were not aligned, we must loop through a few iterations
  // before we may use SIMD functions.

  peel = ( misalgna ? ((ALGN - misalgna) / (sizeof(float))) : 0 );
  peel = (peel > len ? len : peel);

  // Below is the only place i gets initialized! It should never be initialized
  // in the 'for' loops.
  i = 0;

  // Peel off any unaligned beginning for the array.
  for ( ; i < peel; i += 2){
    ar = *aptr;
    ai = *(aptr+1);
    br = *bptr;
    bi = *(bptr+1);

    // Note that 'a' is the input we conjugate

    cr = ar*br + ai*bi;
    ci = ar*bi - ai*br;

    // Store output and increment pointers

    *cptr = cr;
    *(cptr+1) = ci;

    aptr += 2;
    bptr += 2;
    cptr += 2;
  }

  // We use the following to flip sign bits of imaginary parts
  // of the first vector.
  signed_zero = _mm_castsi128_ps( _mm_set1_epi32(0x80000000) );

  // Main loop using SSE; unrolled once.
  for (; i <= len - 2*ALGN_FLT ; i += 2*ALGN_FLT){

      // Load everything into registers
      ymm1 = _mm_load_ps(aptr);
      ymm4 = _mm_load_ps(aptr+ALGN_FLT);
      ymm0 = _mm_load_ps(bptr);
      ymm3 = _mm_load_ps(bptr+ALGN_FLT);

      // First iteration
      ymm2 = _mm_movehdup_ps(ymm1);  // only imag parts (duplicated) of inconj
      ymm1 = _mm_moveldup_ps(ymm1);  // only real parts (duplicated) of inconj
      ymm1 = _mm_mul_ps(ymm1, ymm0); // 4 x [inconj_re * innoconj_re, inconj_re * innoconj_im]
      ymm0 = _mm_shuffle_ps(ymm0, ymm0, 0xB1); // Swap re and im in innoconj
      ymm2 = _mm_xor_ps(ymm2, signed_zero);    // Change sign of inconj imag
      ymm2 = _mm_mul_ps(ymm2, ymm0);    // 4 x [-inconj_im * innoconj_im, -inconj_im * innoconj_re]
      ymm0 = _mm_addsub_ps(ymm1, ymm2); // 4 x [inconj_re * innoconj_re + inconj_im * innoconj_im,
                                        //      inconj_re * innoconj_im - inconj_im * innoconj_re]

      // Second iteration
      ymm5 = _mm_movehdup_ps(ymm4);  // only imag parts (dupl) of inconj
      ymm4 = _mm_moveldup_ps(ymm4);  // only real parts (dupl) of inconj
      ymm4 = _mm_mul_ps(ymm4, ymm3); // 4 x [inconj_re * innoconj_re, inconj_re * innoconj_im]
      ymm3 = _mm_shuffle_ps(ymm3, ymm3, 0xB1); // Swap re and im in innoconj
      ymm5 = _mm_xor_ps(ymm5, signed_zero);    // Change sign of inconj imag
      ymm5 = _mm_mul_ps(ymm5, ymm3);    // 4 x [-inconj_im * innoconj_im, -inconj_im * innoconj_re]
      ymm3 = _mm_addsub_ps(ymm4, ymm5); // 4 x [inconj_re * innoconj_re + inconj_im * innoconj_im,
                                        //      inconj_re * innoconj_im - inconj_im * innoconj_re]

      // Write out results
      _mm_store_ps(cptr, ymm0);
      _mm_store_ps(cptr+ALGN_FLT, ymm3);

      aptr += 2*ALGN_FLT;
      bptr += 2*ALGN_FLT;
      cptr += 2*ALGN_FLT;
  }

  // One last SIMD loop, not unrolled (which we may not always be able to do)
  for (; i <= len - ALGN_FLT; i += ALGN_FLT){
      // Load everything into registers
      ymm1 = _mm_load_ps(aptr);
      ymm0 = _mm_load_ps(bptr);

      // Do the work
      ymm2 = _mm_movehdup_ps(ymm1);  // only imag parts (duplicated) of inconj
      ymm1 = _mm_moveldup_ps(ymm1);  // only real parts (duplicated) of inconj
      ymm1 = _mm_mul_ps(ymm1, ymm0); // 4 x [inconj_re * innoconj_re, inconj_re * innoconj_im]
      ymm0 = _mm_shuffle_ps(ymm0, ymm0, 0xB1); // Swap re and im in innoconj
      ymm2 = _mm_xor_ps(ymm2, signed_zero);    // Change sign of inconj imag
      ymm2 = _mm_mul_ps(ymm2, ymm0);    // 4 x [-inconj_im * innoconj_im, -inconj_im * innoconj_re]
      ymm0 = _mm_addsub_ps(ymm1, ymm2); // 4 x [inconj_re * innoconj_re + inconj_im * innoconj_im,
                                        //      inconj_re * innoconj_im - inconj_im * innoconj_re]

      // Write out results
      _mm_store_ps(cptr, ymm0);

      aptr += ALGN_FLT;
      bptr += ALGN_FLT;
      cptr += ALGN_FLT;
  }

#else
 // If we have no SSE, all we have to do is initialize
 // our loop counter, and the last "cleanup" loop
 // will in fact do all the work.

 i = 0;

#endif

  for ( ; i < len; i += 2){
    ar = *aptr;
    ai = *(aptr+1);
    br = *bptr;
    bi = *(bptr+1);

    // Note that 'a' is the input we conjugate

    cr = ar*br + ai*bi;
    ci = ar*bi - ai*br;

    // Store output and increment pointers

    *cptr = cr;
    *(cptr+1) = ci;

    aptr += 2;
    bptr += 2;
    cptr += 2;
  }

 return;

}

void ccorrf_parallel(std::complex<float> * __restrict inconj,
                     std::complex<float> * __restrict innoconj,
                     std::complex<float> * __restrict out,
                     const int64_t arrlen, const int64_t segsize){

  /*

   This function takes three complex vectors: the complex conjugate of
   the first is elementwise multiplied by the second with the output
   written into the third. The fourth argument is the common length of
   all three vectors, and the fifth argument is  the 'segment' size. It
   is the length of a vector that should be handled by a single thread,
   and should generally be such that three vectors of that length will
   fit in the cache local to a single CPU core (not shared).

  */


  int64_t i, *seglens;
  int nsegs;

  // We will always give 'segsize' length arrays. If segsize does not
  // divide arrlen evenly, the last triple of arrays will be shorter.
  // The next lines enforce that logic: if segsize divides arrlen
  // evenly, the number of segments is their ratio; if not, the number
  // of segments is one more than the floor of that ratio, and the
  // last segment has a shorter length than the others.

  nsegs = ( (arrlen % segsize) ? (arrlen/segsize) + 1 : (arrlen/segsize) );

  seglens = (int64_t *) malloc(nsegs * sizeof(int64_t));
  if (seglens == NULL){
    error(EXIT_FAILURE, ENOMEM, "ccorrf_parallel: could not allocate temporary memory");
  }

  // Setup the segment lengths; only the last might be different.
  // The factor of two is because our inputs and outputs are complex
  // arrays, but ccorrf_simd needs float array arguments.
  for (i = 0; i < nsegs-1; i++){
    seglens[i] = 2*segsize;
  }
  seglens[i] = 2*(arrlen - i*segsize);

#pragma omp parallel for schedule(dynamic,1)
  for (i = 0; i < nsegs; i++){
    ccorrf_simd( (float *) &inconj[i*segsize], (float *) &innoconj[i*segsize],
                 (float *) &out[i*segsize], seglens[i]);
  }

  free(seglens);
  return;
}

"""

corr_simd_code = """
ccorrf_simd(htilde, stilde, qtilde, (int64_t) arrlen);
"""

def correlate_simd(ht, st, qt):
    htilde = _np.array(ht.data, copy = False).view(dtype = float32)
    stilde = _np.array(st.data, copy = False).view(dtype = float32)
    qtilde = _np.array(qt.data, copy = False).view(dtype = float32)
    arrlen = len(htilde)
    inline(corr_simd_code, ['htilde', 'stilde', 'qtilde', 'arrlen'],
           extra_compile_args = ['-march=native -O3 -w'],
           #extra_compile_args = ['-mno-avx -mno-sse2 -mno-sse3 -mno-ssse3 -mno-sse4 -mno-sse4.1 -mno-sse4.2 -mno-sse4a -O2 -w'],
           #extra_compile_args = ['-msse3 -O3 -w'],
           support_code = corr_support, auto_downcast = 1)

corr_parallel_code = """
ccorrf_parallel(htilde, stilde, qtilde, (int64_t) arrlen, (int64_t) segsize);
"""
# We need a segment size (number of complex elements) such that *three* segments
# of that size will fit in the L2 cache. We also want it to be a power of two.
# We are dealing with single-precision complex numbers, which each require 8 bytes.

# Our kernel is written to assume a complex correlation of single-precision vectors,
# so that's all we support here.  Note that we are assuming that the correct target
# is that the vectors should fit in L2 cache.  Figuring out cache topology dynamically
# is a harder problem than we attempt to solve here.

if pycbc.opt.HAVE_GETCONF:
    # Since we need 3 vectors fitting in L2 cache, divide by 3
    # We find the nearest power-of-two that fits, and the length
    # of the single-precision complex array that fits into that size.
    pow2 = int( _np.log( pycbc.opt.LEVEL2_CACHE_SIZE/3.0 )/_np.log(2.0) )
    default_segsize = pow(2, pow2)/ _np.dtype( _np.complex64).itemsize
else:
    # Seems to work for Sandy Bridge/Ivy Bridge/Haswell, for now?
    default_segsize = 8192

def correlate_parallel(ht, st, qt):
    htilde = _np.array(ht.data, copy = False)
    stilde = _np.array(st.data, copy = False)
    qtilde = _np.array(qt.data, copy = False)
    arrlen = len(htilde)
    segsize = default_segsize
    inline(corr_parallel_code, ['htilde', 'stilde', 'qtilde', 'arrlen', 'segsize'],
           extra_compile_args = ['-march=native -O3 -w'] + omp_flags, libraries = omp_libs,
           #extra_compile_args = ['-mno-avx -mno-sse2 -mno-sse3 -mno-ssse3 -mno-sse4 -mno-sse4.1 -mno-sse4.2 -mno-sse4a -O2 -w'],
           #extra_compile_args = ['-msse3 -O3 -w'],
           support_code = corr_support, auto_downcast = 1)


