/*
 *  Copyright (C) 2014 Josh Willis
 *
 *  This program is free software; you can redistribute it and/or modify it
 *  under the terms of the GNU General Public License as published by the
 *  Free Software Foundation; either version 3 of the License, or (at your
 *  option) any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
 *  Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc.,
 *  51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */

/*
 * This module defines several C functions that compute
 * a combined thresholding and time clustering of a complex array using a
 * multithreaded, SIMD vectoried code.
 *
 * There are four C functions defined:
 *
 * parallel_threshold: A multithreaded function that will identify points in a
 *                     timeseries above a user-input threshold.
 *
 * max_simple: A single-threaded function that does NOT use explicit 
 *             SIMD vectorization,
 *             to compute the maximum of a complex time series over a given
 *             window.
 *
 * windowed_max: A single threaded (but not vectorized) function that finds the
 *               locations, norms, and complex values of the maxima in each of
 *               a set of predefined windows in a given complex array. It does
 *               this by calling max_simple on each window with the
 *               appropriate parameters.
 *
 * parallel_thresh_cluster: A multithreaded function that finds the maxima in
 *                          each of a set of contiguous, fixed-length windows
 *                          by calling windowed_max in parallel. It then sweeps
 *                          through the results of that (single-threaded) and
 *                          tests for above threshold, and time-clusters the
 *                          surviving triggers.
 *
 * A user calls only the last function; the other three exist to conveniently
 * compartmentalize SIMD code from OpenMP code.
 */

#ifndef __clang__
    #include <omp.h> 
#endif
#include <stdint.h> // For uint32_t, int64_t
#include <complex>
#include <math.h> // For M_SQRT2
#include <cstdio>
#include <cstring>


/* Rough approx of GCC's error function. */
void error(int status, int errnum, const char *format) {
  fprintf(stderr, format);
  if (status != 0) {
    exit(status);
  }
}

void _parallel_threshold(int64_t N, std::complex<float> * __restrict arr,
                         std::complex<float> * __restrict outv,
                         uint32_t * __restrict outl,
                         uint32_t * __restrict count,
                         const float v) {
        unsigned int num_parallel_regions = 16;
        unsigned int t=0;

        #pragma omp parallel for ordered shared(t)
        for (unsigned int p=0; p<num_parallel_regions; p++){
            unsigned int start  = (N * p) / num_parallel_regions;
            unsigned int end    = (N * (p+1)) / num_parallel_regions;
            unsigned int c = 0;

            for (unsigned int i=start; i<end; i++){
                float r = std::real(arr[i]);
                float im = std::imag(arr[i]);
                if ((r * r + im * im) > v){
                    outl[c+start] = i;
                    outv[c+start] = std::complex<float>(r, im);
                    c++;
                }
            }

            #pragma omp ordered
            {
                t+=c;
            }
            memmove(outl+t-c, outl+start, sizeof(unsigned int)*c);
            memmove(outv+t-c, outv+start, sizeof(std::complex<float>)*c);

        }

        count[0] = t;
}
                        

/*
 * The following maximizes over an interval that can be no longer than
 * the clustering window.  We do *not* assume alignment,
 * because the window length
 * itself may not be a multiple of the alignment.
 */

void _max_simple(float * __restrict inarr, float * __restrict mval,
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


void _windowed_max(std::complex<float> * __restrict inarr,
                   const int64_t arrlen,
                   std::complex<float> * __restrict cvals,
                   float * __restrict norms,
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
   requires *real* arrays, and lengths will be converted where appropriate.

  */

  int64_t i, nwindows;

  nwindows = ( (arrlen % winsize) ? (arrlen/winsize) + 1 : (arrlen/winsize) );

  // Everything but the last window, which may not be full length

  for (i = 0; i < nwindows-1; i++){
    // The factor of 2 multiplying lengths[i] is because max_simd needs its length as a real
    // length, not complex.  But startpts and startoffset are complex values.
    _max_simple((float *) &inarr[i*winsize], (float *) &cvals[i],
                &norms[i], &locs[i], startoffset + i*winsize, 2*winsize);
  }
  // Now the last window (which will be the only window if arrlen <= winzise)
  _max_simple((float *) &inarr[i*winsize], (float *) &cvals[i],
               &norms[i], &locs[i], startoffset + i*winsize, 2*(arrlen - i*winsize));

  return;
}

int _parallel_thresh_cluster(std::complex<float> * __restrict inarr,
                             const uint32_t arrlen,
                             std::complex<float> * __restrict values,
                             uint32_t * __restrict locs,
                             const float thresh, const uint32_t winsize,
                             const uint32_t segsize){

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

  int64_t i, nsegs, nwins_ps, last_arrlen, last_nwins_ps, outlen, true_segsize;
  int64_t *seglens, *mlocs, cnt, s_segsize, s_arrlen, s_winsize;
  float *norms, thr_sqr;
  std::complex<float> *cvals;

  thr_sqr = (thresh * thresh);

  // Signed versions of all our array length inputs, so loop
  // logic and calls to other functions are safe.

  s_segsize = (int64_t) segsize;
  s_arrlen = (int64_t) arrlen;
  s_winsize = (int64_t) winsize;

  // If the length of a window is less than the length of a segment, then the
  // number of windows per segment is segsize/winsize (which rounds *down* if
  // it does not divide evenly). If not, then the number of windows per segment
  // is one, since we will then always give a window-sized chunks to each core,
  // regardless of segsize.

  if (s_winsize < s_segsize) {
    nwins_ps =  s_segsize / s_winsize;
  } else {
    nwins_ps = 1;
  }

  // Now the actual segment size that we will use is the number of windows per
  // segment times the window size; this formula is true whether winsize is
  // larger or smaller than segsize.

  true_segsize = nwins_ps * s_winsize;

  // We need to allocate our temporary arrays that hold the results of clustering
  // over each window. So we just divide the array length by the window size, and
  // if this does not divide evenly, then the last window can be shorter.

  outlen = ( (s_arrlen % s_winsize) ? (s_arrlen/s_winsize) + 1 : (s_arrlen/s_winsize) );

  // If true_segsize divides arrlen evenly, then the number of segments 'nsegs' is that
  // ratio; if not, it is one more than the floor of that ratio and the last segment
  // will be shorter than all of the others.

  nsegs = ( (s_arrlen % true_segsize) ? (s_arrlen/true_segsize) + 1 : (s_arrlen/true_segsize) );

  // The amount of the initial array left over for the last segment could be different

  last_arrlen = s_arrlen - true_segsize * (nsegs -1);

  // Now dynamic allocation.  No reason really to align this memory; it will be parceled
  // out to different cores anyway.

  cvals = (std::complex<float> *) malloc(outlen * sizeof(std::complex<float>) );
  norms = (float *) malloc(outlen * sizeof(float) );
  mlocs = (int64_t *) malloc(outlen * sizeof(int64_t) );

  // The next array allows us to dynamically communicate possibly changed sizes to the
  // many parallel calls to windowed_max:

  seglens = (int64_t *) malloc(nsegs * sizeof(int64_t) );

  // check to see if anything failed
  if ( (cvals == NULL) || (norms == NULL) || (mlocs == NULL) || (seglens == NULL) ){
    error(EXIT_FAILURE, ENOMEM, "Could not allocate temporary memory needed by parallel_thresh_cluster");
  }

  for (i = 0; i < (nsegs-1); i++){
    seglens[i] = true_segsize;
  }
  seglens[i] = last_arrlen;

  // Now the real work, in an OpenMP parallel for loop:
#pragma omp parallel for schedule(dynamic,1)
  for (i = 0; i < nsegs; i++){
    _windowed_max(&inarr[i*true_segsize], seglens[i], &cvals[i*nwins_ps],
                 &norms[i*nwins_ps], &mlocs[i*nwins_ps],
                 s_winsize, i*true_segsize);
  }

  /*

  We have now the requisite maxima over our windows in cvals,
  norms, and mlocs. We want to apply the threshold and cluster over
  time.

  To do this we compare each element to the one before and after it,
  taking special care for the first and last elements.

  */

  cnt = 0;

  // Handle the first element, including the case that
  // the whole array might consist of just one window.

  if (outlen > 1) {
    if ( (norms[0] > norms[1]) && (norms[0] > thr_sqr) ) {
      cnt++;
      values[cnt-1] = cvals[0];
      locs[cnt-1] = (uint32_t) mlocs[0];
    }
  } else {
    if ( norms[0] > thr_sqr ) {
      cnt++;
      values[cnt-1] = cvals[0];
      locs[cnt-1] = (uint32_t) mlocs[0];
    }
  }

  // Now loop through the second through next to last
  // elements, comparing to the elements before and after,
  // and to the threshold.

  for (i = 1; i < (outlen-1); i++) {
    if ( (norms[i] > thr_sqr) && (norms[i] > norms[i-1]) && (norms[i] >= norms[i+1]) ) {
      cnt++;
      values[cnt-1] = cvals[i];
      locs[cnt-1] = (uint32_t) mlocs[i];
    }
  }

  // Finally, handle the last element. We test for whether
  // (outlen > 1), but unlike the first element, we do nothing
  // if (outlen == 1), because that will have been already
  // handled when we looked at the first element.

  if (outlen > 1) {
    if ( (norms[outlen-1] > norms[outlen-2]) && (norms[outlen-1] > thr_sqr) ) {
      cnt++;
      values[cnt-1] = cvals[outlen-1];
      locs[cnt-1] = (uint32_t) mlocs[outlen-1];
    }
  }

  free(cvals);
  free(norms);
  free(mlocs);
  free(seglens);

  return cnt;
}

