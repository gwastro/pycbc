// Copyright (C) 2011 Josh Willis
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at your
// option) any later version.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


//
// =============================================================================
//
//                                   Preamble
//
// =============================================================================
//
// fft plan constructor and destructor prototypes for pycbc

#ifndef FFT_H
#define FFT_H


/* The pointer to the actual plan is void because it may be either CPU or GPU */

#define FFTPLAN_STRUCT_DATA \
  void          *theplan; \
  unsigned long size; \
  int           fwdflag; \


typedef struct {
FFTPLAN_STRUCT_DATA
} fft_real_single_plan;

typedef struct {
FFTPLAN_STRUCT_DATA
} fft_real_double_plan;

typedef struct {
FFTPLAN_STRUCT_DATA
} fft_complex_single_plan;

typedef struct {
FFTPLAN_STRUCT_DATA
} fft_complex_double_plan;


#endif /* FFTCPU_H */
