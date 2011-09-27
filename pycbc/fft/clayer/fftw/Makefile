# Copyright (C) 2011 Josh Willis
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
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
# datavector local makefile for pycbc

RM = /bin/rm
PYTHON = python
SWIG = swig

all: fft_fftw.i fft_fftw.h fft_fftw.c fft_fftw_private.h
	$(SWIG) -python fft_fftw.i
	$(PYTHON) setup.py build_ext --inplace
	cp fft_fftw.py ../../
	cp _fft_fftw.so ../../

clean: 
	$(RM) -rf fft_fftw_wrap.c fft_fftw.py ../../fft_fftw.py \
        fft_fftw.pyc _fft_fftw.so ../../_fft_fftw.so ../../fft_fftw.pyc  build 
 
