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

all: fftcpu.i fftcpu.h fftcpu.c
	$(SWIG) -python fftcpu.i
	$(PYTHON) setup.py build_ext --inplace
	cp fftcpu.py ../../
	cp _fftcpu.so ../../

tidy:
	$(RM) -f *~
clean: tidy
	$(RM) -rf fftcpu_wrap.c fftcpu.py ../../fftcpu.py \
        fftcpu.pyc _fftcpu.so ../../_fftcpu.so ../../fftcpu.pyc  build 
 
