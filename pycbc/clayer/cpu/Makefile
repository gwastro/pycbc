# Copyright (C) 2011 Karsten Wiesner
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

all: pycbccpu.i pycbccpu_types.h pycbccpu_prototypes.h pycbccpu.c
	$(SWIG) -python pycbccpu.i
	$(PYTHON) setup.py build_ext --inplace
	cp pycbccpu.py ../../
	cp _pycbccpu.so ../../

tidy:
	$(RM) -f *~
clean: tidy
	$(RM) -rf pycbccpu_wrap.c pycbccpu.py ../../pycbccpu.py \
        pycbccpu.pyc _pycbccpu.so ../../_pycbccpu.so build 
 
