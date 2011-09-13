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
OPENCL_INC = /usr/local/nvidia/sdk-3.2/OpenCL/common/inc/

all: pycbcopencl.i pycbcopencl_types.h pycbcopencl_prototypes.h pycbcopencl.c gpu_inspiral_gpuutils.h gpu_inspiral_gpuutils.c
	$(SWIG) -python pycbcopencl.i
	$(PYTHON) setup.py build_ext --inplace
	cp pycbcopencl.py ../../
	cp _pycbcopencl.so ../../

tidy:
	$(RM) -f *~

clean: tidy
	$(RM) -rf pycbcopencl_wrap.c pycbcopencl.py ../../pycbcopencl.py \
        pycbcopencl.pyc _pycbcopencl.so ../../_pycbcopencl.so build 
 
