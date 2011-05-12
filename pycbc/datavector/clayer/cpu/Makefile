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

# patching the proxy class definition of the swig-gcc output in order to inherit a swigged class
# from a python base class:
#PATCHSCRIPT = inherit_swigged_class_from_py_base.py
#PARENTCLASS = real_vector_t
#BASECLASSMODULE = datavector_base
#BASECLASS = DataVectorBase

all: datavectorcpu.i datavectorcpu_types.h datavectorcpu_prototypes.h datavectorcpu.c
	$(SWIG) -python datavectorcpu.i
	$(PYTHON) setup.py build_ext --inplace
	#$(PYTHON) $(PATCHSCRIPT) $(PARENTCLASS) $(BASECLASSMODULE) $(BASECLASS) \
	#				< datavectorcpu.py > ../../datavectorcpu.py
	cp datavectorcpu.py ../../
	cp _datavectorcpu.so ../../

tidy:
	$(RM) -f *~
clean: tidy
	$(RM) -rf datavectorcpu_wrap.c datavectorcpu.py datavectorcpu.pyc _datavectorcpu.so build
 
