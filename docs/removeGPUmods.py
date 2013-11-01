#!/usr/bin/env python

# Copyright (C) 2011 Ian W. Harry
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

"""
PLEASE FIX: This code is here because sphinx-apidoc does not have the
capability to exclude specific files from the documentation. The GPU modules
must be removed as we cannot import their stuff. This feature has been
requested in sphinx, remove this code and use it if it gets added.
We may need to write a patch for Sphinx to add this functionality.
"""

import os,sys
import glob

# Can these be named somewhat more obviously!
excludes=['cuda', 'opencl', 'cufft', 'cuda_pyfft', 'cl_pyfft',\
          'pycbc_phenomC_tmplt', 'TaylorF2']

fileList = glob.glob('pycbc.*.rst')

for file in fileList:
    output = []
    fp = open(file,'r')
    addLine = True
    for line in fp:
        if (':mod:' in line) and ('Module' in line):
            for excludeNam in excludes:
                if excludeNam in line:
                    addLine=False
                    break
            else:
                addLine=True
        if addLine:
            output.append(line)
    fp.close()
    fp = open(file,'w')
    fp.writelines(output)
    fp.close()
                 
