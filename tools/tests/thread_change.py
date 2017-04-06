#!/usr/bin/env python

# Copyright (C) 2015 Joshua Willis
#
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



import argparse
from pycbc import scheme
import weave

test_support_code = """
#include <omp.h>
#include <stdio.h>
"""

test_code = """
  int tid, nthreads;

#pragma omp parallel private(nthreads, tid)
  {
   tid = omp_get_thread_num();
   if (tid == 0){
     nthreads = omp_get_num_threads();
     printf("Total number of threads is %d\\n", nthreads);
    }
  }
"""

def print_total_threads():
    weave.inline(test_code, [], extra_compile_args = ['-fopenmp'],
                 support_code = test_support_code, libraries = ['gomp'])

parser = argparse.ArgumentParser(usage='',
    description="Test changing number of threads in pycbc")

scheme.insert_processing_option_group(parser)
opt = parser.parse_args()
scheme.verify_processing_options(opt, parser)
ctx = scheme.from_cli(opt)

print_total_threads()

with ctx:
    print_total_threads()

new_ctx = scheme.CPUScheme(num_threads = 4)

with new_ctx:
    print_total_threads()
