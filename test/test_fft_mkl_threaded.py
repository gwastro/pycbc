# Copyright (C) 2012  Josh Willis, Andrew Miller
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

#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#
"""
These are the unit-tests for the pycbc.fft subpackage, testing only unthreaded
backends for the various schemes.
"""
import platform
import unittest

import pycbc.fft
from pycbc.scheme import CPUScheme
from sys import exit as _exit
from utils import parse_args_cpu_only, simple_exit
from fft_base import _BaseTestFFTClass

parse_args_cpu_only("MKL threaded backend")

# See if we can get set the FFTW backend to 'openmp'; if not, say so and exit.
if 'arm64' in platform.machine():
    print("MKL not supported on arm64, skipping")
    pass
elif not 'mkl' in pycbc.fft.get_backend_names():
    print("MKL does not seem to be available; why isn't it installed?")
    _exit(0)
else:
    # Now set the number of threads to something nontrivial

    # Most of the work is now done in fft_base.

    FFTTestClasses = []
    for num_threads in [2,4,6,8]:
        kdict = {'backends' : ['mkl'],
                 'scheme' : 'cpu',
                 'context' : CPUScheme(num_threads=num_threads)}
        klass = type('FFTW_OpenMP_test',
                     (_BaseTestFFTClass,),kdict)
        klass.__test__ = True
        FFTTestClasses.append(klass)


    # Finally, we create suites and run them

    if __name__ == '__main__':

        suite = unittest.TestSuite()
        for klass in FFTTestClasses:
            suite.addTest(unittest.TestLoader().loadTestsFromTestCase(klass))

        results = unittest.TextTestRunner(verbosity=2).run(suite)
        simple_exit(results)
