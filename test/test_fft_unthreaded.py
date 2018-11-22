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

import pycbc.fft
import unittest
from utils import parse_args_all_schemes, simple_exit
from fft_base import _BaseTestFFTClass

_scheme, _context = parse_args_all_schemes("FFT")

# Most of the work is now done in fft_base.  Below are factories for
# creating a test for each backend of each scheme.

# Get our list of backends:

backends = pycbc.fft.get_backend_names()

FFTTestClasses = []
for backend in backends:
    # This creates, for each backend, a new class derived from
    # both _BaseTestFFTClass and unittest.TestCase, and with
    # the additional property 'self.backend' set to the value
    # of backend.  One such class for each backend is appended
    # to the list
    kdict = {'backends' : [backend], 'scheme' : _scheme, 'context' : _context}
    klass = type('{0}_{1}_test'.format(_scheme,backend),
                 (_BaseTestFFTClass,),kdict)
    FFTTestClasses.append(klass)

# Finally, we create suites and run them

if __name__ == '__main__':

    suite = unittest.TestSuite()
    for klass in FFTTestClasses:
        suite.addTest(unittest.TestLoader().loadTestsFromTestCase(klass))

    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
