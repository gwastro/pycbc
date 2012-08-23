# Copyright (C) 2012  Tito Dal Canton
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
'''
These are the unittests for the pycbc PSD module.
'''

import sys
import os
import tempfile
import pycbc
import pycbc.psd
import unittest
import numpy
import optparse

_parser = optparse.OptionParser()

def _check_scheme(option, opt_str, scheme, parser):
    if scheme == 'cuda' and not pycbc.HAVE_CUDA:
        raise optparse.OptionValueError("CUDA not found")
    if scheme == 'opencl' and not pycbc.HAVE_OPENCL:
        raise optparse.OptionValueError("OpenCL not found")
    setattr(parser.values, option.dest, scheme)

_parser.add_option('--scheme', '-s', action='callback', type='choice',
    choices=('cpu', 'cuda', 'opencl'), default='cpu', dest='scheme',
    callback=_check_scheme,
    help='specifies processing scheme, can be cpu [default], cuda, or opencl')

_parser.add_option('--device-num', '-d', action='store', type='int',
    dest='devicenum', default=0,
    help='specifies a GPU device to use for CUDA or OpenCL, 0 by default')

(_options, _args) = _parser.parse_args()

if _options.scheme == 'cuda':
    _context = pycbc.scheme.CUDAScheme(device_num=_options.devicenum)
elif _options.scheme == 'opencl':
    _context = pycbc.scheme.OpenCLScheme(device_num=_options.devicenum)
elif _options.scheme == 'cpu':
    _context = pycbc.scheme.DefaultScheme()

class TestPSD(unittest.TestCase):
    def setUp(self):
        self.psd_len = 1024
        self.psd_delta_f = 0.1
        self.psd_low_freq_cutoff = 10.
    
    def test_analytical(self):
        with _context:
            psd_list = pycbc.psd.analytical.get_list()
            self.assertTrue(psd_list)
            for psd_name in psd_list:
                psd = pycbc.psd.analytical.from_string(psd_name, self.psd_len,
                                    self.psd_delta_f, self.psd_low_freq_cutoff)
                psd_min = psd.min()
                self.assertTrue(psd_min >= 0,
                                          msg=(psd_name + ': negative values'))
                self.assertTrue(psd.min() < 1e-40,
                                msg=(psd_name + ': unreasonably high minimum'))

    def test_read(self):
        test_data = numpy.zeros((self.psd_len, 2))
        test_data[:, 0] = numpy.linspace(0.,
                           (self.psd_len - 1) * self.psd_delta_f, self.psd_len)
        test_data[:, 1] = numpy.sqrt(test_data[:, 0])
        file_desc, file_name = tempfile.mkstemp()
        os.close(file_desc)
        numpy.savetxt(file_name, test_data)
        test_data[test_data[:, 0] < self.psd_low_freq_cutoff, 1] = 0.
        with _context:
            psd = pycbc.psd.read.from_txt(file_name, self.psd_len,
                                    self.psd_delta_f, self.psd_low_freq_cutoff)
            self.assertAlmostEqual(abs(psd - test_data[:, 1] ** 2).max(), 0)
        os.unlink(file_name)

suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestPSD))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    
    NotImpErrors = 0
    for error in results.errors:
        for errormsg in error:
            if type(errormsg) is str:
                if 'NotImplemented' in errormsg:
                    NotImpErrors +=1
                    break
    if results.wasSuccessful():
        sys.exit(0)
    elif len(results.failures)==0 and len(results.errors)==NotImpErrors:
        sys.exit(1)
    else:
        sys.exit(2)
