# Copyright (C) 2012  Alex Nitz, Josh Willis
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
These are the unittests for the pycbc.waveform module
"""
import unittest
import numpy
from pycbc.types import *
from pycbc.scheme import *
from utils import parse_args_all_schemes, simple_exit

_scheme, _context = parse_args_all_schemes("correlate")

from pycbc.vetoes.chisq_cpu import chisq_accum_bin_numpy
from pycbc.vetoes import chisq_accum_bin, power_chisq_bins, power_chisq
from pycbc.vetoes import power_chisq_at_points_from_precomputed
from pycbc.filter import resample_to_delta_t
from pycbc.catalog import Merger
from pycbc.psd import interpolate, inverse_spectrum_truncation
from pycbc.waveform import get_fd_waveform
from pycbc.filter import matched_filter_core


trusted_accum = chisq_accum_bin_numpy

class TestChisq(unittest.TestCase):
    def setUp(self,*args):
        self.context = _context
        self.scheme = _scheme
        self.tolerance = 1e-6
        xr = numpy.random.uniform(low=-1, high=1.0, size=2**20)
        xi = numpy.random.uniform(low=-1, high=1.0, size=2**20)
        self.x = Array(xr + xi * 1.0j, dtype=complex64)
        self.z = zeros(2**20, dtype=float32)
        for i in range(0, 4):
            trusted_accum(self.z, self.x)

        m = Merger("GW170814")

        ifos = ['H1', 'L1', 'V1']
        data = {}
        psd = {}
        for ifo in ifos:
            # Read in and condition the data and measure PSD
            ts = m.strain(ifo).highpass_fir(15, 512)
            data[ifo] = resample_to_delta_t(ts, 1.0/2048).crop(2, 2)
            p = data[ifo].psd(2)
            p = interpolate(p, data[ifo].delta_f)
            p = inverse_spectrum_truncation(p, int(2 * data[ifo].sample_rate),
                                            low_frequency_cutoff=15.0)
            psd[ifo] = p
        hp, _ = get_fd_waveform(approximant="IMRPhenomD",
                                 mass1=31.36, mass2=31.36,
                                 f_lower=20.0, delta_f=data[ifo].delta_f)
        hp.resize(len(psd[ifo]))

        # For each ifo use this template to calculate the SNR time series
        snr = {}
        snr_unnorm = {}
        norm = {}
        corr = {}
        for ifo in ifos:
            snr_unnorm[ifo], corr[ifo], norm[ifo] = \
                matched_filter_core(hp, data[ifo], psd=psd[ifo],
                                    low_frequency_cutoff=20)
            snr[ifo] = snr_unnorm[ifo] * norm[ifo]

        self.snr = snr
        self.snr_unnorm = snr_unnorm
        self.norm = norm
        self.corr = corr
        self.hp = hp
        self.data = data
        self.psd = psd
        self.ifos = ifos

    def test_accum(self):
        with self.context:
            z = zeros(2**20, dtype=float32)
            for i in range(0, 4):
                chisq_accum_bin(z, self.x)
            self.assertTrue(self.z.almost_equal_elem(z, self.tolerance))

    def test_chisq(self):
        chisq_quick = {}
        chisq_full = {}
        chisq_ref = {}
        for ifo in self.ifos:
            nbins = 26
            dof = nbins * 2 - 2

            bins = power_chisq_bins(self.hp, nbins, self.psd[ifo],
                                    low_frequency_cutoff=20.0)
            chisq = power_chisq_at_points_from_precomputed \
                (self.corr[ifo], self.snr_unnorm[ifo][27402:27492].data,
                 self.norm[ifo], bins, indices=numpy.arange(27402,27492))
            chisq_quick[ifo] = chisq /  dof

            chisq = power_chisq(self.hp, self.data[ifo], nbins,
                                self.psd[ifo], low_frequency_cutoff=20.0)
            chisq_full[ifo] = chisq[27402:27492] / dof
            max_diff = max(abs(chisq_full[ifo] - chisq_quick[ifo]))
            self.assertTrue(max_diff < 1E-5)


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestChisq))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
