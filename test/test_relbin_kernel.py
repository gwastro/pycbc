# Copyright (C) 2026  Alex Nitz
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
These are the unittests for the relative binning likelihood kernels
in pycbc.inference.models.relbin_cpu.

Each kernel computes <d|h> and <h|h> for a bank of marginalization
samples. Each test checks a kernel against direct(), which computes the
same sums from their definition. The kernels sum the same terms in a
different order, so agreement is to about 1e-13 rather than exact.
"""

import os
import unittest

import numpy

from pycbc.inference.models.relbin import setup_bins
from pycbc.inference.models.relbin_cpu import (likelihood_parts,
                                               likelihood_parts_v_pol,
                                               likelihood_parts_v_pol_time,
                                               likelihood_parts_v_time,
                                               likelihood_parts_vector,
                                               likelihood_parts_vectorp,
                                               likelihood_parts_vectort)

# a 32 second segment from 20 Hz, as used for a neutron star binary
FLOW, FHIGH, DF = 20.0, 1024.0, 1. / 32.
NSAMPLE = 64

# relative tolerance against direct()
TOL = 1e-12


def get_seed():
    """ Use PYCBC_VALIDATION_SEED if set, otherwise a random seed
    """
    seed = os.environ.get('PYCBC_VALIDATION_SEED')
    if seed is None:
        seed = numpy.random.randint(0, 2 ** 31)
    return int(seed)


def inspiral(freqs, mchirp):
    """ Stationary phase inspiral with an f^(-7/6) amplitude
    """
    return ((freqs / FLOW) ** (-7. / 6.) * numpy.exp(
        1.0j * 3. / 128. / (numpy.pi * mchirp * 4.9e-6 * freqs) ** (5. / 3.)))


def summary(h1, h2, freqs, psd, bins):
    """ Calculate the summary values for the inner product <h1|h2>
    """
    h12 = numpy.conjugate(h1) * h2 / psd
    a0 = numpy.array([4.0 * DF * h12[lo:hi].sum() for lo, hi in bins])
    a1 = numpy.array([4.0 / (hi - lo)
                      * (h12[lo:hi] * (freqs[lo:hi] - freqs[lo])).sum()
                      for lo, hi in bins])
    return a0, a1


def make_case(seed):
    """ Build one set of kernel arguments
    """
    rng = numpy.random.default_rng(seed)

    fine = numpy.arange(FLOW, FHIGH + DF, DF)
    psd = 1e-46 * ((fine / 100.) ** -4.14 + 1.0 + (fine / 300.) ** 2)

    # bins as the model lays them down for this band
    index = setup_bins(fine, FLOW, FHIGH)
    edges = fine[index]
    bins = list(zip(index[:-1], index[1:]))

    # fiducial waveform, and data holding a signal slightly off it
    h00f = 1e-23 * inspiral(fine, 1.2)
    noise = ((rng.normal(size=fine.size) + 1.0j * rng.normal(size=fine.size))
             * numpy.sqrt(psd))
    data = 1.05 * 1e-23 * inspiral(fine, 1.2005) + noise

    a0, a1 = summary(h00f, data, fine, psd, bins)
    b0, b1 = summary(h00f, h00f, fine, psd, bins)

    # a small second component keeps the polarizations from being a
    # pure phase apart
    cosi = 0.6
    dom = 1e-23 * inspiral(edges, 1.2005)
    sub = 1e-24 * inspiral(edges, 1.35)
    hp = 0.5 * (1.0 + cosi ** 2) * dom + 0.30 * sub
    hc = 1.0j * cosi * dom + 0.22 * numpy.exp(0.7j) * sub

    return dict(freqs=edges, hp=hp, hc=hc, h00=1e-23 * inspiral(edges, 1.2),
                a0=a0, a1=a1, b0=numpy.abs(b0), b1=numpy.abs(b1),
                fp=rng.uniform(-1, 1, NSAMPLE),
                fc=rng.uniform(-1, 1, NSAMPLE),
                dtc=rng.uniform(-0.1, 0.1, NSAMPLE))


def direct(freqs, fp, fc, dtc, *, hp, hc, h00, a0, a1, b0, b1):
    """ Calculate <d|h> and <h|h> directly from their definition
    """
    r = numpy.exp(-2.0j * numpy.pi * dtc * freqs) * (fp * hp + fc * hc) / h00
    x = numpy.abs(r) ** 2
    return (numpy.conjugate((a0 * r[:-1] + a1 * numpy.diff(r)).sum()),
            (b0 * x[:-1] + b1 * numpy.diff(x)).sum())


class TestRelbinKernel(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.seed = get_seed()
        print('using PYCBC_VALIDATION_SEED=%s' % cls.seed)
        cls.case = make_case(cls.seed)
        # with earth rotation fp, fc and the arrival are per bin
        rng = numpy.random.RandomState(cls.seed + 1)
        n = len(cls.case['freqs'])
        cls.band = dict(fp=rng.uniform(-1, 1, n), fc=rng.uniform(-1, 1, n),
                        times=rng.uniform(-0.01, 0.01, n),
                        dtc=rng.uniform(-0.01, 0.01, n),
                        pol=numpy.exp(2.0j * rng.uniform(0, numpy.pi,
                                                         NSAMPLE)))

    def check(self, out, effective, label):
        """ Compare kernel output to direct()

        effective(j) returns the fp, fc and arrival for sample j, either
        as scalars or as per bin arrays.
        """
        c = self.case
        hdv, hhv = out
        err = numpy.empty((NSAMPLE, 2))
        for j in range(NSAMPLE):
            hd, hh = direct(c['freqs'], *effective(j), hp=c['hp'], hc=c['hc'],
                            h00=c['h00'], a0=c['a0'], a1=c['a1'],
                            b0=c['b0'], b1=c['b1'])
            err[j] = (abs(hdv[j] - hd) / abs(hd), abs(hhv[j] - hh) / abs(hh))
        # numpy.max, not the builtin, so a nan does not pass
        worst = err.max()
        print('%s: worst relative difference from the definition %.2e'
              % (label, worst))
        self.assertTrue(worst < TOL)

    def rotated(self, j):
        """ Rotate (fp, fc) by the polarization of sample j
        """
        b = self.band
        cp, sp = b['pol'][j].real, b['pol'][j].imag
        return b['fp'] * cp - b['fc'] * sp, b['fp'] * sp + b['fc'] * cp

    def test_the_scalar_kernel_is_the_definition(self):
        """ Check direct() against the unmodified scalar kernel
        """
        c = self.case
        self.check((numpy.array([likelihood_parts(
                        c['freqs'], c['fp'][j], c['fc'][j], c['dtc'][j],
                        c['hp'], c['hc'], c['h00'], c['a0'], c['a1'],
                        c['b0'], c['b1'])[i] for j in range(NSAMPLE)])
                    for i in (0, 1)),
                   lambda j: (c['fp'][j], c['fc'][j], c['dtc'][j]),
                   'likelihood_parts')

    def test_vector(self):
        """ Sky or time points, response and shift vary per sample
        """
        c = self.case
        self.check(likelihood_parts_vector(
                       c['freqs'], c['fp'], c['fc'], c['dtc'], c['hp'],
                       c['hc'], c['h00'], c['a0'], c['a1'], c['b0'], c['b1']),
                   lambda j: (c['fp'][j], c['fc'][j], c['dtc'][j]),
                   'likelihood_parts_vector')

    def test_vectort(self):
        """ Time points only, so <h|h> is one value for the bank
        """
        c = self.case
        fp, fc = float(c['fp'][0]), float(c['fc'][0])
        out = likelihood_parts_vectort(
            c['freqs'], fp, fc, c['dtc'], c['hp'], c['hc'], c['h00'],
            c['a0'], c['a1'], c['b0'], c['b1'])
        self.assertLess(out[1].max() - out[1].min(),
                        abs(out[1].mean()) * TOL, 'a shift changed <h|h>')
        self.check(out, lambda j: (fp, fc, c['dtc'][j]),
                   'likelihood_parts_vectort')

    def test_vectorp(self):
        """ Polarization points only, the shift is scalar
        """
        c = self.case
        dtc = float(c['dtc'][0])
        self.check(likelihood_parts_vectorp(
                       c['freqs'], c['fp'], c['fc'], dtc, c['hp'], c['hc'],
                       c['h00'], c['a0'], c['a1'], c['b0'], c['b1']),
                   lambda j: (c['fp'][j], c['fc'][j], dtc),
                   'likelihood_parts_vectorp')

    def test_v_pol(self):
        """ Earth rotation, polarization points
        """
        c, b = self.case, self.band
        self.check(likelihood_parts_v_pol(
                       c['freqs'], b['fp'], b['fc'], b['dtc'], b['pol'],
                       c['hp'], c['hc'], c['h00'], c['a0'], c['a1'],
                       c['b0'], c['b1']),
                   lambda j: self.rotated(j) + (b['dtc'],),
                   'likelihood_parts_v_pol')

    def test_v_time(self):
        """ Earth rotation, time points
        """
        c, b = self.case, self.band
        self.check(likelihood_parts_v_time(
                       c['freqs'], b['fp'], b['fc'], b['times'], c['dtc'],
                       c['hp'], c['hc'], c['h00'], c['a0'], c['a1'],
                       c['b0'], c['b1']),
                   lambda j: (b['fp'], b['fc'], b['times'] + c['dtc'][j]),
                   'likelihood_parts_v_time')

    def test_v_pol_time(self):
        """ Earth rotation, polarization and time points
        """
        c, b = self.case, self.band
        self.check(likelihood_parts_v_pol_time(
                       c['freqs'], b['fp'], b['fc'], b['times'], c['dtc'],
                       b['pol'], c['hp'], c['hc'], c['h00'], c['a0'],
                       c['a1'], c['b0'], c['b1']),
                   lambda j: self.rotated(j) + (b['times'] + c['dtc'][j],),
                   'likelihood_parts_v_pol_time')


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestRelbinKernel))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    from utils import simple_exit
    simple_exit(results)
