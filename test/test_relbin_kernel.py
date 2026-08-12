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
"""Tests the vectorized relative binning kernels against the scalar one.

The three vector kernels exist only to do, for a whole bank of samples at
once, what likelihood_parts does for one. They are written differently for
speed: the division by the fiducial waveform is hoisted out of the sample
loop, <h|h> is obtained from a quadratic form rather than accumulated, and
the time shift phase is evaluated in a pass of its own. None of that is
allowed to change the answer, so the property to check is the one that
defines them, that sample by sample they reproduce the scalar kernel.

The rearrangement reassociates floating point sums, so the agreement is
close but not to the last bit; the tolerance here is far tighter than any
reassociation can reach and far looser than bit equality. The scalar kernel
is in turn checked against a direct sum written straight from the
definition, so that a shared mistake could not pass unnoticed.

The inputs are a chirping inspiral seen against a nearby fiducial waveform,
on the bins the model itself would choose for it, with summary weights built
from an analytic noise curve. The two polarizations are given a small extra
component so they are not a pure phase apart: if they were, the cross term
of the quadratic form would vanish and the test would not reach it.
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

# a 32 second segment analysed from 20 Hz, as a neutron star binary is
FLOW, FHIGH, DF = 20.0, 1024.0, 1. / 32.
NSAMPLE = 64

# how closely the vector kernels must track the scalar one, relative
TOL = 1e-12

# the kernels spell pi as a truncated literal rather than taking it from
# math.h, so the reference sum below has to use the same one or it would be
# comparing two slightly different time shifts and disagreeing at 1e-8 for a
# reason that has nothing to do with what is being tested here
KERNEL_PI = 3.141592653


def get_seed():
    """Draw a fresh seed each run, or take the one being reproduced."""
    seed = os.environ.get('PYCBC_VALIDATION_SEED')
    if seed is None:
        seed = numpy.random.randint(0, 2 ** 31)
    return int(seed)


def inspiral(freqs, mchirp):
    """A stationary phase inspiral: an f^(-7/6) amplitude under a phase
    that sweeps thousands of radians across the band, which is what makes
    the time shift in the kernel a hard thing to get right."""
    return ((freqs / FLOW) ** (-7. / 6.) * numpy.exp(
        1.0j * 3. / 128. / (numpy.pi * mchirp * 4.9e-6 * freqs) ** (5. / 3.)))


def summary(h1, h2, freqs, psd, bins):
    """The constant and linear summary weights of <h1|h2> over each bin,
    formed exactly as the relative binning model forms them."""
    h12 = numpy.conjugate(h1) * h2 / psd
    a0 = numpy.array([4.0 * DF * h12[l:h].sum() for l, h in bins])
    a1 = numpy.array([4.0 / (h - l)
                      * (h12[l:h] * (freqs[l:h] - freqs[l])).sum()
                      for l, h in bins])
    return a0, a1


def make_case(seed):
    """Build one realistic argument set for the kernels."""
    rng = numpy.random.default_rng(seed)

    fine = numpy.arange(FLOW, FHIGH + DF, DF)
    psd = 1e-46 * ((fine / 100.) ** -4.14 + 1.0 + (fine / 300.) ** 2)

    # the bins the model would lay down for this band, so the edges are
    # crowded where the phase turns over quickly and sparse where it does not
    index = setup_bins(fine, FLOW, FHIGH)
    edges = fine[index]
    bins = list(zip(index[:-1], index[1:]))

    # the fiducial waveform, and data holding a signal a little off it
    h00f = 1e-23 * inspiral(fine, 1.2)
    noise = ((rng.normal(size=fine.size) + 1.0j * rng.normal(size=fine.size))
             * numpy.sqrt(psd))
    data = 1.05 * 1e-23 * inspiral(fine, 1.2005) + noise

    a0, a1 = summary(h00f, data, fine, psd, bins)
    b0, b1 = summary(h00f, h00f, fine, psd, bins)

    # a dominant quadrupole plus a small second piece, which is what stops
    # the two polarizations from being a pure phase apart
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


def direct(freqs, fp, fc, dtc, hp, hc, h00, a0, a1, b0, b1):
    """<d|h> and <h|h> written straight from the definition, with no
    reference to how the kernels are organised internally."""
    r = numpy.exp(-2.0j * KERNEL_PI * dtc * freqs) * (fp * hp + fc * hc) / h00
    x = numpy.abs(r) ** 2
    return (numpy.conjugate((a0 * r[:-1] + a1 * numpy.diff(r)).sum()),
            (b0 * x[:-1] + b1 * numpy.diff(x)).sum())


class TestRelbinKernel(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.seed = get_seed()
        print('using PYCBC_VALIDATION_SEED=%s' % cls.seed)
        cls.case = make_case(cls.seed)

    def args(self, **override):
        c = dict(self.case, **override)
        return [c[k] for k in ('freqs', 'fp', 'fc', 'dtc', 'hp', 'hc',
                               'h00', 'a0', 'a1', 'b0', 'b1')]

    def scalar(self, j):
        """What the untouched scalar kernel says about sample j."""
        c = self.case
        return likelihood_parts(c['freqs'], c['fp'][j], c['fc'][j],
                                c['dtc'][j], c['hp'], c['hc'], c['h00'],
                                c['a0'], c['a1'], c['b0'], c['b1'])

    def check_against_scalar(self, hdv, hhv, label):
        err = numpy.empty((NSAMPLE, 2))
        for j in range(NSAMPLE):
            hd, hh = self.scalar(j)
            err[j] = (abs(hdv[j] - hd) / abs(hd), abs(hhv[j] - hh) / abs(hh))
        # numpy.max, not the builtin, so that a nan cannot pass unnoticed
        worst = err.max()
        print('%s: worst relative difference from the scalar kernel %.2e'
              % (label, worst))
        self.assertTrue(worst < TOL)

    def test_scalar_kernel_is_the_definition(self):
        """Everything below is measured against the scalar kernel, so first
        confirm that it is itself the sum it claims to be."""
        c = self.case
        for j in (0, NSAMPLE // 2, NSAMPLE - 1):
            hd, hh = self.scalar(j)
            rhd, rhh = direct(c['freqs'], c['fp'][j], c['fc'][j], c['dtc'][j],
                              c['hp'], c['hc'], c['h00'], c['a0'], c['a1'],
                              c['b0'], c['b1'])
            self.assertLess(abs(hd - rhd) / abs(rhd), TOL)
            self.assertLess(abs(hh - rhh) / abs(rhh), TOL)

    def test_vector_matches_scalar(self):
        """Sky or time points: every sample has its own response and its
        own time shift."""
        hdv, hhv = likelihood_parts_vector(*self.args())
        self.check_against_scalar(hdv, hhv, 'likelihood_parts_vector')

    def test_vectort_matches_scalar(self):
        """Time points only: the response is a scalar, so <h|h> is one
        number for the whole bank and must still be the right one."""
        c = self.case
        fp, fc = float(c['fp'][0]), float(c['fc'][0])
        hdv, hhv = likelihood_parts_vectort(
            c['freqs'], fp, fc, c['dtc'], c['hp'], c['hc'], c['h00'],
            c['a0'], c['a1'], c['b0'], c['b1'])
        self.case = dict(c, fp=numpy.full(NSAMPLE, fp),
                         fc=numpy.full(NSAMPLE, fc))
        try:
            self.check_against_scalar(hdv, hhv, 'likelihood_parts_vectort')
        finally:
            self.case = c

    def test_vectorp_matches_scalar(self):
        """Polarization points only: the time shift is a scalar and is
        folded into the waveform ratios ahead of the sample loop."""
        c = self.case
        dtc = float(c['dtc'][0])
        hdv, hhv = likelihood_parts_vectorp(
            c['freqs'], c['fp'], c['fc'], dtc, c['hp'], c['hc'], c['h00'],
            c['a0'], c['a1'], c['b0'], c['b1'])
        self.case = dict(c, dtc=numpy.full(NSAMPLE, dtc))
        try:
            self.check_against_scalar(hdv, hhv, 'likelihood_parts_vectorp')
        finally:
            self.case = c

    def rotating(self):
        """The earth-rotation inputs: response and arrival vary across the
        band, so fp, fc and the arrival time are per bin rather than per
        sample. Only the polarization and the offset are per sample."""
        c = self.case
        n = len(c['freqs'])
        rng = numpy.random.RandomState(self.seed + 1)
        return dict(fpb=rng.uniform(-1, 1, n), fcb=rng.uniform(-1, 1, n),
                    times=rng.uniform(-0.01, 0.01, n),
                    dtcb=rng.uniform(-0.01, 0.01, n),
                    pol=numpy.exp(2.0j * rng.uniform(0, numpy.pi, NSAMPLE)))

    def check_against_direct(self, hdv, hhv, effective, label):
        """effective(j) gives the per-bin fp, fc and arrival for sample j."""
        c = self.case
        err = numpy.empty((NSAMPLE, 2))
        for j in range(NSAMPLE):
            fp, fc, dt = effective(j)
            hd, hh = direct(c['freqs'], fp, fc, dt, c['hp'], c['hc'],
                            c['h00'], c['a0'], c['a1'], c['b0'], c['b1'])
            err[j] = (abs(hdv[j] - hd) / abs(hd), abs(hhv[j] - hh) / abs(hh))
        worst = err.max()
        print('%s: worst relative difference from the definition %.2e'
              % (label, worst))
        self.assertTrue(worst < TOL)

    def test_v_pol_matches_the_definition(self):
        """Earth rotation, polarization points: the arrival varies across
        the band but not between samples, so the whole time shift comes out
        of the sample loop."""
        c, r = self.case, self.rotating()
        hdv, hhv = likelihood_parts_v_pol(
            c['freqs'], r['fpb'], r['fcb'], r['dtcb'], r['pol'],
            c['hp'], c['hc'], c['h00'], c['a0'], c['a1'], c['b0'], c['b1'])

        def effective(j):
            cp, sp = r['pol'][j].real, r['pol'][j].imag
            return (r['fpb'] * cp - r['fcb'] * sp,
                    r['fpb'] * sp + r['fcb'] * cp, r['dtcb'])
        self.check_against_direct(hdv, hhv, effective,
                                  'likelihood_parts_v_pol')

    def test_v_time_matches_the_definition(self):
        """Earth rotation, time points: the response is per bin and does
        not move with the sample, so <h|h> is one number for the bank."""
        c, r = self.case, self.rotating()
        hdv, hhv = likelihood_parts_v_time(
            c['freqs'], r['fpb'], r['fcb'], r['times'], c['dtc'],
            c['hp'], c['hc'], c['h00'], c['a0'], c['a1'], c['b0'], c['b1'])
        self.assertLess(hhv.max() - hhv.min(), abs(hhv.mean()) * TOL,
                        'a time shift changed <h|h>')

        def effective(j):
            return r['fpb'], r['fcb'], r['times'] + c['dtc'][j]
        self.check_against_direct(hdv, hhv, effective,
                                  'likelihood_parts_v_time')

    def test_v_pol_time_matches_the_definition(self):
        """Earth rotation, both: the per-bin arrival is folded in ahead of
        the sample loop and only the per-sample shift is left in it."""
        c, r = self.case, self.rotating()
        hdv, hhv = likelihood_parts_v_pol_time(
            c['freqs'], r['fpb'], r['fcb'], r['times'], c['dtc'], r['pol'],
            c['hp'], c['hc'], c['h00'], c['a0'], c['a1'], c['b0'], c['b1'])

        def effective(j):
            cp, sp = r['pol'][j].real, r['pol'][j].imag
            return (r['fpb'] * cp - r['fcb'] * sp,
                    r['fpb'] * sp + r['fcb'] * cp,
                    r['times'] + c['dtc'][j])
        self.check_against_direct(hdv, hhv, effective,
                                  'likelihood_parts_v_pol_time')

    def test_cross_term_is_reached(self):
        """<h|h> is a quadratic form in the antenna response, and the
        agreement above only means something if the fixture actually
        exercises its cross term. Measure how big that term is: what is
        left of <h|h> once the two pure terms are subtracted off. It is a
        fraction of a percent, which is small but is a trillion times the
        tolerance above, so getting it wrong could not go unnoticed."""
        c = self.case
        one, zero = numpy.ones(NSAMPLE), numpy.zeros(NSAMPLE)
        _, pp = likelihood_parts_vector(*self.args(fp=one, fc=zero))
        _, cc = likelihood_parts_vector(*self.args(fp=zero, fc=one))
        _, hh = likelihood_parts_vector(*self.args())
        cross = hh - c['fp'] ** 2 * pp - c['fc'] ** 2 * cc
        share = numpy.abs(cross / hh).max()
        print('cross term is up to %.2f%% of <h|h> here' % (100 * share))
        self.assertGreater(share, 1e-3)


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestRelbinKernel))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    from utils import simple_exit
    simple_exit(results)
