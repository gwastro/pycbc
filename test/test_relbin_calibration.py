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
"""Tests calibration error support in the relative binning model.

The correction is a spline through a few control points, so it varies far
more slowly with frequency than the waveform does. That is what makes it
affordable here: it can be evaluated at the bin edges alone. These check
that doing so gives the same answer as correcting the waveform at every
frequency of the data, which is what the full resolution model does.
"""

import copy
import unittest
from configparser import ConfigParser

import numpy
from scipy.interpolate import UnivariateSpline
from utils import simple_exit

from pycbc.detector import Detector
from pycbc.distributions import JointDistribution, SinAngle, Uniform
from pycbc.inference import models
from pycbc.noise import noise_from_psd
from pycbc.psd import aLIGOZeroDetHighPower
from pycbc.strain.calibration import CubicSpline, Recalibrate
from pycbc.waveform import get_td_waveform

TC = 1187008882.42840
FLOW, SEGLEN, SRATE = 25., 32, 2048
INJ = dict(mass1=1.42, mass2=1.36, distance=50., inclination=0.5,
           ra=1.7, dec=-0.4, polarization=0.3, coa_phase=1.1)
NPOINT = 4


class TestRelbinCalibration(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        flen = int(SRATE * SEGLEN / 2) + 1
        psd = aLIGOZeroDetHighPower(flen, 1. / SEGLEN, FLOW)
        hp, hc = get_td_waveform(
            approximant='TaylorF2', f_lower=FLOW, delta_t=1. / SRATE,
            **{k: INJ[k] for k in ('mass1', 'mass2', 'distance',
                                   'inclination', 'coa_phase')})
        peak = float(hp.sample_times[numpy.argmax(abs(hp.data) ** 2
                                                  + abs(hc.data) ** 2)])
        cls.data, cls.psds = {}, {}
        seed = 5
        for ifo in ['H1', 'L1']:
            ts = noise_from_psd(int(SEGLEN * SRATE), 1. / SRATE, psd,
                                seed=seed)
            seed += 97
            ts._epoch = TC - SEGLEN / 2
            signal = Detector(ifo).project_wave(
                hp, hc, INJ['ra'], INJ['dec'], INJ['polarization'])
            signal.start_time += TC - peak
            cls.data[ifo] = ts.add_into(signal).to_frequencyseries()
            cls.psds[ifo] = psd

        cls.flow = {ifo: FLOW for ifo in cls.data}
        cls.static = dict(mass1=INJ['mass1'], mass2=INJ['mass2'],
                          f_lower=FLOW, approximant='TaylorF2', tc=TC,
                          ra=INJ['ra'], dec=INJ['dec'],
                          polarization=INJ['polarization'])
        cls.variable = ['distance', 'inclination']
        cls.prior = JointDistribution(
            cls.variable, SinAngle(inclination=None),
            Uniform(distance=(10, 200)))
        cls.point = {'distance': INJ['distance'],
                     'inclination': INJ['inclination']}

    def recalibration(self):
        """Built from a config section, the way a model builds it.

        Reading the section is what decides which class the correction
        ends up being, so the tests go through it rather than naming a
        class themselves.
        """
        cp = ConfigParser()
        cp.add_section('calibration')
        for ifo in self.data:
            cp.set('calibration', '%s_model' % ifo, 'cubic_spline')
            cp.set('calibration', '%s_minimum_frequency' % ifo, str(FLOW))
            cp.set('calibration', '%s_maximum_frequency' % ifo, '1000')
            cp.set('calibration', '%s_n_points' % ifo, str(NPOINT))
        return {ifo: Recalibrate.from_config(cp, ifo, 'calibration')
                for ifo in self.data}

    def calibration_params(self, amplitude=0.05, phase=0.05):
        """A correction big enough to matter, and not flat with frequency."""
        params = {}
        for ifo in self.data:
            for i in range(NPOINT):
                sign = 1 if i % 2 == 0 else -1
                name = '%s_%d' % (ifo.lower(), i)
                params['recalib_amplitude_' + name] = sign * amplitude
                params['recalib_phase_' + name] = sign * phase
        return params

    def relbin(self, recalibration=None, epsilon=0.03):
        return models.Relative(
            list(self.variable), copy.deepcopy(self.data),
            low_frequency_cutoff=self.flow, psds=self.psds,
            static_params=self.static, prior=self.prior,
            fiducial_params={}, epsilon=epsilon,
            recalibration=recalibration)

    def exact(self, recalibration=None):
        return models.MarginalizedPhaseGaussianNoise(
            list(self.variable), copy.deepcopy(self.data),
            low_frequency_cutoff=self.flow, psds=self.psds,
            static_params=self.static, prior=self.prior,
            recalibration=recalibration)

    def test_relbin_matches_full_resolution_with_calibration(self):
        """Correcting at the bin edges must match correcting everywhere."""
        params = dict(self.point, **self.calibration_params())

        exact = self.exact(self.recalibration())
        exact.update(**params)

        model = self.relbin(self.recalibration())
        model.update(**params)

        self.assertAlmostEqual(model.loglr, exact.loglr, delta=0.05)

        # agreeing to 0.05 means nothing unless the correction is worth
        # more than that in the first place
        without = self.exact()
        without.update(**self.point)
        self.assertGreater(abs(exact.loglr - without.loglr), 1.)

    def test_it_holds_for_a_larger_correction(self):
        """A correction well beyond what is expected must still agree."""
        params = dict(self.point,
                      **self.calibration_params(amplitude=0.2, phase=0.2))

        exact = self.exact(self.recalibration())
        exact.update(**params)

        model = self.relbin(self.recalibration())
        model.update(**params)

        self.assertAlmostEqual(model.loglr, exact.loglr, delta=0.05)

    def test_the_matrix_is_dropped_where_the_spline_stops_being_linear(self):
        """The matrix stands in for a spline that is only sometimes linear.

        scipy's spline smooths rather than interpolates, and reduces to a
        linear map of its control values only while its residual
        constraint does not bind. Nothing a calibration model is asked for
        comes near that bound, but the matrix is kept between calls and
        the control values are not, so what is safe on the call that built
        it does not stay safe by itself.
        """
        # ten control points, as a real calibration model has: with four
        # the constraint never binds over any range worth trying
        npoint = 10
        recalib = CubicSpline(minimum_frequency=FLOW, maximum_frequency=1000.,
                              n_points=npoint, ifo_name='H1')
        freqs = numpy.linspace(FLOW, 1000., 512)

        # build the matrix on a correction of the expected size
        small = {'recalib_amplitude_H1_%d' % i: 0.05 * (-1) ** i
                 for i in range(npoint)}
        small.update({'recalib_phase_H1_%d' % i: 0.05 * (-1) ** i
                      for i in range(npoint)})
        ones = numpy.ones(len(freqs))
        recalib.set_params(**small)
        recalib.apply_calibration(ones, frequencies=freqs)
        recalib.apply_calibration(ones, frequencies=freqs)
        self.assertIsNotNone(recalib.spline_basis(freqs),
                             "the matrix was never built")

        # then ask for one far past where the spline is that linear map
        big = {k: 30. * v for k, v in small.items()}
        recalib.set_params(**big)
        factor = recalib.apply_calibration(ones, frequencies=freqs)

        amplitude = [recalib.params['amplitude_H1_%d' % i]
                     for i in range(npoint)]
        phase = [recalib.params['phase_H1_%d' % i] for i in range(npoint)]
        expected = UnivariateSpline(recalib.spline_points, amplitude)(freqs)
        direct = UnivariateSpline(recalib.spline_points, phase)(freqs)
        expected = ((1.0 + expected) * (2.0 + 1j * direct)
                    / (2.0 - 1j * direct))
        self.assertTrue(numpy.allclose(factor, expected, rtol=1e-12),
                        "the matrix answered where the spline is not one")


suite = unittest.TestSuite()
suite.addTest(
    unittest.TestLoader().loadTestsFromTestCase(TestRelbinCalibration))

if __name__ == '__main__':
    from astropy.utils import iers
    iers.conf.auto_download = False
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
