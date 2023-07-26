# Copyright (C) 2021 Alex Nitz
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
These are the unittests for pycbc.inference.models
"""
import unittest
import copy
from utils import simple_exit
from pycbc.catalog import Merger
from pycbc.psd import interpolate, inverse_spectrum_truncation
from pycbc.frame import read_frame
from pycbc.filter import highpass, resample_to_delta_t
from astropy.utils.data import download_file
from pycbc.inference import models
from pycbc.distributions import Uniform, JointDistribution, SinAngle, UniformAngle

class TestModels(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        ###### Get data for references analysis of 170817
        m = Merger("GW170817")
        ifos = ['H1', 'V1', 'L1']
        cls.psds = {}
        cls.data = {}

        for ifo in ifos:
            print("Processing {} data".format(ifo))

            # Download the gravitational wave data for GW170817
            url = "https://dcc.ligo.org/public/0146/P1700349/001/"
            url += "{}-{}1_LOSC_CLN_4_V1-1187007040-2048.gwf"
            fname = download_file(url.format(ifo[0], ifo[0]), cache=True)
            ts = read_frame(fname, "{}:LOSC-STRAIN".format(ifo),
                            start_time=int(m.time - 260),
                            end_time=int(m.time + 40))
            ts = highpass(ts, 15.0)
            ts = resample_to_delta_t(ts, 1.0/2048)
            ts = ts.time_slice(m.time-112, m.time + 16)
            cls.data[ifo] = ts.to_frequencyseries()

            psd = interpolate(ts.psd(4), ts.delta_f)
            psd = inverse_spectrum_truncation(psd, int(4 * psd.sample_rate),
                                              trunc_method='hann',
                                              low_frequency_cutoff=20.0)
            cls.psds[ifo] = psd

        cls.static = {'mass1':1.3757,
                  'mass2':1.3757,
                  'f_lower':20.0,
                  'approximant':"TaylorF2",
                  'polarization':0,
                  'ra': 3.44615914,
                  'dec': -0.40808407,
                  'tc':  1187008882.42840,
                 }

        cls.variable = ['distance',
                    'inclination',
                      ]
        cls.flow = {'H1':25, 'L1':25, 'V1':25}
        inclination_prior = SinAngle(inclination=None)
        distance_prior = Uniform(distance=(10, 100))
        tc_prior = Uniform(tc=(m.time-0.1, m.time+0.1))
        pol = UniformAngle(polarization=None)
        cls.prior = JointDistribution(cls.variable, inclination_prior,
                                       distance_prior)

        # set up for marginalized polarization tests
        cls.static2 = cls.static.copy()
        cls.static2.pop('polarization')
        cls.variable2 = cls.variable + ['polarization']
        cls.prior2 = JointDistribution(cls.variable2, inclination_prior,
                                       distance_prior, pol)

        ###### Expected answers
        # Answer taken from marginalized gaussian model
        cls.q1 = {'distance':42.0, 'inclination':2.5}
        cls.a1 = 541.8235746138382

        # answer taken from brute marginize pol + phase
        cls.a2 = 542.581
        cls.pol_samples = 200

    def test_base_phase_marg(self):
        model = models.MarginalizedPhaseGaussianNoise(
                                self.variable, copy.deepcopy(self.data),
                                low_frequency_cutoff=self.flow,
                                psds=self.psds,
                                static_params=self.static,
                                prior=self.prior,
                                )
        model.update(**self.q1)
        self.assertAlmostEqual(self.a1, model.loglr, delta=1e-3)

    def test_relative_phase_marg(self):
        model = models.Relative(self.variable, copy.deepcopy(self.data),
                                 low_frequency_cutoff=self.flow,
                                 psds = self.psds,
                                 static_params = self.static,
                                 prior = self.prior,
                                 fiducial_params = {'mass1':1.3756},
                                 epsilon = .1,
                                )
        model.update(**self.q1)
        self.assertAlmostEqual(self.a1, model.loglr, delta=0.002)

    def test_single_phase_marg(self):
        model = models.SingleTemplate(
                        self.variable, copy.deepcopy(self.data),
                        low_frequency_cutoff=self.flow,
                        psds = self.psds,
                        static_params = self.static,
                        prior = self.prior,
                        )
        model.update(**self.q1)
        self.assertAlmostEqual(self.a1, model.loglr, delta=0.02)

    def test_single_pol_phase_marg(self):
        model = models.SingleTemplate(
                        self.variable2, copy.deepcopy(self.data),
                        low_frequency_cutoff=self.flow,
                        psds = self.psds,
                        static_params = self.static2,
                        prior = self.prior2,
                        marginalize_vector_samples = 1000,
                        marginalize_vector_params = 'polarization',
                        )
        model.update(**self.q1)
        self.assertAlmostEqual(self.a2, model.loglr, delta=0.04)

    def test_brute_pol_phase_marg(self):
        # Uses the old polarization syntax untill we decide to remove it.
        # Untill then, this also tests that that interface stays working.
        model = models.BruteParallelGaussianMarginalize(
                        self.variable, data=copy.deepcopy(self.data),
                        low_frequency_cutoff=self.flow,
                        psds = self.psds,
                        static_params = self.static2,
                        prior = self.prior,
                        marginalize_phase=400,
                        cores=1,
                        base_model='marginalized_polarization',
                        )
        model.update(**self.q1)
        self.assertAlmostEqual(self.a2, model.loglr, delta=0.002)

suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestModels))

if __name__ == '__main__':
    from astropy.utils import iers
    iers.conf.auto_download = False
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
