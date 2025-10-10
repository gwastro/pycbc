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
import numpy
from scipy import special
from pycbc.catalog import Merger
from pycbc.psd import interpolate, inverse_spectrum_truncation, aLIGOZeroDetHighPower
from pycbc.noise import noise_from_psd
from pycbc.frame import read_frame
from pycbc.filter import highpass, resample_to_delta_t
from astropy.utils.data import download_file
from pycbc.inference import models
from pycbc.distributions import Uniform, JointDistribution, SinAngle, UniformAngle
from pycbc.waveform.waveform import FailedWaveformError

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


class TestWaveformErrors(unittest.TestCase):
    """Tests that models handle no waveform errors correctly."""

    @classmethod
    def setUpClass(cls):
        cls.psds = {}
        cls.data = {}
        # static params for the test
        tc = 1187008882.42840
        flow = 20
        cls.static = {
            'approximant':"IMRPhenomD",
            'mass1': 40.,
            'mass2': 40.,
            'polarization': 0,
            'ra': 3.44615914,
            'dec': -0.40808407,
            'tc': tc,
            'distance': 100.,
            'inclination': 2.5
            }
        cls.variable = ['spin1z', 'f_lower']
        ifos = ['H1', 'L1', 'V1']
        # generate the reference psd
        seglen = 8
        delta_f = 1./seglen
        sample_rate = 4096
        delta_t = 1./sample_rate
        flen = int(sample_rate * seglen / 2) + 1
        psd = aLIGOZeroDetHighPower(flen, delta_f, flow)
        # put non-zero values in the beginning and end of the psd
        # so the gating models will work
        psd[0:int(flow/delta_f+1)] = psd[int(flow/delta_f+1)]
        psd[-2:] = psd[-2]
        seed = 1000
        cls.flow = {'H1': flow, 'L1': flow, 'V1': flow}
        # generate the noise
        for ifo in ifos:
            tsamples = int(seglen * sample_rate)
            ts = noise_from_psd(tsamples, delta_t, psd, seed=seed)
            ts._epoch = cls.static['tc'] - seglen/2
            seed += 1027
            cls.data[ifo] = ts.to_frequencyseries()
            cls.psds[ifo] = psd
        # setup priors
        spin_prior = Uniform(spin1z=(-1., 2.))
        flowbad = 4000.
        flower_prior = Uniform(f_lower=(flow, flowbad+100.))
        pol = UniformAngle(polarization=None)
        cls.prior = JointDistribution(cls.variable, spin_prior, flower_prior)

        # set up for marginalized polarization tests
        cls.static2 = cls.static.copy()
        cls.static2.pop('polarization')
        cls.variable2 = cls.variable + ['polarization']
        cls.prior2 = JointDistribution(cls.variable2, spin_prior, flower_prior,
                                       pol)
        # set up gated parameters
        staticgate = cls.static.copy()
        staticgate['t_gate_start'] = tc - 0.05
        staticgate['t_gate_end'] = tc
        cls.staticgate = staticgate
        # margpol
        staticgate2 = cls.static2.copy()
        staticgate2['t_gate_start'] = tc - 0.05
        staticgate2['t_gate_end'] = tc
        cls.staticgate2 = staticgate2
        # the parameters to test:
        # these parameters should pass
        cls.pass_params = {'spin1z': 0., 'f_lower': flow}
        # these parameters should trigger a NoWaveformError
        cls.nowf_params = {'spin1z': 0., 'f_lower': flowbad}
        # these parameters should cause a FailedWaveformError
        cls.fail_params = {'spin1z': 2., 'f_lower': flow}

    def _run_tests(self, model, check_pass=True, check_nowf=True,
                   check_failed=True, check_raises=True):
        # check that the model works
        if check_pass:
            model.update(**self.pass_params)
            self.assertTrue(numpy.isfinite(model.loglr))
        # check that a no waveform error is caught correctly
        if check_nowf:
            model.update(**self.nowf_params)
            self.assertEqual(model.loglr, -numpy.inf)
        # check that a failed waveform is caught correctly
        if check_failed:
            model.update(**self.fail_params)
            self.assertEqual(model.loglr, -numpy.inf)
        # check that an error is raised if ignore_failed_waveforms is False
        if check_raises:
            model.ignore_failed_waveforms = False
            model.update(**self.fail_params)
            with self.assertRaises((FailedWaveformError, RuntimeError)):
                model.loglr

    def test_base_phase_marg(self):
        model = models.MarginalizedPhaseGaussianNoise(
                                self.variable, copy.deepcopy(self.data),
                                low_frequency_cutoff=self.flow,
                                psds=self.psds,
                                static_params=self.static,
                                prior=self.prior,
                                ignore_failed_waveforms=True)
        self._run_tests(model)

    def test_relative_phase_marg(self):
        model = models.Relative(self.variable, copy.deepcopy(self.data),
                                 low_frequency_cutoff=self.flow,
                                 psds = self.psds,
                                 static_params = self.static,
                                 prior = self.prior,
                                 fiducial_params = {},
                                 #fiducial_params = {'mass1': 40.},
                                 epsilon = .1,
                                 ignore_failed_waveforms=True)
        # relative model doesn't respect flower, so no point in testing nowf
        self._run_tests(model, check_nowf=False)

    def test_brute_pol_phase_marg(self):
        # Uses the old polarization syntax untill we decide to remove it.
        # Untill then, this also tests that that interface stays working.
        model = models.BruteParallelGaussianMarginalize(
                        self.variable, data=copy.deepcopy(self.data),
                        low_frequency_cutoff=self.flow,
                        psds = self.psds,
                        static_params = self.static2,
                        prior = self.prior,
                        marginalize_phase=4,
                        cores=1,
                        base_model='marginalized_polarization',
                        ignore_failed_waveforms=True
                        )
        # we need to do the check raises test separately because the underlying
        # base model's ignore_failed_waveforms needs to be set
        self._run_tests(model, check_raises=False)
        model = models.BruteParallelGaussianMarginalize(
                        self.variable, data=copy.deepcopy(self.data),
                        low_frequency_cutoff=self.flow,
                        psds = self.psds,
                        static_params = self.static2,
                        prior = self.prior,
                        marginalize_phase=4,
                        cores=1,
                        base_model='marginalized_polarization',
                        ignore_failed_waveforms=False
                        )
        self._run_tests(model, check_pass=False, check_nowf=False,
                        check_failed=False, check_raises=True)

    def test_gated_gaussian_noise(self):
        model = models.GatedGaussianNoise(
            self.variable, data=copy.deepcopy(self.data),
            low_frequency_cutoff=self.flow,
            psds=self.psds,
            static_params=self.staticgate,
            prior=self.prior,
            ignore_failed_waveforms=True)
        self._run_tests(model)

    def test_gated_gaussian_margpol(self):
        model = models.GatedGaussianMargPol(
            self.variable, data=copy.deepcopy(self.data),
            low_frequency_cutoff=self.flow,
            psds=self.psds,
            static_params=self.staticgate2,
            prior=self.prior,
            ignore_failed_waveforms=True)
        self._run_tests(model)


class TestMarginalizedPolModels(unittest.TestCase):
    """Tests that marginalized polarization models return the same
    loglikelihood as the unmarginalized model numerically integrated over the
    same set of polarization samples.
    """

    @classmethod
    def setUpClass(cls):
        # Use the same setup as TestWaveformErrors
        cls.psds = {}
        cls.data = {}
        tc = 1187008882.42840
        flow = 20
        cls.static = {
            'approximant': "IMRPhenomD",
            'mass1': 40.,
            'mass2': 40.,
            'ra': 3.44615914,
            'dec': -0.40808407,
            'tc': tc,
            'distance': 100.,
            'inclination': 2.5,
            'f_lower': flow
        }
        cls.marg_variable = ['spin1z']
        cls.orig_variable = cls.marg_variable + ['polarization']
        ifos = ['H1', 'L1', 'V1']
        seglen = 8
        delta_f = 1. / seglen
        sample_rate = 4096
        delta_t = 1. / sample_rate
        flen = int(sample_rate * seglen / 2) + 1
        psd = aLIGOZeroDetHighPower(flen, delta_f, flow)
        psd[0:int(flow / delta_f + 1)] = psd[int(flow / delta_f + 1)]
        psd[-2:] = psd[-2]
        seed = 1000
        cls.flow = {'H1': flow, 'L1': flow, 'V1': flow}
        for ifo in ifos:
            tsamples = int(seglen * sample_rate)
            ts = noise_from_psd(tsamples, delta_t, psd, seed=seed)
            ts._epoch = cls.static['tc'] - seglen / 2
            seed += 1027
            cls.data[ifo] = ts.to_frequencyseries()
            cls.psds[ifo] = psd
        # Set up the static parameters for the gated models
        staticgate = cls.static.copy()
        staticgate['t_gate_start'] = tc - 0.5
        staticgate['t_gate_end'] = tc
        cls.staticgate = staticgate

    def _test_models(self, margpol_model, orig_model, polarization_samples):
        """Tests that manual integration of the loglikelihood over
        polarizations gives the same result as the marginalized model."""
        # Get the loglikelihood from the margpol model
        params = {'spin1z': 0.}
        margpol_model.update(**params)
        margpol_logl = margpol_model.loglikelihood
        # Numerically integrate the loglikelihood over polarizations
        logls = numpy.zeros(len(polarization_samples))
        norm = numpy.log(len(logls))
        for ii,pol in enumerate(polarization_samples):
            params['polarization'] = pol
            orig_model.update(**params)
            logls[ii] = orig_model.loglikelihood
        integrated_logl = special.logsumexp(logls) - norm
        # Assert that the loglikelihoods are close
        self.assertAlmostEqual(
            margpol_logl, integrated_logl, places=2,
            msg=f"Loglikelihood mismatch: MargPol={margpol_logl}, "
                f"Integrated={integrated_logl}"
        ) 

    def test_gaussian_models(self):
        """Tests that the GaussianNoise models are consistent."""
        margpol_model = models.MarginalizedPolarization(
            self.marg_variable, data=copy.deepcopy(self.data),
            low_frequency_cutoff=self.flow,
            psds=self.psds,
            static_params=self.static,
            ignore_failed_waveforms=True
        )
        orig_model = models.GaussianNoise(
            self.orig_variable, data=copy.deepcopy(self.data),
            low_frequency_cutoff=self.flow,
            psds=self.psds,
            static_params=self.static,
            ignore_failed_waveforms=True
        )
        # now test
        polsamples = margpol_model.marginalize_vector_params['polarization']
        self._test_models(margpol_model, orig_model, polsamples)

    def test_gated_models(self):
        """Tests that the Gated models are consistent."""
        # Initialize the GatedGaussianMargPol model
        margpol_model = models.GatedGaussianMargPol(
            self.marg_variable, data=copy.deepcopy(self.data),
            low_frequency_cutoff=self.flow,
            psds=self.psds,
            static_params=self.staticgate,
            ignore_failed_waveforms=True
        )
        # Initialize the GatedGaussianNoise model
        orig_model = models.GatedGaussianNoise(
            self.orig_variable, data=copy.deepcopy(self.data),
            low_frequency_cutoff=self.flow,
            psds=self.psds,
            static_params=self.staticgate,
            ignore_failed_waveforms=True
        )
        # now test
        polsamples = margpol_model.pol
        self._test_models(margpol_model, orig_model, polsamples)


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestModels))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestWaveformErrors))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestMarginalizedPolModels))

if __name__ == '__main__':
    from astropy.utils import iers
    iers.conf.auto_download = False
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
