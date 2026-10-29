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

import copy
import unittest

import numpy
from scipy import special
from utils import simple_exit

from pycbc.catalog import Merger
from pycbc.conversions import mass1_from_mtotal_q, mass2_from_mtotal_q
from pycbc.detector import Detector
from pycbc.distributions import (
    CosAngle,
    JointDistribution,
    SinAngle,
    Uniform,
    UniformAngle,
)
from pycbc.filter import highpass, match, resample_to_delta_t, sigma
from pycbc.frame import read_frame
from pycbc.inference import models
from pycbc.io import get_file
from pycbc.noise import noise_from_psd
from pycbc.psd import aLIGOZeroDetHighPower, interpolate, inverse_spectrum_truncation
from pycbc.waveform import get_fd_waveform
from pycbc.waveform.waveform import FailedWaveformError


class TestModels(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        ###### Get data for references analysis of 170817
        m = Merger("GW170817")
        ifos = ["H1", "V1", "L1"]
        cls.psds = {}
        cls.data = {}

        for ifo in ifos:
            print(f"Processing {ifo} data")

            # Download the gravitational wave data for GW170817
            # Not using dcc.ligo.org because it sometimes hangs in the test
            # suite. But original URL is commented below
            # url = "https://dcc.ligo.org/public/0146/P1700349/001/"
            url = "https://github.com/gwastro/pycbc_data/releases/download/ci-data/"
            url += "{}-{}1_LOSC_CLN_4_V1-1187007040-2048.gwf"
            fname = get_file(url.format(ifo[0], ifo[0]), cache=True)
            ts = read_frame(
                fname,
                f"{ifo}:LOSC-STRAIN",
                start_time=int(m.time - 260),
                end_time=int(m.time + 40),
            )
            ts = highpass(ts, 15.0)
            ts = resample_to_delta_t(ts, 1.0 / 2048)
            ts = ts.time_slice(m.time - 112, m.time + 16)
            cls.data[ifo] = ts.to_frequencyseries()

            psd = interpolate(ts.psd(4), ts.delta_f)
            psd = inverse_spectrum_truncation(
                psd,
                int(4 * psd.sample_rate),
                trunc_method="hann",
                low_frequency_cutoff=20.0,
            )
            cls.psds[ifo] = psd

        cls.static = {
            "mass1": 1.3757,
            "mass2": 1.3757,
            "f_lower": 20.0,
            "approximant": "TaylorF2",
            "polarization": 0,
            "ra": 3.44615914,
            "dec": -0.40808407,
            "tc": 1187008882.42840,
        }

        cls.variable = [
            "distance",
            "inclination",
        ]
        cls.flow = {"H1": 25, "L1": 25, "V1": 25}
        inclination_prior = SinAngle(inclination=None)
        distance_prior = Uniform(distance=(10, 100))
        tc_prior = Uniform(tc=(m.time - 0.1, m.time + 0.1))
        pol = UniformAngle(polarization=None)
        cls.prior = JointDistribution(cls.variable, inclination_prior, distance_prior)

        # set up for marginalized polarization tests
        cls.static2 = cls.static.copy()
        cls.static2.pop("polarization")
        cls.variable2 = cls.variable + ["polarization"]
        cls.prior2 = JointDistribution(
            cls.variable2, inclination_prior, distance_prior, pol
        )

        # set up for gated model tests
        cls.static3 = cls.static.copy()
        cls.static3["t_gate_start"] = cls.static3["tc"]
        cls.static3["t_gate_end"] = cls.static3["tc"] + 0.05
        cls.variable3 = cls.variable
        cls.prior3 = JointDistribution(cls.variable3, inclination_prior, distance_prior)

        ###### Expected answers
        # Answer taken from marginalized gaussian model
        cls.q1 = {"distance": 42.0, "inclination": 2.5}
        cls.a1 = 541.8235746138382

        # answer taken from brute marginize pol + phase
        cls.a2 = 542.581
        cls.pol_samples = 200

        # answer with gate applied, no normalization
        cls.a3 = -1246.1948739646468

    def test_base_phase_marg(self):
        model = models.MarginalizedPhaseGaussianNoise(
            self.variable,
            copy.deepcopy(self.data),
            low_frequency_cutoff=self.flow,
            psds=self.psds,
            static_params=self.static,
            prior=self.prior,
        )
        model.update(**self.q1)
        self.assertAlmostEqual(self.a1, model.loglr, delta=1e-3)

    def test_relative_phase_marg(self):
        model = models.Relative(
            self.variable,
            copy.deepcopy(self.data),
            low_frequency_cutoff=self.flow,
            psds=self.psds,
            static_params=self.static,
            prior=self.prior,
            fiducial_params={"mass1": 1.3756},
            epsilon=0.1,
        )
        model.update(**self.q1)
        self.assertAlmostEqual(self.a1, model.loglr, delta=0.002)

    def test_single_phase_marg(self):
        model = models.SingleTemplate(
            self.variable,
            copy.deepcopy(self.data),
            low_frequency_cutoff=self.flow,
            psds=self.psds,
            static_params=self.static,
            prior=self.prior,
        )
        model.update(**self.q1)
        self.assertAlmostEqual(self.a1, model.loglr, delta=0.02)

    def test_single_pol_phase_marg(self):
        model = models.SingleTemplate(
            self.variable2,
            copy.deepcopy(self.data),
            low_frequency_cutoff=self.flow,
            psds=self.psds,
            static_params=self.static2,
            prior=self.prior2,
            marginalize_vector_samples=1000,
            marginalize_vector_params="polarization",
        )
        model.update(**self.q1)
        self.assertAlmostEqual(self.a2, model.loglr, delta=0.04)

    def test_gated_gaussian_psd_opts(self):
        model_toeplitz = models.GatedGaussianNoise(
            self.variable3,
            copy.deepcopy(self.data),
            low_frequency_cutoff=self.flow,
            psds=self.psds,
            static_params=self.static3,
            prior=self.prior3,
        )
        model_matmul = models.GatedGaussianNoise(
            self.variable3,
            copy.deepcopy(self.data),
            low_frequency_cutoff=self.flow,
            psds=self.psds,
            static_params=self.static3,
            prior=self.prior3,
            paint_method="matmul",
        )
        model_toeplitz.update(**self.q1)
        model_matmul.update(**self.q1)
        # check likelihoods match calculated
        self.assertAlmostEqual(self.a3, model_toeplitz.loglr, delta=0.01)
        self.assertAlmostEqual(self.a3, model_matmul.loglr, delta=0.01)
        # check paint method is being set correctly
        self.assertEqual("toeplitz", model_toeplitz.paint_method)
        self.assertEqual("matmul", model_matmul.paint_method)

    def test_brute_pol_phase_marg(self):
        # Uses the old polarization syntax untill we decide to remove it.
        # Untill then, this also tests that that interface stays working.
        model = models.BruteParallelGaussianMarginalize(
            self.variable,
            data=copy.deepcopy(self.data),
            low_frequency_cutoff=self.flow,
            psds=self.psds,
            static_params=self.static2,
            prior=self.prior,
            marginalize_phase=400,
            cores=1,
            base_model="marginalized_polarization",
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
            "approximant": "IMRPhenomD",
            "mass1": 40.0,
            "mass2": 40.0,
            "polarization": 0,
            "ra": 3.44615914,
            "dec": -0.40808407,
            "tc": tc,
            "distance": 100.0,
            "inclination": 2.5,
        }
        cls.variable = ["spin1z", "f_lower"]
        ifos = ["H1", "L1", "V1"]
        # generate the reference psd
        seglen = 8
        delta_f = 1.0 / seglen
        sample_rate = 4096
        delta_t = 1.0 / sample_rate
        flen = int(sample_rate * seglen / 2) + 1
        psd = aLIGOZeroDetHighPower(flen, delta_f, flow)
        # put non-zero values in the beginning and end of the psd
        # so the gating models will work
        psd[0 : int(flow / delta_f + 1)] = psd[int(flow / delta_f + 1)]
        psd[-2:] = psd[-2]
        seed = 1000
        cls.flow = {"H1": flow, "L1": flow, "V1": flow}
        # generate the noise
        for ifo in ifos:
            tsamples = int(seglen * sample_rate)
            ts = noise_from_psd(tsamples, delta_t, psd, seed=seed)
            ts._epoch = cls.static["tc"] - seglen / 2
            seed += 1027
            cls.data[ifo] = ts.to_frequencyseries()
            cls.psds[ifo] = psd
        # setup priors
        spin_prior = Uniform(spin1z=(-1.0, 2.0))
        flowbad = 4000.0
        flower_prior = Uniform(f_lower=(flow, flowbad + 100.0))
        pol = UniformAngle(polarization=None)
        cls.prior = JointDistribution(cls.variable, spin_prior, flower_prior)

        # set up for marginalized polarization tests
        cls.static2 = cls.static.copy()
        cls.static2.pop("polarization")
        cls.variable2 = cls.variable + ["polarization"]
        cls.prior2 = JointDistribution(cls.variable2, spin_prior, flower_prior, pol)
        # set up gated parameters
        staticgate = cls.static.copy()
        staticgate["t_gate_start"] = tc - 0.05
        staticgate["t_gate_end"] = tc
        cls.staticgate = staticgate
        # margpol
        staticgate2 = cls.static2.copy()
        staticgate2["t_gate_start"] = tc - 0.05
        staticgate2["t_gate_end"] = tc
        cls.staticgate2 = staticgate2
        # the parameters to test:
        # these parameters should pass
        cls.pass_params = {"spin1z": 0.0, "f_lower": flow}
        # these parameters should trigger a NoWaveformError
        cls.nowf_params = {"spin1z": 0.0, "f_lower": flowbad}
        # these parameters should cause a FailedWaveformError
        cls.fail_params = {"spin1z": 2.0, "f_lower": flow}

    def _run_tests(
        self,
        model,
        check_pass=True,
        check_nowf=True,
        check_failed=True,
        check_raises=True,
    ):
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
            self.variable,
            copy.deepcopy(self.data),
            low_frequency_cutoff=self.flow,
            psds=self.psds,
            static_params=self.static,
            prior=self.prior,
            ignore_failed_waveforms=True,
        )
        self._run_tests(model)

    def test_relative_phase_marg(self):
        model = models.Relative(
            self.variable,
            copy.deepcopy(self.data),
            low_frequency_cutoff=self.flow,
            psds=self.psds,
            static_params=self.static,
            prior=self.prior,
            fiducial_params={},
            # fiducial_params = {'mass1': 40.},
            epsilon=0.1,
            ignore_failed_waveforms=True,
        )
        # relative model doesn't respect flower, so no point in testing nowf
        self._run_tests(model, check_nowf=False)

    def test_brute_pol_phase_marg(self):
        # Uses the old polarization syntax untill we decide to remove it.
        # Untill then, this also tests that that interface stays working.
        model = models.BruteParallelGaussianMarginalize(
            self.variable,
            data=copy.deepcopy(self.data),
            low_frequency_cutoff=self.flow,
            psds=self.psds,
            static_params=self.static2,
            prior=self.prior,
            marginalize_phase=4,
            cores=1,
            base_model="marginalized_polarization",
            ignore_failed_waveforms=True,
        )
        # we need to do the check raises test separately because the underlying
        # base model's ignore_failed_waveforms needs to be set
        self._run_tests(model, check_raises=False)
        model = models.BruteParallelGaussianMarginalize(
            self.variable,
            data=copy.deepcopy(self.data),
            low_frequency_cutoff=self.flow,
            psds=self.psds,
            static_params=self.static2,
            prior=self.prior,
            marginalize_phase=4,
            cores=1,
            base_model="marginalized_polarization",
            ignore_failed_waveforms=False,
        )
        self._run_tests(
            model,
            check_pass=False,
            check_nowf=False,
            check_failed=False,
            check_raises=True,
        )

    def test_gated_gaussian_noise(self):
        model = models.GatedGaussianNoise(
            self.variable,
            data=copy.deepcopy(self.data),
            low_frequency_cutoff=self.flow,
            psds=self.psds,
            static_params=self.staticgate,
            prior=self.prior,
            ignore_failed_waveforms=True,
        )
        self._run_tests(model)

    def test_gated_gaussian_margpol(self):
        model = models.GatedGaussianMargPol(
            self.variable,
            data=copy.deepcopy(self.data),
            low_frequency_cutoff=self.flow,
            psds=self.psds,
            static_params=self.staticgate2,
            prior=self.prior,
            ignore_failed_waveforms=True,
        )
        self._run_tests(model)


class TestMarginalizedPolModels(unittest.TestCase):
    """
    Tests that marginalized polarization models return the same
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
            "approximant": "IMRPhenomD",
            "mass1": 40.0,
            "mass2": 40.0,
            "ra": 3.44615914,
            "dec": -0.40808407,
            "tc": tc,
            "distance": 100.0,
            "inclination": 2.5,
            "f_lower": flow,
        }
        cls.marg_variable = ["spin1z"]
        cls.orig_variable = cls.marg_variable + ["polarization"]
        ifos = ["H1", "L1", "V1"]
        seglen = 8
        delta_f = 1.0 / seglen
        sample_rate = 4096
        delta_t = 1.0 / sample_rate
        flen = int(sample_rate * seglen / 2) + 1
        psd = aLIGOZeroDetHighPower(flen, delta_f, flow)
        psd[0 : int(flow / delta_f + 1)] = psd[int(flow / delta_f + 1)]
        psd[-2:] = psd[-2]
        seed = 1000
        cls.flow = {"H1": flow, "L1": flow, "V1": flow}
        for ifo in ifos:
            tsamples = int(seglen * sample_rate)
            ts = noise_from_psd(tsamples, delta_t, psd, seed=seed)
            ts._epoch = cls.static["tc"] - seglen / 2
            seed += 1027
            cls.data[ifo] = ts.to_frequencyseries()
            cls.psds[ifo] = psd
        # Set up the static parameters for the gated models
        staticgate = cls.static.copy()
        staticgate["t_gate_start"] = tc - 0.5
        staticgate["t_gate_end"] = tc
        cls.staticgate = staticgate

    def _test_models(self, margpol_model, orig_model, polarization_samples):
        """
        Tests that manual integration of the loglikelihood over
        polarizations gives the same result as the marginalized model.
        """
        # Get the loglikelihood from the margpol model
        params = {"spin1z": 0.0}
        margpol_model.update(**params)
        margpol_logl = margpol_model.loglikelihood
        # Numerically integrate the loglikelihood over polarizations
        logls = numpy.zeros(len(polarization_samples))
        norm = numpy.log(len(logls))
        for ii, pol in enumerate(polarization_samples):
            params["polarization"] = pol
            orig_model.update(**params)
            logls[ii] = orig_model.loglikelihood
        integrated_logl = special.logsumexp(logls) - norm
        # Assert that the loglikelihoods are close
        self.assertAlmostEqual(
            margpol_logl,
            integrated_logl,
            places=2,
            msg=f"Loglikelihood mismatch: MargPol={margpol_logl}, "
            f"Integrated={integrated_logl}",
        )

    def test_gaussian_models(self):
        """Tests that the GaussianNoise models are consistent."""
        margpol_model = models.MarginalizedPolarization(
            self.marg_variable,
            data=copy.deepcopy(self.data),
            low_frequency_cutoff=self.flow,
            psds=self.psds,
            static_params=self.static,
            ignore_failed_waveforms=True,
        )
        orig_model = models.GaussianNoise(
            self.orig_variable,
            data=copy.deepcopy(self.data),
            low_frequency_cutoff=self.flow,
            psds=self.psds,
            static_params=self.static,
            ignore_failed_waveforms=True,
        )
        # now test
        polsamples = margpol_model.marginalize_vector_params["polarization"]
        self._test_models(margpol_model, orig_model, polsamples)

    def _margpol_model(self, **kwargs):
        return models.MarginalizedPolarization(
            self.marg_variable,
            data=copy.deepcopy(self.data),
            low_frequency_cutoff=self.flow,
            psds=self.psds,
            static_params=self.static,
            ignore_failed_waveforms=True,
            **kwargs,
        )

    def test_margpol_samples(self):
        """
        Tests that the number of polarization samples in the
        MarginalizedPolarization model is set correctly.
        """
        model = self._margpol_model()
        self.assertEqual(len(model.marginalize_vector_params["polarization"]), 1000)
        model = self._margpol_model(marginalize_vector_samples=50)
        self.assertEqual(len(model.marginalize_vector_params["polarization"]), 50)
        # the removed argument is rejected rather than silently ignored
        with self.assertRaises(ValueError):
            self._margpol_model(polarization_samples=60)

    def test_gated_models(self):
        """Tests that the Gated models are consistent."""
        # Initialize the GatedGaussianMargPol model
        margpol_model = models.GatedGaussianMargPol(
            self.marg_variable,
            data=copy.deepcopy(self.data),
            low_frequency_cutoff=self.flow,
            psds=self.psds,
            static_params=self.staticgate,
            ignore_failed_waveforms=True,
        )
        # Initialize the GatedGaussianNoise model
        orig_model = models.GatedGaussianNoise(
            self.orig_variable,
            data=copy.deepcopy(self.data),
            low_frequency_cutoff=self.flow,
            psds=self.psds,
            static_params=self.staticgate,
            ignore_failed_waveforms=True,
        )
        # now test
        polsamples = margpol_model.pol
        self._test_models(margpol_model, orig_model, polsamples)

    def test_marginalize_vector_samples(self):
        """
        Test that marginalized model matches brute force Gaussian marg
        when specifying a non-default number of polarization samples.
        """
        nsamples = 100
        margpol_model = models.MarginalizedPolarization(
            self.marg_variable,
            data=copy.deepcopy(self.data),
            low_frequency_cutoff=self.flow,
            psds=self.psds,
            static_params=self.static,
            ignore_failed_waveforms=True,
            marginalize_vector_samples=nsamples,
        )
        orig_model = models.GaussianNoise(
            self.orig_variable,
            data=copy.deepcopy(self.data),
            low_frequency_cutoff=self.flow,
            psds=self.psds,
            static_params=self.static,
            ignore_failed_waveforms=True,
        )
        polsamples = margpol_model.marginalize_vector_params["polarization"]
        self.assertEqual(len(polsamples), nsamples)
        self.assertEqual(margpol_model.vsamples, nsamples)
        # now test
        self._test_models(margpol_model, orig_model, polsamples)



class TestBruteTotalMassMarginalize(unittest.TestCase):
    """Tests the brute force total mass marginalization model.

    The model rescales a single reference waveform onto a grid of total
    masses instead of regenerating one per grid point, so these tests check
    the rescaling against explicitly regenerated waveforms, and the resulting
    marginalized likelihood against an explicit sum over the same grid.
    """

    @classmethod
    def setUpClass(cls):
        cls.tc = 1187008882.42840
        cls.flow = 20.0
        seglen = 8
        sample_rate = 2048
        delta_f = 1. / seglen
        delta_t = 1. / sample_rate
        flen = int(sample_rate * seglen / 2) + 1

        # grid of total masses, and the reference the waveform is scaled from
        cls.mtotal_grid = [45., 50., 55., 60., 65., 70., 75.]
        cls.fiducial_mtotal = 60.
        cls.injected_mtotal = 55.
        cls.injected_q = 1.4

        psd = aLIGOZeroDetHighPower(flen, delta_f, cls.flow)
        psd[0:int(cls.flow / delta_f + 1)] = psd[int(cls.flow / delta_f + 1)]
        psd[-2:] = psd[-2]

        cls.static = {
            'approximant': "IMRPhenomD",
            'f_lower': cls.flow,
            'spin1z': 0.,
            'spin2z': 0.,
            'ra': 3.44615914,
            'dec': -0.40808407,
            'tc': cls.tc,
            'distance': 500.,
            'inclination': 0.5,
            'polarization': 0.3,
        }
        cls.variable = ['q', 'mtotal']
        cls.ifos = ['H1', 'L1']
        cls.low_frequency_cutoff = {ifo: cls.flow for ifo in cls.ifos}

        # noise, with a signal at the injected total mass added in
        hp, _ = get_fd_waveform(
            approximant=cls.static['approximant'],
            mass1=mass1_from_mtotal_q(cls.injected_mtotal, cls.injected_q),
            mass2=mass2_from_mtotal_q(cls.injected_mtotal, cls.injected_q),
            spin1z=0., spin2z=0., inclination=cls.static['inclination'],
            distance=cls.static['distance'], delta_f=delta_f,
            f_lower=cls.flow)
        hp.resize(flen)

        cls.data = {}
        cls.psds = {}
        seed = 1000
        for ifo in cls.ifos:
            ts = noise_from_psd(int(seglen * sample_rate), delta_t, psd,
                                seed=seed)
            ts._epoch = cls.tc - seglen / 2
            seed += 1027
            fs = ts.to_frequencyseries()
            fs += hp
            cls.data[ifo] = fs
            cls.psds[ifo] = psd

        cls.prior = JointDistribution(
            cls.variable,
            Uniform(q=(1., 2.)),
            Uniform(mtotal=(min(cls.mtotal_grid), max(cls.mtotal_grid))))

    def _marg_model(self, **kwargs):
        opts = dict(
            data=copy.deepcopy(self.data),
            psds=self.psds,
            low_frequency_cutoff=self.low_frequency_cutoff,
            static_params=self.static,
            prior=self.prior,
            base_model='marginalized_time',
            marginalize_mtotal=True,
            mtotal_grid=self.mtotal_grid,
            fiducial_mtotal=self.fiducial_mtotal,
            marginalize_vector_samples=64,
            cores=1,
        )
        opts.update(kwargs)
        return models.BruteTotalMassMarginalize(self.variable, **opts)

    def _base_model(self):
        """The base model, set up to take the component masses directly so
        the waveform can be regenerated at each point of the grid."""
        variable = ['mass1', 'mass2']
        prior = JointDistribution(
            variable,
            Uniform(mass1=(1., 200.)),
            Uniform(mass2=(1., 200.)))
        return models.MarginalizedTime(
            variable, data=copy.deepcopy(self.data), psds=self.psds,
            low_frequency_cutoff=self.low_frequency_cutoff,
            static_params=self.static, prior=prior,
            marginalize_vector_samples=64)

    def _regenerated_loglrs(self, base, q):
        """loglr at each point of the grid, regenerating the waveform."""
        logls = numpy.zeros(len(self.mtotal_grid))
        for ii, mtotal in enumerate(self.mtotal_grid):
            base.update(mass1=mass1_from_mtotal_q(mtotal, q),
                        mass2=mass2_from_mtotal_q(mtotal, q))
            logls[ii] = base.loglr
        return logls

    def test_model_is_registered(self):
        """The model is retrievable by name."""
        self.assertIs(models.get_model('marginalized_mtotal'),
                      models.BruteTotalMassMarginalize)

    def test_waveform_rescaling(self):
        """A rescaled reference waveform matches one generated directly at
        the target total mass."""
        model = self._marg_model()
        q = self.injected_q
        delta_f = self.data[self.ifos[0]].delta_f
        flen = len(self.data[self.ifos[0]])

        # reference generated the way the model does it: low enough in
        # frequency that the lightest grid point still covers the band
        kappa_min = min(self.mtotal_grid) / self.fiducial_mtotal
        href, _ = get_fd_waveform(
            approximant=self.static['approximant'],
            mass1=mass1_from_mtotal_q(self.fiducial_mtotal, q),
            mass2=mass2_from_mtotal_q(self.fiducial_mtotal, q),
            spin1z=0., spin2z=0., inclination=self.static['inclination'],
            distance=self.static['distance'], delta_f=delta_f,
            f_lower=self.flow * kappa_min)
        href.resize(flen)

        for mtotal in self.mtotal_grid:
            scaled, _ = model.scale_waveform(
                href, href, self.fiducial_mtotal, mtotal)
            direct, _ = get_fd_waveform(
                approximant=self.static['approximant'],
                mass1=mass1_from_mtotal_q(mtotal, q),
                mass2=mass2_from_mtotal_q(mtotal, q),
                spin1z=0., spin2z=0.,
                inclination=self.static['inclination'],
                distance=self.static['distance'], delta_f=delta_f,
                f_lower=self.flow)
            direct.resize(flen)

            # the rescaled waveform must agree in both shape and amplitude
            mm, _ = match(scaled, direct, psd=self.psds[self.ifos[0]],
                          low_frequency_cutoff=self.flow)
            self.assertGreater(
                mm, 0.99,
                msg="match between rescaled and regenerated waveform is "
                    f"only {mm} at mtotal={mtotal}")

            s_scaled = sigma(scaled, psd=self.psds[self.ifos[0]],
                             low_frequency_cutoff=self.flow)
            s_direct = sigma(direct, psd=self.psds[self.ifos[0]],
                             low_frequency_cutoff=self.flow)
            self.assertAlmostEqual(
                s_scaled / s_direct, 1.0, places=2,
                msg="amplitude of rescaled waveform is wrong at "
                    f"mtotal={mtotal}")

    def test_loglr_matches_regeneration(self):
        """The marginalized loglr equals an explicit sum over the grid using
        regenerated waveforms."""
        model = self._marg_model()
        base = self._base_model()
        norm = numpy.log(len(self.mtotal_grid))

        for q in [1.0, self.injected_q, 1.9]:
            model.update(q=q, mtotal=self.fiducial_mtotal)
            marg_loglr = model.loglr

            integrated = special.logsumexp(
                self._regenerated_loglrs(base, q)) - norm

            self.assertAlmostEqual(
                marg_loglr, integrated, delta=0.5,
                msg=f"marginalized loglr {marg_loglr} does not match explicit "
                    f"integration {integrated} at q={q}")

    def test_maxl_mtotal_matches_grid(self):
        """The recorded maximum likelihood total mass is the grid point that
        explicit regeneration also picks out."""
        model = self._marg_model()
        base = self._base_model()
        q = self.injected_q

        model.update(q=q, mtotal=self.fiducial_mtotal)
        model.loglr

        logls = self._regenerated_loglrs(base, q)
        expected = self.mtotal_grid[int(numpy.argmax(logls))]
        self.assertEqual(model.current_stats['maxl_mtotal'], expected)

    def test_grid_from_prior_bounds(self):
        """The grid can be built from the mtotal prior instead of given."""
        model = self._marg_model(mtotal_grid=None, mtotal_grid_num=7)
        numpy.testing.assert_allclose(model.mtotal_grid, self.mtotal_grid)
        # and without a number of points there is nothing to build a grid from
        with self.assertRaises(ValueError):
            self._marg_model(mtotal_grid=None, mtotal_grid_num=None)

    def test_marginalization_off(self):
        """With the marginalization off the base model loglr is returned.

        The base model generates its own waveform in this case, so it is set
        up with the component masses as its variable parameters.
        """
        variable = ['mass1', 'mass2']
        prior = JointDistribution(
            variable,
            Uniform(mass1=(1., 200.)),
            Uniform(mass2=(1., 200.)))
        model = models.BruteTotalMassMarginalize(
            variable, data=copy.deepcopy(self.data), psds=self.psds,
            low_frequency_cutoff=self.low_frequency_cutoff,
            static_params=self.static, prior=prior,
            base_model='marginalized_time', marginalize_mtotal=False,
            marginalize_vector_samples=64, cores=1)

        self.assertIsNone(model.mtotal_grid)
        model.update(mass1=mass1_from_mtotal_q(self.injected_mtotal,
                                               self.injected_q),
                     mass2=mass2_from_mtotal_q(self.injected_mtotal,
                                               self.injected_q))
        self.assertTrue(numpy.isfinite(model.loglr))


class TestMarginalizedTimeRefFrame(unittest.TestCase):
    """
    Tests that the marginalized_time model gives the same loglr whether
    tc is sampled in the geocentric frame or in a detector's frame.
    """

    @classmethod
    def setUpClass(cls):
        cls.tc = 1126259462.420
        cls.flow = 20.0
        cls.inj = {
            "ra": 2.2,
            "dec": -1.25,
            "polarization": 1.75,
            "distance": 400.0,
            "mass1": 37.0,
            "mass2": 32.0,
            "inclination": 0.3,
            "coa_phase": 0.0,
        }
        cls.ifos = ["H1", "L1"]
        seglen = 8
        sample_rate = 2048
        tsamples = seglen * sample_rate
        flen = tsamples // 2 + 1
        delta_f = 1.0 / seglen
        psd = aLIGOZeroDetHighPower(flen, delta_f, cls.flow)
        hp, hc = get_fd_waveform(
            approximant="IMRPhenomD", delta_f=delta_f, f_lower=cls.flow, **cls.inj
        )
        hp.resize(flen)
        hc.resize(flen)
        start = int(cls.tc) - seglen / 2
        cls.data = {}
        cls.psds = {}
        for i, ifo in enumerate(cls.ifos):
            ts = noise_from_psd(tsamples, 1.0 / sample_rate, psd, seed=1000 + i)
            ts.start_time = start
            det = Detector(ifo)
            fp, fc = det.antenna_pattern(
                cls.inj["ra"], cls.inj["dec"], cls.inj["polarization"], cls.tc
            )
            tdet = det.arrival_time(cls.tc, cls.inj["ra"], cls.inj["dec"])
            h = (fp * hp + fc * hc).cyclic_time_shift(tdet - start)
            d = ts.to_frequencyseries()
            d += h.astype(d.dtype)
            d._epoch = start
            cls.data[ifo] = d
            cls.psds[ifo] = psd

    def _loglr(self, frame, marg_sky=False):
        """
        Return the mean marginalized loglr at the injected parameters,
        with tc sampled in the given frame.
        """
        numpy.random.seed(0)
        static = dict(self.inj, approximant="IMRPhenomD", f_lower=self.flow)
        static.pop("mass1")
        if frame == "geocentric":
            ref_tc = self.tc
        else:
            static["tc_ref_frame"] = frame
            ref_tc = Detector(frame).arrival_time(
                self.tc, self.inj["ra"], self.inj["dec"]
            )
        variable = ["tc", "mass1"]
        dists = [
            Uniform(tc=(ref_tc - 0.05, ref_tc + 0.05)),
            Uniform(mass1=(20.0, 50.0)),
        ]
        kwargs = {
            "marginalize_vector_params": "tc",
            "marginalize_vector_samples": 2000,
        }
        if marg_sky:
            static.pop("ra")
            static.pop("dec")
            variable += ["ra", "dec"]
            dists += [UniformAngle(ra=None), CosAngle(dec=None)]
            kwargs = {
                "marginalize_vector_params": "tc, ra, dec",
                "marginalize_vector_samples": 5000,
                "marginalize_sky_initial_samples": int(2e5),
            }
        model = models.MarginalizedTime(
            variable,
            copy.deepcopy(self.data),
            dict.fromkeys(self.ifos, self.flow),
            psds=self.psds,
            static_params=static,
            sample_rate=2048,
            prior=JointDistribution(variable, *dists),
            **kwargs,
        )
        loglrs = []
        for _ in range(3):
            model.update(mass1=self.inj["mass1"])
            loglrs.append(model.loglr)
        return numpy.mean(loglrs)

    def test_time_marg_ref_frame(self):
        """
        Test that time marginalization works regardless of the frame times
        are sampled in.
        """
        geo = self._loglr("geocentric")
        # Make sure the signal is actually recovered
        self.assertGreater(geo, 100)
        for ifo in self.ifos:
            det = self._loglr(ifo)
            msg = (
                f"loglr with tc in {ifo} frame ({det}) differs from geocentric ({geo})"
            )
            self.assertAlmostEqual(geo, det, delta=1.0, msg=msg)

    def test_sky_time_marg_ref_frame(self):
        """
        Test that sky/time marginalization works regardless of the frame times
        are sampled in.
        """
        geo = self._loglr("geocentric", marg_sky=True)
        self.assertGreater(geo, 100)
        for ifo in self.ifos:
            det = self._loglr(ifo, marg_sky=True)
            msg = (
                f"loglr with tc in {ifo} frame ({det}) differs from geocentric ({geo})"
            )
            self.assertAlmostEqual(geo, det, delta=5.0, msg=msg)


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestModels))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestWaveformErrors))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestMarginalizedPolModels))
suite.addTest(
    unittest.TestLoader().loadTestsFromTestCase(TestBruteTotalMassMarginalize)
)
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestMarginalizedTimeRefFrame))

if __name__ == "__main__":
    from astropy.utils import iers

    iers.conf.auto_download = False
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
