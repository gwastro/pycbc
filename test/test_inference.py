# Copyright (C) 2017 Christopher M. Biwer
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
These are the unittests for samplers in the pycbc.inference subpackage.
"""
import sys
import pycbc
import unittest
import numpy
from pycbc import distributions
from pycbc.inference import likelihood
from pycbc.inference import sampler
from pycbc.psd import analytical
from pycbc.waveform import generator
from utils import parse_args_cpu_only
from utils import simple_exit

# tests only need to happen on the CPU
parse_args_cpu_only("Samplers")

class TestSamplers(unittest.TestCase):

    def setUp(self, *args):

        # set random seed
        numpy.random.seed(1024)

        # set data parameters
        self.ifos = ["H1", "L1", "V1"]
        self.data_length = 4 # in seconds
        self.sample_rate = 2048 # in Hertz
        self.fdomain_samples = self.data_length * self.sample_rate / 2 + 1
        self.delta_f = 1.0 / self.data_length
        self.fmin = 30.0

        # set an analyitcal PSD for each detector
        psd = analytical.aLIGOZeroDetHighPower(self.fdomain_samples,
                                               self.delta_f, self.fmin)
        self.psds = {ifo : psd for ifo in self.ifos}

        # set parameter to use for generating waveform of test CBC signal
        cbc_test_parameters = (
            ("mass1", 30.0),
            ("mass2", 30.0),
            ("tc", 100.0),
            ("coa_phase", 1.1),
            ("spin1x", 0.0),
            ("spin1y", 0.0),
            ("spin1z", 0.0),
            ("spin2x", 0.0),
            ("spin2y", 0.0),
            ("spin2z", 0.0),
            ("ra", 0.1),
            ("dec", 0.1),
            ("polarization", 0.1),
            ("inclination", 0.1),
            ("distance", 300.0),
        )
        self.parameters, self.values = zip(*cbc_test_parameters)
        self.epoch = dict(cbc_test_parameters)["tc"]

        # get list of evaluators to test
        self.likelihood_evals = [self.get_cbc_fdomain_likelihood_evaluator()]

        # get a set of simulated command line options for sampler
        class Arguments(object):
            ntemps = 2
            nwalkers = 30
            niterations = 4
            update_interval = 2
            nprocesses = 2
        self.opts = Arguments()

    def get_cbc_fdomain_generator(self, approximant="IMRPhenomPv2"):
        """ Returns the waveform generator class for a CBC signal in the
        detector frame.
        """
        waveform_gen = generator.FDomainDetFrameGenerator(
                           generator.FDomainCBCGenerator, self.epoch,
                           variable_args=self.parameters,
                           detectors=self.ifos,
                           delta_f=self.delta_f, f_lower=self.fmin,
                           approximant=approximant)
        return waveform_gen

    def get_prior_evaluator(self, parameters, values):
        """ Returns the prior evaluator class initialized with a set of
        pre-defined distributions for each parameters.
        """
        prior_dists = []
        for param, val in zip(parameters, values):
            if param in ["mass1", "mass2"]:
                dist = distributions.Uniform(**{param : (6, 50)})
            elif param in ["inclination", "dec"]:
                dist = distributions.SinAngle(**{param : None})
            elif param in ["polarization", "ra", "coa_phase"]:
                dist = distributions.Uniform(**{param : (0, 2 * 3.1415)})
            elif param in ["distance"]:
                dist = distributions.UniformRadius(distance=(val - 100,
                                                             val + 300))
            elif param in ["spin1x", "spin1y", "spin1z",
                           "spin2x", "spin2y", "spin2z",]:
                dist = distributions.Uniform(**{param : (-0.1, 0.1)})
            elif param in ["tc"]:
                dist = distributions.Uniform(tc=(val - 0.2, val + 0.2))
            else:
                raise KeyError("Do not recognize parameter %s" % param)
            prior_dists.append(dist)
        return distributions.JointDistribution(parameters, *prior_dists)

    def get_likelihood_evaluator(self, waveform_gen, data, prior_eval):
        """ Returns the likelihood evaluator class.
        """
        likelihood_eval = likelihood.GaussianLikelihood(waveform_gen.variable_args,
                                             waveform_gen, data, self.fmin,
                                             psds=self.psds, prior=prior_eval,
                                             return_meta=False)
        return likelihood_eval

    def get_cbc_fdomain_likelihood_evaluator(self):
        """ Returns the likelihood evalauator class for a CBC signal.
        The data `FrequencySeries` used in the inner product is the signal
        generated from the waveform parameters in `self.setUp`.
        """

        # assert that waveform generator returns a dict of complex
        waveform_gen = self.get_cbc_fdomain_generator()
        signal = waveform_gen.generate(**{p : v
                                          for p, v in zip(self.parameters,
                                                          self.values)})
        assert(isinstance(signal, dict))
        assert(signal.values()[0].dtype == numpy.complex128)

        # assert that prior evaluator class returns a float
        prior_eval = self.get_prior_evaluator(self.parameters, self.values)
        p = prior_eval(**{p : v for p, v in zip(self.parameters,
                                                self.values)})
        assert(p.dtype == numpy.float64)

        # assert that likelihood evaluator returns a float
        # use generated waveform as data
        likelihood_eval = self.get_likelihood_evaluator(waveform_gen,
                                                        signal, prior_eval)
        assert(likelihood_eval(self.values).dtype == numpy.float64)

        return likelihood_eval

    def test_sampler(self):
        """ Runs each sampler for 4 iterations.
        """
        for likelihood_eval in self.likelihood_evals:
            for _, sampler_class in sampler.samplers.iteritems():
                s = sampler_class.from_cli(self.opts, likelihood_eval)
                s.set_p0()
                s.run(self.opts.niterations)

suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestSamplers))

if __name__ == "__main__":
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
