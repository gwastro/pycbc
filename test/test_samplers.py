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
from pycbc.inference import prior
from pycbc.inference import sampler
from pycbc.psd import analytical
from pycbc.waveform import generator
from utils import parse_args_cpu_only
from utils import  simple_exit

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
            ("spin1z", 0.0),
            ("spin2z", 0.0),
            ("ra", 0.1),
            ("dec", 0.1),
            ("polarization", 0.1),
            ("inclination", 0.1),
            ("distance", 300.0),
        )
        self.parameters, self.values = zip(*cbc_test_parameters)
        self.epoch = dict(cbc_test_parameters)["tc"]

    def get_cbc_fdomain_generator(self, approximant="SEOBNRv2_ROM_DoubleSpin"):
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
            elif param in ["spin1z", "spin2z"]:
                dist = distributions.Uniform(**{param : (-1, 1)})
            elif param in ["tc"]:
                dist = distributions.Uniform(tc=(val - 0.2, val + 0.2))
            else:
                raise KeyError("Do not recognize parameter %s" % param)
            prior_dists.append(dist)
        return prior.PriorEvaluator(parameters, *prior_dists)

    def test_cbc_fdomain_generator(self):
        """ Test that we can generate a dict of FrequencySeries waveforms
        for each detector.
        """
        waveform_gen = self.get_cbc_fdomain_generator()
        signal = waveform_gen.generate(**{p : v
                                          for p, v in zip(self.parameters,
                                                          self.values)})
        assert(type(signal) == dict)
        assert(signal.values()[0].dtype == numpy.complex128)

    def test_cbc_prior_evaluator(self):
        """ Test that prior evaluator class returns a float.
        """
        prior_eval = self.get_prior_evaluator(self.parameters, self.values)
        p = prior_eval(**{p : v for p, v in zip(self.parameters,
                                                self.values)})
        assert(p.dtype == numpy.float64)

    def test_cbc_fdomain_sampler(self):

        # get waveform generator in radiation frame
        waveform_gen = self.get_cbc_fdomain_generator()

        # get prior evaluator
        prior_eval = self.get_prior_evaluator(self.parameters, self.values)

        # get likelihood evaluator

        for sampler_class, sampler_name in sampler.samplers.iteritems():
            print sampler_class, sampler_name

suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestSamplers))

if __name__ == "__main__":
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
