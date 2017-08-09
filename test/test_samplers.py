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
from pycbc.inference import sampler
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
        self.f_lower = 20.
        self.delta_f = 1.0 / 64

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

    def get_cbc_fdomain_generator(self, approximant="SEOBNRv2_ROM_DoubleSpin"):
        """ Returns the waveform generator class for a CBC signal in the
        detector frame.
        """
        waveform_gen = generator.FDomainDetFrameGenerator(
                           generator.FDomainCBCGenerator, self.epoch,
                           variable_args=self.parameters,
                           detectors=self.ifos,
                           delta_f=self.delta_f, f_lower=self.f_lower,
                           approximant=approximant)
        return waveform_gen

    def test_cbc_fdomain_generator(self):
        """ Test that we can generate a dict of FrequencySeries waveforms
        for each detector.
        """
        waveform_gen = self.get_cbc_fdomain_generator()
        signal = waveform_gen.generate(**{p : v
                                          for p, v in zip(self.parameters,
                                                          self.values)})

    def test_cbc_fdomain_sampler(self):

        # get waveform generator in radiation frame
        waveform_gen = self.get_cbc_fdomain_generator()

        # get likelihood evaluator

        # get prior evaluator

        for sampler_class, sampler_name in sampler.samplers.iteritems():
            print sampler_class, sampler_name

suite = unittest.TestSuite()
#suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestGenerators))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestSamplers))

if __name__ == "__main__":
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
