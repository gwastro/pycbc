# Copyright (C) 2016 Christopher M. Biwer
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
These are the unittests for the pycbc.inference subpackage
"""
import sys
import pycbc
import unittest
import numpy
from pycbc import distributions
from pycbc import inference
from pycbc import psd as pypsd
from utils import parse_args_cpu_only, simple_exit
from pycbc.waveform import generator

# tests only need to happen on the CPU
parse_args_cpu_only("Inference")

class TestInference(unittest.TestCase):
    def setUp(self, *args):
        numpy.random.seed(1024)

    def test_likelihood_evaluator_init(self):

        # data args
        seglen = 4
        sample_rate = 2048
        N = seglen * sample_rate/2 + 1
        fmin = 30.

        # setup waveform generator and signal parameters
        m1, m2, s1z, s2z, tsig = 38.6, 29.3, 0., 0., 3.1
        ra, dec, pol, dist = 1.37, -1.26, 2.76, 3 * 500.
        gen = generator.FDomainDetFrameGenerator(
                       generator.FDomainCBCGenerator, 0., variable_args=["tc"],
                       detectors=["H1", "L1"], delta_f=1./seglen,
                       f_lower=fmin, approximant="TaylorF2",
                       mass1=m1, mass2=m2, spin1z=s1z, spin2z=s2z, ra=ra,
                       dec=dec, polarization=pol, distance=dist)
        signal = gen.generate(tc=tsig)

        # get PSDs
        psd = pypsd.aLIGOZeroDetHighPower(N, 1. / seglen, 20.)
        psds = {"H1": psd, "L1": psd}

        # get a prior evaluator
        uniform_prior = distributions.Uniform(tc=(tsig - 0.2,tsig + 0.2))
        prior_eval = inference.prior.PriorEvaluator(["tc"], uniform_prior)

        # setup likelihood evaluator
        likelihood_eval = inference.GaussianLikelihood(
                                gen, signal, fmin,
                                psds=psds, prior=prior_eval, return_meta=False)

suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestInference))

if __name__ == "__main__":
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
