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

"""
These are the unittests for distributions in the pycbc.distribtions subpackage.
"""

import itertools
import numpy
import os
import unittest
from pycbc import distributions
from pycbc.inference import entropy
from utils import parse_args_cpu_only
from utils import simple_exit
from pycbc.workflow import WorkflowConfigParser

# distributions to exclude from one-dimensional distribution unit tests
# some of these distributons have their own specific unit test
EXCLUDE_DIST_NAMES = ["fromfile", "arbitrary",
                      "external", "external_func_fromfile",
                      "uniform_solidangle", "uniform_sky",
                      "independent_chip_chieff",
                      "uniform_f0_tau", "fixed_samples"]

# tests only need to happen on the CPU
parse_args_cpu_only("Distributions")

def cartesian(arrays):
    """ Returns a cartesian product from a list of iterables.
    """
    return numpy.array([numpy.array(element)
                        for element in itertools.product(*arrays)])

class TestDistributions(unittest.TestCase):

    def setUp(self):

        # set random seed
        numpy.random.seed(1024)

        # path to example configuration file for testing
        config_path = "/".join([os.path.dirname(os.path.realpath(__file__)),
                                "../examples/distributions/example.ini"])

        # get a set of simulated command line options for
        # configuration file reading
        class Arguments(object):
            config_overrides = []
            config_delete = []
            config_files = [config_path]
        self.opts = Arguments()

        # read configuration files
        self.cp = WorkflowConfigParser.from_cli(self.opts)
        self.variable_args, self.static_args = \
            distributions.read_params_from_config(self.cp)
        self.constraints = distributions.read_constraints_from_config(
            self.cp, static_args=self.static_args)

        # read distributions
        self.dists = distributions.read_distributions_from_config(self.cp)

        # check that all distriubtions will be tested
        for dname in distributions.distribs:
            dclass = distributions.distribs[dname]
            if (not numpy.any([isinstance(dist, dclass)
                               for dist in self.dists])
                    and dname not in EXCLUDE_DIST_NAMES):
                raise ValueError("There is no test for {}".format(dname))

    def test_pdf_rvs(self):
        """ Check the Kullback-Leibler divergence between draws of random
        samples form the distribution and the probability density function
        of the distribution. This implementation only works for
        one dimensional distriubtions.
        """

        # set threshold for KL divergence
        threshold = 0.1

        # number of samples in random draw for test
        n_samples = int(1e6)

        # step size to take in PDF evaluation
        step = 0.1

        # loop over distributions
        for dist in self.dists:
            if dist.name in EXCLUDE_DIST_NAMES:
                continue
            for param in dist.params:

                # get min and max
                hist_min = dist.bounds[param][0]
                hist_max = dist.bounds[param][1]

                # generate some random draws
                samples = dist.rvs(n_samples)[param]

                # get the PDF
                x = numpy.arange(hist_min, hist_max, step)
                pdf = numpy.array([dist.pdf(**{param : xx}) for xx in x])

                # compute the KL divergence and check if below threshold
                kl_val = entropy.kl(samples, pdf, bins=pdf.size, pdf2=True,
                                    hist_min=hist_min, hist_max=hist_max)
                if not (kl_val < threshold):
                    raise ValueError(
                              "Class {} KL divergence is {} which is "
                              "greater than the threshold "
                              "of {}".format(dist.name, kl_val, threshold))

    def test_pdf_logpdf(self):
        """ Checks that the probability density function (PDF) is within some
        tolerance of the natural logarithm of the PDF. This implementation is
        for one dimensional distributions.
        """

        # assign tolerance for element-wise ratio of logarithm of PDF
        tolerance = 0.01

        # step size to take in PDF evaluation
        step = 0.1

        # loop over distributions
        for dist in self.dists:
            if dist.name in EXCLUDE_DIST_NAMES:
                continue
            for param in dist.params:

                # get min and max
                hist_min = dist.bounds[param][0]
                hist_max = dist.bounds[param][1]

                # get the PDF and logarithm of the PDF from the distriubtion
                x = numpy.arange(hist_min, hist_max, step)
                pdf = numpy.array([dist.pdf(**{param : xx}) for xx in x])
                logpdf = numpy.array([dist.logpdf(**{param : xx}) for xx in x])

                # find the logarithm of the PDF
                pdf_log = numpy.log(pdf)

                # see if each element in ratio of these two logarithm of PDF
                # values are within the specified tolerance
                if not numpy.all(abs(1.0 - logpdf / pdf_log) < tolerance):
                    raise ValueError("The PDF and logarithm of the PDF "
                                     "functions for distribution {} "
                                     "do not agree".format(dist.name))

    def test_solid_angle(self):
        """ The uniform solid angle and uniform sky position distributions
        are two independent one-dimensional distributions. This tests checks
        that the two indepdent one-dimensional distributions agree by comparing
        the `rvs`, `pdf`, and `logpdf` functions' output.
        """

        # set tolerance for comparing PDF and logPDF functions
        tolerance = 0.01

        # set threshold for KL divergence test
        threshold = 0.1

        # number of random draws for KL divergence test
        n_samples = int(1e6)

        # create generic angular distributions for test
        sin_dist = distributions.SinAngle(theta=(0, numpy.pi))
        cos_dist = distributions.CosAngle(theta=(-numpy.pi/2.0, numpy.pi/2.0))
        ang_dist = distributions.UniformAngle(theta=(0, numpy.pi*2.0))

        # step size for PDF calculation
        step = 0.1

        # valid range of parameters
        polar_sin = numpy.arange(0, numpy.pi, step)
        polar_cos = numpy.arange(-numpy.pi, numpy.pi, step)
        azimuthal = numpy.arange(0, 2 * numpy.pi, step)

        # get Cartesian product to explore the two-dimensional space
        cart_sin = cartesian([polar_sin, azimuthal])

        # loop over distributions
        for dist in self.dists:
            if dist.name == distributions.UniformSolidAngle.name:
                polar_vals = polar_sin
                polar_dist = sin_dist
            elif dist.name == distributions.UniformSky.name:
                polar_vals = polar_cos
                polar_dist = cos_dist
            else:
                continue

            # check PDF equilvalent
            pdf_1 = numpy.array([dist.pdf(**{dist.polar_angle : p,
                                             dist.azimuthal_angle : a})
                                 for p, a in cart_sin])
            pdf_2 = numpy.array([polar_dist.pdf(**{"theta" : p})
                                 * ang_dist.pdf(**{"theta" : a})
                                 for p, a in cart_sin])

            # Catch and silence warnings here
            with numpy.errstate(invalid="ignore", divide='ignore'):
                if not (numpy.all(numpy.nan_to_num(abs(1.0 - pdf_1 / pdf_2))
                                  < tolerance)):
                    raise ValueError("The {} distribution PDF does not match "
                                     "its component "
                                     "distributions.".format(dist.name))

                # check logarithm of PDF equivalent
                pdf_1 = numpy.array([dist.logpdf(**{dist.polar_angle : p,
                                                    dist.azimuthal_angle : a})
                                     for p, a in cart_sin])
                pdf_2 = numpy.array([polar_dist.logpdf(**{"theta" : p})
                                     + ang_dist.logpdf(**{"theta" : a})
                                     for p, a in cart_sin])
                if not (numpy.all(numpy.nan_to_num(abs(1.0 - pdf_1 / pdf_2))
                           < tolerance)):
                    raise ValueError("The {} distribution PDF does not match "
                                     "its component "
                                     "distriubtions.".format(dist.name))

            # check random draws from polar angle equilvalent
            ang_1 = dist.rvs(n_samples)[dist.polar_angle]
            ang_2 = polar_dist.rvs(n_samples)["theta"]
            kl_val = entropy.kl(ang_1, ang_2, bins=polar_vals.size,
                                hist_min=polar_vals.min(),
                                hist_max=polar_vals.max())
            if not (kl_val < threshold):
                raise ValueError(
                          "Class {} KL divergence is {} which is "
                          "greater than the threshold for polar angle"
                          "of {}".format(dist.name, kl_val, threshold))

            # check random draws from azimuthal angle equilvalent
            ang_1 = dist.rvs(n_samples)[dist.azimuthal_angle]
            ang_2 = ang_dist.rvs(n_samples)["theta"]
            kl_val = entropy.kl(ang_1, ang_2, bins=azimuthal.size,
                                hist_min=azimuthal.min(),
                                hist_max=azimuthal.max())
            if not (kl_val < threshold):
                raise ValueError(
                          "Class {} KL divergence is {} which is "
                          "greater than the threshold for azimuthal angle"
                          "of {}".format(dist.name, kl_val, threshold))

suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestDistributions))

if __name__ == "__main__":
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
