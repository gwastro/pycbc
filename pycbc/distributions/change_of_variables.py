# Copyright (C) 2017  Christopher M. Biwer
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
This modules provides classes for evaluating distributions that include
a change of variables.
"""

import numpy
from pycbc import transforms
from pycbc.distributions import bounded

class ChangeOfVariables(bounded.BoundedDist):
    """ Class for sampling in one set of parameters and using a Jacobian
    transform to another set of parameters to calculate the
    probability density function.
    """

    name = "change_of_variables"

    def __init__(self, sampling_dist, prior_dist, transform):
        self._sampling_dist = sampling_dist
        self._prior_dist = prior_dist
        self._transform = transform

    @property
    def bounds(self):
        return self._sampling_dist.bounds

    @property
    def params(self):
        return self._sampling_dist.params

    def pdf(self, **kwargs):
        prior_dist_params = self._transform.convert(kwargs)
        return self._transform.jacobian(**kwargs) \
                               * self._prior_dist.pdf(**prior_dist_params)

    def logpdf(self, **kwargs):
        prior_dist_params = self._transform.convert(kwargs)
        return self._transform.jacobian(**kwargs) \
                               + self._prior_dist.logpdf(**prior_dist_params)

    @classmethod
    def from_config(cls, cp, section, variable_args):

        # import mapping for Distributions
        # cannot be a top-level import because of circular dependencies
        from pycbc.distributions import distribs

        # a list of restricted options in the configuration file section
        # used by internally by ChangeOfVariables
        restricted_opts = ["parameters"]

        # get name of Distributions and Transformation
        sampling_name = cp.get_opt_tag(section, "sampling-name", variable_args)
        prior_name = cp.get_opt_tag(section, "prior-name", variable_args)
        transform_name = cp.get_opt_tag(
                                    section, "transform-name", variable_args)

        # get options that need to be passed to Distributions
        sampling_opts = {}
        prior_opts = {}
        for opt in cp.options("-".join([section, variable_args])):
            if opt.startswith("sampling-"):
               stripped_opt_name = opt.split("sampling-")[1]
               sampling_opts[stripped_opt_name] = cp.get_opt_tag(
                                                   section, opt, variable_args)
            elif opt.startswith("prior-"):
               stripped_opt_name = opt.split("prior-")[1]
               prior_opts[stripped_opt_name] = cp.get_opt_tag(
                                                   section, opt, variable_args)

        # create a new section in configuration file for new Distributions
        cov_section = "changeofvariable"
        sampling_section = "-".join([cov_section, sampling_opts["parameters"]])
        prior_section = "-".join([cov_section, prior_opts["parameters"]])
        for sec, opts in zip([sampling_section, prior_section],
                             [sampling_opts, prior_opts]):
            cp.add_section(sec)
            sec_opts = [(key, val) for key, val in opts.items()
                        if key not in restricted_opts]
            cp.add_options_to_section(sec, sec_opts)

        # create Distributions
        sampling_dist = bounded.bounded_from_config(
                               distribs[sampling_name], cp,
                               cov_section, sampling_opts["parameters"])
        prior_dist = bounded.bounded_from_config(
                               distribs[prior_name], cp,
                               cov_section, prior_opts["parameters"])

        # create Transformation
        for converter in transforms.all_converters:
            if converter.name == transform_name:
                transform = converter

        return cls(sampling_dist, prior_dist, transform)
