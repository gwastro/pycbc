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
        return self._prior_dist.pdf(**prior_dist_params)

    def logpdf(self, **kwargs):
        prior_dist_params = self._transform.convert(kwargs)
        return self._prior_dist.logpdf(**prior_dist_params)

    @classmethod
    def from_config(cls, cp, section, variable_args):

        # import mapping for Distributions
        # cannot be a top-level import because of circular dependencies
        from pycbc.distributions import distribs

        print variable_args

        sampling_name = cp.get_opt_tag(section, "sampling-name", variable_args)

        cp.add_section("-".join(["cov"])

        sampling_dist = bounded.bounded_from_config(
                                    distribs[sampling_name], cp, section,
                                    variable_args, additional_opts=None)

        # get Distributions
        sampling_dist = distribs[cp.get_opt_tag("prior", key, "")]
        prior_dist = distribs["prior_name"]
        return 0
