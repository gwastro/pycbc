# Copyright (C) 2016  Christopher M. Biwer
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
This modules provides classes and functions for evaluating the prior
for parameter estimation.
"""

import numpy

def flat_prior(param):
    """ Simple flat prior for a single parameter.
    """
    return 0.0

# different prior single parameters can have
prior_distributions = {
    "flat" : flat_prior,
}

class PriorEvaluator(object):
    """
    Callable class that calculates the prior.

    Parameters
    ----------
    variable_args : list
        A list of str that contain the names of the variable parameters.
    params_dist : {'flat'}
        A list of str that contain the names of the function in the prior
        module to use to calculate the prior for that variable.
    params_min : {None, list}
        The lowest acceptable value for parameter.
    params_max : {None, list}
        The largest acceptable value for parameter.
    """

    def __init__(self, variable_args, params_dist, params_min=None,
                 params_max=None):

        # store the names of the variable params
        self.variable_args = variable_args

        # store the distribution for the variable params
        self.params_dist = params_dist

        # store the minimum and maximum acceptable values for a variable param
        # if nothing is specified then use -inf to inf
        if params_min is None:
            self.params_min = numpy.ones(len(self.variable_args)) * -numpy.inf
        else:
            self.params_min = numpy.array(params_min)
        if params_max is None:
            self.params_max = numpy.ones(len(self.variable_args)) * numpy.inf
        else:
            self.params_max = numpy.array(params_max)

    def __call__(self, params):
        """ Evalualate prior for parameters.
        """

        # check if values are valid and if params are outside of min and max
        # acceptable values then return -inf
        params = numpy.array(params)
        if ((params < self.params_min) | (params > self.params_max)).all():
            return -numpy.inf

        # evaluate prior for each parameter
        return sum([prior_distributions[self.params_dist[i]](param) for i,param in enumerate(params)])

