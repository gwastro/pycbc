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
"""
This modules provides classes for evaluating distributions that use a Jacobian
tranformation to sample in a different set of variables.
"""

import numpy
from pycbc import conversions
from pycbc.distributions import bounded
from pycbc.distributions import uniform

class UniformChangeOfVariables(uniform.Uniform):
    """ Class for sampling uniform in one set of parameters (eg. mchirp, q)
    such that its uniform in another set of parameters.

    To add a new transformation, you need to add the input and output
    parameters in the transformation to UniformChangeOfVariables.mappings,
    add a function that computes the jacobian to
    UniformChangeOfVariables.jacobian, add a function that converts from input
    to output parameters to UniformChangeOfVariables.convert, and add a
    function that computes the limits of the output parameters in the
    transformation.
    """

    name = "uniform_change_of_variables"

    def __init__(self, **params):
        _from_params = {p : lims for p,lims in params.iteritems()
                        if not p.startswith("to")}
        super(UniformChangeOfVariables, self).__init__(**_from_params)

        # get set of to-parameters
        self.to_params = set(params["to_params"])

        # get str name for transformation
        self.mapping = self._get_mapping(self.params, self.to_params)

        # get bounds on to-parameters such that it fully encompasses all
        # possible values for from-parameters
        _to_params = self.to_params_bounds()

        # get Uniform distribution for from-parameters
        self.from_dist = uniform.Uniform(**_from_params)

        # get Uniform distribution for to-parameters
        self.to_dist = uniform.Uniform(**_to_params)

    def mass1_mass2_from_mchirp_q_jacobian(**kwargs):
        """ Returns the Jacobian of mass1 and mass2 from chirp mass
        and mass ratio."""
        mass2 = conversions.mass2_from_mchirp_q(kwargs["mchirp"], kwargs["q"])
        return kwargs["mchirp"] / mass2**2

    def mass1_mass2_from_mchirp_q_convert(**kwargs):
        """ Returns a dict with the mass1 and mass2 values from mchirp
        and q."""
        return {"mass1" : conversions.mass1_from_mchirp_q(
                                                kwargs["mchirp"], kwargs["q"]),
                "mass2" : conversions.mass2_from_mchirp_q(
                                                kwargs["mchirp"], kwargs["q"])}

    def mass1_mass2_from_mchirp_q_limit(self, **kwargs):
        """ Returns a dict with the (min, max) bounds for mass1 and mass2
        from mchirp and q."""
        min_mchirp = self.bounds["mchirp"][0]
        max_mchirp = self.bounds["mchirp"][1]
        min_q = self.bounds["q"][0]
        max_q = self.bounds["q"][1]
        min_mass1 = conversions.mass1_from_mchirp_q(min_mchirp, min_q)
        max_mass1 = conversions.mass1_from_mchirp_q(max_mchirp, max_q)
        min_mass2 = conversions.mass2_from_mchirp_q(min_mchirp, max_q)
        max_mass2 = conversions.mass2_from_mchirp_q(max_mchirp, min_q)
        return {"mass1" : (min_mass1, max_mass1),
                "mass2" : (min_mass2, max_mass2)}

    # a dict with a string name as key and (in-parameters, out-parameters)
    # as value
    mappings = {
        "mass1_mass2_from_mchirp_q" : (frozenset(["mchirp", "q"]),
                                       frozenset(["mass1", "mass2"])),
    }

    # a dict with (in-parameters, out-parameters) as key and a function
    # that returns the jacobian as value
    jacobians = {
        "mass1_mass2_from_mchirp_q" : mass1_mass2_from_mchirp_q_jacobian,
    }

    # a dict with (in-parameters, out-parameters) as key and a function
    # that returns the converted in-parameters in the new coordinates
    converts = {
        "mass1_mass2_from_mchirp_q" : mass1_mass2_from_mchirp_q_convert,
    }

    # a dict with (in-parameters, out-parameters) as key and a function
    # that returns the bounds of each out-parameter in a dict
    limits = {
        "mass1_mass2_from_mchirp_q" : mass1_mass2_from_mchirp_q_limit,
    }

    @classmethod
    def _get_mapping(cls, from_params, to_params):
        """ Returns value of dictionary where key is
        (in-parameters, out-parameters)."""
        for key, val in cls.mappings.iteritems():
            from_params_in_mapping = set(val[0])
            to_params_in_mapping = set(val[1])
            if len(from_params_in_mapping.difference(from_params)) == 0 \
                     and len(to_params_in_mapping.difference(to_params)) == 0:
                return key
        raise ValueError("Did not find %s to %s "
                         "transformation" % (str(from_params), str(to_params)))

    def jacobian(self, **kwargs):
        """ Returns Jacobian."""
        return self.jacobians[self.mapping](**kwargs)

    def convert(self, **kwargs):
        """ Returns converted parameters."""
        return self.converts[self.mapping](**kwargs)

    def to_params_bounds(self, **kwargs):
        """ Returns bounds of out-parameters."""
        return self.limits[self.mapping](self, **kwargs)

    def _pdf(self, **kwargs):
        """Returns the pdf at the given values. The keyword arguments must
        contain all of parameters in self's params. Unrecognized arguments are
        ignored.
        """
        jacobian = 1.0 / self.jacobian(**kwargs)
        to_params = self.convert(**kwargs)
        return jacobian * self.to_dist.pdf(**to_params) \
                                                 / self.from_dist.pdf(**kwargs)

    def _logpdf(self, **kwargs):
        """Returns the log of the pdf at the given values. The keyword
        arguments must contain all of parameters in self's params. Unrecognized
        arguments are ignored.
        """
        if kwargs in self:
            return numpy.log(self._pdf(**kwargs))
        else:
            return -numpy.inf

    @classmethod
    def from_config(cls, cp, section, variable_args):
        """ Returns a distribution based on a configuration file. The
        parameters for the distribution are retrieved from the section titled
        "[`section`-`variable_args`]" in the config file.

        Boundary arguments should be provided in the same way as described in
        `get_param_bounds_from_config`. As an example, this uses the jacobian
        for mchirp and q to mass1 and mass2:

        .. code-block:: ini

            [prior-mchirp+q]
            name = uniform_change_of_variables
            min-mchirp = 10.
            max-mchirp = 20.
            min-q = 1.
            max-q = 5.
            to-params = mass1+mass2

        Parameters
        ----------
        cp : pycbc.workflow.WorkflowConfigParser
            A parsed configuration file that contains the distribution
            options.
        section : str
            Name of the section in the configuration file.
        variable_args : str
            The names of the parameters for this distribution, separated by
            `prior.VARARGS_DELIM`. These must appear in the "tag" part
            of the section header.

        Returns
        -------
        UniformChangeOfVariables
            A distribution instance from the pycbc.distribution subpackage.
        """
        to_params = cp.get_opt_tags(section + "-" + variable_args,
                                    "to-params", []).split("+")
        return bounded.bounded_from_config(
                                     cls, cp, section, variable_args,
                                     bounds_required=True, to_params=to_params)

