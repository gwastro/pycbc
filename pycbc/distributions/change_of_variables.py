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

    Examples
    --------
    How to initialize a chirp mass and mass ratio to component mass instance:

    >>> from pycbc import transforms
    >>> from pycbc import distributions
    >>> mchirp_q_dist = distributions.Uniform(mchirp=(7, 40), q=(1, 10))
    >>> mass1_mass2_dist = distributions.Uniform(mass1=(8, 46), mass2=(2, 46))
    >>> cov_dist = distributions.ChangeOfVariables(mchirp_q_dist, mass1_mass2_dist, transforms.MchirpQToMass1Mass2)
    """
    name = "change_of_variables"

    def __init__(self, sampling_dist, prior_dist, transform):
        self._sampling_dist = sampling_dist
        self._prior_dist = prior_dist
        self._transform = transform
        self._bounds = self._sampling_dist._bounds
        self._params = self._sampling_dist._params

    def _pdf(self, **kwargs):
        """Returns the pdf at the given values. The keyword arguments must
        contain all of parameters in self's params. Unrecognized arguments are
        ignored. Any boundary conditions are applied to the values before the
        pdf is evaluated.
        """
        prior_dist_params = self._transform.convert(kwargs)
        return self._transform.jacobian(kwargs) \
                               * self._prior_dist.pdf(**prior_dist_params)

    def _logpdf(self, **kwargs):
        """Returns the log of the pdf at the given values. The keyword
        arguments must contain all of parameters in self's params.
        Unrecognized arguments are ignored. Any boundary conditions are
        applied to the values before the pdf is evaluated.
        """
        prior_dist_params = self._transform.convert(kwargs)
        return self._transform.jacobian(kwargs) \
                               + self._prior_dist.logpdf(**prior_dist_params)

    def rvs(self, size=1, param=None):
        """Gives a set of random values drawn from this distribution.

        Parameters
        ----------
        size : {1, int}
            The number of values to generate; default is 1.
        param : {None, string}
            If provided, will just return values for the given parameter.
            Otherwise, returns random values for each parameter.

        Returns
        -------
        structured array
            The random values in a numpy structured array. If a param was
            specified, the array will only have an element corresponding to the
            given parameter. Otherwise, the array will have an element for each
            parameter in self's params.
        """
        return self._sampling_dist.rvs(size=size, param=param)

    @classmethod
    def from_config(cls, cp, section, variable_args):
        """Returns a distribution based on a configuration file. The parameters
        for the distribution are retrieved from the section titled
        "[`section`-`variable_args`]" in the config file.

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
        ChangeOfVariable
            A ChangeOfVariable distribution instance.
        """

        # import mapping for Distributions
        # cannot be a top-level import because of circular dependencies
        from pycbc.distributions import distribs

        # a list of restricted options in the configuration file section
        # used internally by ChangeOfVariables
        # these options to do get transcribed to the new sections in the
        # configuration file that are created by from_config
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
            if cp.has_section(sec):
                cp.remove_section(sec)
            cp.add_section(sec)
            sec_opts = [(key, val) for key, val in opts.items()
                        if key not in restricted_opts]
            cp.add_options_to_section(sec, sec_opts)

        # create Distributions
        sampling_dist = distribs[sampling_name].from_config(
                           cp, cov_section, sampling_opts["parameters"])
        prior_dist = distribs[prior_name].from_config(
                           cp, cov_section, prior_opts["parameters"])

        # create Transformation
        for converter in transforms.all_converters:
            if converter.name == transform_name:
                transform = converter

        return cls(sampling_dist, prior_dist, transform)

__all__ = ["ChangeOfVariables"]
