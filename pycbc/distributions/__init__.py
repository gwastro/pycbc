# Copyright (C)  2016 Collin Capano, Christopher M. Biwer, Alex Nitz
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
This modules provides classes and functions for drawing and calculating the
probability density function of distributions.
"""
# imports needed for functions below
from pycbc.workflow import ConfigParser as _ConfigParser
from pycbc.distributions import constraints
from pycbc import VARARGS_DELIM as _VARARGS_DELIM

# Promote some classes/functions to the distributions name space
from pycbc.distributions.angular import UniformAngle, SinAngle, CosAngle, \
                                        UniformSolidAngle
from pycbc.distributions.arbitrary import Arbitrary, FromFile
from pycbc.distributions.gaussian import Gaussian
from pycbc.distributions.power_law import UniformPowerLaw, UniformRadius
from pycbc.distributions.sky_location import UniformSky
from pycbc.distributions.uniform import Uniform
from pycbc.distributions.uniform_log import UniformLog10
from pycbc.distributions.spins import IndependentChiPChiEff
from pycbc.distributions.qnm import UniformF0Tau
from pycbc.distributions.joint import JointDistribution

# a dict of all available distributions
distribs = {
    IndependentChiPChiEff.name : IndependentChiPChiEff,
    Arbitrary.name : Arbitrary,
    FromFile.name : FromFile,
    Gaussian.name : Gaussian,
    UniformPowerLaw.name : UniformPowerLaw,
    UniformRadius.name : UniformRadius,
    Uniform.name : Uniform,
    UniformAngle.name : UniformAngle,
    CosAngle.name : CosAngle,
    SinAngle.name : SinAngle,
    UniformSolidAngle.name : UniformSolidAngle,
    UniformSky.name : UniformSky,
    UniformLog10.name : UniformLog10,
    UniformF0Tau.name : UniformF0Tau,
}

def read_distributions_from_config(cp, section="prior"):
    """Returns a list of PyCBC distribution instances for a section in the
    given configuration file.

    Parameters
    ----------
    cp : WorflowConfigParser
        An open config file to read.
    section : {"prior", string}
        Prefix on section names from which to retrieve the distributions.

    Returns
    -------
    list
        A list of the parsed distributions.
    """
    dists = []
    variable_args = []
    for subsection in cp.get_subsections(section):
        name = cp.get_opt_tag(section, "name", subsection)
        dist = distribs[name].from_config(cp, section, subsection)
        if set(dist.params).isdisjoint(variable_args):
            dists.append(dist)
            variable_args += dist.params
        else:
            raise ValueError("Same parameter in more than one distribution.")
    return dists


def _convert_liststring_to_list(lstring):
    """Checks if an argument of the configuration file is a string of a list
    and returns the corresponding list (of strings).

    The argument is considered to be a list if it starts with '[' and ends
    with ']'. List elements should be comma separated. For example, passing
    `'[foo bar, cat]'` will result in `['foo bar', 'cat']` being returned. If
    the argument does not start and end with '[' and ']', the argument will
    just be returned as is.
    """
    if lstring[0]=='[' and lstring[-1]==']':
        lstring = [str(lstring[1:-1].split(',')[n].strip().strip("'"))
                      for n in range(len(lstring[1:-1].split(',')))]
    return lstring


def read_params_from_config(cp, prior_section='prior',
                            vargs_section='variable_args',
                            sargs_section='static_args'):
    """Loads static and variable parameters from a configuration file.

    Parameters
    ----------
    cp : WorkflowConfigParser
        An open config parser to read from.
    prior_section : str, optional
        Check that priors exist in the given section. Default is 'prior.'
    vargs_section : str, optional
        The section to get the parameters that will be varied/need priors
        defined for them. Default is 'variable_args'.
    sargs_section : str, optional
        The section to get the parameters that will remain fixed. Default is
        'static_args'.

    Returns
    -------
    variable_args : list
        The names of the parameters to vary in the PE run.
    static_args : dict
        Dictionary of names -> values giving the parameters to keep fixed.
    """
    # sanity check that each parameter in [variable_args] has a priors section
    variable_args = cp.options(vargs_section)
    subsections = cp.get_subsections(prior_section)
    tags = set([p for tag in subsections for p in tag.split('+')])
    missing_prior = set(variable_args) - tags
    if any(missing_prior):
        raise KeyError("You are missing a priors section in the config file "
                       "for parameter(s): {}".format(', '.join(missing_prior)))
    # sanity check that each parameter with a priors section is in
    # [variable_args]
    missing_variable = tags - set(variable_args)
    if any(missing_variable):
        raise KeyError("Prior section found for parameter(s) {} but not "
                       "listed as variable parameter(s)."
                       .format(', '.join(missing_variable)))
    # get static args
    try:
        static_args = dict([(key, cp.get_opt_tags(sargs_section, key, []))
                           for key in cp.options(sargs_section)])
    except _ConfigParser.NoSectionError:
        static_args = {}
    # try converting values to float
    for key in static_args:
        val = static_args[key]
        try:
            # the following will raise a ValueError if it cannot be cast to
            # float (as we would expect for string arguments)
            static_args[key] = float(val)
        except ValueError:
            # try converting to a list of strings; this function will just
            # return val if it does not begin (end) with [ (])
            static_args[key] = _convert_liststring_to_list(val)
    return variable_args, static_args


def read_constraints_from_config(cp, transforms=None,
                                 constraint_section='constraint'):
    """Loads parameter constraints from a configuration file.

    Parameters
    ----------
    cp : WorkflowConfigParser
        An open config parser to read from.
    transforms : list, optional
        List of transforms to apply to parameters before applying constraints.
    constraint_section : str, optional
        The section to get the constraints from. Default is 'constraint'.

    Returns
    -------
    list
        List of ``Constraint`` objects. Empty if no constraints were provided.
    """
    cons = []
    for subsection in cp.get_subsections(constraint_section):
        name = cp.get_opt_tag(constraint_section, "name", subsection)
        constraint_arg = cp.get_opt_tag(
            constraint_section, "constraint_arg", subsection)
        # get any other keyword arguments
        kwargs = {}
        section = constraint_section + "-" + subsection
        extra_opts = [key for key in cp.options(section)
                      if key not in ["name", "constraint_arg"]]
        for key in extra_opts:
            val = cp.get(section, key)
            if key == "required_parameters":
                val = val.split(_VARARGS_DELIM)
            else:
                try:
                    val = float(val)
                except ValueError:
                    pass
            kwargs[key] = val
        cons.append(constraints.constraints[name](constraint_arg,
                                                  transforms=transforms,
                                                  **kwargs))

    return cons
