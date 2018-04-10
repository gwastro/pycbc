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
from pycbc.distributions.masses import UniformComponentMasses
from pycbc.distributions.spins import IndependentChiPChiEff
from pycbc.distributions.joint import JointDistribution

# a dict of all available distributions
distribs = {
    UniformComponentMasses.name : UniformComponentMasses,
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


def read_args_from_config(cp, section_group=None, prior_section='prior'):
    """Loads static and variable arguments from a configuration file.

    Parameters
    ----------
    cp : WorkflowConfigParser
        An open config parser to read from.
    section_group : {None, str}
        When reading the config file, only read from sections that begin with
        `{section_group}_`. For example, if `section_group='foo'`, the
        variable arguments will be retrieved from section
        `[foo_variable_args]`. If None, no prefix will be appended to section
        names.
    prior_section : str, optional
        Check that priors exist in the given section. Default is 'prior.'

    Returns
    -------
    variable_args : list
        The names of the parameters to vary in the PE run.
    static_args : dict
        Dictionary of names -> values giving the parameters to keep fixed.
    """
    if section_group is not None:
        section_prefix = '{}_'.format(section_group)
    else:
        section_prefix = ''

    # sanity check that each parameter in [variable_args] has a priors section
    variable_args = cp.options("{}variable_args".format(section_prefix))
    subsections = cp.get_subsections("{}{}".format(section_prefix,
                                                   prior_section))
    tags = set([p for tag in subsections for p in tag.split('+')])
    missing_prior = set(variable_args) - tags
    if any(missing_prior):
        raise KeyError("You are missing a priors section in the config file "
                       "for parameter(s): {}".format(', '.join(missing_prior)))

    # get parameters that do not change in sampler
    try:
        static_args = dict([(key,cp.get_opt_tags(
            "{}static_args".format(section_prefix), key, []))
            for key in cp.options("{}static_args".format(section_prefix))])
    except _ConfigParser.NoSectionError:
        static_args = {}
    # try converting values to float
    for key,val in static_args.iteritems():
        try:
            # the following will raise a ValueError if it cannot be cast to
            # float (as we would expect for string arguments)
            static_args[key] = float(val)
        except ValueError:
            # try converting to a list of strings; this function will just
            # return val if it does not begin (end) with [ (])
            static_args[key] = _convert_liststring_to_list(val)

    # get additional constraints to apply in prior
    cons = []
    section = "{}constraint".format(section_prefix)
    for subsection in cp.get_subsections(section):
        name = cp.get_opt_tag(section, "name", subsection)
        constraint_arg = cp.get_opt_tag(section, "constraint_arg", subsection)
        kwargs = {}
        for key in cp.options(section + "-" + subsection):
            if key in ["name", "constraint_arg"]:
                continue
            val = cp.get_opt_tag(section, key, subsection)
            if key == "required_parameters":
                kwargs["required_parameters"] = val.split(_VARARGS_DELIM)
                continue
            try:
                val = float(val)
            except ValueError:
                pass
            kwargs[key] = val
        cons.append(constraints.constraints[name](variable_args,
                                                  constraint_arg, **kwargs))

    return variable_args, static_args, cons
