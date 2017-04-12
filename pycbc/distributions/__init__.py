# Copyright (C)  2016 Collin Capano, Christopher M. Biwer
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

from pycbc.distributions.angular import *
from pycbc.distributions.change_of_variables import *
from pycbc.distributions.from_file import *
from pycbc.distributions.gaussian import *
from pycbc.distributions.power_law import *
from pycbc.distributions.sky_location import *
from pycbc.distributions.uniform import *

# a dict of all available distributions
distribs = {
    FromFile.name : FromFile,
    Gaussian.name : Gaussian,
    UniformRadius.name : UniformRadius,
    Uniform.name : Uniform,
    UniformAngle.name : UniformAngle,
    CosAngle.name : CosAngle,
    SinAngle.name : SinAngle,
    UniformChangeOfVariables.name : UniformChangeOfVariables,
    UniformSolidAngle.name : UniformSolidAngle,
    UniformSky.name : UniformSky,
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
