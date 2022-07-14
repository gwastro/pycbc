# Copyright (C) 2021  Shichao Wu
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
This module provides functions for drawing samples from a standalone .ini file
in a Python script, rather than in the command line.
"""

import numpy as np
from pycbc.types.config import InterpolatingConfigParser
from pycbc import transforms
from pycbc import distributions


def prior_from_config(cp, prior_section='prior'):
    """Loads a prior distribution from the given config file.

    Parameters
    ----------
    cp : pycbc.workflow.WorkflowConfigParser
        The config file to read.
    sections : list of str, optional
        The sections to retrieve the prior from. If ``None`` (the default),
        will look in sections starting with 'prior'.

    Returns
    -------
    distributions.JointDistribution
        The prior distribution.
    """

    # Read variable and static parameters from the config file
    variable_params, static_params = distributions.read_params_from_config(
        cp, prior_section=prior_section, vargs_section='variable_params',
        sargs_section='static_params')

    # Read waveform_transforms to apply to priors from the config file
    if any(cp.get_subsections('waveform_transforms')):
        waveform_transforms = transforms.read_transforms_from_config(
                cp, 'waveform_transforms')
    else:
        waveform_transforms = None

    # Read constraints to apply to priors from the config file
    constraints = distributions.read_constraints_from_config(
        cp, transforms=waveform_transforms, static_args=static_params)

    # Get PyCBC distribution instances for each variable parameter in the
    # config file
    dists = distributions.read_distributions_from_config(cp, prior_section)

    # construct class that will return draws from the prior
    return distributions.JointDistribution(variable_params, *dists,
                                           **{"constraints": constraints})


def draw_samples_from_config(path, num=1, seed=150914):
    r""" Generate sampling points from a standalone .ini file.

    Parameters
    ----------
    path : str
        The path to the .ini file.
    num : int
        The number of samples.
    seed: int
        The random seed for sampling.

    Returns
    --------
    samples : pycbc.io.record.FieldArray
        The parameter values and names of sample(s).

    Examples
    --------
    Draw a sample from the distribution defined in the .ini file:

    >>> import numpy as np
    >>> from pycbc.distributions.utils import draw_samples_from_config

    >>> # A path to the .ini file.
    >>> CONFIG_PATH = "./pycbc_bbh_prior.ini"
    >>> random_seed = np.random.randint(low=0, high=2**32-1)
    >>> sample = draw_samples_from_config(
    >>>          path=CONFIG_PATH, num=1, seed=random_seed)

    >>> # Print all parameters.
    >>> print(sample.fieldnames)
    >>> print(sample)
    >>> # Print a certain parameter, for example 'mass1'.
    >>> print(sample[0]['mass1'])
    """

    np.random.seed(seed)

    # Initialise InterpolatingConfigParser class.
    config_parser = InterpolatingConfigParser()
    # Read the file
    file = open(path, 'r')
    config_parser.read_file(file)
    file.close()

    # Construct class that will draw the samples.
    prior_dists = prior_from_config(cp=config_parser)
    # Draw samples from prior distribution.
    samples = prior_dists.rvs(size=int(num))

    # Apply parameter transformation.
    if any(config_parser.get_subsections('waveform_transforms')):
        waveform_transforms = transforms.read_transforms_from_config(
                config_parser, 'waveform_transforms')
        samples = transforms.apply_transforms(samples, waveform_transforms)

    return samples
