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

    # Get the vairable arguments from the .ini file.
    variable_args, static_args = distributions.read_params_from_config(
        config_parser, prior_section='prior', vargs_section='variable_params',
        sargs_section='static_params')
    constraints = distributions.read_constraints_from_config(
        config_parser, static_args=static_args)

    if any(config_parser.get_subsections('waveform_transforms')):
        waveform_transforms = transforms.read_transforms_from_config(
            config_parser, 'waveform_transforms')
    else:
        waveform_transforms = None

    # Get prior distribution for each variable parameter.
    dists = distributions.read_distributions_from_config(config_parser)

    # Construct class that will draw the samples.
    randomsampler = distributions.JointDistribution(
                                variable_args, *dists,
                                **{"constraints": constraints})

    # Draw samples from prior distribution.
    samples = randomsampler.rvs(size=int(num))

    # Apply parameter transformation.
    if waveform_transforms is not None:
        samples = transforms.apply_transforms(samples, waveform_transforms)
    else:
        pass

    return samples
