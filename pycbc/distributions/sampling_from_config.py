# Copyright (C) 2021  Shichao Wu
#
#
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
# with this program; If not, see <http://www.gnu.org/licenses/>.


import numpy as np
from pycbc.types.config import InterpolatingConfigParser
from pycbc import transforms
import pycbc.distributions as distributions



def draw_samples_from_config(config_path, samples_num=1, seed=150914, **kwds):
    """ Generate sampling points from a standalone .ini file.
    Parameters
    ----------
    config_path : str
        The path of the .ini file.
    samples_num : int
        The number of samples.
    seed: int
        The random seed for sampling.
    Returns
    --------
    samples : list
        The parameter values [{param:value}] of sampling points in the parameter space.
    """

    np.random.seed(seed)

    # Initialise InterpolatingConfigParser class.
    cp = InterpolatingConfigParser()
    # Read the file
    fp = open(config_path,'r')
    cp.read_file(fp)
    fp.close()

    # Get the vairable and static arguments from the .ini file.
    variable_args, static_args = \
        distributions.read_params_from_config(cp, prior_section='prior',
                            vargs_section='variable_params',
                            sargs_section='static_params')
    constraints = distributions.read_constraints_from_config(cp)

    if any(cp.get_subsections('waveform_transforms')):
        waveform_transforms = transforms.read_transforms_from_config(cp,
            'waveform_transforms')
    else:
        waveform_transforms = None
        write_args = variable_args

    # Get prior distribution for each variable parameter.
    dists = distributions.read_distributions_from_config(cp)

    # Construct class that will draw the samples.
    randomsampler = distributions.JointDistribution(variable_args, *dists,
                                **{"constraints" : constraints})

    # Draw samples from prior distribution.
    samples = randomsampler.rvs(size=int(samples_num))

    # Apply parameter transformation.
    if waveform_transforms is not None:
        samples = transforms.apply_transforms(samples, waveform_transforms)
    else:
        pass

    # Pair the value with the parameter name.
    samples_list = []
    for sample_index in range(samples_num):
        sample = {}
        for param_index in samples.fieldnames:
            sample[param_index] = samples[sample_index][param_index]
        samples_list.append(sample)
    
    return samples_list
