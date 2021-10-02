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


import argparse
import numpy as np
from pycbc.workflow import WorkflowConfigParser, configuration
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
    samples : numpy.record
        The parameter values of sampling points in the parameter space.
    """

    np.random.seed(seed)
    parser = argparse.ArgumentParser(description=__doc__,
                formatter_class=argparse.RawDescriptionHelpFormatter)
    configuration.add_workflow_command_line_group(parser)
    parser.add_argument('--dist-section', default='prior',
                        help='What section in the config file to load '
                            'distributions from. Default is prior.')
    parser.add_argument('--variable-params-section', default='variable_params',
                        help='What section in the config file to load the '
                            'parameters to vary. Default is variable_params.')
    parser.add_argument('--static-params-section', default="static_params",
                        help='What section to load the static params from. Default '
                            'is static_params.')
    opts = parser.parse_args()
    # add .ini file path
    opts.config_files = [config_path]
    cp = WorkflowConfigParser.from_cli(opts)

    # get the vairable and static arguments from the .ini file
    variable_params, static_params = distributions.read_params_from_config(cp,
        prior_section=opts.dist_section,
        vargs_section=opts.variable_params_section,
        sargs_section=opts.static_params_section)
    constraints = distributions.read_constraints_from_config(cp)

    if any(cp.get_subsections('waveform_transforms')):
        waveform_transforms = transforms.read_transforms_from_config(cp,
            'waveform_transforms')
    else:
        waveform_transforms = None
        write_args = variable_params

    # get prior distribution for each variable parameter
    dists = distributions.read_distributions_from_config(cp, opts.dist_section)

    # construct class that will draw the samples
    randomsampler = distributions.JointDistribution(variable_params, *dists,
                                **{"constraints" : constraints})

    # draw samples from prior distribution
    samples = randomsampler.rvs(size=int(samples_num))
    
    return samples[0]
