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
import pycbc.distributions as distributions
from pycbc.coordinates import spherical_to_cartesian
from pycbc.cosmology import redshift_from_comoving_volume, distance_from_comoving_volume


def draw_samples_from_config(config_path, samples_num=1, seed=150914, **kwds):
    """ Generate sampling points from a standlone .ini file.
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
    # Add .ini file path.
    opts.config_files = [config_path]
    cp = WorkflowConfigParser.from_cli(opts)

    # Get the vairable and static arguments from the .ini file.
    variable_params, static_params = distributions.read_params_from_config(cp,
        prior_section=opts.dist_section,
        vargs_section=opts.variable_params_section,
        sargs_section=opts.static_params_section)
    constraints = distributions.read_constraints_from_config(cp)

    # Get prior distribution for each variable parameter.
    dists = distributions.read_distributions_from_config(cp, opts.dist_section)

    # Construct class that will draw the samples.
    randomsampler = distributions.JointDistribution(variable_params, *dists,
                                **{"constraints" : constraints})

    # Draw samples from prior distribution.
    samples = randomsampler.rvs(size=int(samples_num))

    # Pair the value with the parameter name.
    samples_list = []
    for sample_index in range(samples_num):
        sample = {}
        for param_index in variable_params:
            sample[param_index] = samples[sample_index][param_index]

        # Maps spherical coordinates to cartesian coordinates.
        sample['spin1x'],sample['spin1y'],sample['spin1z'] = spherical_to_cartesian(
                            sample['spin1_a'],sample['spin1_azimuthal'],sample['spin1_polar'])
        sample['spin2x'],sample['spin2y'],sample['spin2z'] = spherical_to_cartesian(
                            sample['spin2_a'],sample['spin2_azimuthal'],sample['spin2_polar'])

        # Returns the redshift & luminosity distance from the given comoving volume.
        sample['redshift'] = redshift_from_comoving_volume(sample['comoving_volume'])
        sample['distance'] = distance_from_comoving_volume(sample['comoving_volume'])
        
        # Calculate detector-frame masses.
        sample['mass1'] = sample['srcmass1'] * (1+sample['redshift'])
        sample['mass2'] = sample['srcmass2'] * (1+sample['redshift'])

        samples_list.append(sample)
    
    return samples_list
