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
This modules provides a list of implemented samplers for parameter estimation.
"""

from __future__ import absolute_import
# pylint: disable=unused-import
from .base import (initial_dist_from_config, create_new_output_file)
from .emcee import EmceeEnsembleSampler
from .emcee_pt import EmceePTSampler
from .multinest import MultinestSampler

# list of available samplers
samplers = {cls.name: cls for cls in (
    EmceeEnsembleSampler,
    EmceePTSampler,
    MultinestSampler
)}

try:
    from .cpnest import CPNestSampler
    samplers[CPNestSampler.name] = CPNestSampler
except ImportError:
    pass

def load_from_config(cp, model, **kwargs):
    """Loads a sampler from the given config file.

    This looks for a name in the section ``[sampler]`` to determine which
    sampler class to load. That sampler's ``from_config`` is then called.

    Parameters
    ----------
    cp : WorkflowConfigParser
        Config parser to read from.
    model : pycbc.inference.model
        Which model to pass to the sampler.
    \**kwargs :
        All other keyword arguments are passed directly to the sampler's
        ``from_config`` file.

    Returns
    -------
    sampler :
        The initialized sampler.
    """
    name = cp.get('sampler', 'name')
    return samplers[name].from_config(cp, model, **kwargs)
