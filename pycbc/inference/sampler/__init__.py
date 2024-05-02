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
This module provides a list of implemented samplers for parameter estimation.
"""

import logging

# pylint: disable=unused-import
from .base import (initial_dist_from_config, create_new_output_file)
from .multinest import MultinestSampler
from .ultranest import UltranestSampler
from .dummy import DummySampler
from .refine import RefineSampler
from .snowline import SnowlineSampler

# list of available samplers
samplers = {cls.name: cls for cls in (
    MultinestSampler,
    UltranestSampler,
    DummySampler,
    RefineSampler,
    SnowlineSampler,
)}

try:
    from .emcee import EmceeEnsembleSampler
    from .emcee_pt import EmceePTSampler
    samplers[EmceeEnsembleSampler.name] = EmceeEnsembleSampler
    samplers[EmceePTSampler.name] = EmceePTSampler
except ImportError:
    pass

try:
    from .epsie import EpsieSampler
    samplers[EpsieSampler.name] = EpsieSampler
except ImportError:
    pass

try:
    from .ptemcee import PTEmceeSampler
    samplers[PTEmceeSampler.name] = PTEmceeSampler
except ImportError:
    pass

try:
    from .cpnest import CPNestSampler
    samplers[CPNestSampler.name] = CPNestSampler
except ImportError:
    pass

try:
    from .dynesty import DynestySampler
    samplers[DynestySampler.name] = DynestySampler
except ImportError:
    pass

try:
    from .nessai import NessaiSampler
    samplers[NessaiSampler.name] = NessaiSampler
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
    **kwargs :
        All other keyword arguments are passed directly to the sampler's
        ``from_config`` file.

    Returns
    -------
    sampler :
        The initialized sampler.
    """
    if len(model.variable_params) == 0:
        logging.info('No variable params, so assuming Dummy Sampler')
        return DummySampler.from_config(cp, model, **kwargs)

    name = cp.get('sampler', 'name')
    try:
        return samplers[name].from_config(cp, model, **kwargs)
    except KeyError:
        raise ImportError(
            f"No available sampler named {name}. Please check "
            "if the name is correct or the required package "
            "for this sampler is installed correctly. "
            f"Available samplers: {', '.join(list(samplers.keys()))}"
        )
