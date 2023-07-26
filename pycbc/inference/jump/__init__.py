# Copyright (C) 2019  Collin Capano
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
"""Provides custom jump proposals for samplers."""

from .normal import (EpsieNormal, EpsieAdaptiveNormal, EpsieATAdaptiveNormal)
from .bounded_normal import (EpsieBoundedNormal, EpsieAdaptiveBoundedNormal,
                             EpsieATAdaptiveBoundedNormal)
from .angular import (EpsieAngular, EpsieAdaptiveAngular,
                      EpsieATAdaptiveAngular)
from .discrete import (EpsieNormalDiscrete, EpsieBoundedDiscrete,
                       EpsieAdaptiveNormalDiscrete,
                       EpsieAdaptiveBoundedDiscrete)


epsie_proposals = {
    EpsieNormal.name: EpsieNormal,
    EpsieAdaptiveNormal.name: EpsieAdaptiveNormal,
    EpsieATAdaptiveNormal.name: EpsieATAdaptiveNormal,
    EpsieBoundedNormal.name: EpsieBoundedNormal,
    EpsieAdaptiveBoundedNormal.name: EpsieAdaptiveBoundedNormal,
    EpsieATAdaptiveBoundedNormal.name: EpsieATAdaptiveBoundedNormal,
    EpsieAngular.name: EpsieAngular,
    EpsieAdaptiveAngular.name: EpsieAdaptiveAngular,
    EpsieATAdaptiveAngular.name: EpsieATAdaptiveAngular,
    EpsieNormalDiscrete.name: EpsieNormalDiscrete,
    EpsieAdaptiveNormalDiscrete.name: EpsieAdaptiveNormalDiscrete,
    EpsieBoundedDiscrete.name: EpsieBoundedDiscrete,
    EpsieAdaptiveBoundedDiscrete.name: EpsieAdaptiveBoundedDiscrete,
}


def epsie_proposals_from_config(cp, section='jump_proposal'):
    """Loads epsie jump proposals from the given config file.

    This loads jump proposals from sub-sections starting with ``section``
    (default is 'jump_proposal'). The tag part of the sub-sections' headers
    should list the parameters the proposal is to be used for.

    Example::

        [jump_proposal-mtotal+q]
        name = adaptive_normal
        adaptation-duration = 1000
        min-q = 1
        max-q = 8
        min-mtotal = 20
        max-mtotal = 160

        [jump_proposal-spin1_a]
        name = normal

    Parameters
    ----------
    cp : WorkflowConfigParser instance
        The config file to read.
    section : str, optional
        The section name to read jump proposals from. Default is
        ``'jump_proposal'``.

    Returns
    -------
    list :
        List of the proposal instances.
    """
    tags = cp.get_subsections(section)
    proposals = []
    for tag in tags:
        # get the name of the proposal
        name = cp.get_opt_tag(section, "name", tag)
        prop = epsie_proposals[name].from_config(cp, section, tag)
        proposals.append(prop)
    return proposals
