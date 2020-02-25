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

from __future__ import absolute_import


from .normal import (EpsieNormal, EpsieAdaptiveNormal)


epsie_proposals = {
    EpsieNormal.name: EpsieNormal,
    EpsieAdaptiveNormal.name: EpsieAdaptiveNormal
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
    dict :
        Dictionary mapping parameter names to proposal instances.
    """
    tags = cp.get_subsections(section)
    proposals = {}
    for tag in tags:
        # get the name of the proposal
        name = cp.get_opt_tag(section, "name", tag)
        prop = epsie_proposals[name].from_config(cp, section, tag)
        proposals[prop.parameters] = prop
    return proposals
