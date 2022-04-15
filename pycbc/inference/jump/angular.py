# Copyright (C) 2020  Collin Capano
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

"""Jump proposals that use cyclic boundaries on [0, 2pi)."""

from epsie import proposals as epsie_proposals

from .normal import (epsie_from_config, epsie_adaptive_from_config,
                     epsie_at_adaptive_from_config)


class EpsieAngular(epsie_proposals.Angular):
    """Adds ``from_config`` method to epsie's angular proposal."""

    @classmethod
    def from_config(cls, cp, section, tag):
        """Loads a proposal from a config file.

        This calls :py:func:`epsie_from_config` with ``cls`` set to
        :py:class:`epsie.proposals.Angular` and ``with_boundaries`` set
        to False. See that function for details on options that can be read.

        Example::

            [jump_proposal-ra]
            name = angular
            var-ra = 0.01

        Parameters
        ----------
        cp : WorkflowConfigParser instance
            Config file to read from.
        section : str
            The name of the section to look in.
        tag : str
            :py:const:`pycbc.VARARGS_DELIM` separated list of parameter names
            to create proposals for.

        Returns
        -------
        :py:class:`epsie.proposals.Angular`:
            An angular proposal for use with ``epsie`` samplers.
        """
        return epsie_from_config(cls, cp, section, tag, with_boundaries=False)


class EpsieAdaptiveAngular(epsie_proposals.AdaptiveAngular):
    """Adds ``from_config`` method to epsie's adaptive angular proposal."""

    @classmethod
    def from_config(cls, cp, section, tag):
        r"""Loads a proposal from a config file.

        This calls :py:func:`epsie_adaptive_from_config` with ``cls`` set to
        :py:class:`epsie.proposals.AdaptiveBoundedNormal` and
        ``with_boundaries`` set to False (since the boundaries for the angular
        proposals are always :math:`[0, 2\pi)`). See that function
        for details on options that can be read.

        Example::

            [jump_proposal-ra]
            name = adaptive_angular
            adaptation-duration = 1000

        Parameters
        ----------
        cp : WorkflowConfigParser instance
            Config file to read from.
        section : str
            The name of the section to look in.
        tag : str
            :py:const:`pycbc.VARARGS_DELIM` separated list of parameter names
            to create proposals for.

        Returns
        -------
        :py:class:`epsie.proposals.AdaptiveAngular`:
            An adaptive angular proposal for use with ``epsie`` samplers.
        """
        return epsie_adaptive_from_config(cls, cp, section, tag,
                                          with_boundaries=False)


class EpsieATAdaptiveAngular(epsie_proposals.ATAdaptiveAngular):
    """Adds ``from_config`` method to epsie's adaptive angular proposal."""

    @classmethod
    def from_config(cls, cp, section, tag):
        r"""Loads a proposal from a config file.

        This calls :py:func:`epsie_adaptive_from_config` with ``cls`` set to
        :py:class:`epsie.proposals.AdaptiveBoundedNormal` and
        ``with_boundaries`` set to False (since the boundaries for the angular
        proposals are always :math:`[0, 2\pi)`). See that function
        for details on options that can be read.

        Example::

            [jump_proposal-ra]
            name = adaptive_angular_proposal
            adaptation-duration = 1000

        Parameters
        ----------
        cp : WorkflowConfigParser instance
            Config file to read from.
        section : str
            The name of the section to look in.
        tag : str
            :py:const:`pycbc.VARARGS_DELIM` separated list of parameter names
            to create proposals for.

        Returns
        -------
        :py:class:`epsie.proposals.AdaptiveAngularProposal`:
            An adaptive angular proposal for use with ``epsie`` samplers.
        """
        return epsie_at_adaptive_from_config(cls, cp, section, tag,
                                             with_boundaries=False)
