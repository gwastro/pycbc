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

"""Jump proposals that use a bounded normal distribution."""

from epsie import proposals as epsie_proposals

from .normal import (epsie_from_config, epsie_adaptive_from_config)


class EpsieNormalDiscrete(epsie_proposals.NormalDiscrete):
    """Adds ``from_config`` method to epsie's normal discrete proposal."""

    @classmethod
    def from_config(cls, cp, section, tag):
        r"""Loads a proposal from a config file.

        This calls :py:func:`epsie_from_config` with ``cls`` set to
        :py:class:`epsie.proposals.NormalDiscrete` and ``with_boundaries`` set
        to False. See that function for details on options that can be read.

        Example::

            [jump_proposal-index]
            name = discrete

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
        :py:class:`epsie.proposals.BoundedDiscrete`:
            A bounded discrete proposal for use with ``epsie`` samplers.
        """
        return epsie_from_config(cls, cp, section, tag, with_boundaries=False)


class EpsieBoundedDiscrete(epsie_proposals.BoundedDiscrete):
    """Adds ``from_config`` method to epsie's bounded discrete proposal."""

    @classmethod
    def from_config(cls, cp, section, tag):
        r"""Loads a proposal from a config file.

        This calls :py:func:`epsie_from_config` with ``cls`` set to
        :py:class:`epsie.proposals.BoundedDiscrete` and ``with_boundaries`` set
        to True. See that function for details on options that can be read.

        Example::

            [jump_proposal-index]
            name = bounded_discrete
            min-index = 0
            max-index = 19

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
        :py:class:`epsie.proposals.BoundedDiscrete`:
            A bounded discrete proposal for use with ``epsie`` samplers.
        """
        return epsie_from_config(cls, cp, section, tag, with_boundaries=True)


class EpsieAdaptiveNormalDiscrete(epsie_proposals.AdaptiveNormalDiscrete):
    """Adds ``from_config`` method to epsie's adaptive bounded discrete
    proposal."""

    @classmethod
    def from_config(cls, cp, section, tag):
        """Loads a proposal from a config file.

        This calls :py:func:`epsie_adaptive_from_config` with ``cls`` set to
        :py:class:`epsie.proposals.AdaptiveNormalDiscrete`. See that function
        for details on options that can be read.

        Example::

            [jump_proposal-index]
            name = adaptive_normal_discrete
            adaptation-duration = 1000
            min-index = 0
            max-index = 42

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
        :py:class:`epsie.proposals.AdaptiveBoundedDiscrete`:
            An adaptive normal proposal for use with ``epsie`` samplers.
        """
        return epsie_adaptive_from_config(cls, cp, section, tag,
                                          boundary_arg_name='prior_widths')


class EpsieAdaptiveBoundedDiscrete(epsie_proposals.AdaptiveBoundedDiscrete):
    """Adds ``from_config`` method to epsie's adaptive bounded discrete
    proposal."""

    @classmethod
    def from_config(cls, cp, section, tag):
        """Loads a proposal from a config file.

        This calls :py:func:`epsie_adaptive_from_config` with ``cls`` set to
        :py:class:`epsie.proposals.AdaptiveBoundedDiscrete`. See that function
        for details on options that can be read.

        Example::

            [jump_proposal-index]
            name = adaptive_bounded_discrete
            adaptation-duration = 1000
            min-index = 0
            max-index = 42

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
        :py:class:`epsie.proposals.AdaptiveBoundedDiscrete`:
            An adaptive normal proposal for use with ``epsie`` samplers.
        """
        return epsie_adaptive_from_config(cls, cp, section, tag)
