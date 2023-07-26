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

from .normal import (epsie_from_config, epsie_adaptive_from_config,
                     epsie_at_adaptive_from_config)


class EpsieBoundedNormal(epsie_proposals.BoundedNormal):
    """Adds ``from_config`` method to epsie's bounded normal proposal."""

    @classmethod
    def from_config(cls, cp, section, tag):
        r"""Loads a proposal from a config file.

        This calls :py:func:`epsie_from_config` with ``cls`` set to
        :py:class:`epsie.proposals.BoundedNormal` and ``with_boundaries`` set
        to True. See that function for details on options that can be read.

        Example::

            [jump_proposal-mchrip+q]
            name = bounded_normal
            min-q = 1
            max-q = 8
            min-mchirp = 20
            max-mchirp = 80

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
        :py:class:`epsie.proposals.BoundedNormal`:
            A bounded normal proposal for use with ``epsie`` samplers.
        """
        return epsie_from_config(cls, cp, section, tag, with_boundaries=True)


class EpsieAdaptiveBoundedNormal(epsie_proposals.AdaptiveBoundedNormal):
    """Adds ``from_config`` method to epsie's adaptive normal proposal."""

    @classmethod
    def from_config(cls, cp, section, tag):
        """Loads a proposal from a config file.

        This calls :py:func:`epsie_adaptive_from_config` with ``cls`` set to
        :py:class:`epsie.proposals.AdaptiveBoundedNormal`. See that function
        for details on options that can be read.

        Example::

            [jump_proposal-q]
            name = adaptive_bounded_normal
            adaptation-duration = 1000
            min-q = 1
            max-q = 8

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
        :py:class:`epsie.proposals.AdaptiveBoundedNormal`:
            An adaptive normal proposal for use with ``epsie`` samplers.
        """
        return epsie_adaptive_from_config(cls, cp, section, tag)


class EpsieATAdaptiveBoundedNormal(epsie_proposals.ATAdaptiveBoundedNormal):
    """Adds ``from_config`` method to epsie's adaptive bounded proposal."""

    @classmethod
    def from_config(cls, cp, section, tag):
        """Loads a proposal from a config file.

        This calls :py:func:`epsie_adaptive_from_config` with ``cls`` set to
        :py:class:`epsie.proposals.AdaptiveBoundedProposal`. See that function
        for details on options that can be read.

        Example::

            [jump_proposal-q]
            name = adaptive_bounded_proposal
            min-q = 1
            max-q = 8

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
        :py:class:`epsie.proposals.AdaptiveBoundedProposal`:
            An adaptive bounded proposal for use with ``epsie`` samplers.
        """
        return epsie_at_adaptive_from_config(cls, cp, section, tag,
                                             with_boundaries=True)
