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


from __future__ import absolute_import

import numpy

from epsie import proposals as epsie_proposals

from pycbc import VARARGS_DELIM
from pycbc import boundaries
from .normal import (load_opts, get_variance, get_param_boundaries,
                     get_epsie_adaptation_settings)


class EpsieBoundedNormal(epsie_proposals.BoundedNormal):
    """Adds ``from_config`` method to epsie's boundd normal proposal."""

    @classmethod
    def from_config(cls, cp, section, tag):
        r"""Loads a proposal from a config file.

        The section that is read should have the format ``[{section}-{tag}]``,
        where ``{tag}`` is a :py:const:`pycbc.VARARGS_DELIM` separated list
        of the parameters to create the jump proposal for.

        Boundaries must be provied for every must be provided. The syntax
        is ``(min|max)-{param} = float``. Variances for each parameter may also
        be specified, by giving options ``var-{param} = val``. Any parameter
        not specified will use a default variance of :math:`(\Delta p/10)^2`,
        where :math:`\Delta p` is the boundary width for parameter :math:`p`.

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
        :py:class:`epsie.proposals.Normal`:
            A normal proposal for use with ``epsie`` samplers.
        """
        # check that the name matches
        assert cp.get_opt_tag(section, "name", tag) == cls.name, (
            "name in specified section must match mine")
        params, opts = load_opts(cp, seciton, tag, skip=['name'])
        boundaries = get_param_boundaries(params, opts)
        if opts:
            cov = get_variance(params, opts)
        else:
            cov = numpy.array([abs(boundaries[p])/10. for p in params])**2.
        # no other options should remain
        if opts:
            raise ValueError("unrecognized options {}"
                             .format(', '.join(opts.keys())))
        return cls(params, boundaries, cov=cov)


class EpsieAdaptiveBoundedNormal(epsie_proposals.AdaptiveBoundedNormal):
    """Adds ``from_config`` method to epsie's adaptive normal proposal."""

    @classmethod
    def from_config(cls, cp, section, tag):
        """Loads a proposal from a config file.

        The section that is read should have the format ``[{section}-{tag}]``,
        where ``{tag}`` is a :py:const:`pycbc.VARARGS_DELIM` separated list
        of the parameters to create the jump proposal for.

        Options that are read:

        * name : str
            Required. Must match the name of the proposal.
        * adaptation-duration : int
            Required. Sets the ``adaptation_duration``.
        * min-{param} : float
        * max-{param} : float
            Required. Bounds must be provided for every parameter.
        * var-{param} : float
            Optional. Initial variance to use. If not provided, will use a
            default based on the bounds (see
            :py:class:`epsie.proposals.AdaptiveSupport` for details).
        * adaptation-decay : int
            Optional. Sets the ``adaptation_decay``. If not provided, will use
            the class's default.
        * start-iteration : int
            Optional. Sets the ``start_iteration``.If not provided, will use
            the class's default.
        * target-rate : float
            Optional. Sets the ``target_rate``. If not provided, will use
            the class's default.

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
        :py:class:`epsie.proposals.AdaptiveNormal`:
            An adaptive normal proposal for use with ``epsie`` samplers.
        """
        # check that the name matches
        assert cp.get_opt_tag(section, "name", tag) == cls.name, (
            "name in specified section must match mine")
        params, opts = load_opts(cp, seciton, tag, skip=['name'])
        args = {'parameters': params}
        # get the bounds
        args['boundaries'] = get_param_boundaries(params, opts)
        # get the adaptation parameters
        args.update(get_epsie_adaptation_settings(opts))
        # if there are any other options, assume they are for setting the
        # initial standard deviation
        if opts:
            var = get_variance(params, opts)
            args['initial_std'] = var**0.5
            # at this point, there should be no options left
            if opts:
                raise ValueError('unrecognized options {} in section {}'
                                 .format(', '.join(opts.keys()), readsection))
        return cls(**args)
