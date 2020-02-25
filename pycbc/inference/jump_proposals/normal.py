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

"""Jump proposals that use a normal distribution."""

from __future__ import absolute_import

import numpy

from epsie import proposals as epsie_proposals

from pycbc import VARARGS_DELIM
from pycbc import boundaries


class EpsieNormal(epsie_proposals.Normal):
    """Adds ``from_config`` method to epsie's normal proposal."""

    @classmethod
    def from_config(cls, cp, section, tag):
        """Loads a proposal from a config file.

        The section that is read should have the format ``[{section}-{tag}]``,
        where ``{tag}`` is a :py:const:`pycbc.VARARGS_DELIM` separated list
        of the parameters to create the jump proposal for.

        Variances for each parameter may also be specified, by giving options
        ``var-{param} = val``. Any parameter not specified will use a default
        variance of 1.

        Example::

            [jump_proposal-mchrip+q]
            name = normal
            var-q = 0.1

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
        params = tag.split(VARARGS_DELIM)
        # see if any variances were provided
        readsection = '-'.join([section, tag])
        opts = [opt for opt in cp.options(readsection) if opt != 'name']
        varfmt = 'var-{}'
        if opts:
            optvals = {opt: cp.get(readsection, opt) for opt in opts}
            cov = numpy.array([float(optvals.pop(varfmt.format(param), 1.))
                               for param in params])
            # check that there are no unrecognized options
            if optvals:
                raise ValueError("unrecognized options {}"
                                 .format(', '.join(optvals.keys())))
        else:
            cov = None
        return cls(params, cov=cov)


class EpsieAdaptiveNormal(epsie_proposals.AdaptiveNormal):
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
            Required. Bounds must be provided for every parameter. These are
            used to determine the prior widths.
        * adaptation-decay : int
            Optional. Sets the ``adaptation_decay``. If not provided, will use
            the class's default.
        * start-iteration : int
            Optional. Sets the ``start_iteration``.If not provided, will use
            the class's default.
        * target-rate : float
            Optional. Sets the ``target_rate``. If not provided, will use
            the class's default.

        .. note::
           The min and max parameter bounds are only used for setting the width
           of the covariance of the proposal; they are not used as bounds on
           the proposal itself. In other words, it is possible to get proposals
           outside of the given min and max values.

        Example::

            [jump_proposal-mchirp+q]
            name = adaptive_normal
            adaptation-duration = 1000
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
        :py:class:`epsie.proposals.AdaptiveNormal`:
            An adaptive normal proposal for use with ``epsie`` samplers.
        """
        # check that the name matches
        assert cp.get_opt_tag(section, "name", tag) == cls.name, (
            "name in specified section must match mine")
        params = tag.split(VARARGS_DELIM)
        # get options
        readsection = '-'.join([section, tag])
        opts = {opt: cp.get(readsection, opt)
                for opt in cp.options(readsection) if opt != 'name'}
        adaptation_duration = opts.pop('adaptation-duration', None)
        if adaptation_duration is None:
            raise ValueError("No adaptation-duration specified")
        adaptation_duration = int(adaptation_duration)
        # get the bounds
        prior_widths = {}
        for param in params:
            minbound = opts.pop('min-{}'.format(param), None)
            if minbound is None:
                raise ValueError("Must provide a minimum bound for {p}."
                                 "Syntax is min-{p} = val".format(p=param))
            maxbound = opts.pop('max-{}'.format(param), None)
            if maxbound is None:
                raise ValueError("Must provide a maximum bound for {p}."
                                 "Syntax is max-{p} = val".format(p=param))
            prior_widths[param] = boundaries.Bounds(float(minbound),
                                                    float(maxbound))
        # optional args
        optional_args = {}
        adaptation_decay = opts.pop('adaptation-decay', None)
        if adaptation_decay is not None:
            optional_args['adaptation_decay'] = int(adaptation_decay)
        start_iteration = opts.pop('start-iteration', None)
        if start_iteration is not None:
            optional_args['start_iteration'] = int(start_iteration)
        target_rate = opts.pop('target_rate', None)
        if target_rate is not None:
            optional_args['target_rate'] = float(target_rate)
        # check that there are no unrecognized options
        if opts:
            raise ValueError('unrecognized options {} in section {}'
                             .format(', '.join(opts.keys()), readsection))
        return cls(params, prior_widths, adaptation_duration, **optional_args)
