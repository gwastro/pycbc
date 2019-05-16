# Copyright (C) 2018  Collin Capano
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


#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#
"""Provides constructor classes provide support for parallel tempered MCMC
samplers."""

from __future__ import absolute_import

from six import string_types

import numpy
from pycbc.filter import autocorrelation


class MultiTemperedSupport(object):
    """Provides methods for supporting multi-tempered samplers.
    """
    _ntemps = None

    @property
    def ntemps(self):
        """The number of temeratures that are set."""
        return self._ntemps


class MultiTemperedAutocorrSupport(object):
    """Provides class methods for calculating multi-tempered ACFs/ACLs.
    """

    @classmethod
    def compute_acf(cls, filename, start_index=None, end_index=None,
                    per_walker=False, walkers=None, parameters=None,
                    temps=None):
        """Computes the autocorrleation function of the model params in the
        given file.

        By default, parameter values are averaged over all walkers at each
        iteration. The ACF is then calculated over the averaged chain for each
        temperature. An ACF per-walker will be returned instead if
        ``per_walker=True``.

        Parameters
        -----------
        filename : str
            Name of a samples file to compute ACFs for.
        start_index : {None, int}
            The start index to compute the acl from. If None, will try to use
            the number of burn-in iterations in the file; otherwise, will start
            at the first sample.
        end_index : {None, int}
            The end index to compute the acl to. If None, will go to the end
            of the current iteration.
        per_walker : optional, bool
            Return the ACF for each walker separately. Default is False.
        walkers : optional, int or array
            Calculate the ACF using only the given walkers. If None (the
            default) all walkers will be used.
        parameters : optional, str or array
            Calculate the ACF for only the given parameters. If None (the
            default) will calculate the ACF for all of the model params.
        temps : optional, (list of) int or 'all'
            The temperature index (or list of indices) to retrieve. If None
            (the default), the ACF will only be computed for the coldest (= 0)
            temperature chain. To compute an ACF for all temperates pass 'all',
            or a list of all of the temperatures.

        Returns
        -------
        dict :
            Dictionary of arrays giving the ACFs for each parameter. If
            ``per-walker`` is True, the arrays will have shape
            ``ntemps x nwalkers x niterations``. Otherwise, the returned array
            will have shape ``ntemps x niterations``.
        """
        acfs = {}
        with cls._io(filename, 'r') as fp:
            if parameters is None:
                parameters = fp.variable_params
            if isinstance(parameters, string_types):
                parameters = [parameters]
            if isinstance(temps, int):
                temps = [temps]
            elif temps == 'all':
                temps = numpy.arange(fp.ntemps)
            elif temps is None:
                temps = [0]
            for param in parameters:
                subacfs = []
                for tk in temps:
                    if per_walker:
                        # just call myself with a single walker
                        if walkers is None:
                            walkers = numpy.arange(fp.nwalkers)
                        arrays = [cls.compute_acfs(filename,
                                                   start_index=start_index,
                                                   end_index=end_index,
                                                   per_walker=False,
                                                   walkers=ii,
                                                   parameters=param,
                                                   temps=tk)[param][0, :]
                                  for ii in walkers]
                        # we'll stack all of the walker arrays to make a single
                        # nwalkers x niterations array; when these are stacked
                        # below, we'll get a ntemps x nwalkers x niterations
                        # array
                        subacfs.append(numpy.vstack(arrays))
                    else:
                        samples = fp.read_raw_samples(
                            param, thin_start=start_index,
                            thin_interval=1, thin_end=end_index,
                            walkers=walkers, temps=tk, flatten=False)[param]
                        # contract the walker dimension using the mean, and
                        # flatten the (length 1) temp dimension
                        samples = samples.mean(axis=1)[0, :]
                        thisacf = autocorrelation.calculate_acf(
                            samples).numpy()
                        subacfs.append(thisacf)
                # stack the temperatures
                acfs[param] = numpy.stack(subacfs)
        return acfs

    @classmethod
    def compute_acl(cls, filename, start_index=None, end_index=None,
                    min_nsamples=10):
        """Computes the autocorrleation length for all model params and
        temperatures in the given file.

        Parameter values are averaged over all walkers at each iteration and
        temperature.  The ACL is then calculated over the averaged chain.

        Parameters
        -----------
        filename : str
            Name of a samples file to compute ACLs for.
        start_index : {None, int}
            The start index to compute the acl from. If None, will try to use
            the number of burn-in iterations in the file; otherwise, will start
            at the first sample.
        end_index : {None, int}
            The end index to compute the acl to. If None, will go to the end
            of the current iteration.
        min_nsamples : int, optional
            Require a minimum number of samples to compute an ACL. If the
            number of samples per walker is less than this, will just set to
            ``inf``. Default is 10.

        Returns
        -------
        dict
            A dictionary of ntemps-long arrays of the ACLs of each parameter.
        """
        acls = {}
        with cls._io(filename, 'r') as fp:
            if end_index is None:
                end_index = fp.niterations
            tidx = numpy.arange(fp.ntemps)
            for param in fp.variable_params:
                these_acls = numpy.zeros(fp.ntemps)
                for tk in tidx:
                    samples = fp.read_raw_samples(
                        param, thin_start=start_index, thin_interval=1,
                        thin_end=end_index, temps=tk, flatten=False)[param]
                    # contract the walker dimension using the mean, and flatten
                    # the (length 1) temp dimension
                    samples = samples.mean(axis=1)[0, :]
                    if samples.size < min_nsamples:
                        acl = numpy.inf
                    else:
                        acl = autocorrelation.calculate_acl(samples)
                    if acl <= 0:
                        acl = numpy.inf
                    these_acls[tk] = acl
                acls[param] = these_acls
        return acls
