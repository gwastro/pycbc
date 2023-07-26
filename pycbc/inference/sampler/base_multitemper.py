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


import logging
import numpy
import h5py
from pycbc.filter import autocorrelation
from pycbc.inference.io import loadfile


class MultiTemperedSupport(object):
    """Provides methods for supporting multi-tempered samplers.
    """
    _ntemps = None

    @property
    def ntemps(self):
        """The number of temeratures that are set."""
        return self._ntemps

    @staticmethod
    def betas_from_config(cp, section):
        """Loads number of temperatures or betas from a config file.

        This looks in the given section for:

        * ``ntemps`` :
            The number of temperatures to use. Either this, or
            ``inverse-temperatures-file`` must be provided (but not both).
        * ``inverse-temperatures-file`` :
            Path to an hdf file containing the inverse temperatures ("betas")
            to use. The betas will be retrieved from the file's
            ``.attrs['betas']``. Either this or ``ntemps`` must be provided
            (but not both).

        Parameters
        ----------
        cp : WorkflowConfigParser instance
            Config file object to parse.
        section : str
            The name of the section to look in.

        Returns
        -------
        ntemps : int or None
            The number of temperatures to use, if it was provided.
        betas : array
            The array of betas to use, if a inverse-temperatures-file was
            provided.
        """
        if cp.has_option(section, "ntemps") and \
                cp.has_option(section, "inverse-temperatures-file"):
            raise ValueError("Must specify either ntemps or "
                             "inverse-temperatures-file, not both.")
        if cp.has_option(section, "inverse-temperatures-file"):
            # get the path of the file containing inverse temperatures values.
            inverse_temperatures_file = cp.get(section,
                                               "inverse-temperatures-file")
            betas = read_betas_from_hdf(inverse_temperatures_file)
            ntemps = betas.shape[0]
        else:
            # get the number of temperatures
            betas = None
            ntemps = int(cp.get(section, "ntemps"))
        return ntemps, betas


def read_betas_from_hdf(filename):
    """Loads inverse temperatures from the given file.
    """
    # get the path of the file containing inverse temperatures values.
    with h5py.File(filename, "r") as fp:
        try:
            betas = numpy.array(fp.attrs['betas'])
            # betas must be in decending order
            betas = numpy.sort(betas)[::-1]
        except KeyError:
            raise AttributeError("No attribute called betas")
    return betas


#
# =============================================================================
#
#              Functions for computing autocorrelation lengths
#
# =============================================================================
#


def compute_acf(filename, start_index=None, end_index=None,
                chains=None, parameters=None, temps=None):
    """Computes the autocorrleation function for independent MCMC chains with
    parallel tempering.

    Parameters
    -----------
    filename : str
        Name of a samples file to compute ACFs for.
    start_index : int, optional
        The start index to compute the acl from. If None (the default),
        will try to use the burn in iteration for each chain;
        otherwise, will start at the first sample.
    end_index : {None, int}
        The end index to compute the acl to. If None, will go to the end
        of the current iteration.
    chains : optional, int or array
        Calculate the ACF for only the given chains. If None (the
        default) ACFs for all chains will be estimated.
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
        Dictionary parameter name -> ACF arrays. The arrays have shape
        ``ntemps x nchains x niterations``.
    """
    acfs = {}
    with loadfile(filename, 'r') as fp:
        if parameters is None:
            parameters = fp.variable_params
        if isinstance(parameters, str):
            parameters = [parameters]
        temps = _get_temps_idx(fp, temps)
        if chains is None:
            chains = numpy.arange(fp.nchains)
        for param in parameters:
            subacfs = []
            for tk in temps:
                subsubacfs = []
                for ci in chains:
                    samples = fp.read_raw_samples(
                        param, thin_start=start_index, thin_interval=1,
                        thin_end=end_index, chains=ci, temps=tk)[param]
                    thisacf = autocorrelation.calculate_acf(samples).numpy()
                    subsubacfs.append(thisacf)
                # stack the chains
                subacfs.append(subsubacfs)
            # stack the temperatures
            acfs[param] = numpy.stack(subacfs)
    return acfs


def compute_acl(filename, start_index=None, end_index=None,
                min_nsamples=10):
    """Computes the autocorrleation length for independent MCMC chains with
    parallel tempering.

    ACLs are calculated separately for each chain.

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
        A dictionary of ntemps x nchains arrays of the ACLs of each
        parameter.
    """
    # following is a convenience function to calculate the acl for each chain
    # defined here so that we can use map for this below
    def _getacl(si):
        # si: the samples loaded for a specific chain; may have nans in it
        si = si[~numpy.isnan(si)]
        if len(si) < min_nsamples:
            acl = numpy.inf
        else:
            acl = autocorrelation.calculate_acl(si)
        if acl <= 0:
            acl = numpy.inf
        return acl
    acls = {}
    with loadfile(filename, 'r') as fp:
        tidx = numpy.arange(fp.ntemps)
        for param in fp.variable_params:
            these_acls = numpy.zeros((fp.ntemps, fp.nchains))
            for tk in tidx:
                samples = fp.read_raw_samples(
                    param, thin_start=start_index, thin_interval=1,
                    thin_end=end_index, temps=tk, flatten=False)[param]
                # flatten out the temperature
                samples = samples[0, ...]
                # samples now has shape nchains x maxiters
                if samples.shape[-1] < min_nsamples:
                    these_acls[tk, :] = numpy.inf
                else:
                    these_acls[tk, :] = list(map(_getacl, samples))
            acls[param] = these_acls
        # report the mean ACL: take the max over the temps and parameters
        act = acl_from_raw_acls(acls)*fp.thinned_by
        finite = act[numpy.isfinite(act)]
        logging.info("ACTs: min %s, mean (of finite) %s, max %s",
                     str(act.min()),
                     str(finite.mean() if finite.size > 0 else numpy.inf),
                     str(act.max()))
    return acls


def acl_from_raw_acls(acls):
    """Calculates the ACL for one or more chains from a dictionary of ACLs.

    This is for parallel tempered MCMCs in which the chains are independent
    of each other.

    The ACL for each chain is maximized over the temperatures and parameters.

    Parameters
    ----------
    acls : dict
        Dictionary of parameter names -> ntemps x nchains arrays of ACLs (the
        thing returned by :py:func:`compute_acl`).

    Returns
    -------
    array
        The ACL of each chain.
    """
    return numpy.array(list(acls.values())).max(axis=0).max(axis=0)


def ensemble_compute_acf(filename, start_index=None, end_index=None,
                         per_walker=False, walkers=None, parameters=None,
                         temps=None):
    """Computes the autocorrleation function for a parallel tempered, ensemble
    MCMC.

    By default, parameter values are averaged over all walkers at each
    iteration. The ACF is then calculated over the averaged chain for each
    temperature. An ACF per-walker will be returned instead if
    ``per_walker=True``.

    Parameters
    ----------
    filename : str
        Name of a samples file to compute ACFs for.
    start_index : int, optional
        The start index to compute the acl from. If None (the default), will
        try to use the number of burn-in iterations in the file; otherwise,
        will start at the first sample.
    end_index : int, optional
        The end index to compute the acl to. If None (the default), will go to
        the end of the current iteration.
    per_walker : bool, optional
        Return the ACF for each walker separately. Default is False.
    walkers : int or array, optional
        Calculate the ACF using only the given walkers. If None (the
        default) all walkers will be used.
    parameters : str or array, optional
        Calculate the ACF for only the given parameters. If None (the
        default) will calculate the ACF for all of the model params.
    temps : (list of) int or 'all', optional
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
    with loadfile(filename, 'r') as fp:
        if parameters is None:
            parameters = fp.variable_params
        if isinstance(parameters, str):
            parameters = [parameters]
        temps = _get_temps_idx(fp, temps)
        for param in parameters:
            subacfs = []
            for tk in temps:
                if per_walker:
                    # just call myself with a single walker
                    if walkers is None:
                        walkers = numpy.arange(fp.nwalkers)
                    arrays = [ensemble_compute_acf(filename,
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


def ensemble_compute_acl(filename, start_index=None, end_index=None,
                         min_nsamples=10):
    """Computes the autocorrleation length for a parallel tempered, ensemble
    MCMC.

    Parameter values are averaged over all walkers at each iteration and
    temperature.  The ACL is then calculated over the averaged chain.

    Parameters
    -----------
    filename : str
        Name of a samples file to compute ACLs for.
    start_index : int, optional
        The start index to compute the acl from. If None (the default), will
        try to use the number of burn-in iterations in the file; otherwise,
        will start at the first sample.
    end_index : int, optional
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
    with loadfile(filename, 'r') as fp:
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
        maxacl = numpy.array(list(acls.values())).max()
        logging.info("ACT: %s", str(maxacl*fp.thinned_by))
    return acls


def _get_temps_idx(fp, temps):
    """Gets the indices of temperatures to load for computing ACF.
    """
    if isinstance(temps, int):
        temps = [temps]
    elif temps == 'all':
        temps = numpy.arange(fp.ntemps)
    elif temps is None:
        temps = [0]
    return temps
