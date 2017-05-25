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


#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#
"""
This modules provides classes and functions for using the kombine sampler
packages for parameter estimation.
"""

import numpy
from pycbc.inference.sampler_base import BaseMCMCSampler

#
# =============================================================================
#
#                                   Samplers
#
# =============================================================================
#


class KombineSampler(BaseMCMCSampler):
    """This class is used to construct the MCMC sampler from the kombine
    package.

    Parameters
    ----------
    likelihood_evaluator : likelihood class
        An instance of the likelihood class from the
        pycbc.inference.likelihood module.
    nwalkers : int
        Number of walkers to use in sampler.
    transd : bool
        If True, the sampler will operate across parameter spaces using a
        kombine.clustered_kde.TransdimensionalKDE proposal distribution. In
        this mode a masked array with samples in each of the possible sets of
        dimensions must be given for the initial ensemble distribution.
    processes : {None, int}
        Number of processes to use with multiprocessing. If None, all available
        cores are used.
    min_burn_in : {None, int}
        Set the minimum number of burn in iterations to use. If None,
        `burn_in_iterations` will be initialized to `0`.
    update_interval : {None, int}
        Make the sampler update the proposal densities every `update_interval`
        iterations.
    """
    name = "kombine"

    def __init__(self, likelihood_evaluator, nwalkers, transd=False,
                 min_burn_in=None, pool=None, likelihood_call=None,  
                 update_interval=None):

        try:
            import kombine
        except ImportError:
            raise ImportError("kombine is not installed.")

        if likelihood_call is None:
            likelihood_call = likelihood_evaluator

        # construct sampler for use in KombineSampler
        ndim = len(likelihood_evaluator.waveform_generator.variable_args)
        sampler = kombine.Sampler(nwalkers, ndim, likelihood_call,
                                  transd=transd, pool=pool,
                                  processes=pool.count)
        # initialize
        super(KombineSampler, self).__init__(sampler, likelihood_evaluator,
                                             min_burn_in=min_burn_in)
        self._nwalkers = nwalkers
        self.update_interval = update_interval

    @classmethod
    def from_cli(cls, opts, likelihood_evaluator, pool=None, likelihood_call=None):
        """Create an instance of this sampler from the given command-line
        options.

        Parameters
        ----------
        opts : ArgumentParser options
            The options to parse.
        likelihood_evaluator : LikelihoodEvaluator
            The likelihood evaluator to use with the sampler.

        Returns
        -------
        KombineSampler
            A kombine sampler initialized based on the given arguments.
        """
        return cls(likelihood_evaluator, opts.nwalkers, likelihood_call=likelihood_call,
                   min_burn_in=opts.min_burn_in, pool=pool, update_interval=opts.update_interval)

    def run(self, niterations, **kwargs):
        """Advance the sampler for a number of samples.

        Parameters
        ----------
        niterations : int
            Number of samples to get from sampler.

        Returns
        -------
        p : numpy.array
            An array of current walker positions with shape (nwalkers, ndim).
        lnpost : numpy.array
            The list of log posterior probabilities for the walkers at
            positions p, with shape (nwalkers, ndim).
        lnprop : numpy.array
            The list of log proposal densities for the walkers at positions p,
            with shape (nwalkers, ndim).
        """
        blob0 = None
        if self.niterations == 0:
            # first time running, use the initial positions
            p0 = self.p0
            if self.likelihood_evaluator.return_meta:
                blob0 = [self.likelihood_evaluator(p0[wi, :])[1]
                         for wi in range(self.nwalkers)]
        else:
            p0 = None
            # kombine requires blob data to be specified
            blob0 = self._currentblob
        kwargs['blob0'] = blob0
        if 'update_interval' not in kwargs:
            # use the internal update interval
            kwargs['update_interval'] = self.update_interval
        res = self._sampler.run_mcmc(niterations, p0=p0, **kwargs)
        p, lnpost, lnprop = res[0], res[1], res[2]
        # update the positions
        self._pos = p
        if self.likelihood_evaluator.return_meta:
            self._currentblob = self._sampler.blobs[-1]
        return p, lnpost, lnprop

    @property
    def lnpost(self):
        """ Get the natural logarithm of the likelihood as an
        nwalkers x niterations array.
        """
        # kombine returns niterations x nwaklers
        return self._sampler.lnpost.transpose()

    @property
    def chain(self):
        """Get all past samples as an nwalker x niterations x ndim array."""
        # kombine returns niterations x nwalkers x ndim
        return self._sampler.chain.transpose((1, 0, 2))

    def clear_chain(self):
        """Clears the chain and blobs from memory.
        """
        # store the iteration that the clear is occuring on
        self._lastclear = self.niterations
        # kombine stores its chain as niterations x nwalkers x ndim
        current_shape = self._sampler._chain.shape
        new_shape = (0, current_shape[1], current_shape[2])
        if isinstance(self._sampler._chain, numpy.ma.MaskedArray):
            self._sampler._chain = numpy.ma.resize(self._sampler._chain,
                                                   new_shape)
        else:
            self._sampler._chain.resize(new_shape)
        self._sampler.stored_iterations = 0
        # clear the blobs
        self._sampler._blobs = []

    def burn_in(self):
        """Use kombine's `burnin` routine to advance the sampler.

        If a minimum number of burn-in iterations was specified, this will run
        the burn-in until it has advanced at least as many steps as desired.
        The initial positions (p0) must be set prior to running.

        For more details, see `kombine.sampler.burnin`.

        Returns
        -------
        p : numpy.array
            An array of current walker positions with shape (nwalkers, ndim).
        lnpost : numpy.array
            The list of log posterior probabilities for the walkers at
            positions p, with shape (nwalkers, ndim).
        lnprop : numpy.array
            The list of log proposal densities for the walkers at positions p,
            with shape (nwalkers, ndim).
        """
        # check that we haven't already burned in
        if self.pos is not None:
            raise ValueError("burn-in already run")
        # run once
        p0 = self.p0
        if self.likelihood_evaluator.return_meta:
            blob0 = [self.likelihood_evaluator(p0[wi, :])[1]
                     for wi in range(self.nwalkers)]
        else:
            blob0 = None
        res = self._sampler.burnin(self.p0, blob0=blob0)
        p, post, q = res[0], res[1], res[2]
        # continue running until minimum burn in is satistfied
        while self.niterations < self.burn_in_iterations:
            p0 = p
            res = self._sampler.burnin(p0)
            p, post, q = res[0], res[1], res[2]
            # update position
            self._pos = p
            self._currentblob = self._sampler.blobs[-1]
        self.burn_in_iterations = self.niterations
        return p, post, q

    def _write_kde(self, fp, dataset_name, kde):
        """Writes the given kde to the file."""
        shape = kde.data.shape
        try:
            if shape != fp[dataset_name].shape:
                # resize the dataset
                fp[dataset_name].resize(shape)
            fp[dataset_name][:] = kde.data
        except KeyError:
            # dataset doesn't exist yet
            fp.create_dataset(dataset_name, shape,
                              maxshape=(self.nwalkers,
                                        len(self.variable_args)),
                              dtype=float)
            fp[dataset_name][:] = kde.data

    def write_state(self, fp):
        """ Saves the state of the sampler in a file.

        Parameters
        ----------
        fp : InferenceFile
            File to store sampler state.
        """

        # save clustered KDE data
        subgroup = "clustered_kde"
        dataset_name ="/".join([fp.sampler_group, subgroup])
        clustered_kde = self._sampler._kde
        self._write_kde(fp, dataset_name, clustered_kde)
        # metadata
        fp[dataset_name].attrs["nclusters"] = clustered_kde.nclusters
        fp[dataset_name].attrs["assignments"] = clustered_kde._assignments
        fp[dataset_name].attrs["centroids"] = clustered_kde.centroids
        fp[dataset_name].attrs["logweights"] = clustered_kde._logweights
        fp[dataset_name].attrs["mean"] = clustered_kde._mean
        fp[dataset_name].attrs["std"] = clustered_kde._std
        # save individual KDE data
        for i, kde in enumerate(clustered_kde._kdes):
            dataset_name = "/".join([fp.sampler_group, "kde" + str(i)])
            self._write_kde(fp, dataset_name, kde)

    def set_state_from_file(self, fp):
        """ Sets the state of the sampler back to the instance saved in a file.

        Parameters
        ----------
        fp : InferenceFile
            File with sampler state stored.
        """

        try:
            import kombine
        except ImportError:
            raise ImportError("kombine is not installed.")

        # create a ClusteredKDE
        dataset_name = "/".join([fp.sampler_group, "clustered_kde"])
        clustered_kde = kombine.clustered_kde.ClusteredKDE(
                                           fp[dataset_name][:],
                                           fp[dataset_name].attrs["nclusters"])
        clustered_kde._assignments = fp[dataset_name].attrs["assignments"]
        clustered_kde._centroids = fp[dataset_name].attrs["centroids"]
        clustered_kde.logweights = fp[dataset_name].attrs["logweights"]
        clustered_kde._mean = fp[dataset_name].attrs["mean"]
        clustered_kde._std = fp[dataset_name].attrs["std"]

        # add KDEs
        clustered_kde._kdes = []
        for i in range(clustered_kde.nclusters):
            dataset_name = "/".join([fp.sampler_group, "kde" + str(i)])
            clustered_kde._kdes.append(
                                kombine.clustered_kde.KDE(fp[dataset_name][:]))

        # overwrite ClusteredKDE in Sampler instance
        self._sampler._kde = clustered_kde
