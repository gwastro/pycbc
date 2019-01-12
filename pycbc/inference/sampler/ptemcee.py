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


"""
This modules provides classes and functions for using the ptemcee sampler
packages for parameter estimation.
"""

from __future__ import absolute_import

import numpy
import ptemcee
import logging
from pycbc.pool import choose_pool

from .base import BaseSampler
from .base_mcmc import (BaseMCMC, raw_samples_to_dict,
                        get_optional_arg_from_config)
from .base_multitemper import (MultiTemperedSupport,
                               MultiTemperedAutocorrSupport)
from ..burn_in import MultiTemperedMCMCBurnInTests
from pycbc.inference.io import PTEmceeFile
from .. import models


class PTEmceeSampler(MultiTemperedAutocorrSupport, MultiTemperedSupport,
                     BaseMCMC, BaseSampler):
    """This class is used to construct a parallel-tempered MCMC sampler from
    ptemcee.

    Parameters
    ----------
    model : model
        A model from ``pycbc.inference.models``.
    nwalkers : int
        Number of walkers to use in sampler.
    loglikelihood_function : str, optional
        The name of the likelihood function to use. If None (the default),
        will use 'loglikelihood'.
    ntemps : int, optional
        Specify the number of temps to use. Either this, ``Tmax``, or ``betas``
        must be specified.
    Tmax : float, optional
        Specify the maximum temperature to use. This may be used with
        ``ntemps``; see :py:func:`ptemcee.make_ladder` for details. Either
        this, ``ntemps``, or ``betas`` must be specified.
    betas : list of float, optional
        Specify the betas to use. Must be provided if ``ntemps`` and ``Tmax``
        are not given. Will override ``ntemps`` and ``Tmax`` if provided.
    adaptive : bool, optional
        Whether or not to use adaptive temperature levels. Default is True.
    adaptation_lag : int, optional
        See :py:mod:`ptemcee.Sampler` for details. Default is 1000.
    adaptation_time : int, optional
        See :py:mod:`ptemcee.Sampler` for details. Default is 100.
    scale_factor : float, optional
        See :py:mod:`ptemcee.Sampler` for details.
    checkpoint_interval : int, optional
        The number of iterations to run for before checkpointing. Pass None
        (the default) to not do any checkpointing.
    nprocesses : int, optional
        The number of processes to use. Default is 1.
    use_mpi : bool, optional
        Use an MPI for parallelization.
    """
    name = "ptemcee"
    _io = PTEmceeFile
    burn_in_class = MultiTemperedMCMCBurnInTests

    def __init__(self, model, nwalkers, ntemps=None, Tmax=None, betas=None,
                 loglikelihood_function=None,
                 adaptive=True, adaptation_lag=1000, adaptation_time=100,
                 scale_factor=None,
                 checkpoint_interval=None,
                 nprocesses=1, use_mpi=False):

        self.model = model
        ndim = len(model.variable_params)
        # create a wrapper for calling the model
        if loglikelihood_function is None:
            loglikelihood_function = 'loglikelihood'
        # frustratingly, ptemcee does not support blob data, so we have to
        # turn it off
        model_call = models.CallModel(model, loglikelihood_function,
                                      return_all_stats=False)
        # create temperature ladder if needed
        if ntemps is None and Tmax is None and betas is None:
            raise ValueError("must provide either ntemps/Tmax or betas")
        if betas is None:
            betas = ptemcee.make_ladder(ndim, ntemps=ntemps, Tmax=Tmax)
        # construct the keyword arguments to pass; if a kwarg is None, we
        # won't pass it, resulting in ptemcee's defaults being used
        kwargs = {}
        kwargs['adaptive'] = adaptive
        kwargs['betas'] = betas
        kwargs['adaptation_lag'] = adaptation_lag
        kwargs['adaptation_time'] = adaptation_time
        if scale_factor is not None:
            kwargs['scale_factor'] = scale_factor
        # Set up the pool
        if nprocesses > 1:
            # these are used to help paralleize over multiple cores / MPI
            models._global_instance = model_call
            model_call = models._call_global_model
            prior_call = models._call_global_model_logprior
        else:
            prior_call = models.CallModel(model, 'logprior',
                                          return_all_stats=False)
        pool = choose_pool(mpi=use_mpi, processes=nprocesses)
        if pool is not None:
            pool.count = nprocesses
        self._sampler = ptemcee.Sampler(nwalkers=nwalkers, ndim=ndim,
                                        logl=model_call, logp=prior_call,
                                        mapper=pool.map, **kwargs)
        self._nwalkers = nwalkers
        self._ntemps = len(self._sampler.betas)
        self._checkpoint_interval = checkpoint_interval
        # we'll initialize ensemble and chain to None
        self._chain = None
        self._ensemble = None

    @property
    def io(self):
        return self._io

    @property
    def base_shape(self):
        return (self.ntemps, self.nwalkers,)

    @property
    def betas(self):
        # chain betas has shape niterations x ntemps; transpose to
        # ntemps x niterations
        return self._chain.betas.transpose()

    @property
    def starting_betas(self):
        # the initial betas that were used
        return self._sampler.betas

    @property
    def adaptive(self):
        return self._sampler.adaptive

    @property
    def adaptation_lag(self):
        return self._sampler.adaptation_lag

    @property
    def adaptation_time(self):
        return self._sampler.adaptation_time

    @property
    def scale_factor(self):
        return self._sampler.scale_factor

    @property
    def ensemble(self):
        """Returns the current ptemcee ensemble.

        The ensemble stores the current location of and temperatures of the
        walkers. If the ensemble hasn't been setup yet, will set one up
        using p0 for the positions. If set_p0 hasn't been run yet, this will
        result in a ValueError.
        """
        if self._ensemble is None:
            if self._p0 is None:
                raise ValueError("initial positions not set; run set_p0")
            # use the global numpy random state
            rstate = numpy.random.mtrand._rand
            # self._p0 has base_shape x ndim = ntemps x nwalkers x ndim (see
            # BaseMCMC.set_p0). ptemcee's Ensemble expects
            # ntemps x nwalkers x ndim... so we're good
            self._ensemble = self._sampler.ensemble(self._p0, rstate)
        return self._ensemble

    @property
    def _pos(self):
        """Uses the ensemble for the position."""
        # BaseMCMC expects _pos to have shape ntemps x nwalkers x ndim,
        # which is the same shape as ensemble.x
        return self.ensemble.x

    @property
    def chain(self):
        """The current chain of samples in memory.

        The chain is returned as a :py:mod:`ptemcee.chain.Chain` instance. If
        no chain has been created yet (``_chain`` is None), then will create
        a new chain using the current ``ensemble``.
        """
        if self._chain is None:
            # create a chain
            self._chain = ptemcee.chain.Chain(self.ensemble)
        return self._chain

    def clear_samples(self):
        """Clears the chain and blobs from memory.
        """
        # store the iteration that the clear is occuring on
        self._lastclear = self.niterations
        self._itercounter = 0
        # set _chain to None; this will both cause the current chain to
        # get garbage collected, and will cause a new chain to be created
        # the next time self.chain is called
        self._chain = None

    @property
    def samples(self):
        """A dict mapping ``variable_params`` to arrays of samples currently
        in memory.

        The arrays have shape ``ntemps x nwalkers x niterations``.
        """
        # chain.x has shape niterations x ntemps x nwalkers x ndim
        # we'll transpose to ntemps x nwalkers x niterations x ndim
        raw_samples = self._chain.x.transpose((1, 2, 0, 3))
        return raw_samples_to_dict(self, raw_samples)

    @property
    def model_stats(self):
        """Returns the log likelihood ratio and log prior as a dict of arrays.

        The returned array has shape ntemps x nwalkers x niterations.

        Unfortunately, because ``ptemcee`` does not have blob support, this
        will only return the loglikelihood, logprior, and logjacobian,
        regardless of what stats the model can return.
        """
        # log likelihood and prior have shape
        # niterations x ntemps x nwalkers; we'll tranpose to have shape
        # ntemps x nwalkers x niterations
        logl = self._chain.logl.transpose((1, 2, 0))
        logp = self._chain.logP.transpose((1, 2, 0))
        logjacobian = numpy.zeros(logp.size)
        # if different coordinates were used for sampling, get the jacobian
        if self.model.sampling_transforms is not None:
            samples = self.samples
            flattened_samples = {param: arr.ravel()
                                 for param, arr in samples.items()}
            for ii in range(logp.size):
                these_samples = {param: vals[ii]
                                 for param, vals in flattened_samples.items()}
                self.model.update(**these_samples)
                logjacobian[ii] = self.model.logjacobian
        logjacobian = logjacobian.reshape(logp.shape)
        # put the logprior into the variable_params space
        logp -= logjacobian
        return {'loglikelihood': logl, 'logprior': logp,
                'logjacobian': logjacobian}

    def set_state_from_file(self, filename):
        """Sets the state of the sampler back to the instance saved in a file.
        """
        with self.io(filename, 'r') as fp:
            rstate = fp.read_random_state()
            # set the numpy random state
            numpy.random.set_state(rstate)
            # set the ensemble to its last state
            ensemble = self.ensemble
            for attr, val in fp.read_ensemble_attrs().items():
                setattr(ensemble, attr, val)
            ensemble.betas = fp.read_betas(iteration=-1)
            ensemble.time = fp.niterations

    def run_mcmc(self, niterations):
        """Advance the ensemble for a number of samples.

        Parameters
        ----------
        niterations : int
            Number of samples to get from sampler.
        """
        self.chain.run(niterations)

    def write_results(self, filename):
        """Writes samples, model stats, acceptance fraction, and random state
        to the given file.

        Parameters
        -----------
        filename : str
            The file to write to. The file is opened using the ``io`` class
            in an an append state.
        """
        with self.io(filename, 'a') as fp:
            # write samples
            fp.write_samples(self.samples, self.model.variable_params)
            # write stats
            fp.write_samples(self.model_stats)
            # write betas
            fp.write_betas(self.betas)
            # write random state
            fp.write_random_state()
            # write attributes of the ensemble
            fp.write_ensemble_attrs(self.ensemble)

    @classmethod
    def calculate_logevidence(cls, filename, thin_start=None, thin_end=None,
                              thin_interval=None):
        """Calculates the log evidence from the given file.

        This uese ``ptemcee``'s thermodynamic integration.

        Parameters
        ----------
        filename : str
            Name of the file to read the samples from. Should be an
            ``PTEmceeFile``.
        thin_start : int
            Index of the sample to begin returning stats. Default is to read
            stats after burn in. To start from the beginning set thin_start
            to 0.
        thin_interval : int
            Interval to accept every i-th sample. Default is to use the
            `fp.acl`. If `fp.acl` is not set, then use all stats
            (set thin_interval to 1).
        thin_end : int
            Index of the last sample to read. If not given then
            `fp.niterations` is used.

        Returns
        -------
        lnZ : float
            The estimate of log of the evidence.
        dlnZ : float
            The error on the estimate.
        """
        with cls._io(filename, 'r') as fp:
            logls = fp.read_raw_samples(['loglikelihood'],
                                        thin_start=thin_start,
                                        thin_interval=thin_interval,
                                        thin_end=thin_end,
                                        temps='all', flatten=False)
            logls = logls['loglikelihood']
            # we need the betas that were used
            betas = fp.read_betas(thin_start=thin_start,
                                  thin_interval=thin_interval,
                                  thin_end=thin_end)
            # we'll separate betas out by their unique temperatures
            # there's probably a faster way to do this...
            mean_logls = []
            unique_betas = []
            ntemps = betas.shape[0]
            for ti in range(ntemps):
                ubti, idx = numpy.unique(betas[ti,:], return_inverse=True)
                unique_idx = numpy.unique(idx)
                loglsti = logls[ti, :, :]
                for ii in unique_idx:
                    # average over the walkers and iterations with the same
                    # betas
                    getiters = numpy.where(ii == unique_idx)[0]
                    mean_logls.append(loglsti[:, getiters].mean())
                    unique_betas.append(ubti[ii])
        return ptemcee.util.thermodynamic_integration_log_evidence(
            numpy.array(unique_betas), numpy.array(mean_logls))

    @classmethod
    def from_config(cls, cp, model, nprocesses=1, use_mpi=False):
        """Loads the sampler from the given config file.

        The following arguments are searched for, and passed to the specified
        argument when initializing the class:

         * ``nwalkers = INT``: used for the ``nwalkers`` argument
         * ``checkpoint-interval = INT``: used for ``checkpoint_interval``
         * ``logl-function = STR``: used for ``loglikelihood_function``
         * ``ntemps = INT``: used for ``ntemps``
         * ``Tmax = FLOAT``: used for ``Tmax``
         * ``betas = FLOAT1, FLOAT2, [...]``: used for the ``betas``
         * ``no-adaptation =``: sets the ``adaptive`` argument to False
           if listed
         * ``adaptation-time = INT``: sets the ``adaptation_time``
         * ``adaptation-lag = INT``: sets the ``adaptation_lag``
         * ``scale-factor = FLOAT``: sets the ``scale_factor``

        Parameters
        ----------
        cp : ConfigParser
            ConfigParser instance to read.
        model : Model instance
            Model instance to use.
        nprocesses : int, optional
            The number of processors to use. Default is 1.
        use_mpi : bool, optional
            Use MPI pool for parallelization. Default is False.
        """
        section = "sampler"
        # check name
        assert cp.get(section, "name") == cls.name, (
            "name in section [sampler] must match mine")
        # get the number of walkers to use
        nwalkers = int(cp.get(section, "nwalkers"))
        # get the checkpoint interval, if it's specified
        checkpoint_interval = cls.checkpoint_from_config(cp, section)
        # get the loglikelihood function
        logl = get_optional_arg_from_config(cp, section, 'logl-function')
        # get optional args
        optargs = {}
        ntemps = get_optional_arg_from_config(cp, section, 'ntemps', int)
        if ntemps is not None:
            optargs['ntemps'] = ntemps
        tmax = get_optional_arg_from_config(cp, section, 'Tmax', float)
        if tmax is not None:
            optargs['Tmax'] = tmax
        betas = get_optional_arg_from_config(cp, section, 'betas')
        if betas is not None:
            # convert to list
            optargs['betas'] = map(float, betas.split(','))
        optargs['adaptive'] = not cp.has_option(section, 'no-adaptation')
        adaptation_lag = get_optional_arg_from_config(cp, section,
                                                      'adaptation-lag', int)
        if adaptation_lag is not None:
            optargs['adaptation_lag'] = adaptation_lag
        adaptation_time = get_optional_arg_from_config(cp, section,
                                                       'adaptation-time', int)
        if adaptation_time is not None:
            optargs['adaptation_time'] = adaptation_time
        scale_factor = get_optional_arg_from_config(cp, section,
                                                    'scale-factor', float)
        if scale_factor is not None:
            optargs['scale_factor'] = scale_factor
        obj = cls(model, nwalkers, loglikelihood_function=logl,
                  checkpoint_interval=checkpoint_interval,
                  nprocesses=nprocesses, use_mpi=use_mpi,
                  **optargs)
        # set target
        obj.set_target_from_config(cp, section)
        # add burn-in if it's specified
        obj.set_burn_in_from_config(cp)
        return obj

