# Copyright (C) 2016  Collin Capano
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
This modules provides classes and functions for using the emcee_pt sampler
packages for parameter estimation.
"""

from __future__ import absolute_import

import numpy
import emcee
import h5py
import logging
from pycbc.pool import choose_pool

from .base import BaseSampler
from .base_mcmc import (BaseMCMC, raw_samples_to_dict,
                        get_optional_arg_from_config)
from .base_multitemper import (MultiTemperedSupport,
                               MultiTemperedAutocorrSupport)
from ..burn_in import MultiTemperedMCMCBurnInTests
from pycbc.inference.io import EmceePTFile
from .. import models


class EmceePTSampler(MultiTemperedAutocorrSupport, MultiTemperedSupport,
                     BaseMCMC, BaseSampler):
    """This class is used to construct a parallel-tempered MCMC sampler from
    the emcee package's PTSampler.

    Parameters
    ----------
    model : model
        A model from ``pycbc.inference.models``.
    ntemps : int
        Number of temeratures to use in the sampler.
    nwalkers : int
        Number of walkers to use in sampler.
    betas : array
        An array of inverse temperature values to be used in emcee_pt's
        temperature ladder. If not provided, emcee_pt will use the number of
        temperatures and the number of dimensions of the parameter space to
        construct the ladder with geometrically spaced temperatures.
    pool : function with map, Optional
        A provider of a map function that allows a function call to be run
        over multiple sets of arguments and possibly maps them to
        cores/nodes/etc.
    """
    name = "emcee_pt"
    _io = EmceePTFile
    burn_in_class = MultiTemperedMCMCBurnInTests

    def __init__(self, model, ntemps, nwalkers, betas=None,
                 checkpoint_interval=None, checkpoint_signal=None,
                 loglikelihood_function=None,
                 nprocesses=1, use_mpi=False):

        self.model = model

        # create a wrapper for calling the model
        if loglikelihood_function is None:
            loglikelihood_function = 'loglikelihood'
        # frustratingly, emcee_pt does not support blob data, so we have to
        # turn it off
        model_call = models.CallModel(model, loglikelihood_function,
                                      return_all_stats=False)

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

        # construct the sampler: PTSampler needs the likelihood and prior
        # functions separately
        ndim = len(model.variable_params)
        self._sampler = emcee.PTSampler(ntemps, nwalkers, ndim,
                                        model_call, prior_call, pool=pool,
                                        betas=betas)
        self._nwalkers = nwalkers
        self._ntemps = ntemps
        self._checkpoint_interval = checkpoint_interval
        self._checkpoint_signal = checkpoint_signal

    @property
    def io(self):
        return self._io

    @property
    def base_shape(self):
        return (self.ntemps, self.nwalkers,)

    @property
    def betas(self):
        return self._sampler.betas

    @classmethod
    def from_config(cls, cp, model, nprocesses=1, use_mpi=False):
        """
        Loads the sampler from the given config file.

        For generating the temperature ladder to be used by emcee_pt, either
        the number of temperatures (provided by the option 'ntemps'),
        or the path to a file storing inverse temperature values (provided
        under a subsection inverse-temperatures-file) can be loaded from the
        config file. If the latter, the file should be of hdf format, having
        an attribute named 'betas' storing the list of inverse temperature
        values to be provided to emcee_pt. If the former, emcee_pt will
        construct the ladder with "ntemps" geometrically spaced temperatures.
        """
        section = "sampler"
        # check name
        assert cp.get(section, "name") == cls.name, (
            "name in section [sampler] must match mine")
        # get the number of walkers to use
        nwalkers = int(cp.get(section, "nwalkers"))
        if cp.has_option(section, "ntemps") and \
                cp.has_option(section, "inverse-temperatures-file"):
            raise ValueError("Must specify either ntemps or "
                             "inverse-temperatures-file, not both.")
        if cp.has_option(section, "inverse-temperatures-file"):
            # get the path of the file containing inverse temperatures values.
            inverse_temperatures_file = cp.get(section,
                                               "inverse-temperatures-file")
            with h5py.File(inverse_temperatures_file, "r") as fp:
                try:
                    betas = numpy.array(fp.attrs['betas'])
                    ntemps = betas.shape[0]
                except KeyError:
                    raise AttributeError("No attribute called betas")
        else:
            # get the number of temperatures
            betas = None
            ntemps = int(cp.get(section, "ntemps"))
        # get the checkpoint interval, if it's specified
        checkpoint_interval = cls.checkpoint_from_config(cp, section)
        checkpoint_signal = cls.ckpt_signal_from_config(cp, section)
        # get the loglikelihood function
        logl = get_optional_arg_from_config(cp, section, 'logl-function')
        obj = cls(model, ntemps, nwalkers, betas=betas,
                  checkpoint_interval=checkpoint_interval,
                  checkpoint_signal=checkpoint_signal,
                  loglikelihood_function=logl, nprocesses=nprocesses,
                  use_mpi=use_mpi)
        # set target
        obj.set_target_from_config(cp, section)
        # add burn-in if it's specified
        obj.set_burn_in_from_config(cp)
        # set prethin options
        obj.set_thin_interval_from_config(cp, section)
        return obj

    @property
    def samples(self):
        """A dict mapping ``variable_params`` to arrays of samples currently
        in memory.

        The arrays have shape ``ntemps x nwalkers x niterations``.
        """
        # emcee stores samples to it's chain attribute as a
        # nwalker x niterations x ndim array
        raw_samples = self._sampler.chain
        return raw_samples_to_dict(self, raw_samples)

    @property
    def model_stats(self):
        """Returns the log likelihood ratio and log prior as a dict of arrays.

        The returned array has shape ntemps x nwalkers x niterations.

        Unfortunately, because ``emcee_pt`` does not have blob support, this
        will only return the loglikelihood and logprior (with the logjacobian
        set to zero) regardless of what stats the model can return.


        .. warning::
            Since the `logjacobian` is not saved by `emcee_pt`, the `logprior`
            returned here is the log of the prior pdf in the sampling
            coordinate frame rather than the variable params frame. This
            differs from the variable params frame by the log of the Jacobian
            of the transform from one frame to the other. If no sampling
            transforms were used, then the `logprior` is the same.
        """
        # likelihood has shape ntemps x nwalkers x niterations
        logl = self._sampler.lnlikelihood
        # get prior from posterior
        logp = self._sampler.lnprobability - logl
        logjacobian = numpy.zeros(logp.shape)
        return {'loglikelihood': logl, 'logprior': logp,
                'logjacobian': logjacobian}

    def clear_samples(self):
        """Clears the chain and blobs from memory.
        """
        # store the iteration that the clear is occuring on
        self._lastclear = self.niterations
        self._itercounter = 0
        # now clear the chain
        self._sampler.reset()

    def set_state_from_file(self, filename):
        """Sets the state of the sampler back to the instance saved in a file.
        """
        with self.io(filename, 'r') as fp:
            rstate = fp.read_random_state()
        # set the numpy random state
        numpy.random.set_state(rstate)

    def run_mcmc(self, niterations):
        """Advance the ensemble for a number of samples.

        Parameters
        ----------
        niterations : int
            Number of samples to get from sampler.
        """
        pos = self._pos
        if pos is None:
            pos = self._p0
        res = self._sampler.run_mcmc(pos, niterations)
        p, _, _ = res[0], res[1], res[2]
        # update the positions
        self._pos = p

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
            fp.write_samples(self.samples, self.model.variable_params,
                             last_iteration=self.niterations)
            # write stats
            fp.write_samples(self.model_stats, last_iteration=self.niterations)
            # write accpetance
            fp.write_acceptance_fraction(self._sampler.acceptance_fraction)
            # write random state
            fp.write_random_state()

    @classmethod
    def calculate_logevidence(cls, filename, thin_start=None, thin_end=None,
                              thin_interval=None):
        """Calculates the log evidence from the given file using ``emcee_pt``'s
        thermodynamic integration.

        Parameters
        ----------
        filename : str
            Name of the file to read the samples from. Should be an
            ``EmceePTFile``.
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
            betas = fp.betas
            # annoyingly, theromdynaimc integration in PTSampler is an instance
            # method, so we'll implement a dummy one
            ntemps = fp.ntemps
            nwalkers = fp.nwalkers
            ndim = len(fp.variable_params)
        dummy_sampler = emcee.PTSampler(ntemps, nwalkers, ndim, None,
                                        None, betas=betas)
        return dummy_sampler.thermodynamic_integration_log_evidence(
            logls=logls, fburnin=0.)

    def finalize(self):
        """Calculates the log evidence and writes to the checkpoint file.

        The thin start/interval/end for calculating the log evidence are
        retrieved from the checkpoint file's thinning attributes.
        """
        logging.info("Calculating log evidence")
        # get the thinning settings
        with self.io(self.checkpoint_file, 'r') as fp:
            thin_start = fp.thin_start
            thin_interval = fp.thin_interval
            thin_end = fp.thin_end
        # calculate
        logz, dlogz = self.calculate_logevidence(
            self.checkpoint_file, thin_start=thin_start, thin_end=thin_end,
            thin_interval=thin_interval)
        logging.info("log Z, dlog Z: {}, {}".format(logz, dlogz))
        # write to both the checkpoint and backup
        for fn in [self.checkpoint_file, self.backup_file]:
            with self.io(fn, "a") as fp:
                fp.write_logevidence(logz, dlogz)
