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


import numpy
import emcee
import logging
from pycbc.pool import choose_pool

from .base import (BaseSampler, setup_output)
from .base_mcmc import (BaseMCMC, EnsembleSupport, raw_samples_to_dict,
                        get_optional_arg_from_config)
from .base_multitemper import (MultiTemperedSupport,
                               ensemble_compute_acf, ensemble_compute_acl)
from ..burn_in import EnsembleMultiTemperedMCMCBurnInTests
from pycbc.inference.io import EmceePTFile
from .. import models


if emcee.__version__ >= '3.0.0':
    raise ImportError


class EmceePTSampler(MultiTemperedSupport, EnsembleSupport, BaseMCMC,
                     BaseSampler):
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
        temperature ladder. If not provided, ``emcee_pt`` will use the number
        of temperatures and the number of dimensions of the parameter space to
        construct the ladder with geometrically spaced temperatures.
    loglikelihood_function : str, optional
        Set the function to call from the model for the ``loglikelihood``.
        Default is ``loglikelihood``.
    nprocesses : int, optional
        The number of parallel processes to use. Default is 1
        (no paralleliztion).
    use_mpi : bool, optional
        Use MPI for parallelization. Default (False) will use python's
        multiprocessing.
    """
    name = "emcee_pt"
    _io = EmceePTFile
    burn_in_class = EnsembleMultiTemperedMCMCBurnInTests

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

        # these are used to help paralleize over multiple cores / MPI
        models._global_instance = model_call
        model_call = models._call_global_model
        prior_call = models._call_global_model_logprior
        self.pool = choose_pool(mpi=use_mpi, processes=nprocesses)

        # construct the sampler: PTSampler needs the likelihood and prior
        # functions separately
        ndim = len(model.variable_params)
        self._sampler = emcee.PTSampler(ntemps, nwalkers, ndim,
                                        model_call, prior_call, pool=self.pool,
                                        betas=betas)
        self.nwalkers = nwalkers
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

    @staticmethod
    def compute_acf(filename, **kwargs):
        r"""Computes the autocorrelation function.

        Calls :py:func:`base_multitemper.ensemble_compute_acf`; see that
        function for details.

        Parameters
        ----------
        filename : str
            Name of a samples file to compute ACFs for.
        \**kwargs :
            All other keyword arguments are passed to
            :py:func:`base_multitemper.ensemble_compute_acf`.

        Returns
        -------
        dict :
            Dictionary of arrays giving the ACFs for each parameter. If
            ``per-walker=True`` is passed as a keyword argument, the arrays
            will have shape ``ntemps x nwalkers x niterations``. Otherwise, the
            returned array will have shape ``ntemps x niterations``.
        """
        return ensemble_compute_acf(filename, **kwargs)

    @staticmethod
    def compute_acl(filename, **kwargs):
        r"""Computes the autocorrelation length.

        Calls :py:func:`base_multitemper.ensemble_compute_acl`; see that
        function for details.

        Parameters
        -----------
        filename : str
            Name of a samples file to compute ACLs for.
        \**kwargs :
            All other keyword arguments are passed to
            :py:func:`base_multitemper.ensemble_compute_acl`.

        Returns
        -------
        dict
            A dictionary of ntemps-long arrays of the ACLs of each parameter.
        """
        return ensemble_compute_acl(filename, **kwargs)

    @classmethod
    def from_config(cls, cp, model, output_file=None, nprocesses=1,
                    use_mpi=False):
        """Loads the sampler from the given config file.

        The following options are retrieved in the ``[sampler]`` section:

        * ``name`` :
            Required. This must match the samlper's name.
        * ``nwalkers`` :
            Required. The number of walkers to use.
        * ``ntemps`` :
            The number of temperatures to use. Either this, or
            ``inverse-temperatures-file`` must be provided (but not both).
        * ``inverse-temperatures-file`` :
            Path to an hdf file containing the inverse temperatures ("betas")
            to use. The betas will be retrieved from the file's
            ``.attrs['betas']``. Either this or ``ntemps`` must be provided
            (but not both).
        * ``niterations`` :
            The number of iterations to run the sampler for. Either this or
            ``effective-nsamples`` must be provided (but not both).
        * ``effective-nsamples`` :
            Run the sampler until the given number of effective samples are
            obtained. A ``checkpoint-interval`` must also be provided in this
            case. Either this or ``niterations`` must be provided (but not
            both).
        * ``thin-interval`` :
            Thin the samples by the given value before saving to disk. May
            provide this, or ``max-samples-per-chain``, but not both. If
            neither options are provided, will save all samples.
        * ``max-samples-per-chain`` :
            Thin the samples such that the number of samples per chain per
            temperature that are saved to disk never exceeds the given value.
            May provide this, or ``thin-interval``, but not both. If neither
            options are provided, will save all samples.
        * ``checkpoint-interval`` :
            Sets the checkpoint interval to use. Must be provided if using
            ``effective-nsamples``.
        * ``checkpoint-signal`` :
            Set the checkpoint signal, e.g., "USR2". Optional.
        * ``logl-function`` :
            The attribute of the model to use for the loglikelihood. If
            not provided, will default to ``loglikelihood``.

        Settings for burn-in tests are read from ``[sampler-burn_in]``. In
        particular, the ``burn-in-test`` option is used to set the burn in
        tests to perform. See
        :py:func:`MultiTemperedMCMCBurnInTests.from_config` for details. If no
        ``burn-in-test`` is provided, no burn in tests will be carried out.

        Parameters
        ----------
        cp : WorkflowConfigParser instance
            Config file object to parse.
        model : pycbc.inference.model.BaseModel instance
            The model to use.
        output_file : str, optional
            The name of the output file to checkpoint and write results to.
        nprocesses : int, optional
            The number of parallel processes to use. Default is 1.
        use_mpi : bool, optional
            Use MPI for parallelization. Default is False.

        Returns
        -------
        EmceePTSampler :
            The sampler instance.
        """
        section = "sampler"
        # check name
        assert cp.get(section, "name") == cls.name, (
            "name in section [sampler] must match mine")
        # get the number of walkers to use
        nwalkers = int(cp.get(section, "nwalkers"))
        # get the temps/betas
        ntemps, betas = cls.betas_from_config(cp, section)
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
        # Set up the output file
        setup_output(obj, output_file)
        if not obj.new_checkpoint:
            obj.resume_from_checkpoint()
        else:
            obj.set_start_from_config(cp)
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
            fp.write_samples(self.samples,
                             parameters=self.model.variable_params,
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

    def _correctjacobian(self, samples):
        """Corrects the log jacobian values stored on disk.

        Parameters
        ----------
        samples : dict
            Dictionary of the samples.
        """
        # flatten samples for evaluating
        orig_shape = list(samples.values())[0].shape
        flattened_samples = {p: arr.ravel()
                             for p, arr in list(samples.items())}
        # convert to a list of tuples so we can use map function
        params = list(flattened_samples.keys())
        size = flattened_samples[params[0]].size
        logj = numpy.zeros(size)
        for ii in range(size):
            these_samples = {p: flattened_samples[p][ii] for p in params}
            these_samples = self.model.sampling_transforms.apply(these_samples)
            self.model.update(**these_samples)
            logj[ii] = self.model.logjacobian
        return logj.reshape(orig_shape)

    def finalize(self):
        """Calculates the log evidence and writes to the checkpoint file.

        If sampling transforms were used, this also corrects the jacobian
        stored on disk.

        The thin start/interval/end for calculating the log evidence are
        retrieved from the checkpoint file's thinning attributes.
        """
        if self.model.sampling_transforms is not None:
            # fix the lobjacobian values stored on disk
            logging.info("Correcting logjacobian values on disk")
            with self.io(self.checkpoint_file, 'r') as fp:
                samples = fp.read_raw_samples(self.variable_params,
                                              thin_start=0,
                                              thin_interval=1, thin_end=None,
                                              temps='all', flatten=False)
            logjacobian = self._correctjacobian(samples)
            # write them back out
            for fn in [self.checkpoint_file, self.backup_file]:
                with self.io(fn, "a") as fp:
                    fp[fp.samples_group]['logjacobian'][()] = logjacobian
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
