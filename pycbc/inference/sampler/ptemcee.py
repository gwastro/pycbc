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


import shlex
import numpy
import ptemcee
import logging
from pycbc.pool import choose_pool

from .base import (BaseSampler, setup_output)
from .base_mcmc import (BaseMCMC, EnsembleSupport, raw_samples_to_dict,
                        get_optional_arg_from_config)
from .base_multitemper import (read_betas_from_hdf,
                               ensemble_compute_acf, ensemble_compute_acl)
from ..burn_in import EnsembleMultiTemperedMCMCBurnInTests
from pycbc.inference.io import PTEmceeFile
from .. import models


class PTEmceeSampler(EnsembleSupport, BaseMCMC, BaseSampler):
    """This class is used to construct the parallel-tempered ptemcee sampler.

    Parameters
    ----------
    model : model
        A model from ``pycbc.inference.models``.
    nwalkers : int
        Number of walkers to use in sampler.
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
        Whether or not to use adaptive temperature levels. Default is False.
    adaptation_lag : int, optional
        Only used if ``adaptive`` is True; see :py:mod:`ptemcee.Sampler` for
        details. If not provided, will use ``ptemcee``'s default.
    adaptation_time : int, optional
        Only used if ``adaptive`` is True; see :py:mod:`ptemcee.Sampler` for
        details. If not provided, will use ``ptemcee``'s default.
    scale_factor : float, optional
        Scale factor used for the stretch proposal; see
        :py:mod:`ptemcee.Sampler` for details. If not provided, will use
        ``ptemcee``'s default.
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
    name = "ptemcee"
    _io = PTEmceeFile
    burn_in_class = EnsembleMultiTemperedMCMCBurnInTests

    def __init__(self, model, nwalkers, ntemps=None, Tmax=None, betas=None,
                 adaptive=False, adaptation_lag=None, adaptation_time=None,
                 scale_factor=None,
                 loglikelihood_function=None,
                 checkpoint_interval=None, checkpoint_signal=None,
                 nprocesses=1, use_mpi=False):

        self.model = model
        ndim = len(model.variable_params)
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
        if adaptation_lag is not None:
            kwargs['adaptation_lag'] = adaptation_lag
        if adaptation_time is not None:
            kwargs['adaptation_time'] = adaptation_time
        if scale_factor is not None:
            kwargs['scale_factor'] = scale_factor
        # create a wrapper for calling the model
        if loglikelihood_function is None:
            loglikelihood_function = 'loglikelihood'
        # frustratingly, ptemcee does not support blob data, so we have to
        # turn it off
        model_call = models.CallModel(model, loglikelihood_function,
                                      return_all_stats=False)
        # these are used to help paralleize over multiple cores / MPI
        models._global_instance = model_call
        model_call = models._call_global_model
        prior_call = models._call_global_model_logprior
        self.pool = choose_pool(mpi=use_mpi, processes=nprocesses)
        # construct the sampler
        self._sampler = ptemcee.Sampler(nwalkers=nwalkers, ndim=ndim,
                                        logl=model_call, logp=prior_call,
                                        mapper=self.pool.map, **kwargs)
        self.nwalkers = nwalkers
        self._ntemps = ntemps
        self._checkpoint_interval = checkpoint_interval
        self._checkpoint_signal = checkpoint_signal
        # we'll initialize ensemble and chain to None
        self._chain = None
        self._ensemble = None

    @property
    def io(self):
        return self._io

    @property
    def ntemps(self):
        """The number of temeratures that are set."""
        return self._ntemps

    @property
    def base_shape(self):
        return (self.ntemps, self.nwalkers,)

    @property
    def betas(self):
        """Returns the beta history currently in memory."""
        # chain betas has shape niterations x ntemps; transpose to
        # ntemps x niterations
        return self._chain.betas.transpose()

    @property
    def starting_betas(self):
        """Returns the betas that were used at startup."""
        # the initial betas that were used
        return self._sampler.betas

    @property
    def adaptive(self):
        """Whether or not the betas are adapted."""
        return self._sampler.adaptive

    @property
    def adaptation_lag(self):
        """The adaptation lag for the beta evolution."""
        return self._sampler.adaptation_lag

    @property
    def adaptation_time(self):
        """The adaptation time for the beta evolution."""
        return self._sampler.adaptation_time

    @property
    def scale_factor(self):
        """The scale factor used by ptemcee."""
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
        will only return the loglikelihood and logprior (with the logjacobian
        set to zero) regardless of what stats the model can return.


        .. warning::
            Since the ``logjacobian`` is not saved by ``ptemcee``, the
            ``logprior`` returned here is the log of the prior pdf in the
            sampling coordinate frame rather than the variable params frame.
            This differs from the variable params frame by the log of the
            Jacobian of the transform from one frame to the other. If no
            sampling transforms were used, then the ``logprior`` is the same.
        """
        # log likelihood and prior have shape
        # niterations x ntemps x nwalkers; we'll tranpose to have shape
        # ntemps x nwalkers x niterations
        logl = self._chain.logl.transpose((1, 2, 0))
        logp = self._chain.logP.transpose((1, 2, 0))
        logjacobian = numpy.zeros(logp.shape)
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

    @classmethod
    def calculate_logevidence(cls, filename, thin_start=None, thin_end=None,
                              thin_interval=None):
        """Calculates the log evidence from the given file.
        This uses ``ptemcee``'s thermodynamic integration.

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
                ubti, idx = numpy.unique(betas[ti, :], return_inverse=True)
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
        ----------
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

        * ``name = STR`` :
            Required. This must match the sampler's name.
        * ``nwalkers = INT`` :
            Required. The number of walkers to use.
        * ``ntemps = INT`` :
            The number of temperatures to use. This may be used in combination
            with ``Tmax``. Either this, ``Tmax``, ``betas`` or ``betas-file``
            must be provided.
        * ``tmax = FLOAT`` :
            The maximum temperature to use. This may be used in combination
            with ``ntemps``, or alone.
        * ``betas = FLOAT1 FLOAT2 [...]`` :
            Space-separated list of (intial) inverse temperatures ("betas") to
            use. This sets both the number of temperatures and the tmax. A
            ``ValueError`` will be raised if both this and ``ntemps`` or
            ``Tmax`` are provided.
        * ``betas-file = STR`` :
            Path to an hdf file containing the inverse temperatures ("betas")
            to use. The betas will be retrieved from the file's
            ``.attrs['betas']``. A ``ValueError`` will be raised if both this
            and ``betas`` are provided.
        * ``adaptive =`` :
            If provided, temperature adaptation will be turned on.
        * ``adaptation-lag = INT`` :
            The adaptation lag to use (see ptemcee for details).
        * ``adaptation-time = INT`` :
            The adaptation time to use (see ptemcee for details).
        * ``scale-factor = FLOAT`` :
            The scale factor to use for the emcee stretch.
        * ``niterations = INT`` :
            The number of iterations to run the sampler for. Either this or
            ``effective-nsamples`` must be provided (but not both).
        * ``effective-nsamples = INT`` :
            Run the sampler until the given number of effective samples are
            obtained. A ``checkpoint-interval`` must also be provided in this
            case. Either this or ``niterations`` must be provided (but not
            both).
        * ``thin-interval = INT`` :
            Thin the samples by the given value before saving to disk. May
            provide this, or ``max-samples-per-chain``, but not both. If
            neither options are provided, will save all samples.
        * ``max-samples-per-chain = INT`` :
            Thin the samples such that the number of samples per chain per
            temperature that are saved to disk never exceeds the given value.
            May provide this, or ``thin-interval``, but not both. If neither
            options are provided, will save all samples.
        * ``checkpoint-interval = INT`` :
            Sets the checkpoint interval to use. Must be provided if using
            ``effective-nsamples``.
        * ``checkpoint-signal = STR`` :
            Set the checkpoint signal, e.g., "USR2". Optional.
        * ``logl-function = STR`` :
            The attribute of the model to use for the loglikelihood. If
            not provided, will default to ``loglikelihood``.

        Settings for burn-in tests are read from ``[sampler-burn_in]``. In
        particular, the ``burn-in-test`` option is used to set the burn in
        tests to perform. See
        :py:func:`EnsembleMultiTemperedMCMCBurnInTests.from_config` for
        details. If no ``burn-in-test`` is provided, no burn in tests will be
        carried out.

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
        # get the checkpoint interval, if it's specified
        checkpoint_interval = cls.checkpoint_from_config(cp, section)
        checkpoint_signal = cls.ckpt_signal_from_config(cp, section)
        optargs = {}
        # get the temperature level settings
        ntemps = get_optional_arg_from_config(cp, section, 'ntemps', int)
        if ntemps is not None:
            optargs['ntemps'] = ntemps
        tmax = get_optional_arg_from_config(cp, section, 'tmax', float)
        if tmax is not None:
            optargs['Tmax'] = tmax
        betas = get_optional_arg_from_config(cp, section, 'betas')
        if betas is not None:
            # convert to list sorted in descencding order
            betas = numpy.sort(list(map(float, shlex.split(betas))))[::-1]
            optargs['betas'] = betas
        betas_file = get_optional_arg_from_config(cp, section, 'betas-file')
        if betas_file is not None:
            optargs['betas'] = read_betas_from_hdf(betas_file)
        # check for consistency
        if betas is not None and betas_file is not None:
            raise ValueError("provide either betas or betas-file, not both")
        if 'betas' in optargs and (ntemps is not None or tmax is not None):
            raise ValueError("provide either ntemps/tmax or betas/betas-file, "
                             "not both")
        # adaptation parameters
        adaptive = get_optional_arg_from_config(cp, section, 'adaptive')
        if adaptive is not None:
            optargs['adaptive'] = True
        else:
            optargs['adaptive'] = False
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
        # get the loglikelihood function
        logl = get_optional_arg_from_config(cp, section, 'logl-function')
        obj = cls(model, nwalkers,
                  checkpoint_interval=checkpoint_interval,
                  checkpoint_signal=checkpoint_signal,
                  loglikelihood_function=logl, nprocesses=nprocesses,
                  use_mpi=use_mpi, **optargs)
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

    def write_results(self, filename):
        """Writes samples, model stats, acceptance fraction, and random state
        to the given file.

        Parameters
        ----------
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
            # write random state
            fp.write_random_state()
            # write betas
            fp.write_betas(self.betas, last_iteration=self.niterations)
            # write random state
            fp.write_random_state()
            # write attributes of the ensemble
            fp.write_ensemble_attrs(self.ensemble)

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
