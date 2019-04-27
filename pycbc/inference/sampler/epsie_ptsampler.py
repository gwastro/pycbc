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

"""This module provides classes for interacting with epsie samplers.
"""

from __future__ import absolute_import

import epsie

from .base_mcmc import (BaseMCMC, raw_samples_to_dict,
                        get_optional_arg_from_config)
from .base_multitemper import (MultiTemperedSupport,
                               MultiTemperedAutocorrSupport)
from ..burn_in import MultiTemperedMCMCBurnInTests


class EpsiePTSampler(MultiTemperedAutocorrSupport, MultiTemperedSupport,
                     BaseMCMC, BaseSampler):
    """Constructs an MCMC sampler using epsie's parallel-tempered sampler.

    Parameters
    ----------
    model : model
        A model from ``pycbc.inference.models``.
    nwalkers : int
        Number of walkers to use in the sampler.
    ntemps : int, optional
        Number of temperatures to use in the sampler. A geometrically-spaced
        temperature ladder with the gievn number of levels will be constructed
        based on the number of parameters. If not provided, must provide
        ``betas``.
    betas : array, optional
        An array of inverse temperature values to be used in for the
        temperature ladder. If not provided, must provide ``ntemps``.
    proposals : dict, optional
        Dictionary mapping sampling parameter names to proposal classes. Any
        parameters not listed will use the ``default_proposal``. **Note:**
        proposals should be specified for the sampling parameters, not the
        variable parameters.
    default_proposal : an epsie.Proposal class, optional
        The default proposal to use for parameters not in ``proposals``.
        Default is :py:class:`epsie.proposals.Normal`.
    default_proposal_args : dict, optional
        Dictionary of arguments to pass to the default proposal.
    seed : int, optional
        Seed for epsie's random number generator. If None provided, will create
        one.
    checkpoint_interval : int, optional
        Specify the number of iterations to do between checkpoints. If not
        provided, no checkpointin will be done.
    checkpoint_signal : str, optional
        Set the signal to use when checkpointing. For example, 'USR2'.
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
    name = "epsie_ptsampler"
    _io = EpsiePTFile
    burn_in_class = MultiTemperedMCMCBurnInTests

    def __init__(self, model, nchains, ntemps=None, betas=None,
                 proposals=None, default_proposal=None,
                 default_proposal_args=None, seed=None,
                 checkpoint_interval=None, checkpoint_signal=None,
                 loglikelihood_function=None,
                 nprocesses=1, use_mpi=False):

        self.model = model
        # create a wrapper for calling the model
        model_call = _EpsieCallModel(model, loglikelihood_function)
        # Set up the pool
        if nprocesses > 1:
            # these are used to help paralleize over multiple cores / MPI
            models._global_instance = model_call
        pool = choose_pool(mpi=use_mpi, processes=nprocesses)
        if pool is not None:
            pool.count = nprocesses
        # initialize the sampler
        self._sampler = epsie.ParallelTemperedSampler(
            model.sampling_params, model, nwalkers, betas=betas,
            proposals=proposals, default_proposal=default_proposal,
            default_proposal_args=default_proposal_args,
            seed=seed, pool=pool) 
        # set other parameters
        self._nwalkers = nchains
        self._ntemps = ntemps
        self._checkpoint_interval = checkpoint_interval
        self._checkpoint_signal = checkpoint_signal

    @property
    def io(self):
        return self._io

    @property
    def base_shape(self):
        return (self.ntemps, self.nchains,)

    @property
    def nchains(self):
        """Alias for ``nwalkers``."""
        return self._nwalkers

    @property
    def betas(self):
        return self._sampler.betas

    @classmethod
    def from_config(cls, cp, model, nprocesses=1, use_mpi=False):
        """Loads the sampler from the given config file.

        The following options are retrieved in the ``[sampler]`` section:

        * ``name`` :
            (required) must match the samlper's name
        * ``nchains`` :
            (required) the number of chains to use
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
        * ``seed`` :
            The seed to use for epsie's random number generator. If not
            provided, epsie will create one.
        * ``logl-function`` :
            The attribute of the model to use for the loglikelihood. If
            not provided, will default to ``loglikelihood``.

        In addition, the following sections are read:

        * ``[sampler-burn_in]``
            Sets the burn in test and options to do. In particular, the
            ``burn-in-test`` option is used to set the burn in tests to
            perform. See :py:func:`MultiTemperedMCMCBurnInTests.from_config`
            for details. If no ``burn-in-test`` is provided, no burn in tests
            will be carried out.
        * ``[jump_proposal-{params}]``
            Sets the jump proposals and arguments to use for specific
            sampling parameters. See :py:func:`sampler.proposals.from_config`
            for details.
        * ``[jump_proposal-default]``
            Sets the default jump proposal and arguments to use for all other
            parameters. See :py:func:`sampler.proposals.from_config` for
            details. If no default is provided, will use
            :py:class:`epsie.ParallelTemperedSampler`'s default jump proposal.

        .. note::
            Jump proposals should be specified for **sampling parameters**,
            not **variable parameters**.

        Parameters
        ----------
        cp : WorkflowConfigParser instance
            Config file object to parse.
        model : pycbc.inference.model.BaseModel instance
            The model to use.
        nprocesses : int, optional
            The number of parallel processes to use. Default is 1.
        use_mpi : bool, optional
            Use MPI for parallelization. Default is False.

        Returns
        -------
        EpsiePTSampler :
            The sampler instance.
        """
        section = "sampler"
        # check name
        assert cp.get(section, "name") == cls.name, (
            "name in section [sampler] must match mine")
        nchains = int(cp.get(section, "nchains"))
        seed = int(cp.get(section, "seed"))
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
        # initialize
        obj = cls(model, nchains, ntemps=ntemps, betas=betas,
                  proposals=None, default_proposal=None,
                  default_proposal_args=None, seed=seed,
                  checkpoint_interval=checkpoint_interval,
                  checkpoint_signal=checkpoint_signal,
                  loglikelihood_function=logl,
                  nprocesses=nprocesses, use_mpi=use_mpi)
        # set target
        obj.set_target_from_config(cp, section)
        # add burn-in if it's specified
        obj.set_burn_in_from_config(cp)
        # set prethin options
        obj.set_thin_interval_from_config(cp, section)
        return obj


class _EpsieCallModel(object):
    """Model wrapper for epsie.
    
    Allows model to be called like a function. Returns the loglikelihood
    function, logprior, and the model's default stats.
    """

    def __init__(self, model, loglikelihood_function=None):
        self.model = model
        if loglikelihood_function is None:
            loglikelihood_function = 'loglikelihood'
        self.loglikelihood_function = loglikelihood_function

    def __call__(self, **kwargs):
        """Calls update, then calls the loglikelihood and logprior."""
        self.model.update(**kwargs)
        logl = getattr(model, self.loglikelhood_function)
        logp = model.logprior
        return logl, logp, self.model.get_current_stats()
