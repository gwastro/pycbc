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


import numpy

import epsie
from epsie.samplers import ParallelTemperedSampler

# we'll use emcee_pt's default beta ladder for temperature levels
from emcee.ptsampler import default_beta_ladder

from pycbc.pool import choose_pool

from .base import (BaseSampler, setup_output)
from .base_mcmc import (BaseMCMC, get_optional_arg_from_config,
                        nsamples_in_chain)
from .base_multitemper import (MultiTemperedSupport, compute_acf, compute_acl,
                               acl_from_raw_acls)
from ..burn_in import MultiTemperedMCMCBurnInTests
from ..jump import epsie_proposals_from_config
from ..io import EpsieFile
from .. import models


class EpsieSampler(MultiTemperedSupport, BaseMCMC, BaseSampler):
    """Constructs an MCMC sampler using epsie's parallel-tempered sampler.

    Parameters
    ----------
    model : model
        A model from ``pycbc.inference.models``.
    nchains : int
        Number of chains to use in the sampler.
    ntemps : int, optional
        Number of temperatures to use in the sampler. A geometrically-spaced
        temperature ladder with the gievn number of levels will be constructed
        based on the number of parameters. If not provided, must provide
        ``betas``.
    betas : array, optional
        An array of inverse temperature values to be used in for the
        temperature ladder. If not provided, must provide ``ntemps``.
    proposals : list, optional
        List of proposals to use. Any parameters that do not have a proposal
        provided will use the ``default_propsal``. **Note:** proposals should
        be specified for the sampling parameters, not the
        variable parameters.
    default_proposal : an epsie.Proposal class, optional
        The default proposal to use for parameters not in ``proposals``.
        Default is :py:class:`epsie.proposals.Normal`.
    default_proposal_args : dict, optional
        Dictionary of arguments to pass to the default proposal.
    swap_interval : int, optional
        The number of iterations between temperature swaps. Default is 1.
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
    name = "epsie"
    _io = EpsieFile
    burn_in_class = MultiTemperedMCMCBurnInTests

    def __init__(self, model, nchains, ntemps=None, betas=None,
                 proposals=None, default_proposal=None,
                 default_proposal_args=None, seed=None,
                 swap_interval=1,
                 checkpoint_interval=None, checkpoint_signal=None,
                 loglikelihood_function=None,
                 nprocesses=1, use_mpi=False):

        # create the betas if not provided
        if betas is None:
            betas = default_beta_ladder(len(model.variable_params),
                                        ntemps=ntemps)
        self.model = model
        # create a wrapper for calling the model
        model_call = _EpsieCallModel(model, loglikelihood_function)

        # these are used to help paralleize over multiple cores / MPI
        models._global_instance = model_call
        model_call = models._call_global_model

        # Set up the pool
        pool = choose_pool(mpi=use_mpi, processes=nprocesses)

        # initialize the sampler
        self._sampler = ParallelTemperedSampler(
            model.sampling_params, model_call, nchains, betas=betas,
            swap_interval=swap_interval,
            proposals=proposals, default_proposal=default_proposal,
            default_proposal_args=default_proposal_args,
            seed=seed, pool=pool)
        # set other parameters
        self.nchains = nchains
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
    def betas(self):
        """The inverse temperatures being used."""
        return self._sampler.betas

    @property
    def seed(self):
        """The seed used for epsie's random bit generator.

        This is not the same as the seed used for the prior distributions.
        """
        return self._sampler.seed

    @property
    def swap_interval(self):
        """Number of iterations between temperature swaps."""
        return self._sampler.swap_interval

    @staticmethod
    def compute_acf(filename, **kwargs):
        r"""Computes the autocorrelation function.

        Calls :py:func:`base_multitemper.compute_acf`; see that
        function for details.

        Parameters
        ----------
        filename : str
            Name of a samples file to compute ACFs for.
        \**kwargs :
            All other keyword arguments are passed to
            :py:func:`base_multitemper.compute_acf`.

        Returns
        -------
        dict :
            Dictionary of arrays giving the ACFs for each parameter. The arrays
            will have shape ``ntemps x nchains x niterations``.
        """
        return compute_acf(filename, **kwargs)

    @staticmethod
    def compute_acl(filename, **kwargs):
        r"""Computes the autocorrelation length.

        Calls :py:func:`base_multitemper.compute_acl`; see that
        function for details.

        Parameters
        -----------
        filename : str
            Name of a samples file to compute ACLs for.
        \**kwargs :
            All other keyword arguments are passed to
            :py:func:`base_multitemper.compute_acl`.

        Returns
        -------
        dict
            A dictionary of ntemps-long arrays of the ACLs of each parameter.
        """
        return compute_acl(filename, **kwargs)

    @property
    def acl(self):  # pylint: disable=invalid-overridden-method
        """The autocorrelation lengths of the chains.
        """
        return acl_from_raw_acls(self.raw_acls)

    @property
    def effective_nsamples(self):  # pylint: disable=invalid-overridden-method
        """The effective number of samples post burn-in that the sampler has
        acquired so far.
        """
        act = self.act
        if act is None:
            act = numpy.inf
        if self.burn_in is None:
            start_iter = 0
        else:
            start_iter = self.burn_in.burn_in_iteration
        nperchain = nsamples_in_chain(start_iter, act, self.niterations)
        if self.burn_in is not None:
            # ensure that any chain not burned in has zero samples
            nperchain[~self.burn_in.is_burned_in] = 0
            # and that any chain that is burned in has at least one sample
            nperchain[self.burn_in.is_burned_in & (nperchain < 1)] = 1
        return int(nperchain.sum())

    @property
    def samples(self):
        """A dict mapping ``variable_params`` to arrays of samples currently
        in memory.

        The arrays have shape ``ntemps x nchains x niterations``.

        The dictionary also contains sampling parameters.
        """
        samples = epsie.array2dict(self._sampler.positions)
        # apply boundary conditions
        samples = self.model.prior_distribution.apply_boundary_conditions(
            **samples)
        # apply transforms to go to model's variable params space
        if self.model.sampling_transforms is not None:
            samples = self.model.sampling_transforms.apply(
                samples, inverse=True)
        return samples

    @property
    def model_stats(self):
        """A dict mapping the model's ``default_stats`` to arrays of values.

        The arrays have shape ``ntemps x nchains x niterations``.
        """
        return epsie.array2dict(self._sampler.blobs)

    def clear_samples(self):
        """Clears the chain and blobs from memory.
        """
        # store the iteration that the clear is occuring on
        self._lastclear = self.niterations
        self._itercounter = 0
        # now clear the sampler
        self._sampler.clear()

    def set_state_from_file(self, filename):
        """Sets the state of the sampler back to the instance saved in a file.
        """
        with self.io(filename, 'r') as fp:
            # get the numpy state
            numpy_rstate_group = '/'.join([fp.sampler_group,
                                           'numpy_random_state'])
            rstate = fp.read_random_state(group=numpy_rstate_group)
            # set the sampler state for epsie
            self._sampler.set_state_from_checkpoint(fp, path=fp.sampler_group)
        # set the global numpy random state for pycbc
        numpy.random.set_state(rstate)

    def set_p0(self, samples_file=None, prior=None):
        p0 = super(EpsieSampler, self).set_p0(samples_file=samples_file,
                                              prior=prior)
        self._sampler.start_position = p0

    @property
    def pos(self):
        """A dictionary of the current chain positions."""
        # we override BaseMCMC's pos property because this can be directly
        # retrieved from epsie
        return self._sampler.current_positions

    def run_mcmc(self, niterations):
        """Advance the chains for a number of iterations.

        Parameters
        ----------
        niterations : int
            Number of samples to get from sampler.
        """
        self._sampler.run(niterations)

    def write_results(self, filename):
        """Writes samples, model stats, acceptance ratios, and random state
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
            # write accpetance ratio
            acceptance = self._sampler.acceptance
            fp.write_acceptance_ratio(acceptance['acceptance_ratio'],
                                      last_iteration=self.niterations)
            # write temperature data
            if self.ntemps > 1:
                temp_ar = self._sampler.temperature_acceptance
                temp_swaps = self._sampler.temperature_swaps
                fp.write_temperature_data(temp_swaps, temp_ar,
                                          self.swap_interval,
                                          last_iteration=self.niterations)
            # write numpy's global state (for the distributions)
            numpy_rstate_group = '/'.join([fp.sampler_group,
                                           'numpy_random_state'])
            fp.write_random_state(group=numpy_rstate_group)
            # write the sampler's state
            self._sampler.checkpoint(fp, path=fp.sampler_group)

    def finalize(self):
        pass

    @classmethod
    def from_config(cls, cp, model, output_file=None, nprocesses=1,
                    use_mpi=False):
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
        * ``swap-interval`` :
            The number of iterations between temperature swaps. Default is 1.

        Jump proposals must be provided for every sampling
        parameter. These are retrieved from subsections
        ``[jump_proposal-{params}]``, where params is a
        :py:const:`pycbc.VARARGS_DELIM` separated list of parameters the
        proposal should be used for. See
        :py:func:`inference.jump.epsie_proposals_from_config` for
        details.

        .. note::
            Jump proposals should be specified for **sampling parameters**,
            not **variable parameters**.

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
        EpsiePTSampler :
            The sampler instance.
        """
        section = "sampler"
        # check name
        assert cp.get(section, "name") == cls.name, (
            "name in section [sampler] must match mine")
        nchains = int(cp.get(section, "nchains"))
        seed = get_optional_arg_from_config(cp, section, 'seed', dtype=int)
        ntemps, betas = cls.betas_from_config(cp, section)
        # get the swap interval
        swap_interval = get_optional_arg_from_config(cp, section,
                                                     'swap-interval',
                                                     dtype=int)
        if swap_interval is None:
            swap_interval = 1
        # get the checkpoint interval, if it's specified
        checkpoint_interval = cls.checkpoint_from_config(cp, section)
        checkpoint_signal = cls.ckpt_signal_from_config(cp, section)
        # get the loglikelihood function
        logl = get_optional_arg_from_config(cp, section, 'logl-function')
        # get the proposals
        proposals = epsie_proposals_from_config(cp)
        # check that all of the sampling parameters have a specified
        # proposal
        sampling_params = set(model.sampling_params)
        proposal_params = set(param for prop in proposals
                              for param in prop.parameters)
        missing = sampling_params - proposal_params
        if missing:
            raise ValueError("Missing jump proposals for sampling parameters "
                             "{}".format(', '.join(missing)))
        # initialize
        obj = cls(model, nchains,
                  ntemps=ntemps, betas=betas, proposals=proposals,
                  swap_interval=swap_interval, seed=seed,
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
        # Set up the output file
        setup_output(obj, output_file)
        if obj.new_checkpoint:
            obj.set_start_from_config(cp)
        else:
            obj.resume_from_checkpoint()
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
        logp = self.model.logprior
        if logp == -numpy.inf:
            # don't try to call the log likelihood if the prior rules it out
            logl = numpy.nan
        else:
            logl = getattr(self.model, self.loglikelihood_function)
        return logl, logp, self.model.current_stats
