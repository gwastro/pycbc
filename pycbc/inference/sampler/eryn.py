# Copyright (C) 2026 Alexander Nitz
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


"""This module provides classes for using the ``eryn`` sampler for parameter
estimation.

This is a fixed-dimension integration: ``eryn`` is driven as a parallel-
tempered affine-invariant ensemble sampler (a single branch with a single
leaf, no reversible jump), which makes it a close analog of the ``ptemcee``
sampler. The reversible-jump / trans-dimensional capabilities of ``eryn`` are
not exposed here.
"""


import ast
import logging
import shlex

import numpy

# eryn 1.2.6 uses numpy aliases removed in numpy 2.0; restore them.
if not hasattr(numpy, "in1d"):
    numpy.in1d = numpy.isin
if not hasattr(numpy, "trapz"):
    numpy.trapz = numpy.trapezoid

import eryn.moves as eryn_moves
from eryn.ensemble import EnsembleSampler
from eryn.state import State

from pycbc.inference.io import ErynFile
from pycbc.inference.evidence import (
    ladder_thermodynamic_integration, mean_logl_by_temperature)
from pycbc.pool import choose_pool

from .. import models
from ..burn_in import EnsembleMultiTemperedMCMCBurnInTests
from .base import BaseSampler, setup_output
from .base_mcmc import (
    BaseMCMC,
    EnsembleSupport,
    get_optional_arg_from_config,
    raw_samples_to_dict,
)
from .base_multitemper import (
    ensemble_compute_acf,
    ensemble_compute_acl,
    read_betas_from_hdf,
)

# fixed-dimension model: a single branch with a single leaf
BRANCH = "model_0"


class _ErynPrior(object):
    """Wraps a pycbc model's prior as an ``eryn`` prior container.

    ``eryn`` evaluates the prior in the main process by calling ``logpdf`` on
    an array of shape ``(N, ndim)`` (the flattened ``ntemps x nwalkers``
    positions) and expects an ``(N,)`` array of log-prior values back. Each row
    is forwarded to the pycbc model's ``logprior`` via the global model-call
    instance -- the same wiring used by the ptemcee sampler -- so that the
    prior (including any constraints) is exactly the model's prior.

    Because ``eryn`` evaluates the prior serially in the main process (only the
    likelihood is dispatched to the pool), this object does not need to be
    picklable.
    """

    def __init__(self, model):
        self.model = model
        self.ndim = len(model.variable_params)
        # required by eryn's prior setter
        self.key_order = list(range(self.ndim))

    def logpdf(self, x):
        """Returns the log prior for each of the ``N`` rows of ``x``."""
        x = numpy.atleast_2d(x)
        return numpy.array([models._call_global_model_logprior(row)
                            for row in x])

    def rvs(self, size=1):
        """Draws samples from the model prior in the sampling frame.

        This is only used if ``eryn`` needs to generate points itself; for the
        fixed-dimension runs here the initial positions are always provided, so
        this is a safety fallback.
        """
        if isinstance(size, (tuple, list)):
            shape = tuple(int(s) for s in size)
        else:
            shape = (int(size),)
        nsamples = int(numpy.prod(shape))
        draws = self.model.prior_rvs(size=nsamples)
        if self.model.sampling_transforms is not None:
            draws = self.model.sampling_transforms.apply(draws)
        out = numpy.zeros((nsamples, self.ndim))
        for i, param in enumerate(self.model.sampling_params):
            out[:, i] = draws[param]
        return out.reshape(shape + (self.ndim,))


class ErynSampler(EnsembleSupport, BaseMCMC, BaseSampler):
    """This class is used to construct the parallel-tempered eryn sampler.

    Parameters
    ----------
    model : model
        A model from ``pycbc.inference.models``.
    nwalkers : int
        Number of walkers to use in the sampler.
    ntemps : int, optional
        The number of temperatures to use. Either this, ``Tmax``, or ``betas``
        must be specified.
    Tmax : float, optional
        The maximum temperature to use. May be used with ``ntemps``.
    betas : array of float, optional
        The (initial) inverse temperatures to use. Overrides ``ntemps`` and
        ``Tmax`` if provided.
    adaptive : bool, optional
        Whether to use adaptive temperature levels. Default is False. Note that
        the temperature ladder must be frozen for the evidence estimate to be
        valid; ``finalize`` computes the evidence from the betas stored on disk
        regardless.
    moves : list, optional
        A list of ``(eryn.moves.Move, weight)`` tuples giving the proposal
        schedule. If not provided, eryn's default (a single ``StretchMove``)
        is used.
    moves_labels : str, optional
        A human-readable label for the moves used, stored in the output file.
    loglikelihood_function : str, optional
        The attribute of the model to use for the loglikelihood. Default is
        ``loglikelihood``.
    nprocesses : int, optional
        The number of parallel processes to use. Default is 1.
    use_mpi : bool, optional
        Use MPI for parallelization. Default is False.
    """
    name = "eryn"
    _io = ErynFile
    burn_in_class = EnsembleMultiTemperedMCMCBurnInTests

    def __init__(self, model, nwalkers, ntemps=None, Tmax=None, betas=None,
                 adaptive=False, moves=None, moves_labels=None,
                 loglikelihood_function=None,
                 checkpoint_interval=None, checkpoint_signal=None,
                 nprocesses=1, use_mpi=False):
        self.model = model
        ndim = len(model.variable_params)
        self._ndim = ndim
        # build the tempering configuration
        if ntemps is None and Tmax is None and betas is None:
            raise ValueError("must provide either ntemps/Tmax or betas")
        tempering_kwargs = {'adaptive': adaptive}
        if betas is not None:
            tempering_kwargs['betas'] = numpy.asarray(betas, dtype=float)
        else:
            if ntemps is not None:
                tempering_kwargs['ntemps'] = int(ntemps)
            if Tmax is not None:
                tempering_kwargs['Tmax'] = float(Tmax)
        # likelihood is dispatched to the pool via the picklable global
        # instance (as in ptemcee); the prior is evaluated in the main process
        if loglikelihood_function is None:
            loglikelihood_function = 'loglikelihood'
        model_call = models.CallModel(model, loglikelihood_function,
                                      return_all_stats=False)
        models._global_instance = model_call
        loglike_call = models._call_global_model
        self.pool = choose_pool(mpi=use_mpi, processes=nprocesses)
        priors = {BRANCH: _ErynPrior(model)}
        self._sampler = EnsembleSampler(
            nwalkers, ndim, loglike_call, priors,
            tempering_kwargs=tempering_kwargs,
            branch_names=[BRANCH],
            nleaves_max=1, nleaves_min=1,
            moves=moves,
            pool=self.pool,
            vectorize=False,
        )
        self.nwalkers = nwalkers
        self._ntemps = self._sampler.ntemps
        self._starting_betas = numpy.asarray(
            self._sampler.temperature_control.betas, dtype=float).copy()
        self._moves_labels = moves_labels if moves_labels is not None \
            else 'StretchMove'
        self._checkpoint_interval = checkpoint_interval
        self._checkpoint_signal = checkpoint_signal
        # eryn state, carried across checkpoints; seeded on first run_mcmc
        self._state = None

    @property
    def io(self):
        return self._io

    @property
    def ntemps(self):
        """The number of temperatures that are set."""
        return self._ntemps

    @property
    def base_shape(self):
        return (self.ntemps, self.nwalkers,)

    @property
    def starting_betas(self):
        """The betas that were used at startup."""
        return self._starting_betas

    @property
    def adaptive(self):
        """Whether or not the betas are adapted."""
        tc = self._sampler.temperature_control
        return bool(getattr(tc, 'adaptive', False))

    @property
    def moves_labels(self):
        """A label describing the proposal schedule that was used."""
        return self._moves_labels

    @property
    def betas(self):
        """Returns the beta history currently in memory.

        The returned array has shape ``ntemps x niterations``.
        """
        betas = self._sampler.get_betas()
        if betas is None:
            # no history stored (single temperature); use a constant ladder
            betas = numpy.ones((self.niterations, self.ntemps))
        # eryn stores betas as niterations x ntemps
        return numpy.asarray(betas).transpose()

    def _initial_state(self, betas=None):
        """Builds an eryn ``State`` from the current ``_p0`` positions."""
        # _p0 is ntemps x nwalkers x ndim; eryn wants an extra leaf axis
        coords = self._p0[:, :, numpy.newaxis, :]
        if betas is None:
            betas = numpy.asarray(
                self._sampler.temperature_control.betas, dtype=float).copy()
        return State({BRANCH: coords}, betas=betas)

    def run_mcmc(self, niterations):
        """Advance the ensemble by ``niterations`` iterations."""
        if self._state is None:
            self._state = self._initial_state()
        self._state = self._sampler.run_mcmc(
            self._state, niterations, thin_by=1, store=True, progress=False)

    def clear_samples(self):
        """Clears the chain from memory."""
        self._lastclear = self.niterations
        self._itercounter = 0
        # reset eryn's in-memory backend; positions live on self._state and the
        # (possibly adapting) ladder on the sampler, so both survive the reset
        self._sampler.reset(
            ntemps=self.ntemps, branch_names=[BRANCH], nleaves_max=1,
            rj=False, moves=self._sampler.move_keys,
            key_order=self._sampler.key_order)

    @property
    def samples(self):
        """A dict mapping ``variable_params`` to arrays of samples currently
        in memory. The arrays have shape ``ntemps x nwalkers x niterations``.
        """
        # get_chain: niter x ntemps x nwalkers x nleaves(=1) x ndim; drop the
        # leaf axis and transpose to ntemps x nwalkers x niter x ndim
        raw = self._sampler.get_chain()[BRANCH][:, :, :, 0, :]
        raw_samples = raw.transpose((1, 2, 0, 3))
        return raw_samples_to_dict(self, raw_samples)

    @property
    def model_stats(self):
        """Returns the loglikelihood and logprior as a dict of arrays.

        The returned arrays have shape ``ntemps x nwalkers x niterations``.

        Like the ptemcee sampler, eryn does not carry pycbc's blob stats, so
        only the loglikelihood and logprior are returned, with the logjacobian
        set to zero (it is corrected on disk in ``finalize`` if sampling
        transforms were used).
        """
        # niterations x ntemps x nwalkers -> ntemps x nwalkers x niterations
        logl = self._sampler.get_log_like().transpose((1, 2, 0))
        logp = self._sampler.get_log_prior().transpose((1, 2, 0))
        logjacobian = numpy.zeros(logp.shape)
        return {'loglikelihood': logl, 'logprior': logp,
                'logjacobian': logjacobian}

    def set_state_from_file(self, filename):
        """Sets the state of the sampler back to the instance in a file."""
        with self.io(filename, 'r') as fp:
            rstate = fp.read_random_state()
            betas = fp.read_betas(iteration=-1)
        numpy.random.set_state(rstate)
        # restore the ladder and re-seed the state from _p0 (set by set_p0)
        betas = numpy.asarray(betas, dtype=float)
        self._sampler.temperature_control.betas = betas.copy()
        self._state = self._initial_state(betas=betas.copy())

    @staticmethod
    def compute_acf(filename, **kwargs):
        """Computes the autocorrelation function.

        Calls :py:func:`base_multitemper.ensemble_compute_acf`; see that
        function for details.
        """
        return ensemble_compute_acf(filename, **kwargs)

    @staticmethod
    def compute_acl(filename, **kwargs):
        """Computes the autocorrelation length.

        Calls :py:func:`base_multitemper.ensemble_compute_acl`; see that
        function for details.
        """
        return ensemble_compute_acl(filename, **kwargs)

    @staticmethod
    def _parse_moves(cp, section):
        """Parses the (optional) eryn move schedule from the config file.

        The ``moves`` option in the ``[sampler]`` section is a space-separated
        list of ``MoveClass[:weight]`` entries, where ``MoveClass`` is the name
        of a class in :py:mod:`eryn.moves`. Keyword arguments for a move can be
        given in a dedicated ``[eryn_move-MoveClass]`` section; values are
        parsed as Python literals (e.g. ``a = 2.0``). For example::

            [sampler]
            name = eryn
            moves = StretchMove:0.8 GaussianMove:0.2

            [eryn_move-StretchMove]
            a = 2.0

        Returns
        -------
        moves : list of (Move, weight) tuples or None
            The move schedule to pass to eryn, or ``None`` to use eryn's
            default (a single ``StretchMove``).
        label : str
            A human-readable label describing the schedule.
        """
        if not cp.has_option(section, 'moves'):
            return None, 'StretchMove'
        moves = []
        labels = []
        for token in shlex.split(cp.get(section, 'moves')):
            if ':' in token:
                clsname, weight = token.rsplit(':', 1)
                weight = float(weight)
            else:
                clsname, weight = token, 1.0
            try:
                cls = getattr(eryn_moves, clsname)
            except AttributeError:
                raise ValueError(
                    "unrecognized eryn move '{}'; must be a class in "
                    "eryn.moves".format(clsname)) from None
            # read any kwargs for this move from its dedicated section
            move_kwargs = {}
            move_section = 'eryn_move-{}'.format(clsname)
            if cp.has_section(move_section):
                for opt in cp.options(move_section):
                    move_kwargs[opt] = ast.literal_eval(
                        cp.get(move_section, opt))
            moves.append((cls(**move_kwargs), weight))
            labels.append(token)
        return moves, ' '.join(labels)

    @classmethod
    def from_config(cls, cp, model, output_file=None, nprocesses=1,
                    use_mpi=False):
        """Loads the sampler from the given config file.

        The following options are retrieved in the ``[sampler]`` section:

        * ``name = STR`` :
            Required. This must match the sampler's name (``eryn``).
        * ``nwalkers = INT`` :
            Required. The number of walkers to use.
        * ``ntemps = INT`` :
            The number of temperatures to use. May be combined with ``tmax``.
            Either this, ``tmax``, ``betas``, or ``betas-file`` must be given.
        * ``tmax = FLOAT`` :
            The maximum temperature to use.
        * ``betas = FLOAT1 FLOAT2 [...]`` :
            Space-separated list of (initial) inverse temperatures. Mutually
            exclusive with ``ntemps``/``tmax``.
        * ``betas-file = STR`` :
            Path to an hdf file containing the inverse temperatures, retrieved
            from the file's ``.attrs['betas']``. Mutually exclusive with
            ``betas``.
        * ``adaptive =`` :
            If present, temperature adaptation is turned on.
        * ``moves = STR`` :
            Optional. A space-separated list of ``MoveClass[:weight]`` entries
            naming proposals from :py:mod:`eryn.moves`. See
            :py:meth:`_parse_moves`. If not given, eryn's default is used.
        * ``niterations`` / ``effective-nsamples`` :
            The stopping condition (exactly one required).
        * ``thin-interval`` / ``max-samples-per-chain`` :
            Optional pre-thinning (mutually exclusive).
        * ``checkpoint-interval`` / ``checkpoint-signal`` :
            Checkpointing controls. ``checkpoint-interval`` is required when
            using ``effective-nsamples``.
        * ``logl-function = STR`` :
            The model attribute to use for the loglikelihood. Default is
            ``loglikelihood``.

        Settings for burn-in tests are read from ``[sampler-burn_in]``.
        """
        section = "sampler"
        # check name
        assert cp.get(section, "name") == cls.name, (
            "name in section [sampler] must match mine")
        # get the number of walkers to use
        nwalkers = int(cp.get(section, "nwalkers"))
        # checkpointing
        checkpoint_interval = cls.checkpoint_from_config(cp, section)
        checkpoint_signal = cls.ckpt_signal_from_config(cp, section)
        optargs = {}
        # temperature settings
        ntemps = get_optional_arg_from_config(cp, section, 'ntemps', int)
        if ntemps is not None:
            optargs['ntemps'] = ntemps
        tmax = get_optional_arg_from_config(cp, section, 'tmax', float)
        if tmax is not None:
            optargs['Tmax'] = tmax
        betas = get_optional_arg_from_config(cp, section, 'betas')
        if betas is not None:
            betas = numpy.sort(list(map(float, shlex.split(betas))))[::-1]
            optargs['betas'] = betas
        betas_file = get_optional_arg_from_config(cp, section, 'betas-file')
        if betas_file is not None:
            optargs['betas'] = read_betas_from_hdf(betas_file)
        # consistency checks
        if betas is not None and betas_file is not None:
            raise ValueError("provide either betas or betas-file, not both")
        if 'betas' in optargs and (ntemps is not None or tmax is not None):
            raise ValueError("provide either ntemps/tmax or betas/betas-file, "
                             "not both")
        # adaptation
        adaptive = get_optional_arg_from_config(cp, section, 'adaptive')
        optargs['adaptive'] = adaptive is not None
        # proposal (moves) schedule
        moves, moves_labels = cls._parse_moves(cp, section)
        # loglikelihood function
        logl = get_optional_arg_from_config(cp, section, 'logl-function')
        obj = cls(model, nwalkers,
                  moves=moves, moves_labels=moves_labels,
                  checkpoint_interval=checkpoint_interval,
                  checkpoint_signal=checkpoint_signal,
                  loglikelihood_function=logl, nprocesses=nprocesses,
                  use_mpi=use_mpi, **optargs)
        # set target
        obj.set_target_from_config(cp, section)
        # burn-in
        obj.set_burn_in_from_config(cp)
        # pre-thinning
        obj.set_thin_interval_from_config(cp, section)
        # output file
        setup_output(obj, output_file)
        if not obj.new_checkpoint:
            obj.resume_from_checkpoint()
        else:
            obj.set_start_from_config(cp)
        return obj

    def write_results(self, filename):
        """Writes samples, model stats, betas, and random state to a file."""
        with self.io(filename, 'a') as fp:
            fp.write_samples(self.samples,
                             parameters=self.model.variable_params,
                             last_iteration=self.niterations)
            fp.write_samples(self.model_stats, last_iteration=self.niterations)
            fp.write_random_state()
            fp.write_betas(self.betas, last_iteration=self.niterations)

    @classmethod
    def calculate_logevidence(cls, filename, thin_start=None, thin_end=None,
                              thin_interval=None):
        """Calculates the log evidence from the given file.

        This uses thermodynamic integration over the stored inverse
        temperatures, in the same way as the ptemcee sampler.

        Parameters
        ----------
        filename : str
            Name of the file to read the samples from (an ``ErynFile``).
        thin_start : int, optional
            Index of the sample to begin reading stats.
        thin_interval : int, optional
            Interval to accept every i-th sample.
        thin_end : int, optional
            Index of the last sample to read.

        Returns
        -------
        lnZ : float
            The estimate of the log evidence.
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
            betas = fp.read_betas(thin_start=thin_start,
                                  thin_interval=thin_interval,
                                  thin_end=thin_end)
        # average the loglikelihoods within each unique temperature
        return ladder_thermodynamic_integration(
            *mean_logl_by_temperature(logls, betas))

    def finalize(self):
        """Calculates the log evidence and writes it to the output files.

        If sampling transforms were used, this also corrects the logjacobian
        stored on disk.
        """
        if self.model.sampling_transforms is not None:
            logging.info("Correcting logjacobian values on disk")
            with self.io(self.checkpoint_file, 'r') as fp:
                samples = fp.read_raw_samples(self.variable_params,
                                              thin_start=0,
                                              thin_interval=1, thin_end=None,
                                              temps='all', flatten=False)
            logjacobian = self._correctjacobian(samples)
            for fn in [self.checkpoint_file, self.backup_file]:
                with self.io(fn, "a") as fp:
                    fp[fp.samples_group]['logjacobian'][()] = logjacobian
        logging.info("Calculating log evidence")
        with self.io(self.checkpoint_file, 'r') as fp:
            thin_start = fp.thin_start
            thin_interval = fp.thin_interval
            thin_end = fp.thin_end
        logz, dlogz = self.calculate_logevidence(
            self.checkpoint_file, thin_start=thin_start, thin_end=thin_end,
            thin_interval=thin_interval)
        logging.info("log Z, dlog Z: {}, {}".format(logz, dlogz))
        for fn in [self.checkpoint_file, self.backup_file]:
            with self.io(fn, "a") as fp:
                fp.write_logevidence(logz, dlogz)
