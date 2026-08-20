# Copyright (C) 2026  Alex Nitz
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
This modules provides classes and functions for using the nautilus sampler
package for parameter estimation.
"""

import logging
import os

import h5py

from pycbc.inference.io.nautilus import NautilusFile
from pycbc.pool import choose_pool
from .base import (BaseSampler, setup_output)
from .base_mcmc import get_optional_arg_from_config
from .base_cube import setup_calls


def sync_state(src, dst):
    """Mirrors the hdf group ``src`` into the group ``dst``, replacing what
    is already there.
    """
    for name in set(dst) - set(src):
        del dst[name]
    for name, item in src.items():
        if isinstance(item, h5py.Group):
            sync_state(item, dst.require_group(name))
        else:
            if name in dst:
                del dst[name]
            src.copy(name, dst, name=name)
    dst.attrs.update(src.attrs)


#
# =============================================================================
#
#                                   Samplers
#
# =============================================================================
#

class NautilusSampler(BaseSampler):
    """This class is used to construct a nautilus sampler from the nautilus
    package.

    Parameters
    ----------
    model : model
        A model from ``pycbc.inference.models``.
    nlive : int
        Number of live points to use in the sampler.
    nprocesses : int, optional
        Number of parallel processes to use. Default is 1.
    checkpoint_time_interval : float, optional
        Seconds to sample for between checkpoints. Default (None) is to run
        to convergence and checkpoint once at the end.
    """
    name = "nautilus"
    _io = NautilusFile

    def __init__(self, model, nlive, nprocesses=1, use_mpi=False,
                 loglikelihood_function=None, run_kwds=None, extra_kwds=None,
                 checkpoint_time_interval=None):
        import nautilus
        super().__init__(model)
        self.nautilus = nautilus
        self.loglikelihood_call, self.prior_call = setup_calls(
            model, loglikelihood_function=loglikelihood_function,
            copy_prior=True)
        self.pool = choose_pool(mpi=use_mpi, processes=nprocesses)
        self.nlive = nlive
        self.ndim = len(model.sampling_params)
        self.run_kwds = run_kwds or {}
        self.extra_kwds = extra_kwds or {}
        self.checkpoint_time_interval = checkpoint_time_interval
        self._sampler = None

        # nautilus wraps cyclic parameters around the edge of the unit cube
        cyclic = model.prior_distribution.cyclic
        self.periodic = [i for i, param in enumerate(self.variable_params)
                         if param in cyclic]
        for i in self.periodic:
            logging.info('Param: %s will be cyclic', self.variable_params[i])

    def setup_sampler(self, filepath=None):
        """Constructs the underlying sampler, resuming from the nautilus
        state in ``filepath`` if one is given.

        ``filepath`` is None during a run, which stops nautilus writing a
        file of its own; ``checkpoint`` writes its state instead.
        """
        self._sampler = self.nautilus.Sampler(
            self.prior_call, self.loglikelihood_call,
            n_dim=self.ndim, n_live=self.nlive,
            periodic=self.periodic or None, pool=self.pool,
            filepath=filepath, resume=filepath is not None,
            **self.extra_kwds)

    def run(self):
        if self._sampler is None:
            self.setup_sampler()
        if not self.checkpoint_time_interval:
            self._sampler.run(**self.run_kwds)
            return
        # Sample in chunks so that pycbc decides when to checkpoint. run
        # picks up where it left off, and reports whether it stopped because
        # it converged or because it ran into one of its limits.
        run_kwds = self.run_kwds.copy()
        remaining = run_kwds.pop('timeout', float('inf'))
        while remaining > 0:
            n_like = self._sampler.n_like
            interval = min(self.checkpoint_time_interval, remaining)
            done = self._sampler.run(timeout=interval, **run_kwds)
            remaining -= interval
            if done:
                break
            if self._sampler.n_like == n_like:
                logging.info("nautilus is not making progress, stopping")
                break
            logging.info("Checkpointing after %s likelihood calls",
                         self._sampler.n_like)
            self.checkpoint()

    @property
    def io(self):
        return self._io

    @classmethod
    def from_config(cls, cp, model, output_file=None, nprocesses=1,
                    use_mpi=False, loglikelihood_function=None):
        """Loads the sampler from the given config file. Many options are
        directly passed to the underlying nautilus sampler, see the official
        nautilus documentation for more details on these.

        The following options are retrieved in the ``[sampler]`` section:

        * ``name = STR``:
            Required. This must match the sampler's name.
        * ``nlive = INT``:
            Required. The number of live points to use.
        * ``n_update = INT``:
            Likelihood calls between bound updates. Default is ``nlive``.
        * ``enlarge_per_dim = FLOAT``:
            Factor to enlarge the ellipsoid volume by, per dimension.
        * ``n_points_min = INT``:
            Minimum size of each neural network training set.
        * ``split_threshold = FLOAT``:
            Threshold for splitting the multi-ellipsoidal bound.
        * ``n_networks = INT``:
            Number of neural networks per ensemble.
        * ``n_batch = INT``:
            Number of likelihood evaluations performed at once.
        * ``n_like_new_bound = INT``:
            Maximum likelihood calls before a new bound is created.
        * ``seed = INT``:
            Seed for nautilus' random number generator.
        * ``f_live = FLOAT``:
            Maximum fraction of the evidence in the live set before a new
            bound is built.
        * ``n_shell = INT``:
            Minimum number of points per shell.
        * ``n_eff = FLOAT``:
            Minimum effective sample size; the main stopping condition.
        * ``n_like_max = INT``:
            Maximum number of likelihood calls.
        * ``discard_exploration = BOOL``:
            Discard the exploration phase. Recommended by nautilus when
            publishing evidence estimates.
        * ``timeout = FLOAT``:
            Timeout in seconds for the sampling phase.
        * ``checkpoint_time_interval = FLOAT``:
            Seconds to sample for between checkpoints. Default is to run to
            convergence and checkpoint once at the end.
        * ``loglikelihood-function``:
            The model attribute to use for the loglikelihood. Default is
            ``loglikelihood``.
        """
        section = "sampler"
        # check name
        assert cp.get(section, "name") == cls.name, (
            "name in section [sampler] must match mine")
        nlive = int(cp.get(section, "nlive"))
        loglikelihood_function = \
            get_optional_arg_from_config(cp, section, 'loglikelihood-function')

        # optional arguments for nautilus.Sampler
        cargs = {'n_update': int,
                 'enlarge_per_dim': float,
                 'n_points_min': int,
                 'split_threshold': float,
                 'n_networks': int,
                 'n_batch': int,
                 'n_like_new_bound': int,
                 'seed': int,
                 }

        # optional arguments for nautilus.Sampler.run
        rargs = {'f_live': float,
                 'n_shell': int,
                 'n_eff': float,
                 'n_like_max': int,
                 'discard_exploration': bool,
                 'timeout': float,
                 }

        def given(types):
            return {opt: (cp.getboolean(section, opt) if dtype is bool
                          else dtype(cp.get(section, opt)))
                    for opt, dtype in types.items()
                    if cp.has_option(section, opt)}

        checkpoint_time_interval = None
        if cp.has_option(section, 'checkpoint_time_interval'):
            checkpoint_time_interval = float(
                cp.get(section, 'checkpoint_time_interval'))

        obj = cls(model, nlive=nlive, nprocesses=nprocesses, use_mpi=use_mpi,
                  loglikelihood_function=loglikelihood_function,
                  run_kwds=given(rargs), extra_kwds=given(cargs),
                  checkpoint_time_interval=checkpoint_time_interval)
        setup_output(obj, output_file, check_nsamples=False)
        if not obj.new_checkpoint:
            obj.resume_from_checkpoint()
        return obj

    def checkpoint(self):
        """Writes nautilus' state to the checkpoint and backup files.

        Nautilus can only serialize itself to a file of its own, so its state
        is written to a scratch file and copied into the sampler group of our
        files, keeping the layout nautilus wrote so that it reads it back
        itself on resume.
        """
        scratch = self.checkpoint_file + '.nautilus.h5'
        self._sampler.write(scratch, overwrite=True)
        try:
            with h5py.File(scratch, 'r') as state:
                for fn in [self.checkpoint_file, self.backup_file]:
                    with self.io(fn, 'a') as fp:
                        sync_state(state, fp.require_group(self.io.state_path))
        finally:
            os.remove(scratch)
        for fn in [self.checkpoint_file, self.backup_file]:
            self.write_results(fn)

    def resume_from_checkpoint(self):
        """Rebuilds nautilus from the state in the checkpoint file."""
        scratch = self.checkpoint_file + '.nautilus.h5'
        with self.io(self.checkpoint_file, 'r') as fp:
            with h5py.File(scratch, 'w') as state:
                sync_state(fp[self.io.state_path], state)
        try:
            self.setup_sampler(filepath=scratch)
        finally:
            os.remove(scratch)
        # the scratch file is gone; checkpoint writes the state from here on
        self._sampler.filepath = None
        logging.info("Resumed nautilus after %s likelihood calls",
                     self._sampler.n_like)

    def finalize(self):
        logging.info("log Z, dlog Z: %s, %s", self.logz, self.logz_err)
        self.checkpoint()

    @property
    def model_stats(self):
        return {}

    @property
    def samples(self):
        """Returns the raw nested samples and their log weights."""
        points, logwt, loglikelihood = self._sampler.posterior()
        samples = {param: points[:, i]
                   for i, param in enumerate(self.variable_params)}
        # nautilus normalises its log weights to sum to one, whereas the
        # other nested samplers store weights that sum to the evidence
        samples['logwt'] = logwt + self.logz
        samples['loglikelihood'] = loglikelihood
        return samples

    def write_results(self, filename):
        """Writes samples and the log evidence to the given file."""
        with self.io(filename, 'a') as fp:
            fp.write_raw_samples(self.samples)
            fp.write_logevidence(self.logz, self.logz_err)

    @property
    def logz(self):
        """Return the log evidence estimated by nautilus."""
        return self._sampler.log_z

    @property
    def logz_err(self):
        """Return the error on the log evidence.

        Nautilus does not report this directly; its documentation gives the
        statistical uncertainty on log Z as one over the square root of the
        effective sample size.
        """
        return self._sampler.effective_sample_size() ** -0.5
