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


#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#
"""
This modules provides classes and functions for using the emcee sampler
packages for parameter estimation.
"""


import numpy
import emcee
from pycbc.pool import choose_pool

from .base import (BaseSampler, setup_output)
from .base_mcmc import (BaseMCMC, EnsembleSupport,
                        ensemble_compute_acf, ensemble_compute_acl,
                        raw_samples_to_dict,
                        blob_data_to_dict, get_optional_arg_from_config)
from ..burn_in import EnsembleMCMCBurnInTests
from pycbc.inference.io import EmceeFile
from .. import models


#
# =============================================================================
#
#                                   Samplers
#
# =============================================================================
#

if emcee.__version__ >= '3.0.0':
    raise ImportError


class EmceeEnsembleSampler(EnsembleSupport, BaseMCMC, BaseSampler):
    """This class is used to construct an MCMC sampler from the emcee
    package's EnsembleSampler.

    Parameters
    ----------
    model : model
        A model from ``pycbc.inference.models``.
    nwalkers : int
        Number of walkers to use in sampler.
    pool : function with map, Optional
        A provider of a map function that allows a function call to be run
        over multiple sets of arguments and possibly maps them to
        cores/nodes/etc.
    """
    name = "emcee"
    _io = EmceeFile
    burn_in_class = EnsembleMCMCBurnInTests

    def __init__(self, model, nwalkers,
                 checkpoint_interval=None, checkpoint_signal=None,
                 logpost_function=None, nprocesses=1, use_mpi=False):

        self.model = model
        # create a wrapper for calling the model
        if logpost_function is None:
            logpost_function = 'logposterior'
        model_call = models.CallModel(model, logpost_function)

        # these are used to help paralleize over multiple cores / MPI
        models._global_instance = model_call
        model_call = models._call_global_model
        pool = choose_pool(mpi=use_mpi, processes=nprocesses)

        # set up emcee
        self.nwalkers = nwalkers
        ndim = len(model.variable_params)
        self._sampler = emcee.EnsembleSampler(nwalkers, ndim, model_call,
                                              pool=pool)
        # emcee uses it's own internal random number generator; we'll set it
        # to have the same state as the numpy generator
        rstate = numpy.random.get_state()
        self._sampler.random_state = rstate
        self._checkpoint_interval = checkpoint_interval
        self._checkpoint_signal = checkpoint_signal

    @property
    def io(self):
        return self._io

    @property
    def base_shape(self):
        return (self.nwalkers,)

    @property
    def samples(self):
        """A dict mapping ``variable_params`` to arrays of samples currently
        in memory.

        The arrays have shape ``nwalkers x niterations``.
        """
        # emcee stores samples to it's chain attribute as a
        # nwalker x niterations x ndim array
        raw_samples = self._sampler.chain
        return raw_samples_to_dict(self, raw_samples)

    @property
    def model_stats(self):
        """A dict mapping the model's ``default_stats`` to arrays of values.

        The returned array has shape ``nwalkers x niterations``.
        """
        stats = self.model.default_stats
        return blob_data_to_dict(stats, self._sampler.blobs)

    def clear_samples(self):
        """Clears the samples and stats from memory.
        """
        # store the iteration that the clear is occuring on
        self._lastclear = self.niterations
        self._itercounter = 0
        # now clear the chain
        self._sampler.reset()
        self._sampler.clear_blobs()

    def set_state_from_file(self, filename):
        """Sets the state of the sampler back to the instance saved in a file.
        """
        with self.io(filename, 'r') as fp:
            rstate = fp.read_random_state()
        # set the numpy random state
        numpy.random.set_state(rstate)
        # set emcee's generator to the same state
        self._sampler.random_state = rstate

    def run_mcmc(self, niterations):
        """Advance the ensemble for a number of samples.

        Parameters
        ----------
        niterations : int
            Number of iterations to run the sampler for.
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
            fp.write_samples(self.model_stats,
                             last_iteration=self.niterations)
            # write accpetance
            fp.write_acceptance_fraction(self._sampler.acceptance_fraction)
            # write random state
            fp.write_random_state(state=self._sampler.random_state)

    def finalize(self):
        """All data is written by the last checkpoint in the run method, so
        this just passes."""
        pass

    @staticmethod
    def compute_acf(filename, **kwargs):
        r"""Computes the autocorrelation function.

        Calls :py:func:`base_mcmc.ensemble_compute_acf`; see that
        function for details.

        Parameters
        ----------
        filename : str
            Name of a samples file to compute ACFs for.
        \**kwargs :
            All other keyword arguments are passed to
            :py:func:`base_mcmc.ensemble_compute_acf`.

        Returns
        -------
        dict :
            Dictionary of arrays giving the ACFs for each parameter. If
            ``per-walker`` is True, the arrays will have shape
            ``nwalkers x niterations``.
        """
        return ensemble_compute_acf(filename, **kwargs)

    @staticmethod
    def compute_acl(filename, **kwargs):
        r"""Computes the autocorrelation length.

        Calls :py:func:`base_mcmc.ensemble_compute_acl`; see that
        function for details.

        Parameters
        -----------
        filename : str
            Name of a samples file to compute ACLs for.
        \**kwargs :
            All other keyword arguments are passed to
            :py:func:`base_mcmc.ensemble_compute_acf`.

        Returns
        -------
        dict
            A dictionary giving the ACL for each parameter.
        """
        return ensemble_compute_acl(filename, **kwargs)

    @classmethod
    def from_config(cls, cp, model, output_file=None, nprocesses=1,
                    use_mpi=False):
        """Loads the sampler from the given config file."""
        section = "sampler"
        # check name
        assert cp.get(section, "name") == cls.name, (
            "name in section [sampler] must match mine")
        # get the number of walkers to use
        nwalkers = int(cp.get(section, "nwalkers"))
        # get the checkpoint interval, if it's specified
        checkpoint_interval = cls.checkpoint_from_config(cp, section)
        checkpoint_signal = cls.ckpt_signal_from_config(cp, section)
        # get the logpost function
        lnpost = get_optional_arg_from_config(cp, section, 'logpost-function')
        obj = cls(model, nwalkers,
                  checkpoint_interval=checkpoint_interval,
                  checkpoint_signal=checkpoint_signal,
                  logpost_function=lnpost, nprocesses=nprocesses,
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
