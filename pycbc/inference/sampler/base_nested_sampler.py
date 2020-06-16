# Copyright (C) 2020 Prayush Kumar
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
"""Provides constructor classes for NestedSampling samplers."""

from __future__ import (absolute_import, division)

import os
import sys
import resource
import signal
import logging
from abc import (ABCMeta, abstractmethod)

from six import add_metaclass

import dill
from pycbc.inference.io import validate_checkpoint_files
from .base_mcmc import get_optional_arg_from_config

#
# =============================================================================
#
#                              BaseNestedSampler definition
#
# =============================================================================
#


@add_metaclass(ABCMeta)
class BaseNestedSampler(object):
    """Abstract base class that provides methods common to Nested Samplers.

    This is not a sampler class itself. Sampler classes can inherit from this
    along with ``BaseSampler``.

    This class provides ``checkpoint`` and ``resume_from_checkpoint`` methods,
    which are some of the abstract methods required by ``BaseSampler``.
    Through them, it adds generic checkpointing facility that can be inherited
    by individual sampler classes with little development overhead.

    This class introduces the following abstract properties and methods:

    * set_sampler_specific_state_from_file(filename)
        Should set the random state of the sampler using the given filename.
        Called by ``set_initial_conditions``.

    Attributes
    ----------
    is_main_process
    checkpoint_interval
    checkpoint_period
    checkpoint_signal
    checkpoint_pickle
    target_eff_nsamples
    effective_nsamples
    """
    _is_main_process = None
    _checkpoint_interval = None
    _checkpoint_period = None
    _checkpoint_signal = None
    _checkpoint_on_signal = None
    _target_eff_nsamples = None

    @property
    def is_main_process(self):
        """Check if this is the main control process.

        To be used when handling one time tasks
        """
        return self._is_main_process

    @property
    def checkpoint_interval(self):
        """The number of iterations to do between checkpoints."""
        return self._checkpoint_interval

    @property
    def checkpoint_period(self):
        """The number of iterations to do between checkpoints."""
        return self._checkpoint_period

    @property
    def checkpoint_signal(self):
        """The signal to send when checkpointing."""
        return self._checkpoint_signal

    @property
    def checkpoint_pickle(self):
        """The signal to send when checkpointing."""
        # pylint: disable=no-member
        return self.checkpoint_file + '.pkl'

    @staticmethod
    def checkpoint_from_config(cp, section):
        """Gets the checkpoint interval from the given config file.

        This looks for 'checkpoint-interval' in the section.

        Parameters
        ----------
        cp : ConfigParser
            Open config parser to retrieve the argument from.
        section : str
            Name of the section to retrieve from.

        Return
        ------
        int or None :
            The checkpoint interval, if it is in the section.
        """
        return get_optional_arg_from_config(cp, section, 'checkpoint-interval',
                                            dtype=int)

    @staticmethod
    def ckpt_signal_from_config(cp, section):
        """Gets the checkpoint signal from the given config file.

        This looks for 'checkpoint-signal' in the section.

        Parameters
        ----------
        cp : ConfigParser
            Open config parser to retrieve the argument from.
        section : str
            Name of the section to retrieve from.

        Return
        ------
        int or None :
            The checkpoint interval, if it is in the section.
        """
        return get_optional_arg_from_config(cp, section, 'checkpoint-signal',
                                            dtype=str)

    def checkpoint_on_signal(self, signum, frame):
        """Interrupt signal handler to checkpoint the sampler.

        Parameters
        ----------
        signum : int
            Open config parser to retrieve the argument from.
        frame : stack frame object
            Name of the section to retrieve from.
        """
        if self.is_main_process:
            self.checkpoint()
        del frame
        self._sampler.pool.close()
        self._sampler.pool.terminate()
        sys.exit(signum)

    def set_state_from_file(self, filename):
        """Sets the state of the sampler to the instance saved in a pickle file.
        """
        with open(filename, 'r') as fin:
            self._sampler = dill.load(fin)  # pylint: disable=no-member

    @abstractmethod
    def set_sampler_specific_state_from_file(self, filename):
        """Sets the state of the sampler to the instance saved at checkpoint.
        """
        # pylint: disable=empty
        pass

    @abstractmethod
    def getstate(self):
        """Strips unserializable attributes of sampler
        """
        # pylint: disable=empty
        pass

    def checkpoint(self):
        """Checkpoint the sampler.

        Algorithm
        ----------
        - Write out samples and data to regular inference.hdf.*
          files used to checkpoint sampler state
        - Strip un-serializable attributes of sampler object
          (typically multiprocessing-related attrs)
        - Serialize sampler object and write to disk as a pickle
        """
        # Write results collected so far
        logging.info("Checkpoint: start")
        self.finalize()  # pylint: disable=no-member

        # Remove multiprocessing / mpi pool-related attributes
        # as they cannot be pickled
        self._sampler.__getstate__ = self.getstate

        # Pickle the sampler
        logging.info("Checkpoint: picking sampler to %s",
                     self.checkpoint_pickle)
        if dill.pickles(self._sampler):
            tmp_filename = self.checkpoint_pickle + ".tmp"
            # guess for maximum recursion stack depth
            max_rec = 0x100000
            # May segfault without this line. 0x100 is a guess at the size of
            # each stack frame.
            each_stack_size = 0x100
            resource.setrlimit(resource.RLIMIT_STACK,
                               [each_stack_size * max_rec,
                                resource.RLIM_INFINITY])
            sys.setrecursionlimit(max_rec)
            with open(tmp_filename, 'w') as fout:
                dill.dump(self._sampler, fout)  # , recurse=True)
            os.rename(tmp_filename, self.checkpoint_pickle)
        else:
            logging.warning("Checkpoint: could not pickle sampler..!")

        # check validity
        logging.info("Validating checkpoint and backup files")
        is_checkpoint_valid = validate_checkpoint_files(
            self.checkpoint_file, self.backup_file)  # pylint: disable=no-member
        if not is_checkpoint_valid:
            raise IOError("error writing to checkpoint file")
        if self.checkpoint_signal is not None:
            # kill myself with the specified signal
            logging.info("Exiting with SIG%s", self.checkpoint_signal)
            # pylint: disable=use-exec
            exec("os.kill(os.getpid(), signal.SIG{})".format(
                self.checkpoint_signal), {}) in globals()
        logging.info("Checkpoint: completed")

    def resume_from_checkpoint(self):
        """Resume the sampler from the checkpoint file.

        Algorithm
        ----------
        - Read in and deserialize the external sampler object
        - Set attributes in it that were not serializable
        """
        logging.info("Resuming sampler from ckpt: %s", self.checkpoint_pickle)
        # First load the sampler object
        self.set_state_from_file(self.checkpoint_pickle)
        # Set state of sampler object that could not be pickled
        self.set_sampler_specific_state_from_file(self.checkpoint_pickle)
        logging.info("Resuming successfully")
