# Copyright (C) 2023  Alex Nitz
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
This modules provides classes and functions for using the snowline sampler
packages for parameter estimation.
"""

import sys
import logging

from pycbc.inference.io.snowline import SnowlineFile
from pycbc.io.hdf import dump_state
from pycbc.pool import use_mpi
from .base import (BaseSampler, setup_output)
from .base_cube import setup_calls


#
# =============================================================================
#
#                                   Samplers
#
# =============================================================================
#

class SnowlineSampler(BaseSampler):
    """This class is used to construct an Snowline sampler from the snowline
    package.

    Parameters
    ----------
    model : model
        A model from ``pycbc.inference.models``
    """
    name = "snowline"
    _io = SnowlineFile

    def __init__(self, model, **kwargs):
        super().__init__(model)

        import snowline
        log_likelihood_call, prior_call = setup_calls(model, copy_prior=True)

        self._sampler = snowline.ReactiveImportanceSampler(
            list(self.model.variable_params),
            log_likelihood_call,
            transform=prior_call)

        do_mpi, _, rank = use_mpi()
        self.main = (not do_mpi) or (rank == 0)
        self.nlive = 0
        self.ndim = len(self.model.variable_params)
        self.kwargs = kwargs

    def run(self):
        self.result = self._sampler.run(**self.kwargs)
        if not self.main:
            sys.exit(0)
        self._sampler.print_results()

    @property
    def io(self):
        return self._io

    @property
    def niterations(self):
        return self.result['niter']

    @classmethod
    def from_config(cls, cp, model, output_file=None, **kwds):
        """
        Loads the sampler from the given config file.
        """
        skeys = {}
        opts = {'num_global_samples': int,
                'num_gauss_samples': int,
                'max_ncalls': int,
                'min_ess': int,
                'max_improvement_loops': int
                }
        for opt_name in opts:
            if cp.has_option('sampler', opt_name):
                value = cp.get('sampler', opt_name)
                skeys[opt_name] = opts[opt_name](value)
        inst = cls(model, **skeys)

        do_mpi, _, rank = use_mpi()
        if not do_mpi or (rank == 0):
            setup_output(inst, output_file)
        return inst

    def checkpoint(self):
        """ There is currently no checkpointing implemented"""
        pass

    def resume_from_checkpoint(self):
        """ There is currently no checkpointing implemented"""
        pass

    def finalize(self):
        logging.info("Writing samples to files")
        for fn in [self.checkpoint_file, self.backup_file]:
            self.write_results(fn)

    @property
    def model_stats(self):
        return {}

    @property
    def samples(self):
        samples = self.result['samples']
        params = list(self.model.variable_params)
        samples_dict = {p: samples[:, i] for i, p in enumerate(params)}
        return samples_dict

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
            fp.write_samples(self.samples, self.samples.keys())
            # write log evidence
            fp.write_logevidence(self.logz, self.logz_err)

            # write full results
            dump_state(self.result, fp,
                       path='sampler_info',
                       dsetname='presult')

    @property
    def logz(self):
        """Return bayesian evidence estimated by snowline sampler.
        """
        return self.result['logz']

    @property
    def logz_err(self):
        """Return error in bayesian evidence estimated by snowline sampler.
        """
        return self.result['logzerr']
