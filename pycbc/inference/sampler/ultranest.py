# Copyright (C) 2020  Alex Nitz
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
This modules provides classes and functions for using the dynesty sampler
packages for parameter estimation.
"""


from __future__ import absolute_import

import logging

from pycbc.inference.io.ultranest import UltranestFile
from .base import (BaseSampler, setup_output)
from .base_cube import setup_calls


#
# =============================================================================
#
#                                   Samplers
#
# =============================================================================
#

class UltranestSampler(BaseSampler):
    """This class is used to construct an Ultranest sampler from the ultranest
    package.

    Parameters
    ----------
    model : model
        A model from ``pycbc.inference.models``.
    """
    name = "ultranest"
    _io = UltranestFile

    def __init__(self, model, **kwargs):
        super(UltranestSampler, self).__init__(model)

        import ultranest
        log_likelihood_call, prior_call = setup_calls(model, copy_prior=True)

        self._sampler = ultranest.ReactiveNestedSampler(
            list(self.model.variable_params),
            log_likelihood_call,
            prior_call)

        self.nlive = 0
        self.ndim = len(self.model.variable_params)
        self.result = None
        self.kwargs = kwargs  # Keywords for the run method of ultranest

    def run(self):
        self.result = self._sampler.run(**self.kwargs)

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
        opts = {'update_interval_iter_fraction': float,
                'update_interval_ncall': int,
                'log_interval': int,
                'show_status': bool,
                'dlogz': float,
                'dKL': float,
                'frac_remain': float,
                'Lepsilon': float,
                'min_ess': int,
                'max_iters': int,
                'max_ncalls': int,
                'max_num_improvement_loops': int,
                'min_num_live_points': int,
                'cluster_num_live_points:': int}
        for opt_name in opts:
            if cp.has_option('sampler', opt_name):
                value = cp.get('sampler', opt_name)
                skeys[opt_name] = opts[opt_name](value)
        inst = cls(model, **skeys)
        setup_output(inst, output_file)
        return inst

    def checkpoint(self):
        pass

    def resume_from_checkpoint(self):
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
        samples_dict = {p: self.result['samples'][:, i] for p, i in
                        zip(self.model.variable_params, range(self.ndim))}
        return samples_dict

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
            fp.write_samples(self.samples, self.model.variable_params)
            # write log evidence
            fp.write_logevidence(self.logz, self.logz_err)

    @property
    def logz(self):
        """
        return bayesian evidence estimated by
        dynesty sampler
        """
        return self.result['logz']

    @property
    def logz_err(self):
        """
        return error in bayesian evidence estimated by
        dynesty sampler
        """
        return self.result['logzerr']
