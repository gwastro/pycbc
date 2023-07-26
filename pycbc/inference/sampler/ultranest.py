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
This modules provides classes and functions for using the ultranest sampler
packages for parameter estimation.
"""

import sys
import logging
import numpy

from pycbc.inference.io.ultranest import UltranestFile
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

class UltranestSampler(BaseSampler):
    """This class is used to construct an Ultranest sampler from the ultranest
    package.

    Parameters
    ----------
    model : model
        A model from ``pycbc.inference.models``.
    log_dir : str
        Folder where files should be stored for resuming (optional).
    stepsampling : bool
        If false, uses rejection sampling. If true, uses
        hit-and-run sampler, which scales better with dimensionality.
    """
    name = "ultranest"
    _io = UltranestFile

    def __init__(self, model, log_dir=None,
                 stepsampling=False,
                 enable_plots=False,
                 **kwargs):
        super(UltranestSampler, self).__init__(model)

        import ultranest
        log_likelihood_call, prior_call = setup_calls(model, copy_prior=True)

        # Check for cyclic boundaries
        periodic = []
        cyclic = self.model.prior_distribution.cyclic
        for param in self.variable_params:
            if param in cyclic:
                logging.info('Param: %s will be cyclic', param)
                periodic.append(True)
            else:
                periodic.append(False)

        self._sampler = ultranest.ReactiveNestedSampler(
            list(self.model.variable_params),
            log_likelihood_call,
            prior_call, log_dir=log_dir,
            wrapped_params=periodic,
            resume=True)

        if stepsampling:
            import ultranest.stepsampler
            self._sampler.stepsampler = ultranest.stepsampler.RegionBallSliceSampler(
                nsteps=100, adaptive_nsteps='move-distance',
                region_filter=True)

        self.enable_plots = enable_plots
        self.nlive = 0
        self.ndim = len(self.model.variable_params)
        self.result = None
        self.kwargs = kwargs  # Keywords for the run method of ultranest

        do_mpi, _, rank = use_mpi()
        self.main = (not do_mpi) or (rank == 0)

    def run(self):
        self.result = self._sampler.run(**self.kwargs)
        if not self.main:
            sys.exit(0)
        self._sampler.print_results()

        if self.enable_plots:
            self._sampler.plot()

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
                'log_dir': str,
                'stepsampling': bool,
                'enable_plots': bool,
                'max_num_improvement_loops': int,
                'min_num_live_points': int,
                'cluster_num_live_points:': int}
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
        from ultranest.utils import resample_equal

        # we'll do the resampling ourselves so we can pick up
        # additional parameters
        try:  # Remove me on next ultranest release
            wsamples = self.result['weighted_samples']['v']
            weights = self.result['weighted_samples']['w']
            logl = self.result['weighted_samples']['L']
        except KeyError:
            wsamples = self.result['weighted_samples']['points']
            weights = self.result['weighted_samples']['weights']
            logl = self.result['weighted_samples']['logl']

        wsamples = numpy.column_stack((wsamples, logl))
        params = list(self.model.variable_params) + ['loglikelihood']
        samples = resample_equal(wsamples, weights / weights.sum())
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

            # write full ultranest formatted results
            dump_state(self.result, fp,
                       path='sampler_info',
                       dsetname='presult')

    @property
    def logz(self):
        """Return bayesian evidence estimated by ultranest sampler.
        """
        return self.result['logz']

    @property
    def logz_err(self):
        """Return error in bayesian evidence estimated by ultranest sampler.
        """
        return self.result['logzerr']
