# Copyright (C) 2019  Collin Capano, Sumit Kumar
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
This modules provides classes and functions for using the cpnest sampler
packages for parameter estimation.
"""

from __future__ import absolute_import

import logging
import os
import cpnest
import cpnest.model as cpm
from pycbc.inference.io import (CPNestFile, validate_checkpoint_files)
from .base import BaseSampler



#
# =============================================================================
#
#                                   Samplers
#
# =============================================================================
#

class CPNestSampler(BaseSampler):
    """This class is used to construct an CPNest sampler from the cpnest
    package by John Veitch.

    Parameters
    ----------
    model : model
        A model from ``pycbc.inference.models``.
    nlive : int
        Number of live points to use in sampler.
    pool : function with map, Optional
        A provider of a map function that allows a function call to be run
        over multiple sets of arguments and possibly maps them to
        cores/nodes/etc.
    """
    name = "cpnest"
    _io = CPNestFile

    def __init__(self, model, nlive, verbose, nthreads=8, maxmcmc=1000):

        self.model = model
        self.nlive = nlive
        self.maxmcmc = maxmcmc
        self.nthreads = nthreads
        self.verbose = verbose
        # create a wrapper for calling the model
        self.model_call = CPnestModel(model)
        self._sampler = None
        self._nested_samples = None
        self._posterior_samples = None
        self._logz = None
        self._dlogz = None
        self.checkpoint_file = None

    def run(self):
        out_dir = os.path.dirname(os.path.abspath(self.checkpoint_file))
        if self._sampler is None:
            self._sampler = cpnest.CPNest(self.model_call, verbose=1,
                                          output=out_dir,
                                          nthreads=self.nthreads,
                                          nlive=self.nlive,
                                          maxmcmc=self.maxmcmc, resume=True)
        res = self._sampler.run()

    @property
    def io(self):
        return self._io

    @property
    def niterations(self):
        return len(self.samples.values()[0])

    @classmethod
    def from_config(cls, cp, model, nprocesses=1, use_mpi=False):
        """
        Loads the sampler from the given config file.
        """
        section = "sampler"
        # check name
        assert cp.get(section, "name") == cls.name, (
            "name in section [sampler] must match mine")
        # get the number of live points to use
        nlive = int(cp.get(section, "nlive"))
        maxmcmc = int(cp.get(section, "maxmcmc"))
        nthreads = int(cp.get(section, "nthreads"))
        verbose = int(cp.get(section, "verbose"))
        obj = cls(model, nlive=nlive, maxmcmc=maxmcmc,
                  nthreads=nthreads, verbose=verbose)
        return obj

    def checkpoint(self):
        pass

    def finalize(self):
        logz = self._sampler.NS.logZ
        dlogz = 0.1  #######FIXME!!!!!###############
        logging.info("log Z, dlog Z: {}, {}".format(logz, dlogz))
        for fn in [self.checkpoint_file]:
            with self.io(fn, "a") as fp:
                fp.write_logevidence(logz, dlogz)
        logging.info("Writing samples to files")
        for fn in [self.checkpoint_file, self.backup_file]:
            self.write_results(fn)
        logging.info("Validating checkpoint and backup files")
        checkpoint_valid = validate_checkpoint_files(
            self.checkpoint_file, self.backup_file)
        if not checkpoint_valid:
            raise IOError("error writing to checkpoint file")

    @property
    def model_stats(self):
        logl = self._sampler.posterior_samples['logL']
        logp = self._sampler.posterior_samples['logPrior']
        return {'loglikelihood': logl, 'logprior': logp}

    @property
    def samples(self):
        samples_dict = {p: self._sampler.posterior_samples[p] for p in
                        self.posterior_samples.dtype.names}
        return samples_dict

    def set_initial_conditions(self, initial_distribution=None,
                               samples_file=None):
        """Sets up the starting point for the sampler.

        Should also set the sampler's random state.
        """
        pass

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
            # write stats
            fp.write_samples(self.model_stats)
            # write log evidence
            fp.write_logevidence(self._sampler.NS.logZ, 0.1)

    @property
    def nested_samples(self):
        return self._sampler.nested_samples

    @property
    def posterior_samples(self):
        return self._sampler.posterior_samples

    @property
    def logz(self):
        return self._logz

    @property
    def dlogz(self):
        return self._dlogz


class CPnestModel(cpm.Model):
    """
    Class for making PyCBC Inference 'model class'
    compatible with CPNest 'model class'

    Parameters
    ----------
    model : inference.BaseModel instance
             A model instance from pycbc.
    """
    def __init__(self, model):
        if model.sampling_transforms is not None:
            raise ValueError("CPNest does not support sampling transforms")
        self.model = model
        self.names = list(model.sampling_params)
        bounds = {}
        for dist in model.prior_distribution.distributions:
            bounds.update(dist.bounds)
        self.bounds = [bounds[params] for params in self.names]

    def log_likelihood(self, xx):
        """
        Modify the log likelihood which will be passed to CPNest 'model class'
        """
        self.model.update(**xx)
        return self.model.loglikelihood
