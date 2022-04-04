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


import logging
import os
import array
import cpnest
import cpnest.model as cpm
from pycbc.inference.io import (CPNestFile, validate_checkpoint_files)
from .base import (BaseSampler, setup_output)
from .base_mcmc import get_optional_arg_from_config



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

    def __init__(self, model, nlive, maxmcmc=1000, nthreads=1, verbose=1,
                 loglikelihood_function=None):
        self.model = model
        self.nlive = nlive
        self.maxmcmc = maxmcmc
        self.nthreads = nthreads
        self.verbose = verbose
        # create a wrapper for calling the model
        self.model_call = CPNestModel(model, loglikelihood_function)
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
        return len(tuple(self.samples.values())[0])

    @classmethod
    def from_config(cls, cp, model, output_file=None, nprocesses=1,
                    use_mpi=False):
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
        loglikelihood_function = \
            get_optional_arg_from_config(cp, section, 'loglikelihood-function')
        obj = cls(model, nlive=nlive, maxmcmc=maxmcmc, nthreads=nthreads,
                  verbose=verbose,
                  loglikelihood_function=loglikelihood_function)

        setup_output(obj, output_file, check_nsamples=False)
        if not obj.new_checkpoint:
            obj.resume_from_checkpoint()
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
            self.checkpoint_file, self.backup_file, check_nsamples=False)
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

    def resume_from_checkpoint(self):
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


class CPNestModel(cpm.Model):
    """
    Class for making PyCBC Inference 'model class'
    compatible with CPNest 'model class'

    Parameters
    ----------
    model : inference.BaseModel instance
             A model instance from pycbc.
    """
    def __init__(self, model, loglikelihood_function=None):
        if model.sampling_transforms is not None:
            raise ValueError("CPNest does not support sampling transforms")
        self.model = model
        self.names = list(model.sampling_params)
        # set up lohlikelihood_function
        if loglikelihood_function is None:
            loglikelihood_function = 'loglikelihood'
        self.loglikelihood_function = loglikelihood_function
        bounds = {}
        for dist in model.prior_distribution.distributions:
            bounds.update(dist.bounds)
        self.bounds = [bounds[params] for params in self.names]

    def new_point(self):
        point = self.model.prior_rvs()
        return cpm.LivePoint(list(self.model.sampling_params),
                             array.array('d', [point[p] for p in self.model.sampling_params]))

    def log_prior(self,xx):
        self.model.update(**xx)
        return self.model.logprior

    def log_likelihood(self, xx):
        """
        Modify the log likelihood which will be passed to CPNest 'model class'
        """
        self.model.update(**xx)
        return getattr(self.model, self.loglikelihood_function)
