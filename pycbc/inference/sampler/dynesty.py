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
This modules provides classes and functions for using the dynesty sampler
packages for parameter estimation.
"""

from __future__ import absolute_import

import logging
import os
import dynesty
from pycbc.inference.io import (DynestyFile, validate_checkpoint_files)
from .base import BaseSampler



#
# =============================================================================
#
#                                   Samplers
#
# =============================================================================
#

class DynestySampler(BaseSampler):
    """This class is used to construct an Dynesty sampler from the dynesty
    package.

    Parameters
    ----------
    model : model
        A model from ``pycbc.inference.models``.
    nlive : int
        Number of live points to use in sampler.
    dlogz: float
        Tolerance limit to the value of logz (also a convergence criteria)
    pool : function with map, Optional
        A provider of a map function that allows a function call to be run
        over multiple sets of arguments and possibly maps them to
        cores/nodes/etc.
    """
    name = "dynesty"
    _io = DynestyFile

    def __init__(self, model, **kwargs):

        self.model = model
        self.nlive = nlive
        self.dlogz = dlogz
        self.names = model.sampling_params
        self.ndim = len(model.sampling_params)
        self._sampler = None
        self.checkpoint_file = None
        
    def log_likelihood(self,cube):
        params = {p: v for p, v in zip(self.model.variable_params, cube)}
        self.model.update(**params)
        return self.model.loglikelihood
    
    def prior_transform(self,cube):
        prior_dists = self.model.prior_distribution.distributions
        dist_dict = {}
        for dist in prior_dists:
            dist_dict.update({param: dist for param in dist.params})
        for i, param in enumerate(self.model.variable_params):
            cube[i] = dist_dict[param].cdfinv(param, cube[i])
        return cube
    
    def run(self, **kwargs):
        OUTDIR = os.path.dirname(os.path.abspath(self.checkpoint_file))
        if self._sampler is None:
            self._sampler=dynesty.NestedSampler(self.log_likelihood,
                                                self.prior_transform, self.ndim,
                                                nlive=self.nlive,
                                                dlogz=self.dlogz **kwargs)
        res = self._sampler.run_nested()

    @property
    def io(self):
        return self._io

    @property
    def niterations(self):
        return len(tuple(self.samples.values())[0])

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
        dlogz = float(cp.get(section, "dlogz"))
        obj = cls(model, nlive=nlive, dlogz=dlogz)
        return obj

    def checkpoint(self):
        pass

    def finalize(self):
        logz = self._sampler.results.logz[-1:][0]
        dlogz = self._sampler.results.logzerr[-1:][0]
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
        logl = self._sampler.results.logl
        return {'loglikelihood': logl}

    @property
    def samples(self):
        samples_dict = {p: self._sampler.results.samples[:,i] for p,i in
                        zip(self.model.sampling_params,self.ndim)}
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
    def posterior_samples(self):
        return self._sampler.results.samples

    @property
    def logz(self):
        return self._sampler.results.logz[-1:][0]

    @property
    def dlogz(self):
        return self._sampler.results.dlogz[-1:][0]

