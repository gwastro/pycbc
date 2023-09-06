"""
This modules provides class for using the nessai sampler package for parameter
estimation.
"""
import os

import nessai.flowsampler
import nessai.model
import nessai.livepoint
import numpy.lib.recfunctions as rfn

from .base import BaseSampler, setup_output
from .base_mcmc import get_optional_arg_from_config
from ..io import NessaiFile


class NessaiSampler(BaseSampler):
    """Class to construct a FlowSampler from the nessai package."""

    name = "nessai"
    _io = NessaiFile

    def __init__(self, model, nlive, loglikelihood_function, **kwargs):
        super().__init__(model)

        # TODO: add other options
        self.nlive = nlive
        self.model_call = NessaiModel(self.model, loglikelihood_function)

        #TODO: handle multiprocessing

        self._sampler = None
        self._nested_samples = None
        self._posterior_samples = None
        self._logz = None
        self._dlogz = None
        self.checkpoint_file = None

    @property
    def io(self):
        return self._io

    @property
    def model_stats(self):
        return {
            "loglikelihood": self._sampler.posterior_samples["logL"],
            "logprior": self._sampler.posterior_samples["logP"],
        }
    
    @property
    def samples(self):
        return nessai.livepoint.live_points_to_dict(
            self._sampler.posterior_samples,
            self.model.sampling_params,
        )

    def run(self, resume=False):
        out_dir = os.path.dirname(os.path.abspath(self.checkpoint_file))
        if self._sampler is None:
            self._sampler = nessai.flowsampler.FlowSampler(
                self.model_call,
                output=out_dir,
                nlive=self.nlive,
                resume=resume,
            )
        self._sampler.run()

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
        # TODO: add other options
        nlive = int(cp.get(section, "nlive"))

        loglikelihood_function = \
            get_optional_arg_from_config(cp, section, 'loglikelihood-function')
        obj = cls(
            model,
            nlive=nlive,
            loglikelihood_function=loglikelihood_function,
        )

        setup_output(obj, output_file, check_nsamples=False)
        if not obj.new_checkpoint:
            obj.resume_from_checkpoint()
        return obj

    def set_initial_conditions(self, initial_distribution=None,
                               samples_file=None):
        """Sets up the starting point for the sampler.

        Should also set the sampler's random state.
        """
        pass

    def checkpoint(self):
        self._sampler.ns.checkpoint()

    def resume_from_checkpoint(self):
        # TODO: check this works
        self.run(resume=True)

    def finalize(self):
        logz = self._sampler.log_evidence
        dlogz = self._sampler.log_evidence_error

        for fn in [self.checkpoint_file]:
            with self.io(fn, "a") as fp:
                fp.write_logevidence(logz, dlogz)
        for fn in [self.checkpoint_file, self.backup_file]:
            self.write_results(fn)

    def write_results(self, filename):
        with self.io(filename, "a") as fp:
            fp.write_samples(self.samples, self.model.sampling_params)
            fp.write_samples(self.model_stats)
            fp.write_logevidence(
                self._sampler.log_evidence,
                self._sampler.log_evidence_error,
            )


class NessaiModel(nessai.model.Model):
    """Wrapper for PyCBC Inference model class for use with nessai.
    
    Parameters
    ----------
    model : inference.BaseModel instance
        A model instance from PyCBC.
    loglikelihood_function : str
        Name of the log-likelihood method to call.
    """
    def __init__(self, model, loglikelihood_function=None):
        self.model = model
        self.names = list(model.sampling_params)

        # Configure the log-likelihood function
        if loglikelihood_function is None:
            loglikelihood_function = 'loglikelihood'
        self.loglikelihood_function = loglikelihood_function

        # Configure the priors bounds
        bounds = {}
        for dist in model.prior_distribution.distributions:
            bounds.update(**{k: [v.min, v.max] for k, v in dist.bounds.items()})
        self.bounds = bounds

        # Prior and likelihood are not vectorised
        self.vectorised_likelihood = False
        self.vectorised_prior = False

    def to_dict(self, x):
        return {n: x[n].item() for n in self.names}
    
    def to_live_points(self, x):
        """Convert to the structured arrays used by nessai"""
        # TODO: could this be improved?
        return nessai.livepoint.numpy_array_to_live_points(
            rfn.structured_to_unstructured(x),
            self.names,
        )

    def new_point(self, N=1):
        """Draw a new point"""
        return self.to_live_points(self.model.prior_rvs(size=N))

    def new_point_log_prob(self, x):
        """Log-probability for the ``new_point`` method"""
        return self.batch_evaluate_log_prior(x)

    def log_prior(self, x):
        """Compute the log-prior"""
        self.model.update(**self.to_dict(x))
        return self.model.logprior

    def log_likelihood(self, x):
        """Compute the log-likelihood"""
        self.model.update(**self.to_dict(x))
        return getattr(self.model, self.loglikelihood_function)
