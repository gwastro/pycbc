"""
This modules provides class for using the nessai sampler package for parameter
estimation.
"""
import ast
import logging
import os

import nessai.flowsampler
import nessai.model
import nessai.livepoint
import nessai.utils.multiprocessing
import nessai.utils.settings
import numpy.lib.recfunctions as rfn

from .base import BaseSampler, setup_output
from .base_mcmc import get_optional_arg_from_config
from ..io import NessaiFile
from ...pool import choose_pool


class NessaiSampler(BaseSampler):
    """Class to construct a FlowSampler from the nessai package."""

    name = "nessai"
    _io = NessaiFile

    def __init__(
        self,
        model,
        loglikelihood_function,
        nlive=1000,
        nprocesses=1,
        use_mpi=False,
        run_kwds=None,
        extra_kwds=None,
    ):
        super().__init__(model)

        self.nlive = nlive
        self.model_call = NessaiModel(self.model, loglikelihood_function)

        self.extra_kwds = extra_kwds if extra_kwds is not None else {}
        self.run_kwds = run_kwds if run_kwds is not None else {}

        nessai.utils.multiprocessing.initialise_pool_variables(self.model_call)
        self.pool = choose_pool(mpi=use_mpi, processes=nprocesses) 
        self.nprocesses = nprocesses

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
        """The raw nested samples including the corresponding weights"""
        samples = nessai.livepoint.live_points_to_dict(
            self._sampler.nested_samples,
            self.model.sampling_params,
        )
        samples["logwt"] = self._sampler.ns.state.log_posterior_weights
        samples["loglikelihood"] = self._sampler.nested_samples["loglikelihood"]
        return samples

    def run(self):
        out_dir = os.path.join(
            os.path.dirname(os.path.abspath(self.checkpoint_file)),
            "nessai",
        )
        if self._sampler is None:
            self._sampler = nessai.flowsampler.FlowSampler(
                self.model_call,
                output=out_dir,
                pool=self.pool,
                n_pool=self.nprocesses,
                close_pool=False,
                signal_handling=False,
                **self.extra_kwds,
            )
        self._sampler.run(**self.run_kwds)

    @classmethod
    def from_config(
        cls, cp, model, output_file=None, nprocesses=1, use_mpi=False
    ):
        """
        Loads the sampler from the given config file.
        """
        section = "sampler"
        # check name
        assert cp.get(section, "name") == cls.name, (
            "name in section [sampler] must match mine")

        if cp.has_option(section, "importance_nested_sampler"):
            importance_nested_sampler = cp.get(
                section, "importance_nested_sampler",
            )
        else:
            importance_nested_sampler = False

        # Determine all possible keyword arguments that are not hardcoded
        default_kwds, default_run_kwds = nessai.utils.settings.get_all_kwargs(
            importance_nested_sampler=importance_nested_sampler,
            split_kwargs=True,
        )

        # Keyword arguments the user cannot configure via the config
        remove_kwds = [
            "pool", "n_pool", "close_pool", "signal_handling"
        ]

        for kwd in remove_kwds:
            default_kwds.pop(kwd, None)
            default_run_kwds.pop(kwd, None)

        kwds = {}
        run_kwds = {}

        for d_out, d_defaults in zip(
            [kwds, run_kwds], [default_kwds, default_run_kwds]
        ):
            for k in d_defaults.keys():
                if cp.has_option(section, k):
                    d_out[k] = ast.literal_eval(cp.get(section, k))

        # Specified kwds
        ignore_kwds = {"nlive", "name"}
        invalid_kwds = (
            cp[section].keys()
            - set().union(kwds.keys(), run_kwds.keys())
            - ignore_kwds
        )

        if invalid_kwds:
            raise RuntimeError(
                f"Config contains unknown options: {invalid_kwds}"
            )
        logging.info(f"nessai keyword arguments: {kwds}")
        logging.info(f"nessai run keyword arguments: {run_kwds}")

        loglikelihood_function = \
            get_optional_arg_from_config(cp, section, 'loglikelihood-function')

        obj = cls(
            model,
            loglikelihood_function=loglikelihood_function,
            nprocesses=nprocesses,
            use_mpi=use_mpi,
            run_kwds=run_kwds,
            extra_kwds=kwds,
        )

        setup_output(obj, output_file, check_nsamples=False, validate=False)
        return obj

    def set_initial_conditions(
        self,
        initial_distribution=None,
        samples_file=None,
    ):
        """Sets up the starting point for the sampler.

        Should also set the sampler's random state.
        """
        pass

    def checkpoint(self):
        self._sampler.ns.checkpoint()
        for fn in [self.checkpoint_file, self.backup_file]:
            self.write_results(fn)

    def resume_from_checkpoint(self):
        pass

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
            bounds.update(
                **{k: [v.min, v.max] for k, v in dist.bounds.items() if k in self.names}
            )
        self.bounds = bounds
        # Prior and likelihood are not vectorised
        self.vectorised_likelihood = False
        self.vectorised_prior = False
        # Use the pool for computing the prior
        self.parallelise_prior = True

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
