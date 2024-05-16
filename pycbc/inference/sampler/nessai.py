"""
This modules provides class for using the nessai sampler package for parameter
estimation.

Documentation for nessai: https://nessai.readthedocs.io/en/latest/
"""
import ast
import logging
import os

import nessai.flowsampler
import nessai.model
import nessai.livepoint
import nessai.utils.multiprocessing
import nessai.utils.settings
import numpy
import numpy.lib.recfunctions as rfn

from .base import BaseSampler, setup_output
from .base_mcmc import get_optional_arg_from_config
from ..io import NessaiFile, loadfile
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
        self.resume_data = None

    @property
    def io(self):
        return self._io

    @property
    def model_stats(self):
        pass

    @property
    def samples(self):
        """The raw nested samples including the corresponding weights"""
        if self._sampler.ns.nested_samples:
            ns = numpy.array(self._sampler.ns.nested_samples)
            samples = nessai.livepoint.live_points_to_dict(
                ns,
                self.model.sampling_params,
            )
            samples["logwt"] = self._sampler.ns.state.log_posterior_weights
            samples["loglikelihood"] = ns["logL"]
            samples["logprior"] = ns["logP"]
            samples["it"] = ns["it"]
        else:
            samples = {}
        return samples

    def run(self, **kwargs):
        """Run the sampler"""
        default_kwds, default_run_kwds = self.get_default_kwds(
            importance_nested_sampler=self.extra_kwds.get(
                "importance_nested_sampler", False
            )
        )

        extra_kwds = self.extra_kwds.copy()
        run_kwds = self.run_kwds.copy()

        # Output in kwargs takes priority of extra kwds.
        output = kwargs.pop("output", extra_kwds.pop("output", None))
        # If neither have been specified, use the path from the checkpoint file
        if output is None:
            output = os.path.join(
                os.path.dirname(os.path.abspath(self.checkpoint_file)),
                "outdir_nessai",
            )

        if kwargs is not None:
            logging.info("Updating keyword arguments with %s", kwargs)
            extra_kwds.update(
                {k: v for k, v in kwargs.items() if k in default_kwds}
            )
            run_kwds.update(
                {k: v for k, v in kwargs.items() if k in default_run_kwds}
            )

        if self._sampler is None:
            logging.info("Initialising nessai FlowSampler")
            self._sampler = nessai.flowsampler.FlowSampler(
                self.model_call,
                output=output,
                pool=self.pool,
                n_pool=self.nprocesses,
                close_pool=False,
                signal_handling=False,
                resume_data=self.resume_data,
                checkpoint_callback=self.checkpoint_callback,
                **extra_kwds,
            )
        logging.info("Starting sampling with nessai")
        self._sampler.run(**run_kwds)

    @staticmethod
    def get_default_kwds(importance_nested_sampler=False):
        """Return lists of all allowed keyword arguments for nessai.

        Returns
        -------
        default_kwds : list
            List of keyword arguments that can be passed to FlowSampler
        run_kwds: list
            List of keyword arguments that can be passed to FlowSampler.run
        """
        return nessai.utils.settings.get_all_kwargs(
            importance_nested_sampler=importance_nested_sampler,
            split_kwargs=True,
        )

    @classmethod
    def from_config(
        cls, cp, model, output_file=None, nprocesses=1, use_mpi=False
    ):
        """
        Loads the sampler from the given config file.
        """
        section = "sampler"
        # check name
        assert (
            cp.get(section, "name") == cls.name
        ), "name in section [sampler] must match mine"

        if cp.has_option(section, "importance_nested_sampler"):
            importance_nested_sampler = cp.get(
                section,
                "importance_nested_sampler",
            )
        else:
            importance_nested_sampler = False

        # Requires additional development work, see the model class below
        if importance_nested_sampler is True:
            raise NotImplementedError(
                "Importance nested sampler is not currently supported"
            )

        default_kwds, default_run_kwds = cls.get_default_kwds(
            importance_nested_sampler
        )

        # Keyword arguments the user cannot configure via the config
        remove_kwds = [
            "pool",
            "n_pool",
            "close_pool",
            "signal_handling",
            "checkpoint_callback",
        ]

        for kwd in remove_kwds:
            default_kwds.pop(kwd, None)
            default_run_kwds.pop(kwd, None)

        kwds = {}
        run_kwds = {}

        # ast.literal_eval is used here since specifying a dictionary with all
        # various types would be difficult. However, one may wish to revisit
        # this in future, e.g. if evaluating code is a concern.
        for d_out, d_defaults in zip(
            [kwds, run_kwds], [default_kwds, default_run_kwds]
        ):
            for k in d_defaults.keys():
                if cp.has_option(section, k):
                    option = cp.get(section, k)
                    try:
                        # This will fail for e.g. a string with an underscore
                        option = ast.literal_eval(option)
                    except ValueError:
                        pass
                    d_out[k] = option

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
        logging.info("nessai keyword arguments: %s", kwds)
        logging.info("nessai run keyword arguments: %s", run_kwds)

        loglikelihood_function = get_optional_arg_from_config(
            cp, section, "loglikelihood-function"
        )

        obj = cls(
            model,
            loglikelihood_function=loglikelihood_function,
            nprocesses=nprocesses,
            use_mpi=use_mpi,
            run_kwds=run_kwds,
            extra_kwds=kwds,
        )

        # Do not need to check number of samples for a nested sampler
        setup_output(obj, output_file, check_nsamples=False)
        if not obj.new_checkpoint:
            obj.resume_from_checkpoint()
        return obj

    def set_initial_conditions(
        self,
        initial_distribution=None,
        samples_file=None,
    ):
        """Sets up the starting point for the sampler.

        This is not used for nessai.
        """

    def checkpoint_callback(self, state):
        """Callback for checkpointing.

        This will be called periodically by nessai.
        """
        for fn in [self.checkpoint_file, self.backup_file]:
            with self.io(fn, "a") as fp:
                fp.write_pickled_data_into_checkpoint_file(state)
            self.write_results(fn)

    def checkpoint(self):
        """Checkpoint the sampler"""
        self.checkpoint_callback(self._sampler.ns)

    def resume_from_checkpoint(self):
        """Reads the resume data from the checkpoint file."""
        try:
            with loadfile(self.checkpoint_file, "r") as fp:
                self.resume_data = fp.read_pickled_data_from_checkpoint_file()
            logging.info(
                "Found valid checkpoint file: %s", self.checkpoint_file
            )
        except Exception as e:
            logging.info("Failed to load checkpoint file with error: %s", e)

    def finalize(self):
        """Finalize sampling"""
        logz = self._sampler.ns.log_evidence
        dlogz = self._sampler.ns.log_evidence_error
        logging.info("log Z, dlog Z: %s, %s", logz, dlogz)
        self.checkpoint()

    def write_results(self, filename):
        """Write the results to a given file.

        Writes the nested samples, log-evidence and log-evidence error.
        """
        with self.io(filename, "a") as fp:
            fp.write_raw_samples(self.samples)
            fp.write_logevidence(
                self._sampler.ns.log_evidence,
                self._sampler.ns.log_evidence_error,
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
            loglikelihood_function = "loglikelihood"
        self.loglikelihood_function = loglikelihood_function

        # Configure the priors bounds
        bounds = {}
        for dist in model.prior_distribution.distributions:
            bounds.update(
                **{
                    k: [v.min, v.max]
                    for k, v in dist.bounds.items()
                    if k in self.names
                }
            )
        self.bounds = bounds
        # Prior and likelihood are not vectorised
        self.vectorised_likelihood = False
        self.vectorised_prior = False
        # Use the pool for computing the prior
        self.parallelise_prior = True

    def to_dict(self, x):
        """Convert a nessai live point array to a dictionary"""
        return {n: x[n].item() for n in self.names}

    def to_live_points(self, x):
        """Convert to the structured arrays used by nessai"""
        # It is possible this could be made faster
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

    def from_unit_hypercube(self, x):
        """Map from the unit-hypercube to the prior."""
        # Needs to be implemented for importance nested sampler
        # This method is already available in pycbc but the inverse is not
        raise NotImplementedError

    def to_unit_hypercube(self, x):
        """Map to the unit-hypercube to the prior."""
        # Needs to be implemented for importance nested sampler
        raise NotImplementedError
