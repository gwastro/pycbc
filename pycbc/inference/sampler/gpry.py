# Copyright (C) 2025 Jahed Abedi
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
This modules provides classes and functions for using the Fast and robust Bayesian inference
using Gaussian processes with GPry sampler packages for parameter estimation.
"""

from pycbc.distributions.uniform import Uniform
#from .base import BaseSampler  # Import from local base.py
from .base import (BaseSampler, setup_output)
from pycbc.inference.io.base_hdf import BaseInferenceFile  # Import from base.py
from pycbc.inference.io import loadfile
import gpry  # Ensure GPry is installed
from gpry.run import Runner # Correct import path
import os
import inspect
import numpy
import h5py
import datetime
import logging

logger = logging.getLogger(__name__)
FLOOR = -numpy.finfo(float).max  # a finite floor, e.g. -1.797693e+308

class GPrySampler(BaseSampler):
    name = "gpry"  # Unique identifier
    
    def __init__(self, model, gp_config=None):
        super().__init__(model)
        self.gp_config = gp_config or {}

        # Debug: Print model details
        print("\n===== Initializing GPrySampler =====")
        print("Model type:", type(model).__name__)
        print("Model variable parameters:", model.variable_params)
        print("Model prior distributions:", model.prior_distribution)
        print("===============================\n")

        # Do NOT initialize gp_runner here
        self._log_likelihood = None
        self._bounds_dict = None

    def run(self):
        """Run the GPry sampler"""
        # Initialize GPry components here (after output_file is set)
        self._log_likelihood, self._bounds_dict = self._to_gpry_model()
        
        # Ensure bounds are ordered consistently with parameters
        param_names = self.model.variable_params
        bounds = [self._bounds_dict[p] for p in param_names]


        #self.gp_runner = Runner(
        #    model=self._log_likelihood,
        #    bounds=self._bounds_dict,
        #    checkpoint=self.output_file,  # Now accessible
        #    load_checkpoint="overwrite",
        #    options=self.gp_config.get("options", {})
        #)
        # Initialize GPry Runner with ordered bounds
        #self.gp_runner = Runner(
        #    model=self._log_likelihood,
        #    bounds=bounds,  # Pass as ordered list, not dict
        #    checkpoint=self.checkpoint_dir,
        #    load_checkpoint="overwrite",
        #    options=self.gp_config
        #)

        self.gp_runner = Runner(
        loglike=self._log_likelihood,      # <-- use loglike, not model
        bounds=bounds,                     # same as before
        checkpoint=self.checkpoint_dir,    # where to store checkpoints
        load_checkpoint="overwrite",       # resume behavior
        options=self.gp_config             # BO + MC settings bundled here
        )

        #if not hasattr(self.gp_runner, "has_converged"):
        #    raise AttributeError("GPrySampler: GPry's Runner object does not have 'has_converged'.")


        # Proceed with GPry's active learning loop
        #while not self.gp_runner.has_converged:
        #    print("GPry Runner callable methods:", [method for method in dir(self.gp_runner) if callable(getattr(self.gp_runner, method))])
        #    self.gp_runner.run()
        #    self.checkpoint()
        self.gp_runner.run()
        self.checkpoint()

    @classmethod
    def from_config(cls, cp, model, output_file=None, **kwargs):
        mc_opts = {
            "Rminus1_stop": cp.getfloat("mcmc", "Rminus1_stop", fallback=0.05),
            "max_tries":    cp.getint(  "mcmc", "max_tries",    fallback=1000)
        }
        gp_config = {
            "n_initial": cp.getint("gpry", "n_initial", fallback=1000),
            "max_initial": cp.getint("gpry", "max_initial", fallback=1000000),  # Increased
            "acquisition_function": cp.get("gpry", "gp_acquisition", fallback="LogExp"),
            "checkpoint_interval": cp.getint("gpry", "checkpoint_interval", fallback=100),
            "fit_full_every": cp.getint("gpry", "fit_full_every", fallback=10),
            "account_for_inf": None,  # Allow GPry to continue even with infinite likelihoods
            "preprocessing_y": None,  # Disable GPry’s normalization
            "allow_extreme_values": True  # New parameter to avoid filtering
        }
        obj = cls(model, gp_config=gp_config)

        # Ensure output_file is a valid HDF5 file
        if output_file is None:
            output_file = "inference.hdf"

        if not output_file.endswith(".hdf") and not output_file.endswith(".h5"):
            output_file += ".hdf"

        # Ensure PyCBC writes to a file but GPry gets a checkpoint directory
        obj.output_file = output_file
        #obj.checkpoint_dir = output_file + "_checkpoint"  # Separate folder for GPry
        # Ensure GPry gets a proper directory for checkpointing
        if not os.path.exists(output_file + "_checkpoint"):
            os.makedirs(output_file + "_checkpoint")

        obj.checkpoint_dir = output_file + "_checkpoint"

        # Remove mistaken folder
        if os.path.exists(output_file) and os.path.isdir(output_file):
            import shutil
            shutil.rmtree(output_file)

        setup_output(obj, output_file, check_nsamples=False)
        
        # Ensure the output file is properly formatted
        with h5py.File(output_file, "a") as f:
            if "filetype" not in f.attrs:
                print("Warning: 'filetype' attribute missing, adding manually.")
                f.attrs["filetype"] = "gpry"

        return obj


    def _to_gpry_model(self):
        print("DEBUG: self.model =", self.model, "Type:", type(self.model))  # Debugging line
        param_names = self.model.variable_params
        bounds_dict = {}

        # Extract priors for each parameter
        prior_map = {param: dist for dist in self.model.prior_distribution.distributions for param in dist.params}

        # Validate and extract bounds
        for param in param_names:
            prior = prior_map.get(param)
            if not isinstance(prior, Uniform):
                raise ValueError(f"Prior for {param} must be uniform.")
            bounds_dict[param] = [prior._bounds[param].min, prior._bounds[param].max]
        
        print("GPry Model Type:", type(self.model))
        print("GPry Model Attributes:", dir(self.model))
        # Define the log-likelihood function with correct method
        #def log_likelihood(*args, **kwargs):
        #    """Wrapper that explicitly lists all parameters as keyword arguments."""
        #    params = {param: kwargs[param] for param in param_names}
        #    try:
        #        # DEBUG: Print parameters being evaluated
        #        #print(f"\nEvaluating parameters: {params}")
        #
        #        # Set current parameters in the model
        #        #print('DEBUG: self.model=', dir(self.model))
        #        #print("self.model callable methods:", [method for method in dir(self.model) if callable(getattr(self.model, method))])
        #        #self.model._current_params = params
        #        self.model.update(**params)
        #        #self.model.update()  # No arguments; uses current_params
        #        #print('DEBUG: params=',params)
        #        #print('DEBUG: self.model._loglikelihood=',dir(self.model._loglikelihood))
        #        #print('DEBUG: self.model.loglr=',dir(self.model.loglr))
        #        #print('DEBUG: self.model.loglikelihood=', dir(self.model.loglikelihood))
        #        # Calculate log-likelihood
        #        lnl = self.model.loglikelihood  # No arguments; uses current_params
        #        lnl = lnl/100.0
        #        #lnl = self.model.loglr  # No arguments; uses current_params  
        #        # DEBUG: Print result
        #        print(f"Log-Likelihood: {lnl}")
        #        if not numpy.isfinite(lnl):
        #            print(f"Non-finite likelihood at: {params}")
        #        return lnl
        #    except Exception as e:
        #        print(f"Error evaluating {params}: {str(e)}")
        #        return -numpy.inf

        def log_likelihood(*args, **kwargs):
            params = {p: kwargs[p] for p in param_names}
            try:
                # Update model and compute raw log‑likelihood
                self.model.update(**params)
                lnl = self.model.loglikelihood / 100.0
            except Exception as e:
                # On any error, log a warning and return the finite FLOOR
                print(f"Warning: exception for {params}: {e}. Returning FLOOR.")
                return FLOOR

            # Convert to array for vector inputs, then clamp
            lnl = numpy.array(lnl, copy=False)

            # Identify any non-finite entries
            if not numpy.all(numpy.isfinite(lnl)):
                print(f"Warning: non-finite LNL {lnl} at {params}. Clamping to FLOOR.")
                # Clip everything below FLOOR to FLOOR; leave upper bound infinite
                lnl = numpy.clip(lnl, FLOOR, None)  # :contentReference[oaicite:1]{index=1}

            return lnl

        # Dynamically set the function signature to include parameters
        params_str = ", ".join(param_names)
        func_def = f"def _log_likelihood_wrapper({params_str}):\n" \
                   "    return log_likelihood(**locals())"
        namespace = {'log_likelihood': log_likelihood}
        exec(func_def, namespace)
        log_likelihood_wrapper = namespace['_log_likelihood_wrapper']

        # Ensure the wrapper has the correct __annotations__ and __kwdefaults__
        log_likelihood_wrapper.__annotations__ = {p: float for p in param_names}
        log_likelihood_wrapper.__kwdefaults__ = {p: None for p in param_names}

        return log_likelihood_wrapper, bounds_dict


    def set_initial_conditions(self):
        # GPry generates initial points internally; override if needed
        pass

    
    #def checkpoint(self):
    #    """Checkpoint the sampler"""
    #    with self.io(self.checkpoint_file, "a") as fp:
    #        #fp.write_checkpoint(self.gp_runner.state)
    #        self.gp_runner.save_checkpoint()
    #        #fp.write_checkpoint(self.gp_runner.random_state)

    def checkpoint(self):
        """Checkpoint the sampler"""
        if self.gp_runner.checkpoint is None:
            self.gp_runner.checkpoint = self.checkpoint_dir  # Ensure checkpoint is set
        print("Saving GPry checkpoint at:", self.gp_runner.checkpoint)
        self.gp_runner.save_checkpoint()

    def resume_from_checkpoint(self):
        """Resume the sampler from a checkpoint file."""
        # Load the last checkpoint
        #checkpoint_file = self.io.get_checkpoint_file()
        checkpoint_file = self.checkpoint_dir
        if not os.path.exists(checkpoint_file):
            raise ValueError("No checkpoint directory found to resume from.")
        self.gp_runner.read_checkpoint(checkpoint_file)
        if checkpoint_file is None:
            raise ValueError("No checkpoint file found to resume from.")

        # Load the state from the checkpoint
        with self.io(checkpoint_file, "r") as fp:
            state = fp.read_sampler_state()

        # Restore the sampler's state (e.g., GP model, training points)
        self.gp_runner = state["gp_runner"]
        self.iteration = state["iteration"]

    #@property
    #def samples(self):
    #    # Return samples from the GP posterior (e.g., mean/covariance)
    #    return self.gp_runner.get_posterior_samples()


    @property
    #def samples(self):
    #    # If no MC sample has been generated yet, do it now
    #    if not hasattr(self.gp_runner, '_last_mc_samples'):
    #        # delegate to GPry to draw its MC posterior samples
    #        self.gp_runner.generate_mc_sample()
    #    # return the raw dict of arrays that BaseInferenceFile knows how to write
    #    return self.gp_runner.last_mc_samples(copy=False)
    def samples(self):
        # Only run MCMC once; thereafter reuse the stored samples
        if not hasattr(self.gp_runner, '_last_mc_samples'):
            # Pull MCMC options (e.g. Rminus1_stop, max_tries) from gp_config
            mc_opts = self.gp_config.get("mc_options", {})  
            # Delegate to GPry to draw MC posterior samples with your custom settings
            self.gp_runner.generate_mc_sample(add_options=mc_opts)  # :contentReference[oaicite:0]{index=0}

        # Return the dict of arrays for BaseInferenceFile.write_samples()
        return self.gp_runner.last_mc_samples(copy=False)  # :contentReference[oaicite:1]{index=1}

    @property
    def model_stats(self):
        # Retrieve log likelihood and prior from the GP model
        return {
            "log_likelihood": self.gp_runner.log_likelihoods,
            "log_prior": self.gp_runner.log_priors
        }

    @property
    def io(self):
        from pycbc.inference.io.gpry import GPryIO
        return GPryIO


    #def finalize(self):
    #    """Finalize the sampler and ensure PyCBC writes the run end time properly."""
    #    print("Finalizing GPrySampler: Ensuring proper HDF5 output.")
    #
    #    # Use PyCBC's internal method for writing the run end time
    #    with self.io(self.output_file, "a") as fp:
    #        if hasattr(fp, "write_run_end_time"):
    #            fp.write_run_end_time()
    #        else:
    #            print("Warning: 'write_run_end_time' method not found in PyCBC IO.")

    #def finalize(self):
    #    """Finalize the sampler and ensure PyCBC writes the run end time properly."""
    #    print("Finalizing GPrySampler: Ensuring proper HDF5 output.")
    #    with self.io(self.output_file, "a") as fp:
    #        try:
    #            # Try to retrieve the method; if not present, KeyError is caught
    #            write_end = getattr(fp, "write_run_end_time")
    #        except KeyError:
    #            write_end = None
    #
    #        if write_end is not None:
    #            write_end()
    #        else:
    #            print("Warning: 'write_run_end_time' method not found in PyCBC IO.")


    def finalize(self):
        """Finalize the sampler: write end time *and* the posterior."""
        print("Finalizing GPrySampler: Ensuring proper HDF5 output.")

        # -- existing run_end_time logic --
        with self.io(self.output_file, "a") as fp:
            write_end = getattr(fp, "write_run_end_time", None)
            if callable(write_end):
                write_end()
            else:
                print("Warning: 'write_run_end_time' method not found.")

        # -- NEW: write out the samples! --
        with self.io(self.output_file, "a") as fp:
            fp.write_samples(self.samples)
        print(f"Wrote {len(self.samples.get(next(iter(self.samples))))} samples to 'samples' group.")

