from pycbc.inference.io.base_hdf import BaseInferenceFile
from pycbc.io.record import FieldArray
import h5py
import numpy as np
import time
import datetime

class GPryIO(BaseInferenceFile):
    name = "gpry"
    _file_extensions = ["hdf", "h5"]

    def write_resume_point(self):
        """Write minimal data needed to resume a run."""
        # GPry handles checkpoints internally, so this can be a no-op
        pass

    def write_run_start_time(self):
        """Write the start time of the run to the file."""
        self.attrs["run_start_time"] = time.time()

    def write_sampler_metadata(self, sampler):
        """Write metadata specific to GPry sampler."""
        self.attrs["sampler"] = "gpry"
        # Add GPry-specific metadata (e.g., acquisition function, n_initial)
        self.attrs["n_initial"] = sampler.gp_config.get("n_initial", 1000)
        self.attrs["acquisition_function"] = sampler.gp_config.get("acquisition_function", "LogExp")
        self.attrs["checkpoint_interval"] = sampler.gp_config.get("checkpoint_interval", 100)

    #def write_samples(self, samples, **kwargs):
    #    """Write samples to the HDF file (required by PyCBC)."""
    #    if "samples" not in self:
    #        self.create_group("samples")
    #    group = self["samples"]
    #    for param in samples:
    #        if param in group:
    #            del group[param]
    #        group.create_dataset(param, data=samples[param])

    def write_samples(self, samples, **kwargs):
        """
        Write posterior samples into the 'samples' group.
        samples : dict mapping param names → numpy arrays of shape (nsamples,)
        """
        group = self.require_group('samples')
        for param, arr in samples.items():
            # overwrite if it already exists
            if param in group:
                del group[param]
            group.create_dataset(param, data=arr)

    def read_samples(self, parameters=None, **kwargs):
        """Read samples from the HDF file (required by PyCBC)."""
        group = self["samples"]
        if parameters is None:
            parameters = list(group.keys())
        return {param: group[param][:] for param in parameters}

    # Keep existing methods (write_raw_samples, read_raw_samples, etc.)
    def write_raw_samples(self, samples):
        """Write training points from GPry's active learning phase."""
        if "raw_samples" in self:
            del self["raw_samples"]
        group = self.create_group("raw_samples")
        for param in samples:
            group.create_dataset(param, data=samples[param])

    def read_raw_samples(self, **kwargs):
        """Read training points."""
        group = self["raw_samples"]
        return {param: group[param][:] for param in group.keys()}

    def write_posterior_samples(self, samples):
        """Write final posterior samples from GP approximation."""
        if "posterior" in self:
            del self["posterior"]
        group = self.create_group("posterior")
        if isinstance(samples, dict):
            dtype = [(name, 'f8') for name in samples.keys()]
            arr = np.zeros(len(next(iter(samples.values()))), dtype=dtype)
            for name in samples:
                arr[name] = samples[name]
        elif isinstance(samples, FieldArray):
            arr = samples
        else:
            raise ValueError("Unsupported samples type")
        group.create_dataset("samples", data=arr)

    def write_logevidence(self, logz, dlogz):
        """Write log-evidence and its error."""
        self.attrs["logz"] = logz
        self.attrs["dlogz"] = dlogz

    def read_logevidence(self):
        """Read log-evidence and its error."""
        return (self.attrs["logz"], self.attrs["dlogz"])

    def write_checkpoint(self, state):
        """Save GPry's state (training points, GP hyperparameters)."""
        if "checkpoint" in self:
            del self["checkpoint"]
        group = self.create_group("checkpoint")
        # Example: Save training points
        group.create_dataset("training_points", data=state.training_points)
        # Add other state components (log_likelihoods, hyperparameters, etc.)

    def read_checkpoint(self):
        """Load GPry's state from the checkpoint."""
        group = self["checkpoint"]
        state = {
            "training_points": group["training_points"][:],
            # Load other components
        }
        return state

    def write_run_end_time(self):
        # Either write the run end time to the file attributes or log a message.
        self.attrs["run_end_time"] = str(datetime.datetime.now())
