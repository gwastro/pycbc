"""Provides IO for the nessai sampler"""

from .base_nested_sampler import BaseNestedSamplerFile

from ...io.hdf import dump_state, load_state
from .posterior import write_samples_to_file


class NessaiFile(BaseNestedSamplerFile):
    """Class to handle file IO for the ``nessai`` sampler."""

    name = "nessai_file"

    def write_pickled_data_into_checkpoint_file(self, data):
        """Write the pickled data into a checkpoint file"""
        if "sampler_info/saved_state" not in self:
            self.create_group("sampler_info/saved_state")
        dump_state(data, self, path="sampler_info/saved_state")

    def read_pickled_data_from_checkpoint_file(self):
        """Read the pickled data from a checkpoint file"""
        return load_state(self, path="sampler_info/saved_state")

    def write_raw_samples(self, data, parameters=None):
        """Write the nested samples to the file"""
        if "samples" not in self:
            self.create_group("samples")
        write_samples_to_file(
            self, data, parameters=parameters, group="samples"
        )
