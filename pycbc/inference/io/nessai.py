"""Provides IO for the nessai sampler"""
from .base_nested_sampler import BaseNestedSamplerFile


class NessaiFile(BaseNestedSamplerFile):
    """Class to handle file IO for the ``nessai`` sampler."""

    name = 'nessai_file'
