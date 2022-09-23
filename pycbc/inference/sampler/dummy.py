""" Dummy class when no actual sampling is needed, but we may want to do
some reconstruction supported by the likelihood model.
"""

import numpy

from pycbc.inference.io import PosteriorFile
from pycbc.inference import models
from pycbc.pool import choose_pool

from .base import (BaseSampler, setup_output)


def call_reconstruct(iteration):
    """ Accessor to update the global model and call its reconstruction
    routine.
    """
    models._global_instance.update()
    return models._global_instance.reconstruct(seed=iteration)


class DummySampler(BaseSampler):
    """Dummy sampler for not doing sampling

    Parameters
    ----------
    model : Model
        An instance of a model from ``pycbc.inference.models``.
    """
    name = 'dummy'

    def __init__(self, model, *args, nprocesses=1, use_mpi=False,
                 num_samples=1000, **kwargs):
        super().__init__(model, *args)

        models._global_instance = model
        self.num_samples = int(num_samples)
        self.pool = choose_pool(mpi=use_mpi, processes=nprocesses)
        self._samples = {}

    @classmethod
    def from_config(cls, cp, model, output_file=None, nprocesses=1,
                    use_mpi=False):
        """This should initialize the sampler given a config file.
        """
        kwargs = {k: cp.get('sampler', k) for k in cp.options('sampler')}
        obj = cls(model, nprocesses=nprocesses, use_mpi=use_mpi, **kwargs)
        setup_output(obj, output_file, check_nsamples=False, validate=False)
        return obj

    @property
    def samples(self):
        """A dict mapping variable_params to arrays of samples currently
        in memory. The dictionary may also contain sampling_params.

        The sample arrays may have any shape, and may or may not be thinned.
        """
        return self._samples

    @property
    def model_stats(self):
        pass

    def run(self):
        samples = self.pool.map(call_reconstruct,
                                range(self.num_samples))
        self._samples = {k: numpy.array([x[k] for x in samples])
                         for k in samples[0]}

    def finalize(self):
        with self.io(self.checkpoint_file, "a") as fp:
            fp.write_samples(samples=self._samples)

    checkpoint = resume_from_checkpoint = run

    @property
    def io(self):
        """A class that inherits from ``BaseInferenceFile`` to handle IO with
        an hdf file.

        This should be a class, not an instance of class, so that the sampler
        can initialize it when needed.
        """
        return PosteriorFile
