# Copyright (C) 2016  Christopher M. Biwer, Collin Capano
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
Defines the base sampler class to be inherited by all samplers.
"""


from abc import ABCMeta, abstractmethod, abstractproperty
import shutil
import logging
import numpy

from six import add_metaclass

from pycbc import distributions
from pycbc.inference.io import validate_checkpoint_files

#
# =============================================================================
#
#                           Base Sampler definition
#
# =============================================================================
#


@add_metaclass(ABCMeta)
class BaseSampler(object):
    """Abstract base class for all inference samplers.

    All sampler classes must inherit from this class and implement its abstract
    methods.

    Parameters
    ----------
    model : Model
        An instance of a model from ``pycbc.inference.models``.
    """
    name = None

    def __init__(self, model):
        self.model = model
        self.checkpoint_file = None
        self.backup_file = None
        self.checkpoint_valid = None
        self.new_checkpoint = None

    # @classmethod <--uncomment when we move to python 3.3
    @abstractmethod
    def from_config(cls, cp, model, output_file=None, nprocesses=1,
                    use_mpi=False):
        """This should initialize the sampler given a config file.
        """
        pass

    @property
    def variable_params(self):
        """Returns the parameters varied in the model.
        """
        return self.model.variable_params

    @property
    def sampling_params(self):
        """Returns the sampling params used by the model.
        """
        return self.model.sampling_params

    @property
    def static_params(self):
        """Returns the model's fixed parameters.
        """
        return self.model.static_params

    @abstractproperty
    def samples(self):
        """A dict mapping variable_params to arrays of samples currently
        in memory. The dictionary may also contain sampling_params.

        The sample arrays may have any shape, and may or may not be thinned.
        """
        pass

    @abstractproperty
    def model_stats(self):
        """A dict mapping model's metadata fields to arrays of values for
        each sample in ``raw_samples``.

        The arrays may have any shape, and may or may not be thinned.
        """
        pass

    @abstractmethod
    def run(self):
        """This function should run the sampler.

        Any checkpointing should be done internally in this function.
        """
        pass

    @abstractproperty
    def io(self):
        """A class that inherits from ``BaseInferenceFile`` to handle IO with
        an hdf file.

        This should be a class, not an instance of class, so that the sampler
        can initialize it when needed.
        """
        pass

    def model_stats_from_cache(self, samples=None):
        """The model's ``default_stats`` as arrays over the given samples.

        Taken from what the model remembered as it evaluated. A point it did
        not remember is worked out, so the answer is complete either way and
        the cache only decides what it costs.

        Parameters
        ----------
        samples : dict, optional
            The points to look up. Defaults to the sampler's own, which an
            MCMC sampler clears as it checkpoints, so a caller running after
            the fact should read them back from the file instead.
        """
        model = self.model
        model.gather_stats_cache(getattr(self, 'pool', None))
        params = list(model.variable_params)
        samples = self.samples if samples is None else samples
        shape = numpy.asarray(samples[params[0]]).shape
        flat = {p: numpy.asarray(samples[p]).flatten() for p in params}
        names = model.default_stats
        out = numpy.full((flat[params[0]].size, len(names)), numpy.nan)
        missed = 0
        for i in range(out.shape[0]):
            point = {p: flat[p][i] for p in params}
            cached = model.cached_stats(point, names)
            if cached is None:
                # Not remembered here: a sampler that evaluates somewhere we
                # cannot reach, or a bound that dropped it. Work it out, which
                # is what getting these afterwards has always cost.
                model.update(**point)
                model.logposterior  # pylint:disable=pointless-statement
                cached = model.get_current_stats(names)
                missed += 1
            out[i] = cached
        if missed:
            logging.info("Recalculated model stats for %s of %s points",
                         missed, out.shape[0])
        return {name: out[:, j].reshape(shape)
                for j, name in enumerate(names)}

    def write_cached_stats(self):
        """Fill in any stats the sampler did not write for itself.

        The MCMC samplers store them as they checkpoint; the rest drop them.
        Stats already in the file are left alone.
        """
        with self.io(self.checkpoint_file, 'r') as fp:
            group = fp[fp.samples_group]
            missing = [s for s in self.model.default_stats if s not in group]
            if not missing:
                return
            samples = {p: group[p][()] for p in self.model.variable_params}
        stats = self.model_stats_from_cache(samples)
        for fn in [self.checkpoint_file, self.backup_file]:
            with self.io(fn, 'a') as fp:
                for name in missing:
                    fp.write_data(name, stats[name], path=fp.samples_group,
                                  append=False)

    @abstractmethod
    def checkpoint(self):
        """The sampler must have a checkpoint method for dumping raw samples
        and stats to the file type defined by ``io``.
        """
        pass

    @abstractmethod
    def finalize(self):
        """Do any finalization to the samples file before exiting."""
        pass

    @abstractmethod
    def resume_from_checkpoint(self):
        """Resume the sampler from the output file.
        """
        pass

#
# =============================================================================
#
#                           Convenience functions
#
# =============================================================================
#


def setup_output(sampler, output_file, check_nsamples=True, validate=True):
    r"""Sets up the sampler's checkpoint and output files.

    The checkpoint file has the same name as the output file, but with
    ``.checkpoint`` appended to the name. A backup file will also be
    created.

    Parameters
    ----------
    sampler : sampler instance
        Sampler
    output_file : str
        Name of the output file.
    """
    # check for backup file(s)
    checkpoint_file = output_file + '.checkpoint'
    backup_file = output_file + '.bkup'
    # check if we have a good checkpoint and/or backup file
    logging.info("Looking for checkpoint file")
    checkpoint_valid = False
    if validate:
        checkpoint_valid = validate_checkpoint_files(checkpoint_file,
                                                     backup_file,
                                                     check_nsamples)
    # Create a new file if the checkpoint doesn't exist, or if it is
    # corrupted
    sampler.new_checkpoint = False  # keeps track if this is a new file or not
    if not checkpoint_valid:
        logging.info("Checkpoint not found or not valid")
        create_new_output_file(sampler, checkpoint_file)
        # now the checkpoint is valid
        sampler.new_checkpoint = True
        # copy to backup
        shutil.copy(checkpoint_file, backup_file)
    # write the command line, startup
    for fn in [checkpoint_file, backup_file]:
        with sampler.io(fn, "a") as fp:
            fp.write_command_line()

            fp.write_resume_point()
            fp.write_run_start_time()
    # store
    sampler.checkpoint_file = checkpoint_file
    sampler.backup_file = backup_file


def create_new_output_file(sampler, filename, **kwargs):
    r"""Creates a new output file.

    Parameters
    ----------
    sampler : sampler instance
        Sampler
    filename : str
        Name of the file to create.
    \**kwargs :
        All other keyword arguments are passed through to the file's
        ``write_metadata`` function.
    """
    logging.info("Creating file {}".format(filename))
    with sampler.io(filename, "w") as fp:
        # create the samples group and sampler info group
        fp.create_group(fp.samples_group)
        fp.create_group(fp.sampler_group)
        # save the sampler's metadata
        fp.write_sampler_metadata(sampler)


def initial_dist_from_config(cp, variable_params, static_params=None):
    r"""Loads a distribution for the sampler start from the given config file.

    A distribution will only be loaded if the config file has a [initial-\*]
    section(s).

    Parameters
    ----------
    cp : Config parser
        The config parser to try to load from.
    variable_params : list of str
        The variable parameters for the distribution.
    static_params : dict, optional
        The static parameters used to place constraints on the
        distribution.

    Returns
    -------
    JointDistribution or None :
        The initial distribution. If no [initial-\*] section found in the
        config file, will just return None.
    """
    if len(cp.get_subsections("initial")):
        logging.info("Using a different distribution for the starting points "
                     "than the prior.")
        initial_dists = distributions.read_distributions_from_config(
            cp, section="initial")
        constraints = distributions.read_constraints_from_config(
            cp, constraint_section="initial_constraint",
            static_args=static_params)
        init_dist = distributions.JointDistribution(
            variable_params, *initial_dists,
            **{"constraints": constraints})
    else:
        init_dist = None
    return init_dist
