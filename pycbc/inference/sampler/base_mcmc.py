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
"""Provides constructor classes and convenience functions for MCMC samplers."""

from __future__ import absolute_import

from abc import (ABCMeta, abstractmethod, abstractproperty)
import logging
import numpy
from pycbc.workflow import ConfigParser
from pycbc.filter import autocorrelation

from pycbc.io.inference import validate_checkpoint_files

#
# =============================================================================
#
#                              Convenience functions
#
# =============================================================================
#


def raw_samples_to_dict(sampler, raw_samples):
    """Convenience function for converting ND array to a dict of samples.

    The samples are assumed to have dimension
    ``[sampler.base_shape x] niterations x len(sampler.sampling_params)``.

    Parameters
    ----------
    sampler : sampler instance
        An instance of an MCMC sampler.
    raw_samples : array
        The array of samples to convert.

    Returns
    -------
    dict :
        A dictionary mapping the raw samples to the variable params. If the
        sampling params are not the same as the variable params, they will
        also be included. Each array will have shape
        ``[sampler.base_shape x] niterations``.
    """
    sampling_params = sampler.sampling_params
    # convert to dictionary
    samples = {param: raw_samples[..., ii] for
               ii, param in enumerate(sampling_params)}
    # apply boundary conditions
    samples = sampler.model.prior_distribution.apply_boundary_conditions(
        **samples)
    # apply transforms to go to model's variable params space
    if sampler.model.sampling_transforms is not None:
        samples = sampler.model.sampling_transforms.apply(
            samples, inverse=True)
    return samples


def blob_data_to_dict(stat_names, blobs):
    """Converts list of "blobs" to a dictionary of model stats.

    Samplers like ``emcee`` store the extra tuple returned by ``CallModel`` to
    a list called blobs. This is a list of lists of tuples with shape
    niterations x nwalkers x nstats, where nstats is the number of stats
    returned by the model's ``default_stats``. This converts that list to a
    dictionary of arrays keyed by the stat names.

    Parameters
    ----------
    stat_names : list of str
        The list of the stat names.
    blobs : list of list of tuples
        The data to convert.

    Returns
    -------
    dict :
        A dictionary mapping the model's ``default_stats`` to arrays of values.
        Each array will have shape ``nwalkers x niterations``.
    """
    # get the dtypes of each of the stats; we'll just take this from the
    # first iteration and walker
    dtypes = [type(val) for val in blobs[0][0]]
    assert len(stat_names) == len(dtypes), (
        "number of stat names must match length of tuples in the blobs")
    # convert to an array; to ensure that we get the dtypes correct, we'll
    # cast to a structured array
    raw_stats = numpy.array(blobs, dtype=zip(stat_names, dtypes))
    # transpose so that it has shape nwalkers x niterations
    raw_stats = raw_stats.transpose()
    # now return as a dictionary
    return {stat: raw_stats[stat] for stat in stat_names}


def get_optional_arg_from_config(cp, section, arg, dtype=str):
    """Convenience function to retrieve an optional argument from a config
    file.

    Parameters
    ----------
    cp : ConfigParser
        Open config parser to retrieve the argument from.
    section : str
        Name of the section to retrieve from.
    arg : str
        Name of the argument to retrieve.
    dtype : datatype, optional
        Cast the retrieved value (if it exists) to the given datatype. Default
        is ``str``.

    Returns
    -------
    val : None or str
        If the argument is present, the value. Otherwise, None.
    """
    if cp.has_option(section, arg):
        val = dtype(cp.get(section, arg))
    else:
        val = None
    return val


#
# =============================================================================
#
#                              BaseMCMC definition
#
# =============================================================================
#


class BaseMCMC(object):
    """This class provides methods common to MCMCs.

    It is not a sampler class itself. Sampler classes can inherit from this
    along with ``BaseSampler``.

    Attributes
    ----------
    p0 : dict
        A dictionary of the initial position of the walkers. Set by using
        ``set_p0``. If not set yet, a ``ValueError`` is raised when the
        attribute is accessed.
    pos : dict
        A dictionary of the current walker positions. If the sampler hasn't
        been run yet, returns p0.
    """
    __metaclass__ = ABCMeta

    _lastclear = None  # the iteration when samples were cleared from memory
    _itercounter = None  # the number of iterations since the last clear
    _pos = None
    _p0 = None
    _nwalkers = None
    _burn_in = None
    _checkpoint_interval = None
    _target_niterations = None
    _target_eff_nsamples = None

    @abstractproperty
    def base_shape(self):
        """What shape the sampler's samples arrays are in, excluding
        the iterations dimension.

        For example, if a sampler uses 20 walkers and 3 temperatures, this
        would be ``(3, 20)``. If a sampler only uses a single walker and no
        temperatures this would be ``()``.
        """
        pass

    @property
    def nwalkers(self):
        """Get the number of walkers."""
        if self._nwalkers is None:
            raise ValueError("number of walkers not set")
        return self._nwalkers

    @property
    def niterations(self):
        """Get the current number of iterations."""
        itercounter = self._itercounter
        if itercounter is None:
            itercounter = 0
        lastclear = self._lastclear
        if lastclear is None:
            lastclear = 0
        return itercounter + lastclear

    @property
    def checkpoint_interval(self):
        """The number of iterations to do between checkpoints."""
        return self._checkpoint_interval

    @property
    def target_niterations(self):
        """The number of iterations the sampler should run for."""
        return self._target_niterations

    @property
    def target_eff_nsamples(self):
        """The target number of effective samples the sampler should get."""
        return self._target_eff_nsamples

    def set_target(self, niterations=None, eff_nsamples=None):
        """Sets the target niterations/nsamples for the sampler.

        One or the other must be provided, not both.
        """
        if niterations is None and eff_nsamples is None:
            raise ValueError("Must provide a target niterations or "
                             "eff_nsamples")
        if niterations is not None and eff_nsamples is not None:
            raise ValueError("Must provide a target niterations or "
                             "eff_nsamples, not both")
        self._target_niterations = int(niterations) \
            if niterations is not None else None
        self._target_eff_nsamples = int(eff_nsamples) \
            if eff_nsamples is not None else None

    @abstractmethod
    def clear_samples(self):
        """A method to clear samples from memory."""
        pass

    @property
    def pos(self):
        pos = self._pos
        if pos is None:
            return self.p0
        # convert to dict
        pos = {param: self._pos[..., k]
               for (k, param) in enumerate(self.sampling_params)}
        return pos

    @property
    def p0(self):
        """The starting position of the walkers in the sampling param space.

        The returned object is a dict mapping the sampling parameters to the
        values.
        """
        if self._p0 is None:
            raise ValueError("initial positions not set; run set_p0")
        # convert to dict
        p0 = {param: self._p0[..., k]
              for (k, param) in enumerate(self.sampling_params)}
        return p0

    def set_p0(self, samples_file=None, prior=None):
        """Sets the initial position of the walkers.

        Parameters
        ----------
        samples_file : InferenceFile, optional
            If provided, use the last iteration in the given file for the
            starting positions.
        prior : JointDistribution, optional
            Use the given prior to set the initial positions rather than
            ``model``'s prior.

        Returns
        -------
        p0 : dict
            A dictionary maping sampling params to the starting positions.
        """
        # if samples are given then use those as initial positions
        if samples_file is not None:
            with self.io(samples_file, 'r') as fp:
                samples = fp.read_samples(self.variable_params,
                                          iteration=-1)
                # make sure we have the same shape
                assert samples.shape == self.base_shape, (
                       "samples in file {} have shape {}, but I have shape {}".
                       format(samples_file, samples.shape, self.base_shape))
            # transform to sampling parameter space
            if self.model.sampling_transforms is not None:
                samples = self.model.sampling_transforms.apply(samples)
        # draw random samples if samples are not provided
        else:
            nsamples = numpy.prod(self.base_shape)
            samples = self.model.prior_rvs(size=nsamples, prior=prior).reshape(
                self.base_shape)
        # store as ND array with shape [base_shape] x nparams
        ndim = len(self.variable_params)
        p0 = numpy.ones(list(self.base_shape)+[ndim])
        for i, param in enumerate(self.sampling_params):
            p0[..., i] = samples[param]
        self._p0 = p0
        return self.p0

    def set_initial_conditions(self, initial_distribution=None,
                               samples_file=None):
        """Sets the initial starting point for the MCMC.

        If a starting samples file is provided, will also load the random
        state from it.
        """
        self.set_p0(samples_file=samples_file, prior=initial_distribution)
        # if a samples file was provided, use it to set the state of the
        # sampler
        if samples_file is not None:
            self.set_state_from_file(samples_file)

    @abstractmethod
    def set_state_from_file(self, filename):
        """Sets the state of the sampler to the instance saved in a file.
        """
        pass

    def run(self):
        """Runs the sampler."""
        if self.target_eff_nsamples and self.checkpoint_interval is None:
            raise ValueError("A checkpoint interval must be set if "
                             "targetting an effective number of samples")
        # get the starting number of samples:
        # "nsamples" keeps track of the number of samples we've obtained (if
        # target_eff_nsamples is not None, this is the effective number of
        # samples; otherwise, this is the total number of samples).
        # _lastclear is the number of iterations that the file already
        # contains (either due to sampler burn-in, or a previous checkpoint)
        if self.new_checkpoint:
            self._lastclear = 0
        else:
            with self.io(self.checkpoint_file, "r") as fp:
                self._lastclear = fp.niterations
        if self.target_eff_nsamples is not None:
            target_nsamples = self.target_eff_nsamples
            with self.io(self.checkpoint_file, "r") as fp:
                nsamples = fp.effective_nsamples
        elif self.target_niterations is not None:
            # the number of samples is the number of iterations times the
            # number of walkers
            target_nsamples = self.nwalkers * self.target_niterations
            nsamples = self._lastclear * self.nwalkers
        else:
            raise ValueError("must set either target_eff_nsamples or "
                             "target_niterations; see set_target")
        self._itercounter = 0
        # figure out the interval to use
        iterinterval = self.checkpoint_interval
        if iterinterval is None:
            iterinterval = self.target_niterations
        # run sampler until we have the desired number of samples
        while nsamples < target_nsamples:
            # adjust the interval if we would go past the number of iterations
            if self.target_niterations is not None and (
                    self.niterations + iterinterval > self.target_niterations):
                iterinterval = self.target_niterations - self.niterations
            # run sampler and set initial values to None so that sampler
            # picks up from where it left off next call
            logging.info("Running sampler for {} to {} iterations".format(
                self.niterations, self.niterations + iterinterval))
            # run the underlying sampler for the desired interval
            self.run_mcmc(iterinterval)
            # update the itercounter
            self._itercounter = self._itercounter + iterinterval
            # dump the current results
            self.checkpoint()
            # update nsamples for next loop
            if self.target_eff_nsamples is not None:
                nsamples = self.effective_nsamples
                logging.info("Have {} effective samples post burn in".format(
                    nsamples))
            else:
                nsamples += iterinterval * self.nwalkers

    @property
    def burn_in(self):
        """The class for doing burn-in tests (if specified)."""
        return self._burn_in

    def set_burn_in(self, burn_in):
        """Sets the object to use for doing burn-in tests."""
        self._burn_in = burn_in

    @property
    def effective_nsamples(self):
        """The effective number of samples post burn-in that the sampler has
        acquired so far."""
        try:
            acl = numpy.array(self.acls.values()).max()
        except (AttributeError, TypeError):
            acl = numpy.inf
        if self.burn_in is None:
            nperwalker = max(int(self.niterations // acl), 1)
        elif self.burn_in.is_burned_in:
            nperwalker = int(
                (self.niterations - self.burn_in.burn_in_iteration) // acl)
            # after burn in, we always have atleast 1 sample per walker
            nperwalker = max(nperwalker, 1)
        else:
            nperwalker = 0
        return self.nwalkers * nperwalker

    @abstractmethod
    def run_mcmc(self, niterations):
        """Run the MCMC for the given number of iterations."""
        pass

    @abstractmethod
    def write_results(self, filename):
        """Should write all samples currently in memory to the given file."""
        pass

    def checkpoint(self):
        """Dumps current samples to the checkpoint file."""
        # write new samples
        logging.info("Writing samples to files")
        for fn in [self.checkpoint_file, self.backup_file]:
            self.write_results(fn)
            with self.io(fn, "a") as fp:
                # write the current number of iterations
                fp.write_niterations(self.niterations)
        # check for burn in, compute the acls
        self.acls = None
        if self.burn_in is not None:
            logging.info("Updating burn in")
            self.burn_in.evaluate(self.checkpoint_file)
            burn_in_iter = self.burn_in.burn_in_iteration
            logging.info("Is burned in: {}".format(self.burn_in.is_burned_in))
            if self.burn_in.is_burned_in:
                logging.info("Burn-in iteration: {}".format(
                    self.burn_in.burn_in_iteration))
        else:
            burn_in_iter = 0
        # Compute acls; the burn_in test may have calculated an acl and saved
        # it, in which case we don't need to do it again.
        if self.acls is None:
            logging.info("Computing acls")
            self.acls = self.compute_acl(self.checkpoint_file,
                                         start_index=burn_in_iter)
        logging.info("ACL: {}".format(numpy.array(self.acls.values()).max()))
        # write
        for fn in [self.checkpoint_file, self.backup_file]:
            with self.io(fn, "a") as fp:
                if self.burn_in is not None:
                    fp.write_burn_in(self.burn_in)
                if self.acls is not None:
                    fp.write_acls(self.acls)
                # write effective number of samples
                fp.write_effective_nsamples(self.effective_nsamples)
        # check validity
        logging.info("Validating checkpoint and backup files")
        checkpoint_valid = validate_checkpoint_files(
            self.checkpoint_file, self.backup_file)
        if not checkpoint_valid:
            raise IOError("error writing to checkpoint file")
        # clear the in-memory chain to save memory
        logging.info("Clearing samples from memory")
        self.clear_samples()

    @staticmethod
    def checkpoint_from_config(cp, section):
        """Gets the checkpoint interval from the given config file.

        This looks for 'checkpoint-interval' in the section.

        Parameters
        ----------
        cp : ConfigParser
            Open config parser to retrieve the argument from.
        section : str
            Name of the section to retrieve from.

        Return
        ------
        int or None :
            The checkpoint interval, if it is in the section. Otherw
        """
        return get_optional_arg_from_config(cp, section, 'checkpoint-interval',
                                            dtype=int)

    def set_target_from_config(self, cp, section):
        """Sets the target using the given config file.

        This looks for 'niterations' to set the ``target_niterations``, and
        'effective-nsamples' to set the ``target_eff_nsamples``.

        Parameters
        ----------
        cp : ConfigParser
            Open config parser to retrieve the argument from.
        section : str
            Name of the section to retrieve from.
        """
        if cp.has_option(section, "niterations"):
            niterations = int(cp.get(section, "niterations"))
        else:
            niterations = None
        if cp.has_option(section, "effective-nsamples"):
            nsamples = int(cp.get(section, "effective-nsamples"))
        else:
            nsamples = None
        self.set_target(niterations=niterations, eff_nsamples=nsamples)

    def set_burn_in_from_config(self, cp):
        """Sets the burn in class from the given config file.

        If no burn-in section exists in the file, then this just set the
        burn-in class to None.
        """
        try:
            bit = self.burn_in_class.from_config(cp, self)
        except ConfigParser.Error:
            bit = None
        self.set_burn_in(bit)

    @abstractmethod
    def compute_acf(cls, filename, **kwargs):
        """A method to compute the autocorrelation function of samples in the
        given file."""
        pass

    @abstractmethod
    def compute_acl(cls, filename, **kwargs):
        """A method to compute the autocorrelation length of samples in the
        given file."""
        pass


class MCMCAutocorrSupport(object):
    """Provides class methods for calculating ensemble ACFs/ACLs.
    """

    @classmethod
    def compute_acf(cls, filename, start_index=None, end_index=None,
                    per_walker=False, walkers=None, parameters=None):
        """Computes the autocorrleation function of the model params in the
        given file.

        By default, parameter values are averaged over all walkers at each
        iteration. The ACF is then calculated over the averaged chain. An
        ACF per-walker will be returned instead if ``per_walker=True``.

        Parameters
        -----------
        filename : str
            Name of a samples file to compute ACFs for.
        start_index : {None, int}
            The start index to compute the acl from. If None, will try to use
            the number of burn-in iterations in the file; otherwise, will start
            at the first sample.
        end_index : {None, int}
            The end index to compute the acl to. If None, will go to the end
            of the current iteration.
        per_walker : optional, bool
            Return the ACF for each walker separately. Default is False.
        walkers : optional, int or array
            Calculate the ACF using only the given walkers. If None (the
            default) all walkers will be used.
        parameters : optional, str or array
            Calculate the ACF for only the given parameters. If None (the
            default) will calculate the ACF for all of the model params.

        Returns
        -------
        dict :
            Dictionary of arrays giving the ACFs for each parameter. If
            ``per-walker`` is True, the arrays will have shape
            ``nwalkers x niterations``.
        """
        acfs = {}
        with cls._io(filename, 'r') as fp:
            if parameters is None:
                parameters = fp.variable_params
            if isinstance(parameters, str) or isinstance(parameters, unicode):
                parameters = [parameters]
            for param in parameters:
                if per_walker:
                    # just call myself with a single walker
                    if walkers is None:
                        walkers = numpy.arange(fp.nwalkers)
                    arrays = [
                        cls.compute_acf(filename, start_index=start_index,
                                        end_index=end_index,
                                        per_walker=False, walkers=ii,
                                        parameters=param)[param]
                        for ii in walkers]
                    acfs[param] = numpy.vstack(arrays)
                else:
                    samples = fp.read_raw_samples(
                        param, thin_start=start_index, thin_interval=1,
                        thin_end=end_index, walkers=walkers,
                        flatten=False)[param]
                    samples = samples.mean(axis=0)
                    acfs[param] = autocorrelation.calculate_acf(
                        samples).numpy()
        return acfs

    @classmethod
    def compute_acl(cls, filename, start_index=None, end_index=None,
                    min_nsamples=10):
        """Computes the autocorrleation length for all model params in the
        given file.

        Parameter values are averaged over all walkers at each iteration.
        The ACL is then calculated over the averaged chain. If an ACL cannot
        be calculated because there are not enough samples, it will be set
        to ``inf``.

        Parameters
        -----------
        filename : str
            Name of a samples file to compute ACLs for.
        start_index : int, optional
            The start index to compute the acl from. If None, will try to use
            the number of burn-in iterations in the file; otherwise, will start
            at the first sample.
        end_index : int, optional
            The end index to compute the acl to. If None, will go to the end
            of the current iteration.
        min_nsamples : int, optional
            Require a minimum number of samples to compute an ACL. If the
            number of samples per walker is less than this, will just set to
            ``inf``. Default is 10.

        Returns
        -------
        dict
            A dictionary giving the ACL for each parameter.
        """
        acls = {}
        with cls._io(filename, 'r') as fp:
            for param in fp.variable_params:
                samples = fp.read_raw_samples(
                    param, thin_start=start_index, thin_interval=1,
                    thin_end=end_index, flatten=False)[param]
                samples = samples.mean(axis=0)
                # if < min number of samples, just set to inf
                if samples.size < min_nsamples:
                    acl = numpy.inf
                else:
                    acl = autocorrelation.calculate_acl(samples)
                if acl <= 0:
                    acl = numpy.inf
                acls[param] = acl
        return acls
