# Copyright (C) 2016  Christopher M. Biwer
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
This modules provides classes and functions for using different sampler
packages for parameter estimation.
"""

import numpy
from pycbc.io import WaveformArray

#
# =============================================================================
#
#                                   Helper functions
#
# =============================================================================
#
def get_slice(fp, thin_start=None, thin_interval=None, thin_end=None):
    """Formats a slice using the given arguments that can be used to retrieve
    a thinned array from an InferenceFile.

    Parameters
    ----------
    fp : InferenceFile
        An open inference file. The `burn_in_iterations` and acl will try to
        be obtained from the file's attributes.
    thin_start : {None, int}
        The starting index to use. If None, will try to retrieve the
        `burn_in_iterations` from the given file. If no `burn_in_iterations`
        exists, will default to the start of the array.
    thin_interval : {None, int}
        The interval to use. If None, will try to retrieve the acl from the
        given file. If no acl attribute exists, will default to 1.
    thin_end : {None, int}
        The end index to use. If None, will retrieve to the end of the array.

    Returns
    -------
    slice :
        The slice needed.
    """
    # default is to skip burn in samples
    if thin_start is None:
        try:
            thin_start = fp.burn_in_iterations
        except KeyError:
            pass
    # default is to use stored ACL and accept every i-th sample
    if thin_interval is None:
        try:
            thin_interval = int(numpy.ceil(fp.acl))
        except KeyError:
            pass
    return slice(thin_start, thin_end, thin_interval)


#
# =============================================================================
#
#                                   Samplers
#
# =============================================================================
#

class _BaseSampler(object):
    """Base container class for running the inference sampler that will
    generate the posterior distributions.

    Parameters
    ----------
    sampler : sampler class
        An instance of the sampler class from its package.
    """
    name = None

    def __init__(self, likelihood_evaluator):
        self.likelihood_evaluator = likelihood_evaluator
        self.burn_in_iterations = 0

    @property
    def ifos(self):
        """Returns the ifos that were sampled over."""
        return self.likelihood_evaluator.waveform_generator.detector_names

    @property
    def variable_args(self):
        """Returns the variable args used by the likelihood evaluator.
        """
        return self.likelihood_evaluator.waveform_generator.variable_args

    @property
    def niterations(self):
        """Get the current number of iterations."""
        return self.chain.shape[0]

    @property
    def acceptance_fraction(self):
        """This function should return the fraction of walkers that accepted
        each step as an array.
        """
        return NotImplementedError("acceptance_fraction function not set.")

    @property
    def lnpost(self):
        """This function should return the natural logarithm of the likelihood
        as an niterations x nwalker array.
        """
        return NotImplementedError("lnpost function not set.")

    @property
    def chain(self):
        """This function should return the past samples as a
        niterations x nwalker x ndim array.
        """
        return NotImplementedError("chain function not set.")

    def burn_in(self, initial_values):
        """This function should burn in the sampler.
        """
        raise NotImplementedError("This sampler has no burn_in function.")

    def run(self, niterations):
        """This function should run the sampler.
        """
        raise NotImplementedError("run function not set.")

    # write and read functions
    def write_metadata(self, fp):
        """Writes metadata about this sampler to the given file. Metadata is
        written to the file's `attrs`.

        Parameters
        ----------
        fp : InferenceFile
            A file handler to an open inference file.
        """
        fp.attrs['sampler'] = self.name
        fp.attrs['ifos'] = self.ifos
        fp.attrs['variable_args'] = self.variable_args
        fp.attrs["niterations"] = self.niterations


class _BaseMCMCSampler(_BaseSampler):
    """This class is used to construct the MCMC sampler from the kombine-like
    packages.

    Parameters
    ----------
    likelihood_evaluator : likelihood class
        An instance of the likelihood class from the
        pycbc.inference.likelihood module.
    sampler : sampler instance
        An instance of an MCMC sampler similar to kombine or emcee.

    Attributes
    ----------
    sampler :
        The MCMC sampler instance used.
    p0 : nwalkers x ndim array
        The initial position of the walkers. Set by using set_p0. If not set
        yet, a ValueError is raised when the attribute is accessed.
    pos : {None, array}
        An array of the current walker positions.
    """
    name = None
    def __init__(self, sampler, likelihood_evaluator):
        self._sampler = sampler 
        self._pos = None
        self._p0 = None
        self.burn_in_iterations = 0
        # initialize
        super(_BaseMCMCSampler, self).__init__(likelihood_evaluator)

    @property
    def sampler(self):
        return self._sampler

    @property
    def pos(self):
        return self._pos

    def set_p0(self, p0):
        """Sets the initial position of the walkers.

        Parameters
        ----------
        p0 : numpy.array
            An nwalkers x ndim array of initial values for walkers.
        """
        self._p0 = p0

    @property
    def p0(self):
        if self._p0 is None:
            raise ValueError("initial positions not set; run set_p0")
        return self._p0

    @property
    def nwalkers(self):
        """Get the number of walkers."""
        return self.chain.shape[1]

    @property
    def acceptance_fraction(self):
        """Get the fraction of walkers that accepted each step as an arary.
        """
        return self._sampler.acceptance_fraction

    # write and read functions
    def write_metadata(self, fp):
        """Writes metadata about this sampler to the given file. Metadata is
        written to the file's `attrs`.

        Parameters
        ----------
        fp : InferenceFile
            A file handler to an open inference file.
        """
        super(_BaseMCMCSampler, self).write_metadata(fp)
        # add info about walkers, burn in
        fp.attrs["nwalkers"] = self.nwalkers
        fp.attrs['burn_in_iterations'] = self.burn_in_iterations

    def write_chain(self, fp, max_iterations=None):
        """Writes the samples from the current chain to the given file. Results
        are written to: `fp[fp.samples_group/{vararg}/walker{i}]`, where
        `{vararg}` is the name of a variable arg, and `{i}` is the index of
        a walker.

        Parameters
        -----------
        fp : InferenceFile
            A file handler to an open inference file.
        max_iterations : {None, int}
            If samples have not previously been written to the file, a new
            dataset will be created. By default, the size of this dataset will
            be whatever the length of the sampler's chain is at this point. If
            you intend to run more iterations, set this value to that size so
            that the array in the file will be large enough to accomodate
            future data.
        """
        # transpose samples to get an (ndim,nwalkers,niteration) array
        samples = numpy.transpose(self.chain)
        _, nwalkers, niterations = samples.shape

        group = fp.samples_group + '/{name}/walker{wnum}'

        # create an empty array if desired, in case this is the first time
        # writing
        if max_iterations is not None:
            if max_iterations < niterations:
                raise IndexError("The provided max size is less than the "
                    "number of iterations")
            out = numpy.zeros(max_iterations, dtype=samples.dtype)

        # loop over number of dimensions
        for i,dim_name in enumerate(self.variable_args):
            # loop over number of walkers
            for j in range(nwalkers):
                dataset_name = group.format(name=dim_name, wnum=j)
                try:
                    fp[dataset_name][:niterations] = samples[i,j,:]
                except KeyError:
                    # dataset doesn't exist yet, see if a larger array is
                    # desired
                    if max_iterations is not None:
                        out[:niterations] = samples[i,j,:]
                        fp[dataset_name] = out
                    else:
                        fp[dataset_name] = samples[i,j,:]

    def write_lnpost(self, fp, max_iterations=None):
        """Writes the `lnpost`s to the given file. Results are written to:
        `fp[fp.samples_group/lnpost/walker{i}]`, where `{i}` is the index of
        a walker.

        Parameters
        -----------
        fp : InferenceFile
            A file handler to an open inference file.
        max_iterations : {None, int}
            If `lnpost`s have not previously been written to the file, a new
            dataset will be created. By default, the size of this dataset will
            be whatever the length of the sampler's chain is at this point. If
            you intend to run more iterations, set this value to that size so
            that the array in the file will be large enough to accomodate
            future data.
        """
        # lnposts are an (niterations,nwalkers) array
        lnposts = self.lnpost
        niterations, nwalkers = lnposts.shape

        group = fp.samples_group + '/lnpost/walker{wnum}'

        # create an empty array if desired, in case this is the first time
        # writing
        if max_iterations is not None:
            if max_iterations < niterations:
                raise IndexError("The provided max size is less than the "
                    "number of iterations")
            out = numpy.zeros(max_iterations, dtype=lnposts.dtype)

        # loop over number of walkers
        for j in range(nwalkers):
            dataset_name = group.format(wnum=j)
            try:
                fp[dataset_name][:niterations] = lnposts[:,j]
            except KeyError:
                # dataset doesn't exist yet, see if a larger array is
                # desired
                if max_iterations is not None:
                    out[:niterations] = lnposts[:,j]
                    fp[dataset_name] = out
                else:
                    fp[dataset_name] = lnposts[:,j]

    def write_acceptance_fraction(self, fp, max_iterations=None):
        """Write acceptance_fraction data to file. Results are written to
        `fp[acceptance_fraction]`.

        Parameters
        -----------
        fp : InferenceFile
            A file handler to an open inference file.
        max_iterations : {None, int}
            If acceptance fraction have not previously been written to the
            file, a new dataset will be created. By default, the size of this
            dataset will be whatever the length of the sampler's chain is at
            this point. If you intend to run more iterations, set this value
            to that size so that arrays in the file will be large enough to
            accomodate future data.
        """
        dataset_name = "acceptance_fraction"
        acf = self.acceptance_fraction
        try:
            fp[dataset_name][:acf.size] = acf
        except KeyError:
            # dataset doesn't exist yet, see if a larger array is
            # desired
            if max_iterations is not None:
                out = numpy.zeros(max_iterations, dtype=acf.dtype)
                out[:acf.size] = acf
            else:
                out = acf
            fp[dataset_name] = out

    def write_results(self, fp, max_iterations=None):
        """Writes metadata, samples, lnpost, and acceptance fraction to the
        given file. See the write function for each of those for details.

        Parameters
        -----------
        fp : InferenceFile
            A file handler to an open inference file.
        max_iterations : {None, int}
            If results have not previously been written to the
            file, new datasets will be created. By default, the size of these
            datasets will be whatever the length of the sampler's chain is at
            this point. If you intend to run more iterations in the future,
            set this value to that size so that the array in the file will be
            large enough to accomodate future data.
        """
        self.write_metadata(fp)
        self.write_chain(fp, max_iterations=max_iterations)
        self.write_lnpost(fp, max_iterations=max_iterations)
        self.write_acceptance_fraction(fp, max_iterations=max_iterations)

    @staticmethod
    def read_samples(fp, parameters,
             thin_start=None, thin_interval=None, thin_end=None,
             iteration=None,
             walkers=None, flatten=True):
        """Reads samples for the given parameter(s).

        Parameters
        -----------
        fp : InferenceFile
            An open file handler to read the samples from.
        parameters : (list of) strings
            The parameter(s) to retrieve. A parameter can be the name of any
            field in `fp[fp.samples_group]`, a virtual field or method of
            `WaveformArray` (as long as the file contains the necessary fields
            to derive the virtual field or method), and/or a function of
            these.
        thin_start : int
            Index of the sample to begin returning samples. Default is to read
            samples after burn in. To start from the beginning set thin_start
            to 0.
        thin_interval : int
            Interval to accept every i-th sample. Default is to use the
            `fp.acl`. If `fp.acl` is not set, then use all samples
            (set thin_interval to 1).
        thin_end : int
            Index of the last sample to read. If not given then
            `fp.niterations` is used.
        iteration : int
            Get a single iteration. If provided, will override the
            `thin_{start/interval/end}` arguments.
        walkers : {None, (list of) int}
            The walker index (or a list of indices) to retrieve. If None,
            samples from all walkers will be obtained.
        flatten : {True, bool}
            The returned array will be one dimensional, with all desired
            samples from all desired walkers concatenated together. If False,
            the returned array will have dimension requested iterations x
            requested walkers.

        Returns
        -------
        WaveformArray
            Samples for the given parameters, as an instance of a
            WaveformArray.
        """
        # get the names of fields needed for the given parameters
        possible_fields = dict([[str(name), float]
            for name in fp[fp.samples_group].keys()])
        loadfields = WaveformArray.parse_parameters(parameters,
            possible_fields=possible_fields)

        # walkers to load
        if walkers is None:
            walkers = range(fp.nwalkers)
        if isinstance(walkers, int):
            walkers = [walkers]

        # get the slice to use
        if iteration is not None:
            get_index = iteration
        else:
            if thin_end is None:
                # use the number of current iterations
                thin_end = fp.niterations
            get_index = get_slice(fp, thin_start=thin_start, thin_end=thin_end,
                thin_interval=thin_interval)

        # load
        arrays = {}
        group = fp.samples_group + '/{name}/walker{wnum}'
        for name in loadfields:
            these_arrays = [
                    fp[group.format(name=name, wnum=ii)][get_index]
                    for ii in walkers]
            if flatten:
                arrays[name] = numpy.hstack(these_arrays)
            else:
                arrays[name] = numpy.vstack(these_arrays).transpose()
        return WaveformArray.from_kwargs(**arrays)


class KombineSampler(_BaseMCMCSampler):
    """This class is used to construct the MCMC sampler from the kombine
    package.

    Parameters
    ----------
    likelihood_evaluator : likelihood class
        An instance of the likelihood class from the
        pycbc.inference.likelihood module.
    nwalkers : int
        Number of walkers to use in sampler.
    ndim : int
        Number of dimensions in the parameter space. If transd is True this is
        the number of unique dimensions across the parameter spaces.
    transd : bool
        If True, the sampler will operate across parameter spaces using a
        kombine.clustered_kde.TransdimensionalKDE proposal distribution. In
        this mode a masked array with samples in each of the possible sets of
        dimensions must be given for the initial ensemble distribution.
    processes : {None, int}
        Number of processes to use with multiprocessing. If None, all available
        cores are used.
    """
    name = "kombine"

    def __init__(self, likelihood_evaluator, nwalkers=0, ndim=0,
                        transd=False, processes=None):
        try:
            import kombine
        except ImportError:
            raise ImportError("kombine is not installed.")

        # construct sampler for use in KombineSampler
        sampler = kombine.Sampler(nwalkers, ndim, likelihood_evaluator,
                                          transd=transd, processes=processes)
        # initialize
        super(KombineSampler, self).__init__(sampler, likelihood_evaluator)

    def run(self, niterations, **kwargs):
        """Advance the sampler for a number of samples.

        Parameters
        ----------
        niterations : int
            Number of samples to get from sampler.

        Returns
        -------
        p : numpy.array
            An array of current walker positions with shape (nwalkers, ndim).
        lnpost : numpy.array
            The list of log posterior probabilities for the walkers at
            positions p, with shape (nwalkers, ndim).
        lnprop : numpy.array
            The list of log proposal densities for the walkers at positions p,
            with shape (nwalkers, ndim).
        """
        if self.burn_in_iterations == 0:
            # no burn in, use the initial positions
            p0 = self.p0
        else:
            p0 = None
        p, lnpost, lnprop = self._sampler.run_mcmc(niterations, p0=p0,
            **kwargs)
        # update the positions
        self._pos = p
        return p, lnpost, lnprop

    @property
    def lnpost(self):
        """ Get the natural logarithm of the likelihood as an 
        niterations x nwalkers array.
        """
        return self._sampler.lnpost

    @property
    def chain(self):
        """Get all past samples as an niterations x nwalker x ndim array."""
        return self._sampler.chain

    def burn_in(self):
        """Evolve an ensemble until the acceptance rate becomes roughly
        constant. This is done by splitting acceptances in half and checking
        for statistical consistency. This isn't guaranteed to return a fully
        burned-in ensemble, but usually does. The initial positions (p0) must
        be set prior to running.

        Returns
        -------
        p : numpy.array
            An array of current walker positions with shape (nwalkers, ndim).
        lnpost : numpy.array
            The list of log posterior probabilities for the walkers at
            positions p, with shape (nwalkers, ndim).
        lnprop : numpy.array
            The list of log proposal densities for the walkers at positions p,
            with shape (nwalkers, ndim).
        """
        if self.burn_in_iterations == 0:
            res = self._sampler.burnin(self.p0)
            if len(res) == 4:
                p, post, q, _ = res
            else:
                p, post, q = res
            self.burn_in_iterations = self.chain.shape[0]
        else:
            raise ValueError("Burn in has already been performed")
        # update position
        self._pos = p
        return p, post, q


class EmceeEnsembleSampler(_BaseMCMCSampler):
    """This class is used to construct an MCMC sampler from the emcee
    package's EnsembleSampler.

    Parameters
    ----------
    likelihood_evaluator : likelihood class
        An instance of the likelihood class from the
        pycbc.inference.likelihood module.
    nwalkers : int
        Number of walkers to use in sampler.
    ndim : int
        Number of dimensions in the parameter space. If transd is True this is
        the number of unique dimensions across the parameter spaces.
    processes : {None, int}
        Number of processes to use with multiprocessing. If None, all available
        cores are used.
    """

    name = "emcee"

    def __init__(self, likelihood_evaluator, nwalkers=0, ndim=0,
                        processes=None):

        try:
            import emcee
        except ImportError:
            raise ImportError("emcee is not installed.")

        # initialize the pool to use
        if processes == 1:
            pool = None
        else:
            pool = emcee.interruptible_pool.InterruptiblePool(
                processes=processes)

        # construct the sampler
        sampler = emcee.EnsembleSampler(nwalkers, ndim, likelihood_evaluator,
            pool=pool)

        # initialize
        super(EmceeEnsembleSampler, self).__init__(sampler,
            likelihood_evaluator)

    @property
    def lnpost(self):
        """Get the natural logarithm of the likelihood as an 
        niterations x nwalkers array.
        """
        return self._sampler.lnprobability

    @property
    def chain(self):
        """Get all past samples as an niterations x nwalker x ndim array."""
        # emcee returns the chain as nwalker x niterations x ndim, so need to
        # transpose
        return self._sampler.chain.transpose((1,0,2))

    def run(self, niterations, **kwargs):
        """Advance the sampler for a number of samples.

        Parameters
        ----------
        niterations : int
            Number of samples to get from sampler.

        Returns
        -------
        p : numpy.array
            An array of current walker positions with shape (nwalkers, ndim).
        lnpost : numpy.array
            The list of log posterior probabilities for the walkers at
            positions p, with shape (nwalkers, ndim).
        lnprop : numpy.array
            The list of log proposal densities for the walkers at positions p,
            with shape (nwalkers, ndim).
        """
        pos = self._pos
        if pos is None:
            pos = self.p0
        p, lnpost, lnprop = self._sampler.run_mcmc(pos, niterations, **kwargs)
        # update the positions
        self._pos = p
        return p, lnpost, lnprop

    # Emcee defines acceptance fraction differently, so have to override
    # write functions
    def write_acceptance_fraction(self, fp):
        """Write acceptance_fraction data to file. Results are written to
        `fp[acceptance_fraction]`.

        Parameters
        -----------
        fp : InferenceFile
            A file handler to an open inference file.
        """
        dataset_name = "acceptance_fraction"
        try:
            fp[dataset_name][:] = self.acceptance_fraction
        except KeyError:
            # dataset doesn't exist yet, create it
            fp[dataset_name] = self.acceptance_fraction

    def write_results(self, fp, max_iterations=None):
        """Writes metadata, samples, lnpost, and acceptance fraction to the
        given file. See the write function for each of those for details.

        Parameters
        -----------
        fp : InferenceFile
            A file handler to an open inference file.
        max_iterations : {None, int}
            If results have not previously been written to the
            file, new datasets will be created. By default, the size of these
            datasets will be whatever the length of the sampler's chain is at
            this point. If you intend to run more iterations in the future,
            set this value to that size so that the array in the file will be
            large enough to accomodate future data.
        """
        self.write_metadata(fp)
        self.write_chain(fp, max_iterations=max_iterations)
        self.write_lnpost(fp, max_iterations=max_iterations)
        self.write_acceptance_fraction(fp)


samplers = {
    KombineSampler.name : KombineSampler,
    EmceeEnsembleSampler.name : EmceeEnsembleSampler,
}
