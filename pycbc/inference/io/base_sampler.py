# Copyright (C) 2019 Collin Capano
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# self.option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

from __future__ import absolute_import

from .base_hdf import BaseInferenceFile
from .posterior import PosteriorFile

"""Provides abstract base class for all samplers."""

class BaseSamplerFile(BaseInferenceFile):
    """Base HDF class for all samplers.
    
    This adds abstract methods ``write_resume_point``,
    ``write_sampler_metadata``, and ``read_posterior_samples`` to
    :py:class:`BaseInferenceFile`. Also adds a ``write_posterior`` method that
    extracts posterior samples and writes them to a ``PosteriorFile``.
    """

    @abstractmethod
    def write_resume_point(self):
        """Should write the point that a sampler starts up.

        How the resume point is indexed is up to the sampler. For example,
        MCMC samplers use the number of iterations that are stored in the
        checkpoint file.
        """
        pass

    @abstractmethod
    def write_sampler_metadata(self, sampler):
        """This should write the given sampler's metadata to the file.

        This should also include the model's metadata.
        """
        pass

    @abstractmethod
    def read_posterior_samples(self, parameters):
        """Reads posterior samples from the file.

        This is used by ``write_posterior``. It should be a wrapper around
        ``read_samples`` that automatically sets any arguments needed to
        extract a posterior from the raw samples.  Since this is
        sampler-specific, this must be implemented by each sampler file class.

        Parameters
        ----------
        parameters : list of str
            The names of the parameters to read.

        Returns
        -------
        FieldArray :
            The posterior samples, as a 1D ``FieldArray``.
        """
        pass

    def write_posterior(self, filename, parameterdict=None, skip_groups=None):
        """Write posterior samples to a ``PosteriorFile``.

        Parameters
        ----------
        filename : str
            Name of output file to store posterior.
        parameterdict : dict, optional
            Dictionary mapping parameters to read to the names they should
            be called in the posterior file. The keys may be functions of the
            parameters in the samples group. If none provided, will just load
            all of the data sets in the ``samples`` group, writing them out
            with the same name.
        skip_groups : {None, 'all', list of str}, optional
            Do not write the given list of groups to the posterior file. The
            ``samples`` group is always written. This can be used to specify
            wether some or all the other groups, such as ``sampler_info``
            should be written too. May provide ``None``, in which case all of
            the other groups are written; a list of strings given the groups to
            skip; or ``'all'``, in which case only the samples group will be
            written. Default is ``None``.
        """
        # get the posterior samples
        if parameterdict is None:
            # load everything in the samples group
            parameterdict = {p: p for p in self[self.samples_group].keys()}
        samples = self.read_posterior_samples(parameterdict.keys())
        # samples is currently a FieldArray; convert to dict in which the keys
        # are the values of the parameter dict
        samples = {parameterdict[p]: samples[p] for p in parameterdict}
        # now write
        f = PosteriorFile(filename, 'w')
        f.write_samples(samples)
        # Preserve samples group metadata
        for key, val in self[self.samples_group].attrs.items():
            f[f.samples_group].attrs[key] = val
        # Preserve top-level metadata
        for key in self.attrs:
            f.attrs[key] = self.attrs[key]
        # write the other groups
        if skip_groups == 'all':
            skip_groups = [group for group in self.keys()
                           if group != self.samples_group]
        self.copy_info(f, ignore=skip_groups)
        f.close()
