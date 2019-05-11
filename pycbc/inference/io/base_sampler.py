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

    def write_posterior(self, filename, samples, labels=None,
                        skip_groups=None):
        """Write posterior samples to a ``PosteriorFile``.

        Parameters
        ----------
        filename : str
            Name of output file to store posterior.
        samples : FieldArray
            Samples to write. Must be a 1D array.
        labels : dict, optional
            Dictionary mapping parameters to the names they should
            be called in the posterior file. The keys may be functions of the
            fields in ``samples``. If none provided, will write all of the
            fields in ``samples``, with the same name.
        skip_groups : {None, 'all', list of str}, optional
            Do not write the given list of groups to the posterior file. The
            ``samples`` group is always written. This can be used to specify
            wether some or all the other groups, such as ``sampler_info``
            should be written too. May provide ``None``, in which case all of
            the other groups are written; a list of strings given the groups to
            skip; or ``'all'``, in which case only the samples group will be
            written. Default is ``None``.
        """
        if not samples.ndim == 1:
            raise ValueError("samples must be a 1D array")
        # get the posterior samples
        if labels is None:
            # load everything in the samples
            labels = {p: p for p in samples.fieldnames}
        # convert samples to a dict in which the keys
        # are the values of the parameter dict
        samples = {labels[p]: samples[p] for p in labels}
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
