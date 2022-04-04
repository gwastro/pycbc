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

"""Provides abstract base class for all samplers."""


import time
from abc import (ABCMeta, abstractmethod)

from .base_hdf import BaseInferenceFile


class BaseSamplerFile(BaseInferenceFile, metaclass=ABCMeta):
    """Base HDF class for all samplers.

    This adds abstract methods ``write_resume_point`` and
    ``write_sampler_metadata`` to :py:class:`BaseInferenceFile`.
    """

    def write_run_start_time(self):
        """Writes the current (UNIX) time to the file.

        Times are stored as a list in the file's ``attrs``, with name
        ``run_start_time``. If the attrbute already exists, the current time
        is appended. Otherwise, the attribute will be created and time added.
        """
        attrname = "run_start_time"
        try:
            times = self.attrs[attrname].tolist()
        except KeyError:
            times = []
        times.append(time.time())
        self.attrs[attrname] = times

    @property
    def run_start_time(self):
        """The (UNIX) time pycbc inference began running.

        If the run resumed from a checkpoint, the time the last checkpoint
        started is reported.
        """
        return self.attrs['run_start_time'][-1]

    def write_run_end_time(self):
        """"Writes the curent (UNIX) time as the ``run_end_time`` attribute.
        """
        self.attrs["run_end_time"] = time.time()

    @property
    def run_end_time(self):
        """The (UNIX) time pycbc inference finished.
        """
        return self.attrs["run_end_time"]

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

    def update_checkpoint_history(self):
        """Writes a copy of relevant metadata to the file's checkpoint history.

        All data are written to ``sampler_info/checkpoint_history``. If the
        group does not exist yet, it will be created.

        This function writes the current time and the time since the last
        checkpoint to the file. It will also call
        :py:func:`_update_sampler_history` to write sampler-specific history.
        """
        path = '/'.join([self.sampler_group, 'checkpoint_history'])
        try:
            history = self[path]
        except KeyError:
            # assume history doesn't exist yet
            self.create_group(path)
            history = self[path]
        # write the checkpoint time
        current_time = time.time()
        self.write_data('checkpoint_time', current_time, path=path,
                        append=True)
        # get the amount of time since the last checkpoint
        checkpoint_times = history['checkpoint_time'][()]
        if len(checkpoint_times) == 1:
            # this is the first checkpoint, get the run time for comparison
            lasttime = self.run_start_time
        else:
            lasttime = checkpoint_times[-2]
            # if a resume happened since the last checkpoint, use the resume
            # time instad
            if lasttime < self.run_start_time:
                lasttime = self.run_start_time
        self.write_data('checkpoint_dt', current_time-lasttime, path=path,
                        append=True)
        # write any sampler-specific history
        self._update_sampler_history()

    def _update_sampler_history(self):
        """Writes sampler-specific history to the file.

        This function does nothing. Classes that inherit from it may override
        it to add any extra information they would like written. This is
        called by :py:func:`update_checkpoint_history`.
        """
        pass

    def validate(self):
        """Runs a validation test.

        This checks that a samples group exist, and that there are more than
        one sample stored to it.

        Returns
        -------
        bool :
            Whether or not the file is valid as a checkpoint file.
        """
        try:
            group = '{}/{}'.format(self.samples_group, self.variable_params[0])
            checkpoint_valid = self[group].size != 0
        except KeyError:
            checkpoint_valid = False
        return checkpoint_valid
