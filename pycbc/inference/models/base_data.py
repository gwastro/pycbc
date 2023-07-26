# Copyright (C) 2018  Collin Capano
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

"""Base classes for mofdels with data.
"""

import numpy
from abc import (ABCMeta, abstractmethod)
from .base import BaseModel


class BaseDataModel(BaseModel, metaclass=ABCMeta):
    r"""Base class for models that require data and a waveform generator.

    This adds propeties for the log of the likelihood that the data contain
    noise, ``lognl``, and the log likelihood ratio ``loglr``.

    Classes that inherit from this class must define ``_loglr`` and ``_lognl``
    functions, in addition to the ``_loglikelihood`` requirement inherited from
    ``BaseModel``.

    Parameters
    ----------
    variable_params : (tuple of) string(s)
        A tuple of parameter names that will be varied.
    data : dict
        A dictionary of data, in which the keys are the detector names and the
        values are the data.
    recalibration : dict of pycbc.calibration.Recalibrate, optional
        Dictionary of detectors -> recalibration class instances for
        recalibrating data.
    gates : dict of tuples, optional
        Dictionary of detectors -> tuples of specifying gate times. The
        sort of thing returned by `pycbc.gate.gates_from_cli`.
    injection_file : str, optional
        If an injection was added to the data, the name of the injection file
        used. If provided, the injection parameters will be written to
        file when ``write_metadata`` is called.

    \**kwargs :
        All other keyword arguments are passed to ``BaseModel``.


    See ``BaseModel`` for additional attributes and properties.
    """

    def __init__(self, variable_params, data, recalibration=None, gates=None,
                 injection_file=None, no_save_data=False, **kwargs):
        self._data = None
        self.data = data
        self.recalibration = recalibration
        self.no_save_data = no_save_data
        self.gates = gates
        self.injection_file = injection_file
        super(BaseDataModel, self).__init__(variable_params, **kwargs)

    @property
    def data(self):
        """dict: Dictionary mapping detector names to data."""
        return self._data

    @data.setter
    def data(self, data):
        """Store a copy of the data."""
        self._data = {det: d.copy() for (det, d) in data.items()}

    @property
    def _extra_stats(self):
        """Adds ``loglr`` and ``lognl`` to the ``default_stats``."""
        return ['loglr', 'lognl']

    @property
    def lognl(self):
        """The log likelihood of the model assuming the data is noise.

        This will initially try to return the ``current_stats.lognl``.
        If that raises an ``AttributeError``, will call `_lognl`` to
        calculate it and store it to ``current_stats``.
        """
        return self._trytoget('lognl', self._lognl)

    @abstractmethod
    def _lognl(self):
        """Low-level function that calculates the lognl."""
        pass

    @property
    def loglr(self):
        """The log likelihood ratio at the current parameters.

        This will initially try to return the ``current_stats.loglr``.
        If that raises an ``AttributeError``, will call `_loglr`` to
        calculate it and store it to ``current_stats``.
        """
        return self._trytoget('loglr', self._loglr, apply_transforms=True)

    @abstractmethod
    def _loglr(self):
        """Low-level function that calculates the loglr."""
        pass

    @property
    def logplr(self):
        """Returns the log of the prior-weighted likelihood ratio at the
        current parameter values.

        The logprior is calculated first. If the logprior returns ``-inf``
        (possibly indicating a non-physical point), then ``loglr`` is not
        called.
        """
        logp = self.logprior
        if logp == -numpy.inf:
            return logp
        else:
            return logp + self.loglr

    @property
    def detectors(self):
        """list: Returns the detectors used."""
        return list(self._data.keys())

    def write_metadata(self, fp, group=None):
        """Adds data to the metadata that's written.

        Parameters
        ----------
        fp : pycbc.inference.io.BaseInferenceFile instance
            The inference file to write to.
        group : str, optional
            If provided, the metadata will be written to the attrs specified
            by group, i.e., to ``fp[group].attrs``. Otherwise, metadata is
            written to the top-level attrs (``fp.attrs``).
        """
        super().write_metadata(fp, group=group)
        if not self.no_save_data:
            fp.write_stilde(self.data, group=group)
        # save injection parameters
        if self.injection_file is not None:
            fp.write_injections(self.injection_file, group=group)
