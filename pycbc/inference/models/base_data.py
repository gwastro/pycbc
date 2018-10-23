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

"""Base classes for models with data.
"""

import numpy
import logging
from abc import (ABCMeta, abstractmethod)

from pycbc import transforms
from pycbc.waveform import generator

from .base import BaseModel


class BaseDataModel(BaseModel):
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
    waveform_generator : generator class
        A generator class that creates waveforms.
    waveform_transforms : list, optional
        List of transforms to use to go from the variable args to parameters
        understood by the waveform generator.

    \**kwargs :
        All other keyword arguments are passed to ``BaseModel``.

    Attributes
    ----------
    waveform_generator : dict
        The waveform generator that the class was initialized with.
    data : dict
        The data that the class was initialized with.

    Properties
    ----------
    lognl :
        Returns the log likelihood of the noise.
    loglr :
        Returns the log of the likelihood ratio.
    logplr :
        Returns the log of the prior-weighted likelihood ratio.

    See ``BaseModel`` for additional attributes and properties.
    """
    __metaclass__ = ABCMeta

    def __init__(self, variable_params, data, waveform_generator,
                 waveform_transforms=None, **kwargs):
        # we'll store a copy of the data
        self._data = {ifo: d.copy() for (ifo, d) in data.items()}
        self._waveform_generator = waveform_generator
        self._waveform_transforms = waveform_transforms
        super(BaseDataModel, self).__init__(
            variable_params, **kwargs)

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
        return self._trytoget('loglr', self._loglr)

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
    def waveform_generator(self):
        """Returns the waveform generator that was set."""
        return self._waveform_generator

    @property
    def data(self):
        """Returns the data that was set."""
        return self._data

    @property
    def detectors(self):
        """Returns the detectors used."""
        return self._data.keys()

    def _transform_params(self, **params):
        """Adds waveform transforms to parent's ``_transform_params``."""
        params = super(BaseDataModel, self)._transform_params(**params)
        # apply waveform transforms
        if self._waveform_transforms is not None:
            params = transforms.apply_transforms(params,
                                                 self._waveform_transforms,
                                                 inverse=False)
        return params

    @classmethod
    def _init_args_from_config(cls, cp):
        """Adds loading waveform_transforms to parent function.

        For details on parameters, see ``from_config``.
        """
        args = super(BaseDataModel, cls)._init_args_from_config(cp)
        # add waveform transforms to the arguments
        if any(cp.get_subsections('waveform_transforms')):
            logging.info("Loading waveform transforms")
            args['waveform_transforms'] = \
                transforms.read_transforms_from_config(
                    cp, 'waveform_transforms')
        return args

    @classmethod
    def from_config(cls, cp, data=None, delta_f=None, delta_t=None,
                    gates=None, recalibration=None,
                    **kwargs):
        """Initializes an instance of this class from the given config file.

        Parameters
        ----------
        cp : WorkflowConfigParser
            Config file parser to read.
        data : dict
            A dictionary of data, in which the keys are the detector names and
            the values are the data. This is not retrieved from the config
            file, and so must be provided.
        delta_f : float
            The frequency spacing of the data; needed for waveform generation.
        delta_t : float
            The time spacing of the data; needed for time-domain waveform
            generators.
        recalibration : dict of pycbc.calibration.Recalibrate, optional
            Dictionary of detectors -> recalibration class instances for
            recalibrating data.
        gates : dict of tuples, optional
            Dictionary of detectors -> tuples of specifying gate times. The
            sort of thing returned by `pycbc.gate.gates_from_cli`.
        \**kwargs :
            All additional keyword arguments are passed to the class. Any
            provided keyword will over ride what is in the config file.
        """
        if data is None:
            raise ValueError("must provide data")
        args = cls._init_args_from_config(cp)
        args['data'] = data
        args.update(kwargs)

        variable_params = args['variable_params']
        try:
            static_params = args['static_params']
        except KeyError:
            static_params = {}

        # set up waveform generator
        try:
            approximant = static_params['approximant']
        except KeyError:
            raise ValueError("no approximant provided in the static args")
        generator_function = generator.select_waveform_generator(approximant)
        waveform_generator = generator.FDomainDetFrameGenerator(
            generator_function, epoch=data.values()[0].start_time,
            variable_args=variable_params, detectors=data.keys(),
            delta_f=delta_f, delta_t=delta_t,
            recalib=recalibration, gates=gates,
            **static_params)
        args['waveform_generator'] = waveform_generator

        return cls(**args)

    def write_metadata(self, fp):
        """Adds data to the metadata that's written.

        Parameters
        ----------
        fp : pycbc.inference.io.BaseInferenceFile instance
            The inference file to write to.
        """
        super(BaseDataModel, self).write_metadata(fp)
        fp.write_stilde(self.data)
