# Copyright (C) 2018 Colm Talbot
#
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
""" Functions for adding calibration factors to waveform templates.
"""

import numpy as np
from scipy.interpolate import UnivariateSpline
from abc import (ABCMeta, abstractmethod)


class Recalibrate(metaclass=ABCMeta):
    name = None

    def __init__(self, ifo_name):
        self.ifo_name = ifo_name
        self.params = dict()

    @abstractmethod
    def apply_calibration(self, waveform, frequencies=None):
        """Apply calibration model

        This method should be overwritten by subclasses

        Parameters
        ----------
        waveform : FrequencySeries, array, or list of either
            The waveform to be recalibrated. Several waveforms sharing
            frequencies may be given together, and are corrected on one
            evaluation of the model.
        frequencies : array, optional
            The frequencies the waveform is sampled at. Models that can be
            evaluated away from a waveform's own frequencies accept this,
            which lets a likelihood correct a waveform wherever it happens
            to have sampled it. Taken from the waveform if not given.

        Return
        ------
        waveform_adjusted : same as the input
            The recalibrated waveform, or waveforms.
        """
        return

    def set_params(self, prefix='recalib_', **params):
        """Pick the calibration parameters out of a set of sampling
        parameters and keep them for the next correction.

        Parameters
        ----------
        prefix: str
            Prefix for calibration parameter names
        params : dict
            Dictionary of sampling parameters which includes
            calibration parameters.
        """
        self.params.update({
            key[len(prefix):]: params[key]
            for key in params if prefix in key and self.ifo_name in key})

    def map_to_adjust(self, waveform, prefix='recalib_', frequencies=None,
                      **params):
        """Map an input dictionary of sampling parameters to the
        adjust_strain function by filtering the dictionary for the
        calibration parameters, then calling adjust_strain.

        Parameters
        ----------
        waveform : FrequencySeries, array, or list of either
            The waveform to be recalibrated.
        prefix: str
            Prefix for calibration parameter names
        frequencies : array, optional
            The frequencies the waveform is sampled at, if it is not one
            that carries them.
        params : dict
            Dictionary of sampling parameters which includes
            calibration parameters.
        Return
        ------
        waveform_adjusted : FrequencySeries or array
            The recalibrated waveform.
        """
        self.set_params(prefix=prefix, **params)

        return self.apply_calibration(waveform, frequencies=frequencies)

    @classmethod
    def from_config(cls, cp, ifo, section):
        """Read a config file to get calibration options and transfer
        functions which will be used to intialize the model.

        Parameters
        ----------
        cp : WorkflowConfigParser
            An open config file.
        ifo : string
            The detector (H1, L1) for which the calibration model will
            be loaded.
        section : string
            The section name in the config file from which to retrieve
            the calibration options.
        Return
        ------
        instance
            An instance of the class.
        """
        all_params = dict(cp.items(section))
        params = {key[len(ifo)+1:]: all_params[key]
                  for key in all_params if ifo.lower() in key}
        model = params.pop('model')
        params['ifo_name'] = ifo.lower()

        return all_models[model](**params)


class CubicSpline(Recalibrate):
    name = 'cubic_spline'

    def __init__(self, minimum_frequency, maximum_frequency, n_points,
                 ifo_name):
        """
        Cubic spline recalibration

        see https://dcc.ligo.org/LIGO-T1400682/public

        This assumes the spline points follow
        np.logspace(np.log(minimum_frequency), np.log(maximum_frequency),
                    n_points)

        Parameters
        ----------
        minimum_frequency: float
            minimum frequency of spline points
        maximum_frequency: float
            maximum frequency of spline points
        n_points: int
            number of spline points
        """
        Recalibrate.__init__(self, ifo_name=ifo_name)
        minimum_frequency = float(minimum_frequency)
        maximum_frequency = float(maximum_frequency)
        n_points = int(n_points)
        if n_points < 4:
            raise ValueError(
                'Use at least 4 spline points for calibration model')
        self.n_points = n_points
        self.spline_points = np.logspace(np.log10(minimum_frequency),
                                         np.log10(maximum_frequency), n_points)
        self._shape = None
        self._frequencies = None
        self._basis = None

    def waveform_frequencies(self, waveform):
        """The frequencies a waveform is sampled at, kept between calls.

        Asking a series for them builds a new array each time, which costs
        more than the correction does, and a model is normally handed the
        same shape of waveform for the length of a run.
        """
        shape = (len(waveform), float(waveform.delta_f))
        if shape != self._shape:
            self._shape = shape
            self._frequencies = waveform.sample_frequencies.numpy()
            self._basis = None
        return self._frequencies

    def spline_basis(self, frequencies):
        """The two splines reduced to one matrix, or None if that would
        not pay for itself.

        Over the range of control values a calibration correction covers,
        the spline is a linear map of them, so for a fixed set of
        frequencies it becomes a matrix product with the control values,
        which is what makes it cheap enough to apply inside a likelihood.
        The range is not unlimited and the caller checks it; see
        :py:meth:`apply_calibration`.

        Building the matrix costs more than one correction does, so it is
        only built once the same frequencies have been asked for twice.
        Frequencies that change must arrive as a new array, which is what
        rebinning and :py:meth:`waveform_frequencies` both produce.
        """
        if frequencies is not self._frequencies:
            self._shape = None
            self._frequencies = frequencies
            self._basis = None
        elif self._basis is None:
            self._basis = np.column_stack(
                [UnivariateSpline(self.spline_points, unit)(frequencies)
                 for unit in np.eye(self.n_points)])
        return self._basis

    def apply_calibration(self, waveform, frequencies=None):
        """Apply calibration model

        This applies cubic spline calibration to the waveform.

        Parameters
        ----------
        waveform : FrequencySeries, array, or list of either
            The waveform to be recalibrated. Several waveforms sharing
            frequencies may be given together, and are corrected on one
            evaluation of the spline.
        frequencies : array, optional
            The frequencies the waveform is sampled at. Taken from the
            waveform itself if not given.

        Return
        ------
        waveform_adjusted : same as the input
            The recalibrated waveform, or waveforms.
        """
        several = isinstance(waveform, (list, tuple))
        if frequencies is None:
            frequencies = self.waveform_frequencies(
                waveform[0] if several else waveform)

        amplitude_parameters = \
            [self.params['amplitude_{}_{}'.format(self.ifo_name, ii)]
             for ii in range(self.n_points)]
        phase_parameters = \
            [self.params['phase_{}_{}'.format(self.ifo_name, ii)]
             for ii in range(self.n_points)]

        basis = self.spline_basis(frequencies)
        if basis is not None and max(
                np.dot(amplitude_parameters, amplitude_parameters),
                np.dot(phase_parameters, phase_parameters)) > self.n_points:
            # UnivariateSpline smooths rather than interpolates, and is a
            # linear map of the control values only while its residual
            # constraint does not bind. It cannot bind while their sum of
            # squares stays under scipy's default s, which is n_points.
            # A calibration correction is a few percent, some 300 times
            # below that at ten control points, so this is not the path
            # anything real takes; past the bound the matrix and the
            # spline part company by order unity rather than gracefully,
            # so it is checked rather than assumed.
            basis = None
        if basis is None:
            delta_amplitude = UnivariateSpline(
                self.spline_points, amplitude_parameters)(frequencies)
            delta_phase = UnivariateSpline(
                self.spline_points, phase_parameters)(frequencies)
        else:
            delta_amplitude = basis.dot(amplitude_parameters)
            delta_phase = basis.dot(phase_parameters)

        factor = ((1.0 + delta_amplitude) * (2.0 + 1j * delta_phase)
                  / (2.0 - 1j * delta_phase))

        if several:
            return type(waveform)(each * factor for each in waveform)
        return waveform * factor


all_models = {
    CubicSpline.name: CubicSpline
}
