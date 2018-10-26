""" Classes and functions for adjusting strain data.
"""
# Copyright (C) 2015 Ben Lackey, Christopher M. Biwer,
#                    Daniel Finstad, Colm Talbot, Alex Nitz
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

from abc import (ABCMeta, abstractmethod)

import numpy as np
from scipy.interpolate import UnivariateSpline
from pycbc.types import FrequencySeries

class Recalibrate(object):
    """ Base class for modifying calibration """

    __metaclass__ = ABCMeta
    name = None

    def __init__(self, ifo_name):
        self.ifo_name = ifo_name
        self.params = dict()

    @abstractmethod
    def apply_calibration(self, strain):
        """Apply calibration model

        This method should be overwritten by subclasses

        Parameters
        ----------
        strain : FrequencySeries
            The strain to be recalibrated.

        Return
        ------
        strain_adjusted : FrequencySeries
            The recalibrated strain.
        """
        return

    def map_to_adjust(self, strain, prefix='recalib_', **params):
        """Map an input dictionary of sampling parameters to the
        adjust_strain function by filtering the dictionary for the
        calibration parameters, then calling adjust_strain.

        Parameters
        ----------
        strain : FrequencySeries
            The strain to be recalibrated.
        prefix: str
            Prefix for calibration parameter names
        params : dict
            Dictionary of sampling parameters which includes
            calibration parameters.
        Return
        ------
        strain_adjusted : FrequencySeries
            The recalibrated strain.
        """

        self.params.update({
            key[len(prefix):]: params[key]
            for key in params if prefix in key and self.ifo_name in key})

        strain_adjusted = self.apply_calibration(strain)

        return strain_adjusted

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
        params = {key: params[key] for key in params}
        params.pop('model')
        params['ifo_name'] = ifo.lower()

        return cls(**params)


class CubicSpline(Recalibrate):
    """Cubic spline recalibration

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
    name = 'cubic_spline'

    def __init__(self, minimum_frequency, maximum_frequency, n_points,
                 ifo_name):
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

    def apply_calibration(self, strain):
        """Apply calibration model

        This applies cubic spline calibration to the strain.

        Parameters
        ----------
        strain : FrequencySeries
            The strain to be recalibrated.

        Return
        ------
        strain_adjusted : FrequencySeries
            The recalibrated strain.
        """
        amplitude_parameters =\
            [self.params['amplitude_{}_{}'.format(self.ifo_name, ii)]
             for ii in range(self.n_points)]
        amplitude_spline = UnivariateSpline(self.spline_points,
                                            amplitude_parameters)
        delta_amplitude = amplitude_spline(strain.sample_frequencies.numpy())

        phase_parameters =\
            [self.params['phase_{}_{}'.format(self.ifo_name, ii)]
             for ii in range(self.n_points)]
        phase_spline = UnivariateSpline(self.spline_points, phase_parameters)
        delta_phase = phase_spline(strain.sample_frequencies.numpy())

        strain_adjusted = strain * (1.0 + delta_amplitude)\
            * (2.0 + 1j * delta_phase) / (2.0 - 1j * delta_phase)

        return strain_adjusted


class PhysicalModel(object):
    """ Class for adjusting time-varying calibration parameters of given
    strain data.

    Attributes
    ----------
    name : 'physical_model'
        The name of this calibration model.

    Parameters
    ----------
    strain : FrequencySeries
        The strain to be adjusted.
    freq : array
        The frequencies corresponding to the values of c0, d0, a0 in Hertz.
    fc0 : float
        Coupled-cavity (CC) pole at time t0, when c0=c(t0) and a0=a(t0) are
        measured.
    c0 : array
        Initial sensing function at t0 for the frequencies.
    d0 : array
        Digital filter for the frequencies.
    a_tst0 : array
        Initial actuation function for the test mass at t0 for the
        frequencies.
    a_pu0 : array
        Initial actuation function for the penultimate mass at t0 for the
        frequencies.
    fs0 : float
        Initial spring frequency at t0 for the signal recycling cavity.
    qinv0 : float
        Initial inverse quality factor at t0 for the signal recycling
        cavity.
    """

    name = 'physical_model'
    def __init__(self, freq=None, fc0=None, c0=None, d0=None,
                 a_tst0=None, a_pu0=None, fs0=None, qinv0=None):
        self.freq = np.real(freq)
        self.c0 = c0
        self.d0 = d0
        self.a_tst0 = a_tst0
        self.a_pu0 = a_pu0
        self.fc0 = float(fc0)
        self.fs0 = float(fs0)
        self.qinv0 = float(qinv0)

        # initial detuning at time t0
        init_detuning = self.freq**2 / (self.freq**2 - 1.0j * self.freq * \
                                        self.fs0 * self.qinv0 + self.fs0**2)

        # initial open loop gain
        self.g0 = self.c0 * self.d0 * (self.a_tst0 + self.a_pu0)

        # initial response function
        self.r0 = (1.0 + self.g0) / self.c0

        # residual of c0 after factoring out the coupled cavity pole fc0
        self.c_res = self.c0 * (1 + 1.0j * self.freq / self.fc0)/init_detuning

    def update_c(self, fs=None, qinv=None, fc=None, kappa_c=1.0):
        """ Calculate the sensing function c(f,t) given the new parameters
        kappa_c(t), kappa_a(t), f_c(t), fs, and qinv.

        Parameters
        ----------
        fc : float
            Coupled-cavity (CC) pole at time t.
        kappa_c : float
            Scalar correction factor for sensing function at time t.
        fs : float
            Spring frequency for signal recycling cavity.
        qinv : float
            Inverse quality factor for signal recycling cavity.

        Returns
        -------
        c : numpy.array
            The new sensing function c(f,t).
        """
        detuning_term = self.freq**2 / (self.freq**2 - 1.0j *self.freq*fs * \
                                        qinv + fs**2)
        return self.c_res * kappa_c / (1 + 1.0j * self.freq/fc)*detuning_term

    def update_g(self, fs=None, qinv=None, fc=None, kappa_tst_re=1.0,
                 kappa_tst_im=0.0, kappa_pu_re=1.0, kappa_pu_im=0.0,
                 kappa_c=1.0):
        """ Calculate the open loop gain g(f,t) given the new parameters
        kappa_c(t), kappa_a(t), f_c(t), fs, and qinv.

        Parameters
        ----------
        fc : float
            Coupled-cavity (CC) pole at time t.
        kappa_c : float
            Scalar correction factor for sensing function c at time t.
        kappa_tst_re : float
            Real part of scalar correction factor for actuation function
            a_tst0 at time t.
        kappa_pu_re : float
            Real part of scalar correction factor for actuation function
            a_pu0 at time t.
        kappa_tst_im : float
            Imaginary part of scalar correction factor for actuation function
            a_tst0 at time t.
        kappa_pu_im : float
            Imaginary part of scalar correction factor for actuation function
            a_pu0 at time t.
        fs : float
            Spring frequency for signal recycling cavity.
        qinv : float
            Inverse quality factor for signal recycling cavity.

        Returns
        -------
        g : numpy.array
            The new open loop gain g(f,t).
        """
        c = self.update_c(fs=fs, qinv=qinv, fc=fc, kappa_c=kappa_c)
        a_tst = self.a_tst0 * (kappa_tst_re + 1.0j * kappa_tst_im)
        a_pu = self.a_pu0 * (kappa_pu_re + 1.0j * kappa_pu_im)
        return c * self.d0 * (a_tst + a_pu)

    def update_r(self, fs=None, qinv=None, fc=None, kappa_c=1.0,
                 kappa_tst_re=1.0, kappa_tst_im=0.0, kappa_pu_re=1.0,
                 kappa_pu_im=0.0):
        """ Calculate the response function R(f,t) given the new parameters
        kappa_c(t), kappa_a(t), f_c(t), fs, and qinv.

        Parameters
        ----------
        fc : float
            Coupled-cavity (CC) pole at time t.
        kappa_c : float
            Scalar correction factor for sensing function c at time t.
        kappa_tst_re : float
            Real part of scalar correction factor for actuation function
            a_tst0 at time t.
        kappa_pu_re : float
            Real part of scalar correction factor for actuation function
            a_pu0 at time t.
        kappa_tst_im : float
            Imaginary part of scalar correction factor for actuation function
            a_tst0 at time t.
        kappa_pu_im : float
            Imaginary part of scalar correction factor for actuation function
            a_pu0 at time t.
        fs : float
            Spring frequency for signal recycling cavity.
        qinv : float
            Inverse quality factor for signal recycling cavity.

        Returns
        -------
        r : numpy.array
            The new response function r(f,t).
        """
        c = self.update_c(fs=fs, qinv=qinv, fc=fc, kappa_c=kappa_c)
        g = self.update_g(fs=fs, qinv=qinv, fc=fc, kappa_c=kappa_c,
                          kappa_tst_re=kappa_tst_re,
                          kappa_tst_im=kappa_tst_im,
                          kappa_pu_re=kappa_pu_re, kappa_pu_im=kappa_pu_im)
        return (1.0 + g) / c

    def adjust_strain(self, strain, delta_fs=None, delta_qinv=None,
                      delta_fc=None, kappa_c=1.0, kappa_tst_re=1.0,
                      kappa_tst_im=0.0, kappa_pu_re=1.0, kappa_pu_im=0.0):
        """Adjust the FrequencySeries strain by changing the time-dependent
        calibration parameters kappa_c(t), kappa_a(t), f_c(t), fs, and qinv.

        Parameters
        ----------
        strain : FrequencySeries
            The strain data to be adjusted.
        delta_fc : float
            Change in coupled-cavity (CC) pole at time t.
        kappa_c : float
            Scalar correction factor for sensing function c0 at time t.
        kappa_tst_re : float
            Real part of scalar correction factor for actuation function
            A_{tst0} at time t.
        kappa_tst_im : float
            Imaginary part of scalar correction factor for actuation function
            A_tst0 at time t.
        kappa_pu_re : float
            Real part of scalar correction factor for actuation function
            A_{pu0} at time t.
        kappa_pu_im : float
            Imaginary part of scalar correction factor for actuation function
            A_{pu0} at time t.
        fs : float
            Spring frequency for signal recycling cavity.
        qinv : float
            Inverse quality factor for signal recycling cavity.

        Returns
        -------
        strain_adjusted : FrequencySeries
            The adjusted strain.
        """
        fc = self.fc0 + delta_fc if delta_fc else self.fc0
        fs = self.fs0 + delta_fs if delta_fs else self.fs0
        qinv = self.qinv0 + delta_qinv if delta_qinv else self.qinv0

        # calculate adjusted response function
        r_adjusted = self.update_r(fs=fs, qinv=qinv, fc=fc, kappa_c=kappa_c,
                                   kappa_tst_re=kappa_tst_re,
                                   kappa_tst_im=kappa_tst_im,
                                   kappa_pu_re=kappa_pu_re,
                                   kappa_pu_im=kappa_pu_im)

        # calculate error function
        k = r_adjusted / self.r0

        # decompose into amplitude and unwrapped phase
        k_amp = np.abs(k)
        k_phase = np.unwrap(np.angle(k))

        # convert to FrequencySeries by interpolating then resampling
        order = 1
        k_amp_off = UnivariateSpline(self.freq, k_amp, k=order, s=0)
        k_phase_off = UnivariateSpline(self.freq, k_phase, k=order, s=0)
        freq_even = strain.sample_frequencies.numpy()
        k_even_sample = k_amp_off(freq_even) * \
                        np.exp(1.0j * k_phase_off(freq_even))
        strain_adjusted = FrequencySeries(strain.numpy() * \
                                          k_even_sample,
                                          delta_f=strain.delta_f)

        return strain_adjusted

    @classmethod
    def tf_from_file(cls, path, delimiter=" "):
        """Convert the contents of a file with the columns
        [freq, real(h), imag(h)] to a numpy.array with columns
        [freq, real(h)+j*imag(h)].

        Parameters
        ----------
        path : string
        delimiter : {" ", string}

        Return
        ------
        numpy.array
        """
        data = np.loadtxt(path, delimiter=delimiter)
        freq = data[:, 0]
        h = data[:, 1] + 1.0j * data[:, 2]
        return np.array([freq, h]).transpose()

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
            An instance of the Recalibrate class.
        """
        # read transfer functions
        tfs = []
        tf_names = ["a-tst", "a-pu", "c", "d"]
        for tag in ['-'.join([ifo, "transfer-function", name])
                    for name in tf_names]:
            tf_path = cp.get_opt_tag(section, tag, None)
            tfs.append(cls.tf_from_file(tf_path))
        a_tst0 = tfs[0][:, 1]
        a_pu0 = tfs[1][:, 1]
        c0 = tfs[2][:, 1]
        d0 = tfs[3][:, 1]
        freq = tfs[0][:, 0]

        # if upper stage actuation is included, read that in and add it
        # to a_pu0
        uim_tag = '-'.join([ifo, 'transfer-function-a-uim'])
        if cp.has_option(section, uim_tag):
            tf_path = cp.get_opt_tag(section, uim_tag, None)
            a_pu0 += cls.tf_from_file(tf_path)[:, 1]

        # read fc0, fs0, and qinv0
        fc0 = cp.get_opt_tag(section, '-'.join([ifo, "fc0"]), None)
        fs0 = cp.get_opt_tag(section, '-'.join([ifo, "fs0"]), None)
        qinv0 = cp.get_opt_tag(section, '-'.join([ifo, "qinv0"]), None)

        return cls(freq=freq, fc0=fc0, c0=c0, d0=d0, a_tst0=a_tst0,
                   a_pu0=a_pu0, fs0=fs0, qinv0=qinv0)

    def map_to_adjust(self, strain, **params):
        """Map an input dictionary of sampling parameters to the
        adjust_strain function by filtering the dictionary for the
        calibration parameters, then calling adjust_strain.

        Parameters
        ----------
        strain : FrequencySeries
            The strain to be recalibrated.
        params : dict
            Dictionary of sampling parameters which includes
            calibration parameters.

        Return
        ------
        strain_adjusted : FrequencySeries
            The recalibrated strain.
        """
        # calibration param names
        arg_names = ['delta_fs', 'delta_fc', 'delta_qinv', 'kappa_c',
                     'kappa_tst_re', 'kappa_tst_im', 'kappa_pu_re',
                     'kappa_pu_im']
        # calibration param labels as they exist in config files
        arg_labels = [''.join(['calib_', name]) for name in arg_names]
        # default values for calibration params
        default_values = [0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0]
        # make list of calibration param values
        calib_args = []
        for arg, val in zip(arg_labels, default_values):
            if arg in params:
                calib_args.append(params[arg])
            else:
                calib_args.append(val)
        # adjust the strain using calibration param values
        strain_adjusted = self.adjust_strain(strain, delta_fs=calib_args[0],
                              delta_fc=calib_args[1], delta_qinv=calib_args[2],
                              kappa_c=calib_args[3],
                              kappa_tst_re=calib_args[4],
                              kappa_tst_im=calib_args[5],
                              kappa_pu_re=calib_args[6],
                              kappa_pu_im=calib_args[7])
        return strain_adjusted
