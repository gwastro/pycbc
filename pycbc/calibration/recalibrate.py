""" Classes and functions for adjusting strain data.
"""
# Copyright (C) 2015 Ben Lackey, Christopher M. Biwer, Daniel Finstad
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

import numpy
from scipy.interpolate import UnivariateSpline
from pycbc.types import FrequencySeries

class Recalibrate(object):
    """ Class for adjusting time-varying calibration parameters of given
    strain data.

    Parameters
    ----------
    strain : FrequencySeries
        The strain to be adjusted.
    calib_dict : dict
        Dictionary of calibration parameters and their values at reference
        time t0
    """

    def __init__(self, strain, calib_dict=None):

        self.strain = strain
        self.freq = numpy.real(calib_dict["freq"])
        self.c0 = calib_dict["c0"]
        self.d0 = calib_dict["d0"]
        self.a_tst0 = calib_dict["a_tst0"]
        self.a_pu0 = calib_dict["a_pu0"]
        self.fc0 = float(calib_dict["fc0"])
        self.fs0 = float(calib_dict["fs0"])
        self.qinv0 = float(calib_dict["qinv0"])

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

    def adjust_strain(self, params):
        """Adjust the FrequencySeries strain by changing the time-dependent
        calibration parameters kappa_c(t), kappa_a(t), f_c(t), fs, and qinv.
    
        Parameters
        ----------
        fc : float
            Coupled-cavity (CC) pole at time t.
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
        fs = params["calib_fs"] if "calib_fs" in params else self.fs0
        qinv = params["calib_qinv"] if "calib_qinv" in params else self.qinv0
        fc = self.fc0+params["calib_deltafc"] if "calib_deltafc" in params \
             else self.fc0
        kappa_c = params["calib_kappa_c"] if "calib_kappa_c" in params else 1.0
        kappa_tst_re = params["calib_kappa_tst_re"] if "calib_kappa_tst_re" in \
                       params else 1.0
        kappa_tst_im = params["calib_kappa_tst_im"] if "calib_kappa_tst_im" in \
                       params else 0.0
        kappa_pu_re = params["calib_kappa_pu_re"] if "calib_kappa_pu_re" in \
                      params else 1.0
        kappa_pu_im = params["calib_kappa_pu_im"] if "calib_kappa_pu_im" in \
                      params else 0.0

        r_adjusted = self.update_r(fs=fs, qinv=qinv, fc=fc, kappa_c=kappa_c,
                                   kappa_tst_re=kappa_tst_re,
                                   kappa_tst_im=kappa_tst_im,
                                   kappa_pu_re=kappa_pu_re,
                                   kappa_pu_im=kappa_pu_im)

        # calculate error function
        k = r_adjusted / self.r0

        # decompose into amplitude and unwrapped phase
        k_amp = numpy.abs(k)
        k_phase = numpy.unwrap(numpy.angle(k))

        # convert to FrequencySeries by interpolating then resampling
        order = 1
        k_amp_off = UnivariateSpline(self.freq, k_amp, k=order, s=0)
        k_phase_off = UnivariateSpline(self.freq, k_phase, k=order, s=0)
        freq_even = self.strain.sample_frequencies.numpy()
        k_even_sample = k_amp_off(freq_even) * \
                        numpy.exp(1.0j * k_phase_off(freq_even))
        strain_adjusted = FrequencySeries(self.strain.numpy() * \
                                          k_even_sample,
                                          delta_f=self.strain.delta_f)

        return strain_adjusted
