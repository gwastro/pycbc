# Module to hold recalibration stuff used during parameter estimation
# This will draw heavily from pycbc_adjust_strain and cal.py which
# are part of Chris Biwer's pycbc-cal repository

import numpy
from scipy.interpolate import UnivariateSpline
from pycbc.types import FrequencySeries

class Recalibrate:
    """ Class for adjusting time-varying calibration parameters.
    """

    def __init__(self, calib_dict=None):
        """ Initialize the class with transfer functions and calibration 
        parameters for a given epoch that starts at time t0.
        """

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
        self.c_res = self.c0 * (1 + 1.0j * self.freq / self.fc0) / init_detuning

    def update_c(self, fs=None, qinv=None, fc=None, kappa_c=1.0):
        detuning_term = self.freq**2 / (self.freq**2 - 1.0j * self.freq * fs * \
                                        qinv + fs**2)
        return self.c_res * kappa_c / (1 + 1.0j * self.freq/fc) * detuning_term

    def update_g(self, fs=None, qinv=None, fc=None, kappa_tst_re=1.0,
                 kappa_tst_im=0.0, kappa_pu_re=1.0, kappa_pu_im=0.0,
                 kappa_c=1.0):
        c = self.update_c(fs=fs, qinv=qinv, fc=fc, kappa_c=kappa_c)
        a_tst = self.a_tst0 * (kappa_tst_re + 1.0j * kappa_tst_im)
        a_pu = self.a_pu0 * (kappa_pu_re + 1.0j * kappa_pu_im)
        return c * self.d0 * (a_tst + a_pu)

    def update_r(self, fs=None, qinv=None, fc=None, kappa_c=1.0,
                 kappa_tst_re=1.0, kappa_tst_im=0.0, kappa_pu_re=1.0,
                 kappa_pu_im=0.0):
        c = self.update_c(fs=fs, qinv=qinv, fc=fc, kappa_c=kappa_c)
        g = self.update_g(fs=fs, qinv=qinv, fc=fc, kappa_c=kappa_c,
                          kappa_tst_re=kappa_tst_re, kappa_tst_im=kappa_tst_im,
                          kappa_pu_re=kappa_pu_re, kappa_pu_im=kappa_pu_im)
        return (1.0 + g) / c

    def adjust_strain(self, strain, params): #fs=None, qinv=None, fc=None, kappa_c=1.0,
                      #kappa_tst_re=1.0, kappa_tst_im=0.0, kappa_pu_re=1.0,
                      #kappa_pu_im=0.0):
        """Adjust the FrequencySeries strain
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
        freq_even = strain.sample_frequencies.numpy()
        k_even_sample = k_amp_off(freq_even) * \
                        numpy.exp(1.0j * k_phase_off(freq_even))
        strain_adjusted = FrequencySeries(strain.numpy() * \
                                          k_even_sample, delta_f=strain.delta_f)

        return strain_adjusted
