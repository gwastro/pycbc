import numpy as np
import numpy
import math
import os
from typing import Optional, Tuple
import lal
from lal import LIGOTimeGPS
from scipy import signal

import pycbc.psd
import pycbc.waveform
import pycbc.strain
import pycbc.noise

# The GSTLal FIR minimal phase routine
class PSDFirKernel(object):
    def __init__(self):
        self.revplan = None
        self.fwdplan = None
        self.target_phase = None
        self.target_phase_mask = None

    def set_phase(
        self,
        psd: lal.REAL8FrequencySeries,
        f_low: float = 10.0,
        m1: float = 1.4,
        m2: float = 1.4
    ) -> None:
        """
        Compute the phase response of zero-latency whitening filter
        given a reference PSD.
        """
        kernel, latency, sample_rate = self.psd_to_linear_phase_whitening_fir_kernel(psd)
        kernel, phase = self.linear_phase_fir_kernel_to_minimum_phase_whitening_fir_kernel(kernel, sample_rate)

        # get merger model for SNR = 1.
        f_psd = psd.f0 + numpy.arange(len(psd.data.data)) * psd.deltaF
        horizon_distance = HorizonDistance(f_low, f_psd[-1], psd.deltaF, m1, m2)
        f_model, model= horizon_distance(psd, 1.)[1]

        # find the range of frequency bins covered by the merger
        # model
        kmin, kmax = f_psd.searchsorted(f_model[0]), f_psd.searchsorted(f_model[-1]) + 1

        # compute SNR=1 model's (d SNR^2 / df) spectral density
        unit_snr2_density = numpy.zeros_like(phase)
        unit_snr2_density[kmin:kmax] = model / psd.data.data[kmin:kmax]

        # integrate across each frequency bin, converting to
        # snr^2/bin.  NOTE:  this step is here for the record, but
        # is commented out because it has no effect on the result
        # given the renormalization that occurs next.
        #unit_snr2_density *= psd.deltaF

        # take 16th root, then normalize so max=1.  why?  I don't
        # know, just feels good, on the whole.
        unit_snr2_density = unit_snr2_density**(1./16)
        unit_snr2_density /= unit_snr2_density.max()

        # record phase vector and SNR^2 density vector
        self.target_phase = phase
        self.target_phase_mask = unit_snr2_density

    def psd_to_linear_phase_whitening_fir_kernel(
        self,
        psd: lal.REAL8FrequencySeries,
        invert: Optional[bool] = True,
        nyquist: Optional[float] = None
    ) -> Tuple[numpy.ndarray, int, int]:
        """
        Compute an acausal finite impulse-response filter kernel
        from a power spectral density conforming to the LAL
        normalization convention, such that if colored Gaussian
        random noise with the given PSD is fed into an FIR filter
        using the kernel the filter's output will be zero-mean
        unit-variance Gaussian random noise.  The PSD must be
        provided as a lal.REAL8FrequencySeries object.

        The phase response of this filter is 0, just like whitening
        done in the frequency domain.

        Args:
            psd:
                lal.REAL8FrequencySeries, the reference PSD
            invert:
                bool, default true, whether to invert the kernel
            nyquist:
                float, disabled by default, whether to change
                the Nyquist frequency.

        Returns:
            Tuple[numpy.ndarray, int, int], the kernel, latency,
            sample rate pair.  The kernel is a numpy array containing
            the filter kernel, the latency is the filter latency in
            samples and the sample rate is in Hz.  The kernel and
            latency can be used, for example, with gstreamer's stock
            audiofirfilter element.

        """
        #
        # this could be relaxed with some work
        #

        assert psd.f0 == 0.0

        #
        # extract the PSD bins and determine sample rate for kernel
        #

        data = psd.data.data / 2
        sample_rate =  2 * (psd.f0 + (len(data) - 1) * psd.deltaF)

        #
        # remove LAL normalization
        #

        data *= sample_rate

        #
        # change Nyquist frequency if requested.  round to nearest
        # available bin
        #

        if nyquist is not None:
            i = int(round((nyquist - psd.f0) / psd.deltaF))
            assert i < len(data)
            data = data[:i + 1]
            sample_rate = 2 * int(round(psd.f0 + (len(data) - 1) * psd.deltaF))

        #
        # compute the FIR kernel.  it always has an odd number of
        # samples and no DC offset.
        #

        data[0] = data[-1] = 0.0
        if invert:
            data_nonzeros = (data != 0.)
            data[data_nonzeros] = 1./data[data_nonzeros]
        # repack data:  data[0], data[1], 0, data[2], 0, ....
        tmp = numpy.zeros((2 * len(data) - 1,), dtype = data.dtype)
        tmp[len(data)-1:] = data
        #tmp[:len(data)] = data
        data = tmp

        kernel_fseries = lal.CreateCOMPLEX16FrequencySeries(
            name = "double sided psd",
            epoch = LIGOTimeGPS(0),
            f0 = 0.0,
            deltaF = psd.deltaF,
            length = len(data),
            sampleUnits = lal.Unit("strain s")
        )

        kernel_tseries = lal.CreateCOMPLEX16TimeSeries(
            name = "timeseries of whitening kernel",
            epoch = LIGOTimeGPS(0.),
            f0 = 0.,
            deltaT = 1.0 / sample_rate,
            length = len(data),
            sampleUnits = lal.Unit("strain")
        )

        # FIXME check for change in length
        if self.revplan is None:
            self.revplan = lal.CreateReverseCOMPLEX16FFTPlan(len(data), 1)

        kernel_fseries.data.data = numpy.sqrt(data) + 0.j
        lal.COMPLEX16FreqTimeFFT(kernel_tseries, kernel_fseries, self.revplan)
        kernel = kernel_tseries.data.data.real
        kernel = numpy.roll(kernel, (len(data) - 1) // 2) / sample_rate * 2

        #
        # apply a Tukey window whose flat bit is 50% of the kernel.
        # preserve the FIR kernel's square magnitude
        #

        norm_before = numpy.dot(kernel, kernel)
        kernel *= lal.CreateTukeyREAL8Window(len(data), .5).data.data
        kernel *= math.sqrt(norm_before / numpy.dot(kernel, kernel))

        #
        # the kernel's latency
        #

        latency = (len(data) - 1) // 2

        #
        # done
        #

        return kernel, latency, sample_rate

    def linear_phase_fir_kernel_to_minimum_phase_whitening_fir_kernel(
        self,
        linear_phase_kernel: numpy.ndarray,
        sample_rate: int
    ) -> Tuple[numpy.ndarray, numpy.ndarray]:
        """
        Compute the minimum-phase response filter (zero latency)
        associated with a linear-phase response filter (latency
        equal to half the filter length).

        From "Design of Optimal Minimum-Phase Digital FIR Filters
        Using Discrete Hilbert Transforms", IEEE Trans. Signal
        Processing, vol. 48, pp. 1491-1495, May 2000.

        Args:
            linear_phase_kernel:
                numpy.ndarray, the kernel to compute the minimum-phase kernel with
            sample_rate:
                int, the sample rate

        Returns:
            Tuple[numpy.ndarray. numpy.ndarray], the kernel and the phase response.
            The kernel is a numpy array containing the filter kernel. The kernel
            can be used, for example, with gstreamer's stock audiofirfilter element.

        """
        #
        # compute abs of FFT of kernel
        #

        # FIXME check for change in length
        if self.fwdplan is None:
            self.fwdplan = lal.CreateForwardCOMPLEX16FFTPlan(len(linear_phase_kernel), 1)
        if self.revplan is None:
            self.revplan = lal.CreateReverseCOMPLEX16FFTPlan(len(linear_phase_kernel), 1)

        deltaT = 1. / sample_rate
        deltaF = 1. / (len(linear_phase_kernel) * deltaT)
        working_length = len(linear_phase_kernel)

        kernel_tseries = lal.CreateCOMPLEX16TimeSeries(
            name = "timeseries of whitening kernel",
            epoch = LIGOTimeGPS(0.),
            f0 = 0.,
            deltaT = deltaT,
            length = working_length,
            sampleUnits = lal.Unit("strain")
        )
        
        kernel_tseries.data.data = linear_phase_kernel

        absX = lal.CreateCOMPLEX16FrequencySeries(
            name = "absX",
            epoch = LIGOTimeGPS(0),
            f0 = 0.0,
            deltaF = deltaF,
            length = working_length,
            sampleUnits = lal.Unit("strain s")
        )

        logabsX = lal.CreateCOMPLEX16FrequencySeries(
            name = "absX",
            epoch = LIGOTimeGPS(0),
            f0 = 0.0,
            deltaF = deltaF,
            length = working_length,
            sampleUnits = lal.Unit("strain s")
        )

        cepstrum = lal.CreateCOMPLEX16TimeSeries(
            name = "cepstrum",
            epoch = LIGOTimeGPS(0.),
            f0 = 0.,
            deltaT = deltaT,
            length = working_length,
            sampleUnits = lal.Unit("strain")
        )

        theta = lal.CreateCOMPLEX16FrequencySeries(
            name = "theta",
            epoch = LIGOTimeGPS(0),
            f0 = 0.0,
            deltaF = deltaF,
            length = working_length,
            sampleUnits = lal.Unit("strain s")
        )

        min_phase_kernel = lal.CreateCOMPLEX16TimeSeries(
            name = "min phase kernel",
            epoch = LIGOTimeGPS(0.),
            f0 = 0.,
            deltaT = deltaT,
            length = working_length,
            sampleUnits = lal.Unit("strain")
        )

        lal.COMPLEX16TimeFreqFFT(absX, kernel_tseries, self.fwdplan)
        absX.data.data[:] = abs(absX.data.data)

        #
        # compute the cepstrum of the kernel (i.e., the iFFT of the
        # log of the abs of the FFT of the kernel)
        #

        logabsX.data.data[:] = numpy.log(absX.data.data)
        lal.COMPLEX16FreqTimeFFT(cepstrum, logabsX, self.revplan)

        #
        # multiply cepstrum by sgn
        #

        cepstrum.data.data[0] = 0.
        cepstrum.data.data[working_length // 2] = 0.
        cepstrum.data.data[working_length // 2 + 1:] = -cepstrum.data.data[working_length // 2 + 1:]

        #
        # compute theta
        #

        lal.COMPLEX16TimeFreqFFT(theta, cepstrum, self.fwdplan)

        #
        # compute the gain and phase of the zero-phase
        # approximation relative to the original linear-phase
        # filter
        #

        theta_data = theta.data.data[working_length // 2:]
        #gain = numpy.exp(theta_data.real)
        phase = -theta_data.imag

        #
        # apply optional masked phase adjustment
        #

        if self.target_phase is not None:
            # compute phase adjustment for +ve frequencies
            phase_adjustment = (self.target_phase - phase) * self.target_phase_mask

            # combine with phase adjustment for -ve frequencies
            phase_adjustment = numpy.concatenate((phase_adjustment[1:][-1::-1].conj(), phase_adjustment))

            # apply adjustment.  phase adjustment is what we
            # wish to add to the phase.  theta's imaginary
            # component contains the negative of the phase, so
            # we need to add -phase to theta's imaginary
            # component
            theta.data.data += -1.j * phase_adjustment

            # report adjusted phase
            #phase = -theta.data.data[working_length // 2:].imag

        #
        # compute minimum phase kernel
        #

        absX.data.data *= numpy.exp(theta.data.data)
        lal.COMPLEX16FreqTimeFFT(min_phase_kernel, absX, self.revplan)

        kernel = min_phase_kernel.data.data.real

        #
        # this kernel needs to be reversed to follow conventions
        # used with the audiofirfilter and lal_firbank elements
        #

        kernel = kernel[-1::-1]

        #
        # done
        #

        return kernel, phase

def generate_early_warning_psds(psd_path):
    """
    Definitely make this less hardcoded!!
    """
    tlen = 3140000
    sample_rate = 0.2
    length = int(tlen * sample_rate)
    flen = length // 2 + 1
    LISA_A_PSD = pycbc.psd.from_txt(os.path.join(psd_path, 'A_psd_4_smooth.txt'), flen, 1./tlen, 1./tlen, is_asd_file=False)
    psd_kern = PSDFirKernel()

    lisa_a_psd_lal = LISA_A_PSD.lal()

    first_psd_kern, latency, sample_rate = psd_kern.psd_to_linear_phase_whitening_fir_kernel(lisa_a_psd_lal)
    lisa_a_zero_phase_kern, phase = psd_kern.linear_phase_fir_kernel_to_minimum_phase_whitening_fir_kernel(first_psd_kern, sample_rate)
    lisa_a_zero_phase_kern = lisa_a_zero_phase_kern * sample_rate**(1.5) / 2**0.5
    kernel_length = 10000
    lisa_a_zero_phase_kern_cut = lisa_a_zero_phase_kern[-kernel_length:]

    lisa_a_zero_phase_kern_pycbc = pycbc.types.TimeSeries(pycbc.types.zeros(2592000//5), delta_t=5.)
    lisa_a_zero_phase_kern_pycbc.data[-kernel_length:] = lisa_a_zero_phase_kern[-kernel_length:]
    lisa_a_zero_phase_kern_pycbc.data[0] = lisa_a_zero_phase_kern[0] # Do I want to do this??

    lisa_a_zero_phase_kern_pycbc_fd = pycbc.types.FrequencySeries(pycbc.types.zeros(len(lisa_a_zero_phase_kern_pycbc) //2 + 1, dtype=np.complex128), delta_f = 1/2592000)
    pycbc.fft.fft(lisa_a_zero_phase_kern_pycbc, lisa_a_zero_phase_kern_pycbc_fd)
    lisa_a_zero_phase_kern_pycbc_td = lisa_a_zero_phase_kern_pycbc

    LISA_E_PSD = pycbc.psd.from_txt(os.path.join(psd_path, 'E_psd_4_smooth.txt'), flen, 1./tlen, 1./tlen, is_asd_file=False)
    psd_kern = PSDFirKernel()

    lisa_e_psd_lal = LISA_E_PSD.lal()

    first_psd_kern, latency, sample_rate = psd_kern.psd_to_linear_phase_whitening_fir_kernel(lisa_e_psd_lal)
    lisa_e_zero_phase_kern, phase = psd_kern.linear_phase_fir_kernel_to_minimum_phase_whitening_fir_kernel(first_psd_kern, sample_rate)
    lisa_e_zero_phase_kern = lisa_e_zero_phase_kern * sample_rate**(1.5) / 2**0.5
    kernel_length = 10000
    lisa_e_zero_phase_kern_cut = lisa_e_zero_phase_kern[-kernel_length:]

    lisa_e_zero_phase_kern_pycbc = pycbc.types.TimeSeries(pycbc.types.zeros(2592000//5), delta_t=5.)
    lisa_e_zero_phase_kern_pycbc.data[-kernel_length:] = lisa_e_zero_phase_kern[-kernel_length:]
    lisa_e_zero_phase_kern_pycbc.data[0] = lisa_e_zero_phase_kern[0] # Do I want to do this??

    lisa_e_zero_phase_kern_pycbc_fd = pycbc.types.FrequencySeries(pycbc.types.zeros(len(lisa_e_zero_phase_kern_pycbc) //2 + 1, dtype=np.complex128), delta_f = 1/2592000)
    pycbc.fft.fft(lisa_e_zero_phase_kern_pycbc, lisa_e_zero_phase_kern_pycbc_fd)
    lisa_e_zero_phase_kern_pycbc_td = lisa_e_zero_phase_kern_pycbc
    return [lisa_a_zero_phase_kern_pycbc_fd, lisa_a_zero_phase_kern_pycbc_td], [lisa_e_zero_phase_kern_pycbc_fd, lisa_e_zero_phase_kern_pycbc_td]


def generate_data_lisa_ew(
    waveform_params,
    psds_for_datagen,
    psds_for_whitening,
    seed,
    window_length,
    cutoff_time,
    extra_forward_zeroes=0
):
    window = signal.windows.hann(window_length * 2 + 1)[:window_length]

    outs = pycbc.waveform.get_fd_det_waveform(ifos=['LISA_A','LISA_E','LISA_T'], **waveform_params)
    tout_A = outs['LISA_A'].to_timeseries()
    tout_E = outs['LISA_E'].to_timeseries()
    strain_A = pycbc.noise.noise_from_psd(len(tout_A), tout_A.delta_t, psds_for_datagen['LISA_A'], seed=seed+137) * 0.
    strain_E = pycbc.noise.noise_from_psd(len(tout_E), tout_E.delta_t, psds_for_datagen['LISA_E'], seed=seed+13812476) * 0.
    strain_A.start_time += tout_A.start_time
    strain_E.start_time += tout_E.start_time
    strain_A[:] += tout_A[:]
    strain_E[:] += tout_E[:]
    strain_w_A = strain_A[:].copy()
    if extra_forward_zeroes:
        strain_w_A.data[:int(extra_forward_zeroes//5)] = 0
    strain_w_A.data[int(extra_forward_zeroes//5):int(extra_forward_zeroes//5)+window_length] *= window
    strain_w_A.data[-int(cutoff_time//5):] = 0
    strain_w_E = strain_E[:].copy()
    if extra_forward_zeroes:
        strain_w_E.data[:int(extra_forward_zeroes//5)] = 0
    strain_w_E.data[int(extra_forward_zeroes//5):int(extra_forward_zeroes//5)+window_length] *= window
    strain_w_E.data[-int(cutoff_time//5):] = 0

    strain_fout_A = pycbc.strain.strain.execute_cached_fft(
        strain_w_A,
        copy_output=False,
        uid=1236
    )
    strain_fout_A = strain_fout_A * (psds_for_whitening['LISA_A']).conj()
    strain_ww_A = pycbc.strain.strain.execute_cached_ifft(
        strain_fout_A,
        copy_output=False,
        uid=1237
    )
    # FIXME: Might need this!
    strain_ww_A.data[-int(cutoff_time//5):] = 0

    strain_fout_E = pycbc.strain.strain.execute_cached_fft(
        strain_w_E,
        copy_output=False,
        uid=1238
    )
    strain_fout_E = strain_fout_E * (psds_for_whitening['LISA_E']).conj()
    strain_ww_E = pycbc.strain.strain.execute_cached_ifft(
        strain_fout_E,
        copy_output=False,
        uid=1239
    )
    # FIXME: Might need this!
    strain_ww_E.data[-int(cutoff_time//5):] = 0

    return strain_ww_A, strain_ww_E


_WINDOW = None
def generate_waveform_lisa_ew(
    waveform_params,
    psds_for_whitening,
    window_length,
    cutoff_time,
    kernel_length,
    extra_forward_zeroes=0,
):
    """
    
    Parameters
    ----------
    waveform_params: dict
        A dictionary of waveform parameters that will be passed to the waveform
        generator.
    psds_for_whitening: dict[str: FrequencySeries]
        Power spectral denisities for whitening in the frequency-domain.
    window_length : int
        Length (in samples) of time-domain window applied to the start of the
        waveform.
    cutoff_time: float
        Time (in seconds) from the end of the waveform to cutoff.
    kernel_length : int
        Unused.
    extra_foward_zeroes : int
        Time (in seconds) to set to zero at the start of the waveform. If used,
        the window will be applied starting after the zeroes.
    """
    global _WINDOW
    if _WINDOW is None or len(_WINDOW) != window_length:
        _WINDOW = signal.windows.hann(window_length * 2 + 1)[:window_length]
    window = _WINDOW
    
    outs = pycbc.waveform.get_fd_det_waveform(ifos=['LISA_A','LISA_E','LISA_T'], **waveform_params)
    tout_A = pycbc.strain.strain.execute_cached_ifft(
        outs['LISA_A'],
        copy_output=False,
        uid=1234
    )
    tout_E = pycbc.strain.strain.execute_cached_ifft(
        outs['LISA_E'],
        copy_output=False,
        uid=1233
    )
    
    if extra_forward_zeroes:
        tout_A.data[:int(extra_forward_zeroes//5)] = 0
    tout_A.data[int(extra_forward_zeroes//5):int(extra_forward_zeroes//5)+window_length] *= window
    tout_A.data[-int(cutoff_time//5):] = 0
    
    if extra_forward_zeroes:
        tout_E.data[:int(extra_forward_zeroes//5)] = 0
    tout_E.data[int(extra_forward_zeroes//5):int(extra_forward_zeroes//5)+window_length] *= window
    tout_E.data[-int(cutoff_time//5):] = 0

    #print(len(tout))
    #tout.roll(cutoff_time//5)
    #plt.plot(tout)
    fout_A = pycbc.strain.strain.execute_cached_fft(
        tout_A,
        copy_output=True,
        uid=1235
    )
    #print(1./fout.delta_f, 1./lisa_a_zero_phase_kern_pycbc_fd.delta_f)
    fout_A.data[:] = fout_A.data[:] * (psds_for_whitening['LISA_A'].data[:]).conj()
    
    fout_E = pycbc.strain.strain.execute_cached_fft(
        tout_E,
        copy_output=True,
        uid=12350
    )
    #print(1./fout.delta_f, 1./lisa_a_zero_phase_kern_pycbc_fd.delta_f)
    fout_E.data[:] = fout_E.data[:] * (psds_for_whitening['LISA_E'].data[:]).conj()

    fout_ww_A = pycbc.strain.strain.execute_cached_ifft(
        fout_A,
        copy_output=False,
        uid=5237
    )
    # FIXME: Might need this!
    fout_ww_A.data[-int(cutoff_time//5):] = 0

    fout_A = pycbc.strain.strain.execute_cached_fft(
        fout_ww_A,
        copy_output=False,
        uid=5238
    )

    fout_ww_E = pycbc.strain.strain.execute_cached_ifft(
        fout_E,
        copy_output=False,
        uid=5247
    )
    # FIXME: Might need this!
    fout_ww_E.data[-int(cutoff_time//5):] = 0

    fout_E = pycbc.strain.strain.execute_cached_fft(
        fout_ww_E,
        copy_output=False,
        uid=5248
    )

    fouts = {}
    fouts['LISA_A'] = fout_A
    fouts['LISA_E'] = fout_E

    return fouts
