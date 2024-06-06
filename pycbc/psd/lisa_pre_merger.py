import lal
from lal import LIGOTimeGPS
import math
import numpy as np
from typing import Optional, Tuple

import pycbc.fft
import pycbc.psd
import pycbc.types


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
        raise NotImplementedError(
            "`PSDFirKernel.set_phase` is not implemented!"
        )
        kernel, latency, sample_rate = self.psd_to_linear_phase_whitening_fir_kernel(psd)
        kernel, phase = self.linear_phase_fir_kernel_to_minimum_phase_whitening_fir_kernel(kernel, sample_rate)

        # get merger model for SNR = 1.
        f_psd = psd.f0 + np.arange(len(psd.data.data)) * psd.deltaF
        horizon_distance = HorizonDistance(f_low, f_psd[-1], psd.deltaF, m1, m2)
        f_model, model= horizon_distance(psd, 1.)[1]

        # find the range of frequency bins covered by the merger
        # model
        kmin, kmax = f_psd.searchsorted(f_model[0]), f_psd.searchsorted(f_model[-1]) + 1

        # compute SNR=1 model's (d SNR^2 / df) spectral density
        unit_snr2_density = np.zeros_like(phase)
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
        invert: bool = True,
        nyquist: Optional[float] = None
    ) -> Tuple[np.ndarray, int, int]:
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

        Parameters
        ----------
        psd : lal.REAL8FrequencySeries
            The reference PSD.
        invert : bool
            Whether to invert the kernel. Defaults to True.
        nyquist :  float, optional
            Whether to change the Nyquist frequency. Disabled by default.

        Returns
        -------
        numpy.ndarray
            Array containing the filter kernel.
        int
            Filter latency in samples
        int
            The sample rate in Hz.
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
        tmp = np.zeros((2 * len(data) - 1,), dtype = data.dtype)
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

        kernel_fseries.data.data = np.sqrt(data) + 0.j
        lal.COMPLEX16FreqTimeFFT(kernel_tseries, kernel_fseries, self.revplan)
        kernel = kernel_tseries.data.data.real
        kernel = np.roll(kernel, (len(data) - 1) // 2) / sample_rate * 2

        #
        # apply a Tukey window whose flat bit is 50% of the kernel.
        # preserve the FIR kernel's square magnitude
        #

        norm_before = np.dot(kernel, kernel)
        kernel *= lal.CreateTukeyREAL8Window(len(data), .5).data.data
        kernel *= math.sqrt(norm_before / np.dot(kernel, kernel))

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
        linear_phase_kernel: np.ndarray,
        sample_rate: int
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute the minimum-phase response filter (zero latency)
        associated with a linear-phase response filter (latency
        equal to half the filter length).

        From "Design of Optimal Minimum-Phase Digital FIR Filters
        Using Discrete Hilbert Transforms", IEEE Trans. Signal
        Processing, vol. 48, pp. 1491-1495, May 2000.

        Parameters
        -----------
        linear_phase_kernel : numpy.ndarray
            The kernel to compute the minimum-phase kernel with.
        sample_rate : int
            The sample rate in Hz.

        Returns
        -------
        numpy.ndarray
            The filter kernel.
        numpy.ndarray
            The phase response.
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

        logabsX.data.data[:] = np.log(absX.data.data)
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
        #gain = np.exp(theta_data.real)
        phase = -theta_data.imag

        #
        # apply optional masked phase adjustment
        #

        if self.target_phase is not None:
            # compute phase adjustment for +ve frequencies
            phase_adjustment = (self.target_phase - phase) * self.target_phase_mask

            # combine with phase adjustment for -ve frequencies
            phase_adjustment = np.concatenate((phase_adjustment[1:][-1::-1].conj(), phase_adjustment))

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

        absX.data.data *= np.exp(theta.data.data)
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
    

def generate_pre_merger_psds(
    psd_file,
    duration,
    sample_rate,
    kernel_length=17280,
):
    """Generate the time- and frequency-domain pre-merger PSDs
    
    Parameters
    ----------
    psd_file : str
        Path to the PSD file.
    sample_rate : float
        The sample rate.
    duration : float
        Duration in seconds.
    kernel_length : int
        Length of the whitening kernel in samples.

    Returns
    -------
    dict
        A dictionary contain the PSDs. The keys denote frequency-domain (FD)
        and time-domain (TD).
    """
    flen = int(duration * sample_rate) // 2 + 1
    delta_f = 1 / duration
    delta_t = 1 / sample_rate
    td_psd_length = int(duration * sample_rate)

    psd = pycbc.psd.from_txt(
        psd_file, flen, delta_f, delta_f, is_asd_file=False
    )

    psd_kern = PSDFirKernel()

    pdf_lal = psd.lal()

    first_psd_kern, latency, sample_rate = \
        psd_kern.psd_to_linear_phase_whitening_fir_kernel(pdf_lal)
    zero_phase_kern, phase = \
        psd_kern.linear_phase_fir_kernel_to_minimum_phase_whitening_fir_kernel(first_psd_kern, sample_rate)
    zero_phase_kern = zero_phase_kern * sample_rate**(1.5) / 2**0.5

    # Time domain pycbc PSD
    zero_phase_kern_pycbc = pycbc.types.TimeSeries(
        pycbc.types.zeros(td_psd_length),
        delta_t=delta_t,
    )
    zero_phase_kern_pycbc.data[-kernel_length:] = zero_phase_kern[-kernel_length:]
    zero_phase_kern_pycbc.data[0] = zero_phase_kern[0]

    # Frequency domain pycbc PSD
    zero_phase_kern_pycbc_fd = pycbc.types.FrequencySeries(
        pycbc.types.zeros(
            len(zero_phase_kern_pycbc) //2 + 1,
            dtype=np.complex128,
        ),
        delta_f=delta_f
    )
    pycbc.fft.fft(zero_phase_kern_pycbc, zero_phase_kern_pycbc_fd)

    zero_phase_kern_pycbc_td = zero_phase_kern_pycbc
    return {
        "TD": zero_phase_kern_pycbc_td,
        "FD": zero_phase_kern_pycbc_fd,
    }
