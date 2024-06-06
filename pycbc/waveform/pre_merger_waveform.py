from scipy import signal

import pycbc.fft
import pycbc.noise
import pycbc.strain
import pycbc.waveform


def apply_pre_merger_kernel(
    tout,
    whitening_psd,
    window,
    window_length,
    nefz,
    nctf,
    uids,
):
    """Helper function to apply the pre-merger kernel.
    
    Parameters
    ----------
    tout : pycbc.types.TimeSeries
        Time series to apply the kernel to.
    whitening_psd : pycbc.types.FrequencySeries
        PSD for whitening the data in the frequency-domain.
    window : numpy.ndarray
        Window array.
    window_length : int
        Pre-computed length of the window in samples.
    nefz : int
        Number of extra forward zeroes.
    nctf : int
        Number of samples to zero at the end of the data.
    uids : tuple
        Unique UIDs for computing the (i)FFTs. Must be length 2.

    Returns
    -------
    pycbc.types.TimeSeries
        Whitened time series.
    """
    # Zero initial data
    tout.data[:nefz] = 0
    # Apply window
    tout.data[nefz:nefz+window_length] *= window
    # Zero data from cutoff
    tout.data[-nctf:] = 0

    # TD to FD for whitening 
    fout = pycbc.strain.strain.execute_cached_fft(
        tout,
        copy_output=False,
        uid=uids[0],
    )
    # Whiten data
    fout.data[:] = fout.data[:] * (whitening_psd.data[:]).conj()

    # TD to FD to reapply zeroes
    tout_ww = pycbc.strain.strain.execute_cached_ifft(
        fout,
        copy_output=False,
        uid=uids[1],
    )
    # Zero initial data
    tout_ww.data[:nefz] = 0
    # Zero from cutoff
    tout_ww.data[-nctf:] = 0
    return tout_ww


def generate_data_lisa_pre_merger(
    waveform_params,
    psds_for_datagen,
    sample_rate,
    seed=137,
    zero_noise=False,
    no_signal=False,
    duration=None,
):
    """Generate pre-merger LISA data.

    UIDs used for FFTs: 4235(0), 4236(0)
    
    Parameters
    ----------
    waveform_params : dict
        Dictionary of waveform parameters
    psds_for_datagen : dict
        PSDs for data generation.
    sample_rate : float
        Sampling rate in Hz. 
    extra_forward_zeros : float
        Time (in seconds) to set to zero at the start of the waveform. If used,
        the window will be applied starting after the zeroes.
    seed : int
        Random seed used for generating the noise.
    zero_noise : bool
        If true, the noise will be set to zero.
    no_signal : bool
        If true, the signal will not be added to data and only noise will
        be returned.
    duration : float, optional
        If specified, the waveform will be truncated to match the specified
        duration.

    Returns
    -------
    Dict[str: pycbc.types.TimeSeries]
        Dictionary containing the time-domain data for each channel.
    """
    # Generate injection
    outs = pycbc.waveform.get_fd_det_waveform(
        ifos=['LISA_A','LISA_E','LISA_T'],
        **waveform_params
    )

    # Shift waveform so the merger is not at the end of the data
    outs['LISA_A'] = outs['LISA_A'].cyclic_time_shift(-waveform_params['additional_end_data'])
    outs['LISA_E'] = outs['LISA_E'].cyclic_time_shift(-waveform_params['additional_end_data'])

    # FS waveform to TD
    tout_A = outs['LISA_A'].to_timeseries()
    tout_E = outs['LISA_E'].to_timeseries()

    # Generate TD noise from the original PSDs
    strain_w_A = pycbc.noise.noise_from_psd(
        len(tout_A),
        tout_A.delta_t,
        psds_for_datagen['LISA_A'],
        seed=seed,
    )
    strain_w_E = pycbc.noise.noise_from_psd(
        len(tout_E),
        tout_E.delta_t,
        psds_for_datagen['LISA_E'],
        seed=seed + 1,
    )

    # We need to make sure the noise times match the signal
    strain_w_A._epoch = tout_A._epoch
    strain_w_E._epoch = tout_E._epoch

    # If zero noise, set noise to zero
    if zero_noise:
        strain_w_A *= 0.0
        strain_w_E *= 0.0

    # Only add signal if no_signal=False
    if not no_signal:
        strain_w_A[:] += tout_A[:]
        strain_w_E[:] += tout_E[:]

    # If duration is specified, discard the extra data
    if duration is not None:
        if duration > tout_A.duration:
            raise RuntimeError(
                "Specified duration is longer than the generated waveform"
            )
        nkeep = int(duration * sample_rate)
        # New start time will be nkeep sample time
        new_epoch = strain_w_A.sample_times[-nkeep]
        strain_w_A = pycbc.types.TimeSeries(
            strain_w_A.data[-nkeep:],
            delta_t=strain_w_A.delta_t,
        )
        strain_w_E = pycbc.types.TimeSeries(
            strain_w_E.data[-nkeep:],
            delta_t=strain_w_E.delta_t,
        )
        # Set the start time so that the GPS time is still correct
        strain_w_A.start_time = new_epoch
        strain_w_E.start_time = new_epoch
    
    return {
        "LISA_A": strain_w_A,
        "LISA_E": strain_w_E,
    }


def pre_process_data_lisa_pre_merger(
    data,
    sample_rate,
    psds_for_whitening,
    window_length,
    cutoff_time,
    extra_forward_zeroes=0,
):
    """Pre-process the pre-merger data.

    The data is truncated, windowed and whitened.

    data : dict
        Dictionary containing time-domain data. 
    sample_rate : float
        Sampling rate in Hz. 
    psds_for_whitening : dict
        PSDs for whitening.
    window_length : int
        Length of the hann window use to taper the start of the data.
    cutoff_time : float
        Time (in seconds) from the end of the waveform to cutoff.
    extra_forward_zeros : float
        Time (in seconds) to set to zero at the start of the waveform. If used,
        the window will be applied starting after the zeroes.

    Returns
    -------
    Dict[str: pycbc.types.TimeSeries]
        Dictionary containing the time-domain data for each channel.
    """
    # Compute the hann window 
    window = signal.windows.hann(window_length * 2 + 1)[:window_length]

    # Number of samples to zero
    nefz = int(extra_forward_zeroes * sample_rate)
    nctf = int(cutoff_time * sample_rate)

    # Apply pre-merger kernel to both channels
    strain_ww = {}
    strain_ww["LISA_A"] = apply_pre_merger_kernel(
        data["LISA_A"],
        whitening_psd=psds_for_whitening["LISA_A"],
        window=window,
        window_length=window_length,
        nefz=nefz,
        nctf=nctf,
        uids=(4235, 4236),
    )
    strain_ww["LISA_E"] = apply_pre_merger_kernel(
        data["LISA_E"],
        whitening_psd=psds_for_whitening["LISA_E"],
        window=window,
        window_length=window_length,
        nefz=nefz,
        nctf=nctf,
        uids=(42350, 42360),
    )
    return strain_ww



_WINDOW = None
def generate_waveform_lisa_pre_merger(
    waveform_params,
    psds_for_whitening,
    sample_rate,
    window_length,
    cutoff_time,
    extra_forward_zeroes=0,
):
    """Generate a pre-merger LISA waveform.

    UIDs used for FFTs: 1234(0), 1235(0), 1236(0), 1237(0) 
    
    Parameters
    ----------
    waveform_params: dict
        A dictionary of waveform parameters that will be passed to the waveform
        generator.
    psds_for_whitening: dict[str: FrequencySeries]
        Power spectral denisities for whitening in the frequency-domain.
    sample_rate : float
        Sampling rate.
    window_length : int
        Length (in samples) of time-domain window applied to the start of the
        waveform.
    cutoff_time: float
        Time (in seconds) from the end of the waveform to cutoff.
    kernel_length : int
        Unused.
    extra_forward_zeroes : int
        Time (in seconds) to set to zero at the start of the waveform. If used,
        the window will be applied starting after the zeroes.
    """
    global _WINDOW
    if _WINDOW is None or len(_WINDOW) != window_length:
        _WINDOW = signal.windows.hann(window_length * 2 + 1)[:window_length]
    window = _WINDOW
    
    nefz = int(extra_forward_zeroes * sample_rate)
    nctf = int(cutoff_time * sample_rate)

    # FIXME: do we need to generate LISA_T
    outs = pycbc.waveform.get_fd_det_waveform(
        ifos=['LISA_A','LISA_E','LISA_T'], **waveform_params
    )

    # FD waveform to TD so we can apply pre-merger kernel
    tout_A = pycbc.strain.strain.execute_cached_ifft(
        outs["LISA_A"],
        copy_output=False,
        uid=1234,
    )
    tout_E = pycbc.strain.strain.execute_cached_ifft(
        outs["LISA_E"],
        copy_output=False,
        uid=12340,
    )

    # Apply pre-merger kernel
    tout_A_ww = apply_pre_merger_kernel(
        tout_A,
        whitening_psd=psds_for_whitening["LISA_A"],
        window=window,
        window_length=window_length,
        nefz=nefz,
        nctf=nctf,
        uids=(1235, 1236),
    )
    tout_E_ww = apply_pre_merger_kernel(
        tout_E,
        whitening_psd=psds_for_whitening["LISA_E"],
        window=window,
        window_length=window_length,
        nefz=nefz,
        nctf=nctf,
        uids=(12350, 12360),
    )
    
    # Back to FD for search/inference
    fouts_ww = {}
    fouts_ww["LISA_A"] = pycbc.strain.strain.execute_cached_fft(
        tout_A_ww,
        copy_output=False,
        uid=1237,
    )
    fouts_ww["LISA_E"] = pycbc.strain.strain.execute_cached_fft(
        tout_E_ww,
        copy_output=False,
        uid=12370,
    )
    return fouts_ww
