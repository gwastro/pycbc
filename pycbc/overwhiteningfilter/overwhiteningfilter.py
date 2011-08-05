class Spectrum():
    def __init__(self, dt, n):
        """
        Create memory for the power spectrum
        """
        self.__spectrum = new_real_vector(N)
        self.__spectrum_dt = dt
        self.__spectrum_n = n

    def power_spectrum(self, data, psd_type, window_type, 
            block_length, block_overlap):
        """
        psd_type: estimation method (e.g. mean, median, medianmean)
        window_type: time-domain window for fft segments (e.g. hann, hanning)
        block_length: length of time domain segments used to compute average psd
        block_overlap: overlap between time domain segments

        If the desired output n and dt (specified when the class is initialized)
        do not match the resolution created by length, overlap, data.dt and data.n
        then this function will interpolate the computed psd to the desired
        resolution.
        """
        pass

    def inverse(self):
        """
        Compute the inverse of the data in self.__spectrum
        """
        pass

    def truncate(self, length_in_seconds):
        """
        Truncate the PSD to a specified length by converting to the time domain
        and then back to the frequency domain
        """
        pass
