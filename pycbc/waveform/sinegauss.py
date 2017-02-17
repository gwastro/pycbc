""" Generation of sine-Gaussian bursty type things
"""

import pycbc.types
import numpy

def fd_sine_gaussian(amp, quality, central_frequency, fmin, fmax, delta_f):
    """ Generate a Fourier domain sine-Gaussian

    Parameters
    ----------
    amp: float
        Amplitude of the sine-Gaussian
    quality: float
        The quality factor
    central_frequency: float
        The central frequency of the sine-Gaussian
    fmin: float
        The minimum frequency to generate the sine-Gaussian. This determines
        the length of the output vector.
    fmax: float
        The maximum frequency to generate the sine-Gaussian
    delta_f: float
        The size of the frequency step
    
    Returns
    -------
    sg: pycbc.types.Frequencyseries
        A Fourier domain sine-Gaussian
    """
    f = numpy.arange(fmin, fmax, delta_f)
    kmax = int(fmax / delta_f)
    kmin = int(fmin / delta_f)
    tau = quality / 2 / numpy.pi / central_frequency
    A = amp * numpy.pi ** 0.5 / 2 * tau
    d = A * numpy.exp(-(numpy.pi  * tau  * (f - central_frequency))**2.0)
    d *= (1 + numpy.exp(-quality ** 2.0 * f / central_frequency))
    v = numpy.zeros(kmax, dtype=numpy.complex128)
    v[kmin:kmax] = d[:]
    return pycbc.types.FrequencySeries(v, delta_f=delta_f, copy=False)
