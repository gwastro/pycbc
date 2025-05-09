###################################################
Performing FFTs in PyCBC
###################################################

============
Introduction
============

Many applications in gravitational-wave analysis rely on Fast Fourier
Transforms. These are often the dominant computational cost in analyses. PyCBC
needs to balance the requirement that analyses be efficient with ease of use
for end users. To meet this requirement PyCBC provides two different APIs to
do FFTs:

* A function based API, which is easy to use, but not optimized.
* A class-based API, which is a little more involved to use, but allows the use of optimized FFT routines.

These APIs offer access to a number of FFT backends. PyCBC knows how to do FFTs
using the FFTW, MKL and numpy backends, and will enable these if they are
present on your system. By default FFTW will be used, then MKL if FFTW is not
and numpy will be used only if the other 2 are not present. However, you can
override this and choose a specific backend if multiple are available.

When running on GPUs, PyCBC knows how to do CUDA FFTs through the same
interface. 

============================
Using the function based API
============================

The PyCBC function based API offers a simple way to Fourier transform an
input array into an output array. This is done by::

    >>> from pycbc import fft
    >>> fft.fft(input_array, output_array)

Or for an inverse FFT::

    >>> from pycbc import fft
    >>> fft.ifft(input_array, output_array)

To do this requires having the output_array, which would be the Fourier
transform of the input, already existing as an array of 0s. This output array
must be the correct *length*, the correct *type* (complex or real) and the
correct *precision* (double or float precision). It's also worth noting that
fast Fourier transforms are more efficient if their lengths are 2**N where
N is some integer. This becomes a little complicated if doing real->complex
or complex->real transforms; in those cases the longer array should have a
length of 2**N to be efficient.

Here's a few examples::

    >>> import numpy as np
    >>> from pycbc import types
    >>> from pycbc import fft
    >>> inarr = types.Array(np.ones([64], dtype=np.complex64))
    >>> outarr = types.Array(np.zeros([64], dtype=np.complex64))
    >>> fft.fft(inarr, outarr)
    >>> print(outarr)

or (note here the length of the complex array is the length of the real array divided by 2 and then + 1)::

    >>> import numpy as np
    >>> from pycbc import types
    >>> from pycbc import fft
    >>> inarr = types.Array(np.ones([64], dtype=np.float32))
    >>> outarr = types.Array(np.zeros([33], dtype=np.complex64))
    >>> fft.fft(inarr, outarr)
    >>> print(outarr)

or (this one is an inverse FFT)::

    >>> import numpy as np
    >>> from pycbc import types
    >>> from pycbc import fft
    >>> inarr = types.Array(np.ones([33], dtype=np.complex64))
    >>> outarr = types.Array(np.zeros([64], dtype=np.float32))
    >>> fft.ifft(inarr, outarr)
    >>> print(outarr)

This will work the pycbc Timeseries and Frequencyseries as well. Except you
must FFT a TimeSeries to a FrequencySeries or IFFT a FrequencySeries to a
TimeSeries. In this case the time and frequency spacing must also be
consistent. For this reason we provide convenience functions that use the
function API, but figure out these details, and create the output array, for
you. As an example::

    >>> import numpy as np
    >>> from pycbc import types
    >>> inarr = types.TimeSeries(np.ones([64], dtype=np.float64), delta_t=1./64.)
    >>> outarr = inarr.to_frequencyseries()

or::

    >>> import numpy as np
    >>> from pycbc import types
    >>> inarr = types.FrequencySeries(np.ones([33], dtype=np.complex128), delta_f=1.)
    >>> outarr = inarr.to_timeseries()


=========================
Using the class-based API
=========================

The PyCBC class-based API should be used if you care about performance. If you
are performing FFTs many times, with inputs that are the same size each time,
using this will offer significance perfommance improvement. 

Here's how to use this::

    >>> from pycbc import fft
    >>> fft_class = pycbc.fft.FFT(inarr, outarr)
    >>> fft_class.execute()
    >>> outarr *= inarr._delta_t # ONLY IF inarr is a TimeSeries
    >>> outarr *= inarr._delta_f # ONLY IF inarr is a FrequencySeries

Or for an inverse FFT::

    >>> from pycbc import fft
    >>> ifft_class = pycbc.fft.IFFT(inarr, outarr)
    >>> ifft_class.execute()
    >>> outarr *= inarr._delta_t # ONLY IF inarr is a TimeSeries
    >>> outarr *= inarr._delta_f # ONLY IF inarr is a FrequencySeries

The idea would be that the `fft_class` or `ifft_class` would only be created
*once* and the execute command called many times. You would change the contents
of inarr before each call and outarr will update when execute is run. After
creating the FFT class *do not* reassign inarr, but instead set values. So::

    >>> fft_class = pycbc.fft.FFT(inarr, outarr)
    >>> inarr = types.TimeSeries(np.ones([64], dtype=np.float64), delta_t=1./64.)
    >>> fft_class.execute()

would not work! Instead do::

    >>> fft_class = pycbc.fft.FFT(inarr, outarr)
    >>> inarr[:] = np.ones([64], dtype=np.float64)[:]
    >>> fft_class.execute()


===========================
Choosing a specific backend
===========================

If you want to choose a specific backend, you can see what is available with::

    >>> from pycbc.fft import backend_support
    >>> backend_support.get_backend_names()

and then do::

    >>> from pycbc.fft import backend_support
    >>> backend_support.set_backend(['mkl'])

to set a specific backend. Running::

    >>> from pycbc.fft import backend_support
    >>> backend_support.get_backend()

will tell you what you are currently using. You can also use the
MKL `Scheme` to default to using MKL FFTs, instead of FFTW.

====================
Method documentation
====================

.. automodule:: pycbc.fft
    :noindex:
    :members:
