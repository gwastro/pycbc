# Copyright (C) 2012  Tito Dal Canton, Josh Willis
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

"""
Provides a class representing a frequency series.
"""
import os as _os, h5py
from pycbc.types.array import Array, _convert, zeros, _noreal
import lal as _lal
import numpy as _numpy

class FrequencySeries(Array):
    """Models a frequency series consisting of uniformly sampled scalar values.

    Parameters
    ----------
    initial_array : array-like
        Array containing sampled data.
    delta_f : float
        Frequency between consecutive samples in Hertz.
    epoch : {None, lal.LIGOTimeGPS}, optional
        Start time of the associated time domain data in seconds.
    dtype : {None, data-type}, optional
        Sample data type.
    copy : boolean, optional
        If True, samples are copied to a new array.
    """

    def __init__(self, initial_array, delta_f=None, epoch="", dtype=None, copy=True):
        if len(initial_array) < 1:
            raise ValueError('initial_array must contain at least one sample.')
        if delta_f is None:
            try:
                delta_f = initial_array.delta_f
            except AttributeError:
                raise TypeError('must provide either an initial_array with a delta_f attribute, or a value for delta_f')
        if not delta_f > 0:
            raise ValueError('delta_f must be a positive number')
        # We gave a nonsensical default value to epoch so we can test if it's been set.
        # If the user passes in an initial_array that has an 'epoch' attribute and doesn't
        # pass in a value of epoch, then our new object's epoch comes from initial_array.
        # But if the user passed in a value---even 'None'---that will take precedence over
        # anything set in initial_array.  Finally, if the user passes in something without
        # an epoch attribute *and* doesn't pass in a value of epoch, it becomes 'None'
        if not isinstance(epoch,_lal.LIGOTimeGPS):
            if epoch == "":
                if isinstance(initial_array,FrequencySeries):
                    epoch = initial_array._epoch
                else:
                    epoch = _lal.LIGOTimeGPS(0)
            elif epoch is not None:
                try:
                    if isinstance(epoch, _numpy.generic):
                        # In python3 lal LIGOTimeGPS will not work on numpy
                        # types as input. A quick google on how to generically
                        # convert numpy floats/ints to the python equivalent
                        # https://stackoverflow.com/questions/9452775/
                        epoch = _lal.LIGOTimeGPS(epoch.item())
                    else:
                        epoch = _lal.LIGOTimeGPS(epoch)
                except:
                    raise TypeError('epoch must be either None or a lal.LIGOTimeGPS')
        Array.__init__(self, initial_array, dtype=dtype, copy=copy)
        self._delta_f = delta_f
        self._epoch = epoch

    def _return(self, ary):
        return FrequencySeries(ary, self._delta_f, epoch=self._epoch, copy=False)

    def _typecheck(self, other):
        if isinstance(other, FrequencySeries):
            try:
                _numpy.testing.assert_almost_equal(other._delta_f,
                                                   self._delta_f)
            except:
                raise ValueError('different delta_f')
            # consistency of _epoch is not required because we may want
            # to combine frequency series estimated at different times
            # (e.g. PSD estimation)

    def get_delta_f(self):
        """Return frequency between consecutive samples in Hertz.
        """
        return self._delta_f
    delta_f = property(get_delta_f,
                       doc="Frequency between consecutive samples in Hertz.")

    def get_epoch(self):
        """Return frequency series epoch as a LIGOTimeGPS.
        """
        return self._epoch
    epoch = property(get_epoch,
                     doc="Frequency series epoch as a LIGOTimeGPS.")

    def get_sample_frequencies(self):
        """Return an Array containing the sample frequencies.
        """
        return Array(range(len(self))) * self._delta_f
    sample_frequencies = property(get_sample_frequencies,
                                  doc="Array of the sample frequencies.")

    def _getslice(self, index):
        if index.step is not None:
            new_delta_f = self._delta_f * index.step
        else:
            new_delta_f = self._delta_f
        return FrequencySeries(Array._getslice(self, index),
                               delta_f=new_delta_f,
                               epoch=self._epoch,
                               copy=False)

    def at_frequency(self, freq):
        """ Return the value at the specified frequency
        """
        return self[int(freq / self.delta_f)]

    @property
    def start_time(self):
        """Return the start time of this vector
        """
        return self.epoch

    @start_time.setter
    def start_time(self, time):
        """ Set the start time
        """
        self._epoch = _lal.LIGOTimeGPS(time)

    @property
    def end_time(self):
        """Return the end time of this vector
        """
        return self.start_time + self.duration

    @property
    def duration(self):
        """Return the time duration of this vector
        """
        return 1.0 / self.delta_f

    @property
    def delta_t(self):
        """Return the time between samples if this were a time series.
        This assume the time series is even in length!
        """
        return 1.0 / self.sample_rate

    @property
    def sample_rate(self):
        """Return the sample rate this would have in the time domain. This
        assumes even length time series!
        """
        return (len(self) - 1) * self.delta_f * 2.0

    def __eq__(self,other):
        """
        This is the Python special method invoked whenever the '=='
        comparison is used.  It will return true if the data of two
        frequency series are identical, and all of the numeric meta-data
        are identical, irrespective of whether or not the two
        instances live in the same memory (for that comparison, the
        Python statement 'a is b' should be used instead).

        Thus, this method returns 'True' if the types of both 'self'
        and 'other' are identical, as well as their lengths, dtypes,
        epochs, delta_fs and the data in the arrays, element by element.
        It will always do the comparison on the CPU, but will *not* move
        either object to the CPU if it is not already there, nor change
        the scheme of either object. It is possible to compare a CPU
        object to a GPU object, and the comparison should be true if the
        data and meta-data of the two objects are the same.

        Note in particular that this function returns a single boolean,
        and not an array of booleans as Numpy does.  If the numpy
        behavior is instead desired it can be obtained using the numpy()
        method of the PyCBC type to get a numpy instance from each
        object, and invoking '==' on those two instances.

        Parameters
        ----------
        other: another Python object, that should be tested for equality
            with 'self'.

        Returns
        -------
        boolean: 'True' if the types, dtypes, lengths, epochs, delta_fs
            and data of the two objects are each identical.
        """
        if super(FrequencySeries,self).__eq__(other):
            return (self._epoch == other._epoch and self._delta_f == other._delta_f)
        else:
            return False

    def almost_equal_elem(self,other,tol,relative=True,dtol=0.0):
        """
        Compare whether two frequency series are almost equal, element
        by element.

        If the 'relative' parameter is 'True' (the default) then the
        'tol' parameter (which must be positive) is interpreted as a
        relative tolerance, and the comparison returns 'True' only if
        abs(self[i]-other[i]) <= tol*abs(self[i])
        for all elements of the series.

        If 'relative' is 'False', then 'tol' is an absolute tolerance,
        and the comparison is true only if
        abs(self[i]-other[i]) <= tol
        for all elements of the series.

        The method also checks that self.delta_f is within 'dtol' of
        other.delta_f; if 'dtol' has its default value of 0 then exact
        equality between the two is required.

        Other meta-data (type, dtype, length, and epoch) must be exactly
        equal.  If either object's memory lives on the GPU it will be
        copied to the CPU for the comparison, which may be slow. But the
        original object itself will not have its memory relocated nor
        scheme changed.

        Parameters
        ----------
        other: another Python object, that should be tested for
            almost-equality with 'self', element-by-element.
        tol: a non-negative number, the tolerance, which is interpreted
            as either a relative tolerance (the default) or an absolute
            tolerance.
        relative: A boolean, indicating whether 'tol' should be interpreted
            as a relative tolerance (if True, the default if this argument
            is omitted) or as an absolute tolerance (if tol is False).
        dtol: a non-negative number, the tolerance for delta_f. Like 'tol',
            it is interpreted as relative or absolute based on the value of
            'relative'.  This parameter defaults to zero, enforcing exact
            equality between the delta_f values of the two FrequencySeries.

        Returns
        -------
        boolean: 'True' if the data and delta_fs agree within the tolerance,
            as interpreted by the 'relative' keyword, and if the types,
            lengths, dtypes, and epochs are exactly the same.
        """
        # Check that the delta_f tolerance is non-negative; raise an exception
        # if needed.
        if (dtol < 0.0):
            raise ValueError("Tolerance in delta_f cannot be negative")
        if super(FrequencySeries,self).almost_equal_elem(other,tol=tol,relative=relative):
            if relative:
                return (self._epoch == other._epoch and
                        abs(self._delta_f-other._delta_f) <= dtol*self._delta_f)
            else:
                return (self._epoch == other._epoch and
                        abs(self._delta_f-other._delta_f) <= dtol)
        else:
            return False

    def almost_equal_norm(self,other,tol,relative=True,dtol=0.0):
        """
        Compare whether two frequency series are almost equal, normwise.

        If the 'relative' parameter is 'True' (the default) then the
        'tol' parameter (which must be positive) is interpreted as a
        relative tolerance, and the comparison returns 'True' only if
        abs(norm(self-other)) <= tol*abs(norm(self)).

        If 'relative' is 'False', then 'tol' is an absolute tolerance,
        and the comparison is true only if
        abs(norm(self-other)) <= tol

        The method also checks that self.delta_f is within 'dtol' of
        other.delta_f; if 'dtol' has its default value of 0 then exact
        equality between the two is required.

        Other meta-data (type, dtype, length, and epoch) must be exactly
        equal.  If either object's memory lives on the GPU it will be
        copied to the CPU for the comparison, which may be slow. But the
        original object itself will not have its memory relocated nor
        scheme changed.

        Parameters
        ----------
        other: another Python object, that should be tested for
            almost-equality with 'self', based on their norms.
        tol: a non-negative number, the tolerance, which is interpreted
            as either a relative tolerance (the default) or an absolute
            tolerance.
        relative: A boolean, indicating whether 'tol' should be interpreted
            as a relative tolerance (if True, the default if this argument
            is omitted) or as an absolute tolerance (if tol is False).
        dtol: a non-negative number, the tolerance for delta_f. Like 'tol',
            it is interpreted as relative or absolute based on the value of
            'relative'.  This parameter defaults to zero, enforcing exact
            equality between the delta_f values of the two FrequencySeries.

        Returns
        -------
        boolean: 'True' if the data and delta_fs agree within the tolerance,
            as interpreted by the 'relative' keyword, and if the types,
            lengths, dtypes, and epochs are exactly the same.
        """
        # Check that the delta_f tolerance is non-negative; raise an exception
        # if needed.
        if (dtol < 0.0):
            raise ValueError("Tolerance in delta_f cannot be negative")
        if super(FrequencySeries,self).almost_equal_norm(other,tol=tol,relative=relative):
            if relative:
                return (self._epoch == other._epoch and
                        abs(self._delta_f-other._delta_f) <= dtol*self._delta_f)
            else:
                return (self._epoch == other._epoch and
                        abs(self._delta_f-other._delta_f) <= dtol)
        else:
            return False

    @_convert
    def lal(self):
        """Produces a LAL frequency series object equivalent to self.

        Returns
        -------
        lal_data : {lal.*FrequencySeries}
            LAL frequency series object containing the same data as self.
            The actual type depends on the sample's dtype. If the epoch of
            self was 'None', the epoch of the returned LAL object will be
            LIGOTimeGPS(0,0); otherwise, the same as that of self.

        Raises
        ------
        TypeError
            If frequency series is stored in GPU memory.
        """

        lal_data = None
        if self._epoch is None:
            ep = _lal.LIGOTimeGPS(0,0)
        else:
            ep = self._epoch

        if self._data.dtype == _numpy.float32:
            lal_data = _lal.CreateREAL4FrequencySeries("",ep,0,self.delta_f,_lal.SecondUnit,len(self))
        elif self._data.dtype == _numpy.float64:
            lal_data = _lal.CreateREAL8FrequencySeries("",ep,0,self.delta_f,_lal.SecondUnit,len(self))
        elif self._data.dtype == _numpy.complex64:
            lal_data = _lal.CreateCOMPLEX8FrequencySeries("",ep,0,self.delta_f,_lal.SecondUnit,len(self))
        elif self._data.dtype == _numpy.complex128:
            lal_data = _lal.CreateCOMPLEX16FrequencySeries("",ep,0,self.delta_f,_lal.SecondUnit,len(self))

        lal_data.data.data[:] = self.numpy()

        return lal_data

    def save(self, path, group=None, ifo='P1'):
        """
        Save frequency series to a Numpy .npy, hdf, or text file. The first column
        contains the sample frequencies, the second contains the values.
        In the case of a complex frequency series saved as text, the imaginary
        part is written as a third column.  When using hdf format, the data is stored
        as a single vector, along with relevant attributes.

        Parameters
        ----------
        path: string
            Destination file path. Must end with either .hdf, .npy or .txt.

        group: string
            Additional name for internal storage use. Ex. hdf storage uses
            this as the key value.

        Raises
        ------
        ValueError
            If path does not end in .npy or .txt.
        """

        ext = _os.path.splitext(path)[1]
        if ext == '.npy':
            output = _numpy.vstack((self.sample_frequencies.numpy(),
                                    self.numpy())).T
            _numpy.save(path, output)
        elif ext == '.txt':
            if self.kind == 'real':
                output = _numpy.vstack((self.sample_frequencies.numpy(),
                                        self.numpy())).T
            elif self.kind == 'complex':
                output = _numpy.vstack((self.sample_frequencies.numpy(),
                                        self.numpy().real,
                                        self.numpy().imag)).T
            _numpy.savetxt(path, output)
        elif ext == '.xml' or path.endswith('.xml.gz'):
            from pycbc.io.ligolw import make_psd_xmldoc
            from ligo.lw import utils

            if self.kind != 'real':
                raise ValueError('XML only supports real frequency series')
            output = self.lal()
            # When writing in this format we must *not* have the 0 values at
            # frequencies less than flow. To resolve this we set the first
            # non-zero value < flow.
            data_lal = output.data.data
            first_idx = _numpy.argmax(data_lal>0)
            if not first_idx == 0:
                data_lal[:first_idx] = data_lal[first_idx]
            psddict = {ifo: output}
            utils.write_filename(
                make_psd_xmldoc(psddict),
                path,
                compress='auto'
            )
        elif ext == '.hdf':
            key = 'data' if group is None else group
            with h5py.File(path, 'a') as f:
                ds = f.create_dataset(key, data=self.numpy(),
                                      compression='gzip',
                                      compression_opts=9, shuffle=True)
                if self.epoch is not None:
                    ds.attrs['epoch'] = float(self.epoch)
                ds.attrs['delta_f'] = float(self.delta_f)
        else:
            raise ValueError('Path must end with .npy, .txt, .xml, .xml.gz '
                             'or .hdf')

    def to_frequencyseries(self):
        """ Return frequency series """
        return self

    @_noreal
    def to_timeseries(self, delta_t=None):
        """ Return the Fourier transform of this time series.

        Note that this assumes even length time series!

        Parameters
        ----------
        delta_t : {None, float}, optional
            The time resolution of the returned series. By default the
        resolution is determined by length and delta_f of this frequency
        series.

        Returns
        -------
        TimeSeries:
            The inverse fourier transform of this frequency series.
        """
        from pycbc.fft import ifft
        from pycbc.types import TimeSeries, real_same_precision_as
        nat_delta_t =  1.0 / ((len(self)-1)*2) / self.delta_f
        if not delta_t:
            delta_t = nat_delta_t

        # add 0.5 to round integer
        tlen  = int(1.0 / self.delta_f / delta_t + 0.5)
        flen = int(tlen / 2 + 1)

        if flen < len(self):
            raise ValueError("The value of delta_t (%s) would be "
                             "undersampled. Maximum delta_t "
                             "is %s." % (delta_t, nat_delta_t))
        if not delta_t:
            tmp = self
        else:
            tmp = FrequencySeries(zeros(flen, dtype=self.dtype),
                             delta_f=self.delta_f, epoch=self.epoch)
            tmp[:len(self)] = self[:]

        f = TimeSeries(zeros(tlen,
                           dtype=real_same_precision_as(self)),
                           delta_t=delta_t)
        ifft(tmp, f)
        f._delta_t = delta_t
        return f

    @_noreal
    def cyclic_time_shift(self, dt):
        """Shift the data and timestamps by a given number of seconds

        Shift the data and timestamps in the time domain a given number of
        seconds. To just change the time stamps, do ts.start_time += dt.
        The time shift may be smaller than the intrinsic sample rate of the data.
        Note that data will be cycliclly rotated, so if you shift by 2
        seconds, the final 2 seconds of your data will now be at the
        beginning of the data set.

        Parameters
        ----------
        dt : float
            Amount of time to shift the vector.

        Returns
        -------
        data : pycbc.types.FrequencySeries
            The time shifted frequency series.
        """
        from pycbc.waveform import apply_fseries_time_shift
        data = apply_fseries_time_shift(self, dt)
        data.start_time = self.start_time - dt
        return data

    def match(self, other, psd=None,
              low_frequency_cutoff=None, high_frequency_cutoff=None):
        """ Return the match between the two TimeSeries or FrequencySeries.

        Return the match between two waveforms. This is equivalent to the overlap
        maximized over time and phase. By default, the other vector will be
        resized to match self. Beware, this may remove high frequency content or the
        end of the vector.

        Parameters
        ----------
        other : TimeSeries or FrequencySeries
            The input vector containing a waveform.
        psd : Frequency Series
            A power spectral density to weight the overlap.
        low_frequency_cutoff : {None, float}, optional
            The frequency to begin the match.
        high_frequency_cutoff : {None, float}, optional
            The frequency to stop the match.
        index: int
            The number of samples to shift to get the match.

        Returns
        -------
        match: float
        index: int
            The number of samples to shift to get the match.
        """
        from pycbc.types import TimeSeries
        from pycbc.filter import match

        if isinstance(other, TimeSeries):
            if other.duration != self.duration:
                other = other.copy()
                other.resize(int(other.sample_rate * self.duration))

            other = other.to_frequencyseries()

        if len(other) != len(self):
            other = other.copy()
            other.resize(len(self))

        if psd is not None and len(psd) > len(self):
            psd = psd.copy()
            psd.resize(len(self))

        return match(self, other, psd=psd,
                     low_frequency_cutoff=low_frequency_cutoff,
                     high_frequency_cutoff=high_frequency_cutoff)

    def plot(self, **kwds):
        """ Basic plot of this frequency series
        """
        from matplotlib import pyplot

        if self.kind == 'real':
            plot = pyplot.plot(self.sample_frequencies, self, **kwds)
            return plot
        elif self.kind == 'complex':
            plot1 = pyplot.plot(self.sample_frequencies, self.real(), **kwds)
            plot2 = pyplot.plot(self.sample_frequencies, self.imag(), **kwds)
            return plot1, plot2

def load_frequencyseries(path, group=None):
    """Load a FrequencySeries from an HDF5, ASCII or Numpy file. The file type
    is inferred from the file extension, which must be `.hdf`, `.txt` or
    `.npy`.

    For ASCII and Numpy files, the first column of the array is assumed to
    contain the frequency. If the array has two columns, a real frequency
    series is returned. If the array has three columns, the second and third
    ones are assumed to contain the real and imaginary parts of a complex
    frequency series.

    For HDF files, the dataset is assumed to contain the attribute `delta_f`
    giving the frequency resolution in Hz. The attribute `epoch`, if present,
    is taken as the start GPS time (epoch) of the data in the series.

    The default data types will be double precision floating point.

    Parameters
    ----------
    path: string
        Input file path. Must end with either `.npy`, `.txt` or `.hdf`.

    group: string
        Additional name for internal storage use. When reading HDF files, this
        is the path to the HDF dataset to read.

    Raises
    ------
    ValueError
        If the path does not end in a supported extension.
        For Numpy and ASCII input files, this is also raised if the array
        does not have 2 or 3 dimensions.
    """
    ext = _os.path.splitext(path)[1]
    if ext == '.npy':
        data = _numpy.load(path)
    elif ext == '.txt':
        data = _numpy.loadtxt(path)
    elif ext == '.hdf':
        key = 'data' if group is None else group
        with h5py.File(path, 'r') as f:
            data = f[key][:]
            delta_f = f[key].attrs['delta_f']
            epoch = f[key].attrs['epoch'] if 'epoch' in f[key].attrs else None
            series = FrequencySeries(data, delta_f=delta_f, epoch=epoch)
        return series
    else:
        raise ValueError('Path must end with .npy, .hdf, or .txt')

    delta_f = (data[-1][0] - data[0][0]) / (len(data) - 1)
    if data.ndim == 2:
        return FrequencySeries(data[:,1], delta_f=delta_f, epoch=None)
    elif data.ndim == 3:
        return FrequencySeries(data[:,1] + 1j*data[:,2], delta_f=delta_f,
                               epoch=None)

    raise ValueError('File has %s dimensions, cannot convert to FrequencySeries, \
                      must be 2 (real) or 3 (complex)' % data.ndim)
