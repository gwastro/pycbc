# Copyright (C) 2014  Tito Dal Canton, Josh Willis, Alex Nitz
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
Provides a class representing a time series.
"""

import os as _os, h5py
from pycbc.types.array import Array, _convert, complex_same_precision_as, zeros
from pycbc.types.array import _nocomplex
from pycbc.types.frequencyseries import FrequencySeries
import lal as _lal
import numpy as _numpy
from scipy.io.wavfile import write as write_wav


class TimeSeries(Array):
    """Models a time series consisting of uniformly sampled scalar values.

    Parameters
    ----------
    initial_array : array-like
        Array containing sampled data.
    delta_t : float
        Time between consecutive samples in seconds.
    epoch : {None, lal.LIGOTimeGPS}, optional
        Time of the first sample in seconds.
    dtype : {None, data-type}, optional
        Sample data type.
    copy : boolean, optional
        If True, samples are copied to a new array.

    Attributes
    ----------
    delta_t
    duration
    start_time
    end_time
    sample_times
    sample_rate
    """

    def __init__(self, initial_array, delta_t=None, epoch="", dtype=None, copy=True):
        if len(initial_array) < 1:
            raise ValueError('initial_array must contain at least one sample.')
        if delta_t is None:
            try:
                delta_t = initial_array.delta_t
            except AttributeError:
                raise TypeError('must provide either an initial_array with a delta_t attribute, or a value for delta_t')
        if not delta_t > 0:
            raise ValueError('delta_t must be a positive number')
        # We gave a nonsensical default value to epoch so we can test if it's been set.
        # If the user passes in an initial_array that has an 'epoch' attribute and doesn't
        # pass in a value of epoch, then our new object's epoch comes from initial_array.
        # But if the user passed in a value---even 'None'---that will take precedence over
        # anything set in initial_array.  Finally, if the user passes in something without
        # an epoch attribute *and* doesn't pass in a value of epoch, it becomes 'None'
        if not isinstance(epoch,_lal.LIGOTimeGPS):
            if epoch == "":
                if isinstance(initial_array, TimeSeries):
                    epoch = initial_array._epoch
                else:
                    epoch = _lal.LIGOTimeGPS(0)
            elif epoch is not None:
                try: 
                    epoch = _lal.LIGOTimeGPS(epoch)
                except:
                    raise TypeError('epoch must be either None or a lal.LIGOTimeGPS')
        Array.__init__(self, initial_array, dtype=dtype, copy=copy)
        self._delta_t = delta_t
        self._epoch = epoch

    def _return(self, ary):
        return TimeSeries(ary, self._delta_t, epoch=self._epoch, copy=False)

    def _typecheck(self, other):
        if isinstance(other, TimeSeries):
            if other._delta_t != self._delta_t:
                raise ValueError('different delta_t')
            if self._epoch != other._epoch:
                raise ValueError('different epoch')

    def _getslice(self, index):
        # Set the new epoch---note that index.start may also be None
        if self._epoch is None:
            new_epoch = None
        elif index.start is None:
            new_epoch = self._epoch
        else:
            if index.start < 0:
                raise ValueError('Negative start index not supported')
            new_epoch = self._epoch + index.start * self._delta_t

        if index.step is not None:
            new_delta_t = self._delta_t * index.step
        else:
            new_delta_t = self._delta_t
        
        return TimeSeries(Array._getslice(self, index), new_delta_t,
                          new_epoch, copy=False)


    def prepend_zeros(self, num):
        """Prepend num zeros onto the beginning of this TimeSeries. Update also
        epoch to include this prepending.
        """
        self.resize(len(self) + num)
        self.roll(num)
        self._epoch = self._epoch - num * self._delta_t

    def append_zeros(self, num):
        """Append num zeros onto the end of this TimeSeries.
        """
        self.resize(len(self) + num)

    def get_delta_t(self):
        """Return time between consecutive samples in seconds.
        """
        return self._delta_t
    delta_t = property(get_delta_t,
                       doc="Time between consecutive samples in seconds.")

    def get_duration(self):
        """Return duration of time series in seconds.
        """
        return len(self) * self._delta_t
    duration = property(get_duration,
                        doc="Duration of time series in seconds.")
    
    def get_sample_rate(self):
        """Return the sample rate of the time series.
        """
        return int(1.0/self.delta_t)
    sample_rate = property(get_sample_rate,
                           doc="The sample rate of the time series.")
    @property
    def delta_f(self):
        """Return the delta_f this ts would have in the frequency domain
        """
        return 1.0 / self.duration

    @property
    def start_time(self):
        """Return time series start time as a LIGOTimeGPS.
        """
        return self._epoch

    @start_time.setter
    def start_time(self, time):
        """ Set the start time
        """
        self._epoch = _lal.LIGOTimeGPS(time)

    def get_end_time(self):
        """Return time series end time as a LIGOTimeGPS.
        """
        return self._epoch + self.get_duration()
    end_time = property(get_end_time,
                        doc="Time series end time as a LIGOTimeGPS.")

    def get_sample_times(self):
        """Return an Array containing the sample times.
        """
        if self._epoch is None:
            return Array(range(len(self))) * self._delta_t
        else:
            return Array(range(len(self))) * self._delta_t + float(self._epoch)
    sample_times = property(get_sample_times,
                            doc="Array containing the sample times.")

    def __eq__(self,other):
        """
        This is the Python special method invoked whenever the '=='
        comparison is used.  It will return true if the data of two
        time series are identical, and all of the numeric meta-data
        are identical, irrespective of whether or not the two
        instances live in the same memory (for that comparison, the
        Python statement 'a is b' should be used instead).

        Thus, this method returns 'True' if the types of both 'self'
        and 'other' are identical, as well as their lengths, dtypes,
        epochs, delta_ts and the data in the arrays, element by element.
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
        boolean: 'True' if the types, dtypes, lengths, epochs, delta_ts
            and data of the two objects are each identical.
        """
        if super(TimeSeries,self).__eq__(other):
            return (self._epoch == other._epoch and self._delta_t == other._delta_t)
        else:
            return False

    def almost_equal_elem(self,other,tol,relative=True,dtol=0.0):
        """
        Compare whether two time series are almost equal, element
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

        The method also checks that self.delta_t is within 'dtol' of
        other.delta_t; if 'dtol' has its default value of 0 then exact
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
        dtol: a non-negative number, the tolerance for delta_t. Like 'tol',
            it is interpreted as relative or absolute based on the value of
            'relative'.  This parameter defaults to zero, enforcing exact
            equality between the delta_t values of the two TimeSeries.

        Returns
        -------
        boolean: 'True' if the data and delta_ts agree within the tolerance,
            as interpreted by the 'relative' keyword, and if the types,
            lengths, dtypes, and epochs are exactly the same.
        """
        # Check that the delta_t tolerance is non-negative; raise an exception
        # if needed.
        if (dtol < 0.0):
            raise ValueError("Tolerance in delta_t cannot be negative")
        if super(TimeSeries,self).almost_equal_elem(other,tol=tol,relative=relative):
            if relative:
                return (self._epoch == other._epoch and
                        abs(self._delta_t-other._delta_t) <= dtol*self._delta_t)
            else:
                return (self._epoch == other._epoch and
                        abs(self._delta_t-other._delta_t) <= dtol)
        else:
            return False

    def almost_equal_norm(self,other,tol,relative=True,dtol=0.0):
        """
        Compare whether two time series are almost equal, normwise.

        If the 'relative' parameter is 'True' (the default) then the
        'tol' parameter (which must be positive) is interpreted as a
        relative tolerance, and the comparison returns 'True' only if
        abs(norm(self-other)) <= tol*abs(norm(self)).

        If 'relative' is 'False', then 'tol' is an absolute tolerance,
        and the comparison is true only if
        abs(norm(self-other)) <= tol

        The method also checks that self.delta_t is within 'dtol' of
        other.delta_t; if 'dtol' has its default value of 0 then exact
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
        dtol: a non-negative number, the tolerance for delta_t. Like 'tol',
            it is interpreted as relative or absolute based on the value of
            'relative'.  This parameter defaults to zero, enforcing exact
            equality between the delta_t values of the two TimeSeries.

        Returns
        -------
        boolean: 'True' if the data and delta_ts agree within the tolerance,
            as interpreted by the 'relative' keyword, and if the types,
            lengths, dtypes, and epochs are exactly the same.
        """
        # Check that the delta_t tolerance is non-negative; raise an exception
        # if needed.
        if (dtol < 0.0):
            raise ValueError("Tolerance in delta_t cannot be negative")
        if super(TimeSeries,self).almost_equal_norm(other,tol=tol,relative=relative):
            if relative:
                return (self._epoch == other._epoch and
                        abs(self._delta_t-other._delta_t) <= dtol*self._delta_t)
            else:
                return (self._epoch == other._epoch and
                        abs(self._delta_t-other._delta_t) <= dtol)
        else:
            return False

    @_convert
    def lal(self):
        """Produces a LAL time series object equivalent to self.

        Returns
        -------
        lal_data : {lal.*TimeSeries}
            LAL time series object containing the same data as self.
            The actual type depends on the sample's dtype.  If the epoch of
            self is 'None', the epoch of the returned LAL object will be
            LIGOTimeGPS(0,0); otherwise, the same as that of self.

        Raises
        ------
        TypeError
            If time series is stored in GPU memory.
        """
        lal_data = None
        if self._epoch is None:
            ep = _lal.LIGOTimeGPS(0,0)
        else:
            ep = self._epoch

        if self._data.dtype == _numpy.float32:
            lal_data = _lal.CreateREAL4TimeSeries("",ep,0,self.delta_t,_lal.SecondUnit,len(self))
        elif self._data.dtype == _numpy.float64:
            lal_data = _lal.CreateREAL8TimeSeries("",ep,0,self.delta_t,_lal.SecondUnit,len(self))
        elif self._data.dtype == _numpy.complex64:
            lal_data = _lal.CreateCOMPLEX8TimeSeries("",ep,0,self.delta_t,_lal.SecondUnit,len(self))
        elif self._data.dtype == _numpy.complex128:
            lal_data = _lal.CreateCOMPLEX16TimeSeries("",ep,0,self.delta_t,_lal.SecondUnit,len(self))

        lal_data.data.data[:] = self.numpy()

        return lal_data
      
    def save_to_wav(self, file_name):
        """ Save this time series to a wav format audio file.
        
        Parameters
        ----------
        file_name : string
             The output file name
        """
        scaled = _numpy.int16(self.numpy()/max(abs(self)) * 32767)
        write_wav(file_name, self.sample_rate, scaled)

    def save(self, path, group = None):
        """
        Save time series to a Numpy .npy, hdf, or text file. The first column
        contains the sample times, the second contains the values.
        In the case of a complex time series saved as text, the imaginary
        part is written as a third column. When using hdf format, the data is stored
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
            output = _numpy.vstack((self.sample_times.numpy(), self.numpy())).T
            _numpy.save(path, output)
        elif ext == '.txt':
            if self.kind == 'real':
                output = _numpy.vstack((self.sample_times.numpy(),
                                        self.numpy())).T
            elif self.kind == 'complex':
                output = _numpy.vstack((self.sample_times.numpy(),
                                        self.numpy().real,
                                        self.numpy().imag)).T
            _numpy.savetxt(path, output)
        elif ext =='.hdf':
            key = 'data' if group is None else group
            f = h5py.File(path)
            ds = f.create_dataset(key, data=self.numpy(), compression='gzip',
                                  compression_opts=9, shuffle=True)
            ds.attrs['start_time'] = float(self.start_time)
            ds.attrs['delta_t'] = float(self.delta_t)
        else:
            raise ValueError('Path must end with .npy, .txt or .hdf')
                
    @_nocomplex
    def to_frequencyseries(self, delta_f=None):
        """ Return the Fourier transform of this time series
        
        Parameters
        ----------
        delta_f : {None, float}, optional
            The frequency resolution of the returned frequency series. By 
        default the resolution is determined by the duration of the timeseries.
        
        Returns
        -------        
        FrequencySeries: 
            The fourier transform of this time series. 
        """
        from pycbc.fft import fft
        if not delta_f:
            delta_f = 1.0 / self.duration

        # add 0.5 to round integer
        tlen  = int(1.0 / delta_f / self.delta_t + 0.5)
        flen = tlen / 2 + 1

        if tlen < len(self):
            raise ValueError("The value of delta_f (%s) would be "
                             "undersampled. Maximum delta_f "
                             "is %s." % (delta_f, 1.0 / self.duration))
        if not delta_f:
            tmp = self
        else:
            tmp = TimeSeries(zeros(tlen, dtype=self.dtype), 
                             delta_t=self.delta_t, epoch=self.start_time)
            tmp[:len(self)] = self[:]
        
        f = FrequencySeries(zeros(flen, 
                           dtype=complex_same_precision_as(self)),
                           delta_f=delta_f)
        fft(tmp, f)
        return f

def load_timeseries(path, group=None):
    """
    Load a TimeSeries from a .hdf, .txt or .npy file. The
    default data types will be double precision floating point.

    Parameters
    ----------
    path: string
        source file path. Must end with either .npy or .txt.

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
        data = _numpy.load(path)    
    elif ext == '.txt':
        data = _numpy.loadtxt(path)
    elif ext == '.hdf':
        key = 'data' if group is None else group
        f = h5py.File(path)
        data = f[key][:]
        series = TimeSeries(data, delta_t=f[key].attrs['delta_t'],
                                  epoch=f[key].attrs['start_time']) 
        f.close()
        return series
    else:
        raise ValueError('Path must end with .npy, .hdf, or .txt')
        
    if data.ndim == 2:
        delta_t = (data[-1][0] - data[0][0]) / (len(data)-1)
        epoch = _lal.LIGOTimeGPS(data[0][0])
        return TimeSeries(data[:,1], delta_t=delta_t, epoch=epoch)
    elif data.ndim == 3:
        delta_t = (data[-1][0] - data[0][0]) / (len(data)-1)
        epoch = _lal.LIGOTimeGPS(data[0][0])
        return TimeSeries(data[:,1] + 1j*data[:,2],
                          delta_t=delta_t, epoch=epoch)
    else:
        raise ValueError('File has %s dimensions, cannot convert to Array, \
                          must be 2 (real) or 3 (complex)' % data.ndim)
