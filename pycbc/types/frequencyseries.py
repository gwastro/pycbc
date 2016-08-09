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

    Attributes
    ----------
    delta_f
    epoch
    sample_frequencies
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
    delta_f = property(get_delta_f)

    def get_epoch(self):
        """Return frequency series epoch as a LIGOTimeGPS.
        """
        return self._epoch
    epoch = property(get_epoch)

    def get_sample_frequencies(self):
        """Return an Array containing the sample frequencies.
        """
        return Array(range(len(self))) * self._delta_f
    sample_frequencies = property(get_sample_frequencies)

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
            from pylal import series as lalseries
            from glue.ligolw import utils
            assert(self.kind == 'real')
            output = self.lal()
            # When writing in this format we must *not* have the 0 values at
            # frequencies less than flow. To resolve this we set the first
            # non-zero value < flow.
            data_lal = output.data.data
            first_idx = _numpy.argmax(data_lal>0)
            if not first_idx == 0:
                data_lal[:first_idx] = data_lal[first_idx]
            psddict = {ifo: output}
            utils.write_filename(lalseries.make_psd_xmldoc(psddict), path,
                                 gz=path.endswith(".gz"))
        elif ext =='.hdf':
            key = 'data' if group is None else group
            f = h5py.File(path)
            ds = f.create_dataset(key, data=self.numpy(), compression='gzip',
                                  compression_opts=9, shuffle=True)
            ds.attrs['epoch'] = float(self.epoch)
            ds.attrs['delta_f'] = float(self.delta_f)
        else:
            raise ValueError('Path must end with .npy, .txt, .xml, .xml.gz '
                             'or .hdf')

    @_noreal
    def to_timeseries(self, delta_t=None):
        """ Return the Fourier transform of this time series
        
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
        flen = tlen / 2 + 1
        
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
        return f
            

def load_frequencyseries(path, group=None):
    """
    Load a FrequencySeries from a .hdf, .txt or .npy file. The
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
        series = FrequencySeries(data, delta_f=f[key].attrs['delta_f'],
                                       epoch=f[key].attrs['epoch']) 
        f.close()
        return series
    else:
        raise ValueError('Path must end with .npy, .hdf, or .txt')
        
    if data.ndim == 2:
        delta_f = (data[-1][0] - data[0][0]) / (len(data)-1)
        epoch = _lal.LIGOTimeGPS(data[0][0])
        return FrequencySeries(data[:,1], delta_f=delta_f)
    elif data.ndim == 3:
        delta_f = (data[-1][0] - data[0][0]) / (len(data)-1)
        epoch = _lal.LIGOTimeGPS(data[0][0])
        return FrequencySeries(data[:,1] + 1j*data[:,2], delta_f=delta_f)
    else:
        raise ValueError('File has %s dimensions, cannot convert to Array, \
                          must be 2 (real) or 3 (complex)' % data.ndim)
