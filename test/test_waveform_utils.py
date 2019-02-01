import unittest
import numpy

from utils import simple_exit

from pycbc.waveform.utils import apply_fd_time_shift
from pycbc.types import (FrequencySeries, TimeSeries)


class TestFDTimeShift(unittest.TestCase):
    """Tests ``apply_fd_time_shift``.""" 
    def setUp(self):
        # we'll use a sine wave time series to do the testing, with the
        # the segment length such that an interger number of cycles fit, so
        # we don't have to worry about boundary effects
        self.freq = 128
        self.sample_rate = 4096
        self.seglen = 1
        ncycles = self.freq * self.seglen
        t = numpy.linspace(0, ncycles*2*numpy.pi,
                           num=self.sample_rate*self.seglen,
                           endpoint=False)
        self.time_series = TimeSeries(t, delta_t=1./self.sample_rate, epoch=0)
        self.tdsinx = numpy.sin(self.time_series)
        self.fdsinx = self.tdsinx.to_frequencyseries()

    def _shift_and_ifft(self, fdsinx, tshift, fseries=None):
        """Calls apply_fd_time_shift, and iFFTs to the time domain.
        """
        start_time = self.time_series.start_time
        tdshift = apply_fd_time_shift(fdsinx, start_time+tshift,
                                      fseries=fseries)
        if not isinstance(tdshift, FrequencySeries):
            # cast to FrequencySeries so time series will work
            tdshift = FrequencySeries(tdshift, delta_f=fdsinx.delta_f,
                                      epoch=fdsinx.epoch)
        return tdshift.to_timeseries()

    def _test_apply_fd_time_shift(self, fdsinx, fseries=None):
        """Tests ``apply_fd_time_shift`` with the given fdseries.

        If ``fdsinx`` is a FrequencySeries, this will test the shift code
        written in C. Otherwise, this will test the numpy version.
        """
        # shift by -pi/2: should be the same as the cosine
        start_time = self.time_series.start_time
        tshift = 1./(4*self.freq)
        tdshift = self._shift_and_ifft(fdsinx, -tshift, fseries=fseries)
        # check
        comp = numpy.cos(self.time_series)
        self.assertTrue(numpy.isclose(tdshift, comp).all())
        # shift by +pi/2: should be the same as the -cosine
        tdshift = self._shift_and_ifft(fdsinx, tshift, fseries=fseries)
        self.assertTrue(numpy.isclose(tdshift, -1*comp).all())
        # shift by a non-integer fraction of the period; we'll do this by
        # shifting by a prime number times dt / 3
        # forward:
        tshift = 193 * self.time_series.delta_t / 3.
        tdshift = self._shift_and_ifft(fdsinx, tshift, fseries=fseries)
        comp = numpy.sin(self.time_series - 2*numpy.pi*self.freq*tshift)
        self.assertTrue(numpy.isclose(tdshift, comp).all())
        # backward:
        tdshift = self._shift_and_ifft(fdsinx, -tshift, fseries=fseries)
        comp = numpy.sin(self.time_series + 2*numpy.pi*self.freq*tshift)
        self.assertTrue(numpy.isclose(tdshift, comp).all())

    def test_fd_time_shift(self):
        """Applies shifts to fdsinx using C code, and compares the result to
        applying the shift directly in the time domain.
        """
        self._test_apply_fd_time_shift(self.fdsinx)

    def test_fseries_time_shift(self):
        """Applies shifts to fdsinx using numpy code, and compares the result
        to applying the shift directly in the time domain.
        """
        fdsinx = self.fdsinx.numpy()
        # attach things will work needed to make this work
        fdsinx.delta_f = self.fdsinx.delta_f
        fdsinx.epoch = self.fdsinx.epoch
        fdsinx.precision = self.fdsinx.precision
        fseries = self.fdsinx.sample_frequencies.numpy()
        self._test_apply_fd_time_shift(fdsinx, fseries)

suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestFDTimeShift))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
