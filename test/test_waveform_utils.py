import unittest
import numpy

from utils import simple_exit

from pycbc import cosmology
from pycbc.waveform import get_td_waveform, get_fd_waveform
from pycbc.waveform.utils import (
    apply_fd_time_shift,
    redshift_waveform,
)
from pycbc.types import (TimeSeries)


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
        return tdshift.to_timeseries()

    def _test_apply_fd_time_shift(self, fdsinx, fseries=None, atol=1e-8):
        """Tests ``apply_fd_time_shift`` with the given fdseries.

        If ``fdsinx`` is a FrequencySeries, this will test the shift code
        written in C. Otherwise, this will test the numpy version.

        Parameters
        ----------
        fdsinx : FrequencySeries
            The frequency series to shift and test.
        fseires : array, optional
            Array of the sample frequencies of ``fdsinx``. This is only needed
            for the numpy version.
        atol : float, optional
            The absolute tolerance for the comparison test. See
            ``numpy.isclose`` for details.
        """
        # shift by -pi/2: should be the same as the cosine
        tshift = 1./(4*self.freq)
        tdshift = self._shift_and_ifft(fdsinx, -tshift, fseries=fseries)
        # check
        comp = numpy.cos(self.time_series)
        if tdshift.precision == 'single':
            # cast to single
            comp = comp.astype(numpy.float32)
        self.assertTrue(numpy.isclose(tdshift, comp, atol=atol).all())
        # shift by +pi/2: should be the same as the -cosine
        tdshift = self._shift_and_ifft(fdsinx, tshift, fseries=fseries)
        self.assertTrue(numpy.isclose(tdshift, -1*comp, atol=atol).all())
        # shift by a non-integer fraction of the period; we'll do this by
        # shifting by a prime number times dt / 3
        # forward:
        tshift = 193 * self.time_series.delta_t / 3.
        tdshift = self._shift_and_ifft(fdsinx, tshift, fseries=fseries)
        comp = numpy.sin(self.time_series - 2*numpy.pi*self.freq*tshift)
        if tdshift.precision == 'single':
            # cast to single
            comp = comp.astype(numpy.float32)
        self.assertTrue(numpy.isclose(tdshift, comp, atol=atol).all())
        # backward:
        tdshift = self._shift_and_ifft(fdsinx, -tshift, fseries=fseries)
        comp = numpy.sin(self.time_series + 2*numpy.pi*self.freq*tshift)
        if tdshift.precision == 'single':
            # cast to single
            comp = comp.astype(numpy.float32)
        self.assertTrue(numpy.isclose(tdshift, comp, atol=atol).all())

    def test_fd_time_shift(self):
        """Applies shifts to fdsinx using cython code, and compares the
        result to applying the shift directly in the time domain.
        """
        self._test_apply_fd_time_shift(self.fdsinx)

    def test_fd_time_shift32(self):
        """Tests the cython code using single precision.
        """
        # we need to increase the tolerance on isclose
        self._test_apply_fd_time_shift(self.fdsinx.astype(numpy.complex64),
                                       atol=1e-4)

    def test_fseries_time_shift(self):
        """Applies shifts to fdsinx using numpy code, and compares the result
        to applying the shift directly in the time domain.
        """
        fdsinx = self.fdsinx.copy()
        fseries = self.fdsinx.sample_frequencies.numpy()
        self._test_apply_fd_time_shift(fdsinx, fseries)


class TestRedshiftWaveform(unittest.TestCase):
    """Tests ``redshift_waveform`` against detector-frame generation.
    
    Specifically, this tests that a waveform generated in the source frame and
    then redshifted using ``redshift_waveform`` matches a waveform generated
    directly in the detector frame with redshifted masses. This is done for
    both TD and FD waveforms.
    """

    def setUp(self):
        self.srcm1 = 30.0
        self.srcm2 = 20.0
        self.distance = 5672.12 #4000.0
        self.z = cosmology.redshift(self.distance)

        # Detector-frame settings.
        self.flow = 30.0
        self.sample_rate = 4096.0
        self.seglen = 64.0

        # Source-frame settings.
        self.srcflow = self.flow * (1 + self.z)
        self.srcsr = self.sample_rate * (1 + self.z)
        self.srcseglen = self.seglen / (1 + self.z)

        self.detm1 = self.srcm1 * (1 + self.z)
        self.detm2 = self.srcm2 * (1 + self.z)

    def _relative_l2_error(self, test, ref):
        """Returns relative L2 norm of ``test - ref``."""
        return numpy.linalg.norm(test - ref) / numpy.linalg.norm(ref)


    def _check_epochs(self, redshifted_hp, det_hp, err=0.001):
        """Checks that the epochs of the two waveforms are close enough.

        Small differences in the epochs may arise due to floating point
        errors. That may cause a failure when we compute the relative L2
        error, so we'll check that the epochs are close enough.
        """
        isclose = numpy.isclose(redshifted_hp._epoch, det_hp._epoch,
                                      rtol=0., atol=err*det_hp.delta_t)
        self.assertTrue(isclose,
                        msg=f"Epochs differ by more than {err*det_hp.delta_t}:"
                            f" |redshifted - detector epoch| = "
                            f"{abs(redshifted_hp._epoch - det_hp._epoch)}")

    def test_td_redshift_matches_redshifted_masses(self):
        """Redshifting a source-frame TD waveform matches detector-frame TD."""
        src_hp, _ = get_td_waveform(
            approximant='SEOBNRv4',
            mass1=self.srcm1,
            mass2=self.srcm2,
            distance=self.distance,
            delta_t=1.0 / self.srcsr,
            f_lower=self.srcflow,
        )
        redshifted_hp = redshift_waveform(src_hp, self.z)

        det_hp, _ = get_td_waveform(
            approximant='SEOBNRv4',
            mass1=self.detm1,
            mass2=self.detm2,
            distance=self.distance,
            delta_t=1.0 / self.sample_rate,
            f_lower=self.flow,
        )
        self._check_epochs(redshifted_hp, det_hp)
        # if passed, set the redshifted_hp epoch to be the same as the det_hp
        # epoch, so that we can compare the waveforms directly
        redshifted_hp._epoch = det_hp._epoch
        relerr = self._relative_l2_error(redshifted_hp, det_hp)
        self.assertLess(relerr, 1e-3)


    def test_fd_redshift_matches_redshifted_masses(self):
        """Redshifting a source-frame FD waveform matches detector-frame FD."""
        src_hptilde, _ = get_fd_waveform(
            approximant='IMRPhenomXPHM',
            mass1=self.srcm1,
            mass2=self.srcm2,
            distance=self.distance,
            delta_f=1.0 / self.srcseglen,
            f_lower=self.srcflow,
            f_final=self.srcsr / 2.0,
        )
        redshifted_hptilde = redshift_waveform(src_hptilde, self.z)

        det_hptilde, _ = get_fd_waveform(
            approximant='IMRPhenomXPHM',
            mass1=self.detm1,
            mass2=self.detm2,
            distance=self.distance,
            delta_f=1.0 / self.seglen,
            f_lower=self.flow,
            f_final=self.sample_rate / 2.0,
        )

        redshifted_hp = redshifted_hptilde.to_timeseries()
        det_hp = det_hptilde.to_timeseries()
        self._check_epochs(redshifted_hp, det_hp)
        # if passed, set the redshifted_hp epoch to be the same as the det_hp
        # epoch, so that we can compare the waveforms directly
        redshifted_hp._epoch = det_hp._epoch
        relerr = self._relative_l2_error(redshifted_hp, det_hp)
        self.assertLess(relerr, 1e-3)


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestFDTimeShift))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestRedshiftWaveform))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
