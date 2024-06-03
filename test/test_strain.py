"""
These unit tests are for the pycbc.strain.strain module
"""
import numpy
from pycbc.types import TimeSeries
from pycbc.strain.strain import (
    execute_cached_fft,
    execute_cached_ifft,
)
import unittest

from utils import simple_exit


class TestStrain(unittest.TestCase):

    def setUp(self):

        self.rng = numpy.random.default_rng()
        self.td_data = TimeSeries(
            self.rng.normal(size=100), delta_t=0.2, epoch=1123456789.6,
        )
        self.fd_data = self.td_data.to_frequencyseries()
        # Tolerance for float64
        self.tol = 1e-14

    def test_cached_fft(self):
        fd_data = execute_cached_fft(
            self.td_data,
            uid=87651,
            copy_output=True,
        )
        self.assertTrue(
            fd_data.almost_equal_norm(
                self.fd_data, tol=self.tol, dtol=self.tol
            )
        )

    def test_cached_ifft(self):
        td_data = execute_cached_ifft(
            self.fd_data,
            uid=87652,
            copy_output=True,
        )
        self.assertTrue(
            td_data.almost_equal_norm(
                self.td_data, tol=self.tol, dtol=self.tol
            )
        )


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestStrain))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
