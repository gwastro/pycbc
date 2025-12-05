import datetime
import unittest
from astropy.time import Time

from pycbc.live import plotting_utils as pu

class TestStripTime(unittest.TestCase):

    def _expected(self, gps_time):
        # compute expected midnight GPS and date string using astropy
        expected_date = Time(gps_time, format='gps', scale='utc').to_datetime().date()
        midnight_dt = datetime.datetime.combine(expected_date, datetime.time.min)
        expected_midnight_gps = Time(midnight_dt, format='datetime').gps
        return expected_midnight_gps, expected_date.strftime("%Y-%m-%d")

    def test_integer_gps(self):
        gps = 1262304000  # arbitrary known GPS second
        midnight, date_str = pu.strip_time(gps)
        exp_midnight, exp_date = self._expected(gps)
        self.assertEqual(date_str, exp_date)
        # allow small floating point diffs
        self.assertAlmostEqual(midnight, exp_midnight, places=6)

    def test_float_gps(self):
        gps = 1262304000.5
        midnight, date_str = pu.strip_time(gps)
        exp_midnight, exp_date = self._expected(gps)
        self.assertEqual(date_str, exp_date)
        self.assertAlmostEqual(midnight, exp_midnight, places=6)

    def test_midnight_is_before_or_equal_and_within_day(self):
        gps = 1262304000 + 3600 * 12  # midday
        midnight, _ = pu.strip_time(gps)
        # midnight should be <= provided time and within 24 hours difference
        self.assertLessEqual(midnight, gps)
        self.assertLess(gps - midnight, 86400)


if __name__ == '__main__':
    unittest.main()
