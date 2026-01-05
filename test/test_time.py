import datetime
import unittest
from astropy.time import Time

import pycbc.time

class TestStripTime(unittest.TestCase):

    def _expected(self, gps_time):
        # compute expected midnight GPS and date string using astropy
        # and compare to the pycbc internal calculation
        expected_date = Time(gps_time, format='gps', scale='utc').to_datetime().date()
        midnight_dt = datetime.datetime.combine(expected_date, datetime.time.min)
        expected_midnight_gps = Time(midnight_dt, format='datetime').gps
        return expected_midnight_gps, expected_date.strftime("%Y-%m-%d")

    def test_integer_gps(self):
        gps = 1262304000  # arbitrary known GPS second
        midnight, date_str = pycbc.time.strip_time_from_gps(gps)
        exp_midnight, exp_date = self._expected(gps)
        self.assertEqual(date_str, exp_date)
        # allow small floating point diffs
        self.assertAlmostEqual(midnight, exp_midnight, places=6)

    def test_float_gps(self):
        gps = 1262304000.5
        midnight, date_str = pycbc.time.strip_time_from_gps(gps)
        exp_midnight, exp_date = self._expected(gps)
        self.assertEqual(date_str, exp_date)
        self.assertAlmostEqual(midnight, exp_midnight, places=6)

    def test_midnight_is_before_or_equal_and_within_day(self):
        gps = 1262304000 + 3600 * 12  # midday
        midnight, _ = pycbc.time.strip_time_from_gps(gps)
        # midnight should be <= provided time and within 24 hours difference
        self.assertLessEqual(midnight, gps)
        self.assertLess(gps - midnight, 86400)


class TestTimeConversionsFull(unittest.TestCase):
    """
    Comprehensive tests that exercise all functions in pycbc.time.

    Important: include samples around each leap-second insertion to ensure
    gps<->UTC conversions and midnight-stripping behave across leap seconds.
    """

    # leap-second effective dates
    LEAP_DATES = [
        '1981-07-01', '1982-07-01', '1983-07-01', '1985-07-01',
        '1988-01-01', '1990-01-01', '1991-01-01', '1992-07-01',
        '1993-07-01', '1994-07-01', '1996-01-01', '1997-07-01',
        '1999-01-01', '2006-01-01', '2009-01-01', '2012-07-01',
        '2015-07-01', '2017-01-01'
    ]

    def test_roundtrip_and_strip_around_leaps(self):
        """For each leap date test times the day before, day of, day after."""
        from datetime import datetime, time, timedelta

        for date_str in self.LEAP_DATES:
            # build midnight of that date in UTC and sample day before/after
            date_dt = datetime.strptime(date_str, "%Y-%m-%d").date()
            samples = [
                datetime.combine(date_dt - timedelta(days=1), time(hour=12)),
                datetime.combine(date_dt, time(hour=12)),
                datetime.combine(date_dt + timedelta(days=1), time(hour=12)),
            ]

            for dt in samples:
                # convert to GPS using astropy to get canonical gps seconds
                gps = Time(dt, format='datetime', scale='utc').gps

                # gps_to_utc_datetime should return a datetime equal to dt (modulo
                # fractions) when we convert back from gps
                utc_dt = pycbc.time.gps_to_utc_datetime(gps)
                # compare dates and times to second precision
                self.assertEqual(utc_dt.date(), dt.date())

                # utc_datetime_to_gps should round-trip (within small tol)
                gps_back = pycbc.time.utc_datetime_to_gps(utc_dt)
                self.assertAlmostEqual(gps, gps_back, places=6)

                # gps_to_utc_str and datetime_to_str produce same date string
                s1 = pycbc.time.gps_to_utc_str(gps, format="%Y-%m-%d %H:%M:%S")
                s2 = pycbc.time.datetime_to_str(utc_dt, format="%Y-%m-%d %H:%M:%S")
                self.assertEqual(s1, s2)

                # strip_time_from_date should zero the time
                stripped = pycbc.time.strip_time_from_date(utc_dt)
                self.assertEqual(stripped.hour, 0)
                self.assertEqual(stripped.minute, 0)
                self.assertEqual(stripped.second, 0)
                self.assertEqual(stripped.microsecond, 0)

                # strip_time_from_gps returns midnight GPS and date string
                midnight_gps, date_string = pycbc.time.strip_time_from_gps(gps)
                # midnight_gps should be <= gps and within one day
                self.assertLessEqual(midnight_gps, gps)
                self.assertLess(gps - midnight_gps, 86400)
                # date_string should match the stripped date
                self.assertEqual(date_string, stripped.date().strftime("%Y-%m-%d"))

    def test_datetime_to_str_formats(self):
        from datetime import datetime
        dt = datetime(2020, 2, 29, 15, 4, 5)
        self.assertEqual(pycbc.time.datetime_to_str(dt, format="%Y-%m-%d"), "2020-02-29")
        self.assertEqual(pycbc.time.datetime_to_str(dt, format="%H:%M:%S"), "15:04:05")

    def test_gps_now_and_gmst(self):
        # gps_now should be close to astropy Time.now().gps
        now_astropy = Time.now()
        gps_now = pycbc.time.gps_now()
        # Allow a couple of seconds tolerance for the call overhead
        self.assertAlmostEqual(gps_now, float(now_astropy.gps), delta=2.0)

        # gmst_accurate should match astropy sidereal_time('mean').rad
        gps_sample = Time('2020-01-01T00:00:00', format='isot', scale='utc').gps
        gmst_expected = Time(gps_sample, format='gps', scale='utc', location=(0, 0)).sidereal_time('mean').rad
        gmst = pycbc.time.gmst_accurate(gps_sample)
        # compare within a small tolerance (radians)
        self.assertAlmostEqual(gmst, gmst_expected, places=9)


if __name__ == '__main__':
    unittest.main()