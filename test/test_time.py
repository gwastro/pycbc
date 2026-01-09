import datetime as dt
import unittest
from astropy.time import Time
from datetime import datetime, timezone

import pycbc.time
from pycbc.constants import PI

class TestStripTime(unittest.TestCase):

    def _expected(self, gps_time):
        # compute expected midnight GPS and date string using astropy
        # and compare to the pycbc internal calculation
        expected_date = Time(gps_time, format='gps', scale='utc').to_datetime().date()
        midnight_dt = datetime.combine(expected_date, dt.time.min)
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

# Truth conversions for certain times around leap seconds.
# These are converted using lal_tconvert,
# rather than astropy conversions used internally
TRUTH_CONVERSIONS = {
    '1981-06-30 01:02:03': (46746123.0, 5.123104000000),
    '1981-07-01 02:03:04': (46836185.0, 5.407271000000),
    '1982-06-30 04:05:06': (78293107.0, 5.919831000000),
    '1982-07-01 05:06:07': (78383169.0, 6.203998000000),
    '1983-06-30 07:08:09': (109840091.0, 0.433372000000),
    '1983-07-01 08:09:10': (109930153.0, 0.717539000000),
    '1985-06-30 10:11:12': (173009475.0, 1.243135000000),
    '1985-07-01 11:12:13': (173099537.0, 1.527303000000),
    '1987-12-31 13:14:15': (251990059.0, 5.201009000000),
    '1988-01-01 14:15:16': (252080121.0, 5.485177000000),
    '1989-12-31 16:17:18': (315159443.0, 6.010773000000),
    '1990-01-01 17:18:19': (315249505.0, 0.011755000000),
    '1990-12-31 19:20:21': (346706427.0, 0.524314000000),
    '1991-01-01 20:21:22': (346796489.0, 0.808481000000),
    '1992-06-30 22:23:24': (393978211.0, 4.451949000000),
    '1992-07-01 23:24:25': (394068273.0, 4.736116000000),
    '1993-06-30 01:02:03': (425437331.0, 5.124716000000),
    '1993-07-01 02:03:04': (425527393.0, 5.408884000000),
    '1994-06-30 04:05:06': (456984315.0, 5.921443000000),
    '1994-07-01 05:06:07': (457074377.0, 6.205610000000),
    '1995-12-31 07:08:09': (504428899.0, 3.600298000000),
    '1996-01-01 08:09:10': (504518961.0, 3.884466000000),
    '1997-06-30 10:11:12': (551700683.0, 1.244748000000),
    '1997-07-01 11:12:13': (551790745.0, 1.528915000000),
    '1998-12-31 13:14:15': (599145267.0, 5.206788000000),
    '1999-01-01 14:15:16': (599235329.0, 5.490956000000),
    '2005-12-31 16:17:18': (820081051.0, 6.012923000000),
    '2006-01-01 17:18:19': (820171113.0, 0.013905000000),
    '2008-12-31 19:20:21': (914786435.0, 0.535335000000),
    '2009-01-01 20:21:22': (914876497.0, 0.819502000000),
    '2012-06-30 22:23:24': (1025130219.0, 4.454637000000),
    '2012-07-01 23:24:25': (1025220281.0, 4.738804000000),
    '2015-06-30 01:02:03': (1119661339.0, 5.119072000000),
    '2015-07-01 02:03:04': (1119751401.0, 5.403239000000),
    '2016-12-31 04:05:06': (1167192323.0, 2.815130000000),
    '2017-01-01 05:06:07': (1167282385.0, 3.099297000000),
}

class TestTimeConversionsFull(unittest.TestCase):
    """
    Comprehensive tests that exercise all functions in pycbc.time.

    Important: include samples around each leap-second insertion to ensure
    gps<->UTC conversions and midnight-stripping behave across leap seconds.
    """

    def test_roundtrip(self):
        """Use hard-coded authoritative conversions for UTC->GPS and back."""

        for utc_str, (gps, _) in TRUTH_CONVERSIONS.items():
            # gps_to_utc_str should recover the same UTC datetime
            test_str = pycbc.time.gps_to_utc_str(gps)
            self.assertEqual(utc_str, test_str)

            # converting the string to a datetime as then converting to gps
            # should match the hard-coded GPS value
            utc_dt = datetime.strptime(utc_str, "%Y-%m-%d %H:%M:%S")
            gps_back = pycbc.time.utc_datetime_to_gps(utc_dt)
            self.assertAlmostEqual(gps, gps_back, places=6)

            # gps_to_utc_str and datetime_to_str produce same date string
            s1 = pycbc.time.gps_to_utc_str(gps, format="%Y-%m-%d %H:%M:%S")
            s2 = pycbc.time.datetime_to_str(
                utc_dt,
                format="%Y-%m-%d %H:%M:%S"
            )
            self.assertEqual(s1, s2)


    def test_gmst_accurate(self):
        """Calculate GMST values and confirm they match hard-coded truth"""
        tol = 1e-4
        for (gps_val, gmst_expected) in TRUTH_CONVERSIONS.values():
            gmst = pycbc.time.gmst_accurate(gps_val)
            delta = abs((gmst - gmst_expected + PI) % (2*PI) - PI)
            self.assertLess(delta, tol)


    def test_strip_around_leaps(self):
        """Test stripping of time around leap seconds."""

        for utc_str, (gps, _) in TRUTH_CONVERSIONS.items():
            # strip_time_from_date should zero the time
            utc_dt = datetime.strptime(utc_str, "%Y-%m-%d %H:%M:%S")
            stripped = pycbc.time.strip_time_from_date(utc_dt)
            self.assertEqual(stripped.hour, 0)
            self.assertEqual(stripped.minute, 0)
            self.assertEqual(stripped.second, 0)
            self.assertEqual(stripped.microsecond, 0)
            self.assertEqual(
                utc_str.split()[0],
                pycbc.time.datetime_to_str(stripped, format="%Y-%m-%d")
            )

            # strip_time_from_gps returns midnight GPS and date string
            midnight_gps, date_string = pycbc.time.strip_time_from_gps(gps)
            # midnight_gps should be <= gps and within one day
            self.assertLessEqual(midnight_gps, gps)
            self.assertLess(gps - midnight_gps, 86400)
            # date_string should match the date from the hardcoded truth
            self.assertEqual(
                date_string,
                utc_str.split()[0]
            )


    def test_datetime_to_str_formats(self):
        """Confirm conversion from datetime to strings with different formats"""
        dt = datetime(2020, 2, 29, 15, 4, 5)
        self.assertEqual(
            pycbc.time.datetime_to_str(dt, format="%Y-%m-%d"),
            "2020-02-29"
        )
        self.assertEqual(
            pycbc.time.datetime_to_str(dt, format="%H:%M:%S"),
            "15:04:05"
        )


    def test_gps_now(self):
        """confirm gps_now converts to a UTC datetime close to the system UTC
        time"""

        gps_now = pycbc.time.gps_now()
        dt_from_gps = pycbc.time.gps_to_utc_datetime(gps_now)
        # Ensure dt_from_gps is timezone-aware in UTC for a fair comparison
        if dt_from_gps.tzinfo is None:
            dt_from_gps = dt_from_gps.replace(tzinfo=timezone.utc)

        now_dt = datetime.now(timezone.utc)
        # difference in seconds
        delta_secs = abs((dt_from_gps - now_dt).total_seconds())
        self.assertLessEqual(delta_secs, 2.0)



if __name__ == '__main__':
    unittest.main()