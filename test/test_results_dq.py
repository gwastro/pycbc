
import datetime
import unittest

from pycbc import results



class TestGetSummaryPageLink(unittest.TestCase):


    def test_valid_datetime_returns_string(self):
        out = results.get_summary_page_link('H1', datetime.date(2020, 1, 2))
        self.assertIsInstance(out, str)
        self.assertIn('Summary', out)

    def test_valid_date_tuple_returns_string(self):
        out = results.get_summary_page_link('H1', (2020, 1, 2))
        self.assertIsInstance(out, str)
        self.assertIn('Summary', out)

    def test_bad_input_float_raises_type_error(self):
        with self.assertRaises(TypeError) as cm:
            results.get_summary_page_link('H1', 123.456)
        msg = str(cm.exception)
        self.assertIn('utc_time must be a datetime/date or a sequence', msg)
        self.assertIn('got', msg)

    def test_bad_input_short_sequence_raises_type_error(self):
        with self.assertRaises(TypeError):
            results.get_summary_page_link('H1', (2020, 1))

    def test_bad_input_strings_raises_type_error(self):
        with self.assertRaises(TypeError):
            results.get_summary_page_link('H1', ('year', 'month', 'day'))


if __name__ == '__main__':
    unittest.main()