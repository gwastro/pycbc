
import os
import importlib.util
import datetime
import unittest


def load_dq_module():
    # Path to the dq.py file inside the env_pycbc_test_src workspace
    base = os.path.dirname(__file__)
    dq_path = os.path.abspath(os.path.join(base, '..', 'pycbc', 'results', 'dq.py'))
    spec = importlib.util.spec_from_file_location('temp_dq', dq_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TestGetSummaryPageLink(unittest.TestCase):
    def setUp(self):
        self.dq = load_dq_module()

    def test_valid_datetime_returns_string(self):
        out = self.dq.get_summary_page_link('H1', datetime.date(2020, 1, 2))
        self.assertIsInstance(out, str)
        self.assertIn('Summary', out)

    def test_valid_date_tuple_returns_string(self):
        out = self.dq.get_summary_page_link('H1', (2020, 1, 2))
        self.assertIsInstance(out, str)
        self.assertIn('Summary', out)

    def test_bad_input_float_raises_type_error(self):
        with self.assertRaises(TypeError) as cm:
            self.dq.get_summary_page_link('H1', 123.456)
        msg = str(cm.exception)
        self.assertIn('utc_time must be a datetime/date or a sequence', msg)
        self.assertIn('got', msg)

    def test_bad_input_short_sequence_raises_type_error(self):
        with self.assertRaises(TypeError):
            self.dq.get_summary_page_link('H1', (2020, 1))

    def test_bad_input_strings_raises_type_error(self):
        with self.assertRaises(TypeError):
            self.dq.get_summary_page_link('H1', ('year', 'month', 'day'))


if __name__ == '__main__':
    unittest.main()