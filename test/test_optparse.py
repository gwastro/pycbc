# Copyright (C) 2026 Alexander Harvey Nitz
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
Unit tests for robust integer parsing and optparse actions in PyCBC.
"""

import argparse
import unittest
from pycbc.types.optparse import (
    to_int,
    MultiDetOptionAction,
    MultiDetMultiColonOptionAction,
    MultiDetOptionAppendAction,
    DictOptionAction,
    MultiDetDictOptionAction,
    positive_int,
    nonnegative_int,
)
from pycbc.types.config import InterpolatingConfigParser


class TestToInt(unittest.TestCase):
    def test_exact_large_integers(self):
        """Verify arbitrary-precision integer conversion without precision loss."""
        # IEEE 754 float64 boundary 2^53 = 9007199254740992
        # int(float('9007199254740993')) would silently round to 9007199254740992
        val_boundary = 9007199254740993
        self.assertEqual(to_int("9007199254740993"), val_boundary)
        self.assertEqual(to_int(val_boundary), val_boundary)

        # GPS nanoseconds (~1.2e18 >> 2^53)
        val_gps = 1187000000123456789
        self.assertEqual(to_int("1187000000123456789"), val_gps)
        self.assertEqual(to_int(val_gps), val_gps)

        # Large 64-bit random seeds
        val_seed = 1234567890123456789
        self.assertEqual(to_int("1234567890123456789"), val_seed)
        self.assertEqual(to_int(val_seed), val_seed)

        # 64-bit bitmask ((1 << 62) + 5)
        val_mask = 4611686018427387909
        self.assertEqual(to_int("4611686018427387909"), val_mask)

        # Confirm the bug that naive int(float(...)) produces
        self.assertNotEqual(int(float("9007199254740993")), val_boundary)
        self.assertNotEqual(int(float("1187000000123456789")), val_gps)

    def test_whole_number_float_strings(self):
        """Verify float strings representing exact whole numbers convert cleanly."""
        self.assertEqual(to_int("2048.0"), 2048)
        self.assertIsInstance(to_int("2048.0"), int)

        self.assertEqual(to_int("1e3"), 1000)
        self.assertIsInstance(to_int("1e3"), int)

        self.assertEqual(to_int("-5.0"), -5)
        self.assertEqual(to_int("1.0e6"), 1000000)
        self.assertEqual(to_int("0.0"), 0)
        self.assertEqual(to_int("-0.0"), 0)

        # Float objects
        self.assertEqual(to_int(2048.0), 2048)
        self.assertEqual(to_int(0.0), 0)
        self.assertEqual(to_int(-10.0), -10)

    def test_rejection_of_fractional_numbers(self):
        """Verify fractional strings and floats are strictly rejected with ValueError."""
        with self.assertRaises(ValueError):
            to_int("2048.5")

        with self.assertRaises(ValueError):
            to_int("1.5")

        with self.assertRaises(ValueError):
            to_int("-0.7")

        with self.assertRaises(ValueError):
            to_int("1e-3")

        with self.assertRaises(ValueError):
            to_int(2048.5)

        with self.assertRaises(ValueError):
            to_int(1.5)

    def test_rejection_of_invalid_inputs(self):
        """Verify invalid strings, non-numeric values, and infinities raise ValueError."""
        for invalid in ["abc", "", "   ", "2048.0.0", "inf", "-inf", "nan"]:
            with self.assertRaises(ValueError):
                to_int(invalid)

        with self.assertRaises(ValueError):
            to_int(None)


class TestMultiDetOptionAction(unittest.TestCase):
    def test_multidet_int_large_values(self):
        """Verify MultiDetOptionAction parses large integers > 2^53 without precision loss."""
        parser = argparse.ArgumentParser()
        parser.add_argument(
            "--fake-strain-seed",
            type=int,
            nargs="+",
            action=MultiDetOptionAction,
        )
        args = parser.parse_args([
            "--fake-strain-seed",
            "H1:9007199254740993",
            "L1:1187000000123456789",
        ])
        self.assertEqual(args.fake_strain_seed["H1"], 9007199254740993)
        self.assertEqual(args.fake_strain_seed["L1"], 1187000000123456789)

    def test_multidet_int_float_strings(self):
        """Verify MultiDetOptionAction accepts float strings for integer options."""
        parser = argparse.ArgumentParser()
        parser.add_argument(
            "--sample-rate",
            type=int,
            nargs="+",
            action=MultiDetOptionAction,
        )
        args = parser.parse_args([
            "--sample-rate",
            "H1:2048.0",
            "L1:1e3",
        ])
        self.assertEqual(args.sample_rate["H1"], 2048)
        self.assertEqual(args.sample_rate["L1"], 1000)

    def test_multidet_int_global_default(self):
        """Verify MultiDetOptionAction accepts global single-value float strings."""
        parser = argparse.ArgumentParser()
        parser.add_argument(
            "--sample-rate",
            type=int,
            nargs="+",
            action=MultiDetOptionAction,
        )
        args = parser.parse_args(["--sample-rate", "2048.0"])
        self.assertEqual(args.sample_rate["H1"], 2048)
        self.assertEqual(args.sample_rate["V1"], 2048)

    def test_multidet_int_fractional_rejected(self):
        """Verify MultiDetOptionAction rejects fractional inputs for integer options."""
        parser = argparse.ArgumentParser()
        parser.add_argument(
            "--pad-data",
            type=int,
            nargs="+",
            action=MultiDetOptionAction,
        )
        with self.assertRaises(ValueError):
            parser.parse_args(["--pad-data", "H1:1.5"])

        with self.assertRaises(ValueError):
            parser.parse_args(["--pad-data", "2048.5"])


class TestOtherOptparseActions(unittest.TestCase):
    def test_multidet_multi_colon_action(self):
        """Verify MultiDetMultiColonOptionAction with type=int."""
        parser = argparse.ArgumentParser()
        parser.add_argument(
            "--rate",
            type=int,
            nargs="+",
            action=MultiDetMultiColonOptionAction,
        )
        args = parser.parse_args(["--rate", "H1:2048.0", "L1:1e3"])
        self.assertEqual(args.rate["H1"], 2048)
        self.assertEqual(args.rate["L1"], 1000)

        with self.assertRaises(ValueError):
            parser.parse_args(["--rate", "H1:2048.5"])

    def test_multidet_option_append_action(self):
        """Verify MultiDetOptionAppendAction with type=int."""
        parser = argparse.ArgumentParser()
        parser.add_argument(
            "--seed-list",
            type=int,
            nargs="+",
            action=MultiDetOptionAppendAction,
        )
        args = parser.parse_args([
            "--seed-list",
            "H1:1e3",
            "H1:2048.0",
        ])
        self.assertEqual(args.seed_list["H1"], [1000, 2048])

        with self.assertRaises(ValueError):
            parser.parse_args(["--seed-list", "H1:1.5"])

    def test_dict_option_action(self):
        """Verify DictOptionAction with type=int."""
        parser = argparse.ArgumentParser()
        parser.add_argument(
            "--params",
            type=int,
            nargs="+",
            action=DictOptionAction,
        )
        args = parser.parse_args(["--params", "rate:2048.0", "count:1e3"])
        self.assertEqual(args.params["rate"], 2048)
        self.assertEqual(args.params["count"], 1000)

        with self.assertRaises(ValueError):
            parser.parse_args(["--params", "rate:1.5"])

    def test_multidet_dict_option_action(self):
        """Verify MultiDetDictOptionAction with type=int."""
        parser = argparse.ArgumentParser()
        parser.add_argument(
            "--rates",
            type=int,
            nargs="+",
            action=MultiDetDictOptionAction,
        )
        args = parser.parse_args(["--rates", "H1:rate:2048.0"])
        self.assertEqual(args.rates["H1"]["rate"], 2048)

        with self.assertRaises(ValueError):
            parser.parse_args(["--rates", "H1:rate:1.5"])

    def test_positive_and_nonnegative_int(self):
        """Verify positive_int and nonnegative_int accept whole-number float strings."""
        self.assertEqual(positive_int("2048.0"), 2048)
        self.assertEqual(positive_int("1e3"), 1000)
        self.assertEqual(positive_int("9007199254740993"), 9007199254740993)

        self.assertEqual(nonnegative_int("0.0"), 0)
        self.assertEqual(nonnegative_int("0"), 0)
        self.assertEqual(nonnegative_int("1e2"), 100)

        with self.assertRaises(argparse.ArgumentTypeError):
            positive_int("0")

        with self.assertRaises(argparse.ArgumentTypeError):
            positive_int("-5")

        with self.assertRaises(argparse.ArgumentTypeError):
            nonnegative_int("-1")

        with self.assertRaises(argparse.ArgumentTypeError):
            positive_int("2048.5")

        with self.assertRaises(argparse.ArgumentTypeError):
            nonnegative_int("1.5")


class TestInterpolatingConfigParserGetInt(unittest.TestCase):
    def setUp(self):
        self.cp = InterpolatingConfigParser()
        self.ini_text = """
[workflow]
sample-rate = 2048.0
niterations = 1e6
large-seed = 9007199254740993
gps-nano = 1187000000123456789
negative-int = -5.0
fractional-val = 2048.5
"""
        self.cp.read_string(self.ini_text)

    def test_getint_whole_float_strings_and_precision(self):
        """Verify getint safely parses float strings and preserves arbitrary precision."""
        self.assertEqual(self.cp.getint("workflow", "sample-rate"), 2048)
        self.assertIsInstance(self.cp.getint("workflow", "sample-rate"), int)

        self.assertEqual(self.cp.getint("workflow", "niterations"), 1000000)
        self.assertIsInstance(self.cp.getint("workflow", "niterations"), int)

        self.assertEqual(self.cp.getint("workflow", "large-seed"), 9007199254740993)
        self.assertEqual(self.cp.getint("workflow", "gps-nano"), 1187000000123456789)
        self.assertEqual(self.cp.getint("workflow", "negative-int"), -5)

    def test_getint_fractional_rejected(self):
        """Verify getint raises ValueError on fractional strings."""
        with self.assertRaises(ValueError):
            self.cp.getint("workflow", "fractional-val")

    def test_getint_fallback(self):
        """Verify getint returns fallback when option or section is missing."""
        self.assertEqual(self.cp.getint("workflow", "missing-option", fallback=42), 42)
        self.assertEqual(self.cp.getint("missing-section", "missing-option", fallback=99), 99)


if __name__ == "__main__":
    unittest.main()
