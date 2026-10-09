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

"""Unit tests for robust integer parsing and optparse actions in PyCBC."""

import argparse
import os
import sys
import unittest

# Ensure local repository pycbc is imported
_here = os.path.abspath(os.path.dirname(__file__))
_repo_root = os.path.dirname(_here)
if _repo_root not in sys.path:
    sys.path.insert(0, _repo_root)

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
    def test_to_int_integers_and_large_ints(self):
        """Verify basic ints and large integers >= 2^53 without float truncation."""
        # Basic ints
        for val, expected in [(42, 42), ("42", 42), ("-10", -10), (0, 0), ("0", 0)]:
            self.assertEqual(to_int(val), expected)
            self.assertIsInstance(to_int(val), int)

        # Exact large integers >= 2^53 (where float64 loses precision)
        large_vals = [
            ("9007199254740993", 9007199254740993),       # 2^53 + 1
            (1187000000123456789, 1187000000123456789),   # GPS nanosecond timestamp
            ("18446744073709551615", 18446744073709551615),  # 2^64 - 1 bitmask
        ]
        for val, expected in large_vals:
            res = to_int(val)
            self.assertEqual(res, expected)
            self.assertIsInstance(res, int)

    def test_to_int_whole_number_float_strings(self):
        """Verify whole-number float strings and scientific notation convert to int."""
        cases = [("2048.0", 2048), ("1e3", 1000), ("-5.0", -5), ("0.0", 0), ("+100.0", 100), (2048.0, 2048)]
        for s, expected in cases:
            res = to_int(s)
            self.assertEqual(res, expected)
            self.assertIsInstance(res, int)

    def test_to_int_fractional_and_invalid_rejected(self):
        """Verify fractional and invalid inputs raise ValueError."""
        for invalid in ["2048.5", "1.5", "-0.7", 2048.5, 1e-3, "abc", "", "inf", "-inf", "nan", None]:
            with self.assertRaises(ValueError):
                to_int(invalid)


class TestOptparseActions(unittest.TestCase):
    def test_multidet_option_action_int(self):
        """Verify MultiDetOptionAction accepts float strings, large ints, and global defaults."""
        parser = argparse.ArgumentParser()
        parser.add_argument("--rate", type=int, nargs="+", action=MultiDetOptionAction)
        parser.add_argument("--seed", type=int, nargs="+", action=MultiDetOptionAction)

        # Detector-specific values with float strings and large ints >= 2^53
        args = parser.parse_args(["--rate", "H1:2048.0", "L1:1e3", "--seed", "H1:9007199254740993"])
        self.assertEqual(args.rate["H1"], 2048)
        self.assertEqual(args.rate["L1"], 1000)
        self.assertEqual(args.seed["H1"], 9007199254740993)

        # Uniform global default
        args_global = parser.parse_args(["--rate", "2048.0"])
        self.assertEqual(args_global.rate["H1"], 2048)
        self.assertEqual(args_global.rate["V1"], 2048)

    def test_multidet_option_action_rejects_fractional(self):
        """Verify MultiDetOptionAction rejects fractional inputs for integer options."""
        parser = argparse.ArgumentParser()
        parser.add_argument("--rate", type=int, nargs="+", action=MultiDetOptionAction)
        with self.assertRaises(ValueError):
            parser.parse_args(["--rate", "H1:2048.5"])
        with self.assertRaises(ValueError):
            parser.parse_args(["--rate", "2048.5"])

    def test_other_actions_with_int(self):
        """Verify MultiColon, OptionAppend, DictOption, and MultiDetDictOption actions with type=int."""
        parser = argparse.ArgumentParser()
        parser.add_argument("--colon", type=int, nargs="+", action=MultiDetMultiColonOptionAction)
        parser.add_argument("--append", type=int, nargs="+", action=MultiDetOptionAppendAction)
        parser.add_argument("--dict", type=int, nargs="+", action=DictOptionAction)
        parser.add_argument("--mdict", type=int, nargs="+", action=MultiDetDictOptionAction)

        args = parser.parse_args([
            "--colon", "H1:2048.0", "L1:1e3",
            "--append", "H1:1e3", "H1:2048.0",
            "--dict", "rate:2048.0", "count:1e3",
            "--mdict", "H1:rate:2048.0",
        ])
        self.assertEqual(args.colon["H1"], 2048)
        self.assertEqual(args.append["H1"], [1000, 2048])
        self.assertEqual(args.dict["rate"], 2048)
        self.assertEqual(args.mdict["H1"]["rate"], 2048)

        # Rejection of fractionals
        with self.assertRaises(ValueError):
            parser.parse_args(["--colon", "H1:1.5"])
        with self.assertRaises(ValueError):
            parser.parse_args(["--append", "H1:1.5"])
        with self.assertRaises(ValueError):
            parser.parse_args(["--dict", "rate:1.5"])
        with self.assertRaises(ValueError):
            parser.parse_args(["--mdict", "H1:rate:1.5"])

    def test_positive_and_nonnegative_int(self):
        """Verify positive_int and nonnegative_int types accept float strings and reject bad values."""
        self.assertEqual(positive_int("2048.0"), 2048)
        self.assertEqual(positive_int("1e3"), 1000)
        self.assertEqual(positive_int("9007199254740993"), 9007199254740993)

        self.assertEqual(nonnegative_int("0.0"), 0)
        self.assertEqual(nonnegative_int("0"), 0)
        self.assertEqual(nonnegative_int("1e2"), 100)

        for bad in ["0", "-5", "2048.5", "abc"]:
            with self.assertRaises(argparse.ArgumentTypeError):
                positive_int(bad)
        for bad in ["-1", "1.5", "abc"]:
            with self.assertRaises(argparse.ArgumentTypeError):
                nonnegative_int(bad)


class TestConfigParserGetInt(unittest.TestCase):
    def test_getint(self):
        """Verify InterpolatingConfigParser.getint float strings, large ints, and fallback."""
        cp = InterpolatingConfigParser()
        cp.read_string("""
[workflow]
sample-rate = 2048.0
niterations = 1e6
large-seed = 9007199254740993
gps-nano = 1187000000123456789
negative-int = -5.0
fractional = 2048.5
""")
        self.assertEqual(cp.getint("workflow", "sample-rate"), 2048)
        self.assertEqual(cp.getint("workflow", "niterations"), 1000000)
        self.assertEqual(cp.getint("workflow", "large-seed"), 9007199254740993)
        self.assertEqual(cp.getint("workflow", "gps-nano"), 1187000000123456789)
        self.assertEqual(cp.getint("workflow", "negative-int"), -5)
        self.assertEqual(cp.getint("workflow", "missing", fallback=42), 42)

    def test_getint_fractional_rejected(self):
        """Verify getint raises ValueError on fractional string."""
        cp = InterpolatingConfigParser()
        cp.read_string("[workflow]\nfractional = 2048.5\n")
        with self.assertRaises(ValueError):
            cp.getint("workflow", "fractional")


if __name__ == "__main__":
    unittest.main()
