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
import unittest
from pycbc.types.optparse import (
    to_int,
    MultiDetOptionAction,
    MultiDetMultiColonOptionAction,
    MultiDetOptionAppendAction,
    DictOptionAction,
    positive_int,
    nonnegative_int,
)
from pycbc.types.config import InterpolatingConfigParser


class TestToInt(unittest.TestCase):
    def test_to_int_valid(self):
        """Verify float strings and exact arbitrary-precision integers."""
        # Exact large integers >= 2^53 (where float64 loses precision)
        self.assertEqual(to_int("9007199254740993"), 9007199254740993)
        self.assertEqual(to_int(1187000000123456789), 1187000000123456789)
        # Whole-number float strings and scientific notation
        for s, expected in [("2048.0", 2048), ("1e3", 1000), ("-5.0", -5), ("0.0", 0), (2048.0, 2048)]:
            res = to_int(s)
            self.assertEqual(res, expected)
            self.assertIsInstance(res, int)

    def test_to_int_rejected(self):
        """Verify fractional and invalid inputs raise ValueError."""
        for invalid in ["2048.5", "1.5", "-0.7", 2048.5, "abc", "", "inf", "nan", None]:
            with self.assertRaises(ValueError):
                to_int(invalid)


class TestOptparseActions(unittest.TestCase):
    def test_multidet_option_action(self):
        """Verify MultiDetOptionAction accepts float strings, large ints, and rejects fractionals."""
        parser = argparse.ArgumentParser()
        parser.add_argument("--rate", type=int, nargs="+", action=MultiDetOptionAction)
        parser.add_argument("--seed", type=int, nargs="+", action=MultiDetOptionAction)
        args = parser.parse_args(["--rate", "H1:2048.0", "L1:1e3", "--seed", "H1:9007199254740993"])
        self.assertEqual(args.rate["H1"], 2048)
        self.assertEqual(args.rate["L1"], 1000)
        self.assertEqual(args.seed["H1"], 9007199254740993)

        # Global single-value default
        args_global = parser.parse_args(["--rate", "2048.0"])
        self.assertEqual(args_global.rate["H1"], 2048)
        self.assertEqual(args_global.rate["V1"], 2048)

        # Fractional rejected
        with self.assertRaises(ValueError):
            parser.parse_args(["--rate", "H1:2048.5"])

    def test_other_actions_with_int(self):
        """Verify MultiColon, OptionAppend, and DictOption actions with type=int."""
        parser = argparse.ArgumentParser()
        parser.add_argument("--colon", type=int, nargs="+", action=MultiDetMultiColonOptionAction)
        parser.add_argument("--append", type=int, nargs="+", action=MultiDetOptionAppendAction)
        parser.add_argument("--dict", type=int, nargs="+", action=DictOptionAction)

        args = parser.parse_args([
            "--colon", "H1:2048.0", "L1:1e3",
            "--append", "H1:1e3", "H1:2048.0",
            "--dict", "rate:2048.0", "count:1e3",
        ])
        self.assertEqual(args.colon["H1"], 2048)
        self.assertEqual(args.append["H1"], [1000, 2048])
        self.assertEqual(args.dict["rate"], 2048)

        with self.assertRaises(ValueError):
            parser.parse_args(["--colon", "H1:1.5"])
        with self.assertRaises(ValueError):
            parser.parse_args(["--dict", "rate:1.5"])

    def test_positive_and_nonnegative_int(self):
        """Verify positive_int and nonnegative_int types."""
        self.assertEqual(positive_int("2048.0"), 2048)
        self.assertEqual(positive_int("9007199254740993"), 9007199254740993)
        self.assertEqual(nonnegative_int("0.0"), 0)
        self.assertEqual(nonnegative_int("1e2"), 100)

        for bad in ["0", "-5", "2048.5"]:
            with self.assertRaises(argparse.ArgumentTypeError):
                positive_int(bad)
        for bad in ["-1", "1.5"]:
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
fractional = 2048.5
""")
        self.assertEqual(cp.getint("workflow", "sample-rate"), 2048)
        self.assertEqual(cp.getint("workflow", "niterations"), 1000000)
        self.assertEqual(cp.getint("workflow", "large-seed"), 9007199254740993)
        self.assertEqual(cp.getint("workflow", "missing", fallback=42), 42)

        with self.assertRaises(ValueError):
            cp.getint("workflow", "fractional")


if __name__ == "__main__":
    unittest.main()
