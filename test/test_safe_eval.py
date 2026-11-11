"""Tests for the restricted expression evaluator in pycbc.io.safe_eval."""

import math
import os
import unittest

import numpy as np

from pycbc.io.safe_eval import (
    AttributeNamespace,
    UnsafeExpressionError,
    parse_expression,
    safe_eval,
)

# Expressions that would execute code or read files with a plain eval
ESCAPE_PAYLOADS = [
    "__import__('os').system('true')",
    "().__class__.__base__.__subclasses__()",
    "(lambda: 1)()",
    "[x for x in ()]",
    "self._values",
    "f.__globals__",
    "numpy.load('x.npy')",
    "numpy.lib.format",
    "numpy.ctypeslib",
    "arr.tofile('x')",
    "(x := 1)",
    "func(**kw)",
    "open('/etc/passwd')",
]


class TestSafeEval(unittest.TestCase):
    def test_arithmetic_and_numpy(self):
        x = np.array([1.0, 4.0, 9.0])
        np.testing.assert_allclose(
            safe_eval("numpy.sqrt(x) * 2 + 1", {"x": x}), [3.0, 5.0, 7.0]
        )
        self.assertAlmostEqual(safe_eval("math.log10(100) ** 2"), 4.0)
        self.assertEqual(safe_eval("max(a, b)", {"a": 1, "b": 3, "max": max}), 3)

    def test_comparisons_and_logic(self):
        x = np.arange(5)
        np.testing.assert_array_equal(
            safe_eval("(x > 1) & (x < 4)", {"x": x}), [False, False, True, True, False]
        )
        self.assertTrue(safe_eval("1 < a <= 3", {"a": 3}))
        self.assertFalse(safe_eval("1 < a <= 3", {"a": 4}))
        self.assertTrue(safe_eval("a and not b", {"a": True, "b": False}))
        self.assertEqual(safe_eval("a or b", {"a": 0, "b": 7}), 7)

    def test_conditional_and_subscript(self):
        self.assertEqual(
            safe_eval("'A' if m > 4 else 'B'", {"m": 5}),
            "A",
        )
        self.assertEqual(safe_eval("x[1:3]", {"x": [0, 1, 2, 3]}), [1, 2])

    def test_attribute_namespace(self):
        ns = AttributeNamespace({"snr": np.array([5.0, 7.0])})
        np.testing.assert_array_equal(
            safe_eval("self.snr > 6", {"self": ns}), [False, True]
        )
        with self.assertRaises(AttributeError):
            safe_eval("self.chisq", {"self": ns})

    def test_undefined_name(self):
        with self.assertRaises(NameError):
            safe_eval("missing + 1")

    def test_escape_payloads_rejected(self):
        names = {
            "self": AttributeNamespace({}),
            "f": lambda: None,
            "arr": np.zeros(1),
            "func": dict,
            "kw": {},
        }
        for payload in ESCAPE_PAYLOADS:
            # builtins are absent, so unknown names fail with NameError
            with (
                self.subTest(payload=payload),
                self.assertRaises((UnsafeExpressionError, NameError)),
            ):
                safe_eval(payload, names)

    def test_modules_in_namespace_not_traversable(self):
        with self.assertRaises(UnsafeExpressionError):
            safe_eval("os.getcwd()", {"os": os})
        with self.assertRaises(UnsafeExpressionError):
            safe_eval("g('x')", {"g": eval})

    def test_parse_is_cached_and_reusable(self):
        tree = parse_expression("a + 1")
        self.assertIs(tree, parse_expression("a + 1"))
        self.assertEqual(safe_eval(tree, {"a": 1}), 2)
        self.assertEqual(safe_eval(tree, {"a": math.pi}), math.pi + 1)


class TestCallers(unittest.TestCase):
    def test_burn_in_tests(self):
        from pycbc.inference.burn_in import evaluate_tests

        is_burned, _ = evaluate_tests(
            "max_posterior & nacl",
            {"max_posterior": True, "nacl": True},
            {"max_posterior": 10, "nacl": 20},
        )
        self.assertTrue(is_burned)

    def test_distribution_from_string(self):
        from pycbc.distributions import UniformSky, distribution_from_string

        self.assertIsInstance(distribution_from_string("UniformSky()"), UniformSky)
        with self.assertRaises(UnsafeExpressionError):
            distribution_from_string("UniformSky.__init__.__globals__")

    def test_render_function_lookup(self):
        from pycbc.results.render import get_render_function, render_text

        self.assertIs(get_render_function("render_text"), render_text)
        for name in ("os.system", "get_embedded_config", "render_missing"):
            with self.assertRaises(ValueError):
                get_render_function(name)

    def test_fieldarray_expression(self):
        from pycbc.io.record import FieldArray

        arr = FieldArray.from_kwargs(
            mass1=np.array([1.0, 2.0]), mass2=np.array([3.0, 4.0])
        )
        np.testing.assert_allclose(arr["mass1 + mass2"], [4.0, 6.0])
        np.testing.assert_allclose(arr["log(mass2)"], np.log([3.0, 4.0]))
        with self.assertRaises(UnsafeExpressionError):
            arr["mass1 + ().__class__.__base__"]


if __name__ == "__main__":
    unittest.main()
