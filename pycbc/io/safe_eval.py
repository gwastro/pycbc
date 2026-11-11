# Copyright (C) 2026 PyCBC developers
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
"""Restricted evaluation of user-supplied arithmetic expressions.

Expressions such as ``'self.snr > 6'``, ``'mass1 + mass2'`` or
``'numpy.log(x) * 2'`` are supplied on command lines and in configuration
files. Python's :py:func:`eval` would execute arbitrary code from these
strings, and passing ``{"__builtins__": None}`` does not prevent that (an
expression can walk ``().__class__.__base__.__subclasses__()`` to reach
``os``).

This module parses the expression with :py:mod:`ast` and interprets the
tree itself, allowing only arithmetic, comparisons, boolean logic,
conditional expressions, subscripts, attribute access and function calls.
Names and attributes starting with an underscore are rejected, modules other
than numpy / math cannot be traversed, and file-access and introspection
callables cannot be reached.
"""

import ast
import builtins
import cmath
import functools
import logging
import math
import operator
import types

import numpy as np

logger = logging.getLogger("pycbc.io.safe_eval")

__all__ = [
    "AttributeNamespace",
    "UnsafeExpressionError",
    "parse_expression",
    "safe_eval",
]


class UnsafeExpressionError(ValueError):
    """Raised when an expression uses a construct that is not permitted."""


_BINARY_OPS = {
    ast.Add: operator.add,
    ast.Sub: operator.sub,
    ast.Mult: operator.mul,
    ast.Div: operator.truediv,
    ast.FloorDiv: operator.floordiv,
    ast.Mod: operator.mod,
    ast.Pow: operator.pow,
    ast.MatMult: operator.matmul,
    ast.BitAnd: operator.and_,
    ast.BitOr: operator.or_,
    ast.BitXor: operator.xor,
    ast.LShift: operator.lshift,
    ast.RShift: operator.rshift,
}

_UNARY_OPS = {
    ast.UAdd: operator.pos,
    ast.USub: operator.neg,
    ast.Not: operator.not_,
    ast.Invert: operator.invert,
}

_COMPARE_OPS = {
    ast.Eq: operator.eq,
    ast.NotEq: operator.ne,
    ast.Lt: operator.lt,
    ast.LtE: operator.le,
    ast.Gt: operator.gt,
    ast.GtE: operator.ge,
    ast.Is: operator.is_,
    ast.IsNot: operator.is_not,
    ast.In: lambda a, b: a in b,
    ast.NotIn: lambda a, b: a not in b,
}

# Every AST node type that may appear in an accepted expression
_ALLOWED_NODES = (
    ast.Expression,
    ast.Constant,
    ast.Name,
    ast.Load,
    ast.Attribute,
    ast.Call,
    ast.keyword,
    ast.BinOp,
    ast.UnaryOp,
    ast.BoolOp,
    ast.And,
    ast.Or,
    ast.Compare,
    ast.IfExp,
    ast.Subscript,
    ast.Slice,
    ast.Tuple,
    ast.List,
    *_BINARY_OPS,
    *_UNARY_OPS,
    *_COMPARE_OPS,
)

# Modules whose attributes may be looked up inside an expression
_ALLOWED_MODULES = frozenset(
    {"math", "cmath", "numpy", "numpy.linalg", "numpy.fft", "numpy.random"}
)

# Attribute names that give file-system, pickle or low-level access
_FORBIDDEN_ATTRIBUTES = frozenset(
    {
        "ctypes",
        "ctypeslib",
        "DataSource",
        "dump",
        "dumps",
        "f2py",
        "format",
        "fromfile",
        "fromregex",
        "genfromtxt",
        "lib",
        "load",
        "loadtxt",
        "memmap",
        "save",
        "savetxt",
        "savez",
        "savez_compressed",
        "testing",
        "tofile",
    }
)

# Builtins that give code execution or introspection if passed in a namespace
_FORBIDDEN_CALLABLES = frozenset(
    id(getattr(builtins, _name))
    for _name in (
        "eval",
        "exec",
        "compile",
        "open",
        "getattr",
        "setattr",
        "delattr",
        "globals",
        "locals",
        "vars",
        "dir",
        "type",
        "input",
        "breakpoint",
        "__import__",
    )
)


def _check_identifier(name, kind):
    if name.startswith("_"):
        raise UnsafeExpressionError(
            f"{kind} '{name}' is not allowed: names starting with '_' "
            "are not permitted in expressions"
        )
    if kind == "attribute" and name in _FORBIDDEN_ATTRIBUTES:
        raise UnsafeExpressionError(f"attribute '{name}' is not permitted")


def _check_value(value, label):
    """Reject values that would let an expression escape the sandbox."""
    if isinstance(value, types.ModuleType):
        if value.__name__ not in _ALLOWED_MODULES:
            raise UnsafeExpressionError(
                f"access to module '{value.__name__}' via '{label}' is not permitted"
            )
    elif id(value) in _FORBIDDEN_CALLABLES:
        raise UnsafeExpressionError(f"'{label}' is not permitted")
    return value


@functools.lru_cache(maxsize=1024)
def parse_expression(expression):
    """Parse ``expression`` and check that it only uses permitted syntax.

    Parameters
    ----------
    expression : str
        The expression to parse.

    Returns
    -------
    ast.Expression
        The validated syntax tree, which can be passed to :py:func:`safe_eval`.

    Raises
    ------
    SyntaxError
        If ``expression`` is not a valid Python expression.
    UnsafeExpressionError
        If ``expression`` uses a construct that is not permitted.
    """
    tree = ast.parse(str(expression).strip(), mode="eval")
    for node in ast.walk(tree):
        if not isinstance(node, _ALLOWED_NODES):
            raise UnsafeExpressionError(
                f"'{type(node).__name__}' is not allowed in expression {expression!r}"
            )
        if isinstance(node, ast.Name):
            _check_identifier(node.id, "name")
        elif isinstance(node, ast.Attribute):
            _check_identifier(node.attr, "attribute")
        elif isinstance(node, ast.keyword) and node.arg is None:
            raise UnsafeExpressionError("'**' keyword unpacking is not allowed")
    logger.debug("Parsed expression %r", expression)
    return tree


class _Interpreter:
    """Evaluate a validated expression tree against a namespace."""

    def __init__(self, names):
        self.names = names

    def visit(self, node):
        handler = getattr(self, "visit_" + type(node).__name__.lower(), None)
        if handler is None:
            raise UnsafeExpressionError(f"'{type(node).__name__}' is not allowed")
        return handler(node)

    def visit_expression(self, node):
        return self.visit(node.body)

    def visit_constant(self, node):
        return node.value

    def visit_name(self, node):
        try:
            value = self.names[node.id]
        except KeyError:
            raise NameError(f"name '{node.id}' is not defined") from None
        return _check_value(value, node.id)

    def visit_attribute(self, node):
        _check_identifier(node.attr, "attribute")
        base = self.visit(node.value)
        return _check_value(getattr(base, node.attr), node.attr)

    def visit_call(self, node):
        func = self.visit(node.func)
        args = [self.visit(arg) for arg in node.args]
        kwargs = {kw.arg: self.visit(kw.value) for kw in node.keywords}
        return func(*args, **kwargs)

    def visit_binop(self, node):
        return _BINARY_OPS[type(node.op)](self.visit(node.left), self.visit(node.right))

    def visit_unaryop(self, node):
        return _UNARY_OPS[type(node.op)](self.visit(node.operand))

    def visit_boolop(self, node):
        # Same short-circuit semantics as Python's ``and`` / ``or``
        is_and = isinstance(node.op, ast.And)
        result = self.visit(node.values[0])
        for value in node.values[1:]:
            if bool(result) != is_and:
                return result
            result = self.visit(value)
        return result

    def visit_compare(self, node):
        # Chained comparisons stop at the first false result, as in Python
        left = self.visit(node.left)
        outcome = True
        for op, comparator in zip(node.ops, node.comparators, strict=True):
            right = self.visit(comparator)
            outcome = _COMPARE_OPS[type(op)](left, right)
            if len(node.ops) > 1 and not outcome:
                return outcome
            left = right
        return outcome

    def visit_ifexp(self, node):
        branch = node.body if self.visit(node.test) else node.orelse
        return self.visit(branch)

    def visit_subscript(self, node):
        return self.visit(node.value)[self.visit(node.slice)]

    def visit_slice(self, node):
        return slice(
            None if node.lower is None else self.visit(node.lower),
            None if node.upper is None else self.visit(node.upper),
            None if node.step is None else self.visit(node.step),
        )

    def visit_tuple(self, node):
        return tuple(self.visit(elt) for elt in node.elts)

    def visit_list(self, node):
        return [self.visit(elt) for elt in node.elts]


def safe_eval(expression, names=None):
    """Evaluate a restricted arithmetic expression.

    Parameters
    ----------
    expression : str or ast.Expression
        The expression, or a tree returned by :py:func:`parse_expression`.
    names : dict, optional
        Mapping of the names the expression may use to their values.
        ``numpy``, ``np`` and ``math`` are always available unless overridden.

    Returns
    -------
    object
        The value of the expression.
    """
    tree = expression
    if not isinstance(expression, ast.Expression):
        tree = parse_expression(expression)
    namespace = {"numpy": np, "np": np, "math": math, "cmath": cmath}
    if names:
        namespace.update(names)
    return _Interpreter(namespace).visit(tree)


class AttributeNamespace:
    """Expose the entries of a dictionary as read-only attributes.

    This allows expressions written as ``'self.snr > 6'`` to be evaluated
    against a dictionary of arrays without giving access to the real object.
    """

    def __init__(self, values):
        object.__setattr__(self, "_values", dict(values))

    def __getattr__(self, name):
        try:
            return self._values[name]
        except KeyError:
            raise AttributeError(f"no column named '{name}'") from None

    def __setattr__(self, name, value):
        raise AttributeError("AttributeNamespace is read-only")
