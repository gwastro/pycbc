# Copyright (C) 2016  Collin Capano
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


#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#
"""
This modules provides utilities for manipulating parameter boundaries. Namely,
classes are offered that will map values to a specified domain using either
cyclic boundaries or reflected boundaries.
"""

import numpy

class _Bound(float):
    """Adds methods to float for boundary comparisons."""

    name = None

    def larger(self, other):
        """A function to determine whether or not `other` is larger
        than the bound. This raises a NotImplementedError; classes that
        inherit from this must define it.
        """
        raise NotImplementedError("larger function not set")

    def smaller(self, other):
        """A function to determine whether or not `other` is smaller
        than the bound. This raises a NotImplementedError; classes that
        inherit from this must define it.
        """
        raise NotImplementedError("smaller function not set")


class OpenBound(_Bound):
    """Sets larger and smaller functions to be `>` and `<`, respectively."""

    name = 'open'

    def larger(self, other):
        """Returns True if `other` is `>`, False otherwise"""
        return self > other

    def smaller(self, other):
        """Returns True if `other` is `<`, False otherwise."""
        return self < other


class ClosedBound(_Bound):
    """Sets larger and smaller functions to be `>=` and `<=`, respectively."""

    name = 'closed'

    def larger(self, other):
        return self >= other

    def smaller(self, other):
        return self <= other


class ReflectedBound(ClosedBound):
    """Inherits from `ClosedBound`, adding reflection functions."""

    name = 'reflected'

    def reflect(self, value):
        return 2*self - value

    def reflect_left(self, value):
        """Only reflects the value if is > self."""
        if value > self:
            value = self.reflect(value)
        return value

    def reflect_right(self, value):
        """Only reflects the value if is < self."""
        if value < self:
            value = self.reflect(value)
        return value


boundary_types = {
    OpenBound.name: OpenBound,
    ClosedBound.name: ClosedBound,
    ReflectedBound.name: ReflectedBound
}


#
#   Helper functions for applying conditions to boundaries
#

def apply_cyclic(value, bounds):
    """Given a value, applies cyclic boundary conditions between the minimum
    and maximum bounds.

    Parameters
    ----------
    value : float
        The value to apply the cyclic conditions to.
    bounds : Bounds instance
        Boundaries to use for applying cyclic conditions.

    Returns
    -------
    float
        The value after the cyclic bounds are applied.
    """
    return (value - bounds._min) %(bounds._max - bounds._min) + bounds._min

def reflect_well(value, bounds):
    """Given some boundaries, reflects the value until it falls within both
    boundaries. This is done iteratively, reflecting left off of the
    `boundaries.max`, then right off of the `boundaries.min`, etc.

    Parameters
    ----------
    value : float
        The value to apply the reflected boundaries to.
    bounds : Bounds instance
        Boundaries to reflect between. Both `bounds.min` and `bounds.max` must
        be instances of `ReflectedBound`, otherwise an AttributeError is
        raised.

    Returns
    -------
    float
        The value after being reflected between the two bounds.
    """
    while value not in bounds:
        value = bounds._max.reflect_left(value)
        value = bounds._min.reflect_right(value)
    return value


def _pass(value):
    """Just return the given value."""
    return value


#
#   Bounds class
#

class Bounds(object):
    """Creates and stores bounds using the given values.

    The type of boundaries used can be set using the `btype_(min|max)`
    parameters. These arguments set what kind of boundary is used at the
    minimum and maximum bounds. Specifically, if `btype_min` (`btype_max`) is
    set to:

     * "open": the minimum (maximum) boundary will be an instance of
       `OpenBound`. This means that a value must be `>` (`<`) the bound
       for it to be considered within the bounds.
     * "closed": the minimum (maximum) boundary will be an instance of
       `ClosedBound`. This means that a value must be `>=` (`<=`) the bound
       for it to be considered within the bounds.
     * "reflected": the minimum (maximum) boundary will be an isntance of
       `ReflectedBound`. This means that a value will be reflected to the
       right (left) if `apply_conditions` is used on the value. For more
       details see `apply_conditions`.

    If the `cyclic` keyword is set to True, then `apply_conditions` will cause
    values to be wrapped around to the minimum (maximum) bound if the value
    is > (<=) the maximum (minimum) bound. For more details see
    `apply_conditions`.

    Values can be checked whether or not they occur within the bounds using
    `in`; e.g., `6 in bounds`. This is done without applying any boundary
    conditions. To apply conditions, then check whether the value is in
    bounds, use the `contains_conditioned` method.

    The default is for the minimum bound to be "closed" and the maximum bound
    to be "open", i.e., a right-open interval.

    Parameters
    ----------
    min_bound : {-numpy.inf, float}
        The value of the lower bound. Default is `-inf`.
    max_bound : {numpy.inf, float}
        The value of the upper bound. Default is `inf`.
    btype_min : {'closed', string}
        The type of the lower bound; options are "closed", "open", or
        "reflected". Default is "closed".
    btype_min : {'open', string}
        The type of the lower bound; options are "closed", "open", or
        "reflected". Default is "open".
    cyclic : {False, bool}
        Whether or not to make the bounds cyclic; default is False. If True,
        both the minimum and maximum bounds must be finite.

    Examples
    --------
    Create a right-open interval between -1 and 1 and test whether various
    values are within them:
    >>> bounds = Bounds(-1., 1.)
    >>> -1 in bounds
    True
    >>> 0 in bounds
    True
    >>> 1 in bounds
    False

    Create an open interval between -1 and 1 and test the same values:
    >>> bounds = Bounds(-1, 1, btype_min="open")
    >>> -1 in bounds
    False
    >>> 0 in bounds
    True
    >>> 1 in bounds
    False

    Create cyclic bounds between -1 and 1 and plot the effect of conditioning
    on points between -10 and 10:
    >>> bounds = Bounds(-1, 1, cyclic=True)
    >>> x = numpy.linspace(-10, 10, num=1000)
    >>> conditioned_x = bounds.apply_conditions(x)
    >>> fig = pyplot.figure()
    >>> ax = fig.add_subplot(111)
    >>> ax.plot(x, x, c='b', lw=2, label='input')
    >>> ax.plot(conditioned_x, x, c='r', lw=1)
    >>> ax.vlines([-1., 1.], x.min(), x.max(), color='k', linestyle='--')
    >>> ax.set_title('cyclic bounds between x=-1,1')
    >>> fig.show()

    Create a reflected bound at -1 and plot the effect of conditioning:
    >>> bounds = Bounds(-1, 1, btype_min='reflected')
    >>> x = numpy.linspace(-10, 10, num=1000)
    >>> conditioned_x = bounds.apply_conditions(x)
    >>> fig = pyplot.figure()
    >>> ax = fig.add_subplot(111)
    >>> ax.plot(x, x, c='b', lw=2, label='input')
    >>> ax.plot(conditioned_x, x, c='r', lw=1)
    >>> ax.vlines([-1., 1.], x.min(), x.max(), color='k', linestyle='--')
    >>> ax.set_title('reflected right at x=-1')
    >>> fig.show()

    Create a reflected bound at 1 and plot the effect of conditioning:
    >>> bounds = Bounds(-1, 1, btype_max='reflected')
    >>> x = numpy.linspace(-10, 10, num=1000)
    >>> conditioned_x = bounds.apply_conditions(x)
    >>> fig = pyplot.figure()
    >>> ax = fig.add_subplot(111)
    >>> ax.plot(x, x, c='b', lw=2, label='input')
    >>> ax.plot(conditioned_x, x, c='r', lw=1)
    >>> ax.vlines([-1., 1.], x.min(), x.max(), color='k', linestyle='--')
    >>> ax.set_title('reflected left at x=1')
    >>> fig.show()

    Create reflected bounds at -1 and 1 and plot the effect of conditioning:
    >>> bounds = Bounds(-1, 1, btype_min='reflected', btype_max='reflected')
    >>> x = numpy.linspace(-10, 10, num=1000)
    >>> conditioned_x = bounds.apply_conditions(x)
    >>> fig = pyplot.figure()
    >>> ax = fig.add_subplot(111)
    >>> ax.plot(x, x, c='b', lw=2, label='input')
    >>> ax.plot(conditioned_x, x, c='r', lw=1)
    >>> ax.vlines([-1., 1.], x.min(), x.max(), color='k', linestyle='--')
    >>> ax.set_title('reflected betewen x=-1,1')
    >>> fig.show()
    """

    def __init__(self, min_bound=-numpy.inf, max_bound=numpy.inf,
            btype_min='closed', btype_max='open', cyclic=False):
        # check boundary values
        if min_bound >= max_bound:
            raise ValueError("min_bound must be < max_bound")
        if cyclic and not (
                numpy.isfinite(min_bound) and numpy.isfinite(max_bound)):
            raise ValueError("if using cyclic, min and max bounds must both "
                "be finite")
        # store bounds
        try:
            self._min = boundary_types[btype_min](min_bound)
        except KeyError:
            raise ValueError("unrecognized btype_min {}".format(btype_min))
        try:
            self._max = boundary_types[btype_max](max_bound)
        except KeyError:
            raise ValueError("unrecognized btype_max {}".format(btype_max))
        # store cyclic conditions
        self._cyclic = bool(cyclic)
        # store reflection conditions; we'll vectorize them here so that they
        # can be used with arrays
        if self._min.name == 'reflected' and self._max.name == 'reflected':
            self._reflect = numpy.vectorize(self._reflect_well)
            self.reflected = 'well'
        elif self._min.name == 'reflected':
            self._reflect = numpy.vectorize(self._min.reflect_right)
            self.reflected = 'min'
        elif self._max.name == 'reflected':
            self._reflect = numpy.vectorize(self._max.reflect_left)
            self.reflected = 'max'
        else:
            self._reflect = _pass
            self.reflected = False

    def __repr__(self):
        return str(self.__class__)[:-1] + " " + " ".join(
                   map(str, ["min", self._min, "max", self._max,
                             "cyclic", self._cyclic])) + ">"

    @property
    def min(self):
        """_bounds instance: The minimum bound """
        return self._min

    @property
    def max(self):
        """_bounds instance: The maximum bound """
        return self._max

    @property
    def cyclic(self):
        """bool: Whether the bounds are cyclic or not.
        """
        return self._cyclic

    def __getitem__(self, ii):
        if ii == 0:
            return self._min
        elif ii == 1:
            return self._max
        else:
            raise IndexError("index {} out of range".format(ii))

    def __abs__(self):
        return abs(self._max - self._min)

    def __contains__(self, value):
        return self._min.smaller(value) & self._max.larger(value)

    def _reflect_well(self, value):
        """Thin wrapper around `reflect_well` that passes self as the `bounds`.
        """
        return reflect_well(value, self)

    def _apply_cyclic(self, value):
        """Thin wrapper around `apply_cyclic` that passes self as the `bounds`.
        """
        return apply_cyclic(value, self)

    def apply_conditions(self, value):
        """Applies any boundary conditions to the given value.

        The value is manipulated according based on the following conditions:

         * If `self.cyclic` is True then `value` is wrapped around to the
           minimum (maximum) bound if `value` is `>= self.max` (`< self.min`)
           bound. For example, if the minimum and maximum bounds are `0, 2*pi`
           and `value = 5*pi`, then the returned value will be `pi`.
         * If `self.min` is a reflected boundary then `value` will be
           reflected to the right if it is `< self.min`. For example, if
           `self.min = 10` and `value = 3`, then the returned value will be
           17.
         * If `self.max` is a reflected boundary then `value` will be
           reflected to the left if it is `> self.max`. For example, if
           `self.max = 20` and `value = 27`, then the returned value will be
           13.
         * If `self.min` and `self.max` are both reflected boundaries, then
           `value` will be reflected between the two boundaries until it
           falls within the bounds. The first reflection occurs off of the
           maximum boundary. For example, if `self.min = 10`, `self.max =
           20`, and `value = 42`, the returned value will be 18 ( the first
           reflection yields -2, the second 22, and the last 18).
         * If neither bounds are reflected and cyclic is False, then the
           value is just returned as-is.

        Parameters
        ----------
        value : float
            The value to apply the conditions to.

        Returns
        -------
        float
            The value after the conditions are applied; see above for details.
        """
        retval = value
        if self._cyclic:
            retval = apply_cyclic(value, self)
        retval = self._reflect(retval)
        if isinstance(retval, numpy.ndarray) and retval.size == 1:
            try:
                retval = retval[0]
            except IndexError:
                retval = float(retval)
        return retval

    def contains_conditioned(self, value):
        """Runs `apply_conditions` on the given value before testing whether it
        is in bounds. Note that if `cyclic` is True, or both bounds
        are reflected, than this will always return True.

        Parameters
        ----------
        value : float
            The value to test.

        Returns
        -------
        bool
            Whether or not the value is within the bounds after the boundary
            conditions are applied.
        """
        return self.apply_conditions(value) in self
