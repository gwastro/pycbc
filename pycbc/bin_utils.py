from bisect import bisect_right
try:
    from fpconst import PosInf, NegInf
except ImportError:
    # fpconst is not part of the standard library and might not be available
    PosInf = float("+inf")
    NegInf = float("-inf")
import numpy
import math


class Bins(object):

    """
    Parent class for 1-dimensional binnings.

    Not intended to be used directly, but to be subclassed for use in real
    bins classes.
    """
    def __init__(self, minv, maxv, n):
        """
        Initialize a Bins instance.  The three arguments are the
        minimum and maximum of the values spanned by the bins, and
        the number of bins to place between them.  Subclasses may
        require additional arguments, or different arguments
        altogether.
        """
        # convenience code to do some common initialization and
        # input checking
        if not isinstance(n, int):
            raise TypeError(n)
        if n < 1:
            raise ValueError(n)
        if maxv <= minv:
            raise ValueError((minv, maxv))
        self.minv = minv
        self.maxv = maxv
        self.n = n

    def __len__(self):
        return self.n

    def __getitem__(self, x):
        """
        Convert a co-ordinate to a bin index.  The co-ordinate can
        be a single value, or a Python slice instance describing a
        range of values.  If a single value is given, it is mapped
        to the bin index corresponding to that value.  If a slice
        is given, it is converted to a slice whose lower bound is
        the index of the bin in which the slice's lower bound
        falls, and whose upper bound is 1 greater than the index of
        the bin in which the slice's upper bound falls.  Steps are
        not supported in slices.
        """
        if isinstance(x, slice):
            if x.step is not None:
                raise NotImplementedError("step not supported: %s" % repr(x))
            return slice(self[x.start] if x.start is not None
                         else 0, self[x.stop] + 1 if x.stop is not None
                         else len(self))
        raise NotImplementedError

    def __iter__(self):
        """
        If __iter__ does not exist, Python uses __getitem__ with
        range(0) as input to define iteration. This is nonsensical
        for bin objects, so explicitly unsupport iteration.
        """
        raise NotImplementedError

    def lower(self):
        """
        Return an array containing the locations of the lower
        boundaries of the bins.
        """
        raise NotImplementedError

    def centres(self):
        """
        Return an array containing the locations of the bin
        centres.
        """
        raise NotImplementedError

    def upper(self):
        """
        Return an array containing the locations of the upper
        boundaries of the bins.
        """
        raise NotImplementedError


class IrregularBins(Bins):

    """
    Bins with arbitrary, irregular spacing.  We only require strict
    monotonicity of the bin boundaries.  N boundaries define N-1 bins.

    Example:

    >>> x = IrregularBins([0.0, 11.0, 15.0, numpy.inf])
    >>> len(x)
    3
    >>> x[1]
    0
    >>> x[1.5]
    0
    >>> x[13]
    1
    >>> x[25]
    2
    >>> x[4:17]
    slice(0, 3, None)
    >>> IrregularBins([0.0, 15.0, 11.0])
    Traceback (most recent call last):
        ...
    ValueError: non-monotonic boundaries provided
    >>> y = IrregularBins([0.0, 11.0, 15.0, numpy.inf])
    >>> x == y
    True
    """
    def __init__(self, boundaries):
        """
        Initialize a set of custom bins with the bin boundaries.
        This includes all left edges plus the right edge.  The
        boundaries must be monotonic and there must be at least two
        elements.
        """
        # check pre-conditions
        if len(boundaries) < 2:
            raise ValueError("less than two boundaries provided")
        boundaries = tuple(boundaries)
        if any(a > b for a, b in zip(boundaries[:-1], boundaries[1:])):
            raise ValueError("non-monotonic boundaries provided")

        self.boundaries = boundaries
        self.n = len(boundaries) - 1
        self.minv = boundaries[0]
        self.maxv = boundaries[-1]

    def __getitem__(self, x):
        if isinstance(x, slice):
            return super(IrregularBins, self).__getitem__(x)
        if self.minv <= x < self.maxv:
            return bisect_right(self.boundaries, x) - 1
        # special measure-zero edge case
        if x == self.maxv:
            return len(self.boundaries) - 2
        raise IndexError(x)

    def lower(self):
        return numpy.array(self.boundaries[:-1])

    def upper(self):
        return numpy.array(self.boundaries[1:])

    def centres(self):
        return (self.lower() + self.upper()) / 2.0


class LinearBins(Bins):

    """
    Linearly-spaced bins.  There are n bins of equal size, the first
    bin starts on the lower bound and the last bin ends on the upper
    bound inclusively.

    Example:

    >>> x = LinearBins(1.0, 25.0, 3)
    >>> x.lower()
    array([  1.,   9.,  17.])
    >>> x.upper()
    array([  9.,  17.,  25.])
    >>> x.centres()
    array([  5.,  13.,  21.])
    >>> x[1]
    0
    >>> x[1.5]
    0
    >>> x[10]
    1
    >>> x[25]
    2
    >>> x[0:27]
    Traceback (most recent call last):
        ...
    IndexError: 0
    >>> x[1:25]
    slice(0, 3, None)
    >>> x[:25]
    slice(0, 3, None)
    >>> x[10:16.9]
    slice(1, 2, None)
    >>> x[10:17]
    slice(1, 3, None)
    >>> x[10:]
    slice(1, 3, None)
    """
    def __init__(self, minv, maxv, n):
        super(LinearBins, self).__init__(minv, maxv, n)
        self.delta = float(maxv - minv) / n

    def __getitem__(self, x):
        if isinstance(x, slice):
            return super(LinearBins, self).__getitem__(x)
        if self.minv <= x < self.maxv:
            return int(math.floor((x - self.minv) / self.delta))
        if x == self.maxv:
            # special "measure zero" corner case
            return len(self) - 1
        raise IndexError(x)

    def lower(self):
        return numpy.linspace(self.minv, self.maxv - self.delta, len(self))

    def centres(self):
        return numpy.linspace(self.minv + self.delta / 2.,
                              self.maxv - self.delta / 2., len(self))

    def upper(self):
        return numpy.linspace(self.minv + self.delta, self.maxv, len(self))


class LinearPlusOverflowBins(Bins):

    """
    Linearly-spaced bins with overflow at the edges.

    There are n-2 bins of equal size.  The bin 1 starts on the lower bound and
    bin n-2 ends on the upper bound.  Bins 0 and n-1 are overflow going from
    -infinity to the lower bound and from the upper bound to +infinity
    respectively.  Must have n >= 3.

    Example:

    >>> x = LinearPlusOverflowBins(1.0, 25.0, 5)
    >>> x.centres()
    array([-inf,   5.,  13.,  21.,  inf])
    >>> x.lower()
    array([-inf,   1.,   9.,  17.,  25.])
    >>> x.upper()
    array([  1.,   9.,  17.,  25.,  inf])
    >>> x[float("-inf")]
    0
    >>> x[0]
    0
    >>> x[1]
    1
    >>> x[10]
    2
    >>> x[24.99999999]
    3
    >>> x[25]
    4
    >>> x[100]
    4
    >>> x[float("+inf")]
    4
    >>> x[float("-inf"):9]
    slice(0, 3, None)
    >>> x[9:float("+inf")]
    slice(2, 5, None)
    """
    def __init__(self, minv, maxv, n):
        if n < 3:
            raise ValueError("n must be >= 3")
        super(LinearPlusOverflowBins, self).__init__(minv, maxv, n)
        self.delta = float(maxv - minv) / (n - 2)

    def __getitem__(self, x):
        if isinstance(x, slice):
            return super(LinearPlusOverflowBins, self).__getitem__(x)
        if self.minv <= x < self.maxv:
            return int(math.floor((x - self.minv) / self.delta)) + 1
        if x >= self.maxv:
            # +infinity overflow bin
            return len(self) - 1
        if x < self.minv:
            # -infinity overflow bin
            return 0
        raise IndexError(x)

    def lower(self):
        return numpy.concatenate(
            (numpy.array([NegInf]),
             self.minv + self.delta * numpy.arange(len(self) - 2),
             numpy.array([self.maxv]))
        )

    def centres(self):
        return numpy.concatenate(
            (numpy.array([NegInf]),
             self.minv + self.delta * (numpy.arange(len(self) - 2) + 0.5),
             numpy.array([PosInf]))
        )

    def upper(self):
        return numpy.concatenate(
            (numpy.array([self.minv]),
             self.minv + self.delta * (numpy.arange(len(self) - 2) + 1),
             numpy.array([PosInf]))
        )


class LogarithmicBins(Bins):

    """
    Logarithmically-spaced bins.

    There are n bins, each of whose upper and lower bounds differ by the same
    factor.  The first bin starts on the lower bound, and the last bin ends on
    the upper bound inclusively.

    Example:

    >>> x = LogarithmicBins(1.0, 25.0, 3)
    >>> x[1]
    0
    >>> x[5]
    1
    >>> x[25]
    2
    """
    def __init__(self, minv, maxv, n):
        super(LogarithmicBins, self).__init__(minv, maxv, n)
        self.delta = (math.log(maxv) - math.log(minv)) / n

    def __getitem__(self, x):
        if isinstance(x, slice):
            return super(LogarithmicBins, self).__getitem__(x)
        if self.minv <= x < self.maxv:
            return int(math.floor((math.log(x) - math.log(self.minv)) /
                                  self.delta))
        if x == self.maxv:
            # special "measure zero" corner case
            return len(self) - 1
        raise IndexError(x)

    def lower(self):
        return numpy.exp(
            numpy.linspace(math.log(self.minv), math.log(self.maxv) -
                           self.delta, len(self))
        )

    def centres(self):
        return numpy.exp(
            numpy.linspace(math.log(self.minv), math.log(self.maxv) -
                           self.delta, len(self)) + self.delta / 2.
        )

    def upper(self):
        return numpy.exp(
            numpy.linspace(math.log(self.minv) + self.delta,
                           math.log(self.maxv), len(self))
        )


class LogarithmicPlusOverflowBins(Bins):

    """
    Logarithmically-spaced bins plus one bin at each end that goes to
    zero and positive infinity respectively.  There are n-2 bins each
    of whose upper and lower bounds differ by the same factor.  Bin 1
    starts on the lower bound, and bin n-2 ends on the upper bound
    inclusively.  Bins 0 and n-1 are overflow bins extending from 0 to
    the lower bound and from the upper bound to +infinity respectively.
    Must have n >= 3.

    Example:

    >>> x = LogarithmicPlusOverflowBins(1.0, 25.0, 5)
    >>> x[0]
    0
    >>> x[1]
    1
    >>> x[5]
    2
    >>> x[24.999]
    3
    >>> x[25]
    4
    >>> x[100]
    4
    >>> x.lower()
    array([ 0.   ,  1.        ,  2.92401774,  8.54987973, 25.      ])
    >>> x.upper()
    array([ 1.   ,  2.92401774,  8.54987973, 25.        ,       inf])
    >>> x.centres()
    array([ 0.   ,  1.70997595,  5.        , 14.62008869,       inf])
    """
    def __init__(self, minv, maxv, n):
        if n < 3:
            raise ValueError("n must be >= 3")
        super(LogarithmicPlusOverflowBins, self).__init__(minv, maxv, n)
        self.delta = (math.log(maxv) - math.log(minv)) / (n - 2)

    def __getitem__(self, x):
        if isinstance(x, slice):
            return super(LogarithmicPlusOverflowBins, self).__getitem__(x)
        if self.minv <= x < self.maxv:
            return 1 + int(math.floor((math.log(x) - math.log(self.minv)) /
                                      self.delta))
        if x >= self.maxv:
            # infinity overflow bin
            return len(self) - 1
        if x < self.minv:
            # zero overflow bin
            return 0
        raise IndexError(x)

    def lower(self):
        return numpy.concatenate(
            (numpy.array([0.]),
             numpy.exp(numpy.linspace(math.log(self.minv), math.log(self.maxv),
                                      len(self) - 1))
             )
        )

    def centres(self):
        return numpy.concatenate(
            (numpy.array([0.]),
             numpy.exp(numpy.linspace(math.log(self.minv), math.log(self.maxv) -
                       self.delta, len(self) - 2) + self.delta / 2.),
             numpy.array([PosInf])
             )
        )

    def upper(self):
        return numpy.concatenate(
            (numpy.exp(numpy.linspace(math.log(self.minv), math.log(self.maxv),
                       len(self) - 1)),
             numpy.array([PosInf])
             )
        )


class NDBins(tuple):

    """
    Multi-dimensional co-ordinate binning.  An instance of this object
    is used to convert a tuple of co-ordinates into a tuple of bin
    indices.  This can be used to allow the contents of an array object
    to be accessed with real-valued coordinates.

    NDBins is a subclass of the tuple builtin, and is initialized with
    an iterable of instances of subclasses of Bins.  Each Bins subclass
    instance describes the binning to apply in the corresponding
    co-ordinate direction, and the number of them sets the dimensions
    of the binning.

    Example:

    >>> x = NDBins((LinearBins(1, 25, 3), LogarithmicBins(1, 25, 3)))
    >>> x[1, 1]
    (0, 0)
    >>> x[1.5, 1]
    (0, 0)
    >>> x[10, 1]
    (1, 0)
    >>> x[1, 5]
    (0, 1)
    >>> x[1, 1:5]
    (0, slice(0, 2, None))
    >>> x.centres()
    (array([ 5., 13., 21.]), array([ 1.70997595,  5.        , 14.62008869]))

    Note that the co-ordinates to be converted must be a tuple, even if
    it is only a 1-dimensional co-ordinate.
    """
    def __new__(cls, *args):
        new = tuple.__new__(cls, *args)
        new.minv = tuple(b.minv for b in new)
        new.maxv = tuple(b.maxv for b in new)
        new.shape = tuple(len(b) for b in new)
        return new

    def __getitem__(self, coords):
        """
        When coords is a tuple, it is interpreted as an
        N-dimensional co-ordinate which is converted to an N-tuple
        of bin indices by the Bins instances in this object.
        Otherwise coords is interpeted as an index into the tuple,
        and the corresponding Bins instance is returned.

        Example:

        >>> x = NDBins((LinearBins(1, 25, 3), LogarithmicBins(1, 25, 3)))
        >>> x[1, 1]
        (0, 0)
        >>> type(x[1])
        <class 'pylal.rate.LogarithmicBins'>

        When used to convert co-ordinates to bin indices, each
        co-ordinate can be anything the corresponding Bins instance
        will accept.  Note that the co-ordinates to be converted
        must be a tuple, even if it is only a 1-dimensional
        co-ordinate.
        """
        if isinstance(coords, tuple):
            if len(coords) != len(self):
                raise ValueError("dimension mismatch")
            return tuple(map(lambda b, c: b[c], self, coords))
        else:
            return tuple.__getitem__(self, coords)

    def lower(self):
        """
        Return a tuple of arrays, where each array contains the
        locations of the lower boundaries of the bins in the
        corresponding dimension.
        """
        return tuple(b.lower() for b in self)

    def centres(self):
        """
        Return a tuple of arrays, where each array contains the
        locations of the bin centres for the corresponding
        dimension.
        """
        return tuple(b.centres() for b in self)

    def upper(self):
        """
        Return a tuple of arrays, where each array contains the
        locations of the upper boundaries of the bins in the
        corresponding dimension.
        """
        return tuple(b.upper() for b in self)


class BinnedArray(object):

    """
    A convenience wrapper, using the NDBins class to provide access to
    the elements of an array object.  Technical reasons preclude
    providing a subclass of the array object, so the array data is made
    available as the "array" attribute of this class.

    Examples:

    Note that even for 1 dimensional arrays the index must be a tuple.

    >>> x = BinnedArray(NDBins((LinearBins(0, 10, 5),)))
    >>> x.array
    array([ 0.,  0.,  0.,  0.,  0.])
    >>> x[0,] += 1
    >>> x[0.5,] += 1
    >>> x.array
    array([ 2.,  0.,  0.,  0.,  0.])
    >>> x.argmax()
    (1.0,)

    Note the relationship between the binning limits, the bin centres,
    and the co-ordinates of the BinnedArray

    >>> x = BinnedArray(NDBins((LinearBins(-0.5, 1.5, 2), \
    LinearBins(-0.5, 1.5, 2))))
    >>> x.bins.centres()
    (array([ 0.,  1.]), array([ 0.,  1.]))
    >>> x[0, 0] = 0
    >>> x[0, 1] = 1
    >>> x[1, 0] = 2
    >>> x[1, 1] = 4
    >>> x.array
    array([[ 0.,  1.],
           [ 2.,  4.]])
    >>> x[0, 0]
    0.0
    >>> x[0, 1]
    1.0
    >>> x[1, 0]
    2.0
    >>> x[1, 1]
    4.0
    >>> x.argmin()
    (0.0, 0.0)
    >>> x.argmax()
    (1.0, 1.0)
    """
    def __init__(self, bins, array=None, dtype="double"):
        self.bins = bins
        if array is None:
            self.array = numpy.zeros(bins.shape, dtype=dtype)
        else:
            if array.shape != bins.shape:
                raise ValueError("input array and input bins must have the "
                                 "same shape")
            self.array = array

    def __getitem__(self, coords):
        return self.array[self.bins[coords]]

    def __setitem__(self, coords, val):
        self.array[self.bins[coords]] = val

    def __len__(self):
        return len(self.array)

    def copy(self):
        """
        Return a copy of the BinnedArray.  The .bins attribute is
        shared with the original.
        """
        return type(self)(self.bins, self.array.copy())

    def centres(self):
        """
        Return a tuple of arrays containing the bin centres for
        each dimension.
        """
        return self.bins.centres()

    def argmin(self):
        """
        Return the co-ordinates of the bin centre containing the
        minimum value.  Same as numpy.argmin(), converting the
        indexes to bin co-ordinates.
        """
        return tuple(centres[index] for centres, index in
                     zip(self.centres(), numpy.unravel_index(self.array.argmin(),
                                                             self.array.shape)))

    def argmax(self):
        """
        Return the co-ordinates of the bin centre containing the
        maximum value.  Same as numpy.argmax(), converting the
        indexes to bin co-ordinates.
        """
        return tuple(centres[index] for centres, index in
                     zip(self.centres(), numpy.unravel_index(self.array.argmax(),
                                                             self.array.shape)))

    def logregularize(self, epsilon=2**-1074):
        """
        Find bins <= 0, and set them to epsilon, This has the
        effect of allowing the logarithm of the array to be
        evaluated without error.
        """
        self.array[self.array <= 0] = epsilon
        return self


class BinnedRatios(object):

    """
    Like BinnedArray, but provides a numerator array and a denominator
    array.  The incnumerator() method increments a bin in the numerator
    by the given weight, and the incdenominator() method increments a
    bin in the denominator by the given weight.  There are no methods
    provided for setting or decrementing either, but the they are
    accessible as the numerator and denominator attributes, which are
    both BinnedArray objects.
    """
    def __init__(self, bins, dtype="double"):
        self.numerator = BinnedArray(bins, dtype=dtype)
        self.denominator = BinnedArray(bins, dtype=dtype)

    def __getitem__(self, coords):
        return self.numerator[coords] / self.denominator[coords]

    def bins(self):
        return self.numerator.bins

    def incnumerator(self, coords, weight=1):
        """
        Add weight to the numerator bin at coords.
        """
        self.numerator[coords] += weight

    def incdenominator(self, coords, weight=1):
        """
        Add weight to the denominator bin at coords.
        """
        self.denominator[coords] += weight

    def ratio(self):
        """
        Compute and return the array of ratios.
        """
        return self.numerator.array / self.denominator.array

    def regularize(self):
        """
        Find bins in the denominator that are 0, and set them to 1.
        Presumably the corresponding bin in the numerator is also
        0, so this has the effect of allowing the ratio array to be
        evaluated without error, returning zeros in those bins that
        have had no weight added to them.
        """
        self.denominator.array[self.denominator.array == 0] = 1
        return self

    def logregularize(self, epsilon=2**-1074):
        """
        Find bins in the denominator that are 0, and set them to 1,
        while setting the corresponding bin in the numerator to
        float epsilon.  This has the effect of allowing the
        logarithm of the ratio array to be evaluated without error.
        """
        self.numerator.array[self.denominator.array == 0] = epsilon
        self.denominator.array[self.denominator.array == 0] = 1
        return self

    def centres(self):
        """
        Return a tuple of arrays containing the bin centres for
        each dimension.
        """
        return self.numerator.bins.centres()
