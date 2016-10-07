import numpy


class Bins(object):
    """
    Parent class for 1-dimensional binnings. 

    Not intended to be used directly, but to be subclassed for use in real
    bins classes.
    """
    def __init__(self, min, max, n):
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
        if max <= min:
            raise ValueError((min, max))
        self.min = min
        self.max = max
        self.n = n

    def __len__(self):
        return self.n

    def __cmp__(self, other):
        """
        Two binnings are the same if they are instances of the same
        class, have the same lower and upper bounds, and the same
        count of bins.
        """
        if not isinstance(other, type(self)):
            return -1
        return cmp((type(self), self.min, self.max, len(self)), (type(other), other.min, other.max, len(other)))

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
            return slice(self[x.start] if x.start is not None else 0, self[x.stop] + 1 if x.stop is not None else len(self))
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
        self.min = boundaries[0]
        self.max = boundaries[-1]

    def __cmp__(self, other):
        """
        Two binnings are the same if they are instances of the same
        class, and have the same boundaries.
        """
        if not isinstance(other, type(self)):
            return -1
        return cmp(self.boundaries, other.boundaries)

    def __getitem__(self, x):
        if isinstance(x, slice):
            return super(IrregularBins, self).__getitem__(x)
        if self.min <= x < self.max:
            return bisect_right(self.boundaries, x) - 1
        # special measure-zero edge case
        if x == self.max:
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
    def __init__(self, min, max, n):
        super(LinearBins, self).__init__(min, max, n)
        self.delta = float(max - min) / n

    def __getitem__(self, x):
        if isinstance(x, slice):
            return super(LinearBins, self).__getitem__(x)
        if self.min <= x < self.max:
            return int(math.floor((x - self.min) / self.delta))
        if x == self.max:
            # special "measure zero" corner case
            return len(self) - 1
        raise IndexError(x)

    def lower(self):
        return numpy.linspace(self.min, self.max - self.delta, len(self))

    def centres(self):
        return numpy.linspace(self.min + self.delta / 2., self.max - self.delta / 2., len(self))

    def upper(self):
        return numpy.linspace(self.min + self.delta, self.max, len(self))


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
    def __init__(self, min, max, n):
        if n < 3:
            raise ValueError("n must be >= 3")
        super(LinearPlusOverflowBins, self).__init__(min, max, n)
        self.delta = float(max - min) / (n - 2)

    def __getitem__(self, x):
        if isinstance(x, slice):
            return super(LinearPlusOverflowBins, self).__getitem__(x)
        if self.min <= x < self.max:
            return int(math.floor((x - self.min) / self.delta)) + 1
        if x >= self.max:
            # +infinity overflow bin
            return len(self) - 1
        if x < self.min:
            # -infinity overflow bin
            return 0
        raise IndexError(x)

    def lower(self):
        return numpy.concatenate((numpy.array([NegInf]), self.min + self.delta * numpy.arange(len(self) - 2), numpy.array([self.max])))

    def centres(self):
        return numpy.concatenate((numpy.array([NegInf]), self.min + self.delta * (numpy.arange(len(self) - 2) + 0.5), numpy.array([PosInf])))

    def upper(self):
        return numpy.concatenate((numpy.array([self.min]), self.min + self.delta * (numpy.arange(len(self) - 2) + 1), numpy.array([PosInf])))


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
    def __init__(self, min, max, n):
        super(LogarithmicBins, self).__init__(min, max, n)
        self.delta = (math.log(max) - math.log(min)) / n

    def __getitem__(self, x):
        if isinstance(x, slice):
            return super(LogarithmicBins, self).__getitem__(x)
        if self.min <= x < self.max:
            return int(math.floor((math.log(x) - math.log(self.min)) / self.delta))
        if x == self.max:
            # special "measure zero" corner case
            return len(self) - 1
        raise IndexError(x)

    def lower(self):
        return numpy.exp(numpy.linspace(math.log(self.min), math.log(self.max) - self.delta, len(self)))

    def centres(self):
        return numpy.exp(numpy.linspace(math.log(self.min), math.log(self.max) - self.delta, len(self)) + self.delta / 2.)

    def upper(self):
        return numpy.exp(numpy.linspace(math.log(self.min) + self.delta, math.log(self.max), len(self)))


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
    array([  0.        ,   1.        ,   2.92401774,   8.54987973,  25.        ])
    >>> x.upper()
    array([  1.        ,   2.92401774,   8.54987973,  25.        ,          inf])
    >>> x.centres()
    array([  0.        ,   1.70997595,   5.        ,  14.62008869,          inf])
    """
    def __init__(self, min, max, n):
        if n < 3:
            raise ValueError("n must be >= 3")
        super(LogarithmicPlusOverflowBins, self).__init__(min, max, n)
        self.delta = (math.log(max) - math.log(min)) / (n - 2)

    def __getitem__(self, x):
        if isinstance(x, slice):
            return super(LogarithmicPlusOverflowBins, self).__getitem__(x)
        if self.min <= x < self.max:
            return 1 + int(math.floor((math.log(x) - math.log(self.min)) / self.delta))
        if x >= self.max:
            # infinity overflow bin
            return len(self) - 1
        if x < self.min:
            # zero overflow bin
            return 0
        raise IndexError(x)

    def lower(self):
        return numpy.concatenate((numpy.array([0.]), numpy.exp(numpy.linspace(math.log(self.min), math.log(self.max), len(self) - 1))))

    def centres(self):
        return numpy.concatenate((numpy.array([0.]), numpy.exp(numpy.linspace(math.log(self.min), math.log(self.max) - self.delta, len(self) - 2) + self.delta / 2.), numpy.array([PosInf])))

    def upper(self):
        return numpy.concatenate((numpy.exp(numpy.linspace(math.log(self.min), math.log(self.max), len(self) - 1)), numpy.array([PosInf])))


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
	(array([  5.,  13.,  21.]), array([  1.70997595,   5.        ,  14.62008869]))

	Note that the co-ordinates to be converted must be a tuple, even if
	it is only a 1-dimensional co-ordinate.
	"""
	def __new__(cls, *args):
		new = tuple.__new__(cls, *args)
		new.min = tuple(b.min for b in new)
		new.max = tuple(b.max for b in new)
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
