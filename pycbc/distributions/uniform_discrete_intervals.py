import numpy
from pycbc.distributions import bounded

class UniformIntervals(bounded.BoundedDist):
    name = 'uniform_intervals'

    def __init__(self, **params):
         # temporarily suppress numpy divide by 0 warning
        numpy.seterr(divide='ignore')
        self._bounds = {}
        self._stride = {}
        self._interval = {}

        stride_args = [p for p in params if p.endswith('_stride')]
        interval_args = [p for p in params if p.endswith('_interval')]

        self._stride = dict([[p.split("_stride")[0], params.pop(p)] for p in stride_args])
        self._interval = dict([[p.split("_interval")[0], params.pop(p)] for p in interval_args])

        super(UniformIntervals, self).__init__(**params)

        self._lognorm = -sum([numpy.log(abs(bnd[1]-bnd[0]))
                             for bnd in self._bounds.values()])
        self._norm = numpy.exp(self._lognorm)

        missing = set(self._stride.keys()) - set(params.keys())
        print self._stride.keys()
        print set(params.keys())

        if any(missing):
            raise ValueError("stride provided for unknow params {}".format(
                             ', '.join(missing)))
        missing = set(self._interval.keys()) - set(params.keys())
        if any(missing):
            raise ValueError("interval provided for unknow params {}".format(
                             ', '.join(missing)))
        # set default mean/var for params not specified
        self._stride.update(dict([[p, 0.]
            for p in params if p not in self._stride]))
        self._interval.update(dict([[p, 1.]
            for p in params if p not in self._interval]))

        numpy.seterr(divide='warn')

    @property
    def norm(self):
        return self._norm

    @property
    def lognorm(self):
        return self._lognorm

    @property
    def stride(self):
        return self._stride

    @property
    def interval(self):
        return self._interval

    def _pdf(self, **kwargs):
        """Returns the pdf at the given values. The keyword arguments must
        contain all of parameters in self's params. Unrecognized arguments are
        ignored.
        """
        if kwargs in self:
            return self._norm
        else:
            return 0.

    def _logpdf(self, **kwargs):
        """Returns the log of the pdf at the given values. The keyword
        arguments must contain all of parameters in self's params. Unrecognized
        arguments are ignored.
        """
        if kwargs in self:
            return self._lognorm
        else:
            return -numpy.inf


    def rvs(self, size=1, param=None):
        """Gives a set of random values drawn from this distribution.
        Parameters
        ----------
        size : {1, int}
            The number of values to generate; default is 1.
        param : {None, string}
            If provided, will just return values for the given parameter.
            Otherwise, returns random values for each parameter.
        Returns
        -------
        structured array
            The random values in a numpy structured array. If a param was
            specified, the array will only have an element corresponding to the
            given parameter. Otherwise, the array will have an element for each
            parameter in self's params.
        """
        if param is not None:
            dtype = [(param, float)]
        else:
            dtype = [(p, float) for p in self.params]
        arr = numpy.zeros(size, dtype=dtype)
        for (p,_) in dtype:
            dx = self._interval[p] + self._stride[p]
            x = numpy.arange(self._bounds[p][0], self._bounds[p][1],
                             dx)
            arr[p] = numpy.random.uniform(x, x + self._interval[p])

        return arr

    @classmethod
    def from_config(cls, cp, section, variable_args):
        """Returns a distribution based on a configuration file. The parameters
        for the distribution are retrieved from the section titled
        "[`section`-`variable_args`]" in the config file.
        Parameters
        ----------
        cp : pycbc.workflow.WorkflowConfigParser
            A parsed configuration file that contains the distribution
            options.
        section : str
            Name of the section in the configuration file.
        variable_args : str
            The names of the parameters for this distribution, separated by
            ``VARARGS_DELIM``. These must appear in the "tag" part
            of the section header.
        Returns
        -------
        Uniform
            A distribution instance from the pycbc.inference.prior module.
        """
        return super(UniformIntervals, cls).from_config(cp, section, variable_args,
                     bounds_required=True)


__all__ = ['UniformIntervals']
