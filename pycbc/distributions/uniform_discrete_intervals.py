import numpy
from pycbc.distributions import uniform
from pycbc.distributions import bounded

class UniformIntervals(bounded.BoundedDist):
    name = 'uniform_intervals'

    def __init__(self, **params):
         # temporarily suppress numpy divide by 0 warning
        self._stride = {}
        self._size = 1

        stride_args = [p for p in params if p.endswith('_stride')]

        self._stride = dict([[p.split("_stride")[0], params.pop(p)] for p in stride_args])

        uni_distr_obj = super(UniformIntervals, self).__init__(**params)
        missing = set(self._stride.keys()) - set(params.keys())

        if any(missing):
            raise ValueError("stride provided for unknown params {}".format(
                             ', '.join(missing)))
        self._stride.update(dict([[p, 0.]
            for p in params if p not in self._stride]))

        numpy.seterr(divide='ignore')
        self._lognorm = -sum([numpy.log(abs(bnd[1]-bnd[0]))
                                    for bnd in self._bounds.values()])
        self._norm = numpy.exp(self._lognorm)
        numpy.seterr(divide='warn')

### FIXME this needs to be reset if multiple calls to rvs will keep shrinking it
    @property
    def norm(self):
        self._norm = self._size
        return self._norm

    @property
    def lognorm(self):
        return numpy.log(norm())

    @property
    def stride(self):
        return self._stride

    def set_size(self, size=1):
        self._size = size

    def _pdf(self, **kwargs):
        """Returns the pdf at the given values. The keyword arguments must
        contain all of parameters in self's params. Unrecognized arguments are
        ignored.
        """
        for p in self.params:
            width = self._bounds[p][0] + self._bounds[p][1]
            print "width of bounds", width
            print "width of bounds + stride", width + self._stride[p]
            print "input value", kwargs[p]
            cond = kwargs[p] % (width + self._stride[p])
            print "cond", cond
            print cond > width
            if cond > width :
                return 0.

        return self._norm

    def _logpdf(self, size=1, **kwargs):
        """Returns the log of the pdf at the given values. The keyword
        arguments must contain all of parameters in self's params. Unrecognized
        arguments are ignored.
        """
        if kwargs in self:
            return numpy.log(self._pdf(size=size))
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
            x = numpy.arange(0, size, 1)
            b = self.bounds[p][1] + self.stride[p] * x
            a = self.bounds[p][0] + self.stride[p] * x
            arr[p] = numpy.random.uniform(a, b) 

        self.set_size(size)
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
