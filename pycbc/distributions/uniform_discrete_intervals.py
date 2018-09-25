import numpy
from pycbc.distributions import uniform

class UniformDiscreteIntervals(uniform.Uniform):
    """


    """
    name = "uniform_discrete_intervals"

def __init__(self, **params):
    super(UniformDiscreteIntervals, self).__init__(**params)
    self._norm = (1.0 - bnd[0]) / (bnd[1] - bnd[0] + 1.0)
    self._log_norm = numpy.log(self._norm)

def norm(self):
    return self._norm

def log_norm(self):
    return self._log_norm

def rvs(self, size=1, param=None):
    """
    """
    if param is not None:
        dtype = [(param, float)]
    else:
        dtype = [(p, float) for p in self.params]

    arr = numpy.zeros(size, dtype=dtype)
    for (p,_) in dtype:
        n = self._bounds[p][1] - self.bounds[p][0] + 1.0
        arr[p] = numpy.floor(numpy.random.uniform(size=size) * n + self._bounds[p][0] - 1.0) 

    return arr

__all__ = ['UniformDiscreteIntervals']
