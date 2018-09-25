import numpy
from pycbc.distributions import bounded

class UniformDiscreteIntervals(bounded.BoundedDist):
    name = "uniform_discrete_intervals"

def __init__(self, **params):
    super(UniformDiscreteIntervals, self).__init__(**params)

"""
def __init__(self, **params):
    self._bounds = {}
    self._stride = {}
    self._interval = {}
#    self._norm = (1.0 - bnd[0]) / (bnd[1] - bnd[0] + 1.0)
#    self._log_norm = numpy.log(self._norm)
    self._stride = [p for p in params if p.endswith('_stride')]
    self._interval = [p for p in params if p.endswith('_interval')]

#    super(UniformDiscreteIntervals, self).__init__(**params)
    print "Hello World"
    print self._stride
    print self._interval
    print self._bounds[params][0]
    print self._bounds[params][1]
    for p,bnds in self._bounds.items():
        a, b = bnds

#    super(UniformDiscreteIntervals, self).__init__(**params)
#def norm(self):
#    return self._norm

#def log_norm(self):
#    return self._log_norm

@property
def stride(self):
    return self._stride

@property
def interval(self):
    return self._interval

def rvs(self, size=1, param=None):
    if param is not None:
        dtype = [(param, float)]
    else:
        dtype = [(p, float) for p in self.params]

    arr = numpy.zeros(size, dtype=dtype)
    for (p,_) in dtype:
        x = numpy.arange(self._bounds[p][0], self._bounds[p][1],
                         self._stride + self._interval)
        arr[p] = numpy.random.uniform(x, x + self._stride, size=size)
    return arr
"""
__all__ = ['UniformDiscreteIntervals']
