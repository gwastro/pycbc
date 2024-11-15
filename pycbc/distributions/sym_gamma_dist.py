import logging
import numpy
import scipy
from pycbc.distributions import bounded

logger = logging.getLogger('pycbc.distributions.sym_gamma_dist')

class SymGammaDist(bounded.BoundedDist):
    name = 'sym_gamma_dist'
    """
    Samples from a weighted sum of one dimensional symmetric gamma distributions charecterised by a power $n$,
    scale factors $s_{i}$, and weights $w_{i}$
    \[
    p(x) = \sum_{i} w_{i} |x|^{n} exp(-|x|/s_{i})
    \]
    The inverse cdf is calculated by interpolating the cdf over the range of the variable
    The cdf is calculated using regularised lower incomplete gamma function as defined in scipy.special

    Parameters
    -----------
    Example config file

    [prior-x]
    name = sym_gamma_dist
    min-x = -10
    max-x = 10
    scales = 1,2
    weights = 0.5,0.5
    power = 2
    interp_points = 500

    """



    def __init__(self, scales = None, weights = None,
    power = None, interp_points = None, **params):
        super(SymGammaDist, self).__init__(**params)
        self._scales = {}
        self._weights = {}
        self._power = {}
        self._interp_points = {}
        self._norms = {}
        for p, bounds in self._bounds.items():
            if type(scales) == str:
            	self._scales[p] = [float(s) for s in scales.split(",")]
            else:
            	self._scales[p] = [float(scales)]
            if type(weights) == str:
            	self._weights[p] = [float(w) for w in weights.split(",")]
            else:
            	self._weights[p] = [float(weights)]
            self._power[p] = float(power)
            self._interp_points[p] = int(interp_points)
            if len(self._scales[p]) != len(self._weights[p]):
                raise ValueError(f'Unequal number of scales and weights'
                + f'provided for parameter {p}')
            if numpy.sum(self._weights[p]) != float(1):
            	raise ValueError('Weights should add to 1.0')
            self._norms[p] = [0]*len(self._scales[p])
            if numpy.sign(bounds[0]) == numpy.sign(bounds[1]):
                for i in range(len(self._norms[p])):
                    self._norms[p][i] = numpy.abs(scipy.special.factorial(self._power[p])
                                                  *self._scales[p][i]**(self._power[p]+1)
                                                  *(scipy.special.gammainc(self._power[p]+1, numpy.abs(bounds[0]/self._scales[p][i]))
                                                    -scipy.special.gammainc(self._power[p]+1, numpy.abs(bounds[1]/self._scales[p][i]))))
            if numpy.sign(bounds[0]) != numpy.sign(bounds[1]):
                for i in range(len(self._norms[p])):
                    self._norms[p][i] = numpy.abs(scipy.special.factorial(self._power[p])
                                                  *self._scales[p][i]**(self._power[p]+1)
                                                  *(scipy.special.gammainc(self._power[p]+1, numpy.abs(bounds[0]/self._scales[p][i]))
                                                    +scipy.special.gammainc(self._power[p]+1, numpy.abs(bounds[1]/self._scales[p][i]))))

    def c0(self, x, power, scale):
        return scipy.special.gammainc(power+1,numpy.abs(x/scale))*\
        scipy.special.factorial(power)*scale**(power+1)
                    
    def _logpdf(self, **kwargs):
        if kwargs in self:
            lpdf_p = 0
            for p in self._params:
                pdf_i = 0
                for i in range(len(self._scales[p])):
                    pdf_i += (self._weights[p][i]/self._norms[p][i])* \
                        numpy.abs(kwargs[p])**(self._power[p])* \
                        numpy.e**(-numpy.abs(kwargs[p])/self._scales[p][i])
                lpdf_p += numpy.log(pdf_i)
            return lpdf_p
        else:
            return -numpy.inf
    def _pdf(self, **kwargs):
        return numpy.e**(self._logpdf(**kwargs))

    def _cdfinv_param(self, param, value):
        param_min, param_max = self._bounds[param][0], self._bounds[param][1]
        param_array = numpy.linspace(param_min,
                                     param_max,
                                     self._interp_points[param])
        
        if numpy.sign(param_min) != numpy.sign(param_max):
            cdf_array = numpy.zeros(self._interp_points[param])
            for i in range(len(self._scales[param])):
                mask = param_array <= 0
                cdf_array[mask] = cdf_array[mask] + (self._weights[param][i]/self._norms[param][i])* \
                                    (self.c0(numpy.abs(param_min), self._power[param], self._scales[param][i])
                                    -self.c0(numpy.abs(param_array[mask]), self._power[param], self._scales[param][i]))
                cdf_array[~mask] = cdf_array[~mask] + (self._weights[param][i]/self._norms[param][i])* \
                                    (self.c0(param_min, self._power[param], self._scales[param][i])
                                    +self.c0(param_array[~mask], self._power[param], self._scales[param][i]))
            inv_cdf_interpolate = scipy.interpolate.interp1d(cdf_array, param_array)
        
            return inv_cdf_interpolate(value)
        if numpy.sign(param_min) == numpy.sign(param_max):
            cdf_array = numpy.zeros(self._interp_points[param])
            for i in range(len(self._scales[param])):
                cdf_array += (self._weights[param][i]/self._norms[param][i])* \
                                numpy.abs(self.c0(param_min, self._power[param], self._scales[param][i])
                                          -self.c0(param_array, self._power[param], self._scales[param][i]))
            inv_cdf_interpolate = scipy.interpolate.interp1d(cdf_array, param_array)
            return inv_cdf_interpolate(value)
__all__ = ['SymGammaDist']
