import logging
import numpy
import scipy
from pycbc.distributions import bounded
from numpy import sign, abs
from scipy.special import factorial, gammainc
from scipy.interpolate import interp1d

logger = logging.getLogger('pycbc.distributions.sym_gamma_dist')

class SymGammaDist(bounded.BoundedDist):
    r"""Samples from a weighted sum of one dimensional symmetric gamma 
    distributions between the bounds provided.
    The PDF is given by 
    .. math::
        p(x) = \sum_{i} w_{i} |x|^{n} e^{\frac{-|x|}{s_{i}}}
    
    where :math:'w_{i}' are the weights, :math:'n' is the power,
    and :math:'s_{i}' are the scale factors of the individual distributions

    The CDF between the bounds :math:'[a,b]' is given by
    .. math::
        c(r) = \sum_{i} w_{i} \int_{a}^{r} |x|^{n} e^{\frac{-|x|}{s_{i}}}

    The CDF is calculated using a rescaling of scipy's implementation of 
    regularised lower incomplete gamma function gammainc
    :math:'\gamma(n+1,\frac{r}{s_{i}})' as follows

    ..math::
        c(r) = |\gamma(n+1,\frac{r}{s_{i}}) - \sigma(r,b)\gamma(n+1,\frac{b}{s_{i}})|

    where :math:'\sigma(u, v) = 1' if sign(u)=sign(v) else -1.

    The inverse cdf is calculated by interpolating the cdf over the 
    range of the parameter. The number of points over which the 
    interpolation is done is 5000 by default.
    
    Parameters
    -----------
    Bounds : The minimum and maximum range of the parameter. Can be 
        provided as bounded.dist object
    scales : The scale factors for the individual distributions. Must be
        provided as a space separated list
    weights : The weighting factors for the individual distributions. Must
        be provided as a space separated list. Can provide unnormalized
        weights. Defaults to 1 if only one distribution is used.
    power : The power for the gamma distribution.
    interp_points : The number of points over which the CDF is evaluated
        to interpolate the inverse CDF

    Example config file

    [prior-x]
    name = sym_gamma_dist
    min-x = -10
    max-x = 10
    scales = 3 5
    weights = 1 3
    power = 2
    interp_points = 1000
    """

    name = 'sym_gamma_dist'

    def __init__(self,scales=None,weights=None,
        power=None,interp_points=5000,**params):
        super(SymGammaDist, self).__init__(**params)

        if isinstance(scales,str):
            self._scales = [float(s) for s in scales.split()]
        else:
            self._scales = [float(scales)]

        if isinstance(weights,str):
            weight_floats = [float(w) for w in weights.split()]
            self._weights = [w/sum(weight_floats) for w in weight_floats]
        else:
            self._weights = [float(1.0)]

        if len(self._scales) != len(self._weights):
            raise ValueError("Unequal number of scales and weights provided")
        
        self._power = float(power)
        self._interp_points = int(interp_points)
        self._interpolated_invcdf = {}
        self._norms = {}
        for p, bounds in self._bounds.items():
            lower, upper = bounds[0], bounds[1]
            param_array = numpy.linspace(lower,upper,self._interp_points)
            cdf_array = numpy.zeros(len(param_array))
            norms = numpy.zeros(len(self._scales))
            for i in range(len(self._scales)):
                norms[i] = self.integral_gamma(
                    lower,upper,self._power,self._scales[i]
                    )
                cdf_array += (
                    self._weights[i]
                    *self.integral_gamma(
                        lower,param_array,self._power,self._scales[i]
                        )
                        /norms[i]
                    )
            self._interpolated_invcdf[p] = interp1d(cdf_array,param_array)
            self._norms[p] = norms

    def rescaled_gammainc(self,x,power,scale):
        """ Rescales the lower incomplete gamma function
        """
        return gammainc(power+1,abs(x/scale))*\
        factorial(power)*scale**(power+1)

    def integral_gamma(self,lower,upper,power,scale):
        """ The definite integral of the indivdual PDF between the
        limits 'lower' and 'upper'
        """
        lower, upper = numpy.asarray(lower), numpy.asarray(upper)
        rescaled_upper = self.rescaled_gammainc(upper,power,scale)
        rescaled_lower = self.rescaled_gammainc(lower,power,scale)
        sign_factor = numpy.where(sign(lower)==sign(upper),1.0,-1.0)
        return numpy.abs(rescaled_upper - sign_factor*rescaled_lower)
                   
    def _logpdf(self,**kwargs):
        if kwargs in self:
            lpdf_p = 0
            for p in self._params:
                pdf_i = 0
                for i in range(len(self._scales)):
                    pdf_i += (self._weights[i]/self._norms[p][i])* \
                        abs(kwargs[p])**(self._power)* \
                        numpy.e**(-abs(kwargs[p])/self._scales[i])
                lpdf_p += numpy.log(pdf_i)
            return lpdf_p
        else:
            return -numpy.inf
    def _pdf(self,**kwargs):
        return numpy.e**(self._logpdf(**kwargs))

    def _cdfinv_param(self,param,value):
        invcdf = self._interpolated_invcdf[param]
        return invcdf(value)


__all__ = ['SymGammaDist']