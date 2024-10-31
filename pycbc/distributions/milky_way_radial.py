import logging
import numpy
import scipy
from pycbc.distributions import bounded

logger = logging.getLogger('pycbc.distributions.milky_way_radial')

class MilkyWayRadial(bounded.BoundedDist):
    name = 'milky_way_radial'
    '''
    Returns the radius (r-cordinate) in galactocentic frame following an exponential decaydistribution for 3 discs
    Each disc reuires a weighting factor and a scale factor
    
    the config file should be
    [prior-r]
    name = milky_way_height
    min-r =    # Should not be negative
    max-r =
    r_scale1 =
    r_scale2 = 
    r_scale3 = 
    r_weight1 = 
    r_weight2 = 
    r_weight3 =
    
    '''
    def __init__(self, **params):
        self._bounds = {}
        self._scale1 = {}
        self._scale2 = {}
        self._scale3 = {}
        self._weight1 = {}
        self._weight2 = {}
        self._weight3 = {}
        self._norm = {}
        self._lognorm = {}
        self._rarray = {}
        self._cdfarray = {}
        scale1_args = [p for p in params if p.endswith('_scale1')]
        self._scale1 = dict([p[:-7], params.pop(p)] for p in scale1_args)
        scale2_args = [p for p in params if p.endswith('_scale2')]
        self._scale2 = dict([p[:-7], params.pop(p)] for p in scale2_args)
        scale3_args = [p for p in params if p.endswith('_scale3')]
        self._scale3 = dict([p[:-7], params.pop(p)] for p in scale3_args)
        ##
        weight1_args = [p for p in params if p.endswith('_weight1')]
        self._weight1 = dict([p[:-8], params.pop(p)] for p in weight1_args)
        weight2_args = [p for p in params if p.endswith('_weight2')]
        self._weight2 = dict([p[:-8], params.pop(p)] for p in weight2_args)
        weight3_args = [p for p in params if p.endswith('_weight3')]
        self._weight3 = dict([p[:-8], params.pop(p)] for p in weight3_args)
        ## _scale keys removed from the params.
        ## Now pass to bounded.dist
        super(MilkyWayRadial, self).__init__(**params)
        for p, bounds in self._bounds.items():
            scale1 = self._scale1[p]
            scale2 = self._scale2[p]
            scale3 = self._scale3[p]
            weight1 = self._weight1[p]
            weight2 = self._weight2[p]
            weight3 = self._weight3[p]
            r_min, r_max = bounds
            #if r_min < 0:
                #raise ValueError(f'Minimum value of {p} must be greater than or equal to zero ')
            r = numpy.linspace(r_min, r_max, 1000)
            cdf_array = 1. - ((weight1*numpy.e**(-r/scale1)*(1+r/scale1))+
                              (weight2*numpy.e**(-r/scale2)*(1+r/scale2))+
                              (weight3*numpy.e**(-r/scale3)*(1+r/scale3)))
            self._rarray[p] = r
            self._cdfarray[p] = cdf_array
            self._norm[p] = 1.
            self._lognorm[p] = -numpy.inf
    @property
    def norm(self):
        return 1.
    def lognorm(self):
        return numpy.inf
    def _logpdf(self, **kwargs):
        if kwargs in self:
            return sum(
                [ numpy.log(kwargs[p]*((self._weight1[p]*numpy.e**(-kwargs[p]/self._scale1[p])/self._scale1[p]**2)
                                       +(self._weight2[p]*numpy.e**(-kwargs[p]/self._scale2[p])/self._scale2[p]**2)
                                       +(self._weight3[p]*numpy.e**(-kwargs[p]/self._scale3[p])/self._scale3[p]**2))) for p in self._params])
                                       
        else:
        	return -numpy.inf
       
    def _pdf(self, **kwargs):
        return numpy.e**(_logpdf(**kwargs))
    def _cdfinv_param(self, param, value):
        _interpolation = scipy.interpolate.interp1d(self._cdfarray[param], self._rarray[param])
        return _interpolation(value)
    
__all__ = ['MilkWayRadial']
