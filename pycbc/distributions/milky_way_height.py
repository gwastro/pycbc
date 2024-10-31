import logging
import numpy
import scipy
from pycbc.distributions import bounded

logger = logging.getLogger('pycbc.distributions.milky_way_height')

class MilkyWayHeight(bounded.BoundedDist):
    name = 'milky_way_height'
    '''
    Returns the height (z-cordinate) in galactocentic frame following a symmetric exponential decay for 3 discs
    Each disc reuires a weighting factor and a scale factor
    
    the config file should be
    [prior-z]
    name = milky_way_height
    min-z = 
    max-z =
    z_scale1 =
    z_scale2 = 
    z_scale3 = 
    z_weight1 = 
    z_weight2 = 
    z_weight3 =
    
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
        ## Take out the elements in the params dict that ends with _scale by popping it out
        ## so that the remaining elements can be passed to BoundedDist to create a bounded.dist
        ## object we can call bounds from later.
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
        super(MilkyWayHeight, self).__init__(**params)
        ## raising error if scale provided for param not in kwargs
        #missing = set(self._scale.keys()) - set(params.keys())
        #if any(missing):
            #raise ValueError(f"scales provided for unknown params {missing}")
        
        ## Set default scale to 1 if scale not provided for param in kwargs
        #self._scale.update(dict([[p,1.] for p in params if p not in self._scale]))
        
        ## Calculate the norm and lognorm for each parameter and fill in the empty dics we created.
        for p,bounds in self._bounds.items():
            self._norm[p] = 1.
            self._lognorm[p] = -numpy.inf
        
    @property
    def norm(self):
        return 1.
    def lognorm(self):
        return -numpy.inf
        
    def _logpdf(self, **kwargs):
        if kwargs in self:
                return sum([numpy.log((self._weight1[p]/(2*self._scale1[p]*(1 - numpy.e**(self._bounds[p][0]/self._scale1[p]))))*numpy.e**(-numpy.abs(kwargs[p])/self._scale1[p])
                                + (self._weight2[p]/(2*self._scale2[p]*(1 - numpy.e**(self._bounds[p][0]/self._scale2[p]))))*numpy.e**(-numpy.abs(kwargs[p])/self._scale2[p])
                                + (self._weight3[p]/(2*self._scale3[p]*(1 - numpy.e**(self._bounds[p][0]/self._scale3[p]))))*numpy.e**(-numpy.abs(kwargs[p])/self._scale3[p])) for p in self._params])
        else:
            return -numpy.inf
    def _pdf(self, **kwargs):
        return numpy.e**(self._logpdf(**kwargs))
        
    def _cdfinv_param(self, param, value):
        z_min = self._bounds[param][0]
        z_max = self._bounds[param][1]
        z_d1 = self._scale1[param]
        z_d2 = self._scale2[param]
        z_d3 = self._scale3[param]
        w1 = self._weight1[param]
        w2 = self._weight2[param]
        w3 = self._weight3[param]
        beta1 = w1/(2*(1 - numpy.e**(z_min/z_d1)))
        beta2 = w2/(2*(1 - numpy.e**(z_min/z_d2)))
        beta3 = w3/(2*(1 - numpy.e**(z_min/z_d3)))
        param_array = numpy.linspace(z_min, z_max, 500)
        mask = param_array <= 0.0
        _cdf = numpy.zeros(len(param_array))
        _cdf[mask] = (beta1)*(numpy.e**(param_array[mask]/z_d1) - numpy.e**(z_min/z_d1)) + (beta2)*(numpy.e**(param_array[mask]/z_d2) - numpy.e**(z_min/z_d2)) + (beta3)*(numpy.e**(param_array[mask]/z_d3) - numpy.e**(z_min/z_d3))
        _cdf[~mask] = (beta1)*(2 - (numpy.e**(z_min/z_d1) + numpy.e**(-param_array[~mask]/z_d1))) + (beta2)*(2 - (numpy.e**(z_min/z_d2) + numpy.e**(-param_array[~mask]/z_d2))) + (beta3)*(2 - (numpy.e**(z_min/z_d3) + numpy.e**(-param_array[~mask]/z_d3)))
        _cdfinv = scipy.interpolate.interp1d(_cdf, param_array)
        return _cdfinv(value)
    
   # def from_config(cls, cp, section, variable_args):
       # return super(SymExpDecay, cls).from_config(cp, section, variable_args, bounds_required = True)
    
__all__ = ['MilkyWayHeight']
