""" Common utility functions for calculation of likelihoods
"""
import numpy
from scipy.special import logsumexp, i0e
from scipy.integrate import quad

def marginalize_likelihood(sh, hh,
                           phase=False,
                           distance=False):
    """ Return the marginalized likelihood. 
    
    Parameters
    ----------
    sh: complex float or array
    hh: array
    
    Returns
    -------
    loglr: float
        The marginalized loglikehood ratio
    """  
    if isinstance(sh, float):
        clogweights = 0
    else:
        sh = sh.flatten()
        hh = hh.flatten()
        clogweights = numpy.log(len(sh))
        
    vweights = 1
    if distance:
        # brute force distance path
        # scale = dref / dists
        dref, dists, dist_weights = distance
        scale = dref / dists
        
        sh = numpy.multiply.outer(sh, scale) 
        hh = numpy.multiply.outer(hh, scale ** 2.0)
        if len(sh.shape) == 2:
            vweights = numpy.resize(dist_weights, (sh.shape[1], sh.shape[0])).T
        else:
            vweights = dist_weights
    
    if phase:
        vloglr = numpy.log(i0e(abs(sh)))
        vloglr += abs(sh) - 0.5 * hh
    else:
        vloglr = sh.real - 0.5 * hh 
    
    if isinstance(vloglr, float):
        vloglr = float(vloglr)
    else:
        vloglr = float(logsumexp(vloglr, b=vweights)) - clogweights

    return vloglr
