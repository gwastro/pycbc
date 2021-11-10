""" Common utility functions for calculation of likelihoods
"""
import numpy
import tqdm
import logging
from scipy.special import logsumexp, i0e
from scipy.integrate import quad
from scipy.interpolate import RectBivariateSpline
from pycbc.distributions import JointDistribution


class DistMarg(object):
    """Help class to add bookkeeping for distance marginalization"""
    
    def setup_distance_marginalization(self,
         variable_params,
         marginalize_phase=False,
         marginalize_distance=False,
         marginalize_distance_param='distance',
         marginalize_distance_samples=int(1e4),
         marginalize_distance_interpolator=False,
         **kwargs):
             
        self.marginalize_phase = marginalize_phase
        self.distance_marginalization = False     
        if not marginalize_distance:
            return variable_params, kwargs
            
        logging.info('Marginalizing over distance')

        # Take distance out of the variable params since we'll handle it
        # manually now
        variable_params.remove(marginalize_distance_param)
        old_prior = kwargs['prior']
        dists = [d for d in old_prior.distributions \
                 if marginalize_distance_param not in d.params]
        dprior = [d for d in old_prior.distributions \
                 if marginalize_distance_param in d.params][0]
        prior = JointDistribution(variable_params, *dists,
                                            **old_prior.kwargs)
        kwargs['prior'] = prior

        if len(dprior.params) != 1 or not hasattr(dprior, 'bounds'):
            raise ValueError('Distance Marginalization requires a '
                             'univariate and bounded prior')                 

        # Set up distance prior vector and samples
        # (1) prior is using distance
        if dprior.params[0] == 'distance':
            logging.info("Prior is directly on distance, setting up "
                         "%s grid weights", marginalize_distance_samples)
            dmin, dmax = dprior.bounds['distance']
            dist_locs = numpy.linspace(dmin, dmax, 
                                       int(marginalize_distance_samples))
            dist_weights = [dprior.pdf(distance=l) for l in dist_locs]
            dist_weights = numpy.array(dist_weights)

        # (2) prior is univariate and can be converted to distance
        # TODO
        elif marginalize_distance_param != 'distance':
            pass
        else:
            raise ValueError("No prior seems to determine the distance")

        dist_weights /= dist_weights.sum()
        dist_ref = 0.5 * (dmax + dmin)
        self.distance_marginalization = dist_ref / dist_locs, dist_weights 
        self.distance_interpolator = None
        if marginalize_distance_interpolator:
            i = setup_distance_marg_interpolant(self.distance_marginalization,
                                                phase=self.marginalize_phase)
            self.distance_interpolator = i
        kwargs['static_params']['distance'] = dist_ref
        return variable_params, kwargs

    def marginalize_loglr(sh_total, hh_total):
        return marginalize_likelihood(sh_total, hh_total,
                              phase=self.marginalize_phase,
                              interpolator=self.distance_interpolator,
                              distance=self.distance_marginalization)            

def setup_distance_marg_interpolant(dist_marg,
                                    phase=False,
                                    snr_range=(1, 50),
                                    density=(1000, 1000)):
    """ Create the interpolant for distance marginalization
    
    Parameters
    ----------
    dist_marg: tuple of two arrays
        The (dist_loc, dist_weight) tuple which defines the grid
        for integrating over distance
    snr_range: tuple of (float, float)
        Tuple of min, max SNR that the interpolant is expected to work
        for.
    density: tuple of (float, float)
        The number of samples in either dimension of the 2d interpolant
    """
    (dist_rescale, dist_weights) = dist_marg
    logging.info("Interpolator valid for SNRs in", snr_range)
    logging.info("Interpolator using grid", density)
    # approximate maximum shr and hhr values, assuming the true SNR is
    # within the indicated range (and neglecting noise fluctuations)
    snr_min, snr_max = snr_range
    smax = dist_rescale.max()
    smin = dist_rescale.min()
    shr_max = snr_max ** 2.0 / smin
    hhr_max = snr_max ** 2.0 / smin / smin

    shr_min = snr_min ** 2.0 / smax
    hhr_min = snr_min ** 2.0 / smax / smax

    shr = numpy.geomspace(shr_min, shr_max, density[0])
    hhr = numpy.geomspace(hhr_min, hhr_max, density[1])
    lvals = numpy.zeros((len(shr), len(hhr)))
    logging.info('Setup up likelihood interpolator')
    for i, sh in enumerate(tqdm.tqdm(shr)):
        for j, hh in enumerate(hhr):
            lvals[i, j] = marginalize_likelihood(sh, hh,
                                                 distance=dist_marg,
                                                 phase=phase)
    interp = RectBivariateSpline(shr, hhr, lvals)
    
    def interp_wrapper(x, y):
        return interp(x, y, grid=False)
    return interp_wrapper

def marginalize_likelihood(sh, hh,
                           phase=False,
                           distance=False,
                           interpolator=None):
    """ Return the marginalized likelihood
    
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

    if phase:
        sh = abs(sh)
    else:
        sh = sh.real

    vweights = 1
    if interpolator:
        # pre-calculated result for this function
        vloglr = interpolator(sh, hh)
    else:
        #explicit calculation
        if distance:
            # brute force distance path
            dist_rescale, dist_weights = distance
            
            sh = numpy.multiply.outer(sh, dist_rescale) 
            hh = numpy.multiply.outer(hh, dist_rescale ** 2.0)
            if len(sh.shape) == 2:
                vweights = numpy.resize(dist_weights,
                                        (sh.shape[1], sh.shape[0])).T
            else:
                vweights = dist_weights
        
        if phase:
            sh = numpy.log(i0e(sh)) + sh
        else:
            sh = sh.real
            
        # Calculate loglikelihood ratio
        vloglr = sh - 0.5 * hh 
    
    # Do brute-force marginalization if loglr is a vector
    if isinstance(vloglr, float):
        vloglr = float(vloglr)
    else:
        vloglr = float(logsumexp(vloglr, b=vweights)) - clogweights

    return vloglr
