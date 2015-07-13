""" This module contains utilities for calculating search sensitivity
"""
import numpy

def volume_to_distance_with_errors(vol, vol_err):
    """ Return the distance and standard deviation upper and lower bounds
    
    Parameters
    ----------
    volume: float
    volume_dev: float
    
    Returns
    -------
    distance: float
    err_upper: float
    err_lower: float
    
    """
    dist = (vol * 3.0/4.0/numpy.pi) ** (1.0/3.0)
    ehigh = ((vol + vol_err) * 3.0/4.0/numpy.pi) ** (1.0/3.0) - dist
    elow = dist - ((vol - vol_err) * 3.0/4.0/numpy.pi) ** (1.0/3.0)
    return dist, ehigh, elow

def volume_montecarlo(found_d, missed_d, found_mchirp, missed_mchirp,
                      distribution_param, distribution, limits_param,
                      max_param=None, min_param=None):
    """
    Compute the sensitive volume and standard error using a direct Monte Carlo
    integral.  For the result to be useful injections should be made over a
    range of distances D such that sensitive volume due to signals closer than
    D_min is negligible, and efficiency at distances above D_max is negligible

    Parameters
    -----------
    found_d: numpy.ndarray
        The distances of found injections
    missed_d: numpy.ndarray
        The distances of missed injections
    found_mchirp: numpy.ndarray
        Chirp mass of found injections
    missed_mchirp: numpy.ndarray
        Chirp mass of missed injections
    distribution_param: string
        Parameter D of the injections used to generate a distribution over
        distance, may be 'distance', 'chirp_distance".
    distribution: string
        form of the distribution over the parameter, may be 
        'log' (uniform in log D)
        'uniform' (uniform in D)
        'distancesquared' (uniform in D**2)
        'volume' (uniform in D***3)
    limits_param: string
        Parameter Dlim specifying limits inside which injections were made
        may be 'distance', 'chirp distance'
    max_param: float
        maximum value of Dlim out to which injections were made; if None
        the maximum actually injected value will be used
    min_param: float
        minimum value of Dlim at which injections were made; only used for
        log distribution, then if None the minimum actually injected value
        will be used

    Returns
    --------
    volume: float
        Volume estimate
    volume_error: float
        The standard error in the volume
    """
    d_power = {
        'log'             : 3.,
        'uniform'         : 2.,
        'distancesquared' : 1.,
        'volume'          : 0.
    }[distribution]
    mchirp_power = {
            'log'             : 0.,
            'uniform'         : 5. / 6.,
            'distancesquared' : 5. / 3.,
            'volume'          : 15. / 6.
    }[distribution]

    # establish maximum physical distance: first in case of chirp distance distribution
    if limits_param == 'chirp_distance':
        mchirp_standard_bns = 1.4 * (2. ** (-1. / 5.))
        all_mchirp = numpy.concatenate((found_mchirp, missed_mchirp))
        max_mchirp = all_mchirp.max()
        if max_param is not None:
            # use largest actually injected mchirp for conversion
            max_distance = max_param * (max_mchirp / mchirp_standard_bns) ** (5. / 6.)
    elif limits_param == 'distance' and max_param is not None:
        max_distance = max_param
    else: raise NotImplementedError("%s is not a recognized parameter" % limits_param)

    # if no max distance param given, use maximum physical distance actually injected
    if max_param == None:
        max_distance = max(found_d.max(), missed_d.max())

    # volume of sphere
    montecarlo_vtot = (4. / 3.) * numpy.pi * max_distance ** 3.

    # arrays of weights for the MC integral
    if distribution_param == 'distance':
        found_weights = found_d ** d_power
        missed_weights = missed_d ** d_power
    elif distribution_param == 'chirp_distance':
        # weight by a power of mchirp to rescale injection density to the
        # target mass distribution
        found_weights = found_d ** d_power * \
                        found_mchirp ** mchirp_power
        missed_weights = missed_d ** d_power * \
                         missed_mchirp ** mchirp_power
    else:
        raise NotImplementedError("%s is not a recognized distance parameter"
                                                              % distance_param)

    all_weights = numpy.concatenate((found_weights, missed_weights))

    # measured weighted efficiency is w_i for a found inj and 0 for missed
    mc_weight_samples = numpy.concatenate((found_weights, 0*missed_weights))

    # MC integral is volume of sphere * (sum of found weights)/(sum of all weights)
    # over injections covering the sphere
    mc_weight_samples = numpy.concatenate((found_weights, 0 * missed_weights))
    mc_sum = sum(mc_weight_samples)

    if limits_param == 'distance':
        mc_norm = sum(all_weights)
    elif limits_param == 'chirp_distance':
        # if injections are made up to a maximum chirp distance, account for
        # extra missed injections that would occur when injecting up to
        # maximum physical distance : this works out to a 'chirp volume' factor
        mc_norm = sum(all_weights * (max_mchirp / all_mchirp) ** (5. / 2.))

    # take out a constant factor
    mc_prefactor = montecarlo_vtot / mc_norm

    # count the samples
    if limits_param == 'distance':
        Ninj = len(mc_weight_samples)
    elif limits_param == 'chirp_distance':
        # find the total expected number after extending from maximum chirp
        # dist up to maximum physical distance
        if distribution == 'log':
            # only need minimum distance in this one case
            if min_param is not None:
                min_distance = min_param * \
                     (numpy.min(all_mchirp) / mchirp_standard_bns) ** (5. / 6.)
            else:
                min_distance = min(numpy.min(found_d), numpy.min(missed_d))
            logrange = numpy.log(max_distance / min_distance)
            Ninj = len(mc_weight_samples) + (5. / 6.) * \
                             sum(numpy.log(max_mchirp / all_mchirp) / logrange)
        else:
            Ninj = sum((max_mchirp / all_mchirp) ** mchirp_power)

    # sample variance of efficiency: mean of the square - square of the mean
    mc_sample_variance = sum(mc_weight_samples ** 2.) / Ninj - \
                                                          (mc_sum / Ninj) ** 2.

    # return MC integral and its standard deviation; variance of mc_sum scales
    # relative to sample variance by Ninj (Bienayme' rule)
    vol = mc_prefactor * mc_sum
    vol_err = mc_prefactor * (Ninj * mc_sample_variance) ** 0.5
    return vol, vol_err

def volume_binned_pylal(f_dist, m_dist, bins=15):
    """ Compute the sensitive volume using a distanced 
    binned efficiency estimate
    
    Parameters
    -----------
    found_distance: numpy.ndarray
        The distances of found injections
    missed_dsistance: numpy.ndarray
        The distances of missed injections
        
    Returns
    --------
    volume: float
        Volume estimate
    volume_error: float
        The standared error in the volume
    """
    def sims_to_bin(sim):
        return (sim, 0)

    from pylal import rate
    from pylal.imr_utils import compute_search_volume_in_bins, compute_search_efficiency_in_bins
    found = f_dist
    total = numpy.concatenate([f_dist, m_dist])
    ndbins = rate.NDBins([rate.LinearBins(min(total), max(total), bins), rate.LinearBins(0., 1, 1)]) 
    vol, verr = compute_search_volume_in_bins(found, total, ndbins, sims_to_bin)
    return vol.array[0], verr.array[0]

def volume_shell(f_dist, m_dist):
    """ Compute the sensitive volume using sum over spherical shells.
    
    Parameters
    -----------
    found_distance: numpy.ndarray
        The distances of found injections
    missed_dsistance: numpy.ndarray
        The distances of missed injections
        
    Returns
    --------
    volume: float
        Volume estimate
    volume_error: float
        The standared error in the volume
    """
    f_dist.sort()
    m_dist.sort()
    distances = numpy.concatenate([f_dist, m_dist])
    dist_sorting = distances.argsort()
    distances = distances[dist_sorting]
    low = 0
    vol = 0
    vol_err = 0
    for i in range(len(distances)):
        if i == len(distances) - 1:
            break
    
        high = (distances[i+1] + distances[i]) / 2
        bin_width = high - low
        
        if dist_sorting[i] < len(f_dist):
            vol += 4 * numpy.pi * distances[i]**2.0 * bin_width
            vol_err += (4 * numpy.pi * distances[i]**2.0 * bin_width)**2.0
            
        low = high
    vol_err = vol_err ** 0.5
    return vol, vol_err
