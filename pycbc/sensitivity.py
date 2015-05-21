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

def volume_montecarlo(found_d, missed_d, found_mchirp, 
                missed_mchirp, distance_param, 
                distance_distribution):
    """ Compute the sensitive volume and standard error using a direct
    Monte Carlo integral

    Parameters
    -----------
    found_distance: numpy.ndarray
        The distances of found injections
    missed_dsistance: numpy.ndarray
        The distances of missed injections
    found_mchirp : numpy.ndarray
        Chirp mass of found injections
    missed_mchirp : numpy.ndarray
        Chirp mass of missed injections
    distance_param: 
        Parameter of the injections used to generate a distribution over distance.
        -may be 'distance', 'chirp_distance".
    distance_distribution: 
        form of the distribution over the parameter
      - may be 'log' (uniform in log D), 'uniform' (uniform in D), '
        distancesquared' (uniform in D**2),
        'volume' (uniform in D***3)
      - It is assumed that injections were carried out over a range of D such 
      that the sensitive
        volume due to signals at distance below D_min is negligible and the 
        efficiency at distances
        above D_max is negligibly small
        
    Returns
    --------
    volume: float
        Volume estimate
    volume_error: float
        The standared error in the volume
    """
    d_weight_power = {
        'log'             : 3.,
        'uniform'         : 2.,
        'distancesquared' : 1.,
        'volume'          : 0.
    }[distance_distribution]
    
    # if no max distance param given, use maximum physical distance actually injected
    max_distance = missed_d.max()

    # all montecarlo integrals are proportional to the volume out to max_distance
    montecarlo_vtot = (4. / 3.) * numpy.pi * max_distance ** 3.

    # set up arrays of weights for the MC average of efficiency
    if distance_param == 'distance':
        found_weights = found_d ** d_weight_power
        missed_weights = missed_d ** d_weight_power
        
    elif distance_param == 'chirp_distance':
        mchirp_weight_power = {
            'log'             : 0.,
            'uniform'         : 5. / 6.,
            'distancesquared' : 5. / 3.,
            'volume'          : 15. / 6.
        }[distance_distribution]
        
        # for a distribution over dchirp, weights get a power of mchirp to rescale
        # the injection density to the target mass distribution
        found_weights = found_d ** d_weight_power * found_mchirp ** mchirp_weight_power
        missed_weights = missed_d ** d_weight_power * missed_mchirp ** mchirp_weight_power
    else: 
        raise NotImplementedError("%s is not a recognized distance parameter" % distance_param)

    # measured weighted efficiency is w_i for a found inj and 0 for missed
    mc_weight_samples = numpy.concatenate((found_weights, 0*missed_weights))

    # MC integral is total volume of sphere times sum of found weights over sum of all weights
    # Treat (total volume / sum of weights) as a constant prefactor
    mc_norm = sum(found_weights) + sum(missed_weights)
    mc_prefactor = montecarlo_vtot / mc_norm
    mc_sum = sum(mc_weight_samples)

    # Sample variance of injection efficiency: mean of the square - square of the mean
    Ninj = len(mc_weight_samples)
    mc_sample_variance = sum(mc_weight_samples ** 2.) / Ninj - (mc_sum / Ninj) ** 2.
    # Variance of sum over efficiencies scales up with Ninj (Bienayme' rule)
    mc_sum_variance = Ninj * mc_sample_variance

    # return MC integral and its standard deviation
    vol, vol_err = mc_prefactor * mc_sum, mc_prefactor * mc_sum_variance ** 0.5
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
