import numpy
import bisect
from . import bin_utils


def integral_element(mu, pdf):
    '''
    Returns an array of elements of the integrand dP = p(mu) dmu
    for a density p(mu) defined at sample values mu ; samples need
    not be equally spaced.  Uses a simple trapezium rule.
    Number of dP elements is 1 - (number of mu samples).
    '''
    dmu = mu[1:] - mu[:-1]
    bin_mean = (pdf[1:] + pdf[:-1]) / 2.
    return dmu * bin_mean


def normalize_pdf(mu, pofmu):
    """
    Takes a function pofmu defined at rate sample values mu and
    normalizes it to be a suitable pdf. Both mu and pofmu must be
    arrays or lists of the same length.
    """
    if min(pofmu) < 0:
        raise ValueError("Probabilities cannot be negative, don't ask me to "
                         "normalize a function with negative values!")
    if min(mu) < 0:
        raise ValueError("Rates cannot be negative, don't ask me to "
                         "normalize a function over a negative domain!")

    dp = integral_element(mu, pofmu)
    return mu, pofmu/sum(dp)


def compute_upper_limit(mu_in, post, alpha=0.9):
    """
    Returns the upper limit mu_high of confidence level alpha for a
    posterior distribution post on the given parameter mu.
    The posterior need not be normalized.
    """
    if 0 < alpha < 1:
        dp = integral_element(mu_in, post)
        high_idx = bisect.bisect_left(dp.cumsum() / dp.sum(), alpha)
        # if alpha is in (0,1] and post is non-negative, bisect_left
        # will always return an index in the range of mu since
        # post.cumsum()/post.sum() will always begin at 0 and end at 1
        mu_high = mu_in[high_idx]
    elif alpha == 1:
        mu_high = numpy.max(mu_in[post > 0])
    else:
        raise ValueError("Confidence level must be in (0,1].")

    return mu_high


def compute_lower_limit(mu_in, post, alpha=0.9):
    """
    Returns the lower limit mu_low of confidence level alpha for a
    posterior distribution post on the given parameter mu.
    The posterior need not be normalized.
    """
    if 0 < alpha < 1:
        dp = integral_element(mu_in, post)
        low_idx = bisect.bisect_right(dp.cumsum() / dp.sum(), 1 - alpha)
        # if alpha is in [0,1) and post is non-negative, bisect_right
        # will always return an index in the range of mu since
        # post.cumsum()/post.sum() will always begin at 0 and end at 1
        mu_low = mu_in[low_idx]
    elif alpha == 1:
        mu_low = numpy.min(mu_in[post > 0])
    else:
        raise ValueError("Confidence level must be in (0,1].")

    return mu_low


def confidence_interval_min_width(mu, post, alpha=0.9):
    '''
    Returns the minimal-width confidence interval [mu_low, mu_high] of
    confidence level alpha for a posterior distribution post on the parameter mu.
    '''
    if not 0 < alpha < 1:
        raise ValueError("Confidence level must be in (0,1).")

    # choose a step size for the sliding confidence window
    alpha_step = 0.01

    # initialize the lower and upper limits
    mu_low = numpy.min(mu)
    mu_high = numpy.max(mu)

    # find the smallest window (by delta-mu) stepping by dalpha
    for ai in numpy.arange(0, 1 - alpha, alpha_step):
        ml = compute_lower_limit(mu, post, 1 - ai)
        mh = compute_upper_limit(mu, post, alpha + ai)
        if mh - ml < mu_high - mu_low:
            mu_low = ml
            mu_high = mh

    return mu_low, mu_high


def hpd_coverage(mu, pdf, thresh):
    '''
    Integrates a pdf over mu taking only bins where
    the mean over the bin is above a given threshold
    This gives the coverage of the HPD interval for
    the given threshold.
    '''
    dp = integral_element(mu, pdf)
    bin_mean = (pdf[1:] + pdf[:-1]) / 2.

    return dp[bin_mean > thresh].sum()


def hpd_threshold(mu_in, post, alpha, tol):
    '''
    For a PDF post over samples mu_in, find a density
    threshold such that the region having higher density
    has coverage of at least alpha, and less than alpha
    plus a given tolerance.
    '''
    norm_post = normalize_pdf(mu_in, post)
    # initialize bisection search
    p_minus = 0.0
    p_plus = max(post)
    while abs(hpd_coverage(mu_in, norm_post, p_minus) -
              hpd_coverage(mu_in, norm_post, p_plus)) >= tol:
        p_test = (p_minus + p_plus) / 2.
        if hpd_coverage(mu_in, post, p_test) >= alpha:
            # test value was too low or just right
            p_minus = p_test
        else:
            # test value was too high
            p_plus = p_test
    # p_minus never goes above the required threshold and p_plus never goes below
    # thus on exiting p_minus is at or below the required threshold and the
    # difference in coverage is within tolerance
    return p_minus


def hpd_credible_interval(mu_in, post, alpha=0.9, tolerance=1e-3):
    '''
    Returns the minimum and maximum rate values of the HPD
    (Highest Posterior Density) credible interval for a posterior
    post defined at the sample values mu_in.  Samples need not be
    uniformly spaced and posterior need not be normalized.

    Will not return a correct credible interval if the posterior
    is multimodal and the correct interval is not contiguous;
    in this case will over-cover by including the whole range from
    minimum to maximum mu.
    '''
    if alpha == 1:
        nonzero_samples = mu_in[post > 0]
        mu_low = numpy.min(nonzero_samples)
        mu_high = numpy.max(nonzero_samples)
    elif 0 < alpha < 1:
        # determine the highest PDF for which the region with
        # higher density has sufficient coverage
        pthresh = hpd_threshold(mu_in, post, alpha, tol=tolerance)
        samples_over_threshold = mu_in[post > pthresh]
        mu_low = numpy.min(samples_over_threshold)
        mu_high = numpy.max(samples_over_threshold)

    return mu_low, mu_high


# Following functions are for the old pylal volume vs mass calculations
# These were replaced by 'imr_utils' functions now contained in sensitivity.py
# and bin_utils.py


def integrate_efficiency(dbins, eff, err=0, logbins=False):

    if logbins:
        logd = numpy.log(dbins)
        dlogd = logd[1:] - logd[:-1]
        # use log midpoint of bins
        dreps = numpy.exp((numpy.log(dbins[1:]) + numpy.log(dbins[:-1])) / 2.)
        vol = numpy.sum(4.*numpy.pi * dreps**3. * eff * dlogd)
        # propagate errors in eff to errors in v
        verr = numpy.sqrt(
            numpy.sum((4.*numpy.pi * dreps**3. * err * dlogd)**2.)
        )
    else:
        dd = dbins[1:] - dbins[:-1]
        dreps = (dbins[1:] + dbins[:-1]) / 2.
        vol = numpy.sum(4. * numpy.pi * dreps**2. * eff * dd)
        # propagate errors
        verr = numpy.sqrt(numpy.sum((4.*numpy.pi * dreps**2. * err * dd)**2.))

    return vol, verr


def compute_efficiency(f_dist, m_dist, dbins):
    '''
    Compute the efficiency as a function of distance for the given sets of found
    and missed injection distances.
    Note that injections that do not fit into any dbin get lost :(
    '''
    efficiency = numpy.zeros(len(dbins) - 1)
    error = numpy.zeros(len(dbins) - 1)
    for j, dlow in enumerate(dbins[:-1]):
        dhigh = dbins[j + 1]
        found = numpy.sum((dlow <= f_dist) * (f_dist < dhigh))
        missed = numpy.sum((dlow <= m_dist) * (m_dist < dhigh))
        if found+missed == 0:
            # avoid divide by 0 in empty bins
            missed = 1.
        efficiency[j] = float(found) / (found + missed)
        error[j] = numpy.sqrt(efficiency[j] * (1 - efficiency[j]) /
                              (found + missed))

    return efficiency, error


def mean_efficiency_volume(found, missed, dbins):

    if len(found) == 0:
        # no efficiency here
        return numpy.zeros(len(dbins) - 1), numpy.zeros(len(dbins) - 1), 0, 0

    # only need distances
    f_dist = numpy.array([l.distance for l in found])
    m_dist = numpy.array([l.distance for l in missed])

    # compute the efficiency and its variance
    eff, err = compute_efficiency(f_dist, m_dist, dbins)
    vol, verr = integrate_efficiency(dbins, eff, err)

    return eff, err, vol, verr


def filter_injections_by_mass(injs, mbins, bin_num, bin_type, bin_num2=None):
    '''
    For a given set of injections (sim_inspiral rows), return the subset
    of injections that fall within the given mass range.
    '''
    if bin_type == "Mass1_Mass2":
        m1bins = numpy.concatenate((mbins.lower()[0],
                                    numpy.array([mbins.upper()[0][-1]])))
        m1lo = m1bins[bin_num]
        m1hi = m1bins[bin_num + 1]
        m2bins = numpy.concatenate((mbins.lower()[1],
                                    numpy.array([mbins.upper()[1][-1]])))
        m2lo = m2bins[bin_num2]
        m2hi = m2bins[bin_num2 + 1]
        newinjs = [l for l in injs if
                   ((m1lo <= l.mass1 < m1hi and m2lo <= l.mass2 < m2hi) or
                    (m1lo <= l.mass2 < m1hi and m2lo <= l.mass1 < m2hi))]
        return newinjs

    mbins = numpy.concatenate((mbins.lower()[0],
                               numpy.array([mbins.upper()[0][-1]])))
    mlow = mbins[bin_num]
    mhigh = mbins[bin_num + 1]
    if bin_type == "Chirp_Mass":
        newinjs = [l for l in injs if (mlow <= l.mchirp < mhigh)]
    elif bin_type == "Total_Mass":
        newinjs = [l for l in injs if (mlow <= l.mass1 + l.mass2 < mhigh)]
    elif bin_type == "Component_Mass":
        # here it is assumed that m2 is fixed
        newinjs = [l for l in injs if (mlow <= l.mass1 < mhigh)]
    elif bin_type == "BNS_BBH":
        if bin_num in [0, 2]:
            # BNS/BBH case
            newinjs = [l for l in injs if
                       (mlow <= l.mass1 < mhigh and mlow <= l.mass2 < mhigh)]
        else:
            # NSBH
            newinjs = [l for l in injs if (mbins[0] <= l.mass1 < mbins[1] and
                                           mbins[2] <= l.mass2 < mbins[3])]
            # BHNS
            newinjs += [l for l in injs if (mbins[0] <= l.mass2 < mbins[1] and
                                            mbins[2] <= l.mass1 < mbins[3])]

    return newinjs


def compute_volume_vs_mass(found, missed, mass_bins, bin_type, dbins=None):
    """
    Compute the average luminosity an experiment was sensitive to

    Assumes that luminosity is uniformly distributed in space.
    Input is the sets of found and missed injections.
    """
    # mean and std estimate for luminosity
    volArray = bin_utils.BinnedArray(mass_bins)
    vol2Array = bin_utils.BinnedArray(mass_bins)

    # found/missed stats
    foundArray = bin_utils.BinnedArray(mass_bins)
    missedArray = bin_utils.BinnedArray(mass_bins)

    # compute the mean luminosity in each mass bin
    effvmass = []
    errvmass = []
    # 2D case first
    if bin_type == "Mass1_Mass2":
        for j, mc1 in enumerate(mass_bins.centres()[0]):
            for k, mc2 in enumerate(mass_bins.centres()[1]):
                newfound = filter_injections_by_mass(
                                              found, mass_bins, j, bin_type, k)
                newmissed = filter_injections_by_mass(
                                             missed, mass_bins, j, bin_type, k)

                foundArray[(mc1, mc2)] = len(newfound)
                missedArray[(mc1, mc2)] = len(newmissed)

                # compute the volume using this injection set
                meaneff, efferr, meanvol, volerr = mean_efficiency_volume(
                                                    newfound, newmissed, dbins)
                effvmass.append(meaneff)
                errvmass.append(efferr)
                volArray[(mc1, mc2)] = meanvol
                vol2Array[(mc1, mc2)] = volerr

        return volArray, vol2Array, foundArray, missedArray, effvmass, errvmass

    for j, mc in enumerate(mass_bins.centres()[0]):

        # filter out injections not in this mass bin
        newfound = filter_injections_by_mass(found, mass_bins, j, bin_type)
        newmissed = filter_injections_by_mass(missed, mass_bins, j, bin_type)

        foundArray[(mc, )] = len(newfound)
        missedArray[(mc, )] = len(newmissed)

        # compute the volume using this injection set
        meaneff, efferr, meanvol, volerr = mean_efficiency_volume(
                                                    newfound, newmissed, dbins)
        effvmass.append(meaneff)
        errvmass.append(efferr)
        volArray[(mc, )] = meanvol
        vol2Array[(mc, )] = volerr

    return volArray, vol2Array, foundArray, missedArray, effvmass, errvmass
