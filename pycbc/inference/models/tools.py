""" Common utility functions for calculation of likelihoods
"""

import logging
from distutils.util import strtobool

import numpy
import numpy.random
import tqdm

from scipy.special import logsumexp, i0e
from scipy.interpolate import RectBivariateSpline, interp1d
from pycbc.distributions import JointDistribution


def str_to_tuple(sval, ftype):
    """ Convenience parsing to convert str to tuple"""
    return tuple(ftype(x) for x in sval.split(','))


def str_to_bool(sval):
    """ Ensure value is a bool if it can be converted """
    if isinstance(sval, str):
        return strtobool(sval)
    return sval


class DistMarg():
    """Help class to add bookkeeping for distance marginalization"""

    marginalize_phase = None
    distance_marginalization = None
    distance_interpolator = None

    def setup_distance_marginalization(self,
                                       variable_params,
                                       marginalize_phase=False,
                                       marginalize_distance=False,
                                       marginalize_distance_param='distance',
                                       marginalize_distance_samples=int(1e4),
                                       marginalize_distance_interpolator=False,
                                       marginalize_distance_snr_range=None,
                                       marginalize_distance_density=None,
                                       marginalize_vector=False,
                                       marginalize_vector_params=None,
                                       **kwargs):
        """ Setup the model for use with distance marginalization

        This function sets up precalculations for distance / phase
        marginalization. For distance margininalization it modifies the
        model to internally remove distance as a parameter.

        Parameters
        ----------
        variable_params: list of strings
            The set of variable parameters
        marginalize_phase: bool, False
            Do analytic marginalization (appopriate only for 22 mode waveforms)
        marginalize_distance: bool, False
            Marginalize over distance
        marginalize_distance_param: str
            Name of the parameter that is used to determine the distance.
            This might be 'distance' or a parameter which can be converted
            to distance by a provided univariate transformation.
        marginalize_distance_interpolator: bool
            Use a pre-calculated interpolating function for the distance
            marginalized likelihood.
        marginalize_distance_snr_range: tuple of floats, (1, 50)
            The SNR range for the interpolating function to be defined in.
            If a sampler goes outside this range, the logl will be returned
            as -numpy.inf.
        marginalize_distance_density: tuple of intes, (1000, 1000)
            The dimensions of the interpolation grid over (sh, hh).

        Returns
        -------
        variable_params: list of strings
            Set of variable params (missing distance-related parameter).
        kwags: dict
            The keyword arguments to the model initialization, may be modified
            from the original set by this function.
        """
        self.marginalize_vector = marginalize_vector
        self.marginalize_vector_params = marginalize_vector_params

        self.reconstructing_distance = False
        self.marginalize_phase = str_to_bool(marginalize_phase)
        self.distance_marginalization = False
        self.distance_interpolator = None

        marginalize_distance = str_to_bool(marginalize_distance)
        self.marginalize_distance = marginalize_distance
        if not marginalize_distance:
            return variable_params, kwargs

        if isinstance(marginalize_distance_snr_range, str):
            marginalize_distance_snr_range = \
                str_to_tuple(marginalize_distance_snr_range, float)

        if isinstance(marginalize_distance_density, str):
            marginalize_distance_density = \
                str_to_tuple(marginalize_distance_density, int)

        logging.info('Marginalizing over distance')

        # Take distance out of the variable params since we'll handle it
        # manually now
        variable_params.remove(marginalize_distance_param)
        old_prior = kwargs['prior']
        dists = [d for d in old_prior.distributions
                 if marginalize_distance_param not in d.params]
        dprior = [d for d in old_prior.distributions
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
        elif marginalize_distance_param != 'distance':
            waveform_transforms = kwargs['waveform_transforms']
            pname = dprior.params[0]
            logging.info("Settings up transform,  prior is in terms of"
                         " %s", pname)
            wtrans = [d for d in waveform_transforms
                      if 'distance' not in d.outputs]
            if len(wtrans) == 0:
                wtrans = None
            kwargs['waveform_transforms'] = wtrans
            dtrans = [d for d in waveform_transforms
                      if 'distance' in d.outputs][0]
            v = dprior.rvs(int(1e8))
            d = dtrans.transform({pname: v[pname]})['distance']
            d.sort()
            cdf = numpy.arange(1, len(d)+1) / len(d)
            i = interp1d(d, cdf)
            dmin, dmax = d.min(), d.max()
            logging.info('Distance range %s-%s', dmin, dmax)
            x = numpy.linspace(dmin, dmax,
                               int(marginalize_distance_samples) + 1)
            xl, xr = x[:-1], x[1:]
            dist_locs = 0.5 * (xr + xl)
            dist_weights = i(xr) - i(xl)
        else:
            raise ValueError("No prior seems to determine the distance")

        dist_weights /= dist_weights.sum()
        dist_ref = 0.5 * (dmax + dmin)
        self.dist_locs = dist_locs
        self.distance_marginalization = dist_ref / dist_locs, dist_weights
        self.distance_interpolator = None
        if str_to_bool(marginalize_distance_interpolator):
            setup_args = {}
            if marginalize_distance_snr_range:
                setup_args['snr_range'] = marginalize_distance_snr_range
            if marginalize_distance_density:
                setup_args['density'] = marginalize_distance_density
            i = setup_distance_marg_interpolant(self.distance_marginalization,
                                                phase=self.marginalize_phase,
                                                **setup_args)
            self.distance_interpolator = i
        kwargs['static_params']['distance'] = dist_ref
        return variable_params, kwargs

    def marginalize_loglr(self, sh_total, hh_total,
                          skip_vector=False, return_peak=False):
        """ Return the marginal likelihood

        Parameters
        -----------
        sh_total: float or ndarray
            The total <s|h> inner product summed over detectors
        hh_total: float or ndarray
            The total <h|h> inner product summed over detectors
        skip_vector: bool, False
            If true, and input is a vector, do not marginalize over that
            vector, instead return the likelihood values as a vector.
        """
        interpolator = self.distance_interpolator
        if self.reconstructing_distance:
            skip_vector = True
            interpolator = None

        return marginalize_likelihood(sh_total, hh_total,
                                      phase=self.marginalize_phase,
                                      interpolator=interpolator,
                                      distance=self.distance_marginalization,
                                      skip_vector=skip_vector,
                                      return_peak=return_peak)

    def reconstruct(self):
        """ Reconstruct the distance or vectored marginalized parameter
        of this class.
        """
        rec = {}

        def draw_sample(params, loglr):
            x = numpy.random.uniform()
            cdf = numpy.exp(loglr).cumsum()
            cdf /= cdf[-1]
            xl = numpy.searchsorted(cdf, x)
            return params[xl], xl

        if self.distance_marginalization:
            # turn off interpolator and set to return vector output
            self.reconstructing_distance = True

        loglr_full = self.loglr
        if self.marginalize_vector and self.distance_marginalization:
            loglr = logsumexp(loglr_full, axis=0)
        else:
            loglr = loglr_full

        if self.distance_marginalization:
            # call likelihood to get vector output
            _, weights = self.distance_marginalization
            loglr += numpy.log(weights)

            # draw distance sample
            dist, xl = draw_sample(self.dist_locs, loglr)
            rec['distance'] = dist

        if self.marginalize_vector:
            if self.distance_marginalization:
                # logl along for the selected distance
                vlr = loglr_full[:, xl]
            else:
                vlr = loglr

            vec_param, _ = draw_sample(self.marginalize_vector_params, vlr)
            rec[self.marginalize_vector] = vec_param

        self.reconstructing_distance = False
        return rec


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

    Returns
    -------
    interp: function
        Function which returns the precalculated likelihood for a given
        inner product sh/hh.
    """
    dist_rescale, _ = dist_marg
    logging.info("Interpolator valid for SNRs in %s", snr_range)
    logging.info("Interpolator using grid %s", density)
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

    def interp_wrapper(x, y, bounds_check=True):
        k = None
        if bounds_check:
            if isinstance(x, float):
                if x > shr_max or x < shr_min or y > hhr_max or y < hhr_min:
                    return -numpy.inf
            else:
                k = (x > shr_max) | (x < shr_min)
                k = k | (y > hhr_max) | (y < hhr_min)

        v = interp(x, y, grid=False)
        if k is not None:
            v[k] = -numpy.inf
        return v
    return interp_wrapper


def marginalize_likelihood(sh, hh,
                           phase=False,
                           distance=False,
                           skip_vector=False,
                           interpolator=None,
                           return_peak=False,
                           ):
    """ Return the marginalized likelihood

    Parameters
    ----------
    sh: complex float or numpy.ndarray
        The data-template inner product
    hh: complex float or numpy.ndarray
        The template-template inner product
    phase: bool, False
        Enable phase marginalization. Only use if orbital phase can be related
        to just a single overall phase (e.g. not true for waveform with
        sub-dominant modes)
    skip_vector: bool, False
        Don't apply marginalization of vector component of input (i.e. leave
        as vector).
    interpolator: function, None
        If provided, internal calculation is skipped in favor of a
        precalculated interpolating function which takes in sh/hh
        and returns the likelihood.

    Returns
    -------
    loglr: float
        The marginalized loglikehood ratio
    """
    if isinstance(hh, float):
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

        if skip_vector:
            return vloglr
    else:
        # explicit calculation
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

    if return_peak:
        maxv = vloglr.argmax()
        maxl = vloglr[maxv]

    # Do brute-force marginalization if loglr is a vector
    if isinstance(vloglr, float):
        vloglr = float(vloglr)
    elif not skip_vector:
        vloglr = float(logsumexp(vloglr, b=vweights)) - clogweights

    if return_peak:
        return vloglr, maxv, maxl
    return vloglr
