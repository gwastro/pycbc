""" Common utility functions for calculation of likelihoods
"""

import logging
import warnings
from distutils.util import strtobool

import numpy
import numpy.random
import tqdm

from scipy.special import logsumexp, i0e
from scipy.interpolate import RectBivariateSpline, interp1d
from pycbc.distributions import JointDistribution

from pycbc.detector import Detector


# Earth radius in seconds
EARTH_RADIUS = 0.031


def str_to_tuple(sval, ftype):
    """ Convenience parsing to convert str to tuple"""
    if sval is None:
        return ()
    return tuple(ftype(x.strip(' ')) for x in sval.split(','))


def str_to_bool(sval):
    """ Ensure value is a bool if it can be converted """
    if isinstance(sval, str):
        return strtobool(sval)
    return sval


def draw_sample(loglr, size=None):
    """ Draw a random index from a 1-d vector with loglr weights
    """
    if size:
        x = numpy.random.uniform(size=size)
    else:
        x = numpy.random.uniform()
    loglr = loglr - loglr.max()
    cdf = numpy.exp(loglr).cumsum()
    cdf /= cdf[-1]
    xl = numpy.searchsorted(cdf, x)
    return xl


class DistMarg():
    """Help class to add bookkeeping for distance marginalization"""

    marginalize_phase = None
    distance_marginalization = None
    distance_interpolator = None

    def setup_marginalization(self,
                              variable_params,
                              marginalize_phase=False,
                              marginalize_distance=False,
                              marginalize_distance_param='distance',
                              marginalize_distance_samples=int(1e4),
                              marginalize_distance_interpolator=False,
                              marginalize_distance_snr_range=None,
                              marginalize_distance_density=None,
                              marginalize_vector_params=None,
                              marginalize_vector_samples=1e3,
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
        def pop_prior(param):
            variable_params.remove(param)
            old_prior = kwargs['prior']

            dists = [d for d in old_prior.distributions
                     if param not in d.params]
            dprior = [d for d in old_prior.distributions
                      if param in d.params][0]
            prior = JointDistribution(variable_params,
                                      *dists, **old_prior.kwargs)
            kwargs['prior'] = prior
            return dprior

        self.reconstruct_phase = False
        self.reconstruct_distance = False
        self.reconstruct_vector = False

        # Handle any requested parameter vector / brute force marginalizations
        self.marginalize_vector_params = {}
        self.marginalized_vector_priors = {}
        self.vsamples = int(marginalize_vector_samples)

        for param in str_to_tuple(marginalize_vector_params, str):
            logging.info('Marginalizing over %s, %s points from prior',
                         param, self.vsamples)
            self.marginalized_vector_priors[param] = pop_prior(param)

        # Remove in the future, backwards compatibility
        if 'polarization_samples' in kwargs:
            warnings.warn("use marginalize_vector_samples rather "
                          "than 'polarization_samples'", DeprecationWarning)
            pol_uniform = numpy.linspace(0, numpy.pi * 2.0, self.vsamples)
            self.marginalize_vector_params['polarization'] = pol_uniform
            self.vsamples = int(kwargs['polarization_samples'])
            kwargs.pop('polarization_samples')

        self.reset_vector_params()

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
        dprior = pop_prior(marginalize_distance_param)

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

    def reset_vector_params(self):
        """ Redraw vector params from their priors
        """
        for param in self.marginalized_vector_priors:
            vprior = self.marginalized_vector_priors[param]
            values = vprior.rvs(self.vsamples)[param]
            self.marginalize_vector_params[param] = values

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
        return_complex = False
        distance = self.distance_marginalization

        if self.reconstruct_vector:
            skip_vector = True

        if self.reconstruct_distance:
            interpolator = None
            skip_vector = True

        if self.reconstruct_phase:
            interpolator = None
            distance = False
            skip_vector = True
            return_complex = True

        return marginalize_likelihood(sh_total, hh_total,
                                      logw=self.marginalize_vector_weights,
                                      phase=self.marginalize_phase,
                                      interpolator=interpolator,
                                      distance=distance,
                                      skip_vector=skip_vector,
                                      return_complex=return_complex,
                                      return_peak=return_peak)

    def snr_draw(self, snrs):
        """ Improve the monte-carlo vector marginalization using the SNR time
        series of each detector
        """
        p = self.current_params

        if (not numpy.isscalar(p['tc']) and
            'tc' in self.marginalized_vector_priors and
            not ('ra' in self.marginalized_vector_priors
                 or 'dec' in self.marginalized_vector_priors)):
            self.draw_times(snrs)
        elif (not numpy.isscalar(p['tc']) and
              'tc' in self.marginalized_vector_priors and
              'ra' in self.marginalized_vector_priors and
              'dec' in self.marginalized_vector_priors):
            self.draw_sky_times(snrs)
        else:
            # OK, we couldn't do anything with the requested monte-carlo
            # marginalizations.
            pass

    def draw_times(self, snrs):
        """ Draw times consistent with the incoherent network SNR

        Parameters
        ----------
        snrs: dist of TimeSeries
        """
        if not hasattr(self, 'tinfo'):
            # determine the rough time offsets for this sky location
            tcmin, tcmax = self.marginalized_vector_priors['tc'].bounds['tc']
            tcave = (tcmax + tcmin) / 2.0
            ifos = list(snrs.keys())
            if hasattr(self, 'keep_ifos'):
                ifos = self.keep_ifos
            d = {ifo: Detector(ifo, reference_time=tcave) for ifo in ifos}
            self.tinfo = tcmin, tcmax, tcave, ifos, d

        tcmin, tcmax, tcave, ifos, d = self.tinfo

        # Determine the weights for the valid time range
        ra = self._current_params['ra']
        dec = self._current_params['dec']

        # Determine the common valid time range
        iref = ifos[0]
        dref = d[iref]
        dt = dref.time_delay_from_earth_center(ra, dec, tcave)

        starts = []
        ends = []

        tmin, tmax = tcmin - dt, tcmax + dt
        if hasattr(self, 'tstart'):
            tmin = self.tstart[iref]
            tmax = self.tend[iref]

        starts.append(max(tmin, snrs[iref].start_time))
        ends.append(min(tmax, snrs[iref].end_time))

        idels = {}
        for ifo in ifos[1:]:
            dti = d[ifo].time_delay_from_detector(dref, ra, dec, tcave)
            idel = round(dti / snrs[iref].delta_t) * snrs[iref].delta_t
            idels[ifo] = idel

            starts.append(snrs[ifo].start_time - idel)
            ends.append(snrs[ifo].end_time - idel)
        start = max(starts)
        end = min(ends)
        if end <= start:
            return

        # get the weights
        snr = snrs[iref].time_slice(start, end, mode='nearest')
        logweight = snr.squared_norm().numpy()
        for ifo in ifos[1:]:
            idel = idels[ifo]
            snrv = snrs[ifo].time_slice(snr.start_time + idel,
                                        snr.end_time + idel,
                                        mode='nearest')
            logweight += snrv.squared_norm().numpy()
        logweight /= 2.0

        # Draw proportional to the incoherent likelihood
        # Draw first which time sample
        tci = draw_sample(logweight, size=self.vsamples)
        # Second draw a subsample size offset so that all times are covered
        tct = numpy.random.uniform(-snr.delta_t / 2.0,
                                   snr.delta_t / 2.0,
                                   size=self.vsamples)
        tc = tct + tci * snr.delta_t + float(snr.start_time) - dt

        # Update the current proposed times and the marginalization values
        self.marginalize_vector_params['tc'] = tc
        self._current_params['tc'] = tc

        # Update the importance weights for each vector sample
        logw = self.marginalize_vector_weights - logweight[tci]
        self.marginalize_vector_weights = logw - logsumexp(logw)

    def draw_sky_times(self, snrs):
        """ Draw ra, dec, and tc together using SNR timeseries to determine
        monte-carlo weights.
        """
        # First setup
        # precalculate dense sky grid and make dict and or array of the results
        ifos = list(snrs.keys())
        if hasattr(self, 'keep_ifos'):
            ifos = self.keep_ifos
        ikey = ''.join(ifos)

        # No good SNR peaks, go with prior draw
        if len(ifos) == 0:
            return

        def make_init():
            logging.info('pregenerating sky pointings')
            size = int(1e6)
            logging.info('drawing samples')
            ra = self.marginalized_vector_priors['ra'].rvs(size=size)['ra']
            dec = self.marginalized_vector_priors['dec'].rvs(size=size)['dec']
            tcmin, tcmax = self.marginalized_vector_priors['tc'].bounds['tc']
            tcave = (tcmax + tcmin) / 2.0
            d = {ifo: Detector(ifo, reference_time=tcave) for ifo in ifos}

            # What data structure to hold times? Dict of offset -> list?
            logging.info('sorting into time delay dict')
            dts = []
            for i in range(len(ifos) - 1):
                dt = d[ifos[0]].time_delay_from_detector(d[ifos[i+1]],
                                                         ra, dec, tcave)
                dt = numpy.rint(dt / snrs[ifos[0]].delta_t)
                dts.append(dt)

            dtc = d[ifos[0]].time_delay_from_earth_center(ra, dec, tcave)

            dmap = {}
            for i, t in enumerate(tqdm.tqdm(zip(*dts))):
                v = ra[i], dec[i], dtc[i]
                if t not in dmap:
                    dmap[t] = []
                dmap[t].append(v)

            if len(ifos) == 1:
                dmap[()] = list(zip(ra, dec, dtc))
            return ifos, dmap, tcmin, tcmax

        if not hasattr(self, 'tinfo'):
            self.tinfo = {}

        if ikey not in self.tinfo:
            self.tinfo[ikey] = make_init()

        ifos, dmap, tcmin, tcmax = self.tinfo[ikey]

        # draw times from each snr time series
        # Is it worth doing this if some detector has low SNR?
        sref = None
        iref = None
        idx = []
        dx = []
        mcweight = None
        for ifo in ifos:
            snr = snrs[ifo]
            tmin, tmax = tcmin - EARTH_RADIUS, tcmax + EARTH_RADIUS
            if hasattr(self, 'tstart'):
                tmin = self.tstart[ifo]
                tmax = self.tend[ifo]

            start = max(tmin, snrs[ifo].start_time)
            end = min(tmax, snrs[ifo].end_time)

            snr = snr.time_slice(start, end, mode='nearest')

            w = snr.squared_norm().numpy() / 2.0
            i = draw_sample(w, size=self.vsamples)

            if sref is not None:
                mcweight -= w[i]
                delt = float(snr.start_time - sref.start_time)
                i += round(delt / sref.delta_t)
                dx.append(iref - i)
            else:
                sref = snr
                iref = i
                mcweight = -w[i]

            idx.append(i)

        # check if delay is in dict, if not, throw out
        ra = []
        dec = []
        dtc = []
        ti = []
        wi = []
        rand = numpy.random.uniform(0, 1, size=self.vsamples)
        for i in range(self.vsamples):
            t = tuple(x[i] for x in dx)
            if t in dmap:
                randi = int(rand[i] * (len(dmap[t])))
                xra, xdec, xdtc = dmap[t][randi]
                ra.append(xra)
                dec.append(xdec)
                dtc.append(xdtc)
                wi.append(len(dmap[t]))
                ti.append(i)

        # If we had really poor efficiency at finding a point, we should
        # give up and just use the original random draws
        if len(ra) < 0.05 * self.vsamples:
            return

        # fill back to fixed size with repeat samples
        # sample order is random, so this should be OK statistically

        ra = numpy.resize(numpy.array(ra), self.vsamples)
        dec = numpy.resize(numpy.array(dec), self.vsamples)
        dtc = numpy.resize(numpy.array(dtc), self.vsamples)
        ti = numpy.resize(numpy.array(ti, dtype=int), self.vsamples)
        wi = numpy.resize(numpy.array(wi), self.vsamples)

        # Second draw a subsample size offset so that all times are covered
        tct = numpy.random.uniform(-snr.delta_t / 2.0,
                                   snr.delta_t / 2.0,
                                   size=len(ti))

        tc = tct + iref[ti] * snr.delta_t + float(sref.start_time) - dtc

        # Update the current proposed times and the marginalization values
        self.marginalize_vector_params['tc'] = tc
        self.marginalize_vector_params['ra'] = ra
        self.marginalize_vector_params['dec'] = dec
        self._current_params.update(self.marginalize_vector_params)

        # Update the importance weights for each vector sample
        logw = self.marginalize_vector_weights + mcweight[ti] + numpy.log(wi)
        self.marginalize_vector_weights = logw - logsumexp(logw)

    def setup_peak_lock(self,
                        sample_rate=4096,
                        snrs=None,
                        peak_lock_snr=None,
                        peak_lock_ratio=1e4,
                        peak_lock_region=4,
                        **kwargs):
        """ Determine where to constrain marginalization based on
        the observed reference SNR peaks.

        Parameters
        ----------
        sample_rate : float
            The SNR sample rate
        snrs : Dict of SNR time series
            Either provide this or the model needs a function
            to get the reference SNRs.
        peak_lock_snr: float
            The minimum SNR to bother restricting from the prior range
        peak_lock_ratio: float
            The likelihood ratio (not log) relative to the peak to
            act as a threshold bounding region.
        peak_lock_region: int
            Number of samples to inclue beyond the strict region
            determined by the relative likelihood
        """

        if 'tc' not in self.marginalized_vector_priors:
            return

        tcmin, tcmax = self.marginalized_vector_priors['tc'].bounds['tc']
        tstart = tcmin - EARTH_RADIUS
        tmax = tcmax - tcmin + EARTH_RADIUS * 2.0
        num_samples = int(tmax * sample_rate)
        self.tstart = {ifo: tstart for ifo in self.data}
        self.num_samples = {ifo: num_samples for ifo in self.data}

        if snrs is None:
            if not hasattr(self, 'ref_snr'):
                raise ValueError("Model didn't have a reference SNR!")
            snrs = self.ref_snr

        # Restrict the time range for constructing SNR time series
        # to identifiable peaks
        if peak_lock_snr is not None:
            peak_lock_snr = float(peak_lock_snr)
            peak_lock_ratio = float(peak_lock_ratio)
            peak_lock_region = int(peak_lock_region)

            for ifo in snrs:
                z = snrs[ifo].time_slice(tstart, tstart + tmax, mode='nearest')
                peak_snr, imax = z.abs_max_loc()
                times = z.sample_times
                peak_time = times[imax]

                logging.info('%s: Max Ref SNR Peak of %s at %s',
                             ifo, peak_snr, peak_time)

                if peak_snr > peak_lock_snr:
                    target = peak_snr ** 2.0 / 2.0 - numpy.log(peak_lock_ratio)
                    target = (target * 2.0) ** 0.5

                    region = numpy.where(abs(z) > target)[0]
                    ts = times[region[0]] - peak_lock_region / sample_rate
                    te = times[region[-1]] + peak_lock_region / sample_rate
                    self.tstart[ifo] = ts
                    self.num_samples[ifo] = int((te - ts) * sample_rate)

            # Check times are commensurate with each other
            for ifo in snrs:
                ts = self.tstart[ifo]
                te = ts + self.num_samples[ifo] / sample_rate

                for ifo2 in snrs:
                    if ifo == ifo2:
                        continue
                    ts2 = self.tstart[ifo2]
                    te2 = ts2 + self.num_samples[ifo2] / sample_rate
                    det = Detector(ifo)
                    dt = Detector(ifo2).light_travel_time_to_detector(det)

                    ts = max(ts, ts2 - dt)
                    te = min(te, te2 + dt)

                self.tstart[ifo] = ts
                self.num_samples[ifo] = int((te - ts) * sample_rate) + 1
                logging.info('%s: use region %s-%s, %s points',
                             ifo, ts, te, self.num_samples[ifo])

        self.tend = self.tstart.copy()
        for ifo in snrs:
            self.tend[ifo] += self.num_samples[ifo] / sample_rate

    def draw_ifos(self, snrs, peak_snr_threshold=4.0, log=True, **kwargs):
        """ Helper utility to determine which ifos we should use based on the
        reference SNR time series.
        """
        if 'tc' not in self.marginalized_vector_priors:
            return

        tcmin, tcmax = self.marginalized_vector_priors['tc'].bounds['tc']
        ifos = list(snrs.keys())
        keep_ifos = []
        psnrs = []
        for ifo in ifos:
            snr = snrs[ifo]
            start = max(tcmin - EARTH_RADIUS, snr.start_time)
            end = min(tcmax + EARTH_RADIUS, snr.end_time)
            snr = snr.time_slice(start, end, mode='nearest')
            psnr = abs(snr).max()
            if psnr > peak_snr_threshold:
                keep_ifos.append(ifo)
            psnrs.append(psnr)

        if log:
            logging.info("Ifos used for SNR based draws:"
                         " %s, snrs: %s, peak_snr_threshold=%s",
                         keep_ifos, psnrs, peak_snr_threshold)

        self.keep_ifos = keep_ifos
        return keep_ifos

    @property
    def current_params(self):
        """ The current parameters

        If a parameter has been vector marginalized, the likelihood should
        expect an array for the given parameter. This allows transparent
        vectorization for many models.
        """
        params = self._current_params
        for k in self.marginalize_vector_params:
            if k not in params:
                params[k] = self.marginalize_vector_params[k]
        self.marginalize_vector_weights = - numpy.log(self.vsamples)
        return params

    def reconstruct(self, rec=None, seed=None):
        """ Reconstruct the distance or vectored marginalized parameter
        of this class.
        """
        if seed:
            numpy.random.seed(seed)

        if rec is None:
            rec = {}

        def get_loglr():
            p = self.current_params.copy()
            p.update(rec)
            self.update(**p)
            return self.loglr

        if self.marginalize_vector_params:
            logging.info('Reconstruct vector')
            self.reconstruct_vector = True
            self.reset_vector_params()
            loglr = get_loglr()
            xl = draw_sample(loglr + self.marginalize_vector_weights)
            for k in self.marginalize_vector_params:
                rec[k] = self.marginalize_vector_params[k][xl]
            self.reconstruct_vector = False

        if self.distance_marginalization:
            logging.info('Reconstruct distance')
            # call likelihood to get vector output
            self.reconstruct_distance = True
            _, weights = self.distance_marginalization
            loglr = get_loglr()
            xl = draw_sample(loglr + numpy.log(weights))
            rec['distance'] = self.dist_locs[xl]
            self.reconstruct_distance = False

        if self.marginalize_phase:
            logging.info('Reconstruct phase')
            self.reconstruct_phase = True
            s, h = get_loglr()
            phasev = numpy.linspace(0, numpy.pi*2.0, int(1e4))
            # This assumes that the template was conjugated in inner products
            loglr = (numpy.exp(-2.0j * phasev) * s).real + h
            xl = draw_sample(loglr)
            rec['coa_phase'] = phasev[xl]
            self.reconstruct_phase = False

        rec['loglr'] = loglr[xl]
        rec['loglikelihood'] = self.lognl + rec['loglr']
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
                           logw=None,
                           phase=False,
                           distance=False,
                           skip_vector=False,
                           interpolator=None,
                           return_peak=False,
                           return_complex=False,
                           ):
    """ Return the marginalized likelihood.

    Apply various marginalizations to the data, including phase, distance,
    and brute-force vector marginalizations. Several options relate
    to how the distance marginalization is approximated and others allow for
    special return products to aid in parameter reconstruction.

    Parameters
    ----------
    sh: complex float or numpy.ndarray
        The data-template inner product
    hh: complex float or numpy.ndarray
        The template-template inner product
    logw:
        log weighting factors if vector marginalization is used, if not
        given, each sample is assumed to be equally weighted
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
    return_peak: bool, False
        Return the peak likelihood and index if using passing an array as
        input in addition to the marginalized over the array likelihood.
    return_complex: bool, False
        Return the sh / hh data products before applying phase marginalization.
        This option is intended to aid in reconstucting phase marginalization
        and is unlikely to be useful for other purposes.

    Returns
    -------
    loglr: float
        The marginalized loglikehood ratio
    """
    if distance and not interpolator and not numpy.isscalar(sh):
        raise ValueError("Cannot do vector marginalization "
                         "and distance at the same time")

    if logw is None:
        if isinstance(hh, float):
            logw = 0
        else:
            logw = -numpy.log(len(sh))

    if return_complex:
        pass
    elif phase:
        sh = abs(sh)
    else:
        sh = sh.real

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
            sh = sh * dist_rescale
            hh = hh * dist_rescale ** 2.0
            logw = numpy.log(dist_weights)

        if return_complex:
            return sh, -0.5 * hh

        # Apply the phase marginalization
        if phase:
            sh = numpy.log(i0e(sh)) + sh

        # Calculate loglikelihood ratio
        vloglr = sh - 0.5 * hh

    if return_peak:
        maxv = vloglr.argmax()
        maxl = vloglr[maxv]

    # Do brute-force marginalization if loglr is a vector
    if isinstance(vloglr, float):
        vloglr = float(vloglr)
    elif not skip_vector:
        vloglr = float(logsumexp(vloglr, b=numpy.exp(logw)))

    if return_peak:
        return vloglr, maxv, maxl
    return vloglr
