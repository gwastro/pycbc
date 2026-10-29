# Copyright (C) 2013 Ian W. Harry
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

import logging

import numpy

from pycbc.tmpltbank.lambda_mapping import generate_mapping

logger = logging.getLogger(__name__)


def determine_eigen_directions(
    metric_params, preserve_moments=False, vary_fmax=False, vary_density=None
):
    """
    This function will calculate the coordinate transfomations that are needed
    to rotate from a coordinate system described by the various Lambda
    components in the frequency expansion, to a coordinate system where the
    metric is Cartesian.

    Parameters
    ----------
    metric_params : MetricParameters instance
        Structure holding all the options for construction of the metric.
    preserve_moments : boolean, optional (default False)
        Currently only used for debugging.
        If this is given then if the moments structure is already set
        within metric_params then they will not be recalculated.
    vary_fmax : boolean, optional (default False)
        If set to False the metric and rotations are calculated once, for the
        full range of frequency [f_low,f_upper).
        If set to True the metric and rotations are calculated multiple times,
        for frequency ranges [f_low,f_low + i*vary_density), where i starts at
        1 and runs up until f_low + (i+1)*vary_density > f_upper.
        Thus values greater than f_upper are *not* computed.
        The calculation for the full range [f_low,f_upper) is also done.
    vary_density : float, optional
        If vary_fmax is True, this will be used in computing the frequency
        ranges as described for vary_fmax.

    Returns
    -------
    metric_params : MetricParameters instance
        Structure holding all the options for construction of the metric.
        **THIS FUNCTION ONLY RETURNS THE CLASS**
        The following will be **added** to this structure
    metric_params.evals : Dictionary of numpy.array
        Each entry in the dictionary corresponds to the different frequency
        ranges described in vary_fmax. If vary_fmax = False, the only entry
        will be f_upper, this corresponds to integrals in [f_low,f_upper). This
        entry is always present. Each other entry will use floats as keys to
        the dictionary. These floats give the upper frequency cutoff when it is
        varying.
        Each numpy.array contains the eigenvalues which, with the eigenvectors
        in evecs, are needed to rotate the
        coordinate system to one in which the metric is the identity matrix.
    metric_params.evecs : Dictionary of numpy.matrix
        Each entry in the dictionary is as described under evals.
        Each numpy.matrix contains the eigenvectors which, with the eigenvalues
        in evals, are needed to rotate the
        coordinate system to one in which the metric is the identity matrix.
    metric_params.metric : Dictionary of numpy.matrix
        Each entry in the dictionary is as described under evals.
        Each numpy.matrix contains the metric of the parameter space in the
        Lambda_i coordinate system.
    metric_params.moments : Moments structure
        See the structure documentation for a description of this. This
        contains the result of all the integrals used in computing the metrics
        above. It can be used for additional metric calculations.

    """
    evals = {}
    evecs = {}
    metric = {}
    unmax_metric = {}

    # First step is to get the moments needed to calculate the metric
    if not (metric_params.moments and preserve_moments):
        get_moments(metric_params, vary_fmax=vary_fmax, vary_density=vary_density)

    # What values are going to be in the moments
    # J7 is the normalization factor so it *MUST* be present
    list = metric_params.moments["J7"].keys()

    # We start looping over every item in the list of metrics
    for item in list:
        # Here we convert the moments into a form easier to use here
        Js = {}  # noqa: N806
        for i in range(-7, 18):
            Js[i] = metric_params.moments[f"J{i}"][item]

        log_js = {}
        for i in range(-1, 18):
            log_js[i] = metric_params.moments[f"log{i}"][item]

        loglog_js = {}
        for i in range(-1, 18):
            loglog_js[i] = metric_params.moments[f"loglog{i}"][item]

        logloglog_js = {}
        for i in range(-1, 18):
            logloglog_js[i] = metric_params.moments[f"logloglog{i}"][item]

        loglogloglog_js = {}
        for i in range(-1, 18):
            loglogloglog_js[i] = metric_params.moments[f"loglogloglog{i}"][item]

        mapping = generate_mapping(metric_params.pn_order)

        # Calculate the metric
        gs, unmax_metric_curr = calculate_metric(
            Js, log_js, loglog_js, logloglog_js, loglogloglog_js, mapping
        )
        metric[item] = gs
        unmax_metric[item] = unmax_metric_curr

        # And the eigenvalues
        evals[item], evecs[item] = numpy.linalg.eig(gs)

        # Numerical error can lead to small negative eigenvalues.
        for i in range(len(evals[item])):
            if evals[item][i] < 0:
                # Due to numerical imprecision the very small eigenvalues can
                # be negative. Make these positive.
                evals[item][i] = -evals[item][i]
            if evecs[item][i, i] < 0:
                # We demand a convention that all diagonal terms in the matrix
                # of eigenvalues are positive.
                # This is done to help visualization of the spaces (increasing
                # mchirp always goes the same way)
                evecs[item][:, i] = -evecs[item][:, i]

    metric_params.evals = evals
    metric_params.evecs = evecs
    metric_params.metric = metric
    metric_params.time_unprojected_metric = unmax_metric

    return metric_params


def get_moments(metric_params, vary_fmax=False, vary_density=None):
    """
    This function will calculate the various integrals (moments) that are
    needed to compute the metric used in template bank placement and
    coincidence.

    Parameters
    ----------
    metric_params : MetricParameters instance
        Structure holding all the options for construction of the metric.
    vary_fmax : boolean, optional (default False)
        If set to False the metric and rotations are calculated once, for the
        full range of frequency [f_low,f_upper).
        If set to True the metric and rotations are calculated multiple times,
        for frequency ranges [f_low,f_low + i*vary_density), where i starts at
        1 and runs up until f_low + (i+1)*vary_density > f_upper.
        Thus values greater than f_upper are *not* computed.
        The calculation for the full range [f_low,f_upper) is also done.
    vary_density : float, optional
        If vary_fmax is True, this will be used in computing the frequency
        ranges as described for vary_fmax.

    Returns
    -------
    None : None
        **THIS FUNCTION RETURNS NOTHING**
        The following will be **added** to the metric_params structure
    metric_params.moments : Moments structure
        This contains the result of all the integrals used in computing the
        metrics above. It can be used for additional metric calculations.
        This is composed of two compound
        dictionaries. The first entry indicates which moment is being
        calculated and the second entry indicates the upper frequency cutoff
        that was used.

        In all cases x = f/f0.

        For the first entries the options are:

        moments['J%d' %(i)][f_cutoff]
        This stores the integral of
        x**((-i)/3.) * delta X / PSD(x)

        moments['log%d' %(i)][f_cutoff]
        This stores the integral of
        (numpy.log(x**(1./3.))) x**((-i)/3.) * delta X / PSD(x)

        moments['loglog%d' %(i)][f_cutoff]
        This stores the integral of
        (numpy.log(x**(1./3.)))**2 x**((-i)/3.) * delta X / PSD(x)

        moments['loglog%d' %(i)][f_cutoff]
        This stores the integral of
        (numpy.log(x**(1./3.)))**3 x**((-i)/3.) * delta X / PSD(x)

        moments['loglog%d' %(i)][f_cutoff]
        This stores the integral of
        (numpy.log(x**(1./3.)))**4 x**((-i)/3.) * delta X / PSD(x)

        The second entry stores the frequency cutoff used when computing
        the integral. See description of the vary_fmax option above.

        All of these values are nomralized by a factor of

        x**((-7)/3.) * delta X / PSD(x)

        The normalization factor can be obtained in

        moments['I7'][f_cutoff]

    """
    # NOTE: Unless the TaylorR2F4 metric is used the log^3 and log^4 terms are
    # not needed. As this calculation is not too slow compared to bank
    # placement we just do this anyway.

    psd_amp = metric_params.psd.data
    psd_f = numpy.arange(len(psd_amp), dtype=float) * metric_params.delta_f
    new_f, new_amp = interpolate_psd(psd_f, psd_amp, metric_params.delta_f)

    # Need I7 first as this is the normalization factor
    funct = lambda x, f0: 1
    I7 = calculate_moment(  # noqa: N806
        new_f,
        new_amp,
        metric_params.f_low,
        metric_params.f_upper,
        metric_params.f0,
        funct,
        vary_fmax=vary_fmax,
        vary_density=vary_density,
    )

    # Do all the J moments
    moments = {}
    moments["I7"] = I7
    for i in range(-7, 18):
        funct = lambda x, f0: x ** ((-i + 7) / 3.0)
        moments[f"J{i}"] = calculate_moment(
            new_f,
            new_amp,
            metric_params.f_low,
            metric_params.f_upper,
            metric_params.f0,
            funct,
            norm=I7,
            vary_fmax=vary_fmax,
            vary_density=vary_density,
        )

    # Do the logx multiplied by some power terms
    for i in range(-1, 18):
        funct = lambda x, f0: (
            (numpy.log((x * f0) ** (1.0 / 3.0))) * x ** ((-i + 7) / 3.0)
        )
        moments[f"log{i}"] = calculate_moment(
            new_f,
            new_amp,
            metric_params.f_low,
            metric_params.f_upper,
            metric_params.f0,
            funct,
            norm=I7,
            vary_fmax=vary_fmax,
            vary_density=vary_density,
        )

    # Do the loglog term
    for i in range(-1, 18):
        funct = lambda x, f0: (
            (numpy.log((x * f0) ** (1.0 / 3.0))) ** 2 * x ** ((-i + 7) / 3.0)
        )
        moments[f"loglog{i}"] = calculate_moment(
            new_f,
            new_amp,
            metric_params.f_low,
            metric_params.f_upper,
            metric_params.f0,
            funct,
            norm=I7,
            vary_fmax=vary_fmax,
            vary_density=vary_density,
        )

    # Do the logloglog term
    for i in range(-1, 18):
        funct = lambda x, f0: (
            (numpy.log((x * f0) ** (1.0 / 3.0))) ** 3 * x ** ((-i + 7) / 3.0)
        )
        moments[f"logloglog{i}"] = calculate_moment(
            new_f,
            new_amp,
            metric_params.f_low,
            metric_params.f_upper,
            metric_params.f0,
            funct,
            norm=I7,
            vary_fmax=vary_fmax,
            vary_density=vary_density,
        )

    # Do the logloglog term
    for i in range(-1, 18):
        funct = lambda x, f0: (
            (numpy.log((x * f0) ** (1.0 / 3.0))) ** 4 * x ** ((-i + 7) / 3.0)
        )
        moments[f"loglogloglog{i}"] = calculate_moment(
            new_f,
            new_amp,
            metric_params.f_low,
            metric_params.f_upper,
            metric_params.f0,
            funct,
            norm=I7,
            vary_fmax=vary_fmax,
            vary_density=vary_density,
        )

    metric_params.moments = moments


def interpolate_psd(psd_f, psd_amp, delta_f):
    """
    Function to interpolate a PSD to a different value of delta_f. Uses linear
    interpolation.

    Parameters
    ----------
    psd_f : numpy.array or list or similar
        List of the frequencies contained within the PSD.
    psd_amp : numpy.array or list or similar
        List of the PSD values at the frequencies in psd_f.
    delta_f : float
        Value of delta_f to interpolate the PSD to.

    Returns
    -------
    new_psd_f : numpy.array
       Array of the frequencies contained within the interpolated PSD
    new_psd_amp : numpy.array
       Array of the interpolated PSD values at the frequencies in new_psd_f.

    """
    # In some cases this will be a no-op. I thought about removing this, but
    # this function can take unequally sampled PSDs and it is difficult to
    # check for this. As this function runs quickly anyway (compared to the
    # moment calculation) I decided to always interpolate.

    new_psd_f = []
    new_psd_amp = []
    fcurr = psd_f[0]

    for i in range(len(psd_f) - 1):
        f_low = psd_f[i]
        f_high = psd_f[i + 1]
        amp_low = psd_amp[i]
        amp_high = psd_amp[i + 1]
        while 1:
            if fcurr > f_high:
                break
            new_psd_f.append(fcurr)
            gradient = (amp_high - amp_low) / (f_high - f_low)
            f_diff = fcurr - f_low
            new_psd_amp.append(amp_low + f_diff * gradient)
            fcurr = fcurr + delta_f
    return numpy.asarray(new_psd_f), numpy.asarray(new_psd_amp)


def calculate_moment(
    psd_f, psd_amp, fmin, fmax, f0, funct, norm=None, vary_fmax=False, vary_density=None
):
    r"""
    Function for calculating one of the integrals used to construct a template
    bank placement metric. The integral calculated will be

    \int funct(x) * (psd_x)**(-7./3.) * delta_x / PSD(x)

    where x = f / f0. The lower frequency cutoff is given by fmin, see
    the parameters below for details on how the upper frequency cutoff is
    chosen

    Parameters
    ----------
    psd_f : numpy.array
       numpy array holding the set of evenly spaced frequencies used in the PSD
    psd_amp : numpy.array
       numpy array holding the PSD values corresponding to the psd_f
       frequencies
    fmin : float
        The lower frequency cutoff used in the calculation of the integrals
        used to obtain the metric.
    fmax : float
        The upper frequency cutoff used in the calculation of the integrals
        used to obtain the metric. This can be varied (see the vary_fmax
        option below).
    f0 : float
        This is an arbitrary scaling factor introduced to avoid the potential
        for numerical overflow when calculating this. Generally the default
        value (70) is safe here.
    funct : Lambda function
        The function to use when computing the integral as described above.
    norm : Dictionary of floats
        If given then moment[f_cutoff] will be divided by norm[f_cutoff]
    vary_fmax : boolean, optional (default False)
        If set to False the metric and rotations are calculated once, for the
        full range of frequency [f_low,f_upper).
        If set to True the metric and rotations are calculated multiple times,
        for frequency ranges [f_low,f_low + i*vary_density), where i starts at
        1 and runs up until f_low + (i+1)*vary_density > f_upper.
        Thus values greater than f_upper are *not* computed.
        The calculation for the full range [f_low,f_upper) is also done.
    vary_density : float, optional
        If vary_fmax is True, this will be used in computing the frequency
        ranges as described for vary_fmax.

    Returns
    -------
    moment : Dictionary of floats
        moment[f_cutoff] will store the value of the moment at the frequency
        cutoff given by f_cutoff.

    """
    # Must ensure delta_f in psd_f is constant
    psd_x = psd_f / f0
    deltax = psd_x[1] - psd_x[0]

    mask = numpy.logical_and(psd_f > fmin, psd_f < fmax)
    psdf_red = psd_f[mask]
    comps_red = (
        psd_x[mask] ** (-7.0 / 3.0) * funct(psd_x[mask], f0) * deltax / psd_amp[mask]
    )
    moment = {}
    moment[fmax] = comps_red.sum()
    if norm:
        moment[fmax] = moment[fmax] / norm[fmax]
    if vary_fmax:
        for t_fmax in numpy.arange(fmin + vary_density, fmax, vary_density):
            moment[t_fmax] = comps_red[psdf_red < t_fmax].sum()
            if norm:
                moment[t_fmax] = moment[t_fmax] / norm[t_fmax]
    return moment


def calculate_metric(Js, log_js, loglog_js, logloglog_js, loglogloglog_js, mapping):  # noqa: N803
    """
    This function will take the various integrals calculated by get_moments and
    convert this into a metric for the appropriate parameter space.

    Parameters
    ----------
    Js : Dictionary
        The list of (log^0 x) * x**(-i/3) integrals computed by get_moments()
        The index is Js[i]
    log_js : Dictionary
        The list of (log^1 x) * x**(-i/3) integrals computed by get_moments()
        The index is log_js[i]
    loglog_js : Dictionary
        The list of (log^2 x) * x**(-i/3) integrals computed by get_moments()
        The index is loglog_js[i]
    logloglog_js : Dictionary
        The list of (log^3 x) * x**(-i/3) integrals computed by get_moments()
        The index is logloglog_js[i]
    loglogloglog_js : Dictionary
        The list of (log^4 x) * x**(-i/3) integrals computed by get_moments()
        The index is loglogloglog_js[i]
    mapping : dictionary
        Used to identify which Lambda components are active in this parameter
        space and map these to entries in the metric matrix.

    Returns
    -------
    metric : numpy.matrix
        The resulting metric.

    """
    # How many dimensions in the parameter space?
    max_len = len(mapping.keys())

    metric = numpy.zeros(shape=(max_len, max_len), dtype=float)
    unmax_metric = numpy.zeros(shape=(max_len + 1, max_len + 1), dtype=float)

    for i in range(16):
        for j in range(16):
            calculate_metric_comp(
                metric,
                unmax_metric,
                i,
                j,
                Js,
                log_js,
                loglog_js,
                logloglog_js,
                loglogloglog_js,
                mapping,
            )
    return metric, unmax_metric


def calculate_metric_comp(
    gs,
    unmax_metric,
    i,
    j,
    Js,  # noqa: N803
    log_js,
    loglog_js,
    logloglog_js,
    loglogloglog_js,
    mapping,
):
    """
    Used to compute part of the metric. Only call this from within
    calculate_metric(). Please see the documentation for that function.
    """
    # Time term in unmax_metric. Note that these terms are recomputed a bunch
    # of time, but this cost is insignificant compared to computing the moments
    unmax_metric[-1, -1] = Js[1] - Js[4] * Js[4]

    # Normal terms
    if f"Lambda{i}" in mapping and f"Lambda{j}" in mapping:
        gammaij = Js[17 - i - j] - Js[12 - i] * Js[12 - j]
        gamma0i = Js[9 - i] - Js[4] * Js[12 - i]
        gamma0j = Js[9 - j] - Js[4] * Js[12 - j]
        gs[mapping[f"Lambda{i}"], mapping[f"Lambda{j}"]] = 0.5 * (
            gammaij - gamma0i * gamma0j / (Js[1] - Js[4] * Js[4])
        )
        unmax_metric[mapping[f"Lambda{i}"], -1] = gamma0i
        unmax_metric[-1, mapping[f"Lambda{j}"]] = gamma0j
        unmax_metric[mapping[f"Lambda{i}"], mapping[f"Lambda{j}"]] = gammaij
    # Normal,log cross terms
    if f"Lambda{i}" in mapping and f"LogLambda{j}" in mapping:
        gammaij = log_js[17 - i - j] - log_js[12 - j] * Js[12 - i]
        gamma0i = Js[9 - i] - Js[4] * Js[12 - i]
        gamma0j = log_js[9 - j] - log_js[12 - j] * Js[4]
        gs[mapping[f"Lambda{i}"], mapping[f"LogLambda{j}"]] = gs[
            mapping[f"LogLambda{j}"], mapping[f"Lambda{i}"]
        ] = 0.5 * (gammaij - gamma0i * gamma0j / (Js[1] - Js[4] * Js[4]))
        unmax_metric[mapping[f"Lambda{i}"], -1] = gamma0i
        unmax_metric[-1, mapping[f"Lambda{i}"]] = gamma0i
        unmax_metric[-1, mapping[f"LogLambda{j}"]] = gamma0j
        unmax_metric[mapping[f"LogLambda{j}"], -1] = gamma0j
        unmax_metric[mapping[f"Lambda{i}"], mapping[f"LogLambda{j}"]] = gammaij
        unmax_metric[mapping[f"LogLambda{j}"], mapping[f"Lambda{i}"]] = gammaij
    # Log,log terms
    if f"LogLambda{i}" in mapping and f"LogLambda{j}" in mapping:
        gammaij = loglog_js[17 - i - j] - log_js[12 - j] * log_js[12 - i]
        gamma0i = log_js[9 - i] - Js[4] * log_js[12 - i]
        gamma0j = log_js[9 - j] - log_js[12 - j] * Js[4]
        gs[mapping[f"LogLambda{i}"], mapping[f"LogLambda{j}"]] = 0.5 * (
            gammaij - gamma0i * gamma0j / (Js[1] - Js[4] * Js[4])
        )
        unmax_metric[mapping[f"LogLambda{i}"], -1] = gamma0i
        unmax_metric[-1, mapping[f"LogLambda{j}"]] = gamma0j
        unmax_metric[mapping[f"LogLambda{i}"], mapping[f"LogLambda{j}"]] = gammaij

    # Normal,loglog cross terms
    if f"Lambda{i}" in mapping and f"LogLogLambda{j}" in mapping:
        gammaij = loglog_js[17 - i - j] - loglog_js[12 - j] * Js[12 - i]
        gamma0i = Js[9 - i] - Js[4] * Js[12 - i]
        gamma0j = loglog_js[9 - j] - loglog_js[12 - j] * Js[4]
        gs[mapping[f"Lambda{i}"], mapping[f"LogLogLambda{j}"]] = gs[
            mapping[f"LogLogLambda{j}"], mapping[f"Lambda{i}"]
        ] = 0.5 * (gammaij - gamma0i * gamma0j / (Js[1] - Js[4] * Js[4]))
        unmax_metric[mapping[f"Lambda{i}"], -1] = gamma0i
        unmax_metric[-1, mapping[f"Lambda{i}"]] = gamma0i
        unmax_metric[-1, mapping[f"LogLogLambda{j}"]] = gamma0j
        unmax_metric[mapping[f"LogLogLambda{j}"], -1] = gamma0j
        unmax_metric[mapping[f"Lambda{i}"], mapping[f"LogLogLambda{j}"]] = gammaij
        unmax_metric[mapping[f"LogLogLambda{j}"], mapping[f"Lambda{i}"]] = gammaij

    # log,loglog cross terms
    if f"LogLambda{i}" in mapping and f"LogLogLambda{j}" in mapping:
        gammaij = logloglog_js[17 - i - j] - loglog_js[12 - j] * log_js[12 - i]
        gamma0i = log_js[9 - i] - Js[4] * log_js[12 - i]
        gamma0j = loglog_js[9 - j] - loglog_js[12 - j] * Js[4]
        gs[mapping[f"LogLambda{i}"], mapping[f"LogLogLambda{j}"]] = gs[
            mapping[f"LogLogLambda{j}"], mapping[f"LogLambda{i}"]
        ] = 0.5 * (gammaij - gamma0i * gamma0j / (Js[1] - Js[4] * Js[4]))
        unmax_metric[mapping[f"LogLambda{i}"], -1] = gamma0i
        unmax_metric[-1, mapping[f"LogLambda{i}"]] = gamma0i
        unmax_metric[-1, mapping[f"LogLogLambda{j}"]] = gamma0j
        unmax_metric[mapping[f"LogLogLambda{j}"], -1] = gamma0j
        unmax_metric[mapping[f"LogLambda{i}"], mapping[f"LogLogLambda{j}"]] = gammaij
        unmax_metric[mapping[f"LogLogLambda{j}"], mapping[f"LogLambda{i}"]] = gammaij

    # Loglog,loglog terms
    if f"LogLogLambda{i}" in mapping and f"LogLogLambda{j}" in mapping:
        gammaij = loglogloglog_js[17 - i - j] - loglog_js[12 - j] * loglog_js[12 - i]
        gamma0i = loglog_js[9 - i] - Js[4] * loglog_js[12 - i]
        gamma0j = loglog_js[9 - j] - loglog_js[12 - j] * Js[4]
        gs[mapping[f"LogLogLambda{i}"], mapping[f"LogLogLambda{j}"]] = 0.5 * (
            gammaij - gamma0i * gamma0j / (Js[1] - Js[4] * Js[4])
        )
        unmax_metric[mapping[f"LogLogLambda{i}"], -1] = gamma0i
        unmax_metric[-1, mapping[f"LogLogLambda{j}"]] = gamma0j
        unmax_metric[mapping[f"LogLogLambda{i}"], mapping[f"LogLogLambda{j}"]] = gammaij
