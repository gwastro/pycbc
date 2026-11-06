import logging

import numpy

from pycbc.tmpltbank.coord_utils import get_cov_params

logger = logging.getLogger(__name__)


# best_xis is unused, but kept so that existing calls keep working
def get_physical_covaried_masses(
    xis,
    best_masses,
    best_xis,  # noqa: ARG001
    req_match,
    mass_range_params,
    metric_params,
    f_upper,
    give_up_thresh=5000,
):
    """
    This function takes the position of a point in the xi parameter space and
    iteratively finds a close point in the physical coordinate space (masses
    and spins).

    Parameters
    ----------
    xis : list or array
        Desired position of the point in the xi space. If only N values are
        provided and the xi space's dimension is larger then it is assumed that
        *any* value in the remaining xi coordinates is acceptable.
    best_masses : list
        Contains [total_mass, eta, spin1z, spin2z]. Is a physical position
        mapped to xi coordinates in best_xis that is close to the desired point.
        This is aimed to give the code a starting point.
    best_xis : list
        Contains the position of best_masses in the xi coordinate system.
    req_match : float
        Desired maximum mismatch between xis and the obtained point. If a point
        is found with mismatch < req_match immediately stop and return that
        point. A point with this mismatch will not always be found.
    mass_range_params : MassRangeParameters instance
        Instance holding all the details of mass ranges and spin ranges.
    metric_params : MetricParameters instance
        Structure holding all the options for construction of the metric
        and the eigenvalues, eigenvectors and covariance matrix
        needed to manipulate the space.
    f_upper : float
        The value of f_upper that was used when obtaining the xi_i
        coordinates. This lets us know how to rotate potential physical points
        into the correct xi_i space. This must be a key in metric_params.evals,
        metric_params.evecs and metric_params.evecs_cv
        (ie. we must know how to do the transformation for
        the given value of f_upper)
    give_up_thresh : int, optional (default = 5000)
        The program will try this many iterations. If no close matching point
        has been found after this it will give up.

    Returns
    -------
    mass1 : float
        The heavier mass of the obtained point.
    mass2 : float
        The smaller mass of the obtained point
    spin1z : float
        The heavier bodies spin of the obtained point.
    spin2z : float
        The smaller bodies spin of the obtained point.
    count : int
        How many iterations it took to find the point. For debugging.
    mismatch : float
        The mismatch between the obtained point and the input xis.
    new_xis : list
        The position of the point in the xi space

    """
    # TUNABLE PARAMETERS GO HERE!
    # This states how far apart to scatter test points in the first proposal
    orig_scale_factor = 1

    # Set up
    xi_size = len(xis)
    scale_factor = orig_scale_factor
    best_chirpmass = best_masses[0] * (best_masses[1]) ** (3.0 / 5.0)
    count = 0
    un_fixed_count = 0
    curr_dist = 100000000000000000
    while 1:
        # If we are a long way away we use larger jumps
        if count:
            if curr_dist > 1 and scale_factor == orig_scale_factor:
                scale_factor = orig_scale_factor * 10
        # Get a set of test points with mass -> xi mappings
        totmass, eta, spin1z, spin2z, mass1, mass2, new_xis = get_mass_distribution(
            [best_chirpmass, best_masses[1], best_masses[2], best_masses[3]],
            scale_factor,
            mass_range_params,
            metric_params,
            f_upper,
        )
        c_dist = (new_xis[0] - xis[0]) ** 2
        for j in range(1, xi_size):
            c_dist += (new_xis[j] - xis[j]) ** 2
        if c_dist.min() < req_match:
            idx = c_dist.argmin()
            scale_factor = orig_scale_factor
            new_xis_list = [new_xis[ldx][idx] for ldx in range(len(new_xis))]
            return (
                mass1[idx],
                mass2[idx],
                spin1z[idx],
                spin2z[idx],
                count,
                c_dist.min(),
                new_xis_list,
            )
        if c_dist.min() < curr_dist:
            idx = c_dist.argmin()
            best_masses[0] = totmass[idx]
            best_masses[1] = eta[idx]
            best_masses[2] = spin1z[idx]
            best_masses[3] = spin2z[idx]
            best_chirpmass = best_masses[0] * (best_masses[1]) ** (3.0 / 5.0)
            curr_dist = c_dist.min()
            un_fixed_count = 0
            scale_factor = orig_scale_factor
        count += 1
        un_fixed_count += 1
        if un_fixed_count > give_up_thresh:
            # Stop at this point
            diff = (best_masses[0] * best_masses[0] * (1 - 4 * best_masses[1])) ** 0.5
            mass1 = (best_masses[0] + diff) / 2.0
            mass2 = (best_masses[0] - diff) / 2.0
            new_xis_list = [new_xis[ldx][0] for ldx in range(len(new_xis))]
            return (
                mass1,
                mass2,
                best_masses[2],
                best_masses[3],
                count,
                curr_dist,
                new_xis_list,
            )
        if not un_fixed_count % 100:
            scale_factor *= 2
        if scale_factor > 64:
            scale_factor = 1
    # Shouldn't be here!
    raise RuntimeError


# spin2z_jump_fac is unused, but kept so that existing calls keep working
def get_mass_distribution(
    best_masses,
    scale_factor,
    mass_range_params,
    metric_params,
    f_upper,
    num_jump_points=100,
    chirp_mass_jump_fac=0.0001,
    eta_jump_fac=0.01,
    spin1z_jump_fac=0.01,
    spin2z_jump_fac=0.01,  # noqa: ARG001
):
    """
    Given a set of masses, this function will create a set of points nearby
    in the mass space and map these to the xi space.

    Parameters
    ----------
    best_masses : list
        Contains [ChirpMass, eta, spin1z, spin2z]. Points will be placed around
        tjos
    scale_factor : float
        This parameter describes the radius away from best_masses that points
        will be placed in.
    mass_range_params : MassRangeParameters instance
        Instance holding all the details of mass ranges and spin ranges.
    metric_params : MetricParameters instance
        Structure holding all the options for construction of the metric
        and the eigenvalues, eigenvectors and covariance matrix
        needed to manipulate the space.
    f_upper : float
        The value of f_upper that was used when obtaining the xi_i
        coordinates. This lets us know how to rotate potential physical points
        into the correct xi_i space. This must be a key in metric_params.evals,
        metric_params.evecs and metric_params.evecs_cv
        (ie. we must know how to do the transformation for
        the given value of f_upper)
    num_jump_points : int, optional (default = 100)
        The number of points that will be generated every iteration
    chirp_mass_jump_fac : float, optional (default=0.0001)
        The jump points will be chosen with fractional variation in chirpMass
        up to this multiplied by scale_factor.
    eta_jump_fac : float, optional (default=0.01)
        The jump points will be chosen with fractional variation in eta
        up to this multiplied by scale_factor.
    spin1z_jump_fac : float, optional (default=0.01)
        The jump points will be chosen with absolute variation in spin1z up to
        this multiplied by scale_factor.
    spin2z_jump_fac : float, optional (default=0.01)
        The jump points will be chosen with absolute variation in spin2z up to
        this multiplied by scale_factor.

    Returns
    -------
    Totmass : numpy.array
        Total mass of the resulting points
    Eta : numpy.array
        Symmetric mass ratio of the resulting points
    Spin1z : numpy.array
        Spin of the heavier body of the resulting points
    Spin2z : numpy.array
        Spin of the smaller body of the resulting points
    Diff : numpy.array
        Mass1 - Mass2 of the resulting points
    Mass1 : numpy.array
        Mass1 (mass of heavier body) of the resulting points
    Mass2 : numpy.array
        Mass2 (mass of smaller body) of the resulting points
    new_xis : list of numpy.array
        Position of points in the xi coordinates

    """
    # FIXME: It would be better if rejected values could be drawn from the
    # full possible mass/spin distribution. However speed in this function is
    # a major factor and must be considered.
    best_chirpmass = best_masses[0]
    best_eta = best_masses[1]
    best_spin1z = best_masses[2]
    best_spin2z = best_masses[3]

    # Firstly choose a set of values for masses and spins
    chirpmass = best_chirpmass * (
        1
        - (numpy.random.random(num_jump_points) - 0.5)
        * chirp_mass_jump_fac
        * scale_factor
    )
    eta_range = mass_range_params.max_eta - mass_range_params.min_eta
    curr_jump_fac = eta_jump_fac * scale_factor
    curr_jump_fac = min(curr_jump_fac, eta_range)
    eta = best_eta * (1 - (numpy.random.random(num_jump_points) - 0.5) * curr_jump_fac)

    max_spin_mag = max(
        mass_range_params.max_ns_spin_mag, mass_range_params.max_bh_spin_mag
    )
    min_spin_mag = min(
        mass_range_params.max_ns_spin_mag, mass_range_params.max_bh_spin_mag
    )
    # Note that these two are cranged by spinxzFac, *not* spinxzFac/spinxz
    curr_jump_fac = spin1z_jump_fac * scale_factor
    curr_jump_fac = min(curr_jump_fac, max_spin_mag)

    # Actually set the new spin trial points
    if mass_range_params.nsbh_flag or (max_spin_mag == min_spin_mag):
        curr_spin_1z_jump_fac = curr_jump_fac
        curr_spin_2z_jump_fac = curr_jump_fac
        # Check spins aren't going to be unphysical
        if curr_jump_fac > mass_range_params.max_bh_spin_mag:
            curr_spin_1z_jump_fac = mass_range_params.max_bh_spin_mag
        if curr_jump_fac > mass_range_params.max_ns_spin_mag:
            curr_spin_2z_jump_fac = mass_range_params.max_ns_spin_mag
        spin1z = best_spin1z + (
            (numpy.random.random(num_jump_points) - 0.5) * curr_spin_1z_jump_fac
        )
        spin2z = best_spin2z + (
            (numpy.random.random(num_jump_points) - 0.5) * curr_spin_2z_jump_fac
        )
    else:
        # If max_ns_spin_mag is very low (0) and max_bh_spin_mag is high we can
        # find it hard to place any points. So mix these when
        # masses are swapping between the NS and BH.
        curr_spin_bh_jump_fac = curr_jump_fac
        curr_spin_ns_jump_fac = curr_jump_fac
        # Check spins aren't going to be unphysical
        if curr_jump_fac > mass_range_params.max_bh_spin_mag:
            curr_spin_bh_jump_fac = mass_range_params.max_bh_spin_mag
        if curr_jump_fac > mass_range_params.max_ns_spin_mag:
            curr_spin_ns_jump_fac = mass_range_params.max_ns_spin_mag
        spin1z = numpy.zeros(num_jump_points, dtype=float)
        spin2z = numpy.zeros(num_jump_points, dtype=float)
        split_point = int(num_jump_points / 2)
        # So set the first half to be at least within the BH range and the
        # second half to be at least within the NS range
        spin1z[:split_point] = best_spin1z + (
            (numpy.random.random(split_point) - 0.5) * curr_spin_bh_jump_fac
        )
        spin1z[split_point:] = best_spin1z + (
            (numpy.random.random(num_jump_points - split_point) - 0.5)
            * curr_spin_ns_jump_fac
        )
        spin2z[:split_point] = best_spin2z + (
            (numpy.random.random(split_point) - 0.5) * curr_spin_bh_jump_fac
        )
        spin2z[split_point:] = best_spin2z + (
            (numpy.random.random(num_jump_points - split_point) - 0.5)
            * curr_spin_ns_jump_fac
        )

    # Point[0] is always set to the original point
    chirpmass[0] = best_chirpmass
    eta[0] = best_eta
    spin1z[0] = best_spin1z
    spin2z[0] = best_spin2z

    # Remove points where eta becomes unphysical
    eta[eta > mass_range_params.max_eta] = mass_range_params.max_eta
    if mass_range_params.min_eta:
        eta[eta < mass_range_params.min_eta] = mass_range_params.min_eta
    else:
        eta[eta < 0.0001] = 0.0001

    # Total mass, masses and mass diff
    totmass = chirpmass / (eta ** (3.0 / 5.0))
    diff = (totmass * totmass * (1 - 4 * eta)) ** 0.5
    mass1 = (totmass + diff) / 2.0
    mass2 = (totmass - diff) / 2.0

    # Check the validity of the spin values
    # Do the first spin

    if max_spin_mag == 0:
        # Shortcut if non-spinning
        pass
    elif mass_range_params.nsbh_flag or (max_spin_mag == min_spin_mag):
        # Simple case where I don't have to worry about correlation with mass
        numploga = abs(spin1z) > mass_range_params.max_bh_spin_mag
        spin1z[numploga] = 0
    else:
        # Do have to consider masses
        boundary_mass = mass_range_params.ns_bh_boundary_mass
        numploga1 = numpy.logical_and(
            mass1 >= boundary_mass, abs(spin1z) <= mass_range_params.max_bh_spin_mag
        )
        numploga2 = numpy.logical_and(
            mass1 < boundary_mass, abs(spin1z) <= mass_range_params.max_ns_spin_mag
        )
        numploga = numpy.logical_or(numploga1, numploga2)
        numploga = numpy.logical_not(numploga)
        spin1z[numploga] = 0

    # Same for the second spin

    if max_spin_mag == 0:
        # Shortcut if non-spinning
        pass
    elif mass_range_params.nsbh_flag or (max_spin_mag == min_spin_mag):
        numplogb = abs(spin2z) > mass_range_params.max_ns_spin_mag
        spin2z[numplogb] = 0
    else:
        # Do have to consider masses
        boundary_mass = mass_range_params.ns_bh_boundary_mass
        numplogb1 = numpy.logical_and(
            mass2 >= boundary_mass, abs(spin2z) <= mass_range_params.max_bh_spin_mag
        )
        numplogb2 = numpy.logical_and(
            mass2 < boundary_mass, abs(spin2z) <= mass_range_params.max_ns_spin_mag
        )
        numplogb = numpy.logical_or(numplogb1, numplogb2)
        numplogb = numpy.logical_not(numplogb)
        spin2z[numplogb] = 0

    if (max_spin_mag) and (numploga[0] or numplogb[0]):
        raise ValueError("Cannot remove the guide point!")

    # And remove points where the individual masses are outside of the physical
    # range. Or the total masses are.
    # These "removed" points will have metric distances that will be much, much
    # larger than any thresholds used in the functions in brute_force_utils.py
    # and will always be rejected. An unphysical value cannot be used as it
    # would result in unphysical metric distances and cause failures.
    totmass[mass1 < mass_range_params.min_mass1 * 0.9999] = 0.0001
    totmass[mass1 > mass_range_params.max_mass1 * 1.0001] = 0.0001
    totmass[mass2 < mass_range_params.min_mass2 * 0.9999] = 0.0001
    totmass[mass2 > mass_range_params.max_mass2 * 1.0001] = 0.0001
    # There is some numerical error which can push this a bit higher. We do
    # *not* want to reject the initial guide point. This error comes from
    # Masses -> totmass, eta -> masses conversion, we will have points pushing
    # onto the boudaries of the space.
    totmass[totmass > mass_range_params.max_tot_mass * 1.0001] = 0.0001
    totmass[totmass < mass_range_params.min_tot_mass * 0.9999] = 0.0001
    if mass_range_params.max_chirp_mass:
        totmass[chirpmass > mass_range_params.max_chirp_mass * 1.0001] = 0.0001
    if mass_range_params.min_chirp_mass:
        totmass[chirpmass < mass_range_params.min_chirp_mass * 0.9999] = 0.0001

    if totmass[0] < 0.00011:
        raise ValueError("Cannot remove the guide point!")

    mass1[totmass < 0.00011] = 0.0001
    mass2[totmass < 0.00011] = 0.0001

    # Then map to xis
    new_xis = get_cov_params(mass1, mass2, spin1z, spin2z, metric_params, f_upper)
    return totmass, eta, spin1z, spin2z, mass1, mass2, new_xis


def stack_xi_direction_brute(
    xis,
    best_masses,
    best_xis,
    direction_num,
    req_match,
    mass_range_params,
    metric_params,
    f_upper,
    scale_factor=0.8,
    num_iterations=3000,
):
    """
    This function is used to assess the depth of the xi_space in a specified
    dimension at a specified point in the higher dimensions. It does this by
    iteratively throwing points at the space to find maxima and minima.

    Parameters
    ----------
    xis : list or array
        Position in the xi space at which to assess the depth. This can be only
        a subset of the higher dimensions than that being sampled.
    best_masses : list
        Contains [total_mass, eta, spin1z, spin2z]. Is a physical position
        mapped to xi coordinates in best_xis that is close to the xis point.
        This is aimed to give the code a starting point.
    best_xis : list
        Contains the position of best_masses in the xi coordinate system.
    direction_num : int
        The dimension that you want to assess the depth of (0 = 1, 1 = 2 ...)
    req_match : float
        When considering points to assess the depth with, only consider points
        with a mismatch that is smaller than this with xis.
    mass_range_params : MassRangeParameters instance
        Instance holding all the details of mass ranges and spin ranges.
    metric_params : MetricParameters instance
        Structure holding all the options for construction of the metric
        and the eigenvalues, eigenvectors and covariance matrix
        needed to manipulate the space.
    f_upper : float
        The value of f_upper that was used when obtaining the xi_i
        coordinates. This lets us know how to rotate potential physical points
        into the correct xi_i space. This must be a key in metric_params.evals,
        metric_params.evecs and metric_params.evecs_cv
        (ie. we must know how to do the transformation for
        the given value of f_upper)
    scale_factor : float, optional (default = 0.8)
        The value of the scale factor that is used when calling
        pycbc.tmpltbank.get_mass_distribution.
    num_iterations : int, optional (default = 3000)
        The number of times to make calls to get_mass_distribution when
        assessing the maximum/minimum of this parameter space. Making this
        smaller makes the code faster, but at the cost of accuracy.

    Returns
    -------
    xi_min : float
        The minimal value of the specified dimension at the specified point in
        parameter space.
    xi_max : float
       The maximal value of the specified dimension at the specified point in
        parameter space.

    """
    # Find minimum
    ximin = find_xi_extrema_brute(
        xis,
        best_masses,
        best_xis,
        direction_num,
        req_match,
        mass_range_params,
        metric_params,
        f_upper,
        find_minimum=True,
        scale_factor=scale_factor,
        num_iterations=num_iterations,
    )

    # Find maximum
    ximax = find_xi_extrema_brute(
        xis,
        best_masses,
        best_xis,
        direction_num,
        req_match,
        mass_range_params,
        metric_params,
        f_upper,
        find_minimum=False,
        scale_factor=scale_factor,
        num_iterations=num_iterations,
    )

    return ximin, ximax


# best_xis is unused, but kept so that existing calls keep working
def find_xi_extrema_brute(
    xis,
    best_masses,
    best_xis,  # noqa: ARG001
    direction_num,
    req_match,
    mass_range_params,
    metric_params,
    f_upper,
    find_minimum=False,
    scale_factor=0.8,
    num_iterations=3000,
):
    """
    This function is used to find the largest or smallest value of the xi
    space in a specified
    dimension at a specified point in the higher dimensions. It does this by
    iteratively throwing points at the space to find extrema.

    Parameters
    ----------
    xis : list or array
        Position in the xi space at which to assess the depth. This can be only
        a subset of the higher dimensions than that being sampled.
    best_masses : list
        Contains [total_mass, eta, spin1z, spin2z]. Is a physical position
        mapped to xi coordinates in best_xis that is close to the xis point.
        This is aimed to give the code a starting point.
    best_xis : list
        Contains the position of best_masses in the xi coordinate system.
    direction_num : int
        The dimension that you want to assess the depth of (0 = 1, 1 = 2 ...)
    req_match : float
        When considering points to assess the depth with, only consider points
        with a mismatch that is smaller than this with xis.
    mass_range_params : MassRangeParameters instance
        Instance holding all the details of mass ranges and spin ranges.
    metric_params : MetricParameters instance
        Structure holding all the options for construction of the metric
        and the eigenvalues, eigenvectors and covariance matrix
        needed to manipulate the space.
    f_upper : float
        The value of f_upper that was used when obtaining the xi_i
        coordinates. This lets us know how to rotate potential physical points
        into the correct xi_i space. This must be a key in metric_params.evals,
        metric_params.evecs and metric_params.evecs_cv
        (ie. we must know how to do the transformation for
        the given value of f_upper)
    find_minimum : boolean, optional (default = False)
        If True, find the minimum value of the xi direction. If False find the
        maximum value.
    scale_factor : float, optional (default = 0.8)
        The value of the scale factor that is used when calling
        pycbc.tmpltbank.get_mass_distribution.
    num_iterations : int, optional (default = 3000)
        The number of times to make calls to get_mass_distribution when
        assessing the maximum/minimum of this parameter space. Making this
        smaller makes the code faster, but at the cost of accuracy.

    Returns
    -------
    xi_extent : float
        The extremal value of the specified dimension at the specified point in
        parameter space.

    """
    # Setup
    xi_size = len(xis)
    best_chirpmass = best_masses[0] * (best_masses[1]) ** (3.0 / 5.0)
    if find_minimum:
        xiextrema = 10000000000
    else:
        xiextrema = -100000000000

    for _ in range(num_iterations):
        # Evaluate extrema of the xi direction specified
        totmass, eta, spin1z, spin2z, _, _, new_xis = get_mass_distribution(
            [best_chirpmass, best_masses[1], best_masses[2], best_masses[3]],
            scale_factor,
            mass_range_params,
            metric_params,
            f_upper,
        )
        c_dist = (new_xis[0] - xis[0]) ** 2
        for j in range(1, xi_size):
            c_dist += (new_xis[j] - xis[j]) ** 2
        red_c_dist = c_dist[c_dist < req_match]
        if len(red_c_dist):
            if not find_minimum:
                new_xis[direction_num][c_dist > req_match] = -10000000
                curr_xi_extrema = (new_xis[direction_num]).max()
                idx = (new_xis[direction_num]).argmax()
            else:
                new_xis[direction_num][c_dist > req_match] = 10000000
                curr_xi_extrema = (new_xis[direction_num]).min()
                idx = (new_xis[direction_num]).argmin()
            if ((not find_minimum) and (curr_xi_extrema > xiextrema)) or (
                find_minimum and (curr_xi_extrema < xiextrema)
            ):
                xiextrema = curr_xi_extrema
                best_masses[0] = totmass[idx]
                best_masses[1] = eta[idx]
                best_masses[2] = spin1z[idx]
                best_masses[3] = spin2z[idx]
                best_chirpmass = best_masses[0] * (best_masses[1]) ** (3.0 / 5.0)
    return xiextrema
