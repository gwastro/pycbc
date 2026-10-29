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

from pycbc import conversions, pnutils
from pycbc.neutron_stars import load_ns_sequence
from pycbc.tmpltbank.lambda_mapping import get_chirp_params

logger = logging.getLogger(__name__)


def estimate_mass_range(
    num_points, mass_range_params, metric_params, f_upper, covary=True
):
    """
    This function will generate a large set of points with random masses and
    spins (using pycbc.tmpltbank.get_random_mass) and translate these points
    into the xi_i coordinate system for the given upper frequency cutoff.

    Parameters
    ----------
    num_points : int
        Number of systems to simulate
    mass_range_params : MassRangeParameters instance
        Instance holding all the details of mass ranges and spin ranges.
    metric_params : MetricParameters instance
        Structure holding all the options for construction of the metric
        and the eigenvalues, eigenvectors and covariance matrix
        needed to manipulate the space.
    f_upper : float
        The value of f_upper to use when getting the mu coordinates from the
        lambda coordinates. This must be a key in metric_params.evals and
        metric_params.evecs (ie. we must know how to do the transformation for
        the given value of f_upper). It also must be a key in
        metric_params.evecs_cv if covary=True.
    covary : boolean, optional (default = True)
        If this is given then evecs_cv will be used to rotate from the Cartesian
        coordinate system into the principal coordinate direction (xi_i). If
        not given then points in the original Cartesian coordinates are
        returned.


    Returns
    -------
    xis : numpy.array
        A list of the positions of each point in the xi_i coordinate system.

    """
    vals_set = get_random_mass(num_points, mass_range_params)
    mass1 = vals_set[0]
    mass2 = vals_set[1]
    spin1z = vals_set[2]
    spin2z = vals_set[3]
    if covary:
        lambdas = get_cov_params(mass1, mass2, spin1z, spin2z, metric_params, f_upper)
    else:
        lambdas = get_conv_params(mass1, mass2, spin1z, spin2z, metric_params, f_upper)

    return numpy.array(lambdas)


def get_random_mass_point_particles(num_points, mass_range_params):
    """
    This function will generate a large set of points within the chosen mass
    and spin space. It will also return the corresponding PN spin coefficients
    for ease of use later (though these may be removed at some future point).

    Parameters
    ----------
    num_points : int
        Number of systems to simulate
    mass_range_params : MassRangeParameters instance
        Instance holding all the details of mass ranges and spin ranges.

    Returns
    -------
    mass1 : float
        Mass of heavier body.
    mass2 : float
        Mass of lighter body.
    spin1z : float
        Spin of body 1.
    spin2z : float
        Spin of body 2.

    """
    # WARNING: We expect mass1 > mass2 ALWAYS

    # First we choose the total masses from a unifrom distribution in mass
    # to the -5/3. power.
    mass = numpy.random.random(num_points) * (
        mass_range_params.min_tot_mass ** (-5.0 / 3.0)
        - mass_range_params.max_tot_mass ** (-5.0 / 3.0)
    ) + mass_range_params.max_tot_mass ** (-5.0 / 3.0)
    mass = mass ** (-3.0 / 5.0)

    # Next we choose the mass ratios, this will take different limits based on
    # the value of total mass
    maxmass2 = numpy.minimum(mass / 2.0, mass_range_params.max_mass2)
    minmass1 = numpy.maximum(mass_range_params.min_mass1, mass / 2.0)
    mineta = numpy.maximum(
        mass_range_params.min_comp_mass
        * (mass - mass_range_params.min_comp_mass)
        / (mass * mass),
        mass_range_params.max_comp_mass
        * (mass - mass_range_params.max_comp_mass)
        / (mass * mass),
    )
    # Note that mineta is a numpy.array because mineta depends on the total
    # mass. Therefore this is not precomputed in the mass_range_params instance
    if mass_range_params.min_eta:
        mineta = numpy.maximum(mass_range_params.min_eta, mineta)
    # Eta also restricted by chirp mass restrictions
    if mass_range_params.min_chirp_mass:
        eta_val_at_min_chirp = mass_range_params.min_chirp_mass / mass
        eta_val_at_min_chirp = eta_val_at_min_chirp ** (5.0 / 3.0)
        mineta = numpy.maximum(mineta, eta_val_at_min_chirp)

    maxeta = numpy.minimum(
        mass_range_params.max_eta, maxmass2 * (mass - maxmass2) / (mass * mass)
    )
    maxeta = numpy.minimum(maxeta, minmass1 * (mass - minmass1) / (mass * mass))
    # max eta also affected by chirp mass restrictions
    if mass_range_params.max_chirp_mass:
        eta_val_at_max_chirp = mass_range_params.max_chirp_mass / mass
        eta_val_at_max_chirp = eta_val_at_max_chirp ** (5.0 / 3.0)
        maxeta = numpy.minimum(maxeta, eta_val_at_max_chirp)

    if (maxeta < mineta).any():
        err_msg = "ERROR: Maximum eta is smaller than minimum eta!!"
        raise ValueError(err_msg)
    eta = numpy.random.random(num_points) * (maxeta - mineta) + mineta

    # Also calculate the component masses; mass1 > mass2
    diff = (mass * mass * (1 - 4 * eta)) ** 0.5
    mass1 = (mass + diff) / 2.0
    mass2 = (mass - diff) / 2.0
    # Check the masses are where we want them to be (allowing some floating
    # point rounding error).
    if (mass1 > mass_range_params.max_mass1 * 1.001).any() or (
        mass1 < mass_range_params.min_mass1 * 0.999
    ).any():
        err_msg = "Mass1 is not within the specified mass range."
        raise ValueError(err_msg)
    if (mass2 > mass_range_params.max_mass2 * 1.001).any() or (
        mass2 < mass_range_params.min_mass2 * 0.999
    ).any():
        err_msg = "Mass2 is not within the specified mass range."
        raise ValueError(err_msg)

    # Next up is the spins. First check if we have non-zero spins
    if (
        mass_range_params.max_ns_spin_mag == 0
        and mass_range_params.max_bh_spin_mag == 0
    ):
        spin1z = numpy.zeros(num_points, dtype=float)
        spin2z = numpy.zeros(num_points, dtype=float)
    elif mass_range_params.nsbh_flag:
        # Spin 1 first
        mspin = numpy.zeros(len(mass1))
        mspin += mass_range_params.max_bh_spin_mag
        spin1z = (2 * numpy.random.random(num_points) - 1) * mspin
        # Then spin2
        mspin = numpy.zeros(len(mass2))
        mspin += mass_range_params.max_ns_spin_mag
        spin2z = (2 * numpy.random.random(num_points) - 1) * mspin
    else:
        boundary_mass = mass_range_params.ns_bh_boundary_mass
        # Spin 1 first
        mspin = numpy.zeros(len(mass1))
        mspin += mass_range_params.max_ns_spin_mag
        mspin[mass1 > boundary_mass] = mass_range_params.max_bh_spin_mag
        spin1z = (2 * numpy.random.random(num_points) - 1) * mspin
        # Then spin 2
        mspin = numpy.zeros(len(mass2))
        mspin += mass_range_params.max_ns_spin_mag
        mspin[mass2 > boundary_mass] = mass_range_params.max_bh_spin_mag
        spin2z = (2 * numpy.random.random(num_points) - 1) * mspin

    return mass1, mass2, spin1z, spin2z


def get_random_mass(num_points, mass_range_params, eos="2H"):
    """
    This function will generate a large set of points within the chosen mass
    and spin space, and with the desired minimum remnant disk mass (this applies
    to NS-BH systems only). It will also return the corresponding PN spin
    coefficients for ease of use later (though these may be removed at some
    future point).

    Parameters
    ----------
    num_points : int
        Number of systems to simulate
    mass_range_params : MassRangeParameters instance
        Instance holding all the details of mass ranges and spin ranges.
    eos : string
        Name of equation of state of neutron star.

    Returns
    -------
    mass1 : float
        Mass of heavier body.
    mass2 : float
        Mass of lighter body.
    spin1z : float
        Spin of body 1.
    spin2z : float
        Spin of body 2.

    """
    # WARNING: We expect mass1 > mass2 ALWAYS

    # Check if EM contraints are required, i.e. if the systems must produce
    # a minimum remnant disk mass.  If this is not the case, proceed treating
    # the systems as point particle binaries
    if mass_range_params.remnant_mass_threshold is None:
        mass1, mass2, spin1z, spin2z = get_random_mass_point_particles(
            num_points, mass_range_params
        )
    # otherwise, load EOS dependent data, generate the EM constraint
    # (i.e. compute the minimum symmetric mass ratio needed to
    # generate a given remnant disk mass as a function of the NS
    # mass and the BH spin along z) and then proceed by accepting
    # only systems that can yield (at least) the desired remnant
    # disk mass and that pass the mass and spin range cuts.
    else:
        max_ns_g_mass = load_ns_sequence(mass_range_params.ns_eos)[1]
        boundary_mass = mass_range_params.ns_bh_boundary_mass
        if max_ns_g_mass < boundary_mass:
            warn_msg = "WARNING: "
            warn_msg += f"Option of ns-bh-boundary-mass is  {boundary_mass} "
            warn_msg += "which is higher than the maximum NS gravitational "
            warn_msg += "mass admitted by the EOS that was prescribed "
            warn_msg += f"({max_ns_g_mass}). "
            warn_msg += "The code will proceed using the latter value "
            warn_msg += "as the boundary mass."
            logger.warning(warn_msg)
            boundary_mass = max_ns_g_mass

        # Empty arrays to store points that pass all cuts
        mass1_out = []
        mass2_out = []
        spin1z_out = []
        spin2z_out = []

        # As the EM cut can remove several randomly generated
        # binaries, track the number of accepted points that pass
        # all cuts and stop only once enough of them are generated
        num_points_found = 0
        while num_points_found < num_points:
            # Generate the random points within the required mass
            # and spin cuts
            mass1, mass2, spin1z, spin2z = get_random_mass_point_particles(
                num_points - num_points_found, mass_range_params
            )

            # Now proceed with cutting out EM dim systems
            # Use a logical mask to track points that do not correspond to
            # BBHs. The remaining points will be BNSs and NSBHs.
            # Further down, EM-dim NSBHs will also be removed.
            mask_not_bbh = numpy.zeros(len(mass1), dtype=bool)

            # Keep a point if:
            # 1) the secondary object is a not a BH (mass2 < boundary mass)
            #    [Store masses and spins of non BBH]
            mask_not_bbh[mass2 < boundary_mass] = True
            mass1_not_bbh = mass1[mask_not_bbh]
            mass2_not_bbh = mass2[mask_not_bbh]
            spin1z_not_bbh = spin1z[mask_not_bbh]
            spin2z_not_bbh = spin2z[mask_not_bbh]
            # 2) and if the primary mass is a NS (i.e., it is a BNS), or...
            mask_nsbh = numpy.zeros(len(mass1_not_bbh), dtype=bool)
            #    [mask_nsbh identifies NSBH systems]
            mask_nsbh[mass1_not_bbh > boundary_mass] = True
            #    [mask_bns identifies BNS systems]
            mask_bns = ~mask_nsbh
            #    [Store masses and spins of BNSs]
            mass1_bns = mass1_not_bbh[mask_bns]
            mass2_bns = mass2_not_bbh[mask_bns]
            spin1z_bns = spin1z_not_bbh[mask_bns]
            spin2z_bns = spin2z_not_bbh[mask_bns]
            # 3) ...it is an NS-BH with remnant mass greater than the threshold
            #    required to have a counterpart
            #    [Store masses and spins of all NSBHs]
            mass1_nsbh = mass1_not_bbh[mask_nsbh]
            mass2_nsbh = mass2_not_bbh[mask_nsbh]
            spin1z_nsbh = spin1z_not_bbh[mask_nsbh]
            spin2z_nsbh = spin2z_not_bbh[mask_nsbh]
            #    [Store etas of all NSBHs]
            eta_nsbh = conversions.eta_from_mass1_mass2(mass1_nsbh, mass2_nsbh)
            #    [mask_bright_nsbh will identify NSBH systems with high enough
            #     threshold mass]
            mask_bright_nsbh = numpy.zeros(len(mass1_nsbh), dtype=bool)
            if eta_nsbh.size != 0:
                remnant = conversions.remnant_mass_from_mass1_mass2_cartesian_spin_eos(
                    mass1_nsbh,
                    mass2_nsbh,
                    spin1x=0.0,
                    spin1y=0.0,
                    spin1z=spin1z_nsbh,
                    eos=eos,
                )
                mask_bright_nsbh[remnant > mass_range_params.remnant_mass_threshold] = (
                    True
                )

            # Keep only points that correspond to binaries that can produce an
            # EM counterpart (i.e., BNSs and EM-bright NSBHs) and add their
            # properties to the pile of accpeted points to output
            mass1_out = numpy.concatenate(
                (mass1_out, mass1_bns, mass1_nsbh[mask_bright_nsbh])
            )
            mass2_out = numpy.concatenate(
                (mass2_out, mass2_bns, mass2_nsbh[mask_bright_nsbh])
            )
            spin1z_out = numpy.concatenate(
                (spin1z_out, spin1z_bns, spin1z_nsbh[mask_bright_nsbh])
            )
            spin2z_out = numpy.concatenate(
                (spin2z_out, spin2z_bns, spin2z_nsbh[mask_bright_nsbh])
            )

            # Number of points that survived all cuts
            num_points_found = len(mass1_out)

        # Ready to go
        mass1 = mass1_out
        mass2 = mass2_out
        spin1z = spin1z_out
        spin2z = spin2z_out

    return mass1, mass2, spin1z, spin2z


def get_cov_params(
    mass1,
    mass2,
    spin1z,
    spin2z,
    metric_params,
    f_upper,
    lambda1=None,
    lambda2=None,
    quadparam1=None,
    quadparam2=None,
):
    """
    Function to convert between masses and spins and locations in the xi
    parameter space. Xi = Cartesian metric and rotated to principal components.

    Parameters
    ----------
    mass1 : float
        Mass of heavier body.
    mass2 : float
        Mass of lighter body.
    spin1z : float
        Spin of body 1.
    spin2z : float
        Spin of body 2.
    metric_params : MetricParameters instance
        Structure holding all the options for construction of the metric
        and the eigenvalues, eigenvectors and covariance matrix
        needed to manipulate the space.
    f_upper : float
        The value of f_upper to use when getting the mu coordinates from the
        lambda coordinates. This must be a key in metric_params.evals,
        metric_params.evecs and metric_params.evecs_cv
        (ie. we must know how to do the transformation for
        the given value of f_upper)

    Returns
    -------
    xis : list of floats or numpy.arrays
        Position of the system(s) in the xi coordinate system

    """
    # Do this by doing masses - > lambdas -> mus
    mus = get_conv_params(
        mass1,
        mass2,
        spin1z,
        spin2z,
        metric_params,
        f_upper,
        lambda1=lambda1,
        lambda2=lambda2,
        quadparam1=quadparam1,
        quadparam2=quadparam2,
    )
    # and then mus -> xis
    xis = get_covaried_params(mus, metric_params.evecs_cv[f_upper])
    return xis


def get_conv_params(
    mass1,
    mass2,
    spin1z,
    spin2z,
    metric_params,
    f_upper,
    lambda1=None,
    lambda2=None,
    quadparam1=None,
    quadparam2=None,
):
    """
    Function to convert between masses and spins and locations in the mu
    parameter space. Mu = Cartesian metric, but not principal components.

    Parameters
    ----------
    mass1 : float
        Mass of heavier body.
    mass2 : float
        Mass of lighter body.
    spin1z : float
        Spin of body 1.
    spin2z : float
        Spin of body 2.
    metric_params : MetricParameters instance
        Structure holding all the options for construction of the metric
        and the eigenvalues, eigenvectors and covariance matrix
        needed to manipulate the space.
    f_upper : float
        The value of f_upper to use when getting the mu coordinates from the
        lambda coordinates. This must be a key in metric_params.evals and
        metric_params.evecs (ie. we must know how to do the transformation for
        the given value of f_upper)

    Returns
    -------
    mus : list of floats or numpy.arrays
        Position of the system(s) in the mu coordinate system

    """
    # Do this by masses -> lambdas
    lambdas = get_chirp_params(
        mass1,
        mass2,
        spin1z,
        spin2z,
        metric_params.f0,
        metric_params.pn_order,
        lambda1=lambda1,
        lambda2=lambda2,
        quadparam1=quadparam1,
        quadparam2=quadparam2,
    )
    # and lambdas -> mus
    mus = get_mu_params(lambdas, metric_params, f_upper)
    return mus


def get_mu_params(lambdas, metric_params, f_upper):
    """
    Function to rotate from the lambda coefficients into position in the mu
    coordinate system. Mu = Cartesian metric, but not principal components.

    Parameters
    ----------
    lambdas : list of floats or numpy.arrays
        Position of the system(s) in the lambda coefficients
    metric_params : MetricParameters instance
        Structure holding all the options for construction of the metric
        and the eigenvalues, eigenvectors and covariance matrix
        needed to manipulate the space.
    f_upper : float
        The value of f_upper to use when getting the mu coordinates from the
        lambda coordinates. This must be a key in metric_params.evals and
        metric_params.evecs (ie. we must know how to do the transformation for
        the given value of f_upper)

    Returns
    -------
    mus : list of floats or numpy.arrays
        Position of the system(s) in the mu coordinate system

    """
    lambdas = numpy.asarray(lambdas)
    # If original inputs were floats we need to make this a 2D array
    if len(lambdas.shape) == 1:
        resize_needed = True
        lambdas = lambdas[:, None]
    else:
        resize_needed = False

    evecs = metric_params.evecs[f_upper]
    evals = metric_params.evals[f_upper]

    evecs = numpy.asarray(evecs)

    mus = ((lambdas.T).dot(evecs)).T
    mus = mus * numpy.sqrt(evals)[:, None]

    if resize_needed:
        mus = numpy.ndarray.flatten(mus)

    return mus


def get_covaried_params(mus, evecs_cv):
    """
    Function to rotate from position(s) in the mu_i coordinate system into the
    position(s) in the xi_i coordinate system

    Parameters
    ----------
    mus : list of floats or numpy.arrays
        Position of the system(s) in the mu coordinate system
    evecs_cv : numpy.matrix
        This matrix is used to perform the rotation to the xi_i
        coordinate system.

    Returns
    -------
    xis : list of floats or numpy.arrays
        Position of the system(s) in the xi coordinate system

    """
    mus = numpy.asarray(mus)
    # If original inputs were floats we need to make this a 2D array
    if len(mus.shape) == 1:
        resize_needed = True
        mus = mus[:, None]
    else:
        resize_needed = False

    xis = ((mus.T).dot(evecs_cv)).T

    if resize_needed:
        xis = numpy.ndarray.flatten(xis)

    return xis


def rotate_vector(evecs, old_vector, rescale_factor, index):
    """
    Function to find the position of the system(s) in one of the xi_i or mu_i
    directions.

    Parameters
    ----------
    evecs : numpy.matrix
        Matrix of the eigenvectors of the metric in lambda_i coordinates. Used
        to rotate to a Cartesian coordinate system.
    old_vector : list of floats or numpy.arrays
        The position of the system(s) in the original coordinates
    rescale_factor : float
        Scaling factor to apply to resulting position(s)
    index : int
        The index of the final coordinate system that is being computed. Ie.
        if we are going from mu_i -> xi_j, this will give j.

    Returns
    -------
    positions : float or numpy.array
        Position of the point(s) in the resulting coordinate.

    """
    temp = 0
    for i in range(len(evecs)):
        temp += (evecs[i, index] * rescale_factor) * old_vector[i]
    return temp


def get_point_distance(point1, point2, metric_params, f_upper):
    """
    Function to calculate the mismatch between two points, supplied in terms
    of the masses and spins, using the xi_i parameter space metric to
    approximate the mismatch of the two points. Can also take one of the points
    as an array of points and return an array of mismatches (but only one can
    be an array!)

    point1 : List of floats or numpy.arrays
        point1[0] contains the mass(es) of the heaviest body(ies).
        point1[1] contains the mass(es) of the smallest body(ies).
        point1[2] contains the spin(es) of the heaviest body(ies).
        point1[3] contains the spin(es) of the smallest body(ies).
    point2 : List of floats
        point2[0] contains the mass of the heaviest body.
        point2[1] contains the mass of the smallest body.
        point2[2] contains the spin of the heaviest body.
        point2[3] contains the spin of the smallest body.
    metric_params : MetricParameters instance
        Structure holding all the options for construction of the metric
        and the eigenvalues, eigenvectors and covariance matrix
        needed to manipulate the space.
    f_upper : float
        The value of f_upper to use when getting the mu coordinates from the
        lambda coordinates. This must be a key in metric_params.evals,
        metric_params.evecs and metric_params.evecs_cv
        (ie. we must know how to do the transformation for
        the given value of f_upper)

    Returns
    -------
    dist : float or numpy.array
        Distance between the point2 and all points in point1
    xis1 : List of floats or numpy.arrays
        Position of the input point1(s) in the xi_i parameter space
    xis2 : List of floats
        Position of the input point2 in the xi_i parameter space

    """
    a_mass1 = point1[0]
    a_mass2 = point1[1]
    a_spin1 = point1[2]
    a_spin2 = point1[3]

    b_mass1 = point2[0]
    b_mass2 = point2[1]
    b_spin1 = point2[2]
    b_spin2 = point2[3]

    a_xis = get_cov_params(a_mass1, a_mass2, a_spin1, a_spin2, metric_params, f_upper)

    b_xis = get_cov_params(b_mass1, b_mass2, b_spin1, b_spin2, metric_params, f_upper)

    dist = (a_xis[0] - b_xis[0]) ** 2
    for i in range(1, len(a_xis)):
        dist += (a_xis[i] - b_xis[i]) ** 2

    return dist, a_xis, b_xis


def calc_point_dist(vs, entry):
    r"""
    This function is used to determine the distance between two points.

    Parameters
    ----------
    vs : list or numpy.array or similar
        An array of point 1's position in the \chi_i coordinate system
    entry : list or numpy.array or similar
        An array of point 2's position in the \chi_i coordinate system
    minimal_mismatch : float
        The minimal mismatch allowed between the points

    Returns
    -------
    val : float
        The metric distance between the two points.

    """
    chi_diffs = vs - entry
    val = ((chi_diffs) * (chi_diffs)).sum()
    return val


def test_point_dist(point_1_chis, point_2_chis, distance_threshold):
    r"""
    This function tests if the difference between two points in the chi
    parameter space is less than a distance threshold. Returns True if it is
    and False if it is not.

    Parameters
    ----------
    point_1_chis : numpy.array
        An array of point 1's position in the \chi_i coordinate system
    point_2_chis : numpy.array
        An array of point 2's position in the \chi_i coordinate system
    distance_threshold : float
        The distance threshold to use.

    """
    return calc_point_dist(point_1_chis, point_2_chis) < distance_threshold


def calc_point_dist_vary(
    mus1, f_upper1, mus2, f_upper2, f_map, norm_map, minimal_mismatch
):
    r"""
    Function to determine if two points, with differing upper frequency cutoffs
    have a mismatch < minimal_mismatch for *both* upper frequency cutoffs.

    Parameters
    ----------
    mus1 : List of numpy arrays
        mus1[i] will give the array of point 1's position in the \chi_j
        coordinate system. The i element corresponds to varying values of the
        upper frequency cutoff. f_map is used to map between i and actual
        frequencies
    f_upper1 : float
        The upper frequency cutoff of point 1.
    mus2 : List of numpy arrays
        mus2[i] will give the array of point 2's position in the \chi_j
        coordinate system. The i element corresponds to varying values of the
        upper frequency cutoff. f_map is used to map between i and actual
        frequencies
    f_upper2 : float
        The upper frequency cutoff of point 2.
    f_map : dictionary
        f_map[f_upper] will give the index needed to get the \chi_j coordinates
        in the two sets of mus
    norm_map : dictionary
        norm_map[f_upper] will give the relative frequency domain template
        amplitude (sigma) at the given value of f_upper.
    minimal_mismatch
        The minimal mismatch allowed between the points

    Returns
    -------
    Boolean
        True if the points have a mismatch < minimal_mismatch
        False if the points have a mismatch > minimal_mismatch

    """
    f_upper = min(f_upper1, f_upper2)
    f_other = max(f_upper1, f_upper2)
    idx = f_map[f_upper]
    vecs1 = mus1[idx]
    vecs2 = mus2[idx]
    val = ((vecs1 - vecs2) * (vecs1 - vecs2)).sum()
    if val > minimal_mismatch:
        return False
    # Reduce match to account for normalization.
    norm_fac = norm_map[f_upper] / norm_map[f_other]
    val = 1 - (1 - val) * norm_fac
    return val < minimal_mismatch


def find_max_and_min_frequencies(name, mass_range_params, freqs):
    """
    ADD DOCS
    """
    cutoff_fns = pnutils.named_frequency_cutoffs
    if name not in cutoff_fns.keys():
        err_msg = f"{name} not recognized as a valid cutoff frequency choice."
        err_msg += "Recognized choices: " + " ".join(cutoff_fns.keys())
        raise ValueError(err_msg)

    # Can I do this quickly?
    total_mass_approxs = {
        "SchwarzISCO": pnutils.f_schwarz_isco,
        "LightRing": pnutils.f_light_ring,
        "ERD": pnutils.f_erd,
    }

    if name in total_mass_approxs:
        # This can be done quickly if the cutoff only depends on total mass
        # Assumes that lower total mass = higher cutoff frequency
        upper_f_cutoff = total_mass_approxs[name](mass_range_params.min_tot_mass)
        lower_f_cutoff = total_mass_approxs[name](mass_range_params.max_tot_mass)
    else:
        # Do this numerically
        # FIXME: Is 1000000 the right choice? I think so, but just highlighting
        mass1, mass2, spin1z, spin2z = get_random_mass(1000000, mass_range_params)
        mass_dict = {}
        mass_dict["mass1"] = mass1
        mass_dict["mass2"] = mass2
        mass_dict["spin1z"] = spin1z
        mass_dict["spin2z"] = spin2z
        tmp_freqs = cutoff_fns[name](mass_dict)
        upper_f_cutoff = tmp_freqs.max()
        lower_f_cutoff = tmp_freqs.min()

    cutoffs = numpy.array([lower_f_cutoff, upper_f_cutoff])
    if lower_f_cutoff < freqs.min():
        warn_msg = "WARNING: "
        warn_msg += f"Lowest frequency cutoff is {lower_f_cutoff} Hz "
        warn_msg += "which is lower than the lowest frequency calculated "
        warn_msg += f"for the metric: {freqs.min()} Hz. "
        warn_msg += "Distances for these waveforms will be calculated at "
        warn_msg += "the lowest available metric frequency."
        logger.warning(warn_msg)
    if upper_f_cutoff > freqs.max():
        warn_msg = "WARNING: "
        warn_msg += f"Highest frequency cutoff is {upper_f_cutoff} Hz "
        warn_msg += "which is larger than the highest frequency calculated "
        warn_msg += f"for the metric: {freqs.max()} Hz. "
        warn_msg += "Distances for these waveforms will be calculated at "
        warn_msg += "the largest available metric frequency."
        logger.warning(warn_msg)
    return find_closest_calculated_frequencies(cutoffs, freqs)


def return_nearest_cutoff(name, mass_dict, freqs):
    """
    Given an array of total mass values and an (ascending) list of
    frequencies, this will calculate the specified cutoff formula for each
    mtotal and return the nearest frequency to each cutoff from the input
    list.
    Currently only supports cutoffs that are functions of the total mass
    and no other parameters (SchwarzISCO, LightRing, ERD)

    Parameters
    ----------
    name : string
        Name of the cutoff formula to be approximated
    mass_dict : Dictionary where the keys are used to call the functions
        returned by tmpltbank.named_frequency_cutoffs. The values can be
        numpy arrays or single values.
    freqs : list of floats
        A list of frequencies (must be sorted ascending)

    Returns
    -------
    numpy.array
        The frequencies closest to the cutoff for each value of totmass.

    """
    # A bypass for the redundant case
    if len(freqs) == 1:
        return numpy.zeros(len(mass_dict["m1"]), dtype=float) + freqs[0]
    cutoff_fns = pnutils.named_frequency_cutoffs
    if name not in cutoff_fns.keys():
        err_msg = f"{name} not recognized as a valid cutoff frequency choice."
        err_msg += "Recognized choices: " + " ".join(cutoff_fns.keys())
        raise ValueError(err_msg)
    f_cutoff = cutoff_fns[name](mass_dict)
    return find_closest_calculated_frequencies(f_cutoff, freqs)


def find_closest_calculated_frequencies(input_freqs, metric_freqs):
    """
    Given a value (or array) of input frequencies find the closest values in
    the list of frequencies calculated in the metric.

    Parameters
    ----------
    input_freqs : numpy.array or float
        The frequency(ies) that you want to find the closest value in
        metric_freqs
    metric_freqs : numpy.array
        The list of frequencies calculated by the metric

    Returns
    -------
    output_freqs : numpy.array or float
        The list of closest values to input_freqs for which the metric was
        computed

    """
    try:
        ref_ev = numpy.zeros(len(input_freqs), dtype=float)
    except TypeError:
        ref_ev = numpy.zeros(1, dtype=float)
        input_freqs = numpy.array([input_freqs])

    if len(metric_freqs) == 1:
        ref_ev[:] = metric_freqs[0]
        return ref_ev

    # FIXME: This seems complicated for what is a simple operation. Is there
    #        a simpler *and* faster way of doing this?
    # NOTE: This function assumes a sorted list of frequencies
    # NOTE: totmass and f_cutoff are both numpy arrays as this function is
    #       designed so that the cutoff can be calculated for many systems
    #       simulataneously
    for i in range(len(metric_freqs)):
        if i == 0:
            # If frequency is lower than halfway between the first two entries
            # use the first (lowest) value
            logic_arr = input_freqs < ((metric_freqs[0] + metric_freqs[1]) / 2.0)
        elif i == (len(metric_freqs) - 1):
            # If frequency is larger than halfway between the last two entries
            # use the last (highest) value
            logic_arr = input_freqs > ((metric_freqs[-2] + metric_freqs[-1]) / 2.0)
        else:
            # For frequencies within the range in freqs, check which points
            # should use the frequency corresponding to index i.
            logic_arr_a = input_freqs > ((metric_freqs[i - 1] + metric_freqs[i]) / 2.0)
            logic_arr_b = input_freqs < ((metric_freqs[i] + metric_freqs[i + 1]) / 2.0)
            logic_arr = numpy.logical_and(logic_arr_a, logic_arr_b)
        if logic_arr.any():
            ref_ev[logic_arr] = metric_freqs[i]
    return ref_ev


def outspiral_loop(max_dist):
    """
    Return a list of points that will loop outwards in a 2D lattice in terms
    of distance from a central point. So if max_dist=2 this will be [0,0], [0,1],
    [0,-1],[1,0],[-1,0],[1,1] .... This is useful when you want to loop over
    a number of bins, but want to start in the center and work outwards.
    """
    # Create a 2D lattice of all points
    x_vals, y_vals = numpy.meshgrid(
        numpy.arange(-max_dist, max_dist + 1), numpy.arange(-max_dist, max_dist + 1)
    )

    # Flatten it
    x_vals = numpy.ndarray.flatten(x_vals)
    y_vals = numpy.ndarray.flatten(y_vals)

    # Force to an integer
    x_vals = numpy.array(x_vals, dtype=int)
    y_vals = numpy.array(y_vals, dtype=int)

    # Calculate distances
    dist = numpy.sqrt(x_vals**2 + y_vals**2)

    # Combine back into an array
    out_arr = numpy.array([x_vals, y_vals, dist])

    # And order correctly
    sorted_out_arr = out_arr[:, out_arr[2].argsort()]

    return sorted_out_arr[:2, :].T
