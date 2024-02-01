# Copyright (C) 2023 Arthur Tolley, Gareth Cabourn Davies
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

#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#
""" This module contains functions for optimizing the signal-to-noise ratio
of triggers produced by PyCBC Live. Also contained within this module are the
command line arguments required and options group for the SNR optimization.
This module is primarily used in the pycbc_optimize_snr program.
"""

import time
import logging
import types
import numpy
from scipy.optimize import differential_evolution, shgo
from pycbc import (
    DYN_RANGE_FAC, waveform
)
from pycbc.types import zeros
import pycbc.waveform.bank
from pycbc.filter import matched_filter_core
import pycbc.conversions as cv

try:
    import pyswarms as ps
except:
    ps = None

# Set a minimum mass for points tried in optimization allowing for
# minimal slop relative to the lightest template
MIN_CPT_MASS = 0.99

# Set a large maximum total mass
MAX_MTOTAL = 500.

Nfeval = 0
start_time = time.time()


def callback_func(Xi, convergence=0):
    global Nfeval
    logging.info("Currently at %d %s", Nfeval, convergence)
    # Time out if the optimization takes longer than 6 minutes
    if (time.time() - start_time) > 360:
        return True
    Nfeval += 1


def compute_network_snr_core(v, data, coinc_times, ifos, flen, approximant,
                             flow, f_end, delta_f, sample_rate, raise_err=False):
    """
    Compute network SNR as a function over mchirp, eta, and two aligned
    spin components, stored in that order in the sequence v.

    Parameters
    ----------
    v : list
        A list containing the input values for mchirp, eta, and spin
        components.
    data : dict
        A dictionary containing keys of ifos ('H1', 'L1') and
        values of the frequency series data for those ifos
    coinc_times : dict
        A dictionary containing the coincidence times for the network.
    ifos : list
        A list of the ifos, e.g. ['H1', 'L1']
    flen : float
        The length of the data.
    approximant : str
        The approximant used for the waveform model.
    flow : float
        The lower frequency bound.
    f_end : float
        The upper frequency bound.
    delta_f : float
        The frequency spacing.
    sample_rate : float
        The sampling rate of the data.
    raise_err : bool, optional
        A flag indicating whether to raise an error if an exception
        occurs during the computation. Defaults to False.

    Returns
    -------
    network_snr : float
        The computed network SNR (Signal-to-Noise Ratio) value.
    snr_series_dict : dict
        A dictionary containing the snr timeseries from each ifo.

    """
    distance = 1.0 / DYN_RANGE_FAC
    mtotal = cv.mtotal_from_mchirp_eta(v[0], v[1])
    mass1 = cv.mass1_from_mtotal_eta(mtotal, v[1])
    mass2 = cv.mass2_from_mtotal_eta(mtotal, v[1])

    # enforce broadly accepted search space boundaries
    if mass1 < MIN_CPT_MASS or mass2 < MIN_CPT_MASS or mtotal > MAX_MTOTAL:
        return -numpy.inf, {}

    try:
        htilde = waveform.get_waveform_filter(
                zeros(flen, dtype=numpy.complex64),
                approximant=approximant,
                mass1=mass1, mass2=mass2, spin1z=v[2], spin2z=v[3],
                f_lower=flow, f_final=f_end, delta_f=delta_f,
                delta_t=1./sample_rate, distance=distance)
    except RuntimeError:
        if raise_err:
            raise
        # Assume a failure in the waveform approximant
        # due to the choice of parameters and carry on
        return -numpy.inf, {}

    if not hasattr(htilde, 'params'):
        htilde.params = dict(mass1=mass1, mass2=mass2,
                             spin1z=v[2], spin2z=v[3])
    if not hasattr(htilde, 'end_idx'):
        htilde.end_idx = int(f_end / htilde.delta_f)
    htilde.approximant = approximant
    htilde.sigmasq = types.MethodType(pycbc.waveform.bank.sigma_cached,
                                      htilde)
    htilde.min_f_lower = flow
    htilde.end_frequency = f_end
    htilde.f_lower = flow
    network_snrsq = 0
    snr_series_dict = {}
    for ifo in ifos:
        sigmasq = htilde.sigmasq(data[ifo].psd)
        snr, _, norm = matched_filter_core(htilde, data[ifo],
                                           h_norm=sigmasq)
        duration = 0.095
        half_dur_samples = int(snr.sample_rate * duration / 2)
        onsource_idx = (float(coinc_times[ifo] - snr.start_time) *
                        snr.sample_rate)
        onsource_idx = int(round(onsource_idx))
        onsource_slice = slice(onsource_idx - half_dur_samples,
                               onsource_idx + half_dur_samples + 1)
        snr_series = snr[onsource_slice] * norm
        snr_series_dict[ifo] = snr * norm
        snr_series_dict['sigmasq_' + ifo] = sigmasq
        network_snrsq += max(abs(snr_series._data)) ** 2.

    return network_snrsq ** 0.5, snr_series_dict


def compute_minus_network_snr(v, *argv):
    if len(argv) == 1:
        argv = argv[0]
    nsnr, _ = compute_network_snr_core(v, *argv)
    logging.debug('snr: %s', nsnr)
    return -nsnr


def compute_minus_network_snr_pso(v, *argv, **kwargs):
    argv = kwargs['args']
    nsnr_array = numpy.array([
        compute_network_snr_core(v_i, *argv)[0]
        for v_i in v])
    return -nsnr_array


def optimize_di(bounds, cli_args, extra_args, initial_point):
    # Convert from dict to array with parameters in a given order
    bounds = numpy.array([
        bounds['mchirp'],
        bounds['eta'],
        bounds['spin1z'],
        bounds['spin2z']
    ])
    # Initialize the population with random values within specified bounds
    population = numpy.random.uniform(
        bounds[:, 0],
        bounds[:, 1],
        size=(int(cli_args.snr_opt_di_popsize), len(bounds))
    )
    if cli_args.snr_opt_include_candidate:
        # add the initial point to the population
        population = numpy.concatenate((population[:-1],
                                        initial_point))
    logging.debug('Initial population: %s', population)

    results = differential_evolution(
        compute_minus_network_snr,
        bounds,
        maxiter=int(cli_args.snr_opt_di_maxiter),
        workers=(cli_args.cores or -1),
        popsize=int(cli_args.snr_opt_di_popsize),
        mutation=(0.5, 1),
        recombination=0.7,
        callback=callback_func,
        args=extra_args,
        init=population
    )
    return results.x


def optimize_shgo(bounds, cli_args, extra_args, initial_point): # pylint: disable=unused-argument
    bounds = [
        bounds['mchirp'],
        bounds['eta'],
        bounds['spin1z'],
        bounds['spin2z']
    ]
    results = shgo(
        compute_minus_network_snr,
        bounds=bounds,
        args=extra_args,
        iters=cli_args.snr_opt_shgo_iters,
        n=cli_args.snr_opt_shgo_samples,
        sampling_method="sobol"
    )
    return results.x


def optimize_pso(bounds, cli_args, extra_args, initial_point):
    options = {
        'c1': cli_args.snr_opt_pso_c1,
        'c2': cli_args.snr_opt_pso_c2,
        'w': cli_args.snr_opt_pso_w
    }
    min_bounds = numpy.array([
        bounds['mchirp'][0],
        bounds['eta'][0],
        bounds['spin1z'][0],
        bounds['spin2z'][0]
    ])
    max_bounds = numpy.array([
        bounds['mchirp'][1],
        bounds['eta'][1],
        bounds['spin1z'][1],
        bounds['spin2z'][1]
    ])

    # Initialize the population with random values within specified bounds
    population = numpy.random.uniform(
        min_bounds,
        max_bounds,
        size=(int(cli_args.snr_opt_pso_particles), len(bounds))
    )

    if cli_args.snr_opt_include_candidate:
        # add the initial point to the population
        population = numpy.concatenate((population[:-1],
                                        initial_point))
    logging.debug('Initial population: %s', population)

    optimizer = ps.single.GlobalBestPSO(
        n_particles=int(cli_args.snr_opt_pso_particles),
        dimensions=4,
        options=options,
        bounds=(min_bounds, max_bounds),
        init_pos=population
    )
    _, results = optimizer.optimize(
        compute_minus_network_snr_pso,
        iters=int(cli_args.snr_opt_pso_iters),
        n_processes=cli_args.cores,
        args=extra_args
    )
    return results


optimize_funcs = {
    'differential_evolution': optimize_di,
    'shgo': optimize_shgo,
    'pso': optimize_pso
}

# The following sets the default values of the options, but allows us to check
# if the option has been given on the command line

# For each optimizer, we have a dictionary of the options, its help
# message and default value

option_dict = {
    'differential_evolution': {
        'maxiter': ('The maximum number of generations over which the entire '
                    'population is evolved.', 50),
        'popsize': ('A multiplier for setting the total population size.',
                    100),
    },
    'shgo': {
        'samples': ('Number of sampling points used in the construction of '
                    'the simplicial complex.', 76),
        'iters': ('Number of iterations used in the construction of the '
                  'simplicial complex.', 3),
    },
    'pso': {
        'iters': ('Number of iterations used in the particle swarm '
                  'optimization.', 5),
        'particles': ('Number of particles used in the swarm.', 250),
        'c1': ('The hyperparameter c1: the cognitive parameter.', 0.5),
        'c2': ('The hyperparameter c2: the social parameter.', 2.0),
        'w': ('The hyperparameter w: the inertia parameter.', 0.01),
    }
}

def insert_snr_optimizer_options(parser):
    opt_opt_group = parser.add_argument_group("SNR optimizer configuration "
                                              "options.")
    # Option to choose which optimizer to use:
    optimizer_choices = sorted(list(option_dict.keys()))
    opt_opt_group.add_argument('--snr-opt-method',
        default='differential_evolution',
        choices=optimizer_choices,
        help='SNR Optimizer choices: ' + ', '.join(optimizer_choices))

    # Add the generic options
    opt_opt_group.add_argument('--snr-opt-include-candidate',
        action='store_true',
        help='Include parameters of the candidate event in the initialized '
             'array for the optimizer. Only relevant for --optimizer pso or '
             'differential_evolution')
    opt_opt_group.add_argument('--snr-opt-seed',
        default='42',
        help='Seed to supply to the random generation of initial array to '
             'pass to the optimizer. Only relevant for --optimizer pso or '
             'differential_evolution. Set to ''random'' for a random seed')

    # For each optimizer, add the possible options
    for optimizer, option_subdict in option_dict.items():
        optimizer_name = optimizer.replace('_', '-')
        if optimizer_name == 'differential-evolution':
            optimizer_name = 'di'
        for opt_name, opt_help_default in option_subdict.items():
            option_name = f"--snr-opt-{optimizer_name}-{opt_name}"
            opt_opt_group.add_argument(option_name,
                type=float,
                help=f'Only relevant for --optimizer {optimizer}: ' +
                     opt_help_default[0] +
                     f' Default = {opt_help_default[1]}')


def check_snr_optimizer_options(args, parser):
    """
    Deal with default options and required parameters given optimizer option
    """
    options = {}
    options['differential_evolution'] = [args.snr_opt_di_maxiter,
                                         args.snr_opt_di_popsize]
    options['shgo'] = [args.snr_opt_shgo_samples, args.snr_opt_shgo_iters]
    options['pso'] = [args.snr_opt_pso_iters, args.snr_opt_pso_particles,
                      args.snr_opt_pso_c1, args.snr_opt_pso_c2,
                      args.snr_opt_pso_w]

    if args.snr_opt_method == 'pso' and ps is None:
        parser.error('You need to install pyswarms to use the pso optimizer.')

    # Check all the options are suitable for the chosen optimizer
    for k in options.keys():
        if args.snr_opt_method == k:
            continue
        if any(options[k]):
            parser.error("Argument has been supplied which is not suitable " +
                         f"for the optimizer given ({args.snr_opt_method})")

    # Give the arguments the default values according to the dictionary
    optimizer_name = args.snr_opt_method.replace('_', '-')
    if optimizer_name == 'differential-evolution':
        optimizer_name = 'di'
    for key, value in option_dict[args.snr_opt_method].items():
        key_name = f'snr_opt_{optimizer_name}_{key}'
        if not getattr(args, key_name):
            setattr(args, key_name, value[1])


def args_to_string(args):
    """
    Convert the supplied arguments for SNR optimization config into
    a string - this is to be used when running subprocesses
    """
    argstr = f'--snr-opt-method {args.snr_opt_method} '
    optimizer_name = args.snr_opt_method.replace('_', '-')
    if optimizer_name == 'differential-evolution':
        optimizer_name = 'di'
    for opt in option_dict[args.snr_opt_method]:
        key_name = f'snr_opt_{optimizer_name}_{opt}'
        option_value = getattr(args, key_name)
        # If the option is not given, don't pass it and use default
        if option_value is None:
            continue
        option_fullname = f'--snr-opt-{optimizer_name}-{opt}'
        argstr += f'{option_fullname} {option_value} '

    return argstr
