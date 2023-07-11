import argparse, numpy
import time
import logging
import types
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
    if (time.time() - start_time) > 360:
        return True
    Nfeval += 1


def compute_network_snr_core(v, *argv, raise_err=False):
    """
    Compute network SNR as a function over mchirp, eta and two aligned
    spin components, stored in that order in the sequence v
    """
    data = argv[0]
    coinc_times = argv[1]
    ifos = argv[2]
    flen = argv[3]
    approximant = argv[4]
    flow = argv[5]
    f_end = argv[6]
    delta_f = argv[7]
    sample_rate = argv[8]
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
        # assume a failure in the waveform approximant
        # due to the choice of parameters
        else:
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
        onsource_idx = float(coinc_times[ifo] - snr.start_time) * snr.sample_rate
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
    return -nsnr


def compute_minus_network_snr_pso(v, *argv, **kwargs):
    argv = kwargs['args']
    nsnr_array = numpy.array([compute_network_snr_core(v_i, *argv)[0] for v_i in v])
    return -nsnr_array


def normalize_initial_point(initial_point, bounds):
    return (initial_point - bounds[:,0]) / (bounds[:,1] - bounds[:,0])

def optimize_di(bounds, cli_args, extra_args, initial_point):
    bounds = numpy.array([
        bounds['mchirp'],
        bounds['eta'],
        bounds['spin1z'],
        bounds['spin2z']
    ])


    # Currently only implemented for random seed initial array
    rng = numpy.random.mtrand._rand
    population_shape=(cli_args.di_popsize, 4)
    population = rng.uniform(size=population_shape)
    if cli_args.include_candidate_in_optimizer:
        # Re-normalize the initial point into the correct range
        point_init = normalize_initial_point(initial_point, bounds)
        # add the initial point to the population
        population = numpy.concatenate((population[:-1], point_init))

    results = differential_evolution(
        compute_minus_network_snr,
        bounds,
        maxiter=cli_args.di_maxiter,
        workers=(cli_args.cores or -1),
        popsize=cli_args.di_popsize,
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
        iters=args.shgo_iters,
        n=args.shgo_samples,
        sampling_method="sobol"
    )
    return results.x

def normalize_population(population, min_bounds, max_bounds):
    norm_pop = min_bounds + population * (max_bounds - min_bounds)

    return norm_pop

def optimize_pso(bounds, cli_args, extra_args, initial_point):
    options = {
        'c1': cli_args.pso_c1,
        'c2': cli_args.pso_c2,
        'w': cli_args.pso_w
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

    # Manually generate the initial points, this is the same as the default
    # method, but allows us to make some modifications
    population = numpy.random.uniform(
        low=0.0, high=1.0, size=(cli_args.pso_particles, 4)
    )
    population = normalize_population(population, min_bounds, max_bounds)

    if cli_args.include_candidate_in_optimizer:
        # add the initial point to the population
        population = numpy.concatenate((population[:-1],
                                        initial_point))

    optimizer = ps.single.GlobalBestPSO(
        n_particles=cli_args.pso_particles,
        dimensions=4,
        options=options,
        bounds=(min_bounds, max_bounds),
        init_pos=population
    )
    _, results = optimizer.optimize(
        compute_minus_network_snr_pso,
        iters=cli_args.pso_iters,
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

# TODO: Add the options for the initial array
#parser.add_argument('--include-candidate-in-optimizer', action='store_true',
#             help='Include parameters of the candidate event in the '
#             'initialised array for the optimizer. Only relevant for '
#             '--optimizer pso or differential_evolution')
# parser.add_argument('--seed', type=int, default=42,
#             help='Seed to supply to the random generation of initial '
#             'array to pass to the optimizer. Only relevant for '
#             '--optimizer pso or differential_evolution')

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
    },
}

def insert_snr_optimizer_options(parser):
    opt_opt_group = parser.add_argument_group(
        "SNR optimizer "
        "configuration options."
    )
    # Option to choose which optimizer to use:
    optimizer_choices = sorted(list(option_dict.keys()))
    opt_opt_group.add_argument('--optimizer',
        type=str,
        default='differential_evolution',
        choices=optimizer_choices,
        help='The optimizer to use, choose from '
              ', '.join(optimizer_choices))

    # For each optimizer, add the possible options
    for optimizer, option_subdict in option_dict.items():
        optimizer_name = optimizer.replace('_','-')
        for opt_name, opt_help_default in option_subdict.items():
            option_name = f"--{optimizer_name}-{opt_name}"
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
    options['differential_evolution'] = [args.differential_evolution_maxiter,
                                         args.differential_evolution_popsize]
    options['shgo'] = [args.shgo_samples, args.shgo_iters]
    options['pso'] = [args.pso_particles, args.pso_c1, args.pso_c2, args.pso_w]

    if args.optimizer == 'pso' and ps == None:
        parser.error('You need to install pyswarms to use the pso optimizer.')

    # Have any options been given for a different optimizer?
    # If so, raise a warning
    for k in options.keys():
        if args.optimizer == k: continue
        if any(options[k]):
            parser.error("Argument has been supplied which is not suitable " +
                         f"for the optimizer given ({args.optimizer})")

    # Give the arguments the default values according to the dictionary
    for key, value in option_dict[args.optimizer].items():
        key_name = f'{args.optimizer}_{key}'
        if not getattr(args, key_name):
            setattr(args, key_name, value[1])

def args_to_string(args):
    """
    Convert the supplied arguments for SNR optimization config into
    a string - this is to be used when running subprocesses
    """
    argstr = f' --optimizer {args.optimizer} '
    for opt in option_dict[args.optimizer]:
        option_fullname = f"--{args.optimizer}-{opt}"
        key_name = f'{args.optimizer}_{opt}'
        option_value = getattr(args, key_name)
        argstr += f'{option_fullname} {option_value} '
    return argstr
