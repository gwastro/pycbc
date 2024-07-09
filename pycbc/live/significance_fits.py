"""
Functions for defining the live significance fits
"""

import logging
import h5py
import numpy

logger = logging.getLogger('pycbc.live.single_fits')


def add_live_significance_trigger_pruning_options(parser):
    """
    Add options used for pruning in live singles significance fits
    """
    pruning_group = parser.add_argument_group("Trigger pruning")
    pruning_group.add_argument(
        "--prune-loudest",
        type=int,
        help="Maximum number of loudest trigger clusters to "
             "remove from each bin."
    )
    pruning_group.add_argument(
        "--prune-window",
        type=float,
        help="Window (seconds) either side of the --prune-loudest "
             "loudest triggers in each duration bin to remove."
    )
    pruning_group.add_argument(
        "--prune-stat-threshold",
        type=float,
        help="Minimum statistic value to consider a "
             "trigger for pruning."
    )


def verify_live_significance_trigger_pruning_options(args, parser):
    """
    Verify options used for pruning in live singles significance fits
    """
    # Pruning options are mutually required or not needed
    prune_options = [args.prune_loudest, args.prune_window,
                     args.prune_stat_threshold]

    if any(prune_options) and not all(prune_options):
        parser.error("Require all or none of --prune-loudest, "
                     "--prune-window and --prune-stat-threshold")


def add_live_significance_duration_bin_options(parser):
    """
    Add options used to calculate duration bin edges in live
    singles significance fits
    """
    durbin_group = parser.add_argument_group('Duration Bins')
    durbin_group.add_argument(
        "--duration-bin-edges",
        nargs='+',
        type=float,
        help="Durations to use for bin edges. "
             "Use if specifying exact bin edges, "
             "Not compatible with --duration-bin-start, "
             "--duration-bin-end and --num-duration-bins"
    )
    durbin_group.add_argument(
        "--duration-bin-start",
        type=float,
        help="Shortest duration to use for duration bins."
             "Not compatible with --duration-bins, requires "
             "--duration-bin-end and --num-duration-bins."
    )
    durbin_group.add_argument(
        "--duration-bin-end", type=float,
        help="Longest duration to use for duration bins."
    )
    durbin_group.add_argument(
        "--duration-from-bank",
        help="Path to the template bank file to get max/min "
             "durations from."
    )
    durbin_group.add_argument(
        "--num-duration-bins",
        type=int,
        help="How many template duration bins to split the bank "
             "into before fitting."
    )
    durbin_group.add_argument(
        "--duration-bin-spacing",
        choices=['linear', 'log'],
        default='log',
        help="How to set spacing for bank split "
             "if using --num-duration-bins and "
             "--duration-bin-start + --duration-bin-end "
             "or --duration-from-bank."
    )


def verify_live_significance_duration_bin_options(args, parser):
    """
    Verify options used to calculate duration bin edges in live
    singles significance fits
    """
    # Check the bin options
    if args.duration_bin_edges:
        if (args.duration_bin_start or args.duration_bin_end or
                args.duration_from_bank or args.num_duration_bins):
            parser.error("Cannot use --duration-bin-edges with "
                         "--duration-bin-start, --duration-bin-end, "
                         "--duration-from-bank or --num-duration-bins.")
    else:
        if not args.num_duration_bins:
            parser.error("--num-duration-bins must be set if not using "
                         "--duration-bin-edges.")
        if not ((args.duration_bin_start and args.duration_bin_end) or
                args.duration_from_bank):
            parser.error("--duration-bin-start & --duration-bin-end or "
                         "--duration-from-bank must be set if not using "
                         "--duration-bin-edges.")
    if args.duration_bin_end and \
            args.duration_bin_end <= args.duration_bin_start:
        parser.error("--duration-bin-end must be greater than "
                     "--duration-bin-start, got "
                     f"{args.duration_bin_end} and {args.duration_bin_start}")


def duration_bins_from_cli(args):
    """Create the duration bins from CLI options.
    """
    if args.duration_bin_edges:
        # direct bin specification
        return numpy.array(args.duration_bin_edges)
    # calculate bins from min/max and number
    min_dur = args.duration_bin_start
    max_dur = args.duration_bin_end
    if args.duration_from_bank:
        # read min/max duration directly from the bank itself
        with h5py.File(args.duration_from_bank, 'r') as bank_file:
            temp_durs = bank_file['template_duration'][:]
        min_dur, max_dur = min(temp_durs), max(temp_durs)
    if args.duration_bin_spacing == 'log':
        return numpy.logspace(
            numpy.log10(min_dur),
            numpy.log10(max_dur),
            args.num_duration_bins + 1
        )
    if args.duration_bin_spacing == 'linear':
        return numpy.linspace(
            min_dur,
            max_dur,
            args.num_duration_bins + 1
        )
    raise RuntimeError("Invalid duration bin specification")


__all__ = [
    'add_live_significance_trigger_pruning_options',
    'verify_live_significance_trigger_pruning_options',
    'add_live_significance_duration_bin_options',
    'verify_live_significance_duration_bin_options',
    'duration_bins_from_cli',
]
