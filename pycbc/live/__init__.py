"""
This packages contains modules to help with pycbc live running
"""

from .significance_fits import (
    add_live_significance_duration_bin_options,
    add_live_significance_trigger_pruning_options,
    duration_bins_from_cli,
    verify_live_significance_duration_bin_options,
    verify_live_significance_trigger_pruning_options,
)
from .snr_optimizer import (
    MAX_MTOTAL,
    MIN_CPT_MASS,
    Nfeval,
    callback_func,
    check_snr_optimizer_options,
    compute_minus_network_snr,
    compute_minus_network_snr_pso,
    compute_network_snr_core,
    insert_snr_optimizer_options,
    optimize_di,
    optimize_funcs,
    optimize_pso,
    optimize_shgo,
    option_dict,
    start_time,
)
from .supervision import (
    dict_to_args,
    ensure_directories,
    mail_volunteers_error,
    run_and_error,
    symlink,
    wait_for_utc_time,
)

__all__ = [
    "MAX_MTOTAL",
    "MIN_CPT_MASS",
    "Nfeval",
    "add_live_significance_duration_bin_options",
    "add_live_significance_trigger_pruning_options",
    "callback_func",
    "check_snr_optimizer_options",
    "compute_minus_network_snr",
    "compute_minus_network_snr_pso",
    "compute_network_snr_core",
    "dict_to_args",
    "duration_bins_from_cli",
    "ensure_directories",
    "insert_snr_optimizer_options",
    "mail_volunteers_error",
    "optimize_di",
    "optimize_funcs",
    "optimize_pso",
    "optimize_shgo",
    "option_dict",
    "run_and_error",
    "start_time",
    "symlink",
    "verify_live_significance_duration_bin_options",
    "verify_live_significance_trigger_pruning_options",
    "wait_for_utc_time",
]
