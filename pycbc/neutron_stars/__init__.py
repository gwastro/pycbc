import os.path

# Setup the directory with the NS equilibrium sequence(s)
NS_DATA_DIRECTORY = os.path.join(os.path.dirname(__file__), "ns_data")
NS_SEQUENCES = [
    f.replace("equil_", "").replace(".dat", "")
    for f in os.listdir(NS_DATA_DIRECTORY)
    if f.endswith(".dat")
]
# These import NS_DATA_DIRECTORY and NS_SEQUENCES from this package, so
# they must come after them
from pycbc.neutron_stars.eos_utils import (  # noqa: E402
    foucart18,
    initialize_eos,
    interp_grav_mass_to_baryon_mass,
    interp_grav_mass_to_compactness,
    load_ns_sequence,
)
from pycbc.neutron_stars.pg_isso_solver import (  # noqa: E402
    isco_solution,
    isso_eq_at_pole,
    isso_eq_at_pole_dr,
    isso_eq_at_pole_dr2,
    pg_isso_eq,
    pg_isso_eq_dr,
    pg_isso_eq_dr2,
    solve_pg_isso,
)

__all__ = [
    "NS_DATA_DIRECTORY",
    "NS_SEQUENCES",
    "foucart18",
    "initialize_eos",
    "interp_grav_mass_to_baryon_mass",
    "interp_grav_mass_to_compactness",
    "isco_solution",
    "isso_eq_at_pole",
    "isso_eq_at_pole_dr",
    "isso_eq_at_pole_dr2",
    "load_ns_sequence",
    "pg_isso_eq",
    "pg_isso_eq_dr",
    "pg_isso_eq_dr2",
    "solve_pg_isso",
]
