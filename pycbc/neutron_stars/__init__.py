import os.path
# Setup the directory with the NS equilibrium sequence(s)
NS_DATA_DIRECTORY = os.path.join(
    os.path.dirname(__file__), 'ns_data')
NS_SEQUENCES = [
    f.replace('equil_', '').replace('.dat', '')
    for f in os.listdir(NS_DATA_DIRECTORY) if f.endswith('.dat')]
from pycbc.neutron_stars.eos_utils import *
from pycbc.neutron_stars.pg_isso_solver import *
