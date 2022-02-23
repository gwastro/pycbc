import os.path
# Setup the directory with the NS equilibrium sequence(s)
NS_SEQUENCE_FILE_DIRECTORY = os.path.join(
    os.path.dirname(__file__), 'ns_sequences')
NS_SEQUENCES = [
    f.replace('equil_', '').replace('.dat', '')
    for f in os.listdir(NS_SEQUENCE_FILE_DIRECTORY) if f.endswith('.dat')]
from pycbc.neutron_stars.ns_functions import *
