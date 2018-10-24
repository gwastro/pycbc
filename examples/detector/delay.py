from pycbc.detector import Detector
from astropy.utils import iers

# Make sure the documentation can be built without an internet connection
iers.conf.auto_download = False

# The source of the gravitational waves
right_ascension = 0.7
declination = -0.5

# Reference location will be the Hanford detector
# see the `time_delay_from_earth_center` method to use use geocentric time
# as the reference
dref = Detector("H1")

# Time in GPS seconds that the GW passes
time = 100000000

# Time that the GW will (or has) passed through the given detector
for ifo in ["H1", "L1", "V1"]:
    d = Detector(ifo)
    dt = d.time_delay_from_detector(dref, right_ascension, declination, time)
    st = "GW passed through {} {} seconds relative to passing by Hanford"
    print(st.format(ifo, dt))
