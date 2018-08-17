from pycbc.detector import Detector, get_available_detectors

# We can list the available detectors. This gives their detector abbreviation
# along with a longer name. Note that some of these are not physical detectors
# but may be useful for testing or study purposes

for abv, long_name in get_available_detectors():
    d = Detector(abv)

    # Note that units are all in radians
    print("{} {} Latitude {} Longitude {}".format(long_name, abv,
                                                  d.latitude,
                                                  d.longitude))
