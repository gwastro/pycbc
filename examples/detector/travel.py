from pycbc.detector import Detector

for ifo1 in ['H1', 'L1', 'V1']:
    for ifo2 in ['H1', 'L1', 'V1']:
        dt = Detector(ifo1).light_travel_time_to_detector(Detector(ifo2))
        print("Direct Time from {} to {} is {} seconds".format(ifo1, ifo2, dt))
