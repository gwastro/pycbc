import fstat
import math
import numpy
from pycbc.waveform import get_td_waveform
import pycbc.detector

def zeros_generator():
	##############
	# WAVEFORM GENERATION
	#####################
	# generate the waveform
	f_low = 30
	sample_rate = 16384
	dist = 3e7
	inc = math.pi/6
	pol = 0.5 # math.pi/4
	phase = 0.5

	original_a = fstat.params_to_a(dist, numpy.cos(inc), pol, phase)
	print("Starting F-stat A values")
	print original_a

	#Generate a waveform
	hp, hc = get_td_waveform(approximant="TaylorT4",
			mass1=1.4,
			mass2=1.4,
			f_lower=f_low,
			delta_t=1.0/sample_rate,
			distance=dist,
			inclination=inc,
			coa_phase=phase)
	#ceil(x) is the smallest integer larger than x
	newlen = int(2**(math.ceil(math.log(len(hp), 2))))
	hp.resize(newlen)
	hc.resize(newlen)

	#########################
	# Find strain
	#########################

	ifoList = ['H1','L1']

	#Make up sky location and time
	dec = 0.5
	ra = 0.5
	t_gps = 1000

	#Get the antenna functions
	ifo_dictionary = {}
	for key in ifoList:
	    ifo_dictionary[key] = pycbc.detector.Detector(key)

	Fp_dictionary = {}
	Fc_dictionary = {}
	strain = {}
	for key in ifoList:
	    Fp_dictionary[key], Fc_dictionary[key] = ifo_dictionary[key].antenna_pattern(ra, dec, pol, t_gps)
	    strain[key] = hp * Fp_dictionary[key] + hc * Fc_dictionary[key]
	print("Initial parameters: D = %.2f Mpc; i = %.2f; psi = %.2f; phi = %.2f"
	    % (dist, inc, pol, phase))
        return strain
