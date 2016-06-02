#! /bin/bash

#
# Injection parameters
#
INJ_OUTPUT=test_injection.xml.gz
INJ_APPROX=SEOBNRv2
TRIGGER_TIME=1126259462.42
# masses are observed, not source
MASS1=37.
MASS2=32.
RA=2.21535724066
DEC=-1.23649695537
INC=0
DISTANCE=400000 # in kpc
INJ_F_MIN=28.

#
# Parameters for MCMC
#
# the amount of time simulated; the start,
# stop times of the segment will be the trigger
# time +/-SEGLEN, and the df used will = 1/(2*SEGLEN)
OUTPUT=ringdown_example.hdf
SEGLEN=2
IFOS="H1 L1"
STRAIN="H1:aLIGOZeroDetHighPower L1:aLIGOZeroDetHighPower"
SAMPLE_RATE=1024
F_MIN=30.
N_WALKERS=500
N_ITERATIONS=1000
PROCESSING_SCHEME=cpu
CONFIG_PATH=ringdown_example.ini


#
# Do not edit below this
#
TRIGGER_TIME_INT=${TRIGGER_TIME%.*}

GPS_START_TIME=$((${TRIGGER_TIME_INT} - ${SEGLEN}))
GPS_END_TIME=$((${TRIGGER_TIME_INT} + ${SEGLEN}))

echo "Creating injection file"
lalapps_inspinj \
    --output ${INJ_OUTPUT} \
    --seed 1000 \
    --f-lower ${INJ_F_MIN} \
    --waveform ${INJ_APPROX} \
    --amp-order 7 \
    --gps-start-time ${TRIGGER_TIME} \
    --gps-end-time ${TRIGGER_TIME} \
    --time-step 1 \
    --t-distr fixed \
    --l-distr fixed \
    --longitude ${RA} \
    --latitude ${DEC} \
    --d-distr uniform \
    --min-distance ${DISTANCE} \
    --max-distance ${DISTANCE} \
    --i-distr fixed \
    --fixed-inc ${INC} \
    --m-distr fixMasses \
    --fixed-mass1 ${MASS1} \
    --fixed-mass2 ${MASS2} \
    --disable-spin

echo "Running MCMC"
echo "============================"
pycbc_mcmc --verbose \
    --instruments ${IFOS} \
    --gps-start-time ${GPS_START_TIME} \
    --gps-end-time ${GPS_END_TIME} \
    --psd-model ${STRAIN} \
    --fake-strain ${STRAIN} \
    --sample-rate ${SAMPLE_RATE} \
    --low-frequency-cutoff ${F_MIN} \
    --processing-scheme ${PROCESSING_SCHEME} \
    --sampler kombine \
    --likelihood-evaluator gaussian \
    --nwalkers ${N_WALKERS} \
    --niterations ${N_ITERATIONS} \
    --config-file ${CONFIG_PATH} \
    --output-file ${OUTPUT}
