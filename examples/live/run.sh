#!/bin/bash

# example/test of running PyCBC Live on simulated data

export OMP_NUM_THREADS=4
export HDF5_USE_FILE_LOCKING="FALSE"

# Get the current time in UTC format
current_time=$(date -u +"%Y-%m-%d %H:%M:%S")

# GPS time started on January 6, 1980 (Unix timestamp: 315964800)
gps_epoch_timestamp=315964800

# Convert the current time to seconds since the GPS epoch
gps_start_time=$(($(date -u -d "$current_time" +"%s") - $gps_epoch_timestamp))
gps_end_time=$((gps_start_time + 512))

echo "GPS Time: $gps_start_time"
echo "Current UTC Time: $current_time"
echo "GPS End Time: $gps_end_time"

f_min=17

CONF_DIR=/home/pycbc.live/analysis/prod/o4/full_bandwidth/
#--bank-file ${CONF_DIR}/bank/O4_DESIGN_OPT_FLOW_HYBRID_BANK_O3_CONFIG.hdf \

# make phase-time-amplitude histogram files, if needed

if [[ ! -f statHL.hdf ]]
then
    echo -e "\\n\\n>> [`date`] Making phase-time-amplitude files"

    bash ../search/stats.sh
else
    echo -e "\\n\\n>> [`date`] Pre-existing phase-time-amplitude files found"
fi

# delete old outputs if they exist
rm -rf ./output

echo -e "\\n\\n>> [`date`] Running PyCBC Live"

mpirun \
    -hostfile mpi_hosts.txt \
    -n 4 \
    -ppn 1 \
    \
python -m mpi4py `which pycbc_live` \
--bank-file bank_3k.hdf \
--sample-rate 2048 \
--enable-bank-start-frequency \
--low-frequency-cutoff ${f_min} \
--max-length 256 \
--approximant "SPAtmplt:mtotal<4" "SEOBNRv4_ROM:else" \
--chisq-bins "0.72*get_freq('fSEOBNRv4Peak',params.mass1,params.mass2,params.spin1z,params.spin2z)**0.7" \
--snr-abort-threshold 500 \
--snr-threshold 4.5 \
--newsnr-threshold 4.5 \
--max-triggers-in-batch 30 \
--store-loudest-index 50 \
--analysis-chunk 8 \
--autogating-threshold 50 \
--autogating-pad 0.5 \
--autogating-cluster 1 \
--autogating-width 0.25 \
--autogating-taper 0.25 \
--highpass-frequency 13 \
--highpass-bandwidth 5 \
--highpass-reduction 200 \
--psd-samples 30 \
--max-psd-abort-distance 600 \
--min-psd-abort-distance 68 \
--psd-abort-difference .15 \
--psd-recalculate-difference .01 \
--psd-inverse-length 3.5 \
--psd-segment-length 4 \
--trim-padding .5 \
--store-psd \
--increment-update-cache \
    H1:/dev/shm/kafka/H1_O3ReplayMDC \
    L1:/dev/shm/kafka/L1_O3ReplayMDC \
--frame-src \
    H1:/dev/shm/kafka/H1_O3ReplayMDC/* \
    L1:/dev/shm/kafka/L1_O3ReplayMDC/* \
--frame-read-timeout 50 \
--channel-name \
    H1:GDS-CALIB_STRAIN_INJ1_O3Replay \
    L1:GDS-CALIB_STRAIN_INJ1_O3Replay \
--state-channel \
    H1:GDS-CALIB_STATE_VECTOR \
    L1:GDS-CALIB_STATE_VECTOR \
--processing-scheme cpu:4 \
--fftw-measure-level 0 \
--fftw-input-float-wisdom-file ${CONF_DIR}/cit/fftw_wisdom \
--fftw-threads-backend openmp \
--increment 8 \
--max-batch-size 16777216 \
--output-path output \
--day-hour-output-prefix \
--sngl-ranking newsnr_sgveto_psdvar_threshold \
--ranking-statistic phasetd \
--statistic-files \
    statHL.hdf \
--sgchisq-snr-threshold 4 \
--sgchisq-locations "mtotal>40:20-30,20-45,20-60,20-75,20-90,20-105,20-120" \
--enable-background-estimation \
--background-ifar-limit 100 \
--timeslide-interval 0.1 \
--pvalue-combination-livetime 0.0005 \
--ifar-double-followup-threshold 0.0001 \
--ifar-upload-threshold 0.0002 \
--round-start-time 4 \
--start-time $gps_start_time \
--end-time $gps_end_time \
--src-class-mchirp-to-delta 0.01 \
--src-class-eff-to-lum-distance 0.74899 \
--src-class-lum-distance-to-delta -0.51557 -0.32195 \
--enable-profiling 1 \
--psd-variation \
--verbose

#--sngl-ranking newsnr_sgveto_psdvar_threshold \
#--ranking-statistic phasetd_exp_fit_fgbg_bbh_norm \
#    statHV.hdf \
#    statLV.hdf \
#    statHLV.hdf \
#    V1:Hrec_hoft_16384Hz_INJ1_O3Replay \
#    V1:DQ_ANALYSIS_STATE_VECTOR \
