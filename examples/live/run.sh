#!/bin/bash

# example/test of running PyCBC Live on simulated data

set -e

export OMP_NUM_THREADS=4
export HDF5_USE_FILE_LOCKING="FALSE"

gps_start_time=1272790000
gps_end_time=1272790512
f_min=18


# test if there is a template bank. If not, make one

if [[ ! -f template_bank.hdf ]]
then
    echo -e "\\n\\n>> [`date`] Making template bank"
    curl \
        --remote-name \
        --silent \
        --show-error \
        https://raw.githubusercontent.com/gwastro/pycbc-config/710dbfd3590bd93d7679d7822da59fcb6b6fac0f/O2/bank/H1L1-HYPERBANK_SEOBNRv4v2_VARFLOW_THORNE-1163174417-604800.xml.gz

    pycbc_coinc_bank2hdf \
        --bank-file H1L1-HYPERBANK_SEOBNRv4v2_VARFLOW_THORNE-1163174417-604800.xml.gz \
        --output-file template_bank_full.hdf

    rm -f H1L1-HYPERBANK_SEOBNRv4v2_VARFLOW_THORNE-1163174417-604800.xml.gz

    pycbc_hdf5_splitbank \
        --bank-file template_bank_full.hdf \
        --output-prefix template_bank_ \
        --random-sort \
        --random-seed 831486 \
        --templates-per-bank 50

    mv template_bank_0.hdf template_bank.hdf
    rm -f template_bank_*.hdf
else
    echo -e "\\n\\n>> [`date`] Pre-existing template bank found"
fi


# test if there is a single fits file. If not, make a representative one
if [[ ! -f single_trigger_fits.hdf ]]
then
    echo -e "\\n\\n>> [`date`] Making single fits file"
    python make_singles_fits_file.py
else
    echo -e "\\n\\n>> [`date`] Pre-existing single fits file found"
fi


# test if there is a injection file.
# If not, make one and delete any existing strain

if [[ -f injections.hdf ]]
then
    echo -e "\\n\\n>> [`date`] Pre-existing injections found"
else
    echo -e "\\n\\n>> [`date`] Generating injections"

    rm -rf ./strain

    ./generate_injections.py
fi


# test if strain files exist. If they don't, make them

if [[ ! -d ./strain ]]
then
    echo -e "\\n\\n>> [`date`] Generating simulated strain"

    function simulate_strain { # detector PSD_model random_seed
        mkdir -p strain/$1

        out_path="strain/$1/$1-SIMULATED_STRAIN-{start}-{duration}.gwf"

        pycbc_condition_strain \
            --fake-strain $2 \
            --fake-strain-seed $3 \
            --output-strain-file $out_path \
            --gps-start-time $gps_start_time \
            --gps-end-time $4 \
            --sample-rate 16384 \
            --low-frequency-cutoff 10 \
            --channel-name $1:SIMULATED_STRAIN \
            --frame-duration 32 \
            --injection-file injections.hdf
    }
    # L1 ends 32s later, so that we can inject in single-detector time
    simulate_strain H1 aLIGOMidLowSensitivityP1200087 1234 $((gps_end_time - 32))
    simulate_strain L1 aLIGOMidLowSensitivityP1200087 2345 $gps_end_time
    simulate_strain V1 AdVEarlyLowSensitivityP1200087 3456 $((gps_end_time - 32))

else
    echo -e "\\n\\n>> [`date`] Pre-existing strain data found"
fi


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
-host localhost,localhost \
-n 2 \
--bind-to none \
 -x PYTHONPATH -x LD_LIBRARY_PATH -x OMP_NUM_THREADS -x VIRTUAL_ENV -x PATH -x HDF5_USE_FILE_LOCKING \
\
python -m mpi4py `which pycbc_live` \
--bank-file template_bank.hdf \
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
--min-psd-abort-distance 20 \
--psd-abort-difference .15 \
--psd-recalculate-difference .01 \
--psd-inverse-length 3.5 \
--psd-segment-length 4 \
--trim-padding .5 \
--store-psd \
--increment-update-cache \
    H1:strain/H1 \
    L1:strain/L1 \
    V1:strain/V1 \
--frame-src \
    H1:strain/H1/* \
    L1:strain/L1/* \
    V1:strain/V1/* \
--frame-read-timeout 10 \
--channel-name \
    H1:SIMULATED_STRAIN \
    L1:SIMULATED_STRAIN \
    V1:SIMULATED_STRAIN \
--processing-scheme cpu:4 \
--fftw-measure-level 0 \
--fftw-threads-backend openmp \
--increment 8 \
--max-batch-size 16777216 \
--output-path output \
--day-hour-output-prefix \
--sngl-ranking newsnr_sgveto \
--ranking-statistic phasetd \
--statistic-files statHL.hdf statHV.hdf statLV.hdf \
--sgchisq-snr-threshold 4 \
--sgchisq-locations "mtotal>40:20-30,20-45,20-60,20-75,20-90,20-105,20-120" \
--enable-background-estimation \
--background-ifar-limit 100 \
--timeslide-interval 0.1 \
--pvalue-combination-livetime 0.0005 \
--ifar-double-followup-threshold 0.0001 \
--ifar-upload-threshold 0.0001 \
--round-start-time 4 \
--start-time $gps_start_time \
--end-time $gps_end_time \
--src-class-mchirp-to-delta 0.01 \
--src-class-eff-to-lum-distance 0.74899 \
--src-class-lum-distance-to-delta -0.51557 -0.32195 \
--run-snr-optimization \
--enable-single-detector-background \
--single-newsnr-threshold 9 \
--single-duration-threshold 7 \
--single-reduced-chisq-threshold 2 \
--single-fit-file single_trigger_fits.hdf \
--sngl-ifar-est-dist conservative \
--verbose

# note that, at this point, some SNR optimization processes may still be
# running, so the checks below may ignore their results

# cat the logs of pycbc_optimize_snr so we can check them
for opt_snr_log in `find output -type f -name optimize_snr.log | sort`
do
    echo -e "\\n\\n>> [`date`] Showing log of SNR optimizer, ${opt_snr_log}"
    cat ${opt_snr_log}
done

echo -e "\\n\\n>> [`date`] Checking results"
./check_results.py \
    --gps-start ${gps_start_time} \
    --gps-end ${gps_end_time} \
    --f-min ${f_min} \
    --bank template_bank.hdf \
    --injections injections.hdf \
    --detectors H1 L1 V1

echo -e "\\n\\n>> [`date`] Running Bayestar"
for XMLFIL in `find output -type f -name \*.xml\* | sort`
do
    pushd `dirname ${XMLFIL}`
    bayestar-localize-coincs --f-low ${f_min} `basename ${XMLFIL}` `basename ${XMLFIL}`
    test -f 0.fits
    popd
done
