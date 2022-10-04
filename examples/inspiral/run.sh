#!/bin/bash -e

function inspiral_run {
echo -e "\\n\\n>> [`date`] Running pycbc inspiral $1:$3 with $2 threads"
# Uncomment if you want profile information
#python -m cProfile -o output_$1_$2_$3.pstats
`which pycbc_inspiral` \
--frame-files DATA_FILE.gwf \
--sample-rate 2048 \
--sgchisq-snr-threshold 6.0 \
--sgchisq-locations "mtotal>40:20-30,20-45,20-60,20-75,20-90,20-105,20-120" \
--segment-end-pad 16 \
--low-frequency-cutoff 30 \
--pad-data 8 \
--cluster-window 1 \
--cluster-function symmetric \
--injection-window 4.5 \
--segment-start-pad 112 \
--psd-segment-stride 8 \
--psd-inverse-length 16 \
--filter-inj-only \
--psd-segment-length 16 \
--snr-threshold 5.5 \
--segment-length 256 \
--autogating-width 0.25 \
--autogating-threshold 100 \
--autogating-cluster 0.5 \
--autogating-taper 0.25 \
--newsnr-threshold 5 \
--psd-estimation median \
--strain-high-pass 20 \
--order -1 \
--chisq-bins "1.75*(get_freq('fSEOBNRv2Peak',params.mass1,params.mass2,params.spin1z,params.spin2z)-60.)**0.5" \
--channel-name H1:LOSC-STRAIN \
--gps-start-time 1126259078 \
--gps-end-time 1126259846 \
--output H1-INSPIRAL_$1_$2_$3-OUT.hdf \
--approximant "SPAtmplt:mtotal<4" "SEOBNRv4_ROM:mtotal<20" "SEOBNRv2_ROM_DoubleSpin:else" \
--fft-backends $1 \
--processing-scheme cpu:$2 \
--fftw-threads-backend $3 \
--use-compressed-waveforms \
--bank-file COMPRESSED_BANK.hdf \
#--verbose 2> inspiral_$1_$2_$3.log
# Uncomment above if you want logging

# Uncomment for profile pngs
#gprof2dot -f pstats output_$1_$2_$3.pstats | dot -Tpng -o $1_$2_$3.png

echo -e "\\n\\n>> [`date`] test for GW150914"
python ./check_GW150914_detection.py H1-INSPIRAL_$1_$2_$3-OUT.hdf
}


# Test pycbc inspiral by running over GW150914 with a limited template bank
echo -e "\\n\\n>> [`date`] Getting template bank"
wget -nv -nc https://github.com/gwastro/pycbc-config/raw/master/test/inspiral/SMALLER_BANK_FOR_GW150914.hdf

echo -e "\\n\\n>> [`date`] Compressing template bank"
pycbc_compress_bank --bank-file SMALLER_BANK_FOR_GW150914.hdf --output COMPRESSED_BANK.hdf --sample-rate 4096 --segment-length 256 --compression-algorithm mchirp --psd-model aLIGOZeroDetHighPower --low-frequency-cutoff 30 --approximant "SEOBNRv4_ROM"

echo -e "\\n\\n>> [`date`] Creating data file"
pycbc_condition_strain \
--frame-type LOSC \
--sample-rate 2048 \
--pad-data 8 \
--autogating-width 0.25 \
--autogating-threshold 100 \
--autogating-cluster 0.5 \
--autogating-taper 0.25 \
--strain-high-pass 10 \
--channel-name H1:LOSC-STRAIN \
--gps-start-time 1126258578 \
--gps-end-time 1126259946 \
--output-strain-file DATA_FILE.gwf \


start=`date +%s`
inspiral_run fftw 16 openmp
end=`date +%s`
runtime_fftw_openmp_16=$((end-start))

start=`date +%s`
inspiral_run fftw 16 pthreads
end=`date +%s`
runtime_fftw_pthreads_16=$((end-start))

# In all cases below the last argument is irrelevant
start=`date +%s`
inspiral_run fftw 1 openmp
end=`date +%s`
runtime_fftw_1=$((end-start))

start=`date +%s`
inspiral_run mkl 16 openmp
end=`date +%s`
runtime_mkl_16=$((end-start))

start=`date +%s`
inspiral_run mkl 1 openmp
end=`date +%s`
runtime_mkl_1=$((end-start))

# Numpy doesn't include threading. We just want it to run, it won't be fast!
start=`date +%s`
inspiral_run numpy 1 openmp
end=`date +%s`
runtime_numpy_1=$((end-start))

echo "RUN TIMES"
echo "FFTW 16 threads with openmp " ${runtime_fftw_openmp_16}
echo "FFTW 16 threads with pthreads" ${runtime_fftw_pthreads_16}
echo "FFTW single core" ${runtime_fftw_1}
echo "MKL 16 threads" ${runtime_mkl_16}
echo "MKL 1 thread" ${runtime_mkl_1}
echo "Numpy 1 thread" ${runtime_numpy_1}
