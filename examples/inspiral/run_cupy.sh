
# Select bank file: use argument or default to small bank
BANK_FILE=${1:-BANK_SPLIT0.hdf}

# Uncomment if you want profile information

#nvprof --print-callstack --log-file profout.txt \

#/home/ian.harry/.conda/envs/env_lisa_premerger/bin/../nsight-compute/2024.1.1/ncu -o profile \

#python -m cProfile -o output_cupy.pstats \

#/home/ian.harry/nsight-systems-2024.7.1/bin/nsys profile \
#  --trace cuda,osrt,nvtx \
#  --cuda-memory-usage true \
#  --force-overwrite true \
#  --output profile_run_v1 \
#  --python-sampling=true \
#`which python` -m nvtx -- ~/.conda/envs/env_lisa_premerger/bin/pycbc_inspiral \
#python -m cProfile -o output_cupy.pstats `which pycbc_inspiral` \
#time ~/.conda/envs/env_lisa_premerger/bin/pycbc_inspiral \

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
--output H1-INSPIRAL_GPU-OUT.hdf \
--approximant "SPAtmplt" \
--processing-scheme cupy \
--bank-file ${BANK_FILE} \
--gpu-batch-size 128 \
--chisq-snr-threshold 5.25 \
--verbose
#--verbose 2> inspiral_$1_$2_$3.log
# Uncomment above if you want logging

# Uncomment for profile pngs
#gprof2dot -f pstats output_$1_$2_$3.pstats | dot -Tpng -o $1_$2_$3.png
