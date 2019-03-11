#!/bin/bash
# Test pycbc inspiral by running over GW150914 with a limited template bank
echo -e "\\n\\n>> [`date`] Getting template bank"
wget -nc https://github.com/gwastro/pycbc-config/raw/master/test/inspiral/H1L1-SBANK_FOR_GW150914ER10.xml.gz

echo -e "\\n\\n>> [`date`] Running pycbc inspiral"
pycbc_inspiral \
--frame-type LOSC \
--sample-rate 2048 \
--sgchisq-snr-threshold 6.0 \
--sgchisq-locations "mtotal>40:20-30,20-45,20-60,20-75,20-90,20-105,20-120" \
--segment-end-pad 16 \
--cluster-method window \
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
--output H1-INSPIRAL-OUT.hdf \
--verbose \
--approximant "SPAtmplt:mtotal<4" "SEOBNRv4_ROM:mtotal<20" "SEOBNRv2_ROM_DoubleSpin:else" \
--bank-file H1L1-SBANK_FOR_GW150914ER10.xml.gz 2> inspiral.log

cat inspiral.log | head -10
cat inspiral.log | tail -20

echo -e "\\n\\n>> [`date`] test for GW150914"
python ./check_GW150914_detection.py H1-INSPIRAL-OUT.hdf
