#!/bin/bash

set -v

INST=${HOME}/inst
source ${INST}/etc/lal-user-env.sh
export LD_LIBRARY_PATH=${INST}/lib:${INST}/lib64
export PKG_CONFIG_PATH=${INST}/lib/pkgconfig:${PKG_CONFIG_PATH}
export PATH=/usr/lib/ccache:${PATH}:${INST}/bin

# Using python setup.py test has two issues:
#     Some tests fail for reasons not necessarily related to PyCBC
#     Setup.py seems to returns 0 even when tests fail
# So we rather run specific tests manually

RESULT=0

python test/test_array_lal.py
test $? -ne 0 && RESULT=1

python test/test_array.py
test $? -ne 0 && RESULT=1

#python test/test_autochisq.py
#test $? -ne 0 && RESULT=1

python test/test_chisq.py
test $? -ne 0 && RESULT=1

python test/test_conversions.py
test $? -ne 0 && RESULT=1

python test/test_correlate.py
test $? -ne 0 && RESULT=1

#python test/test_fft_unthreaded.py
#test $? -ne 0 && RESULT=1

python test/test_frame.py
test $? -ne 0 && RESULT=1

python test/test_frequencyseries.py
test $? -ne 0 && RESULT=1

python test/test_injection.py
test $? -ne 0 && RESULT=1

python test/test_matchedfilter.py
test $? -ne 0 && RESULT=1

python test/test_pnutils.py
test $? -ne 0 && RESULT=1

python test/test_psd.py
test $? -ne 0 && RESULT=1

python test/test_resample.py
test $? -ne 0 && RESULT=1

#python test/test_schemes.py
#test $? -ne 0 && RESULT=1

python test/test_threshold.py
test $? -ne 0 && RESULT=1

python test/test_timeseries.py
test $? -ne 0 && RESULT=1

python test/test_tmpltbank.py
test $? -ne 0 && RESULT=1

python test/test_spatmplt.py
test $? -ne 0 && RESULT=1

python test/test_inference.py
test $? -ne 0 && RESULT=1

# check that all executables that do not require
# special environments can return a help message
for prog in `find ${PATH//:/ } -maxdepth 1 -name 'pycbc*' -print | egrep -v '(pycbc_fit_sngl_trigs|pycbc_live|pycbc_live_nagios_monitor|pycbc_make_grb_summary_page|pycbc_make_offline_grb_workflow|pycbc_mvsc_get_features|pycbc_upload_xml_to_gracedb)'`
do
    echo "Checking $prog --help"
    $prog --help > /dev/null
    test $? -ne 0 && RESULT=1
done

exit ${RESULT}
