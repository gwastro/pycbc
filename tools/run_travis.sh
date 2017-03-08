#!/bin/bash

set -v

BUILD=${HOME}/build
BUILDDIRNAME="pycbc-build"
PYCBC="$BUILD/$BUILDDIRNAME"
PYTHON_PREFIX="$PYCBC"
ENVIRONMENT="$PYCBC/environment"
PREFIX="$ENVIRONMENT"
PATH="$PREFIX/bin:$PYTHON_PREFIX/bin:$PATH"
export LD_LIBRARY_PATH="$PREFIX/lib:$PREFIX/bin:$PYTHON_PREFIX/lib:/usr/local/lib:$LD_LIBRARY_PATH"
export PKG_CONFIG_PATH="$PREFIX/lib/pkgconfig:$PYTHON_PREFIX/lib/pkgconfig:/usr/local/lib/pkgconfig:$PKG_CONFIG_PATH"
source ${BUILD}/pycbc-build/environment/etc/lalsuite-user-env.sh
source ${BUILD}/pycbc-build/environment/bin/activate

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
