#!/bin/bash

echo -e "\\n\\n>> [`date`] Starting PyCBC test suite"

LOG_FILE=$(mktemp -t pycbc-test-log.XXXXXXXXXX)
echo -e "\\n\\n>> [`date`] writing test log to $LOG_FILE"

function exit_on_error {
    echo "--- Error or interrupt ----------------------------------------" >&4
    if [ -f $LOG_FILE ] ; then
        cat $LOG_FILE >&4
    fi
    echo "---------------------------------------------------------------" >&4
exit 1
}
trap exit_on_error ERR INT

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

# make a copy of stdin and stdout and close them
exec 3>&1-
exec 4>&2-

# open stdout as $LOG_FILE file for read and write.
exec 1<>$LOG_FILE

# redirect stderr to stdout
exec 2>&1

set -v

RESULT=0

# Using python setup.py test has two issues:
#     Some tests fail for reasons not necessarily related to PyCBC
#     Setup.py seems to returns 0 even when tests fail
# So we rather run specific tests manually
for prog in `find test -name '*.py' -print | egrep -v '(autochisq|bankveto|fft|schemes|long)'`
do 
    echo -e "\\n\\n>> [`date`] running unit test for $prog" >&3
    python $prog
    test $? -ne 0 && RESULT=1
done

# check that all executables that do not require
# special environments can return a help message
for prog in `find ${PATH//:/ } -maxdepth 1 -name 'pycbc*' -print | egrep -v '(pycbc_fit_sngl_trigs|pycbc_live|pycbc_live_nagios_monitor|pycbc_make_grb_summary_page|pycbc_make_offline_grb_workflow|pycbc_mvsc_get_features|pycbc_upload_xml_to_gracedb)'`
do
    echo -e "\\n\\n>> [`date`] running $prog --help" >&3
    $prog --help
    test $? -ne 0 && RESULT=1
done

exit ${RESULT}
