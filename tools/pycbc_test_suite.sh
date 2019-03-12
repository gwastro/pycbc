#!/bin/bash

echo -e "\\n>> [`date`] Starting PyCBC test suite"

LOG_FILE=$(mktemp -t pycbc-test-log.XXXXXXXXXX)

RESULT=0

# Using python setup.py test has two issues:
#     Some tests fail for reasons not necessarily related to PyCBC
#     Setup.py seems to returns 0 even when tests fail
# So we rather run specific tests manually
#
# test_psd disabled because analytic files not included in lalsuite wheel

for prog in `find test -name '*.py' -print | egrep -v '(long|lalsim|test_waveform|test_psd)'`
do
    echo -e ">> [`date`] running unit test for $prog"
    python $prog &> $LOG_FILE
    if test $? -ne 0 ; then
        RESULT=1
        echo -e "    FAILED!"
        echo -e "---------------------------------------------------------"
        cat $LOG_FILE
        echo -e "---------------------------------------------------------"
    else
        echo -e "    Pass."
    fi
done

# check that all executables that do not require
# special environments can return a help message
for prog in `find ${PATH//:/ } -maxdepth 1 -name 'pycbc*' -print 2>/dev/null | egrep -v '(pycbc_live_nagios_monitor|pycbc_make_grb_summary_page|pycbc_make_offline_grb_workflow|pycbc_mvsc_get_features|pycbc_upload_xml_to_gracedb)'`
do
    echo -e ">> [`date`] running $prog --help"
    $prog --help &> $LOG_FILE
    if test $? -ne 0 ; then
        RESULT=1
        echo -e "    FAILED!"
        echo -e "---------------------------------------------------------"
        cat $LOG_FILE
        echo -e "---------------------------------------------------------"
    else
        echo -e "    Pass."
    fi
done

#run pycbc inspiral test
pushd examples/inspiral
bash -e run.sh
if test $? -ne 0 ; then
    RESULT=1
    echo -e "    FAILED!"
    echo -e "---------------------------------------------------------"
else
    echo -e "    Pass."
fi
popd

# Run Inference Scripts
pushd examples/inference/analytic-normal2d
bash -e run.sh
if test $? -ne 0 ; then
    RESULT=1
    echo -e "    FAILED!"
    echo -e "---------------------------------------------------------"
else
    echo -e "    Pass."
fi
popd

echo -e "\\n>> [`date`] Building documentation"

python setup.py build_gh_pages &> $LOG_FILE
if test $? -ne 0 ; then
    echo -e "    FAILED!"
    echo -e "---------------------------------------------------------"
    cat $LOG_FILE
    echo -e "---------------------------------------------------------"
    RESULT=1
fi

exit ${RESULT}
