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
for prog in `find ${PATH//:/ } -maxdepth 1 -name 'pycbc*' -print 2>/dev/null | egrep -v '(pycbc_live_nagios_monitor|pycbc_make_grb_summary_page|pycbc_make_offline_grb_workflow|pycbc_mvsc_get_features|pycbc_upload_xml_to_gracedb|pycbc_make_skymap)'`
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

# run PyCBC Live test if running in Python 3
if [ "x$PYTHON3" = "xTRUE" ]
then
    pushd examples/live
    bash -e run.sh
    if test $? -ne 0 ; then
        RESULT=1
        echo -e "    FAILED!"
        echo -e "---------------------------------------------------------"
    else
        echo -e "    Pass."
    fi
    popd
fi

# Run Inference Scripts
## Run inference on 2D-normal analytic likelihood function
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

## Run inference on BBH example; this will also run
## a test of create_injections
pushd examples/inference/bbh-injection
bash -e make_injection.sh
if test $? -ne 0 ; then
    RESULT=1
    echo -e "    FAILED!"
    echo -e "---------------------------------------------------------"
else
    echo -e "    Pass."
fi
# now run inference
bash -e run_test.sh
if test $? -ne 0 ; then
    RESULT=1
    echo -e "    FAILED!"
    echo -e "---------------------------------------------------------"
else
    echo -e "    Pass."
fi
popd

## Run inference on GW150914 data
pushd examples/inference/gw150914
bash -e run_test.sh
if test $? -ne 0 ; then
    RESULT=1
    echo -e "    FAILED!"
    echo -e "---------------------------------------------------------"
else
    echo -e "    Pass."
fi
popd

## Run inference using single template model
pushd examples/inference/single
bash -e get.sh
bash -e run.sh
if test $? -ne 0 ; then
    RESULT=1
    echo -e "    FAILED!"
    echo -e "---------------------------------------------------------"
else
    echo -e "    Pass."
fi
popd

## Run inference using relative model
pushd examples/inference/relative
bash -e get.sh
bash -e run.sh
if test $? -ne 0 ; then
    RESULT=1
    echo -e "    FAILED!"
    echo -e "---------------------------------------------------------"
else
    echo -e "    Pass."
fi
popd

## Run inference samplers
pushd examples/inference/samplers
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

python setup.py build_gh_pages
if test $? -ne 0 ; then
    echo -e "    FAILED!"
    echo -e "---------------------------------------------------------"
    RESULT=1
fi

exit ${RESULT}
