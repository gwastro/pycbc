#!/bin/bash

echo -e "\\n>> [`date`] Starting PyCBC test suite"

PYTHON_VERSION=`python -c 'import sys; print(sys.version_info.major)'`
echo -e "\\n>> [`date`] Python Major Version:" $PYTHON_VERSION
PYTHON_MINOR_VERSION=`python -c 'import sys; print(sys.version_info.minor)'`
echo -e "\\n>> [`date`] Python Minor Version:" $PYTHON_MINOR_VERSION

# This will work from anywhere within the pycbc directory
this_script_dir=`dirname -- "$( readlink -f -- "$0"; )"`
cd $this_script_dir
cd ..

LOG_FILE=$(mktemp -t pycbc-test-log.XXXXXXXXXX)

RESULT=0
cat_output=true

function test_result {
    if test $? -ne 0 ; then
        RESULT=1
        echo -e "    FAILED!"
        if $cat_output ; then
            echo -e "---------------------------------------------------------"
            cat $LOG_FILE
            echo -e "---------------------------------------------------------"
        fi
    else
        echo -e "    Pass"
    fi
}

if [ "$PYCBC_TEST_TYPE" = "unittest" ] || [ -z ${PYCBC_TEST_TYPE+x} ]; then
    for prog in `find test -name '*.py' -print | egrep -v '(long|lalsim|test_waveform)'`
    do
        prog_short=`echo $prog | rev | cut -d"/" -f1 | rev`
        echo -e ">> [`date`] running unit test for $prog_short"
        python $prog &> $LOG_FILE
        test_result
    done
fi

if [ "$PYCBC_TEST_TYPE" = "help" ] || [ -z ${PYCBC_TEST_TYPE+x} ]; then
    # check that all executables that do not require
    # special environments can return a help message
    for prog in `find ${PATH//:/ } -maxdepth 1 -name 'pycbc*' -print 2>/dev/null | egrep -v '(pycbc_live_nagios_monitor|pycbc_mvsc_get_features)' | sort | uniq`
    do
        echo -e ">> [`date`] running $prog --help"
        $prog --help &> $LOG_FILE
        test_result
        if [[ `echo $prog | egrep '(pycbc_copy_output_map|pycbc_submit_dax|pycbc_stageout_failed_workflow)'` ]] ; then
            continue
        fi
        echo -e ">> [`date`] running $prog --version"
        $prog --version &> $LOG_FILE
        test_result
    done
    # also check that --version with increased modifiers works for one executable
    echo -e ">> [`date`] running pycbc_inspiral --version with modifiers"
    for modifier in "" 0 1 2 3
    do
        echo -e ">> [`date`] running pycbc_inspiral --version ${modifier}"
        pycbc_inspiral --version ${modifier} &> $LOG_FILE
        test_result
    done
fi

cat_output=false

if [ "$PYCBC_TEST_TYPE" = "search" ] || [ -z ${PYCBC_TEST_TYPE+x} ]; then
    # run pycbc inspiral test
    pushd examples/inspiral
    bash -e run.sh
    test_result
    popd

    # run a quick bank placement example
    pushd examples/tmpltbank
    bash -e testNonspin2.sh
    test_result
    popd

    # run PyCBC Live test
    pushd examples/live
    bash -e run.sh
    test_result
    popd

    # run pycbc_multi_inspiral (PyGRB) test
    pushd examples/multi_inspiral
    bash -e run.sh
    test_result
    popd
fi

if [ "$PYCBC_TEST_TYPE" = "inference" ] || [ -z ${PYCBC_TEST_TYPE+x} ]; then
    # Run Inference Scripts
    ## Run inference on 2D-normal analytic likelihood function
    pushd examples/inference/analytic-normal2d
    bash -e run.sh
    test_result
    popd

    ## Run inference on BBH example; this will also run
    ## a test of create_injections
    pushd examples/inference/bbh-injection
    bash -e make_injection.sh
    test_result
    # now run inference
    bash -e run_test.sh
    test_result
    popd

    ## Run inference on GW150914 data
    pushd examples/inference/gw150914
    bash -e run_test.sh
    test_result
    popd

    ## Run inference using single template model
    pushd examples/inference/single
    bash -e get.sh
    bash -e run.sh
    test_result
    popd

    ## Run inference using relative model
    pushd examples/inference/relative
    bash -e get.sh
    bash -e run.sh
    test_result
    popd

    ## Run inference using the hierarchical model
    pushd examples/inference/hierarchical
    bash -e run_test.sh
    test_result
    popd

    ## Run inference samplers
    pushd examples/inference/samplers
    bash -e run.sh
    test_result
    popd

    ## Run pycbc_make_skymap example
    pushd examples/make_skymap
    bash -e simulated_data.sh
    test_result
    popd
fi

if [ "$PYCBC_TEST_TYPE" = "docs" ] || [ -z ${PYCBC_TEST_TYPE+x} ]; then
    echo -e "\\n>> [`date`] Building documentation"

    python setup.py build_gh_pages
    test_result
fi

exit ${RESULT}
