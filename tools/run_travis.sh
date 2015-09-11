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

python test/test_correlate.py
test $? -ne 0 && RESULT=1

#python test/test_fft_unthreaded.py
#test $? -ne 0 && RESULT=1

#python test/test_frame.py
#test $? -ne 0 && RESULT=1

python test/test_frequencyseries.py
test $? -ne 0 && RESULT=1

#python test/test_injection.py
#test $? -ne 0 && RESULT=1

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

# check for trivial failures of important executables

function test_exec_help {
    $1 --help > /dev/null
    test $? -ne 0 && RESULT=1
}

test_exec_help pycbc_inspiral
test_exec_help pycbc_geom_nonspinbank
test_exec_help pycbc_aligned_stoch_bank
test_exec_help pycbc_geom_aligned_bank
test_exec_help pycbc_splitbank
test_exec_help pycbc_strip_injections
test_exec_help pycbc_coinc_bank2hdf
test_exec_help pycbc_coinc_mergetrigs
test_exec_help pycbc_coinc_findtrigs
test_exec_help pycbc_coinc_hdfinjfind
test_exec_help pycbc_coinc_statmap
test_exec_help pycbc_coinc_statmap_inj
test_exec_help pycbc_foreground_censor
test_exec_help pycbc_plot_singles_vs_params
test_exec_help pycbc_page_snrchi
test_exec_help pycbc_page_coinc_snrchi
test_exec_help pycbc_page_sensitivity
test_exec_help pycbc_page_foreground
test_exec_help pycbc_page_foundmissed
test_exec_help pycbc_page_ifar
test_exec_help pycbc_page_injtable
test_exec_help pycbc_page_segments
#test_exec_help pycbc_page_segplot
#test_exec_help pycbc_page_segtable
test_exec_help pycbc_page_snrifar
#test_exec_help pycbc_page_vetotable
test_exec_help pycbc_plot_bank_bins
test_exec_help pycbc_plot_hist
test_exec_help pycbc_calculate_psd
test_exec_help pycbc_plot_psd_file
test_exec_help pycbc_plot_range
test_exec_help pycbc_average_psd
#test_exec_help pycbc_make_html_page

exit ${RESULT}
