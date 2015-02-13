#! /bin/bash

rm -rf /home/cbiwer/public_html/pycbc_plotting/970012743-970099143/*

./pycbc_make_html_page --analysis-title "\"PyCBC Coincidence Analysis\"" --analysis-subtitle "\"Christopher M. Biwer\"" --template-dir "./templates/" --template-name "base.html" --output-path "/home/cbiwer/public_html/pycbc_plotting/970012743-970099143/" --plots-dir "./rundir/970012743-970099143/results_plots/" 

