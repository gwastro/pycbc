###################################################################
PyCBC calibrate data documentation (``pycbc_calibrate_data``)
###################################################################

===================
Introduction
===================

This page gives details on how to use the `pycbc_calibrate_data` executable
available in the pyCBC examples subdirectory. The purpose of this example is to
calibrate h(t) data the same as in the front end of the interferometer.

The executable `pycbc_calibrate_data` will read in the `DARM_CTRL` and `DARM_ERR`
channels, and apply the filterbanks described in a foton filter file to the data.
This example assumes the `OAF` model is being used.

=================================================
Usage Example
=================================================

---------------------
Get data
---------------------

This must be done at the LLO or LHO clusters at the moment, since they have the foton and ROOT modules correctly installed. First you need to get the location of the frame files. An example looks like::

   cd examples/cal/calibrate_data/
   gw_data_find --observatory L --type L1_R --gps-start-time 1102660100 --gps-end-time 1102660200 --url-type file --lal-cache > L1-DATAFIND_FRAMES-1102660100-100.lcf

---------------------
Calibrate the data
---------------------

Then you will calibrate the data. In order to calibrate the data you will need to give it a foton filter file, one is included in this example subdirectory at `./filter_files/L1OAF_8379.txt`, and the frame cache that was just wtitten by `gw_data_find`. An example looks like::

   ./pycbc_calibrate_data --gps-start-time 1102660100 --gps-end-time 1102660200 --frame-cache ./L1-DATAFIND_FRAMES-1102660100-100.lcf --filter-file ./filter_files/L1OAF_8379.txt --output-file L1-RECALIBRATE_DATA-1102660100-100.txt

---------------------
Diagnostic plots
---------------------

There are two executables in this example subdirectory to make diagnostic plots. You can plot the timeseries in the frame files against the calibrated data. An example looks like::

    mkdir -p ~/public_html/tmp/cal_dev/
    ./pycbc_plot_txt L1-RECALIBRATE_DATA-1102660100-100.txt ~/public_html/tmp/cal_dev/calibrate_oaf_cal_darm_dq_5.png
    ./pycbc_plot_frame L1:OAF-CAL_DARM_DQ ./L1-DATAFIND_FRAMES-1102660100-100.lcf 1102660195 1102660200 ~/public_html/tmp/cal_dev/cal_dev_frame.png

The arguments for `pycbc_plot_txt` are the `pycbc_calibrate_data` data file and the output location for the plot. The arguments for `pycbc_plot_frame` are the channel name, frame cache, start time, end time, and output location for the plot.
