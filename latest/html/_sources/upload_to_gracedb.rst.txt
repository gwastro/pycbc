################################################
Uploading triggers to gracedb
################################################

============
Introduction
============

This page describes how to upload your results to gracedb via the code pycbc_upload_xml_to_gracedb.

==========================
Requirements
==========================

To be able to run this code in production you should use the pre-packaged binaries. These will come with all necessary dependancies.

To run without the pre-packaged binaries you need:

* A version of PyCBC installed
* A version of ligo.gracedb client installed (pip install ligo-gracedb)

You also need to be authorized on the gracedb end to do the upload. The initial list of authorized people was given from the list of people on the front page of the PyCBC documentation under "PyCBC Contributors". If you are not on this list, or only recently added, please contact Branson Stephens in advance of box opening, to ensure that you are authorized.

============================
Finding input files
============================

To run this code requires some input files from your run. These can be found as follows where RUN_DIR is the place where you ran pycbc_submit_dax::

   TRIG_FILE=${RUN_DIR}/results/7._result/H1L1-PAGE_FOREGROUND_XMLLOUDEST-1130754617-1019335.xml

   PSD_FILES=${RUN_DIR}/psds/*MERGE_PSDS*hdf

============================
Running the command
============================

To run the upload command::

   pycbc_upload_xml_to_gracedb --input-file ${TRIG_FILE} --psd-files ${PSD_FILES}

Adding the --testing flag will upload events to the TEST group of gracedb, and can be used for testing purposes.

**WARNING**, if running with config files from before November 15 2015, the TRIG_FILE will have *many* triggers in it. Uploading hundreds of triggers to gracedb will be slow and a bad idea so the process has been to edit the TRIG_FILE by hand to restrict to the first 20 events in the coinc_inspiral and coinc_event tables and then uploading.
