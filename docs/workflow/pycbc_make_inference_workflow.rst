############################################################################
``pycbc_make_inference_workflow``: A parameter estimation workflow generator
############################################################################

===============
Introduction
===============

The executable ``pycbc_make_inference_workflow`` is a tool used to analyse data

===========================
Workflow configuration file
===========================

If you want to analyse a different set of ifos, you will have to edit some additional options::

    [workflow]
    end-time = 1137254417
    start-time = 1133173817
    file-retention-level = all_triggers
    h1-channel-name = H1:DCS-CALIB_STRAIN_C01
    l1-channel-name = L1:DCS-CALIB_STRAIN_C01

    [workflow-ifos]
    h1 =
    l1 =

    [workflow-datafind]
    datafind-l1-frame-type = L1_HOFT_C01
    datafind-method = AT_RUNTIME_SINGLE_FRAMES
    datafind-check-segment-gaps = raise_error
    datafind-h1-frame-type = H1_HOFT_C01
    datafind-check-frames-exist = raise_error
    datafind-check-segment-summary = warn

    [workflow-inference]
    num-events = 1
    data-seconds-before-trigger = 220
    data-seconds-after-trigger = 1828

    [executables]
    inference = ${which:pycbc_inference}
    inference_acf = ${which:pycbc_inference_plot_acf}
    inference_acl = ${which:pycbc_inference_plot_acl}
    inference_posterior = ${which:pycbc_inference_plot_posterior}
    inference_prior = ${which:pycbc_inference_plot_prior}
    inference_rate = ${which:pycbc_inference_plot_acceptance_rate}
    inference_samples = ${which:pycbc_inference_plot_samples}
    inference_table = ${which:pycbc_inference_table_summary}
    plot_spectrum = ${which:pycbc_plot_psd_file}

    [datafind]
    urltype = file

    [inference]
    sample-rate = 2048
    low-frequency-cutoff = 40
    strain-high-pass = 30
    pad-data = 8
    psd-estimation = median
    psd-segment-length = 16
    psd-segment-stride = 8
    psd-inverse-length = 16
    processing-scheme = mkl
    sampler = kombine
    likelihood-evaluator = gaussian
    skip-burn-in =
    nwalkers = 10
    niterations = 320

    [pegasus_profile-inference]
    condor|universe = local
    condor|request_memory = 150000

    [inference_acf]
    ymax = 1.1
    ymin = -1.1

    [inference_acl]

    [inference_posterior]

    [inference_prior]

    [inference_rate]

    [inference_samples]

    [inference_table]

    [plot_spectrum]

============================
Inference configuration file
============================

=====================
Generate the workflow
=====================

When you are ready, you can generate the workflow. Here is an example::

    # remove proxy from env
    unset X509_USER_PROXY

    # name of the workflow
    WORKFLOW_NAME="r1"

    # path to output dir
    OUTPUT_DIR=output

    # input configuration files
    CONFIG_PATH=workflow.ini
    INFERENCE_CONFIG_PATH=inference.ini

If you want to run on the loudest triggers from a PyCBC coincident search workflow then run::

    # run workflow generator on triggers from workflow
    pycbc_make_inference_workflow --workflow-name ${WORKFLOW_NAME} \
        --config-files ${CONFIG_PATH} \
        --inference-config-file ${INFERENCE_CONFIG_PATH} \
        --output-dir ${OUTPUT_DIR} \
        --output-file ${WORKFLOW_NAME}.dax \
        --output-map ${OUTPUT_MAP_PATH} \
        --bank-file ${BANK_PATH} \
        --statmap-file ${STATMAP_PATH} \
        --single-detector-triggers ${SNGL_H1_PATHS} ${SNGL_L1_PATHS}

Where ``${BANK_FILE}`` is the path to the template bank HDF file, ``${STATMAP_FILE}`` is the path to the combined statmap HDF file, and ``${SNGL_H1_PATHS}`` and ``${SNGL_L1_PATHS}`` are the paths to the merged single-detector HDF files.

Else you can run from a specific GPS end time with the ``--gps-end-time`` option like::

    # run workflow generator on specific GPS end time
    pycbc_make_inference_workflow --workflow-name ${WORKFLOW_NAME} \
        --config-files ${CONFIG_PATH} \
        --inference-config-file ${INFERENCE_CONFIG_PATH} \
        --output-dir ${OUTPUT_DIR} \
        --output-file ${WORKFLOW_NAME}.dax \
        --output-map ${OUTPUT_MAP_PATH} \
        --gps-end-time ${GPS_END_TIME}

Where ``${GPS_END_TIME}`` is the GPS end time of the trigger.

Finally plan and submit the workflow with::

    # submit workflow
    pycbc_submit_dax --dax ${WORKFLOW_NAME}.dax \
        --accounting-group ligo.dev.o2.cbc.explore.test \
        --no-create-proxy

