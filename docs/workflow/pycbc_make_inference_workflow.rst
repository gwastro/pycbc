############################################################################
``pycbc_make_inference_workflow``: A parameter estimation workflow generator
############################################################################

===============
Introduction
===============

The executable ``pycbc_make_inference_workflow`` is a workflow generator to setup a parameter estimation analysis.

===========================
Workflow configuration file
===========================

A simple workflow configuration file::

    [workflow]
    ; basic information used by the workflow generator
    file-retention-level = all_triggers
    h1-channel-name = H1:DCS-CALIB_STRAIN_C02
    l1-channel-name = L1:DCS-CALIB_STRAIN_C02

    [workflow-ifos]
    ; the IFOs to analyze
    h1 =
    l1 =

    [workflow-datafind]
    ; how the workflow generator should get frame data
    datafind-h1-frame-type = H1_HOFT_C02
    datafind-l1-frame-type = L1_HOFT_C02
    datafind-method = AT_RUNTIME_SINGLE_FRAMES
    datafind-check-segment-gaps = raise_error
    datafind-check-frames-exist = raise_error
    datafind-check-segment-summary = warn

    [workflow-inference]
    ; how the workflow generator should setup inference nodes
    num-events = 1

    [executables]
    ; paths to executables to use in workflow
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
    ; datafind options
    urltype = file

    [inference]
    ; command line options use --help for more information
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
    ; pegasus profile for inference nodes
    condor|universe = local
    condor|request_memory = 150000

    [inference_acf]
    ; command line options use --help for more information
    ymax = 1.1
    ymin = -1.1

    [inference_acl]
    ; command line options use --help for more information

    [inference_posterior]
    ; command line options use --help for more information

    [inference_prior]
    ; command line options use --help for more information

    [inference_rate]
    ; command line options use --help for more information

    [inference_samples]
    ; command line options use --help for more information

    [inference_table]
    ; command line options use --help for more information

    [plot_spectrum]
    ; command line options use --help for more information

============================
Inference configuration file
============================

You will also need a configuration file with sections that tells ``pycbc_inference`` how to construct the priors. A simple inference configuration file is::

    [variable_args]
    ; parameters to vary in inference sampler
    tc =
    mass1 =
    mass2 =
    distance =
    coa_phase =
    inclination =
    ra =
    dec =
    polarization =

    [labels]
    ; LaTeX expressions to use in HTML and plotting executables
    tc = $t_{c}$
    mass1 = $m_{1}$
    mass2 = $m_{2}$
    distance = $d$
    coa_phase = $\phi_{c}$
    inclination = $\iota$
    ra = $\alpha$
    dec = $\delta$

    [static_args]
    ; parameters that do not vary in inference sampler
    approximant = TaylorF2
    f_lower = 40.0

    [prior-tc]
    ; how to construct prior distribution
    name = uniform
    min-tc = 1137215767.92
    max-tc = 1137215768.04

    [prior-mass1]
    ; how to construct prior distribution
    name = uniform
    min-mass1 = 1.3
    max-mass1 = 10.0

    [prior-mass2]
    ; how to construct prior distribution
    name = uniform
    min-mass2 = 1.3
    max-mass2 = 10.0

    [prior-distance]
    ; how to construct prior distribution
    name = uniform
    min-distance = 30.0
    max-distance = 100.0

    [prior-coa_phase]
    ; how to construct prior distribution
    name = uniform_angle
    ; uniform_angle defaults to [0,2pi), so we
    ; don't need to specify anything here

    [prior-inclination]
    ; how to construct prior distribution
    name = uniform_angle
    ; inclination between 0 and pi
    min-inclination = 0
    max-inclination = 1

    [prior-ra+dec]
    ; how to construct prior distribution
    name = uniform_sky

    [prior-polarization]
    ; how to construct prior distribution
    name = uniform_angle

If you want to use another variable parameter in the inference sampler then add its name to ``[variable_args]`` and add a prior section like shown above.

=====================
Generate the workflow
=====================

To generate a workflow you will need your configuration files. We set the following enviroment variables for this example::

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
        --config-overrides workflow:start-time:${WORKFLOW_START_TIME} \
                           workflow:end-time:$((${WORKFLOW_END_TIME} \
                           workflow-inference:data-seconds-before-trigger:1024 \
                           workflow-inference:data-seconds-after-trigger:1024

Where ``${BANK_FILE}`` is the path to the template bank HDF file, ``${STATMAP_FILE}`` is the path to the combined statmap HDF file, ``${SNGL_H1_PATHS}`` and ``${SNGL_L1_PATHS}`` are the paths to the merged single-detector HDF files,  and ``${WORKFLOW_START_TIME}`` and ``${WORKFLOW_END_TIME}`` are the start and end time of the coincidence workflow.

Else you can run from a specific GPS end time with the ``--gps-end-time`` option like::

    # run workflow generator on specific GPS end time
    pycbc_make_inference_workflow --workflow-name ${WORKFLOW_NAME} \
        --config-files ${CONFIG_PATH} \
        --inference-config-file ${INFERENCE_CONFIG_PATH} \
        --output-dir ${OUTPUT_DIR} \
        --output-file ${WORKFLOW_NAME}.dax \
        --output-map ${OUTPUT_MAP_PATH} \
        --gps-end-time ${GPS_END_TIME} \
        --config-overrides workflow:start-time:$((${GPS_END_TIME}-1024)) \
                           workflow:end-time:$((${GPS_END_TIME}+1024)) \
                           workflow-inference:data-seconds-before-trigger:1024 \
                           workflow-inference:data-seconds-after-trigger:1024

Where ``${GPS_END_TIME}`` is the GPS end time of the trigger.

=============================
Plan and execute the workflow
=============================

Finally plan and submit the workflow with::

    # submit workflow
    pycbc_submit_dax --dax ${WORKFLOW_NAME}.dax \
        --accounting-group ligo.dev.o2.cbc.explore.test \
        --no-create-proxy

