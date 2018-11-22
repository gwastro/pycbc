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
    datafind-check-segment-summary = no_test

    [workflow-inference]
    ; how the workflow generator should setup inference nodes
    num-events = 1
    plot-1d-mass = mass1 mass2 mchirp q
    plot-1d-orientation = ra dec tc polarization inclination coa_phase
    plot-1d-distance = distance redshift

    [executables]
    ; paths to executables to use in workflow
    inference = ${which:pycbc_inference}
    inference_posterior = ${which:pycbc_inference_plot_posterior}
    inference_prior = ${which:pycbc_inference_plot_prior}
    inference_rate = ${which:pycbc_inference_plot_acceptance_rate}
    inference_samples = ${which:pycbc_inference_plot_samples}
    inference_table = ${which:pycbc_inference_table_summary}
    plot_spectrum = ${which:pycbc_plot_psd_file}
    results_page = ${which:pycbc_make_html_page}

    [datafind]
    ; datafind options
    urltype = file

    [inference]
    ; command line options use --help for more information
    sample-rate = 2048
    low-frequency-cutoff = 30
    strain-high-pass = 20
    pad-data = 8
    psd-estimation = median
    psd-segment-length = 16
    psd-segment-stride = 8
    psd-inverse-length = 16
    processing-scheme = mkl
    sampler = kombine
    likelihood-evaluator = gaussian
    nwalkers = 500
    niterations = 100000
    save-psd =

    [pegasus_profile-inference]
    ; pegasus profile for inference nodes
    condor|request_memory = 20G
    condor|request_cpus = 12

    [inference_posterior]
    ; command line options use --help for more information
    plot-density =
    plot-contours =
    plot-marginal =

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

    [results_page]
    ; command line options use --help for more information
    analysis-title = "PyCBC Inference Test"

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

    [static_args]
    ; parameters that do not vary in inference sampler
    approximant = SEOBNRv2_ROM_DoubleSpin
    f_lower = 28.0

    [prior-tc]
    ; how to construct prior distribution
    name = uniform
    min-tc = 1126259462.2
    max-tc = 1126259462.6

    [prior-mass1]
    ; how to construct prior distribution
    name = uniform
    min-mass1 = 10.
    max-mass1 = 80.

    [prior-mass2]
    ; how to construct prior distribution
    name = uniform
    min-mass2 = 10.
    max-mass2 = 80.

    [prior-distance]
    ; how to construct prior distribution
    name = uniform
    min-distance = 10
    max-distance = 500

    [prior-coa_phase]
    ; how to construct prior distribution
    name = uniform_angle
    ; uniform_angle defaults to [0,2pi), so we
    ; don't need to specify anything here

    [prior-inclination]
    ; how to construct prior distribution
    name = sin_angle

    [prior-ra+dec]
    ; how to construct prior distribution
    name = uniform_sky

    [prior-polarization]
    ; how to construct prior distribution
    name = uniform_angle

A simple configuration file for parameter estimation on the ringdown is::

    [variable_args]
    ; parameters to vary in inference sampler
    tc =
    f_0 =
    tau =
    amp =
    phi =

    [labels]
    ; LaTeX expressions to use in HTML and plotting executables
    tc = $t_c$
    f_0 = $f_0$
    tau = $\tau$
    amp = $A$
    phi = $\phi_0$

    [static_args]
    ; parameters that do not vary in inference sampler
    approximant = FdQNM
    ra = 2.21535724066
    dec = -1.23649695537
    polarization = 0.
    f_lower = 28.0
    f_final = 512

    [prior-tc]
    ; how to construct prior distribution
    name = uniform
    min-tc = 1126259462.4
    max-tc = 1126259462.5

    [prior-f_0]
    ; how to construct prior distribution
    name = uniform
    min-f_0 = 200.
    max-f_0 = 300.

    [prior-tau]
    ; how to construct prior distribution
    name = uniform
    min-tau = 0.0008
    max-tau = 0.020

    [prior-amp]
    ; how to construct prior distribution
    name = uniform
    min-amp = 0
    max-amp = 1e-20

    [prior-phi]
    ; how to construct prior distribution
    name = uniform
    min-phi = 0
    max-phi = 6.283185307179586

If you want to use another variable parameter in the inference sampler then add its name to ``[variable_args]`` and add a prior section like shown above.

=====================
Generate the workflow
=====================

To generate a workflow you will need your configuration files. We set the following enviroment variables for this example::

    # name of the workflow
    WORKFLOW_NAME="r1"

    # path to output dir
    OUTPUT_DIR=output

    # input configuration files
    CONFIG_PATH=workflow.ini
    INFERENCE_CONFIG_PATH=inference.ini

Specify a directory to save the HTML pages::

    # directory that will be populated with HTML pages
    HTML_DIR=${HOME}/public_html/inference_test

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
                           workflow:end-time:${WORKFLOW_END_TIME} \
                           workflow-inference:data-seconds-before-trigger:8 \
                           workflow-inference:data-seconds-after-trigger:8 \
                           results_page:output-path:${HTML_DIR} \
                           results_page:analysis-subtitle:${WORKFLOW_NAME}

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
        --config-overrides workflow:start-time:$((${GPS_END_TIME}-2)) \
                           workflow:end-time:$((${GPS_END_TIME}+2)) \
                           workflow-inference:data-seconds-before-trigger:2 \
                           workflow-inference:data-seconds-after-trigger:2 \
                           inference:psd-start-time:$((${GPS_END_TIME}-300)) \
                           inference:psd-end-time:$((${GPS_END_TIME}+1748)) \
                           results_page:output-path:${HTML_DIR} \
                           results_page:analysis-subtitle:${WORKFLOW_NAME}


Where ``${GPS_END_TIME}`` is the GPS end time of the trigger.

For the CBC example above define the environment variables ``GPS_END_TIME=1126259462`` and ``OUTPUT_MAP_PATH=output.map``. 

=============================
Plan and execute the workflow
=============================

Finally plan and submit the workflow with::

    # submit workflow
    pycbc_submit_dax --dax ${WORKFLOW_NAME}.dax \
        --accounting-group ligo.dev.o2.cbc.explore.test

