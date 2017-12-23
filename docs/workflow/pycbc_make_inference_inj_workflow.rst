############################################################################
``pycbc_make_inference_workflow``: A parameter estimation workflow generator
############################################################################

===============
Introduction
===============

The executable ``pycbc_make_inference_inj_workflow`` is a workflow generator to setup a parameter estimation analysis.

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

    [workflow-inference]
    ; how the workflow generator should setup inference nodes
    num-injections = 3
    plot-group-mass = mass1 mass2 mchirp q
    plot-group-orientation =  inclination polarization ra dec
    plot-group-distance = distance redshift
    plot-group-time = tc coa_phase

    [executables]
    ; paths to executables to use in workflow
    create_injections = ${which:pycbc_create_injections}
    inference = ${which:pycbc_inference}
    inference_intervals = ${which:pycbc_inference_plot_inj_intervals}
    inference_posterior = ${which:pycbc_inference_plot_posterior}
    inference_rate = ${which:pycbc_inference_plot_acceptance_rate}
    inference_recovery = ${which:pycbc_inference_plot_inj_recovery}
    inference_samples = ${which:pycbc_inference_plot_samples}
    inference_table = ${which:pycbc_inference_table_summary}
    results_page = ${which:pycbc_make_html_page}

    [create_injections]
    ; command line options use --help for more information
    ninjections = 1
    dist-section = prior

    [inference]
    ; command line options use --help for more information
    processing-scheme = mkl
    sampler = kombine
    likelihood-evaluator = gaussian
    nwalkers = 100
    n-independent-samples = 10
    checkpoint-interval = 10
    nprocesses = 24
    fake-strain = aLIGOZeroDetHighPower
    psd-model = aLIGOZeroDetHighPower
    pad-data = 8
    strain-high-pass = 20
    sample-rate = 1024
    low-frequency-cutoff = 30
    update-interval = 5
    config-overrides = static_args:approximant:TaylorF2

    [pegasus_profile-inference]
    ; pegasus profile for inference nodes
    condor|request_memory = 20G
    condor|request_cpus = 24

    [inference_intervals]
    ; command line options use --help for more information

    [inference_posterior]
    ; command line options use --help for more information
    plot-scatter =
    plot-contours =
    plot-marginal =
    z-arg = logposterior

    [inference_rate]
    ; command line options use --help for more information

    [inference_recovery]
    ; command line options use --help for more information

    [inference_samples]
    ; command line options use --help for more information

    [inference_table]
    ; command line options use --help for more information

    [results_page]
    ; command line options use --help for more information
    analysis-title = "PyCBC Inference Test"

Use the ``ninjections`` option in the ``[workflow-inference]`` section to set the number of injections in the analysis.

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

If you want to run with a test likelihood function use::

    # option for using test likelihood functions
    DATA_TYPE=analytical

Otherwise if you want to run with simulated data use::

    # option for using simulated data
    DATA_TYPE=simulated_data

If you want to run on the loudest triggers from a PyCBC coincident search workflow then run::

    # run workflow generator on simulated data
    pycbc_make_inference_inj_workflow \
        --workflow-name ${WORKFLOW_NAME} \
        --data-type ${DATA_TYPE} \
        --output-dir output \
        --output-file ${WORKFLOW_NAME}.dax \
        --inference-config-file ${INFERENCE_CONFIG_PATH} \
        --config-files ${CONFIG_PATH} \
        --config-overrides results_page:output-path:${HTML_DIR} \
                           workflow:start-time:${GPS_START_TIME} \
                           workflow:end-time:${GPS_END_TIME}

Where ``${GPS_START_TIME}`` and ``${GPS_END_TIME}`` are the GPS times of data to read.

=============================
Plan and execute the workflow
=============================

Finally plan and submit the workflow with::

    # submit workflow
    pycbc_submit_dax --dax ${WORKFLOW_NAME}.dax \
        --accounting-group ligo.dev.o2.cbc.explore.test

