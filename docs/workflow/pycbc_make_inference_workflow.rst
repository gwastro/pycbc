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
    h1-channel-name = H1:GDS-FAKE_STRAIN
    l1-channel-name = L1:OAF-CAL_DARM_DQ

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

Where ``--gps-end-time`` is the GPS end time of the trigger.

Finally plan and submit the workflow with::

    # submit workflow
    pycbc_submit_dax --dax ${WORKFLOW_NAME}.dax \
        --accounting-group ligo.dev.o2.cbc.explore.test \
        --no-create-proxy

