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

A sample workflow configuration file:

.. literalinclude:: ../../examples/workflow/inference/workflow_config.ini
    :language: ini

:download:`Download <../../examples/workflow/inference/workflow_config.ini>`

============================
Inference configuration file
============================

You will also need a configuration file with sections that tells ``pycbc_inference`` how to construct the priors. A sample inference configuration file is:

.. literalinclude:: ../../examples/workflow/inference/inference.ini
    :language: ini

:download:`Download <../../examples/workflow/inference/inference.ini>`

A sample configuration file for parameter estimation on the ringdown is:

.. literalinclude:: ../../examples/workflow/inference/ringdown_inference.ini
    :language: ini

:download:`Download <../../examples/workflow/inference/ringdown_inference.ini>`

If you want to use another variable parameter in the inference sampler then add its name to ``[variable_args]`` and add a prior section like shown above.

=====================
Generate the workflow
=====================

To generate a workflow you will need your configuration files. If you want to run on the loudest triggers from a PyCBC coincident search workflow then run:

.. literalinclude:: ../../examples/workflow/inference/run_pycbc_make_inference_workflow.sh
   :language: bash

:download:`Download <../../examples/workflow/inference/run_pycbc_make_inference_workflow.sh>`

Where ``${BANK_FILE}`` is the path to the template bank HDF file, ``${STATMAP_FILE}`` is the path to the combined statmap HDF file, ``${SNGL_H1_PATHS}`` and ``${SNGL_L1_PATHS}`` are the paths to the merged single-detector HDF files,  and ``${WORKFLOW_START_TIME}`` and ``${WORKFLOW_END_TIME}`` are the start and end time of the coincidence workflow.

Else you can run from a specific GPS end time with the ``--gps-end-time`` option like:

.. literalinclude:: ../../examples/workflow/inference/run_pycbc_make_inference_workflow_2.sh
   :language: bash

:download:`Download <../../examples/workflow/inference/run_pycbc_make_inference_workflow_2.sh>`

Where ``${GPS_END_TIME}`` is the GPS end time of the trigger.

For the CBC example above define the environment variables ``GPS_END_TIME=1126259462`` and ``OUTPUT_MAP_PATH=output.map``. 

=============================
Plan and execute the workflow
=============================

If you are on LDG, you need to define an accounting group. Plan and submit the workflow with::

    # submit workflow
    cd ${OUTPUT_DIR}
    pycbc_submit_dax --dax ${WORKFLOW_NAME}.dax \
        --no-grid \
        --enable-shared-filesystem \
        --accounting-group ${ACCOUNTING_GROUP}

Where ``${ACCOUNTING_GROUP}`` is the appropriate tag for your workflow.

