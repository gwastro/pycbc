###############################################################################################
``pycbc_make_inference_inj_workflow``: A parameter estimation workflow generator for injections
###############################################################################################

===============
Introduction
===============

The executable ``pycbc_make_inference_inj_workflow`` is a workflow generator to setup a parameter estimation analysis.


===========================
Workflow configuration file
===========================

A sample workflow configuration file:

.. literalinclude:: ../../examples/workflow/inference/inj_workflow_config.ini
   :language: ini

:download:`Download <../../examples/workflow/inference/inj_workflow_config.ini>`


============================
Inference configuration file
============================
A sample inference configuration file:

.. literalinclude:: ../../examples/workflow/inference/inference.ini
   :language: ini

:download:`Download <../../examples/workflow/inference/inference.ini>`

=====================
Generate the workflow
=====================

To generate a workflow you will need your configuration files. Generate the workflow using following example run script:

.. literalinclude:: ../../examples/workflow/inference/run_pycbc_make_inference_inj_workflow.sh
   :language: bash

:download:`Download <../../examples/workflow/inference/run_pycbc_make_inference_inj_workflow.sh>`


=============================
Plan and execute the workflow
=============================

If you are on LDG, you need to define an accounting group. Finally plan and submit the workflow with:

::

    # submit workflow
    cd ${output_dir}
    pycbc_submit_dax --dax ${WORKFLOW_NAME}.dax \
        --no-grid \
        --enable-shared-filesystem \
        --accounting-group ${ACCOUNTING_GROUP}

