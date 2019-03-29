################################################################################
``pycbc_make_inference_inj_workflow``: A parameter estimation workflow generator
################################################################################

===============
Introduction
===============

The executable ``pycbc_make_inference_inj_workflow`` is a workflow generator to setup a parameter estimation analysis.


===========================
Workflow configuration file
===========================

A sample workflow configuration file::
.. literalinclude:: ../examples/workflow/inference_inj/workflow_config.ini

:download:`Download <../examples/workflow/inference_inj/workflow_config.ini>`


============================
Inference configuration file
============================
A sample inference configuration file::
.. literalinclude:: ../examples/workflow/inference_inj/inference.ini

:download:`Download <../examples/workflow/inference_inj/inference.ini>`

=====================
Generate the workflow
=====================

To generate a workflow you will need your configuration files. Generate the workflow using following example run script::
A sample inference configuration file::
.. literalinclude:: ../examples/workflow/inference_inj/run_pycbc_make_inference_inj_workflow.sh

:download:`Download <../examples/workflow/inference_inj/run_pycbc_make_inference_inj_workflow.sh>`


=============================
Plan and execute the workflow
=============================

Finally plan and submit the workflow with::

    # submit workflow
    cd ${OUT_DIR}
    pycbc_submit_dax --dax ${WORKFLOW_NAME}.dax \
        --accounting-group ${ACCOUNTING_GROUP}

