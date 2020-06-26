###############################################################################################
``pycbc_make_inference_inj_workflow``: A parameter estimation workflow generator for injections
###############################################################################################

===============
Introduction
===============

The executable ``pycbc_make_inference_inj_workflow`` is a workflow generator to
setup a parameter estimation analysis on one or more simulated signals.
Optionally, it can also run a percentile-percentile on the injections it
analyzed.

The workflow is very similar to the standard inference workflow created by
`pycbc_make_inference_workflow <pycbc_make_inference_workflow>`_. The main
differences are:

 * Rather than providing one or more ``[event-{label}]`` sections in the
   workflow config file, you provide a single ``[workflow-inference]`` section.
   The syntax for this section is very similar to the ``[event]`` section(s) in
   the standard workflow, as it sets the configuration files that are used
   by ``pycbc_inference``. The difference is that the same settings are used
   for all injections.
 * When you create the workflow, you either pass it a ``--num-injections``
   or a ``--injection-file``. If the former, the workflow will draw the
   specified number of injections from the prior given to ``pycbc_inference``
   and analyze them. If the latter, the workflow will analyze the injections
   specified in the given injection file. The file must be an HDF file;
   see ``pycbc_create_injections`` for details. In either case, each injection
   is treated as an independent event, with its own summary section in the
   results page.
 * You may optionally have the workflow do a percentile-percentile test on
   the injections. You do this by adding the necessary executables and
   corresponding sections to the ``workflow_config.ini`` file. See the example
   below for details. If a percentile-percentile test is done, the results
   page will have an additional tab that gives a summary of the PP test on
   all of the parameters, as well as PP plots and plots of injected versus
   recoverd values.
 * It is recommend (though not required) that you add
   ``plot-injection-parameters`` to the ``[plot_posterior]`` and
   ``[plot_posterior_summary]`` sections. Doing so will cause redlines to
   be plotted at the injected parameter values on the posterior plots, so
   that you may visually inspect how well the injected values are recovered.
   This may also require providing an ``injection-samples-map`` argument.
   See the example file below for details.

In the `standard workflow <pycbc_make_inference_workflow>`_ we used two
workflow configuration files, a ``workflow_config.ini`` and an ``events.ini``.
For the injection workflow, we can use the same ``workflow_config.ini``; we
just need to setup an ``injections_config.ini`` to add the needed sections
and arguments for the injections workflow.

In the example below, we demonstrate how to use the injections workflow
using the same prior and sampler settings as given in the
`standard workflow <pycbc_make_inference_workflow>`_ example.


========================================
Example: BBH injections with ``dynesty``
========================================

In this example we use the same prior and sampler settings as the example
of analyzing GW150914 and GW170814 in the
`pycbc_make_inference_workflow <pycbc_make_inference_workflow>`_
documentation. We will analyze 10 injections, and do a percentile-percentile
test on them. (This is only as an example. To do a full PP test, we recommend
using at least 100 injections.)

-------------------------------------
Get the inference configuration files
-------------------------------------

We can use the same
:download:`prior <../../examples/inference/priors/bbh-uniform_comoving_volume.ini>`,
:download:`model <../../examples/inference/models/marginalized_phase.ini>`,
and :download:`sampler <../../examples/inference/samplers/dynesty.ini>`
configuration files as used in the  
`pycbc_make_inference_workflow <pycbc_make_inference_workflow>`_ example.
However, instead of analyzing O1 or O2 data, we will create fake Gaussian
noise. To do that, we will use the
:download:`data.ini <../../examples/inference/bbh-injection/data.ini>` file
used for the `BBH simulation example <../inference/examples/bbh>`_.

-------------------------------------
Setup the workflow configuration file
-------------------------------------

As discussed above, we can use the same :download:`workflow configuration file <../../examples/workflow/inference/gw150914_gw170814-dynesty/workflow_config.ini>` as used in
the ``dynesty`` example in the standard workflow. We need to create
an ``injections_config.ini`` file to go along with the ``workflow_config.ini``:

.. literalinclude:: ../../examples/workflow/inference/bbh_inj-dynesty/injections_config.ini
   :language: ini

:download:`Download <../../examples/workflow/inference/bbh_inj-dynesty/injections_config.ini>`

---------------------
Generate the workflow
---------------------

Assuming that you have downloaded all of the configuration files to the
same directory, you can generate the workflow by running the following script:

.. literalinclude:: ../../examples/workflow/inference/bbh_inj-dynesty/create_inj_workflow.sh
   :language: bash

:download:`Download <../../examples/workflow/inference/bbh_inj-dynesty/create_inj_workflow.sh>`

Note that you need to set the ``HTML_DIR`` before running. This tells the
workflow where to save the results page when done. You can also change
``WORKFLOW_NAME`` if you like.

You should also change the ``SEED`` everytime you create a different workflow.
This sets the seed that is passed to ``pycbc_inference`` (you set it here
because it will be incremented for every ``pycbc_inference`` job that will be
run in the workflow).

After the workflow has finished it will have created a directory named
``${WORKFLOW_NAME}-output``. This contains the ``dax`` and all necessary files
to run the workflow.

-----------------------------
Plan and execute the workflow
-----------------------------

Change directory into the ``${WORKFLOW_NAME}-output`` directory::

    cd ${WORKFLOW_NAME}-output

If you are on the ATLAS cluster (at AEI Hannover) or on an LDG cluster, you
need to define an accounting group tag (talk to your cluster admins if you do
not know what this is). Once you know what accounting-group tag to use, plan
and submit the workflow with::

    # submit workflow
    pycbc_submit_dax --dax ${WORKFLOW_NAME}.dax \
        --no-grid \
        --enable-shared-filesystem \
        --accounting-group ${ACCOUNTING_GROUP}

Here, ``${ACCOUNTING_GROUP}`` is the appropriate tag for your workflow.

Once it is running, you can monitor the status of the workflow by running
``./status`` from within the ``${WORKFLOW_NAME}-output`` directory. If your
workflow fails for any reason, you can see what caused the failure by running
``./debug``. If you need to stop the workflow at any point, run ``./stop``.
To resume a workflow, run ``./start``. If the ``pycbc_inference`` jobs were
still running, and they had checkpointed, they will resume from their last
checkpoint upon restart.

------------
Results page
------------

When the workflow has completed successfully it will write out the results
page to the directory you specified in the ``create_inj_workflow.sh`` script.
You can see what the result page will look like `here <https://www.atlas.aei.uni-hannover.de/~work-cdcapano/scratch/inference_workflow_docs/bbh_injections-dynesty/>`_.
