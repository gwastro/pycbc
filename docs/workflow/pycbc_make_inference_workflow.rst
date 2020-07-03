############################################################################
``pycbc_make_inference_workflow``: A parameter estimation workflow generator
############################################################################

===============
Introduction
===============

The executable ``pycbc_make_inference_workflow`` is a workflow generator to
setup a parameter estimation analysis. It can be setup to run on one or more
events at once. For each event, the workflow:

 #. Runs ``pycbc_inference``. If desired, you can run multiple independent
    instances of ``pycbc_inference`` on the same event.
 #. Extracts a posterior file using ``pycbc_inference_extract_samples``. If
    multiple instances of ``pycbc_inference`` were run on the same event, the
    samples from all of the runs will be combined into a single posterior file.
    You can also have derived parameters written out to the posterior file.
 #. Makes various posterior plots and tables. The prior is also plotted. If
    you are analyzing gravitational-wave data, a plot of power spectral density
    (PSD) used for each event is also created.
 #. If you are working in a Python 3.x environment you can optionally have
    the workflow produce a skymap for each event (this requires ``ligo.skymap``
    to be installed).
 #. Optionally creates sampler-dependent diagnostic plots.
 #. Generates a results html page that gathers all of the results.

The workflow generator requires a configuration file that tells it what plots
to make, what parameters to produce posteriors for, which events to analyze,
and any other settings to use for the various executables that are run.

For each event, one or more inference configuration files (the file(s) passed
to ``pycbc_inference``) must also be provided. These are separate from the
workflow configuration file, as they describe how to analyze each event.  You
tell the workflow how many events to analyze and which inference configuration
files to use for each event via ``[event-{label}]`` sections in the workflow
configuration file. Here, ``{label}`` is a unique label for each event.

To illustrate how to setup and use a workflow, below we provide an example
of how to setup the workflow to analyze two binary black hole events at once
-- GW150914 and GW170814.


================================================
Example: GW150914 and GW170814 with ``emcee_pt``
================================================

In this example we setup a workflow to analyze GW150914 and GW170814 using
``emcee_pt``. We will use a prior that is uniform in comoving volume and
uniform in source masses. As we will be using the ``IMRPhenomPv2`` waveform
approximant, we will use the ``marginalized_phase`` Gaussian noise model.

This workflow will produce a results page that looks like the example
`here <https://www.atlas.aei.uni-hannover.de/~work-cdcapano/scratch/inference_workflow_docs/inference/inference-gw150914_gw170814/>`_.

The inference configuration files we will use can all be found in the pycbc
``examples`` directory. Below, we provide instructions on what files need
to be downloaded, and how to setup and run the workflow.


-------------------------------------
Get the inference configuration files
-------------------------------------

We need the configuration files for ``pycbc_inference``. These define the
prior, model, sampler, and data to use for each event.

**The prior:**

.. literalinclude:: ../../examples/inference/priors/bbh-uniform_comoving_volume.ini 
   :language: ini

:download:`Download <../../examples/inference/priors/bbh-uniform_comoving_volume.ini>`

**The model:**

.. literalinclude:: ../../examples/inference/models/marginalized_phase.ini
   :language: ini

:download:`Download <../../examples/inference/models/marginalized_phase.ini>`

**The sampler:**

.. literalinclude:: ../../examples/inference/samplers/emcee_pt-srcmasses_comoving_volume.ini 
   :language: ini

:download:`Download <../../examples/inference/samplers/emcee_pt-srcmasses_comoving_volume.ini>`

**The data:** We also need configuration files for the data. Since GW150914
occured during O1 while GW170814 occurred during O2, we need both the standard
O1 and O2 files:

.. literalinclude:: ../../examples/inference/data/o1.ini
   :language: ini

:download:`Download <../../examples/inference/data/o1.ini>`

.. literalinclude:: ../../examples/inference/data/o2.ini
   :language: ini

:download:`Download <../../examples/inference/data/o2.ini>`


-------------------------------------
Setup the workflow configuration file
-------------------------------------

As discussed above, the workflow configuration file specifes what events to
analyze, what programs to run, and what settings to use for those programs.
Since the same general workflow settings can be used for different classes of
events, here we have split the workflow configuration file into two separate
files, ``events.ini`` and ``workflow_config.ini``. The former specifies what
events we are analyzing in this run, while the latter specifies all of the
other settings. As we will see below, we can simply provide these two files to
``pycbc_make_inference_workflow``'s ``--config-file`` argument; it will
automatically combine them into a single file.

The events:

.. literalinclude:: ../../examples/workflow/inference/gw150914_gw170814-emcee_pt/events.ini 
   :language: ini

:download:`Download <../../examples/workflow/inference/gw150914_gw170814-emcee_pt/events.ini>`

The rest of the configuration file:

.. literalinclude:: ../../examples/workflow/inference/gw150914_gw170814-emcee_pt/workflow_config.ini 
   :language: ini

:download:`Download <../../examples/workflow/inference/gw150914_gw170814-emcee_pt/workflow_config.ini>`

**Notes**:

 * Since the ``[executables]`` section contains entries for
   ``create_fits_file`` and ``plot_skymap``, the workflow will try to create
   sky maps. **This requires a Python 3.x environment and** ``ligo.skymap``
   **to be installed.** If you have not installed ``ligo.skymap`` yet, do so by
   running::

        pip install ligo.skymap

 * If you do not want to create sky maps, or are running a Python 2.7
   environment, you can turn this off by simply commenting out or removing
   ``create_fits_file`` and ``plot_skymap`` from the ``[executables]`` section.

 * The number of cores that will be used by ``pycbc_inference`` is set by the
   ``nprocesses`` argument in the ``[inference]`` section. You should set this
   to the number of cores you expect to be able to get on your cluster. In the
   configuration presented here, we are limited to shared memory cores. (It
   is possible to run using MPI in order to parallelize over a larger number
   of cores, but that requires special condor settings that must be implemented
   by your cluster admins. That is outside the scope of these instructions.)

 * Notice that the number of processes that ``pycbc_inference`` will use is
   referenced by the ``condor|request_cpus`` argument in the
   ``[pegasus_profile-inference]`` section. This argurment is what tells
   condor how many cores to assign to the job, and so sets the actual number
   of resources ``pycbc_inference`` will get. Generally, you want this to
   be the same as what is fed to ``pycbc_inference``'s ``nprocesses``
   option.

The ``workflow_config.ini`` file can be used with any of the MCMC samplers when
analyzing a gravitational wave that involves the parameters mentioned in the
file. If you wanted to analyze other binary black holes, you could use this
same file, simply changing the ``events.ini`` file to point to the events
you want to analyze.


---------------------
Generate the workflow
---------------------

Assuming that you have downloaded all of the configuration files to the
same directory, you can generate the workflow by running the following script:

.. literalinclude:: ../../examples/workflow/inference/gw150914_gw170814-emcee_pt/create_workflow.sh 
   :language: bash

:download:`Download <../../examples/workflow/inference/gw150914_gw170814-emcee_pt/create_workflow.sh>`

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
page to the directory you specified in the ``create_workflow.sh`` script.
You can see what the result page will look like `here <https://www.atlas.aei.uni-hannover.de/~work-cdcapano/scratch/inference_workflow_docs/inference/inference-gw150914_gw170814/>`_.


===============================================
Example: GW150914 and GW170814 with ``dynesty``
===============================================

In this example, we repeat the above analysis, but using the `dynesty`
sampler. We can use the same
:download:`prior <../../examples/inference/priors/bbh-uniform_comoving_volume.ini>`,
:download:`model <../../examples/inference/models/marginalized_phase.ini>`,
and :download:`o1 <../../examples/inference/data/o1.ini>` and 
:download:`o2 <../../examples/inference/data/o2.ini>` inference configuration
files as above. New files that we need are:

 * The sampler configuration file for ``dynesty``:

.. literalinclude:: ../../examples/inference/samplers/dynesty.ini 
   :language: ini

:download:`Download <../../examples/inference/samplers/dynesty.ini>`

 * An ``events`` file which uses ``dynesty``:

.. literalinclude:: ../../examples/workflow/inference/gw150914_gw170814-dynesty/events.ini
   :language: ini

:download:`Download <../../examples/workflow/inference/gw150914_gw170814-dynesty/events.ini>`

Note that here, we are not running ``pycbc_inference`` multiple times. This is
because a single run of ``dynesty`` with the settings we are using (2000 live
points) produces a large number of (O(10 000)) samples.

We also need a slightly different
:download:`workflow configuration file <../../examples/workflow/inference/gw150914_gw170814-dynesty/workflow_config.ini>`. The only difference from the workflow configuration file from the one above
is that the diagnostic plot executable have been removed
(``plot_acceptance_rate`` and ``plot_samples``). This is because these
diagnostics do not work for ``dynesty``, a nested sampler. As above, **set the
nprocesses argument in the** ``[inference]`` **section to the number of cores that
works for your cluster.***

Note that we could have run both the ``emcee_pt`` analysis, above, and the
``dynesty`` analysis together in a single workflow. However, to do so, we would
need to remove any diagnostic plots that are unique to each sampler.

Once you have downloaded the necessary files, create the workflow and launch
it using the same ``create_workflow.sh`` script and ``pycbc_submit_dax``
commands as above, making sure to change the ``WORKFLOW_NAME`` and ``SEED``.

This will produce a results page that looks like the example
`here <https://www.atlas.aei.uni-hannover.de/~work-cdcapano/scratch/inference_workflow_docs/inference/inference-dynesty-gw150914_gw170814/>`_.
