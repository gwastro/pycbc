######################################################
Workflow Generator for Doing Faithfulness Comparisons
######################################################

============
Introduction
============

This page describes how to use ``pycbc_make_faithsim_workflow``, the
workflow generator within PyCBC for comparing the faithfulness of pairs of
waveform approximants.

The workflow:

 #. Generates a population of injections with ``lalapps_inspinj`` or
    ``pycbc_create_injections``.
 #. Splits the injection set into a number of smaller pieces with
    whichever of ``pycbc_split_inspinj``/``pycbc_hdf_splitinj``/
    ``pycbc_splitbank``/``pycbc_hdf5_splitbank`` is configured.
 #. Runs ``pycbc_faithsim`` on each piece, computing the match between the
    two configured approximants for every injection.
 #. Collects the per-piece results into a single results file with
    ``pycbc_faithsim_collect_results``.
 #. Generates summary plots with ``pycbc_faithsim_plots``.

Like all PyCBC workflow generators, this creates a Pegasus DAX that is
planned and run using HTCondor. A Condor pool (or a submission point to
one, such as an LDG cluster) is required to actually run the workflow;
this page only covers generating and submitting it.

.. note::

    ``pycbc_make_faithsim_workflow`` replaces the older
    ``pycbc_make_faithsim`` script, which has been removed. If you have
    configuration files written for the old script, note that the section
    names and options below are different and the old files are not
    directly compatible.

==========================
How to generate a workflow
==========================

------------------------------------
Creating a configuration (.ini) file
------------------------------------

All the choices when setting up a faithsim workflow are contained in one
or more configuration (.ini) files, which are merged together when passed
via ``--config-files``.

Below is a complete, working example, defaulting to HDF for the injection
set (LIGOLW XML works equally well; use ``lalapps_inspinj`` and
``pycbc_split_inspinj``/``pycbc_splitbank`` in ``[executables]`` instead
to work with XML).

.. literalinclude:: ../examples/faith/faithsim_workflow_config.ini

The important sections are:

 #. ``[workflow]``

    Must contain ``start-time``/``end-time``, which set the GPS time
    range over which injections are placed (spaced according to
    ``time-step`` in ``[injection]`` below); together with ``time-step``
    they determine how many injections are generated.

 #. ``[injection]``

    Options for generating the injection set, passed to
    ``lalapps_inspinj`` or ``pycbc_create_injections`` (selected by the
    ``injection`` entry in ``[executables]``).

    With ``pycbc_create_injections`` (the default, needed for HDF
    injections), the parameter distributions themselves are not given
    here: ``config-files`` instead points at a separate, standard PyCBC
    distributions config file (``[variable_params]``/``[prior-XXX]``/
    ``[static_params]`` sections; see :ref:`inference` for the format).
    Below is a working example.

    .. literalinclude:: ../examples/faith/injection_priors.ini

    ``config-files`` must be an absolute path: the workflow generator
    changes into ``--output-dir`` before resolving it. Do not give a
    value for ``tc``; the workflow sets injection times itself via
    ``time-step``/``time-window``. Note that the approximants actually
    compared are set in ``[pycbc_faithsim]`` below, not here.

    With ``lalapps_inspinj`` (XML injections only), options are instead
    given directly in ``[injection]`` following that program's own
    command-line options (``min-mass1``, ``i-distr``, ``waveform``, etc.).

 #. ``[executables]``

    The paths to the executables used by the workflow: ``injection``,
    ``injsplit``, ``pycbc_faithsim``, ``pycbc_faithsim_collect_results``
    and ``pycbc_faithsim_plots``. The ``${which:program-name}`` syntax
    shown above resolves each path via the shell ``PATH``, so it is
    normally not necessary to hard-code full paths.

 #. ``[splitbank]``

    Contains ``num_banks``, the number of pieces the injection set is
    split into (by whichever executable is configured as ``injsplit``).
    Each piece is processed by an independent ``pycbc_faithsim`` job, so
    this option is directly proportional to the number of jobs (and
    inversely proportional to the running time of each).

 #. ``[pycbc_faithsim]``

    The options passed to the ``pycbc_faithsim`` executable: the PSD, and
    the two approximants (``waveform1-approximant``/
    ``waveform2-approximant``) and their starting frequencies to be
    compared, along with the filtering settings.

 #. ``[pycbc_faithsim_plots]`` and ``[pycbc_faithsim_plots-XXX]``

    ``[pycbc_faithsim_plots]`` sets options common to all summary plots
    (e.g. ``colormap``). Each ``[pycbc_faithsim_plots-XXX]`` subsection
    adds one more plot of the results, giving the parameters to use for the
    x-, y- and (optionally) z-axes.

 #. ``[pegasus_profile]``

    Sets Condor/Pegasus properties for the workflow. Faithsim workflows
    running on LDG clusters must include an ``accounting-group`` value
    valid for your account; the value must be chosen according to the
    `Accounting information web page <https://ldas-gridmon.ligo.caltech.edu/ldg_accounting/>`_.

------------------------
Generating the workflow
------------------------

Once a configuration file has been made, create a workspace directory and
place the file (and this script) into it. Running the following will
generate and submit the workflow:

.. literalinclude:: ../examples/faith/run_workflow.sh
   :language: bash

Alternatively, generate the workflow without submitting it by omitting
``--submit-now``:

.. code-block:: bash

    pycbc_make_faithsim_workflow \
        --workflow-name my_faithsim_workflow \
        --config-files faithsim_workflow_config.ini \
        --output-dir output

This will create the ``output`` directory containing the Pegasus DAX and
supporting files, along with ``status``/``debug``/``stop``/``start``
helper scripts once the workflow has been planned.

-------------------------
Understanding the results
-------------------------

The main result of the workflow is a single results file, produced by
``pycbc_faithsim_collect_results``. This is a whitespace-separated ASCII
file with one row per injection, giving the match, overlap, time offset,
and the recovered/injected parameters.

Summary plots (as configured in ``[pycbc_faithsim_plots-XXX]``) are also
generated automatically.

-------------------------------
Running pycbc_faithsim directly
-------------------------------

For a single, ad-hoc comparison (e.g. testing one pair of approximants
against a handful of injections) that doesn't need a full workflow,
``pycbc_faithsim`` can be run directly on a pre-made parameter file.

.. literalinclude:: ../examples/faith/run.sh
   :language: bash

------------------------------
Notes on recent fixes
------------------------------

A number of issues affecting this workflow have recently been fixed:

* ``pycbc_faithsim`` and ``pycbc_faithsim_collect_results`` can now read
  HDF parameter/injection files, in addition to LIGOLW XML.
* The workflow generator previously crashed immediately due to a
  variable-name typo, and separately never created its own output
  directory; both are fixed.
* ``pycbc_faithsim_collect_results`` previously crashed due to a
  ``NameError`` and two incorrect argument names; these are fixed, and it
  had never previously run successfully.
