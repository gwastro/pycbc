#################################################################
Workflow Generator for Verifying Template Banks (Bank Verifier)
#################################################################

.. _bank_verifier:

============
Introduction
============

``pycbc_make_bank_verifier_workflow`` is a workflow generator that verifies
the ability of a template bank to recover sets of injections, by measuring
the fitting factor (effectualness) achieved by the bank across the
parameter space it was designed for. It is the replacement for the older
``pycbc_make_banksim`` script, which has been removed.

Like all PyCBC workflow generators, this creates a Pegasus DAX that is
planned and run using HTCondor. A Condor pool (or a submission point to
one, such as an LDG cluster) is required to actually run the workflow;
this page only covers generating and submitting it.

The workflow verifies a bank against two kinds of injection sets, both
configured in the same way:

* **Point-injection sets**: small, targeted populations (e.g. drawn from a
  narrow region around a single mass pair) used to check recovery at
  specific, interesting points in parameter space.
* **Broad-injection sets**: larger populations covering a wide region of
  parameter space, used to check overall bank coverage.

For each configured injection set, the workflow:

 #. Splits the template bank into pieces with ``pycbc_splitbank`` or
    ``pycbc_hdf5_splitbank``.
 #. Generates the injection set with ``lalapps_inspinj`` or
    ``pycbc_create_injections``.
 #. Splits the injection set into pieces with ``pycbc_split_inspinj`` or
    ``pycbc_hdf5_splitbank``.
 #. Computes the fitting factor of every injection against every bank
    piece with ``pycbc_banksim``.
 #. Combines the per-bank-piece results for each injection with
    ``pycbc_banksim_combine_banks``.
 #. Combines the per-injection results into one file per injection set
    with ``pycbc_banksim_match_combine``.

Once every injection set has been processed, the workflow:

 #. Builds a summary table of the point-injection sets with
    ``pycbc_banksim_table_point_injs``.
 #. Generates fitting-factor plots for the point- and broad-injection sets
    with ``pycbc_banksim_plot_fitting_factors`` and
    ``pycbc_banksim_plot_eff_fitting_factor``.
 #. Assembles an HTML results page.

If no point-injection sets are configured, steps that only apply to them
are skipped rather than planned with no input (see :ref:`bank_verifier_notes`
below).

.. note::

    The structure of this page follows an earlier community-contributed
    draft written by GitHub user MPillas:
    `pycbc_make_bank_verifier_workflow.rst
    <https://github.com/MPillas/banksim-documentation/blob/main/pycbc_make_bank_verifier_workflow.rst>`_.
    Thank you for that write-up. The content here has been revised to
    match, and has been verified against, the current implementation.

==========================
How to generate a workflow
==========================

------------------------------------
Creating a configuration (.ini) file
------------------------------------

All the choices when setting up a bank verifier workflow are contained in
one or more configuration (.ini) files, which are merged together when
passed via ``--config-files``. Splitting the configuration across several
files is useful when different injection sets need different parameter
distributions (see below); a single file is enough for simple cases.

Below is a complete, working example, defaulting to HDF for both the
template bank and the injection sets (LIGOLW XML works equally well for
both; use ``pycbc_splitbank``/``pycbc_split_inspinj`` instead of
``pycbc_hdf5_splitbank``/``pycbc_hdf_splitinj`` in ``[executables]`` to
work with XML instead).

.. literalinclude:: ../examples/banksim/bank_verifier_config.ini

The important sections are:

 #. ``[workflow]``

    Must contain ``input-bank``, the path to the template bank to verify
    (a ``sngl_inspiral`` table, in either LIGOLW XML or HDF format), and
    ``start-time``/``end-time``, required when generating injections with
    ``pycbc_create_injections`` (not otherwise used). ``input-bank`` must
    be an absolute path: the workflow generator changes into
    ``--output-dir`` before resolving it.

 #. ``[injection]``

    Options for generating the injection sets, passed to
    ``lalapps_inspinj`` or ``pycbc_create_injections`` (selected by the
    ``injection`` entry in ``[executables]``). These apply to every
    injection set unless overridden in a tagged section named
    ``[injection-TAG]``, where ``TAG`` is one of the names given in
    ``[workflow-pointinjs]``/``[workflow-broadinjs]`` below.

    With ``pycbc_create_injections`` (the default, needed for HDF
    injections), the parameter distributions themselves are not given
    here: ``config-files`` instead points at a separate, standard PyCBC
    distributions config file (``[variable_params]``/``[prior-XXX]``/
    ``[static_params]`` sections; see :ref:`inference` for the format).
    Below is a working example.

    .. literalinclude:: ../examples/banksim/injection_priors.ini

    This is how different injection sets (e.g. distinct mass regions) are
    given different parameter distributions: give each set its own
    ``[injection-TAG]`` section overriding ``config-files`` to point at a
    different distributions file. ``config-files`` must be an absolute
    path, for the same reason as ``input-bank`` above. Do not give a value
    for ``tc`` in the distributions file; the workflow sets the injection
    times itself via ``time-step``/``time-window``.

    With ``lalapps_inspinj`` (XML injections only), options are instead
    given directly in ``[injection]``/``[injection-TAG]`` following that
    program's own command-line options (``min-mass1``, ``i-distr``,
    ``waveform``, etc.)

 #. ``[executables]``

    The paths to every executable used by the workflow. The
    ``${which:program-name}`` syntax shown above resolves each path via
    the shell ``PATH``.

 #. ``[workflow-pointinjs]`` and ``[workflow-broadinjs]``

    One ``tag = number-of-injections`` entry per named injection set. Any
    number of sets (including zero) may be given in either section.

 #. ``[workflow-splittable-shortinjbanksplit]`` and
    ``[workflow-splittable-broadinjbanksplit]``

    Control how the template bank is split up for the point- and
    broad-injection banksim jobs respectively. ``splittable-exe-tag``
    names the entry in ``[executables]`` to use for splitting (dispatched
    by executable basename to ``pycbc_splitbank`` or
    ``pycbc_hdf5_splitbank``); ``splittable-num-banks`` sets how many
    pieces the bank is split into. Increase ``splittable-num-banks`` if
    jobs are running for too long.

 #. ``[workflow-splittable-shortinjs]`` and
    ``[workflow-splittable-broadinjs]``

    The same, but controlling how each injection set is split up before
    being distributed to banksim jobs.

 #. ``[banksim]``

    The options passed to every ``pycbc_banksim`` job: the PSD, the
    template and signal approximants and their starting frequencies, and
    the filtering settings. See ``pycbc_banksim --help`` for the full set
    of options, including ``--mchirp-window``/``--tau0-window`` to
    restrict which templates are checked against each injection.

    As of this release, ``--signal-file`` (and ``--template-file``) accept
    either an HDF or LIGOLW XML file, so ``splittable-exe-tag`` may point
    at ``pycbc_hdf5_splitbank`` for either the injection or the bank
    splitting step.

 #. ``[banksim_plot_eff_fitting_fac-XXX]`` and
    ``[banksim_plot_fitting_factors-XXX]``

    Each subsection adds one more summary plot: a mass1-vs-mass2 effective
    fitting-factor plot, or a fitting-factor distribution plot,
    respectively (the quantities plotted are fixed, not configurable).
    Useful options include ``plot-title``/``plot-caption``/``log-axes``/
    ``log-colorbar``/``filter-injections``; see
    ``pycbc_banksim_plot_eff_fitting_factor --help``/
    ``pycbc_banksim_plot_fitting_factors --help`` for the full set.

 #. ``[results_page]``

    Must contain ``output-path``, the directory (relative to the
    workflow's output directory) that the HTML results page is written
    into.

 #. ``[pegasus_profile]``

    Sets Condor/Pegasus properties for the workflow. Bank verifier
    workflows running on LDG clusters must include an ``accounting-group``
    valid for your account; the value must be chosen according to the
    `Accounting information web page <https://ldas-gridmon.ligo.caltech.edu/ldg_accounting/>`_.
    This is not required for non-LVK users.

------------------------------
Generating and submitting
------------------------------

Once a configuration file has been made, create a workspace directory and
place the file (and this script) into it. Running the following will
generate and submit the workflow:

.. literalinclude:: ../examples/banksim/run_verifier.sh
   :language: bash

Set ``output-dir`` to a suitable location before running, and change
``WORKFLOW_NAME`` if desired. Once planned, Pegasus creates
``status``/``debug``/``stop``/``start`` helper scripts in the output
directory: ``./status`` reports progress, ``./debug`` helps diagnose
failed jobs, and ``./stop``/``./start`` halt and resume the running
workflow.

-------------
Results page
-------------

The results page is written to the ``output-path`` given in
``[results_page]``, inside the workflow's output directory. It includes
the summary table and fitting-factor plots for the point-injection sets,
and the fitting-factor plots for the broad-injection sets.

.. _bank_verifier_notes:

------------------------------
Notes on recent fixes (#4931)
------------------------------

A few issues affecting the bank verifier workflow and ``pycbc_banksim``
have recently been fixed:

* ``pycbc_banksim`` can now read injections from an HDF file via
  ``--signal-file``, in addition to LIGOLW XML, so injection sets no
  longer need to be split to XML.
* If no point-injection sets are configured (``[workflow-pointinjs]`` is
  empty), the summary table and eff. fitting factor plotting jobs are no
  longer planned, rather than being planned with no input files and
  failing when run.
* When waveform generation fails inside ``pycbc_banksim``, the offending
  approximant and parameters are now logged before the job exits, to make
  it easier to identify problematic points in a large bank or injection
  set without needing to reproduce the failure interactively.
* ``pycbc_banksim`` previously crashed with ``AttributeError`` whenever it
  tapered a waveform for an injection that had a ``taper`` value set (a
  common setting for real production injection sets); this is fixed.
