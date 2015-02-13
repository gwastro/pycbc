.. _workflowpostprocprepmod:

###############################################
The workflow postprocessing preparation module
###############################################

=============
Introduction
=============

The postprocessing preparation module is used to prepare the output of the coincidence stage for the postprocessing (ie. calculation of trigger significance and of rate statements). For the case of pipedown-style running this involves combining together all the trigger files adding in the injection and segment tables, performing injection finding and some clustering of coincident triggers.

The return from this module is a list of pycbc workflow File
objects corresponding to the output files to be directly used in the post-processing.

======
Usage
======

Using this module requires a number of things

* A configuration file (or files) containing the information needed to tell this module how to do the postprocessing preparation.
* An initialized instance of the pycbc Workflow class, containing the ConfigParser.
* A FileList returned by the coincidence module containing the triggers.
* Other FileLists that may be needed by different methods, as described below.

This module is then called according to

.. autofunction:: pycbc.workflow.setup_postprocessing_preparation
   :noindex:

-------------------------
Configuration file setup
-------------------------

Here we describe the options given in the configuration file used in the
workflow that will be needed in this section

$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
[workflow-postprocprep] section
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

The configuration file must have a [workflow-postprocprep] section, which is used to provide instructions to the workflow module on how to set up the postprocessing preparation stage. The first option to choose and provide is

* postprocprep-method = VALUE

The choices here and their description are as described below

* PIPEDOWN_WORKFLOW - This will prepare the output files for pipedown-style postprocessing. This involves combining all triggers, injection sets, segments into a single file, performing injection finding and performing clustering.

.. autofunction:: pycbc.workflow.setup_postprocprep_pipedown_workflow
   :noindex:

With the PIPEDOWN_WORKFLOW method the following options apply

* postprocprep-combiner1-exe=NAME
* postprocprep-combiner2-exe=NAME
* postprocprep-cluster-exe=NAME
* postprocprep-injfind-exe=NAME

these specify which executables will be used for each step. These are described more fully below.

* GSTLAL_POSTPROCPREP - This will perform the gstlal-style post-processing. This involves computing likelihood files, from there FAP and then some plotting routines.

.. autofunction:: pycbc.workflow.setup_postprocprep_gstlal_workflow
   :noindex:


With the GSTLAL_POSTPROCPREP method the following options apply

* postprocprep-runsqlite-exe=NAME
* postprocprep-ligolwsqlite-exe=NAME
* postprocprep-inspinjfind-exe=NAME
* postprocprep-sqltoxml-exe=NAME
* postprocprep-picklehor-exe=NAME
* postprocprep-combllhood-exe=NAME
* postprocprep-genranking-exe=NAME
* postprocprep-compllhood-exe=NAME
* postprocprep-marglikelihood-exe=NAME
* postprocprep-fargstlal-exe=NAME
* postprocprep-plotsummary-exe=NAME
* postprocprep-plotsensitivity-exe=NAME
* postprocprep-plotbackground-exe=NAME
* postprocprep-summarypage-exe=NAME

these specify which executables will be used for each step. These are described 
more fully below.

$$$$$$$$$$$$$$$
[executables]
$$$$$$$$$$$$$$$

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
For the PIPEDOWN_WORKFLOW method
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Executables required by this module are provided in the [executables] section. Any executable names specified in the [workflow-postprocprep] section must appear here. For instance if the [workflow-postprocprep] section reads

* postprocprep-combiner1-exe=pycbcsqlite
* postprocprep-combiner2-exe=pycbcsqlite
* postprocprep-cluster-exe=clustercoincs
* postprocprep-injfind-exe=databaseinjfind

you would need to have

* pycbcsqlite = ${which:pycbc_sqlite_simplify}
* clustercoincs = ${which:ligolw_cbc_cluster_coincs}
* databaseinjfind = ${which:ligolw_dbinjfind}

in the [executables] section.

Sections, in this case [pycbcsqlite], [clustercoincs] and [databaseinjfind], will be used to specify the constant command line options that are sent to all jobs with the corresponding exe name. How to set up the [{exe_name}] section, and which executables are currently supported is discussed below.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
For the GSTLAL_POSTPROCPREP method
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Executables required by this module are provided in the [executables] section. Any executable names specified in the [workflow-postprocprep] section must appear here. For instance if the [workflow-postprocprep] section reads

* ADD

you would need to have

* ADD

in the [executables] section.

Sections will be used to specify the constant command line options that are sent to all jobs with the corresponding exe name. How to set up the [{exe_name}] section, and which executables are currently supported is discussed below.


----------------------------------------------------------------------------
Supported post-processing preparation codes and instructions for using them
----------------------------------------------------------------------------

$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
For the PIPEDOWN_WORKFLOW method
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

The following coincidence codes are currently supported in the workflow module

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
postprocprep-combiner1-exe
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* pycbc_sqlite_simplify

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
postprocprep-combiner2-exe
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* pycbc_sqlite_simplify

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
postprocprep-cluster-exe
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ligolw_cbc_cluster_coincs

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
postprocprep-injfind-exe
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ligolw_dbinjfind

$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
For the GSTLAL_POSTPROCPREP method
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
postprocprep-runsqlite-exe
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ADD

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
postprocprep-ligolwsqlite-exe
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ADD


^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
postprocprep-inspinjfind-exe
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ADD


^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
postprocprep-sqltoxml-exe
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ADD


^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
postprocprep-picklehor-exe
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ADD

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
postprocprep-combllhood-exe
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ADD

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
postprocprep-genranking-exe
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ADD

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
postprocprep-compllhood-exe
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ADD

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
postprocprep-marglikelihood-exe
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ADD

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
postprocprep-fargstlal-exe
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ADD

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
postprocprep-plotsummary-exe
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ADD

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
postprocprep-plotsensitivity-exe
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ADD

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
postprocprep-plotbackground-exe
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ADD

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
postprocprep-summarypage-exe
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ADD


---------------------------
Instructions for each code
---------------------------

Adding a new executable is not too hard, please ask a developer for some pointers on how to do this if you want to add a new code. This may require a new method if the new code uses a different method of generaring coincidences or time sliding.

The instructions for currently supported codes is given below

$$$$$$$$$$$$$$$$$$$$$$$
pycbc_sqlite_simplify
$$$$$$$$$$$$$$$$$$$$$$$

This code is responsible for combining together multiple xml or sql files into a single database. It also includes the code in dbsimplify for removing repeated entries to avoid the resulting database from getting too huge.

.. command-output:: pycbc_sqlite_simplify --help

Of these options the workflow module will automatically add the following, which are unique for each job. **DO NOT ADD THESE OPTIONS IN THE CONFIGURATION FILE**.

* --injection-file
* --simulation-tag
* --output-file
* INPUT_FILE(s) (positional argument)

$$$$$$$$$$$$$$$$$$$$$$$$$$
ligolw_cbc_cluster_coincs
$$$$$$$$$$$$$$$$$$$$$$$$$$

This code is used to perform a clustering stage on a set of coincident triggers.

.. command-output:: ligolw_cbc_cluster_coincs --help

Of these options the workflow module will automatically add the following, which are unique for each job. **DO NOT ADD THESE OPTIONS IN THE CONFIGURATION FILE**.

* --input
* --output

$$$$$$$$$$$$$$$$$$$$$$$$$$
ligolw_dbinjfind
$$$$$$$$$$$$$$$$$$$$$$$$$$

This code is used to perform "injection finding" - it associates injection entries in the sim_inspiral table, with coincident triggers in the coinc_inspiral table.

.. command-output:: ligolw_dbinjfind --help

Of these options the workflow module will automatically add the following, which are unique for each job. **DO NOT ADD THESE OPTIONS IN THE CONFIGURATION FILE**.

* --input
* --output
