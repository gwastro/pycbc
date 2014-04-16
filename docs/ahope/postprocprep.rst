.. _ahopepostprocprepmod:

############################################
The ahope postprocessing preperation module
############################################

=============
Introduction
=============

The postprocessing preparation module is used to prepare the output of the coincidence stage for the postprocessing (ie. calculation of trigger significance and of rate statements). For the case of pipedown-style running this involves combining together all the trigger files adding in the injection and segment tables, performing injection finding and some clustering of coincident triggers.

The return from this module is a list of AhopeOutFile
objects corresponding to the output files to be directly used in the post-processing.

======
Usage
======

Using this module requires a number of things

* A configuration file (or files) containing the information needed to tell this module how to do the postprocessing preparation.
* An initialized instance of the ahope workflow class, containing the ConfigParser.
* An AhopeFileList returned by the coincidence module containing the triggers.
* Other AhopeFileLists that may be needed by different methods, as described below.

This module is then called according to

.. autofunction:: pycbc.ahope.setup_postprocessing_preperation
   :noindex:

-------------------------
Configuration file setup
-------------------------

Here we describe the options given in the configuration file used in the ahope
workflow that will be needed in this section

$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
[ahope-postprocprep] section
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

The configuration file must have an [ahope-postprocprep] section, which is used to provide instructions to ahope on how to set up the postprocessing preparation stage. The first option to choose and provide is

* postprocprep-method = VALUE

The choices here and their description are as described below

* PIPEDOWN_AHOPE - This will prepare the output files for pipedown-style postprocessing. This involves combining all triggers, injection sets, segments into a single file, performing injection finding and performing clustering.

Currently only one option, but others can be added. The subfunctions used are described here

.. autofunction:: pycbc.ahope.setup_postprocprep_pipedown_ahope
   :noindex:

With the PIPEDOWN_AHOPE method the following options apply

* postprocprep-combiner1-exe=NAME
* postprocprep-combiner2-exe=NAME
* postprocprep-cluster-exe=NAME
* postprocprep-injfind-exe=NAME

these specify which executables will be used for each step. These are described more fully below.

$$$$$$$$$$$$$$$
[executables]
$$$$$$$$$$$$$$$

Executables required by this module are provided in the [executables] section. Any executable names specified in the [ahope-postprocprep] section must appear here. For instance if the [ahope-postprocprep] section reads

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

----------------------------------------------------------------------------
Supported post-processing preparation codes and instructions for using them
----------------------------------------------------------------------------

The following coincidence codes are currently supported in ahope

$$$$$$$$$$$$$$$$$$$$$$$$$$$
postprocprep-combiner1-exe
$$$$$$$$$$$$$$$$$$$$$$$$$$$

* pycbc_sqlite_simplify

$$$$$$$$$$$$$$$$$$$$$$$$$$$$
postprocprep-combiner2-exe
$$$$$$$$$$$$$$$$$$$$$$$$$$$$

* pycbc_sqlite_simplify

$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
postprocprep-cluster-exe
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

* ligolw_cbc_cluster_coincs

$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
postprocprep-injfind-exe
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

* ligolw_dbinjfind

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

Of these options ahope will automatically add the following, which are unique for each job. **DO NOT ADD THESE OPTIONS IN THE CONFIGURATION FILE**.

* --injection-file
* --simulation-tag
* --output-file
* INPUT_FILE(s) (positional argument)

$$$$$$$$$$$$$$$$$$$$$$$$$$
ligolw_cbc_cluster_coincs
$$$$$$$$$$$$$$$$$$$$$$$$$$

This code is used to perform a clustering stage on a set of coincident triggers.

.. command-output:: ligolw_cbc_cluster_coincs --help

Of these options ahope will automatically add the following, which are unique for each job. **DO NOT ADD THESE OPTIONS IN THE CONFIGURATION FILE**.

* --input
* --output

$$$$$$$$$$$$$$$$$$$$$$$$$$
ligolw_dbinjfind
$$$$$$$$$$$$$$$$$$$$$$$$$$

This code is used to perform "injection finding" - it associates injection entries in the sim_inspiral table, with coincident triggers in the coinc_inspiral table.

.. command-output:: ligolw_dbinjfind --help

Of these options ahope will automatically add the following, which are unique for each job. **DO NOT ADD THESE OPTIONS IN THE CONFIGURATION FILE**.

* --input
* --output
