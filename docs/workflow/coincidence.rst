.. _workflowcoincmod:

##################################
The workflow coincidence module
##################################

=============
Introduction
=============

The coincidence module of the workflow module is responsible for
performing coincidence
between single detector triggers to determine if any of these triggers are
coincident and potential gravitational wave candidates. This can be run on
foreground triggers, but can also be run with time shifts included between
the sets of detectors to assess the background rate of coincident triggers.

This module should be capable of performing different "flavours" of coincidence.
The standard method is ethinca coincidence over masses and time between sets
of triggers, however other methods include "exact mass" coincidence, and
possibly coincidence algorithms including the effects of spin.

The return from the coincidence module is a list of pycbc.workflow.File
objects corresponding to the files holding the lists of coincidence triggers.

======
Usage
======

Using this module requires a number of things

* A configuration file (or files) containing the information needed to tell this module how to determine coincident GW triggers.
* An initialized instance of the Workflow class, containing the ConfigParser.
* A FileList returned by the matchedfilter module containing the single detector triggers to perform coincidence over.
* A FileList of the data quality veto files if they are going to be included during the coincidence stage.
* A FileList of the time slides that will be performed. If you are doing no time slides a "no time slides" time slide file still needs to be provided.

This module is then called according to

.. autofunction:: pycbc.workflow.setup_coincidence_workflow
   :noindex:

-------------------------
Configuration file setup
-------------------------

Here we describe the options given in the configuration file used in the
workflow that will be needed in this section

$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
[workflow-coincidence] section
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

The configuration file must have a [workflow-coincidence] section, which is used to provide instructions to the workflow module on how to set up the coincidence stage. The first option to choose and provide is

* coincidence-method = VALUE

The choices here and their description are as described below

* WORKFLOW_DISCRETE_SLIDES - This will perform coincidence at the discrete time shifts provided in the time slide file(s). Multiple time slide tags can be provided to run the code on multiple time slide files. If you only want to obtain foreground coincidence triggers simply provide a time slide file with only a single, no time shift, entry. This is the same method as used in ihope runs.

Currently only one option, but others can be added. The subfunctions used are described here

.. autofunction:: pycbc.workflow.setup_coincidence_workflow_ligolw_thinca
   :noindex:

* maximum-extent = EXTENT_IN_SECONDS (default=3600). The maximum amount of zerolag or time-shifted analysis time a single coincidence job will process.  

If running with ethinca coincidence, the computational cost of each coincidence job will scale with the time extent, thus if jobs are evicted before finishing it may help to rerun with a smaller value, the number of jobs will scale up inversely.

If exact-match coincidence is used, then the coincidence jobs can be parallelized to reduce the memory footprint. This requires that the inspiral jobs also use a split bank. To parallelize the coincidence, add under [workflow-coincidence]

* coincidence-exact-match-parallelize =

To reduce output file size, it can be useful to put clustering jobs directly
after the coincidence jobs. To do this add

* coincidence-post-cluster =

If it is desired to use the gstlal-likelihood-based post-processing, which requires exact-match, then the option

* coincidence-write-likelihood =

Can be added, which instructs the code to write gstlal-style likelihood files,
and will convert output to gstlal format. At the moment this means that the
chisq column is replaced with chisq / dof (hopefully an agreement can be made
so that the two write the *same* format files).

$$$$$$$$$$$$$$$
[executables]
$$$$$$$$$$$$$$$

Executables required by this module are provided in the [executables] section. Currently the executables must called llwadd and thinca, but this will be changed soon to be more flexible

llwadd = ${which:ligolw_add}
thinca = ${which:ligolw_sstinca}

A section, in this case [llwadd] and [thinca], will be used to specify the constant command line options that are sent to all llwadd/thinca jobs. How to set up the [{exe_name}] section, and which executables are currently supported is discussed below.

-----------------------------------------------------------------------
Supported coincidence codes and instructions for using them
-----------------------------------------------------------------------

The following coincidence codes are currently supported

* ligolw_add followed by ligolw_sstinca

Adding a new executable is not too hard, please ask a developer for some pointers on how to do this if you want to add a new code. This may require a new method if the new code uses a different method of generaring coincidences or time sliding.

$$$$$$$$$$$$$$$
ligolw_add
$$$$$$$$$$$$$$$

This code simply adds together the input files in preparation for ligolw_sstinca. We plan to have a front end combining this functionality soon.

.. command-output:: ligolw_add --help

Of these options the workflow module will automatically add the following, which are unique for each job. **DO NOT ADD THESE OPTIONS IN THE CONFIGURATION FILE**.

* --input-cache (this will not be used and should not be used)
* --output
* Positional arguments supplying the input files

$$$$$$$$$$$$$$$$$$$
ligolw_cbc_sstinca
$$$$$$$$$$$$$$$$$$$

This code is responsible for identifying coincidences and applying and vetoes that are present and supplied.

.. command-output:: ligolw_cbc_sstinca --help

Of these options the workflow will automatically add the following, which are unique for each job. **DO NOT ADD THESE OPTIONS IN THE CONFIGURATION FILE**.

* Positional arguments supplying the input file
* --vetoes-name
* --coinc-end-time-segment
