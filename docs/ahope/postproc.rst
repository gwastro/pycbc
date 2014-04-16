.. _ahopepostprocmod:

############################################
The ahope postprocessing module
############################################

=============
Introduction
=============

The postprocessing module of ahope is used to assess the significance of all triggers identified at earlier stages of the pipeline and to compute any rate statements that can be made from the analysis results.

The return from this module is a list of AhopeOutFile objects corresponding to the files containing the trigger significances and rate statements.

======
Usage
======

Using this module requires a number of things

* A configuration file (or files) containing the information needed to tell this module how to do the postprocessing preparation.
* An initialized instance of the ahope workflow class, containing the ConfigParser.
* An AhopeFileList containing the analysis results.

This module is then called according to

.. autofunction:: pycbc.ahope.setup_postprocessing
   :noindex:

-------------------------
Configuration file setup
-------------------------

Here we describe the options given in the configuration file used in the ahope
workflow that will be needed in this section

$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
[ahope-postproc] section
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

The configuration file must have an [ahope-postproc] section, which is used to provide instructions to ahope on how to set up the post-processinge stage. The first option to choose and provide is

* postproc-method = VALUE

The choices here and their description are as described below

* PIPEDOWN_AHOPE - This will perform pipedown style post-processing, using the same methods as was used by pipedown to assess trigger significance and compute rate statements (currently rate statements are not calculated)

Currently only one option, but others can be added. The subfunctions used are described here

.. autofunction:: pycbc.ahope.setup_postproc_pipedown_ahope
   :noindex:

With the PIPEDOWN_AHOPE method the following options apply

* postproc-computedurations-exe=NAME
* postproc-cfar-exe=NAME

these specify which executables will be used for each step. These are described more fully below.

$$$$$$$$$$$$$$$
[executables]
$$$$$$$$$$$$$$$

Executables required by this module are provided in the [executables] section. Any executable names specified in the [ahope-postproc] section must appear here. For instance if the [ahope-postproc] section reads

* postproc-computedurations-exe=computedurs
* postproc-cfar-exe=pycbccfar

you would need to have

* computedurs = ${which:pycbc_compute_durations}
* pycbccfar = ${which:pycbc_calculate_far}

in the [executables] section.

Sections, in this case [computedurs] and [pycbccfar], will be used to specify the constant command line options that are sent to all jobs with the corresponding exe name. How to set up the [{exe_name}] section, and which executables are currently supported is discussed below.

----------------------------------------------------------------------------
Supported post-processing codes and instructions for using them
----------------------------------------------------------------------------

The following coincidence codes are currently supported in ahope

$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
postproc-computedurations-exe
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

* pycbc_compute_durations

$$$$$$$$$$$$$$$$$$$$$$$$$$$$
postproc-cfar-exe
$$$$$$$$$$$$$$$$$$$$$$$$$$$$

* pycbc_calculate_far

---------------------------
Instructions for each code
---------------------------

Adding a new executable is not too hard, please ask a developer for some pointers on how to do this if you want to add a new code. This may require a new method if the new code uses a different method of generaring coincidences or time sliding.

The instructions for currently supported codes is given below

$$$$$$$$$$$$$$$$$$$$$$$$
pycbc_compute_durations
$$$$$$$$$$$$$$$$$$$$$$$$

This code is responsible for determining how how time was analysed, in coincidence, and after vetoes have been applied. This is then needed to assess the false alarm rates.

.. command-output:: pycbc_compute_durations --help

Of these options ahope will automatically add the following, which are unique for each job. **DO NOT ADD THESE OPTIONS IN THE CONFIGURATION FILE**.

* --input
* --output

$$$$$$$$$$$$$$$$$$$$$$$$$$
pycbc_calculate_far
$$$$$$$$$$$$$$$$$$$$$$$$$$

This code will calculate the false alarm rate associated with each trigger. It will calculate an "uncombined" false alarm rate, for each triggers mass bin and coincidence type, and a "combined" false alarm rate, after trials factor from the number of mass bins * coincident types has been accounted for.

.. command-output:: pycbc_calculate_far --help

Of these options ahope will automatically add the following, which are unique for each job. **DO NOT ADD THESE OPTIONS IN THE CONFIGURATION FILE**.

* --input
* --output
