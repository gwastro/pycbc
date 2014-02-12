.. _ahopetimeslidesmod:

########################################
The ahope time-slide generation module
########################################

==============
Introduction
==============

The time-slide generation module of ahope is responsible for generating time
slide input files, which are used to tell analysis codes which time slides they
need to perform when estimating the background using the time-slides method.
An example of where this is needed is when using ligolw_thinca a time slide xml file is required as input, even if no time slides are being performed (in that case the time slide file corresponds to no-lags only).

Currently the only option in this module is to have the time slide files generated at runtime. This is because the interface to ligolw_tisi is not compatible with ConfigParser. This introduces a number of issues/hacks in this module which need resolving (as mentioned on the to-do list).

The return of the time-slide generation module is an AhopeFileList of the time-slide files generated/supplied to this module. It is possible to call this module multiple times by using the tags key-word argument, but we have not foreseen
a case where this would be needed. One call should generate all needed time-slide input files.

=======
Usage
=======

Using this module requires a number of things

* A configuration file (or files) containing the information needed to tell this module how to generate (or gather) the template banks (described below).
* An initialized instance of the ahope workflow class, containing the ConfigParser.
* A list of science segments, which are only used to identify start and end times for file naming and because AhopeFiles require associated valid segment times.

This module is then called according to

.. autofunction:: pycbc.ahope.setup_timeslides_workflow
   :noindex:

Configuration file setup
-------------------------

Here we describe the options given in the configuration file used in the ahope
workflow that will be needed in this section

$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
[ahope-timeslides] section
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

The configuration file must have an [ahope-timeslides] section, which is used to
tell the workflow how to construct (or gather) the time slide files. The first option to choose and provide is

* timeslides-method = VALUE


The choices here and their description are as described below

* AT_RUNTIME - The time slide files are generated within the call to this module.
* Action items to add other methods (IN_WORKFLOW and PREGENERATED) can be found in the to-do list, but some of this requires the ligolw_tisi issues to be addressed first.

With only one current option, this module has no sub-modules. This may change in the future.

Currently no additional options are required in the [ahope-injections] section.

$$$$$$$$$$$$$$$
[executables]
$$$$$$$$$$$$$$$

In this section you need to supply the executable that will be used to generate the time slide files. This is done in the [executables] section by adding something like::

timeslides = /path/to/ligolw_tisi

the option, in this case 'timeslides', will be used to specify the constant command line options that are sent to all ligolw_tisi jobs. The tag 'timeslides' can be changed, it is currently supplied as a key-word argument when calling the function. 'ligolw_tisi' is the default.

**PLEASE NOTE THE ISSUES WHEN USING ligolw_tisi AS REPORTED BELOW**

----------------------------------------------------------------
Supported injection executables and instructions for using them
----------------------------------------------------------------

The following time slide generation executables are currently supported in ahope

* ligolw_tisi

Adding a new executable is not too hard, please ask a developer for some pointers on how to do this if you want to add a new code.

$$$$$$$$$$$$$$$$$
ligolw_tisi
$$$$$$$$$$$$$$$$$

ligolw_tisi is the pylal ligolw python code that was responsible for making time slide tables in ihope. The help message for ligolw_tisi follow

.. command-output:: ligolw_tisi --help

At the moment these options cannot be supplied in the ConfigParser because of the repeated use of the --instrument option. Instead the following options need to be given in a [tisi] section

If

inspiral-num-slides = VALUE

is in the section then add that to the command line call. Otherwise look for

X1-slide-start = INT
X1-slide-end = INT
X1-slide-step = INT

where X1 loops over all active ifos. This is translated to --instrument X1=slide_start:slide_end:slide_step

Then if

remove-zero-lag

is present this will be added to the command line call.

The output file argument is automatically added.

**THIS NEEDS FIXING AND MAKING COMPATIBLE WITH OTHER AHOPE EXE CALLS**

An example of a ligolw_tisi call, for time slides, is given below

.. code-block:: bash
   ligolw_tisi -v --inspiral-num-slides 100 --remove-zero-lag /home/spxiwh/lscsoft_git/src/pycbc/examples/ahope/weekly_ahope/961585543-961671943/time_slide_files/H1L1-TIMESLIDES_SLIDES-961585543-86400.xml.gz

and for zero-lag only

.. code-block:: bash
   ligolw_tisi -v -i H1=0:0:0 -i L1=0:0:0 /home/spxiwh/lscsoft_git/src/pycbc/examples/ahope/weekly_ahope/961585543-961671943/time_slide_files/H1L1-TIMESLIDES_ZEROLAG-961585543-86400.xml.gz

============================================
:mod:`pycbc.ahope.timeslides_utils` Module
============================================

This is complete documentation of this module's code

.. automodule:: pycbc.ahope.timeslides_utils
    :noindex:
    :members:
    :undoc-members:
    :show-inheritance:




