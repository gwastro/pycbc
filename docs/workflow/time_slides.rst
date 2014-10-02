.. _workflowtimeslidesmod:

##########################################
The workflow time-slide generation module
##########################################

==============
Introduction
==============

The time-slide generation sub-module of pycbc's workflow module
is responsible for generating time
slide input files, which are used to tell analysis codes which time slides they
need to perform when estimating the background using the time-slides method.
An example of where this is needed is when using ligolw_thinca a time slide xml file is required as input, even if no time slides are being performed (in that case the time slide file corresponds to no-lags only). The module can either
generate time slide files at runtime, within the workflow or take pregenerated
input files.

The return of the time-slide generation module is an pycbc FileList of the time-slide files generated/supplied to this module. It is possible to call this module multiple times by using the tags key-word argument, but we have not foreseen
a case where this would be needed. One call should generate all needed time-slide input files.

=======
Usage
=======

Using this module requires a number of things

* A configuration file (or files) containing the information needed to tell this module how to generate (or gather) the template banks (described below).
* An initialized instance of the pycbc Workflow class, containing the ConfigParser.

This module is then called according to

.. autofunction:: pycbc.workflow.setup_timeslides_workflow
   :noindex:

Configuration file setup
-------------------------

Here we describe the options given in the configuration file used in the
workflow that will be needed in this section

$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
[workflow-timeslides] section
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

The configuration file must have a [workflow-timeslides] section, which is used to
tell the workflow how to construct (or gather) the time slide files. The first option to choose and provide is

* timeslides-method = VALUE

The choices here and their description are as described below

* IN_WORKFLOW - The time slide file generation will be added as jobs in the workflow and will be generated after submission of the workflow.
* AT_RUNTIME - The time slide jobs will be run directly within this module.
* PREGENERATED - The time slide files will be supplied as pregenerated files.

For all cases (even PREGENERATED) the following option is needed:

* timeslides-exe = VALUE - This tag is used to locate the timeslides exe and the options needed for that exe. It will also identify which time slides are being done, see below for details.

When using PREGENERATED the additional option is needed:

* timeslides-pregenerated-file = PATH - Specifies the location of the pregenerated time slides file.

$$$$$$$$$$$$$$$
[executables]
$$$$$$$$$$$$$$$

**PLEASE READ THIS EVEN IF USING method=PREGENERATED**

In this section, if not using PREGENERATED you need to supply the executable that will be used to generate the time slide files. This is done in the [executables] section by adding something like:

timeslides = /path/to/pycbc_timeslides

the option, in this case 'timeslides', will be used to specify the constant command line options that are sent to all pycbc_timeslides jobs. This is done in the [timeslides] section and the specific options for each executable is discussed below. The tag 'timeslides' is taken from the workflow-timeslides section of the configuration file, as described above.

**IMPORTANT** the time slides tag name, (e.g. 'timeslides'), is used to identify which time slide files are to be produced. The module will search for modules called [timeslides-XXXX] and generate one time slide file for every sub-section using XXXX as a tag.

When using the IN_WORKFLOW or AT_RUNTIME options, value pairs in [timeslides-XXXX] are also used to specify options specific to that
particular time slide job. These options will be directly supplied to the
executable being used to generate the file.

**PLEASE NOTE:** When using PREGENERATED the time slide names are still identified by looking in the "tisi-XXXX" sections (a path in [executables] is not needed). These sections will be empty in this case. Also note that as XXXX is used as a tag you can use [workflow-timeslides-XXXX] section to supply unique pregenerated time slide files. You can also use this to supply a mix of pregenerated, and not-pregenerated time slide files.

----------------------------------------------------------------
Supported time slideexecutables and instructions for using them
----------------------------------------------------------------

The following time slide generation executables are currently supported

* pycbc_timeslides

Adding a new executable is not too hard, please ask a developer for some pointers on how to do this if you want to add a new code.

$$$$$$$$$$$$$$$$$
pycbc_timeslides
$$$$$$$$$$$$$$$$$

pycbc_timeslides is the pycbc python code that is almost identical to the code used for making time slide tables in ihope. The help message for pycbc_timeslides follows

.. command-output:: pycbc_timeslides --help

An example of a pycbc_timeslides call, for time slides, is given below

.. code-block:: bash

   /home/spxiwh/lscsoft_git/executables_master/bin/pycbc_timeslides --output-files /home/spxiwh/lscsoft_git/src/pycbc/examples/ahope/weekly_ahope/961585543-961671944/time_slide_files/H1L1V1-TISI_SLIDES-961585543-86401.xml.gz --remove-zero-lag --inspiral-num-slides 100:H1=0,L1=5,V1=10

and for zero-lag only

.. code-block:: bash

   /home/spxiwh/lscsoft_git/executables_master/bin/pycbc_timeslides --output-files /home/spxiwh/lscsoft_git/src/pycbc/examples/ahope/weekly_ahope/961585543-961671944/time_slide_files/H1L1V1-TISI_ZEROLAG-961585543-86401.xml.gz --tisi-slides H1=0:0:0 L1=0:0:0 V1=0:0:0

The following options are added by the workflow module and **must not** be given in the configuration file

* --output-files

============================================
:mod:`pycbc.workflow.timeslides` Module
============================================

This is complete documentation of this module's code

.. automodule:: pycbc.workflow.timeslides
    :noindex:
    :members:
    :undoc-members:
    :show-inheritance:

