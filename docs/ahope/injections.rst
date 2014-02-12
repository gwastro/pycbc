.. _ahopeinjectionsmod:

######################################
The ahope injection generation module
######################################

==============
Introduction
==============

The injection generation module of ahope is responsible for generating or
gathering the injection files that will be used to assess the performance of a
pipeline to detect simulated signals. It is foreseen that this could be run
either by generating injection files at run time, by generating injection files
within the workflow or by using pregenerated injection files.

Currently the only option is to add jobs that generate injection files within
the workflow.

The return of the injection generation module is an AhopeFileList of the injection files generated/supplied to this module. It is possible to call this module
multiple times by using the tags key-word argument, but we have not foreseen
a case where this would be needed.

======
Usage
======

Using this module requires a number of things

* A configuration file (or files) containing the information needed to tell this module how to generate (or gather) the template banks (described below).
* An initialized instance of the ahope workflow class, containing the ConfigParser.
* A list of science segments, which are only used to identify start and end times of the injection files.

The module is then called according to

.. autofunction:: pycbc.ahope.setup_injection_workflow
   :noindex:

-------------------------
Configuration file setup
-------------------------

Here we describe the options given in the configuration file used in the ahope
workflow that will be needed in this section

$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
[ahope-injections] section
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

The configuration file must have an [ahope-injections] section, which is used to
tell the workflow how to construct (or gather) the injection files. The first option to choose and provide is

* injections-method = VALUE

The choices here and their description are as described below

* IN_WORKFLOW - The injection file generation will be added as jobs in the workflow and will be generated after submission of the workflow.
* Action items to add other methods (AT_RUNTIME and PREGNERATED) can be found in the to-do list.

With only one current option, this module has no sub-modules. This may change in the future.

Currently no additional options are required in the [ahope-injections] section.

$$$$$$$$$$$$$$
[executables]
$$$$$$$$$$$$$$

In this section you need to supply the executable that will be used to generate the injection files. This is done in the [executables] section by adding something like::

injections = /path/to/lalapps_inspinj

the option, in this case 'injections', will be used to specify the constant command line options that are sent to all lalapps_inspinj jobs. The tag 'injections' can be changed, it is currently supplied as a key-word argument when calling the function. 'injections' is the default, and must be used if you want to use pipedown.

----------------------------------------------------------------
Supported injection executables and instructions for using them
----------------------------------------------------------------

The following injection generation executables are currently supported in ahope

* lalapps_inspinj

Adding a new executable is not too hard, please ask a developer for some pointers on how to do this if you want to add a new code.

$$$$$$$$$$$$$$$$$
lalapps_inspinj
$$$$$$$$$$$$$$$$$

lalapps_inspinj is the legacy C-code that has been used to generate injection files for gravitational-wave data analysis since the dawn of time. The help message for lalapps_inspinj follows

.. command-output:: lalapps_inspinj --help

Of these options ahope will automatically add the following. **DO NOT ADD THESE OPTIONS IN THE CONFIGURATION FILE**.

* --gps-start-time
* --gps-end-time
* --output

All other options must be provided in the configuration file. Here is an example of a lalapps_inspinj call.

.. code-block:: bash
   lalapps_inspinj --gps-end-time 961671943 --i-distr uniform --max-mass1 3.1 --max-mass2 3.1 --m-distr componentMass --disable-spin  --min-mtotal 2.0 --output /home/spxiwh/lscsoft_git/src/pycbc/examples/ahope/weekly_ahope/961585543-961671943/inj_files/HL-INJECTIONS_BNSLININJ_2134-961585543-86400.xml --max-mtotal 6.2 --waveform TaylorT4threePointFivePN --time-interval 300 --time-step 837.155 --min-mass2 1.0 --f-lower 30 --l-distr random --min-mass1 1.0 --min-distance 1000 --gps-start-time 961585543 --d-distr uniform --max-distance 60000

==========================================
:mod:`pycbc.ahope.injection_utils` Module
==========================================

This is complete documentation of this module's code

.. automodule:: pycbc.ahope.injection_utils
    :noindex:
    :members:
    :undoc-members:
    :show-inheritance:

