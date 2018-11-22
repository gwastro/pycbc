.. _workflowconfigparsermod:

#########################################################################
Pycbc's workflow module configuration file(s) and command line interface
#########################################################################

============
Introduction
============

The workflow module at its core is designed to be flexible and allow the user
to do what they
want to create the pipeline that they want to run. One of the ways to allow
this is by having a, sometimes large, configuration file that serves two
purposes

* Tell the workflow planner, how to run the various stages specified in
  the top-level workflow script.
* Specify, as completely as possible, all command line options that will be
  sent to every executable that is run in the pipeline. Tags are used to
  identify options sent a subset of jobs, as described more fully later.

The idea is that the only input that a user *needs* is the configuration file.
However, it may often be useful for certain options, such as user-specific
locations and analysis start/end times, to be supplied on the command line.
To allow this we allow a method by which configuration file options can be
supplied, or overridden, on the command line.

Ihope used similar .ini files in every analysis. However, it was noted that these files grew huge and it becomes difficult for a novice to understand which options can be safely changed and which ones to leave well alone. It is also difficult so see which options are going to which job, inspiral.c for example looks for options in > 10 places and it isn't clear where those places are.

To attempt to solve this the workflow module has a number of features

* Multiple configuration files: You can now supply multiple configuration files to, for e.g. identify a file containing only injection generation parameters, which a user may want to change often. It is even possible to have sections split across files, so one could have a configuration file of key options, ones that might be changed, and another file of "leave alone" options.
* Direct command line options: In the workflow module command line options are not drawn from obscure sections, they correspond one-to-one with the executables. Options in the [inspiral] section will be sent to the inspiral executable and *only* to the inspiral executable.
* Combined sections: To avoid the issue of specifiying common options repeatedly we have allowed the ability of combined sections. So if you have two executables with a large set of shared options you can specify a [exe1&exe2] section to provide the shared options and [exe1] and [exe2] sections to supply the individual options. One can also use the [sharedoptions-NAME] sections to acheive the same thing.
* Interpolation: As in configparser 3.0+ we have the ability to specify an option in one place and use an interpolation string to also provide it in other places, this is described below.
* Tags/subsections: In some cases options may only need to be sent to certain jobs, or you may want to call individual modules multiple times and do different things. To accomodate this the workflow module includes a tagging (or subsections) system to provide options to only a subset of jobs, or to a specific call to a module. For example, options in [inspiral] are sent to all inspiral jobs, options in [inspiral-h1] would be sent to inspiral jobs running only on h1 data.
* Executable expanding: The workflow module includes macros to enable the user to more easily specify executable paths. For example $(which:exe1} will be expanded to the location of exe1 in the users path automatically.

most of these features will be applied directly after reading in the configuration file. The workflow module will then dump the parser configuration back to disk so the user/reviewer can more easily see what the analysis is actually doing.

In this page we describe the layout of the workflow module .ini configuration file and what the various sections mean, how they are used, and how an ini file should be set out. 

**NOTE**: A number of features that have been put in here, are available in the python 3.X version of ConfigParser. In addition this version also has a duplicate option check. In python 2.X if I do::

    [inspiral]
    detect-gravitational-waves = True
    LOTS OF GARBAGE
    detect-gravitational-waves = False

it will set the value to False, and proceed happily. THERE IS NO WAY TO CATCH THIS! There is a python 2.X backport of this new version, it is available in pypi, but not in macports. It would be good to pick up this new version and have some of these features available natively.

======================================================================
Supplying the config file on the command line and overriding options
======================================================================

The workflow module only uses two command line options, one to specify the
configuration files
and one to specify and overriding options. First the config files:

* --config-files FILE1 [FILE2 FILE3 ....]

where FILEX corresponds to the configuration files. Second the overriding options:

* --config-overrides section1:option1:value1 [section2:option2:value2 ...]

These specify options that should be *added* to the config files, or if already present *overwritten*. The section, option and value refer to the section option and value to be added. If the section doesn't already exist in the configuration file it will be added. In some cases you will want to supply an option without a value. This can be done with either

section:option:

or

section:option

--------------
Example
--------------

Here is an example of running a workflow from the command line::

  python weekly_ahope.py --config-files weekly_ahope.ini pipedown.ini inj.ini --config-overrides workflow:start-time:${GPS_START_TIME} workflow:end-time:${GPS_END_TIME}

Here the analysis start and end times are being overriden with values from the user's environment.

=======================================
Global options - the [workflow] section
=======================================

The [workflow] section and [workflow-XXX] subsections should appear at the top of a configuration file.

The [workflow] section and [workflow-XXX] subsections of the configuration file are used to store options that the workflow module uses to make decisions on what paths to take when deciding how to construct the workflow. Options in here are *not* going to end up supplied to any executable on the command line.

The [workflow] section must contain two entries

* start-time=START
* end-time=END

which are used to tell the workflow that is only to consider times in [START,END) for analysis. These will often be supplied as override options directly on the command line.

Another optional entry in the [workflow] section, that we recommend be used is the:

* file-retention-level = all_files

entry. This can take one of 4 values: "all_files", "all_triggers", "merged_triggers" or "results". These specify how many files produced during the workflow should be stored after the workflow finishes. With "all_files", which is the default value, everything produced in the workflow will be stored. With "results" only the critical result files are stored. "all_triggers" and "merged_triggers" store some subset of the full set of files. Defining whether a file should be stored under each of these levels is the job of the Executable class, which carries a current_retention_level attribute (one of executable.INTERMEDIATE_PRODUCT, executable.ALL_TRIGGERS, executable.MERGED_TRIGGERS or executable.FINAL_RESULT). When building workflows one can set this atrribute when creating executable instances to set under what conditions a file should be stored.

It is okay to store other *important and widely used* values in here. You might often see cases where channel names are given here as these are sent to a number of codes on the command line, and it is easier to refer to them here, at the very top of the .ini file, so that the user can more easily see and change such values.

---------------------------------
[workflow-XXX] subsections
---------------------------------

Each module that you use when setting up your workflow will need an [workflow-XXX] subsection. The name of the subsection and the particular options needed can be found in each module's documentation page.

If you want to call any module more than once you will need to use the workflow module's tagging system. As an example let's say I want to call the template bank module twice, once to set up a pycbc template bank and once to set up a SVD template bank. I could then create [workflow-tmpltbank-pycbc] and [workflow-tmpltbank-svd] sections to provide options that are unique to each tag. I could also use [exename-pycbc] and [exename-svd] sections if the two methods are using the same executable, but need different options. In both cases options in [workflow-tmpltbank] and [exename] would be used for *both* tags. (If the two codes were using different executables then [exename1] and [exename2] sections would suffice.)

An example of where this section might be used is in the template bank stage where one can either run with a pre-generated bank or generate banks within the workflow. This information would be provided in this section.

----------------------
Requirements
----------------------

The [workflow] section in *every* .ini file should contain a link to this page to see what options are needed.

The [workflow-XXX] sections in *every* .ini file should start with a link to that module's documentation to see what options/values are relevant for that section.

----------------------
Example
----------------------

Here is an example of the [workflow] section of a .ini file::

  [workflow]
  ; https://ldas-jobs.ligo.caltech.edu/~cbc/docs/pycbc/workflow/initialization.html
  ; provides details of how to set up a pycbc workflow configuration .ini file
  h1-channel-name = H1:LDAS-STRAIN
  l1-channel-name = L1:LDAS-STRAIN
  ;h2-channel-name = H2:LDAS-STRAIN
  workflow-html-basedir = /home/spxiwh/public_html/workflow/development/weekly_ahope/test

  [workflow-ifos]
  ; This is the list of ifos to analyse
  h1 =
  l1 =

  [workflow-datafind]
  ; See https://ldas-jobs.ligo.caltech.edu/~cbc/docs/pycbc/workflow/datafind.html
  datafind-method = AT_RUNTIME_SINGLE_FRAMES
  datafind-h1-frame-type = H1_LDAS_C02_L2
  datafind-l1-frame-type = L1_LDAS_C02_L2
  ;datafind-h2-frame-type = H2_LDAS_C02_L2
  datafind-check-segment-gaps = update_times
  datafind-check-frames-exist = raise_error
  datafind-check-segment-summary = no_test
  ; Set this to sepcify the datafind server. If this is not set the code will
  ; use the value in ${LIGO_DATAFIND_SERVER}
  ;datafind-ligo-datafind-server = ""

  [workflow-segments]
  ; See https://ldas-jobs.ligo.caltech.edu/~cbc/docs/pycbc/workflow/segments.html
  ; PIPEDOWN demands we use AT_RUNTIME
  segments-method = AT_RUNTIME
  segments-H1-science-name = H1:DMT-SCIENCE:4
  segments-L1-science-name = L1:DMT-SCIENCE:4
  ;segments-V1-science-name = V1:ITF_SCIENCEMODE:6
  segments-database-url = https://segdb.ligo.caltech.edu
  segments-veto-definer-url = https://www.lsc-group.phys.uwm.edu/ligovirgo/cbc/public/segments/S6/H1L1V1-S6_CBC_LOWMASS_B_OFFLINE-937473702-0.xml
  segments-veto-categories = 2,3,4
  segments-minimum-segment-length = 2000
  segments-generate-coincident-segments =

  [workflow-tmpltbank]
  ; See https://ldas-jobs.ligo.caltech.edu/~cbc/docs/pycbc/workflow/template_bank.html
  tmpltbank-method=WORKFLOW_INDEPENDENT_IFOS
  ; Remove the option below to disable linking with matchedfilter_utils
  tmpltbank-link-to-matchedfltr=

  [workflow-injections]
  ; See https://ldas-jobs.ligo.caltech.edu/~cbc/docs/pycbc/workflow/injections.html
  injections-method=IN_WORKFLOW

  [workflow-timeslides]
  ; See https://ldas-jobs.ligo.caltech.edu/~cbc/docs/pycbc/workflow/time_slides.html
  timeslides-method=AT_RUNTIME


=================================================
Executable locations - the [executables] section
=================================================

This section should contain the names of each of the executables that will be used in the workflow and their locations. 

-------------------
executable macros
-------------------

The following macros can be used **only** within this section to automatically fill in full path names

$$$$$$$$$$$$$$$$$$$$$
which(executable)
$$$$$$$$$$$$$$$$$$$$$

In the following example tmpltbank's value will be replaced with the output of which(lalapps_tmpltbank)::

  [executables]
  tmpltbank = ${which:lalapps_tmpltbank}
  inspiral = /full/path/to/lalapps_inspiral

------------------
Requirements
------------------

All executables used in the workflow should be supplied in this section, and *only* in this section.

-------------------
Example
-------------------

Here is an example of the [executables] section of a pycbc workflow .ini file::

  [executables]
  tmpltbank         = /home/cbc/opt/s6b/ab577e4e5dad14e46fce511cffdb04917836ba36/bin/lalapps_tmpltbank
  inspiral          = /home/cbc/opt/s6b/ab577e4e5dad14e46fce511cffdb04917836ba36/bin/lalapps_inspiral
  inspinj           = /home/cbc/opt/s6b/ab577e4e5dad14e46fce511cffdb04917836ba36/bin/lalapps_inspinj
  thinca            = ${which:ligolw_thinca}


===================
Executable options
===================

For each of the executables in the [executables] section, options for that executable should be listed under the section corresponding to that executable. Options in the [tmpltbank] section are sent to lalapps_tmpltbank, options in the [inspiral] section are sent to lalapps_inspiral etc.

It is possible to have more than one [tmpltbank] section, ConfigParser will simply combine them together when reading in. Therefore '''important options''' and '''options that a novice user might want to change''' could be supplied in a first [tmpltbank] section near the top of the .ini file. This section could be commented accordingly. The modules documentation page should also include instructions for each of the supported executables (usually the code's own help message). Options that are not so important and ones that a novice user would not want to change could be placed in a second [tmpltbank] section at the bottom of the ini file, this section would be labelled accordingly and also contain a link to documentation for that executable.

Some options are only sent to a subset of jobs using a given executable. For example those running on H1 data. Options like these will be provided in sections labelled [executable_name-subset_tag]. So for the H1 example the section would be called [tmpltbank-H1]. As well as obeying the rules above these section must clearly state ''which'' jobs will be sent those options. This can also be used when calling a section multiple times with different tags. Nested tags are not supported (ie [tmpltbank-H1-pycbc])

Some options need to be sent to more than one executable, for example the channel names are used by any code that reads the data. Such sections should be given as the combination of executable names separated by the & token. So options sent to tmpltbank '''and''' inspiral would go in a section called [tmpltbank&inspiral]. The code parsing the .ini file will automatically separate and duplicate these options in memory. All of the above rules apply. If I want to send an option to all tmpltbank and inspiral jobs running on H1 data, I might do something like [tmpltbank-H1&inspiral-H1].

If an option is given in more than one section (ie. if I specify --time-window 0.5 in [inspiral] and --time-window 1.0 in another [inspiral] or [inspiral&tmpltbank] or [inspiral-H1] the code will throw an error. Specifying --time-window 1.0 in [inspiral-H1] and --time-window 0.5 in [inspiral-L1] is valid as long as the subset of H1 jobs and the subset of L1 jobs do not overlap.

If a particular code (let's say inspiral) wants to use an option supplied in the [workflow] section (for e.g. the channel names) it can do this by using::

  [inspiral-h1]
  channel-name = ${workflow|h1-channel}

  [inspiral-l1]
  channel-name = ${workflow|l1-channel}

  [inspiral-v1]
  channel-name = ${workflow|v1-channel}

Similar macros can be added as needed, but these should be limited to avoid namespace confusion. 

------------------------------------
Example complete workflow .ini file
------------------------------------

Please see individual workflow documentation pages for some examples of complete .ini files and example workflows.

========================
[sharedoptions] section
========================

An alternative to the [exe1&exe2] section, especially when options are split
well into groups of options, is to use the [sharedoptions] section. An example of this follows::

  [sharedoptions]
  massranges = exe1,exe2,exe3-mass
  metric = exe1,exe2-range,exe3-metric, exe5

  [sharedoptions-massranges]
  min-mass1 = 2.0
  max-mass1 = 48.0
  min-mass2 = 2.0
  max-mass2 = 48.0
  max-total-mass = 4.2
  min-total-mass = 4.0
  max-eta = 0.25
  max-ns-spin-mag = 0.9899
  max-bh-spin-mag = 0.9899

  [sharedoptions-metric]
  pn-order = threePointFivePN
  f0 = 70.0
  f-low = 30.0
  f-upper = 1100.0
  delta-f = 0.01

This will ensure that all options in [sharedoptions-massranges] are added to the [exe1], [exe2] and [exe3-mass] sections. All options in [sharedoptions-metric] are added to [exe1], [exe2-range], [exe-metric] and [exe5].


=====================
Code documentation
=====================

The parsing of .ini files and command line parsing is done from within the pycbc.workflow.configuration module. The functions in this module are shown below

--------------------------------------------
:mod:`pycbc.workflow.configuration` Module
--------------------------------------------

.. automodule:: pycbc.workflow.configuration
    :noindex:
    :members:
    :undoc-members:
    :show-inheritance:

