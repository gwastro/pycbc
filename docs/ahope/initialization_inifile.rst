###################################################
The ahope .ini configuration file
###################################################

============
Introduction
============

One of the key inputs to ahope is the ahope.ini configuration file. This file contains all the configuration options to all of the code that will be used in the ahope pipeline. With a few exceptions (start as --gps-start, --gps-end, --output-file and similar options that change for every job), every command line option that a job is run with is inherited from this file.

Ihope's .ini files were used in every analysis. However, it has been noted that these files have grown huge and it is difficult for a novice to understand which options can be safely changed and which a novice would never want to go near. It is also difficult so see which options are going to which job, inspiral.c for example looks for options in > 10 places and it isn't clear where that is.

In this page we describe the layout of the ahope .ini configuration file and what the various sections mean, how they are used, and how an ini file should be set out. 

###################
Features to fit somewhere else on the page
###################

- aHope should be provided a nicely formatted ini file as input. After parsing the .ini file, aHope will dump the parsed .ini file (unordered and free of comments) back to the analysis directory. This file would be more useful for an expert to sanity check the workflow than the input file is.
- The ability to supply secondary ini files (as many as needed) within the ini file should be given. When parsing the ini file, aHope would begin by reading these in and concatenating them with the original configuration file.

###################
Global options - the [ahope] section
###################

The [ahope] section of the configuration file should be used to store options that are used when constructing an ahope workflow that are not supplied on the command line.

An example of where this section might be used is in the template bank stage where one can either run with a pre-generated bank or generate banks within the workflow. This information would be provided in this section.

If it proves necessary ahope "subsections" can be created with names like [ahope-segments] or [ahope-datafind] .....

----------------------
Requirements
----------------------

**Every** option given in the [ahope] section of the ini file should be preceeded by a comment explaining what that option does, and give links to the documentation page where these options are described if appropriate.

----------------------
Example
----------------------

Here is an example of the [ahope] section of a .ini file::

  [ahope]
  ; This should give a link for all the options that can be given in this section and what they do
  ; LINK HERE

  ; This is the lists of instruments to analyse
  analyse-ifos = H1,L1,V1

  ; Location of the veto-definer file to be used for generating veto segments
  veto-def-file = https://www.lsc-group.phys.uwm.edu/ligovirgo/cbc/public/segments/S6/H1L1V1-S6_CBC_LOWMASS_B_OFFLINE-937473702-0.xml
  ; This is the name of the science segments for H1,L1 and V1
  h1-science = H1:DMT-SCIENCE:4
  l1-science = L1:DMT-SCIENCE:4
  v1-science = V1:ITF_SCIENCEMODE:6
  ; This is the name of the channels to look at for H1,L1 and V1
  h1-channel = H1:LDAS-STRAIN
  l1-channel = L1:LDAS-STRAIN
  v1-channel = V1:HrecV3
  ; This is the frame type to use when finding the data
  h1-frame-type = H1_LDAS_C02_L2
  l1-frame-type = L1_LDAS_C02_L2
  v1-frame-type = V1_HrecV3
  ; veto categories to analyze (starting after cat1 which is always applied)
  veto-categories = 2,3,4
  ; Location of the segment database to get segments from
  segment-url = https://segdb.ligo.caltech.edu
  ; The following option is used to tell ahope to use a pre-generated template bank and gives the location
  tmplt-bank = /PATH/TO/TEMPLATE/BANK/template_bank.xml
  ; EXAMPLE TO BE IMPROVED WHEN AN INI FILE ACTUALLY EXISTS!

###################
Executable locations - the [executables] section
###################

This section should contain the names of each of the executables that will be used in the ahope workflow and their locations. 

-------------------
Example
-------------------

Here is an example of the [executables] section of an ahope .ini file::

  [executables]
  tmpltbank         = /home/cbc/opt/s6b/ab577e4e5dad14e46fce511cffdb04917836ba36/bin/lalapps_tmpltbank
  inspiral          = /home/cbc/opt/s6b/ab577e4e5dad14e46fce511cffdb04917836ba36/bin/lalapps_inspiral
  inspinj           = /home/cbc/opt/s6b/ab577e4e5dad14e46fce511cffdb04917836ba36/bin/lalapps_inspinj

###################
Executable options
###################

For each of the executables in the [executables] section, options for that executable should be listed under the section corresponding to that executable. Options in the [tmpltbank] section are sent to lalapps_tmpltbank, options in the [inspiral] section are sent to lalapps_inspiral etc.

It is possible to have more than one [tmpltbank] section, ConfigParser will simply combine them together when reading in. Therefore '''important options''' and '''options that a novice user might want to change''' will be supplied in a first [tmpltbank] section near the top of the .ini file. This section will be labelled accordingly and '''all options should be described'''. A link to documentation for that executable should also be given. Options that are not so important and ones that a novice user would not want to change will be placed in a second [tmpltbank] section at the bottom of the ini file, this section would be labelled accordingly and also contain a link to documentation for that executable.

Some options are only sent to a subset of jobs using a given executable. For example those running on H1 data. Options like these will be provided in sections labelled [executable_name-subset_tag]. So for the H1 example the section would be called [tmpltbank-H1]. As well as obeying the rules above these section must clearly state ''which'' jobs will be sent those options.

Some options need to be sent to more than one executable, for example the channel names are used by any code that reads the data. Such sections should be given as the combination of executable names separated by the & token. So options sent to tmpltbank '''and''' inspiral would go in a section called [tmpltbank&inspiral]. The code parsing the .ini file will automatically separate and duplicate these options in memory. All of the above rules apply. If I want to send an option to all tmpltbank and inspiral jobs running on H1 data, I might do something like [tmpltbank-H1&inspiral-H1].

If an option is given in more than one section (ie. if I specify --time-window 0.5 in [inspiral] and --time-window 1.0 in another [inspiral] or [inspiral&tmpltbank] or [inspiral-H1] the code will throw an error. Specifying --time-window 1.0 in [inspiral-H1] and --time-window 0.5 in [inspiral-L1] is valid as long as the subset of H1 jobs and the subset of L1 jobs do not overlap.

If a particular code (let's say inspiral) wants to use an option supplied in the [ahope] section (for e.g. the channel names) it can do this by using::

  [inspiral-h1]
  channel-name = MACRO_GET(ahope,h1-channel)

  [inspiral-l1]
  channel-name = MACRO_GET(ahope,l1-channel)

  [inspiral-v1]
  channel-name = MACRO_GET(ahope,v1-channel)

Similar macros can be added as needed, but these should be limited to avoid namespace confusion. This also means that no value supplied to an option can begin with MACRO (unless we add some functionality to allow this to be overridden (value = MACRO_MACRO(MACRO_SPOTTY_WOT) might be interpreted to value = MACRO_SPOTTY_WOT ). '''Any feature like this must be clearly documented'''

---------------------
Example complete ahope .ini file
---------------------

Provided here is an example of a complete ahope.ini file. ''File to be added once we start putting one together.''

#####################
Code documentation
#####################

The parsing of ahope .ini files is done from within the pycbc.ahope.configparserutils module. The functions in this module are documented below

.. autofunction:: pycbc.ahope.configparserutils
