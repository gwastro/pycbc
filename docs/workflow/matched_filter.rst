.. _workflowinspiralmod:

###################################
The workflow matched-filter module
###################################

=============
Introduction
=============

The matched-filter section of pycbc's workflow module
is responsible for matched-filtering the
data against the template bank(s) from the template bank section and generating
a list of "triggers" for each interferometer. These triggers should be a list
of any event where the signal to noise ratio and any signal consistency test
are such that that point should be sent forward to check for coincidence in
other ifos.

Any single-ifo signal consistency tests (ie. chi-squared tests etc.) should
be computed in this section and stored within the lists of triggers. The
workflow
does not make any specifications on the output format, but obviously code in
the next stages of the workflow must know how to process that input.

The matched-filtering section should be as independent of the other stages of
the workflow as possible. This means that we don't require the data read in by
matched-filter jobs to match that read in by template bank jobs (however this
may be desirable in some cases, so **should** be possible where sensible).
Options should also not be hardcoded (so there are no cases where an option
that gets sent to a template bank job also gets sent to a matched-filter job
without any way to stop that). However, it is possible to duplicate options
where this is desireable (see :ref:`workflowconfigparsermod`).

The return from the matched-filter section of the workflow module
is a list of File
objects corresponding to each actual file (one for each job) that will be 
generated within the workflow. This will be the **only** thing that will be 
passed from the matched-filter section to the future sections.

======
Usage
======

Using this module requires a number of things

* A configuration file (or files) containing the information needed to tell this module how to generate GW triggers.
* An initialized instance of the workflow class, containing the ConfigParser.
* A list of segments to be analysed by this module.
* A FileList returned by the templatebank module containing the template banks available for use.
* A FileList returned by the datafind module containing the frames that contain the data that will be used to make the template banks.
* If desired an injection file for assessing sensitivity to simulated signals.
 
The module is then called according to

.. autofunction:: pycbc.workflow.setup_matchedfltr_workflow
   :noindex:

-------------------------
Configuration file setup
-------------------------

Here we describe the options given in the configuration file used in the 
workflow that will be needed in this section

$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
[workflow-matchedfilter] section
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

The configuration file must have an [workflow-matchedfilter] section, which is used to tell the workflow how to construct (or gather) the template banks. The first option to choose and provide is

* matchedfilter-method = VALUE

The choices here and their description are as described below

* WORKFLOW_INDEPENDENT_IFOS - Matched-filter trigger files will be generated within the workflow. These banks will be made to cover only short (normally ~ 2000s) of data to reflect PSD changes over time and will be independent and distinct for each analysed interferometer. This uses the setup_matchedfltr_dax_generated sub-module.

Currently only one option, but others can be added. The subfunctions used are described here

.. autofunction:: pycbc.workflow.setup_matchedfltr_dax_generated
      :noindex:

When using the setup_matchedfltr_dax_generated sub-module the following additional options apply in the [workflow-matchedfilter] section:

* matchedfilter-link-to-tmpltbank - OPTIONAL. If this is given the workflow module will attempt to ensure a one-to-one correspondence between template banks and matched-filter outputs. This may not work in all cases and should be considered an option to be used for comparing with ihope output.
* matchedfilter-compatibility-mode - OPTIONAL. If this is given the workflow module will tile the matched-filter jobs in the same way as inspiral_hipe used to. This requires the link option above and that the template bank and matched-filtering jobs are reading the same amount of data in each job.
* max-analysis-segments = (*NOT* used for lalapps_inspiral) - REQUIRED. The maximum number of analysis segments to analyze within a single inspiral job. Note that triggers may not be produced for the entire span of time. 
* min-analysis-segments = (*NOT* used for lalapps_inspiral) - REQUIRED. The minimum number of analysis segments to analyze within a single inspiral job. This may be the same as the maximum.

$$$$$$$$$$$$$$$
[executables]
$$$$$$$$$$$$$$$

inspiral = /path/to/inspiral_exec

A section, in this case [inspiral], will be used to specify the constant command line options that are sent to all inspiral jobs. How to set up the [{exe_name}] section, and which executables are currently supported is discussed below.

-----------------------------------------------------------------------
Supported inspiral trigger generators and instructions for using them
-----------------------------------------------------------------------

The following inspiral trigger generators are currently supported in pycbc's
workflow module

* lalapps_inspiral
* pycbc_inspiral

Adding a new executable is not too hard, please ask a developer for some pointers on how to do this if you want to add a new code.

$$$$$$$$$$$$$$$$$$$$$$$
lalapps_inspiral_ahope
$$$$$$$$$$$$$$$$$$$$$$$

Lalapps_inspiral is the legacy C-code that has been used for years to find gravitational-wave triggers in  It is a little inflexible in terms of output file names.

lalapps_inspiral is supported in the workflow module via a wrapper script lalapps_inspiral_ahope, this allows us to specify all the frame files and the output file name directly. Documentation for the lalapps_inspiral command line arguments can be found at http://software.ligo.org/docs/lalsuite/lalapps/inspiral_8c.html

Of these options the workflow module or the wrapper script will automatically add the following, which are unique for each job. **DO NOT ADD THESE OPTIONS IN THE CONFIGURATION FILE**.

* --gps-start-time
* --gps-end-time
* --trig-start-time
* --trig-end-time
* --frame-cache
* --user-tag
* --ifo-tag

All other options must be provided in the configuration file. Here is an example of a lalapps_inspiral call.

.. code-block:: bash

   lalapps_inspiral --do-rsq-veto  --trig-end-time 971614817 --enable-rsq-veto  --dynamic-range-exponent 69.0 --autochisq-stride 2 --bank-file datafind/L1-TMPLTBANK_19-971612833-2048.xml.gz --high-pass-order 8 --strain-high-pass-order 8 --ifo-tag FULL_DATA --user-tag 19 --gps-end-time 971614881 --calibrated-data real_8 --channel-name L1:LDAS-STRAIN --snr-threshold 5.5 --cluster-method template --number-of-segments 15 --trig-start-time 971613852 --enable-high-pass 30.0 --gps-start-time 971612833 --enable-filter-inj-only  --maximization-interval 30 --high-pass-attenuation 0.1 --chisq-bins 2 --inverse-spec-length 16 --rsq-veto-threshold 15.0 --segment-length 1048576 --low-frequency-cutoff 40.0 --pad-data 8 --autochisq-two-sided  --sample-rate 4096 --chisq-threshold 10.0 --rsq-veto-max-snr 12.0 --resample-filter ldas --strain-high-pass-atten 0.1 --strain-high-pass-freq 30 --bank-veto-time-freq  --segment-overlap 524288 --frame-cache datafind/L1-DATAFIND-968556757-3058132.lcf --chisq-delta 0.2 --bank-veto-subbank-size 20 --approximant FindChirpSP --rsq-veto-time-thresh 0.0002 --write-compress  --autochisq-length 100 --enable-output  --rsq-veto-window 6.0 --order threePointFivePN --spectrum-type median 


$$$$$$$$$$$$$$$$$$$$$$$$
pycbc_inspiral
$$$$$$$$$$$$$$$$$$$$$$$$

pycbc_inspiral is pycbc's inspiral matched-filtering program. Designed as a replacement and improvement of lalapps_inspiral. The help message of pycbc_inspiral follows:

.. command-output:: pycbc_inspiral --help

Of these options the workflow module will automatically add the following, which are unique fo
r each job. **DO NOT ADD THESE OPTIONS IN THE CONFIGURATION FILE**.

* --gps-start-time
* --gps-end-time
* --frame-cache
* --output

All other options must be provided in the configuration file. Here is an example of a pycbc_inspiral call.

.. code-block:: bash

   pycbc_inspiral --trig-end-time 961592867 --verbose  --cluster-method window --bank-filetmpltbank/L1-TMPLTBANK_01-961591486-1382.xml.gz --gps-end-time 961592884 --channel-name L1:LDAS-STRAIN --processing-scheme cuda --snr-threshold 5.5 --psd-estimation median --trig-start-time 961591534 --gps-start-time 961590836 --chisq-bins 16 --segment-end-pad 16 --segment-length 2048 --low-frequency-cutoff 15 --pad-data 8 --cluster-window 1 --sample-rate 4096 --segment-start-pad 650 --psd-segment-stride 32 --psd-inverse-length 16 --psd-segment-length 64 --frame-cache datafind/L1-DATAFIND-961585543-7349.lcf --approximant SPAtmplt --output inspiral/L1-INSPIRAL_1-961591534-1333.xml.gz --strain-high-pass 30 --order 7 


===========================================
:mod:`pycbc.workflow.matched_filter` Module
===========================================

This is complete documentation of this module's code

.. automodule:: pycbc.workflow.matched_filter
    :noindex:
    :members:
    :undoc-members:
    :show-inheritance:
