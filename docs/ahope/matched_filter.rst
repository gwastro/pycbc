.. _ahopeinspiralmod:

###############################
The ahope matched-filter module
###############################

=============
Introduction
=============

The matched-filter section of ahope is responsible for matched-filtering the
data against the template bank(s) from the template bank section and generating
a list of "triggers" for each interferometer. These triggers should be a list
of any event where the signal to noise ratio and any signal consistency test
are such that that point should be sent forward to check for coincidence in
other ifos.

Any single-ifo signal consistency tests (ie. chi-squared tests etc.) should
be computed in this section and stored within the lists of triggers. Ahope
does not make any specifications on the output format, but obviously code in
the next stages of the workflow must know how to process that input.

The matched-filtering section should be as independent of the other stages of
the workflow as possible. This means that we don't require the data read in by
matched-filter jobs to match that read in by template bank jobs (however this
may be desireable in some cases, so **should** be possible where sensible).
Options should also not be hardcoded (so there are no cases where an option
that gets sent to a template bank job also gets sent to a matched-filter job
without any way to stop that). However, it is possible to duplicate options
where this is desireable (see :ref:`ahopeconfigparsermod`).

The return from the matched-filter section of ahope is a list of AhopeOutFile
objects corresponding to each actual file (one for each job) that will be 
generated within the workflow. This will be the **only** thing that will be 
passed from the matched-filter section to the future sections.

======
Usage
======

Using this module requires a number of things

* A configuration file (or files) containing the information needed to tell this module how to generate GW triggers.
* An initialized instance of the ahope workflow class, containing the ConfigParser.
* A list of segments to be analysed by this module.
* An AhopeFileList returned by the templatebank module containing the template banks
available for use.
 
The module is then called according to

.. autofunction:: pycbc.ahope.setup_matchedfltr_workflow
   :noindex:

-------------------------
Configuration file setup
-------------------------

Here we describe the options given in the configuration file used in the ahope
workflow that will be needed in this section

$$$$$$$$$$$$$$$$$$$$$$$$$$
[ahope-inspiral] section
$$$$$$$$$$$$$$$$$$$$$$$$$$

The configuration file requires an [ahope-inspiral] section only in the case
where an option is required by the intended executables. See below for options
required by each executable.

The following options apply only to a specific executable.

* analysis-length = LENGTH_IN_SECONDS (pycbc_inspiral only) - REQUIRED. The amount of time in seconds that will be matched-filtered. Note that triggers may not be produced for the entire span of time. 

$$$$$$$$$$$$$$$
[executables]
$$$$$$$$$$$$$$$

inspiral = /path/to/inspiral_exec

A section, in this case [inspiral], will be used to specify the constant command line options that are sent to all inspiral jobs. How to set up the [{exe_name}] section, and which executables are currently supported is discussed below.

---------------------------------------------------------------------
Supported inspiral trigger generators and instructions for using them
---------------------------------------------------------------------

The following inspiral trigger generators are currently supported in ahope

* lalapps_inspiral
* pycbc_inspiral

Adding a new executable is not too hard, please ask a developer for some pointers on how to do this if you want to add a new code.

$$$$$$$$$$$$$$$$$$
lalapps_inspiral
$$$$$$$$$$$$$$$$$$

Lalapps_inspiral is the legacy C-code that has been used for years to find gravitational-wave triggers in  It is a little inflexible in terms of output file names.

lalapps_inspiral is supported in ahope via a wrapper script lalapps_inspiral_ahope, this allows us to specify all the frame files and the output file name directly.

.. command-output:: lalapps_inspiral --help

Of these options ahope or the wrapper script will automatically add the following, which are unique for each job. **DO NOT ADD THESE OPTIONS IN THE CONFIGURATION FILE**.

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

Of these options ahope will automatically add the following, which are unique fo
r each job. **DO NOT ADD THESE OPTIONS IN THE CONFIGURATION FILE**.

* --gps-start-time
* --gps-end-time
* --frame-cache
* --output

All other options must be provided in the configuration file. Here is an example of a pycbc_inspiral call.

.. code-block:: bash

   pycbc_inspiral --trig-end-time 961592867 --verbose  --cluster-method window --bank-filetmpltbank/L1-TMPLTBANK_01-961591486-1382.xml.gz --gps-end-time 961592884 --channel-name L1:LDAS-STRAIN --processing-scheme cuda --snr-threshold 5.5 --psd-estimation median --trig-start-time 961591534 --gps-start-time 961590836 --chisq-bins 16 --segment-end-pad 16 --segment-length 2048 --low-frequency-cutoff 15 --pad-data 8 --cluster-window 1 --sample-rate 4096 --segment-start-pad 650 --psd-segment-stride 32 --psd-inverse-length 16 --psd-segment-length 64 --frame-cache datafind/L1-DATAFIND-961585543-7349.lcf --approximant SPAtmplt --output inspiral/L1-INSPIRAL_1-961591534-1333.xml.gz --strain-high-pass 30 --order 7 


==========================================
:mod:`pycbc.ahope.matchedfltr_utils` Module
==========================================

This is complete documentation of this module's code

.. automodule:: pycbc.ahope.matchedfltr_utils
    :noindex:
    :members:
    :undoc-members:
    :show-inheritance:
