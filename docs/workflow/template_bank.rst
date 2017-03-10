.. _workflowtmpltbankmod:

#############################################
The workflow template bank generation module
#############################################

=============
Introduction
=============

The template bank section of the pycbc workflow module is responsible for gathering/generating the banks of waveforms that will be used to matched-filter against the data.

It can run in a number of different modes

- With a pre-generated template bank, which can be specified for all ifos, or as one bank for each ifo.
- By generating unique and independent template banks for each ifo that are regenerated every approx. 2000s
- Generating template banks that do not vary over the workflow, can be the same template bank in each ifo or different ones.

The template bank module, by default, is independent of other modules, though it
is possible to ensure that there is a one-to-one correspondence between template banks and matched-filter outputs. 
Improvements over ihope:

- The template bank analysis chunks may be longer or shorter than the inspiral chunks
- Options sent to one job can be sent to the other or not, as desired.  No hardcoding of options that are **always** sent to both the matched-filter and template bank stages.

Like other modules, the template bank module returns a pycbc FileList of the template bank files generated/supplied to this module. It is possible to make multiple calls to the template bank module in the same module by using the tags kwarg when calling this.

======
Usage
======

Using this module requires a number of things

* A configuration file (or files) containing the information needed to tell this module how to generate (or gather) the template banks (described below).
* An initialized instance of the pycbc Workflow class, containing the ConfigParser.
* A list of segments to be analysed by this module.
* A FileList returned by the datafind module containing the frames that contain the data that will be used to make the template banks. (If using a pre-supplied PSD, or some other use-case that does not require reading data this can be set to None).

The module is then called according to

.. autofunction:: pycbc.workflow.setup_tmpltbank_workflow
   :noindex:

-------------------------
Configuration file setup
-------------------------

Here we describe the options given in the configuration file used in the
workflow that will be needed in this section

$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
[workflow-tmpltbank] section
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

The configuration file must have a [workflow-tmpltbank] section,
which is used to
tell the workflow how to construct (or gather) the template banks. The first option to choose and provide is

* tmpltbank-method = VALUE

The choices here and their description are as described below

* PREGENERATED_BANK - A pregenerated template bank is supplied and this should be used when performing matched-filtering for all ifos and all times. This uses the setup_tmpltbank_pregenerated sub-module.
* WORKFLOW_INDEPENDENT_IFOS - Template banks will be generated within the workflow. These banks will be made to cover only short (normally ~ 2000s) of data to reflect PSD changes over time and will be independent and distinct for each analysed interferometer. This uses the setup_tmpltbank_dax_generated sub-module.
* WORKFLOW_INDEPENDENT_IFOS_NODATA - Template banks will be generated within the workflow. There will be one bank for each ifo, which will cover all times. No data frames will be used when constructing the workflow (ie. using a design PSD or similar). This uses the setup_tmpltbank_without_frames sub-module.
* WORKFLOW_NO_IFO_VARIATION_NODATA - As WORKFLOW_INDEPENDENT_IFOS_NODATA except only one template bank will be generated that is valid for all ifos. This uses the setup_tmpltbank_without_frames sub-module.

Each of these options will describe which subfunction to use. These are described here

.. autofunction:: pycbc.workflow.setup_tmpltbank_pregenerated
   :noindex:

.. autofunction:: pycbc.workflow.setup_tmpltbank_dax_generated
   :noindex:

.. autofunction:: pycbc.workflow.setup_tmpltbank_without_frames
   :noindex:

When using the setup_tmpltbank_pregenerated sub-module the following additional options apply in the [workflow-tmpltbank] section.

* tmpltbank-pregenerated-bank = PATH - OPTIONAL. This is the location of the pre-generated bank that is to be used for all ifos.
* tmpltbank-pregenerated-bank-[ifo] = PATH = OPTIONAL. This is the location of the pre-generated bank that is to be used for all ifos. Either this option must be given for all active ifos or tmpltbank-pregenerated-bank must be given. If both are supplied this option tmpltbank-pregenerated-bank-[ifo] will be ignored. [ifo] should be replaced by the name of each ifo in lower case. For example tmpltbank-pregenerated-bank-l1 or tmpltbank-pregenerated-bank-v1.

When using the setup_tmpltbank_without_frames sub-module the following additional options apply in the [workflow-tmpltbank] section:

* tmpltbank-write-psd-file - OPTIONAL. If given the template bank code will also write a file containing the PSD. Currently only supported in pycbc template bank codes.

When using the setup_tmpltbank_dax_generated sub-module the following additional options apply in the [workflow-templtbank] section.

* tmpltbank-link-to-matchedfltr - OPTIONAL. If this is given the workflow module will attempt to ensure a one-to-one correspondence between template banks and matched-filter outputs. This may not work in all cases and should be considered an option to be used for comparing with ihope output.
* tmpltbank-compatibility-mode - OPTIONAL. If this is given the workflow module will tile the template bank jobs in the same way as inspiral_hipe used to. This requires the link option above and that the template bank and matched-filtering jobs are reading the same amount of data in each job.
* tmpltbank-write-psd-file - OPTIONAL. If given the template bank code will also write a file containing the PSD. Currently only supported in pycbc template bank codes.

The following options apply only when using setup_tmpltbank_dax_generated and not using lalapps_tmpltbank

* analysis-length = LENGTH_IN_SECONDS (not used when lalapps_tmpltbank is the executable) - REQUIRED. The amount of time in seconds that will be used for frames for PSD generation.

$$$$$$$$$$$$$$$
[executables]
$$$$$$$$$$$$$$$

If using setup_tmpltbank_pregenerated then no executables are needed, you have already generated the template bank.

If using setup_tmpltbank_dax_generated or setup_tmpltbank_without_frames then you need to supply the template bank executable. This is done in the [executables] section by adding something like:

tmpltbank = /path/to/lalapps_tmpltbank

the option, in this case 'tmpltbank', will be used to specify the constant command line options that are sent to all tmpltbank jobs. Currently this value is hardcoded to tmpltbank, but we plan to change this to allow multiple tmpltbank executables to be used in a single workflow. How to set up the [{exe_name}] section, and which executables are currently supported is discussed below.

If template banks are being generated separately for each ifo then sections named [tmpltbank-${ifo}] (for e.g. [tmpltbank-h1] or [tmpltbank-v1]) can be used to specify options that should only be supplied when running on that ifo.

-------------------------------------------------------------
Supported template bank exes and instructions for using them
-------------------------------------------------------------

The following template bank executables are currently supported in the workflow module

* lalapps_tmpltbank
* pycbc_geom_nonspinbank

Adding a new executable is not too hard, please ask a developer for some pointers on how to do this if you want to add a new code.

Also see :ref:`Pycbc.Tmpltbank <tmpltbankmod>` for a description/introduction to pycbc template bank code and papers describing these codes, lalapps_tmpltbank and sBank.

$$$$$$$$$$$$$$$$$$$$$$$$
lalapps_tmpltbank_ahope
$$$$$$$$$$$$$$$$$$$$$$$$

Lalapps_tmpltbank is the legacy C-code that has been used to generate template banks for gravitational-wave data analysis since the dawn of time. It is a little inflexible in terms of output file names. We recommend using the newer pycbc_geom_nonspinbank if possible.

lalapps_tmpltbank is supported in pycbc's workflow module via a wrapper script lalapps_tmpltbank_ahope, this allows us to specify all the frame files and the output file name directly.

The help message for lalapps_tmpltbank is available at http://software.ligo.org/docs/lalsuite/lalapps/tmpltbank_8c.html

Of these options the workflow module and/or the wrapper script will automatically add the following, which are unique for each job. **DO NOT ADD THESE OPTIONS IN THE CONFIGURATION FILE**.

* --gps-start-time
* --gps-end-time
* --frame-cache
* --user-tag
* --ifo-tag

All other options must be provided in the configuration file. Here is an example of a lalapps_tmpltbank call.

.. code-block:: bash

   lalapps_tmpltbank --grid-spacing Hexagonal --dynamic-range-exponent 69.0 --minimal-match 0.97 --high-pass-order 8 --strain-high-pass-order 8 --maximum-mass 25.0 --gps-end-time 961587599 --calibrated-data real_8 --channel-name H1:LDAS-STRAIN --space Tau0Tau3 --number-of-segments 15 --enable-high-pass 30.0 --gps-start-time 961585551 --high-pass-attenuation 0.1 --num-freq-cutoffs 1 --segment-length 1048576 --low-frequency-cutoff 40.0 --pad-data 8 --min-high-freq-cutoff SchwarzISCO --sample-rate 4096 --high-frequency-cutoff 2048.0 --resample-filter ldas --strain-high-pass-atten 0.1 --strain-high-pass-freq 30 --max-total-mass 25.0 --frame-cache /home/spxiwh/lscsoft_git/src/pycbc/examples/ahope/weekly_ahope/961585543-961671943/datafind/H1-DATAFIND-961585543-86400.lcf --disable-compute-moments  --max-high-freq-cutoff SchwarzISCO --approximant TaylorF2 --write-compress  --minimum-mass 1.0 --order twoPN --spectrum-type median

$$$$$$$$$$$$$$$$$$$$$$$$
pycbc_geom_nonspinbank
$$$$$$$$$$$$$$$$$$$$$$$$

pycbc_geom_nonspinbank is pycbc's non-spinning template bank generator. Designed as a replacement and improvement of lalapps_tmpltbank. See :ref:`tmpltbankmod` for documentation of this code. The help message of pycbc_geom_nonspinbank follows:

.. command-output:: pycbc_geom_nonspinbank --help

Of these options the workflow module will automatically add the following, which are unique fo
r each job. **DO NOT ADD THESE OPTIONS IN THE CONFIGURATION FILE**.

* --gps-start-time
* --gps-end-time
* --frame-cache
* --output-file

All other options must be provided in the configuration file. Here is an example of a pycbc_geom_nonspinbank call.

.. code-block:: bash

   pycbc_geom_nonspinbank --pn-order twoPN --f0 40 --f-low 40 --f-upper 2048 --delta-f 0.01 --min-match 0.97 --min-mass1 2.0 --min-mass2 2.0 --max-mass1 3. --max-mass2 3. --verbose --output-file testNonSpin.xml --calculate-ethinca-metric --filter-cutoff SchwarzISCO --psd-estimation median --psd-segment-length 256 --psd-segment-stride 128 --psd-inverse-length 8 --gps-start-time 900000033 --gps-end-time 900002081 --strain-high-pass 30 --pad-data 8 --sample-rate 4096 --frame-cache cache/H-H1_NINJA2_G1000176_EARLY_RECOLORED_CACHE-900000024-10653.lcf --channel-name H1:LDAS-STRAIN --max-total-mass 5.5 --min-total-mass 4.5

==========================================
:mod:`pycbc.workflow.tmpltbank` Module
==========================================

This is complete documentation of this module's code

.. automodule:: pycbc.workflow.tmpltbank
    :noindex:
    :members:
    :undoc-members:
    :show-inheritance:
