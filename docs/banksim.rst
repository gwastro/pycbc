################################################################
Calculating the Effectualness (Fitting Factor) of Template Banks
################################################################

.. _banksim:

===================
Introduction
===================

This page describes how to use the ``banksim`` facilities within PyCBC. 
The ``banksim`` tools calculate the matches, maximized over a set of templates,
for a list of injections waveforms to measure the effectualness (fitting
factor) of a template bank. 

The purpose of this test is to allow the user to investigate the construction of new template banks, as well as act as a sanity check of the template bank generation methodology and code. Therefore the tests run over the same parameter ranges used to generate the bank and using the same sensitivity curve. The tests described here may not be optimal or exhaustive, but should be used to gain confidence that a bank is recovering signals it is designed to recover at an acceptable match.

-----------------------------
Creating a configuration file
-----------------------------

All the choices when setting up a banksim are contained
in a single configuration file. 

Below is an example.

.. literalinclude:: ../examples/banksim/banksim_simple.ini

There are four sections that must be present [inspinj], [executables], [workflow],
and [banksim].

 #. inspinj

    This section sets the paramaters of all of the injection waveforms. 
    The arguments in the configuration file are fed directly to the  
    lalapps_inspinj program to create an injection file.
    The same arguments are available, and the same restrictions apply. 
    The number of injections can be set by using the gps start and end time
    options along with the time step.
    Note, however, that the waveform name is required but does not 
    determine the actual approximants that will be compared. That is set in the [banksim] 
    section. 

 #. executables

    This section lists the location of the pycbc_banksim script. Make note
    that the script is copied to the executables folder
    and that is the version that will be used. 

 #. workflow 

    This section has options that configure the workflow.
    The required options are 'log-path', 'bank-file',
    'injections-per-job', and 'templates-per-job'. The
    'log-path' specifies the directory to store
    condor log files. 'bank-file' sets the template bank
    over which to maximize matches. It must be either 
    a sngl or sim inspiral table in xml format. 
    'injections-per-job' as its name suggests determines
    the maximum number of injections that each job 
    has to calculate fitting factors for. 

    The injection
    file generated from the [inspinj] section is split
    into smaller pieces to satisfy this requirement.
    Note that this option has a direct effect on the memory
    requirements of each banksim job, as each injection
    is pregenerated at the beginning of the program.
    
    The 'templates-per-job' will cause the given template
    bank to be split into smaller portions. This option
    is directly proportional to the running time of 
    each job.

    An optional value 'use-gpus' can be set. This will 
    set up the workflow to choose condor nodes with GPUS
    and sets the options so that the banksim program
    will use the GPU for accelerated processing. Note that
    the default is to treat all results from a GPU as
    unreliable. As such, each job is automatically run
    twice. The results are compared and only kept if 
    they equivelant. Only the GPUs on SUGAR and ATLAS
    are supported at this time.

    Bank simulations running on LDG clusters must include
    the 'accounting-group' option in the workflow section.
    The value must be choosen according to the
    `Accounting information web page <https://ldas-gridmon.ligo.caltech.edu/ldg_accounting/>`_.

 #. banksim

    This section corresponds to the arguments sent to the
    banksim executable. The notable exeption is that the
    correct flag for GPU support will be set if the 'use-gpus'
    option is set in the [workflow] section. The actual
    signal and template approximants, along with their
    PN order paramters (if relevant), are set here. Note that
    the option filter-buffer-length must be set to a value
    greater than the duration of the longest generated
    approximant.

------------------------
Generating the workflow
------------------------

Once a configuration file as been made, create a 
workspace directory and place the file into it. 
Running the following command will generate a dag
that will submit the required jobs. 

.. code-block:: bash

    pycbc_make_banksim --conf YOUR_INI_FILE.ini

The workflow can then be submitted by running the
generated shell script. 

.. code-block:: bash

    sh submit.sh

-------------------------
Understanding the results
-------------------------

The main results of the banksim is a single file called
'results.dat'. This is a space separated ASCII file.

Early (incomplete) results can be generated at any time
by executing the following script.

.. code-block:: bash

   sh partial_results.sh

Some basic plots are also generated automatically and 
placed into the 'plots' folder. 

The pycbc_banksim_plots script located in the 
scripts folder is an example of 
how to read the results file. 

=================================================
Validating template banks for production analysis
=================================================

To validate the uberbanks used in LIGO searches, we the BNS, NSBH and BBH regions, with separate banksim runs. Therefore there will be some overlap between the signals tested. For technical reasons, it is also convenient to split the NSBH and BBH tests up into a run with signals below a total mass of 50 and signals with a total mass above 50. 

We propose to select test signals from mass distributions that a flat in component masses in the respective regions; NS masses between 1 and 3 and BH masses between 2 and 99, with a total mass limit of 100. In addition, we select aligned spin magnitudes uniform in the respective regions; -0.05 to 0.05 for NS and -0.99 to 0.99 for BH. 

We propose to test with 10,000 injection signals in each of the BNS, NSBH and BBH regions, for a total of 30,000. This number is much less than the total number of templates in the bank. 

We propose to use SEOBNRv2 as the test signals, even though the uberbank uses TaylorF2 and SEOBNRv2_ROM_DoubleSpin templates for recovery. This is because we believe that SEOBNRv2 is a more accurate waveform than either TaylorF2 or the ROMs.

---------
Procedure
---------

The bank generation can be verified using the pycbc_banksim code. To run this follow the instructions for running the banksim code. An example ini file to run the NSBH banksim for total masses below 50 is given here

.. literalinclude:: ../examples/banksim/nsbh_below50.ini

To run this you will need to change the banksim option to your local version of pycbc_banksim, the log-path option to a suitable location for your log files on your cluster, the locations of the bank and noise curve and possibly whatever processing_scheme is best on your cluster (mkl works on Atlas with /opt/intel/2015/intel.sh sourced).

Injected spins are up to 0.99, not 0.9895 and the injections are uniform in component mass from 1 to 50 and uniform in spin magnitude (so it contains some highly spinning BNS). Injections are generated from 25Hz but matches are calculated from 30Hz, this gives the signal some "burn-in" time. Source location l-distr is random over the sky and inclination i-distr is uniformly distributed over arccos(i) - although this should not matter for aligned signals.

----------
Evaluation
----------

A stochastic placement method (like sbank) will not be able to guarantee that all points in parameter space are covered at better than 0.97 fitting factor. A convenient measure of the success of the bank generation is if the bank is able to recover 99% of injected signals using the same parameters and templates as the bank is designed for with a fitting factor of 0.97 or better. Further requirements might be that there should be no fitting factors with matches less than 0.95 or that the fitting factors below 0.97 should not be clustered in a particular part of parameter space. To cover all source groups we can run such tests separately for simulated BNS, NSBH and BBH signals when testing a bank that covers all three parameter ranges.

While such tests do not guarantee that the bank will successfully recover all possible signals in the parameter region (for example due different sensitivites in the two detectors, different waveform approximants, precession effects, tidal deformation and disruption etc.) these tests do indicate with a reasonable level of confidence that the template generation has been successful at what it was designed to do.

------------
Known issues
------------

The coverage of the high-mass (>70) and anti-aligned (<-0.5) NSBH region is known to be sparse in some versions.

The mchirp-window size may need to be changed if it is too tight. This is particularly a problem at higher masses.

If speed is an issue, the banksims can be sped up by reducing the number of injection signals, using ROMs instead of SEOBNRv2 as injection signals, reducing the signal-sample-rate or tightening the mchirp-window. Code is being developed to do this dynamically.

The option total-mass-divide is needed to replicate the uberbank switching from using TaylorF2 below total mass of 4 to using ROMs above. This may not exist on current master of pycbc_banksim.

