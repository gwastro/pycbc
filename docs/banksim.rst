###################################################################
Dag Generator for Doing Banksims
###################################################################

===================
Introduction
===================

This page describes how to use the banksim facilities within PyCBC. 
A banksim calculates the matches, maximized over a set of templates,
for a list of injections waveforms. 

================================================
How to generate a workflow to do a large banksim
================================================

------------------------
Creating an INI file
------------------------

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

-----------------------
Workflow considerations
-----------------------
*ADD SOME HERE*

-------------------------
Example config files
-------------------------
*ADD SOME HERE*
