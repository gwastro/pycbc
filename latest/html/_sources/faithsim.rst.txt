################################################
Dag Generator for Doing Faithfulness Comparisons
################################################

============
Introduction
============

This page describes how to use the faithfulness dag generator within
PyCBC.

==========================
How to generate a workflow
==========================

------------------------------------
Creating a configuration (.ini) file
------------------------------------

All the choices when setting up a faithsim are contained
in a single configuration file. 

Below is an example.

.. literalinclude:: ../examples/faith/faithsim_simple.ini

There are four sections that must be present [inspinj], [executables], [workflow], and [faithsim-XXX].

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

    This section lists the location of the pycbc_faithsim script. Make note
    that the script is copied to the executables folder
    and that is the version that will be used. 

 #. workflow 

    This section has options that configure the workflow.
    The required options are 'log-path' and 'templates-per-job'. The
    'log-path' specifies the directory to store
    condor log files. The 'templates-per-job' section determines
    how many faithfulness calculations each job will do. The 
    injection file is split into smaller portions to match this
    restriction. This option
    is directly proportional to the running time of 
    each job.

    Faith simulations running on LDG clusters must include the
    'accounting-group' option in the workflow section. The value
    must be choosen according to the
    `Accounting information web page <https://ldas-gridmon.ligo.caltech.edu/ldg_accounting/>`_.

 #. faithsim-XXX
    Multiple sections with a a name of the form 'faithsim-USER_STRING' can exist.
    The generator will create jobs that correspond to each of these sections 
    and each will generate an independent results file
    file labeled with the same USER_STRING.

    These sections corresponds to the arguments sent to the
    banksim executable. The two approximants to compare, along with their
    PN order paramters (if relevant), are set here. Note that
    the option filter-waveform-length must be set to a value
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

    pycbc_make_faithsim --conf YOUR_INI_FILE.ini

The workflow can then be submitted by running the
generated shell script. 

.. code-block:: bash

    sh submit.sh

-------------------------
Understanding the results
-------------------------

The main results of the faithsim are result files, one for each
faithsim section in the file. These are whitespace separated ASCII files.

Some basic plots are also generated automatically and 
placed into the 'plots' folder. 

The pycbc_faithsim_plots script located in the 
scripts folder is an example of 
how to read the results files. 

--------------------
Example config files
--------------------
*ADD SOME HERE*

Contact Alex Nitz for some more detailed examples of configuration files.
