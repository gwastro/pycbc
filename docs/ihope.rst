###################################################
Using pycbc_inspiral within ihope
###################################################

========================================
Introduction
========================================
These instructions are aimed at allowing the use of pycbc_inspiral within
standard ihope workflows. This is accomplished using a wrapper script that
emulates the interface of lalapps_inspiral. The standard ihope instructions
should be followed with modifications made to the ini file as noted below.

* `ihope instructions <https://www.lsc-group.phys.uwm.edu/ligovirgo/cbcnote/InspiralPipelineDevelopment/101119150805InspiralPipelineDocumentationHow_to_run_ihope_with_the_ligolw_thinca_single_stage_pipeline>`_

========================================
Ensure that you have the latest packages
========================================

Follow the instructions to install lalsuite and pycbc. 

* :doc:`install`
* `LALSuite <https://www.lsc-group.phys.uwm.edu/daswg/projects/lalsuite.html>`_ (with swig bindings enabled)


========================================
Ensure that pycbc is sourced
========================================
PyCBC creates a user environement source script. 

Either directly run

>>> source /path/to/pycbc/install/etc/pycbc-user-env.sh

or add it to your .bash_profile.

========================================
Modifications to the ini file
========================================

* Change the inspiral executable from lalapps_inspiral to pycbc_legacy_inspiral

    The following::

        [condor]
        ....
        inspiral          = /path/to/lalsuite/installation/bin/lalapps_inspiral
        ...

    is replaced by::

        [condor]
        ...
        inspiral          = /path/to/pycbc/installation/bin/pycbc_legacy_inspiral
        ...

* Add an pycbc flag to the [inspiral] section

    This indicates to the workflow generator that the inspiral jobs are in fact
    pycbc_inspiral::

        [inspiral]
        pycbc = 
        ...

* Modify the injection approximant names

    The injection routines within pycbc inpspiral are thin wrappers around functions
    in lalsimulation. As such, the names of approximants must match those available
    within lalsimulation. The phase order of the waveform should also not be appended
    to approximant name.
    
    For example::
    
        [bnsloginj]
        f-lower = 30
        waveform = TaylorT4threePointFivePN
        d-distr = log10
        min-distance = 5000
        
    becomes::

        [bnsloginj]
        f-lower = 30
        waveform = TaylorT4
        d-distr = log10
        min-distance = 5000
    
        
