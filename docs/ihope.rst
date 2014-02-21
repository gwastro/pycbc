###################################################
Using pycbc executables within ihope
###################################################

========================================
Introduction
========================================
These instructions are aimed at allowing the use of a few pycbc executables,
currently pycbc_inspiral and pycbc_geom_nonspinbank, within standard ihope
workflows. This is accomplished using a wrapper script that emulates the
interface of the appropriate lalapps executable. The standard ihope instructions
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

* To use pycbc_inspiral rather than lalapps_inspiral, change the inspiral executable from
  lalapps_inspiral to pycbc_legacy_inspiral:

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


  Similarly, to use pycbc_geom_nonspinbank rather than lalapps_tmpltbank, change the
  tmpltbank executable from lalapps_tmpltbank to pycbc_legacy_tmpltbank

    The following::

        [condor]
        ....
        tmpltbank          = /path/to/lalsuite/installation/bin/lalapps_tmpltbank
        ...

    is replaced by::

        [condor]
        ...
        tmpltbank          = /path/to/pycbc/installation/bin/pycbc_legacy_tmpltbank
        ...


* Add a pycbc flag to the appropriate section

    You must also indicate to the workflow generator that the corresponding executables
    are in fact pycbc executables.  This is done by adding a pycbc flag to the appropriate
    section of the ini file; that will be [inspiral] if you are using pycbc_legacy_inspiral
    and [tmpltbank] if you are using pycbc_legacy_tmpltbank; the example below is for the
    inspiral executable::

        [inspiral]
        pycbc = 
        ...

* For pycbc_legacy_inspiral, modify the injection approximant names

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
    
* For pycbc_legacy_tmpltbank, modify the post-Newtonian order of the bank

    This step is not strictly necessary, since you may use the same order as
    lalapps_tmpltbank.  However, one reason you might want to run using
    pycbc_legacy_tmpltbank is the ability to place banks at a higher PN order.
    If this is the case, then you may for example use::

        [tmpltbank]
        ...
        order = threePointFivePN
        ...
        
    to place the bank templates using 3.5 PN order
