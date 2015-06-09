.. _configuration_files:

###########################
Configuration File Sets
###########################

===========================
The example set
===========================

For instructional purposes, we provide an example configuration set that creates a three detector H1L1V1 search on S6 data. It is broken into three pieces, the core search, the injection sets, 
and the postprocessing.

There are two options for the core search. 

* example_pycbc.ini
 Use PyCBC executables for template placement and matched-filtering.
* example_lalapps.ini
 Use Lalapps executables for template placement and matched-filtering.

The post-processing and injections are available in the following ini files. 

* example_pipedown.ini
* example_inj.ini

=====================================
latest review configuration file
=====================================
The ini file is configured to reproduce the latest reviewed coincident pipeline configuration.

* http://atlas1.atlas.aei.uni-hannover.de/~cbiwer/ahope_review/pycbc_review/s6d/twoweek_pycbc_tmpltbank/970012743-971622087/
  This was the latest review run using both pycbc_inspiral and pycbc_geom_nonspin_tmpltbank.

* ``review/s6d_twoweek_review.ini``

===========================
daily configuration files
===========================

----------------------------------
``daily/example_daily_lalapps.ini``
----------------------------------
This ini file is configured to run on ER5 data.

* ``daily/example_daily_lalapps.ini``

----------------------------------
``daily/example_daily_pycbc_zpk.ini``
----------------------------------
This ini file is configured to run on LLO lock data from 2014. It uses the zero-pole-gain (ZPK)
filtering module to dewhiten the channel ``OAF-CAL_DARM_DQ``.

* ``daily/example_daily_pycbc_zpk.ini``

===========================
sngl configuration files
===========================

----------------------------------
``sngl/example_sngl.ini``
----------------------------------
This ini file is configured to run on ER6 data from 2014. It uses the zero-pole-gain (ZPK)
filtering module to dewhiten the channel ``OAF-CAL_DARM_DQ``.

* ``sngl/example_sngl.ini``

