.. _configuration_files:

===========================
Configuration File Sets
===========================

----------------------
The example set
----------------------

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

----------------------------------
s6d_pycbc_review_twoweek_v1.ini
----------------------------------
The ini file is configured to reproduce the two-stage S6D results here but with a single-stage configuration.

* https://sugar-jobs.phy.syr.edu/~spxiwh/S6/all_sky/badgamma_reruns/968543943-971622087/
 This was the run that obtained the single-detector triggers for the Big Dog
 Plot in S6D.

* s6d_pycbc_review_twoweek_v1.ini

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

----------------------------------
``sngl/example_sngl.ini``
----------------------------------
This ini file is configured to run on ER6 data from 2014. It uses the zero-pole-gain (ZPK)
filtering module to dewhiten the channel ``OAF-CAL_DARM_DQ``.

* ``sngl/example_sngl.ini``

