.. _tmpltbankmod:

###################################################################
PyCBC template bank generation documentation (``pycbc.tmpltbank``)
###################################################################

===================
Introduction
===================

This page gives details on how to use the various bank generation codes
and bank generation workflows
available in pyCBC and specifically the pyCBC.tmpltbank module.
We also describe the sBank code, currently in lalinspiral.
These codes currently consist of 

* :ref:`A program to generate a non-spinning template bank (pycbc_geom_nonspinbank) <tmpltbank_nonspinbank>`
* :ref:`A program to generate an aligned-spin template bank using a geometrical lattice algorithm (pycbc_geom_aligned_bank) <tmpltbank_alignedgeombank>`
* :ref:`Two programs to generate an aligned-spin template bank using a stochastic placement algorithm. This is done either using a metric to approximate distances, or by computing explicit matches from waveforms. (pycbc_aligned_stoch_bank, lalapps_cbc_sbank) <tmpltbank_alignedstochbank>`
* :ref:`A workflow generator for generating large template banks using the lalapps_cbc_sbank executables, but allowing the parallel use of multiple cores to optimize the generation of large template banks. (pycbc_create_sbank_workflow) <tmpltbank_sbankworkflow>`
* :ref:`A workflow generator for generating the "uberbank". This consists of a chain of calls to the geometric aligned-spin lattice code and the sbank workflow generator to create a template bank suitable for analysing Advanced LIGO and Advanced Virgo data. (pycbc_create_uberbank_workflow) <tmpltbank_uberbankworkflow>`

Each of these codes and workflow generators is described in turn and is accompanied by some examples of how to run the code, some relevant background and references and the options described from the codes' help messages.

The output of all these codes is currently a sngl_inspiral_table XML object.
However, if people want to output in HDF or other formats it should be
easy to write a back-end module to do this. The code may not fill every
necessary column of the sngl_inspiral table, please let me know if something
is missing and it will be added!

**NOTE: Aligned-spin template generation is not at the point where it can be considered "plug and play". Here we give some examples for distinct situations, but these options may not be suitable in other situations, for example if changing the mass ranges, or PSD, or frequency ranges etc. Template banks should always be checked (see :ref:`pycbc_banksim <banksim>`) and you should consider if you have placed too many templates ("overcovering" the parameter space), or if the placement has left gaps in the coverage ("undercovering" the parameter space). There are numerous things that can be tuned in the workflow generators, we try to describe them here, but it is often good to ask questions if you are unsure!**

**In some cases recommendations are given for which codes to use where. It
should be noted that these are personal recommendations of the page maintainer
(Ian), and there is always a strong possibility that Ian, and the
recommendations he makes, are wrong. Where multiple methods exist, it is
advised that the user try out the different methods and decide for themselves
which method is best.**

.. _tmpltbank_nonspinbank:

=================================================
Non spin bank placement: pycbc_geom_nonspinbank
=================================================

---------------------
Overview
---------------------

pycbc_geom_nonspinbank is designed as a drop in replacement for 
lalapps_tmpltbank. It uses a globally flat coordinate system, which enables it
to easily place banks at 3.5PN order (and higher order as these become available).

Mapping back from the new coordinate system to masses and spins makes the bank
placement slower than lalapps_tmpltbank. However, it doesn't waste stupid
amounts of time calculating the ethinca metric so it is faster than tmpltbank
if the ethinca components are being computed.

A metric-based approach is not very reliable at high masses, but we have found
that in general this just results in overcoverage, which is not normally a problem as most templates are found at lower masses. One might want to consider the stochastic code sbank if placing a high-mass only non-spinning template bank.

**Ian's recommendation: Use this code instead of lalapps_tmpltbank if you want
a non-spinning bank. If only placing high-mass templates consider using lalapps_cbc_sbank, described below.**

----------------
Background
----------------

The methods, specifically the metric and the lattice algorithm, used in
lalapps_tmpltbank are described in:

* Sathyaprakash, B. S. and Dhurandhar, S. V., Phys.Rev. D44 3819  (1991)
* Poisson, Eric and Will, Clifford M., Phys.Rev. D52 848 (1995)
* Balasubramanian, R. et al., Phys.Rev. D53 3033 (1996)
* Owen, Benjamin J., Phys.Rev. D53 6749 (1996)
* Owen, Benjamin J. and Sathyaprakash, B. S., Phys.Rev. D60 022002 (1998)
* Babak, S. et al., Class.Quant.Grav. 23, 5477 (2006)
* Cokelaer, Thomas Phys.Rev. D76 102004 (2007)
* Babak, S. et al., Phys.Rev. D87 024033 (2013)

Our approach still uses the hexagonal lattice algorithm but instead of using
the chirp times as coordinates uses the xi_i coordinate system described in:

* Brown, Duncan A. et al. Phys.Rev. D86 084017 (2012)
* Ohme, Frank et al. Phys.Rev. D88 024002 (2013)

This allows us to easily use higher order terms in the metric placement as we
do not have to worry about correlations between tau 0 and 3 and the other
coordinates. Our coordinates are independent.

The issue is that unlike in the lalapps_tmpltbank approach it is not trivial
to map back from our coordinates to masses. Currently we use a monte-carlo
approach to find mass positions fitting the desired coordinates. If we could
implement some form of 2D interpolation this would probably be much faster.
If anyone wants to investigate this I can provide lots of data!

-------------------
Some examples
-------------------

This first example computes a non-spinning template bank using a supplied ASD text file. It also calculates the ethinca metric components for a filter with ISCO cutoff

.. code-block:: bash

    pycbc_geom_nonspinbank --pn-order twoPN --f0 40 --f-low 40 --f-upper 2048 --delta-f 0.01 --min-match 0.97 --min-mass1 2.0 --min-mass2 2.0 --max-mass1 3. --max-mass2 3. --verbose --output-file testNonSpin.xml --calculate-ethinca-metric --filter-cutoff SchwarzISCO --asd-file ZERO_DET_high_P.txt

Here is a second example that computes a non-spinning template bank from a frame
cache file. This needs a bunch more options to read data and estimate the PSD
from that.

.. code-block:: bash

    pycbc_geom_nonspinbank --pn-order threePointFivePN --f0 40 --f-low 40 --f-upper 2048 --delta-f 0.01 --min-match 0.97 --min-mass1 2.0 --min-mass2 2.0 --max-mass1 3. --max-mass2 3. --verbose --output-file testNonSpin.xml --calculate-ethinca-metric --filter-cutoff SchwarzISCO --psd-estimation median --psd-segment-length 256 --psd-segment-stride 128 --psd-inverse-length 8 --gps-start-time 900000033 --gps-end-time 900002081 --strain-high-pass 30 --pad-data 8 --sample-rate 4096 --frame-cache cache/H-H1_NINJA2_G1000176_EARLY_RECOLORED_CACHE-900000024-10653.lcf --channel-name H1:LDAS-STRAIN --max-total-mass 5.5 --min-total-mass 4.5

--------------------------
Command line options
--------------------------

The command line options read as follows

.. command-output:: pycbc_geom_nonspinbank --help

Some notes on these options:

* The value of f0 generally doesn't matter (so just use 70 if unsure). However if ethinca is being calculated f0 and f-low must be equal.
* Choose f-upper wisely! If your signals coalesce at <1000Hz, you might want to use a lower value ... although of course in this regime this inspiral-only metric will lose validity.
* A delta-f value of 1/256 will certainly be accurate enough, a larger value will cause the code to run faster.
* Using binary values for f-upper, delta-f and sample-rate is not needed for the metric calculation (which is a direct integral), but will provide a noticeable speed difference if reading data to compute a PSD.

.. _tmpltbank_alignedgeombank:

=============================================================
Aligned-spin geometric placement: pycbc_geom_aligned_bank
=============================================================

---------------------
Overview
---------------------

pycbc_geom_aligned_bank is a code for placing of aligned-spin banks using the
TaylorF2 metric.

This code will begin by generating a dax (a pegasus workflow, from which a DAG
is created), which the user submits, that will
generate the final template bank. This code will not calculate ethinca, as
ethinca is not defined for spinning waveforms. However, if a new coincidence
scheme is developed it could be used in that. This code forms part of the
"uberbank", described later, that is being used in Advanced LIGO searches.

**Ian's recommendation: This code will produce optimal template banks, if the
metric can be trusted (ie. BNS). For NSBH and BBH the merger matters. Here the
stochastic approaches described below could be more appropriate for your
parameter space. It is recommended to try these, or a combination of the two
as in the uberbank workflow. Note that the merger matters less in
ZDHP than in O1, so this approach will work better for NSBH banks in 2018 than
it does in 2015.**

If any problems are encountered please email ian.harry@ligo.org.

----------------
Background
----------------

This code build upon the methods used to calculate the metric used in
lalapps_tmpltbank. For details of how this code works see:

* Brown et al., Phys.Rev. D86 (2012) 084017
* Harry et al., Phys.Rev. D89 (2014) 024010

As with our non-spinning generator we place a lattice in the flat xi_i
coordinate system and then map the points back to physical coordinates. In
doing so we remove points that are too far away from the physical space, and
push points that are outside the space, but overlap with it, back into the
physical space.

This is all done with monte-carlo techniques, which is why this code requires
a dag. If there were a more efficient way to map xi_i coordinates to physical
values it will **massively** speed this code up. Please contact ian.harry@ligo.org if you want to play around with this.


--------------------
Some examples
--------------------

The command line arguments given to this code are very similar to those given
to the non-spinning code (not surprising as they actually use the same modules).

Here is one example for making an aligned-spin template bank using a pre-generated ASD file.

.. code-block:: bash

    pycbc_geom_aligned_bank --pn-order threePointFivePN --f0 70 --f-low 15 --f-upper 1000 --delta-f 0.01 --min-match 0.97 --min-mass1 2.5 --min-mass2 2.5 --max-mass1 3. --max-mass2 3. --verbose --max-ns-spin-mag 0.05 --max-bh-spin-mag 0.05 --output-file testAligned.xml --split-bank-num 5 --asd-file /home/spxiwh/aLIGO/BBH_template_banks/asd-T1200307v4.txt --intermediate-data-file intermediate.hdf --metadata-file metadata.xml

and then submitted with something like:

.. code-block:: bash

    pycbc_submit_dax --dax bank_gen.dax --accounting-group ligo.dev.o2.cbc.bbh.pycbcoffline

run

.. code-block:: bash

    pycbc_submit_dax --help

for more details on options for submitting an abstract workflow to condor.

--------------------------
Command line options
--------------------------

The command line options read as follows

.. command-output:: pycbc_geom_aligned_bank --help

As before some notes on these options:

Some notes on these options:

* The value of f0 generally doesn't matter (so just use 70 if unsure). However if ethinca is being calculated f0 and f-low must be equal.
* Choose f-upper wisely! If your signals coalesce at <1000Hz, you might want to use a lower value ... although of course in this regime this inspiral-only metric will lose validity.
* A delta-f value of 1/256 will certainly be accurate enough, a larger value will cause the code to run faster.
* Using binary values for f-upper, delta-f and sample-rate is not needed for the metric calculation (which is a direct integral), but will provide a noticeable speed difference if reading data to compute a PSD.
* For the purposes of spin calculations a NS is considered to be anything with mass < 3 solar masses and a BH is considered to be anything with mass > 3 solar masses. But this can be overridden with the `--ns-bh-boundary-mass` option.
* To avoid generating a NSBH bank where you have a wall of triggers with NS mass = 3 and spins up to the black hole maximum use the nsbh-flag. This will ensure that only neutron-star--black-hole systems are generated.
* If a remnant-mass-threshold is specified the code sets up the bank so that it does not target NS-BH systems that cannot produce a remnant disk mass that exceeds the threshold: this is used to remove EM dim NS-BH binaries from the target parameter space, as discussed in Pannarale and Ohme, ApJL 7, 5 (2014).  Please contact francesco.pannarale@ligo.org if you want to know more about this.

.. _tmpltbank_alignedstochbank:

===============================================================================
Aligned-spin stochastic placement
===============================================================================

lalapps_cbc_sbank and pycbc_aligned_stoch_bank

------------------
Introduction
------------------

pycbc_aligned_stoch_bank is an aligned-spin bank generator that uses a
stochastic algorithm and the TaylorF2 (or TaylorR2F4) metric to create a bank
of templates.

lalapps_cbc_sbank is an aligned-spin bank generator that uses a stochastic algorithm to create a bank of templates. It can use a metric to approximate distances (currently available metrics: TaylorF2 single spin approximation). It can also compute banks using direct matches between waveforms instead of a metric. It can also be used to generate precessing banks (and of course non-spin banks), although computational cost of precessing banks makes this difficult; here the workflow code is necessary and one must spend some time constructing a good configuration file.

Comparisons, advantages and disadvantages of the two codes are discussed below.

When using a metric these codes do not require the processing power of the
geometric lattice algorithm and do not produce a dag.
Run the command, wait, get bank. When using the direct match method in sbank
computational cost is much higher. For this purpose a workflow generator is
also supplied to parallelize generation.

**NOTE: The question "Do I need to use the sbank workflow or just call sbank
directly?" is an important one. As a rough rule of thumb if your template bank
will result in less than 20000 templates (assuming quick-to-generate waveforms)
and you use the optimization options, then a single sbank job should be
sufficient. If expecting more than this then expect to need the workflow
generator. Many of the optimizations used to speed up sbank do not work with
time-domain waveforms, which are often much slower to generate than
frequency-domain ones. sbank can run with time-domain waveforms, but it will
generally be too slow to generate anything but the smallest template banks
(approximately 1000 templates). We strongly recommend not using time-domain
waveforms for bank generation!**

------------------
Background
------------------

As opposed to geometric algorithms where the metric is used to determine a
lattice of equally placed points, a stochastic algorithm works by creating
a large set of randomly placed points and then using the metric or actually
computing overlaps to remove points that are too close to each other.

This algorithm has the benefit that it can work in any parameter space, you
do not require a flat metric. However, it does require more templates to cover
a space than a geometric approach. Additionally the computational cost will
increase as a factor of the number of points in the bank to a power between 2 
and 3 - the exact number depends on how well you are able to optimize the 
parameter space and not match **every** point in the bank with every seed point.
Nevertheless it has been found that the computational cost of generating
aligned-spin banks for Advanced LIGO, using direct-match algorithms is not
too large, and this approach forms a large part of the "uberbank" construction
described below.

The stochastic bank code can calculate the ethinca metric, but **only** if both 
component maximum spins are set to 0.0, i.e. a stochastic non-spinning bank. One can also generate only the time component of the metric if doing exact-match coincidence.

Stochastic placement was first proposed in terms of LISA searches in

* Harry et al. Class.Quant.Grav. 25 (2008) 184027
* Babak Class.Quant.Grav. 25 (2008) 195011

(the former using a metric, and the latter using direct match). Then described in more detail in

* Harry et al. Phys.Rev. D80 (2009) 104014 
* Manca and Vallisneri Phys.Rev. D81 (2010) 024004

Recently stochastic placement has been explored for LIGO searches in

* Ajith et al., Phys.Rev. D89 (2014) 084041
* Privitera et al., Phys.Rev. D89 (2014) 024003
* Harry et al., Phys.Rev. D89 (2014) 024010
* Capano et al., Phys.Rev. D93 (2016) 124007

pycbc_aligned_stoch_bank follows the method described in Harry et al. (2009)
and calculates matches using a metric (in this case the F2 metric).
lalapps_cbc_sbank can do stochastic template bank generation with or
without a metric. In the absence of a metric it uses the method introduced in
Babak (2008) and used in Privitera (2013) to compute distances
between points by generating both waveforms and calculating an explicit overlap.

-------------------------------------------------------------
Ian's recommendation: Which stochastic code should I use?
-------------------------------------------------------------

Okay so there are two stochastic codes, and a lot of overlap between the two
codes. We are aware that this is not optimal and hope to combine the features
of these two codes into a single front-end to get the best of both worlds.

In the mean-time here is how I see the breakdown of the two codes:

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
I want to use the F2 metric
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

If you want to use the F2 metric used in the geometric code
then I recommend to use pycbc_aligned_stoch_bank. Here I have found the
code to be faster
than lalapps_cbc_sbank using the F2 reduced spin metric because the pycbc
code is using
a Cartesian parameter space and is therefore able to greatly reduce the
number of matches that are calculated. (Effectively lalapps_cbc_sbank only
limits the number of matches calculated by chirp mass differences,
pycbc_aligned_stoch_bank also uses the second direction (some combination
of mass ratio and the spins).
The pycbc metric also incorporates the effect of both spins.

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
I want to use another metric, ie. an IMR metric
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

If you have and want to use a different metric, such as the IMRPhenomX
metric, where the approach used in the geometric bank cannot be used to
create a Cartesian coordinate system, then use lalapps_cbc_sbank.

This code is more flexible, and is not dependent on the assumptions that are
used in the geometric code. If your metric is not already added then please
contact a developer for help in integrating this.

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
I want to use direct match
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

Use lalapps_cbc_sbank. Be aware though that this is often difficult to run in a
single process and may require the dag generator described below.

--------------------
Some examples
--------------------

The command line arguments given to this code are very similar to those given
to the non-spinning code and the aligned-spin geometric code (they use the same modules).

Here is one example reading data from a frame cache, as in the non-spinning example

.. code-block:: bash

    pycbc_aligned_stoch_bank --pn-order threePointFivePN --f0 60 --f-low 30 -V --delta-f 0.01 --min-match 0.97 --min-mass1 2.5 --max-mass1 3 --min-mass2 2.5 --max-mass2 3 --max-ns-spin-mag 0.05 --max-bh-spin-mag 0.05 --nsbh-flag --psd-estimation median --psd-segment-length 256 --psd-segment-stride 128 --psd-inverse-length 8 --gps-start-time 900000033 --gps-end-time 900002081 --strain-high-pass 30 --pad-data 8 --sample-rate 4096 --frame-cache cache/H-H1_NINJA2_G1000176_EARLY_RECOLORED_CACHE-900000024-10653.lcf --channel-name H1:LDAS-STRAIN --verbose --output-file testStoch.xml --num-seeds 2000000 --f-upper 2000

Another example where we read the PSD from a supplied ASD file:

.. code-block:: bash

    pycbc_aligned_stoch_bank --pn-order threePointFivePN --f0 60 --f-low 30 --delta-f 0.1 --min-match 0.97 --min-mass1 2.5 --max-mass1 3 --min-mass2 2.5 --max-mass2 3 --max-ns-spin-mag 0.05 --max-bh-spin-mag 0.05 --nsbh-flag --verbose --asd-file ZERO_DET_high_P.txt --num-seeds 10000 --output-file "testStoch.xml" --f-upper 2000

And a third example where we use a batshit PN order.

.. code-block:: bash

    pycbc_aligned_stoch_bank --pn-order onePN --f0 60 --f-low 30 --delta-f 0.1 --min-match 0.97 --min-mass1 2.5 --max-mass1 3 --min-mass2 2.5 --max-mass2 3 --max-ns-spin-mag 0.05 --max-bh-spin-mag 0.05 --nsbh-flag --verbose --asd-file ZERO_DET_high_P.txt --num-seeds 10000 --output-file "testStoch.xml" --f-upper 2000

lalapps_cbc_sbank provides example in the help text, see below.

-----------------------------------------------
Command line options: pycbc_aligned_stoch_bank
-----------------------------------------------

The command line options read as follows

.. command-output:: pycbc_aligned_stoch_bank --help

As before some notes on these options:

Some notes on these options:

* The value of f0 generally doesn't matter (so just use 70 if unsure). However if ethinca is being calculated f0 and f-low must be equal.
* Choose f-upper wisely! If your signals coalesce at <1000Hz, you might want to use a lower value ... although of course in this regime this inspiral-only metric will lose validity.
* A delta-f value of 1/256 will certainly be accurate enough, a larger value will cause the code to run faster.
* For the purposes of spin calculations a NS is considered to be anything with mass < 3 solar masses and a BH is considered to be anything with mass > 3 solar masses. But this can be overridden with the `--ns-bh-boundary-mass` option.
* To avoid generating a NSBH bank where you have a wall of triggers with NS mass = 3 and spins up to the black hole maximum use the nsbh-flag. This will ensure that only neutron-star--black-hole systems are generated.
* --num-seeds is the termination condition. The code will throw NUM_SEEDS points at the parameter space, and then filter out those that are too close to each other. If this value is too low, the bank will not converge, if it is too high, the code will take longer to run.
* --num-failed-cutoff can be used as an alternative termination condition. Here the code runs until NUM_FAILED_CUTOFF points have been **consecutively** rejected and will then stop.
* --vary-fupper will allow the code to vary the upper frequency cutoff across the parameter space. Currently this feature is **only** available in the stochastic code. 
* If a remnant-mass-threshold is specified the code sets up the bank so that it does not target NS-BH systems that cannot produce a remnant disk mass that exceeds the threshold: this is used to remove EM dim NS-BH binaries from the target parameter space, as discussed in Pannarale and Ohme, ApJL 7, 5 (2014).  Please contact francesco.pannarale@ligo.org if you want to know more about this.

-----------------------------------------------
Command line options: lalapps_cbc_sbank
-----------------------------------------------

The program lalapps_cbc_sbank is part of lalsuite and documentation and command line options are described in the documentation at http://software.ligo.org/docs/lalsuite/lalapps/lalapps__cbc__sbank_8py_source.html

Some notes on the command line options:

* Using the `--iterative-match-df-max 8` option will greatly speed up the code by generating waveforms with a large frequency step and then iteratively lowering it until the calculated matches converge.
* Using the `--cache-waveforms` option will greatly speed up the code by avoiding generating a waveform more than once. With this enabled memory usage can become a concern. If memory usage becomes high consider using the workflow generator.

.. _tmpltbank_sbankworkflow:

------------------------------------------------
Sbank workflow generator
------------------------------------------------

In the case where one sbank job will not do the job you can use
pycbc_create_sbank_workflow to parallelize your template bank placement
needs.

The command line options for this code read as follows

.. command-output:: pycbc_create_sbank_workflow --help

The main thing here is that a configuration file is supplied. This
configuration file supplies basically all options about the template bank
generation. One can read about the general setup of such files
:ref:`at the following page <workflowconfigparsermod>`, but we provide
some examples here. First we give a simple example, with some options commented
out, but still shown and explained.

.. literalinclude:: resources/sbank_example1.ini
   :language: ini

Then a more complex example, where we exercise more fine-grained control over
the parallel job and use different settings during each cycle of the generator.

.. literalinclude:: resources/sbank_example2.ini
   :language: ini

**PLEASE NOTE: These examples are intended to be illustrative of the workflow
generator and not intended to be used verbatim, expecting this to work for all
examples. You will want to tune these options based around whatever examples
you are running. When doing that you want to balance job run time, memory
usage, wall-clock time, avoiding too many templates at the readder stage.
Please ask for help if you have a specific example and want some help to
optimize.**

.. _tmpltbank_uberbankworkflow:

===================================================
Hybrid approaches: the best of both worlds
===================================================

We have found that in many cases it makes sense to combine the techniques of
sbank with the geometric lattic algorithm to produce a template bank. This
technique was used in Advanced LIGO's first observing run to produce what was
called the "uberbank". The references for this are the following:

* Abbott et al., Phys.Rev. D93 (2016) 122003
* Capano et al., Phys.Rev. D93 (2016) 124007

The uberbank construction process is now written up in a single workflow,
pycbc_create_uberbank_workflow. The current setup of this workflow is as
follows:

* Run in parallel one sbank workflow and one geometric workflow
* Combine the two outputs together (this assumes the two do not overlap)
* Run a further sbank workflow to fill in the remaining space and any holes

The idea to this 3-part workflow is it allows us to first focus on covering
well both the BNS and BBH regions of parameter space, before filling in the
rest with potentially different settings. There is a lot of scope for editing
and reordering this workflow, potentially adding, or removing various of the
stages.

The command-line help for the workflow generator is as follows:

.. command-output:: pycbc_create_uberbank_workflow --help

As with the sbank workflow generator the main bulk of the configuration is the
configuration file, which is provided on the command line. This configuration
file contains options for all 3 stages of the workflow, and so is a little more
involved than the sbank example given above. Here we provide a particularly
detailed configuration for generating a bank equivalent to the uberbank on
O1-like data.

.. literalinclude:: resources/uberbank_example1.ini
   :language: ini

==========================
The module's source code
==========================

Follow the following link to view the documentation (and through that the source code) of the pycbc.tmpltbank module:

.. toctree::
    pycbc.tmpltbank

The code also calls into the PSD module whose documentation is here:

.. toctree::
    pycbc.psd

and the data generation/reading routines, which are here:

.. toctree::
    pycbc.noise

