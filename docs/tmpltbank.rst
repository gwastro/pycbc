.. _tmpltbankmod:

###################################################################
PyCBC template bank generation documentation (``pycbc.tmpltbank``)
###################################################################

===================
Introduction
===================

This page gives details on how to use the various bank generation codes
available in pyCBC, the pyCBC.tmpltbank module. In addition we describe the
sBank code. These codes currently consist of 

* :ref:`A program to generate a non-spinning template bank (pycbc_geom_nonspinbank) <tmpltbank_nonspinbank>`
* :ref:`A program to generate an aligned-spin template bank using a geometrical lattice algorithm (pycbc_geom_aligned_bank) <tmpltbank_alignedgeombank>`
* :ref:`Two programs to generate an aligned-spin template bank using a stochastic placement algorithm. This is done either using a metric to approximate distances, or by computing explicit matches from waveforms. (pycbc_aligned_stoch_bank, lalapps_cbc_sbank) <tmpltbank_alignedstochbank>`

Each of these codes is described in turn and is accompanied by some examples of how to run the code, some relevant background and references and the options described from the codes' help messages.

The output of all these codes is currently a sngl_inspiral_table XML object.
However, if people want to output in SQL or other formats it should be
easy to write a back-end module to do this. The code may not fill every
necessary column of the sngl_inspiral table, please let me know if something
is missing and it will be added!

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


**Ian's recommendation: Use this code instead of lalapps_tmpltbank if you want
a non-spinning bank.**

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

.. _tmpltbank_alignedgeombank:

=============================================================
Aligned-spin geometric placement: pycbc_geom_aligned_bank
=============================================================

---------------------
Overview
---------------------

pycbc_geom_aligned_bank is a code for placing of aligned-spin banks using the
TaylorF2 metric (or TaylorR2F4, which is an approximation of TaylorT4).

This code will begin by generating a dag, which the user submits, that will
generate the final template bank. This code will not calculate ethinca, as
ethinca is not defined for spinning waveforms. However, if a new coincidence
scheme is developed it could be used in that. This code has already been used
in two published works, has been run through the GRB code, MBTA and in
ahope development work ongoing at AEI.

**Ian's recommendation: This code will produce optimal template banks, if the
metric can be trusted (ie. BNS). For NSBH and BBH the merger matters. Here the
stochastic approaches described below could be more appropriate for your
parameter space. It is recommended to try these, and maybe even look into
the hybrid approaches described below. Note that the merger matters less in
ZDHP than in O1, so this approach works better for NSBH banks in 2018 than it
will in 2015.**

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

    pycbc_geom_aligned_bank --pn-order threePointFivePN --f0 70 --f-low 15 --f-upper 1000 --delta-f 0.01 --min-match 0.97 --min-mass1 2.5 --min-mass2 2.5 --max-mass1 3. --max-mass2 3. --verbose --log-path /usr1/spxiwh/logs --max-ns-spin-mag 0.05 --max-bh-spin-mag 0.05 --output-file testAligned.xml --split-bank-num 5 --asd-file ZERO_DET_high_P.txt

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
* For the purposes of spin calculations a NS is considered to be anything with mass < 3 solar masses
* And a BH is considered to be anything with mass > 3 solar masses.
* To avoid generating a NSBH bank where you have a wall of triggers with NS mass = 3 and spins up to the black hole maximum use the nsbh-flag. This will ensure that only neutron-star--black-hole systems are generated.

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

lalapps_cbc_sbank is an aligned-spin bank generator that uses a stochastic algorithm to create a bank of templates. It can use a metric to approximate distances (currently available metrics: TaylorF2 single spin approximation). It can also compute banks using direct matches between waveforms instead of a metric. In theory this could also be used to generate precessing banks (and of course non-spin banks), although computational cost of precessing banks may be way too much to make this practical.

Comparisons, advantages and disadvantages of the two codes are discussed below.

When using a metric these codes do not require the processing power of the
geometric lattice algorithm and does no produce a dag.
Run the command, wait, get bank. When using the direct match method in sbank
computational cost is much higher. For this purpose a dag generator is also
supplied to parallelize generation.

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
Nevertheless it has been found that the computational cost is not prohibitive
for aLIGO-size stochastically-generated spinning banks, when a metric can be
used. 

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

pycbc_aligned_stoch_bank follows the method described in Harry et al. (2009)
and calculates matches using a metric (in this case the F2 or R2F4 metrics).
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
I want to use the F2 or R2F4 metric
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

If you want to use one of the metrics that can be used in the geometric code
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

Use lalapps_cbc_sbank. Be aware though that this is difficult to run in a
single process and will require the dag generator described below. However,
this will produce banks that have significantly overlapping coverage regions
as the parallelized jobs cannot talk to each other. Choose the number of
parallel jobs wisely! It is possible to have 2x-3x more templates than would be
generated from running the code in a single process.

We are also working on a hybrid method where the metric approach is used to
create a first-pass of a template bank and then the direct match is used to
fill the holes. This has the advantage of significantly reducing the
overcoverage from using the dag generator.

--------------------
Some examples
--------------------

The command line arguments given to this code are very similar to those given
to the non-spinning code and the aligned-spin geometric code (they use the same modules).

Here is one example reading data from a frame cache, as in the non-spinning example

.. code-block:: bash

    pycbc_aligned_stoch_bank -v --pn-order threePointFivePN --f0 60 --f-low 30 -V --delta-f 0.01 --min-match 0.97 --min-mass1 2.5 --max-mass1 3 --min-mass2 2.5 --max-mass2 3 --max-ns-spin-mag 0.05 --max-bh-spin-mag 0.05 --nsbh-flag --psd-estimation median --psd-segment-length 256 --psd-segment-stride 128 --psd-inverse-length 8 --gps-start-time 900000033 --gps-end-time 900002081 --strain-high-pass 30 --pad-data 8 --sample-rate 4096 --frame-cache cache/H-H1_NINJA2_G1000176_EARLY_RECOLORED_CACHE-900000024-10653.lcf --channel-name H1:LDAS-STRAIN --verbose --output-file testStoch.xml --num-seeds 2000000

Another example where we read the PSD from a supplied ASD file:

.. code-block:: bash

    pycbc_aligned_stoch_bank -v --pn-order threePointFivePN --f0 60 --f-low 30 --delta-f 0.1 --min-match 0.97 --min-mass1 2.5 --max-mass1 3 --min-mass2 2.5 --max-mass2 3 --max-ns-spin-mag 0.05 --max-bh-spin-mag 0.05 --nsbh-flag --verbose --asd-file ZERO_DET_high_P.txt --num-seeds 10000 --output-file "testStoch.xml"

And a third example where we use a batshit PN order.

.. code-block:: bash

    pycbc_aligned_stoch_bank -v --pn-order onePN --f0 60 --f-low 30 --delta-f 0.1 --min-match 0.97 --min-mass1 2.5 --max-mass1 3 --min-mass2 2.5 --max-mass2 3 --max-ns-spin-mag 0.05 --max-bh-spin-mag 0.05 --nsbh-flag --verbose --asd-file ZERO_DET_high_P.txt --num-seeds 10000 --output-file "testStoch.xml"

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
* For the purposes of spin calculations a NS is considered to be anything with mass < 3 solar masses
* And a BH is considered to be anything with mass > 3 solar masses.
* To avoid generating a NSBH bank where you have a wall of triggers with NS mass = 3 and spins up to the black hole maximum use the nsbh-flag. This will ensure that only neutron-star--black-hole systems are generated.
* --num-seeds is the termination condition. The code will throw NUM_SEEDS points at the parameter space, and then filter out those that are too close to each other. If this value is too low, the bank will not converge, if it is too high, the code will take longer to run.
* --num-failed-cutoff can be used as an alternative termination condition. Here the code runs until NUM_FAILED_CUTOFF points have been **consecutively** rejected and will then stop.
* --vary-fupper will allow the code to vary the upper frequency cutoff across the parameter space. Currently this feature is **only** available in the stochastic code. 

-----------------------------------------------
Command line options: lalapps_cbc_sbank
-----------------------------------------------

The command line options read as follows

.. command-output:: lalapps_cbc_sbank --help

Some notes on these options:

* Anything to highlight Steve?

-----------------------------------------------
Command line options: lalapps_cbc_sbank_pipe
-----------------------------------------------

The command line options read as follows

.. command-output:: lalapps_cbc_sbank_pipe --help

Some notes on these options:

* Anything to highlight Steve?

===================================================
Hybrid approaches: the best of both worlds
===================================================

We are currently looking into hybrid bank construction techniques. This is
where one of the methods used above is used to create a preliminary bank and
then a second method is run on that output to fill in any holes that might 
exist in parameter space.

Documentation on hybrid approaches will be added soon.

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

