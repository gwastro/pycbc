###################################################################
Verifying the construction of banks with pycbc_banksim
###################################################################

===================
Introduction
===================

This page describes how to use the banksim facilities within PyCBC to test a previously built template bank. A banksim calculates the matches, maximized over a set of templates, for a list of injections waveforms. 

The purpose of this test is to act as a sanity check of the template bank generation methodology and code. Therefore the tests run over the same parameter ranges used to generate the bank and using the same sensitivity curve. The tests described here may not be optimal or exhaustive, but should be used to gain confidence that a bank is recovering signals it is designed to recover at an acceptable match.

===========================================
Testing the uberbank
===========================================

To test the PyCBC uberbank that contains templates for BNS, NSBH and BBH in one template bank, we propose to test separately the BNS, NSBH and BBH regions, with separate banksim runs. Therefore there will be some overlap between the signals tested. For technical reasons, it is also convenient to split the NSBH and BBH tests up into a run with signals below a total mass of 50 and signals with a total mass above 50. 

We propose to select test signals from mass distributions that a flat in component masses in the respective regions; NS masses between 1 and 3 and BH masses between 2 and 99, with a total mass limit of 100. In addition, we select aligned spin magnitudes uniform in the respective regions; -0.05 to 0.05 for NS and -0.99 to 0.99 for BH. 

We propose to test with 10,000 injection signals in each of the BNS, NSBH and BBH regions, for a total of 30,000. This number is much less than the total number of templates in the bank. 

We propose to use SEOBNRv2 as the test signals, even though the uberbank uses TaylorF2 and SEOBNRv2_ROM_DoubleSpin templates for recovery. This is because we believe that SEOBNRv2 is a more accurate waveform than either TaylorF2 or the ROMs.

================================================
Procedure
================================================

The bank generation can be verified using the pycbc_banksim code. To run this follow the instructions for running the banksim code. An example ini file to run the NSBH banksim for total masses below 50 is given in examples/banksim/nsbh_below50.ini

.. literalinclude:: ../examples/banksim/nsbh_below50.ini

To run this you will need to change the banksim option to your local version of pycbc_banksim, the log-path option to a suitable location for your log files on your cluster, the locations of the bank and noise curve and possibly whatever processing_scheme is best on your cluster (mkl works on Atlas with /opt/intel/2015/intel.sh sourced).

Injected spins are up to 0.99, not 0.9895 and the injections are uniform in component mass from 1 to 50 and uniform in spin magnitude (so it contains some highly spinning BNS). Injections are generated from 25Hz but matches are calculated from 30Hz, this gives the signal some "burn-in" time. Source location l-distr is random over the sky and inclination i-distr is uniformly distributed over arccos(i) - although this should not matter for aligned signals.

========================================
Evaluation
========================================

A stochastic placement method (like sbank) will not be able to guarantee that all points in parameter space are covered at better than 0.97 fitting factor. A convenient measure of the success of the bank generation is if the bank is able to recover 99% of injected signals using the same parameters and templates as the bank is designed for with a fitting factor of 0.97 or better. Further requirements might be that there should be no fitting factors with matches less than 0.95 or that the fitting factors below 0.97 should not be clustered in a particular part of parameter space. To cover all source groups we can run such tests separately for simulated BNS, NSBH and BBH signals when testing a bank that covers all three parameter ranges.

While such tests do not guarantee that the bank will successfully recover all possible signals in the parameter region (for example due to non-Gaussian noise, different sensitivites in the two detectors, different waveform approximants, precession effects, tidal deformation and disruption etc.) these tests do indicate with a reasonable level of confidence that the template generation has been successful at what it was designed to do.

=================================================
Known issues
=================================================

The coverage of the high-mass (>70) and anti-aligned (<-0.5) NSBH region is known to be sparse in some versions when using SEOBNRv2 signals and SEOBNRv2_ROM_DoubleSpin templates.

The mchirp-window size may need to be changed if it is too tight. This is particularly a problem at higher masses.

If speed is an issue, the banksims can be sped up by reducing the number of injection signals, using ROMs instead of SEOBNRv2 as injection signals, reducing the signal-sample-rate or tightening the mchirp-window. Code is being developed to do this dynamically.

To replicate the behaviour of the uberbank, that switches the template approximant from TaylorF2 to SEOBNRv2_ROM_DoubleSpin at a total mass of 4, the options total-mass-divide and highmass-approximant are needed in the pycbc_banksim code. These options do not exist on older versions (Sep 2015) of pycbc_banksim.

