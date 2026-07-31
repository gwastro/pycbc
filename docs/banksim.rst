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

--------------------------------------------
Running a banksim workflow (recommended)
--------------------------------------------

For anything beyond a single, quick check, the recommended way to run a
banksim is via the ``pycbc_make_bank_verifier_workflow`` workflow
generator, which splits the work across a Condor pool, combines the
results, and produces summary plots and an HTML results page
automatically. See :ref:`bank_verifier` for full details on configuring
and running it.

.. note::

    ``pycbc_make_bank_verifier_workflow`` replaces the older
    ``pycbc_make_banksim`` script, which has been removed. Configuration
    files written for the old script use different section and option
    names and are not directly compatible; see :ref:`bank_verifier` for
    the current format.

----------------------------------
Running pycbc_banksim directly
----------------------------------

For a single, ad-hoc comparison (e.g. testing a small number of templates
against a handful of injections), ``pycbc_banksim`` can be run directly
without going through a workflow.

.. literalinclude:: ../examples/banksim/run.sh
   :language: bash

The ``--template-file`` and ``--signal-file`` options each accept either a
LIGOLW XML file (a ``sngl_inspiral`` table for templates, a
``sim_inspiral``/``sngl_inspiral`` table for signals) or an HDF file, as
read by ``pycbc.waveform.bank.TemplateBank``. See ``pycbc_banksim --help``
for the full set of options.

The main result of running ``pycbc_banksim`` is a single, whitespace
separated ASCII match file, with one row per injection giving the maximum
match found in the bank, the bank file and template index that achieved
it, and the signal's normalization.

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

The bank generation can be verified using the ``pycbc_banksim`` code, run
via a ``pycbc_make_bank_verifier_workflow`` workflow (see
:ref:`bank_verifier`). Each of the BNS/NSBH/BBH mass regions above can be
configured as its own named broad-injection set (a ``[workflow-broadinjs]``
entry), with its own ``[injection-TAG]`` section giving that region's mass
and spin ranges, all pointed at the same ``input-bank``.

To run this you will need to point ``[executables] banksim`` (and the
other executables) at your local install, set an ``accounting-group``
suitable for your cluster in ``[pegasus_profile]``, and give the location
of the bank and noise curve in ``[banksim]``.

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
