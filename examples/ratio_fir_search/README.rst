#################################################
Hierarchical FIR/Ratio Matched Filtering Example
#################################################

.. contents::

=================================================
Introduction
=================================================

This example demonstrates the hierarchical "FIR/ratio" matched-filtering
method (``pycbc_fir_bank`` and ``pycbc_inspiral_fir``), a way to search a
densely-spaced ("fine") template bank without matched-filtering every fine
template directly against the data.

The idea: build a much sparser "coarse" bank, matched-filter only the
coarse templates, and for every fine template design a short real FIR
filter that approximates the frequency-domain ratio between it and its
best-matching coarse template. At search time, applying that filter to the
coarse template's SNR time series reconstructs an approximation of the
fine template's own SNR, without ever running a full matched filter for
it. This can be significantly cheaper than filtering every fine template
directly, when many fine templates share a similar frequency-domain shape.

This example is self-contained (it uses simulated Gaussian noise, not real
strain data) and runs in about a minute. The example directory is::

  examples/ratio_fir_search

=================================================
Running the example
=================================================

Run everything in order with::

  bash run_example.sh

which runs, in order:

#. ``make_banks.sh`` -- builds a loose "coarse" bank (``bank_coarse.hdf``)
   and a tight "fine" bank (``bank_full.hdf``) with ``pycbc_brute_bank``.
#. ``make_fir_bank.sh`` -- builds the FIR ratio bank (``fir_full_bank.hdf``)
   from the two banks above with ``pycbc_fir_bank``.
#. ``run_reference_search.sh`` -- a plain ``pycbc_inspiral`` run over the
   fine bank (``reference.hdf``), used only as a check.
#. ``run_fir_search.sh`` -- the FIR/ratio search itself, over the same data
   and bank, with ``pycbc_inspiral_fir`` (``fir.hdf``).
#. ``compare_triggers.py`` -- matches triggers between the two outputs by
   time and template, and reports the standard deviation of their SNR
   difference.

=================================================
Interpreting the comparison
=================================================

``compare_triggers.py`` compares the measured SNR-difference standard
deviation against the value predicted from the FIR bank's ``--min-match``
target (``sqrt(2 * (1 - min_match))``, the standard mismatch-induced SNR
fluctuation bound), read directly from ``fir_full_bank.hdf`` so it always
reflects whatever ``make_fir_bank.sh`` was actually run with.

The measured value is expected to be several times larger than this
mismatch bound. This isn't a bug: the FIR/ratio reconstruction runs in
single precision for performance, and a FIR filter whose frequency
response has a wide dynamic range (which depends on how different the
coarse and fine templates are) loses some accuracy in that single-precision
arithmetic beyond what the bank's own double-precision construction-time
verification predicts. ``compare_triggers.py`` allows a configurable
margin (``--allowed-margin-factor``, default 4x) above the mismatch
prediction for this reason.

=================================================
Tuning the FIR bank
=================================================

``make_fir_bank.sh``'s options that most affect the FIR bank's accuracy
and search-time cost:

``--min-match``
  Target overlap between each fine template and its FIR-filter
  reconstruction. Higher generally means more taps (and so more
  search-time compute) per template.
``--n-taps``, ``--max-taps``, ``--tap-step``
  The tap-count escalation ladder: start at ``--n-taps``, and if
  ``--min-match`` isn't reached, add ``--tap-step`` taps at a time up to
  ``--max-taps``.
``--regularization-type``, ``--regularization-magnitude``
  Controls the least-squares fit's regularization (ridge, derivative, or
  composite), trading fit accuracy for smoother/smaller filter
  coefficients.
``--ratio-cap``
  Frequency bins where the fine/coarse amplitude ratio exceeds this value
  are excluded from the fit, to avoid numerical instability from
  near-zero coarse-template amplitudes.
