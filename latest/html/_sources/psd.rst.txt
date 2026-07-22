###################################################
Handling PSDs
###################################################

=====================================
Reading / Saving a PSD from a file
=====================================

If you have a PSD in a two column, space separated, (frequency strain), you can
load this into PyCBC.

.. plot:: ../examples/psd/read.py
   :include-source:

.. _Analytic PSDs from lalsimulation:

==============================================
Generating an Analytic PSD from lalsimulation
==============================================

A certain number of PSDs are built into lalsimulation, which you'll be able
to access through PyCBC. Below we show how to see which ones are available, 
and demonstrate how to generate one.

.. plot:: ../examples/psd/analytic.py
   :include-source:

The PSDs from lalsimulation are computed at the required frequencies by
interpolating a fixed set of samples; if the required frequencies fall
outside of the range of the known samples no warnings will be raised,
and (meaningless) extrapolated values will be returned.

Therefore, users should check the validity range of the PSD they are 
using within lalsimulation, lest they get incorrect results.

====================================
Estimating the PSD of a time series
====================================

.. plot:: ../examples/psd/estimate.py
   :include-source:

