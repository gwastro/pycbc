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

==============================================
Generating an Analytic PSD from lalsimulation
==============================================

A certain number of PSDs are built into lalsimulation, which you'll be able
to access through PyCBC. Below we show how to see which ones are available, 
and demonstrate how to generate one.

.. plot:: ../examples/psd/analytic.py
   :include-source:

====================================
Estimating the PSD of a time series
====================================

.. plot:: ../examples/psd/estimate.py
   :include-source:

