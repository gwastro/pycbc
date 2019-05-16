#############################################################
Query times of valid data, hardware injections, and more.
#############################################################

Information about the state of the LIGO/Virgo data can be queried from
resources at the GWOSC and from LVC-proprietary DQSEGDB through the
functions of the :py:mod:`pycbc.dq` module. We will outline a few of common
tasks below.

=================================================
Determine the times an instrument has valid data
=================================================

.. plot:: ../examples/dataquality/on.py
   :include-source:

====================================
Finding times of hardware injections
====================================

.. plot:: ../examples/dataquality/hwinj.py
   :include-source:

========================
What flags can I query?
========================

A list of many of the flags which can be quiered is `available here <https://www.gw-openscience.org/archive/dataset/O1/>`_. Instead, just give the
raw name such as "DATA" instead of "H1_DATA".

There are two additional types of flags which can be queried. These are
the negation of the flags like "NO_CBC_HW_INJ". The flag "CBC_HW_INJ" would
give times where there *is* a hardware injection instead of times when
there isn't. Similarly if you use "CBC_CAT2_VETO" instead of "CBC_CAT2" you
will get the times that are adversely affected instead of the times that
are not.
