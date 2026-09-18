###################################################
Waveforms
###################################################

=====================================
What waveforms can I generate?
=====================================

.. literalinclude:: ../examples/waveform/what_waveform.py
.. command-output:: python ../examples/waveform/what_waveform.py


=====================================
Plotting Time Domain Waveforms 
=====================================

.. plot:: ../examples/waveform/plot_waveform.py
   :include-source:

==============================================
Generating one waveform in multiple detectors
==============================================

.. plot:: ../examples/waveform/plot_detwaveform.py
   :include-source:


===============================================
Selecting which modes to include
===============================================
Gravitational waves can be decomposed into a set
of modes. Some approximants only calculate the dominant
l=2, m=2 mode, while others included higher-order modes. These
often, but not always, include 'HM' in the name. The modes
present in the output polarizations can be selected for these waveforms
as demonstrated below. By default, all modes that the waveform model
supports are typically generated.

.. plot:: ../examples/waveform/higher_modes.py
   :include-source:

=======================================
Calculating the match between waveforms
=======================================

.. literalinclude:: ../examples/waveform/match_waveform.py
.. command-output:: python ../examples/waveform/match_waveform.py

================================================
Plotting a TD and FD waveform together in the TD
================================================

.. plot:: ../examples/waveform/plot_fd_td.py
   :include-source:
   
================================================
Plotting GW phase and amplitude of TD waveform
================================================

.. plot:: ../examples/waveform/plot_phase.py
   :include-source:

================================================
Plotting frequency evolution of TD waveform
================================================

.. plot:: ../examples/waveform/plot_freq.py
   :include-source:
   
=====================================
Adding a custom waveform
=====================================

You can also add your own custom waveform and make it available through 
the waveform interface standard. You can directly call the code as below
or if you include it in a python package, :ref:`PyCBC can directly detect it! <waveform_plugin>`

.. plot:: ../examples/waveform/add_waveform.py
   :include-source:
   
