.. _inference_example_lisa_smbhb_inj:

---------------------------------------------
LISA SMBHB injection and parameter estimation
---------------------------------------------

This example shows how to use PyCBC for time-domain LISA TDI noise generation and 
supermassive black hole binaries (SMBHB) signal injection. This one is similar to 
:ref:`LISA parameter estimation for simulated SMBHB from LDC example
<inference_example_lisa_smbhb_ldc>`, the main difference is we generate our own mock data 
in this example. In order to to that, we use 
`LISA TDI PSD module <https://pycbc.org/pycbc/latest/html/pycbc.psd.html#module-pycbc.psd.analytical_space>`_ 
to generate the stationary and Gaussian noise for each TDI channel in the time domain, then we use 
`waveform injection module <https://pycbc.org/pycbc/latest/html/_modules/pycbc/waveform/waveform.html#get_td_det_waveform_from_fd_det>`_ 
to add the simulated signal into the simulated noise.

First, we use the following configuration file to define the parameters of our SMBHB injection, we use the 
same parameters from the SMBHB signal in :ref:`LISA parameter estimation for simulated SMBHB from LDC example
<inference_example_lisa_smbhb_ldc>`:

.. literalinclude:: ../../../examples/inference/lisa_smbhb_inj/injection_smbhb.ini
   :language: ini

:download:`Download <../../../examples/inference/lisa_smbhb_inj/injection_smbhb.ini>`

Then we run the following bash script to create a .hdf file that contains same information:

.. literalinclude:: ../../../examples/inference/lisa_smbhb_inj/injection_smbhb.sh
   :language: bash

:download:`Download <../../../examples/inference/lisa_smbhb_inj/injection_smbhb.sh>`

Here, we use a similar configuration file for parameter estimation, we also use 
:py:class:`Relative <pycbc.inference.models.relbin.Relative>` model. We also just 
set chirp mass, mass ratio and tc as variable parameters, `tc`, `eclipticlongitude`, `eclipticlatitude` 
and `polarization` are defined in the LISA frame:

.. literalinclude:: ../../../examples/inference/lisa_smbhb_inj/lisa_smbhb_relbin.ini
   :language: ini

:download:`Download <../../../examples/inference/lisa_smbhb_inj/lisa_smbhb_relbin.ini>`


Now run:

.. literalinclude:: ../../../examples/inference/lisa_smbhb_inj/run.sh
   :language: bash

:download:`Download <../../../examples/inference/lisa_smbhb_inj/run.sh>`

To plot the posterior distribution after the last iteration, you can run the following script:

.. literalinclude:: ../../../examples/inference/lisa_smbhb_inj/plot.sh
   :language: bash

:download:`Download <../../../examples/inference/lisa_smbhb_inj/plot.sh>`

In this example it will create the following plot:

.. image:: ../../_include/lisa_smbhb_mass_tc.png
   :scale: 60
   :align: center

The scatter points show each walker's position after the last iteration. The
points are colored by the SNR at that point, with the 50th and 90th
percentile contours drawn. The red lines represent the true parameters of injected signal.
