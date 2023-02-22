----------------------------------------------
LISA parameter estimation for simulated SMBHB
----------------------------------------------

This example shows how to use PyCBC for parameter estimation of supermassive black hole binaries (SMBHB) 
in LISA mock data. The `data <https://zenodo.org/record/7497853>`_ are generated from 
`LISA Data Challenge 2a: Sangria <https://lisa-ldc.lal.in2p3.fr/challenge2a>`_, 
and `BBHx <https://github.com/mikekatz04/BBHx>`_ package is used to generate the ``IMRPhenomD`` template and calculate 
the corresponding TDI response for LISA. Relative binning (heterodyned likelihood) 
is used during sampling to speed up the computation of likelihood functions. Before doing parameter estimation, 
you need to install `BBHx <https://github.com/mikekatz04/BBHx>`_ and `the corresponding PyCBC waveform plugin <https://github.com/ConWea/BBHX-waveform-model>`_, 
please click the corresponding link to see the detailed description of the installation.

First, we create the following configuration file, here we just set chirp mass, mass ratio and tc as variable parameters, 
`tc`, `eclipticlongitude`, `eclipticlatitude` and `polarization` are defined in the LISA frame:

.. literalinclude:: ../../../examples/inference/lisa_smbhb/lisa_smbhb_relbin.ini
   :language: ini

:download:`Download <../../../examples/inference/lisa_smbhb/lisa_smbhb_relbin.ini>`

By setting the model name to ``relative`` we are using
:py:class:`Relative <pycbc.inference.models.relbin.Relative>` model.

In this simple example, we do the parameter estimation for the first SMBHB signal in the LDC Sangria dataset 
(you can also run parameter estimation for other SMBHB signals by choosing appropriate prior range),
we need download the data first (`MBHB_params_v2_LISA_frame.pkl` contains all the true parameters):

.. literalinclude:: ../../../examples/inference/lisa_smbhb/get.sh
   :language: bash

:download:`Download <../../../examples/inference/lisa_smbhb/get.sh>`

Now run:

.. literalinclude:: ../../../examples/inference/lisa_smbhb/run.sh
   :language: bash

:download:`Download <../../../examples/inference/lisa_smbhb/run.sh>`

This will run the ``dynesty`` sampler. When it is done, you will have a file called
``lisa_smbhb.hdf`` which contains the results. It should take about three minutes to
run.

To plot the posterior distribution after the last iteration, you can run the following simplified script:

.. literalinclude:: ../../../examples/inference/lisa_smbhb/plot.sh
   :language: bash

:download:`Download <../../../examples/inference/lisa_smbhb/plot.sh>`

Or you can run the advanced one:

.. literalinclude:: ../../../examples/inference/lisa_smbhb/advanced_plot.py
   :language: python

:download:`Download <../../../examples/inference/lisa_smbhb/advanced_plot.py>`

You can modify this advanced plot script to generate the posterior of any SMBHB signals in the LDC Sangria dataset. 
In this example it will create the following plot:

.. image:: ../../_include/lisa_smbhb_mass_tc_0.png
   :scale: 60
   :align: center

The scatter points show each walker's position after the last iteration. The
points are colored by the SNR at that point, with the 50th and 90th
percentile contours drawn. The red lines represent the true parameters of injected signal.
