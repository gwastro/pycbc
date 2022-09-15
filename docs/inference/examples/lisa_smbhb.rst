--------------------------------------------------
LISA parameter estimation for SMBHB by using BBHx
--------------------------------------------------

This example shows how to use PyCBC for parameter estimation of supermassive black hole binaries (SMBHB) 
in LISA data. The `data <https://zenodo.org/record/7078835>` are generated from 
`LISA Data Challenge 2a: Sangria <https://lisa-ldc.lal.in2p3.fr/challenge2a>`_, 
and `BBHx <https://github.com/mikekatz04/BBHx>` package is used to generate the ``IMRPhenomD`` template and calculate 
the corresponding TDI response for LISA. :ref:`Relative binning (heterodyned likelihood)<relative>`
is used during sampling to speed up the computation of likelihood functions. Before doing parameter estimation, 
you need to install BBHx and `the corresponding PyCBC waveform plugin <https://github.com/ConWea/BBHX-waveform-model>`
(in the future, this waveform plugin will be a part of BBHx).

First, we create the following configuration file, here we just set chirp mass, mass ratio and tc as variable parameters:

.. literalinclude:: ../../../examples/inference/lisa_smbhb/lisa_smbhb_relbin.ini
   :language: ini

:download:`Download <../../../examples/inference/lisa_smbhb/lisa_smbhb_relbin.ini>`

In this simple example, we do the parameter estimation for the first SMBHB signal in the LDC Sangria dataset,
we need download the data first:

.. literalinclude:: ../../../examples/inference/lisa_smbhb/get.sh
   :language: bash

:download:`Download <../../../examples/inference/lisa_smbhb/get.sh>`

By setting the model name to ``relative`` we are using
:py:class:`Relative <pycbc.inference.models.relbin.Relative>`.

Now run:

.. literalinclude:: ../../../examples/inference/lisa_smbhb/run.sh
   :language: bash

:download:`Download <../../../examples/inference/lisa_smbhb/run.sh>`

This will run the ``dynesty`` sampler. When it is done, you will have a file called
``lisa_smbhb.hdf`` which contains the results. It should take about three minutes to
run.

To plot the posterior distribution after the last iteration, run:

.. literalinclude:: ../../../examples/inference/lisa_smbhb/plot.py
   :language: python

:download:`Download <../../../examples/inference/lisa_smbhb/plot.py>`

This will create the following plot:

.. image:: ../../_include/lisa_smbhb_mass_tc_0.png
   :scale: 30
   :align: center

The scatter points show each walker's position after the last iteration. The
points are colored by the SNR at that point, with the 50th and 90th
percentile contours drawn.
