-----------------------------------
Using the single template model
-----------------------------------

The single template model is useful for when you know all the intrinsic parameters
of your signal (say the masses, spins, etc of a merger). In this case, we don't
need to recalculate the waveform model to sample different possible extrinsic 
parameters (i.e. distance, sky location, inclination). This can greatly
speed up the calculation of the likelihood. To use this model
you provide the intrinsic parameters as fixed arguments as in the configuration
file below. 

This example demonstrates using the ``single_template`` model with the
``emcee_pt`` sampler. First, we create the following configuration file:

.. literalinclude:: ../../../examples/inference/single/single.ini
   :language: ini

:download:`Download <../../../examples/inference/single/single.ini>`

For this example, we'll need to download gravitational-wave data for GW170817:

.. literalinclude:: ../../../examples/inference/single/get.sh
   :language: bash

By setting the model name to ``single_template`` we are using
:py:class:`SingleTemplate <pycbc.inference.models.single_template.SingleTemplate>`.

Now run:

.. literalinclude:: ../../../examples/inference/single/run.sh
   :language: bash
   
:download:`Download <../../../examples/inference/single/run.sh>`

This will run the ``emcee_pt`` sampler. When it is done, you will have a file called
``single.hdf`` which contains the results. It should take about a minute or two to
run.

To plot the posterior distribution after the last iteration, run:

.. literalinclude:: ../../../examples/inference/single/plot.sh
   :language: bash

:download:`Download <../../../examples/inference/single/plot.sh>`

This will create the following plot:

.. image:: ../../_include/single.png
   :scale: 30
   :align: center

The scatter points show each walker's position after the last iteration. The
points are colored by the log likelihood at that point, with the 50th and 90th
percentile contours drawn.
