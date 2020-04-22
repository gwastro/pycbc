-----------------------------------
Using the relative model
-----------------------------------

The relative model is useful for when you know the parameters of your signal
(say the masses, spins, etc of a merger) that are near the peak of the
likelihood. In this case, we don't need to calculate the likelihood at the full
frequency resolution in order to accurately sample the neighborhood of the
maximum likelihood. This can greatly speed up the calculation of the likelihood.
To use this model you provide the near-peak parameters as fixed arguments as in
the configuration file below.

This example demonstrates using the ``relative`` model with the
``emcee_pt`` sampler. First, we create the following configuration file:

.. literalinclude:: ../../../examples/inference/relative/relative.ini
   :language: ini

:download:`Download <../../../examples/inference/relative/relative.ini>`

For this example, we'll need to download gravitational-wave data for GW170817:

.. literalinclude:: ../../../examples/inference/relative/get.sh
   :language: bash

By setting the model name to ``relative`` we are using
:py:class:`Relative <pycbc.inference.models.relbin.Relative>`.

Now run:

.. literalinclude:: ../../../examples/inference/relative/run.sh
   :language: bash

:download:`Download <../../../examples/inference/relative/run.sh>`

This will run the ``emcee_pt`` sampler. When it is done, you will have a file called
``relative.hdf`` which contains the results. It should take about a minute or two to
run.

To plot the posterior distribution after the last iteration, run:

.. literalinclude:: ../../../examples/inference/relative/plot.sh
   :language: bash

:download:`Download <../../../examples/inference/relative/plot.sh>`

This will create the following plot:

.. image:: ../../_include/relative.png
   :scale: 30
   :align: center

The scatter points show each walker's position after the last iteration. The
points are colored by the signal-to-noise ratio at that point.
