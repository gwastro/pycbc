.. _marginalized_time_example:

######################################
Marginalized time model
######################################

-------------------------------
Example with GW150914
-------------------------------

The marginalized time model can handle a broad array of waveform models which
it generates at full resolution like the ``gaussian_noise`` model. However,
it is optimized to enable marginalization of time in addition to marginalization
over sky location, polarization, and overall phase (valid if the waveform
approximant is dominant mode only).

This example demonstrates using the ``marginalized_time`` model with the
``dynesty`` sampler in a configuration designed to run in couple minutes on a
laptop. Actual sampling will only occur over the component masses and
inclination. The remaining parameters are marginalized over, but will
be reconstructed after running the parameter estimation by a follow-up
script.

First, we create the following configuration file:

.. dropdown:: Configuration File
    :animate: fade-in-slide-down

    .. literalinclude:: ../../../examples/inference/margtime/margtime.ini
       :language: ini

    :download:`Download <../../../examples/inference/margtime/margtime.ini>`

For this example, we'll need to download gravitational-wave data for GW150914:

.. literalinclude:: ../../../examples/inference/margtime/get.sh
   :language: bash

By setting the model name to ``marginalized_time`` we are using
:py:class:`MarginalizedTime <pycbc.inference.models.marginalized_gaussian_noise.MarginalizedTime>`.

Now run the following script. Note that after the parameter inference is run, we reconstruct
the marginalized parameters by using the `pycbc_inference_model_stats` script with the options
as follows.

.. literalinclude:: ../../../examples/inference/margtime/run.sh
   :language: bash

:download:`Download <../../../examples/inference/margtime/run.sh>`

This will run the ``dynesty`` sampler. When it is done, you will have a file called
``demarg_150914.hdf`` which contains the results. It should take just a few minutes to run.

This will create the following plot:

.. image:: ../../_include/demarg_150914.png
    :width: 400
    :align: center

The scatter points show position of different posterior samples.
