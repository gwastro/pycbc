.. _single_template_examples:

######################################
Using the single template model
######################################

-------------------------------
Quickstart example
-------------------------------

The single template model is useful for when you know all the intrinsic parameters
of your signal (say the masses, spins, etc of a merger). In this case, we don't
need to recalculate the waveform model to sample different possible extrinsic
parameters (i.e. distance, sky location, inclination). This can greatly
speed up the calculation of the likelihood. To use this model
you provide the intrinsic parameters as fixed arguments as in the configuration
file below.

This example demonstrates using the ``single_template`` model with the
``dynesty`` sampler. First, we create the following configuration file:

.. dropdown:: Configuration File
    :animate: fade-in-slide-down

    .. literalinclude:: ../../../examples/inference/single/single_simple.ini
       :language: ini

    :download:`Download <../../../examples/inference/single/single_simple.ini>`

For this example, we'll need to download gravitational-wave data for GW170817:

.. literalinclude:: ../../../examples/inference/single/get.sh
   :language: bash

By setting the model name to ``single_template`` we are using
:py:class:`SingleTemplate <pycbc.inference.models.single_template.SingleTemplate>`.

Now run:

.. literalinclude:: ../../../examples/inference/single/run.sh
   :language: bash

:download:`Download <../../../examples/inference/single/run.sh>`

This will run the ``dynesty`` sampler. When it is done, you will have a file called
``single.hdf`` which contains the results. It should take about a minute or two to
run.

To plot the posterior distribution, run:

.. literalinclude:: ../../../examples/inference/single/plot.sh
    :language: bash

:download:`Download <../../../examples/inference/single/plot.sh>`

This will create the following plot:

.. image:: ../../_include/single.png
    :width: 400
    :align: center

The scatter points show position of different posterior samples. The
points are colored by the log likelihood at that point, with the 50th and 90th
percentile contours drawn.

------------------------------------------------------------
Marginalization subset of parameters
------------------------------------------------------------
The single template model supports marginalization over its parameters. This
can greatly speed up parameter estimation. The
marginalized parameters can also be recovered as in the example below which
extends the previous example to include all the parameters the model supports.

In this example, the sampler will only explore inclination, whereas all other
parameters are either numerically or analytically marginalized over. Note
that the marginalization still takes into account the priors you choose. This
includes if you have placed boundaries or constraints.

First, you'll need the configuration file.

.. dropdown:: Configuration File
    :animate: fade-in-slide-down

    .. literalinclude:: ../../../examples/inference/single/single.ini
       :language: ini

    :download:`Download <../../../examples/inference/single/single.ini>`

Run this script to use this configuration file:

.. literalinclude:: ../../../examples/inference/single/run_marg.sh
   :language: bash

:download:`Download <../../../examples/inference/single/run_marg.sh>`

This will create the following plots:

Before demarginalization:

.. image:: ../../_include/single_marg.png
    :width: 400
    :align: center

After demarginalization:

.. image:: ../../_include/single_demarg.png
   :width: 400
   :align: center

-----------------------------------------------------------
Marginalization over all parameters
-----------------------------------------------------------
All parameters that the single_template model supports can be marginalized over
to do so, you would just configure the model section to include each parameter
(minus any you choose to be static). PyCBC Inference will simply reconstruct
the missing parameters rather than performing any sampling over the parameter
space. The number of samples to take from the posterior is configurable.

First, you'll need the configuration file. Note that in the `sampler` section
you can set the number of samples to draw from the posterior.

.. dropdown:: Configuration File
    :animate: fade-in-slide-down

    .. literalinclude:: ../../../examples/inference/single/single_instant.ini
       :language: ini

    :download:`Download <../../../examples/inference/single/single.ini>`

Run this script to use this configuration file:

.. literalinclude:: ../../../examples/inference/single/run_instant.sh
   :language: bash

:download:`Download <../../../examples/inference/single/run_marg.sh>`

This will create the following plot:

After demarginalization:

.. image:: ../../_include/single_instant.png
   :width: 400
   :align: center

---------------------------------------------------
Abitrary sampling coordinates with nested samplers
---------------------------------------------------

The single template model also supports marginalization over the polarization
angle by numerical sampling. The following example features arbitrary sampling
coordinates with nested samplers
using the ``fixed_samples`` distribution. Here we sample in the time delay
space rather than sky location directly. This functionality is generic
to pycbc inference and could be used with other models or samplers.

.. dropdown:: Configuration File
    :animate: fade-in-slide-down

    .. literalinclude:: ../../../examples/inference/single/single_adv.ini
        :language: ini
