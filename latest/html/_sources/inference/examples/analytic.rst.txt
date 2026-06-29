.. _inference_example_analytic:

-----------------------------------
Running on an analytic distribution
-----------------------------------

Several analytic distributions are available to run tests on. These can be run
quickly on a laptop to check that a sampler is working properly.

This example demonstrates how to sample a 2D normal distribution with the
``emcee`` sampler. First, we create the following configuration file:

.. literalinclude:: ../../../examples/inference/analytic-normal2d/normal2d.ini
   :language: ini

:download:`Download <../../../examples/inference/analytic-normal2d/normal2d.ini>`

By setting the model name to ``test_normal`` we are using
:py:class:`TestNormal <pycbc.inference.models.analytic.TestNormal>`.
The number of dimensions of the distribution is set by the number of
``variable_params``. The names of the parameters do not matter, just that just
that the prior sections use the same names.

Now run:

.. literalinclude:: ../../../examples/inference/analytic-normal2d/run.sh
   :language: bash

:download:`Download <../../../examples/inference/analytic-normal2d/run.sh>`

This will run the ``emcee`` sampler on the 2D analytic normal distribution with
5000 walkers for 100 iterations. When it is done, you will have a file called
``normal2d.hdf`` which contains the results. It should take about a minute to
run. If you have a computer with more cores, you can increase the
parallelization by changing the ``nprocesses`` argument.

To plot the posterior distribution after the last iteration, run:

.. literalinclude:: ../../../examples/inference/analytic-normal2d/plot.sh
   :language: bash

:download:`Download <../../../examples/inference/analytic-normal2d/plot.sh>`

This will create the following plot:

.. image:: ../../../examples/inference/analytic-normal2d/posterior-normal2d.png
   :scale: 30
   :align: center

The scatter points show each walker's position after the last iteration. The
points are colored by the log likelihood at that point, with the 50th and 90th
percentile contours drawn.

See below for more information about using ``pycbc_inference_plot_posterior``.

To make a movie showing how the walkers evolved, run:

.. literalinclude:: ../../../examples/inference/analytic-normal2d/make_movie.sh
   :language: bash

:download:`Download <../../../examples/inference/analytic-normal2d/make_movie.sh>`

.. note::
   You need ``ffmpeg`` installed for the mp4 to be created.

:ref:`For more information on pycbc_inference_plot_movie <inference_make_movie>`.
