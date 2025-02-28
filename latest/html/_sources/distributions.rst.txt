###################################################
Using PyCBC Distributions from PyCBC Inference
###################################################

The aim of this page is to demonstrate some simple uses of the distributions available from the distributions.py module available at :py:mod:`pycbc.distributions`.

============================================================
Generating samples in a Python script by using the .ini file
============================================================

This example shows how to generate samples in a Python script by using a standalone .ini :download:`file <../examples/distributions/pycbc_bbh_prior.ini>`, 
in this example, we can also learn that when we draw samples from source-frame masses and comoving volume uniformly, the distributions of detector-frame masses are not uniform.

.. plot:: ../examples/distributions/sampling_from_config_example.py
   :include-source:

=========================================
Making Mass Distributions in M1 and M2
=========================================

Here we will demonstrate how to make different mass populations of binaries. This example shows uniform and gaussian mass distributions.

.. plot:: ../examples/distributions/mass_examples.py
   :include-source:

====================================================
Generating mchirp and q from uniform mass1 and mass2
====================================================

This example shows chirp mass and mass ratio samples drawn from uniform mass1 and mass2 with boundaries given by chirp mass and mass ratio. The first row of the figures below are the chirp mass and mass ratio samples, and the same set of samples converted to mass1 and mass2, respectively. The second row are the marginalized distribution for chirp mass and mass ratio, with comparisons for analytical probability density function (PDF).

.. plot:: ../examples/distributions/mchirp_q_from_uniform_m1m2_example.py
   :include-source:

========================================================
Sky Location Distribution as Spin Distribution Example 
========================================================

Here we can make a distribution of spins of unit length with equal distribution in x, y, and z to sample the surface of a unit sphere.

.. plot:: ../examples/distributions/spin_examples.py
   :include-source:

Using much of the same code we can also show how this sampling accurately covers the surface of a unit sphere.

.. plot:: ../examples/distributions/spin_spatial_distr_example.py
   :include-source:
