###################################################
Using PyCBC Distributions from PyCBC Inference
###################################################

The aim of this page is to demonstrate some simple uses of the distributions available from the distributions.py module available at :py:mod:`pycbc.distributions`.

=========================================
Making Mass Distributions in M1 and M2
=========================================

Here we will demonstrate how to make different mass populations of binaries. This example shows uniform and gaussian mass distributions.

.. plot:: ../examples/distributions/mass_examples.py
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
