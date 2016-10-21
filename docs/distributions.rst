###################################################
Using PyCBC Distributions from PyCBC Inference
###################################################

The aim of this page is to demonstrate some simple uses of the distributions available from the distributions.py module available at <https://github.com/ligo-cbc/pycbc/blob/master/pycbc/inference/distributions.py>.

=========================================
Making Mass Distributions in M1 and M2
=========================================

.. literalinclude:: ../examples/inference/massExamples.py
.. command-output:: python ../examples/inference/massExamples.py

========================================================
Sky Location Distribution as Spin Distribution Example 
========================================================
Here we can make a distribution of spins of unit length with equal distribution in x, y, and z to sample the surface of a unit sphere.
.. literalinclude:: ../examples/inference/spinExamples.py
.. command-output:: python ../examples/inference/spinExamples.py

Using much of the same code we can also show how this sampling accurately covers the surface of a unit sphere.
.. literalinclude:: ../examples/inferences/spinSpatialDistrExample.py
.. command-output:: python ../examples/inference/spinSpatialExample.py
