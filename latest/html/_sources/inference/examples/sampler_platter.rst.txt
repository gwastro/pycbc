.. _inference_example_samplers:

------------------------------------------------------
Trying out different samplers
------------------------------------------------------

This page shows basic configurations for some of the samplers that PyCBC 
Inference supports. Note that these samplers often have many additional
options which you can learn apart from their respective library documentation
pages. There very likely are more optimal settings than the example ones here
for your problem as these are configured to run quickly.

We'll use a very simple analytic model to test each sampler. The following is
the configuration to set up a unimodal gaussian likelihood. We'll have each
sampler try to fit this.

.. literalinclude:: ../../../examples/inference/samplers/simp.ini
   :language: ini
   
Each sampler needs nees a slightly different configuration. Below are basic
configurations which can analyze this simple problem.

===================================================
`Emcee <https://emcee.readthedocs.io/en/v2.2.1/>`_
===================================================

.. literalinclude:: ../../../examples/inference/samplers/emcee_stub.ini
   :language: ini
   
======================================================================
`Parallel-tempered Emcee <https://emcee.readthedocs.io/en/v2.2.1/>`_
======================================================================
   
.. literalinclude:: ../../../examples/inference/samplers/emcee_pt_stub.ini
   :language: ini
   
===================================================
`PTEmcee <https://github.com/willvousden/ptemcee>`_
===================================================
   
.. literalinclude:: ../../../examples/inference/samplers/ptemcee_stub.ini
   :language: ini
   
===============================================
`Dynesty <https://dynesty.readthedocs.io/>`_
===============================================
   
.. literalinclude:: ../../../examples/inference/samplers/dynesty_stub.ini
   :language: ini
   
============================================================
`Ultranest <https://johannesbuchner.github.io/UltraNest/>`_
============================================================
   
.. literalinclude:: ../../../examples/inference/samplers/ultranest_stub.ini
   :language: ini
   
===============================================
`Epsie <https://github.com/cdcapano/epsie>`_
===============================================
   
.. literalinclude:: ../../../examples/inference/samplers/epsie_stub.ini
   :language: ini
   
The following are also supported, but require either python3 support only (cpnest)
or an external package (multinest).
   
================================================
`cpnest <https://github.com/johnveitch/cpnest>`_
================================================
   
.. literalinclude:: ../../../examples/inference/samplers/cpnest_stub.ini
   :language: ini
   
=============================================================
`Multinest <https://johannesbuchner.github.io/PyMultiNest/>`_
=============================================================
   
.. literalinclude:: ../../../examples/inference/samplers/multinest_stub.ini
   :language: ini
   
============================================================
`Snowline <https://johannesbuchner.github.io/snowline/>`_
============================================================
   
.. literalinclude:: ../../../examples/inference/samplers/snowline_stub.ini
   :language: ini
   
   
============================================================
`nessai <https://github.com/mj-will/nessai>`_
============================================================
   
.. literalinclude:: ../../../examples/inference/samplers/nessai_stub.ini
   :language: ini
   
If we run these samplers, we create the following plot:

.. image:: ../../_include/sample.png
   :scale: 30
   :align: center
