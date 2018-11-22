###################################################
Catalog of Observed Gravitational-wave Mergers
###################################################

Information about the growing catalog of gravitational-wave mergers can be 
accessed through the :py:mod:`pycbc.catalog` package.

===========================================
Which mergers do we have information about?
===========================================

.. literalinclude:: ../examples/catalog/what.py
.. command-output:: python ../examples/catalog/what.py

==============================================
Plotting some key statistics from the catalog
==============================================

.. plot:: ../examples/catalog/stat.py
   :include-source:

==============================================
Accessing data around each event
==============================================

Data around each event can also be easily accessed for any detector.

.. plot:: ../examples/catalog/data.py
   :include-source:
