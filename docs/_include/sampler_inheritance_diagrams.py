# Creates RST for the sampler inheritance diagrams
from __future__ import print_function
from pycbc.inference.sampler import samplers

fname = 'sampler_inheritance_diagrams.rst'

tmplt = """.. _inheritance-{name}:

* ``{name}``:

.. inheritance-diagram:: {module}.{clsname}
   :parts: 3

|

"""
fp = open(fname, 'w')
for sampler, cls in sorted(samplers.items()):
    out = tmplt.format(name=sampler, clsname=cls.__name__,
                       module=cls.__module__)
    print(out, file=fp)

fp.close()
