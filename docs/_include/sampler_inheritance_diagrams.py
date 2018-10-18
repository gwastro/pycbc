# Creates RST for the sampler inheritance diagrams
from __future__ import print_function
from pycbc.inference.sampler import samplers

print("CREATING SAMPLER DIAGRAM")

fname = 'sampler_inheritance_diagrams.rst'

tmplt = """.. _inheritance-{}:

* ``{}``:

.. inheritance-diagram:: pycbc.inference.sampler.{}
   :parts: 3

|

"""
fp = open(fname, 'w')

for sampler, cls in sorted(samplers.items()):
    print(tmplt.format(sampler, sampler, cls.__name__), file=fp)

fp.close()
