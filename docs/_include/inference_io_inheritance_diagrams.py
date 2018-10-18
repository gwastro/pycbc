# Creates RST for the sampler inheritance diagrams
from __future__ import print_function
import inspect
from pycbc.inference.io import filetypes

fname = 'inference_io_inheritance_diagrams.rst'

def get_topclasses(cls):
    """Gets the base classes that are in pycbc."""
    bases = [c for c in inspect.getmro(cls)
             if c.__module__.startswith('pycbc') and c != cls]
    return ', '.join(['{}.{}'.format(c.__module__, c.__name__) for c in bases])

tmplt = """.. _inheritance-io-{name}:

* ``{name}``:

.. inheritance-diagram:: {module}.{clsname}
   :parts: 3
   :top-classes: pycbc.inference.io.base_hdf.BaseInferenceFile

|

"""
fp = open(fname, 'w')

for ftname, cls in sorted(filetypes.items()):
    # get the parents
    topclasses = get_topclasses(cls)
    out = tmplt.format(name=ftname, clsname=cls.__name__,
                       module=cls.__module__)
    print(out, file=fp)

fp.close()
