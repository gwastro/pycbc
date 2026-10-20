import pycbc.catalog

c = pycbc.catalog.Catalog(source='gwtc-2')

# Names of mergers in the catalog
print(c.names)

# Approximate GPS time of the mergers
print([c[m].time for m in c])
