import pycbc.catalog

c = pycbc.catalog.Catalog()

# Names of mergers in the catalog
print(c.names)

# Approximate GPS time of the mergers
print([c[m].time for m in c])
