import pycbc.catalog, pylab

c = pycbc.catalog.Catalog()
mchirp, elow, ehigh = c.median1d('mchirp', return_errors=True)
spin = c.median1d('chi_eff')

pylab.errorbar(mchirp, spin, xerr=[-elow, ehigh], fmt='o', markersize=7)
pylab.xlabel('Chirp Mass')
pylab.xscale('log')
pylab.ylabel('Effective Spin')
pylab.show()
