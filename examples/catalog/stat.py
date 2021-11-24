import matplotlib.pyplot as pp
import pycbc.catalog


c = pycbc.catalog.Catalog(source='gwtc-2')
mchirp, elow, ehigh = c.median1d('mchirp', return_errors=True)
spin = c.median1d('chi_eff')

pp.errorbar(mchirp, spin, xerr=[-elow, ehigh], fmt='o', markersize=7)
pp.xlabel('Chirp Mass')
pp.xscale('log')
pp.ylabel('Effective Spin')
pp.show()
