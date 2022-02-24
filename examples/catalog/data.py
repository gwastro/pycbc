import matplotlib.pyplot as pp
import pycbc.catalog


m = pycbc.catalog.Merger("GW170817", source='gwtc-1')

fig, axs = pp.subplots(2, 1, sharex=True, sharey=True)
for ifo, ax in zip(["L1", "H1"], axs):
    pp.sca(ax)
    pp.title(ifo)
    # Retreive data around the BNS merger
    ts = m.strain(ifo).time_slice(m.time - 15, m.time + 6)

    # Whiten the data with a 4s filter
    white = ts.whiten(4, 4)

    times, freqs, power = white.qtransform(.01, logfsteps=200,
                                        qrange=(110, 110),
                                        frange=(20, 512))
    pp.pcolormesh(times, freqs, power**0.5, vmax=5)

pp.yscale('log')
pp.ylabel("Frequency (Hz)")
pp.xlabel("Time (s)")
pp.show()
