import pycbc.catalog, pylab

m = pycbc.catalog.Merger("GW170817")

# Retreive data from the Livingston observatory around the BNS merger
ts = m.strain("L1").time_slice(m.time - 30, m.time + 6)

# Whiten the data with a 4s filter
white = ts.whiten(4, 4)

times, freqs, power = white.qtransform(.01, logfsteps=200,
                                    qrange=(110, 110),
                                    frange=(20, 512))
pylab.pcolormesh(times, freqs, power**0.5)
pylab.yscale('log')
pylab.ylabel("Frequency (Hz)")
pylab.xlabel("Time (s)")
pylab.show()
