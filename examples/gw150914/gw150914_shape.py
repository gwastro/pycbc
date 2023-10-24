import matplotlib.pyplot as pp
from pycbc.filter import highpass_fir, lowpass_fir
from pycbc.psd import welch, interpolate
from pycbc.catalog import Merger


for ifo in ['H1', 'L1']:
    # Read data and remove low frequency content
    h1 = Merger("GW150914").strain(ifo)
    h1 = highpass_fir(h1, 15, 8)

    # Calculate the noise spectrum
    psd = interpolate(welch(h1), 1.0 / h1.duration)

    # whiten
    white_strain = (h1.to_frequencyseries() / psd ** 0.5).to_timeseries()

    # remove some of the high and low
    smooth = highpass_fir(white_strain, 35, 8)
    smooth = lowpass_fir(smooth, 300, 8)

    # time shift and flip L1
    if ifo == 'L1':
        smooth *= -1
        smooth.roll(int(.007 / smooth.delta_t))

    pp.plot(smooth.sample_times, smooth, label=ifo)

pp.legend()
pp.xlim(1126259462.21, 1126259462.45)
pp.ylim(-150, 150)
pp.ylabel('Smoothed-Whitened Strain')
pp.grid()
pp.xlabel('GPS Time (s)')
pp.show()
