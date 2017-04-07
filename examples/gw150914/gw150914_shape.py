from pycbc.frame import read_frame
from pycbc.filter import highpass_fir, lowpass_fir, matched_filter
from pycbc.waveform import get_fd_waveform 
from pycbc.psd import welch, interpolate
import urllib
import pylab

for ifo in ['H1', 'L1']:
    # Read data and remove low frequency content
    fname = '%s-%s_LOSC_4_V2-1126259446-32.gwf' % (ifo[0], ifo)
    url = "https://losc.ligo.org/s/events/GW150914/" + fname
    urllib.urlretrieve(url, filename=fname)
    h1 = read_frame(fname, '%s:LOSC-STRAIN' % ifo)
    h1 = highpass_fir(h1, 15, 8)

    # Calculate the noise spectrum
    psd = interpolate(welch(h1), 1.0 / 32)

    # whiten
    white_strain = (h1.to_frequencyseries() / psd ** 0.5 * psd.delta_f).to_timeseries()

    # remove some of the high and low
    smooth = highpass_fir(white_strain, 35, 8)
    smooth = lowpass_fir(white_strain, 300, 8)

    # time shift and flip L1
    if ifo == 'L1':
        smooth *= -1
        smooth.roll(int(.007 / smooth.delta_t))

    pylab.plot(smooth.sample_times, smooth, label=ifo)

pylab.legend()
pylab.xlim(1126259462.21, 1126259462.45)
pylab.ylabel('Smoothed-Whitened Strain')
pylab.grid()
pylab.ylim(-5, 5)
pylab.xlabel('GPS Time (s)')
pylab.show()


