from pycbc.frame import read_frame
from pycbc.filter import highpass_fir, lowpass_fir
from pycbc.psd import welch, interpolate
from pycbc.types import TimeSeries
try:
    from urllib.request import urlretrieve
except ImportError:  # python < 3
    from urllib import urlretrieve

# Read data and remove low frequency content
fname = 'H-H1_LOSC_4_V2-1126259446-32.gwf'
url = "https://www.gw-openscience.org/GW150914data/" + fname
urlretrieve(url, filename=fname)
h1 = highpass_fir(read_frame(fname, 'H1:LOSC-STRAIN'), 15.0, 8)

# Calculate the noise spectrum and whiten
psd = interpolate(welch(h1), 1.0 / 32)
white_strain = (h1.to_frequencyseries() / psd ** 0.5 * psd.delta_f).to_timeseries()

# remove some of the high and low frequencies
smooth = highpass_fir(white_strain, 25, 8)
smooth = lowpass_fir(white_strain, 250, 8)

#strech out and shift the frequency upwards to aid human hearing
fdata = smooth.to_frequencyseries()
fdata.roll(int(1200 / fdata.delta_f))
smooth = TimeSeries(fdata.to_timeseries(), delta_t=1.0/1024)

#Take slice around signal
smooth = smooth[len(smooth)//2 - 1500:len(smooth)//2 + 3000]
smooth.save_to_wav('gw150914_h1_chirp.wav')

