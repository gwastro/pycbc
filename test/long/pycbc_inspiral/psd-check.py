from pycbc.frame import read_frame
from pycbc.filter import highpass, resample_to_delta_t
import pycbc
from pycbc.types import float32, FrequencySeries, Array, TimeSeries
import matplotlib
matplotlib.use('gtk')
import pylab
import numpy
import lal
import pycbc.psd

def lal_psd_estimate(numPoints, chan):
    """ This function simulates the psd estimation as done in inspiral.c 
    """
    chan = chan.lal()
    df = 4096/ numPoints
    avgSpecParams = lal.AverageSpectrumParams()
    spec = lal.CreateREAL4FrequencySeries("", chan.epoch, 0, df, lal.lalSecondUnit, numPoints/2+1)
    avgSpecParams.method = lal.useMedian;
    avgSpecParams.plan =  lal.CreateREAL4FFTPlan(numPoints, 1, 1) 
    avgSpecParams.overlap = numPoints / 2;
    avgSpecParams.window = lal.CreateHannREAL4Window(numPoints);
    lal.REAL4AverageSpectrum(spec, chan, avgSpecParams)
    return FrequencySeries(spec.data.data, delta_f=spec.deltaF)

pad_data = 8
sample_rate=4096
duration=2048
start_time=1026019572
seg_len = 256*4096
stride = 128*4096
nsegs = 15
delta_f = 1.0 / 256

# Check that we are using the same hann window
w_pycbc = Array(numpy.hanning(seg_len).astype(float32))
w_lal = lal.CreateHannREAL4Window(seg_len)
print "hann props pycbc", w_pycbc.sum(), w_pycbc.squared_norm().sum()
print "hann props lal", w_lal.sum, w_lal.sumofsquares

# Check that the median bias is the same
print "%s segments" % nsegs
print "BIAS pycbc", pycbc.psd.median_bias(nsegs)
print "BIAS lal", lal.RngMedBias(nsegs)

# Check the psd norm
print "PYCBC psd norm", 2.0 / float(sample_rate) / (w_pycbc.squared_norm().sum()) / pycbc.psd.median_bias(nsegs)

# Same strain preparation for lal and pycbc psd estimation
strain = read_frame("LER2.lcf", "H1:FAKE-STRAIN", start_time=start_time-pad_data, duration=duration+pad_data*2)
strain = highpass(strain, frequency=30)
strain *= pycbc.DYN_RANGE_FAC
strain = TimeSeries(strain, delta_t=strain.delta_t, epoch=strain.start_time, dtype=float32)
strain = resample_to_delta_t(strain, 1.0/sample_rate)
strain = strain[int(pad_data*1.0/strain.delta_t):len(strain)-int(1.0/strain.delta_t*pad_data)]

psd_pycbc = pycbc.psd.welch(strain, seg_len=seg_len, seg_stride=stride)
psd_lal = lal_psd_estimate(seg_len, strain)

print psd_pycbc[-1]
print psd_lal[-1]

pylab.figure(1)
pylab.loglog(psd_pycbc.sample_frequencies.numpy(), psd_pycbc.numpy(), label="PyCBC")
pylab.loglog(psd_lal.sample_frequencies.numpy(), psd_lal.numpy(), label="LAL")
pylab.xlabel("Frequency (Hz)")
pylab.xlim(40, 2000)
pylab.ylabel(" (Strain/rhz * DYN_RANGE_FAC)^2 " ) 
pylab.legend()
pylab.savefig("psdloglog.png")

pylab.figure(2)
reldif = (psd_pycbc-psd_lal)/psd_lal
pylab.plot(reldif.sample_frequencies.numpy(), reldif.numpy())
pylab.xlabel("Frequency (Hz)")
pylab.ylabel("Relative Difference") 
pylab.savefig("psdreldif.png")

