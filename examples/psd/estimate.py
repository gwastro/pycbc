import pycbc.noise
import pycbc.psd
import pylab

# generate some colored gaussian noise
flow = 30.0
delta_f = 1.0 / 16
flen = int(2048 / delta_f) + 1
psd = pycbc.psd.aLIGOZeroDetHighPower(flen, delta_f, flow)

### Generate 128 seconds of noise at 4096 Hz
delta_t = 1.0 / 4096
tsamples = int(128 / delta_t)
ts = pycbc.noise.noise_from_psd(tsamples, delta_t, psd, seed=127)

# Estimate the PSD
## We'll choose 4 seconds PSD samples that are overlapped 50 %
seg_len = int(4 / delta_t)
seg_stride = seg_len / 2
estimated_psd = pycbc.psd.welch(ts,
                      seg_len=seg_len,
                      seg_stride=seg_stride)

pylab.loglog(estimated_psd.sample_frequencies, estimated_psd, label='estimate')
pylab.loglog(psd.sample_frequencies, psd, linewidth=3, label='known psd')
pylab.xlim(xmin=flow, xmax=2000)
pylab.ylim(1e-47, 1e-45)
pylab.legend()
pylab.grid()
pylab.show()
