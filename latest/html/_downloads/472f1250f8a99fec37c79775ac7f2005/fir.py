# Apply an FIR filter. The algorithm is written for high performance so if you
# have a large number of taps, it will resort to a FFT based implementation
# under the hood.
import pycbc.types
import pycbc.filter.resample

# Reference time series
ts = pycbc.types.TimeSeries([-1, 1, -1, 1, -1], delta_t=1.0)

# May also be a numpy array
coeff = pycbc.types.Array([1.0, 0, 1.0])

ts_filtered = pycbc.filter.resample.lfilter(coeff, ts)

# If you want to have a zero phase filter provide a symmetric set of coefficients
# The time delay will be compensated for.

ts_filtered2 = pycbc.filter.resample.fir_zero_filter(coeff, ts)

