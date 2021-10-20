# Let's see when the detector was active.
# Note that all units are in seconds and that the +
from pycbc import dq
import pylab

start_time = 1126051217
end_time = start_time + 10000000

# Get times that the Livingston detecot has CBC injections into the data
segs = dq.query_flag('L1', 'CBC_HW_INJ', start_time, end_time)

pylab.figure(figsize=[10, 2])
for seg in segs:
    start, end = seg
    pylab.axvspan(start, end, color='blue')

pylab.xlabel('Time (s)')
pylab.show()

