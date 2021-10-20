# Let's see when the detector was active.
# Note that all units are in seconds and that the +
from pycbc import dq
import pylab

start_time = 1126051217
end_time = start_time + 100000

# Get times that the Hanford detector has data
hsegs = dq.query_flag('H1', 'DATA', start_time, end_time)

# Get times that the Livingston detector has data
lsegs = dq.query_flag('L1', 'DATA', start_time, end_time)

pylab.figure(figsize=[10,2])
for seg in lsegs:
    start, end = seg
    pylab.axvspan(start, end, color='green', ymin=0.1, ymax=0.4)

for seg in hsegs:
    start, end = seg
    pylab.axvspan(start, end, color='red', ymin=0.6, ymax=0.9)

pylab.xlabel('Time (s)')
pylab.show()
