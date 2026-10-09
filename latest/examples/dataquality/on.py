"""This example shows how to determine when a detector is active."""

import matplotlib.pyplot as pp
from pycbc import dq
from pycbc.results import ifo_color


start_time = 1126051217
end_time = start_time + 100000

# Get times that the Hanford detector has data
hsegs = dq.query_flag('H1', 'DATA', start_time, end_time)

# Get times that the Livingston detector has data
lsegs = dq.query_flag('L1', 'DATA', start_time, end_time)

pp.figure(figsize=[10,2])
for seg in lsegs:
    start, end = seg
    pp.axvspan(start, end, color=ifo_color('L1'), ymin=0.1, ymax=0.4)

for seg in hsegs:
    start, end = seg
    pp.axvspan(start, end, color=ifo_color('H1'), ymin=0.6, ymax=0.9)

pp.xlabel('Time (s)')
pp.show()
