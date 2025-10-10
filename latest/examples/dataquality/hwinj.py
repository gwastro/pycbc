"""This example shows how to determine when a CBC hardware injection is present
in the data from a detector.
"""

import matplotlib.pyplot as pp
from pycbc import dq


start_time = 1126051217
end_time = start_time + 10000000

# Get times that the Livingston detector has CBC injections into the data
segs = dq.query_flag('L1', 'CBC_HW_INJ', start_time, end_time)

pp.figure(figsize=[10, 2])
for seg in segs:
    start, end = seg
    pp.axvspan(start, end, color='blue')

pp.xlabel('Time (s)')
pp.show()

