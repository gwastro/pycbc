"""
"""

from pycbc.array import Array

class TimeSeries(Array):
    def __init__(self,initial_array, time_step, gps_start_time = 0 , dtype=None, copy=True):
        Array.__init__(self,initial_array,dtype=None,copy=True)

        self.time_step = time_step
        self.gps_start_time = gps_start_time

    # This is required so that an operation between two TimeSeries returns a TimeSeries and not an Array
    def _return(self,ary):
        return TimeSeries(ary,self.time_step,gps_start_time=self.gps_start_time, copy=False)

    # This is where you can add typechecking (dt , epoch , etc) to all of the operations that work with another Array type
    def _typecheck(self,other):
        if not isinstance(other,Array):
            return NotImplemented
