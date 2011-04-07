from datavector_base import *

class DataVectorOpenClGlobal(DataVectorBase):
    
    def __init__(self, element_type=0, length=0):
        print "DataVectorOpenClGlobal.__init__ called"
        
        self.data_vector= [length] #:real_vector_t(length, true) 
        
        super(DataVectorOpenClGlobal, self).__init__(element_type, length, 
              4, 'gpu_opencl_global_memory')
    