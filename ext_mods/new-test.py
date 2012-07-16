import testlal
from numpy import *

myplan=testlal.lal.CreateForwardCOMPLEX8FFTPlan(4,0)

a=array([0,0,0,0],dtype=complex64)
b=array([0,0,0,0],dtype=complex64)
b.fill(1+1j)

testlal.XLALCOMPLEX8VectorFFT(a,b,myplan)
