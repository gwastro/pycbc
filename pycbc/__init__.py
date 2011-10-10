"""
PyCBC
=====

Provides
  1. A collection of classes for processing objects for gravitational wave 
     search and analysis
  2. processing objects can be implemented at different processing 
     architectures: On the Host CPU (C/C++), on a GPU (OpenCl, CUDA)
  3. Bindings to other GW search and analysis frameworks (Lalsuite)

"""

__author__ = 'Karsten Wiesner <karsten.wiesner@ligo.org>'
__all__ = ["datavecstim_opencl, datavecterm_cpu, datavector, fft, 
            highpassfilter, injector, matchedfilter, overwhiteningfilter,
            resampler, singledetectorevent, straindata, templatebank"]
