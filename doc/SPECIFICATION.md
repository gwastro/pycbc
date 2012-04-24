# The pyCBC specification                                                      

[TOC]
\tableofcontents

This document is the pyCBC specification, the root of everything.
In principle this is a summarized version of the accepted and
implemented PEPs (pyCBC Enhancement Proposals).

## About pyCBC ###
### The goal of pyCBC ###

- Create a flexible, extensible sw production code for CBC analysis that can be 
  released for the public
- Unix small tools philosophy contains small units that make it simple to do 
  simple things and that makes it easy to assemble complex analysis codes
- Enables simple, easy and transparent access for various many-core architectures
  like GPUs
- Ultimately become the CBC data analysis tool of the 'advanced era'

### Requirements for pyCBC ###

- It is a CBC only package
- High level layer will be written in python only
- Low level layer could be C/C++/OpenCL/CUDA
- Modular, looselycoupled modules
- Swig will be used for the binding
- Supported platforms are that of the LSC reference platforms
- High abstraction allows rapid prototyping (like writing a formula) at the top 
  layer (user frontend)
- Design of basic operations that can be used within the bigger boxes (reusability)
- Allows speed optimization by the use of GPUs and upcoming many-core hardware (APUs)
- System decides on optimal implementation in terms of a heterogeneous compute 
  architecture (CPU/GPU) automatically at runtime
- Parallelism should be hidden for the user (transparent heterogeneous computing)
- Shall run on CPU-only systems (even without CUDA or OpenCl installed)
- Easy, straightforward installation
- Test driven design
- Module verification against golden reference implementations
- Code quality is good enough for reviews by the public or government institutes
- Bottom up design
- Distributed PEP-style teamwork


## General considerations ##
This section contains the definitions, agreements, guidelines which are more 
connected to the pyCBC framework itself rather than to an actual piece of
code.

### Definitions, terms, expressions ###

Thorough this document the words SHALL, MUST, MAY are to be understood as 
it is defined in  [IETF RFC 2119](http://www.ietf.org/rfc/rfc2119.txt). 

### Tools for development and communications
- **Repository**: The pyCBC source base is stored in the 
 [GIT repository of UWM](https://sugwg-git.phy.syr.edu/cgit/pycbc/).
- **Bugs, PEPs**:  Bugs are to be submitted, tracked and handled in the 
  [RedMine bug tracker](https://bugs.ligo.org/redmine/projects/pycbc-production/).
- **Mailing list**: Send an email to: pycbc+cbc@gravity.phys.uwm.edu
- **Project home page**: http://www.gravity.phy.syr.edu/pycbc/

### The development cycles

The pyCBC project is a common LIGO - Virgo community effort. Anyone who accepts the 
development rules listed below is welcome to join. Any feature or change 
request should be submitted in the form of a PEP (pyCBC Enhancement Proposal)
strictly following the format defined by the PEP template
[PEP #312](https://bugs.ligo.org/redmine/issues/312).

### The pyCBC build system ###
The build system for pyCBC is defined in [PEP #313](https://bugs.ligo.org/redmine/issues/313). 
Accordingly, the build and test 
infrastructure for pyCBC should be the Metronome software used by the NMI build and test lab. 
pyCBC will use Metronome to perform a nightly test and build of the software on the 
LIGO Data Grid and EMI/EGI reference platforms. The nightly test will also include running 
unit tests for the library and higher-level testing of pyCBC programs.
Each night, a NMI/Metronome server should perform the following tests of pyCBC:

- Build and test the software from the source git repository.
- Create a distribution tar ball of the git repository.
- Build and test the software from the source tar ball.
- Ensure that it is possible to build the documentation both from the git repository and 
the distribution tar ball.
- Python 2 to 3 compliance testing as appropriate.

"Testing" should include running all unit tests, and any system-level tests that are needed 
for higher-level  programs.

The build and test should be performed on the LIGO Data Grid and EMI/EGI reference operating
systems using an up-to-date LIGO Data Grid and EGI/EMI environment. Errors should be 
reported to the pyCBC project manager for assignment to the responsible person.

### Batch systems and pyCBC ###

As defined in [PEP #325](https://bugs.ligo.org/redmine/issues/325)the pyCBC software pieces 
must be absolutely batch system independent. Each and every line of the code must not 
explicitely or implicetly assume or refeer to any specific type of batch-sytem backend.

- pyCBC executables should be idempotent. In particular, they should not modify in-place 
  their input data to create their output data.
- pyCBC executables should read their input data and write their output data to the 
  directory they are executed from and should not assume any file-system layout. It is the 
  responsibility of the workflow management system to ensure that the input data is present 
  in the directory when pyCBC executables are run and that output data is subsequently 
  staged to the appropriate place
- pyCBC executables should take their parameters and arguments from the command line as 
  getopt_long() compatible arguments, i.e. in the form --option argument. They should 
  not read parameters from configuration files.    
- pyCBC executables should store sufficient metadata in their output files to allow a user 
  to reproduce the analysis that generated the output file.

### The pyCBC documentation tool ###

As defined by [PEP #318](https://bugs.ligo.org/redmine/issues/318) the Doxygen documentation 
system has been chosen to be the source code documentation tool. 
The developers has to  follow only one of the posibble way of source code documentation 
supported by Doxygen.  The commonly agreed form will be  subject of another PEP.

### Style convention for pyCBC code ###

The pyCBC code must follow the style conventions defined by
[Python PEP 8](http://www.python.org/dev/peps/pep-0008/).

### Python versioning for pyCBC ###

Accordinig to [PEP #314](https://bugs.ligo.org/redmine/issues/314) all PyCBC code should 
be written in Python 2.6, as this is the current (or should be) /usr/bin/python on 
LDG machines. Moreover, all such python should be written in the "good" style as specified here:
http://wiki.python.org/moin/Python2orPython3

Some things of particular note are:

- Use new-style classes (inherit ultimately from object) rather than old.
- Use xrange instead of range
- Do not use the print statement (it is acceptable to use the print function).

Other points of particular concern to PyCBC may be added here later as they become apparent.

At this point, it is not part of the specification that 2to3 be run on python code, 
and produce error-free Python 3 code. But that may change in the future.

## The pyCBC specification ##

### The Array Class ###

pycbc.array.Array(object, dtype=None, copy=True)

where,
- object: An array-like object as specified by NumPy, this also includes instances of an 
  underlying data type as described  in section 3 or an instance of the PYCBC Array class 
  itself. This object is used to populate the data of the  array.

- dtype: A NumPy style dtype that describes the type of encapsulated data (float32,compex64, etc)
- copy:

This defines whether the object is copied to instantiate the array or is simply referenced. If copy is false, 
new data is not created, and so the context is ignored. The default is to copy the given object.

- `__len__(self)`
- `__str__(self)`
- `__mul__(self,other)`  Returns a new array that is the element-wise multiplication of other and self.
- `__rmul__(self,other)` Returns a new array that is the element-wise multiplication of other and self.
- `__imul__(self,other)` Does an mulitplication of other and self in place. No new memory is created.
- `__add__(self,other)`  Returns a new array that is the element-wise addition of other and self.
- `__radd__(self,other)` Returns a new array that is the element-wise addition of other and self.
- `__iadd__(self,other)` Does an addition of other and self in place. No new memory is created.
- `__div__(self,other)`  Returns a new array that is the element-wise division of self by other.
- `__rdiv__(self,other)` Returns a new array that is the element-wise division of other by self.
- `__idiv__(self,other)` Does the division of self by other in place. No new memory is created.
- `__sub__(self,other)`  Returns a new array that is the element-wise subtraction of other from self.
- `__rsub__(self,other)` Returns a new array that is the element-wise subtraction of self from other.
- `__isub__(self,other)` Does a subtraction of other from self in place. No new memory is created.
- `__pow__(self,other)`  Returns a new array containing the exponentiation of each element by other.
- `__abs__(self)`        Returns a new array containing the the absolute value of each element.
- `real(self)`           Returns a new array containing the real part of the array, or the original if it is 
only real.
- `imag(self)` Returns a new array containing the imaginary part of the array or None if the array is real.
- `conj(self)` Returns a new array that contains the complex conjugate of the original.
- `sum(self)`  Returns the sum of the array as a cpu scaler.
- `dot(self,other)` Returns the dot product of self and other as a cpu scaler.
- `ptr(self)`  Returns a pointer to the data memory block that is suitable for passing to a c-function.


Allowable types for other Other can either be an instance of Array or a numeric type (int, float ,complex, etc).

Additional restrictions that may be lifted at some later point:

-The underlying Array data type (Numpy, PyCUDA, PyOpenCL) for other must match self

The following private members also exist.

- `_data` A pointer to the underlying data type for this instance
- `_scheme` A pointer to the most recent processing state of the object.
