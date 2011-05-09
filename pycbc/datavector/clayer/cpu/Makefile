RM = /bin/rm
PYTHON = python
SWIG = swig

# patching the proxy class definition of the swig-gcc output in order to inherit a swigged class
# from a python base class:
#PATCHSCRIPT = inherit_swigged_class_from_py_base.py
#PARENTCLASS = real_vector_t
#BASECLASSMODULE = datavector_base
#BASECLASS = DataVectorBase

all: datavectorcpu.i datavectorcpu_types.h datavectorcpu_prototypes.h datavectorcpu.c
	$(SWIG) -python datavectorcpu.i
	$(PYTHON) setup.py build_ext --inplace
	#$(PYTHON) $(PATCHSCRIPT) $(PARENTCLASS) $(BASECLASSMODULE) $(BASECLASS) \
	#				< datavectorcpu.py > ../../datavectorcpu.py
	cp datavectorcpu.py ../../
	cp _datavectorcpu.so ../../

tidy:
	$(RM) -f *~
clean: tidy
	$(RM) -rf datavectorcpu_wrap.c datavectorcpu.py datavectorcpu.pyc _datavectorcpu.so build
 
