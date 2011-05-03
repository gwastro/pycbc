RM = /bin/rm
PYTHON = python
SWIG = swig

all: datavectorcpu.i datavectorcpu_types.h datavectorcpu_prototypes.h datavectorcpu.c
	$(SWIG) -python datavectorcpu.i
	$(PYTHON) setup.py build_ext --inplace
	cp datavectorcpu.py ../../
	cp _datavectorcpu.so ../../

tidy:
	$(RM) -f *~
clean: tidy
	$(RM) -rf datavectorcpu_wrap.c datavectorcpu.py datavectorcpu.pyc _datavectorcpu.so build
 
