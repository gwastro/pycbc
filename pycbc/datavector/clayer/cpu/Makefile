RM = /bin/rm
PYTHON = /usr/bin/python
SWIG = /usr/bin/swig

all: datavectorcpu.i datavectorcpu.h datavectorcpu.c
	$(SWIG) -python datavectorcpu.i
	$(PYTHON) setup.py build_ext --inplace
tidy:
	$(RM) -f *~
clean: tidy
	$(RM) -rf datavectorcpu_wrap.c datavectorcpu.py datavectorcpu.pyc _datavectorcpu.so /build
 
