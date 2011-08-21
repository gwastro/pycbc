RM = /bin/rm
PYTHON = python
SWIG = swig

all: straindataopencl.i straindataopencl_types.h straindataopencl_prototypes.h straindataopencl.c
	$(SWIG) -python straindataopencl.i
	$(PYTHON) setup.py build_ext --inplace
	cp straindataopencl.py ../../
	cp _straindataopencl.so ../../

tidy:
	$(RM) -f *~
clean: tidy
	$(RM) -rf straindataopencl_wrap.c straindataopencl.py ../../straindataopencl.py \
        straindataopencl.pyc _straindataopencl.so ../../_straindataopencl.so build
 
