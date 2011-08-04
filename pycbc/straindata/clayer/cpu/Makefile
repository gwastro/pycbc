RM = /bin/rm
PYTHON = python
SWIG = swig

all: straindatacpu.i straindatacpu_types.h straindatacpu_prototypes.h straindatacpu.c
	$(SWIG) -python straindatacpu.i
	$(PYTHON) setup.py build_ext --inplace
	cp straindatacpu.py ../../
	cp _straindatacpu.so ../../

tidy:
	$(RM) -f *~
clean: tidy
	$(RM) -rf straindatacpu_wrap.c straindatacpu.py ../../straindatacpu.py \
        straindatacpu.pyc _straindatacpu.so ../../_straindatacpu.so build
 
