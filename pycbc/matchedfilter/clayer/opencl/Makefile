RM = /bin/rm
PYTHON = python
SWIG = swig

all: matchedfilteropencl.i matchedfilteropencl_types.h matchedfilteropencl_prototypes.h matchedfilteropencl.c
	$(SWIG) -python matchedfilteropencl.i
	$(PYTHON) setup.py build_ext --inplace
	cp matchedfilteropencl.py ../../
	cp _matchedfilteropencl.so ../../

tidy:
	$(RM) -f *~
clean: tidy
	$(RM) -rf matchedfilteropencl_wrap.c matchedfilteropencl.py ../../matchedfilteropencl.py \
        matchedfilteropencl.pyc _matchedfilteropencl.so ../../_matchedfilteropencl.so build
 
