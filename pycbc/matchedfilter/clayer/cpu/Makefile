RM = /bin/rm
PYTHON = python
SWIG = swig

all: matchedfiltercpu.i matchedfiltercpu_types.h matchedfiltercpu_prototypes.h matchedfiltercpu.c
	$(SWIG) -python matchedfiltercpu.i
	$(PYTHON) setup.py build_ext --inplace
	cp matchedfiltercpu.py ../../
	cp _matchedfiltercpu.so ../../

tidy:
	$(RM) -f *~
clean: tidy
	$(RM) -rf matchedfiltercpu_wrap.c matchedfiltercpu.py matchedfiltercpu.pyc _matchedfiltercpu.so build
 
