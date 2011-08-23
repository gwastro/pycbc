RM = /bin/rm
PYTHON = python
SWIG = swig

all: templatebankcpu.i templatebankcpu_prototypes.h templatebankcpu.c
	$(SWIG) -python templatebankcpu.i
	$(PYTHON) setup.py build_ext --inplace
	cp templatebankcpu.py ../../
	cp _templatebankcpu.so ../../

tidy:
	$(RM) -f *~
clean: tidy
	$(RM) -rf templatebankcpu_wrap.c templatebankcpu.py ../../templatebankcpu.py \
        templatebankcpu.pyc _templatebankcpu.so ../../_templatebankcpu.so build
 
