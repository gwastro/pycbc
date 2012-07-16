PKG_CONFIG = pkg-config
LALINCLUDE = $(shell $(PKG_CONFIG) --cflags lal)
#LALLIBS = $(shell $(PKG_CONFIG) --libs lal)
RM = /bin/rm
PYTHON = python
#SWIG = swig
SWIG = /home/jwillis/local/swig/bin/swig
CFLAGS += -std=gnu99

all: testlal

testlal_wrap.c testlal.py: testlal.i
	$(SWIG) -python -O -builtin $(LALINCLUDE) testlal.i

testlal: testlal_wrap.c testlal.py
	$(PYTHON) setup.py build_ext --inplace

clean:
	$(RM) -rf build testlal.pyc _testlal.so

distclean: clean
	$(RM) -rf testlal_wrap.c testlal.py