#!/bin/bash
export PKG_CONFIG_PATH=$HOME/dev/root/lalsuite/branches/newswig/lib/pkgconfig
export LD_LIBRARY_PATH=$HOME/dev/root/lalsuite/branches/newswig/lib
export PYTHONPATH=$HOME/dev/root/pycbc/branches/master/lib/python2.7/site-packages:$HOME/dev/root/lalsuite/branches/newswig/lib/python2.7/site-packages:$PYTHONPATH
export CFLAGS="-std=gnu99"