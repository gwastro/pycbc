import os, distutils.sysconfig, sys, os.path
import scipy.misc
import scipy, fnmatch

def find(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result

# Handle stupid scipy imports (fixed in scipy 0.13)
scipy.factorial = scipy.misc.factorial

# handle python shared libary name
basedir = sys._MEIPASS
print("Setting up scipy for temporary basedir %s" % basedir)
python_lib = find('libpython*.so.*', basedir)[0]
python_lib_dest = python_lib.split('.so')[0] + '.so'
os.symlink(python_lib, python_lib_dest)

print("LD_LIBRARY_PATH=%s" % os.environ['LD_LIBRARY_PATH'])
