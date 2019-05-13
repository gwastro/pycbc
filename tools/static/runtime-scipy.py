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

# add the libaries to the include and library directories
import weave.inline_tools
weave.inline_tools.inline_real = weave.inline_tools.inline
def inline(code,arg_names=None,local_dict=None, global_dict=None,
           force=0,
           compiler='',
           verbose=0,
           support_code=None,
           headers=None,
           customize=None,
           type_converters=None,
           auto_downcast=1,
           newarr_converter=0,
           **kw):
           if arg_names is None:
                      arg_names = []
           if headers is None:
                      headers = []
           kw['library_dirs'] = [basedir]
           kw['include_dirs'] = [os.path.join(basedir, 'include')]
           call_frame = sys._getframe().f_back
           if local_dict is None:
               local_dict = call_frame.f_locals
           if global_dict is None:
               global_dict = call_frame.f_globals

           return weave.inline_tools.inline_real(code,
                   arg_names, local_dict=local_dict, global_dict=global_dict,
                   force=force,
                   compiler=compiler,
                   verbose=verbose,
                   support_code=support_code,
                   headers=headers,
                   customize=customize,
                   type_converters=type_converters,
                   auto_downcast=auto_downcast,
                   newarr_converter=newarr_converter,
                   **kw)
weave.inline_tools.inline = inline
weave.inline = inline



print("LD_LIBRARY_PATH=%s" % os.environ['LD_LIBRARY_PATH'])
