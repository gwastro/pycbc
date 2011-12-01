Welcome to PYCBC

More information can be found at the website:
https://www.lsc-group.phys.uwm.edu/ligovirgo/cbcnote/InspiralPipelineDevelopment

When installing pycbc:
python setup.py install --prefix=/path/to/install/directory

To install and bypass the unit tests:
python setup.py install --force --prefix=/path/to/install/directory

When installing with opencl extensions:
python setup.py config --with-opencl=/path/to/opencl/toolkit/ --prefix=/path/to/install/directory

After the install to run the unit tests:
python setup.py test

To clean the pycbc src tree:
python setup.py clean

To clean and and remove the build tree:
python setup.py clean --all

The most current development version is always available from our
git repository:
https://sugwg-git.phy.syr.edu/git/pycbc/
