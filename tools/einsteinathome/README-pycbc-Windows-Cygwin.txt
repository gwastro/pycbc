1. Setup a (virtual) machine with Win7 64 Bit (I used 2 cores, 3GB RAM, 180GB HD)
The name of the user account should not contain blanks/spaces

2. Install Cygwin 64 Bit Version from Cygwin.com (https://cygwin.com/setup-x86_64.exe)
I used the mirror http://ftp.inf.tu-dresden.de . Install (required):
- Archive: unzip, zip, hdf5
- Database: libmysqlclient-devel, libmysqlclient18
- Devel: gcc-core, gcc-g++, gcc-gfortran, autoconf (2.5+wrapper), automake (1.15+wrapper), make, pkg-config, binutils, cmake, libtool, python, git, openssl-devel, swig, libhdf5-devel, libhdf5_10, libhdf5cpp
- Libs: libjpeg-devel, libpng-devel, libpcre-devel, libpq-devel, libgsl-devel, libgsl0
- Web: wget
- Net: openssh
- Science: gsl
Recommended: Editors: emacs, Admin: cygrunsrv, shutdown
NOTE: unlike in e.g. Debian apt, -devel packages don't include / require the corresponding standard/runtime packages,
      so you need to select both separately (e.g gsl, hdf5)

3. After Cygwin installation finished, open a Cygwin shell and run a few commands:
python -m ensurepip
pip install --upgrade pip
pip install virtualenv
git config --user.email "your@email.adr"
git config --user.name "Your Name"

Note that for some git commands in the script you need to configure a (syntactically) valid email address.

4. Download and run the build script:
wget https://gitmaster.atlas.aei.uni-hannover.de/einsteinathome/lalsuite/blobs/raw/eah_cbc/lalapps/src/pulsar/EinsteinAtHome/PyCBC_EinsteinAtHome/pycbc_build_etch.sh
chmod +x pycbc_build_etch.sh
./pycbc_build_etch.sh



Some additional Hints for convenience:

1. Autologin: In a VM running on your local machine with only local networking login usually just holds you up.
Open a terminal ("cmd.exe") and type "control userpasswords2" to disable this.

2. If you really need case-sensitive filenames in Cygwin (so far I didn't), set registry value
HKLM\SYSTEM\CurrentControlSet\Control\Session Manager\kernel\obcaseinsensitive
to 0 and reboot the (virtual) machine.

3. Setup ssh server for remote command-line access:
- open a Cygwin terminal as Admin (right-click -> run as administrator)
- in that shell, run "ssh-host-config" and follow the instructions. You will probably run "cygrunsrv -S sshd" at the end
- Make sure the Windows firewall permits connection to port 22 from outside

4. if you get errors from fork of the type "Resource temporarily unavailable", do
- stop "Cygwin sshd" service (serch for "services")
- close all locally open Cygwin windows/terminals
- with Windows Explorer, navigate to C:\Cygwin64\bin, right-click on "ash" and "Run as Administrator"
- type "/bin/rebaseall"
- close the shell ("exit") and start "Cygwin sshd" again
