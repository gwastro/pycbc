Name:           python-pycbc
Summary:        Python tools for CBC data analysis
Version:        0.1
Release:        1
License:        GPL
Group:          Development/Libraries
Source:         %{name}-%{version}.tar.gz
Url:            http://www.lsc-group.phys.uwm.edu/daswg/projects/pycbc.html
BuildRoot:      %{_tmppath}/%{name}-%{version}-root
Requires:       python python-decorator python-pylal glue glue-segments lal lal-python lalframe lalframe-python lalsimulation lalsimulation-python lalinspiral lalinspiral-python numpy scipy
BuildRequires:  python-devel lal-devel lalmetaio-devel lalframe-devel lalsimulation-devel lalinspiral-devel numpy pkgconfig
%description
PyCBC is a python toolkit for analysis of data from gravitational-wave
laser interferometer detectors with the goal of detecting and studying
signals from compact binary coalescences (CBCs).

%package common
Summary: Common files
Group: Development/Libraries
Requires: python
%description common
This is the common package, to be complemented by one of the
-nogpu, -cuda or -opencl packages.

%package nogpu
Summary: CPU-only version
Group: Development/Libraries
Requires: python-pycbc-common
%description nogpu
Version supporting CPU computation only.

#%package cuda
#Summary: CUDA version
#Group: Development/Libraries
#Requires: python-pycbc-common
#%description cuda
#Version supporting GPU computation via CUDA.

#%package opencl
#Summary: OpenCL version
#Group: Development/Libraries
#Requires: python-pycbc-common
#%description opencl
#Version supporting GPU computation via OpenCL.

%prep
%setup

%build
CFLAGS="%{optflags}" %{__python} setup.py build

%install
rm -rf %{buildroot}
%{__python} setup.py install \
    --skip-build \
    --root=%{buildroot}

%files common
%defattr(-,root,root)
%exclude /usr/etc/
/usr/

%files nogpu
