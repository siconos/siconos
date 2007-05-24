Summary: Siconos/Front-End provides python and scilab interfaces for Siconos/Kernel
Name: siconos-front-end
Version: 2.0.1
Release: 1
License: GNU LGPL
Group: Application/Math
URL: http://gforge.inria.fr/projects/siconos
Source0: Siconos-Front-End-v%{version}.tgz
BuildRoot: %{_tmppath}/%{name}-v%{version}-%{release}-root

BuildPreReq: autoconf, automake, gcc, gcc-gfortran, doxygen, atlas, atlas-devel, cppunit, siconos-numerics, swig,scilab 
Requires: atlas, siconos-numerics, siconos-kernel, python, numpy, gnuplot, scilab

%description 
The present package, Siconos/Front-End provides python and scilab
interfaces for Siconos/Kernel.  More details are available on software
documentation pages, Doc/Devel directory of the current distribution,
in "design" chapter, or on http://siconos.gforge.inria.fr/

%define component Front-End
%define gdocs GeneratedDocs
%define docs %{_datadir}/doc/%{name}-%{version}

%prep
%setup -q -c
mkdir -p %{gdocs}/%{component}
mkdir -p %{gdocs}/Tags
pushd %{component}
./autogen.sh

%build
pushd %{component}
%{configure} --enable-scilab
%{__make}

%install
pushd %{component}
rm -rf %{buildroot}
make DESTDIR=%{buildroot} PYTHON_PATH=%{buildroot}/usr/site-package PYTHON_DIR=%{buildroot}/usr/site-package install
mkdir -p %{buildroot}%{docs}/%{component}
%{__install} AUTHORS COPYING ChangeLog NEWS README %{buildroot}%{docs}
popd
pushd %{gdocs}
tar cvf - %{component} | (cd %{buildroot}%{docs} ; tar xvf -)

%clean
rm -rf %{buildroot}

%files
%defattr(-,root,root,-)
%doc %{docs}
%{_bindir}
%{_usr}/site-package
%exclude %{_usrsrc}/debug

%changelog
* Fri May 11 2007 Maurice Bremond <Maurice.Bremond@inria.fr> - 2.0.0-1
- initial rpm
