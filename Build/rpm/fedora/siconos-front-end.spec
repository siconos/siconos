Summary: Siconos/Front-End provides python and scilab interfaces for Siconos/Kernel
Name: siconos-front-end
Version: 2.1.1
Release: 1
License: GNU LGPL
Group: Application/Math
URL: http://gforge.inria.fr/projects/siconos
Source0: %{name}-%{version}.tar.gz
BuildRoot: %{_tmppath}/%{name}-%{version}-%{release}-root

BuildPreReq: autoconf, automake, gcc, gcc-gfortran, doxygen, atlas, atlas-devel, cppunit, siconos-numerics, swig,scilab 
Requires: atlas, siconos-numerics, siconos-kernel, python, numpy, gnuplot, scilab

%description 
The present package, Siconos/Front-End provides python and scilab
interfaces for Siconos/Kernel.  More details are available on software
documentation pages, Doc/Devel directory of the current distribution,
in "design" chapter, or on http://siconos.gforge.inria.fr/

%define component Front-End
%define namev %{name}-%{version}
%define gdocs GeneratedDocs
%define docs %{_datadir}/doc/siconos-%{version}

%prep
%setup -q -c
mkdir -p %{gdocs}/%{component}
mkdir -p %{gdocs}/Tags

%build
pushd %{namev}
%{configure} --enable-scilab
%{__make}

%install
pushd %{namev}
rm -rf %{buildroot}
make DESTDIR=%{buildroot} PYTHON_PATH=%{buildroot}/usr/site-package PYTHON_DIR=%{buildroot}/usr/site-package install
mkdir -p %{buildroot}%{docs}/%{component}
mkdir -p %{buildroot}%{docs}/html
mkdir -p %{buildroot}%{docs}/pdf
%{__install} AUTHORS COPYING ChangeLog NEWS README %{buildroot}%{docs}/%{component}
pdfs=`find . -name \*.pdf`; [ x"$pdfs" = x ] || cp $pdfs %{buildroot}%{docs}/pdf
popd
pushd %{gdocs}
rm -rf `find . -name latex -type d`
pdfs=`find . -name \*.pdf`; [ x"$pdfs" = x ] || cp $pdfs %{buildroot}%{docs}/pdf
cp -r %{component} %{buildroot}%{docs}/html

%clean
rm -rf %{buildroot}

%files
%defattr(-,root,root,-)
%doc %{docs}
%doc %{_mandir}/*/*
%{_bindir}
%{_usr}/site-package
%exclude %{_usrsrc}/debug

%changelog
* Fri May 11 2007 Maurice Bremond <Maurice.Bremond@inria.fr> - 2.0.0-1
- initial rpm
