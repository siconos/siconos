Summary: Siconos/Kernel, is dedicated to the modeling and the simulation of NSDS
Name: siconos-kernel
Version: 2.1.0
Release: 1.fc6
License: GNU LGPL
Group: Development/Libraries
URL: http://gforge.inria.fr/projects/siconos
Source0: %{name}-%{version}.tar.gz
BuildRoot: %{_tmppath}/%{name}-%{version}-%{release}-root

BuildPreReq: autoconf, automake, gcc, doxygen, siconos-numerics, atlas-devel, boost-devel
Requires: atlas, boost, siconos-numerics

%description 
The present package, Siconos/Kernel, is dedicated to the modeling and
the simulation of NSDS, with high level description and numerical
solving strategies.  It relies on Siconos/Numerics which provides
low-level algorithms to compute basic well-identified problems.

%define component Kernel
%define namev %{name}-%{version}
%define gdocs GeneratedDocs
%define docs %{_datadir}/doc/siconos-%{version}

%prep
%setup -q -c
mkdir -p %{gdocs}/%{component}
mkdir -p %{gdocs}/Tags

%build
pushd %{namev}
%{configure} --enable-cppunit
%{__make}
%{__make} doc
%{__make} check

%install
pushd %{namev}
rm -rf %{buildroot}
make DESTDIR=%{buildroot} install
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
%{_libdir}
%{_includedir}
%{_bindir}
%{_datadir}/%{name}
%exclude %{_libdir}/debug

%changelog
* Fri May 11 2007 Maurice Bremond <Maurice.Bremond@inria.fr> - 2.0.0-1
- initial rpm
