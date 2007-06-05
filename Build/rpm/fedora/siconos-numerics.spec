Summary: The Siconos/Numerics Package is dedicated to low-level algorithms
Name: siconos-numerics
Version: 2.0.1
Release: 1
License: GNU LGPL
Group: Development/Libraries
URL: http://gforge.inria.fr/projects/siconos
Source0: %{name}-%{version}.tar.gz
BuildRoot: %{_tmppath}/%{name}-%{version}-%{release}-root
BuildPreReq: autoconf, automake, gcc, gcc-gfortran, doxygen, atlas, atlas-devel, cppunit
Requires: atlas

%description 
Siconos is a program dedicated to modeling, simulation and control of
non smooth dynamical systems. The present package, Siconos/Numerics
provides low-level algorithms to compute basic well-identified
problems.

%define component Numerics
%define namev %{name}-%{version}
%define gdocs GeneratedDocs
%define docs %{_datadir}/doc/siconos-%{version}

%prep
%setup -q -c
mkdir -p %{gdocs}/%{component}
mkdir -p %{gdocs}/Tags
pushd %{namev}

%build
pushd %{namev}
%{configure}
%{__make}
%{__make} doc
%{__make} check

%install
pushd %{namev}
rm -rf %{buildroot}
make DESTDIR=%{buildroot} install
mkdir -p %{buildroot}%{docs}/%{component}
mkdir -p %{buildroot}%{docs}/html
%{__install} AUTHORS COPYING ChangeLog NEWS README %{buildroot}%{docs}/%{component}
popd
pushd %{gdocs}
cp -r %{component} %{buildroot}%{docs}/html

%clean
rm -rf %{buildroot}

%files
%defattr(-,root,root,-)
%doc %{docs}
%{_libdir}
%{_includedir}
%exclude %{_libdir}/debug

%changelog
* Fri May 11 2007 Maurice Bremond <Maurice.Bremond@inria.fr> - 2.0.0-1
- initial rpm
