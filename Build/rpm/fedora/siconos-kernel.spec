Summary: Siconos/Kernel, is dedicated to the modeling and the simulation of NSDS
Name: siconos-kernel
Version: 2.0.1
Release: 1
License: GNU LGPL
Group: Development/Libraries
URL: http://gforge.inria.fr/projects/siconos
Source0: Siconos-Kernel-v%{version}.tgz
BuildRoot: %{_tmppath}/%{name}-v%{version}-%{release}-root

BuildPreReq: autoconf, automake, gcc, gcc-gfortran, doxygen, atlas, atlas-devel, cppunit, siconos-numerics
Requires: atlas, siconos-numerics

%description 
The present package, Siconos/Kernel, is dedicated to the modeling and
the simulation of NSDS, with high level description and numerical
solving strategies.  It relies on Siconos/Numerics which provides
low-level algorithms to compute basic well-identified problems.

%define component Kernel
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
%{configure}
%{__make}
%{__make} doc
#%{__make} check

%install
pushd %{component}
rm -rf %{buildroot}
make DESTDIR=%{buildroot} install
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
%{_libdir}
%{_includedir}
%{_bindir}
%{_datadir}/%{name}
%exclude %{_libdir}/debug

%changelog
* Fri May 11 2007 Maurice Bremond <Maurice.Bremond@inria.fr> - 2.0.0-1
- initial rpm
