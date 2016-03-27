Name:           siconos
Version:        3.9.0
Release:        2%{?dist}
Summary:        Simulation platform dedicated to non-smooth dynamical systems

License:        GPLv3
URL:            http://siconos.gforge.inria.fr
#Source0:        https://github.com/siconos/siconos/archive/master.tar.gz
Source0:        master.tar.gz

BuildRequires:  gcc-gfortran
BuildRequires:  gcc-c++
BuildRequires:	boost-devel
BuildRequires:	boost-filesystem
BuildRequires:	boost-serialization
BuildRequires:	gmp-devel
BuildRequires:	cppunit-devel
BuildRequires:	hdf5-devel
BuildRequires:	openblas-devel
BuildRequires:	lapack-devel
BuildRequires:  cmake

#BuildRequires:	ptscotch-openmpi
#BuildRequires:	ptscotch-openmpi-devel
#BuildRequires:	blacs-openmpi-devel
#BuildRequires:	openmpi-devel
#BuildRequires:	scalapack-openmpi
#BuildRequires:	MUMPS-openmpi-devel
#BuildRequires:	environment-modules

BuildRequires: MUMPS-devel

BuildRequires:  suitesparse-devel
BuildRequires:  python-devel
BuildRequires:  numpy
BuildRequires:  swig
BuildRequires:  lpsolve-devel
BuildRequires:  bullet-devel

Requires:       gcc-c++
Requires:       boost-devel
Requires:       h5py
Requires:       scipy
Requires:       python-lxml
Requires:       python-matplotlib
Requires:       cmake
Requires:       gmp-devel
Requires:       vtk

#Requires:	openmpi%{?_isa}

%if 0%{?fedora}
#BuildRequires: rpm-mpi-hooks
# for devel
#Requires: rpm-mpi-hooks
%endif

%description
Siconos is a simulation platform dedicated to non-smooth dynamical systems

%prep
%autosetup -n %{name}-master
#%setup -q


%build
#export CC=mpicc
#export CXX=mpicxx
#export FC=mpif90
#export F77=mpif77

#%{_openmpi_load}

mkdir build
cd build
%cmake -DFORCE_SKIP_RPATH=1 -DNO_RUNTIME_BUILD_DEP=1 -DWITH_MUMPS=1 -DWITH_HDF5=1 -DWITH_BULLET=1 -DIDONTWANTMPI=1 -DWITH_UMFPACK=1 ..
make %{?_smp_mflags}


%install
rm -rf $RPM_BUILD_ROOT
cd build
%make_install


%files
#%doc README
#%license LICENSE COPYING
%{_bindir}/*
%{_libdir}/*
%{python2_sitelib}/*
%{_includedir}/%{name}/
%{_datadir}/%{name}/


%post -p /sbin/ldconfig

%postun -p /sbin/ldconfig


%changelog
* Sun Mar 27 2016 Olivier Huber <oli.huber@gmail.com> - 3.9.0-2
- revdump for update dependencies

* Sat Mar 26 2016 Olivier Huber <oli.huberATgmail.com> - 3.9.0-1
- initial test
