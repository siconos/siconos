#!/bin/bash

set -x

MXE_PREFIX="$HOME"
CMAKE_INSTALL_PREFIX="$HOME/siconos-windows-VC/install"
SICONOS_SOURCES="$HOME/s-VC"
BUILD_DIR="$HOME/siconos-windows-VC/build"
LIBS_PREFIX="$HOME/libs"
GCC_VER_MXE="4.8.0"
GCC_VER="4.7.0"

if [ -z $1 ]; then
  REL_TYPE="Debug" # or Release
else
  REL_TYPE=$1
fi

build_siconos() {
	rm -rf ${BUILD_DIR}/*
	rm -rf ${CMAKE_INSTALL_PREFIX}/*
	cd ${BUILD_DIR}


	cmake "${SICONOS_SOURCES}/Build" \
	-DCMAKE_TOOLCHAIN_FILE=$HOME/siconos-windows-VC/toolchain.cmake \
	-DON_DASHBOARD=1 \
	-DMODE=Continuous \
	-DWITH_IO=0 \
	-DWITH_MECHANICS=1 \
	-DWITH_CONTROL=1 \
	-DWITH_EXAMPLES=1  \
	-DCMAKE_INSTALL_PREFIX="${CMAKE_INSTALL_PREFIX}" \
	-DINSTALL_COMMAND="make;-ik;install" \
	-DFROM_REPO=0 \
	-DCTEST_OPTIONS="-j2;-V" \
	-DBUILD_TYPE=$REL_TYPE \
        -DCROSSCOMPILING_LINUX_TO_WINDOWS=1 \
	-Dproject_CMAKE_ARGS="-DCMAKE_TOOLCHAIN_FILE=$HOME/siconos-windows-VC/toolchain.cmake;-DCROSSCOMPILING_LINUX_TO_WINDOWS=1;-DLAPACK_INCLUDE_DIRS=$HOME/mxe/usr/i686-pc-mingw32/include-VC;-DLAPACK_HEADER=lapacke.h;-DHAS_LAPACKE=1;-DBLAS_INCLUDE_DIRS=$HOME/mxe/usr/i686-pc-mingw32/include-VC;-DLAPACK_LIBRARIES=$HOME/mxe/usr/i686-pc-mingw32/lib/liblapack.a___$HOME/mxe/usr/i686-pc-mingw32/lib/liblapacke.a;-DEXTRA_LAPACK_LIB=$HOME/libs/siconosCompat.a;-DBLAS_LIBRARIES=$HOME/mxe/usr/i686-pc-mingw32/lib/libblas.a___${WINEPREFIX}/drive_c/MinGW/lib/gcc/mingw32/${GCC_VER}/libgfortran.dll.a___$HOME/mxe/usr/i686-pc-mingw32/lib/libcblas.a;-DCMAKE_INSTALL_PREFIX=$HOME/siconos-windows-VC/install;-DSiconosNumerics_LIBRARY=${CMAKE_INSTALL_PREFIX}/lib/libSiconosNumerics.dll.a;-DSiconosNumerics_INCLUDE_DIRS=${CMAKE_INSTALL_PREFIX}/include/Siconos/Numerics/;-DSiconosKernel_INCLUDE_DIRS=${CMAKE_INSTALL_PREFIX}/include/Siconos/Kernel/;-DSiconosKernel_LIBRARY=${CMAKE_INSTALL_PREFIX}/lib/libSiconosKernel.dll.a;-DSWIG_DIR=/usr/share/swig/2.0.9/;-DPYTHON_LIBRARIES=${LIBS_PREFIX}/python/python27.lib;-DPYTHON_LIBRARY=${LIBS_PREFIX}/python/python27.lib;-DPYTHON_INCLUDE_DIR=${LIBS_PREFIX}/python/;-DPYTHON_NUMPY_INCLUDE_DIR=${LIBS_PREFIX}/numpy/PLATLIB/numpy/core/include;-DBOOST_ROOT=${MXE_PREFIX}/mxe/usr/i686-pc-mingw32/;-DBoost_INCLUDE_DIR=$HOME/mxe/usr/i686-pc-mingw32/include-VC/;-DGMP_LIBRARY=$HOME/mxe/usr/i686-pc-mingw32/lib/libgmp.a;-DGMP_LIBRARIES=$HOME/mxe/usr/i686-pc-mingw32/lib/libgmp.a___$HOME/mxe/usr/lib/gcc/i686-pc-mingw32/4.8.0/libgcc.a;-DGMP_INCLUDE_DIRS=$HOME/mxe/usr/i686-pc-mingw32/include-VC;-DGMP_INCLUDE_DIR=$HOME/mxe/usr/i686-pc-mingw32/include-VC;-DCPPUNIT_INCLUDE_DIR=$HOME/mxe/usr/i686-pc-mingw32/include-VC/cppunit;-DCPPUNIT_LIBRARY=$HOME/libs/cppunit-1.12.1/src/cppunit/${REL_TYPE}/cppunit.lib;-DSiconosKernel_EXE_DIR=${CMAKE_INSTALL_PREFIX}/bin"

	CXXFLAGS='/DBUILD_AS_CPP /DADD_ /DHAVE_LAPACK_CONFIG_H /DLAPACK_COMPLEX_STRUCTURE /DBOOST_ALL_NO_LIB' make VERBOSE=1 -ik
}

. $HOME/msvc/env.sh

#dumb access to VS compiler, the first execution is always so slow ...
cl.exe
build_siconos
