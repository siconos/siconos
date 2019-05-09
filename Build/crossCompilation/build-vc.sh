#!/bin/bash


MXE_PREFIX="$HOME"
CMAKE_INSTALL_PREFIX="$HOME/siconos-windows-VC/install"
SICONOS_SOURCES="$HOME/workspace/crosscompilation"
BUILD_DIR="$HOME/siconos-windows-VC/build"
LIBS_PREFIX="$HOME/libs"
GCC_VER_MXE="4.8.0"
GCC_VER="4.7.0"

REL_TYPE="Debug" # or Release

# need fixes in Siconos first
build_component() {
	rm -rf ${BUILD_DIR}/$1
	mkdir -p ${BUILD_DIR}/$1
	cd ${BUILD_DIR}/$1
	cmake -DBUILD_STATIC_LIBS=ON \
		-DCMAKE_TOOLCHAIN_FILE=${MXE_PREFIX}/mxe/usr/i686-pc-mingw32/share/cmake/mxe-conf.cmake \
		-DCROSSCOMPILING_LINUX_TO_WINDOWS=1 \
		-DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX} \
		${SICONOS_SOURCES}/$1
}

make_component() {
	make ExperimentalStart
#	make ExperimentalConfigure # not working right now -- xhub
	make ExperimentalBuild
	make ExperimentalTest
	make ExperimentalSubmit
}

build_siconos() {
	rm -rf ${BUILD_DIR}/*
	rm -rf ${CMAKE_INSTALL_PREFIX}/*
  build_Numerics_VS
	build_Kernel_VS
  build_FE_VS
}

build_Numerics_VS() {
	mkdir -p "${BUILD_DIR}/Numerics"
	cd "${BUILD_DIR}/Numerics"

	CXXFLAGS='/DBUILD_AS_CPP /DADD_ /DHAVE_LAPACK_CONFIG_H /DLAPACK_COMPLEX_STRUCTURE' cmake \
	-DCMAKE_TOOLCHAIN_FILE=$HOME/siconos-windows-VC/toolchain.cmake \
	-DCROSSCOMPILING_LINUX_TO_WINDOWS=1 \
	-DLAPACK_LIBRARIES="$HOME/mxe/usr/i686-pc-mingw32/lib/liblapack.a;$HOME/mxe/usr/i686-pc-mingw32/lib/liblapacke.a" \
        -DLAPACK_INCLUDE_DIRS=$HOME/mxe/usr/i686-pc-mingw32/include-VC \
	-DLAPACK_HEADER="lapacke.h" \
	-DHAS_LAPACKE=1 \
        -DBLAS_INCLUDE_DIRS=$HOME/mxe/usr/i686-pc-mingw32/include-VC \
	-DBLAS_LIBRARIES="$HOME/mxe/usr/i686-pc-mingw32/lib/libblas.a;${WINEPREFIX}/drive_c/MinGW/lib/gcc/mingw32/${GCC_VER}/libgfortran.dll.a;$HOME/mxe/usr/i686-pc-mingw32/lib/libcblas.a" \
	-DCMAKE_INSTALL_PREFIX=$HOME/siconos-windows-VC/install \
	-DCMAKE_BUILD_TYPE=${REL_TYPE} \
	${SICONOS_SOURCES}/Numerics
	make -k -j2
	make ExperimentalTest
	make -i install -j2
	make ExperimentalSubmit
}

build_Kernel_VS() {
	mkdir -p "${BUILD_DIR}/Kernel"
	cd "${BUILD_DIR}/Kernel"
	CXXFLAGS='/DADD_ /DHAVE_LAPACK_CONFIG_H /DLAPACK_COMPLEX_STRUCTURE /DBOOST_ALL_NO_LIB'  cmake \
		-DCMAKE_TOOLCHAIN_FILE=$HOME/siconos-windows-VC/toolchain.cmake \
		-DCROSSCOMPILING_LINUX_TO_WINDOWS=1 \
		-DLAPACK_LIBRARIES="$HOME/mxe/usr/i686-pc-mingw32/lib/liblapack.a;$HOME/mxe/usr/i686-pc-mingw32/lib/liblapacke.a" \
        	-DLAPACK_INCLUDE_DIRS=$HOME/mxe/usr/i686-pc-mingw32/include-VC \
		-DLAPACK_HEADER="lapacke.h" \
		-DHAS_LAPACKE=1 \
	        -DBLAS_INCLUDE_DIRS=$HOME/mxe/usr/i686-pc-mingw32/include-VC \
		-DBLAS_LIBRARIES="$HOME/mxe/usr/i686-pc-mingw32/lib/libblas.a;${WINEPREFIX}/drive_c/MinGW/lib/gcc/mingw32/${GCC_VER}/libgfortran.dll.a;$HOME/mxe/usr/i686-pc-mingw32/lib/libcblas.a" \
		-DCMAKE_INSTALL_PREFIX=$HOME/siconos-windows-VC/install \
		-DSiconosNumerics_LIBRARY="${CMAKE_INSTALL_PREFIX}/lib/libSiconosNumerics.dll.a" \
		-DSiconosNumerics_INCLUDE_DIRS="${CMAKE_INSTALL_PREFIX}/include/Siconos/Numerics/" \
		-DGMP_LIBRARY=$HOME/mxe/usr/i686-pc-mingw32/lib/libgmp.a \
		-DGMP_LIBRARIES=$HOME/mxe/usr/i686-pc-mingw32/lib/libgmp.a \
		-DGMP_INCLUDE_DIRS=$HOME/mxe/usr/i686-pc-mingw32/include-VC \
		-DGMP_INCLUDE_DIR=$HOME/mxe/usr/i686-pc-mingw32/include-VC \
		-DCPPUNIT_LIBRARIES=$HOME/libs/cppunit-1.12.1/src/cppunit/${REL_TYPE}/cppunit.lib \
		-DCPPUNIT_INCLUDE_DIR=$HOME/mxe/usr/i686-pc-mingw32/include-VC/cppunit \
		-DCPPUNIT_INCLUDE_DIRS=$HOME/mxe/usr/i686-pc-mingw32/include-VC/cppunit \
		-DCPPUNIT_LIBRARY=$HOME/libs/cppunit-1.12.1/src/cppunit/${REL_TYPE}/cppunit.lib \
		-DBoost_INCLUDE_DIR=$HOME//mxe/usr/i686-pc-mingw32/include-VC/ \
		-DCMAKE_BUILD_TYPE=${REL_TYPE} \
		${SICONOS_SOURCES}/Kernel/
	# VS is way too verbose right now
	make -k -j3 2>&1 | egrep '( error |Building|Built|Linking)'
	make ExperimentalTest
	make -i install 2>&1 | egrep '( error |Building|Built|Linking)'
	make ExperimentalSubmit
}

build_FE_VS() {
	mkdir -p "${BUILD_DIR}/Front-End"
	cd "${BUILD_DIR}/Front-End"
	CXXFLAGS='/DADD_ /DHAVE_LAPACK_CONFIG_H /DLAPACK_COMPLEX_STRUCTURE'  cmake \
		-DCMAKE_TOOLCHAIN_FILE=$HOME/siconos-windows-VC/toolchain.cmake \
		-DCROSSCOMPILING_LINUX_TO_WINDOWS=1 \
        	-DLAPACK_INCLUDE_DIRS=$HOME/mxe/usr/i686-pc-mingw32/include-VC \
		-DLAPACK_HEADER="lapacke.h" \
		-DHAS_LAPACKE=1 \
	        -DBLAS_INCLUDE_DIRS=$HOME/mxe/usr/i686-pc-mingw32/include-VC \
		-DLAPACK_LIBRARIES="$HOME/mxe/usr/i686-pc-mingw32/lib/liblapack.a;$HOME/mxe/usr/i686-pc-mingw32/lib/liblapacke.a" \
		-DBLAS_LIBRARIES="$HOME/mxe/usr/i686-pc-mingw32/lib/libblas.a;${WINEPREFIX}/drive_c/MinGW/lib/gcc/mingw32/${GCC_VER}/libgfortran.dll.a;$HOME/mxe/usr/i686-pc-mingw32/lib/libcblas.a" \
		-DCMAKE_INSTALL_PREFIX=$HOME/siconos-windows-VC/install \
		-DSiconosNumerics_LIBRARY="${CMAKE_INSTALL_PREFIX}/lib/libSiconosNumerics.dll.a" \
		-DSiconosNumerics_INCLUDE_DIRS="${CMAKE_INSTALL_PREFIX}/include/Siconos/Numerics/" \
		-DSiconosKernel_INCLUDE_DIRS="${CMAKE_INSTALL_PREFIX}/include/Siconos/Kernel/" \
		-DSiconosKernel_LIBRARY="${CMAKE_INSTALL_PREFIX}/lib/libSiconosKernel.dll.a" \
		-DSWIG_DIR="/usr/share/swig/2.0.9/" \
		-DPYTHON_LIBRARIES="${LIBS_PREFIX}/python/python27.lib" \
		-DPYTHON_LIBRARY="${LIBS_PREFIX}/python/python27.lib" \
		-DPYTHON_INCLUDE_DIR="${LIBS_PREFIX}/python/" \
		-DPYTHON_NUMPY_INCLUDE_DIR="${LIBS_PREFIX}/numpy/PLATLIB/numpy/core/include" \
		-DBOOST_ROOT="${MXE_PREFIX}/mxe/usr/i686-pc-mingw32/" \
		-DBoost_INCLUDE_DIR=$HOME/mxe/usr/i686-pc-mingw32/include-VC/ \
		-DGMP_LIBRARY=$HOME/mxe/usr/i686-pc-mingw32/lib/libgmp.a \
		-DGMP_LIBRARIES=$HOME/mxe/usr/i686-pc-mingw32/lib/libgmp.a \
		-DGMP_INCLUDE_DIRS=$HOME/mxe/usr/i686-pc-mingw32/include-VC \
		-DGMP_INCLUDE_DIR=$HOME/mxe/usr/i686-pc-mingw32/include-VC \
		-DCMAKE_BUILD_TYPE=Release \
		${SICONOS_SOURCES}/Front-End

	make -j3
	make ExperimentalTest
	make -i install
	make ExperimentalSubmit

}

. $HOME/msvc/env.sh

#dumb access to VS compiler, the first execution is always so slow ...
cl.exe
build_siconos
