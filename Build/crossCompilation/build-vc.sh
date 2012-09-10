#!/bin/sh

MXE_PREFIX="/scratch/Olivier/mingw32/"
INSTALL_PREFIX="/scratch/Olivier/siconos-windows-VC/install"
SICONOS_SOURCES="$HOME/siconos"
BUILD_DIR="/scratch/Olivier/siconos-windows-VC/build"

REL_TYPE="Debug" # or Release

# need fixes in Siconos first
build_component() {
	rm -rf ${BUILD_DIR}/$1
	mkdir -p ${BUILD_DIR}/$1
	cd ${BUILD_DIR}/$1
	cmake -DBUILD_STATIC_LIBS=ON \
		-DCMAKE_TOOLCHAIN_FILE=${MXE_PREFIX}/mxe/usr/i686-pc-mingw32/share/cmake/mxe-conf.cmake \
		-DCROSSCOMPILING_LINUX_TO_WINDOWS=1 \
		-DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX} \
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
	rm -rf ${INSTALL_PREFIX}/*
  build_Numerics_VS
	build_Kernel_VS
  build_FE_VS
}

build_Numerics_VS() {
	mkdir -p "${BUILD_DIR}/Numerics"
	cd "${BUILD_DIR}/Numerics"

	CXXFLAGS='/DBUILD_AS_CPP' cmake \
	-DCMAKE_TOOLCHAIN_FILE=/scratch/Olivier/siconos-windows-VC/toolchain.cmake \
	-DCROSSCOMPILING_LINUX_TO_WINDOWS=1 \
	-DLAPACK_LIBRARIES=/scratch/Olivier/mingw32/mxe/usr/i686-pc-mingw32/lib/liblapack.lib \
	-DBLAS_LIBRARIES="/scratch/Olivier/mingw32/mxe/usr/i686-pc-mingw32/lib/libblas.lib;${WINEPREFIX}/drive_c/MinGW/lib/gcc/mingw32/4.7.0/libgfortran.dll.a" \
	-DCMAKE_INSTALL_PREFIX=/scratch/Olivier/siconos-windows-VC/install \
	-DCMAKE_BUILD_TYPE=${REL_TYPE} \
	${SICONOS_SOURCES}/Numerics
	make -k -j5
	make ExperimentalTest
	make -i install
	make ExperimentalSubmit
}

build_Kernel_VS() {
	mkdir -p "${BUILD_DIR}/Kernel"
	cd "${BUILD_DIR}/Kernel"
  cmake \
		-DCMAKE_TOOLCHAIN_FILE=/scratch/Olivier/siconos-windows-VC/toolchain.cmake \
		-DCROSSCOMPILING_LINUX_TO_WINDOWS=1 \
		-DLAPACK_LIBRARIES=/scratch/Olivier/mingw32/mxe/usr/i686-pc-mingw32/lib/liblapack.lib \
		-DBLAS_LIBRARIES="/scratch/Olivier/mingw32/mxe/usr/i686-pc-mingw32/lib/libblas.lib;${WINEPREFIX}/drive_c/MinGW/lib/gcc/mingw32/4.7.0/libgfortran.dll.a" \
		-DCMAKE_INSTALL_PREFIX=/scratch/Olivier/siconos-windows-VC/install \
		-DSiconosNumerics_FOUND="${INSTALL_PREFIX}/lib/libSiconosNumerics.dll.a" \
		-DSiconosNumerics_INCLUDE_DIRS="${INSTALL_PREFIX}/include/Siconos/Numerics/" \
		-DLIBXML2_LIBRARIES="/scratch/Olivier/mingw32/mxe/usr/i686-pc-mingw32/lib/libxml2.a;/scratch/Olivier/mingw32/mxe/usr/i686-pc-mingw32/lib/libiconv.a;/scratch/Olivier/mingw32/mxe/usr/i686-pc-mingw32/lib/libz.a;/scratch/Olivier/mingw32/mxe/usr/i686-pc-mingw32/lib/libws2_32.a;/scratch/Olivier/mingw32/mxe/usr/lib/gcc/i686-pc-mingw32/4.7.0/libgcc.a" \
		-DLIBXML2_INCLUDE_DIR=${MXE_PREFIX}/mxe/usr/i686-pc-mingw32/include-VC/libxml2 \
		-DGMP_FOUND=/scratch/Olivier/mingw32/mxe/usr/i686-pc-mingw32/lib/libgmp.a \
		-DGMP_INCLUDE_DIRS=/scratch/Olivier/mingw32/mxe/usr/i686-pc-mingw32/include-VC \
		-DCPPUNIT_FOUND=/scratch/Olivier/mingw32/cppunit-1.12.1/src/cppunit/${REL_TYPE}/cppunit.lib \
		-DCPPUNIT_INCLUDE_DIR=/scratch/Olivier/mingw32/mxe/usr/i686-pc-mingw32/include-VC/cppunit \
		-DCPPUNIT_LIBRARIES=/scratch/Olivier/mingw32/cppunit-1.12.1/src/cppunit/${REL_TYPE}/cppunit.lib \
		-DBoost_INCLUDE_DIR=/scratch/Olivier/mingw32/boost_1_50_0 \
		-DCMAKE_BUILD_TYPE=${REL_TYPE} \
		${SICONOS_SOURCES}/Kernel/
	# VS is way too verbose right now
	make -k -j5 2>&1 | egrep '( error |Building|Built|Linking)'
	make ExperimentalTest
	make -i install
	make ExperimentalSubmit
}

build_FE_VS() {
	mkdir -p "${BUILD_DIR}/Front-End"
	cd "${BUILD_DIR}/Front-End"
  cmake \
		-DCMAKE_TOOLCHAIN_FILE=/scratch/Olivier/siconos-windows-VC/toolchain.cmake \
		-DCROSSCOMPILING_LINUX_TO_WINDOWS=1 \
		-DLAPACK_LIBRARIES=/scratch/Olivier/mingw32/mxe/usr/i686-pc-mingw32/lib/liblapack.lib \
		-DBLAS_LIBRARIES="/scratch/Olivier/mingw32/mxe/usr/i686-pc-mingw32/lib/libblas.lib;${WINEPREFIX}/drive_c/MinGW/lib/gcc/mingw32/4.7.0/libgfortran.dll.a" \
		-DCMAKE_INSTALL_PREFIX=/scratch/Olivier/siconos-windows-VC/install \
		-DSiconosNumerics_FOUND="${INSTALL_PREFIX}/lib/libSiconosNumerics.dll.a" \
		-DSiconosNumerics_INCLUDE_DIRS="${INSTALL_PREFIX}/include/Siconos/Numerics/" \
		-DSiconosKernel_INCLUDE_DIRS="${INSTALL_PREFIX}/include/Siconos/Kernel/" \
		-DSiconosKernel_FOUND="${INSTALL_PREFIX}/lib/libSiconosKernel.dll.a" \
		-DSWIG_DIR="/usr/share/swig/2.0.7/" \
		-DPYTHON_LIBRARIES="${MXE_PREFIX}/python/python27.lib" \
		-DPYTHON_INCLUDE_DIR="${MXE_PREFIX}/python/" \
		-DPYTHON_NUMPY_INCLUDE_DIR="${MXE_PREFIX}/numpy/PLATLIB/numpy/core/include" \
		-DBOOST_ROOT="${MXE_PREFIX}/mxe/usr/i686-pc-mingw32/" \
		-DBoost_INCLUDE_DIR=/scratch/Olivier/mingw32/boost_1_50_0 \
		-DLIBXML2_LIBRARIES="/scratch/Olivier/mingw32/mxe/usr/i686-pc-mingw32/lib/libxml2.a;/scratch/Olivier/mingw32/mxe/usr/i686-pc-mingw32/lib/libiconv.a;/scratch/Olivier/mingw32/mxe/usr/i686-pc-mingw32/lib/libz.a;/scratch/Olivier/mingw32/mxe/usr/i686-pc-mingw32/lib/libws2_32.a;/scratch/Olivier/mingw32/mxe/usr/lib/gcc/i686-pc-mingw32/4.7.0/libgcc.a" \
		-DLIBXML2_INCLUDE_DIR=${MXE_PREFIX}/mxe/usr/i686-pc-mingw32/include-VC/libxml2 \
		-DGMP_FOUND=/scratch/Olivier/mingw32/mxe/usr/i686-pc-mingw32/lib/libgmp.a \
		-DGMP_INCLUDE_DIRS=/scratch/Olivier/mingw32/mxe/usr/i686-pc-mingw32/include-VC \
		-DCMAKE_BUILD_TYPE=Release \
		${SICONOS_SOURCES}/Front-End

	make -j5
	make ExperimentalTest
	make -i install
	make ExperimentalSubmit

}

source /scratch/Olivier/mingw32/msvc/env.sh
build_siconos

cd ${SICONOS_SOURCES}
git pull
CURRENT_SHA=`git rev-parse HEAD`
if [ -f /tmp/.sha-siconos-windows.old ]; then
	OLD_SHA=`cat /tmp/.sha-siconos-windows.old`
else
	OLD_SHA=""
fi

if [ ${CURRENT_SHA} != ${OLD_SHA} ]; then
	build_siconos
	echo ${CURRENT_SHA} > /tmp/.sha-siconos-windows.old
fi
