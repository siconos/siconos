#!/bin/sh

MXE_PREFIX="/scratch/Olivier/mingw32/"
INSTALL_PREFIX="/scratch/Olivier/siconos-windows/install"
SICONOS_SOURCES="$HOME/siconos"
BUILD_DIR="/scratch/Olivier/siconos-windows/build"

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
	mkdir -p "${BUILD_DIR}/Numerics"
	cd "${BUILD_DIR}/Numerics"

	WINEARCH=win64

	wineserver -p
	wine mspdbsrv.exe -start -shutdowntime -1 >/dev/null 2>&1 &

	CFLAGS='-U__STRICT_ANSI__' cmake \
		-DCMAKE_TOOLCHAIN_FILE="${MXE_PREFIX}/mxe/usr/i686-pc-mingw32/share/cmake/mxe-conf.cmake" \
		-DCROSSCOMPILING_LINUX_TO_WINDOWS=1 \
		-DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}" \
		"${SICONOS_SOURCES}/Numerics"
	make_component
	make -i install

	mkdir -p "${BUILD_DIR}/Kernel"
	cd "${BUILD_DIR}/Kernel"
	cmake \
		-DCMAKE_TOOLCHAIN_FILE="${MXE_PREFIX}/mxe/usr/i686-pc-mingw32/share/cmake/mxe-conf.cmake" \
		-DCROSSCOMPILING_LINUX_TO_WINDOWS=1 \
		-DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}" \
		-DSiconosNumerics_INCLUDE_DIRS="${INSTALL_PREFIX}/Siconos/Numerics/" \
		-DSiconosNumerics_FOUND="${INSTALL_PREFIX}/lib/libSiconosNumerics.dll" \
		"${SICONOS_SOURCES}/Kernel"
	make_component
	make -i install

	mkdir -p "${BUILD_DIR}/Front-End"
	cd "${BUILD_DIR}/Front-End"
	cmake \
		-DCMAKE_TOOLCHAIN_FILE="${MXE_PREFIX}/mxe/usr/i686-pc-mingw32/share/cmake/mxe-conf.cmake" \
		-DCROSSCOMPILING_LINUX_TO_WINDOWS=1 \
		-DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}" \
		-DSiconosNumerics_INCLUDE_DIRS="${INSTALL_PREFIX}/include/Siconos/Numerics/" \
		-DSiconosNumerics_FOUND="${INSTALL_PREFIX}/lib/libSiconosNumerics.dll" \
		-DSiconosKernel_INCLUDE_DIRS="${INSTALL_PREFIX}/include/Siconos/Kernel/" \
		-DSiconosKernel_FOUND="${INSTALL_PREFIX}/lib/libSiconosKernel.dll" \
		-DSWIG_DIR="/usr/share/swig/2.0.7/" \
		-DPYTHON_LIBRARIES="${MXE_PREFIX}/python/python27.dll" \
		-DPYTHON_INCLUDE_DIR="${MXE_PREFIX}/python/" \
		-DPYTHON_NUMPY_INCLUDE_DIR="${MXE_PREFIX}/numpy/PLATLIB/numpy/core/include" \
		"${SICONOS_SOURCES}/Front-End"
	make_component
	make -i install

}

check_git() {
cd ${SICONOS_SOURCES}
git pull
CURRENT_SHA=`git rev-parse HEAD`
if [ -f /tmp/.sha-siconos-windows.old ]; then
	OLD_SHA=`cat /tmp/.sha-siconos-windows.old`
else
	OLD_SHA=""
fi

if [ "${CURRENT_SHA}" != "${OLD_SHA}" ]; then
	build_siconos
	echo ${CURRENT_SHA} > /tmp/.sha-siconos-windows.old
fi
}

build_siconos
