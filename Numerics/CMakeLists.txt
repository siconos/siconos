# Set minimum version for cmake
cmake_minimum_required(VERSION 2.8.7)

# MACOSX_RPATH required from cmake 3.0.
# See http://www.kitware.com/blog/home/post/510
cmake_policy(SET CMP0042 OLD)
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

ENABLE_LANGUAGE(C)
ENABLE_LANGUAGE(CXX)
ENABLE_LANGUAGE(Fortran)

INCLUDE(SiconosProject)

SET(PACKAGE_DESCRIPTION_SUMMARY 
  "Siconos/Numerics provides low-level algorithms (NSDS field).")

SET(PACKAGE_DESCRIPTION 
  "Siconos is a program dedicated to modeling, simulation and control of
Non Smooth Dynamical Systems (NSDS). The present package,
Siconos/Numerics provides low-level algorithms to compute basic
well-identified problems.")

IF(PIPOL_IMAGE MATCHES "etch")
  SET(DEBIAN_PACKAGE_DEPENDS 
    "g++, gfortran, atlas3-base, atlas3-base-dev, atlas3-headers")
ENDIF(PIPOL_IMAGE MATCHES "etch")

IF(PIPOL_IMAGE MATCHES "lenny|ubuntu")
  SET(DEBIAN_PACKAGE_DEPENDS
    "g++, gfortran, libatlas-base-dev, libatlas-headers")
ENDIF(PIPOL_IMAGE MATCHES "lenny|ubuntu")  

IF(FAT_BINARY)
  SET (CMAKE_OSX_ARCHITECTURES "i386;x86_64" )
ENDIF()

SICONOS_PROJECT(Numerics)

# ============= Summary =============
# To have a full print of variables, use -DPRINT_ENV=ON
message(STATUS "\n============================================ Summary ============================================")
message(STATUS "${PROJECT_NAME} version ${VERSION} is now ready for compilation and installation.")
message(STATUS "To proceed run 'make' and 'make install' and optionaly 'make test'.")
message(STATUS " C++ Compiler : ${CMAKE_CXX_COMPILER}")
message(STATUS " C Compiler : ${CMAKE_C_COMPILER}")
message(STATUS " Fortran Compiler : ${CMAKE_Fortran_COMPILER}")
message(STATUS " Compilation mode is : ${CMAKE_BUILD_TYPE}")
message(STATUS " Code Sources are in : ${CMAKE_SOURCE_DIR}")
message(STATUS " Blas is ${WITH_BLAS}. Libraries : ${BLAS_LIBRARIES}. Header : ${BLAS_INCLUDE_DIRS}/${BLAS_HEADER}.")
message(STATUS " Lapack is ${WITH_LAPACK}. Libraries : ${LAPACK_LIBRARIES}. Header : ${LAPACK_INCLUDE_DIRS}/${LAPACK_HEADER}.")
message(STATUS " The project will be installed in ${CMAKE_INSTALL_PREFIX}")
message(STATUS "To get more information about dependencies, config or else, ")
message(STATUS "check CMakeCache.txt file or re-run cmake with -DPRINT_ENV=ON.")
message(STATUS "=================================================================================================\n")
