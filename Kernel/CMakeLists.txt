CMAKE_MINIMUM_REQUIRED(VERSION 2.6.3)

SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

INCLUDE(SiconosProject)

ENABLE_LANGUAGE(C)
ENABLE_LANGUAGE(CXX)
ENABLE_LANGUAGE(Fortran) # see lapack macros
 
SET(PACKAGE_DESCRIPTION_SUMMARY 
  "Siconos/Kernel is dedicated to the modeling and the simulation of NSDS.")

SET(PACKAGE_DESCRIPTION 
  "The present package, Siconos/Kernel, is dedicated to the modeling and
the simulation of Non Smooth Dynamical Systems (NSDS), with high level
description and numerical solving strategies.  It relies on
Siconos/Numerics which provides low-level algorithms to compute basic
well-identified problems.")


IF(PIPOL_IMAGE MATCHES "etch")
  SET(DEBIAN_PACKAGE_DEPENDS 
    "g++, gfortran, atlas3-base, atlas3-base-dev, atlas3-headers, libxml2-dev, libboost-dev, libboost-graph-dev, libgfortran1-dev, libgmp3-dev, siconos-numerics")
ENDIF(PIPOL_IMAGE MATCHES "etch")

IF(PIPOL_IMAGE MATCHES "lenny\\|ubuntu")
  SET(DEBIAN_PACKAGE_DEPENDS
    "g++, gfortran, libatlas-base-dev, libatlas-headers, libboost-dev, libboost-graph-dev, libgmp3-dev, libxml2-dev, siconos-numerics")
ENDIF(PIPOL_IMAGE MATCHES "lenny\\|ubuntu")  

#IF(APPLE)
#  SET (CMAKE_OSX_ARCHITECTURES "i386;x86_64" )
#ENDIF(APPLE)

SICONOS_PROJECT(Kernel 3 4 0)

# To have a full print of variables, use -DPRINT_ENV=ON
message(STATUS "\n============================================ Summary ============================================")
message(STATUS "${PROJECT_NAME} version ${VERSION} is now ready for compulation and installation.")
message(STATUS "To proceed run 'make' and 'make install' and optionaly 'make test'.")
message(STATUS " C++ Compiler : ${CMAKE_CXX_COMPILER}")
message(STATUS " C Compiler : ${CMAKE_C_COMPILER}")
message(STATUS " Fortran Compiler : ${CMAKE_Fortran_COMPILER}")
message(STATUS " Compilation mode is : ${CMAKE_BUILD_TYPE}")
message(STATUS " Code Sources are in : ${CMAKE_SOURCE_DIR}")
message(STATUS " The project will be installed in ${CMAKE_INSTALL_PREFIX}")
message(STATUS " The Numerics library used is : ${SiconosNumerics_FOUND}")
message(STATUS " To get more information about dependencies, config or else, ")
message(STATUS "check CMakeCache.txt file or re-run cmake with -DPRINT_ENV=ON.")
message(STATUS "=================================================================================================\n")
