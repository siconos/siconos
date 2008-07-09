CMAKE_MINIMUM_REQUIRED(VERSION 2.4)

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

SET(DEBIAN_PACKAGE_DEPENDS 
  "g++, gfortran, atlas3-base, atlas3-base-dev, atlas3-headers, libxml2-dev, doxygen, siconos-numerics")

SICONOS_PROJECT(Kernel 3 0 0)



