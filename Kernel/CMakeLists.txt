CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# CMP0011 policy before first include
# others policies are set in SiconosProject
IF(POLICY CMP0011)
  CMAKE_POLICY(SET CMP0011 NEW)
ENDIF(POLICY CMP0011)
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

SICONOS_PROJECT(Kernel 3 4 0)
