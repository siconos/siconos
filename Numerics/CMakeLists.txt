CMAKE_MINIMUM_REQUIRED(VERSION 2.4)

SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
  
if(COMMAND cmake_policy)
  cmake_policy(VERSION 2.6.0)
  cmake_policy(SET CMP0011 NEW)
endif(COMMAND cmake_policy)

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


SICONOS_PROJECT(Numerics 3 0 0)


  

