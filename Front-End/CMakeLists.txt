CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

ENABLE_LANGUAGE(C)
ENABLE_LANGUAGE(CXX)
ENABLE_LANGUAGE(Fortran)

# CMP0011 policy before first include
# others policies are set in SiconosProject
IF(POLICY CMP0011)
  CMAKE_POLICY(SET CMP0011 NEW)
ENDIF(POLICY CMP0011)

INCLUDE(SiconosProject)

SET(PACKAGE_DESCRIPTION_SUMMARY 
  "Siconos/Front-End provides Python bindings to the Siconos library.")

SET(PACKAGE_DESCRIPTION 
  "Siconos is a program dedicated to modeling, simulation and control of
Non Smooth Dynamical Systems (NSDS). The present package,
Siconos/Front-End provides Python bindings.")

SET(DEBIAN_PACKAGE_DEPENDS 
  "g++, python, libboost-dev")


SICONOS_PROJECT(FrontEnd 3 0 0)
