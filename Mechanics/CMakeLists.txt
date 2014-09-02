# Set minimum version for cmake
cmake_minimum_required(VERSION 2.8.7)

# MACOSX_RPATH required from cmake 3.0.
# See http://www.kitware.com/blog/home/post/510
cmake_policy(SET CMP0042 OLD)
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

INCLUDE(SiconosProject)

ENABLE_LANGUAGE(C)
ENABLE_LANGUAGE(CXX)
ENABLE_LANGUAGE(Fortran) # see lapack macros

SET(PACKAGE_DESCRIPTION_SUMMARY
  "Siconos/Mechanics is dedicated to the modeling and the simulation of nonsmooth mechanical systems")

SET(PACKAGE_DESCRIPTION
  "The present package, Siconos/Mechanics, is dedicated to the modeling and
the simulation of Non Smooth Mechanical Systems")

SICONOS_PROJECT(Mechanics)

