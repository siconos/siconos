# Set minimum version for cmake
cmake_minimum_required(VERSION 2.8.7)

# MACOSX_RPATH required from cmake 3.0.
# See http://www.kitware.com/blog/home/post/510
if(POLICY CMP0042)
  cmake_policy(SET CMP0042 OLD)
endif()
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

include(SiconosTools)
get_uninitialized_vars(CMAKE_ARGS_FROM_CACHE)

option(WITH_OCC "compilation with OpenCascade Bindings. Default = OFF" OFF)
option(WITH_MECHANISMS "compilation with OpenCascade Bindings. Default = OFF" OFF)
option(WITH_BULLET "compilation with Bullet Bindings. Default = OFF" OFF)

include(SiconosProject)

enable_language(C)
enable_language(CXX)
enable_language(Fortran) # see lapack macros

set(PACKAGE_DESCRIPTION_SUMMARY
  "Siconos/Mechanics is dedicated to the modeling and the simulation of nonsmooth mechanical systems")

set(PACKAGE_DESCRIPTION
  "The present package, Siconos/Mechanics, is dedicated to the modeling and
the simulation of Non Smooth Mechanical Systems")

siconos_project(Mechanics)
