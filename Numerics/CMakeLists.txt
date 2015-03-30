# Set minimum version for cmake
cmake_minimum_required(VERSION 2.8.7)

# MACOSX_RPATH required from cmake 3.0.
# See http://www.kitware.com/blog/home/post/510
if(POLICY CMP0042)
  cmake_policy(SET CMP0042 OLD)
endif()
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

enable_language(C)
enable_language(CXX)
enable_language(Fortran)

include(SiconosProject)

set(PACKAGE_DESCRIPTION_SUMMARY 
  "Siconos/Numerics provides low-level algorithms (NSDS field).")

set(PACKAGE_DESCRIPTION 
  "Siconos is a program dedicated to modeling, simulation and control of
Non Smooth Dynamical Systems (NSDS). The present package,
Siconos/Numerics provides low-level algorithms to compute basic
well-identified problems.")

if(FAT_BINARY)
  set(CMAKE_OSX_ARCHITECTURES "i386;x86_64" )
endif()

siconos_project(Numerics)
