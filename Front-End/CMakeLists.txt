# Set minimum version for cmake
cmake_minimum_required(VERSION 2.8.7)

# MACOSX_RPATH required from cmake 3.0.
# See http://www.kitware.com/blog/home/post/510
cmake_policy(SET CMP0042 OLD)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

ENABLE_LANGUAGE(C)
ENABLE_LANGUAGE(CXX)

INCLUDE(SiconosProject)

SET(PACKAGE_DESCRIPTION_SUMMARY 
  "Siconos/Front-End provides Python bindings to the Siconos library.")

SET(PACKAGE_DESCRIPTION 
  "Siconos is a program dedicated to modeling, simulation and control of
Non Smooth Dynamical Systems (NSDS). The present package,
Siconos/Front-End provides Python bindings.")

SET(DEBIAN_PACKAGE_DEPENDS 
  "g++, python, libboost-dev")

SICONOS_PROJECT(FrontEnd)

