CMAKE_MINIMUM_REQUIRED(VERSION 2.4)

SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

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


SICONOS_PROJECT(FrontEnd 3 0 0)
