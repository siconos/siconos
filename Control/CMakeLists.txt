CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)

SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

INCLUDE(SiconosProject)

ENABLE_LANGUAGE(CXX)

SET(PACKAGE_DESCRIPTION_SUMMARY 
  "Siconos/Control is dedicated to the control of dynamical systems")

SET(PACKAGE_DESCRIPTION 
  "The present package, Siconos/Control, is dedicated to the control
of dynamical systems. It is build on top of Siconos/Kernel, where all the
modeling is done.")

SICONOS_PROJECT(Control)

