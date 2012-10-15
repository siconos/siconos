CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)

SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
  
ENABLE_LANGUAGE(CXX)

INCLUDE(SiconosProject)

SET(PACKAGE_DESCRIPTION_SUMMARY 
  "Siconos/IO provides input/output services to Siconos/Kernel")

SET(PACKAGE_DESCRIPTION 
  "Siconos is a program dedicated to modeling, simulation and control of
Non Smooth Dynamical Systems (NSDS). The present package,
Siconos/IO provides some services for data input/output")

SICONOS_PROJECT(IO 3 5 0)
