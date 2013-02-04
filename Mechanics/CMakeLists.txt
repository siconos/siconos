CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)

SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

INCLUDE(SiconosProject)

ENABLE_LANGUAGE(C)
ENABLE_LANGUAGE(CXX)
ENABLE_LANGUAGE(Fortran) # see lapack macros
 
SET(PACKAGE_DESCRIPTION_SUMMARY 
  "Siconos/Mechanics is dedicated to the modeling and the simulation of nonsmooth mechanical systems")

SET(PACKAGE_DESCRIPTION 
  "The present package, Siconos/Mechanics, is dedicated to the modeling and
the simulation of Non Smooth Mechanical Systems")

#IF(APPLE)
#  SET (CMAKE_OSX_ARCHITECTURES "i386;x86_64" )
#ENDIF(APPLE)

SICONOS_PROJECT(Mechanics)

