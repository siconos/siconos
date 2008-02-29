SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

ENABLE_LANGUAGE(C)
ENABLE_LANGUAGE(CXX)
ENABLE_LANGUAGE(Fortran)
FIND_PACKAGE(ATLAS)

INCLUDE(SiconosProject)

SICONOS_PROJECT(Numerics 2 2 1 
  "The Siconos/Numerics Package is dedicated to low-level algorithms")

