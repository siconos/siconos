SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

INCLUDE(SiconosProject)

ENABLE_LANGUAGE(C)
ENABLE_LANGUAGE(CXX)
ENABLE_LANGUAGE(Fortran) # see lapack macros

SICONOS_PROJECT(Kernel 2 2 1 
  "Siconos/Kernel is dedicated to the modeling and the simulation of NSDS"

  "The present package, Siconos/Kernel, is dedicated to the modeling
and the simulation of Non Smooth Dynamical Systems (NSDS), with high
level description and numerical solving strategies.  It relies on
Siconos/Numerics which provides low-level algorithms to compute basic
well-identified problems.")

