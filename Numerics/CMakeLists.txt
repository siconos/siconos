SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

ENABLE_LANGUAGE(C)
ENABLE_LANGUAGE(CXX)
ENABLE_LANGUAGE(Fortran)

INCLUDE(SiconosProject)

SICONOS_PROJECT(Numerics 2 2 1 
  "Siconos/Numerics provides low-level algorithms (NSDS field)."

  "Siconos is a program dedicated to modeling, simulation and control
of Non Smooth Dynamical Systems (NSDS). The present package,
Siconos/Numerics provides low-level algorithms to compute basic
well-identified problems.")

