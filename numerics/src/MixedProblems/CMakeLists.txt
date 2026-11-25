# SPDX-FileCopyrightText: 2026 INRIA
#
# SPDX-License-Identifier: Apache-2.0
#
# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.

# Siconos numerics MixedProblems library

set(MixedProblems_src
  MixedProblem.c
  MixedProblem_NSGS.c
)

# Build object library
add_library(MixedProblems OBJECT ${MixedProblems_src})
target_include_directories(MixedProblems PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../tools
  ${CMAKE_CURRENT_SOURCE_DIR}/../tools/internal
  ${CMAKE_CURRENT_SOURCE_DIR}/../LinearAlgebra
  ${CMAKE_CURRENT_SOURCE_DIR}/../FrictionContact
  ${CMAKE_CURRENT_SOURCE_DIR}/../RollingFrictionContact
  ${CMAKE_CURRENT_SOURCE_DIR}/../Plasticity
  ${CMAKE_CURRENT_SOURCE_DIR}/../NonSmoothSolvers
)

# Add to the main numerics library
if(BUILD_SHARED_LIBS)
  set_property(TARGET MixedProblems PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()

# Installation of headers
install(FILES
  MixedProblem.h
  MixedProblems_options.h
  MixedProblem_NSGS.h
  DESTINATION include/siconos/numerics
)
