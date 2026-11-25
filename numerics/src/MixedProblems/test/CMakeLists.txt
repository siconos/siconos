# SPDX-FileCopyrightText: INRIA
#
# SPDX-License-Identifier: Apache-2.0
#
# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.

# CMakeLists.txt for MixedProblems tests

# Simple test with M = Identity, q = -1
add_executable(mixed_problem_simple_test mixed_problem_simple_test.c)
target_link_libraries(mixed_problem_simple_test PRIVATE numerics)
target_include_directories(mixed_problem_simple_test PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/..
  ${CMAKE_SOURCE_DIR}/numerics/src/MixedProblems
  ${CMAKE_SOURCE_DIR}/numerics/src/FrictionContact
  ${CMAKE_SOURCE_DIR}/numerics/src/RollingFrictionContact
  ${CMAKE_SOURCE_DIR}/numerics/src/Plasticity
  ${CMAKE_SOURCE_DIR}/numerics/src/core
)

# Add test
add_test(NAME mixed_problem_simple_test
  COMMAND mixed_problem_simple_test
)
