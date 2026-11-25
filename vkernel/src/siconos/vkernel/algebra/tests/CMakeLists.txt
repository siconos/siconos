# SPDX-FileCopyrightText: INRIA
#
# SPDX-License-Identifier: Apache-2.0
#
# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.

add_executable(ct_eigen ct_eigen.cpp)
add_executable(ct_numerics ct_numerics.cpp)

set_target_properties(ct_eigen PROPERTIES COMPILE_FLAGS "${flags}")# -fconcepts-diagnostics-depth=10")
set_target_properties(ct_numerics PROPERTIES COMPILE_FLAGS "${flags}")# -fconcepts-diagnostics-depth=10")
target_link_libraries(ct_eigen Eigen3::Eigen)
target_link_libraries(ct_numerics Eigen3::Eigen numerics)
