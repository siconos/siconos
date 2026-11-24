# SPDX-FileCopyrightText: 2026 INRIA
#
# SPDX-License-Identifier: Apache-2.0
#
# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.

if(WITH_vkernel_TESTING)
  if (WITH_GPU)
    set(COMPACTNSEARCH cuNSearch)
  else()
    set(COMPACTNSEARCH CompactNSearch)
  endif()

  # The ctest tests of this directory are registered in
  # vkernel/vkernel_tests.cmake (begin_tests/new_test). This file only
  # keeps the targets that are not ctest tests.

  # the comparisons of Moreau-Jean with Anitescu-Tasora and DEM: CSV tables,
  # the target vkernel_study_time_integrators (vkernel/docs/notes) runs
  # it and makes the tables and the plots
  add_executable(nscd_comparisons nscd_comparisons.cpp)
  set_target_properties(nscd_comparisons PROPERTIES COMPILE_FLAGS "${flags}")
  target_link_libraries(nscd_comparisons vkernel numerics Eigen3::Eigen pthread ${COMPACTNSEARCH})

  # bouncing_disk_sf and tnonos are not built.
endif()
