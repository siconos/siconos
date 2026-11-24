# SPDX-FileCopyrightText: 2026 INRIA
#
# SPDX-License-Identifier: Apache-2.0
#
# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.

option(WITH_vkernel_BENCHMARKS "Build the vkernel benchmarks" OFF)

if(WITH_vkernel_BENCHMARKS)
  add_executable(polytope_narrowphase_scale polytope_narrowphase_scale.cpp)
  # same source twice: static "shape" declaration vs per-kind shape sets
  add_executable(disks_shape_static disks_shape_dispatch.cpp)
  add_executable(disks_shape_per_kind disks_shape_dispatch.cpp)
  target_compile_definitions(disks_shape_per_kind PRIVATE PER_KIND)
  # the relation and index sets on the configuration of the python module
  add_executable(mixed_relations mixed_relations.cpp)
  #add_executable(tnonos tnonos.cpp)

  if (WITH_GPU)
    set(COMPACTNSEARCH cuNSearch)
  else()
    set(COMPACTNSEARCH CompactNSearch)
  endif()


  # the parent scope (vkernel/CMakeLists.txt) already sets flags with
  # --std=c++23 -- the vkernel sources need std::ranges::views::zip -- and
  # -fconcepts-diagnostics-depth=0: they must not be overridden here.
  set(gdwarf4)
  check_cxx_compiler_flag("-gdwarf-4" gdwarf4) # clang / perf annotate
  if(gdwarf4)
    set(flags "${flags} -gdwarf-4")
  endif()


  set_target_properties(polytope_narrowphase_scale PROPERTIES COMPILE_FLAGS "${flags}")# -fconcepts-diagnostics-depth=10")

  target_link_libraries(polytope_narrowphase_scale PUBLIC Eigen3::Eigen numerics pthread ${COMPACTNSEARCH})
  foreach(bench disks_shape_static disks_shape_per_kind mixed_relations)
    set_target_properties(${bench} PROPERTIES COMPILE_FLAGS "${flags}")
    target_link_libraries(${bench} PUBLIC Eigen3::Eigen numerics pthread ${COMPACTNSEARCH})
  endforeach()
endif()
