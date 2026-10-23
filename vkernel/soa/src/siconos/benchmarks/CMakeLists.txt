option(WITH_vkernel_BENCHMARKS "Build the vkernel benchmarks" OFF)

if(WITH_vkernel_BENCHMARKS)
  # The two legacy benchmarks (ColumnOfBeadsTS, column_of_beads_ts) were
  # deleted: they had never been configured (this directory used to override
  # the c++ standard) and no longer compile against the current storage API
  # (undeclared symbols / missing siconos/utils/print.hpp). They are kept in
  # git history if anyone wants to revive them.
  add_executable(polytope_narrowphase_scale polytope_narrowphase_scale.cpp)
  # same source twice: static "shape" declaration vs per-kind shape sets
  add_executable(disks_shape_static disks_shape_dispatch.cpp)
  add_executable(disks_shape_per_kind disks_shape_dispatch.cpp)
  target_compile_definitions(disks_shape_per_kind PRIVATE PER_KIND)
  #add_executable(tnonos tnonos.cpp)

  if (WITH_GPU)
    set(COMPACTNSEARCH cuNSearch)
  else()
    set(COMPACTNSEARCH CompactNSearch)
  endif()


  # the parent scope (vkernel/soa/CMakeLists.txt) already sets flags with
  # --std=c++23 -- the vkernel sources need std::ranges::views::zip -- and
  # -fconcepts-diagnostics-depth=0. This directory used to override them
  # with a stale --std=c++20, which is why it had never been configured.
  set(gdwarf4)
  check_cxx_compiler_flag("-gdwarf-4" gdwarf4) # clang / perf annotate
  if(gdwarf4)
    set(flags "${flags} -gdwarf-4")
  endif()


  set_target_properties(polytope_narrowphase_scale PROPERTIES COMPILE_FLAGS "${flags}")# -fconcepts-diagnostics-depth=10")

  # the Siconos:: aliases below never existed as targets: this directory had
  # never been configured, so the stale names went unnoticed
  target_link_libraries(polytope_narrowphase_scale PUBLIC Eigen3::Eigen numerics pthread ${COMPACTNSEARCH})
  foreach(bench disks_shape_static disks_shape_per_kind)
    set_target_properties(${bench} PROPERTIES COMPILE_FLAGS "${flags}")
    target_link_libraries(${bench} PUBLIC Eigen3::Eigen numerics pthread ${COMPACTNSEARCH})
  endforeach()
endif()
