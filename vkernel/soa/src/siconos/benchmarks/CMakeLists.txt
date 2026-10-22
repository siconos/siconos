option(WITH_vkernel_BENCHMARKS "Build the vkernel benchmarks" OFF)

if(WITH_vkernel_BENCHMARKS)
  # ColumnOfBeadsTS (legacy kernel) and column_of_beads_ts (vkernel) do not
  # compile: the first uses undeclared symbols (H, ...), the second includes
  # a nonexistent siconos/utils/print.hpp. They had never been configured
  # (this directory used to override the c++ standard), so they are kept out
  # of the build until someone revives them.
  #add_executable(ColumnOfBeadsTS  ColumnOfBeadsTS.cpp)
  #add_executable(column_of_beads_ts column_of_beads_ts.cpp)
  add_executable(polytope_narrowphase_scale polytope_narrowphase_scale.cpp)
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
endif()
