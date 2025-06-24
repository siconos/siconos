if(WITH_vkernel_BENCHMARKS)
  add_executable(ColumnOfBeadsTS  ColumnOfBeadsTS.cpp)
  add_executable(column_of_beads_ts column_of_beads_ts.cpp)
  #add_executable(tnonos tnonos.cpp)


  set(flags "-Wall --std=c++20")
  set(fconcepts_diagnostics)
  check_cxx_compiler_flag("-fconcepts-diagnostics-depth=0" fconcepts_diagnostics)
  if(fconcepts_diagnostics)
    set(flags "${flags} -fconcepts-diagnostics-depth=0")
  endif()
  set(gdwarf4)
  check_cxx_compiler_flag("-gdwarf-4" gdwarf4) # clang / perf annotate
  if(gdwarf4)
    set(flags "${flags} -gdwarf-4")
  endif()


  set_target_properties(ColumnOfBeadsTS PROPERTIES COMPILE_FLAGS "${flags}")# -fconcepts-diagnostics-depth=10")
  set_target_properties(column_of_beads_ts PROPERTIES COMPILE_FLAGS "${flags}")# -fconcepts-diagnostics-depth=10")

  target_link_libraries(column_of_beads_ts PUBLIC Eigen3::Eigen Siconos::numerics pthread fmt)
  target_link_libraries(ColumnOfBeadsTS PUBLIC Siconos::numerics Siconos::kernel pthread fmt)
endif()
