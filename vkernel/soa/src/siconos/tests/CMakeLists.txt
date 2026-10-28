if(WITH_vkernel_TESTING)
  if (WITH_GPU)
    set(COMPACTNSEARCH cuNSearch)
  else()
    set(COMPACTNSEARCH CompactNSearch)
  endif()

  # The ctest tests of this directory are registered in
  # vkernel/soa/vkernel_tests.cmake (begin_tests/new_test). This file only
  # keeps the targets that are not ctest tests.

  add_executable(registry_remove registry_remove.cpp)
  set_target_properties(registry_remove PROPERTIES COMPILE_FLAGS "${flags}")
  target_link_libraries(registry_remove Eigen3::Eigen numerics pthread)
  add_test(NAME vkernel_registry_remove COMMAND registry_remove
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

  add_executable(fem_t3 fem_t3.cpp)
  set_target_properties(fem_t3 PROPERTIES COMPILE_FLAGS "${flags}")

  target_link_libraries(fem_t3 Eigen3::Eigen numerics kernel mechanics pthread ${COMPACTNSEARCH})

  # not yet cppunit suites
  foreach(test nslaws_multi_group polydisperse_disks)
    add_executable(${test} ${test}.cpp)
    set_target_properties(${test} PROPERTIES COMPILE_FLAGS "${flags}")
    target_link_libraries(${test} Eigen3::Eigen numerics pthread ${COMPACTNSEARCH})
  endforeach()

  # fem_t3 is not registered: its mesh square_200.msh and reference
  # T3_square_200.ref are not in the repository.
  # bouncing_disk_sf and tnonos are not built today.
endif()
