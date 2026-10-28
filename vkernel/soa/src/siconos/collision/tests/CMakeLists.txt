if(WITH_vkernel_TESTING)
  if (WITH_GPU)
    set(COMPACTNSEARCH cuNSearch)
  else()
    set(COMPACTNSEARCH CompactNSearch)
  endif()

  # The ctest tests of this directory are registered in
  # vkernel/soa/vkernel_tests.cmake (begin_tests/new_test). This file only
  # keeps the targets that are not ctest tests: the relation/index sets
  # oracles, which build the same test twice and compare their output.

  # oracle of the relation and index sets (simul/relation_sets.hpp,
  # simul/index_sets.hpp): the same test without the sets (_sets0) and with
  # them, assertions enabled (_sets1), must print the same output. Covers
  # the space filter (bouncing_disks), the cxd rebind and the removals
  # (churn_polytopes), polytopes, disks and FEM together (mixed_bodies) and
  # the rt ct interactions (polytope_on_fem) and the sort of the
  # interactions (sort_interactions).
  foreach(test bouncing_disks churn_polytopes mixed_bodies polytope_on_fem
      sort_interactions)
    foreach(sets 0 1)
      add_executable(${test}_sets${sets} ${test}.cpp ${SIMPLE_TEST_MAIN})
      set_target_properties(${test}_sets${sets} PROPERTIES COMPILE_FLAGS "${flags}")
      target_compile_definitions(${test}_sets${sets} PRIVATE
        VKERNEL_RELATION_SETS=${sets} VKERNEL_INDEX_SETS=${sets})
      target_link_libraries(${test}_sets${sets} Eigen3::Eigen numerics pthread
        ${COMPACTNSEARCH} CPPUNIT::CPPUNIT)
    endforeach()
    target_compile_options(${test}_sets1 PRIVATE -UNDEBUG)
    # a test that checks a reference reads it from its working directory
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test}_sets_oracle)
    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${test}.ref)
      configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${test}.ref
        ${CMAKE_CURRENT_BINARY_DIR}/${test}_sets_oracle/${test}.ref COPYONLY)
    endif()
    add_test(NAME vkernel_${test}_sets_oracle
      COMMAND ${CMAKE_COMMAND}
        -DA=$<TARGET_FILE:${test}_sets0> -DB=$<TARGET_FILE:${test}_sets1>
        -DWD=${CMAKE_CURRENT_BINARY_DIR}/${test}_sets_oracle
        -P ${CMAKE_CURRENT_SOURCE_DIR}/compare_outputs.cmake)
    # a multithreaded neighbor search creates the interactions in any order
    set_tests_properties(vkernel_${test}_sets_oracle PROPERTIES
      ENVIRONMENT OMP_NUM_THREADS=1)
  endforeach()

endif()
