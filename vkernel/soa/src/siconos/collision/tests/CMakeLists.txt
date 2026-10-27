if(WITH_vkernel_TESTING)
  if (WITH_GPU)
    set(COMPACTNSEARCH cuNSearch)
  else()
    set(COMPACTNSEARCH CompactNSearch)
  endif()

  add_executable(cxd_batch_smoke cxd_batch_smoke.cpp)
  add_executable(polytope_cxd polytope_cxd.cpp)
  add_executable(bouncing_polytopes bouncing_polytopes.cpp)
  add_executable(churn_polytopes churn_polytopes.cpp)
  add_executable(polytope_fsegment polytope_fsegment.cpp)
  add_executable(polytopemesh polytopemesh.cpp)
  add_executable(polytopefdisk polytopefdisk.cpp)
  add_executable(polytopedisk polytopedisk.cpp)
  add_executable(rigid_shapes rigid_shapes.cpp)
  add_executable(mixed_bodies mixed_bodies.cpp)
  add_executable(polytope_on_fem polytope_on_fem.cpp)
  add_executable(bouncing_disks bouncing_disks.cpp)
  add_executable(fixed_disks fixed_disks.cpp)
  add_executable(io_contact_points io_contact_points.cpp)

  foreach(t cxd_batch_smoke polytope_cxd bouncing_polytopes churn_polytopes
            polytope_fsegment polytopemesh polytopefdisk polytopedisk
            rigid_shapes mixed_bodies polytope_on_fem bouncing_disks
            fixed_disks io_contact_points)
    set_target_properties(${t} PROPERTIES COMPILE_FLAGS "${flags}")
  endforeach()

  target_link_libraries(cxd_batch_smoke Eigen3::Eigen numerics pthread)
  target_link_libraries(polytope_cxd Eigen3::Eigen numerics pthread)
  target_link_libraries(bouncing_polytopes Eigen3::Eigen numerics pthread ${COMPACTNSEARCH})
  target_link_libraries(churn_polytopes Eigen3::Eigen numerics pthread ${COMPACTNSEARCH})
  target_link_libraries(polytope_fsegment Eigen3::Eigen numerics pthread ${COMPACTNSEARCH})
  target_link_libraries(polytopemesh Eigen3::Eigen numerics pthread ${COMPACTNSEARCH})
  target_link_libraries(polytopefdisk Eigen3::Eigen numerics pthread ${COMPACTNSEARCH})
  target_link_libraries(polytopedisk Eigen3::Eigen numerics pthread ${COMPACTNSEARCH})
  target_link_libraries(rigid_shapes Eigen3::Eigen numerics pthread CPPUNIT::CPPUNIT)
  target_link_libraries(mixed_bodies Eigen3::Eigen numerics pthread ${COMPACTNSEARCH})
  target_link_libraries(polytope_on_fem Eigen3::Eigen numerics pthread ${COMPACTNSEARCH})
  target_link_libraries(bouncing_disks Eigen3::Eigen numerics pthread ${COMPACTNSEARCH})
  target_link_libraries(fixed_disks Eigen3::Eigen numerics pthread ${COMPACTNSEARCH})
  target_link_libraries(io_contact_points Eigen3::Eigen numerics pthread ${COMPACTNSEARCH})

  enable_testing()
  foreach(test cxd_batch_smoke polytope_cxd bouncing_polytopes churn_polytopes
             polytope_fsegment polytopemesh polytopefdisk polytopedisk
             rigid_shapes mixed_bodies polytope_on_fem bouncing_disks
             fixed_disks io_contact_points)
    add_test(NAME vkernel_${test} COMMAND ${test}
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  endforeach()

  # topology::sort_interactions (storage::permute) on both storage backends:
  # the registry, and the compile time database (_static_db)
  foreach(backend registry static_db)
    set(t sort_interactions_${backend})
    add_executable(${t} sort_interactions.cpp)
    set_target_properties(${t} PROPERTIES COMPILE_FLAGS "${flags}")
    target_compile_options(${t} PRIVATE -UNDEBUG)
    if(backend STREQUAL static_db)
      target_compile_definitions(${t} PRIVATE SORT_STATIC_DB)
    endif()
    target_link_libraries(${t} Eigen3::Eigen numerics pthread ${COMPACTNSEARCH})
    add_test(NAME vkernel_${t} COMMAND ${t}
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    # the neighbor search order, hence the interactions, must be reproducible
    set_tests_properties(vkernel_${t} PROPERTIES ENVIRONMENT OMP_NUM_THREADS=1)
  endforeach()

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
      add_executable(${test}_sets${sets} ${test}.cpp)
      set_target_properties(${test}_sets${sets} PROPERTIES COMPILE_FLAGS "${flags}")
      target_compile_definitions(${test}_sets${sets} PRIVATE
        VKERNEL_RELATION_SETS=${sets} VKERNEL_INDEX_SETS=${sets})
      target_link_libraries(${test}_sets${sets} Eigen3::Eigen numerics pthread
        ${COMPACTNSEARCH})
    endforeach()
    target_compile_options(${test}_sets1 PRIVATE -UNDEBUG)
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
