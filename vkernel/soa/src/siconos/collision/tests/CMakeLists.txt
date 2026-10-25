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

endif()
