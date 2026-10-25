if(WITH_vkernel_TESTING)
  if (WITH_GPU)
    set(COMPACTNSEARCH cuNSearch)
  else()
    set(COMPACTNSEARCH CompactNSearch)
  endif()
  add_executable(bouncing_ball bouncing_ball.cpp)
  add_executable(bouncing_ball_registry bouncing_ball_registry.cpp)
  add_executable(bouncing_ball_split bouncing_ball_split.cpp bouncing_ball_split_sim.cpp)
  add_executable(registry_remove registry_remove.cpp)
  add_executable(unbounded_bouncing_ball unbounded_bouncing_ball.cpp)
  add_executable(bouncing_balls bouncing_balls.cpp)
  add_executable(bouncing_disk bouncing_disk.cpp)
  add_executable(bouncing_disk_line bouncing_disk_line.cpp)
  #add_executable(bouncing_disk_sf bouncing_disk_sf.cpp)
  add_executable(nslaws_multi_group nslaws_multi_group.cpp)
  add_executable(polydisperse_disks polydisperse_disks.cpp)
  #add_executable(tnonos tnonos.cpp)

  add_executable(fem_t3 fem_t3.cpp)

  set_target_properties(bouncing_ball_registry PROPERTIES COMPILE_FLAGS "${flags}")
  set_target_properties(bouncing_ball_split PROPERTIES COMPILE_FLAGS "${flags}")
  set_target_properties(registry_remove PROPERTIES COMPILE_FLAGS "${flags}")
  set_target_properties(bouncing_ball PROPERTIES COMPILE_FLAGS "${flags}")# -fconcepts-diagnostics-depth=10")
  set_target_properties(unbounded_bouncing_ball PROPERTIES COMPILE_FLAGS "${flags}")# -fconcepts-diagnostics-depth=10")
  set_target_properties(bouncing_balls PROPERTIES COMPILE_FLAGS "${flags}")# -fconcepts-diagnostics-depth=10")
  set_target_properties(bouncing_disk PROPERTIES COMPILE_FLAGS "${flags}")# -fconcepts-diagnostics-depth=10")
  #set_target_properties(bouncing_disk_sf PROPERTIES COMPILE_FLAGS "${flags}")# -fconcepts-diagnostics-depth=10")
  set_target_properties(bouncing_disk_line PROPERTIES COMPILE_FLAGS "${flags}")# -fconcepts-diagnostics-depth=10")
  #set_target_properties(bouncing_disk_sf PROPERTIES COMPILE_FLAGS "${flags}")# -fconcepts-diagnostics-depth=10")
  set_target_properties(nslaws_multi_group PROPERTIES COMPILE_FLAGS "${flags}")
  set_target_properties(polydisperse_disks PROPERTIES COMPILE_FLAGS "${flags}")
  set_target_properties(fem_t3 PROPERTIES COMPILE_FLAGS "${flags}")

  target_link_libraries(bouncing_ball  Eigen3::Eigen numerics pthread)
  target_link_libraries(bouncing_ball_registry  Eigen3::Eigen numerics pthread)
  target_link_libraries(bouncing_ball_split  Eigen3::Eigen numerics pthread)
  target_link_libraries(registry_remove  Eigen3::Eigen numerics pthread)
  target_link_libraries(unbounded_bouncing_ball  Eigen3::Eigen numerics pthread)
  target_link_libraries(bouncing_balls Eigen3::Eigen numerics pthread)
  target_link_libraries(bouncing_disk Eigen3::Eigen numerics pthread ${COMPACTNSEARCH})
  #target_link_libraries(bouncing_disk_sf Eigen3::Eigen numerics pthread fmt CompactNSearch)
  target_link_libraries(bouncing_disk_line Eigen3::Eigen numerics pthread ${COMPACTNSEARCH})
  target_link_libraries(nslaws_multi_group Eigen3::Eigen numerics pthread ${COMPACTNSEARCH})
  target_link_libraries(polydisperse_disks Eigen3::Eigen numerics pthread ${COMPACTNSEARCH})

  target_link_libraries(fem_t3 Eigen3::Eigen numerics kernel mechanics pthread ${COMPACTNSEARCH})

  # ctest --test-dir <build>/vkernel/soa/src/siconos/tests
  # fem_t3 is not registered: its mesh square_200.msh and reference
  # T3_square_200.ref are not in the repository.
  enable_testing()
  foreach(test bouncing_ball bouncing_ball_registry bouncing_ball_split
      registry_remove unbounded_bouncing_ball bouncing_balls bouncing_disk
      bouncing_disk_line)
    add_test(NAME vkernel_${test} COMMAND ${test}
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  endforeach()

endif()
