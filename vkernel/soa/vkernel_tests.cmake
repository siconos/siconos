include(tools4tests)

if(WITH_TESTING)
  add_custom_target(vkernel-tests echo "Start vkernel tests")

  begin_tests(src/siconos/simul/tests DEPS
              "numerics;kernel;mechanics;CPPUNIT::CPPUNIT")
  new_test(SOURCES RtOsiTest.cpp ${SIMPLE_TEST_MAIN})
  # end_tests()

  begin_tests(src/siconos/storage/tests DEPS "vkernel;CPPUNIT::CPPUNIT")
  new_test(SOURCES StorageTest.cpp ${SIMPLE_TEST_MAIN})
  # the schedules (storage/schedule.hpp) and the Moreau-Jean time step as text
  new_test(NAME vkernel_schedule SOURCES ScheduleTest.cpp ${SIMPLE_TEST_MAIN} DEPS
           "numerics;kernel;mechanics")
  new_test(SOURCES t1.cpp ${SIMPLE_TEST_MAIN} DEPS "numerics;kernel;mechanics")
  # end_tests()

  if(WITH_GPU)
    set(COMPACTNSEARCH cuNSearch)
  else()
    set(COMPACTNSEARCH CompactNSearch)
  endif()

  # Tests that produce a trajectory also compare it to a committed reference
  # file (tests/ref_check.hpp) at a relative Frobenius error of 1e-8.
  # bouncing_ball's reference is the tutorial's BouncingBallTS.ref (external
  # truth); the others lock the current behaviour against regression.
  begin_tests(src/siconos/tests DEPS
              "vkernel;numerics;CPPUNIT::CPPUNIT;Eigen3::Eigen;pthread")
  new_test(
    NAME
    vkernel_bouncing_ball
    SOURCES
    bouncing_ball.cpp
    ${SIMPLE_TEST_MAIN}
    DATA
    bouncing_ball.ref)
  new_test(
    NAME
    vkernel_bouncing_ball_registry
    SOURCES
    bouncing_ball_registry.cpp
    ${SIMPLE_TEST_MAIN}
    DATA
    bouncing_ball_registry.ref)
  new_test(
    NAME
    vkernel_bouncing_ball_split
    SOURCES
    bouncing_ball_split.cpp
    bouncing_ball_split_sim.cpp
    ${SIMPLE_TEST_MAIN}
    DATA
    bouncing_ball_split.ref)
  new_test(
    NAME
    vkernel_unbounded_bouncing_ball
    SOURCES
    unbounded_bouncing_ball.cpp
    ${SIMPLE_TEST_MAIN}
    DATA
    unbounded_bouncing_ball.ref)
  new_test(
    NAME
    vkernel_bouncing_balls
    SOURCES
    bouncing_balls.cpp
    ${SIMPLE_TEST_MAIN}
    DATA
    bouncing_balls.ref)
  new_test(
    NAME
    vkernel_bouncing_disk
    SOURCES
    bouncing_disk.cpp
    ${SIMPLE_TEST_MAIN}
    DATA
    bouncing_disk.ref
    DEPS
    "${COMPACTNSEARCH}")
  new_test(
    NAME
    vkernel_bouncing_disk_line
    SOURCES
    bouncing_disk_line.cpp
    ${SIMPLE_TEST_MAIN}
    DATA
    bouncing_disk_line.ref
    DEPS
    "${COMPACTNSEARCH}")
  new_test(NAME vkernel_registry_remove SOURCES registry_remove.cpp
           ${SIMPLE_TEST_MAIN})
  new_test(NAME vkernel_nslaws_multi_group SOURCES nslaws_multi_group.cpp
           ${SIMPLE_TEST_MAIN} DEPS "${COMPACTNSEARCH}")
  new_test(NAME vkernel_polydisperse_disks SOURCES polydisperse_disks.cpp
           ${SIMPLE_TEST_MAIN} DEPS "${COMPACTNSEARCH}")

  begin_tests(
    src/siconos/collision/tests DEPS
    "vkernel;numerics;CPPUNIT::CPPUNIT;Eigen3::Eigen;pthread;${COMPACTNSEARCH}")
  new_test(
    NAME
    vkernel_bouncing_disks
    SOURCES
    bouncing_disks.cpp
    ${SIMPLE_TEST_MAIN}
    DATA
    bouncing_disks.ref)
  new_test(
    NAME
    vkernel_fixed_disks
    SOURCES
    fixed_disks.cpp
    ${SIMPLE_TEST_MAIN}
    DATA
    fixed_disks.ref)
  new_test(NAME vkernel_cxd_batch_smoke SOURCES cxd_batch_smoke.cpp
           ${SIMPLE_TEST_MAIN})
  new_test(NAME vkernel_polytope_cxd SOURCES polytope_cxd.cpp
           ${SIMPLE_TEST_MAIN})
  new_test(NAME vkernel_bouncing_polytopes SOURCES bouncing_polytopes.cpp
           ${SIMPLE_TEST_MAIN})
  new_test(NAME vkernel_churn_polytopes SOURCES churn_polytopes.cpp
           ${SIMPLE_TEST_MAIN})
  new_test(NAME vkernel_polytope_fsegment SOURCES polytope_fsegment.cpp
           ${SIMPLE_TEST_MAIN})
  new_test(NAME vkernel_polytopemesh SOURCES polytopemesh.cpp
           ${SIMPLE_TEST_MAIN})
  new_test(NAME vkernel_polytopefdisk SOURCES polytopefdisk.cpp
           ${SIMPLE_TEST_MAIN})
  new_test(NAME vkernel_polytopedisk SOURCES polytopedisk.cpp
           ${SIMPLE_TEST_MAIN})
  new_test(NAME vkernel_rigid_shapes SOURCES rigid_shapes.cpp
           ${SIMPLE_TEST_MAIN})
  new_test(NAME vkernel_mixed_bodies SOURCES mixed_bodies.cpp
           ${SIMPLE_TEST_MAIN})
  new_test(NAME vkernel_polytope_on_fem SOURCES polytope_on_fem.cpp
           ${SIMPLE_TEST_MAIN})
  new_test(NAME vkernel_io_contact_points SOURCES io_contact_points.cpp
           ${SIMPLE_TEST_MAIN})

  # the same deposit run twice with a multithreaded neighbor search must give
  # bitwise identical states
  new_test(NAME vkernel_reproducible_disks SOURCES reproducible_disks.cpp
           ${SIMPLE_TEST_MAIN})
  set_property(TEST vkernel_reproducible_disks APPEND PROPERTY ENVIRONMENT
                                                           OMP_NUM_THREADS=4)

  # topology::sort_interactions (storage::permute) on both storage backends: the
  # registry, and the compile time database (_static_db). The neighbor search
  # order, hence the interactions, must be reproducible (the neighbor lists are
  # sorted, whatever the number of threads), and the assertions of the library
  # headers must stay enabled in a Release build.
  new_test(NAME vkernel_sort_interactions_registry SOURCES
           sort_interactions.cpp ${SIMPLE_TEST_MAIN})
  new_test(
    NAME
    vkernel_sort_interactions_static_db
    SOURCES
    sort_interactions.cpp
    ${SIMPLE_TEST_MAIN}
    COMPILE_DEFINITIONS
    SORT_STATIC_DB)

  # the tests were written against the vkernel/soa flags, which include the
  # C++23 standard they need for std::print
  set(vkernel_ref_tests
      vkernel_bouncing_ball
      vkernel_bouncing_ball_registry
      vkernel_bouncing_ball_split
      vkernel_unbounded_bouncing_ball
      vkernel_bouncing_balls
      vkernel_bouncing_disk
      vkernel_bouncing_disk_line
      vkernel_bouncing_disks
      vkernel_fixed_disks
      vkernel_registry_remove
      vkernel_nslaws_multi_group
      vkernel_polydisperse_disks
      vkernel_cxd_batch_smoke
      vkernel_polytope_cxd
      vkernel_bouncing_polytopes
      vkernel_churn_polytopes
      vkernel_polytope_fsegment
      vkernel_polytopemesh
      vkernel_polytopefdisk
      vkernel_polytopedisk
      vkernel_rigid_shapes
      vkernel_mixed_bodies
      vkernel_reproducible_disks
      vkernel_polytope_on_fem
      vkernel_io_contact_points
      vkernel_sort_interactions_registry
      vkernel_sort_interactions_static_db)
  foreach(t IN LISTS vkernel_ref_tests)
    set_target_properties(${t} PROPERTIES COMPILE_FLAGS "${flags}")
  endforeach()
  set_target_properties(
    vkernel_sort_interactions_registry vkernel_sort_interactions_static_db
    PROPERTIES COMPILE_FLAGS "${flags} -UNDEBUG")

  set_target_properties(RtOsiTest PROPERTIES COMPILE_FLAGS "-Wall --std=c++23")
  set_target_properties(t1 PROPERTIES COMPILE_FLAGS "-Wall --std=c++23")
endif()
