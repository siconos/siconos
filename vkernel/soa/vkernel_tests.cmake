include(tools4tests)

if(WITH_TESTING)
  add_custom_target(vkernel-tests echo "Start vkernel tests")

  begin_tests(src/siconos/simul/tests DEPS
              "numerics;kernel;mechanics;CPPUNIT::CPPUNIT")
  new_test(SOURCES RtOsiTest.cpp ${SIMPLE_TEST_MAIN})
  # end_tests()

  begin_tests(src/siconos/storage/tests DEPS "vkernel;CPPUNIT::CPPUNIT")
  new_test(SOURCES StorageTest.cpp ${SIMPLE_TEST_MAIN})
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

  # the legacy targets were built with the vkernel/soa flags, which include the
  # C++23 standard the tests need for std::print
  set(vkernel_ref_tests
      vkernel_bouncing_ball
      vkernel_bouncing_ball_registry
      vkernel_bouncing_ball_split
      vkernel_unbounded_bouncing_ball
      vkernel_bouncing_balls
      vkernel_bouncing_disk
      vkernel_bouncing_disk_line
      vkernel_bouncing_disks
      vkernel_fixed_disks)
  foreach(t IN LISTS vkernel_ref_tests)
    set_target_properties(${t} PROPERTIES COMPILE_FLAGS "${flags}")
  endforeach()

  set_target_properties(RtOsiTest PROPERTIES COMPILE_FLAGS "-Wall --std=c++23")
  set_target_properties(t1 PROPERTIES COMPILE_FLAGS "-Wall --std=c++23")
endif()
