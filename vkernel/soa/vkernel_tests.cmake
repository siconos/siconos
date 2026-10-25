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

  set_target_properties(RtOsiTest PROPERTIES COMPILE_FLAGS "-Wall --std=c++23")
  set_target_properties(t1 PROPERTIES COMPILE_FLAGS "-Wall --std=c++23")
endif()
