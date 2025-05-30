include(tools4tests)

if(WITH_TESTING)
  add_custom_target(vkernel-tests echo "Start vkernel tests")

  begin_tests(src/siconos/simul/tests DEPS "numerics;kernel;mechanics;CPPUNIT::CPPUNIT")
  new_test(SOURCES RtOsiTest.cpp ${SIMPLE_TEST_MAIN})

  #end_tests()
endif()
