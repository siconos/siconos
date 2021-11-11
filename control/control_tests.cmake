include(tools4tests)

if(WITH_TESTING)
  # ----  Control tests ----
  begin_tests(src/tests DEPS "numerics;kernel;CPPUNIT::CPPUNIT")

  if(HAS_FORTRAN)
    new_test(SOURCES PIDTest.cpp ${SIMPLE_TEST_MAIN})
    new_test(SOURCES SMCTest.cpp ${SIMPLE_TEST_MAIN})
    new_test(SOURCES ObserverTest.cpp ${SIMPLE_TEST_MAIN})
    new_test(SOURCES TwistingTest.cpp ${SIMPLE_TEST_MAIN})
  endif()
  
endif()
