include(tools4tests)

if(WITH_${COMPONENT}_TESTING)
  find_package(CPPUNIT REQUIRED)
  
  set(SIMPLE_TEST_MAIN ${CMAKE_SOURCE_DIR}/kernel/src/utils/SiconosMemory/test/TestMain.cpp)

  # ----  Control tests ----
  begin_tests(src/tests DEPS "numerics;kernel;CPPUNIT::CPPUNIT")

  if(HAS_FORTRAN)
    new_test(SOURCES PIDTest.cpp ${SIMPLE_TEST_MAIN})
    new_test(SOURCES SMCTest.cpp ${SIMPLE_TEST_MAIN})
    new_test(SOURCES ObserverTest.cpp ${SIMPLE_TEST_MAIN})
    new_test(SOURCES TwistingTest.cpp ${SIMPLE_TEST_MAIN})
  endif()
  
endif()
