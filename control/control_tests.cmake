include(tools4tests)

set(TEST_WRAP FALSE)

if(WITH_${COMPONENT}_TESTING)
  # We don't use COMPILE_WITH since we don't want to link cppunit with the
  # kernel library
  find_package(CppUnit REQUIRED)
  set(TEST_LIBS ${TEST_LIBS} ${CPPUNIT_LIBRARIES})
  set(TEST_INCLUDE_DIR ${TEST_INCLUDE_DIR} ${CPPUNIT_INCLUDE_DIR})
  
  # the main test driver
  SET(TEST_MAIN src/tests/TestMain.cpp)
  
  # Simulation tests
  BEGIN_TEST(src/tests)
  
  NEW_TEST(tests)
  
  IF(HAS_FORTRAN)
    NEW_TEST(ControlTests PIDTest.cpp SMCTest.cpp ObserverTest.cpp TwistingTest.cpp)
  ENDIF(HAS_FORTRAN)

  END_TEST()
  
endif()
