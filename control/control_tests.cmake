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
  
  IF(WITH_FORTRAN)
    NEW_TEST(PIDTest.cpp SMCTest.cpp ObserverTest.cpp)
  ENDIF()
  END_TEST()
  
endif()
