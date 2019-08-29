include(tools4tests)

if(WITH_${COMPONENT}_TESTING)
  find_package(CPPUNIT REQUIRED)
  set(TEST_LIBS ${TEST_LIBS} ${CPPUNIT_LIBRARIES})
  set(TEST_INCLUDE_DIR ${TEST_INCLUDE_DIR} ${CPPUNIT_INCLUDE_DIR})

  # the main test driver
  set(SIMPLE_TEST_MAIN ${CMAKE_SOURCE_DIR}/kernel/src/utils/SiconosMemory/test/TestMain.cpp)

  # ---- Collision/native tests ----
  begin_tests(src/collision/native/test DEPS "numerics;kernel;CPPUNIT::CPPUNIT")
  new_test(SOURCES MultiBodyTest.cpp ${SIMPLE_TEST_MAIN})

  if(WITH_BULLET)
    begin_tests(src/collision/bullet/test DEPS "numerics;kernel;CPPUNIT::CPPUNIT")
    new_test(SOURCES  ContactTest.cpp ${SIMPLE_TEST_MAIN})
    new_test(SOURCES  Contact2dTest.cpp ${SIMPLE_TEST_MAIN})
  endif()
  
  if(WITH_OCE)
    begin_tests(src/occ/test DEPS "numerics;kernel;CPPUNIT::CPPUNIT")
    new_test(SOURCES  OccTest.cpp ${SIMPLE_TEST_MAIN})
  endif()

endif()
