include(tools4tests)

if(WITH_${COMPONENT}_TESTING)
  # We don't use COMPILE_WITH since we don't want to link cppunit with the
  # kernel library
  find_package(CppUnit REQUIRED)
  set(TEST_LIBS ${TEST_LIBS} ${CPPUNIT_LIBRARIES})
  set(TEST_INCLUDE_DIR ${TEST_INCLUDE_DIR} ${CPPUNIT_INCLUDE_DIR})

  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/src/test/result.ref src/test/result.ref COPYONLY)
  # the main test driver
  SET(TEST_MAIN src/test/TestMain.cpp)

  BEGIN_TEST(src/test)

  SET(TEST_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/src)

  IF(WITH_SERIALIZATION)
    NEW_TEST(ioTests BasicTest.cpp KernelTest.cpp)
  ENDIF()

  END_TEST(test)

endif()
