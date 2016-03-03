include(tools4tests)

if(WITH_${COMPONENT}_TESTING)

  # We don't use COMPILE_WITH since we don't want to link cppunit with the
  # kernel library
  find_package(CppUnit REQUIRED)
  set(TEST_LIBS ${TEST_LIBS} ${CPPUNIT_LIBRARIES})
  set(TEST_INCLUDE_DIR ${TEST_INCLUDE_DIR} ${CPPUNIT_INCLUDE_DIR})

  # the main test driver
  SET(TEST_MAIN src/contactDetection/basicBroadphase/test/TestMain.cpp)

  BEGIN_TEST(src/contactDetection/basicBroadphase/test)
  NEW_TEST(testMultiBody MultiBodyTest.cpp)
  END_TEST()

  BEGIN_TEST(src/proposed)
  NEW_TEST(testContact testContact.cpp)
  END_TEST()

  IF(WITH_OCC)
    BEGIN_TEST(src/occ/test)
    NEW_TEST(testOcc OccTest.cpp)
    END_TEST()
  ENDIF()

endif()
