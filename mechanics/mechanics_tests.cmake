include(tools4tests)

if(WITH_${COMPONENT}_TESTING)

  # We don't use COMPILE_WITH since we don't want to link cppunit with the
  # kernel library
  find_package(CppUnit REQUIRED)
  set(TEST_LIBS ${TEST_LIBS} ${CPPUNIT_LIBRARIES})
  set(TEST_INCLUDE_DIR ${TEST_INCLUDE_DIR} ${CPPUNIT_INCLUDE_DIR})

  # the main test driver
  SET(TEST_MAIN src/collision/native/test/TestMain.cpp)

  BEGIN_TEST(src/collision/native/test)
  NEW_TEST(testMultiBody MultiBodyTest.cpp)
  END_TEST()

  BEGIN_TEST(src/collision/bullet/test)
  NEW_TEST(testContact ContactTest.cpp)
  END_TEST()

  IF(WITH_OCC)
    BEGIN_TEST(src/occ/test)
    NEW_TEST(testOcc OccTest.cpp)
    END_TEST()
  ENDIF()
  
endif()
