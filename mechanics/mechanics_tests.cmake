include(tools4tests)

if(WITH_TESTING)

  # ---- Collision/native tests ----
  begin_tests(src/collision/native/test DEPS "numerics;kernel;CPPUNIT::CPPUNIT")
  new_test(SOURCES MultiBodyTest.cpp ${SIMPLE_TEST_MAIN})

  if(WITH_BULLET)
    begin_tests(src/collision/bullet/test DEPS "numerics;kernel;CPPUNIT::CPPUNIT")
    new_test(SOURCES  ContactTest.cpp ${SIMPLE_TEST_MAIN})
    new_test(SOURCES  Contact2dTest.cpp ${SIMPLE_TEST_MAIN})
  endif()
  
  if(WITH_OCE)
    begin_tests(src/occ/test DEPS "numerics;kernel;CPPUNIT::CPPUNIT;${OCE_LIBRARIES}")
    new_test(SOURCES  OccTest.cpp ${SIMPLE_TEST_MAIN})
  endif()

endif()
