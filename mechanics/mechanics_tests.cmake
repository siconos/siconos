include(tools4tests)

if(WITH_TESTING)

  # ---- Collision/native tests ----
  begin_tests(src/collision/native/test DEPS "numerics;kernel;CPPUNIT::CPPUNIT")
  new_test(SOURCES MultiBodyTest.cpp ${SIMPLE_TEST_MAIN})

  if(SICONOS_HAS_BULLET)
    begin_tests(src/collision/bullet/test DEPS "numerics;kernel;CPPUNIT::CPPUNIT")
    new_test(SOURCES  ContactTest.cpp ${SIMPLE_TEST_MAIN})
    new_test(SOURCES  Contact2dTest.cpp ${SIMPLE_TEST_MAIN})
  endif()
  
  if(WITH_OpenCASCADE)
    begin_tests(src/occ/test DEPS "numerics;kernel;CPPUNIT::CPPUNIT;OpenCASCADE::OpenCASCADE")
    new_test(SOURCES  OccTest.cpp ${SIMPLE_TEST_MAIN})
  endif()

endif()
