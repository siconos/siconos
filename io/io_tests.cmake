include(tools4tests)

if(WITH_TESTING)
  
  add_custom_target(io-tests echo "Start io tests")  
  if(WITH_SERIALIZATION)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/test/result.ref src/test/result.ref COPYONLY)
    begin_tests(src/test DEPS "CPPUNIT::CPPUNIT")
    new_test(SOURCES BasicTest.cpp ${SIMPLE_TEST_MAIN})
    new_test(SOURCES KernelTest.cpp ${SIMPLE_TEST_MAIN})
  endif()
  
endif()
