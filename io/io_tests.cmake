include(tools4tests)

if(WITH_${COMPONENT}_TESTING)
  if(WITH_SERIALIZATION) # Some extra boost components are required
    find_package(Boost 1.61 COMPONENTS serialization filesystem REQUIRED)
  endif()
  find_package(CPPUNIT REQUIRED)

  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/test/result.ref src/test/result.ref COPYONLY)
  # the main test driver
  set(SIMPLE_TEST_MAIN ${CMAKE_SOURCE_DIR}/kernel/src/utils/SiconosMemory/test/TestMain.cpp)


  begin_tests(src/test DEPS "CPPUNIT::CPPUNIT")
  
  if(WITH_SERIALIZATION)
    new_test_1(SOURCES BasicTest.cpp ${SIMPLE_TEST_MAIN})
    new_test_1(SOURCES KernelTest.cpp ${SIMPLE_TEST_MAIN})
  endif()
  
endif()
