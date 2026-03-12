include(tools4tests)

if(WITH_TESTING)
  add_custom_target(control-tests echo "Start control tests")

  # ----  Control tests ----
  begin_tests(src/tests DEPS "numerics;kernel;CPPUNIT::CPPUNIT")

  if(HAS_FORTRAN)
    if(LPSOLVE_FOUND)
      # HAS_EXTREME_POINT_ALGO must be added to compile def if lpsolve has been found.
      list(APPEND compile_defs HAS_EXTREME_POINT_ALGO)
    endif()

    new_test(SOURCES PIDTest.cpp ${SIMPLE_TEST_MAIN})
    new_test(SOURCES SMCTest.cpp ${SIMPLE_TEST_MAIN} COMPILE_DEFINITIONS compile_defs)
    new_test(SOURCES ObserverTest.cpp ${SIMPLE_TEST_MAIN})
    new_test(SOURCES TwistingTest.cpp ${SIMPLE_TEST_MAIN} COMPILE_DEFINITIONS compile_defs)
  endif()
  
endif()
