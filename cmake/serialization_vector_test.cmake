
macro(TEST_SERIALIZATION_VECTOR_BUG)
  if(NOT "${SERIALIZATION_VECTOR_TEST}")
    configure_file(cmake/serialization_vector_test.cpp
      ${CMAKE_BINARY_DIR}/CMakeFiles COPYONLY)
    TRY_RUN(SERIALIZATION_VECTOR_TEST_RUN SERIALIZATION_VECTOR_TEST_COMPILE
      ${CMAKE_BINARY_DIR}/CMakeFiles
      ${CMAKE_BINARY_DIR}/CMakeFiles/serialization_vector_test.cpp
      LINK_LIBRARIES "${Boost_SERIALIZATION_LIBRARY_DEBUG}"
      CMAKE_FLAGS
        "-DINCLUDE_DIRECTORIES=${Boost_INCLUDE_DIR}"
        "-DLINK_DIRECTORIES=${Boost_LIBRARY_DIR_DEBUG}")
    if("${SERIALIZATION_VECTOR_TEST_RUN}" EQUAL 0)
      set(SERIALIZATION_VECTOR_TEST TRUE CACHE BOOL
        "True if passed test for BOOST serialization vector bug")
    endif()
  endif()
  if(NOT "${SERIALIZATION_VECTOR_TEST}")
    message(FATAL_ERROR "Boost serialization vector bug detected.  (Known bug in BOOST 1.58, see https://svn.boost.org/trac/boost/ticket/11612)  Please use BOOST 1.59 or later, or disable WITH_SERIALIZATION.")
  endif()
endmacro(TEST_SERIALIZATION_VECTOR_BUG)
