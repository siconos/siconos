# http://www.cmake.org/pipermail/cmake/2006-October/011446.html
#
# Find the CppUnit includes and library
#
# This module defines
# CPPUNIT_INCLUDE_DIR, where to find tiff.h, etc.
# CPPUNIT_LIBRARIES, the libraries to link against to use CppUnit.
# CPPUNIT_FOUND, If false, do not try to use CppUnit.

# also defined, but not for general use are
# CPPUNIT_LIBRARY, where to find the CppUnit library.
# CPPUNIT_DEBUG_LIBRARY, where to find the CppUnit library in debug mode.
INCLUDE(FindPackageHandleStandardArgs)

FIND_PATH(CPPUNIT_INCLUDE_DIR cppunit/TestCase.h
  /usr/local/include
  /usr/include
)

# With Win32, important to have both
IF(WIN32)
  FIND_LIBRARY(CPPUNIT_LIBRARY cppunit_dll
    ${CPPUNIT_INCLUDE_DIR}/../lib
    /usr/local/lib
    /usr/lib)
  FIND_LIBRARY(CPPUNIT_DEBUG_LIBRARY cppunitd_dll
    ${CPPUNIT_INCLUDE_DIR}/../lib
    /usr/local/lib
    /usr/lib)
ELSE(WIN32)
  # On unix system, debug and release have the same name
  FIND_LIBRARY(CPPUNIT_LIBRARY cppunit
    ${CPPUNIT_INCLUDE_DIR}/../lib
    /usr/local/lib
    /usr/lib)
  FIND_LIBRARY(CPPUNIT_DEBUG_LIBRARY cppunit
    ${CPPUNIT_INCLUDE_DIR}/../lib
    /usr/local/lib
    /usr/lib)
ENDIF(WIN32)

IF(CPPUNIT_LIBRARY OR CPPUNIT_DEBUG_LIBRARY)
  IF(NOT CPPUNIT_INCLUDE_DIR)
    GET_FILENAME_COMPONENT(CPPUNIT_LIB_PATH ${CPPUNIT_LIBRARY} PATH)
    FIND_PATH(CPPUNIT_INCLUDE_DIR cppunit/TestCase.h
    ${CPPUNIT_LIB_PATH}/../include/)
  ENDIF()
ENDIF()

FIND_PACKAGE_HANDLE_STANDARD_ARGS(CPPUNIT
  REQUIRED_VARS CPPUNIT_LIBRARY)

IF(CPPUNIT_INCLUDE_DIR)
  IF(CPPUNIT_LIBRARY)
    SET(CPPUNIT_LIBRARIES ${CPPUNIT_LIBRARY} ${CMAKE_DL_LIBS})
    SET(CPPUNIT_DEBUG_LIBRARIES ${CPPUNIT_DEBUG_LIBRARY}
      ${CMAKE_DL_LIBS})
  ENDIF(CPPUNIT_LIBRARY)
ENDIF(CPPUNIT_INCLUDE_DIR)

IF(NOT CPPUNIT_LIBRARY)
  IF(CppUnit_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR
      "Required Cppunit library not found. Please specify library location in CPPUNIT_LIBRARY")
  ENDIF()

ENDIF(NOT CPPUNIT_LIBRARY)
