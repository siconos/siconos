#
# Here comes IO
#
PROJECT(SiconosIO CXX)

#
# VERSION defined in top-level CMakeLists.txt
#
SET(SiconosIO_VERSION "${VERSION}")

#
# Project directories 
#
SET(SiconosIO_DIRS
  .)

#
# Installation layout
#
SET(SiconosIO_INSTALL_INCLUDE_DIR include/Siconos/IO)



# config.h
INSTALL(FILES ${CMAKE_BINARY_DIR}/IOConfig.h DESTINATION include/Siconos/IO)

# boost patch
INSTALL(FILES 
  ${CMAKE_SOURCE_DIR}/src/boost/serialization/unordered_collections_load_imp.hpp
  ${CMAKE_SOURCE_DIR}/src/boost/serialization/unordered_collections_save_imp.hpp
  ${CMAKE_SOURCE_DIR}/src/boost/serialization/unordered_map.hpp
  ${CMAKE_SOURCE_DIR}/src/boost/serialization/unordered_set.hpp
  DESTINATION include/Siconos/IO/boost/serialization)

IF(WITH_TESTING)
  COMPILE_WITH(CppUnit REQUIRED)
ENDIF(WITH_TESTING)

COMPILE_WITH(SiconosNumerics REQUIRED)
COMPILE_WITH(SiconosKernel REQUIRED)
COMPILE_WITH(LibXml2 REQUIRED)
COMPILE_WITH(VTK)
IF(VTK_FOUND)
  MESSAGE(STATUS "Found vtk-${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION}")
  SET(HAVE_VTK TRUE)
  SET(SiconosIO_DIRS ${SiconosIO_DIRS} VTK)
  IF(VTK_USE_FILE)
    INCLUDE(${VTK_USE_FILE})
  ENDIF()
  SET(SiconosIO_LINK_LIBRARIES "${SiconosIO_LINK_LIBRARIES};vtkIO;vtkCommon;vtkGraphics;vtkRendering")
ENDIF()


FIND_PACKAGE(Boost COMPONENTS serialization REQUIRED)

FIND_PACKAGE(Boost COMPONENTS filesystem REQUIRED)

COMPILE_WITH(Boost REQUIRED)

# remove debug & optimized inside LINK_LIBRARIES...
FOREACH(_L ${SiconosIO_LINK_LIBRARIES})
  IF(NOT ${_L} STREQUAL "debug")
    IF(NOT ${_L} STREQUAL "optimized")
      LIST(APPEND _LL ${_L})
    ENDIF()
  ENDIF()
ENDFOREACH()
SET(SiconosIO_LINK_LIBRARIES ${_LL})

SET(SiconosIO_LINKER_LANGUAGE CXX)

# do not skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 

# the RPATH to be used when installing
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# don't add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

LIBRARY_PROJECT_SETUP()

SET(TEST_WRAP TRUE)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/test/result.ref test/result.ref COPYONLY)

IF(WITH_TESTING)

  ENABLE_TESTING()


  # the main test driver
  SET(TEST_MAIN test/TestMain.cpp)

  BEGIN_TEST(test)
  NEW_TEST(ioTests BasicTest.cpp KernelTest.cpp)
  END_TEST(test)

ENDIF(WITH_TESTING)

