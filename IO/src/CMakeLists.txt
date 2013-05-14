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

COMPILE_WITH(SiconosMechanics)
IF(SiconosMechanics_FOUND)
  SET(HAVE_SICONOS_MECHANICS TRUE)
ENDIF(SiconosMechanics_FOUND)

COMPILE_WITH(LibXml2 REQUIRED)
COMPILE_WITH(GMP REQUIRED)


IF(WITH_VTK)
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
ENDIF()
compile_with(BLAS REQUIRED)
compile_with(LAPACK REQUIRED)

FIND_PACKAGE(Boost 1.47 COMPONENTS serialization REQUIRED)

FIND_PACKAGE(Boost 1.47 COMPONENTS filesystem REQUIRED)

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

# header generation
ADD_CUSTOM_COMMAND(
  OUTPUT SiconosKernelGenerated.hpp
  DEPENDS ${SiconosKernel_INCLUDE_DIRS}/SiconosKernel.hpp
  COMMAND ${CMAKE_SOURCE_DIR}/tools/builder.py --target=Kernel
  -I${SiconosNumerics_INCLUDE_DIRS} -I${SiconosKernel_INCLUDE_DIRS}
  -I${LIBXML2_INCLUDE_DIR})

IF(SiconosMechanics_FOUND)
  ADD_CUSTOM_COMMAND(
    OUTPUT SiconosMechanicsGenerated.hpp
    DEPENDS ${SiconosMechanics_INCLUDE_DIRS}/SpaceFilter.hpp
    COMMAND ${CMAKE_SOURCE_DIR}/tools/builder.py --target=Mechanics
    -I${SiconosNumerics_INCLUDE_DIRS} -I${SiconosKernel_INCLUDE_DIRS} -I${SiconosMechanics_INCLUDE_DIRS}
    -I${LIBXML2_INCLUDE_DIR})
ENDIF()

# set library target
LIBRARY_PROJECT_SETUP()

# add dependencies for header generation
IF(SiconosMechanics_FOUND)
  ADD_CUSTOM_TARGET(SerializersGeneration ALL
    DEPENDS 
    SiconosKernelGenerated.hpp 
    SiconosMechanicsGenerated.hpp)
ELSE()
  ADD_CUSTOM_TARGET(SerializersGeneration ALL
    DEPENDS SiconosKernelGenerated.hpp)
ENDIF()

IF(BUILD_SHARED_LIBS)
  ADD_DEPENDENCIES(SiconosIO_shared SerializersGeneration)
ENDIF()
IF(BUILD_STATIC_LIBS)
  ADD_DEPENDENCIES(SiconosIO_static SerializersGeneration)
ENDIF()

# tests
IF(WITH_TESTING)

  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/test/result.ref test/result.ref COPYONLY)

  ENABLE_TESTING()

  SET(LD_LIBRARY_PATH
    ".:${CMAKE_CURRENT_BINARY_DIR}/plugin:${CMAKE_CURRENT_BINARY_DIR}/plugin/test")
  SET(PATH_FOR_PLUGIN
    ".\;${CMAKE_CURRENT_BINARY_DIR}/plugin\;${CMAKE_CURRENT_BINARY_DIR}/plugin/test")
  SET(SICONOS_LIB_PATH "${CMAKE_INSTALL_PREFIX}/lib${LIB_SUFFIX}")

  # the main test driver
  SET(TEST_MAIN test/TestMain.cpp)

  BEGIN_TEST(test)

  SET(TEST_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR})

  NEW_TEST(ioTests BasicTest.cpp KernelTest.cpp)
  END_TEST(test)

ENDIF(WITH_TESTING)

# print all features
FEATURE_SUMMARY(WHAT ALL)
