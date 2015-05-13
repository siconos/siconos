#
# Here comes IO
#
PROJECT(SiconosIO CXX)

OPTION(WITH_BULLET "compilation with Bullet Bindings. Default = OFF" OFF)
OPTION(WITH_OCC "compilation with Occ Bindings. Default = OFF" OFF)
OPTION(WITH_MECHANISMS "compilation with Siconos/Mechanisms Bindings. Default = OFF" OFF)

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

# Generated
INSTALL(FILES ${CMAKE_BINARY_DIR}/src/SiconosFullGenerated.hpp DESTINATION include/Siconos/IO)

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

COMPILE_WITH(SiconosControl)
IF(SiconosControl_FOUND)
  SET(HAVE_SICONOS_CONTROL TRUE)
ENDIF(SiconosControl_FOUND)


COMPILE_WITH(GMP REQUIRED)

# contact detection
IF(WITH_BULLET)
  COMPILE_WITH(Bullet REQUIRED)
  IF(BULLET_FOUND)
    SET(HAVE_BULLET TRUE)
    MESSAGE( STATUS "Bullet include dirs : ${BULLET_INCLUDE_DIRS}" )
    SET(SiconosMechanics_DIRS ${SiconosMechanics_DIRS} contactDetection/bullet)
    IF(BULLET_USE_DOUBLE_PRECISION)
      APPEND_CXX_FLAGS("-DBT_USE_DOUBLE_PRECISION")
    ENDIF(BULLET_USE_DOUBLE_PRECISION)
  ENDIF(BULLET_FOUND)
ENDIF(WITH_BULLET)

#occ 
IF(WITH_OCC)
  SET(HAVE_OCC TRUE)
  MESSAGE( STATUS "Compilation with Occ.")
ENDIF()

IF(WITH_MECHANISMS)
  SET(HAVE_MECHANISMS TRUE)
  MESSAGE( STATUS "Compilation with Occ.")
ENDIF()

compile_with(BLAS REQUIRED)
compile_with(LAPACK REQUIRED)

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

FIND_PACKAGE(Boost 1.47 COMPONENTS serialization REQUIRED)

# cmake 2.8.12.2 (ubuntu 14.04) need this
SET(SiconosIO_LINK_LIBRARIES  "${SiconosIO_LINK_LIBRARIES};${Boost_LIBRARIES}")

FIND_PACKAGE(Boost 1.47 COMPONENTS filesystem REQUIRED)

# cmake 2.8.12.2 (ubuntu 14.04) need this
SET(SiconosIO_LINK_LIBRARIES  "${SiconosIO_LINK_LIBRARIES};${Boost_LIBRARIES}")

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
SET(GENERATED_TARGETS --targets=Kernel)
SET(GENERATED_HEADERS SiconosFullGenerated.hpp)
SET(GENERATED_DEPENDENCIES ${SiconosKernel_INCLUDE_DIRS}/SiconosKernel.hpp)
SET(GENERATED_INCLUDES -I${SiconosNumerics_INCLUDE_DIRS}
  -I${SiconosKernel_INCLUDE_DIRS})
IF(SiconosMechanics_FOUND)
  SET(GENERATED_TARGETS "${GENERATED_TARGETS},Mechanics")
  SET(GENERATED_DEPENDENCIES ${GENERATED_DEPENDENCIES}
    ${SiconosMechanics_INCLUDE_DIRS}/SpaceFilter.hpp)
  SET(GENERATED_INCLUDES ${GENERATED_INCLUDES} -I${SiconosMechanics_INCLUDE_DIRS})
ENDIF()
IF(SiconosControl_FOUND)
  SET(GENERATED_TARGETS "${GENERATED_TARGETS},Control")
  SET(GENERATED_DEPENDENCIES ${GENERATED_DEPENDENCIES}
    ${SiconosControl_INCLUDE_DIRS}/SiconosControl.hpp)
  SET(GENERATED_INCLUDES ${GENERATED_INCLUDES} -I${SiconosControl_INCLUDE_DIRS})
ENDIF(SiconosControl_FOUND)


ADD_CUSTOM_COMMAND(
  OUTPUT ${GENERATED_HEADERS}
  DEPENDS ${GENERATED_DEPENDENCIES}
  COMMAND ${CMAKE_SOURCE_DIR}/tools/builder.py ${GENERATED_TARGETS}
  ${GENERATED_INCLUDES})

# set library target
LIBRARY_PROJECT_SETUP()

# add dependencies for header generation
ADD_CUSTOM_TARGET(SerializersGeneration ALL
  DEPENDS 
  SiconosFullGenerated.hpp)

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

#end
CLOSE_PROJECT()

