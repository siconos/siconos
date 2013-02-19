#
# Here comes Mechanics
#
PROJECT(SiconosMechanics CXX)

OPTION(WITH_MECHANISMS "compilation with OpenCascade Bindings. Default = OFF" OFF)
OPTION(WITH_BULLET "compilation with Bullet Bindings. Default = OFF" OFF)

#
# VERSION defined in top-level CMakeLists.txt
#
SET(SiconosMechanics_VERSION "${VERSION}")

#
# Project directories 
#
SET(SiconosMechanics_DIRS
  bodies
  joints
  contactDetection/basicBroadphase)

IF(WITH_MECHANISMS)
  SET(OCE_TOOLKITS "TKernel"  "TKMath" "TKService" "TKV3d" "TKV2d" "TKBRep" "TKIGES" "TKSTL" "TKVRML" "TKSTEP" "TKSTEPAttr" "TKSTEP209" "TKSTEPBase" "TKShapeSchema" "TKGeomBase" "TKGeomAlgo" "TKG3d" "TKG2d" "TKXSBase" "TKPShape" "TKShHealing" "TKHLR" "TKTopAlgo" "TKMesh" "TKPrim" "TKCDF" "TKBool" "TKBO" "TKFillet" "TKOffset")
  FIND_PACKAGE(OpenCASCADE REQUIRED COMPONENTS ${OCE_TOOLKITS})
  COMPILE_WITH(OpenCASCADE)

  IF(OpenCASCADE_FOUND) 
    message(STATUS "OpenCASCADE_INCLUDE_DIR = " ${OpenCASCADE_INCLUDE_DIR})
    message(STATUS "OpenCASCADE_LIBRARIES = " ${OpenCASCADE_LIBRARIES})
    message(STATUS "OpenCASCADE_LINK_DIRECTORY = " ${OpenCASCADE_LINK_DIRECTORY})
    include_directories(${OpenCASCADE_INCLUDE_DIR})
  ENDIF()

  SET(SiconosMechanics_DIRS ${SiconosMechanics_DIRS}
    mechanisms/CADMBTB
    mechanisms/MBTB
    mechanisms/TIMERMBTB)

ENDIF(WITH_MECHANISMS)



#
# Installation layout
#
SET(SiconosMechanics_INSTALL_INCLUDE_DIR include/Siconos/Mechanics)



# config.h
INSTALL(FILES ${CMAKE_BINARY_DIR}/MechanicsConfig.h DESTINATION include/Siconos/Mechanics)

IF(WITH_TESTING)
  COMPILE_WITH(CppUnit REQUIRED)
ENDIF(WITH_TESTING)

# contact detection
IF(WITH_BULLET)
  COMPILE_WITH(Bullet REQUIRED)
  IF(BULLET_FOUND)
    SET(HAVE_BULLET TRUE)
    MESSAGE( STATUS " Bullet include dirs : ${BULLET_INCLUDE_DIRS}" )
    SET(SiconosMechanics_DIRS ${SiconosMechanics_DIRS} contactDetection/bullet)
    IF(BULLET_USE_DOUBLE_PRECISION)
      APPEND_FLAGS(CMAKE_CXX_FLAGS "-DBT_USE_DOUBLE_PRECISION")
    ENDIF(BULLET_USE_DOUBLE_PRECISION)
  ENDIF(BULLET_FOUND)
ENDIF(WITH_BULLET)

COMPILE_WITH(SiconosNumerics REQUIRED)
COMPILE_WITH(SiconosKernel REQUIRED)

COMPILE_WITH(BLAS REQUIRED)
COMPILE_WITH(LAPACK REQUIRED)

COMPILE_WITH(Boost REQUIRED)

# remove debug & optimized inside LINK_LIBRARIES...
FOREACH(_L ${SiconosMechanics_LINK_LIBRARIES})
  IF(NOT ${_L} STREQUAL "debug")
    IF(NOT ${_L} STREQUAL "optimized")
      LIST(APPEND _LL ${_L})
    ENDIF()
  ENDIF()
ENDFOREACH()
SET(SiconosMechanics_LINK_LIBRARIES ${_LL})

SET(SiconosMechanics_LINKER_LANGUAGE CXX)

LIBRARY_PROJECT_SETUP()


IF(WITH_TESTING)

  ENABLE_TESTING()

  # the main test driver
  SET(TEST_MAIN contactDetection/basicBroadphase/test/TestMain.cpp)

  BEGIN_TEST(contactDetection/basicBroadphase/test)
  NEW_TEST(testMultiBody MultiBodyTest.cpp)
  END_TEST()

ENDIF(WITH_TESTING)

# print all features
FEATURE_SUMMARY(WHAT ALL)
