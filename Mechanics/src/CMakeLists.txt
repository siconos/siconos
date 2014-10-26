#
# Here comes Mechanics
#
PROJECT(SiconosMechanics CXX)

OPTION(WITH_OCC "compilation with OpenCascade Bindings. Default = OFF" OFF)
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
  .
  bodies
  joints
  contactDetection/basicBroadphase)

IF(WITH_OCC)
  COMPILE_WITH(OCE)
  LIST(REMOVE_ITEM SiconosMechanics_LINK_LIBRARIES DRAWEXE)
  SET(SiconosMechanics_DIRS ${SiconosMechanics_DIRS} occ)
ENDIF()

IF(WITH_MECHANISMS)
  # SET(OCE_TOOLKITS "TKernel"  "TKMath" "TKService" "TKV3d"  "TKBRep" "TKIGES" "TKSTL" "TKVRML" "TKSTEP" "TKSTEPAttr" "TKSTEP209" "TKSTEPBase" "TKShapeSchema" "TKGeomBase" "TKGeomAlgo" "TKG3d" "TKG2d" "TKXSBase" "TKPShape" "TKShHealing" "TKHLR" "TKTopAlgo" "TKMesh" "TKPrim" "TKCDF" "TKBool" "TKBO" "TKFillet" "TKOffset")
  
  # message(STATUS "Searching for OCE ....")
  # find_package(OCE COMPONENTS ${OCE_TOOLKITS})
  
  # if(OCE_FOUND)
  #   message(STATUS "Found OCE version ${OCE_VERSION}")
  #   if(NOT OCE_ALL_FOUND)
  #     set(OCE_FOUND false)
  #     message(WARNING "Ignoring OCE installation due to missing toolkit(s): ${OCE_MISSING_TOOLKITS}")
  #   endif(NOT OCE_ALL_FOUND)
  # endif(OCE_FOUND)
  
  # if(OCE_FOUND)
  #   # Include files reside in ${OCE_INCLUDE_DIRS};
  #   include_directories(${OCE_INCLUDE_DIRS})
  #   LINK_LIBRARIES(CADMBTB ${OCE_LIBRARIES})
  #   # We do not need library path, they will be automatically imported.
  # else(OCE_FOUND)
  #   # OCE not found; either it is not found and user
  #   # has to set OCE_DIR to the directory containing
  #   # OCEConfig.cmake, or OCE is not installed and we
  #   # try to find OpenCascade files.
  #   message(STATUS "OCE not found.  Try to find OpenCascade files.")
    
    

  #   FIND_PACKAGE(OpenCASCADE REQUIRED COMPONENTS ${OCE_TOOLKITS})
  #   COMPILE_WITH(OpenCASCADE)
    
  #   IF(OpenCASCADE_FOUND) 
  #     message(STATUS "OpenCASCADE_INCLUDE_DIR = " ${OpenCASCADE_INCLUDE_DIR})
  #     message(STATUS "OpenCASCADE_LIBRARIES = " ${OpenCASCADE_LIBRARIES})
  #     message(STATUS "OpenCASCADE_LINK_DIRECTORY = " ${OpenCASCADE_LINK_DIRECTORY})
  #     include_directories(${OpenCASCADE_INCLUDE_DIR})
  #   ELSE(OpenCASCADE_FOUND)
  #     MESSAGE(STATUS "OpenCascade Libraries not found in standard paths.")
  #   ENDIF(OpenCASCADE_FOUND)
    
  # endif(OCE_FOUND)


  add_subdirectory(mechanisms)
  # SET(SiconosMechanics_DIRS ${SiconosMechanics_DIRS}
  #   mechanisms/CADMBTB
  #   mechanisms/MBTB
  #   mechanisms/TIMERMBTB
  #   mechanisms/frontEnd)

ENDIF(WITH_MECHANISMS)

SET(SICONOS_LIB_PATH "${CMAKE_INSTALL_PREFIX}/lib${LIB_SUFFIX}")

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
      APPEND_CXX_FLAGS("-DBT_USE_DOUBLE_PRECISION")
    ENDIF(BULLET_USE_DOUBLE_PRECISION)
  ENDIF(BULLET_FOUND)
ENDIF(WITH_BULLET)

INCLUDE(SiconosCompilerFlags)

compile_with(BLAS REQUIRED)
compile_with(LAPACK REQUIRED)

COMPILE_WITH(SiconosNumerics REQUIRED)
COMPILE_WITH(SiconosKernel REQUIRED)

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

IF(CMAKE_SYSTEM_NAME MATCHES "Windows")
  STRING(REGEX REPLACE "\\\\" "/" ENV_PATH "$ENV{PATH}")
  STRING(REGEX REPLACE "\;" "\\\;" ENV_PATH "${ENV_PATH}")
ENDIF()

LIBRARY_PROJECT_SETUP()


IF(WITH_TESTING)

  ENABLE_TESTING()

  # the main test driver
  SET(TEST_MAIN contactDetection/basicBroadphase/test/TestMain.cpp)

  BEGIN_TEST(contactDetection/basicBroadphase/test)
  NEW_TEST(testMultiBody MultiBodyTest.cpp)
  END_TEST()

  IF(WITH_OCC)
    BEGIN_TEST(occ/test)
    NEW_TEST(testOcc OccTest.cpp)
    END_TEST()
  ENDIF()

ENDIF(WITH_TESTING)

# print all features
FEATURE_SUMMARY(WHAT ALL)
