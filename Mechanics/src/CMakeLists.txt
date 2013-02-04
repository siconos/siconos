#
# Here comes Mechanics
#
PROJECT(SiconosMechanics CXX)

#
# VERSION defined in top-level CMakeLists.txt
#
SET(SiconosMechanics_VERSION "${VERSION}")

#
# Project directories 
#
SET(SiconosMechanics_DIRS
  multiBodyTools
  multiBodyTools/2D
  multiBodyTools/3D
  multiBodyTools/joints)

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
    SET(SiconosKernel_DIRS ${SiconosKernel_DIRS} multiBodyTools/contactDetection/bullet)
    IF(BULLET_USE_DOUBLE_PRECISION)
      APPEND_FLAGS(CMAKE_CXX_FLAGS "-DBT_USE_DOUBLE_PRECISION")
    ENDIF(BULLET_USE_DOUBLE_PRECISION)
  ENDIF(BULLET_FOUND)
ENDIF(WITH_BULLET)

COMPILE_WITH(SiconosNumerics REQUIRED)
COMPILE_WITH(SiconosKernel REQUIRED)

COMPILE_WITH(LibXml2 REQUIRED)
# This is needed to overcome some shortcomings in the FindXml2 Module
# When we are crosscompiling to windows, we want to link external lib statically
# But libxml2 archive needs zlib, iconv and an obscure W32api lib. pkg-config
# gives this information, but it is not taken into account in the cmake Module
IF(CROSSCOMPILING_LINUX_TO_WINDOWS)
  FOREACH(_EXTRA_LIB ${PC_LIBXML_LIBRARIES})
    LIST(APPEND ${PROJECT_NAME}_LINK_LIBRARIES ${_EXTRA_LIB})
  ENDFOREACH()
ENDIF(CROSSCOMPILING_LINUX_TO_WINDOWS)

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
  SET(TEST_MAIN multiBodyTools/test/TestMain.cpp)

  BEGIN_TEST(test)
  NEW_TEST(testMultiBody MultiBodytest.cpp)
  END_TEST(test)

ENDIF(WITH_TESTING)

# print all features
FEATURE_SUMMARY(WHAT ALL)
