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

IF(WITH_TESTING)
  COMPILE_WITH(CppUnit REQUIRED)
ENDIF(WITH_TESTING)

FIND_PACKAGE(Boost COMPONENTS serialization)

COMPILE_WITH(Boost REQUIRED)

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

IF(WITH_TESTING)

  # the main test driver
  SET(TEST_MAIN Test/TestMain.cpp)

  BEGIN_TEST(Test)
  NEW_TEST(whole WholeTest.cpp)
  END_TEST(Test)

ENDIF(WITH_TESTING)
