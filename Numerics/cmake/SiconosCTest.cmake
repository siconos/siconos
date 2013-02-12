
IF(BUILDNAME_OPTIONS)
  SET(BUILDNAME "${BUILDNAME}-${BUILDNAME_OPTIONS}")
ENDIF(BUILDNAME_OPTIONS)

IF(CMAKE_BUILD_TYPE)
  SET(BUILDNAME "${BUILDNAME}-${CMAKE_BUILD_TYPE}")
ENDIF(CMAKE_BUILD_TYPE)

IF(PIPOL_IMAGE)
  SET(BUILDNAME "${BUILDNAME}-${PIPOL_IMAGE_NAME}")
  SET(SITE ${PIPOL_SITE})
ELSE(PIPOL_IMAGE)
  SET(BUILDNAME "${BUILDNAME}-${CMAKE_SYSTEM_NAME}")
  IF(CROSSCOMPILING_LINUX_TO_WINDOWS)
    SET(BUILDNAME "${BUILDNAME}-CrossCompilingFromLinux")
  ENDIF()
ENDIF(PIPOL_IMAGE)
  
# Tests coverage (taken from ViSp)
    
#
# Note: all of this is done with a recent cmake version (>2.6.0) with:
# cmake -DCMAKE_BUILD_TYPE=Profile
#
IF(WITH_TESTS_COVERAGE)
  # Add build options for test coverage. Currently coverage is only supported
  # on gcc compiler
  # Because using -fprofile-arcs with shared lib can cause problems like:
  # hidden symbol `__bb_init_func', we add this option only for static
  # library build
  SET(BUILD_SHARED_LIBS)
  SET(CMAKE_BUILD_TYPE Debug)
  CHECK_CXX_ACCEPTS_FLAG(-ftest-coverage CXX_HAVE_FTEST_COVERAGE)
  CHECK_CXX_ACCEPTS_FLAG(-fprofile-arcs CXX_HAVE_PROFILE_ARCS)
  CHECK_C_COMPILER_FLAG(-ftest-coverage C_HAVE_FTEST_COVERAGE)
  CHECK_C_COMPILER_FLAG(-fprofile-arcs C_HAVE_PROFILE_ARCS)
  IF(CXX_HAVE_FTEST_COVERAGE AND CXX_HAVE_PROFILE_ARCS)
    MESSAGE("Adding test coverage flags to CXX compiler : -ftest-coverage -fprofile-arcs")
    SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -ftest-coverage -fprofile-arcs")
    SET (CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -fprofile-arcs -ftest-coverage")
    SET (CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -fprofile-arcs -ftest-coverage")
  ENDIF(CXX_HAVE_FTEST_COVERAGE AND CXX_HAVE_PROFILE_ARCS)
  
  IF(C_HAVE_FTEST_COVERAGE)
    MESSAGE("Adding test coverage flags to C compiler : -ftest-coverage")
    SET(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -ftest-coverage -fprofile-arcs")
    SET (CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -fprofile-arcs -ftest-coverage")
    SET (CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -fprofile-arcs -ftest-coverage")
  ENDIF(C_HAVE_FTEST_COVERAGE)
  
ENDIF(WITH_TESTS_COVERAGE)
