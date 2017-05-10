
IF(BUILDNAME_OPTIONS)
  SET(BUILDNAME "${BUILDNAME}-${BUILDNAME_OPTIONS}")
ENDIF(BUILDNAME_OPTIONS)

IF(CMAKE_BUILD_TYPE)
  SET(BUILDNAME "${BUILDNAME}-${CMAKE_BUILD_TYPE}")
ENDIF(CMAKE_BUILD_TYPE)

FIND_PACKAGE(Git)

EXECUTE_PROCESS(COMMAND
      ${GIT_EXECUTABLE} rev-parse --is-inside-work-tree
      OUTPUT_VARIABLE SOURCE_IS_GIT_REPO
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})


IF(SOURCE_IS_GIT_REPO MATCHES "^true\n")

  MESSAGE(STATUS "In a Git repo")

  EXECUTE_PROCESS(COMMAND
    ${GIT_EXECUTABLE} log -n 1 --pretty=format:%h
    OUTPUT_VARIABLE SOURCE_ABBREV_GIT_SHA1
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

  MESSAGE(STATUS "ABBREV_GIT_SHA1 :: ${SOURCE_ABBREV_GIT_SHA1}")

  # buildname must stay static to have track possibility on CDash.
  # SET(BUILDNAME "${BUILDNAME}-${SOURCE_ABBREV_GIT_SHA1}")

  # Jenkins always works in detached head mode, so we have to work a bit more to get the branch name
  # sorry this does not work for my pull request
  # http://cdash-bipop.inrialpes.fr/viewConfigure.php?buildid=25783
  EXECUTE_PROCESS(COMMAND
    ${GIT_EXECUTABLE} branch -r --contains ${SOURCE_ABBREV_GIT_SHA1}
    OUTPUT_VARIABLE SOURCE_GIT_BRANCHES
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

  IF(SOURCE_GIT_BRANCHES)
    STRING(REGEX REPLACE "^.*origin/" "" SOURCE_GIT_BRANCHES2 ${SOURCE_GIT_BRANCHES})
    STRING(REGEX REPLACE "\n.*$" "" SOURCE_GIT_BRANCH ${SOURCE_GIT_BRANCHES2})
    
    MESSAGE(STATUS "GIT_BRANCH :: ${SOURCE_GIT_BRANCH}")
  ENDIF()

  # just imagine a world without Jenkins...
  IF(NOT SOURCE_GIT_BRANCH)
    EXECUTE_PROCESS(COMMAND
      ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
      OUTPUT_VARIABLE SOURCE_GIT_BRANCH
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    MESSAGE(STATUS "GIT BRANCH      :: ${SOURCE_GIT_BRANCH}")
    
    IF(NOT SOURCE_GIT_BRANCH MATCHES "^master$")
      SET(BUILDNAME "${BUILDNAME}-${SOURCE_GIT_BRANCH}")
    ENDIF()
  ENDIF()
ELSE()
  MESSAGE(STATUS "Not in a Git repo")
ENDIF()

set(BUILDNAME "${BUILDNAME}-${CMAKE_SYSTEM_NAME}-${CMAKE_HOST_SYSTEM_PROCESSOR}")
if(CROSSCOMPILING_LINUX_TO_WINDOWS)
  set(BUILDNAME "${BUILDNAME}-CrossCompilingFromLinuxToWindows")
endif()
  
# Tests coverage (taken from ViSp)
    
#
# Note: all of this is done with a recent cmake version (>2.6.0) with:
# cmake -DCMAKE_BUILD_TYPE=Profile
#
#
IF(WITH_TESTS_COVERAGE)
  # Add build options for test coverage. Currently coverage is only supported
  # on gcc compiler
  # Because using -fprofile-arcs with shared lib can cause problems like:
  # hidden symbol `__bb_init_func', we add this option only for static
  # library build
  SET(BUILD_SHARED_LIBS)
  SET(CMAKE_BUILD_TYPE Debug)
  INCLUDE(TestCXXAcceptsFlag)
  INCLUDE(CheckCCompilerFlag)
  IF(CMAKE_CXX_COMPILER_WORKS)
    CHECK_CXX_ACCEPTS_FLAG(-ftest-coverage CXX_HAVE_FTEST_COVERAGE)
    CHECK_CXX_ACCEPTS_FLAG(-fprofile-arcs CXX_HAVE_PROFILE_ARCS)
  ENDIF(CMAKE_CXX_COMPILER_WORKS)
  IF(CMAKE_C_COMPILER_WORKS)
    CHECK_C_COMPILER_FLAG(-ftest-coverage C_HAVE_FTEST_COVERAGE)
    CHECK_C_COMPILER_FLAG(-fprofile-arcs C_HAVE_PROFILE_ARCS)
  ENDIF(CMAKE_C_COMPILER_WORKS)
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
