#----------------------------------------------------------------
# Ctest driver : command lines common to all
# ctest drivers.
#
# Usage, add at the end of your driver file :
# include(<path-to>/ctest_common.cmake)
# See for instance ctest_driver_install_siconos.cmake
# ---------------------------------------------------------------

# --- Configure setup ---
# - Top level build directory -
# If not specified : current dir.
if(NOT CTEST_BINARY_DIRECTORY)
  set(CTEST_BINARY_DIRECTORY .)
endif()

# Current testing model. Priority: 
# Nightly -> set by scheduler on gricad-gitlab
# Continuous -> set in .gitlab-ci.yml
# Experimental : default
if(NOT model)
  set(model Experimental)
endif()

if(NOT CTEST_CMAKE_GENERATOR)
  set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
endif()

# -- Query host system information --
#include(cmake_host_system_information)
cmake_host_system_information(RESULT hostname QUERY HOSTNAME)
cmake_host_system_information(RESULT fqdn QUERY FQDN)
cmake_host_system_information(RESULT NP QUERY NUMBER_OF_LOGICAL_CORES)
set(CTEST_BUILD_FLAGS -j${NP})

if(${CMAKE_VERSION} VERSION_GREATER "3.10.3") 
  cmake_host_system_information(RESULT osname QUERY OS_NAME)
  cmake_host_system_information(RESULT osrelease QUERY OS_RELEASE)
  cmake_host_system_information(RESULT osversion QUERY OS_VERSION)
  cmake_host_system_information(RESULT osplatform QUERY OS_PLATFORM)
else()
  set(osname ${CMAKE_SYSTEM_NAME})
  set(osversion ${CMAKE_SYSTEM_VERSION})
  set(osplatform ${CMAKE_SYSTEM_PROCESSOR})
endif()

# With gitlab-ci, runner name is too long and useless ...
string(FIND ${hostname} "runner-" on_ci) 
if(on_ci GREATER -1)
  set(hostname "gitlab-ci runner on $ENV{CI_RUNNER_DESCRIPTION}")
endif()

# Host description
if(NOT OSNAME)
  set(OSNAME ${osname}) # Use -DOSNAME=docker_image name on CI
endif()
if(NOT CTEST_SITE)
  set(CTEST_SITE "${OSNAME} ${osrelease}, ${osplatform}, ${hostname}")
endif()

ctest_start(${model})

# Set CTEST_CONFIGURE_COMMAND to cmake followed by siconos options 
set(CTEST_CONFIGURE_COMMAND ${CMAKE_COMMAND})
foreach(option ${SICONOS_CMAKE_OPTIONS})
  set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} ${option}")
endforeach()
set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} ${CTEST_SOURCE_DIRECTORY}")

message("\n\n=============== Start ctest_configure =============== ")
message("- Configure command line :\n ${CTEST_CONFIGURE_COMMAND}\n")
if(CDASH_SUBMIT)
  message("- Results will be submitted to cdash server, ${CTEST_DROP_SITE}.")
else()
  message("- Results won't be submitted to a cdash server.\n")
endif()

if(${CMAKE_VERSION} VERSION_GREATER "3.6.3") 
  ctest_configure(
    RETURN_VALUE CONFIGURE_RESULT
    CAPTURE_CMAKE_ERROR CONFIGURE_STATUS
    QUIET
    )
else()
  ctest_configure(
    RETURN_VALUE CONFIGURE_RESULT
    QUIET)
  set(CONFIGURE_STATUS ${CONFIGURE_RESULT})
endif()

message("=============== End of ctest_configure =============== ")
message("------> Configure status/result : ${CONFIGURE_STATUS}/${CONFIGURE_RESULT}")
if(NOT CONFIGURE_STATUS EQUAL 0 OR NOT CONFIGURE_RESULT EQUAL 0)
  if(CDASH_SUBMIT)
    ctest_submit(PARTS Configure)
  endif()
  message(FATAL_ERROR "\n\n *** Configure (cmake) process failed *** \n\n")
endif()

# --- Build ---

if(NOT CTEST_BUILD_CONFIGURATION)
  set(CTEST_BUILD_CONFIGURATION "Profiling")
endif()

message("\n\n=============== Start ctest_build =============== ")

if(${CMAKE_VERSION} VERSION_GREATER "3.6.3") 
  ctest_build(
      PROJECT_NAME ${current_project}
      CAPTURE_CMAKE_ERROR BUILD_STATUS
      RETURN_VALUE BUILD_RESULT
      QUIET
      )
else()
  ctest_build(
      PROJECT_NAME ${current_project}
      RETURN_VALUE BUILD_RESULT
      QUIET
      )
    set(BUILD_STATUS ${BUILD_RESULT})
  endif()
message("=============== End of ctest_build =============== ")
message("------> Build status/result : ${BUILD_STATUS}/${BUILD_RESULT}")
if(NOT BUILD_STATUS EQUAL 0 OR NOT BUILD_RESULT EQUAL 0)
  if(CDASH_SUBMIT)
    ctest_submit(PARTS Configure Build)
  endif()
  message(FATAL_ERROR " *** Build (make) process failed *** ")
endif()


# -- Tests --
message("\n\n=============== Start ctest_test (nbprocs = ${NP}) =============== ")
if(${CMAKE_VERSION} VERSION_GREATER "3.6.3") 
  ctest_test(
    PARALLEL_LEVEL NP
    CAPTURE_CMAKE_ERROR TEST_STATUS
    SCHEDULE_RANDOM ON
    RETURN_VALUE TEST_RESULT
    QUIET
    )
else()
  ctest_test(
    PARALLEL_LEVEL NP
    RETURN_VALUE TEST_STATUS
    SCHEDULE_RANDOM ON
    QUIET
    )
  set(TEST_STATUS ${TEST_RESULT})
endif()
message("=============== End of ctest_test =============== ")
message("------> Test status/result : ${TEST_STATUS}/${TEST_RESULT}")
# error status check later, we try to submit even if tests failed.

# -- memory check -- Skip this to 'enlight' submit process, since cdash inria is overbooked ...
# if(CTEST_BUILD_CONFIGURATION MATCHES "Profiling")
#   find_program(CTEST_MEMORYCHECK_COMMAND NAMES valgrind)
#   set(CTEST_MEMORYCHECK_COMMAND_OPTIONS "--quiet --leak-check=full --show-leak-kinds=definite,possible --track-origins=yes --error-limit=no --gen-suppressions=all") 
#   set(CTEST_MEMORYCHECK_COMMAND_OPTIONS "--quiet --leak-check=full --show-reachable=yes --error-limit=no --gen-suppressions=all") 
#   ctest_memcheck(PARALLEL_LEVEL NP QUIET)
# endif()
if(NOT CDASH_SUBMIT)
  return()
endif()

# -- Submission to cdash --
message("\n\n=============== Start ctest_submit =============== ")
ctest_submit(
  CAPTURE_CMAKE_ERROR  SUBMISSION_STATUS
  RETRY_COUNT 4 # Retry 4 times, if submission failed ...)
  RETRY_DELAY 1 # seconds
  )
message("=============== End of ctest_test =============== ")

# ============= Summary =============
message(STATUS "\n============================================ Summary ============================================")
message(STATUS "CTest process for ${current_project} has ended.")
message(STATUS "Ctest model is: ${model}")
message(STATUS "Ctest executed on sources directory : ${CTEST_SOURCE_DIRECTORY}")
message(STATUS "Build name (cdash) : ${CTEST_BUILD_NAME}")
message(STATUS "Site (cdash) : ${CTEST_SITE}")
message(STATUS "=================================================================================================\n")

# tests failed?
if(NOT TEST_STATUS EQUAL 0 OR NOT TEST_RESULT EQUAL 0)
  message(FATAL_ERROR " *** test failure *** ")
endif()

# -- Submission failed? --
if(NOT SUBMISSION_STATUS EQUAL 0)
  message(WARNING " *** submission failure *** ")
endif()
