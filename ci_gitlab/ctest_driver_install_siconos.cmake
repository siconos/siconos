#-----------------------------------------------
# Ctest driver for siconos install.
# Target : continuous integration on gitlab-ci,
# aims at providing a proper install of siconos for a given configuration.
#
# Input variables :
# - SICONOS_INSTALL_DIR : where to install siconos. Default : /home/install-siconos
# - USER_FILE : user option file used by cmake to configure siconos. Default : siconos_conf.cmake.
#   Warning : always searched in siconos-tutorials/ci directory.
#   using common commands (lsb_release ...)
# ----------------------------------------------


message("--- Start conf for siconos ctest pipeline.")

# ============= setup  ================

# -- CI_PROJECT_DIR is a required environment variable --
# --> set by default for gitlab-ci, even inside the docker container
# --> unknown in docker container run with travis/siconos pipeline.

if(DEFINED ENV{TRAVIS})
  if($ENV{TRAVIS} STREQUAL true)
    set(CI_TRAVIS ON)
    set(ENV{CI_PROJECT_DIR} ${CTEST_SOURCE_DIRECTORY})
  endif()
elseif(DEFINED ENV{GITLAB_CI})
  if($ENV{GITLAB_CI} STREQUAL true)
    set(CI_GITLAB ON)
  endif()
endif()
  
if(NOT DEFINED ENV{CI_PROJECT_DIR} )
  message(FATAL_ERROR "Please set env variable CI_PROJECT_DIR to siconos sources directory (git repo).")
endif()

# -- Definition of all variables required for ctest --
include($ENV{CI_PROJECT_DIR}/ci_gitlab/ctest_tools.cmake)
if(CI_TRAVIS)
  list(APPEND CMAKE_MODULE_PATH ${CTEST_SOURCE_DIRECTORY}/ci_travis/cmake)
  list(APPEND CMAKE_MODULE_PATH ${CTEST_SOURCE_DIRECTORY}/ci_travis/config)
  list(APPEND CMAKE_MODULE_PATH ${CTEST_SOURCE_DIRECTORY}/ci_travis)
  include(Tools)

  # -- Get config --
  # i.e. set extra options/values (cmake -Doption=value ...)
  # either from file default.cmake or
  # from file CI_CONFIG.cmake
  # --> may set SICONOS_CMAKE_OPTIONS
  # --> may set DSICONOS_COMPONENTS
  # Rq : For Travis CI, we include cmake files while for gitlab CI we use
  # siconos user option file. Todo: one way to rule them all?

  if(CI_CONFIG)
    string(REPLACE "," ";" CI_CONFIG_LIST ${CI_CONFIG})
    foreach(_CI IN LISTS CI_CONFIG_LIST)
      include(${_CI})
    endforeach()
  else()
    set(CI_CONFIG default)
    include(${CI_CONFIG})
  endif()
endif()

# - Source dir and path to siconos install
if(NOT CTEST_SOURCE_DIRECTORY)
  set(CTEST_SOURCE_DIRECTORY $ENV{CI_PROJECT_DIR})
endif()

# - Top level build directory -
# If not specified : current dir.
if(NOT CTEST_BINARY_DIRECTORY)
  set(CTEST_BINARY_DIRECTORY .)
endif()

# Install dir (used as CMAKE_INSTALL_PREFIX)
if(NOT SICONOS_INSTALL_DIR)
  set(SICONOS_INSTALL_DIR /home/install-siconos/)
endif()
# Build name (for cdash)
if(NOT CTEST_BUILD_NAME)
  set_cdash_build_name()
endif()

if(USER_FILE)
  list(APPEND SICONOS_CMAKE_OPTIONS -DUSER_OPTIONS_FILE=${USER_FILE})
endif()

list(APPEND SICONOS_CMAKE_OPTIONS -DCMAKE_INSTALL_PREFIX=${SICONOS_INSTALL_DIR})
list(APPEND SICONOS_CMAKE_OPTIONS -DWITH_GIT=ON) # required to generate siconos-commit.txt to tag cdash build in the examples.

if(DEFINED ENV{OCE_INSTALL}) # set if oce has been installed using oce repo, in install_oce.sh
  message("Search oce in $ENV{OCE_INSTALL}.")
  list(APPEND SICONOS_CMAKE_OPTIONS -DOCE_DIR=$ENV{OCE_INSTALL})
endif()

# Parallel build only for siconos_install. For examples it leads to: warning: jobserver unavailable: using -j1. Add `+' to parent make rule.
#set(CTEST_MEMORYCHECK_SUPPRESSIONS_FILE ${CTEST_SOURCE_DIRECTORY}/cmake/valgrind.supp)

if(NOT CTEST_CMAKE_GENERATOR)
  set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
endif()

if(NOT CTEST_SITE)
  set_site_name()
endif()

if(NOT CTEST_BUILD_CONFIGURATION)
  set(CTEST_BUILD_CONFIGURATION "Release")
endif()

  
# Write a note file for cdash server.
# Content :
# - info. regarding the runner, the system
# - siconos config (user option file)
write_notes()

# =============  Run ctest steps ================
# Either one by one (to split ci jobs) if CTEST_MODE=Configure, Build, Test or
# all in a row if CTEST_MODE=all.
# Submit : only after test phase except if conf or build failed.

# - Configure -- 
if(${CTEST_MODE} STREQUAL "Configure" OR ${CTEST_MODE} STREQUAL "all")

  # Current testing model. Priority: 
  # Nightly -> set by scheduler on gricad-gitlab
  # Continuous -> set in .gitlab-ci.yml
  # Experimental : default
  if(NOT model)
    set(model Experimental)
  endif()

  ctest_start(${model})

  # Set CTEST_CONFIGURE_COMMAND to cmake followed by siconos options 
  set(CTEST_CONFIGURE_COMMAND ${CMAKE_COMMAND})
  foreach(option IN LISTS SICONOS_CMAKE_OPTIONS)
    set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} ${option}")
  endforeach()
  set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} ${CTEST_SOURCE_DIRECTORY}")

  message("\n\n=============== Start ctest_configure =============== ")
  message("- Configure command line :\n ${CTEST_CONFIGURE_COMMAND}\n")

  ctest_configure(
    RETURN_VALUE _RESULT
    CAPTURE_CMAKE_ERROR _STATUS
    #QUIET
    )
  post_ctest(PHASE Configure)
endif()
 
# - Build -
if(${CTEST_MODE} STREQUAL "Build" OR ${CTEST_MODE} STREQUAL "all")

  if(${CTEST_MODE} STREQUAL "Build")
    ctest_start(APPEND) # Restart from existing (configure step) cdash config
  endif()
  # --- Build ---

  message("\n\n=============== Start ctest_build =============== ")

  cmake_host_system_information(RESULT NP QUERY NUMBER_OF_LOGICAL_CORES)
  if(NOT ALLOW_PARALLEL_BUILD)
    set(NP 1)
  endif()
  ctest_build(
    FLAGS -j${NP}
    CAPTURE_CMAKE_ERROR _STATUS
    RETURN_VALUE _RESULT
    #QUIET if quiet, travis failed because of missing outputs during a long time ...
    )
  post_ctest(PHASE Build)
endif()

# - Test -
if(${CTEST_MODE} STREQUAL "Test" OR ${CTEST_MODE} STREQUAL "all")
  # -- Tests --
  
  if(${CTEST_MODE} STREQUAL "Test")
    ctest_start(APPEND)
  endif()
  message("\n\n=============== Start ctest_test (nbprocs = ${NP}) =============== ")
  ctest_test(
    #PARALLEL_LEVEL NP
    CAPTURE_CMAKE_ERROR _STATUS
    #SCHEDULE_RANDOM ON
    RETURN_VALUE _RESULT
    # QUIET
    )
  post_ctest(PHASE Test)

  if(WITH_MEMCHECK AND CTEST_COVERAGE_COMMAND)
    #find_program(CTEST_COVERAGE_COMMAND NAMES gcov)
    # set(CTEST_COVERAGE_COMMAND gcov)

    ctest_coverage(
      CAPTURE_CMAKE_ERROR COVERAGE_STATUS
      RETURN_VALUE COVERAGE_RESULT  
      )
  endif()

  if(WITH_MEMCHECK AND CTEST_MEMORYCHECK_COMMAND)
    ctest_memcheck()
  endif()
  
  if(CDASH_SUBMIT)
    ctest_submit(
      RETURN_VALUE RETURN_STATUS
      CAPTURE_CMAKE_ERROR SUBMISSION_STATUS
      )
    if(NOT SUBMISSION_STATUS EQUAL 0)
      message(WARNING " *** submission failure *** ")
    endif()
  endif()

  # -- memory check -- Skip this to 'enlight' submit process, since cdash inria is overbooked ...
  # if(CTEST_BUILD_CONFIGURATION MATCHES "Profiling")
  #   find_program(CTEST_MEMORYCHECK_COMMAND NAMES valgrind)
  #   set(CTEST_MEMORYCHECK_COMMAND_OPTIONS "--quiet --leak-check=full --show-leak-kinds=definite,possible --track-origins=yes --error-limit=no --gen-suppressions=all") 
  #   set(CTEST_MEMORYCHECK_COMMAND_OPTIONS "--quiet --leak-check=full --show-reachable=yes --error-limit=no --gen-suppressions=all") 
  #   ctest_memcheck(PARALLEL_LEVEL NP QUIET)
  # endif()

endif()


# ============= Summary =============
message(STATUS "\n============================================ Summary ============================================")
message(STATUS "CTest process for Siconos has ended.")
message(STATUS "Ctest model is: ${model}")
message(STATUS "Ctest mode was: ${CTEST_MODE}")
message(STATUS "Ctest executed on sources directory : ${CTEST_SOURCE_DIRECTORY}")
message(STATUS "CTEST_BINARY_DIRECTORY is: ${CTEST_BINARY_DIRECTORY}")
message(STATUS "CTEST_BUILD_CONFIGURATION is: ${CTEST_BUILD_CONFIGURATION}")
if(CDASH_SUBMIT)
  message(STATUS "Cdash server name: ${CTEST_DROP_SITE}/${CTEST_DROP_LOCATION}.")
  message(STATUS "Cdash build name: ${CTEST_BUILD_NAME}")
  message(STATUS "Cdash Site name: ${CTEST_SITE}")
endif()
message(STATUS "=================================================================================================\n")
