#-----------------------------------------------
# Ctest driver for siconos install.
# Target : continuous integration on gitlab-ci,
# aims at providing a proper install of siconos for a given configuration.
#
# Input variables :
# - SICONOS_INSTALL_DIR : where to install siconos. Default : ../install-siconos
# - USER_FILE : user option file used by cmake to configure siconos. Default : siconos_conf.cmake.
#   Warning : always searched in siconos-tutorials/ci directory.
# - OSNAME : host system name (used to qualify cdash build). If not set, try to catch info
#   using common commands (lsb_release ...)
# ----------------------------------------------

# used as CMAKE_INSTALL_PREFIX
if(NOT SICONOS_INSTALL_DIR)
  set(SICONOS_INSTALL_DIR ../install-siconos/)
endif()

# -- job : build and install siconos --
message("--- Start conf for siconos install.")
# - Source dir and path to siconos install
if(NOT CTEST_SOURCE_DIRECTORY)
  set(CTEST_SOURCE_DIRECTORY ..)# $ENV{CI_PROJECT_DIR})
endif()

# Build name (for cdash)
if(NOT CTEST_BUILD_NAME)
  # Get hash for commit of current version of Siconos
  # Saved by CI in CI_COMMIT_SHORT_SHA.
  include(${CTEST_SOURCE_DIRECTORY}/cmake/SiconosVersion.cmake)
  set(CTEST_BUILD_NAME "Siconos (${SICONOS_VERSION}-devel, branch/commit=$ENV{TRAVIS_BRANCH}/$ENV{TRAVIS_COMMIT})")
  if(EXTRA_NAME)
    set(CTEST_BUILD_NAME "${CTEST_BUILD_NAME} - ${EXTRA_NAME}.")
  endif()
endif()

set(CMAKE_MODULE_PATH ${CTEST_SOURCE_DIRECTORY}/ci_travis/cmake;${CTEST_SOURCE_DIRECTORY}/ci_travis/config;${CTEST_SOURCE_DIRECTORY}/ci_travis)
include(Tools)

# -- Get config --
# i.e. set extra options/values (cmake -Doption=value ...)
# either from file default.cmake or
# from file CI_CONFIG.cmake
# --> may set SICONOS_CMAKE_OPTIONS
# --> may set DSICONOS_COMPONENTS
if(CI_CONFIG)
  string(REPLACE "," ";" CI_CONFIG_LIST ${CI_CONFIG})
  foreach(_CI IN LISTS CI_CONFIG_LIST)
    include(${_CI})
  endforeach()
else()
  set(CI_CONFIG default)
  include(${CI_CONFIG})
endif()

list(APPEND SICONOS_CMAKE_OPTIONS -DCMAKE_INSTALL_PREFIX=${SICONOS_INSTALL_DIR})

set(current_project siconos_install)

include(${CTEST_SOURCE_DIRECTORY}/ci_gitlab/ctest_common.cmake)

# ============= Summary =============
message(STATUS "\n============================================ Summary ============================================")
message(STATUS "CTest process for siconos-build has ended.")
message(STATUS "Ctest model is: ${model}")
message(STATUS "Ctest executed on sources directory : ${CTEST_SOURCE_DIRECTORY}")
message(STATUS "CTEST_BINARY_DIRECTORY is: ${CTEST_BINARY_DIRECTORY}")
message(STATUS "CTEST_BUILD_CONFIGURATION is: ${CTEST_BUILD_CONFIGURATION}")
message(STATUS "Build name (cdash) : ${CTEST_BUILD_NAME}")
message(STATUS "Site (cdash) : ${CTEST_SITE}")
message(STATUS "=================================================================================================\n")


