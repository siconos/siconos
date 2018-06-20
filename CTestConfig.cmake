## This file should be placed in the root directory of your project.
## Then modify the CMakeLists.txt file in the root directory of your
## project to incorporate the testing dashboard.
##
## # The following are required to submit to the CDash dashboard:
##   ENABLE_TESTING()
##   INCLUDE(CTest)

set(CTEST_PROJECT_NAME "Siconos")
set(CTEST_NIGHTLY_START_TIME "20:00:00 CET")

# fedora fails with https  
#    Error message was: libcurl was built with SSL disabled, https: not supported!
#    Problems when submitting via HTTP
# CMake Warning at /home/maurice/src/git/siconos/CI/machinery/CTestDriver.cmake:172 (message):
#    *** submission failure ***
# Note: option --system-curl for cmake build is not sufficient.

set(CTEST_DROP_METHOD "http")
set(CTEST_DROP_SITE "cdash-tripop.inrialpes.fr")
set(CTEST_DROP_LOCATION "/submit.php?project=Siconos")
set(CTEST_DROP_SITE_CDASH TRUE)
if(BUILD_NAME)
  set(CTEST_BUILD_NAME Siconos-${BUILD_NAME})
endif()
