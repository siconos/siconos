## This file should be placed in the root directory of your project.
## Then modify the CMakeLists.txt file in the root directory of your
## project to incorporate the testing dashboard.
## # The following are required to uses Dart and the Cdash dashboard
##   ENABLE_TESTING()
##   INCLUDE(CTest)
set(CTEST_PROJECT_NAME "Siconos")
set(CTEST_NIGHTLY_START_TIME "20:00:00 CET")

set(CTEST_DROP_METHOD "http")
set(CTEST_DROP_SITE "cdash-bipop.inrialpes.fr")
set(CTEST_DROP_LOCATION "/submit.php?project=Siconos")
set(CTEST_DROP_SITE_CDASH TRUE)
if(BUILD_NAME)
  set(CTEST_BUILD_NAME Siconos-${BUILD_NAME})
endif()
  
