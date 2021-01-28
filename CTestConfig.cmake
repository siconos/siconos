## This file should be placed in the root directory of your project.
## Then modify the CMakeLists.txt file in the root directory of your
## project to incorporate the testing dashboard.
##
## # The following are required to submit to the CDash dashboard:
##   ENABLE_TESTING()
##   INCLUDE(CTest)

set(CTEST_PROJECT_NAME "Siconos")
set(CTEST_NIGHTLY_START_TIME "20:00:00 CET")
set(CTEST_DROP_METHOD "http")
# -- Drop site = on Nova VM --
set(CTEST_DROP_SITE "siconos-web.univ-grenoble-alpes.fr:8080")
set(CTEST_DROP_LOCATION "/submit.php?project=siconos-dashboard")
#set(CTEST_DROP_SITE "my.cdash.org")
#set(CTEST_DROP_LOCATION "/submit.php?project=sicosandbox")
set(CTEST_DROP_SITE_CDASH TRUE)
