# hack, on opensuse 42.3 cxx1 is needed and there is a compilation failure
# with using std::isnan after including <cmath>
include(default)
set_option(SICONOS_STD_ISNAN_ALREADY_HERE_AND_I_DO_NOT_KNOW_WHY ON)
