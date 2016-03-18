# from default, test solvers with mumps
include(CI/config/default.cmake)
set_option(WITH_CXX OFF)
set_option(COMPONENTS "externals;numerics")
