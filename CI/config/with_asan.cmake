# from default, test solvers with sanitizer
include(config/default.cmake)
set_option(USE_SANITIZER asan)
