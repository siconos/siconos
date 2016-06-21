# from default, with sanitizer
# no python with clang + sanitizer
include(config/without_python.cmake)
set_option(USE_SANITIZER asan)
