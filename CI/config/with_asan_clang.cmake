# from default, with sanitizer
# no python with clang + sanitizer
include(without_python)
set_option(USE_SANITIZER asan)
