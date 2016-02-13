# from default, with sanitizer
# no python with clang + sanitizer
include(CI/config/without_python)
set_option(USE_SANITIZER asan)
